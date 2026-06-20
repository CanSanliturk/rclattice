"""Lattice area calibration (D16): match static deflection + modal periods.

Calibrates lattice strut areas so the lattice reproduces a target static response and the
first N modal periods. Areas are parameterized in GROUPS (orthogonal vs diagonal struts on
the regular grid); a single uniform area is the degenerate 1-group case. The target is either
the continuum reference model or user-supplied periods.

This module orchestrates analysis runs, so it (transitively) uses the OpenSees backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

from .builders import build_continuum, build_lattice
from .opensees import run_modal, run_static
from .problem import Problem

# Default residual weights (w_static, w_T1, w_higher_modes): static + fundamental are firm
# targets (both flexural), higher modes are soft so the fit doesn't chase what an axial lattice
# can't represent (D16).
DEFAULT_WEIGHTS = (1.0, 1.0, 0.3)


@dataclass
class CalibrationTargets:
    """What the lattice is matched to: a static response scalar + a list of periods."""

    static: float
    periods: list[float]


@dataclass
class CalibrationResult:
    areas: dict[str, float]          # {"orthogonal": A, "diagonal": A}
    area_fn: Callable[[float], float]
    success: bool
    residuals: list[float]
    rms: float


def _load_dof(problem: Problem) -> int:
    """The DOF index of the dominant load direction (selector-agnostic)."""
    return int(np.argmax(np.abs(problem.loads[0].total)))


def static_response(model, dof: int) -> float:
    """Mean displacement of the loaded nodes in the load direction (the static metric).

    Reads the loaded node ids from `model.loads`, so it works for any support/load selector
    (edge or box) — the builder has already resolved which nodes carry the load.
    """
    res = run_static(model)
    if res["ok"] != 0:
        raise RuntimeError(f"static analysis failed (rc={res['ok']})")
    loaded = [ld.node for ld in model.loads]
    return float(np.mean([res["disps"][nid][dof] for nid in loaded]))


def orientation_area_fn(a_orthogonal: float, a_diagonal: float, mesh_size: float) -> Callable[[float], float]:
    """Area as a function of strut length: short (orthogonal ~ s) vs long (diagonal ~ s*sqrt2)."""
    threshold = 1.2 * mesh_size  # between s and s*sqrt(2) ~= 1.414 s
    return lambda length: a_orthogonal if length <= threshold else a_diagonal


def continuum_targets(
    problem: Problem,
    mesh_size: float,
    *,
    n_modes: int = 3,
    plane: str = "PlaneStress",
) -> CalibrationTargets:
    """Compute calibration targets from the continuum reference model."""
    dof = _load_dof(problem)
    model, _edges = build_continuum(problem, mesh_size, plane=plane)
    static = static_response(model, dof)
    periods = run_modal(model, n_modes)["periods"]
    return CalibrationTargets(static=static, periods=list(periods)[:n_modes])


def _lattice_response(
    problem: Problem,
    mesh_size: float,
    area_fn: Callable[[float], float],
    n_modes: int,
    horizon: float,
    dof: int,
) -> tuple[float, list[float]]:
    model, _edges = build_lattice(problem, mesh_size, horizon=horizon, strut_area=area_fn)
    static = static_response(model, dof)
    periods = run_modal(model, n_modes)["periods"]
    return static, list(periods)


def _residuals(
    static: float,
    periods: list[float],
    targets: CalibrationTargets,
    weights: tuple[float, float, float],
) -> list[float]:
    """Normalized residuals: weights = (w_static, w_T1, w_higher). Mode 0 (T1) uses w_T1,
    all higher modes use w_higher."""
    w_static, w_t1, w_higher = weights
    n = min(len(periods), len(targets.periods))
    r = [w_static * (static - targets.static) / targets.static]
    for k in range(n):
        w = w_t1 if k == 0 else w_higher
        r.append(w * (periods[k] - targets.periods[k]) / targets.periods[k])
    return r


def nominal_area(problem: Problem, mesh_size: float) -> float:
    """Physically-motivated reference strut area: thickness x grid spacing (D16)."""
    return problem.domain.thickness * mesh_size


def calibrate_lattice(
    problem: Problem,
    mesh_size: float,
    *,
    targets: CalibrationTargets,
    n_modes: int = 3,
    horizon: float = 1.5,
    weights: tuple[float, float, float] = DEFAULT_WEIGHTS,
    area_bounds: tuple[float, float] | None = None,
) -> CalibrationResult:
    """Fit (orthogonal, diagonal) strut areas to the targets via bounded nonlinear least squares.

    `area_bounds` clamps both group areas to a physical range; default is
    (1e-3, 3.0) x nominal_area so the fit stays physical instead of chasing higher modes with
    a near-rigid diagonal (D16). Higher-mode residual error is expected and documented.
    """
    dof = _load_dof(problem)
    nom = nominal_area(problem, mesh_size)
    lo, hi = area_bounds if area_bounds is not None else (1e-3 * nom, 3.0 * nom)
    x0 = np.clip([0.34 * nom, 0.34 * nom], lo, hi)  # start near the static-calibrated scalar

    def fun(params: np.ndarray) -> list[float]:
        area_fn = orientation_area_fn(params[0], params[1], mesh_size)
        static, periods = _lattice_response(problem, mesh_size, area_fn, n_modes, horizon, dof)
        return _residuals(static, periods, targets, weights)

    sol = least_squares(fun, x0, bounds=([lo, lo], [hi, hi]))
    residuals = list(sol.fun)
    return CalibrationResult(
        areas={"orthogonal": float(sol.x[0]), "diagonal": float(sol.x[1])},
        area_fn=orientation_area_fn(sol.x[0], sol.x[1], mesh_size),
        success=bool(sol.success),
        residuals=residuals,
        rms=float(np.sqrt(np.mean(np.square(residuals)))),
    )


def combined_rms(
    problem: Problem,
    mesh_size: float,
    area_fn: Callable[[float], float],
    targets: CalibrationTargets,
    *,
    n_modes: int = 3,
    horizon: float = 1.5,
    weights: tuple[float, float, float] = DEFAULT_WEIGHTS,
) -> float:
    """Diagnostic: weighted RMS of the static + period residuals for a given area_fn."""
    dof = _load_dof(problem)
    static, periods = _lattice_response(problem, mesh_size, area_fn, n_modes, horizon, dof)
    r = _residuals(static, periods, targets, weights)
    return float(np.sqrt(np.mean(np.square(r))))
