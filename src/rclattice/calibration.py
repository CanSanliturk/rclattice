"""Lattice area calibration — two independent methods.

1. `calibrate_lattice` (D16) — a STRUCTURAL match: run static + modal analyses of the whole
   structure and fit strut areas in GROUPS (orthogonal vs diagonal) against a reference model's
   tip deflection and periods. Answers "what areas make THIS structure behave like the reference".

2. `energy_balance_area` (D47) — Aydin's (2017) overlapping-lattice HOMOGENIZATION: equate the
   continuum's stored elastic energy under an affine strain field to the lattice's, and solve for
   one uniform EA. Answers "what area makes the lattice an equivalent elastic CONTINUUM". Needs no
   FE solve at all and no reference model — it is a closed-form sum over the strut list.

Method 1 orchestrates analysis runs, so it (transitively) uses the OpenSees backend; method 2 is
pure geometry and does not.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

from .builders import build_continuum, build_lattice
from .mesh import connect_horizon, mesh_rectangle_nodes
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


# --- Aydin (2017) overlapping-lattice energy balance (D47) -------------------------------------
#
# Aydin, B.B. (2017) "Overlapping lattice modeling for concrete fracture simulations using
# sequentially linear analysis", MSc thesis, METU — Section 2.2, Eqs 2.1-2.3.
#
# Impose an AFFINE displacement field u = F.x on the lattice. Because the field is affine every
# strut's elongation is exact (no FE solve): for a strut vector v = x_j - x_i of length L,
#   elongation dL = (F.v).v / L      axial strain e = dL / L = (F.v).v / L^2
# With EA = 1 the axial force is N = e and the strut strain energy is N^2 L / 2 (Eq 2.2). Summing
# over all struts and equating to the continuum's stored energy for the same field (Eq 2.1) gives
# the one uniform EA that makes the lattice store the right energy (Eq 2.3):
#
#   EA_t = W_continuum / W_lattice(EA=1)
#
# Note EA_t is proportional to E, so the calibrated AREA A_t = EA_t / E is INDEPENDENT of E — one
# area serves every material zone of a model as long as nu, thickness and grid spacing are uniform.


@dataclass
class EnergyBalanceResult:
    """Outcome of Aydin's energy balance, with the isotropy diagnostics it implies.

    `area` is the calibrated uniform strut area (from the normal-strain balance at `nu`). The
    remaining fields expose what the method cannot control: a horizon-strut lattice is an isotropic
    solid only at ONE Poisson ratio (`nu_consistent`), so matching E under normal strain and
    matching G under shear give different areas unless `nu == nu_consistent`.
    """

    area: float             # A_t — the calibrated uniform strut area (from the eps_x balance)
    EA: float               # E * A_t
    nu: float               # Poisson ratio used in the continuum energy (Eq 2.1)
    area_x: float           # A_t from the eps_x balance
    area_y: float           # A_t from the eps_y balance (equals area_x for an isotropic grid)
    area_shear: float       # A_t implied by matching G under pure shear instead
    nu_consistent: float    # the nu at which the normal and shear balances agree
    anisotropy: float       # |area_x - area_y| / area_x — grid directional bias
    isotropy_error: float   # |area_shear - area_x| / area_x at `nu` — cost of the pinned nu
    n_nodes: int
    n_struts: int


def _affine_strut_energy(coords: np.ndarray, pairs, F: np.ndarray) -> float:
    """Total lattice strain energy under the affine field u = F.x with EA = 1 (Aydin Eq 2.2)."""
    c = np.asarray(coords, dtype=float)
    idx = np.asarray(pairs, dtype=int).reshape(-1, 2)   # reshape keeps an empty strut list valid
    v = c[idx[:, 1]] - c[idx[:, 0]]
    L = np.linalg.norm(v, axis=1)
    strain = np.einsum("ij,ij->i", v @ F.T, v) / (L * L)   # (F.v).v / L^2
    return float(np.sum(strain * strain * L) / 2.0)


def energy_balance_area(
    coords: np.ndarray,
    pairs,
    *,
    E: float,
    nu: float,
    thickness: float,
    area_inplane: float,
    strain: float = 1e-3,
) -> EnergyBalanceResult:
    """Aydin's energy-balance calibration: the uniform strut area matching the elastic continuum.

    `coords`/`pairs` are the lattice node coordinates and strut index pairs (as produced by
    `mesh.mesh_rectangle_grid` + `mesh.connect_horizon`); `area_inplane` is the in-plane area of
    the region they cover (Aydin's `A`, i.e. `volume / thickness`). Plane stress.

    Three affine fields are applied: uniaxial `eps_x` (with eps_y = 0), uniaxial `eps_y`, and pure
    shear. The FIRST sets the returned `area`; the other two are diagnostics — `eps_y` measures the
    grid's directional bias, and the shear field reveals the lattice's own Poisson ratio.

    `strain` is arbitrary (the balance is a ratio of two quadratic forms, so it cancels exactly);
    it is exposed only so a caller can confirm that.
    """
    e = float(strain)
    Fx = np.array([[e, 0.0], [0.0, 0.0]])
    Fy = np.array([[0.0, 0.0], [0.0, e]])
    Fs = np.array([[0.0, e / 2.0], [e / 2.0, 0.0]])   # pure shear, engineering gamma = e

    w_x = _affine_strut_energy(coords, pairs, Fx)
    w_y = _affine_strut_energy(coords, pairs, Fy)
    w_s = _affine_strut_energy(coords, pairs, Fs)
    if min(w_x, w_y, w_s) <= 0.0:
        raise ValueError("lattice stores no energy in one or more directions — the strut set is "
                         "degenerate (check mesh_size / horizon)")

    vol = thickness * area_inplane
    # Eq 2.1, plane stress: uniaxial strain with the transverse strain restrained.
    cont_normal = E * e * e * vol / (2.0 * (1.0 - nu * nu))
    # Pure shear: W = G*gamma^2/2 per unit volume, G = E / (2(1+nu)).
    cont_shear = (E / (2.0 * (1.0 + nu))) * e * e * vol / 2.0

    ea_x, ea_y = cont_normal / w_x, cont_normal / w_y
    ea_s = cont_shear / w_s

    # The normal and shear balances agree only at one nu. Setting them equal and cancelling the
    # common (1+nu) factor leaves a closed form:  nu = 1 - 2 * W_shear / W_normal.
    nu_consistent = 1.0 - 2.0 * w_s / w_x

    return EnergyBalanceResult(
        area=ea_x / E,
        EA=ea_x,
        nu=float(nu),
        area_x=ea_x / E,
        area_y=ea_y / E,
        area_shear=ea_s / E,
        nu_consistent=float(nu_consistent),
        anisotropy=abs(ea_x - ea_y) / ea_x,
        isotropy_error=abs(ea_s - ea_x) / ea_x,
        n_nodes=len(coords),
        n_struts=len(pairs),
    )


def energy_balance_rectangle(
    length: float,
    height: float,
    mesh_size: float,
    *,
    E: float,
    nu: float,
    thickness: float,
    horizon: float = 1.5,
    strain: float = 1e-3,
) -> EnergyBalanceResult:
    """`energy_balance_area` on a freshly-meshed rectangular patch of the same grid and horizon.

    The convenience entry point: calibrate against a representative block of the member (Aydin
    applies the balance over the specimen geometry itself). The result is a material-level
    property of the (mesh_size, horizon, nu, thickness) lattice, so it transfers to any model
    built on the same grid — including compound domains and multi-zone models.
    """
    coords = mesh_rectangle_nodes(length, height, mesh_size)
    pairs = connect_horizon(coords, mesh_size, horizon)
    return energy_balance_area(coords, pairs, E=E, nu=nu, thickness=thickness,
                               area_inplane=length * height, strain=strain)
