"""Modal calibration checks (D16): mass/eigen sanity, physical bounds, improvement."""

import numpy as np

from rclattice import ConcreteGrade, EdgeLoad, EdgeSupport, Problem, RectangleDomain
from rclattice.builders import build_continuum
from rclattice.calibration import (
    calibrate_lattice,
    combined_rms,
    continuum_targets,
    nominal_area,
    orientation_area_fn,
)
from rclattice.opensees import run_modal

MESH = 0.1


def _problem() -> Problem:
    return Problem(
        ndm=2,
        ndf=2,
        domain=RectangleDomain(length=1.0, height=0.2, thickness=0.1),
        material=ConcreteGrade(name="C30", E=30e9, nu=0.2, rho=2400.0),
        supports=[EdgeSupport(edge="xmin", fix=(1, 1))],
        loads=[EdgeLoad(edge="xmax", total=(0.0, -10e3))],
    )


def test_continuum_modal_periods_positive_and_sorted():
    cont, _ = build_continuum(_problem(), MESH)
    periods = run_modal(cont, 4)["periods"]
    assert all(np.isfinite(periods)) and all(p > 0 for p in periods)
    assert periods == sorted(periods, reverse=True)  # ascending eigenvalue -> descending period


def test_calibration_stays_within_physical_bounds_and_improves():
    problem = _problem()
    targets = continuum_targets(problem, MESH, n_modes=3)
    nom = nominal_area(problem, MESH)
    lo, hi = 1e-3 * nom, 3.0 * nom

    result = calibrate_lattice(problem, MESH, targets=targets, n_modes=3)
    assert result.success
    for area in result.areas.values():
        assert lo <= area <= hi + 1e-12  # physical bounds respected

    # area-group fit should not be worse than the static-only scalar under the same objective
    scalar = orientation_area_fn(0.34 * nom, 0.34 * nom, MESH)
    assert result.rms <= combined_rms(problem, MESH, scalar, targets, n_modes=3) + 1e-9
