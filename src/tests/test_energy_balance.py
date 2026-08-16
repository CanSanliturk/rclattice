"""Aydin (2017) overlapping-lattice energy-balance calibration (D47).

The balance has an exact closed form for the interior of a regular horizon=1.5 grid, which is what
most of these tests check against. Per unit cell (spacing d, one x-orthogonal + one y-orthogonal +
two diagonal struts) under a uniform eps_x with eps_y restrained:

    lattice energy   = eps^2 * d * (1/2 + sqrt(2)/4)        [y-orthogonal struts carry nothing]
    continuum energy = E * eps^2 * w * d^2 / (2(1 - nu^2))  [Aydin Eq. 2.1, plane stress]

so  A_t / (w * d) = [1 / (2(1 - nu^2))] / (1/2 + sqrt(2)/4).

A finite patch sits slightly below that (its free edges carry proportionally more x-orthogonal
struts), approaching it as the patch grows.
"""

import math

import numpy as np
import pytest

from rclattice.calibration import energy_balance_area, energy_balance_rectangle
from rclattice.mesh import connect_horizon, mesh_rectangle_nodes

E, NU, THK, MESH = 30000.0, 1.0 / 3.0, 200.0, 50.0
_CELL = 0.5 + math.sqrt(2) / 4.0


def _closed_form(nu: float) -> float:
    """Interior-lattice A_t / (thickness * mesh) for horizon = 1.5."""
    return (1.0 / (2.0 * (1.0 - nu * nu))) / _CELL


def _balance(length, height, mesh=MESH, *, nu=NU, e=30000.0, horizon=1.5):
    return energy_balance_rectangle(length, height, mesh, E=e, nu=nu, thickness=THK,
                                    horizon=horizon)


@pytest.mark.parametrize("nu", [1.0 / 3.0, 0.2, 0.0])
def test_matches_closed_form_and_converges_to_it(nu):
    """A big patch lands within 1% of the interior closed form, and refines toward it."""
    small = _balance(1000.0, 1000.0, nu=nu).area / (THK * MESH)
    large = _balance(3000.0, 3000.0, nu=nu).area / (THK * MESH)
    exact = _closed_form(nu)
    assert small < large < exact                      # finite patch is soft, monotone toward exact
    assert abs(large - exact) / exact < 0.01


def test_area_is_independent_of_modulus():
    """EA_t is proportional to E, so the calibrated AREA is not — one area serves every zone."""
    a1 = _balance(1000.0, 3000.0, e=16100.0).area
    a2 = _balance(1000.0, 3000.0, e=200000.0).area
    assert a1 == pytest.approx(a2, rel=1e-12)


def test_area_is_independent_of_the_probing_strain():
    """The balance is a ratio of two quadratic forms, so the probe strain cancels exactly."""
    coords = mesh_rectangle_nodes(1000.0, 1000.0, MESH)
    pairs = connect_horizon(coords, MESH, 1.5)
    kw = dict(E=E, nu=NU, thickness=THK, area_inplane=1000.0 * 1000.0)
    a1 = energy_balance_area(coords, pairs, strain=1e-6, **kw).area
    a2 = energy_balance_area(coords, pairs, strain=12.5, **kw).area
    assert a1 == pytest.approx(a2, rel=1e-12)


def test_calibration_constant_is_mesh_objective():
    """A_t / (thickness * mesh) is a property of the lattice, not of the discretisation.

    It is not perfectly constant: a coarser patch has proportionally more free edge, so the
    residual mesh dependence is the closing of that boundary-layer deficit. Across a 4x refinement
    it moves ~1.5% and always TOWARD the interior closed form, which is the objectivity that
    matters — the calibration is not chasing the discretisation.
    """
    ratios = [_balance(1000.0, 3000.0, m).area / (THK * m) for m in (100.0, 50.0, 25.0)]
    assert max(ratios) - min(ratios) < 0.02 * max(ratios)     # within 2% across a 4x refinement
    assert ratios == sorted(ratios)                            # monotone toward the interior value
    assert max(ratios) < _closed_form(NU)                      # and never overshoots it


def test_lattice_has_its_own_poisson_ratio_independent_of_the_one_assumed():
    """The horizon=1.5 grid is isotropic only near nu ~ 0.18, whatever nu the balance is given."""
    got = [_balance(1000.0, 3000.0, nu=nu).nu_consistent for nu in (0.0, 0.2, 1.0 / 3.0)]
    assert got == pytest.approx([got[0]] * 3, rel=1e-12)
    assert 0.15 < got[0] < 0.20
    # ... and the penalty for assuming a different nu is smallest at that value.
    assert (_balance(1000.0, 3000.0, nu=got[0]).isotropy_error
            < _balance(1000.0, 3000.0, nu=1.0 / 3.0).isotropy_error)


def test_wider_horizon_pushes_the_lattice_toward_a_higher_poisson_ratio():
    """More strut directions change the homogenised solid — h=3.01 sits far above h=1.5."""
    assert _balance(1000.0, 1000.0, horizon=3.01).nu_consistent > \
           _balance(1000.0, 1000.0, horizon=1.5).nu_consistent + 0.1


def test_regular_grid_is_directionally_unbiased():
    """A square patch stores the same energy under eps_x and eps_y."""
    res = _balance(1000.0, 1000.0)
    assert res.area_x == pytest.approx(res.area_y, rel=1e-12)
    assert res.anisotropy == pytest.approx(0.0, abs=1e-12)


def test_degenerate_strut_set_is_rejected():
    """A lattice with no struts stores no energy — that must raise, not divide by zero."""
    coords = mesh_rectangle_nodes(1000.0, 1000.0, MESH)
    with pytest.raises(ValueError, match="degenerate"):
        energy_balance_area(np.asarray(coords), [], E=E, nu=NU, thickness=THK, area_inplane=1e6)
