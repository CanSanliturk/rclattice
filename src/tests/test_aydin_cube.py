"""Aydin, Binici & Tuncay (2021) tension-only lattice: perturbation, softening algebra, material."""

from __future__ import annotations

import numpy as np
import pytest

from rclattice.materials import aydin_lattice_softening, concrete_lattice_aydin
from rclattice.mesh import connect_horizon, mesh_rectangle_nodes, perturb_nodes
from rclattice.problem import ConcreteGrade

FC = 20.0
GRADE = ConcreteGrade("aydin", E=4700.0 * FC**0.5, nu=0.2, rho=2.4e-9,
                      fc=FC, epsc0=2.0 * FC / (4700.0 * FC**0.5), ft=0.35 * FC**0.5)


def test_perturbation_respects_the_boundary():
    """Boundary nodes stay on their face; corners do not move at all."""
    c0 = mesh_rectangle_nodes(100.0, 100.0, 10.0)
    c = perturb_nodes(c0, 10.0, 0.08, seed=3)

    for axis, value in ((1, 0.0), (1, 100.0), (0, 0.0), (0, 100.0)):
        on = np.isclose(c0[:, axis], value)
        assert np.allclose(c[on, axis], value), "an edge node left its face"

    corner = np.isclose(c0[:, 0] % 100.0, 0.0) & np.isclose(c0[:, 1] % 100.0, 0.0)
    assert np.allclose(c[corner], c0[corner]), "a corner moved"

    moved = np.linalg.norm(c - c0, axis=1)
    assert moved.max() <= 0.08 * 10.0 + 1e-12
    assert moved.max() > 0.5 * 0.08 * 10.0, "perturbation is suspiciously small"
    assert np.allclose(perturb_nodes(c0, 10.0, 0.0), c0), "Rmax/d = 0 must be the plain grid"


def test_perturbation_changes_horizon_connectivity():
    """Why the example freezes topology on the unperturbed grid (D60)."""
    c0 = mesh_rectangle_nodes(100.0, 100.0, 10.0)
    c = perturb_nodes(c0, 10.0, 0.08, seed=7)
    assert len(connect_horizon(c, 10.0, 1.5)) != len(connect_horizon(c0, 10.0, 1.5))


def test_softening_tail_dissipates_the_target_fracture_energy():
    """`a2`/`a3` are solved so L * area-under-the-tail == Gf, for any strut length."""
    ft, E, Gf = GRADE.ft, GRADE.E, 0.05
    a1, b1, b2 = 1.5, 0.6, 0.2
    for length in (5.0, 10.0, 14.142, 20.0):
        a2, a3 = aydin_lattice_softening(ft, E, length, Gf=Gf, a1=a1, b1=b1, b2=b2)
        eps_cr = ft / E
        eps = np.array([1.0, a1, a2, a3]) * eps_cr
        sig = np.array([1.0, b1, b2, 0.0]) * ft
        assert length * float(np.trapezoid(sig, eps)) == pytest.approx(Gf, rel=1e-9)
        assert a3 == pytest.approx(5.0 * a2)


def test_softening_refuses_a_grid_too_coarse_for_its_fracture_energy():
    with pytest.raises(ValueError, match="too long for Gf"):
        aydin_lattice_softening(GRADE.ft, GRADE.E, 5000.0, Gf=0.05)


def test_material_has_no_compressive_strength():
    """The paper's whole idea: compression never peaks, with or without the RSM knee."""
    for rsm in (True, False):
        mat = concrete_lattice_aydin(GRADE, 1, 10.0, Gf=0.05, rsm=rsm)
        assert mat.mtype == "ElasticMultiLinear"
        args = list(mat.args)
        strains = np.array(args[1:args.index("-stress")], dtype=float)
        stresses = np.array(args[args.index("-stress") + 1:], dtype=float)

        neg = strains < 0.0
        assert np.all(np.diff(stresses[neg][::-1]) < 0.0), "compression must rise without limit"
        assert abs(stresses[neg]).max() > 20.0 * FC, "compression is capped somewhere"

        # tension: peaks at ft, returns to zero, and stays there
        assert stresses.max() == pytest.approx(GRADE.ft)
        assert stresses[-1] == 0.0 and stresses[-2] == 0.0
        assert np.all(np.diff(strains) > 0.0), "ElasticMultiLinear needs increasing strains"


def test_rsm_knee_softens_compression_at_a_third_of_epsc0():
    mat = concrete_lattice_aydin(GRADE, 1, 10.0, Gf=0.05, rsm=True, alpha=1 / 3, beta=0.4)
    args = list(mat.args)
    strains = np.array(args[1:args.index("-stress")], dtype=float)
    stresses = np.array(args[args.index("-stress") + 1:], dtype=float)
    knee = -GRADE.epsc0 / 3.0
    i = int(np.argmin(np.abs(strains - knee)))
    assert strains[i] == pytest.approx(knee)
    assert stresses[i] / strains[i] == pytest.approx(GRADE.E, rel=1e-9)      # slope E below
    post = (stresses[i] - stresses[i - 1]) / (strains[i] - strains[i - 1])
    assert post == pytest.approx(0.4 * GRADE.E, rel=1e-9)                    # slope 0.4E beyond


def test_divergence_filter_keeps_a_clean_run_and_catches_a_spike():
    """A quasi-static record can fall fast but cannot RISE fast (D61).

    The filter exists because a diverged realization of the 200 mm cube recorded a 139 MPa spike on
    a lattice whose struts have no compressive strength, reported `converged = True`, and poisoned
    the mean it was averaged into.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "examples" / "compression_cube"))
    from aydin_approach import MAX_RISE_PER_STEP, first_divergence

    fc = 15.1
    eps = np.linspace(0.0, 5e-3, 40_000)
    rise = np.minimum(eps / 1.9e-3, 1.0)
    decay = np.where(eps < 1.9e-3, 1.0, np.exp(-(eps - 1.9e-3) / 5e-4))
    clean = 12.0 * rise * decay

    assert first_divergence(clean, fc) == len(clean), "a clean run must not be truncated"
    assert float(np.diff(clean).max()) < 0.01 * MAX_RISE_PER_STEP * fc, \
        "the genuine per-step rise must sit far below the threshold"

    spiked = clean.copy()
    spiked[25_000] = 139.0
    assert first_divergence(spiked, fc) == 25_000

    # A fast DROP is an instability, which is the mechanism being modelled — never truncate it.
    dropped = clean.copy()
    dropped[30_000:] = 0.0
    assert first_divergence(dropped, fc) == len(dropped)
