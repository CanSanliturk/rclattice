"""Optional bond elements (D72) — Aydin, Tuncay & Binici (2019), J. Struct. Eng. 145(9): 04019091.

Perfect bond stays the default; passing `bond_material` to `build_lattice_rc` switches the
reinforcement onto its own nodes, tied back to the concrete by a ring of bond struts. Most of these are pure assembly (no OpenSees solve): they check the topology the paper describes,
the mass bookkeeping the extra nodes force, and the shape of the bond law.

THE ASSEMBLY TESTS ALL PASSED WHILE THE FEATURE WAS BROKEN, which is the lesson recorded here.
Counts, lengths, masses and law shape were all correct and the resulting model was still 2.19x
STIFFER than perfect bond — physically impossible, since bond can only add flexibility. What was
missing was a test of what the assembly DOES rather than what it contains, so
`test_rigid_bond_reduces_to_perfect_bond` is that test; it needs a solve and it currently xfails.
"""


import numpy as np
import pytest

from rclattice import (
    BoxLoad,
    BoxSupport,
    ConcreteGrade,
    Problem,
    Rebar,
    RectangleDomain,
    SteelGrade,
    bond_elastic_brittle,
    build_lattice_rc,
)
from rclattice.materials import concrete_uniaxial_regularized

MESH = 100.0
CONC = ConcreteGrade(name="c", E=25000.0, nu=0.2, rho=2.4e-9, fc=28.0, ft=1.85,
                     epsc0=2 * 28.0 / 25000.0, fcu=5.6, epsU=0.01)
STEEL = SteelGrade(name="s", fy=360.0, E0=200000.0)


def _build(**extra):
    if "bond_material" in extra:
        extra.setdefault("i_accept_the_known_bond_defect", True)
    prob = Problem(ndm=2, ndf=2,
                   domain=RectangleDomain(length=400.0, height=400.0, thickness=120.0),
                   material=CONC,
                   supports=(BoxSupport(box=(-1, 401, -1, 1), fix=(1, 1)),),
                   loads=(BoxLoad(box=(-1, 401, 399, 401), total=(1000.0, 0.0)),))
    bars = (Rebar(path=[(200.0, 0.0), (200.0, 400.0)], area=100.0, steel=STEEL),
            Rebar(path=[(0.0, 200.0), (400.0, 200.0)], area=100.0, steel=STEEL, role="stirrup"))
    model, _ = build_lattice_rc(
        prob, MESH, rebars=bars, strut_area=1000.0,
        material_for=lambda z, L: concrete_uniaxial_regularized(CONC, 0, L, Gf=0.075),
        zone_of=lambda x, y: "w", **extra)
    return model


def _kinds(model):
    out: dict[str, int] = {}
    for e in model.elements:
        out[e.kind] = out.get(e.kind, 0) + 1
    return out


def _length(model, e):
    a = np.array(model.nodes[e.nodes[0]].coords)
    b = np.array(model.nodes[e.nodes[1]].coords)
    return float(np.linalg.norm(b - a))


def test_default_is_still_perfect_bond():
    """No `bond_material` -> no extra nodes, no bond elements: D13 behaviour is untouched."""
    m = _build()
    assert len(m.nodes) == 25            # the 5x5 grid alone
    assert "bond" not in _kinds(m)


def test_steel_gets_its_own_nodes():
    """Two 5-node bars -> 10 steel nodes on top of the grid, so steel and concrete are independent."""
    m = _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0))
    assert len(m.nodes) == 25 + 10
    # every rebar strut now runs between steel nodes, never on a concrete one
    steel_ids = set(range(26, 36))
    for e in m.elements:
        if e.kind in ("longitudinal", "stirrup"):
            assert set(e.nodes) <= steel_ids


def test_bond_ring_matches_the_horizon_and_skips_the_coincident_node():
    """Each steel node links to the concrete nodes within the horizon EXCEPT the one it sits on.

    On this grid that is 8 links for an interior bar node and 5 for one on an edge: 34 per bar.
    A link to the coincident node would be a zero-length truss, so the shortest is one grid step.
    """
    m = _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0))
    bonds = [e for e in m.elements if e.kind == "bond"]
    assert len(bonds) == 68
    lengths = sorted({round(_length(m, e), 6) for e in bonds})
    assert lengths == [MESH, pytest.approx(MESH * np.sqrt(2))]


def test_bond_horizon_is_separable_from_the_strut_horizon():
    """The paper warns bond stiffens and strengthens as its horizon grows, so it is its own knob."""
    near = _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0), bond_horizon=1.0)
    far = _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0), bond_horizon=1.5)
    assert _kinds(near)["bond"] < _kinds(far)["bond"]
    assert _kinds(near)["concrete"] == _kinds(far)["concrete"]   # struts unaffected


def test_mass_is_conserved_when_steel_nodes_are_added():
    """Steel nodes have no tributary area; they take a share of their host's, they do not add mass."""
    plain = sum(v[0] for v in _build().masses.values())
    bonded = _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0))
    assert sum(v[0] for v in bonded.masses.values()) == pytest.approx(plain)
    assert all(v[0] > 0.0 for v in bonded.masses.values())       # no zero-mass DOF


def test_mass_share_that_would_starve_a_concrete_node_is_rejected():
    """Bars cross at (200, 200): two steel nodes share one host, so 0.5 each takes all of it."""
    with pytest.raises(ValueError, match="bond_mass_share"):
        _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0), bond_mass_share=0.5)


def test_bond_law_is_brittle_with_a_residual_plateau():
    """Fig. 1(c): up to ft at eps_cr, a near-vertical drop to a*ft, then flat — symmetric."""
    mat = bond_elastic_brittle(CONC, 7, residual=0.7)
    assert mat.mtype == "ElasticMultiLinear"
    args = list(mat.args)
    strains = np.array(args[args.index("-strain") + 1:args.index("-stress")], dtype=float)
    stresses = np.array(args[args.index("-stress") + 1:], dtype=float)
    assert np.all(np.diff(strains) > 0)                       # ElasticMultiLinear requires this
    assert stresses.max() == pytest.approx(CONC.ft)
    assert stresses[-1] == pytest.approx(0.7 * CONC.ft)       # residual plateau, not zero
    assert strains[np.argmax(stresses)] == pytest.approx(CONC.ft / CONC.E)
    assert np.allclose(stresses, -stresses[::-1])             # symmetric in slip direction


def test_zero_residual_is_rejected():
    """A fully brittle ring would leave the steel node free once it cracks."""
    with pytest.raises(ValueError, match="residual"):
        bond_elastic_brittle(CONC, 1, residual=0.0)


def test_bond_is_refused_without_the_explicit_opt_in():
    """The known defect must not be reachable by accident (D72)."""
    with pytest.raises(ValueError, match="KNOWN DEFECT"):
        _build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0),
               i_accept_the_known_bond_defect=False)


def _tip_stiffness(model, target: float = 0.02) -> float:
    """Tip stiffness of a built model under a small elastic push."""
    from rclattice.builders import select_nodes
    from rclattice.model import Load
    from rclattice.opensees import run_pushover

    top = select_nodes(model, (-1, 401, 399, 401))
    ctrl = select_nodes(model, (-1, 1, 399, 401))[0]
    base = select_nodes(model, (-1, 401, -1, 1))
    res = run_pushover(model, lateral_loads=[Load(n, (1.0 / len(top), 0.0)) for n in top],
                       control_node=ctrl, control_dof=1, dU=target / 20.0, target=target,
                       base_nodes=base)
    assert res["converged"]
    return res["shear"][-1] / res["disp"][-1]


@pytest.mark.xfail(strict=True, reason="D72: the horizon-star bond topology cannot satisfy this "
                                       "limit at any stiffness; it needs the coincident-node "
                                       "zeroLength interface spring")
def test_rigid_bond_reduces_to_perfect_bond():
    """THE test this feature is defined by, and the one whose absence let the defect ship.

    Bond introduces relative slip between steel and concrete, so it can only ever make a model
    SOFTER than shared nodes. Stiffen the bond and the bonded model must approach the perfect-bond
    one FROM BELOW. The horizon star instead diverges upward without limit, because its links are a
    parallel load path through the concrete rather than an interface.
    """
    perfect = _tip_stiffness(_build())
    stiff = _tip_stiffness(_build(bond_material=lambda z, L: bond_elastic_brittle(CONC, 0),
                                  bond_area=1e3 * CONC.ft))
    assert stiff <= perfect * 1.001, (
        f"bond made the model STIFFER: {stiff / perfect:.3f}x perfect bond")
