"""Stage 2 RC lattice: nonlinear material mappings, rebar struts, nonlinear pushover (D19)."""

import numpy as np
import pytest

from rclattice.builders import build_lattice_rc, select_nodes
from rclattice.materials import (
    concrete_uniaxial_elastic,
    concrete_uniaxial_nonlinear,
    concrete_uniaxial_regularized,
    steel_uniaxial,
)
from rclattice.model import Load
from rclattice.opensees import run_gravity, run_pushover
from rclattice.problem import (
    BoxLoad,
    BoxSupport,
    ConcreteGrade,
    Problem,
    Rebar,
    SteelGrade,
    portal_frame,
)
from rclattice.reinforcement import rebar_node_chain

CORE = ConcreteGrade("core", E=4030.0, nu=0.2, rho=2.25e-7, fc=6.0, epsc0=0.004, fcu=5.0, epsU=0.014)
COVER = ConcreteGrade("cover", E=4030.0, nu=0.2, rho=2.25e-7, fc=5.0, epsc0=0.002, fcu=0.0, epsU=0.006)
BEAMC = ConcreteGrade("beam", E=4030.0, nu=0.2, rho=2.25e-7)
STEEL = SteelGrade("rebar", fy=60.0, E0=30000.0, b=0.01)


def test_concrete02_mapping_signs_and_defaults():
    mat = concrete_uniaxial_nonlinear(CORE, 7)
    assert mat.mtype == "Concrete02" and mat.id == 7
    fpc, epsc0, fpcu, epsU, lam, ft, ets = mat.args
    assert (fpc, epsc0, fpcu, epsU) == (-6.0, -0.004, -5.0, -0.014)  # compression negative
    assert lam == 0.1
    assert ft == pytest.approx(0.6) and ets == pytest.approx(403.0)  # ~0.1*fc, ~0.1*E defaults


def test_concrete02_requires_nonlinear_params():
    with pytest.raises(ValueError):
        concrete_uniaxial_nonlinear(BEAMC, 1)  # no fc/epsc0/... -> error


def test_steel02_mapping():
    mat = steel_uniaxial(STEEL, 3)
    assert mat.mtype == "Steel02"
    assert mat.args == (60.0, 30000.0, 0.01, 18.0, 0.925, 0.15)


def test_regularized_softening_is_gentler_for_shorter_struts():
    """Crack-band regularization (D20): a shorter strut gets a steeper tension slope and a
    gentler (larger) compression crushing strain than a longer one — both length-dependent."""
    short = concrete_uniaxial_regularized(CORE, 1, 1.5)
    long = concrete_uniaxial_regularized(CORE, 1, 3.0)
    ets_s, ets_l = short.args[6], long.args[6]
    epsU_s, epsU_l = -short.args[3], -long.args[3]
    assert ets_s < ets_l                 # Ets ~ ft^2*L/(2Gf): grows with length
    assert epsU_s > epsU_l               # epsU ~ epsc0 + c/L: shorter strut crushes more gently
    assert epsU_s >= CORE.epsU           # never steeper than the nominal grade


def test_rebar_node_chain_on_grid():
    coords = np.array([(x, y) for x in (0.0, 1.5, 3.0) for y in (0.0, 1.5, 3.0)])
    chain = rebar_node_chain(coords, [(1.5, 0.0), (1.5, 3.0)])  # vertical line x=1.5
    xs = {coords[i][0] for i in chain}
    ys = [coords[i][1] for i in chain]
    assert xs == {1.5}                       # only the x=1.5 column
    assert ys == sorted(ys) and len(ys) == 3  # ordered, all three nodes


def test_rebar_chain_raises_when_unaligned():
    coords = np.array([(0.0, 0.0), (1.5, 0.0)])
    with pytest.raises(ValueError):
        rebar_node_chain(coords, [(0.7, 0.0), (0.7, 3.0)])  # no node on the path


# --- a small RC frame for the structural smoke tests ------------------------
H, SPAN, COL, BEAM, THK = 12.0, 12.0, 6.0, 6.0, 6.0
MESH, EPS = 1.5, 1e-6


def _rc_model():
    half = COL / 2.0
    domain = portal_frame(height=H, span=SPAN, col_depth=COL, beam_depth=BEAM, thickness=THK)
    problem = Problem(
        ndm=2, ndf=2, domain=domain, material=CORE,
        supports=[BoxSupport((-half - 1, SPAN + half + 1, -EPS, EPS), (1, 1))],
        loads=[BoxLoad((-half, half, H - EPS, H + EPS), (0.0, -20.0)),
               BoxLoad((SPAN - half, SPAN + half, H - EPS, H + EPS), (0.0, -20.0))],
    )

    def zone_of(x, y):
        if y >= H - EPS:
            return "beam"
        c = 0.0 if x < SPAN / 2.0 else SPAN
        return "cover" if abs(x - c) >= half - 1.5 - EPS else "core"

    def material_for(zone, length):
        if zone == "beam":
            return concrete_uniaxial_elastic(BEAMC, 0)
        return concrete_uniaxial_regularized(CORE if zone == "core" else COVER, 0, length)

    rebars = tuple(
        Rebar([(c + dx, 0.0), (c + dx, H)], area, STEEL)
        for c in (0.0, SPAN) for dx, area in ((-1.5, 0.9), (0.0, 0.6), (1.5, 0.9))
    )
    return build_lattice_rc(problem, MESH, material_for=material_for, zone_of=zone_of,
                            rebars=rebars, strut_area=2.0)


def test_build_lattice_rc_has_zone_and_steel_materials():
    model, _ = _rc_model()
    mtypes = [m.mtype for m in model.uniaxial_materials]
    assert mtypes.count("Concrete02") >= 2   # core + cover, per strut length (regularized, D20)
    assert "Elastic" in mtypes               # elastic beam
    assert "Steel02" in mtypes               # rebar
    # steel struts are vertical and sit on the rebar lines (x = +-1.5, 0 about each column)
    steel_tag = next(m.id for m in model.uniaxial_materials if m.mtype == "Steel02")
    steel_elems = [e for e in model.elements if e.args[1] == steel_tag]
    assert steel_elems
    for e in steel_elems:
        (x0, _), (x1, _) = (model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords)
        assert abs(x0 - x1) < 1e-9           # vertical bars


def test_nonlinear_pushover_runs_and_yields():
    model, _ = _rc_model()
    assert run_gravity(model)["ok"] == 0
    ctrl = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(model, (-COL, SPAN + COL, -EPS, EPS))
    ids = select_nodes(model, (-COL, SPAN + COL, H - EPS, H + BEAM + EPS))
    lat = [Load(n, (10.0 / len(ids), 0.0)) for n in ids]
    res = run_pushover(model, lateral_loads=lat, control_node=ctrl, control_dof=1,
                       dU=0.1, target=1.0, base_nodes=base)
    assert len(res["disp"]) > 2
    assert max(res["shear"]) > 0.0           # develops lateral resistance
    # softens relative to the initial tangent (nonlinear, not a straight line)
    k0 = (res["shear"][1] - res["shear"][0]) / (res["disp"][1] - res["disp"][0])
    kf = (res["shear"][-1] - res["shear"][-2]) / (res["disp"][-1] - res["disp"][-2])
    assert kf < k0
