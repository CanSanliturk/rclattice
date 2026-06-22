"""D29: nonlinear 2D continuum reference — ASDConcrete3D+PlaneStress mapping + build_continuum_rc.

Structural (no-solve) tests, mirroring test_rc.py's assembly checks; the live convergence of the
continuum pushover is covered by the example run, not the unit suite (it is heavy)."""

import pytest

from rclattice.builders import build_continuum_rc
from rclattice.materials import concrete_nd_nonlinear, steel_uniaxial_elastic
from rclattice.problem import ConcreteGrade, Problem, Rebar, RectangleDomain, SteelGrade

CORE = ConcreteGrade("core", E=4030.0, nu=0.2, rho=2.25e-7, fc=6.0, epsc0=0.004, fcu=4.0, epsU=0.030, ft=0.6)
COVER = ConcreteGrade("cover", E=4030.0, nu=0.2, rho=2.25e-7, fc=5.0, epsc0=0.002, fcu=0.0, epsU=0.006, ft=0.5)
BEAMC = ConcreteGrade("beam", E=4030.0, nu=0.2, rho=2.25e-7)  # elastic-only (no nonlinear params)
STEEL = SteelGrade("rebar", fy=60.0, E0=30000.0, b=0.01)


def _flagged(args, flag):
    """The run of numeric values following `flag` in an ASDConcrete3D arg tuple (until the next flag)."""
    out = []
    for v in args[args.index(flag) + 1:]:
        if isinstance(v, str):
            break
        out.append(v)
    return out


def test_concrete_nd_nonlinear_pair_and_peaks():
    base, wrapper = concrete_nd_nonlinear(CORE, 3, 4, lch=1.5)
    assert base.mtype == "ASDConcrete3D" and base.id == 3
    assert wrapper.mtype == "PlaneStress" and wrapper.id == 4 and wrapper.args == (3,)  # wraps the 3D base
    assert base.args[0] == CORE.E and base.args[1] == CORE.nu
    assert max(_flagged(base.args, "-Ts")) == pytest.approx(CORE.ft)   # tensile peak = ft
    assert max(_flagged(base.args, "-Cs")) == pytest.approx(CORE.fc)   # compressive peak = fc
    for flag in ("-Te", "-Ts", "-Td", "-Ce", "-Cs", "-Cd"):
        assert len(_flagged(base.args, flag)) >= 3                     # full curve present


def test_nd_nonlinear_crack_band_is_length_regularized():
    """Shorter element (lch) → larger ultimate tensile/crushing strains (crack-band; D20 analog).

    Uses the column's fracture energies so neither branch clamps to the grade's nominal epsU."""
    short = concrete_nd_nonlinear(CORE, 1, 2, lch=1.5, Gf=4.0e-3, Gfc=1.5)[0]
    long = concrete_nd_nonlinear(CORE, 1, 2, lch=3.0, Gf=4.0e-3, Gfc=1.5)[0]
    assert _flagged(short.args, "-Te")[-1] > _flagged(long.args, "-Te")[-1]   # eps_tu ~ 2Gf/(ft*lch)
    assert _flagged(short.args, "-Ce")[-1] > _flagged(long.args, "-Ce")[-1]   # eps_cu ~ 2Gfc/((..)*lch)


def test_concrete_nd_nonlinear_requires_nonlinear_params():
    with pytest.raises(ValueError):
        concrete_nd_nonlinear(BEAMC, 1, 2, lch=1.5)  # no fc/epsc0/... -> error


# --- a small RC continuum column for the assembly smoke test ----------------
W, H, THK = 6.0, 6.0, 3.0
MESH, EPS, COVER_T = 1.5, 1e-6, 1.5


def _continuum_model():
    half = W / 2.0
    domain = RectangleDomain(length=W, height=H, thickness=THK, origin=(-half, 0.0))
    problem = Problem(ndm=2, ndf=2, domain=domain, material=CORE, supports=[], loads=[])

    def zone_of(x, _y):
        return "cover" if abs(x) >= half - COVER_T - EPS else "core"

    def nd_material_for(zone):
        return concrete_nd_nonlinear(CORE if zone == "core" else COVER, 0, 0, lch=MESH)

    rebars = (Rebar([(0.0, 0.0), (0.0, H)], 1.2, STEEL),)  # one central longitudinal bar
    return build_continuum_rc(problem, MESH, nd_material_for=nd_material_for, zone_of=zone_of,
                              rebars=rebars, rebar_material=steel_uniaxial_elastic)


def test_build_continuum_rc_assembly():
    model, _ = _continuum_model()
    quads = [e for e in model.elements if e.etype == "quad"]
    steel = [e for e in model.elements if e.etype == "Truss"]
    assert quads and steel
    # two zones (core, cover) -> two nD-material pairs (base + PlaneStress wrapper) = 4 nd materials
    assert [m.mtype for m in model.nd_materials].count("ASDConcrete3D") == 2
    assert [m.mtype for m in model.nd_materials].count("PlaneStress") == 2
    # each quad references a PlaneStress wrapper tag (the last element arg)
    wrappers = {m.id for m in model.nd_materials if m.mtype == "PlaneStress"}
    assert all(q.args[-1] in wrappers for q in quads)
    # rebar struts are vertical (on the x=0 line) and use the elastic steel material
    assert "Elastic" in [m.mtype for m in model.uniaxial_materials]
    for e in steel:
        (x0, _), (x1, _) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        assert abs(x0 - x1) < 1e-9 and abs(x0) < 1e-9   # vertical, on x=0


def test_continuum_rc_nd_and_uniaxial_tags_independent():
    """nDMaterial and uniaxialMaterial tag namespaces are separate, so both may start at 1."""
    model, _ = _continuum_model()
    assert min(m.id for m in model.nd_materials) == 1
    assert min(m.id for m in model.uniaxial_materials) == 1
