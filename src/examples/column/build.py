"""Model builders + stiffness calibration for the RC cantilever column (shared by every entry
script).

Two model families, built from the same `specimen`:
  - fiber `forceBeamColumn` reference (`beamcolumn_reference` / `beamcolumn_reference_linear`);
  - RC lattice (`rc_lattice` / `rc_lattice_linear`).
Each comes in a NONLINEAR flavour (Concrete02/Steel02, length-regularized struts) and a LINEAR
flavour (Elastic everything). `calibrate_area` matches the plain elastic lattice's strut area to the
beam-column K0; `calibrate_area_linear` refines that for the full RC-topology elastic lattice (rebar
adds a fixed parallel stiffness, so K0 is not proportional to area) via a secant root-find on
`lattice_k0`.
"""

from __future__ import annotations

from rclattice.builders import build_lattice, build_lattice_rc, select_nodes
from rclattice.materials import (
    concrete_uniaxial_elastic,
    concrete_uniaxial_nonlinear,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
)
from rclattice.opensees import run_beamcolumn_cantilever, run_pushover

from specimen import (
    CORE, COVER_C, DU, EPS, H, HORIZON, MESH, P, STEEL, TARGET, W,
    column_problem, lateral_loads, rebars, zone_of,
)


# --- nonlinear models (Concrete02 / Steel02) --------------------------------

def make_material_for(regularize: bool, Gf: float, Gfc: float):
    """Build the per-strut concrete material factory passed to `build_lattice_rc`.

    `regularize=True` -> crack-band / fracture-energy length-regularized Concrete02 (D20):
    each strut's softening is scaled by its length (`Gf` tension, `Gfc` compression) so the
    dissipated energy is mesh-objective and small struts stay convergent past yield.
    `regularize=False` -> plain Concrete02 straight from the grade (no length regularization);
    `Gf`/`Gfc` are then unused. The grade (confined CORE vs unconfined COVER_C) is still picked
    per `zone` either way.
    """
    def material_for(zone: str, length: float):
        grade = CORE if zone == "core" else COVER_C
        if regularize:
            return concrete_uniaxial_regularized(grade, 0, length, Gf=Gf, Gfc=Gfc)
        return concrete_uniaxial_nonlinear(grade, 0)
    return material_for


def beamcolumn_reference() -> dict:
    """Force-based fiber beam-column, material-matched to the lattice (same grade-level Concrete02
    core/cover + Steel02). The lattice's length regularization is strut-specific (no fiber-section
    analog), so the fiber section uses the plain grades."""
    materials = (concrete_uniaxial_nonlinear(CORE, 1),
                 concrete_uniaxial_nonlinear(COVER_C, 2),
                 steel_uniaxial(STEEL, 3))
    return run_beamcolumn_cantilever(height=H, P=P, dU=DU, target=TARGET, materials=materials)


def calibrate_area(k_bc: float):
    """Strut area giving the elastic lattice the same initial lateral stiffness as `k_bc`.
    Returns (area, control_node, base_nodes)."""
    lat0, _ = build_lattice(column_problem(CORE), MESH, strut_area=1.0, horizon=HORIZON)
    ctrl = select_nodes(lat0, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(lat0, (-W, W, -EPS, EPS))
    r = run_pushover(lat0, lateral_loads=lateral_loads(lat0), control_node=ctrl,
                     control_dof=1, dU=DU, target=DU, base_nodes=base)
    return k_bc / (r["shear"][-1] / r["disp"][-1]), ctrl, base


def rc_lattice(regularize: bool, Gf: float, Gfc: float, area: float):
    """The calibrated RC lattice column model (corotTruss struts + rebar)."""
    model, _ = build_lattice_rc(column_problem(CORE), MESH,
                                material_for=make_material_for(regularize, Gf, Gfc),
                                zone_of=zone_of, rebars=rebars(), strut_area=area,
                                horizon=HORIZON,
                                strut_element="corotTruss")
    return model


# --- linear models (Elastic everything) -------------------------------------

def beamcolumn_reference_linear() -> dict:
    """The fiber `forceBeamColumn` cantilever with the EXACT benchmark section but Elastic fibers:
    concrete E (core + cover) and steel E0. Linear, so its pushover is a straight line of slope K0
    (= the transformed-section lateral stiffness, with gravity P-Delta)."""
    materials = (concrete_uniaxial_elastic(CORE, 1),
                 concrete_uniaxial_elastic(COVER_C, 2),
                 steel_uniaxial_elastic(STEEL, 3))
    return run_beamcolumn_cantilever(height=H, P=P, dU=DU, target=TARGET, materials=materials)


def rc_lattice_linear(area: float):
    """The SAME RC lattice topology as the nonlinear case (concrete struts + longitudinal bars +
    stirrups, corotTruss) but fully linear: Elastic concrete struts (E per zone — both grades share
    E=4030) and Elastic rebar (steel E0)."""
    model, _ = build_lattice_rc(
        column_problem(CORE), MESH,
        material_for=lambda zone, length: concrete_uniaxial_elastic(CORE if zone == "core" else COVER_C, 0),
        zone_of=zone_of, rebars=rebars(), strut_area=area, horizon=HORIZON,
        strut_element="corotTruss", rebar_material=steel_uniaxial_elastic,
    )
    return model


def lattice_k0(area: float, ctrl: int, base: list[int]) -> float:
    """Initial lateral stiffness K0 of the elastic RC lattice at concrete strut area `area`:
    gravity (held) then one tiny DisplacementControl step; K0 = base shear / control disp.
    Linear, so a single small step is exact."""
    model = rc_lattice_linear(area)
    r = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                     control_dof=1, dU=1e-3, target=1e-3, base_nodes=base)
    return r["shear"][1] / r["disp"][1]


def calibrate_area_linear(k_bc: float, *, tol: float = 5e-4, max_iter: int = 15):
    """Concrete strut area giving the FULL elastic RC lattice (concrete + rebar + stirrups) the same
    initial lateral stiffness as `k_bc`.

    Unlike the plain-lattice calibration (where K scales linearly with strut area, so one solve
    suffices), the rebar adds a fixed parallel stiffness, so K0(area) is monotone but not
    proportional — a short secant root-find on K0(area) - k_bc is used. The plain-lattice area
    (from `calibrate_area`, rebar ignored) seeds it and supplies the control/base node ids.
    Returns (area, control_node, base_nodes)."""
    a_plain, ctrl, base = calibrate_area(k_bc)   # plain-elastic-lattice seed + node ids
    a0, a1 = a_plain, 0.85 * a_plain             # rebar adds stiffness -> true area < plain area
    f0 = lattice_k0(a0, ctrl, base) - k_bc
    f1 = lattice_k0(a1, ctrl, base) - k_bc
    for _ in range(max_iter):
        if abs(f1) <= tol * k_bc or f1 == f0:
            break
        a2 = max(a1 - f1 * (a1 - a0) / (f1 - f0), 1e-4)   # secant step, keep area positive
        a0, f0, a1 = a1, f1, a2
        f1 = lattice_k0(a1, ctrl, base) - k_bc
    return a1, ctrl, base
