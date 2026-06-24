"""Model builders + stiffness calibration for the RC portal frame (shared by every entry script).

Two model families, built from the same `specimen` (the portal frame = two cantilever columns + a
thinner beam):
  - fiber `forceBeamColumn` FRAME reference (`beamcolumn_reference` / `beamcolumn_reference_linear`,
    via `run_beamcolumn_frame`) — the 1D idealization;
  - RC lattice (`rc_lattice` / `rc_lattice_linear`).
Each comes in a NONLINEAR flavour (Concrete02/Steel02, length-regularized struts) and a LINEAR
flavour (Elastic everything). A 2D plane-stress continuum frame is the alternative reference
(`continuum_reference`, D29). `calibrate_area` matches the plain elastic lattice's strut area to the
reference K0; `calibrate_area_linear` refines that for the full RC-topology elastic lattice (rebar
adds a fixed parallel stiffness, so K0 is not proportional to area) via a secant root-find on
`lattice_k0`. Mirrors `examples/column/build.py`.
"""

from __future__ import annotations

from rclattice.builders import build_continuum_rc, build_lattice, build_lattice_rc
from rclattice.calibration import calibrate_lattice, continuum_targets
from rclattice.materials import (
    concrete_nd_nonlinear,
    concrete_uniaxial_elastic,
    concrete_uniaxial_nonlinear,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
)
from rclattice.opensees import run_beamcolumn_frame, run_dynamic, run_pushover

from specimen import (
    BEAM, CORE, COVER_C, DU, GF, GFC, H, HORIZON, MESH, P, SPAN, STEEL, TARGET, THK,
    add_axial_mass, control_base_nodes, frame_problem, lateral_loads, longitudinal_rebars,
    rebars, zone_of,
)


# --- nonlinear models (Concrete02 / Steel02) --------------------------------

def _grade(zone: str):
    """Map a zone name (column "core"/"cover" or beam "beam_core"/"beam_cover") to its concrete
    grade. Columns and the beam share the SAME confined CORE / unconfined COVER_C grades; the beam
    zones are kept distinct only so the beam concrete law can be toggled elastic/nonlinear."""
    return COVER_C if zone in ("cover", "beam_cover") else CORE


def _beam_fiber_materials(beam_nonlinear: bool):
    """The fiber-frame beam's concrete materials for `run_beamcolumn_frame[_dynamic]`'s
    `beam_materials`: None (reuse the columns' nonlinear core/cover) when `beam_nonlinear`, else an
    Elastic (core, cover) pair so the fiber beam matches the lattice's elastic beam (the default)."""
    if beam_nonlinear:
        return None
    return (concrete_uniaxial_elastic(CORE, 0), concrete_uniaxial_elastic(COVER_C, 0))


def make_material_for(regularize: bool, Gf: float, Gfc: float, beam_nonlinear: bool = False):
    """Build the per-strut concrete material factory passed to `build_lattice_rc`.

    The columns are ALWAYS nonlinear Concrete02 (`regularize=True` -> crack-band length-regularized,
    D20; else plain). The thin BEAM defaults to ELASTIC concrete (`beam_nonlinear=False`) because a
    softening axial-strut beam forms a local lattice mechanism the static pushover cannot trace; pass
    `beam_nonlinear=True` to give the beam the SAME nonlinear Concrete02 as the columns (stable under
    the transient dynamic runs). The grade (CORE/COVER_C) is picked per `zone` either way."""
    def material_for(zone: str, length: float):
        grade = _grade(zone)
        if zone.startswith("beam") and not beam_nonlinear:
            return concrete_uniaxial_elastic(grade, 0)
        if regularize:
            return concrete_uniaxial_regularized(grade, 0, length, Gf=Gf, Gfc=Gfc)
        return concrete_uniaxial_nonlinear(grade, 0)
    return material_for


def beamcolumn_reference(beam_nonlinear: bool = False) -> dict:
    """Fiber `forceBeamColumn` portal-frame pushover, material-matched to the lattice (the SAME
    grade-level Concrete02 core/cover + Steel02 trio). The columns are nonlinear; the beam concrete
    is Elastic by default (matching the lattice's default elastic beam) or nonlinear when
    `beam_nonlinear`. The lattice's length regularization is strut-specific (no fiber-section
    analog), so the fiber sections use the plain grades."""
    materials = (concrete_uniaxial_nonlinear(CORE, 1),
                 concrete_uniaxial_nonlinear(COVER_C, 2),
                 steel_uniaxial(STEEL, 3))
    return run_beamcolumn_frame(height=H, span=SPAN, beam_depth=BEAM, beam_width=THK, P=P,
                                dU=DU, target=TARGET, materials=materials,
                                beam_materials=_beam_fiber_materials(beam_nonlinear))


def _continuum_model(Gf: float = GF, Gfc: float = GFC):
    """Build the RC continuum frame (D29): structured quads with nonlinear nD concrete (ASDConcrete3D
    + PlaneStress, the SAME core/cover grades as the lattice, crack-band regularized to the quad size,
    PURE-DAMAGE hysteresis to match Concrete02 — D30) + the longitudinal bars (columns 3+2+3 + beam
    top/bottom) as steel struts on shared nodes (no stirrups; the 2D continuum supplies the
    lateral/shear path). Returns (model, control_node, base_nodes). Shared by the static reference,
    the K0 probe, and the dynamic reference."""
    def nd_material_for(zone: str):
        return concrete_nd_nonlinear(_grade(zone), 0, 0, lch=MESH, Gf=Gf, Gfc=Gfc)

    model, _ = build_continuum_rc(frame_problem(CORE), MESH, nd_material_for=nd_material_for,
                                  zone_of=zone_of, rebars=longitudinal_rebars())
    ctrl, base = control_base_nodes(model)
    return model, ctrl, base


def continuum_reference(*, dU: float = DU, target: float = TARGET, Gf: float = GF, Gfc: float = GFC) -> dict:
    """2D plane-stress continuum frame pushover, material-matched to the lattice (D29) — the
    alternative reference. Captures the 2D load-spreading / diagonal action a 1D fiber frame cannot,
    so it is the apples-to-apples reference. Heavier than the fiber frame (several-thousand
    ASDConcrete3D quads). Returns {"disp", "shear", "converged"} (shear = base reaction sum)."""
    model, ctrl, base = _continuum_model(Gf, Gfc)
    return run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                        control_dof=1, dU=dU, target=target, base_nodes=base)


def continuum_k0(*, Gf: float = GF, Gfc: float = GFC) -> float:
    """Initial lateral stiffness of the RC continuum frame (one tiny pushover step) — the calibration
    target for the dynamic --reference continuum, mirroring how the fiber-frame K0 is taken (D30)."""
    model, ctrl, base = _continuum_model(Gf, Gfc)
    r = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                     control_dof=1, dU=1e-3, target=1e-3, base_nodes=base)
    return r["shear"][1] / r["disp"][1]


def continuum_dynamic(*, accel, dt_record: float, scale: float, top_mass: float, dt: float = 0.01,
                      Gf: float = GF, Gfc: float = GFC) -> dict:
    """Seismic time-history of the RC continuum frame (D30) — the dynamic counterpart of
    `continuum_reference`. Bakes the axial-load tributary seismic mass `top_mass` (per column) onto
    each column's top nodes (on top of the builder's lumped self-mass, exactly as the lattice does),
    then runs the SAME scaled UniformExcitation via `run_dynamic`."""
    model, ctrl, base = _continuum_model(Gf, Gfc)
    add_axial_mass(model, top_mass)
    return run_dynamic(model, accel=accel, dt_record=dt_record, scale=scale, control_node=ctrl,
                       control_dof=1, base_nodes=base, dt=dt)


def make_reference(name: str, beam_nonlinear: bool = False) -> dict:
    """Dispatch the verification reference by name (mirrors `make_excitation`): "beamcolumn" (fast
    fiber `forceBeamColumn` frame, beam concrete elastic/nonlinear per `beam_nonlinear`) or
    "continuum" (2D plane-stress quads, D29 — the continuum is stable so its beam is ALWAYS nonlinear,
    regardless of `beam_nonlinear`). Returns its pushover curve."""
    if name == "beamcolumn":
        return beamcolumn_reference(beam_nonlinear)
    if name == "continuum":
        return continuum_reference()
    raise ValueError(f"unknown reference {name!r} (expected 'beamcolumn' or 'continuum')")


def calibrate_area(k_ref: float, *, horizon: float = HORIZON):
    """SCALAR calibration: one uniform strut area giving the elastic lattice frame the same initial
    lateral stiffness as `k_ref`. Returns (area, control_node, base_nodes). `horizon` sets the strut
    connectivity (larger = more redundant bracing, D31)."""
    lat0, _ = build_lattice(frame_problem(CORE), MESH, strut_area=1.0, horizon=horizon)
    ctrl, base = control_base_nodes(lat0)
    r = run_pushover(lat0, lateral_loads=lateral_loads(lat0), control_node=ctrl,
                     control_dof=1, dU=DU, target=DU, base_nodes=base)
    return k_ref / (r["shear"][-1] / r["disp"][-1]), ctrl, base


def calibrate_groups(*, n_modes: int = 3, horizon: float = HORIZON):
    """STRONG 2-group calibration (D16): fit orthogonal (length ~ s) and diagonal (length ~ s*sqrt2)
    strut areas, by bounded least squares, to the CONTINUUM's static deflection + first `n_modes`
    modal periods. Returns (area_fn, control_node, base_nodes, CalibrationResult). area_fn maps strut
    length -> area, accepted directly by `build_lattice_rc(strut_area=...)`."""
    problem = frame_problem(CORE)
    targets = continuum_targets(problem, MESH, n_modes=n_modes)
    result = calibrate_lattice(problem, MESH, targets=targets, horizon=horizon, n_modes=n_modes)
    lat0, _ = build_lattice(problem, MESH, strut_area=result.area_fn, horizon=horizon)
    ctrl, base = control_base_nodes(lat0)
    return result.area_fn, ctrl, base, result


def rc_lattice(regularize: bool, Gf: float, Gfc: float, area, *, beam_nonlinear: bool = False,
               horizon: float = HORIZON):
    """The calibrated RC lattice frame model (corotTruss struts + rebar). Columns are nonlinear
    Concrete02; the thin beam concrete is Elastic by default or nonlinear when `beam_nonlinear` (see
    `make_material_for`). `area` is a uniform float or a length->area callable; `horizon` sets strut
    connectivity (larger = more redundant bracing against the post-peak mechanism, D31)."""
    model, _ = build_lattice_rc(frame_problem(CORE), MESH,
                                material_for=make_material_for(regularize, Gf, Gfc, beam_nonlinear),
                                zone_of=zone_of, rebars=rebars(), strut_area=area,
                                horizon=horizon, strut_element="corotTruss")
    return model


# --- linear models (Elastic everything) -------------------------------------

def beamcolumn_reference_linear() -> dict:
    """The fiber `forceBeamColumn` frame with the EXACT sections but Elastic fibers (concrete E core +
    cover, steel E0). Linear, so its pushover is a straight line of slope K0 (the transformed-section
    lateral stiffness with gravity P-Delta)."""
    materials = (concrete_uniaxial_elastic(CORE, 1),
                 concrete_uniaxial_elastic(COVER_C, 2),
                 steel_uniaxial_elastic(STEEL, 3))
    return run_beamcolumn_frame(height=H, span=SPAN, beam_depth=BEAM, beam_width=THK, P=P,
                                dU=DU, target=TARGET, materials=materials)


def rc_lattice_linear(area: float):
    """The SAME RC lattice topology as the nonlinear case (concrete struts + longitudinal bars +
    stirrups, corotTruss) but fully linear: Elastic concrete struts (E per zone) and Elastic rebar
    (steel E0)."""
    model, _ = build_lattice_rc(
        frame_problem(CORE), MESH,
        material_for=lambda zone, length: concrete_uniaxial_elastic(_grade(zone), 0),
        zone_of=zone_of, rebars=rebars(), strut_area=area, horizon=HORIZON,
        strut_element="corotTruss", rebar_material=steel_uniaxial_elastic,
    )
    return model


def lattice_k0(area: float, ctrl: int, base: list[int]) -> float:
    """Initial lateral stiffness K0 of the elastic RC lattice frame at concrete strut area `area`:
    gravity (held) then one tiny DisplacementControl step; K0 = base shear / control disp. Linear,
    so a single small step is exact."""
    model = rc_lattice_linear(area)
    r = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                     control_dof=1, dU=1e-3, target=1e-3, base_nodes=base)
    return r["shear"][1] / r["disp"][1]


def calibrate_area_linear(k_ref: float, *, tol: float = 5e-4, max_iter: int = 15):
    """Concrete strut area giving the FULL elastic RC lattice frame (concrete + rebar + stirrups) the
    same initial lateral stiffness as `k_ref`.

    The rebar adds a fixed parallel stiffness, so K0(area) is monotone but not proportional — a short
    secant root-find on K0(area) - k_ref is used, seeded by the plain-lattice area (rebar ignored)
    from `calibrate_area`, which also supplies the control/base node ids. Returns (area, control_node,
    base_nodes)."""
    a_plain, ctrl, base = calibrate_area(k_ref)   # plain-elastic-lattice seed + node ids
    a0, a1 = a_plain, 0.85 * a_plain              # rebar adds stiffness -> true area < plain area
    f0 = lattice_k0(a0, ctrl, base) - k_ref
    f1 = lattice_k0(a1, ctrl, base) - k_ref
    for _ in range(max_iter):
        if abs(f1) <= tol * k_ref or f1 == f0:
            break
        a2 = max(a1 - f1 * (a1 - a0) / (f1 - f0), 1e-4)   # secant step, keep area positive
        a0, f0, a1 = a1, f1, a2
        f1 = lattice_k0(a1, ctrl, base) - k_ref
    return a1, ctrl, base
