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

from rclattice import viz
from rclattice.builders import build_continuum_rc, build_lattice, build_lattice_rc, select_nodes
from rclattice.calibration import calibrate_lattice, continuum_targets
from rclattice.materials import (
    concrete_nd_elastic_planestress,
    concrete_nd_nonlinear,
    concrete_uniaxial_elastic,
    concrete_uniaxial_nonlinear,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
)
from rclattice.opensees import (
    run_beamcolumn_cantilever, run_beamcolumn_modal, run_dynamic, run_modal, run_pushover,
)

from specimen import (
    CORE, COVER_C, DU, EPS, GF, GFC, H, HORIZON, MESH, P, STEEL, TARGET, W,
    column_problem, lateral_loads, longitudinal_rebars, rebars, zone_of,
)

N_MODES = 3   # first N modal periods drawn + tabulated in the calibration figure (D35)


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


def _continuum_model(Gf: float = GF, Gfc: float = GFC):
    """Build the RC continuum column (D29): structured quads with nonlinear nD concrete (ASDConcrete3D
    + PlaneStress, the SAME core/cover grades, crack-band regularized to the quad size, PURE-DAMAGE
    hysteresis to match Concrete02 — D30) + the three longitudinal bars as steel struts on shared
    nodes (no stirrups; the 2D continuum supplies the lateral/shear path the lattice gets from them).
    Returns (model, control_node, base_nodes). Shared by the static reference, the K0 probe, and the
    dynamic reference."""
    def nd_material_for(zone: str):
        grade = CORE if zone == "core" else COVER_C
        return concrete_nd_nonlinear(grade, 0, 0, lch=MESH, Gf=Gf, Gfc=Gfc)

    model, _ = build_continuum_rc(column_problem(CORE), MESH, nd_material_for=nd_material_for,
                                  zone_of=zone_of, rebars=longitudinal_rebars())
    ctrl = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(model, (-W, W, -EPS, EPS))
    return model, ctrl, base


def continuum_reference(*, dU: float = 0.1, target: float = TARGET, Gf: float = GF, Gfc: float = GFC) -> dict:
    """2D plane-stress continuum pushover, material-matched to the lattice (D29) — the alternative
    reference. Captures the 2D load-spreading / diagonal action a 1D fiber section cannot, so it is
    the apples-to-apples reference. Heavier than the beam-column (1536 ASDConcrete3D quads) → coarser
    `dU`. Returns the pushover curve {"disp", "shear", "converged"} (shear = base reaction sum)."""
    model, ctrl, base = _continuum_model(Gf, Gfc)
    return run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                        control_dof=1, dU=dU, target=target, base_nodes=base)


def continuum_k0(*, Gf: float = GF, Gfc: float = GFC) -> float:
    """Initial lateral stiffness of the RC continuum column (one tiny pushover step) — the calibration
    target for the dynamic --reference continuum, mirroring how the beam-column K0 is taken (D30)."""
    model, ctrl, base = _continuum_model(Gf, Gfc)
    r = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                     control_dof=1, dU=1e-3, target=1e-3, base_nodes=base)
    return r["shear"][1] / r["disp"][1]


def continuum_dynamic(*, accel, dt_record: float, scale: float, top_mass: float, dt: float = 0.01,
                      Gf: float = GF, Gfc: float = GFC) -> dict:
    """Seismic time-history of the RC continuum column (D30) — the dynamic counterpart of
    `continuum_reference`, the reference for the lattice's nonlinear dynamic run. Bakes the axial-load
    tributary seismic mass `top_mass` onto the top nodes (on top of the builder's lumped self-mass,
    exactly as the lattice does), then runs the SAME scaled UniformExcitation via `run_dynamic`.
    Returns the `_transient_uniform_excitation` dict (t/disp/shear histories + peaks) plus periods."""
    model, ctrl, base = _continuum_model(Gf, Gfc)
    top = select_nodes(model, (-W, W, H - EPS, H + EPS))
    for nid in top:
        mx, my = model.masses[nid]
        model.masses[nid] = (mx + top_mass / len(top), my + top_mass / len(top))
    return run_dynamic(model, accel=accel, dt_record=dt_record, scale=scale, control_node=ctrl,
                       control_dof=1, base_nodes=base, dt=dt)


def make_reference(name: str) -> dict:
    """Dispatch the verification reference by name (mirrors `make_excitation`): "beamcolumn" (fast
    fiber `forceBeamColumn`) or "continuum" (2D plane-stress quads, D29). Returns its pushover curve."""
    if name == "beamcolumn":
        return beamcolumn_reference()
    if name == "continuum":
        return continuum_reference()
    raise ValueError(f"unknown reference {name!r} (expected 'beamcolumn' or 'continuum')")


def modal_calibration_figure(*, reference: str, lattice_model, label: str, caption: str,
                             savepath: str, linear: bool = False, n_modes: int = N_MODES,
                             Gf: float = GF, Gfc: float = GFC) -> tuple[list, list]:
    """Calibration output (D35): draw the first `n_modes` mode shapes of the calibrated lattice and
    the SELECTED reference, with a periods table underneath, to `savepath`.

    The reference is whatever `--reference` picked: the 2D RC continuum (`run_modal` on
    `_continuum_model`) or the subdivided fiber beam-column (`run_beamcolumn_modal`). Both are put on
    a mass-consistent footing with the lattice so the periods are directly comparable — the continuum
    shares the builder tributary mass by construction (D16); the beam-column is given the lattice's
    total self-mass. The lattice is the as-built calibrated model (self-mass only, no seismic top
    mass — the modal basis the calibration itself uses). The fit target stays the continuum (D16), so
    the table's column is a plain `Δ vs reference`, not a calibration residual. Returns
    (reference_periods, lattice_periods)."""
    lat = run_modal(lattice_model, n_modes)
    lat_total = sum(m[0] for m in lattice_model.masses.values())   # x-direction total mass

    if reference == "continuum":
        ref_model = (_continuum_model_linear() if linear else _continuum_model(Gf, Gfc))[0]
        ref = run_modal(ref_model, n_modes)
    else:  # beamcolumn — subdivided fiber stick, total mass matched to the lattice
        mats = ((concrete_uniaxial_elastic(CORE, 1), concrete_uniaxial_elastic(COVER_C, 2),
                 steel_uniaxial_elastic(STEEL, 3)) if linear else
                (concrete_uniaxial_nonlinear(CORE, 1), concrete_uniaxial_nonlinear(COVER_C, 2),
                 steel_uniaxial(STEEL, 3)))
        ref = run_beamcolumn_modal(height=H, materials=mats, nelem=int(round(H / MESH)),
                                   self_mass=lat_total, num_modes=n_modes)
        ref_model = ref["model"]

    rows = []
    for i in range(n_modes):
        tr, tl = ref["periods"][i], lat["periods"][i]
        d = 100.0 * (tl - tr) / tr if tr not in (0.0, float("inf")) else float("nan")
        rows.append((i + 1, tr, tl, d))

    viz.figure_modal_calibration(
        {"model": ref_model, "shapes": ref["shapes"], "label": label, "color": "C3"},
        {"model": lattice_model, "shapes": lat["shapes"], "label": "RC lattice", "color": "C0"},
        rows, caption=caption, savepath=savepath,
        title=f"Modal calibration: {label} vs RC lattice (first {n_modes} modes)",
    )
    return ref["periods"], lat["periods"]


def calibrate_area(k_bc: float, *, horizon: float = HORIZON):
    """SCALAR calibration: one uniform strut area giving the elastic lattice the same initial lateral
    stiffness as `k_bc`. Matches a single number (K0). Returns (area, control_node, base_nodes).
    `horizon` sets the strut connectivity (larger = more redundant bracing, D31)."""
    lat0, _ = build_lattice(column_problem(CORE), MESH, strut_area=1.0, horizon=horizon)
    ctrl = select_nodes(lat0, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(lat0, (-W, W, -EPS, EPS))
    r = run_pushover(lat0, lateral_loads=lateral_loads(lat0), control_node=ctrl,
                     control_dof=1, dU=DU, target=DU, base_nodes=base)
    return k_bc / (r["shear"][-1] / r["disp"][-1]), ctrl, base


def calibrate_groups(*, n_modes: int = 3, horizon: float = HORIZON):
    """STRONG 2-group calibration (D16, `calibration.py`): fit orthogonal (length ~ s) and diagonal
    (length ~ s*sqrt2) strut areas, by bounded least squares, to the CONTINUUM's static deflection +
    first `n_modes` modal periods. Two area groups tune the effective stiffness AND the Poisson/shear
    coupling (the diagonal/orthogonal ratio) — not just the single K0 scalar — so the lattice's elastic
    stress field tracks the continuum better (trims the diagonal-strut overstrength). The plain
    concrete skeletons are matched (no rebar, like `calibrate_area`); the rebar is then added equally
    to both. Returns (area_fn, control_node, base_nodes, CalibrationResult). area_fn maps strut length
    -> area, accepted directly by `build_lattice_rc(strut_area=...)`."""
    problem = column_problem(CORE)
    targets = continuum_targets(problem, MESH, n_modes=n_modes)
    result = calibrate_lattice(problem, MESH, targets=targets, horizon=horizon, n_modes=n_modes)
    lat0, _ = build_lattice(problem, MESH, strut_area=result.area_fn, horizon=horizon)
    ctrl = select_nodes(lat0, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(lat0, (-W, W, -EPS, EPS))
    return result.area_fn, ctrl, base, result


def rc_lattice(regularize: bool, Gf: float, Gfc: float, area, *, horizon: float = HORIZON):
    """The calibrated RC lattice column model (corotTruss struts + rebar). `area` is a uniform float
    or a length->area callable; `horizon` sets strut connectivity (larger = more redundant bracing
    against the post-peak mechanism, D31)."""
    model, _ = build_lattice_rc(column_problem(CORE), MESH,
                                material_for=make_material_for(regularize, Gf, Gfc),
                                zone_of=zone_of, rebars=rebars(), strut_area=area,
                                horizon=horizon,
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


# --- linear continuum reference (ElasticIsotropic + PlaneStress quads) ------

def _continuum_model_linear():
    """Build the LINEAR elastic RC continuum column: ElasticIsotropic+PlaneStress plane-stress quads
    + elastic rebar struts on the shared nodes. The linear analog of ``_continuum_model``; used by
    ``--reference continuum`` in the linear pushover and dynamic scripts."""
    def nd_material_for(zone: str) -> tuple:
        grade = CORE if zone == "core" else COVER_C
        return concrete_nd_elastic_planestress(grade, 0)   # tags assigned by the builder

    model, _ = build_continuum_rc(column_problem(CORE), MESH, nd_material_for=nd_material_for,
                                  zone_of=zone_of, rebars=longitudinal_rebars(),
                                  rebar_material=steel_uniaxial_elastic)
    ctrl = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(model, (-W, W, -EPS, EPS))
    return model, ctrl, base


def continuum_reference_linear(*, dU: float = 0.1, target: float = TARGET) -> dict:
    """LINEAR elastic 2D plane-stress continuum pushover — the continuum reference for the linear
    pushover study. Same mesh and geometry as the nonlinear continuum but with ElasticIsotropic
    quads so the curve is a straight line of slope K0."""
    model, ctrl, base = _continuum_model_linear()
    return run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                        control_dof=1, dU=dU, target=target, base_nodes=base)


def continuum_k0_linear() -> float:
    """Initial lateral stiffness of the LINEAR elastic RC continuum column (one tiny pushover step)
    — the calibration target when ``--reference continuum`` is used in the linear studies."""
    model, ctrl, base = _continuum_model_linear()
    r = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                     control_dof=1, dU=1e-3, target=1e-3, base_nodes=base)
    return r["shear"][1] / r["disp"][1]


def continuum_dynamic_linear(*, accel, dt_record: float, scale: float,
                              top_mass: float, dt: float = 0.01) -> dict:
    """Seismic time-history of the LINEAR elastic RC continuum column — the dynamic continuum
    reference for ``--reference continuum`` in ``dynamic_linear.py``. Bakes the axial-load
    tributary seismic mass ``top_mass`` onto the top nodes (same as the lattice) then runs the
    scaled UniformExcitation via ``run_dynamic``."""
    model, ctrl, base = _continuum_model_linear()
    top = select_nodes(model, (-W, W, H - EPS, H + EPS))
    for nid in top:
        mx, my = model.masses[nid]
        model.masses[nid] = (mx + top_mass / len(top), my + top_mass / len(top))
    return run_dynamic(model, accel=accel, dt_record=dt_record, scale=scale, control_node=ctrl,
                       control_dof=1, base_nodes=base, dt=dt)


def make_reference_linear(name: str) -> dict:
    """Dispatch the LINEAR verification reference by name: "beamcolumn" or "continuum"."""
    if name == "beamcolumn":
        return beamcolumn_reference_linear()
    if name == "continuum":
        return continuum_reference_linear()
    raise ValueError(f"unknown reference {name!r} (expected 'beamcolumn' or 'continuum')")
