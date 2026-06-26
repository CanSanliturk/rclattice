"""OpenSees backend — the ONLY module allowed to import openseespy (D8).

Translates a generic FE `Model` into OpenSees commands and runs analyses: a linear static
step (`run_static`), an eigen analysis (`run_modal`), and the staged RC-frame pushover
machinery (`run_gravity` LoadControl + `run_pushover` DisplacementControl with base-shear
recording, D18/D19). The pushover runners are written to also carry the Stage-2 nonlinear
case (Newton stepping with step reduction); Stage 1 exercises them elastically.
"""

from __future__ import annotations

import math

import openseespy.opensees as ops

from .model import Model


def build(model: Model) -> None:
    """Emit `model` into the current OpenSees domain (after ops.wipe / model setup)."""
    ops.model("basic", "-ndm", model.ndm, "-ndf", model.ndf)

    for node in model.nodes.values():
        ops.node(node.id, *node.coords)
    for nid, mvals in model.masses.items():
        ops.mass(nid, *mvals)
    for mat in model.uniaxial_materials:
        ops.uniaxialMaterial(mat.mtype, mat.id, *mat.args)
    for mat in model.nd_materials:
        ops.nDMaterial(mat.mtype, mat.id, *mat.args)
    for el in model.elements:
        ops.element(el.etype, el.id, *el.nodes, *el.args)
    for sup in model.supports:
        ops.fix(sup.node, *sup.fix)


def run_static(model: Model) -> dict:
    """Build `model` and run one linear static load step.

    Returns {"ok": <0 = success>, "disps": {node_id: [u...]}}.
    """
    ops.wipe()
    build(model)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for load in model.loads:
        ops.load(load.node, *load.values)

    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    ok = ops.analyze(1)

    disps = {nid: ops.nodeDisp(nid) for nid in model.nodes}
    return {"ok": ok, "disps": disps}


def run_modal(model: Model, num_modes: int) -> dict:
    """Build `model` and run an eigenvalue analysis (D16).

    Mass must be assigned on `model` (builders do this from density). Returns
    {"eigenvalues": [...], "periods": [Ti...]} sorted ascending by eigenvalue. A
    non-positive eigenvalue (rigid-body / mechanism / spurious mode) yields period inf and is
    a signal that the lattice is under-constrained.
    """
    if not model.masses:
        raise ValueError("run_modal requires nodal mass on the model (none assigned)")
    ops.wipe()
    build(model)
    try:
        eigenvalues = ops.eigen(num_modes)  # default (-genBandArpack): fast for a few modes
    except Exception:
        eigenvalues = ops.eigen("-fullGenLapack", num_modes)  # robust fallback (small models)
    periods = [2.0 * math.pi / math.sqrt(lam) if lam > 0.0 else math.inf for lam in eigenvalues]
    shapes = [
        {nid: ops.nodeEigenvector(nid, mode) for nid in model.nodes}
        for mode in range(1, num_modes + 1)
    ]
    return {"eigenvalues": list(eigenvalues), "periods": periods, "shapes": shapes}


# --- staged pushover machinery (D18/D19) ------------------------------------

def _gravity_loads(model: Model, gravity_loads) -> list:
    """The vertical (gravity) load case: explicit list if given, else the model's loads."""
    return list(gravity_loads) if gravity_loads is not None else list(model.loads)


def _apply_gravity(model: Model, loads, nsteps: int, tol: float) -> int:
    """Set up + run the constant gravity case as LoadControl pattern 1. Returns the rc.

    The solver (system/numberer/constraints) must already be configured. Leaves the gravity
    pattern in the domain; caller does `ops.loadConst` to hold it constant into the pushover.
    """
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for ld in loads:
        ops.load(ld.node, *ld.values)
    ops.test("NormDispIncr", tol, 100)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    return ops.analyze(nsteps)


def run_gravity(model: Model, *, gravity_loads=None, nsteps: int = 10, tol: float = 1e-8) -> dict:
    """Build `model` and apply the gravity loads in `nsteps` LoadControl increments (D18).

    Standalone (wipes + builds + analyzes). `gravity_loads` is a list of `Load`; if omitted the
    model's own loads are used as the gravity case. Returns {"ok", "disps"}.
    """
    ops.wipe()
    build(model)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ok = _apply_gravity(model, _gravity_loads(model, gravity_loads), nsteps, tol)
    disps = {nid: ops.nodeDisp(nid) for nid in model.nodes}
    return {"ok": ok, "disps": disps}


def run_pushover(
    model: Model,
    *,
    lateral_loads,
    control_node: int,
    control_dof: int,
    dU: float,
    target: float,
    gravity_loads=None,
    gravity_steps: int = 10,
    base_nodes=None,
    tol: float = 1e-5,
    max_iter: int = 100,
    element_groups=None,
) -> dict:
    """Gravity (constant) → DisplacementControl pushover, recording base shear (D18/D19).

    Sequence (mirrors the OpenSees RCFrameGravity → RCFramePushOver benchmark):
      1. apply `gravity_loads` (or `model.loads`) as a constant pattern;
      2. add `lateral_loads` (a list of `Load`) as the reference lateral pattern;
      3. step `control_node`'s `control_dof` by `dU` (signed by `target`) up to `target`,
         summing horizontal base reactions into the base shear at each converged step.

    Base shear is `-sum(reaction[control_dof])` over `base_nodes` (defaults to the supported
    nodes) — the structure's resistance, positive for a positive push. A failed step is retried
    with finer sub-steps and a stronger algorithm (Newton → KrylovNewton) before giving up — the
    Stage-2 nonlinear lattice (Concrete02 softening) needs this. `tol`/`max_iter` set the
    NormDispIncr test (1e-6 is the practical tolerance for the softening lattice).

    `element_groups` (optional) is a force-decomposition probe for diagnostics: a dict mapping a
    label to a list of `(element_id, dof, coef)`. At each converged step it records, per label, the
    sum of `eleForce(element_id)[dof] * coef` — a GLOBAL nodal force component the element exerts
    (so the corotational/P-Δ geometry is already baked in, unlike a precomputed direction cosine).
    With `dof` = the horizontal index at an element's node-on-one-side-of-a-cut and `coef`=1, the
    per-label sums reconcile to the base shear; pairing the vertical index with `coef`=that node's x
    gives the overturning-moment share. Lets a caller attribute base shear/overturning to element
    categories (vertical vs diagonal struts, concrete vs rebar) without any ops.* calls of its own.

    Returns {"ok", "converged", "disp": [...], "shear": [...], "control_node", "base_nodes"}, plus
    "groups": {label: [...]} when `element_groups` is given.
    """
    ops.wipe()
    build(model)
    base = list(base_nodes) if base_nodes is not None else [s.node for s in model.supports]

    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")

    grav = _gravity_loads(model, gravity_loads)
    if grav:
        if _apply_gravity(model, grav, gravity_steps, tol) != 0:
            return {"ok": -1, "converged": False, "stage": "gravity",
                    "disp": [], "shear": [], "control_node": control_node, "base_nodes": base}
        ops.loadConst("-time", 0.0)

    # lateral reference pattern (its magnitude is just a shape — DisplacementControl drives it)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for ld in lateral_loads:
        ops.load(ld.node, *ld.values)

    ops.test("NormDispIncr", tol, max_iter)
    ops.algorithm("Newton")  # Stage 1 elastic; ModifiedNewton is an option for Stage 2
    ops.analysis("Static")

    sign = 1.0 if target >= 0.0 else -1.0
    du = sign * abs(dU)
    ops.integrator("DisplacementControl", control_node, control_dof, du)

    def control_disp() -> float:
        return ops.nodeDisp(control_node, control_dof)

    def base_shear() -> float:
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base)

    groups: dict[str, list[float]] = {label: [] for label in (element_groups or {})}
    disp: list[float] = []
    shear: list[float] = []

    def record() -> None:
        disp.append(control_disp())
        shear.append(base_shear())
        for label, members in (element_groups or {}).items():
            groups[label].append(sum(ops.eleForce(eid)[dof] * coef for eid, dof, coef in members))

    record()
    converged = True
    goal = abs(target)
    while sign * control_disp() < goal - 1e-9 * (goal + 1.0):
        if ops.analyze(1) == 0:
            record()
            continue
        # retry the increment with finer sub-steps and stronger algorithms (softening, D20)
        sub_ok = False
        for algo, args in (("Newton", ()), ("KrylovNewton", ()), ("NewtonLineSearch", ("-type", "Bisection"))):
            ops.algorithm(algo, *args)
            for nsub in (5, 20, 50):
                ops.integrator("DisplacementControl", control_node, control_dof, du / nsub)
                if ops.analyze(nsub) == 0:
                    record()
                    sub_ok = True
                    break
            if sub_ok:
                break
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", control_node, control_dof, du)  # restore full step
        if not sub_ok:
            converged = False
            break

    result = {"ok": 0 if converged else -1, "converged": converged, "disp": disp,
              "shear": shear, "control_node": control_node, "base_nodes": base}
    if element_groups is not None:
        result["groups"] = groups
    return result


def run_pushover_dynamic(
    model: Model,
    *,
    control_node: int,
    control_dof: int,
    target: float,
    drive_nodes=None,
    gravity_loads=None,
    gravity_steps: int = 10,
    base_nodes=None,
    periods_to_target: float = 12.0,
    steps_per_period: int = 40,
    damping_ratio: float = 0.6,
    tol: float = 1e-5,
    max_iter: int = 50,
) -> dict:
    """Dynamic-relaxation pushover: a quasi-static TRANSIENT solve that rides through limit
    points / local snap-backs the static Newton can't (D22, user-selected).

    After the constant-gravity stage, the lateral displacement is *imposed* (ramped) on
    `drive_nodes` (default [control_node]) via single-point constraints, and the model is marched
    with Newmark + heavy Rayleigh damping slowly enough (`periods_to_target` fundamental periods
    to reach `target`) that inertia stays negligible — so the recorded base shear vs control
    displacement is the quasi-static pushover, but the mass/damping regularize the instability.
    Returns {"converged", "disp", "shear", "control_node", "base_nodes"}.
    """
    if not model.masses:
        raise ValueError("run_pushover_dynamic requires nodal mass on the model")
    ops.wipe()
    build(model)
    base = list(base_nodes) if base_nodes is not None else [s.node for s in model.supports]
    drive = list(drive_nodes) if drive_nodes is not None else [control_node]

    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")

    grav = _gravity_loads(model, gravity_loads)
    if grav:
        if _apply_gravity(model, grav, gravity_steps, tol) != 0:
            return {"converged": False, "stage": "gravity", "disp": [], "shear": [],
                    "control_node": control_node, "base_nodes": base}
        ops.loadConst("-time", 0.0)

    # fundamental frequency sets the (quasi-static) loading rate, step, and Rayleigh damping
    lam1 = ops.eigen("-fullGenLapack", 1)[0]
    w1 = math.sqrt(lam1) if lam1 > 0 else 1.0
    T1 = 2.0 * math.pi / w1
    total_time = periods_to_target * T1
    dt = T1 / steps_per_period
    rate = target / total_time
    ops.rayleigh(damping_ratio * w1, damping_ratio / w1, 0.0, 0.0)  # ~damping_ratio at w1

    # impose the ramped lateral displacement on the drive nodes (disp = rate * t)
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    for nid in drive:
        ops.sp(nid, control_dof, rate)

    ops.test("NormDispIncr", tol, max_iter)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    def cdisp():
        return ops.nodeDisp(control_node, control_dof)

    def base_shear():
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base)

    disp = [cdisp()]; shear = [base_shear()]
    converged = True
    nsteps = int(round(total_time / dt))
    for _ in range(nsteps):
        print(f"Step {_}", flush=True)
        if ops.analyze(1, dt) != 0:
            sub_ok = False
            for algo in ("KrylovNewton", "NewtonLineSearch"):
                ops.algorithm(algo)
                if ops.analyze(10, dt / 10.0) == 0:
                    sub_ok = True
                    break
            ops.algorithm("Newton")
            if not sub_ok:
                converged = False
                break
        disp.append(cdisp()); shear.append(base_shear())
        if cdisp() >= target - 1e-9:
            break
    return {"converged": converged, "disp": disp, "shear": shear,
            "control_node": control_node, "base_nodes": base}


# --- nonlinear seismic time-history (UniformExcitation) ---------------------

def _transient_uniform_excitation(*, accel, dt_record: float, scale: float, dof: int,
                                  control_node: int, control_dof: int, base_nodes,
                                  dt: float, zeta: float = 0.05, modes: tuple[int, int] = (1, 2),
                                  tol: float = 1e-6, max_iter: int = 50) -> dict:
    """March a UniformExcitation transient on the already-built+massed+gravity-held domain. Sets
    `zeta` MODAL damping at the first max(`modes`) modes internally (D33), then applies the
    ground-acceleration record `accel` (a sequence at spacing `dt_record`) along `dof`, multiplied by
    `scale` (fold in g and the intensity scale factor). HHT-alpha (numerical damping); on a
    non-converged step, retry with finer sub-steps + stronger algorithms. Records relative control
    displacement (= drift, UniformExcitation is in coords) and base shear `-sum reaction[control_dof]`
    at every step. Returns time/disp/shear histories + peak |disp|, residual disp, peak |shear|, and
    the modal `periods`."""
    ops.timeSeries("Path", 90, "-dt", dt_record, "-values", *accel, "-factor", scale)
    ops.pattern("UniformExcitation", 90, dof, "-accel", 90)
    # Drop the gravity stage's Static analysis: while it exists OpenSees REFUSES to set a transient
    # integrator ("can't set transient integrator in static analysis") and silently falls back to the
    # default Newmark, so the HHT numerical damping below is ignored. wipeAnalysis clears only the
    # analysis aggregation (the domain, gravity loadConst, and mass all persist), so the
    # constraints/numberer/system must be re-declared before the transient integrator is honored.
    ops.wipeAnalysis()
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("BandGeneral")
    # Modal damping (D33): assign `zeta` to the first max(modes) modes, formed from the
    # mass-orthonormal eigenvectors on the gravity-held tangent. Unlike stiffness-proportional
    # Rayleigh (a1*K) it has NO term riding the committed tangent, so it does NOT spike the base
    # reaction when elements crack/yield at high velocity (the D28 artifact). Higher (uncomputed)
    # modes get no modal damping; the HHT below dissipates that residual high-frequency content.
    # modalDamping is stored on the domain, so the transient integrator created next uses it.
    nmode = max(modes)
    try:
        lams = ops.eigen(nmode)
    except Exception:
        lams = ops.eigen("-fullGenLapack", nmode)
    ops.modalDamping(zeta)
    w = [math.sqrt(l) if l > 0 else math.inf for l in lams]
    periods = [2.0 * math.pi / wk if wk > 0 else math.inf for wk in w]
    ops.test("NormDispIncr", tol, max_iter)
    ops.algorithm("Newton")
    ops.integrator("HHT", 0.7)   # numerical damping (alpha<1): aids convergence + dissipates the
    ops.analysis("Transient")    # high-frequency content modal damping leaves on the uncomputed
    #                              higher modes (still negligible at the structural period)

    def cdisp() -> float:
        return ops.nodeDisp(control_node, control_dof)

    def base_shear() -> float:
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base_nodes)

    total_time = (len(accel) - 1) * dt_record
    nsteps = int(round(total_time / dt))
    t = [0.0]; disp = [cdisp()]; shear = [base_shear()]
    converged = True
    algos = (("Newton", ()), ("KrylovNewton", ()),
             ("NewtonLineSearch", ("-type", "Bisection")), ("ModifiedNewton", ()))
    for k in range(nsteps):
        if ops.analyze(1, dt) != 0:
            ok = False
            ops.test("NormDispIncr", tol * 10, max_iter * 4)   # modest relax + try harder for this step
            for algo, a in algos:
                ops.algorithm(algo, *a)
                for nsub in (10, 50, 200):
                    if ops.analyze(nsub, dt / nsub) == 0:       # finer sub-steps, same dt advanced
                        ok = True
                        break
                if ok:
                    break
            ops.test("NormDispIncr", tol, max_iter)
            ops.algorithm("Newton")
            if not ok:
                converged = False
                break
        t.append((k + 1) * dt); disp.append(cdisp()); shear.append(base_shear())
    ipk = max(range(len(disp)), key=lambda m: abs(disp[m]))
    return {"t": t, "disp": disp, "shear": shear, "converged": converged, "periods": periods,
            "peak_disp": disp[ipk], "peak_time": t[ipk], "residual_disp": disp[-1],
            "peak_shear": max(shear, key=abs), "control_node": control_node,
            "base_nodes": list(base_nodes)}


def run_dynamic(model: Model, *, accel, dt_record: float, scale: float, control_node: int,
                control_dof: int, base_nodes, extra_mass=None, damping_ratio: float = 0.05,
                modes: tuple[int, int] = (1, 2), dt: float = 0.01, gravity_loads=None,
                gravity_steps: int = 10, tol: float = 1e-6) -> dict:
    """Nonlinear seismic time-history of a generic Model (the lattice) under UniformExcitation.

    Sequence: build → optionally ADD `extra_mass` (dict node_id -> mass, applied to both
    translational DOFs, e.g. tributary axial-load mass at the top) on top of the builder's lumped
    self-mass → apply gravity (constant, held with loadConst) → modal damping (`damping_ratio` at
    `modes`) → UniformExcitation transient of the scaled `accel` record along `control_dof`.
    Returns the _transient_uniform_excitation dict (which carries the modal "periods")."""
    if not model.masses and not extra_mass:
        raise ValueError("run_dynamic requires nodal mass on the model")
    ops.wipe()
    build(model)
    for nid, m in (extra_mass or {}).items():
        mx, my = model.masses.get(nid, (0.0, 0.0))
        ops.mass(nid, mx + m, my + m)   # ops.mass overwrites -> emit the summed value

    base = list(base_nodes) if base_nodes is not None else [s.node for s in model.supports]
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")

    grav = _gravity_loads(model, gravity_loads)
    if grav:
        if _apply_gravity(model, grav, gravity_steps, tol) != 0:
            return {"converged": False, "stage": "gravity", "t": [], "disp": [], "shear": [],
                    "control_node": control_node, "base_nodes": base}
        ops.loadConst("-time", 0.0)

    res = _transient_uniform_excitation(accel=accel, dt_record=dt_record, scale=scale,
                                        dof=control_dof, control_node=control_node,
                                        control_dof=control_dof, base_nodes=base, dt=dt,
                                        zeta=damping_ratio, modes=modes, tol=tol)
    return res


def _rc_fiber_section(sec_tag: int, *, materials=None) -> None:
    """The benchmark RC column fiber section (15 wide z x 24 deep y, 1.5 cover): core (mat 1) /
    cover (mat 2) / steel (mat 3), 3+2+3 bars @0.6 in^2. Defines materials 1-3 too.

    Default trio is the benchmark Concrete01 core / Concrete01 cover / Steel01. Pass `materials`
    as a triple of backend-agnostic `UniaxialMaterial` specs (core, cover, steel) — each carrying
    `.mtype` + `.args` — to emit those instead (e.g. material-match the lattice column's
    Concrete02/Steel02 grades; the section geometry is identical either way)."""
    if materials is None:
        ops.uniaxialMaterial("Concrete01", 1, -6.0, -0.004, -5.0, -0.014)  # confined core
        ops.uniaxialMaterial("Concrete01", 2, -5.0, -0.002, 0.0, -0.006)   # unconfined cover
        ops.uniaxialMaterial("Steel01", 3, 60.0, 30000.0, 0.01)
    else:
        core_mat, cover_mat, steel_mat = materials
        ops.uniaxialMaterial(core_mat.mtype, 1, *core_mat.args)   # confined core
        ops.uniaxialMaterial(cover_mat.mtype, 2, *cover_mat.args)  # unconfined cover
        ops.uniaxialMaterial(steel_mat.mtype, 3, *steel_mat.args)
    y1, z1, cover, As = 12.0, 7.5, 1.5, 0.60
    ops.section("Fiber", sec_tag)
    ops.patch("quad", 1, 1, 10, -y1 + cover, z1 - cover, -y1 + cover, -z1 + cover,
              y1 - cover, -z1 + cover, y1 - cover, z1 - cover)               # core
    ops.patch("quad", 2, 1, 1, -y1, z1, -y1, -z1, -y1 + cover, -z1 + cover, -y1 + cover, z1 - cover)
    ops.patch("quad", 2, 1, 1, y1 - cover, z1 - cover, y1 - cover, -z1 + cover, y1, -z1, y1, z1)
    ops.patch("quad", 2, 1, 1, -y1 + cover, z1 - cover, -y1 + cover, z1, y1 - cover, z1, y1 - cover, z1 - cover)
    ops.patch("quad", 2, 1, 1, -y1 + cover, -z1, -y1 + cover, -z1 + cover, y1 - cover, -z1 + cover, y1 - cover, -z1)
    ops.layer("straight", 3, 3, As, -y1 + cover, z1 - cover, -y1 + cover, -z1 + cover)  # 3 bars y=-10.5
    ops.layer("straight", 3, 2, As, 0.0, z1 - cover, 0.0, -z1 + cover)                  # 2 bars y=0
    ops.layer("straight", 3, 3, As, y1 - cover, z1 - cover, y1 - cover, -z1 + cover)    # 3 bars y=+10.5


def _displacement_pushover(control_node, control_dof, base_nodes, dU, target):
    """Robust DisplacementControl pushover loop on the already-set-up domain (gravity held
    constant, lateral pattern applied). Steps control_node/control_dof to target, recording
    (disp, base shear). On a failed step: finer sub-steps + stronger algorithms. Returns
    {"disp", "shear", "converged"}."""
    def base_shear():
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base_nodes)

    def cdisp():
        return ops.nodeDisp(control_node, control_dof)

    ops.integrator("DisplacementControl", control_node, control_dof, dU)
    disp = [cdisp()]; shear = [base_shear()]
    converged = True
    while cdisp() < target - 1e-9:
        if ops.analyze(1) == 0:
            disp.append(cdisp()); shear.append(base_shear())
            continue
        sub_ok = False
        for algo in ("Newton", "KrylovNewton"):
            ops.algorithm(algo)
            for nsub in (5, 20, 100):
                ops.integrator("DisplacementControl", control_node, control_dof, dU / nsub)
                if ops.analyze(nsub) == 0:
                    disp.append(cdisp()); shear.append(base_shear())
                    sub_ok = True
                    break
            if sub_ok:
                break
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", control_node, control_dof, dU)
        if not sub_ok:
            converged = False
            break
    return {"converged": converged, "disp": disp, "shear": shear}


def run_beamcolumn_cantilever(*, height: float = 144.0, P: float = 180.0,
                              dU: float = 0.05, target: float = 10.0, materials=None) -> dict:
    """Single RC cantilever column via a force-based fiber `forceBeamColumn` — the verification
    reference for the lattice column (D12 single-member check). Fixed base (node 1), free top
    (node 2) carrying constant axial `P` (down) then a DisplacementControl lateral pushover of
    the top in X. Same 15x24 fiber section as the frame benchmark, P-Delta transform, 5 Lobatto
    points. `materials` (optional core/cover/steel `UniaxialMaterial` triple) overrides the default
    Concrete01/Steel01 section trio — pass the lattice's grade materials to material-match it.
    Returns the pushover curve {"disp", "shear", "converged"} (shear = base reaction)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, 0.0, height)
    ops.fix(1, 1, 1, 1)
    _rc_fiber_section(1, materials=materials)
    ops.geomTransf("PDelta", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -P, 0.0)
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    if ops.analyze(10) != 0:
        return {"converged": False, "disp": [], "shear": [], "stage": "gravity"}
    ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
    ops.load(2, 1.0, 0.0, 0.0)
    ops.test("NormDispIncr", 1e-6, 1000); ops.algorithm("Newton")
    ops.analysis("Static")
    return _displacement_pushover(2, 1, [1], dU, target)


def run_beamcolumn_dynamic(*, height: float, P: float, materials, nelem: int,
                           self_mass: float, top_mass: float, accel, dt_record: float,
                           scale: float, damping_ratio: float = 0.05,
                           modes: tuple[int, int] = (1, 2), dt: float = 0.01,
                           tol: float = 1e-6) -> dict:
    """Nonlinear seismic time-history of the fiber `forceBeamColumn` column — the reference for the
    lattice (matched material/mass). The column is SUBDIVIDED into `nelem` equal force-based fiber
    elements (nodes at the same heights as the lattice rows) so `self_mass` distributes by tributary
    length exactly like the lattice's lumped self-mass; `top_mass` (= P/g) is added at the free top.
    Fixed base, constant axial `P` (held), 5% modal damping at `modes`, then UniformExcitation of the
    scaled `accel` record in X. Returns the _transient_uniform_excitation dict (carrying "periods")."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nnode = nelem + 1
    seg = height / nelem
    for i in range(1, nnode + 1):
        ops.node(i, 0.0, seg * (i - 1))
    ops.fix(1, 1, 1, 1)
    for i in range(1, nnode + 1):
        trib = seg if 1 < i < nnode else seg / 2.0          # half-tributary at the ends
        m = self_mass * trib / height + (top_mass if i == nnode else 0.0)
        ops.mass(i, m, m, 0.0)                               # translational only (no rotary inertia)

    _rc_fiber_section(1, materials=materials)
    ops.geomTransf("PDelta", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    for e in range(1, nelem + 1):
        # force-based, consistent with the static/pushover references; the column is still
        # subdivided into `nelem` elements so `self_mass` distributes by tributary length like
        # the lattice's lumped self-mass (the subdivision also keeps the per-element span short,
        # which helps the element-level state determination converge under dynamic increments).
        ops.element("forceBeamColumn", e, e, e + 1, 1, 1)

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(nnode, 0.0, -P, 0.0)
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    if ops.analyze(10) != 0:
        return {"converged": False, "stage": "gravity", "t": [], "disp": [], "shear": [],
                "control_node": nnode, "base_nodes": [1]}
    ops.loadConst("-time", 0.0)

    res = _transient_uniform_excitation(accel=accel, dt_record=dt_record, scale=scale, dof=1,
                                        control_node=nnode, control_dof=1, base_nodes=[1],
                                        dt=dt, zeta=damping_ratio, modes=modes, tol=tol)
    return res


def run_beamcolumn_modal(*, height: float, materials, nelem: int, self_mass: float,
                         num_modes: int) -> dict:
    """First `num_modes` modes of the subdivided fiber `forceBeamColumn` cantilever — the modal
    counterpart of `run_beamcolumn_dynamic`, for the calibration mode-shape report (D35).

    Same discretization (nelem segments, nodes at the lattice row heights) and tributary `self_mass`
    as the dynamic run, so it is MASS-CONSISTENT with the lattice/continuum (equal total mass → the
    periods are directly comparable, D16). No top mass and no gravity, so the eigen sees the
    initial-tangent stiffness at the undeformed state — exactly like the lattice/continuum modal in
    `run_modal`. Returns {"periods", "shapes" (per mode, {node: [ux, uy, rot]}), "model"} where
    `model` is a lightweight drawing `Model` (vertical line of nodes + segments) for the visualizer."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nnode = nelem + 1
    seg = height / nelem
    for i in range(1, nnode + 1):
        ops.node(i, 0.0, seg * (i - 1))
    ops.fix(1, 1, 1, 1)
    for i in range(1, nnode + 1):
        trib = seg if 1 < i < nnode else seg / 2.0          # half-tributary at the ends
        m = self_mass * trib / height
        ops.mass(i, m, m, 0.0)                               # translational only (no rotary inertia)

    _rc_fiber_section(1, materials=materials)
    ops.geomTransf("PDelta", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    for e in range(1, nelem + 1):
        ops.element("forceBeamColumn", e, e, e + 1, 1, 1)

    try:
        eigenvalues = ops.eigen(num_modes)                  # genBandArpack: lowest modes (massless rot DOFs ok)
    except Exception:
        eigenvalues = ops.eigen("-fullGenLapack", num_modes)
    periods = [2.0 * math.pi / math.sqrt(lam) if lam > 0.0 else math.inf for lam in eigenvalues]
    shapes = [
        {i: ops.nodeEigenvector(i, mode) for i in range(1, nnode + 1)}
        for mode in range(1, num_modes + 1)
    ]

    model = Model(ndm=2, ndf=2)                              # drawing-only model (lines along the height)
    for i in range(1, nnode + 1):
        model.add_node(i, (0.0, seg * (i - 1)))
    for e in range(1, nelem + 1):
        model.add_element(e, "line", (e, e + 1), kind="concrete")
    return {"periods": list(periods), "shapes": shapes, "model": model}


def run_benchmark_rc_frame(*, dU: float = 0.1, target: float = 15.0, gravity_P: float = 180.0) -> dict:
    """The original OpenSees RC frame pushover (RCFrameGravity -> RCFramePushOver), kip-in.

    The verification reference for the lattice: 1-bay/1-storey, RC fiber `forceBeamColumn`
    columns (15x24 section, Concrete01 core/cover, Steel01, 3+2+3 bars @0.6 in^2), elastic beam,
    gravity P then DisplacementControl pushover of node 3 (DOF 1). Returns the pushover curve
    {"disp", "shear", "converged"} (shear = base reaction sum, the lattice's comparison target).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, 360.0, 0.0); ops.node(3, 0.0, 144.0); ops.node(4, 360.0, 144.0)
    ops.fix(1, 1, 1, 1); ops.fix(2, 1, 1, 1)

    _rc_fiber_section(1)

    ops.geomTransf("PDelta", 1)  # columns
    ops.geomTransf("Linear", 2)  # beam
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element("forceBeamColumn", 1, 1, 3, 1, 1)
    ops.element("forceBeamColumn", 2, 2, 4, 1, 1)
    ops.element("elasticBeamColumn", 3, 3, 4, 360.0, 4030.0, 8640.0, 2)  # A, E, Iz

    # gravity (constant)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, -gravity_P, 0.0); ops.load(4, 0.0, -gravity_P, 0.0)
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    if ops.analyze(10) != 0:
        return {"converged": False, "disp": [], "shear": [], "stage": "gravity"}
    ops.loadConst("-time", 0.0)

    # pushover
    ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
    ops.load(3, 10.0, 0.0, 0.0); ops.load(4, 10.0, 0.0, 0.0)
    ops.test("NormDispIncr", 1e-6, 1000); ops.algorithm("Newton")
    ops.analysis("Static")
    return _displacement_pushover(3, 1, [1, 2], dU, target)


def _rc_beam_fiber_section(sec_tag: int, *, depth: float, width: float, top_area: float,
                           bot_area: float, cover: float = 1.5, nbar: int = 3,
                           core_mat: int = 1, cover_mat: int = 2, steel_mat: int = 3) -> None:
    """RC beam fiber section: `width` wide (z) x `depth` deep (y), `cover` cover ring, with a
    confined core (`core_mat`) inside an unconfined cover ring (`cover_mat`) and `nbar` longitudinal
    bars in a TOP layer (total `top_area`) and a BOTTOM layer (total `bot_area`), both of `steel_mat`.

    Mirrors `_rc_fiber_section` (the column section) but for the thinner beam and a top/bottom bar
    layout. It REFERENCES materials `core_mat`/`cover_mat`/`steel_mat` (defined by the caller, e.g.
    via `_rc_fiber_section` on the column section) — it does NOT redefine them, so a frame can share
    one Concrete02 core / cover / Steel02 trio between its columns and its beam (the lattice does the
    same via the grade -> material mapping)."""
    y1, z1 = depth / 2.0, width / 2.0
    ops.section("Fiber", sec_tag)
    ops.patch("quad", core_mat, 1, 10, -y1 + cover, z1 - cover, -y1 + cover, -z1 + cover,
              y1 - cover, -z1 + cover, y1 - cover, z1 - cover)                 # core
    ops.patch("quad", cover_mat, 1, 1, -y1, z1, -y1, -z1, -y1 + cover, -z1 + cover, -y1 + cover, z1 - cover)
    ops.patch("quad", cover_mat, 1, 1, y1 - cover, z1 - cover, y1 - cover, -z1 + cover, y1, -z1, y1, z1)
    ops.patch("quad", cover_mat, 1, 1, -y1 + cover, z1 - cover, -y1 + cover, z1, y1 - cover, z1, y1 - cover, z1 - cover)
    ops.patch("quad", cover_mat, 1, 1, -y1 + cover, -z1, -y1 + cover, -z1 + cover, y1 - cover, -z1 + cover, y1 - cover, -z1)
    ops.layer("straight", steel_mat, nbar, top_area / nbar, y1 - cover, z1 - cover, y1 - cover, -z1 + cover)  # top
    ops.layer("straight", steel_mat, nbar, bot_area / nbar, -y1 + cover, z1 - cover, -y1 + cover, -z1 + cover)  # bottom


def _beam_section(sec_tag: int, *, depth: float, width: float, top_area: float, bot_area: float,
                  beam_materials=None) -> None:
    """Emit the beam fiber section. By default the beam concrete REUSES the column section's
    core/cover materials (tags 1/2); pass `beam_materials` (a (core, cover) `UniaxialMaterial` pair,
    e.g. Elastic) to give the beam its OWN concrete law (emitted as tags 4/5) while keeping the
    shared Steel02 (tag 3) — used to model the thin beam elastically (the stable default) without
    touching the columns."""
    if beam_materials is None:
        _rc_beam_fiber_section(sec_tag, depth=depth, width=width, top_area=top_area, bot_area=bot_area)
        return
    bcore, bcover = beam_materials
    ops.uniaxialMaterial(bcore.mtype, 4, *bcore.args)
    ops.uniaxialMaterial(bcover.mtype, 5, *bcover.args)
    _rc_beam_fiber_section(sec_tag, depth=depth, width=width, top_area=top_area, bot_area=bot_area,
                           core_mat=4, cover_mat=5, steel_mat=3)


def run_beamcolumn_frame(*, height: float = 144.0, span: float = 144.0, beam_depth: float = 18.0,
                         beam_width: float = 15.0, P: float = 180.0, dU: float = 0.1,
                         target: float = 15.0, materials=None, beam_materials=None,
                         beam_top_area: float = 1.8, beam_bot_area: float = 1.8) -> dict:
    """Fiber `forceBeamColumn` portal-frame pushover — the 1D verification reference for the RC
    lattice frame (the frame analog of `run_beamcolumn_cantilever`). One bay, one storey:
    two RC fiber columns (the 15x24 section, `_rc_fiber_section`, P-Delta) and a thinner RC fiber
    beam (`beam_width` x `beam_depth`, top+bottom bars) sharing the SAME Concrete02 core/cover +
    Steel02 trio (`materials`, exactly as the lattice/continuum share grades). `beam_materials`
    (optional (core, cover) pair) gives the beam its own concrete law (e.g. Elastic, the stable
    default) while keeping the shared steel. Fixed bases, constant gravity `P` at each column top,
    then a DisplacementControl lateral pushover of the top-left joint. Returns the pushover curve
    {"disp", "shear", "converged"} (shear = base reaction sum)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, span, 0.0)            # bases
    ops.node(3, 0.0, height); ops.node(4, span, height)      # beam-column joints
    ops.fix(1, 1, 1, 1); ops.fix(2, 1, 1, 1)

    _rc_fiber_section(1, materials=materials)                 # column section (defines mats 1/2/3)
    _beam_section(2, depth=beam_depth, width=beam_width, top_area=beam_top_area,
                  bot_area=beam_bot_area, beam_materials=beam_materials)      # beam section
    ops.geomTransf("PDelta", 1)                               # columns (P-Delta under gravity)
    ops.geomTransf("Linear", 2)                               # beam
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.beamIntegration("Lobatto", 2, 2, 5)
    ops.element("forceBeamColumn", 1, 1, 3, 1, 1)            # left column
    ops.element("forceBeamColumn", 2, 2, 4, 1, 1)            # right column
    ops.element("forceBeamColumn", 3, 3, 4, 2, 2)            # beam

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, -P, 0.0); ops.load(4, 0.0, -P, 0.0)     # gravity at the column tops
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    if ops.analyze(10) != 0:
        return {"converged": False, "disp": [], "shear": [], "stage": "gravity"}
    ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
    ops.load(3, 1.0, 0.0, 0.0); ops.load(4, 1.0, 0.0, 0.0)   # lateral reference pattern (+X)
    ops.test("NormDispIncr", 1e-6, 1000); ops.algorithm("Newton")
    ops.analysis("Static")
    return _displacement_pushover(3, 1, [1, 2], dU, target)


def run_beamcolumn_frame_dynamic(*, height: float, span: float, beam_depth: float, beam_width: float,
                                 P: float, materials, ncol: int, nbeam: int, self_mass_col: float,
                                 self_mass_beam: float, top_mass: float, accel, dt_record: float,
                                 scale: float, beam_materials=None, beam_top_area: float = 1.8,
                                 beam_bot_area: float = 1.8, damping_ratio: float = 0.05,
                                 modes: tuple[int, int] = (1, 2), dt: float = 0.01,
                                 tol: float = 1e-6) -> dict:
    """Nonlinear seismic time-history of the fiber `forceBeamColumn` portal frame — the dynamic
    reference for the lattice frame (the frame analog of `run_beamcolumn_dynamic`). Each column is
    SUBDIVIDED into `ncol` force-based fiber elements and the beam into `nbeam`, so `self_mass_col`
    (per column) and `self_mass_beam` distribute by tributary length like the lattice's lumped
    self-mass; the axial-load tributary mass `top_mass` (= P/g, per column) is lumped at each
    beam-column joint. `beam_materials` (optional (core, cover) pair) gives the beam its own concrete
    law (e.g. Elastic) while keeping the shared steel. Fixed bases, constant gravity `P` per column
    (held), 5% modal damping at `modes`, then a UniformExcitation of the scaled `accel` record in X.
    Returns the `_transient_uniform_excitation` dict (carrying the modal "periods")."""
    from collections import defaultdict

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    col_seg, beam_seg = height / ncol, span / nbeam

    for i in range(ncol + 1):                                 # left column nodes 1 .. ncol+1
        ops.node(1 + i, 0.0, col_seg * i)
    joint_L = ncol + 1
    rbase = ncol + 2
    for i in range(ncol + 1):                                 # right column nodes rbase .. rbase+ncol
        ops.node(rbase + i, span, col_seg * i)
    joint_R = rbase + ncol
    beam_nodes = [joint_L]                                    # beam: joint_L, interior, joint_R
    nid = joint_R + 1
    for k in range(1, nbeam):
        ops.node(nid, beam_seg * k, height)
        beam_nodes.append(nid)
        nid += 1
    beam_nodes.append(joint_R)
    ops.fix(1, 1, 1, 1); ops.fix(rbase, 1, 1, 1)

    # accumulate nodal mass (a joint receives column + beam tributary + the axial top mass), emit once
    m: dict[int, float] = defaultdict(float)
    for col_base in (1, rbase):
        for i in range(ncol + 1):
            trib = col_seg if 0 < i < ncol else col_seg / 2.0
            m[col_base + i] += self_mass_col * trib / height
    for idx, n in enumerate(beam_nodes):
        trib = beam_seg if 0 < idx < len(beam_nodes) - 1 else beam_seg / 2.0
        m[n] += self_mass_beam * trib / span
    m[joint_L] += top_mass; m[joint_R] += top_mass
    for n, mv in m.items():
        ops.mass(n, mv, mv, 0.0)                              # translational only (no rotary inertia)

    _rc_fiber_section(1, materials=materials)
    _beam_section(2, depth=beam_depth, width=beam_width, top_area=beam_top_area,
                  bot_area=beam_bot_area, beam_materials=beam_materials)
    ops.geomTransf("PDelta", 1); ops.geomTransf("Linear", 2)
    ops.beamIntegration("Lobatto", 1, 1, 5); ops.beamIntegration("Lobatto", 2, 2, 5)
    eid = 1
    for col_base in (1, rbase):                               # force-based columns (subdivided)
        for i in range(ncol):
            ops.element("forceBeamColumn", eid, col_base + i, col_base + i + 1, 1, 1)
            eid += 1
    for a, b in zip(beam_nodes, beam_nodes[1:]):             # force-based beam (subdivided)
        ops.element("forceBeamColumn", eid, a, b, 2, 2)
        eid += 1

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(joint_L, 0.0, -P, 0.0); ops.load(joint_R, 0.0, -P, 0.0)
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    if ops.analyze(10) != 0:
        return {"converged": False, "stage": "gravity", "t": [], "disp": [], "shear": [],
                "control_node": joint_L, "base_nodes": [1, rbase]}
    ops.loadConst("-time", 0.0)

    return _transient_uniform_excitation(accel=accel, dt_record=dt_record, scale=scale, dof=1,
                                         control_node=joint_L, control_dof=1, base_nodes=[1, rbase],
                                         dt=dt, zeta=damping_ratio, modes=modes, tol=tol)


def run_beamcolumn_frame_modal(*, height: float, span: float, beam_depth: float, beam_width: float,
                               materials, ncol: int, nbeam: int, self_mass_col: float,
                               self_mass_beam: float, num_modes: int, beam_materials=None,
                               beam_top_area: float = 1.8, beam_bot_area: float = 1.8) -> dict:
    """First `num_modes` modes of the subdivided fiber portal frame — the modal counterpart of
    `run_beamcolumn_frame_dynamic` for the calibration mode-shape report (D35), and the frame analog of
    `run_beamcolumn_modal`. Same discretization (each column into `ncol`, the beam into `nbeam`) and
    tributary `self_mass_col`/`self_mass_beam`, so it is mass-consistent with the lattice/continuum
    frame; no top mass and no gravity → the eigen sees the initial-tangent stiffness. Returns
    {"periods", "shapes" (per mode, {node: [ux, uy, rot]}), "model"} where `model` is a lightweight
    drawing `Model` (the two columns + the beam as line segments) for the visualizer."""
    from collections import defaultdict

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    col_seg, beam_seg = height / ncol, span / nbeam
    for i in range(ncol + 1):                                 # left column nodes 1 .. ncol+1
        ops.node(1 + i, 0.0, col_seg * i)
    joint_L = ncol + 1
    rbase = ncol + 2
    for i in range(ncol + 1):                                 # right column nodes rbase .. rbase+ncol
        ops.node(rbase + i, span, col_seg * i)
    joint_R = rbase + ncol
    beam_nodes = [joint_L]                                    # beam: joint_L, interior, joint_R
    nid = joint_R + 1
    for k in range(1, nbeam):
        ops.node(nid, beam_seg * k, height)
        beam_nodes.append(nid)
        nid += 1
    beam_nodes.append(joint_R)
    ops.fix(1, 1, 1, 1); ops.fix(rbase, 1, 1, 1)

    m: dict[int, float] = defaultdict(float)                  # tributary self-mass (no top mass)
    for col_base in (1, rbase):
        for i in range(ncol + 1):
            trib = col_seg if 0 < i < ncol else col_seg / 2.0
            m[col_base + i] += self_mass_col * trib / height
    for idx, n in enumerate(beam_nodes):
        trib = beam_seg if 0 < idx < len(beam_nodes) - 1 else beam_seg / 2.0
        m[n] += self_mass_beam * trib / span
    for n, mv in m.items():
        ops.mass(n, mv, mv, 0.0)                              # translational only (no rotary inertia)

    _rc_fiber_section(1, materials=materials)
    _beam_section(2, depth=beam_depth, width=beam_width, top_area=beam_top_area,
                  bot_area=beam_bot_area, beam_materials=beam_materials)
    ops.geomTransf("PDelta", 1); ops.geomTransf("Linear", 2)
    ops.beamIntegration("Lobatto", 1, 1, 5); ops.beamIntegration("Lobatto", 2, 2, 5)
    eid = 1
    segments: list[tuple[int, int]] = []
    for col_base in (1, rbase):                               # columns (subdivided)
        for i in range(ncol):
            ops.element("forceBeamColumn", eid, col_base + i, col_base + i + 1, 1, 1)
            segments.append((col_base + i, col_base + i + 1)); eid += 1
    for a, b in zip(beam_nodes, beam_nodes[1:]):             # beam (subdivided)
        ops.element("forceBeamColumn", eid, a, b, 2, 2)
        segments.append((a, b)); eid += 1

    try:
        eigenvalues = ops.eigen(num_modes)
    except Exception:
        eigenvalues = ops.eigen("-fullGenLapack", num_modes)
    periods = [2.0 * math.pi / math.sqrt(lam) if lam > 0.0 else math.inf for lam in eigenvalues]
    all_nodes = (list(range(1, joint_L + 1)) + list(range(rbase, joint_R + 1)) + beam_nodes[1:-1])
    shapes = [
        {n: ops.nodeEigenvector(n, mode) for n in all_nodes}
        for mode in range(1, num_modes + 1)
    ]

    model = Model(ndm=2, ndf=2)                              # drawing-only model (columns + beam lines)
    for i in range(ncol + 1):
        model.add_node(1 + i, (0.0, col_seg * i))
    for i in range(ncol + 1):
        model.add_node(rbase + i, (span, col_seg * i))
    for k in range(1, nbeam):
        model.add_node(joint_R + k, (beam_seg * k, height))
    for j, (a, b) in enumerate(segments, start=1):
        model.add_element(j, "line", (a, b), kind="concrete")
    return {"periods": list(periods), "shapes": shapes, "model": model}
