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
    algorithm: "tuple[str, ...]" = ("Newton",),
    element_groups=None,
    capture: bool = False,
    node_history=None,
    node_history_every: int = 1,
    equal_dof=None,
) -> dict:
    """Gravity (constant) → DisplacementControl pushover, recording base shear (D18/D19).

    Sequence (mirrors the OpenSees RCFrameGravity → RCFramePushOver benchmark):
      1. apply `gravity_loads` (or `model.loads`) as a constant pattern;
      2. add `lateral_loads` (a list of `Load`) as the reference lateral pattern;
      3. step `control_node`'s `control_dof` by `dU` (signed by `target`) up to `target`,
         summing horizontal base reactions into the base shear at each converged step.

    Base shear is `-sum(reaction[control_dof])` over `base_nodes` (defaults to the supported
    nodes) — the structure's resistance, positive for a positive push. A failed step is retried
    with finer sub-steps and a stronger algorithm (the `algorithm` primary → KrylovNewton →
    NewtonLineSearch) before giving up — the Stage-2 nonlinear lattice (Concrete02 softening) needs
    this. `tol`/`max_iter` set the NormDispIncr test (1e-6 is the practical tolerance for the
    softening lattice).

    `algorithm` is the PRIMARY solution algorithm (an `ops.algorithm(*algorithm)` arg tuple), used
    for normal stepping and as the first rung of the retry ladder; it defaults to `("Newton",)`. For
    a plastic/softening pushover pass `("ModifiedNewton", "-initial")` — iterating on the constant
    INITIAL (elastic) stiffness never re-forms the singular/negative cracked tangent that trips full
    Newton, so it is the robust default there.

    `element_groups` (optional) is a force-decomposition probe for diagnostics: a dict mapping a
    label to a list of `(element_id, dof, coef)`. At each converged step it records, per label, the
    sum of `eleForce(element_id)[dof] * coef` — a GLOBAL nodal force component the element exerts
    (so the corotational/P-Δ geometry is already baked in, unlike a precomputed direction cosine).
    With `dof` = the horizontal index at an element's node-on-one-side-of-a-cut and `coef`=1, the
    per-label sums reconcile to the base shear; pairing the vertical index with `coef`=that node's x
    gives the overturning-moment share. Lets a caller attribute base shear/overturning to element
    categories (vertical vs diagonal struts, concrete vs rebar) without any ops.* calls of its own.

    `capture` and `equal_dof` mirror `run_pushover_dynamic` so the two solvers are interchangeable
    for a caller: `capture=True` also returns the full nodal displacement field at the last step
    ("disps_final") and at peak |base shear| ("disps_peak"), which a damage figure needs;
    `equal_dof` is a list of `(retained, constrained, dof)` emitted as `ops.equalDOF` (a rigid
    loading platen, say), applied before the gravity stage so it holds throughout.

    `node_history=(node_ids, dof)` follows those nodes' displacement component `dof` (1-based)
    through the run, sampled every `node_history_every` recorded steps. `dof` may be one 1-based component or a SEQUENCE of them; several dofs are recorded dof-major,
    so `values[k*len(nodes):(k+1)*len(nodes)]` is the k-th dof and a one-dof probe is unchanged (D68).
    It is the gauge-length probe
    a STRAIN PROFILE needs: two aligned node rows a fixed distance apart give
    `(u_top - u_bottom)/gauge` — exactly what a vertical LVDT pair on a wall face measures. Returned
    as "node_history" = {"nodes", "dof", "index", "values"}, where `index` holds positions into
    `disp`/`shear`, so each profile can be tied to the drift it was taken at.

    Returns {"ok", "converged", "disp": [...], "shear": [...], "control_node", "base_nodes"}, plus
    "groups": {label: [...]} when `element_groups` is given, "disps_final"/"disps_peak" when
    `capture` is set, and "node_history" when `node_history` is given.
    """
    ops.wipe()
    build(model)
    base = list(base_nodes) if base_nodes is not None else [s.node for s in model.supports]

    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")

    # Multi-point constraints before the gravity stage, so they hold for the whole analysis.
    for retained, constrained, cdof in (equal_dof or ()):
        ops.equalDOF(int(retained), int(constrained), int(cdof))

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
    ops.algorithm(*algorithm)  # primary solver (default Newton; ModifiedNewton -initial for softening)
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

    snap = {"peak": None, "peak_abs": -1.0}
    probe = _node_probe(node_history)
    nh = ({"nodes": probe[0], "dof": probe[1][0], "dofs": list(probe[1]),
           "index": [], "values": []} if probe else None)

    def record() -> None:
        disp.append(control_disp())
        shear.append(base_shear())
        for label, members in (element_groups or {}).items():
            groups[label].append(sum(ops.eleForce(eid)[dof] * coef for eid, dof, coef in members))
        if capture and abs(shear[-1]) > snap["peak_abs"]:
            snap["peak_abs"] = abs(shear[-1])
            snap["peak"] = _nodal_disps(model)
        if nh is not None and (len(disp) - 1) % max(1, node_history_every) == 0:
            nh["index"].append(len(disp) - 1)
            nh["values"].append(_probe_sample(nh))

    record()
    converged = True
    goal = abs(target)
    while sign * control_disp() < goal - 1e-9 * (goal + 1.0):
        if ops.analyze(1) == 0:
            record()
            continue
        # retry the increment with finer sub-steps and stronger algorithms (softening, D20)
        sub_ok = False
        for spec in (algorithm, ("KrylovNewton",), ("NewtonLineSearch", "-type", "Bisection")):
            ops.algorithm(*spec)
            for nsub in (5, 20, 50):
                ops.integrator("DisplacementControl", control_node, control_dof, du / nsub)
                if ops.analyze(nsub) == 0:
                    record()
                    sub_ok = True
                    break
            if sub_ok:
                break
        ops.algorithm(*algorithm)
        ops.integrator("DisplacementControl", control_node, control_dof, du)  # restore full step
        if not sub_ok:
            converged = False
            break

    result = {"ok": 0 if converged else -1, "converged": converged, "disp": disp,
              "shear": shear, "control_node": control_node, "base_nodes": base}
    if element_groups is not None:
        result["groups"] = groups
    if capture:
        result["disps_final"] = _nodal_disps(model)
        result["disps_peak"] = snap["peak"]
    if nh is not None:
        result["node_history"] = nh
    return result


def cyclic_protocol(peaks, *, cycles_per_level: int = 2, start_positive: bool = True) -> list[float]:
    """Expand target drift/displacement amplitudes into a reversed-cyclic displacement history.

    `peaks` are the amplitudes (same unit as the control DOF); each is repeated
    `cycles_per_level` times as a full push-pull cycle returning to zero, e.g. peaks=[1, 2] with
    2 cycles gives [+1,-1,+1,-1, +2,-2,+2,-2, 0]. `cycles_per_level` may also be a per-level list
    (the ACI 374.2R protocol used for these wall tests runs three cycles up to yield and two
    after). The history always ends back at zero.
    """
    counts = ([cycles_per_level] * len(peaks) if isinstance(cycles_per_level, int)
              else list(cycles_per_level))
    if len(counts) != len(peaks):
        raise ValueError(f"cycles_per_level has {len(counts)} entries for {len(peaks)} peaks")
    sign = 1.0 if start_positive else -1.0
    history: list[float] = []
    for peak, n in zip(peaks, counts):
        for _ in range(int(n)):
            history.extend((sign * abs(peak), -sign * abs(peak)))
    history.append(0.0)
    return history


def run_cyclic(
    model: Model,
    *,
    lateral_loads,
    control_node: int,
    control_dof: int,
    history,
    dU: float,
    gravity_loads=None,
    gravity_steps: int = 10,
    base_nodes=None,
    tol: float = 1e-5,
    max_iter: int = 100,
    algorithm: "tuple[str, ...]" = ("ModifiedNewton", "-initial"),
    node_history=None,
    node_history_every: int = 1,
    wipe: bool = True,
) -> dict:
    """Gravity (constant) → reversed-cyclic DisplacementControl, recording base shear (D48).

    The cyclic sibling of `run_pushover`: same gravity-then-lateral-pattern setup, but instead of
    marching to a single target it walks the control DOF to each displacement in `history` in turn,
    reversing the increment sign as needed. Build `history` with `cyclic_protocol`.

    `algorithm` defaults to `("ModifiedNewton", "-initial")` — a reversed-cyclic run on a cracking
    lattice re-forms a softening (possibly negative) tangent at every reversal, and iterating on the
    constant initial stiffness is what survives that (the D44 finding, and doubly true here).
    Failed increments fall back to finer sub-steps and stronger algorithms exactly as in
    `run_pushover`; on exhaustion the run STOPS and returns what it traced with `converged=False`,
    so a partial hysteresis is still usable and honestly labelled.

    `node_history=(node_ids, dof)` follows those nodes' displacement component `dof` (1-based),
    sampled every `node_history_every` steps AND at every reversal. `dof` may be one 1-based component or a SEQUENCE of them; several dofs are recorded dof-major,
    so `values[k*len(nodes):(k+1)*len(nodes)]` is the k-th dof and a one-dof probe is unchanged (D68).
    It is the gauge-length probe a strain
    profile needs (D63). Reversals are always captured because the loop tips are the instants the
    profile is wanted at, and a stride alone would miss them by up to a full stride of drive.

    Returns {"ok", "converged", "disp", "shear", "control_node", "base_nodes", "reached"} where
    `reached` is the number of history targets completed, plus "node_history" =
    {"nodes", "dof", "index", "values"} when `node_history` is given (`index` holds positions into
    `disp`/`shear`).
    """
    # WIPE AND BUILD, like every other runner here — `build()` emits into the CURRENT domain and
    # does not clear it (its docstring says "after ops.wipe"), so without the wipe a second
    # run_cyclic in one process appends to the previous model, OpenSees rejects the repeated node
    # tags, and the solve fails somewhere else entirely. This was the only model-building runner
    # missing it.
    #
    # `wipe=False` HANDS THE DOMAIN TO THE CALLER, and skips the build with it. The case that
    # matters is gravity applied SEPARATELY: `run_gravity` is standalone (it wipes, builds, and
    # leaves the model built and loaded), so wiping here would discard that state while rebuilding
    # on top of it would collide on node tags — which is why that pairing never worked. With
    # wipe=False neither happens and the domain is analysed as it stands; pass `gravity_loads=None`
    # too, or gravity is applied a second time.
    #
    # The ORDINARY path remains this runner's own `gravity_loads`, applied after the build and held
    # constant by the `loadConst` below.
    if wipe:
        ops.wipe()
        build(model)
    elif not ops.getNodeTags():
        raise RuntimeError(
            "run_cyclic(wipe=False) analyses the domain the CALLER has already built (e.g. via "
            "run_gravity), but the domain is empty. Either build it first or leave wipe=True.")
    base = list(base_nodes) if base_nodes is not None else [s.node for s in model.supports]

    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")

    # Under `wipe=False` the caller owns the domain AND its loads, so gravity is NOT applied here:
    # `_gravity_loads(model, None)` falls back to the MODEL's own loads, which after a separate
    # `run_gravity` would apply gravity a second time and collide on the time-series tag. What is
    # still needed is to freeze whatever the caller applied, so the displacement control below
    # starts from a held gravity state — that is the `loadConst` this branch keeps.
    if not wipe:
        ops.loadConst("-time", 0.0)
    else:
        grav = _gravity_loads(model, gravity_loads)
        if grav:
            if _apply_gravity(model, grav, gravity_steps, tol) != 0:
                return {"ok": -1, "converged": False, "stage": "gravity", "disp": [], "shear": [],
                        "control_node": control_node, "base_nodes": base, "reached": 0}
            ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for ld in lateral_loads:
        ops.load(ld.node, *ld.values)

    ops.test("NormDispIncr", tol, max_iter)
    ops.algorithm(*algorithm)
    ops.analysis("Static")

    def control_disp() -> float:
        return ops.nodeDisp(control_node, control_dof)

    def base_shear() -> float:
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base)

    probe = _node_probe(node_history)
    nh = ({"nodes": probe[0], "dof": probe[1][0], "dofs": list(probe[1]),
           "index": [], "values": []} if probe else None)
    disp: list[float] = []
    shear: list[float] = []

    def record() -> None:
        """Append one converged step. Also samples the node probe — on the stride, and always at a
        REVERSAL, since the loop tips are the instants a strain profile is wanted at and a stride
        alone would miss them by up to `node_history_every` steps of drive."""
        disp.append(control_disp())
        shear.append(base_shear())
        if nh is None:
            return
        turned = (len(disp) >= 3
                  and (disp[-1] - disp[-2]) * (disp[-2] - disp[-3]) < 0.0)
        if turned or (len(disp) - 1) % max(1, node_history_every) == 0:
            nh["index"].append(len(disp) - 1)
            nh["values"].append(_probe_sample(nh))

    record()
    converged, reached = True, 0
    for goal in history:
        sign = 1.0 if goal >= control_disp() else -1.0
        du = sign * abs(dU)
        ops.integrator("DisplacementControl", control_node, control_dof, du)
        while sign * (goal - control_disp()) > 1e-9 * (abs(goal) + 1.0):
            if ops.analyze(1) == 0:
                record()
                continue
            sub_ok = False
            for spec in (algorithm, ("KrylovNewton",), ("NewtonLineSearch", "-type", "Bisection")):
                ops.algorithm(*spec)
                for nsub in (5, 20, 50):
                    ops.integrator("DisplacementControl", control_node, control_dof, du / nsub)
                    if ops.analyze(nsub) == 0:
                        record()
                        sub_ok = True
                        break
                if sub_ok:
                    break
            ops.algorithm(*algorithm)
            ops.integrator("DisplacementControl", control_node, control_dof, du)
            if not sub_ok:
                converged = False
                break
        if not converged:
            break
        reached += 1

    out = {"ok": 0 if converged else -1, "converged": converged, "disp": disp, "shear": shear,
           "control_node": control_node, "base_nodes": base, "reached": reached}
    if nh is not None:
        out["node_history"] = nh
    return out


def run_pushover_arclength(
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
    tol: float = 1e-6,
    max_iter: int = 50,
    max_steps: int = 800,
    alpha: float = 0.0,
    stop_shear_frac: float = 0.02,
) -> dict:
    """Gravity (constant) → cylindrical ARC-LENGTH pushover, recording base shear.

    Traces the full force–deformation response INCLUDING snapback — the descending branch that
    ordinary DisplacementControl cannot follow, because the control displacement itself must
    *decrease* while the load drops (a localizing crack band opening while the rest of the specimen
    unloads elastically). Arc-length constrains the displacement-increment NORM and lets the load
    factor float, so it traverses limit points AND snapback. `alpha=0` gives Crisfield's cylindrical
    form (the load term dropped from the constraint) — the robust choice for softening.

    Sequence (mirrors `run_pushover`):
      1. apply `gravity_loads` (or `model.loads`) as a constant, held pattern;
      2. add `lateral_loads` (a list of `Load`) as the reference pattern the arc-length scales;
      3. take ONE DisplacementControl *probe* step of `dU` on `control_node`/`control_dof` — this
         fixes a well-scaled arc length (= the resulting displacement-increment norm) and the initial
         loading direction (the probe is still on the ascending branch, so the tangent is positive
         definite and arc-length then continues forward);
      4. march with `ArcLength` up to `max_steps`, recording (control disp, base shear) per step.

    Stops when `|control disp|` reaches `target`, the base shear falls below `stop_shear_frac` of the
    running peak (crack essentially fully open — a clean end for a softening coupon), or a step fails
    after arc-length reduction. On a failed step the arc length is cut (and stronger algorithms
    tried) before giving up. Base shear = `-Σ reaction[control_dof]` over `base_nodes` (defaults to
    the supported nodes). Returns {"ok","converged","disp","shear","control_node","base_nodes"}.
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
            return {"ok": -1, "converged": False, "stage": "gravity", "disp": [], "shear": [],
                    "control_node": control_node, "base_nodes": base}
        ops.loadConst("-time", 0.0)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for ld in lateral_loads:
        ops.load(ld.node, *ld.values)

    ops.test("NormDispIncr", tol, max_iter)
    ops.algorithm("Newton")
    ops.analysis("Static")

    def cdisp() -> float:
        return ops.nodeDisp(control_node, control_dof)

    def base_shear() -> float:
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base)

    def disp_norm() -> float:
        """L2 norm of the full nodal-displacement vector (= the increment norm for a from-zero step)."""
        total = 0.0
        for nid in model.nodes:
            for d in range(1, model.ndf + 1):
                u = ops.nodeDisp(nid, d)
                total += u * u
        return math.sqrt(total)

    disp = [cdisp()]
    shear = [base_shear()]

    # probe step (DisplacementControl) — sets the arc-length scale and the loading direction
    sign = 1.0 if target >= 0.0 else -1.0
    ops.integrator("DisplacementControl", control_node, control_dof, sign * abs(dU))
    if ops.analyze(1) != 0:
        return {"ok": -1, "converged": False, "stage": "probe", "disp": disp, "shear": shear,
                "control_node": control_node, "base_nodes": base}
    disp.append(cdisp())
    shear.append(base_shear())

    s = disp_norm()   # arc length = the probe step's displacement-increment norm (from zero)
    ops.integrator("ArcLength", s, alpha)

    converged = True
    goal = abs(target)
    for _ in range(max_steps):
        if ops.analyze(1) != 0:
            sub_ok = False
            for algo, args in (("Newton", ()), ("KrylovNewton", ()),
                               ("NewtonLineSearch", ("-type", "Bisection"))):
                ops.algorithm(algo, *args)
                for factor in (0.5, 0.25, 0.1):
                    ops.integrator("ArcLength", s * factor, alpha)
                    if ops.analyze(1) == 0:
                        sub_ok = True
                        break
                if sub_ok:
                    break
            ops.algorithm("Newton")
            ops.integrator("ArcLength", s, alpha)   # restore the full arc length
            if not sub_ok:
                converged = False
                break
        disp.append(cdisp())
        shear.append(base_shear())
        if sign * cdisp() >= goal - 1e-12:
            break
        peak = max(abs(v) for v in shear)
        if peak > 0.0 and abs(shear[-1]) < stop_shear_frac * peak and sign * cdisp() > 0.0:
            break

    return {"ok": 0 if converged else -1, "converged": converged, "disp": disp, "shear": shear,
            "control_node": control_node, "base_nodes": base}


def _nodal_disps(model: Model) -> dict:
    """Snapshot of every node's displacement vector, {node_id: [u1..u_ndf]} (viz's `disp` format)."""
    return {nid: [ops.nodeDisp(nid, d) for d in range(1, model.ndf + 1)] for nid in model.nodes}


def nodal_displacements(model: Model) -> dict:
    """The current displacement field, `{node_id: [u1..u_ndf]}` — public form of the snapshot.

    `capture=True` on the runners returns this at two instants (peak and final), but only once the
    run RETURNS. A run that is killed, diverges or is aborted mid-way leaves nothing behind, which
    is exactly when a damage picture is most wanted. Calling this from a `progress` callback lets a
    caller dump the field periodically, so an interrupted run still has a last known state.

    Feed the result straight to `viz.strut_strains` / `viz.figure_damage`; it needs no recorder.
    Costs one `ops.nodeDisp` per node per DOF, so call it on a progress tick, not every step.
    """
    return _nodal_disps(model)


def _node_probe(node_history):
    """Normalize the `node_history` argument to `(nodes, dofs)` or None.

    `node_history` is a `(node_ids, dof)` pair naming a handful of nodes whose displacement is to be
    followed through the analysis — the gauge-length probe a strain profile needs (D63). Unlike
    `capture`, which snapshots the WHOLE field at two instants, this follows a few nodes through the
    whole history, which is what a physical LVDT does.

    `dof` is either ONE 1-based component (the D63 vertical gauge: one dof, one value per node) or a
    SEQUENCE of them (D68: a diagonal gauge needs both u_x and u_y at each corner, because a chord
    length between two moving nodes is not a function of either component alone). Multiple dofs are
    stored dof-MAJOR — see `_probe_sample` — so a one-dof probe records exactly what it always did
    and every existing reader keeps working unchanged.
    """
    if node_history is None:
        return None
    nodes, dof = node_history
    dofs = (int(dof),) if isinstance(dof, (int, float)) else tuple(int(d) for d in dof)
    if not dofs:
        raise ValueError("node_history: at least one dof is required")
    return [int(n) for n in nodes], dofs


def _probe_sample(nh) -> list:
    """One probe sample: dof-MAJOR blocks, so `values[k*len(nodes):(k+1)*len(nodes)]` is dof k.

    Block-major rather than interleaved on purpose: with a single dof the record is byte-identical
    to the pre-D68 one-dof format, so `specimen.segment_strain`-style row slicing is unaffected, and
    with several dofs each component is still one contiguous slice.
    """
    return [ops.nodeDisp(n, d) for d in nh["dofs"] for n in nh["nodes"]]


# Integrators that march WITHOUT forming or factorizing a tangent. They are conditionally stable —
# dt must sit below 2/w_max, set by the stiffest/lightest element, not by T1 — but each step costs a
# force recovery instead of a factorized Newton solve, which is the trade that makes them worth having
# on a model that will not converge implicitly (D74).
EXPLICIT_INTEGRATORS = frozenset({
    "CentralDifference", "CentralDifferenceAlternative", "CentralDifferenceNoDamping",
    "NewmarkExplicit", "HHTExplicit", "AlphaOS", "AlphaOSGeneralized", "KRAlphaExplicit",
})


def is_explicit(integrator) -> bool:
    return bool(integrator) and str(integrator[0]) in EXPLICIT_INTEGRATORS


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
    rate: float | None = None,
    steps_per_period: int = 40,
    damping_ratio: float = 0.6,
    tol: float = 1e-5,
    max_iter: int = 50,
    capture: bool = False,
    quasi_static: bool = False,
    node_history=None,
    node_history_every: int = 1,
    progress=None,
    progress_every: int = 200,
    integrator=("Newmark", 0.5, 0.25),
    element_groups=None,
    equal_dof=None,
    algorithm: "tuple[str, ...]" = ("Newton",),
) -> dict:
    """Dynamic-relaxation pushover: a quasi-static TRANSIENT solve that rides through limit
    points / local snap-backs the static Newton can't (D22, user-selected).

    After the constant-gravity stage, the lateral displacement is *imposed* (ramped) on
    `drive_nodes` (default [control_node]) via single-point constraints, and the model is marched
    with Newmark + heavy Rayleigh damping slowly enough (`periods_to_target` fundamental periods
    to reach `target`) that inertia stays negligible — so the recorded base shear vs control
    displacement is the quasi-static pushover, but the mass/damping regularize the instability.

    `target` may be NEGATIVE to drive the control DOF backwards — a COMPRESSION pushover, say.
    The drive speed is always a positive magnitude; its sign rides on the imposed-displacement
    constraint, and the stop test compares magnitudes (matching `run_pushover`, D53).

    `capture=True` additionally returns the FULL nodal displacement field at the last step
    ("disps_final") and at the peak |base shear| ("disps_peak"), which is what a damage/crack-pattern
    figure needs. Off by default: it copies the displacement vector on every new peak, which is
    wasted work for a run that only wants the curve.

    `element_groups` is the force-decomposition probe, with the same contract as `run_pushover`: a
    dict mapping a label to a list of `(element_id, dof, coef)`. At each step it records, per label,
    `sum(eleForce(element_id)[dof] * coef)` — a GLOBAL nodal force component the element exerts, so
    corotational geometry is already baked in. Summing the forces the struts crossing one cut exert
    on the nodes below it attributes the carried load to element categories (vertical vs diagonal
    struts, say), which is how a "who is actually carrying this?" question gets answered.

    `algorithm` is the PRIMARY solution algorithm (an `ops.algorithm(*algorithm)` arg tuple),
    matching `run_pushover`. It is used for normal stepping and as the first rung of the retry
    ladder (primary -> KrylovNewton -> NewtonLineSearch, at a tenth of the step). Plain `Newton` is
    the default and is usually right here: unlike the static case, the mass and damping keep the
    effective tangent positive-definite through cracking, so there is no need to iterate on the
    initial stiffness.

    `quasi_static=True` MEASURES whether the drive is slow enough, instead of assuming it. Global
    equilibrium in the drive direction reads `S_base + S_drive = -sum_free (M.a + C.v)`, where each
    S is the `-sum(reaction)` over that node set: the internal forces cancel globally, so whatever
    the two constrained sets do NOT balance is exactly the inertial + damping force riding in the
    recorded base shear. It costs one extra reaction sum over the drive nodes per step and returns
    that residual as "dynamic"; `max|dynamic| / max|shear|` is the contamination of the result.
    It assumes `base_nodes` and `drive_nodes` between them cover EVERY constrained DOF in the drive
    direction (true whenever the supports are the base row) — any other restraint would leak into
    the residual and read as contamination that is not there.

    `equal_dof` ties degrees of freedom together before the analysis: a list of
    `(retained_node, constrained_node, dof)` emitted as `ops.equalDOF`. Use it for a rigid loading
    platen — tie a whole face's vertical DOF to one node and drive only that node, so the face is
    forced to stay flat instead of each node being told independently where to go. Requires the
    `Transformation` constraint handler, which this runner already sets.

    `node_history=(node_ids, dof)` follows those nodes' displacement component `dof` (1-based),
    `integrator` is the time integrator as an `ops.integrator(*args)` tuple. The default
    `("Newmark", 0.5, 0.25)` is average-acceleration: unconditionally stable for a LINEAR system and
    exactly energy-conserving, which means it sustains a spurious high-frequency mode indefinitely
    rather than bleeding it away. `run_cyclic_dynamic` instead uses `("HHT", 0.7)`, whose alpha < 1
    adds numerical damping aimed squarely at the high-frequency response a cracking lattice
    generates. On a model whose struts are very stiff relative to their tributary mass — so that
    local node-pair modes sit far above T1 — the two integrators can behave completely differently
    through cracking, and a screen run on one certifies nothing about the other (D69).

    `progress(step, nsteps, disp, shear)` is called every `progress_every` steps, mirroring
    `run_cyclic_dynamic`. Without it a dynamic pushover is silent for its whole duration, which on a
    redirected log means no way to tell a slow run from a stuck one.

    `node_history=(node_ids, dof)` is
    sampled every `node_history_every` recorded steps — the gauge-length probe a STRAIN PROFILE
    needs (D63): two aligned node rows a fixed distance apart give `(u_top - u_bottom)/gauge`,
    which is what a vertical LVDT pair on a wall face measures. Unlike `capture`, which snapshots
    the whole field at two instants, this follows a few nodes through the whole run.

    Returns {"converged", "disp", "shear", "control_node", "base_nodes", "rate", "T1"}, plus
    "disps_final"/"disps_peak" when `capture` is set, "dynamic" when `quasi_static` is set,
    "node_history" = {"nodes", "dof", "index", "values"} when `node_history` is given (`index`
    holds positions into `disp`/`shear`), and "groups" when `element_groups` is given.
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

    # Drop the gravity stage's Static analysis first: while it exists OpenSees REFUSES to set a
    # transient integrator and silently falls back to a default, so the Newmark below is ignored.
    # wipeAnalysis clears only the analysis aggregation, so the solver must be re-declared (D49).
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")

    # Multi-point constraints go on AFTER wipeAnalysis (which clears the analysis aggregation, not
    # the domain) and BEFORE the eigen solve, so the tied DOFs are reflected in T1 and the damping.
    for retained, constrained, cdof in (equal_dof or ()):
        ops.equalDOF(int(retained), int(constrained), int(cdof))

    # fundamental frequency sets the (quasi-static) loading rate, step, and Rayleigh damping
    try:                                     # the default solver is far faster than fullGenLapack
        lam1 = ops.eigen(1)[0]
    except Exception:
        lam1 = ops.eigen("-fullGenLapack", 1)[0]
    w1 = math.sqrt(lam1) if lam1 > 0 else 1.0
    T1 = 2.0 * math.pi / w1
    dt = T1 / steps_per_period
    # PREFER `rate`. Deriving the speed from `periods_to_target` scales it with the TARGET, so the
    # same `periods_to_target` drives a larger pushover proportionally faster — and the recorded
    # base shear is a sum of reactions that includes inertial and DAMPING (C*v) terms, both of which
    # grow with speed. Calibrate `rate` once against a static run and it transfers between targets.
    drive_rate = float(rate) if rate is not None else abs(target) / (periods_to_target * T1)
    if drive_rate <= 0.0:
        raise ValueError(f"drive rate must be positive, got {drive_rate}")
    # `drive_rate` is a SPEED (positive); the direction rides on the sp() coefficient below, so a
    # negative `target` (e.g. a compression pushover) drives backwards instead of stopping at once.
    sign = 1.0 if target >= 0.0 else -1.0
    total_time = abs(target) / drive_rate
    # Rayleigh on the INITIAL stiffness (3rd arg), never the current tangent (2nd): a betaK term
    # rides the committed tangent, so when struts crack at speed the damping force explodes and
    # diverges the solve (D49). betaKinit gives the same damping from a matrix that never collapses.
    # UNDER AN EXPLICIT INTEGRATOR THE betaKinit TERM MUST GO. Stiffness-proportional damping adds
    # xi = betaK*w_max/2 at the highest mode, and the stable step shrinks by (sqrt(1+xi^2) - xi) —
    # with these lattices that is a factor of ~90, which would make an explicit run slower than the
    # implicit one it is meant to replace. Mass-proportional damping costs nothing in stability.
    if is_explicit(integrator):
        ops.rayleigh(damping_ratio * w1, 0.0, 0.0, 0.0)
    else:
        ops.rayleigh(damping_ratio * w1, 0.0, damping_ratio / w1, 0.0)

    # impose the ramped lateral displacement on the drive nodes (disp = sign * rate * t)
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    for nid in drive:
        ops.sp(nid, control_dof, sign * drive_rate)

    if is_explicit(integrator):
        # A diagonal (lumped-mass) system is the point of an explicit march: no factorization. The
        # algorithm is Linear because there is nothing to iterate — the step is an explicit update.
        ops.system("Diagonal")
        ops.test("NormDispIncr", tol, 1)
        ops.algorithm("Linear")
    else:
        ops.test("NormDispIncr", tol, max_iter)
        ops.algorithm(*algorithm)
    ops.integrator(*integrator)
    ops.analysis("Transient")

    def cdisp():
        return ops.nodeDisp(control_node, control_dof)

    def base_shear():
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base)

    def dynamic_residual():
        """Inertia + damping riding in the recorded shear. Call right after `base_shear` — it
        reuses that `ops.reactions()` rather than recomputing the whole domain."""
        return shear[-1] - sum(ops.nodeReaction(n)[control_dof - 1] for n in drive)

    groups: dict[str, list[float]] = {label: [] for label in (element_groups or {})}

    def record_groups() -> None:
        for label, members in (element_groups or {}).items():
            groups[label].append(sum(ops.eleForce(eid)[dof] * coef for eid, dof, coef in members))

    probe = _node_probe(node_history)
    nh = ({"nodes": probe[0], "dof": probe[1][0], "dofs": list(probe[1]),
           "index": [], "values": []} if probe else None)

    def record_probe() -> None:
        if nh is not None and (len(disp) - 1) % max(1, node_history_every) == 0:
            nh["index"].append(len(disp) - 1)
            nh["values"].append(_probe_sample(nh))

    disp = [cdisp()]; shear = [base_shear()]
    dynamic = [dynamic_residual()] if quasi_static else []
    record_probe()
    record_groups()
    converged = True
    nsteps = int(round(total_time / dt))
    goal = abs(target)
    disps_peak, peak_abs = (_nodal_disps(model) if capture else None), abs(shear[0])
    explicit = is_explicit(integrator)
    for _istep in range(nsteps):
        if ops.analyze(1, dt) != 0:
            if explicit:
                # An explicit step does not "fail to converge" — it goes unstable, and retrying with
                # a different algorithm cannot help. The fix is a smaller dt, which is the caller's.
                converged = False
                break
            sub_ok = False
            for algo in (("KrylovNewton",), ("NewtonLineSearch",)):
                ops.algorithm(*algo)
                if ops.analyze(10, dt / 10.0) == 0:
                    sub_ok = True
                    break
            ops.algorithm(*algorithm)
            if not sub_ok:
                converged = False
                break
        disp.append(cdisp()); shear.append(base_shear())
        if quasi_static:
            dynamic.append(dynamic_residual())
        record_probe()
        record_groups()
        if capture and abs(shear[-1]) > peak_abs:
            peak_abs = abs(shear[-1])
            disps_peak = _nodal_disps(model)
        if progress is not None and (_istep + 1) % max(1, progress_every) == 0:
            progress(_istep + 1, nsteps, disp[-1], shear[-1])
        if sign * cdisp() >= goal - 1e-9:
            break
    out = {"integrator": tuple(integrator),          # PROVENANCE (D70), see run_cyclic_dynamic
           "converged": converged, "disp": disp, "shear": shear,
           "control_node": control_node, "base_nodes": base, "rate": drive_rate, "T1": T1,
           # The marched step. `disp`/`shear` are sampled every step, so `dt` is what turns a
           # sample index into a TIME — needed to smooth a recorded series over a physical window
           # (a released crack or a ruptured bar rings locally at periods far below T1, D92).
           "dt": dt}
    if quasi_static:
        out["dynamic"] = dynamic
    if nh is not None:
        out["node_history"] = nh
    if capture:
        out["disps_final"] = _nodal_disps(model)
        out["disps_peak"] = disps_peak
    if element_groups is not None:
        out["groups"] = groups
    return out


def run_cyclic_dynamic(
    model: Model,
    *,
    control_node: int,
    control_dof: int,
    history,
    drive_nodes=None,
    gravity_loads=None,
    gravity_steps: int = 10,
    base_nodes=None,
    periods_to_peak: float = 48.0,
    rate: float | None = None,
    steps_per_period: int = 30,
    damping_ratio: float = 0.8,
    damping_modes: int = 5,
    tol: float = 1e-5,
    max_iter: int = 50,
    quasi_static: bool = False,
    node_history=None,
    node_history_every: int = 1,
    capture: bool = False,
    progress=None,
    progress_every: int = 200,
    integrator=("HHT", 0.7),
) -> dict:
    """Dynamic-relaxation REVERSED-CYCLIC analysis (D49): the cyclic sibling of
    `run_pushover_dynamic`.

    A static path-follower cannot cross the softening instability of a cracking lattice (D38/D44/D45);
    dynamic relaxation rides through it because mass and damping regularize the singular tangent
    (D46). This drives the whole reversing history rather than a single ramp.

    Mechanism: the displacement history is carried by a `Path` time series and imposed on
    `drive_nodes` through single-point constraints, so one continuous transient solve traces every
    cycle — no pattern teardown at reversals, and the imposed displacement stays continuous across
    them. The model is marched with `integrator` (HHT(0.7) by default) plus heavy Rayleigh damping.

    QUASI-STATIC RATE — the parameter that decides whether the answer means anything. The recorded
    base shear is a sum of REACTIONS, which in a transient solve include inertial and damping
    contributions; it is the quasi-static cyclic response only if the drive is slow enough that
    those are negligible. Give `rate` (displacement units per second) directly, or leave it None to
    derive it as `max|history| / (periods_to_peak * T1)`. Either way it is held CONSTANT for every
    segment, so a larger cycle simply takes proportionally longer.

    PREFER `rate`. Deriving from `periods_to_peak` scales the speed with the protocol's amplitude,
    so the SAME `periods_to_peak` drives a 4%-drift protocol four times faster than a 1% one and
    silently contaminates the larger run. `rate` is amplitude-independent and transfers between
    protocols. Calibrate it once by running a small-amplitude cycle both statically and
    dynamically: when the peak shears agree, the rate is slow enough.

    `quasi_static=True` MEASURES that, instead of leaving it to the calibration argument above.
    Global equilibrium in the drive direction reads `S_base + S_drive = -sum_free (M.a + C.v)`,
    where each S is the `-sum(reaction)` over that node set: internal forces cancel globally, so
    whatever the two constrained sets do NOT balance is exactly the inertial + damping force riding
    in the recorded base shear. (The horizontal external load is zero — gravity is vertical and
    held — so nothing else enters the sum.) It costs one extra reaction sum over the drive nodes
    per step, reusing the `ops.reactions()` the base shear already triggered, and returns the
    residual as "dynamic". It assumes `base_nodes` and `drive_nodes` between them cover EVERY
    constrained DOF in the drive direction (true whenever the supports are the base row) — any other
    restraint would leak into the residual and read as contamination that is not there.

    READ IT RELATIVELY, NOT ABSOLUTELY (D62). The residual is the instantaneous unbalance at the
    COMMITTED state, and the HHT integrator below leaves a large numerical-dissipation term there
    that does NOT bias the recorded shear: on an elastic wall this runner reads a residual of twice
    the base shear while reproducing the static answer to 0.24%, where the Newmark
    `run_pushover_dynamic` reads 27% for a 2.6% error. So `max|dynamic|/max|shear|` compares two
    runs OF THE SAME RUNNER (it is linear in `rate`, and falls with `damping_ratio` once cracking
    starts); it is not a portable error estimate. For an absolute check, drive an ELASTIC model and
    compare the recorded shear against the static solution.

    `integrator` is the time integrator as an `ops.integrator(*args)` tuple, defaulting to the
    implicit `("HHT", 0.7)`. Pass `("CentralDifference",)` for an EXPLICIT march (D74): ~45x cheaper
    per step against ~20x more steps, so ~2.7x faster end to end, and it needs no tangent, which is
    what lets it walk through the softening a Newton solve has to converge past. THREE THINGS CHANGE
    AUTOMATICALLY under it — the `betaKinit` damping term is dropped (it would shrink the stable
    step ~90x), the solver becomes Diagonal + Linear, and the implicit sub-step rescue ladder is
    skipped (a failed explicit step is instability, not non-convergence). THE CALLER MUST STILL SIZE `steps_per_period`
    FROM `builders.critical_time_step`, NOT from T1: an explicit march is only conditionally stable,
    and a step sized off T1 is typically 15-20x too large, which diverges rather than converging
    slowly.

    Cost scales as `total_path * periods_to_peak * steps_per_period / max|history|`, so a full
    protocol is tens of thousands of steps — hours, not minutes. `progress(i, nsteps, disp, shear)`
    is called every `progress_every` steps if given, so a long run can report where it is.

    `node_history=(node_ids, dof)` follows those nodes' displacement component `dof` (1-based) —
    the gauge-length probe a STRAIN PROFILE needs (D63): two aligned node rows a fixed distance
    apart give `(u_top - u_bottom)/gauge`, which is what a vertical LVDT pair on a wall face
    measures. Sampled every `node_history_every` steps AND at every reversal, because the loop tips
    are the instants a profile is wanted at and on a protocol this long a bare stride would miss
    them by up to a full stride of drive. Keep the stride coarse: a 4% protocol is ~1.8M steps, so
    `every=1` would hold tens of millions of floats.

    Returns {"converged", "disp", "shear", "control_node", "base_nodes", "steps", "rate", "T1"},
    plus "dynamic" when `quasi_static` is set and "node_history" =
    {"nodes", "dof", "index", "values"} when `node_history` is given (`index` holds positions into
    `disp`/`shear`).
    """
    if not model.masses:
        raise ValueError("run_cyclic_dynamic requires nodal mass on the model")
    hist = [float(h) for h in history]
    if not hist:
        raise ValueError("history is empty")
    peak = max(abs(h) for h in hist)
    if peak <= 0.0:
        raise ValueError("history must contain a non-zero amplitude")

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
                    "control_node": control_node, "base_nodes": base, "steps": 0}
        ops.loadConst("-time", 0.0)

    # Drop the gravity stage's Static analysis. While it exists OpenSees REFUSES to set a transient
    # integrator ("can't set transient integrator in static analysis") and silently falls back to a
    # default, so the integrator chosen below would be ignored. wipeAnalysis clears only the analysis
    # aggregation — domain, gravity loadConst and mass persist — so the solver must be re-declared.
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")

    try:                                     # the default solver is far faster than fullGenLapack
        lams = ops.eigen(max(1, damping_modes))
    except Exception:
        lams = ops.eigen("-fullGenLapack", max(1, damping_modes))
    w1 = math.sqrt(lams[0]) if lams[0] > 0 else 1.0
    T1 = 2.0 * math.pi / w1
    drive_rate = float(rate) if rate is not None else peak / (periods_to_peak * T1)
    if drive_rate <= 0.0:
        raise ValueError(f"drive rate must be positive, got {drive_rate}")
    dt = T1 / steps_per_period
    explicit = is_explicit(integrator)

    # Rayleigh damping on the INITIAL stiffness (3rd arg), never the current tangent (2nd arg).
    # A `betaK` term rides the COMMITTED tangent, so when struts crack at speed the damping force
    # explodes — the D28 artifact, and in a cracking lattice it diverges the solve outright rather
    # than merely spiking the reaction. `betaKinit` gives the same relaxation damping from a matrix
    # that never collapses. (Modal damping, D33's answer for the seismic runner, is NOT usable here:
    # dynamic relaxation needs damping near critical, and modalDamping at that level is pathological
    # — it fails to converge from the first step.)
    # UNDER AN EXPLICIT INTEGRATOR THE betaKinit TERM MUST GO (D74). Stiffness-proportional
    # damping shrinks the stable step by roughly (sqrt(1+xi^2) - xi) at the highest mode, which on
    # this lattice is a ~90x reduction — it fails as a run that never finishes, not as an error.
    if explicit:
        ops.rayleigh(damping_ratio * w1, 0.0, 0.0, 0.0)
    else:
        ops.rayleigh(damping_ratio * w1, 0.0, damping_ratio / w1, 0.0)

    # Piecewise-linear displacement path at constant speed `drive_rate`; the Path series then IS the
    # imposed displacement, so sp() carries a unit factor.
    times, values, cur = [0.0], [0.0], 0.0
    for goal in hist:
        times.append(times[-1] + abs(goal - cur) / drive_rate)
        values.append(goal)
        cur = goal
    total_time = times[-1]

    ops.timeSeries("Path", 3, "-time", *times, "-values", *values)
    ops.pattern("Plain", 3, 3)
    for nid in drive:
        ops.sp(nid, control_dof, 1.0)

    if explicit:
        # Nothing to iterate: the step is an explicit update, so Linear + a diagonal (lumped-mass)
        # system is both correct and ~45x cheaper per step than a Newton solve (D74).
        ops.system("Diagonal")
        ops.test("NormDispIncr", tol, max_iter)
        ops.algorithm("Linear")
    else:
        ops.test("NormDispIncr", tol, max_iter)
        ops.algorithm("Newton")
    ops.integrator(*integrator)  # HHT(0.7) by default: numerical damping aids convergence through
    ops.analysis("Transient")    # cracking and dissipates high-frequency content

    def cdisp() -> float:
        return ops.nodeDisp(control_node, control_dof)

    def base_shear() -> float:
        ops.reactions()
        return -sum(ops.nodeReaction(n)[control_dof - 1] for n in base)

    def dynamic_residual() -> float:
        """Inertia + damping riding in the recorded shear. Call right after `base_shear` — it
        reuses that `ops.reactions()` rather than recomputing the whole domain."""
        return shear[-1] - sum(ops.nodeReaction(n)[control_dof - 1] for n in drive)

    probe = _node_probe(node_history)
    nh = ({"nodes": probe[0], "dof": probe[1][0], "dofs": list(probe[1]),
           "index": [], "values": []} if probe else None)

    snap = {"peak": None, "peak_abs": -1.0}

    def record_probe() -> None:
        """Sample the gauge on the stride, and always at a REVERSAL (the loop tips)."""
        if nh is None:
            return
        turned = (len(disp) >= 3
                  and (disp[-1] - disp[-2]) * (disp[-2] - disp[-3]) < 0.0)
        if turned or (len(disp) - 1) % max(1, node_history_every) == 0:
            nh["index"].append(len(disp) - 1)
            nh["values"].append(_probe_sample(nh))

    disp, shear = [cdisp()], [base_shear()]
    dynamic = [dynamic_residual()] if quasi_static else []
    record_probe()
    nsteps = max(1, int(round(total_time / dt)))
    converged = True
    for i in range(nsteps):
        if ops.analyze(1, dt) != 0:
            if explicit:
                # An explicit step does not "fail to converge" — it goes unstable, and retrying with
                # a different algorithm cannot help. The fix is a smaller dt, which is the caller's.
                # Worse, the ladder below ends by SETTING Newton, so without this gate one failed
                # step silently converts the rest of the march to an implicit solve on a Diagonal
                # system, with nothing in the output to say so (D80 item 6).
                converged = False
                break
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
        disp.append(cdisp())
        shear.append(base_shear())
        if quasi_static:
            dynamic.append(dynamic_residual())
        # AFTER the append: testing `shear[-1]` before it holds this step's value snapshots the
        # displacement field of step i against the shear of step i-1, so the "peak" field is one
        # step stale — which matters, since D78 diagnosed the Aldemir tear off exactly this field.
        if capture and abs(shear[-1]) > snap["peak_abs"]:
            snap["peak_abs"] = abs(shear[-1])
            snap["peak"] = _nodal_disps(model)
        record_probe()
        if progress is not None and (i + 1) % max(1, progress_every) == 0:
            progress(i + 1, nsteps, disp[-1], shear[-1])

    if capture:
        out_capture = {"disps_final": _nodal_disps(model), "disps_peak": snap["peak"]}
    # PROVENANCE (D70): the integrator is REPORTED, not assumed — `tuple(integrator)`, never a
    # literal. This runner defaults to HHT(0.7) while `run_pushover_dynamic` defaults to Newmark,
    # and a day was lost to screening one with the other's solver because neither logged which it
    # used. A hardcoded literal here was the same failure one level down: both September cyclic
    # runs were marched with CentralDifference and their data.json said HHT (D80 item 6).
    out = {"integrator": tuple(integrator),
           "converged": converged, "disp": disp, "shear": shear, "control_node": control_node,
           "base_nodes": base, "steps": len(disp) - 1, "rate": drive_rate, "T1": T1,
           "dt": dt}          # the marched step — see run_pushover_dynamic (D92)
    if quasi_static:
        out["dynamic"] = dynamic
    if nh is not None:
        out["node_history"] = nh
    if capture:
        out.update(out_capture)
    return out


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


def run_dispbeamcolumn_tension(*, length: float, width: float, depth: float, material,
                               dU: float, target: float, n_int: int = 3) -> dict:
    """Single DISPLACEMENT-based fiber beam-column in uniaxial tension — the reference for the plain
    concrete cube lattice (a material coupon: one member, two nodes). Node 1 fixed at the origin,
    node 2 at (0, length) pulled axially (+Y) under DisplacementControl; no gravity. The fiber
    section is the cube face (`width` in local z x `depth` in local y) of one plain-concrete
    `material` (a backend-agnostic `UniaxialMaterial` spec carrying `.mtype` + `.args`).

    Why displacement-based: a `dispBeamColumn` interpolates the axial displacement LINEARLY, so the
    axial strain is CONSTANT along the element — the whole element carries eps = u/length uniformly
    (no integration-point localization, unlike a force-based element), and the axial response is
    exactly A*sigma(eps). That makes it the clean 1D reference for the concrete's uniaxial tension
    law. Returns the pushover curve {"disp", "shear", "converged"} where shear = the axial base
    reaction (tension positive)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, 0.0, length)
    ops.fix(1, 1, 1, 1)

    ops.uniaxialMaterial(material.mtype, 1, *material.args)
    y1, z1 = depth / 2.0, width / 2.0
    ops.section("Fiber", 1)
    ops.patch("quad", 1, 4, 4, -y1, z1, -y1, -z1, y1, -z1, y1, z1)   # homogeneous plain-concrete face
    ops.geomTransf("Linear", 1)                                       # small strain (tension), no P-Delta
    ops.beamIntegration("Legendre", 1, 1, n_int)
    ops.element("dispBeamColumn", 1, 1, 2, 1, 1)

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 1.0, 0.0)   # axial reference pull (+Y); DisplacementControl sets the magnitude
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 1000); ops.algorithm("Newton")
    ops.analysis("Static")
    return _displacement_pushover(2, 2, [1], dU, target)


def run_beamcolumn_dynamic(*, height: float, P: float, materials,
                           self_mass: float, top_mass: float, accel, dt_record: float,
                           scale: float, damping_ratio: float = 0.05,
                           modes: tuple[int, int] = (1, 2), dt: float = 0.01,
                           tol: float = 1e-6) -> dict:
    """Nonlinear seismic time-history of the fiber `forceBeamColumn` column — the reference for the
    lattice (matched material/mass). The column is a SINGLE force-based fiber element (fixed base
    node 1, free top node 2, 5 Lobatto points), matching the single-element pushover reference
    (`run_beamcolumn_cantilever`) and `single_beamcolumn.py`. `self_mass` is carried as DISTRIBUTED
    element mass (`-mass self_mass/height`) and `top_mass` (= P/g) is lumped at the free top. Fixed
    base, constant axial `P` (held), 5% modal damping at `modes`, then UniformExcitation of the
    scaled `accel` record in X. Returns the _transient_uniform_excitation dict (carrying "periods")."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, 0.0, height)
    ops.fix(1, 1, 1, 1)
    ops.mass(2, top_mass, top_mass, 0.0)                     # axial-load tributary seismic mass (top)

    _rc_fiber_section(1, materials=materials)
    ops.geomTransf("PDelta", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1, "-mass", self_mass / height)  # distributed self-mass

    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(2, 0.0, -P, 0.0)
    ops.system("BandGeneral"); ops.numberer("RCM"); ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-8, 100); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    if ops.analyze(10) != 0:
        return {"converged": False, "stage": "gravity", "t": [], "disp": [], "shear": [],
                "control_node": 2, "base_nodes": [1]}
    ops.loadConst("-time", 0.0)

    res = _transient_uniform_excitation(accel=accel, dt_record=dt_record, scale=scale, dof=1,
                                        control_node=2, control_dof=1, base_nodes=[1],
                                        dt=dt, zeta=damping_ratio, modes=modes, tol=tol)
    return res


def run_beamcolumn_modal(*, height: float, materials, self_mass: float,
                         num_modes: int) -> dict:
    """First `num_modes` modes of the SINGLE-element fiber `forceBeamColumn` cantilever — the modal
    counterpart of `run_beamcolumn_dynamic`, for the calibration mode-shape report (D35).

    One force-based fiber element (fixed base node 1, free top node 2, 5 Lobatto points), matching
    the dynamic reference. `self_mass` is carried as DISTRIBUTED element mass (`-mass`) so the total
    mass matches the lattice/continuum (the periods stay comparable, D16). No top mass and no gravity,
    so the eigen sees the initial-tangent stiffness at the undeformed state — like the lattice/
    continuum modal in `run_modal`. NB: a single element resolves only the first flexural mode well;
    higher requested modes are coarse (the accepted single-element trade-off). Returns
    {"periods", "shapes" (per mode, {node: [ux, uy, rot]}), "model"} where `model` is a lightweight
    drawing `Model` (two nodes + one segment) for the visualizer."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0); ops.node(2, 0.0, height)
    ops.fix(1, 1, 1, 1)

    _rc_fiber_section(1, materials=materials)
    ops.geomTransf("PDelta", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1, "-mass", self_mass / height)  # distributed self-mass

    try:
        eigenvalues = ops.eigen(num_modes)                  # genBandArpack: lowest modes (massless rot DOFs ok)
    except Exception:
        eigenvalues = ops.eigen("-fullGenLapack", num_modes)
    periods = [2.0 * math.pi / math.sqrt(lam) if lam > 0.0 else math.inf for lam in eigenvalues]
    shapes = [
        {i: ops.nodeEigenvector(i, mode) for i in range(1, 3)}
        for mode in range(1, num_modes + 1)
    ]

    model = Model(ndm=2, ndf=2)                              # drawing-only model (one line along the height)
    model.add_node(1, (0.0, 0.0)); model.add_node(2, (0.0, height))
    model.add_element(1, "line", (1, 2), kind="concrete")
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
