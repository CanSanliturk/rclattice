# TASK: Solve an RC Frame Pushover with Lattices

Handoff document so a fresh session can implement this end-to-end. Read
[CLAUDE.md](CLAUDE.md) and [DECISIONS.md](DECISIONS.md) first for project conventions and the
settled decisions (D1–D17); this file is the spec for the next decision (D18) and its
implementation.

**Goal:** reproduce the classic OpenSees RC frame pushover
(https://openseespydoc.readthedocs.io/en/latest/src/RCFramePushOver.html, which imports
`RCFrameGravity`) using the project's **lattice** modelling — and compare against the
continuum builder and/or the original beam-column result.

> **Standing rule (memory `ask-on-mechanical-decisions`): the user decides every
> mechanical/structural/modelling point. Do NOT assume. The four forks in §6 are still OPEN —
> ask the user (AskUserQuestion) and record as D18 before implementing.**

---

## 1. Current project state (what already works)

Environment: native arm64, Python 3.12 venv at `src/.venv` (uv). Run everything from `src/`:
`uv run python examples/<x>.py`, `uv run --with pytest pytest tests/`. Deps: openseespy 3.8,
gmsh, numpy, scipy, matplotlib, pillow.

Package `src/rclattice/` (all backend-agnostic except `opensees.py`):
- `problem.py` — `Problem(ndm, ndf, domain, material, supports, loads)`; `RectangleDomain`,
  `CompoundRectangles` + `portal_frame(...)` (compound rectangles, joints merged);
  `ConcreteGrade(name, E, nu, rho, fc?, ft?, Gf?)`; selectors `EdgeSupport/EdgeLoad` and
  `BoxSupport/BoxLoad`.
- `mesh.py` — `mesh_rectangle_grid` (structured quads), `mesh_compound_rectangles` (per-rect
  structured mesh, coincident nodes merged), `connect_horizon` (peridynamics horizon struts,
  default 1.5×mesh_size, dedup).
- `model.py` — generic FE `Model`: `Node`, `Element(etype, nodes, args)`, `UniaxialMaterial`,
  `NDMaterial`, `Support`, `Load`, `masses`.
- `materials.py` — grade→OpenSees mapping. **Elastic only so far**:
  `concrete_uniaxial_elastic` (→ `Elastic`), `concrete_nd_elastic` (→ `ElasticIsotropic`).
- `builders.py` — `build_lattice` (horizon truss struts, `strut_area` float|callable, lumped
  tributary mass), `build_continuum` (2D quads, PlaneStress/Strain). Both share the node grid;
  apply supports/loads via edge or box selectors; return `(Model, edges)`.
- `opensees.py` — ONLY module importing openseespy. `build(model)`, `run_static(model)`
  (linear static, returns `{ok, disps}`), `run_modal(model, n)` (eigen → `{eigenvalues,
  periods, shapes}`).
- `calibration.py` — `continuum_targets`, `calibrate_lattice` (orthogonal/diagonal strut
  areas, bounded scipy least_squares, static + first-N periods), `static_response`,
  `combined_rms`, `nominal_area`.
- `viz.py` — matplotlib (Agg): `figure_static`, `figure_modes`, `animate_modes`, `draw_model`.

Examples (in `src/examples/`, outputs to `examples/output/`): `cantilever_hello`,
`cantilever_verify` (elastic lattice↔continuum), `cantilever_calibrate` (modal calibration),
`cantilever_visualize` (shows higher-mode mismatch), `frame_visualize` (portal frame, static +
modes side-by-side + GIF; lattice↔continuum periods agree ~1–3%). Tests: `test_horizon`,
`test_verification`, `test_calibration`, `test_frame` — **9 passing**.

Key established facts: axial lattice has a near-fixed Poisson ratio and can't fully match
higher modes of a solid (documented limitation); flexural frames are matched far better.
`ndm=2, ndf=2` everywhere so far (translational DOFs only — lattice is axial truss).

---

## 2. The benchmark problem (units: kip, in, sec)

`ndm=2, ndf=3` in the original (beam-column with rotations). 1 bay × 1 storey.

Nodes: 1 (0,0), 2 (360,0) — **both fully fixed** (1,1,1); 3 (0,144), 4 (360,144).
- Columns: elem 1 (1→3), elem 2 (2→4). `forceBeamColumn`, P-Delta transf, 5 Lobatto points,
  **RC fiber section** (tag 1). Height 144 in.
- Beam: elem 3 (3→4). `elasticBeamColumn`, A=360, E=4030, Iz=8640, Linear transf. Span 360 in.
  (Beam is ELASTIC in the benchmark.)

Materials:
- Concrete01 core (confined), tag 1: fpc −6.0, epsc0 −0.004, fpcu −5.0, epsU −0.014 (ksi).
- Concrete01 cover (unconfined), tag 2: fpc −5.0, epsc0 −0.002, fpcu 0.0, epsU −0.006.
- Steel01, tag 3: fy 60.0, E0 30000.0, b 0.01.

Column fiber section (tag 1): **15 wide (z) × 24 deep (y), cover 1.5 in.** Core concrete patch
(mat 1) over the confined region; cover concrete patches (mat 2) on the 1.5 in ring.
Reinforcing steel (mat 3), **As = 0.60 in²/bar**: 3 bars at y=+10.5, 2 bars at y=0, 3 bars at
y=−10.5 (i.e. top/mid/bottom layers in the section depth direction).

Loads / analysis:
- Gravity: `load(3, 0, −180, 0)` and `load(4, 0, −180, 0)` kip; LoadControl, 10 steps; Newton;
  system BandGeneral, constraints Transformation, numberer RCM, test NormDispIncr 1e-12.
- Pushover: reference lateral H=10 kip at nodes 3 & 4 (DOF 1); `integrator
  DisplacementControl 3 1 dU` with dU=0.1 in to target 15.0 in; ModifiedNewton (initial
  stiffness); test NormDispIncr 1e-12, 1000 iters.

---

## 3. Lattice idealization of the benchmark

Reuse the compound-rectangle frame (D17) in **kip-in**:
- The frame is in the global X–Y plane. The column's in-plane bending depth is **24 in** (the
  section's `y`/depth), so each **column is a 2D region 24 in wide in X × 144 in tall in Y**,
  out-of-plane thickness **15 in**. Left column centerline at x=0 → x∈[−12,12]; right at x=360
  → x∈[348,372]. Beam at top, span 360, with a chosen depth (see fork §6).
- **Rebar mapping (2D):** longitudinal bars at section depth y=±10.5 map to **vertical steel
  lines at x = ±10.5** in each column (extreme fibers in X), running the full 144 in height;
  the mid layer (y=0, 2 bars) maps to x=0. In 2D we lump out-of-plane bars: x=±10.5 layers =
  3×0.6 = **1.8 in² each**; x=0 layer = 2×0.6 = **1.2 in²**.
- **Mesh size = 1.5 in** is the sweet spot: column depth 24/1.5 = 16 cells → nodes land exactly
  at x=−12,−10.5,…,10.5,12, so rebar lines coincide with lattice nodes and the 1.5 in cover is
  exactly one cell. (Story 144/1.5 = 96 cells; this is a biggish model — see §8 perf.)
- Supports: fix the two column-base regions (BoxSupport at y≈0). Note `ndf=2` (no rotational
  fixity needed; bases pinned-but-the-region-width provides rotational restraint via the
  lattice). Control node: the top-left joint (≈ node 3, at (0,144)); control DOF = X.
- Base shear = sum of horizontal reactions at base nodes (need reaction extraction).
- Pushover curve = base shear vs. control (roof) X-displacement.

---

## 4. New capabilities required (file-by-file plan)

1. **Nonlinear materials (`materials.py`, `problem.py`)**
   - Extend `ConcreteGrade` (or add a `ConcreteNonlinear` spec) with Concrete01/02 params:
     fpc, epsc0, fpcu, epsU, and tensile ft (+ Ets) for Concrete02. Add a `SteelGrade(fy, E0,
     b, ...)` (Steel01/Steel02).
   - Add mappings: `concrete_uniaxial_nonlinear(grade, tag)` →
     `UniaxialMaterial(tag, "Concrete02", (...))`, `steel_uniaxial(grade, tag)` →
     `UniaxialMaterial(tag, "Steel02"/"Steel01", (...))`. (Concrete02 arg order:
     `Fc, epsc0, Fcu, epsU, lambda, ft, Ets` — verify against openseespy docs.)

2. **Reinforcement (the deferred D13) — `problem.py`, `builders.py`, new `reinforcement.py`**
   - `Rebar(path, area, steel_grade)` where `path` is a polyline (list of (x,y)); for the
     frame, vertical lines at x=±10.5 and x=0 from y=0 to y=144, plus beam bars if the beam is
     RC (fork §6).
   - In `build_lattice`: snap rebar paths to existing lattice nodes (mesh aligned so nodes lie
     on the path), then create **steel truss struts** connecting consecutive on-path nodes
     (shared nodes ⇒ perfect bond, D5). Assign the steel area & material.
   - (Continuum embedded bars via gmsh `embed` is the D13 plan but NOT needed for the lattice
     pushover; defer unless comparing to a continuum RC model.)

3. **Nonlinear analysis runners (`opensees.py`)**
   - `run_gravity(model, *, nsteps=10)` — LoadControl, Newton, NormDispIncr test; apply the
     gravity loads (vertical). Keep the gravity pattern loaded (constant) into the pushover.
   - `run_pushover(model, *, control_node, control_dof, dU, target, ...)` — second `Plain`
     pattern with the lateral reference load; `integrator DisplacementControl`; ModifiedNewton;
     step to target; at each step record control displacement and **base shear** (sum of
     `nodeReaction` X at base nodes — call `ops.reactions()` first). Return the pushover curve
     (arrays disp, shear) + convergence info. Handle non-convergence gracefully (reduce step).
   - Reaction extraction helper; ensure base nodes are tracked (builder should expose them).

4. **Units** — everything is unit-agnostic; just feed kip-in numbers. (ρ only matters for
   modal; pushover is static so mass can be omitted/zero — but `run_modal` requires mass.)

5. **Visualization (`viz.py`)** — add a pushover-curve plot (base shear vs drift), overlay
   lattice vs continuum vs benchmark; optionally animate the pushover deformation.

6. **Example `examples/frame_pushover.py`** — define the RC frame Problem, build lattice (+
   rebar), run gravity then pushover, plot the curve, compare.

7. **Tests** — nonlinear material mapping smoke test; rebar struts created on the right nodes;
   gravity converges; pushover produces a monotonic-ish curve with a yield plateau.

---

## 5. Recommended staging

**Stage 1 — geometry + analysis machinery, ELASTIC (de-risk the pipeline):**
- Build the RC-frame geometry in kip-in (compound rectangles, columns 24×144 thk 15, beam at
  top), elastic lattice + continuum.
- Implement `run_gravity` (LoadControl) and `run_pushover` (DisplacementControl + base-shear
  recording) and the pushover-curve plot. Run elastically (linear pushover = straight line) to
  prove control node, reactions, stepping, and the curve plumbing.

**Stage 2 — nonlinear + reinforcement (the physics):**
- Add Concrete02 + Steel02 mappings and the `Rebar` feature; build the RC lattice; run
  gravity→pushover; compare the curve to the benchmark and/or continuum. Calibrate strut areas
  (elastic, existing tooling) first so initial stiffness matches.

---

## 6. OPEN DECISIONS — ask the user, record as D18 (recommendations in **bold**)

1. **Staging** — *user already leaned here via this doc; confirm:* **Stage it (Stage 1 then
   Stage 2).** Alt: full nonlinear in one go.
2. **Beam modelling** — **Elastic beam region (match benchmark; hinging in columns).** Alt:
   nonlinear RC beam region (needs beam rebar detailing).
3. **Concrete constitutive for lattice struts** — **Concrete02 with tensile strength +
   softening** (truss tension stability; keep benchmark fc=6/5 ksi). Alts: Concrete01
   compression-only (unstable in truss), Concrete02 + fracture-energy regularization (most
   rigorous).
4. **Steel + rebar layers** — **Steel02, all three layers** (x=±10.5 @1.8 in², x=0 @1.2 in²).
   Alts: Steel01 (exact benchmark match), outer layers only.

(The user answered these with "." = deferred/no-decision; they remain OPEN. Get explicit
answers before Stage 2 material/rebar work. Stage 1 can proceed without them.)

---

## 7. Verification targets

- Stage 1: linear pushover is a straight line; its slope = elastic lateral stiffness; sanity-
  check against the continuum and a hand estimate.
- Stage 2: pushover curve (base shear vs roof drift) should show elastic branch → yield
  plateau. Compare initial stiffness, yield/peak base shear, and post-yield slope to the
  benchmark beam-column result and/or a continuum RC model. Exact match is NOT expected (axial
  lattice + Concrete02 vs fiber beam-column); aim for the right shape and capacity ballpark,
  and document differences (this is verification, not replication).

---

## 8. Gotchas / pitfalls

- **Truss tension instability:** compression-only concrete ⇒ mechanisms. Use Concrete02 (fork
  3). Even so, watch for non-convergence; have step-reduction in `run_pushover`.
- **ndf mismatch:** benchmark is ndf=3 (rotations); our lattice/continuum are ndf=2. Fixed
  *region* bases provide rotational restraint through the lattice width — fine, but the control
  "node" is a joint region; pick a representative node (top-left corner or joint centroid) and
  state it.
- **Base shear:** call `ops.reactions()` before `ops.nodeReaction`; sum X-reactions over ALL
  base nodes (the base is a region of nodes, not one node).
- **Mass for modal vs none for pushover:** `run_modal` needs mass; pushover is static. Keep
  mass assignment (harmless in static) or guard.
- **Model size:** mesh 1.5 in → columns 16×96 nodes each + beam ⇒ thousands of nodes/struts.
  Static nonlinear is fine; if slow, coarsen to mesh 3.0 in (nodes at x=±10.5 need 24/3=8 cells
  → nodes at −12,−9,−6,−3,0,3,6,9,12 — that MISSES ±10.5! so 1.5 in, or 0.75 in, is required to
  hit ±10.5; do not coarsen below alignment). Alternatively relocate rebar to the nearest grid
  line and document the small offset.
- **Quad CCW / element args:** continuum quads already handled; for any new element types
  verify arg order against the compiled openseespy (probe with `strings` on the .so as in the
  env-setup history).
- **Units:** stay in kip-in throughout this task (do not mix SI).

---

## 9. Quick start (fresh session)

```bash
cd src
uv sync                                   # restore env
uv run --with pytest pytest tests/        # 9 should pass
uv run python examples/frame_visualize.py # see the working frame viz baseline
```
Then implement Stage 1 (§5), wiring `run_gravity`/`run_pushover` into `opensees.py` and a new
`examples/frame_pushover.py`. Confirm the §6 decisions with the user before Stage 2.
