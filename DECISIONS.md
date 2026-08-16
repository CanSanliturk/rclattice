# Decision Log

A running log of key technical decisions for `rclattice`. CLAUDE (and humans) should read
this for context before making changes, and append a new entry whenever a significant
decision is made or an old one is reversed. Newest entries at the bottom. Do not rewrite
history — if a decision changes, add a new entry that supersedes the old one and note it.

Entry format:
- **ID / date** — the decision (one line)
- **Why** — rationale
- **Alternatives rejected** — and why
- **Status** — accepted / superseded by Dxx

---

### D1 — 2026-06-14 — OpenSees backend is `openseespy` 3.8, run natively on arm64
- **Why:** `openseespy` 3.8 is a pure-python shim depending on `openseespymac`, which ships a
  native `macosx_13_0_arm64` wheel. Verified: imports and solves correctly on this machine
  (macOS 15, Apple Silicon). No Rosetta needed.
- **Alternatives rejected:** `opensees` (xara) — arm64 wheels only exist for Python 3.10 and
  3.13, not 3.11/3.12. Rosetta 2 / x86_64 conda env — obsolete workaround now that native
  arm64 wheels exist.
- **Status:** accepted

### D2 — 2026-06-14 — Tooling: uv + native arm64, Python 3.12
- **Why:** User wants uv/venv rather than conda. Native arm64 works (see D1), so no conda
  needed. uv chosen as primary env/dependency manager. Verified with CPython 3.12.13.
- **Alternatives rejected:** conda/Rosetta x86 env (unnecessary). `requires-python` is pinned
  to `>=3.12` because `openseespylinux` 3.8 requires Python >=3.12, so the universal resolver
  rejects anything lower; 3.12 is also the developed/verified target.
- **Status:** accepted

### D3 — 2026-06-14 — Dimension-agnostic core (2D and 3D equally)
- **Why:** User wants both 2D and 3D supported from the start. Nodes carry an `ndm`-length
  coordinate; dimensionality threaded explicitly through the model (`ndm`, `ndf`).
- **Alternatives rejected:** 2D-first then extend (faster but risks baking in 2D
  assumptions).
- **Status:** accepted

### D4 — 2026-06-14 — Lattice struts are axial trusses only
- **Why:** Classic concrete-lattice approach; nodes need only translational DOFs
  (`ndf = ndm`). Struts map to OpenSees `truss`/`corotTruss` with a uniaxial material.
- **Alternatives rejected:** Beam/frame lattice (struts carrying moment) — more DOFs and
  complexity than needed for the intended studies.
- **Caveat:** pure axial-truss lattices can form kinematic mechanisms if not triangulated;
  Delaunay-edge lattices mitigate this. Suspect this first on singular-stiffness errors.
- **Status:** accepted

### D5 — 2026-06-14 — Reinforcement uses shared, perfectly-bonded nodes
- **Why:** Simplest realistic coupling — rebar elements share (or snap onto) lattice nodes,
  giving perfect bond.
- **Alternatives rejected:** bond-slip via zeroLength interface springs (more realistic,
  more complex) — deferred. Interface to be designed so slip can be added later without
  reworking the core.
- **Status:** accepted

### D6 — 2026-06-14 — Mesher backend is gmsh (for NODE generation)
- **Why:** Don't hand-roll meshing math. gmsh has a mature Python API, native arm64 wheel
  (verified: meshed a box into nodes/tets), and handles complex geometry. gmsh places the
  lattice **nodes** (transfinite/structured grid for the regular pattern, Delaunay for
  irregular).
- **Alternatives rejected:** hand-rolled Delaunay (unreliable, reinventing the wheel).
  `scipy.spatial` — viable lighter alternative, kept as a possible second backend behind a
  backend-agnostic mesher interface, but gmsh is the implementation.
- **Status:** accepted. NOTE: the strut-connectivity half of the original D6 ("element edges
  -> struts") is **superseded by D9** — struts now come from a horizon rule, not mesh edges.

### D9 — 2026-06-14 — Strut connectivity via a peridynamics-style horizon
- **Why:** Connect every pair of nodes whose separation is <= `horizon * mesh_size`
  (default `horizon = 1.5`), deduplicated so only ONE element exists between any two nodes
  (reinforcement may be the exception later). On a regular grid this captures orthogonal
  neighbors (`s`) and diagonals (`s*sqrt(2) ~= 1.41 s`) but not the next ring (`2 s`), so the
  lattice is naturally triangulated and can carry shear/bending — no separate bracing step
  needed. Generalizes identically to irregular node clouds.
- **Supersedes:** the "element edges -> struts" rule from D6.
- **Status:** accepted

### D10 — 2026-06-14 — gmsh is the single node source for both patterns
- **Why:** With struts now horizon-based (D9), nodes still come exclusively from gmsh for
  both regular (transfinite/structured) and irregular (Delaunay) patterns. Keeps one node
  source, consistent with D6.
- **Alternatives rejected:** direct grid generation for the regular pattern (two node
  sources to maintain).
- **Status:** accepted

### D11 — 2026-06-14 — `src/` is the project root; flat package + hatchling
- **Why:** All development files (package, tests, examples) and project config
  (pyproject.toml, uv.lock, requirements.txt, .venv) live under `src/`, which is therefore
  the uv/build project root. CLAUDE.md and DECISIONS.md stay one level above it at the repo
  root. Package is flat at `src/rclattice/`; build backend is hatchling; `rclattice` is
  installed editable via `uv sync`.
- **Implication:** run `uv` / build commands from inside `src/`.
- **Alternatives rejected:** config at repo root (user preferred a clean root); nested
  `src/src` layout (redundant naming); setuptools (hatchling chosen).
- **Status:** accepted

### D12 — 2026-06-14 — "Define once, build many": shared Problem + multiple builders
- **Why:** Lattice is the main aim, but it must be *verified/calibrated* against continuum
  and (for single members) beam-column models of the SAME specimen. So define the specimen
  ONCE in a backend-agnostic `Problem` (geometry + reinforcement + physical materials +
  BCs + loads), and translate it with independent builders:
  1. **LatticeBuilder** — gmsh nodes + horizon struts (D9); uniaxial concrete struts + rebar
     struts on shared nodes. (Current `model.py`/`mesh.py`/`opensees.py` are the seed of this.)
  2. **ContinuumBuilder** — 3D solids / 2D plane elements (and shells) with nD concrete.
  3. **BeamColumnBuilder** — single prismatic members via `forceBeamColumn` + fiber section.
  All builders consume the same Problem, so their results are directly comparable. The
  continuum model is the reference for calibrating the lattice.
- **Alternatives rejected:** fully independent per-model definitions (no guarantee they
  represent the same specimen — defeats verification).
- **Status:** accepted. Verified on this OpenSees build: all required elements/materials
  for the three builders are compiled in (solids stdBrick/SSPbrick/tets, quads, shells
  LayeredShell/ShellMITC4, forceBeamColumn, ASDConcrete3D, Concrete02, Steel02).

### D13 — 2026-06-14 — Reinforcement is defined once as free 3D curves; discrete & embedded
- **Why:** Define rebar at the Problem level as free geometric curves (polyline path +
  diameter + steel grade) — exact, arbitrary placement. Each builder consumes the SAME
  definition: BeamColumn -> section fibers at the bar locations; Lattice -> rebar struts on
  shared nodes; Continuum -> **discrete embedded bars** using gmsh's `embed` to make the
  mesh conform to the rebar path, then truss/beam rebar elements sharing the solid/quad
  nodes (perfect bond, D5). Reuses the gmsh + shared-node machinery (D6/D8/D10).
- **Alternatives rejected:** smeared reinforcement (simpler but loses discrete bar
  locations and breaks like-for-like comparison with the lattice). May be added later as an
  option for slabs/walls.
- **Status:** accepted

### D14 — 2026-06-14 — 2D continuum supports BOTH plane-stress and shell elements
- **Why:** Plane-stress quads (`quad`/`SSPquad`/`tri31`, ndf=2) are the like-for-like
  continuum check for a planar (ndf=2) lattice. Shells (`ShellMITC4`/`ShellNLDKGQ` +
  `LayeredShell`, ndf=6) are needed when a "2D member" is really a thin wall/slab with
  out-of-plane bending. Support both, selectable per problem.
- **Alternatives rejected:** plane-stress only (can't model thin-plate bending); shell only
  (not DOF-matched to a planar lattice).
- **Status:** accepted

### D15 — 2026-06-14 — Physical material grades map to per-builder OpenSees materials
- **Why:** Materials are defined ONCE as physical grades (concrete: fc, ft, E, Gf, ...;
  steel: fy, E, b, ...). A mapping layer emits the appropriate OpenSees material per builder
  / element context. This mapping layer is also where lattice **calibration** lives (e.g.
  fracture-energy regularization of strut softening by strut length/area so the lattice
  assembly reproduces the continuum).
- **First material bundle:** uniaxial **Concrete02** (struts & fibers) + nD **ASDConcrete3D**
  (solids/shells) + **Steel02** (rebar/fibers). All confirmed present in this build.
- **Alternatives rejected (for now):** Concrete04/ConcreteCM, PlasticDamageConcrete3d,
  ReinforcingSteel — all available, kept as future options behind the grade->material map.
- **Status:** accepted

### D16 — 2026-06-15 — Lattice area calibration via static + modal periods (multi-group)
- **Why:** Improve lattice calibration beyond the single static-deflection scalar (D-area in
  cantilever_verify) by ALSO matching modal periods. A single uniform area can only match one
  quantity (all periods scale together as `1/sqrt(area)`), so use multiple **area groups**
  (on the regular grid: orthogonal struts length `s` vs diagonal `s*sqrt2`) — tuning their
  ratio adjusts stiffness anisotropy and partially addresses the axial-lattice Poisson limit.
  The existing single-scalar/static method is the degenerate 1-group, static-only case and is
  kept "on top of".
- **Mechanics:**
  - Add **density `rho`** to `ConcreteGrade`. Mass is assigned by **tributary volume × rho**
    (geometry-based, so independent of strut area → area tunes stiffness K only, mass M
    fixed). Lattice and continuum get the **same total mass**.
  - `run_modal(model, n)` via `ops.eigen` returns periods `T_i = 2*pi/sqrt(lambda_i)`.
  - Calibration fits the area-group parameters by **weighted least squares over the static tip
    deflection + the first N modal periods** (N default 3, configurable).
- **Target is selectable:** match the **continuum reference** (inside the verification
  framework) OR **user-supplied target periods** (e.g. experimental).
- **Implemented mechanics:** lumped tributary mass on both models (equal total); mode matching
  by ascending period order; optimizer = `scipy.optimize.least_squares`.
- **Status:** implemented (`calibration.py`, `run_modal`, mass in builders).
- **Open finding (2026-06-15):** for the cantilever, T1 is flexural and nearly redundant with
  the static target (both match easily); higher modes (T2+) diverge. An UNBOUNDED orthogonal/
  diagonal area fit lowers RMS ~26% but drives the diagonal area to ~1.6 m^2 (near-rigid,
  unphysical) chasing higher modes. => calibration needs physical conditioning (area bounds /
  regularization / mode weighting), and matching higher modes with an axial lattice may be
  fundamentally limited (a documented lattice limitation, not a bug). Decision pending.
- **Resolution:** physical area bounds (default 1e-3..3x nominal `t*s`) + weighting (static + T1
  firm, higher modes soft). Implemented in `calibrate_lattice`.

### D17 — 2026-06-15 — Compound-rectangle frame geometry + box selectors + matplotlib viz
- **Why:** Provide a visual lattice-vs-continuum comparison (static + mode shapes) on a
  one-bay/one-storey portal frame.
- **Geometry:** frames are built as **compound axis-aligned rectangles** (`portal_frame` =>
  2 columns + beam) meshed structured per rectangle with **coincident nodes merged** at the
  joints (`mesh_compound_rectangles`). Keeps the validated regular-grid lattice; member dims
  should be multiples of mesh size. `RectangleDomain` is the single-rectangle case.
- **Supports/loads:** generalized from edge-based to **box selectors** (`BoxSupport`/`BoxLoad`,
  axis-aligned region) alongside the existing `Edge*`. Frame scenario: fixed column bases +
  lateral (sway) load distributed over the beam. Calibration's static metric now reads the
  loaded node ids from the built model, so it is selector-agnostic.
- **Visualizer:** `viz.py` (matplotlib, Agg) draws deformed lattice (lines) / continuum (quads)
  side-by-side for static and per-mode (with periods), plus an animated-modes GIF. Pure
  drawing — takes (Model, disp dict); does not import openseespy. `run_modal` now also returns
  nodal **mode shapes**. Deps added: matplotlib, pillow.
- **Result:** static deflected shapes match; with static-scalar area calibration the frame
  periods agree within ~1-3% (continuum [68.0, 22.3, 10.8] ms vs lattice [68.2, 22.6, 11.1] ms)
  — a flexural frame is captured far better by the axial lattice than the solid cantilever.
- **Status:** implemented (`examples/frame_visualize.py`, tests in `test_frame.py`).

### D18 — 2026-06-16 — Solve the RC frame pushover benchmark with lattices (PENDING)
- **Goal:** reproduce the OpenSees RC frame pushover (RCFramePushOver + RCFrameGravity) using
  the lattice builder; compare to continuum / the original beam-column. Full spec, lattice
  idealization, file-by-file plan, staging, and gotchas are in
  [TASK_RC_FRAME_PUSHOVER.md](TASK_RC_FRAME_PUSHOVER.md).
- **Documented (so a new session can continue):** benchmark = 1-bay/1-storey, kip-in, columns
  RC fiber `forceBeamColumn` (15×24 section, Concrete01 core/cover, Steel01, 3+2+3 bars @0.6
  in²), elastic beam, gravity P=180 kip + DisplacementControl pushover to 15 in. Lattice plan:
  columns as 2D regions 24 in (in-plane) × 144 tall, rebar as vertical steel struts at x=±10.5
  & 0, mesh 1.5 in (aligns nodes to bars). New capabilities needed: nonlinear materials
  (Concrete02 + Steel02 mappings), reinforcement (D13: rebar struts on shared nodes), nonlinear
  runners (`run_gravity` LoadControl + `run_pushover` DisplacementControl with base-shear
  recording), reaction extraction, pushover-curve plot.
- **OPEN forks (user deferred with "." — must be answered before Stage 2):**
  (1) staging [recommend: stage it — elastic machinery first, then nonlinear+rebar];
  (2) beam [recommend: elastic beam, match benchmark];
  (3) lattice concrete law [recommend: Concrete02 with tension+softening — Concrete01 is
  unstable in truss struts]; (4) steel + rebar layers [recommend: Steel02, all three layers].
- **Status:** PENDING — documented, not implemented. Next session: confirm forks, do Stage 1.

### D19 — 2026-06-16 — RC frame pushover forks resolved (D18 unblocked); Stage 1 underway
- **Decision (user, AskUserQuestion, all four = recommended option):**
  (1) **Stage it** — Stage 1 builds the RC-frame geometry (kip-in) + `run_gravity`/`run_pushover`
  machinery + pushover-curve plot, ELASTIC, to de-risk control node / reactions / stepping;
  Stage 2 adds the nonlinear physics. (2) **Elastic beam** — match the benchmark; hinging in the
  columns. (3) **Concrete02 with tensile strength + softening** for the lattice struts (truss
  tension stability; keep benchmark fc=6/5 ksi). (4) **Steel02, all three rebar layers** —
  vertical steel struts at x=±10.5 (1.8 in² each) and x=0 (1.2 in²), perfect bond on shared
  nodes.
- **Why:** the recommended path keeps a clean like-for-like comparison with the OpenSees
  beam-column benchmark while avoiding the known axial-truss instabilities (compression-only
  Concrete01 → mechanisms).
- **Status:** accepted. Supersedes the "PENDING / forks open" status of D18.
- **Stage 1 done (2026-06-16):** elastic pipeline landed. `opensees.py` gains `run_gravity`
  (LoadControl) and `run_pushover` (gravity-constant → DisplacementControl, base shear =
  `-Σ reaction[control_dof]` over the base nodes, with step-halving on non-convergence);
  `builders.select_nodes` (post-build box query) picks control/base/load nodes; `viz.figure_pushover`
  plots the curve. [examples/frame_pushover.py](src/examples/frame_pushover.py) builds the
  benchmark frame in kip-in (E=4030, cols 24×144 thk 15, elastic beam), runs lattice + continuum
  pushovers, and calibrates strut area to the continuum lateral stiffness — they match (271.06 vs
  271.10 kip/in, both converge to 15 in). Found: a portal **splays** under gravity (antisymmetric
  drift at column tops) — expected, not a bug. 4 new tests in `test_pushover.py` (13 total pass).
  Next: Stage 2 (Concrete02 + Steel02 mappings, `Rebar` struts on shared nodes, nonlinear runs).

### D20 — 2026-06-16 — Lattice concrete softening uses fracture-energy regularization (supersedes D19 fork 3)
- **Context:** Stage 2 of the RC pushover. The D19 fork-3 choice (plain Concrete02 + tension/
  softening) makes the axial-truss lattice go **unstable shortly after first yield** — it cannot
  reach the 15-in target. Diagnostics (single 24×144 column, calibrated strut area): plain
  Concrete02 fails at 0.20 in; gentler tension softening reaches only ~0.66 in; **near-flat
  tension (Ets≈E/1000) still fails at 0.68 in** — so tension softening rate is NOT the binding
  constraint. The instability is driven by **compression** softening (Concrete02 descending
  branch; the cover crushing to fcu=0 over a short strain) creating negative tangent stiffness →
  singular stiffness / mechanism (the D4 truss-lattice caveat, TASK §8 top gotcha). Measured peak
  base shear (~30 kip on the cantilever column) is a sane capacity — the failure is in sustaining
  the post-yield DRIFT, not in strength.
- **Decision (user, AskUserQuestion):** adopt **fracture-energy regularization** (fork-3's
  rejected alternative). Regularize the Concrete02 softening of each strut by its **length** L
  (crack-band / Bažant), so dissipation is mesh-objective and, for the small struts (L≈1.5 and
  diagonal ≈2.12), the post-peak is gentle enough to stay convergent:
  - tension softening slope `Ets = ft²·L/(2·Gf)` (capped below E to avoid snap-back);
  - compression crushing strain `epsU = epsc0 + 2·Gfc/((fc+fcu)·L)` (gentle for small L; Gfc the
    compression fracture energy, ~250·Gf by default).
  Regularization depends on L only (energy per unit cross-section area), so it is independent of
  the elastic strut-area calibration. On the regular grid there are only two strut lengths, so
  each nonlinear zone yields just two materials (orthogonal + diagonal).
- **Status:** accepted; supersedes the lattice-concrete part of D19 fork 3. Implementation:
  `materials.concrete_uniaxial_regularized(grade, tag, length, *, Gf, Gfc)`; `build_lattice_rc`
  assigns concrete materials per (zone, length).

### D24 — 2026-06-18 — Linear-elastic dynamic verification (lattice vs fiber column)
- **Context:** The nonlinear seismic comparison (`column_dynamic.py`) tracked the fiber beam-column
  far worse than the pushover did, because hysteresis / energy dissipation / residual drift depend
  on constitutive details that differ between axial struts and a fiber section. To isolate the
  problem, take a step back and verify the *linear* dynamic equivalence first.
- **Decision (user, AskUserQuestion):** add linear-material siblings of the column pushover and
  dynamic comparisons. Both models are kept fully `Elastic` (concrete E, steel E0) on the SAME
  specimen as the nonlinear scripts:
  - reference = the SAME fiber `forceBeamColumn` section, but with Elastic fibers (reuse
    `run_beamcolumn_*` with an elastic materials triple — chosen over a plain `elasticBeamColumn`);
  - lattice = the SAME RC topology (concrete struts + longitudinal bars + stirrups, corotTruss),
    but Elastic concrete + Elastic rebar (`build_lattice_rc` gains a `rebar_material` factory so the
    rebar can be elastic; new `materials.steel_uniaxial_elastic`).
  - new scripts `examples/column_pushover_linear.py` + `column_dynamic_linear.py` (the nonlinear
    ones are untouched).
- **Rationale / finding:** with linear materials two models can agree only through mass, lateral
  stiffness (→ period) and damping. The rebar adds a fixed parallel stiffness, so K0(area) is
  monotone but NOT proportional — `calibrate_area_linear` does a short secant root-find on the FULL
  elastic lattice K0 (seeded by the plain-lattice area) instead of the single-solve scaling used
  when there is no rebar. With K0 matched the linear dynamic agrees closely: pushover lines overlap,
  and the El Centro response tracks (T1 0.515 s lat vs 0.520 s BC; peak roof drift 1.02 vs 1.00 in
  at the SAME instant t=2.35 s). Linearity also means the record scale to a target drift is exact in
  one correction (no iteration).
- **Status:** accepted. Note (not changed here): `_transient_uniform_excitation` sets the HHT
  integrator while the stale gravity Static analysis is still active, so OpenSees warns and silently
  falls back to the default Newmark integrator — pre-existing, affects the nonlinear path too; left
  as-is to avoid perturbing existing nonlinear results. The linear comparison is unaffected because
  both models go through the identical code path / integrator.

### D25 — 2026-06-18 — Sinusoidal (resonant) base excitation as the default dynamic input
- **Context:** For simplicity (and to reliably push the column into the nonlinear range in the
  nonlinear seismic file), replace the real El Centro record with a harmonic base acceleration in
  both dynamic comparisons.
- **Decision (user, AskUserQuestion):**
  - **Frequency = resonant at T1.** `a_g(t) = A·sin(2π·t/T1)`, T1 = the structure's fundamental
    period from the modal step (so the sine record is built AFTER `run_modal`, using the lattice's
    own T1). Largest response for a given amplitude; under softening the column detunes as it yields
    (realistic).
  - **Amplitude = auto-tuned to a target peak roof drift** (on the cheap fiber column, then the same
    scale drives the lattice). The sine is unit-amplitude in g and the existing `scale = G·intensity`
    machinery scales it, exactly like the record.
  - **El Centro kept selectable** via `--excitation sine|elcentro` (default `sine`); shared
    `sine_excitation` / `make_excitation` live in `column_dynamic.py`, imported by the linear file.
    Output files are suffixed by excitation (`column_dynamic_<exc>.png`).
- **Tuner (important):** the naive proportional drift-tuning loop DIVERGES for the nonlinear
  resonant case — sustained resonance ratchets a softening column to collapse (drift jumps
  1.4 in → 100s of in across a tiny intensity band), so the fixed-point collapsed the scale toward
  0 (back to elastic). Replaced with a divergence-guarded **bisection** (`tune_intensity`, shared):
  any run that is non-converged or exceeds `cap_factor·target` (or the target) is the upper bracket;
  the largest converged ≤target intensity is the answer. If the target sits past the largest stable
  drift it returns the closest stable (still nonlinear) response instead of chasing collapse. Linear
  path stays a one-shot (response ∝ intensity).
- **Result:** linear sine = clean resonance buildup, lattice within ~7% of the fiber column (lattice
  sits exactly on its own T1, the BC slightly off-resonance). Nonlinear sine = settles ~1.44 in
  (past yield), both converge, **roof drifts match (1.44 vs 1.42 in)**; base shear still differs
  (BC reaction spike ~174 kip vs lattice ~39 kip) — the known nonlinear gap to investigate next.
- **Status:** accepted.

### D26 — 2026-06-18 — Column examples refactored into an `examples/column/` package
- **Context:** The five `column_*.py` example scripts imported specimen/builder/excitation code from
  each other (e.g. `column_dynamic_linear` ← `column_dynamic` + `column_pushover` +
  `column_pushover_linear`), an increasingly tangled web.
- **Decision (user request, behaviour-preserving refactor):** group the column studies under
  `examples/column/` with three shared modules and five thin entry scripts:
  - `specimen.py` — geometry, grades, constants, `column_problem`/`zone_of`/`rebars`/`lateral_loads`;
  - `build.py` — model builders + calibration (nonlinear AND linear): `beamcolumn_reference[_linear]`,
    `rc_lattice[_linear]`, `calibrate_area[_linear]`, `lattice_k0` (was the private `_lattice_k0`);
  - `excitation.py` — `load_elcentro`, `sine_excitation`, `make_excitation`, `tune_intensity`
    (+ `G`, `RHO`, `RECORD`, `N_CYCLES`);
  - entry scripts `pushover.py`, `pushover_linear.py`, `pushover_dynamic.py`, `dynamic.py`,
    `dynamic_linear.py` (each: imports from the shared modules + its own `main()`/CLI).
  Run as e.g. `python examples/column/pushover.py` (sibling imports resolve via sys.path[0]).
- **Output tree mirrors the source:** figures now go to `examples/output/column/` (kept the same PNG
  filenames). `OUT`/`RECORD` resolve via `__file__.parent.parent` so output→`examples/output/column`
  and input→`examples/data` regardless of cwd; `mkdir(parents=True)`.
- **Verification:** pure code MOVE, no logic change. Confirmed identical output before/after — the two
  static pushovers diff IDENTICAL; pushover_linear / dynamic (sine) / dynamic_linear (sine) reproduce
  every printed metric (incl. the full `tune_intensity` bisection trace) to the digit. The five
  original `column_*.py` (and their stale root output PNGs) were then removed.
- **Status:** accepted; supersedes the file locations referenced in D24/D25 (those entries are left
  as the historical record).

### D27 — 2026-06-19 — Frame examples refactored into an `examples/frame/` package
- **Context:** Mirror the D26 column refactor for the frame examples (`frame_pushover.py`,
  `frame_pushover_rc.py`, `frame_visualize.py`).
- **Decision (user request, behaviour-preserving):** group them under `examples/frame/`:
  - `specimen.py` — the SHARED kip-in benchmark specimen (geometry constants, grades
    `FRAME_ELASTIC`/`CORE`/`COVER_C`/`BEAM_C`/`STEEL`, `frame_problem`, `lateral_loads`, `zone_of`,
    `rebars`, `OUT`). Used by both pushover scripts.
  - `pushover.py` (Stage 1 elastic) and `pushover_rc.py` (Stage 2 nonlinear RC) — thin entry scripts;
    each keeps its own `MESH` (6.0 vs 1.5) and its own calibration (Stage 1 inline at TARGET, Stage 2
    a one-step helper — deliberately NOT merged, to keep results bit-identical).
  - `visualize.py` — the SI modal study is a DIFFERENT specimen (3x4 m portal, C30, SI), so it stays
    self-contained and shares only `OUT`. No shared `build.py`/`excitation.py` (unlike the column
    package) because the frame scripts share only the specimen, not builders/calibration/input.
- **Output tree mirrors the source:** figures now go to `examples/output/frame/` (same PNG/GIF
  filenames); `OUT`/`RECORD`-style paths via `__file__.parent.parent`; `mkdir(parents=True)`.
- **Verification:** pure code MOVE. All three scripts reproduce every printed metric IDENTICALLY vs
  the pre-refactor baseline (numbers diffed; only the relocated `saved...` path lines differ):
  pushover K=271.1 kip/in, area 39.6 in²; pushover_rc lattice peak V 121.5 kip (drift 0.61,
  conv=False) vs benchmark 140.4 kip; visualize periods 0.068/0.022/0.011 s. The original
  `frame_*.py` were removed by the user.
- **Status:** accepted; supersedes the file locations in D17/D19 (left as the historical record).

### D7 — 2026-06-14 — Support both regular-grid and irregular lattice patterns
- **Why:** Regular grids are predictable and easy to validate; irregular (Delaunay) lattices
  avoid mesh-induced directional bias in cracking/fracture studies. Selectable per model.
- **Alternatives rejected:** single pattern only (less flexible).
- **Status:** accepted

### D8 — 2026-06-14 — OpenSees is an isolated backend
- **Why:** Only `rclattice/opensees/` may import `openseespy`. Keeps the internal domain
  model testable without OpenSees and allows mocking/swapping the analysis engine. OpenSees
  integer tags are assigned in the translation layer, not stored on domain objects.
- **Alternatives rejected:** sprinkling `ops.*` calls through the domain code (untestable,
  tightly coupled).
- **Status:** accepted

### D28 — 2026-06-19 — Dynamic fiber column uses `forceBeamColumn`; transient runner actually applies HHT (latent bug fixed)
- **Context:** Reviewing the nonlinear beam-column reference used to verify/calibrate the
  lattice. Two issues surfaced in `opensees.py`.
- **Decision 1 (user request) — force-based dynamic column.** `run_beamcolumn_dynamic` now
  emits `forceBeamColumn` (was `dispBeamColumn`), so ALL beam-column references — static
  cantilever, RC frame benchmark, and seismic time-history — are force-based and consistent.
  The column stays subdivided into `nelem` elements (needed so `self_mass` distributes by
  tributary length like the lattice's lumped mass; the short spans also help the element-level
  state determination converge under dynamic increments).
- **Decision 2 (bug fix) — HHT was silently disabled in every transient run.** In
  `_transient_uniform_excitation`, the integrator/analysis chain was rebuilt while the gravity
  stage's **Static** analysis object still existed. OpenSees then refused the transient
  integrator ("can't set transient integrator in static analysis") and fell back to the default
  Newmark with NO numerical damping — so the intended `HHT` was never active in any dynamic
  study (`run_dynamic`, column sine/elcentro, etc.). Fix: `ops.wipeAnalysis()` before
  re-declaring constraints/numberer/system then `integrator("HHT", 0.7)` + `analysis("Transient")`
  (wipeAnalysis clears only the analysis aggregation; the domain, gravity `loadConst`, mass, and
  Rayleigh all persist). α lowered 0.8→0.7 for slightly stronger high-frequency dissipation,
  negligible at the structural period.
- **Effect:** With HHT now genuinely applied, the El Centro nonlinear lattice-vs-fiber **roof-drift**
  agreement tightened from ~19% to ~5% (peak 1.18 vs 1.24 in; small residuals; both converged).
  All prior dynamic results shift (they had been plain-Newmark) and should be re-baselined.
- **Known-open:** HHT did NOT remove the fiber column's **base-shear** spikes (~636 kip vs the
  lattice's clean ~37 kip). Diagnosed as the stiffness-proportional Rayleigh damping force
  (`a1·K_committed·v`) spiking when fiber integration points crack/yield at high velocity — NOT
  inertia (both models share mass + input, only the fiber column spikes) and NOT integrator
  high-frequency modes (HHT didn't touch it). A clean base-shear comparison needs a damping-model
  change (mass-proportional-only or modal damping on the reference) or report-only low-pass — a
  mechanical decision deferred to the user, not yet made.
- **Status:** accepted

### D29 — 2026-06-20 — 2D plane-stress continuum as a selectable verification reference for the RC column
- **Context:** Extend the column study (`examples/column/`) with a *2D-member* reference for
  calibration/result comparison, with "an appropriate material model that matches the lattice's
  material". This is the ContinuumBuilder role (D12/D14): plane-stress `quad` (ndf=2) is the
  like-for-like continuum for the planar lattice, sharing the same structured node grid.
- **Decision (user, AskUserQuestion):**
  (1) **Material = ASDConcrete3D + PlaneStress wrapper** (the D15-planned nD material), curves built
  from the SAME `ConcreteGrade` as the struts and **crack-band regularized to the quad characteristic
  length** `lch` (Gf tension / Gfc compression) — the continuum analog of the lattice's
  length-regularized Concrete02 (D20). The "match" is necessarily at the GRADE level: a uniaxial strut
  law and a multiaxial continuum law are different objects; both consume the same physical grade.
  (2) **Reinforcement = longitudinal bars only** as steel truss struts on shared quad nodes (perfect
  bond, D5/D13) — the 2D continuum itself supplies the lateral/shear path the lattice gets from
  stirrups, so the stirrups would double-count.
  (3) **Reference is arg-selectable** like the dynamic excitations: `--reference {beamcolumn,continuum}`
  picks ONE reference; the lattice strut area is calibrated to THAT reference's K0; the plot shows
  lattice vs the chosen reference (`column_pushover_<reference>.png`). **Pushover first; dynamic later.**
- **Why the continuum is the apples-to-apples reference:** the fiber `forceBeamColumn` is 1D and has
  no analog for the lattice's diagonal-strut / 2D load-spreading (the documented ~17% overstrength;
  see `pushover.diagnose`). A 2D continuum captures the same 2D action, so lattice-vs-continuum
  isolates "is the lattice right?" from "is a 1D idealization missing 2D action?".
- **Mechanics / implementation:**
  - `materials.concrete_nd_nonlinear(grade, base_tag, wrapper_tag, *, lch, Gf, Gfc)` → (ASDConcrete3D,
    PlaneStress) pair. Tension `-Te/-Ts/-Td`: elastic to ft at eps_cr=ft/E, linear softening to a small
    residual at eps_tu=eps_cr+2·Gf/(ft·lch). Compression `-Ce/-Cs/-Cd`: elastic→peak fc at epsc0→soften
    to fcu at eps_cu=max(epsc0+2·Gfc/((fc+fcu)·lch), epsU). Damage d=clip(1−σ/(E·ε),0,0.95) (isotropic-
    damage decomposition). **Verified on a single-quad coupon:** reproduces ft and fc; ASDConcrete3D
    follows the given (Te,Ts) curve with NO hidden auto-regularization, so the pre-regularization is
    the only regularization (no double-counting) as long as the quads are ≈ `lch`.
  - `builders.build_continuum_rc(problem, mesh_size, *, nd_material_for, zone_of, rebars, plane, …)`:
    per-zone nonlinear quads (one nD-material pair per zone — the structured grid has a single quad
    size, so no per-length split as in the lattice) + longitudinal steel struts on shared nodes. Same
    grid / mass / supports plumbing as `build_lattice_rc`. nDMaterial and uniaxialMaterial tags are in
    separate namespaces.
  - `examples/column/build.py`: `continuum_reference()` + `make_reference(name)`; `specimen.py` gains
    `longitudinal_rebars()` (and `rebars()` stays bit-identical). `pushover.py`: `--reference` flag.
- **Result (kip, in):** continuum K0=80.3 kip/in, peak V=39.1 kip (smooth plateau, reached 9.28 in);
  lattice calibrated to that K0 (strut area 16.26 in²) peaks V=39.8 kip (reached 5.50 in) — **initial
  stiffness matched, peak shear within ~2%**, post-peak plateau within ~10–15% (lattice settles ~34 kip,
  noisier). Both `conv=False` past their limit points (normal softening RC). For contrast the 1D
  beam-column reference is K0=61.5, peak 31.4 — the lattice and continuum agree with EACH OTHER far
  better than either does with the 1D beam-column, confirming the lattice "overstrength" is genuine 2D
  action, not an artifact. (A single-quad coupon check is what caught an early compression-curve
  overshoot — the elastic point must reach 0.4·fc at strain 0.4·fc/E, not 0.4·epsc0.)
- **Cost / known-open:** the continuum is heavy (1536 ASDConcrete3D quads, ~6 s/step) → a coarser
  `dU=0.1`; ASDConcrete3D `-implex` would speed/robustify it but introduces error and needs controlled
  steps — left off (a future option). The reference K0 is the first-step secant (like the beam-column),
  so it is `dU`-sensitive because concrete cracks at very low drift. Pre-existing, unrelated:
  `tests/test_rc.py::test_nonlinear_pushover_runs_and_yields` fails on the clean baseline too (portal-
  frame `build_lattice_rc` path, untouched here).
- **Status:** accepted (pushover). A dynamic continuum reference (the seismic siblings) is deferred to
  next, per the user's staging.

### D30 — 2026-06-21 — Dynamic continuum reference; ASDConcrete3D hysteresis matched to Concrete02 by PURE DAMAGE
- **Context:** Extend D29 to the seismic time-history (`dynamic.py`). For a monotonic pushover only the
  backbone matters (D29 matched it at the grade level), but a DYNAMIC run is governed by the cyclic
  HYSTERESIS — and the lattice strut law (uniaxial Concrete02, Kent-Park-Yassin unloading) and the
  continuum law (ASDConcrete3D, a continuum damage model) are different constitutive theories that
  cannot be made identical. The question (user): how to configure the nD material so its hysteresis
  matches Concrete02's.
- **Finding (single-quad cyclic coupon vs a Concrete02 truss):** the ASDConcrete3D damage value sets
  the damage↔plasticity split — it back-computes plastic strain as `eps_p = eps − σ/((1−d)·E)`. With
  `d = 1 − σ/(E·ε)` (PURE isotropic damage, unload toward the origin, no residual) the loops match
  Concrete02 best: dissipated energy ~108%, comparable residual strain. Introducing plasticity (lower
  `d`, residual strain) OVERSHOOTS — ~134% dissipation and larger residual at +30%. The unloading
  SHAPE still differs (ASD unloads ~linearly down the damaged secant; Concrete02 curves), an
  irreducible theory difference, but the energetics (what governs seismic response) are close.
- **Decision (user, AskUserQuestion):** (1) **pure damage** is the dynamic-hysteresis configuration —
  exposed as `concrete_nd_nonlinear(..., plastic_frac=0.0)` (documented; >0 scales damage down to add
  plasticity). The monotonic backbone is independent of `plastic_frac`, so the D29 pushover is
  unchanged. (2) **Wire the continuum dynamic reference, sine only** — a full continuum time-history is
  heavy (~1.5 s/step over 1536 ASDConcrete3D quads; the lattice is ~1.35 s/step too — 7902 corotTruss
  struts are not cheap), so El Centro (~3000+ steps) is impractical and guarded off.
- **Mechanics / implementation:**
  - `materials.concrete_nd_nonlinear` gains `plastic_frac` (scales the damage table by `1−plastic_frac`).
  - `build.py`: `_continuum_model` factored out; `continuum_k0()` (one tiny continuum pushover step) and
    `continuum_dynamic(accel, dt_record, scale, top_mass, dt, …)` (bakes the axial-load tributary mass
    on the top nodes, then `run_dynamic`).
  - `dynamic.py`: `--reference {beamcolumn,continuum}`. For the continuum the lattice is calibrated to
    the **continuum K0** (so lattice and continuum share mass+stiffness+T1 → the comparison isolates the
    Concrete02-vs-ASD hysteresis), the sine cycles are capped (`CONTINUUM_CYCLES=8`), and the intensity
    is still tuned on the **cheap fiber column** (the only affordable proxy; its K0 differs from the
    matched lattice/continuum so the achieved drift is approximate — fine, both get the identical scaled
    record). El Centro + continuum is rejected at the CLI. Output `column_dynamic_<reference>_<exc>.png`
    (+ `_hyst.png`) — note this renames the beam-column dynamic outputs too (was `column_dynamic_<exc>`).
- **Result (kip, in; resonant sine at T1=0.452 s, 8 cycles, scale 0.168×, both `conv=True` through all
  8 cycles):** roof-drift and base-shear **histories track closely**; peak drift lattice 1.12 vs
  continuum 0.85 in, residual +0.26 vs +0.31 in, peak shear 43.3 vs 39.2 kip (~10–25%). The hysteresis
  **envelopes coincide** (both ±40 kip) but the **loop shapes differ** — the continuum's are fuller/
  rounder, the lattice's more pinched — exactly the coupon's irreducible Concrete02-vs-damage difference.
- **Status:** accepted. Same pre-existing unrelated test failure as D29 noted; 25 unit tests pass.
  Possible follow-ups (user): tune the continuum intensity on a continuum-matched cheap proxy for an
  exact drift target; ASDConcrete3D `-implex` to cut the ~1.5 s/step.

### D32 — 2026-06-22 — Opt-in analysis-model drawings (`--draw`) with reinforcement styled by kind
- **Decision:** every analysis entry script (column + frame) gains an OPT-IN `--draw` flag
  (`action="store_true"`, OFF by default) that, in addition to the result plots, saves a drawing of
  the analysis model(s) it builds. In the lattice drawing, reinforcement is styled distinctly from the
  concrete skeleton — a different colour AND a heavier lineweight — and **longitudinal bars vs
  stirrups/ties get different colours** (per user instruction). Default by off so the existing fast runs
  are unchanged.
- **Mechanism:** a backend-agnostic per-element drawing tag `Element.kind` (NOT seen by the OpenSees
  layer): "concrete" (default strut), "quad" (continuum), or a reinforcement role. `Rebar` gains
  `role` ("longitudinal" default / "stirrup"), threaded by the builders into the rebar struts' `kind`
  (`build_lattice_rc`, `build_continuum_rc`); the column specimen tags its hoop ties `role="stirrup"`.
  `viz.figure_model(panels)` / `draw_model_kinds` render the undeformed model, grouping line elements by
  `kind` (concrete grey/thin; longitudinal C3/red; stirrup C2/green; generic rebar C1 — rebar weights
  kept light so they don't swamp the skeleton) with a legend, and continuum quads as light polygons.
  `figure_model` lays out ONE ROW PER MODEL × THREE COLUMNS — concrete skeleton only / reinforcement
  only / combined (`draw_model_kinds(which=...)`), with a shared per-model extent so the panels align.
  Each script saves `<result-stem>_model.png`; the heavy dynamic scripts draw up front (before the
  time-history). For the continuum-reference column runs the continuum model is also drawn (via
  `build._continuum_model`).
- **Output layout:** figures are organised `output/<modeltype>/<analysis>/<reference>/<file>` — e.g.
  `output/column/dynamic/beamcolumn/…`, `output/column/pushover/continuum/…`,
  `output/frame/pushover_rc/beamcolumn/…`. `OUT` (from each `specimen.py`) is the model-type dir;
  each script sets `outdir = OUT / <analysis> / <reference>` and `mkdir(parents=True)`. The reference
  is the CLI `--reference` for the column pushover/dynamic, and fixed otherwise (linear siblings +
  `pushover_rc` → `beamcolumn`; elastic frame `pushover`/`visualize` → `continuum`). Filenames are
  unchanged (args still encoded in the name), so the reference appears in both path and filename.
- **Status:** accepted. Additive only — `kind`/`role` carry defaults and do not change the FE assembly
  numerically (test suite numerically unchanged: 25 pass + the same one pre-existing nonlinear-pushover
  failure from the in-progress materials WIP). Verified end-to-end on `column/pushover_linear.py --draw`.

### D34 — 2026-06-24 — Frame examples rebuilt as the cantilever column + a thinner beam; full column-package mirror
- **Context (user request):** make the frame examples mirror the column package (D26/D27 structure) —
  "whatever we have for the cantilever column, the portal frame should have too" — and **constitute the
  frame from the cantilever column** (`examples/column`) plus a **thinner connecting beam** that uses
  the **same concrete grades and steel grade**. Supersedes the D27 frame example file set.
- **Geometry / reinforcement (user, AskUserQuestion):** one-bay one-storey portal frame, kip-in.
  Columns = the EXACT cantilever column (24-deep × 144-tall, THK 15, confined CORE + unconfined
  COVER_C Concrete02, 3+2+3 Steel02 longitudinal + horizontal stirrups). Beam = thinner **18-in**
  deep × 15 wide, **span 144 in** (centerline-to-centerline), SAME CORE/COVER + STEEL grades;
  **top+bottom 1.8 in² longitudinal layers + vertical stirrup ties** along the clear span. Fixed bases,
  gravity P=180 kip/column at the joints, lateral pattern over the beam level.
- **Mirror (scope = full mirror + replace old):** `examples/frame/` now holds `specimen.py` (shared),
  `build.py` (references {beamcolumn fiber frame, 2D continuum frame} + scalar/2-group calibration +
  lattice + linear variants), `excitation.py` (self-contained copy of the column's), and the four entry
  scripts `pushover.py` / `pushover_linear.py` / `dynamic.py` / `dynamic_linear.py` (selectable
  `--reference` exactly like the column). The old Stage-1 `pushover.py` and Stage-2 `pushover_rc.py`
  were removed; the separate SI modal `visualize.py` is kept as-is.
- **New backend (`opensees.py`, the only `ops.*` module):** `_rc_beam_fiber_section` + `_beam_section`,
  `run_beamcolumn_frame` (pushover) and `run_beamcolumn_frame_dynamic` (members subdivided so self-mass
  distributes like the lattice) — the fiber `forceBeamColumn` portal-frame reference (the frame analog
  of `run_beamcolumn_cantilever` / `_dynamic`). The columns reuse the 15×24 `_rc_fiber_section`; the
  beam gets a 15×18 top/bottom-bar fiber section sharing the same core/cover/steel materials. (The old
  hard-coded `run_benchmark_rc_frame` is left in place but is now unused by the examples.)
- **Beam concrete model (user, AskUserQuestion) — the key fork:** the thin **nonlinear** beam softens
  into a **local lattice mechanism the static Newton pushover cannot trace** (dies ~0.2 in / ~50 kip
  regardless of grade ductility CORE-vs-COVER or `horizon` — the SAME phenomenon as the pre-existing
  failing `test_nonlinear_pushover_runs_and_yields`). The **dynamic** transient is fine with the fully
  nonlinear beam (inertia + HHT regularize the softening — verified it rides past the static death
  point). Decision: **DEFAULT = elastic beam concrete** (columns + ALL rebar nonlinear; beam concrete
  linear-elastic at the same grade E), with an **opt-in `--nonlinear-beam`** for the SAME Concrete02
  beam, in **BOTH pushover and dynamic**. Threaded via distinct beam zones (`zone_of` → "beam_core" /
  "beam_cover") into `make_material_for` (lattice) and `_beam_fiber_materials` (fiber reference, emitted
  as separate mat tags 4/5 keeping shared Steel02). The 2D **continuum** reference keeps the beam
  nonlinear always (the plane-stress quads are stable). The default elastic beam matches the project's
  prior frame precedent (D19 fork-2) and the OpenSees RCFrame benchmark.
- **Results (kip, in; default elastic beam, scalar calibration to the fiber frame):** `pushover_linear`
  K0 match −0.03% (lattice 301.16 vs fiber 301.26 kip/in). `pushover` (nonlinear columns, elastic beam)
  lattice peak V **133.4 vs fiber-frame 131.8 kip (~1.2%)**, lattice `conv=False` past its limit point
  at 0.79 in (column-base sway mechanism, same conv-past-limit behaviour as the column lattice).
  `dynamic_linear` (resonant sine) T1 0.359 vs 0.353 s, peak drift 1.03 vs 1.00 in, peak V 297 vs
  302 kip (~1.6%), both `conv=True`. Nonlinear `dynamic` runs and converges with both the default
  elastic beam and `--nonlinear-beam` (transient is stable either way).
- **Status:** accepted. Test suite unchanged (25 pass + the same one pre-existing nonlinear-pushover
  failure, which is itself this thin-nonlinear-beam-frame instability). Supersedes the frame file set
  in D27 (left as the historical record).

---

## Reconstructed entries (logged 2026-06-26 during a docs↔code audit)

These decisions are **referenced in the source code** (and were implemented) but had **no entry in
this log**. They are reconstructed from the code and its comments; the *rationale* below is inferred,
so **confirm/correct** it. (Also: **D21 is referenced nowhere** in the code or docs — either an
unlogged decision or a skipped number; please confirm.)

### D22 — ~2026-06-17 — Keep the lattice convergent past the limit point: `corotTruss` + residual compression plateau + dynamic-relaxation pushover *(reconstructed)*
- **Decision:** three coupled measures so the RC lattice can be pushed well past first yield without
  the axial-truss mechanism terminating the run (the D4 caveat; the D20 length-regularization alone
  was not enough):
  1. struts use **`corotTruss`** (not `Truss`) so geometry/P-Δ is carried consistently at large drift
     (the P-Δ effect the beam-column reference shows). `build_lattice_rc(strut_element=...)`, set to
     `"corotTruss"` by the column/frame builders.
  2. a **residual compression plateau**: `concrete_uniaxial_regularized(residual_ratio=0.2)` floors
     `fpcu` at `residual_ratio*fc`, so a crushed strut holds a positive flat residual instead of
     dropping to zero tangent — removing the local zero-stiffness mechanism that otherwise ends the
     pushover just past yield (`materials.py`).
  3. **`run_pushover_dynamic`** (user-selected): a quasi-static *transient* (Newmark + heavy Rayleigh,
     imposed ramp on the drive nodes) that rides through limit points / local snap-backs the static
     Newton can't (`opensees.py`).
- **Why (inferred):** sustaining post-yield *drift* (not strength) was the binding problem; these
  three address geometry, the crushed-strut zero-tangent, and solver robustness respectively.
- **Status:** accepted (reconstructed). `corotTruss` + the residual plateau are the defaults; the
  dynamic-relaxation pushover is an opt-in runner.

### D23 — ~2026-06-17 — Confinement: confined-core vs unconfined-cover grades + transverse stirrup/tie struts *(reconstructed)*
- **Decision:** model transverse confinement explicitly in the RC specimens (column and frame):
  - two concrete grades — a **confined `CORE`** (high crushing strain `epsU` + residual `fcu`,
    Mander-style) inside an **unconfined `COVER`** (no residual) — selected per `zone_of(x,y)`;
  - **transverse stirrup/tie steel struts** (`STIRRUP_AREA`, one horizontal tie per mesh row across
    the core), `role="stirrup"`, giving the lattice a **non-softening lateral/shear path** so the
    confined core does not disintegrate into a local mechanism past yield.
  (`examples/column/specimen.py`, `examples/frame/specimen.py`.)
- **Why (inferred):** the unconfined softening law alone lets the core crush into a mechanism; the
  ties + ductile core grade reproduce confinement so the lattice sustains post-yield drift. The 2D
  continuum reference omits stirrups (it supplies that lateral path itself — D29).
- **Status:** accepted (reconstructed).

### D31 — ~2026-06-21 — `horizon` is a tunable knob for redundant bracing against the post-peak mechanism *(reconstructed)*
- **Decision:** expose the strut-connectivity **`horizon`** (D9) as a study parameter (CLI `--horizon`,
  default unchanged at **1.5**). A larger horizon connects each node to more neighbours → more
  redundant triangulation/bracing, which helps the lattice resist forming the post-peak local
  mechanism. Threaded through `calibrate_area` / `rc_lattice` / `build_lattice_rc`.
- **Why (inferred):** a tuning lever, complementary to D22, for the same post-yield-stability problem;
  kept at 1.5 by default so existing results are unchanged.
- **Status:** accepted (reconstructed).

### D33 — ~2026-06-23 — Transient runners use modal damping (replacing stiffness-proportional Rayleigh); resolves the D28 base-shear spike *(reconstructed)*
- **Decision:** `_transient_uniform_excitation` (and the `single_beamcolumn` run) now assign **modal
  damping** (`ops.modalDamping(zeta)` over the first `max(modes)` mass-orthonormal eigenvectors on the
  gravity-held tangent) instead of stiffness-proportional Rayleigh. HHT-α still dissipates the
  uncomputed higher-mode content.
- **Why:** stiffness-proportional Rayleigh (`a1·K_committed·v`) **spiked the base reaction** whenever
  fiber/strut integration points cracked or yielded at high velocity — the **D28 known-open**
  base-shear-spike artifact (fiber column ~636 kip vs lattice ~37 kip). Modal damping has no term
  riding the committed tangent, so the spike disappears and lattice↔reference base-shear histories
  become comparable. **This resolves the mechanical decision D28 left open.**
- **Status:** accepted (reconstructed). Supersedes the Rayleigh damping noted in D24/D28 (those
  entries left as the historical record).

### D35 — 2026-06-26 — Calibration output: first-N mode-shape figure (reference vs lattice) + periods table
- **Decision:** every calibrating run now emits a **modal calibration figure** as calibration output —
  the first `N` (default 3, the existing `n_modes`) **mode shapes** drawn for the SELECTED reference
  (top row) and the lattice (bottom row), columns = modes 1..N, with a **periods table underneath**
  (per mode: `T_ref`, `T_lattice`, `Δ vs reference %`). One figure per run; **no CSV** (the figure
  carries everything — table + a caption with the calibration summary). Scope = **always**: scalar +
  groups, nonlinear + linear, column + frame (8 entry scripts).
- **Methodology (settled with the user, in order):**
  - *Content:* reference + lattice + per-mode Δ. **No MAC** — shapes are compared visually in the
    panels; periods quantitatively in the table. Pairing is **ascending order** (so for a cantilever
    the 3rd "mode" is typically the axial mode in BOTH models — visible and self-explanatory).
  - *Reference = whatever `--reference` selected* ("report whatever is used as reference"): the 2D RC
    continuum (`run_modal` on `_continuum_model`) or the **subdivided fiber** beam-column / portal
    frame. The single-element pushover reference can't give a real bending mode, so two new modal
    runners mirror the *dynamic* discretization: **`run_beamcolumn_modal`** (cantilever) and
    **`run_beamcolumn_frame_modal`** (portal frame) — nodes at the lattice row/column heights, lumped
    tributary mass, **no gravity / no top mass** → eigen on the initial-tangent stiffness.
  - *Fit target stays the continuum (D16, user's call).* So the table column is a plain **`Δ vs
    reference`**, NOT a calibration residual; the calibration *quality* (fitted areas / RMS / success,
    or scalar area + K0 match) lives in the figure **caption**. This also sidesteps the plain-continuum
    (groups fit target) vs RC-continuum (`--reference continuum`) mismatch honestly.
  - *Mass consistency* (required for comparable periods, D16): equal total mass — the continuum shares
    the builder tributary mass by construction; the beam-column is given the lattice's total self-mass
    (column) / the same geometry·ρ per-member self-mass (frame). The lattice modal uses the **as-built
    self-mass** (the dynamic scripts call it BEFORE adding the seismic top mass `P/g`), matching the
    modal basis the calibration itself uses.
- **Mechanics added:** `opensees.run_beamcolumn_modal` / `run_beamcolumn_frame_modal` (eigen +
  per-node eigenvectors + a lightweight line `Model` for drawing); `viz.figure_modal_calibration`
  (gridspec mode-shape grid + a `matplotlib` table axes + caption) and `viz.sign_fix` (deterministic
  per-mode sign — reference and lattice are on different meshes, so `align_sign`'s shared-node-id
  cross-alignment can't be used); `modal_calibration_figure` in `examples/column/build.py` and
  `examples/frame/build.py`, wired into all pushover/dynamic (+ linear) entry scripts.
- **Result (column):** T1 tracks the K0 calibration (~4% vs beam-column scalar; ~6% vs continuum);
  higher modes diverge more for the beam-column (axial-lattice limitation, D16) and much less for the
  continuum (~5%). Frame: sway T1 ~5%. The known higher-mode divergence is now *shown*, not hidden.
- **Status:** implemented + verified (column pushover beamcolumn/continuum, column/frame linear,
  dynamic in-flow). Backend independence preserved (figure I/O lives in the example layer; `viz`
  imports no openseespy).

### D36 — 2026-06-26 — Column beam-column reference is a SINGLE `forceBeamColumn` element everywhere (dynamic + modal), not subdivided
- **Decision:** in the `examples/column/` scripts, when `--reference beamcolumn`, the column is
  represented by **one** force-based fiber `forceBeamColumn` element (fixed base node 1, free top
  node 2, 5 Lobatto points) across **all** analyses — pushover, dynamic, and the D35 modal figure.
  The pushover reference (`run_beamcolumn_cantilever`) was already single-element; this brings the
  *dynamic* (`run_beamcolumn_dynamic`) and *modal* (`run_beamcolumn_modal`) references into line.
  Self-mass is now carried as **distributed element mass** (`-mass self_mass/height`) and the
  axial-load seismic mass `P/g` is lumped at the top — identical to `single_beamcolumn.py`.
- **Why:** the user wants the column reference idealized as a single member (as in
  `single_beamcolumn.py`), not a ~96-element stick. A single force-based element with 5 Lobatto
  points already captures the fiber-section nonlinearity; the subdivision was only there to spread
  lumped nodal mass and shorten per-element spans. `single_beamcolumn.py` proves the single element
  converges under dynamic increments (same `_transient_uniform_excitation` sub-stepping/algorithm
  fallback), so the subdivision is unnecessary.
- **Supersedes (partially):** the column half of D30/D35's "subdivide the fiber reference to the
  lattice row heights." The **frame** references (`run_beamcolumn_frame_dynamic` /
  `run_beamcolumn_frame_modal`) are unchanged — still subdivided.
- **Trade-off accepted (with the user):** a single 2-node element resolves only the **first
  flexural mode** well. In the D35 modal figure the beam-column's modes 2–3 are now axial/rotational
  artifacts (mode 3's massless rotational DOF gives `T→0`). The user chose "single element
  everywhere" knowing this; `N_MODES` stays 3 (the continuum reference still resolves 3 real modes).
  Requesting `≥` the available DOFs makes ARPACK fall back to `fullGenLapack` (instant on the tiny
  system; the existing `try/except` already handles it — only cosmetic OpenSees stderr).
- **Mechanics:** `run_beamcolumn_dynamic` / `run_beamcolumn_modal` lost their `nelem` parameter and
  the nodal-mass distribution loop; callers (`dynamic.py`, `dynamic_linear.py`, `build.py`) drop the
  `nelem` argument (and the now-unused `MESH` import). The modal runner's drawing `Model` is now two
  nodes + one segment.
- **Status:** implemented + verified (single-element seismic reference converges; modal figure
  renders; test suite unchanged at 25 pass + the one pre-existing known-fail). Backend independence
  preserved.

### D37 — 2026-07-14 — Plain-concrete cube uniaxial-tension study (`examples/cube/`): lattice vs a single displacement-based beam-column
- **Context (user request):** the simplest possible lattice check — a 0.1 m plain-concrete cube
  (Concrete02) modelled two ways and compared under UNIAXIAL TENSION (a material coupon), mirroring
  the structure of `examples/column/` but much thinner.
- **Decisions (user, AskUserQuestion):**
  (1) **2D plane-stress square represents the cube** — a 100×100 mm face with 100 mm out-of-plane
  thickness (the section perpendicular to the pull is L×THK = the cube face). The lattice/continuum
  mesher is 2D-only (gmsh plane surface); a true 3D lattice was rejected as a large new capability
  (no 3D mesher/builders) out of scope for "similar to examples/column". Uniaxial tension is fully
  captured by the 2D model.
  (2) **Fully fixed base** — the bottom edge is fixed in X and Y (like the column), the top edge is
  pulled up (+Y) under DisplacementControl. (The frictionless-platen alternative — a purer 1D
  uniaxial stress — was offered and declined.)
  (3) **Plain Concrete02 on both sides** — NO fracture-energy regularization (D20) and no residual
  plateau/corotTruss/dynamic-relaxation (D22). The calibration target is the INITIAL STIFFNESS only;
  the post-peak softening is left plain (so the truss lattice snaps at the crack, as expected).
- **Reference = a single DISPLACEMENT-based fiber `dispBeamColumn`** (new `opensees.run_dispbeamcolumn_tension`),
  NOT the force-based `forceBeamColumn` used elsewhere. A displacement-based element interpolates the
  axial displacement linearly → CONSTANT axial strain along the element → the whole element carries
  eps=u/L uniformly with no integration-point localization, so the axial response is exactly
  A·σ(eps): the clean 1D coupon. Fiber section = the cube face (one plain-concrete material), fixed
  base, axial DisplacementControl pull, Linear geomTransf (small strain, no P-Δ).
- **Units:** SI in **N–mm** (edge = 100 mm, E = 30000 N/mm² = 30 GPa, stresses in MPa) so cracking
  displacements are O(0.01 mm) — well conditioned for the default NormDispIncr tolerances (pure-metre
  units make them ~1e-5 and fight the solver). Grade: generic C30 (fc 30, ft 3, epsc0 0.002, fcu 6,
  epsU 0.006).
- **Files:** `examples/cube/{specimen.py,build.py,pushover.py}` (thin package, sibling imports like
  the column/frame packages; output → `examples/output/cube/pushover/`). `calibrate_area` matches the
  elastic lattice's axial K0 to the beam-column's (K linear in strut area → one tiny tension solve);
  `cube_lattice` builds the plain-Concrete02 strut lattice (`build_lattice_rc`, one trivial zone, no
  rebar, `Truss` struts). Output is the coupon **stress–strain** curve (σ = base reaction / L·THK,
  ε = top disp / L).
- **Result (N, mm; MESH 10 → 121 nodes, 420 struts):** dispBeamColumn K0 = 3.0e6 N/mm (= E·A/L),
  peak = 30000 N = **3.00 MPa = ft**, full softening traced to 0.08 mm. Lattice: strut area calibrated
  to **631.2 mm²** → K0 matched exactly (elastic branches overlap); peak **2.88 MPa (~4% below ft)**,
  then `conv=False` at the cracking strain (u≈0.01 mm) — the plain-Concrete02 truss lattice snaps at
  the crack band (displacement-controlled snapback), the expected consequence of decision (3).
- **Status:** accepted. Additive backend change (new `run_dispbeamcolumn_tension`); test suite
  unchanged (25 pass + the same one pre-existing known-fail). Follow-up available if the lattice
  descending branch is wanted: D20 regularization + D22 residual plateau / corotTruss / dynamic-
  relaxation pushover (deliberately off here to honour the plain-Concrete02 choice). See D38 (the
  descending-branch follow-up: arc-length integrated, but the plain-Concrete02 crack is a mechanism,
  so only dynamic relaxation traces it).

### D38 — 2026-07-14 — Cube descending branch: arc-length integrated, but the plain-Concrete02 tension crack is a MECHANISM (singular), not a snapback → only dynamic relaxation traces it
- **Context (user request):** "integrate arclength properly" so the cube lattice (D37) can follow the
  DESCENDING branch of its uniaxial-tension force–deformation curve — which plain DisplacementControl
  cannot (it dies at the peak, D37).
- **What was built:** `opensees.run_pushover_arclength` — a proper cylindrical arc-length (Crisfield,
  `alpha=0` so the load term drops from the constraint) pushover: a DisplacementControl PROBE step
  fixes a well-scaled arc length (= the probe's displacement-increment norm) and the loading
  direction, then `integrator ArcLength` marches with adaptive arc reduction + stronger-algorithm
  fallback, stopping on a displacement target / a base-shear-fraction floor / a failed step. Same
  `run_pushover`-style interface (gravity → lateral pattern → base-shear recording). Wired into
  `examples/cube/pushover.py` as `--solver {static,arclength,dynamic}` (default `static`).
- **KEY FINDING (empirically established, evidence in-session):** arc-length does NOT trace this
  branch, and the reason is fundamental, not a tuning failure. Tests on the calibrated cube lattice:
  - global arc-length (`ArcLength` alpha=0/1, `MinUnbalDispNorm`): **oscillates/stalls at the peak**;
  - single-node indirect displacement control (control a crack-band node), **± a triggering
    imperfection** (weaken the base band): **fails at the peak**;
  - even with fracture-energy REGULARIZATION (gentler softening, D20), arc-length and DisplacementControl
    **fail at the SAME step** (the delayed limit point) — arc-length gains nothing.
  Root cause: plain **Concrete02 tension softens to ZERO stiffness**, so once a horizontal band cracks
  the block above it is a **kinematic mechanism** and the tangent stiffness is **SINGULAR** (rank
  deficient), not merely negative. A snapback (negative-definite but non-singular tangent) is what
  arc-length is for; a singular tangent / mechanism is not something ANY static path-follower can
  cross. This is the D4 truss-lattice-mechanism caveat, in tension (where — unlike D22's compression
  residual plateau — Concrete02 leaves no residual stiffness floor).
- **What DOES trace it — dynamic relaxation** (`run_pushover_dynamic`, D22): a quasi-static transient
  whose inertia keeps the system well-posed THROUGH the mechanism. On the cube it traces the full
  pre- + post-peak branch, `conv=True` to u=0.12 mm: elastic rise (K0 matched) → peak (dynamic
  OVERSHOOT to ~3.5 MPa vs ft 3.0, from inertia + the finite ramp rate + diagonal redistribution) →
  the crack snap → damped ringing about zero → zero (fully cracked). The snap is genuinely dynamic
  (a brittle unreinforced tension crack IS an unstable event — why real direct-tension tests need
  servo CMOD control), so the post-peak carries dynamic artifacts; that is honest, not a bug.
- **Decision:** keep `run_pushover_arclength` (a correct, general, reusable snapback solver — it will
  matter for problems whose softening retains a non-singular tangent, e.g. residual-plateau /
  reinforced cases), expose all three solvers, and **document that `--solver dynamic` is the one that
  follows the cube's descending branch**. Did NOT silently regularize or add a tension residual to
  force arc-length through — that would change the D37 plain-Concrete02 physics the user chose.
- **Options recorded for a future turn (if a cleaner STATIC descending branch is wanted):** (a) give
  the cracked struts a small residual tension stiffness (parallel Elastic strut / a tension floor) so
  the tangent stays non-singular → arc-length then traces a real snapback; (b) fracture-energy
  regularization to reduce brittleness (still a mechanism at full crack, but a longer stable branch
  first); (c) a dedicated local/energy-dissipation arc-length (crack-opening-constrained), not in
  OpenSees's built-ins.
- **Status:** accepted. Additive (`run_pushover_arclength`; removed a stray `print("Step …")` debug
  line in `run_pushover_dynamic`); test suite unchanged (25 pass + the same one pre-existing known-fail).

### D39 — 2026-07-14 — Cube: the lattice FOLLOWS the reference once its struts are fracture-energy-regularized (localization fix); default flipped to regularized + dynamic
- **Context (user):** "I can't see that lattice model follow[s the] reference model." On the D37/D38
  figure the two descending branches diverge completely — the reference (single element) softens
  gently from ft to 0 over global strain 0→1.1e-3, but the lattice drops near-vertically at ~3.5e-4.
- **Diagnosis (three distinct causes, established with in-session evidence):**
  1. **Localization / mesh-size (the visible "doesn't follow"):** the single-element reference smears
     its softening over the WHOLE 100 mm element; the lattice crack localizes into ONE ~10 mm strut
     row, so it reaches zero stress at ~1/10th the global strain. This is the crack-band problem, and
     the fix is fracture-energy regularization (D20).
  2. **Mechanism / solver (D38):** the cracked band → 0 stiffness → singular tangent, so only dynamic
     relaxation traces the branch.
  3. **Diagonal overstrength (strength, NOT stiffness):** matching K0 does NOT match tensile
     strength — the diagonal struts add a parallel tensile load path, so the lattice PEAK overshoots
     ft (≈4.5 MPa regularized vs 3.0; slowing the dynamic loading to 60 periods did NOT reduce it →
     it is real overstrength, not inertia). **User decision: "You don't need to match peaks, it's
     ok"** — so cause 3 is left uncorrected (a strut-strength calibration is deferred).
- **Decision:** regularize the lattice struts to the reference's IMPLIED fracture energy so the two
  descending branches align. The reference uses plain Concrete02 with default `Ets = 0.1*E` smeared
  over L, i.e. `Gf_ref = ft^2 * L / (2 * 0.1*E)` (= 0.15 N/mm for this grade); the lattice struts are
  crack-band regularized (`concrete_uniaxial_regularized`, D20) to that SAME `Gf`, so the crack
  opening at failure `w = 2Gf/ft` is band-width-independent and both models reach zero at the same
  global strain (≈1.1e-3). Verified: the regularized lattice, traced with dynamic relaxation, now
  overlaps the reference elastic branch (K0), softens over the same strain range, and is a CLEAN
  curve (gentle softening ⇒ the dynamic snap/ringing of the plain case is gone). Only the peak sits
  higher (cause 3, accepted).
- **Defaults flipped (`examples/cube/pushover.py`):** the lattice is now **regularized by default**
  (`--plain` recovers the D37 K0-only baseline) and the default solver is **`dynamic`** (the only one
  that follows the branch; `static`/`arclength` still offered and still die at the mechanism, D38).
  `build.cube_lattice(area, *, regularize=True, Gf=REF_GF)`; `REF_GF` computed from the grade + L.
  Output stem `cube_tension_<solver>[_plain].png`.
- **Status:** accepted. Example-layer change only (`build.py`, `pushover.py`); backend untouched since
  D38; test suite unchanged (25 pass + the same one pre-existing known-fail). Deferred (user may want
  later): a strut-tensile-strength calibration so the lattice PEAK also matches (cause 3).

### D40 — 2026-07-14 — Cube peak-strength correction via strut-area reduction by a saved `PEAK_RATIO` (matches peak, sacrifices K0)
- **Context (user):** address the D39 diagonal-overstrength peak gap by "sav[ing the] peak mismatch
  ratio as a constant, and use it to modify [the] lattice area by reducing the area [by] the same
  ratio as the peak mismatch."
- **Decision:** add `build.PEAK_RATIO = 1.50` (the measured lattice_peak/reference_peak ≈ 4.50/3.00
  on the default regularized+dynamic run) and have `cube_lattice(area, *, correct_peak=True)` divide
  the K0-calibrated strut area by `PEAK_RATIO` before building. The lattice peak scales linearly with
  area, so peak → 4.50/1.50 = **3.00 MPa = reference** (verified: 29997 N ≈ 30000 N).
- **Explicit trade-off (flagged to the user, then done as asked):** strut area scales stiffness AND
  strength together, so a single area scalar CANNOT match both. Reducing area to fix the peak also
  **lowers K0 by the same factor** → the lattice elastic branch is now softer (effective E ≈ E/1.5)
  and peaks at strain ≈4e-4 instead of the reference's 1e-4; the elastic branches no longer overlap.
  This is the D39 stiffness match traded for a strength match — the user's call.
- **Caveat:** `PEAK_RATIO` is the REGULARIZED lattice's overstrength; for `--plain` (overstrength
  ≈1.17) it over-corrects (peak → 2.34 MPa). `--plain` is only the baseline, so left as-is (documented
  on the constant). `correct_peak=False` disables the reduction.
- **The other option (NOT taken, recorded):** reduce the strut TENSILE STRENGTH (`ft`) instead of the
  area — that lowers the peak while keeping the area, so K0 stays matched AND the peak matches (the two
  are then calibrated by two independent knobs: area→K0, ft→strength). Offer for a later turn if the
  user wants both matched.
- **Status:** accepted (default later flipped OFF by D41 — the correction is now opt-in). Example-layer
  only (`build.py`, `pushover.py`); backend untouched; test suite unchanged (25 pass + the same one
  pre-existing known-fail).

### D41 — 2026-07-14 — `PEAK_RATIO` peak correction made a THREE-WAY selectable mode (default off; area OR tensile-strength)
- **Context (user):** "Keep [the] peak mismatch ratio, but make [it] use[d] as follows: by default it
  is not used; if intended to modify area, modify the area; if intended to modify tensile [strength],
  modify [the] material's tensile strength." Turn the D40 always-on area reduction into an opt-in with
  two mechanisms.
- **Decision:** `cube_lattice(..., peak_correction: str = "none")` (CLI `--peak-correction`), three modes:
  - **"none" (DEFAULT)** — `PEAK_RATIO` unused: K0 matched, peak overshoots (~4.50 MPa), and the
    regularized descending branches still follow (the D39 state). Flips the D40 default OFF.
  - **"area"** (D40) — divide the calibrated strut area by `PEAK_RATIO`. Peak → 3.00 MPa, but area
    scales stiffness AND strength together so K0 drops ÷PEAK_RATIO (elastic branch no longer overlaps).
  - **"strength"** — reduce the concrete grade's `ft` by `PEAK_RATIO` for the STRUT materials
    (`dataclasses.replace(CONCRETE, ft=CONCRETE.ft/PEAK_RATIO)`), leaving area/`Gf` alone. Peak → 3.10
    MPa (≈ reference) with the AREA (hence K0) UNCHANGED — the elastic branch STILL overlaps AND the
    peak matches: two targets met by two independent knobs (area→K0, ft→strength). The D40 "not-taken"
    option, now available.
- **Why strength keeps K0 while area doesn't:** peak and K0 both scale ∝ area (one knob → can't hit two
  targets), but `ft` sets strength independently of the elastic stiffness (E, area), so it is the right
  lever when BOTH the initial slope and the peak must match. Trade-off of "strength": with `Gf` fixed the
  smaller `ft` widens the crack-opening at failure (`w=2Gf/ft`), so the softening tail runs slightly past
  the reference (zero at ≈1.2e-3 vs 1.1e-3) and the pre-peak is more rounded (diagonals yielding
  progressively) — minor/expected.
- **Verified (regularized + dynamic):** none → 4.50 MPa; area → 3.00 (K0 ÷1.5); strength → 3.10 (K0
  unchanged, elastic branch overlaps). Output stem gains `_<mode>` for area/strength.
- **Status:** accepted. Example-layer only (`build.py`, `pushover.py`); backend untouched; test suite
  unchanged (25 pass + the same one pre-existing known-fail).

### D42 — 2026-07-21 — Cube reset to a clean LINEAR-ELASTIC baseline (`examples/cube/elastic.py`): scalar-area lattice vs the 1D dispBeamColumn, E matched before any nonlinearity
- **Context (user):** "Results do not match. Let's start from scratch, and let's try to match linear
  elastic result first." Step back from the nonlinear peak/softening mismatch (D39–D41) and
  re-establish a VERIFIED elastic foundation before any cracking is layered on.
- **Diagnosis (linear-elastic, established with in-session evidence):**
  1. The elastic dispBeamColumn reference gives effective E = 30000 MPa exactly (= E·A/L).
  2. A lattice with PHYSICALLY-sized struts (area = tributary MESH·THK = 1000 mm²) is E = 47529 MPa
     = **1.584× too stiff** — every strut (orthogonal + diagonal) contributes axial stiffness in
     parallel, so the sum overshoots E·A/L.
  3. **Root cause = the diagonal struts, not the boundary condition:** a frictionless (roller) base
     gives 1.540 vs the fixed base's 1.584 (~3% — BC barely matters); an orthogonal-only lattice
     (horizon 1.2) is a SINGULAR mechanism (the D4 caveat) that won't solve; adding diagonals
     (horizon 2.1) drives it to 3.636×. The diagonals are structurally REQUIRED and are the source
     of the over-stiffness.
  4. This 1.58× elastic over-stiffness is the SAME phenomenon as the ~1.50× nonlinear peak
     overstrength (`PEAK_RATIO`, D39–D41). A single uniform strut area can cancel it at ONE point
     only (stiffness OR strength, not both) — which is exactly why the nonlinear study needed a
     separate `PEAK_RATIO` knob.
- **Decision (user, AskUserQuestion — "Beam-column E, one area"):** keep the simplest, most defensible
  baseline — reference = the 1D **elastic** dispBeamColumn (target = E only; a 1D coupon has no
  Poisson to match), lattice = a SINGLE uniform strut area, scalar-calibrated so the assembly's
  K0 = E·A/L. The rigorous alternatives offered and DECLINED for step 1: a 2D plane-stress-continuum
  reference with two orthogonal/diagonal strut-area families to also match Poisson ν.
- **What was built:** `examples/cube/elastic.py` — a clean, self-contained elastic-ONLY entry point
  (reuses the backend-agnostic `specimen.py`; carries NONE of the Concrete02 / regularization /
  dynamic-relaxation machinery). Elastic dispBeamColumn reference (`run_dispbeamcolumn_tension` +
  `concrete_uniaxial_elastic`), scalar `calibrate_area` (one tiny tension solve — K linear in area),
  elastic-strut lattice (`build_lattice`), stress–strain overlay via `viz.figure_pushover`. Output →
  `examples/output/cube/elastic/cube_tension_elastic.png`.
- **Result (N, mm; MESH 10 → 121 nodes, 420 struts):** calibrated strut area A = **631.20 mm²** →
  lattice effective E = 30000.0 MPa, secant-E spread over the whole pull = 0.000 MPa (perfectly
  linear), **E_lattice / E_ref = 1.00000 (0.000 % error)**. The two stress–strain lines overlap
  exactly on σ = E·ε.
- **Status:** accepted. ADDITIVE (new `examples/cube/elastic.py`); the D37–D41 nonlinear study
  (`pushover.py`, `build.py`) is LEFT IN PLACE (not deleted — "from scratch" read as "start the
  rebuild from a verified elastic baseline," reversible if a full prune is later wanted). Backend
  untouched; test suite unaffected (additive file, not imported by tests). Next: layer nonlinearity
  back on top of this baseline step by step (Concrete02 elastic branch first, then cracking/softening).

### D43 — 2026-07-21 — Cube elastic baseline extended with a `--material concrete02` mode: Concrete02's ELASTIC BRANCH reproduces E, with a slight sub-cracking curvature from transverse-strut compression
- **Context (user):** "Can you use Concrete02 for that elastic example?" — the D42 next step: swap the
  linear `Elastic` material for the real `Concrete02` law but stay in its elastic branch.
- **Decision:** `elastic.py main(material=...)` / CLI `--material {elastic,concrete02}` (default
  `elastic` = the D42 pure-elastic baseline, unchanged). `concrete02` builds Concrete02 struts
  (`build_lattice_rc`, one zone, `Truss`) and a Concrete02 dispBeamColumn reference, and pulls only to
  **0.98·eps_cr** (eps_cr = ft/E = 1e-4 here) so BOTH models stay on the elastic branch (Concrete02
  tension is linear to ft at eps_cr, then softens). Scalar `calibrate_area` is unchanged — it uses the
  elastic builder, and Concrete02's INITIAL tangent is the same E, so the same area (631.20 mm²)
  matches K0 in either mode. Effective E is reported from the INITIAL tangent (step 1), robust to any
  late-branch curvature. Output stem gains `_concrete02`.
- **Result (concrete02, N–mm):** initial-tangent **E_lattice/E_ref = 0.99992 (0.008 % error)** — the
  elastic modulus still matches. But the lattice curve is NOT perfectly straight: secant-E spread over
  the pull ≈ 1243 MPa (~4 %), and it bends slightly BELOW the reference near eps_cr. **Cause
  (constitutive, not a modelling error — the Elastic mode is dead straight, 0.000 spread):** Concrete02's
  COMPRESSION branch is parabolic (nonlinear from the origin), and the transverse/diagonal struts go
  into compression under Poisson contraction, so the assembly softens a touch as the pull grows. The 1D
  dispBeamColumn reference (single axial fiber, no transverse struts, linear tension to ft) cannot show
  this, so it stays straight → the small lattice-vs-reference divergence near ft. The vertical struts
  themselves stay linear (they crack only AT eps_cr); the diagonals see ~0.35·eps (far from their own
  cracking), so the curvature is transverse-compression, not tension.
- **Status:** accepted. Example-layer only (`elastic.py`); D42 default behaviour unchanged; backend and
  tests untouched. The elastic-modulus match (the step goal) holds; the mild sub-cracking curvature is
  a real Concrete02 feature and previews the next step (tracing THROUGH eps_cr into cracking/softening).

### D44 — 2026-07-21 — Cube PLASTIC pushover (`examples/cube/plastic.py`) on the calibrated area; `run_pushover` gains an `algorithm` knob, default primary = ModifiedNewton -initial for the softening pass
- **Context (user):** "With the adjusted lattice area after calibration, perform plastic pushover
  analysis. For [the] numerical part, use modified newton initial by default." The D42/D43 next step:
  trace the SAME calibrated-area Concrete02 lattice PAST cracking, statically.
- **Backend change (`run_pushover`, backward-compatible):** new `algorithm: tuple = ("Newton",)` — the
  PRIMARY solution algorithm (`ops.algorithm(*algorithm)`), used for normal stepping, as the FIRST rung
  of the failed-step retry ladder (then KrylovNewton → NewtonLineSearch), and on restore. Default
  `("Newton",)` reproduces the prior behaviour EXACTLY (verified: 25 pass + the same one pre-existing
  known-fail `test_rc.py::test_nonlinear_pushover_runs_and_yields`, D34 — unchanged). Callers pass
  `("ModifiedNewton", "-initial")` for a softening pass (iterate on the constant INITIAL elastic
  stiffness, never re-forming the singular cracked tangent). Replaced the old hard-coded
  `ops.algorithm("Newton")` primary/restore and the Newton-first retry tuple.
- **New script `examples/cube/plastic.py`** (reuses `elastic.py`'s `calibrate_area` + `cube_lattice` +
  `beamcolumn_reference` — DRY): plain (unregularized) Concrete02 struts at the elastic K0-calibrated
  area **631.20 mm² (NOT re-tuned for the nonlinear range)**, DisplacementControl tension pull to
  TARGET, **ModifiedNewton -initial by default**; `--full-newton` for contrast. Reference = the
  Concrete02 dispBeamColumn (single element → constant strain → traces the whole backbone to zero, D37).
- **KEY RESULT — the algorithm is NOT the limiter; the crack is a mechanism (confirms D38):** BOTH
  ModifiedNewton -initial AND full Newton trace the entire elastic branch (overlapping the reference)
  and die at the EXACT SAME point — the first crack, `eps = eps_cr = 1.0e-4`, peak **2.88 MPa**
  (≈ ft, ~4 % under from the D43 transverse-compression curvature), `u = 0.010 mm`, `conv=False`. The
  descending branch is not traced. Identical failure point for both algorithms ⇒ the plain-Concrete02
  tension crack is a GLOBAL singular tangent (a kinematic mechanism / brittle snapback), not an
  algorithmic convergence difficulty — so NO static path-follower crosses it, not even initial-stiffness
  iteration (K0 is non-singular, but it cannot track the load-factor snap as the band softens to zero;
  the load factor blows up to ~±30000 at the crack). This is exactly the D38 conclusion, now reproduced
  with the robust initial-stiffness algorithm the user asked for.
- **Decision:** keep **ModifiedNewton -initial as the default** for the plastic pushover (the user's
  request; and it IS the correct robust default for softening in general — it simply cannot beat a
  genuine mechanism in THIS problem). Report the static death at the crack honestly rather than forcing
  a branch. Tracing the full descending branch still needs (deferred, unchanged from D38/D39): crack-band
  regularization (gentler softening) + dynamic relaxation (`run_pushover_dynamic`) — the only combination
  that follows it. Plain Concrete02 + static is a mechanism at the crack, full stop.
- **Status:** accepted. Backend change backward-compatible (verified); new example file `plastic.py`
  (reuses `elastic.py`). Output → `examples/output/cube/plastic/cube_tension_plastic[_newton].png`.

### D45 — 2026-07-21 — Cube plastic pushover: crack-band regularization made optional (`--regularize`/`--gf`); it defeats localization (traces past first-crack) but not the mechanism (static still can't reach zero)
- **Context (user):** "Implement crack band regularization but explain me how it works first, and make
  it optional via cli args." Add the D20 length-regularization to the clean `plastic.py` as a toggle.
- **How it works (explained to the user):** a softening crack LOCALIZES into one band of width `h` (=
  strut length). The dissipated energy per crack AREA is `g_f * h` with `g_f = ft^2/(2*Ets)` (area of the
  linear softening triangle). Keeping the σ–ε curve fixed makes that energy scale with the mesh (smaller
  `h` → more brittle). Crack-band regularization holds the material fracture energy `Gf` constant by
  sizing the softening slope from `h`: **`Ets = ft^2*h/(2*Gf)`** ⇒ crack opening `w = 2*Gf/ft`,
  band-width-independent. Cube numbers: reference (h=100) `Ets`=3000, `Gf`=0.15; PLAIN strut (h=10) keeps
  `Ets`=3000 so `Gf`=0.015 (10× too little → drops ~10× too steeply); REGULARIZED strut (h=10) `Ets`=300
  (10× gentler) → `Gf`=0.15, matched. The existing `materials.concrete_uniaxial_regularized` (D20) is
  the exact implementation.
- **Implementation (example-layer only):** `plastic.cube_lattice(area, *, regularize, Gf)` — `material_for`
  returns `concrete_uniaxial_regularized(CONCRETE, 0, length, Gf, Gfc=250*Gf, residual_ratio=0.0)` when
  on, else plain `concrete_uniaxial_nonlinear`. CLI `--regularize` (default OFF = the D44 plain baseline)
  and `--gf` (default `REF_GF = ft^2*L/(2*0.1*E) = 0.150 N/mm`, the reference's implied Gf). Output stem
  gains `_reg`. No backend change (the D44 `algorithm` knob already there).
- **Result (ModifiedNewton -initial, N–mm):** PLAIN dies AT first crack — 2.88 MPa, `eps = 1.0*eps_cr`,
  conv=False (D44). REGULARIZED traces PAST first crack — up to 3.46 MPa at `eps = 1.6*eps_cr`, conv=False
  — but dies while still ASCENDING the post-cracking OVERSTRENGTH climb (verticals cracked, diagonals
  redistributing; the D39 lattice peak is ~4.5 MPa), before its own peak or ANY descending branch. So
  regularization defeats LOCALIZATION (static now survives first-crack and enters redistribution) but NOT
  the MECHANISM (as more struts crack the band still goes singular). The full aligned descending branch
  still needs regularization + dynamic relaxation (`run_pushover_dynamic`), exactly D38/D39 — deferred,
  the clear next step. The claimed branch ALIGNMENT is a regularized+dynamic property (D39); the
  regularized+STATIC run only demonstrates the localization fix + the overstrength onset.
- **Status:** accepted. Example-layer only (`plastic.py`); backend/tests untouched this turn (D44 already
  verified 25 pass + the one known-fail). Output → `.../cube_tension_plastic_reg[_newton].png`.

### D46 — 2026-07-21 — Cube plastic pushover: `--solver dynamic` (dynamic relaxation) added; regularized + dynamic traces the FULL descending branch and the lattice FOLLOWS the reference (staged rebuild of D39)
- **Context (user):** "Go" — take the offered next step: add the dynamic-relaxation solver so the
  (regularized) lattice can trace the full descending branch through the crack mechanism.
- **Implementation (example-layer only; no backend change — `run_pushover_dynamic` exists since D22):**
  `plastic.main(solver=...)` / CLI `--solver {static,dynamic}` (default `static`), dispatched by a new
  `_lattice_pushover` helper. `dynamic` calls `run_pushover_dynamic` driving the whole TOP EDGE (imposed
  ramp) with the cube-tuned params `periods_to_target=12, steps_per_period=30, damping_ratio=0.8`
  (mirrors the D38 cube run). `--full-newton` is now documented as static-only. Output stem gains
  `_dynamic`. Static path verified unchanged after the refactor (plain 2.88 MPa/1.0·eps_cr, regularized
  3.46 MPa/1.6·eps_cr — identical to D44/D45).
- **Result (dynamic, N–mm; both reach the full target u=0.12, eps=12·eps_cr, conv=True):**
  - **regularized + dynamic** — elastic branch overlaps the reference (K0), peak **4.50 MPa** (the diagonal
    overstrength, accepted per D39), then a CLEAN descending branch spanning the reference's strain range
    (both heading to ~0 near eps≈1.1–1.2e-3). The lattice now FOLLOWS the reference — the original goal
    ("I can't see the lattice follow the reference"). This cleanly reproduces the D39 default on the
    staged foundation (elastic baseline → plastic → regularization → dynamic).
  - **plain + dynamic** — full branch too, peak **3.51 MPa**, brittle snap + dynamic ringing (the D38
    plain-Concrete02 result); regularization is what makes the descending branch smooth and aligned.
  - Confirms the two-part fix: regularization (localization) + dynamic relaxation (mechanism) together are
    necessary and sufficient for the lattice to trace and follow the reference's softening branch.
- **Status:** accepted. Example-layer only (`plastic.py`); backend and tests untouched (`run_pushover_dynamic`
  reused as-is). Output → `.../cube_tension_plastic[_reg]_dynamic.png`. The cube staged rebuild (D42→D46)
  is complete: elastic match → Concrete02 elastic branch → plain static (mechanism) → regularization
  (localization fix) → dynamic (full aligned branch). Peak overstrength (~1.5×) left uncorrected per D39
  (the D40/D41 area/strength peak corrections remain available if wanted).
