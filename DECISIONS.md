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
