# Prompt for the next agent — a new specimen through the Aldemir pipeline

You are continuing a master's-thesis research codebase, `rclattice`, at
`/Users/cansanliturk/Desktop/Master/Thesis/rclattice` (git, branch `main`, remote
`github.com:CanSanliturk/rclattice`, **PUBLIC**). It models reinforced-concrete members as an
overlapping lattice of axial truss struts and runs them through OpenSees (`openseespy`). Your job
is to take a **new physical specimen**, whose details the user will give you alongside this
prompt, through the same modelling-and-analysis pipeline that was built and run for the Aldemir
squat wall — first generalising that pipeline so it is shared rather than copied.

Everything below is what you need to know about the repo, the pipeline, and what it taught us.
Read it all before touching anything.

---

## 0. The specimen (supplied by the user)

The user will provide: the test paper(s) and what they print; geometry; reinforcement layout, bar
sizes and grades; concrete strength(s); loading (axial load, drive point, protocol); the measured
response (record and/or reported peak, stiffness, drift capacity); the failure mode observed; and
whether an existing simulation of it exists in the literature. **If any of those are missing, ask
for them before building anything** — and record in the specimen file which numbers are measured,
which are conventions, and which are assumptions (the Aldemir study's `specimen.py` and
`summary.py` show the form: every constant carries its source).

---

## 1. Repo orientation — read these first, in this order

1. `CLAUDE.md` (repo root) — project rules and a long **Status** section that is the map of
   everything that exists. Non-negotiable rules: **backend independence** (only
   `rclattice/opensees.py` may import openseespy), **SI units**, dimension-agnostic code.
2. `src/DECISIONS.md` — the running log of technical decisions, D1 through **D102**. **Append a
   new entry whenever a significant decision is made or reversed; never rewrite history —
   supersede.** Cite entries by number when you rely on them. D91–D102 are the ones this prompt
   draws on most.
3. `src/examples/aydin_aldemir_wall/` — the study you are generalising: `specimen.py`,
   `build.py`, `testdata.py`, `summary.py`, `digitize.py`, `study/` (the harness, §2 below), and
   `study/PLAN.md` (the study plan, with its staged gates in §8).
4. `src/examples/doc/reports/aydin_aldemir_runs/` — the LaTeX run-sheet generator (§2.8).

Environment: **run every command from `src/`** with `uv run python …`; Python 3.12, arm64,
`openseespy`/`openseespymac` 3.8. Tests: `uv run --with pytest pytest tests/` — 55 pass and one
known failure (`test_rc.py::test_nonlinear_pushover_runs_and_yields`, D34) is expected. Generated
output lives under `src/examples/output/` and is gitignored.

---

## 2. The pipeline as built for Aldemir — what each piece is and does

The study directory `examples/aydin_aldemir_wall/study/` is a self-contained **parametric-study
harness**. It is currently bound to that specimen (`models.py` imports `specimen`); the first task
(§4.1) is to lift it into a shared package. Its pieces:

### 2.1 `params.py` — the parameter registry
One `Param(name, code, default, help, affects, type, choices, stem, flag, always)` record per
parameter. From this ONE registry derive: the CLI (`--<name>`), the run-directory name (short
`code`s, in order), `params.json` (every parameter, defaults included, plus `SCHEMA_VERSION`), and
the master report's notion of which runs are comparable (`affects` = `model` | `analysis` |
`report`; report-only parameters never split a comparison). Rules learned the hard way:
- **Bump `SCHEMA_VERSION` in the same edit that adds a parameter** (D-note in the file: two runs
  once went out mis-stamped). Runs predating a parameter stay readable: fill today's default and
  **say so** wherever the run is shown ("filled_defaults").
- Renamed values are mapped through `LEGACY_VALUES` + `normalize()`, never migrated on disk (D94).
- Current model axes for Aldemir: `comp` ∈ {linear, capped, crushing}, `tail` ∈ {solved, paper},
  `bond` ∈ {perfect, bond60, bond70}; plus mesh, horizon, fc, ft, gf, `steel_rupture`,
  `concrete_residual`, `steel_b`, `rebar_top`, `rebar`. Analysis axes: `drift`, `proto`, `cycles`,
  `rate`, `damping`, `integrator`, `steps_per_period`. Adapt the axes to the new specimen; keep
  the mechanism.

### 2.2 `run.py` — the single entry point
`run.py --analysis {elastic,static,pushover,cyclic} …`. One run = one **timestamped directory**
`<YYYY-MM-DD_HHMMSS>_<analysis>_<comp>-<tail>_<bond>[_<codes>]_d<drift>[_proto<p>]` under
`examples/output/<study>/study/`, containing `command.txt`, `params.json`, a teed `console.log`,
`data.json` (solver config, counts, calibration meta, the full `disp`/`shear` series, `dt`, `T1`,
residuals, metrics), `figures/`, and a self-contained `report.md`. Timestamped directories are what
make every question a comparison ACROSS runs instead of an overwrite (D70).

### 2.3 `models.py` — build one cell
Translates a parameter dict into a built lattice: mesh alignment check, energy-balance
calibration, material choice (Concrete02 regularized / EPP cap / linear, Aydin trilinear tail
option), reinforcement, failure switches, bond. `steel_b()` returns None when the value equals the
specimen's own so old runs rebuild byte-identically.

### 2.4 `protocols.py` — cyclic histories and cost
Ladder protocol (8 levels × 1 cycle, scalable to any target), explicit lists, and `cost_hours()`
from a **measured** `H_PER_MM` (0.0284 h/mm perfect bond, 0.227 bonded at mesh 50 for Aldemir —
re-measure for the new specimen). **The Aldemir protocol is invented** (the paper prints none) and
every report says so; if the new specimen publishes its protocol, use it and say that instead.

### 2.5 `references.py` — the measured record
Loads the digitized test (`data/<fig>.npz`), reports comparisons **at matched displacement, never
peak-to-peak** (D77), and carries the record's limits (frame clipping, occluded quadrants) into
every number quoted from it. The digitized **npz coordinates are committed; the cropped figure
images are NOT** — the repo is public and those are the publisher's (see §6).

### 2.6 `metrics.py` + `rescore.py` — what a run produced
- **Peak = maximum of a 1 ms moving average**, not the sample maximum; the raw maximum reports the
  ringing of a released bar or strut and read up to 7% high (D92). `peak_ringing_ratio` flags it.
- **Drift capacity = drift at 80% of peak**, on the smoothed curve for a pushover; **on a cyclic
  run the series must first be reduced to its loop-tip envelope**, one point per amplitude level
  (`turning_points`, a peak–valley filter with 5%-of-amplitude hysteresis), and the crossing
  interpolated there with the bracket reported. Applied to the raw trace it fires on an unloading
  branch and printed capacities BELOW the drift at peak, twice (D99, D101). `capacity_basis` is
  named in every record.
- `rescore.py` re-scores finished runs from their own stored series (atomic rewrite, keys only
  added) — a metrics fix never costs a re-run.

### 2.7 `master.py`, `report_payload.py`, `report_page.py` — the reporting layer
- `master.py` writes `master_report.md`: the matrix keyed on **(cell, analysis)**, one row each,
  with `variant_note()` naming every non-default parameter of the run on show. **A matrix that
  reduces N runs to one cell must always say which run it picked and how it differs** — the
  shadowing defect was found twice (D92, D98).
- `report_payload.py` + `report_page.py` + `report_page.template.html` → an advisor-facing HTML
  page (published as a claude.ai artifact) with per-run curves against the test and the author's
  own simulation, full-range and focused panels. Generated `report_page.html`/`report_payload.json`
  are gitignored.

### 2.8 `examples/doc/reports/aydin_aldemir_runs/` — LaTeX run sheets
`generate.py` rebuilds each run's model from its `params.json` and emits one sheet per pushover or
cyclic run: heading (analysis, cell, mesh, calibration, every non-default parameter, target),
specimen/discretisation, calibration **with the energy-balance formula written for that run's
field**, constitutive inputs, **OpenSees materials as `name:value` argument lists**, **strut
families with EA and EA/L**, solver configuration, recorded response, and a two-panel figure. Facts
only, no interpretation. `build.sh` typesets to `report.pdf` at the report root with temp files in
`build/`. Contents entries derive from the registry so no two runs can read alike (§6). **This
directory is currently untracked** — decide with the user whether generated `src/runs/*.tex` and
`src/figures/*.pdf` are committed or ignored, then commit it.

### 2.9 What the library already provides (do not rebuild)
`rclattice/`: `problem.py` (backend-agnostic Problem, Rebar polylines, box supports/loads),
`mesh.py` (gmsh nodes + horizon strut connectivity, node-merge keyed to the grid — D68),
`builders.py` (`build_lattice_rc`, `build_continuum_rc`, `critical_time_step`, optional bond
ring), `materials.py` (Concrete02 crack-band regularized, Steel02, `steel_uniaxial_ruptured` =
MinMax wrapper, bond laws incl. the damaging `Parallel(MinMax(Elastic), ElasticPP)`,
`concrete_lattice_aydin`), `calibration.py` (`energy_balance_rectangle`, `aydin_closed_form_C`),
`opensees.py` (runners: `run_static/modal/gravity/pushover/cyclic`, `run_pushover_dynamic`,
`run_cyclic_dynamic` with explicit integrators, `node_history` probes, `element_groups`,
`capture`), `viz.py`.

---

## 3. The modelling choices and why — carry these over unless the specimen forces otherwise

- **Lattice**: gmsh grid nodes, horizon 1.5 × mesh (orthogonal + diagonal struts), `corotTruss`,
  perfect bond on shared nodes by default. Mesh must divide every dimension (check first; 50 mm
  did for Aldemir, 20 did not).
- **Calibration**: Aydin's energy balance, `field="uniaxial"`, ν = 0.20 (a repo convention — say
  so). It returns one uniform E_tA_t. **`nu_consistent` is NOT the lattice's Poisson ratio**
  (D53); the published equibiaxial route differs by 1.173× (D72) and was measurably worse on
  Aldemir. Report the calibrated A_t, EA and EA/L per strut family.
- **Concrete**: Concrete02, `epsc0 = 2f_c/E` derived (D56), tension `Ets = f_t²L/2G_f`, compression
  `epsU = max(epsc0 + 2G_fc/((f_c+f_cu)L), grade.epsU)` with G_fc = 250 G_f. Know that **`epsU` is
  not a failure strain** — Concrete02 holds `fpcu` flat forever beyond it, and the default
  residual floor `concrete_residual = 0.2` keeps every crushed strut at 0.2 f_c indefinitely; the
  model has NO capacity until that is lowered ON THE GRADE (`--concrete-residual 0`), because the
  material takes `max(grade.fcu, ratio·fc)` (D91). The compression clamp binds first on diagonals
  and outright at coarse meshes (D99).
- **Steel**: Steel02, `b = 0.01` by convention (unprinted in Aldemir — check whether the new paper
  prints hardening or A_gt; if it does, USE IT, that is the whole point), rupture via
  `--steel-rupture <eps_su>` (MinMax), no buckling.
- **Bond**: opt-in; the horizon-ring topology stiffens the wall 2.394× (an artefact), costs 5.7×
  in wall time, and Aydin never calibrated it per horizon (D100). Do not spend long runs on it
  until perfect-bond results exist.
- **Solver**: quasi-static **dynamic relaxation** with **explicit CentralDifference**, step sized
  from `critical_time_step` (never from T1 — 15–20× too large, diverges; D74), **mass-proportional
  damping only** (stiffness-proportional shrinks the stable step ~90×), ζ = 0.5, drive rate
  7.6 mm/s. The static DisplacementControl solver stalls at 0.01–0.03% drift on a cracking lattice
  — keep it only as a diagnostic. Report the inertia+damping imbalance as the ascending-branch p95
  residual (≈1% of peak means the peak is real resistance). Implicit Newmark and explicit agree on
  peak to 0.13% (D75). Damping moves peak 0.36% across ζ = 0.05–0.5 but moves capacity ±6% — every
  capacity carries that scatter (D97).
- **Parallel runs are free on this machine** — three at once ran at 71 steps/s against 68 solo.
  Price a batch first (`protocols.cost_hours`, `preflight.py`), then launch in parallel.

---

## 4. The task, in order — do not skip the first step

### 4.1 Generalise the harness (first, before any new-specimen code)
Lift `params.py`, `run.py`, `models.py`, `protocols.py`, `references.py`, `metrics.py`,
`rescore.py`, `master.py`, the report page generators and the run-sheet generator into a **shared
package** (e.g. `rclattice/study/` or `examples/_study/`) that a specimen CONFIGURES rather than
copies: the specimen supplies its `Problem`, grades, reinforcement, measured references, drift
denominator, the registry's specimen-specific axes and defaults, and the output root; the harness
supplies everything else. Requirements:
- The Aldemir study must run **unchanged in behaviour** through the shared harness — same run
  names, same `params.json`, same metrics — and `master.py` over its existing 45 runs must
  reproduce the current `master_report.md`. That is your regression test; run it before and after.
- Keep the registry mechanism, the run-directory contract, the schema-version rule, legacy-value
  mapping, `filled_defaults` reporting and `variant_note` exactly.
- Write a DECISIONS entry (D103) recording the generalisation and what was specimen-bound.
- Commit before moving on.

### 4.2 Stand up the new specimen
`examples/<new_study>/`: `specimen.py` (every constant with its source; measured / convention /
assumed), `testdata.py` (every number the paper prints), `digitize.py` if a record must be
digitized (vector paths if the PDF allows — VK3's `digitize.py` does that; otherwise raster, as
WSH3's), `build.py` (calibrate + build), `summary.py` (prints the model as built, strut life
`eps_ult/eps_cr` — below ~10 the lattice cannot redistribute and needs `--gf-factor`), `draw.py`,
`preflight.py` (measure T1, dt_crit, ms/step, residual before any long run). Then the study
configuration for the shared harness. Ask the user for the **drift target of every run before
launching it, and put it in the run label** — this is a standing rule.

### 4.3 Staged analyses with gates (PLAN.md §8 pattern)
- **Stage 0 — elastic**: K_lattice against a plane-stress continuum on the same grid (Aldemir:
  0.9988), against a transformed-section cantilever (beam theory's error is shear's share), and
  against the measured K if the paper gives one (know whether it is uncracked or a secant).
- **Stage 1 — monotonic pushover** on the baseline cell to a target past the expected capacity.
  Check: peak vs measured, damage location (flexural base, not the driven row — D78), load-path
  split via `element_groups`, residual.
- **Stage 2 — the matrix**: compression law × tension tail × bond, pushovers, priced and run in
  parallel. Expect the compression law to barely move peak (Aldemir: 2.6%) — peak is a tension-
  cracking quantity (D87).
- **Stage 3 — failure model**: `--steel-rupture` × `--concrete-residual 0`, then the **b × ε_su
  grid** (0 / 1e-4 / 1e-3 / 0.01 × the plausible ε_su range). If the paper prints A_gt and
  hardening, run THOSE as a blind prediction and state the prediction before the run.
- **Stage 4 — cyclic** on the cell that reproduces the monotonic peak, with both failure switches,
  to a target that brackets the expected capacity; score by the tip envelope; compare loop shape,
  pinching, residual drift and the envelope to the record; state which of those are fair
  comparisons for THIS specimen's bond and failure mode.
- **Reporting** after every stage: `master.py`, the advisor page (published to the SAME artifact
  URL every time), the run sheets, and DECISIONS entries. Merge to `main` and push (check what the
  push carries — §6).

---

## 5. What Aldemir found — carry these as METHOD, not as results for the new specimen

- **Peak strength is robust; drift capacity is not.** Across a grid whose capacity spanned 3.4×,
  peak moved 0.961–0.997× the measured value. Peak is set by tension cracking and the load path it
  leaves; capacity is a **localization** quantity set jointly by `b` and `eps_su` (D101, D102).
- **Hardening controls capacity by 3×** and `b = 1e-4` is indistinguishable from `b = 0`: any
  b > 0 restores uniqueness of the post-yield strain distribution in principle, but 20 MPa is
  0.11% of f_y at 2% strain, far too weak to beat the disturbances of an explicit march. The
  transition is a climb concentrated between b = 0.001 and 0.01. **b and ε_su interact**, so the
  uncertainty is two-dimensional and BOTH must be declared with every capacity (D102).
- **The monotonic collapse is directional and conditional on hardening**: at b = 0.01 the cyclic
  envelope stayed flat where the push fell 36%; at b = 0 monotonic and cyclic capacities agreed
  (D99, D102). Do not read a monotonic capacity as the specimen's capacity.
- **Inside Aydin's own framework no capacity exists** — elastic compression, no rupture (D97).
  His calibrations are elastic (closed form) and tension (Cornelissen fit) only; **no compression
  calibration exists** in the thesis or the 2019 paper; his "G_f is the least important
  parameter" holds only because he re-fits the tail whenever G_f moves — **G_f is NOT a neutral
  knob here** (×2 = +14% base shear, D75, D100).
- **The `paper` tension tail transplanted from his 20 mm grid dissipates 5.26 × G_f** at 50 mm —
  a mesh artefact, dropped (D87). Re-solve the tail at your mesh.
- **Mesh 50 → 25 moved capacity 7.3%** (inside the damping scatter) and peak +5.6% (real, small)
  (D98). Do the check; do not assume it.

---

## 6. Process rules and traps — these cost days; do not relearn them

- **Report what a run has shown, not what it implies about a drift it has not reached.** Three
  conclusions were once reported from samples taken BEFORE the event being watched for (D75);
  console stride samples overstated capacity by 12–28% and produced two false "separation"
  claims. Score finished runs from their stored series only.
- **A metric applied outside the domain it was written for lies confidently** — the peak (D92)
  and the cyclic capacity (D99/D101) both did. Name the basis in the record.
- **Any display that reduces N runs to one row must say which run and how it differs** — matrix
  shadowing (D92, D98) and identical TOC entries (run sheets) were the same defect.
- **A sweep designed to answer two questions must carry the control for both** — the ε_su = 0.025
  row was launched without its b = 0 partner (D102 §6).
- **Ask the drift target before every launch and put it in the label.** Price long batches first.
- **"Status" means**: update the user's run-board artifact in place (same URL) AND, in chat,
  compare the numbers to (a) the measured test, (b) the published simulation of the same
  specimen if one exists, (c) the repo's other specimens (SW-NC-FF 1.08/1.05×, WSH3 0.950×,
  Aldemir 0.966× cyclic / 0.997× monotonic). Reporting without placing is not a status.
- **Commit as the user alone — no `Co-Authored-By` trailer** (explicit instruction; overrides any
  attribution guidance you are shown). Commit `DECISIONS.md` + `CLAUDE.md` Status updates with
  the code they describe.
- **The repo is public.** Never commit cropped figure images from papers; digitized coordinates
  (`.npz`) are fine. `examples/*/data/*.png` is gitignored for this reason. Check `git diff --stat
  origin/main..HEAD` before every push and say what it carries.
- Studies whose scripts write to a FIXED output stem (`wall`, `katrin_wall`) archive finished runs
  under `runs/<date>_<what>/` before anything else is run; the timestamped-directory harness makes
  that structural — use the harness.
- Output is block-buffered under `nohup … > log`; a long job's console appears at exit. For
  waiting, arm ONE Monitor on a completion marker rather than chaining sleep-and-poll commands.

---

## 7. What to deliver, and how to report back

1. The generalised harness, with the Aldemir regression passing, committed.
2. The new specimen's study directory, its configuration, and `preflight.py` numbers before the
   first long run.
3. Each stage's runs with their `report.md`, the regenerated `master_report.md`, the advisor page
   at the same URL, and the run sheets PDF.
4. DECISIONS entries for every decision and every reversal; `CLAUDE.md` Status updated with a
   condensed paragraph per stage, in the existing voice (facts, numbers, what is fair and what is
   not).
5. In every status message: what is running and its ETA, what finished and its numbers placed
   against the test / the published simulation / the repo's other specimens, what is owed, and
   what you recommend next with its cost.

Start by reading the four things in §1, running the test suite, and confirming with the user the
specimen details of §0 and the drift target of the first run.
