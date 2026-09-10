# CLAUDE.md

Guidance for working in this repository.

> **Read [DECISIONS.md](DECISIONS.md) first.** It is the running log of key technical
> decisions and their rationale. Consult it before making design changes, and **append a new
> entry whenever a significant decision is made or reversed** (never rewrite history —
> supersede instead).
>
> **Project state:** the RC **column** and **portal-frame** studies are implemented — pushover +
> dynamic, nonlinear + linear, lattice vs fiber-`forceBeamColumn` / 2D-continuum references (see the
> **Status** section below and DECISIONS.md through D35). The original RC-frame pushover benchmark
> (D18) is implemented — see D18/D19 and the D34 frame rebuild.

## Project

`rclattice` — a Python library to model **reinforced concrete (RC) members and structures**
in 2D and 3D using the **lattice modelling technique**, and to run structural analysis
through **OpenSees** (via the `openseespy` package).

The workflow is:

1. Define geometry of members (beams, columns, slabs, walls) in 2D or 3D.
2. **Mesh** the geometry into a lattice: generate nodes and connect them with axial struts.
3. Add **reinforcement** (rebar) as elements tied to the lattice nodes.
4. Assign **materials**, **boundary conditions**, and **loads**.
5. **Translate** the internal model into an OpenSees model and **run the analysis**.
6. Parse OpenSees recorder output back into internal result objects for post-processing.

## Core design principles

- **Backend independence.** The internal domain model (geometry, mesh, members,
  reinforcement, materials, loads, BCs) must NOT import or depend on `openseespy`.
  OpenSees is a *backend*. All `ops.*` calls live only in the `rclattice/opensees/`
  layer. This keeps the model testable without OpenSees and allows mocking/swapping the
  analysis engine.
- **Dimension-agnostic.** Code supports 2D and 3D equally from the start. Nodes carry a
  coordinate of length `ndm` (2 or 3). Avoid hard-coding 2D or 3D assumptions; thread the
  dimensionality through explicitly (e.g. a `Model(ndm=..., ndf=...)`).
- **Build internally, simply, then emit.** Create nodes, struts, rebar, etc. as plain,
  well-typed internal objects first. Only at analysis time do we walk the model and emit
  OpenSees commands.
- **Define once, build many (D12).** A specimen is defined ONCE in a backend-agnostic
  `Problem` (geometry + reinforcement + physical material grades + BCs + loads). Independent
  *builders* translate that same Problem into different OpenSees idealizations so results are
  directly comparable. Lattice is the main aim; continuum/beam-column exist to verify and
  calibrate it.
- **Keep it simple.** Prefer small, composable, readable pieces over premature generality.

## Key modelling decisions (settled)

- **Lattice element type: axial truss struts only.** Lattice struts are axial members
  (OpenSees `truss` / `corotTruss` with a uniaxial material). Nodes therefore need only
  **translational DOFs** (`ndf = ndm`: 2 in 2D, 3 in 3D). No rotational DOFs by default.
- **Reinforcement coupling: shared, perfectly-bonded nodes.** Rebar elements share lattice
  nodes (or nodes snapped onto the rebar path). Perfect bond is the default. Design the
  reinforcement interface so bond-slip (zeroLength interface springs) can be added later
  without reworking the core.
- **Caveat — truss lattice stability.** A pure axial-truss lattice (especially in 3D) can
  form kinematic mechanisms if the lattice topology is not sufficiently braced/triangulated.
  When designing the mesher and when debugging singular-stiffness / non-convergence errors,
  consider lattice connectivity and restraint as a likely cause. (A Delaunay-edge lattice is
  naturally triangulated/tetrahedralized, which mitigates this.)
- **Node generation: gmsh.** Do NOT hand-roll meshing math. The mesher uses **gmsh** (Python
  API, native arm64 wheel) to place **nodes** for the member geometry — transfinite/structured
  for the regular pattern, Delaunay for irregular. gmsh is the single node source (D6/D10).
- **Strut connectivity: peridynamics-style horizon (D9).** Struts are NOT taken from mesh
  element edges. Instead, connect every node pair whose separation is `<= horizon * mesh_size`
  (default `horizon = 1.5`), deduplicated so only ONE element exists between any two nodes
  (reinforcement may be the exception later). On a grid this captures orthogonal (`s`) and
  diagonal (`s*sqrt2`) neighbours but not the next ring (`2s`), so the lattice is naturally
  triangulated and carries shear/bending without a separate bracing step.
- **Lattice patterns: regular grid AND irregular.** Support both, selectable per model:
  - *Regular grid* — structured node layout / transfinite mesh; predictable, easy to
    validate.
  - *Irregular* — gmsh's Delaunay meshing (optionally with controlled/min spacing); avoids
    mesh-induced directional bias in cracking/fracture studies.

## Modelling backends / builders (D12–D15)

One backend-agnostic `Problem` (geometry, reinforcement, material grades, BCs, loads) is
translated by three builders into OpenSees:

- **LatticeBuilder** (main aim) — gmsh nodes + horizon struts; uniaxial concrete struts and
  rebar struts on shared nodes.
- **ContinuumBuilder** (verification reference) — 3D solids (`stdBrick`/`SSPbrick`/tets); 2D
  **both** plane-stress quads (`quad`/`SSPquad`/`tri31`, ndf=2, like-for-like with a planar
  lattice) and shells (`ShellMITC4`/`ShellNLDKGQ` + `LayeredShell`, ndf=6, for thin
  walls/slabs) — selectable (D14).
- **BeamColumnBuilder** (single members) — `forceBeamColumn` + fiber section.

Cross-cutting rules:

- **Reinforcement (D13):** defined once as free 3D curves (polyline + diameter + steel grade).
  Builders consume the same definition — fibers (beam-column), rebar struts on shared nodes
  (lattice), and **discrete embedded bars** in continuum via gmsh `embed` (mesh conforms to
  the rebar path; rebar elements share solid/quad nodes → perfect bond, D5).
- **Materials (D15):** physical grades (concrete `fc, ft, E, Gf…`; steel `fy, E, b…`) map to
  per-builder OpenSees materials. First bundle: **Concrete02** (uniaxial struts/fibers),
  **ASDConcrete3D** (nD solids/shells), **Steel02** (rebar). This mapping layer is also where
  lattice **calibration** lives (fracture-energy regularization of strut softening by strut
  length/area so the lattice assembly matches the continuum). All listed materials/elements
  are confirmed compiled into this `openseespymac` build.

## Repository layout

`src/` is the project root (uv/build run from there); CLAUDE.md and DECISIONS.md sit one
level above it (D11).

```
<repo>/
  CLAUDE.md, DECISIONS.md         # docs stay at repo root
  src/                            # PROJECT ROOT (run uv from here)
    pyproject.toml, uv.lock, requirements.txt, .venv/
    rclattice/                    # the package (flat layout, hatchling, editable install)
      __init__.py
      problem.py                  # backend-agnostic Problem (geometry, grades, Rebar, BCs, loads)
      materials.py                # grade -> OpenSees materials: Elastic / Concrete02(+regularized) / ASDConcrete3D / Steel02 (D15/D20/D29)
                                  #   + bond_elastic_brittle: Aydin's BOND law, elastic-brittle with a
                                  #   residual plateau a=0.7 (D72)
                                  #   + concrete_lattice_aydin: TENSION-ONLY ElasticMultiLinear, his
                                  #   Fig. 1 backbone, tail solved from Gf per strut (D60)
      reinforcement.py            # map a Rebar polyline onto the lattice node chain (D13)
      mesh.py                     # gmsh grid (nodes + quads) + horizon strut connectivity
                                  #   + perturb_nodes: Aydin's random grid perturbation (D60).
                                  #   Compound-rectangle node merge is keyed to the GRID, not to a
                                  #   fixed 9-decimal round (D68) — see the CAUTION in Status
      builders.py                 # build_lattice[_rc] / build_continuum[_rc] -> FE Model + lumped mass (D12/D16/D29);
                                  #   critical_time_step() sizes an EXPLICIT run: dt < 2/w_max set by
                                  #   the stiffest/lightest ELEMENT, not T1 (D74)
                                  #   build_lattice_rc takes grid=/pairs= overrides for a perturbed lattice (D60)
                                  #   and OPTIONAL bond_material= (D72): separate steel nodes + a ring of
                                  #   bond struts instead of shared nodes. Default is still perfect bond
      calibration.py              # TWO calibration routes: structural fit to static+modal (D16, scipy)
                                  #   and Aydin's OLM energy balance -> uniform strut area (D47, no FE solve).
                                  #   The energy balance itself now has TWO affine fields (D72):
                                  #   field="uniaxial" (thesis, default) / "equibiaxial" (published 2019
                                  #   Appendix); aydin_closed_form_C reproduces its C=0.621 / 0.102
      viz.py                      # matplotlib (Agg): deformed shapes, modes, pushover, time-history, --draw models, modal-calibration figure (D17/D32/D35)
      model.py                    # generic FE Model (Node/Element[+kind]/Uniaxial+NDMaterial/mass/...)
      opensees.py                 # ONLY module importing openseespy: static/modal/gravity/pushover/
                                  #   EXPLICIT integrators supported too (D74): is_explicit() gates a
                                  #   Diagonal system + Linear algorithm, MASS-PROPORTIONAL damping
                                  #   only (betaKinit would shrink the stable step ~90x) and no
                                  #   sub-step rescue. 45x cheaper/step, 2.7x faster end to end
                                  #   cyclic/dynamic + fiber beam-column refs (pushover/dynamic/modal,
                                  #   cantilever + frame; D35). Cyclic: cyclic_protocol + run_cyclic
                                  #   (D48) and run_cyclic_dynamic (dynamic relaxation, D49).
                                  #   node_history takes ONE dof or a SEQUENCE, stored dof-major (D68)
    tests/                        # pytest (horizon, verification, calibration, rc, pushover, frame, continuum_rc)
    examples/
      cube/                       # plain-concrete 0.1 m cube: uniaxial-tension coupon, lattice vs single dispBeamColumn (D37); clean staged rebuild — elastic baseline (elastic.py, D42/D43) then plastic pushover (plastic.py, ModifiedNewton -initial, D44)
      compression_cube/           # plain-concrete 200 mm cube in uniaxial COMPRESSION between
                                  #   smooth platens (rollers + one pin), Aydin-calibrated, traced by
                                  #   very slow dynamic relaxation (D53-D56). Deliberately minimal:
                                  #   ONE Concrete02 grade at fc = 15.1 MPa, standard Truss struts.
                                  #   specimen/build/pushover; --rigid-platen ties the top face to
                                  #   one vertical DOF; damage + load-path + force-deformation figures.
                                  #   PLUS aydin_approach.py (D61): the SAME cube modelled Aydin's
                                  #   way instead — tension-only struts, perturbed grid, corotTruss —
                                  #   so approach is the only variable; specimen/build/pushover are
                                  #   untouched by it
      column/                     # RC cantilever-column studies (pushover/dynamic, nonlinear/linear, single-element BC)
      frame/                      # portal-frame studies (column + thinner beam; pushover/dynamic +/- linear; SI modal visualize)
      wall/                       # "Kutay's wall" — flexure-controlled RC shear wall on its pedestal
                                  #   (Sahinkaya et al. 2025 SW-NC-FF), calibrated by Aydin's OLM energy
                                  #   balance (D47). specimen/build/draw + elastic.py (D47), pushover.py,
                                  #   cyclic.py (reversed-cyclic to 4% drift, D48/D49), digitize.py
                                  #   (test loops from paper Fig. 14b -> data/*.npz, D51), replot.py and
                                  #   compare.py (redraw from saved JSON — never re-run the analysis),
                                  #   gauge.py (base vertical-strain profile across the section, D63)
      katrin_wall/                # "Katrin's wall" — WSH3 (Dazio, Beyer & Bachmann 2009). TEST DATA
                                  #   (D65): digitize.py (Fig. 7 hysteresis, grayscale panels ->
                                  #   data/wsh*_fig7.npz + protocol-anchored backbone), testdata.py
                                  #   (every number the paper prints, incl. the loading protocol).
                                  #   FULL kutay-wall pipeline (D66): specimen/build/draw/summary +
                                  #   elastic, pushover, cyclic, gauge, replot, compare — plus
                                  #   preflight.py (measure T1/ms-per-step before a long run) and an
                                  #   LVDT CHAIN giving curvature over height (Fig. 11c) and the
                                  #   paper's own extrapolated base curvature (Fig. 15d)
      aydin_aldemir_wall/         # "Aldemir's wall" — the SQUATTEST wall in the repo (3000x2250x210,
                                  #   ratio 0.75), no axial load, designed to yield in SHEAR. Tested by
                                  #   Aldemir, Binici & Canbay (2017) ACI SJ 114(2); simulated by Aydin,
                                  #   Tuncay & Binici (2019) JSE 145(9):04019091, whose PDF sits here.
                                  #   The specimen that paper predicts WORST (+21% peak force), which is
                                  #   the point of it. specimen/testdata/build/draw/summary only so far.
                                  #   CAUTION: the paper gives TWO incompatible geometries; summary.py
                                  #   START AT replica/../HANDOFF.md — see D72-D75. THICKNESS 120->210
                                  #   (D73); replica/ holds a reconstruction of AYDIN'S OWN model,
                                  #   which is how the error was measured
      vk3_wall/                   # "VK3" — squat wall-type BRIDGE PIER (Bimschas 2010, IBK
                                  #   Bericht 326, ETH Zurich, Ch. 5). The first SHEAR-relevant
                                  #   specimen: Lv/lw = 2.20 but rho_sw = 0.08%, shear 20-22% of
                                  #   displacement, combined shear+axial failure at 1.59% drift
                                  #   (D68). Full pipeline (specimen/testdata/build/draw/summary +
                                  #   elastic/pushover/cyclic/gauge/replot/compare/preflight).
                                  #   digitize.py reads VECTOR paths, not pixels -> ordered loops
                                  #   + the measured drive history. --draw on every analysis script
      doc/                        # LaTeX reports: build.sh, shared/ preamble, reports/<study>/report.tex
                                  #   (column, kutay_wall, katrin_wall); PDFs -> doc/out/ (not committed)
```

Working today: the shared `Problem` + material mapping (elastic + nonlinear) + builders
(`build_lattice[_rc]`, `build_continuum[_rc]`, 2D plane-stress) + static/modal/gravity/pushover/
dynamic runners, with the elastic lattice-vs-continuum verification and the nonlinear RC
column/frame studies (see Status). `model.py` is a GENERIC FE assembly (Element with
`etype`/`nodes`/`args`/`kind`), not truss-specific. Target structure to grow into incrementally (D12):

```
rclattice/
  problem/        # backend-agnostic Problem: geometry, reinforcement curves, grades, BCs, loads
  materials/      # physical grades + grade->OpenSees mapping (per builder) + lattice calibration
  mesh/           # gmsh node/element generation + horizon connectivity + rebar embedding
  builders/       # lattice, continuum (solid/quad/shell), beamcolumn -> internal analysis model
  opensees/       # ONLY place importing openseespy: analysis model -> ops.* + runners
  results/        # parse recorders -> result objects; verification/comparison helpers
```
(Add directories as features land, not all at once.)

## Environment

- **Platform:** Apple Silicon (arm64), macOS 15. Run **native arm64** — no Rosetta / x86
  conda environment is needed for OpenSees 3.8.
- **Why native works:** `openseespy` 3.8 is a pure-python shim that depends on
  `openseespymac`, which ships a native `macosx_13_0_arm64` wheel (requires macOS 13+).
- **Do NOT use the `opensees` (xara) package** — its arm64 wheels only cover Python 3.10
  and 3.13, not 3.11/3.12. Stick with `openseespy`.
- **Package/env manager:** prefer **uv**; plain `venv` + `pip` also works. **Target Python
  3.12** (native arm64). `requires-python` is `>=3.12` (lower fails universal resolution —
  see D2).
- **Dependencies live in [src/pyproject.toml](src/pyproject.toml)** (`[project.dependencies]`),
  locked in `src/uv.lock`. `src/requirements.txt` is an auto-generated export for non-uv
  consumers — do not hand-edit it; regenerate (from `src/`) with
  `uv export --no-hashes --no-emit-project -o requirements.txt`. Add/remove deps via
  `uv add <pkg>` / `uv remove <pkg>` (keeps pyproject + lock in sync).
- **All `uv` / build / test commands run from `src/`** (the project root).

### Setup (uv) — run from `src/`

```bash
cd src
uv sync                        # creates .venv + editable-installs rclattice
source .venv/bin/activate
uv run python examples/column/pushover_linear.py   # vertical slice (linear RC column)
uv run --with pytest pytest tests/             # tests
```

### Setup (venv fallback) — run from `src/`

```bash
cd src
python3 -m venv .venv          # ensure the python is arm64
source .venv/bin/activate
pip install -e .               # editable install from src/pyproject.toml (pulls deps)
```

### Sanity check

```bash
python -c "import openseespy.opensees as ops; ops.wipe(); print('OpenSees OK')"
```

## Conventions

- OpenSees integer tags (nodes, elements, materials) are emitted by the `opensees.py`
  translation layer. Current convention: domain objects carry their own integer ids and the
  backend reuses them directly as OpenSees tags. Domain code must not call `ops.*`.
- Use SI units consistently throughout (document the unit system in code/docstrings);
  OpenSees is unit-agnostic, so consistency is the user's responsibility.
- Type hints throughout; prefer dataclasses for the internal domain objects.

## Status

Built incrementally; the elastic verification slice and the nonlinear RC studies work.
- Generic FE `model.py`, gmsh `mesh.py` (nodes + quads + horizon struts), backend `opensees.py`
  (the only `ops.*` module: static / modal / gravity / pushover / dynamic + fiber beam-column refs).
- Shared `Problem` (`problem.py`, incl. `Rebar`) + grade->material mapping (`materials.py`:
  Elastic, Concrete02 plain + length-regularized D20, ASDConcrete3D+PlaneStress D29/D30, Steel02)
  + builders (`build_lattice[_rc]`, `build_continuum[_rc]`, 2D plane-stress/strain).
- Modal calibration (D16): density-based lumped tributary mass, `run_modal` (ops.eigen), and
  `calibration.py` fitting orthogonal/diagonal strut areas (bounded, scipy least_squares) to
  the static deflection + first N periods.
- Calibration output figure (D35): every calibrating run saves a modal figure — the first N mode
  shapes of the SELECTED reference vs the lattice (reference follows `--reference`: continuum via
  `run_modal`, or the subdivided fiber `run_beamcolumn_modal` / `run_beamcolumn_frame_modal`) with a
  per-mode periods table (`T_ref`, `T_lattice`, `Δ vs reference`) underneath. `viz.figure_modal_calibration`;
  `examples/{column,frame}/build.py:modal_calibration_figure`. Equal-mass, self-mass basis (D16).
- Frame + visualizer (D17): compound-rectangle geometry (`portal_frame`, joints merged) and
  box-based supports/loads (`BoxSupport`/`BoxLoad`); matplotlib `viz.py` renders deformed
  lattice/continuum side-by-side for static + mode shapes + an animated GIF.
  [frame/visualize.py](src/examples/frame/visualize.py) — frame periods agree ~1-3%
  lattice-vs-continuum.
- Pushover machinery (D18/D19): `opensees.run_gravity` (LoadControl) and `opensees.run_pushover`
  (gravity-constant → DisplacementControl, base shear = base-reaction sum, step-halving fallback);
  `builders.select_nodes` (post-build box node query); `viz.figure_pushover`. Large-drift extras:
  `corotTruss` struts + a residual compression plateau + `run_pushover_dynamic` (dynamic relaxation,
  D22). Transient `run_dynamic` (UniformExcitation) uses **modal damping** (D33, replacing
  stiffness-proportional Rayleigh — fixes the D28 base-shear spike).

- RC column studies (D19–D33, `examples/column/`): nonlinear RC lattice pushover + dynamic, with
  linear-material siblings, all calibrated to a reference. Confined-core vs cover grades + stirrup
  ties (D23); strut connectivity tunable via `--horizon` (D31); opt-in `--draw` analysis-model
  figures with reinforcement styled by kind (D32); a `single_beamcolumn.py` single-element fiber-BC
  reference variant. Materials done: uniaxial Concrete02
  (length-regularized struts, D20), Steel02 rebar, and nD **ASDConcrete3D + PlaneStress** for the
  continuum (D29). Reinforcement done: `Rebar` struts on shared nodes (D13). Column pushover compares
  the lattice to a **selectable reference** — `--reference {beamcolumn,continuum}` (D29 pushover, D30
  dynamic): the fiber `forceBeamColumn` (1D) or the 2D plane-stress continuum (`build_continuum_rc`,
  material-matched at the grade level). Pushover: lattice↔continuum agree ~1–2% on peak shear (both
  capture the 2D/diagonal action the 1D beam-column lacks). Dynamic (D30): the continuum's ASDConcrete3D
  is configured as **pure damage** (`plastic_frac=0`) — the closest match to the strut Concrete02's
  hysteresis (coupon-verified); seismic histories track closely, loop *shapes* differ slightly (the
  irreducible damage-vs-Concrete02 difference). Continuum dynamic is sine-only (heavy ~1.5 s/step).

- RC frame (D34, `examples/frame/`): the portal frame rebuilt as the cantilever column + a thinner
  (18-in) beam sharing the same grades — a full column-package mirror (pushover/dynamic,
  nonlinear/linear, selectable `--reference`). Beam concrete defaults **elastic** (the thin nonlinear
  beam forms a local mechanism the static pushover can't trace), opt-in `--nonlinear-beam`;
  `visualize.py` keeps the separate SI modal lattice-vs-continuum study.

- RC shear wall (D47, `examples/wall/`): the Sahinkaya et al. (2025) SW-NC-FF nonconforming
  flexure-controlled wall (3000x1000x200) standing on its 1900x600 pedestal, three concrete casts as
  zones, two reinforcement curtains collapsed onto in-plane bar lines. **Calibrated by Aydin's (2017)
  overlapping-lattice energy balance** — a homogenization, not a structural fit: equate continuum and
  lattice stored energy under an affine strain field and solve for one uniform `EA`. No FE solve, no
  reference model, no optimiser. The calibrated AREA is E-independent (one value serves every zone)
  and mesh-objective; `build.report_calibration` prints the EA it gives each zone (D62), since a
  truss's stiffness is EA/L and the balance returns EA and divides by E only to make it transferable
  — 9.729e7 / 1.595e8 / 5.656e8 N for test-unit / upper / pedestal at mesh 50. Key correction to the cited method: with a single `EA_t`, `nu_consistent` (the
  nu at which the normal-strain and shear calibration routes agree) is **0.18** for horizon=1.5, not
  Aydin's stated 1/3 (which overshoots shear stiffness 19%); `EnergyBalanceResult` reports it.
  **CAUTION (D53): `nu_consistent` is NOT the lattice's Poisson ratio** — that is `nu_effective`
  = 0.41, and the lattice is cubic-symmetric rather than isotropic at any nu (`cubic_anisotropy`
  = 1.38). Earlier wording here and in `examples/wall/build.py` called the lattice "isotropic at
  nu = 0.18"; that justification is wrong, though the wall's results stand, since flexure is governed
  by the C11 behaviour the balance matches exactly. Elastic result matches a transformed-section cantilever to ~1%
  (clamped base); the pedestal accounts for a further 6.5%. NOTE the specimen's real mechanism is
  rocking from PLAIN-bar debonding, which a perfect-bond lattice cannot represent (see the scope
  caveat below).

- RC shear wall, NONLINEAR cyclic (D48–D52, `examples/wall/cyclic.py`): the same Aydin-calibrated
  strut areas driven through the test's own protocol (paper Fig. 9) to 4% drift. A static
  path-follower stalls near 0.3% drift on a cracking lattice, so the run uses **dynamic relaxation**
  (`run_cyclic_dynamic`, D49) — drive rate set DIRECTLY in mm/s, Rayleigh damping on `betaKinit`
  never `betaK`, and `ops.wipeAnalysis()` before the transient integrator. Full 4% run at the old
  7.6 mm/s: 463,534 steps, ~7.1 h, converged.
  **Drive settings are now MEASURED, not assumed (D62/D64):** `specimen.QUASI_STATIC_RATE = 7.6`
  mm/s (the published value) and `DAMPING_RATIO = 0.2`, shared by pushover and cyclic so the two
  stay comparable. `quasi_static=True` on both dynamic runners returns the inertia + damping the
  drive and base reactions fail to balance (`S_base + S_drive = -sum_free(M.a + C.v)`).
  **DAMPING is the knob that matters, not rate (D64):** at 7.6 mm/s the steady residual is 10.3% of
  peak shear at zeta = 0.8, 12.4% at 0.05 and **4.0% at 0.2** — the contamination is crack-release
  RINGING, which damping suppresses. Slowing the drive buys nothing once damping is right (4.0 /
  12.5 / 4.6% at rates 7.6 / 3.8 / 2.0 is scatter, not a trend); what it DOES shrink cleanly is the
  one-off START-UP transient and, in a cyclic run, the transient at each of the 49 reversals (22.0 /
  10.5 / 5.3 kN), so `--rate 2.0` is for a loop-SHAPE study — at ~27 h against ~7 h. The pushover
  reports start-up and steady separately; the cyclic run cannot (every reversal is another velocity
  step). Read the residual RELATIVELY: HHT(0.7) inflates it ~60x against Newmark without biasing the
  shear (elastic wall: residual 2x the shear while the static answer is reproduced to 0.24%).
  `--rate` previously defaulted to None and derived the speed from `--periods 48`, which scales with
  amplitude (~20-79 mm/s); it now defaults to the constant. CAUTION: these are single runs and the
  peak shear itself scatters ~±3% (219.1 / 231.9 / 221.7 kN over the three rates) with which struts
  crack in which step.
  **Base strain gauge (D63):** two aligned node rows 250 mm apart at the wall-pedestal face give
  `eps(x) = (uy_top - uy_bottom)/250` across the 1000 mm width — the same construction as the test's
  vertical LVDT line. New backend probe `node_history=(nodes, dof)` on all four path-following
  runners (stride PLUS every reversal, so loop tips are exact); `examples/wall/gauge.py` reduces +
  draws it, `viz.figure_strain_profile` plots in microstrain with a least-squares line whose slope
  is the curvature and whose zero crossing is the neutral axis (a direct test of plane sections).
  Static pushover and dynamic cyclic agree exactly at +0.10% drift (NA -188.0 vs -188.1 mm). Saved
  into each run's `_data.json`, redrawn by `replot.py --levels 0.3 1.0 2.0` without re-running. `digitize.py` extracts the
  measured LOOPS from paper Fig. 14b (a point cloud, not an ordered path — no per-cycle energy);
  `compare.py` overlays model vs test and derives the test backbone as a loop-tip envelope.
  **Results:** peak base shear model +227/−221 kN vs digitized test +210/−212 → **1.08x / 1.05x**;
  envelopes agree within 11% over 0.5–2.0% drift. Two systematic gaps, both predicted by the
  idealization: the model peaks at **0.36% drift vs the test's 1.00%** and over-strengths 1.32x at
  0.3% (perfect bond → bars strain more per unit drift → earlier yield, stiffer after cracking), and
  it shows **no post-peak degradation** (1.34x at 4% drift; `Steel02` has no bar buckling, and the
  test's degradation is bar-buckling + core crushing). NOTE these numbers were produced at zeta =
  0.8 (D62/D64); the RATE is unchanged at the published 7.6, so only the damping differs now — the
  0.3%-drift pushover peak moves 223.1 -> 219.1 kN, inside the ~±3% run-to-run scatter. **SCOPE: strength, stiffness and damage
  location are fair comparisons; drift capacity, pinching, self-centering and energy dissipation are
  NOT** — 74% of the test's drift at 2% is rocking on debonded plain bars, absent by construction.
  Written up in `examples/doc/reports/kutay_wall/report.tex`.

- Plain-concrete COMPRESSION cube (D53-D56, `examples/compression_cube/`): a 200 mm cube between
  SMOOTH platens (base on rollers + one ux pin; top face driven), Aydin-calibrated, Concrete02,
  traced by very slow dynamic relaxation at 0.5 mm/s (~52 s/run). **Deliberately minimal after D56**:
  ONE grade at fc = 15.1 MPa (the wall's test-unit concrete), standard small-displacement `Truss`
  struts, mesh 20 mm, horizon 1.5 (3.01 via `--horizon`). `epsc0` is DERIVED as `2*fc/E` because
  Concrete02's compressive tangent is `2*fc/epsc0` whatever the grade's `E` says — quoting a
  published epsc0 instead silently mismatches the calibration by 15%.
  Exists because Aydin explicitly excludes compression (his Sec. 2.2 p. 30). **Findings.**
  (1) Eq. 2.1 balances a CONFINED energy, so it pins `E/(1-nu^2)` (verified to 1.000) — an
  unconfined cube therefore reads **0.88*E**, mesh-independently, and `nu_consistent` (0.18) is NOT
  the lattice's Poisson ratio (`nu_effective` = 0.41; cubic-symmetric, `cubic_anisotropy` = 1.38).
  (2) Peak is **0.634*fc** at 0.998*epsc0 — a ratio unchanged when fc went 30 -> 15.1, so it is
  geometric. **`--peak-correction strength` (D71) closes that gap on purpose:** it multiplies the
  strut grade's fc by a factor DERIVED from the strut list (1 / the verticals' 0.6341 share of the
  gross face = 1.5771 at horizon 1.5 / mesh 20), leaving the calibrated area — and hence K0 —
  untouched, the compression twin of the tension cube's D41 knob. `epsc0` scales with fc because
  Concrete02's tangent is 2fc/epsc0 and holding it would multiply every strut's stiffness by the
  same factor (tangent preservation verified to 1.000000), so THE PRICE is that the peak STRAIN
  moves by 1.577 as well — one parameter cannot fix stress and hold strain. `ft` is untouched, so
  splitting still starts at the same absolute strain, which against the enlarged epsc0 is
  RELATIVELY earlier (12.7% -> 8.0%): the peak is corrected, the mechanism is not. Under correction
  `capacity_accounting`'s `verticals/continuum` reads 1.0000 as its own self-check. NOT valid at
  horizon 3.01, where the diagonals still hold 22% at peak. NOT YET RUN — the prediction to test is
  `peak_over_fc` = 1.00. **Why (D55):** the lattice has MORE material than the continuum (crossing struts 1.79x
  the face), and compatibility alone would give 1.05x the continuum; transverse cracking at 12.7% of
  epsc0 then sheds the inclined path (diagonals 28% -> 0%), leaving only the vertical struts' 0.634
  of the face. The section is short of a LOAD PATH, not of material.
  (3) `--rigid-platen` (equalDOF on the top face) changes **nothing** — peak, damage and load split
  are identical, because the default already prescribes the same ramp to every top node; a printed
  `top_face_uy_spread` measures this (0.000e+00 mm for the default vs 1.2e-05 mm with the tie —
  identical prescribed motions are exact, the MPC only holds to solver tolerance).
  (4) `--solver static` (D57) adds a DisplacementControl pushover on a uniform top-face traction as
  the cheap alternative. **Whether it works depends entirely on the upper boundary condition (D59).**
  With the FREE (traction) face it stalls while still ascending at 0.1520% strain and stays there at
  max_iter 100 / 1000 / 5000 / **20000** — a genuinely singular tangent, so its maximum is a lower
  bound, not a capacity. With `--rigid-platen --max-iter 1000` it traces the FULL curve including the
  descending branch in 0.3 s, agreeing with the dynamic run (peak 9.85 vs 9.57 MPa, degradation
  74.3% vs 73.5%). A traction boundary lets a softening column shorten unopposed (face warps by
  11.5% of the shortening); tying the face forces redistribution and keeps the tangent non-singular,
  leaving only slow convergence. So `--rigid-platen` decides whether a static pushover is possible
  at all.
  (5) `--algorithm "<spec>"` (D58) selects the solution algorithm on either solver (all 19 OpenSees
  algorithms are accepted by this build). It barely matters statically: 15 of 19 stall at the SAME
  0.1520% strain, which is the proof the barrier is a mechanism rather than a convergence failure.
  Dynamically all complete, and KrylovNewton is ~25% faster than the Newton default. Defaults are
  per-solver: Newton (dynamic), ModifiedNewton -initial (static). Note it takes ONE QUOTED STRING —
  argparse would claim `-initial` as a flag of its own.
  New backend: `run_pushover_dynamic` accepts NEGATIVE targets, and both it and `run_pushover` carry
  `capture` (nodal-displacement snapshots), `element_groups` (force decomposition) and `equal_dof`
  (MPCs), so the two solvers are interchangeable;
  `viz.strut_strains` + `viz.figure_damage` draw the crack/crush pattern on a symlog strain field.
  **Run outputs follow VK3 (D70):** `specimen.run_dir`/`find_run` give every run its OWN timestamped
  directory under `examples/output/compression_cube/runs/` — `<stamp>_pushover_<solver>_<free|rigid>`
  and `<stamp>_aydin_rmax<r>` — each with `command.txt`, `console.log`, `data.json` and generically
  named figures (stress-strain, force-deformation, load-path, damage). The point here is NOT
  protecting a long run (52 s) but that every question on this cube is a comparison ACROSS runs
  (platen / solver / algorithm / horizon / Rmax/d), which a fixed stem overwrites. `console.log` is
  teed by `run_dir` itself, no shell wrapper — and openseespy's own WARNINGs land in it.
  `aydin_approach.py` finds its Concrete02 baseline via `find_run("pushover")` and names the run it
  overlaid.

- BRIDGE PIER VK3 (D68, `examples/vk3_wall/`): Test Unit VK3 of Bimschas (2010), IBK Bericht 326,
  ETH Zurich (doi:10.3929/ethz-a-006237119) — a three-unit campaign on squat wall-type bridge piers (1500x350, L_v = 3300, L_v/l_w = 2.20, 42phi14, phi6@200, N = 1300 kN).
  **The first SHEAR-relevant specimen in the repo**: nearly WSH3's aspect ratio but 3x less transverse
  steel (rho_sw 0.08% vs 0.25%), shear 20-22% of top displacement vs ~12%, and a combined SHEAR +
  AXIAL-LOAD failure at 1.59% drift. The chapter's own reading of that failure is lattice-shaped —
  the base compression zone that SUPPORTED the diagonal strut was destroyed, a load-path collapse
  `element_groups` measures directly (cf. D55). Bond is fair (continuous deformed bars, no splice).
  Mesh 50 mm needs NO geometric rounding at all — every dimension is a whole number of cells, so the
  shear span is exact. **f_c = 34 MPa is the ONLY concrete number the chapter gives**: E_c (SIA 262),
  f_t (EC2 at f_ck = f_cm-8) and G_f (MC90) are conventions, flagged in `summary.py`, and there is no
  M_cr to pin f_t against the way WSH3's paper allowed (D67). Strut life 14.2 — budget `--gf-factor 2`.
  **k0 = 60 kN/mm is NOT an elastic target** (a secant to first yield, 0.24x uncracked).
  `digitize.py` reads Figs. 5.10-5.13 as VECTOR paths via `pdftocairo -svg`, not as a raster: exact
  coordinates, no text/curve separation problem, and the loops stay ORDERED (so energy dissipation and
  per-cycle degradation are recoverable, unlike the WSH3 cloud). Axes calibrated on TICKS not the
  frame (the frame is +-110mm x +-950kN against +-100/+-800 labels — 19% error if assumed), self-checked
  against the independent base-moment axis to 1.0000. It recovers what the chapter never prints: peak
  shear **+891/-876 kN at ~0.95% drift** (vs predicted F_n = 851), the measured ELASTIC amplitudes
  (0.85/2.10/5.08 mm, which k0 mispredicts 3x) giving an uncracked stiffness of ~174 kN/mm, and the
  13 small intermediate cycles. `cyclic.py` therefore defaults to `--protocol measured`.
  NEW INSTRUMENT: the test's four diagonal string pots, reproduced as `specimen.panel_shear` +
  `deformation_components` -> the shear/flexure/base-crack split of Fig. 5.19-right, validated exactly
  on synthetic fields first. SCOPE: strength, stiffness, the shear split, loop shape and damage
  location are fair; anything at/past failure is not — the hoops' 90-degree hooks could OPEN after
  spalling and the model's never can, so expect over-retention beyond ~1.27% drift. No FE run launched.

- **CAUTION (D68), affects `katrin_wall` results:** `mesh_compound_rectangles` merged coincident nodes
  on `round(coords, 9)`, which fails on gmsh's 4e-9 float noise along a shared edge; the survivors were
  then joined by ~4e-9 mm struts whose EA/L is ~10 orders above every real strut — near-rigid links
  that wreck conditioning. Counts before the fix: **katrin_wall 109, vk3_wall 50, wall 0**. The merge
  is now keyed to the grid and `connect_horizon` refuses pairs below `1e-4*mesh`. D66/D67's WSH3
  numbers were produced WITH those links, so D67's attribution of its 0.065%-drift divergence purely
  to strut brittleness may be incomplete. MEASURED since: an ELASTIC A/B gives K_fixed/K_prefix =
  0.9861, so the links added 1.4% stiffness. The NONLINEAR effect is unquantified and the elastic
  number is a weak proxy — at L ~ 4e-9 the regularized tension branch is a FLAT PLASTIC PLATEAU that
  never fails (Ets = f_t^2*L/2Gf ~ 0), so those were ~109 spurious plastic ties at the wall-foundation
  interface, not rigid links. The WSH3 1.02% run nevertheless COMPLETED on the pre-fix mesher
  (143,443 steps, 21.7 h, peak +431/-434 kN = 0.950 of test, base curvature 0.99 of test) and is
  archived under `examples/output/katrin_wall/runs/2026-08-25_cyclic_1p02_gf2/`.
  **`wall` and `katrin_wall` scripts write to a FIXED stem, so completed runs there are archived
  under `examples/output/<study>/runs/<date>_<what>/` BEFORE anything else is run — see
  runs/README.md. `vk3_wall` (D68) and `compression_cube` (D70) instead give every run its own
  timestamped directory, so the archiving is structural rather than a rule to remember.**

Tests: `src/tests/` (horizon, verification, calibration, energy_balance, bond, rc, pushover, frame,
continuum_rc, aydin_cube) — 42 pass and one known-failing (`test_rc.py::test_nonlinear_pushover_runs_and_yields`,
the thin-nonlinear-beam lattice instability, D34).

- Aydin's COMPRESSION cube (D60, `examples/aydin_cube/`): the 100 mm cube of Aydin, Binici & Tuncay
  (2021) — a DIFFERENT source from the 2017 thesis the rest of the repo calibrates against, and the
  one that supplies the compression method the thesis excludes. Struts have **no compressive
  strength**: linear elastic forever (+ an RSM knee at `epsc0/3` dropping to `0.4*Ec`), trilinear
  tension softening, and the cube's strength EMERGES as the stability loss of columns that indirect
  tension has split free. Requires all three of: uncapped compression, a **perturbed grid**
  (`Rmax/d`, the one calibration parameter), and **corotTruss** (the failure is geometric — D53
  measured 6*fc with no failure on small-displacement struts). Results at the paper's own fitted
  `Rmax/d`, 5 realizations, no refitting: **LR 0.924*fc** (CoV 6.7%, peak at 1.00*epsc0),
  **NR 1.031*fc** (CoV 10.7%) — both inside his 10% criterion — and **1.51*fc at `--rmax 0`**,
  reproducing the uniform-grid "vertical locking" he reports as ~2*fc. `ElasticMultiLinear` not
  `HystereticSM`: the latter's plastic unloading puts a cracked strut at -10/-44 MPa at zero strain
  and would destroy the mechanism (the cost is no damage memory — NOT for cyclic work).
  **Key relation to `compression_cube`:** there `fc` is an INPUT and the answer is 0.634*fc; here it
  is an OUTPUT. Perturbation only ever LOWERS strength, so it cannot be bolted onto Concrete02 to
  reach fc — the geometry calibration and the tension-only law are not separable. Full list in
  `examples/aydin_cube/DIFFERENCES.md`.

- OUR cube, Aydin's approach (D61, `examples/compression_cube/aydin_approach.py`): the controlled
  swap. `specimen.py` / `build.py` / `pushover.py` are UNTOUCHED — same 200 mm cube, same fc = 15.1
  grade, same smooth platens, same energy-balance area — while the six modelling choices of D60 are
  swapped together (tension-only struts, perturbed grid, corotTruss, realizations, fc as an output).
  `Rmax/d` calibrated for THIS cube by the Fig. 4 loop: 0.050 -> 0.062 -> 0.093 -> **0.075**; the
  production set of 5 there gives **0.943*fc** (-5.7%, CoV 11.3%), peak at 1.36*epsc0. Two results the swap proves rather than argues:
  (1) **the elastic calibration is shared** — both approaches hit the closed-form modulus target
  (1.008 and 0.996) and the Aydin run independently reproduces the D53 confined-modulus finding
  (0.879*E unconfined), so strength differences are differences of APPROACH;
  (2) **the load-path collapse is topological, not constitutive** — swapping the whole constitutive
  law still sheds the inclined struts (28%->0% becomes 0.4% at peak), testing what D55 argued.
  Uniform grid under this approach **locks at >=1.86*fc and never turns over**, against Concrete02's
  0.634*fc capacity on the same grid. New `first_divergence` guard: one realization recorded a
  139 MPa peak on struts with no compressive strength while reporting `converged = True`, and
  poisoned a mean — diverged records are truncated and discarded, and LOCKING is reported separately
  from divergence.

- RC shear wall #2, TEST DATA ONLY (D65, `examples/katrin_wall/`): WSH3 of Dazio, Beyer & Bachmann
  (2009) — a half-scale slender cantilever (2000x150, L_v = 4560, ratio 2.28) whose DEFORMED bars stay
  bonded, so perfect bond is fair here and drift capacity is a fair comparison (unlike SW-NC-FF, where
  74% of the drift was rocking on debonded plain bars). The model was traced number-by-number back to
  the paper and MATCHES; three things were added to `summary.py` rather than changed: **`f_t` is 17%
  high** (the paper's own M_cr = 527 kN.m pins f_ctm = 2.97 MPa — EC2 at the CHARACTERISTIC f_ck =
  f'_c - 8, not at the mean — against `specimen.FT` = 3.46, so the model's M_cr is 575, +9%; flagged,
  NOT changed), the second concrete cast and strain penetration added to "not modelled", and the
  `nu_consistent` "lattice Poisson ratio" mislabel fixed per D53. `digitize.py` extracts Fig. 7 from
  a GRAYSCALE 3x2 panel raster (no colour to separate text from data — erosion glyph-seeds, axis
  rules cut BEFORE the component pass so the ductility tick stubs detach, then a |V| <= 1.03*V_max
  physical clip): WSH3 comes out at **+453/-449 kN vs Table 5's 454** and **+-93.2 vs 92.4 mm**, inside
  one pixel. The reduced curve is a backbone anchored on the PROTOCOL (amplitudes are known multiples
  of delta_y = 15.4), not an envelope of the cloud. `testdata.py` holds every printed number with its
  source. NOTE for the coming cyclic run: the test ran at **1.2-3.6 mm/MINUTE**, ~130x slower than
  SW-NC-FF's 7.6 mm/s, so damping is the only usable knob (D64); and WSH3 **never met the paper's own
  20%-drop failure criterion**, so its 2.04% drift is a lower bound. Out of reach either way: bar
  buckling and fracture, which are what ended the test at 1.79%.

- WSH3 ANALYSIS PIPELINE (D66, `examples/katrin_wall/`): the SW-NC-FF scripts ported file-for-file,
  with the differences that the specimen forces. **Bond runs the other way**: WSH3's DEFORMED bars
  stay bonded (flexure + ~12% shear, Fig. 9b), so loop SHAPE, energy dissipation, residual drift and
  drift capacity are all FAIR here — the cyclic run is the deliverable, the pushover only a
  diagnostic. Out of reach: bar buckling from 1.70% drift and the corner-bar rupture at 1.79% that
  ended the test. Confinement gets NO separate grade on purpose — the hoops are already in-plane
  rebar struts that restrain the compressed boundary, so a Mander grade would double-count.
  `epsc0` derived as 2fc/Ec (D56). **The static solver stalls at 0.067% drift** (SW-NC-FF's managed
  0.3%) because this wall cracks at ~0.025%, so `cyclic.py` defaults to `--solver dynamic`.
  **The test ran at 1.2-3.6 mm/MINUTE** — any usable drive is ~130x that, so the licence is the
  measured residual; 7.6 mm/s kept (same ABSOLUTE speed as SW-NC-FF, and this wall's shear is 2x, so
  the same contamination reads half). **50 mm is the coarsest legal mesh** (gcd of the model
  dimensions). Elastic: `K_lattice/K_transformed` = 0.930, same as SW-NC-FF; shear share 13% vs the
  measured 12%. `preflight.py` MEASURED T1 = 19.5 ms, 276 ms/step, residual 16.4 kN: to 1.02% drift
  is ~143k steps / **11 h**, full 2.03% ~518k / **40 h** (lower bounds). Post-processing validated on
  a synthetic response first. No cyclic run committed yet.

- ALDEMIR'S SQUAT WALL + two library features from its source paper (D72, `examples/aydin_aldemir_wall/`):
  the 2019 JSE paper (Aydin, Tuncay & Binici, 145(9):04019091 — the MISSING MIDDLE of the trilogy, and
  the only one of the three that models *reinforced* concrete) read in full. It **independently confirms
  D53**: it states the horizon-1.5d lattice's Poisson ratio runs "between 0.26 and 0.42 depending on the
  rotation of the loading axis" and is "about 0.33" at 3.01d, and rotating the loading axis through our
  own `nu_effective` gives **0.266-0.408** and **0.314-0.353** — outside corroboration that
  `nu_consistent` ~ 0.18 is not a Poisson ratio. Two features follow, both OPT-IN, neither changing any
  existing study:
  (1) **a second energy-balance field** — `energy_balance_*(field="equibiaxial")` is the paper's own
  published route (equal stresses, continuum side `E e^2/(1-nu)`), against the thesis default
  (`"uniaxial"`, transverse strain restrained, `E e^2/(2(1-nu^2))`). `aydin_closed_form_C` reproduces
  its published **C = 0.621 / 0.102** to 0.05%. THE TWO ROUTES DISAGREE **1.173x** at nu = 0.20,
  horizon 1.5 — and every truss stiffness scales linearly with EA, so that is the factor a wall's
  global stiffness would move by. The repo's SW-NC-FF `EA = 9.729e7` survives only by a
  NEAR-CANCELLATION (0.973x the paper's closed form at ITS nu = 1/3): field and nu differences almost
  offset. `EnergyBalanceResult` now carries `field` and `area_equibiaxial` so both are always visible.
  (2) **optional bond elements** — `build_lattice_rc(bond_material=...)`: each `Rebar` gets its OWN
  steel nodes, tied to the concrete by a ring of bond struts within `bond_horizon*mesh` (8 at horizon
  1.5, the coincident node excluded), with `materials.bond_elastic_brittle` = his Fig. 1(c) law
  (elastic to `ft`, brittle drop to `a*ft` = 0.7, flat plateau, symmetric). The residual is
  STRUCTURAL — a steel node is held only by its ring — so `residual=0` is rejected. `bond_mass_share`
  moves mass from host to steel nodes conserving the total (a zero-mass DOF is singular for the
  dynamic runners). COST: 2,806 -> 5,424 nodes and 13,471 -> 34,325 elements on this wall
  (CORRECTED 2026-09-05 by a build; the earlier 5,454 / 13,501 / 34,595 are 30/30/270 higher and
  their PROVENANCE IS UNEXPLAINED — 'bars run to the top' was the obvious guess and was tested
  and refuted: that option costs +60 elements, since the top bar line sits TWO cells below the
  top face at mesh 50).
  **The specimen** (Aldemir, Binici & Canbay 2017, ACI SJ 114(2):395-406, via the 2019 paper — the
  2017 paper is NOT in the repo): 3000x2250x120, **aspect ratio 0.75**, NO axial load, Ø8@100 both
  ways with 3 bars/position (rho = 1.257%), f_c 28 / f_t 1.85 (TS 500) / E_c 24,870 (ACI) / f_y 360.
  Mesh 50 divides everything (the paper's d = 20 does not: 2250/20 = 112.5). Measured K = 1,038 kN/mm
  and F = 964 kN; their lattice gave **1.21x and 1.38x on force** at the two horizons — the worst of
  their six specimens, which is why this is the target.
  **GEOMETRY: 3000 wide x 2250 high x 210 THICK (D73, 2026-08-29).** The panel was user-confirmed
  2026-08-28; the THICKNESS was corrected 120 -> 210 a day later and everything before that date is
  43% too thin. 120 dimensions ONE panel of a precast "double wall"; the load-carrying section is the
  whole stack. Found by building `replica/` — Aydin's OWN model, his panel/mesh/calibration/materials
  — which reproduces his Table 2 counts EXACTLY (20,385 nodes, 80,684 struts) and then misses his
  K_sim = 943.16 kN/mm by 1.743; lattice stiffness is exactly linear in thickness, so that factor IS
  the thickness. Five checks agree: the replica at 210 returns 0.979 of his K_sim unfitted;
  210 = 50+100+60 of the Fig. 10(a) plan stack; the uncracked section finally sits 1.31x ABOVE the
  measured 1,038 kN/mm (at 120 it was 0.75x, impossible); and V_flex/F_exp = 1.00 (rho falls from an
  overstated 1.257% to the true 0.718%). Fig. 10(a) is to scale and gives the panel. FIG. 4(f) IS NOT
  A GEOMETRY AT ALL (D82, author 2026-09-05): it is a DETAIL VIEW of the bottom 1500 mm of his model,
  so its "1500" is a crop height — which also dissolves the "not to scale" reading, since a crop is
  what produces an aspect matching none of its labels. Inverting Table 2's counts gives a 150x134
  grid = 3000x2680, matching neither printed geometry — that inversion is what made the thickness
  measurable. **The inversion is unique only UP TO TRANSPOSITION**: 2680x3000 fits both integers
  equally, Fig. 4(f)'s horizontal 2680 points at it, and `replica/` assumes the other. OPEN — it
  would move every published replica ratio.
  **ELASTIC (mesh 50, t=120 numbers, ratios unaffected by thickness):** K_lattice/K_continuum =
  **1.0018**, i.e. the lattice reproduces a plane-stress continuum on the same grid to 0.18%; the
  gap against a cantilever formula is BEAM THEORY's error, since shear carries **58.2%** of the
  flexibility at aspect 0.75. Mesh-objective: A_t halves exactly 50 -> 25 with K moving <1%.
  The published `field="equibiaxial"` route is measurably WORSE here (0.79 vs 0.92 against the
  continuum), so `"uniaxial"` stays the default as a measurement.
  **PUSHOVER: static stalls at 0.026% drift; the quasi-static dynamic run at t=210 peaks at
  903.4 kN = 0.94 of the measured 964 kN** (t=120 gave 591 kN = 0.61) with the first genuine
  descending branch. Residual is now saved as a SERIES: ascending-branch p95 = **1.6% of peak**, so
  the peak is only lightly contaminated — the 49% headline is the collapse transient.
  **THE PREMATURE COLLAPSE IS REAL, NOT NUMERICAL, AND STILL OPEN AFTER TEN EXPERIMENTS (D74/D75).**
  Best model (mesh 50, t=210, BOND, explicit): **peak 1,026.6 kN = 1.065x the measured 963.6** —
  closer than Aydin's own 1.21x overshoot — with an ascending-branch residual of 1.51%, so the peak
  is real resistance. But it loses its load path at **0.237% drift against the test's ~1%**.
  ONLY TWO THINGS EVER MOVED IT: thickness 120->210 (peak x1.53, drift x1.9) and bond (x1.14, x1.21).
  The other eight — damping, Gf, ringing, topology, mesh, constitutive law, reinforcement presence
  and INTEGRATION SCHEME — changed timing only. Implicit Newmark and explicit CentralDifference agree
  on peak to 1.0013 and on collapse drift to 0.004%, which closed the "it's a Newton failure" route.
  What remains unpinned is all BOND PARAMETERS the paper never prints (he gives only a = 0.7).
  `replica/` without bond reaches 0.42 of his published peak; **bond is now WIRED into the replica
  (D76: `run.py --bond`, `preflight.py --bond`) and NOT YET RUN.** The links cost 20,385 nodes /
  88,324 elements -> 28,081 / 149,802, his Table 2 concrete counts unchanged. The area ratio 0.01 is
  the PARENT's, calibrated at mesh 50; at the replica's mesh 20 the failure slip is 0.149 mm, on the
  floor of the Model Code band, and the ratio is NOT calibrated here — `--elastic --bond` against
  `--elastic` is the check to run first.
  Notable traps: **Gf is not a neutral knob** (x2 = +14.2% base shear — a caveat that also lands on
  D67's WSH3 result); **the trilinear backbone made everything worse** (1,510 solver failures,
  `ElasticMultiLinear` is path-independent so struts flip branches — Aydin runs it EXPLICITLY);
  **`figure_damage`'s `eps_crush` wants a POSITIVE magnitude** (a `< --2.3e-03` legend is the tell).
  **PROCESS:** three conclusions were reported from samples taken BEFORE the event being watched
  for, on a collapse already measured at ~30 steps against 250-step sampling. Report what a run has
  shown, not what it implies about a drift it has not reached.

  **BOND FIXES 2026-09-05 (D83/D84/D85), all four before-any-further-run items of D80 closed.**
  Supports now propagate to coincident STEEL nodes (30 unanchored base bar nodes → 0), and
  `select_nodes(kind="any"|"concrete"|"steel")` splits the sets: drive = concrete only (the actuator
  loads concrete; driving the duplicate pins zero slip at the row under study), base reaction set =
  everything (the anchored steel carries part of the shear). `strut_groups` gives bond links their
  own group, so D72's "rebar 0.7% of the overturning moment" was a probe artefact. `run_cyclic_dynamic`
  now gates the sub-step rescue on `explicit` (one failed step used to convert the rest of the march
  to Newton silently), reports `tuple(integrator)` instead of a literal, and snapshots `disps_peak`
  against the right step. `--rebar-to-top` ported to the parent (D78/D84). **Every bonded result
  before 2026-09-05 carries the anchorage defect** — D72's "2.19x stiffer than transformed section"
  in particular needs remeasuring.

  **STUDY HARNESS + STAGE 0/1/2 (D83-D86, `examples/aydin_aldemir_wall/study/`).** The parametric
  study of PLAN.md is built: a parameter REGISTRY (`params.py`) from which the CLI, run-directory
  names, `params.json` and the master report all derive; `run.py` (elastic/static/pushover/cyclic),
  `models.py`, `protocols.py`, `references.py`, `report.py`, `master.py`, `carrier_probe.py`.
  Cyclic amplitudes are an ARGUMENT — `--proto ladder --drift 0.012` scales the eight-level shape to
  any target, or pass an explicit list. Explicit runs SIZE THEMSELVES from `critical_time_step`
  (728 steps/period here, independently reproducing the hand-tuned 750). **STAGE 1 PASSED:** with
  bars to the top the wall traces to the full 0.30% target, converged, peak **943.8 kN = 0.979x the
  measured 963.6** at 0.266% drift, ascending-branch residual 1.1%, and damage stays at the
  flexural base (worst 1% at y = 25-1325, top 150 mm holds 3.3%) instead of tearing at the driven
  row — D78's diagnosis confirmed on the parent. **STAGE 2 PASSED:** K_lattice/K_continuum =
  **0.9988**. Two library fixes fell out: `concrete_lattice_aydin(fc_cap=)` gives the EPP
  compression branch, and **`run_cyclic` never wiped its domain** (the only model-building runner
  that did not) — now wiped by default, with `wipe=False` handing the domain, the build AND the
  gravity stage to the caller so `run_gravity` -> `run_cyclic` works for the first time.
  **STAGES 3-4, NOBOND HALF (D87).** Static cells all stall while still ascending at 0.0098-0.0120%
  drift. Quasi-static: **the COMPRESSION LAW barely matters and the TENSION TAIL dominates.** Three
  different compression branches (Concrete02 softening / linear-forever / EPP cap) span 2.6% on peak
  base shear — 943.8 / 944.2 / 968.3 kN, 0.98-1.005 of the measured — because only 2-19 of 10,905
  struts ever pass compressive yield; the peak is set by TENSION cracking and the load path it
  leaves (cf. D55/D78). The `paper` tail adds +22.6% for BOTH laws but dissipates a measured
  **5.26 x Gf**, so it is a mesh artefact of transplanting his 20 mm a2/a3 fit, not a physical
  alternative. CAUTION: `eppcomp`/solved peaked at 0.2997% of a 0.3000% target, still ascending, so
  its 1.005 is a lower bound. **BOND WIRED** (`wall_lattice(bond_law="aydin")` = PLAN §4's full-area
  residual-plateau law; the parent's own `--bond` keeps D75's force-slip form). Its artefact,
  measured for the first time with anchorage correct: **K_bond/K_perfect = 2.394**, so bonded
  stiffness comparisons measure the ring, not bond. Bond also costs **5.7x**, not the planned 3.4x
  (dt_crit 10.0 -> 4.5 us on light steel nodes: 248,088 steps vs 111,017, on 2.57x the elements).

  **BOND HISTORY — the 2026-08-28 disable is SUPERSEDED.** It was broken and disabled that day
  (D72); re-enabled 2026-08-29 by user instruction once the force-slip law fixed it (D75), and
  extended to `replica/` 2026-08-30 by user instruction (D76). The library guard stands: passing
  `bond_material` still raises unless `i_accept_the_known_bond_defect=True`, because the horizon-ring
  topology cannot reduce to perfect bond in the limit — calibration BOUNDS that artefact (~1% of
  elastic stiffness), it does not remove it. A coincident-node `zeroLength` spring would, and is
  unbuilt.

  **FAILURE MODEL, BOND LAW AND THE LIMITS OF THE REPLICATION (D91-D97, 2026-09-07/09).** Stock,
  this model cannot fail; two switches change that — `--steel-rupture` (MinMax on Steel02) and
  `--concrete-residual 0.0` (D22's floor lifted, ON THE GRADE or it does nothing). With them the
  peak lands at **0.997-0.998x the measured 963.6 kN** across five cells, against **Aydin's own
  1.208x overshoot** on this specimen. Drift capacity is **eps_su-dominated below 0.05** and
  crushing-capped above it: 0.591% (eps 0.025) / 1.033% (0.05) / 1.155% (crushing alone) / 1.374%
  (fracture alone), each **+/-6%** (D97 damping scatter).
  **BUT CAPACITY IS OURS, NOT HIS (D97):** with Aydin's own elastic compression law the wall never
  falls to 80% of peak through 3% drift and the rupture switch **never engages** — identical to the
  no-rupture twin to 0.000%. Crushing is the only reason this model has a capacity at all, and he
  excludes it. Peak strength is unaffected (a tension-cracking result, D87).
  **METRICS (D92):** `peak_shear` was a sample maximum and reported ringing on failure runs;
  `study/metrics.py` adds a 1 ms smoothed peak, a plateau width and an 80%-drop `drift_capacity`,
  and `study/rescore.py` re-scores finished runs from their own series. Sparse console samples
  overstate capacity by 12-28% — never read capacity off a log.
  **BOND (D95/D96/D97):** the law was `ElasticMultiLinear` and therefore REVERSIBLE — a broken bond
  healed on every reversal, dissipating exactly 0.0000 J per cycle. Aydin's word is "brittle", so
  `materials.bond_elastic_brittle_damaging` = `Parallel(MinMax(Elastic), ElasticPP)` makes the drop
  irreversible with an IDENTICAL envelope (verified to 0.49% on 21,244 elements) and is exposed as
  `--bond-damage`. Fig. 1(c) draws origin-oriented unloading on the CONCRETE branch and nothing on
  the Bond line, so the cyclic bond rule is OURS (ours dissipates 1.96x an origin-oriented one).
  `a` has TWO values, both his — published **0.7** (now the default) and **0.6** told to this
  project 2026-09-05 — and they differ by **1.166x on wall strength**, so no bonded result crosses
  that axis. Damping 5% (his) vs 50% (ours) moves peak by 0.36%: immaterial.
  **NAMING (D94):** the `comp` axis is `linear`/`capped`/`crushing` and `bond` is
  `perfect`/`bond60`/`bond70` — `nobond` used to mean PERFECT bond, the opposite of how it read.
  Legacy values are accepted forever via `params.LEGACY_VALUES`; nothing was migrated.

  **MESH OBJECTIVITY (D98):** the whole capacity sweep was measured at mesh 50. Rebuilt at **mesh 25**
  (13,531 -> 48,662 elements, 2.3 h -> 10.7 h) the same cell gives capacity **1.1084% vs 1.0328%** —
  a **7.3% shift, the same size as D97's +/-6% damping scatter** — so the ceiling is a property of the
  MODEL, not the discretization, and the 0.591/1.033/1.155/1.374% ordering stands (but quote no
  capacity better than +/-7%). PEAK does move for real: **+5.6%** (961.0 -> 1,015.1 kN = 0.997 ->
  **1.053x** measured), mean 1.0515 at matched drift. Note Aydin ran at 20 mm, FINER than either, where
  our trend says stronger — so mesh does not explain his 1.208x. Crack-band regularization (D20)
  regularizes DISSIPATION, not collapse drift; the latter is merely insensitive.
  **CYCLIC vs MONOTONIC (D99):** the controlled twin of the 1.5% push, cyclic, 8-level ladder, 11.9 h.
  **THE PUSH DEGRADES AND THE CYCLING DOES NOT** — loop tips 924/930/926/915 kN from 0.45% to 1.5%
  against a monotonic twin that falls 951 -> 933 -> 741 -> 609, i.e. cyc/mono **1.25x at 1.125% and
  1.50x at 1.5%**. Peak **931.1 kN = 0.966x** measured (push 0.997x), capacity **> 1.5%, not
  reached**. So the monotonic collapse is DIRECTIONAL and D97's 0.591/1.033/1.155/1.374% sweep is a
  MONOTONIC sweep. TWO TRAPS: `drift_capacity` is invalid on a cyclic series (it fires on an
  unloading branch — it printed 0.6978%, BELOW the drift at peak; reduce to the tip envelope first),
  and `epsU` is NOT a failure strain — Concrete02 holds `fpcu` flat forever past it, `epsU` is
  9-13x `epsc0`, and a strut at 2x `epsc0` still holds 91.6% of fc. The compression clamp
  `max(regularized, grade.epsU)` binds first on DIAGONALS: 8.01 vs the 8.0 floor at the default
  residual (harmless at residual 0, binds outright at mesh 100). New `--steel-b` (schema v5):
  hardening was hardcoded at 0.01, unprinted like eps_su, and worth 1.27x f_y at eps_su = 0.05.
  **WHAT AYDIN ACTUALLY CALIBRATES (D100, read from the 2019 PDF):** exactly two things — elastic
  `Et*At` in closed form, and TENSION fitted to Cornelissen (1986) via his Fig. 2 loop, after which
  the specimens are BLIND predictions. **No compression calibration exists anywhere** ("Concrete in
  compression is assumed to be elastic"); the only one in the trilogy is the 2021 cube's GEOMETRIC
  `Rmax/d`. He also states he never calibrated bond per horizon, so stop hunting those numbers. And
  his "Gf is the least important parameter" does NOT contradict D75's +14.2%: he RE-FITS a1/a2/a3
  whenever Gf moves, so the comparison was never like-for-like — which also unsettles D67's WSH3
  `--gf-factor 2`.
  **REPORTING (D98, second instance of D92's defect):** the matrix keys on comp x tail x bond x
  analysis, so any OTHER parameter silently substitutes — mesh 25 advertised 1.053x in the cell whose
  headline is 0.997x. `master.variant_note()` now names what is non-default about the run shown
  (amber chip on the advisor page). A matrix that reduces N runs to one cell must say which it picked.

Not yet: the aydin_aldemir_wall replica bond run (staged: `preflight.py --bond --explicit`, then
`run.py --elastic --bond`, then `run.py --drift 0.0025 --bond --explicit`), and its
cyclic/gauge/replot/compare scripts; the WSH3 cyclic run itself (staged: `cyclic.py --drift 0.0102 --gf-factor 2`, and see the
D68 caution above); the VK3 runs (staged: `preflight.py`, then `cyclic.py --gf-factor 2 --draw`); gmsh `embed` discrete bars in the continuum; 3D solids, shells, a packaged BeamColumnBuilder
(fiber refs live as `run_beamcolumn_*` in `opensees.py`); a dedicated results layer.
