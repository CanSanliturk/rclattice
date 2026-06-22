# CLAUDE.md

Guidance for working in this repository.

> **Read [DECISIONS.md](DECISIONS.md) first.** It is the running log of key technical
> decisions and their rationale. Consult it before making design changes, and **append a new
> entry whenever a significant decision is made or reversed** (never rewrite history —
> supersede instead).
>
> **Next planned work:** implementing the RC frame pushover benchmark with lattices — full
> handoff spec in [TASK_RC_FRAME_PUSHOVER.md](TASK_RC_FRAME_PUSHOVER.md) (decision D18, pending).

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
      problem.py                  # backend-agnostic Problem (geometry, grades, BCs, loads)
      materials.py                # grade -> per-builder OpenSees material mapping (D15)
      mesh.py                     # gmsh grid (nodes + quads) + horizon strut connectivity
      builders.py                 # build_lattice / build_continuum -> FE Model + lumped mass (D12/D16)
      calibration.py              # lattice area calibration: static + modal periods (D16, scipy)
      viz.py                      # matplotlib visualizer (Agg): deformed lattice/continuum, modes (D17)
      model.py                    # generic FE Model (Node/Element/Uniaxial+NDMaterial/mass/...)
      opensees.py                 # ONLY module that imports openseespy: run_static + run_modal
    tests/                        # pytest (test_horizon, test_verification)
    examples/                     # cantilever_hello.py, cantilever_verify.py
```

Working today: the shared `Problem` + material mapping + `build_lattice`/`build_continuum`
(2D plane-stress quads) + a linear-static runner, with an elastic lattice-vs-continuum
verification. `model.py` is now a GENERIC FE assembly (Element with `etype`/`nodes`/`args`),
not truss-specific. Target structure to grow into incrementally (D12):

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
uv run python examples/cantilever_hello.py     # hello-world vertical slice
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

Early scaffolding, built incrementally. Working so far:
- Generic FE `model.py`, gmsh `mesh.py` (nodes + quads + horizon struts), linear-static
  `opensees.py`.
- Shared `Problem` (`problem.py`) + grade->material mapping (`materials.py`) + `builders.py`
  (`build_lattice`, `build_continuum` 2D plane-stress/strain).
- [cantilever_hello.py](src/examples/cantilever_hello.py): lattice cantilever (33 nodes, 92
  struts, tip ~ -5.2 mm at A=1e-3).
- [cantilever_verify.py](src/examples/cantilever_verify.py): elastic lattice-vs-continuum on
  one Problem; single-scalar strut-area calibration matches the continuum tip deflection.
- Modal calibration (D16): density-based lumped tributary mass, `run_modal` (ops.eigen), and
  `calibration.py` fitting orthogonal/diagonal strut areas (bounded, scipy least_squares) to
  the static deflection + first N periods. [cantilever_calibrate.py](src/examples/cantilever_calibrate.py).
- Frame + visualizer (D17): compound-rectangle geometry (`portal_frame`, joints merged) and
  box-based supports/loads (`BoxSupport`/`BoxLoad`); matplotlib `viz.py` renders deformed
  lattice/continuum side-by-side for static + mode shapes + an animated GIF.
  [frame/visualize.py](src/examples/frame/visualize.py) — frame periods agree ~1-3%
  lattice-vs-continuum.
- RC frame pushover, Stage 1 (D18/D19, ELASTIC pipeline): `opensees.run_gravity` (LoadControl)
  and `opensees.run_pushover` (gravity-constant → DisplacementControl, base shear = sum of base
  reactions, step-halving fallback); `builders.select_nodes` (post-build box node query);
  `viz.figure_pushover` (base-shear-vs-drift plot). [frame/pushover.py](src/examples/frame/pushover.py):
  benchmark frame in kip-in, lattice strut area calibrated to the continuum lateral stiffness
  (match ~271 kip/in, both converge to 15 in). 13 tests in [src/tests/](src/tests/), all passing.

- RC column studies (D19–D29, `examples/column/`): nonlinear RC lattice pushover + dynamic, with
  linear-material siblings, all calibrated to a reference. Materials done: uniaxial Concrete02
  (length-regularized struts, D20), Steel02 rebar, and nD **ASDConcrete3D + PlaneStress** for the
  continuum (D29). Reinforcement done: `Rebar` struts on shared nodes (D13). Column pushover compares
  the lattice to a **selectable reference** — `--reference {beamcolumn,continuum}` (D29 pushover, D30
  dynamic): the fiber `forceBeamColumn` (1D) or the 2D plane-stress continuum (`build_continuum_rc`,
  material-matched at the grade level). Pushover: lattice↔continuum agree ~1–2% on peak shear (both
  capture the 2D/diagonal action the 1D beam-column lacks). Dynamic (D30): the continuum's ASDConcrete3D
  is configured as **pure damage** (`plastic_frac=0`) — the closest match to the strut Concrete02's
  hysteresis (coupon-verified); seismic histories track closely, loop *shapes* differ slightly (the
  irreducible damage-vs-Concrete02 difference). Continuum dynamic is sine-only (heavy ~1.5 s/step).

Not yet: gmsh `embed` discrete bars in the continuum; 3D solids, shells, BeamColumnBuilder; a dedicated
results layer.
