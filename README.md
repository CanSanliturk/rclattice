# rclattice

A Python library for modelling **reinforced concrete (RC) members and structures** in 2D using the **lattice modelling technique**, with structural analysis powered by [OpenSees](https://opensees.berkeley.edu/) via the `openseespy` package.

## Overview

The lattice approach discretises a concrete body into a network of axial struts connected at nodes. Nonlinear uniaxial material laws (tension softening, compression crushing) are assigned to each strut so that cracking and failure emerge naturally from the strut network rather than being prescribed by element formulations. Reinforcement bars are added as additional struts sharing nodes with the concrete lattice, giving a direct representation of bond and load transfer.

`rclattice` implements this workflow end-to-end:

1. **Define** the geometry (beams, columns, walls) and RC detailing (rebar, material grades).
2. **Mesh** with [gmsh](https://gmsh.info/) — regular structured grid or Delaunay irregular.
3. **Connect** struts via a peridynamics-style horizon rule (all node pairs within `horizon × mesh_size`).
4. **Build** the OpenSees model: either a lattice (main aim), a continuum quad mesh (verification), or a fibre beam-column (calibration reference).
5. **Run** static, pushover, gravity + lateral, or nonlinear dynamic analysis.
6. **Calibrate** strut cross-sections to match a reference continuum response (static deflection + modal periods).
7. **Visualise** deformed shapes, pushover curves, mode shapes, and animated GIFs.

## Features

- **Backend-agnostic domain model** — geometry, mesh, reinforcement, materials, BCs, and loads are defined without any `openseespy` dependency; OpenSees calls are isolated to `rclattice/opensees.py`.
- **Three analysis builders from one `Problem`** — lattice, continuum (plane-stress quads), and fibre beam-column are built from the same backend-agnostic specimen definition so results are directly comparable.
- **Lattice calibration** — `scipy` least-squares fitting of orthogonal and diagonal strut areas to match a continuum static deflection and the first *N* natural periods.
- **Nonlinear RC column** — confined/unconfined concrete zones (Concrete02 with fracture-energy regularisation), longitudinal + stirrup reinforcement (Steel02), axial pre-load, pushover and nonlinear time-history analysis.
- **Portal-frame pushover** — gravity + lateral displacement control, base-shear extraction, and a pushover curve plot.
- **Visualiser** — matplotlib-based renderer for deformed lattice/continuum meshes, mode shapes, and animated mode-shape GIFs.

## Installation

**Requires Python 3.12+ on macOS (arm64), Linux, or Windows.**

### With uv (recommended)

```bash
cd src
uv sync                        # creates .venv and editable-installs rclattice
source .venv/bin/activate
```

### With pip / venv

```bash
cd src
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Sanity check

```bash
python -c "import openseespy.opensees as ops; ops.wipe(); print('OpenSees OK')"
```

> **macOS (Apple Silicon):** `openseespy` 3.8 ships a native `macosx_13_0_arm64` wheel via the `openseespymac` package — no Rosetta or x86 environment needed.

## Quick Start

### RC Column — pushover

```bash
cd src
uv run python examples/column/pushover.py
```

Runs a nonlinear pushover of a 12 ft × 24 in RC cantilever column (Concrete02 core + cover, Steel02 rebar) under 180 kip axial load and monotonic lateral displacement. Plots base-shear vs drift.

### RC Column — nonlinear time-history (El Centro)

```bash
uv run python examples/column/dynamic.py
```

Applies the 1940 El Centro ground motion to the same column via `UniformExcitation`. Plots displacement–time history and base-shear hysteresis.

### Portal Frame — pushover

```bash
uv run python examples/frame/pushover.py
```

Builds a two-bay portal frame, calibrates the lattice to the continuum lateral stiffness, then runs a gravity + lateral pushover.

### Frame — visualiser

```bash
uv run python examples/frame/visualize.py
```

Renders deformed lattice and continuum meshes side-by-side and saves animated mode-shape GIFs.

## Project Structure

```
.
├── CLAUDE.md                   # Codebase guidance for Claude Code
├── DECISIONS.md                # Running log of key technical decisions
├── TASK_RC_FRAME_PUSHOVER.md   # Handoff spec for the RC frame pushover task
└── src/
    ├── pyproject.toml          # Project metadata and dependencies
    ├── requirements.txt        # Pinned dependencies (uv export)
    ├── uv.lock                 # uv lockfile
    ├── rclattice/              # The library
    │   ├── problem.py          # Backend-agnostic Problem (geometry, grades, BCs, loads)
    │   ├── materials.py        # Grade → OpenSees material mapping + calibration
    │   ├── mesh.py             # gmsh node/quad generation + horizon strut connectivity
    │   ├── builders.py         # build_lattice / build_continuum / build_lattice_rc
    │   ├── calibration.py      # Static + modal calibration (scipy least_squares)
    │   ├── model.py            # Generic FE model (Node, Element, Material, …)
    │   ├── opensees.py         # Only module that imports openseespy; analysis runners
    │   ├── reinforcement.py    # Rebar struts on shared lattice nodes
    │   └── viz.py              # matplotlib visualiser (deformed mesh, modes, GIFs)
    ├── examples/
    │   ├── column/             # RC cantilever column studies
    │   │   ├── specimen.py     # Shared geometry / grades / rebar definition
    │   │   ├── build.py        # Lattice + fibre beam-column builders + calibration
    │   │   ├── pushover.py     # Nonlinear pushover
    │   │   ├── pushover_linear.py
    │   │   ├── dynamic.py      # Nonlinear time-history (El Centro + sine)
    │   │   ├── dynamic_linear.py
    │   │   └── excitation.py   # Ground-motion / sine loading helpers
    │   ├── frame/              # Portal frame studies
    │   │   ├── specimen.py
    │   │   ├── pushover.py
    │   │   ├── pushover_rc.py
    │   │   └── visualize.py
    │   ├── data/
    │   │   └── elcentro.at2    # 1940 El Centro NS ground motion record
    │   └── doc/
    │       ├── column_report.tex
    │       └── column_report.pdf
    └── tests/                  # pytest suite (13 tests)
```

## Dependencies

| Package | Role |
|---|---|
| `openseespy` ≥ 3.8 | OpenSees finite element engine |
| `gmsh` ≥ 4.13 | Node placement and mesh generation |
| `numpy` ≥ 2.0 | Array operations |
| `scipy` ≥ 1.17 | Least-squares calibration, ODE integration |
| `matplotlib` ≥ 3.11 | Visualisation and plots |
| `pillow` ≥ 12.2 | Animated GIF export |

## Status

| Feature | Status |
|---|---|
| Generic FE model + gmsh mesher | ✅ |
| Horizon strut connectivity | ✅ |
| Lattice + continuum builders (2D) | ✅ |
| Static + modal calibration | ✅ |
| Portal frame pushover (elastic) | ✅ |
| Nonlinear RC column (Concrete02 + Steel02) | ✅ |
| Nonlinear dynamic (time-history) | ✅ |
| RC frame nonlinear pushover | 🔄 In progress |
| 3D solids / shells | ⬜ Planned |
| Beam-column builder | ⬜ Planned |

## Licence

This project was developed as part of course CE7016. No licence is currently applied — all rights reserved.
