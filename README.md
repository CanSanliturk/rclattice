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
6. **Calibrate** strut cross-sections to match a reference response (static deflection + modal periods, or initial lateral stiffness).
7. **Visualise** deformed shapes, pushover curves, time-histories, mode shapes, and animated GIFs.

## Features

- **Backend-agnostic domain model** — geometry, mesh, reinforcement, materials, BCs, and loads are defined without any `openseespy` dependency; OpenSees calls are isolated to `rclattice/opensees.py`.
- **Three analysis builders from one `Problem`** — a lattice (main aim), a 2D continuum (plane-stress quads, elastic or nonlinear ASDConcrete3D), and a fibre `forceBeamColumn` are built from the same backend-agnostic specimen so results are directly comparable; the column/frame studies pick the verification reference at the CLI (`--reference beamcolumn|continuum`).
- **Lattice calibration** — `scipy` least-squares fitting of orthogonal and diagonal strut areas to match a continuum static deflection and the first *N* natural periods, plus a scalar/secant fit to a target initial lateral stiffness.
- **Nonlinear RC column** — confined/unconfined concrete zones (Concrete02 with fracture-energy regularisation), longitudinal + stirrup reinforcement (Steel02), axial pre-load, pushover and nonlinear time-history analysis, each with a fully linear sibling.
- **RC portal frame** — built from the same RC column plus a thinner connecting beam; pushover and nonlinear time-history mirroring the column studies (gravity + lateral control, base-shear extraction, pushover/hysteresis plots). Beam concrete is elastic by default, with an opt-in nonlinear beam.
- **Visualiser** — matplotlib renderer for deformed lattice/continuum meshes, mode shapes (+ animated GIFs), pushover curves, seismic time-histories and hysteresis loops, plus opt-in `--draw` analysis-model figures with reinforcement styled by role (longitudinal vs stirrup).

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

All example scripts accept `--draw` (save the analysis-model figure) and write outputs under `examples/output/`.

### RC Column — pushover

```bash
cd src
uv run python examples/column/pushover.py
```

Runs a nonlinear pushover of a 12 ft × 24 in RC cantilever column (Concrete02 core + cover, Steel02 rebar) under 180 kip axial load and monotonic lateral displacement, plotting base-shear vs drift against a selectable reference (`--reference {beamcolumn,continuum}`).

### RC Column — nonlinear time-history

```bash
uv run python examples/column/dynamic.py
```

Runs a nonlinear seismic time-history of the same column — by default a resonant sine (`--excitation elcentro` for the 1940 El Centro NS record) via `UniformExcitation` with modal damping — plotting drift/base-shear histories and the hysteresis loop. `--reference {beamcolumn,continuum}` selects the comparison model.

### Portal Frame — pushover

```bash
uv run python examples/frame/pushover.py
```

Builds a one-bay, one-storey RC portal frame (the cantilever column plus a thinner 18 in beam sharing the same grades), calibrates the lattice, then runs a gravity + lateral pushover against the selected reference. The beam concrete is elastic by default (`--nonlinear-beam` makes it Concrete02).

### Frame — visualiser

```bash
uv run python examples/frame/visualize.py
```

Renders deformed lattice and continuum meshes side-by-side and saves animated mode-shape GIFs (a separate SI-unit modal study).

## Project Structure

```
.
├── CLAUDE.md                   # Codebase guidance for Claude Code
├── DECISIONS.md                # Running log of key technical decisions (D1–D34)
└── src/
    ├── pyproject.toml          # Project metadata and dependencies
    ├── requirements.txt        # Pinned dependencies (uv export)
    ├── uv.lock                 # uv lockfile
    ├── rclattice/              # The library
    │   ├── problem.py          # Backend-agnostic Problem (geometry, grades, Rebar, BCs, loads)
    │   ├── materials.py        # Grade → OpenSees materials (Concrete02 / ASDConcrete3D / Steel02)
    │   ├── mesh.py             # gmsh node/quad generation + horizon strut connectivity
    │   ├── builders.py         # build_lattice[_rc] / build_continuum[_rc] / select_nodes
    │   ├── calibration.py      # Static + modal calibration (scipy least_squares)
    │   ├── model.py            # Generic FE model (Node, Element, Material, …)
    │   ├── opensees.py         # Only module that imports openseespy; analysis runners
    │   ├── reinforcement.py    # Rebar struts on shared lattice nodes
    │   └── viz.py              # matplotlib visualiser (deformed mesh, modes, curves, GIFs)
    ├── examples/
    │   ├── column/             # RC cantilever column studies
    │   │   ├── specimen.py     # Shared geometry / grades / rebar definition
    │   │   ├── build.py        # Lattice + fibre-BC + continuum builders + calibration
    │   │   ├── excitation.py   # Ground-motion / resonant-sine helpers
    │   │   ├── pushover.py     # Nonlinear pushover (--reference beamcolumn|continuum)
    │   │   ├── pushover_linear.py
    │   │   ├── dynamic.py      # Nonlinear time-history (sine / El Centro)
    │   │   ├── dynamic_linear.py
    │   │   └── single_beamcolumn.py  # vs a single-element fibre beam-column
    │   ├── frame/              # Portal frame studies (column + thinner beam)
    │   │   ├── specimen.py
    │   │   ├── build.py
    │   │   ├── excitation.py
    │   │   ├── pushover.py
    │   │   ├── pushover_linear.py
    │   │   ├── dynamic.py
    │   │   ├── dynamic_linear.py
    │   │   └── visualize.py    # SI-unit modal lattice-vs-continuum study + GIFs
    │   ├── data/
    │   │   └── elcentro.at2    # 1940 El Centro NS ground motion record
    │   └── doc/                # Standalone column report (LaTeX) + model figures
    │       ├── column_report.tex
    │       ├── make_model_figures.py
    │       └── column_report.pdf
    └── tests/                  # pytest suite (26 tests; 1 known-failing)
```

## Dependencies

| Package | Role |
|---|---|
| `openseespy` ≥ 3.8 | OpenSees finite element engine |
| `gmsh` ≥ 4.13 | Node placement and mesh generation |
| `numpy` ≥ 2.0 | Array operations |
| `scipy` ≥ 1.17 | Least-squares / secant calibration |
| `matplotlib` ≥ 3.11 | Visualisation and plots |
| `pillow` ≥ 12.2 | Animated GIF export |

## Status

| Feature | Status |
|---|---|
| Generic FE model + gmsh mesher | ✅ |
| Horizon strut connectivity | ✅ |
| Lattice + continuum builders (2D) | ✅ |
| Static + modal calibration | ✅ |
| Nonlinear RC column (Concrete02 + Steel02) | ✅ |
| Nonlinear dynamic time-history (modal damping) | ✅ |
| 2D continuum RC reference (ASDConcrete3D) | ✅ |
| RC portal frame (pushover + dynamic) | ✅ |
| gmsh-embedded discrete bars in continuum | ⬜ Planned |
| 3D solids / shells | ⬜ Planned |
| Packaged beam-column builder | ⬜ Planned |

## Licence

This project was developed as part of course CE7016. No licence is currently applied — all rights reserved.
