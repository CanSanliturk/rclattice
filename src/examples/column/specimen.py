"""Shared specimen definition for the single RC cantilever column studies (pushover / dynamic,
nonlinear / linear).

This is the backend-agnostic specimen: geometry, material grades, reinforcement layout, the
`Problem` (supports + axial load), zoning, and the lateral-load shape — everything that is common to
every column model in this package. The model BUILDERS (lattice / fiber beam-column + calibration)
live in `build.py`; the dynamic input lives in `excitation.py`. Units: kip, in.

`OUT` mirrors this package's location under the output tree: figures land in
`examples/output/column/` (the source dir name appended to `examples/output`).
"""

from __future__ import annotations

from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, Rebar, RectangleDomain, SteelGrade

H, W, THK = 144.0, 24.0, 15.0     # height, section depth (X, in-plane), width (out-of-plane)
COVER = 1.5
P, H_REF = 180.0, 10.0            # axial load (down), reference lateral
DU, TARGET = 0.05, 10.0
MESH, EPS = 1.5, 1e-6
E_C = 4030.0
GF, GFC = 4.0e-3, 1.5             # crack-band fracture energies (kip, in), D20
STIRRUP_AREA = 0.3                # transverse (hoop) steel area per horizontal tie strut (D23)
OUT = Path(__file__).resolve().parent.parent / "output" / "column"
HORIZON = 1.5 #3.1

# CORE is the CONFINED concrete (D23): high crushing strain epsU + residual fcu so the core
# stays ductile post-cracking, mimicking transverse confinement (Mander). COVER is unconfined.
CORE = ConcreteGrade("core", E=E_C, nu=0.2, rho=2.25e-7, fc=6.0, epsc0=0.004, fcu=4.0, epsU=0.030, ft=0.6)
COVER_C = ConcreteGrade("cover", E=E_C, nu=0.2, rho=2.25e-7, fc=5.0, epsc0=0.002, fcu=0.0, epsU=0.006, ft=0.5)
STEEL = SteelGrade("rebar", fy=60.0, E0=30000.0, b=0.01)


def column_problem(material: ConcreteGrade) -> Problem:
    half = W / 2.0
    domain = RectangleDomain(length=W, height=H, thickness=THK, origin=(-half, 0.0))
    supports = [BoxSupport(box=(-half - 1.0, half + 1.0, -EPS, EPS), fix=(1, 1))]  # fixed base
    loads = [BoxLoad(box=(-half, half, H - EPS, H + EPS), total=(0.0, -P))]         # axial at top
    return Problem(ndm=2, ndf=2, domain=domain, material=material, supports=supports, loads=loads)


def zone_of(x: float, y: float) -> str:
    return "cover" if abs(x) >= W / 2.0 - COVER - EPS else "core"


def longitudinal_rebars() -> tuple[Rebar, ...]:
    """The three vertical longitudinal bars (x = -10.5, 0, +10.5), as steel struts on shared nodes.
    Shared by the lattice (with stirrups added) and the continuum reference (longitudinal only, D29)."""
    return tuple(Rebar([(dx, 0.0), (dx, H)], area, STEEL)
                 for dx, area in ((-10.5, 1.8), (0.0, 1.2), (10.5, 1.8)))


def rebars() -> tuple[Rebar, ...]:
    """Longitudinal bars (vertical) + transverse stirrup ties (horizontal across the core, one
    per mesh row): the transverse steel gives the lattice a non-softening lateral/shear path so
    the confined core does not disintegrate into a local mechanism past yield (D23)."""
    nrows = int(round(H / MESH))
    stirrups = tuple(Rebar([(-10.5, round(i * MESH, 6)), (10.5, round(i * MESH, 6))], STIRRUP_AREA, STEEL,
                           role="stirrup")
                     for i in range(nrows + 1))
    return longitudinal_rebars() + stirrups


def lateral_loads(model) -> list[Load]:
    ids = select_nodes(model, (-W, W, H - EPS, H + EPS))
    return [Load(nid, (H_REF / len(ids), 0.0)) for nid in ids]
