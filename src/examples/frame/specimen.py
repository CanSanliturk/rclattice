"""Shared specimen for the RC frame benchmark studies (Stage-1 elastic + Stage-2 nonlinear).

Reproduces the OpenSees RCFrameGravity -> RCFramePushOver benchmark (1 bay x 1 storey, kip-in):
geometry, gravity/lateral loading, material grades, reinforcement layout and zoning. Shared by
`pushover.py` (Stage 1, elastic) and `pushover_rc.py` (Stage 2, nonlinear RC); `visualize.py` is a
separate SI modal study and defines its own (different) frame.

`OUT` mirrors this package's location under the output tree: figures land in
`examples/output/frame/`. Units: kip, in.
"""

from __future__ import annotations

from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, Rebar, SteelGrade, portal_frame

# --- benchmark geometry / loading (kip, in) ---------------------------------
H, SPAN = 144.0, 360.0          # storey height, bay span
COL, BEAM, THK = 24.0, 24.0, 15.0   # in-plane depths (X cols, Y beam), out-of-plane width
COVER = 1.5
P_COL = 180.0                   # gravity load per column top (kip, downward)
H_REF = 20.0                    # reference lateral load total (= 10 + 10 kip)
DU, TARGET = 0.1, 15.0          # DisplacementControl step / target roof drift (in)
EPS = 1e-6
E_C = 4030.0                    # concrete modulus (ksi)
GF, GFC = 4.0e-3, 1.5           # crack-band fracture energies (kip, in) for Stage-2 softening (D20)
OUT = Path(__file__).resolve().parent.parent / "output" / "frame"

# Stage-1 elastic grade (only E/nu/rho are used in the elastic pipeline). Stage-2 grades carry the
# Concrete02 backbone params; the beam stays elastic (fork 2).
FRAME_ELASTIC = ConcreteGrade(name="RC-frame", E=E_C, nu=0.2, rho=2.25e-7, fc=6.0)
CORE = ConcreteGrade("core", E=E_C, nu=0.2, rho=2.25e-7, fc=6.0, epsc0=0.004, fcu=5.0, epsU=0.014, ft=0.6)
COVER_C = ConcreteGrade("cover", E=E_C, nu=0.2, rho=2.25e-7, fc=5.0, epsc0=0.002, fcu=0.0, epsU=0.006, ft=0.5)
BEAM_C = ConcreteGrade("beam", E=E_C, nu=0.2, rho=2.25e-7)
STEEL = SteelGrade("rebar", fy=60.0, E0=30000.0, b=0.01)


def frame_problem(material: ConcreteGrade) -> Problem:
    domain = portal_frame(height=H, span=SPAN, col_depth=COL, beam_depth=BEAM, thickness=THK)
    half = COL / 2.0
    supports = [BoxSupport(box=(-half - 1.0, SPAN + half + 1.0, -EPS, EPS), fix=(1, 1))]  # fixed bases
    loads = [  # gravity: P_COL down on each column top
        BoxLoad(box=(-half, half, H - EPS, H + EPS), total=(0.0, -P_COL)),
        BoxLoad(box=(SPAN - half, SPAN + half, H - EPS, H + EPS), total=(0.0, -P_COL)),
    ]
    return Problem(ndm=2, ndf=2, domain=domain, material=material, supports=supports, loads=loads)


def lateral_loads(model) -> list[Load]:
    """Reference lateral pattern: H_REF distributed over the beam-level nodes, +X."""
    ids = select_nodes(model, (-COL, SPAN + COL, H - EPS, H + BEAM + EPS))
    return [Load(nid, (H_REF / len(ids), 0.0)) for nid in ids]


def zone_of(x: float, y: float) -> str:
    """Concrete zone of a strut midpoint: elastic beam above the storey; column cover on the
    outer 1.5-in ring of each column's depth; confined core inside."""
    if y >= H - EPS:
        return "beam"
    center = 0.0 if x < SPAN / 2.0 else SPAN
    return "cover" if abs(x - center) >= COL / 2.0 - COVER - EPS else "core"


def rebars() -> tuple[Rebar, ...]:
    bars = []
    for center in (0.0, SPAN):
        for dx, area in ((-10.5, 1.8), (0.0, 1.2), (10.5, 1.8)):
            bars.append(Rebar(path=[(center + dx, 0.0), (center + dx, H)], area=area, steel=STEEL))
    return tuple(bars)
