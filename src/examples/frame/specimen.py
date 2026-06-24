"""Shared specimen for the RC portal-frame studies (pushover / dynamic, nonlinear / linear).

The frame is CONSTITUTED FROM the single RC cantilever column (examples/column): two of that exact
column (24-in deep x 144-tall, confined CORE + unconfined COVER Concrete02, 3+2+3 Steel02
longitudinal bars + horizontal stirrups) connected at the top by a THINNER 18-in beam that uses the
SAME concrete grades and steel grade — top + bottom longitudinal bars and vertical stirrups. The
result is a one-bay, one-storey portal frame. Everything backend-agnostic lives here: geometry,
material grades, reinforcement layout, the `Problem` (fixed bases + axial gravity per column), zoning,
and the lateral-load shape. The model BUILDERS (lattice / fiber beam-column frame / continuum +
calibration) live in `build.py`; the dynamic input lives in `excitation.py`. Units: kip, in.

`OUT` mirrors this package's location under the output tree: figures land in
`examples/output/frame/` (the source dir name appended to `examples/output`).
"""

from __future__ import annotations

from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, Rebar, SteelGrade, portal_frame

# --- geometry / loading (kip, in) -------------------------------------------
H, SPAN = 144.0, 144.0           # storey height (= column height), bay span (col centerline to centerline)
COL, BEAM, THK = 24.0, 18.0, 15.0  # column depth (X), thinner beam depth (Y), out-of-plane width (both)
COVER = 1.5
P, H_REF = 180.0, 20.0           # axial gravity per column (down), reference lateral total (pattern only)
DU, TARGET = 0.1, 15.0           # DisplacementControl step / target roof drift (in)
MESH, EPS = 1.5, 1e-6
HORIZON = 1.5
E_C = 4030.0
GF, GFC = 4.0e-3, 1.5            # crack-band fracture energies (kip, in), D20
STIRRUP_AREA = 0.3              # transverse (stirrup/tie) steel area per tie strut (D23)
BEAM_TOP_AREA = BEAM_BOT_AREA = 1.8  # beam top/bottom longitudinal layer area (= column corner bars)
OUT = Path(__file__).resolve().parent.parent / "output" / "frame"

# Same grades as the cantilever column: CORE is the CONFINED concrete (high crushing strain epsU +
# residual fcu so the core stays ductile post-cracking), COVER_C is unconfined. STEEL is shared by
# the columns' bars/stirrups and the beam's bars/stirrups.
CORE = ConcreteGrade("core", E=E_C, nu=0.2, rho=2.25e-7, fc=6.0, epsc0=0.004, fcu=4.0, epsU=0.030, ft=0.6)
COVER_C = ConcreteGrade("cover", E=E_C, nu=0.2, rho=2.25e-7, fc=5.0, epsc0=0.002, fcu=0.0, epsU=0.006, ft=0.5)
STEEL = SteelGrade("rebar", fy=60.0, E0=30000.0, b=0.01)

_BEAM_TOP = H + BEAM - COVER     # y of the beam's top longitudinal layer  (= 160.5)
_BEAM_BOT = H + COVER            # y of the beam's bottom longitudinal layer (= 145.5)
_BEAM_X0, _BEAM_X1 = -COL / 2.0 + COVER, SPAN + COL / 2.0 - COVER  # beam bar span inside the cover


def frame_problem(material: ConcreteGrade) -> Problem:
    """The portal-frame Problem: two columns (centerlines x=0 and x=SPAN) + a thinner beam on top,
    fixed bases, constant axial gravity P at each column top (beam-column joint)."""
    domain = portal_frame(height=H, span=SPAN, col_depth=COL, beam_depth=BEAM, thickness=THK)
    half = COL / 2.0
    supports = [BoxSupport(box=(-half - 1.0, SPAN + half + 1.0, -EPS, EPS), fix=(1, 1))]  # fixed bases
    loads = [  # gravity: P down on each column top
        BoxLoad(box=(-half, half, H - EPS, H + EPS), total=(0.0, -P)),
        BoxLoad(box=(SPAN - half, SPAN + half, H - EPS, H + EPS), total=(0.0, -P)),
    ]
    return Problem(ndm=2, ndf=2, domain=domain, material=material, supports=supports, loads=loads)


def zone_of(x: float, y: float) -> str:
    """Concrete zone of a strut midpoint / quad centroid. In the columns the cover is the outer
    1.5-in ring of the section depth (X); in the beam it is the outer 1.5-in ring of the (vertical)
    beam depth (Y). Both columns and the beam use the SAME confined CORE / unconfined COVER grades —
    but the beam zones are returned DISTINCTLY ("beam_core"/"beam_cover") so a builder can model the
    beam concrete elastically (the stable default) or nonlinearly (opt-in) independent of the columns
    (the thin nonlinear beam softens into a local lattice mechanism the static pushover can't trace)."""
    if y >= H - EPS:  # beam region (above the storey)
        return "beam_cover" if (y - H <= COVER + EPS or H + BEAM - y <= COVER + EPS) else "beam_core"
    center = 0.0 if x < SPAN / 2.0 else SPAN
    return "cover" if abs(x - center) >= COL / 2.0 - COVER - EPS else "core"


def _column_longitudinal() -> list[Rebar]:
    """The two columns' vertical 3+2+3 longitudinal bars (x = center +- 10.5 and center), full
    height, identical to the cantilever column's longitudinal layout."""
    bars = []
    for center in (0.0, SPAN):
        for dx, area in ((-10.5, 1.8), (0.0, 1.2), (10.5, 1.8)):
            bars.append(Rebar([(center + dx, 0.0), (center + dx, H)], area, STEEL))
    return bars


def _beam_longitudinal() -> list[Rebar]:
    """The beam's top + bottom horizontal longitudinal layers (1.8 in^2 each), running the full
    beam length inside the cover, on shared lattice/continuum nodes (perfect bond)."""
    return [
        Rebar([(_BEAM_X0, _BEAM_TOP), (_BEAM_X1, _BEAM_TOP)], BEAM_TOP_AREA, STEEL),
        Rebar([(_BEAM_X0, _BEAM_BOT), (_BEAM_X1, _BEAM_BOT)], BEAM_BOT_AREA, STEEL),
    ]


def longitudinal_rebars() -> tuple[Rebar, ...]:
    """All longitudinal bars (columns 3+2+3 + beam top/bottom), no stirrups — shared by the lattice
    (which then adds stirrups) and the 2D continuum reference (longitudinal only: the continuum
    itself supplies the lateral/shear path the lattice gets from stirrups, D29)."""
    return tuple(_column_longitudinal() + _beam_longitudinal())


def rebars() -> tuple[Rebar, ...]:
    """Longitudinal bars + transverse stirrups (the lattice reinforcement). Column stirrups: a
    horizontal tie across each column core, one per mesh row. Beam stirrups: a vertical tie across
    the beam core (bottom bar to top bar), one per mesh column along the clear span. The transverse
    steel gives the lattice a non-softening lateral/shear path past yield (D23)."""
    stirrups: list[Rebar] = []
    nrows = int(round(H / MESH))
    for center in (0.0, SPAN):                       # column ties (horizontal, across the core)
        for i in range(nrows + 1):
            y = round(i * MESH, 6)
            stirrups.append(Rebar([(center - 10.5, y), (center + 10.5, y)], STIRRUP_AREA, STEEL,
                                  role="stirrup"))
    j0, j1 = int(round(COL / 2.0 / MESH)), int(round((SPAN - COL / 2.0) / MESH))  # beam clear span
    for j in range(j0, j1 + 1):                      # beam ties (vertical, across the beam core)
        x = round(j * MESH, 6)
        stirrups.append(Rebar([(x, _BEAM_BOT), (x, _BEAM_TOP)], STIRRUP_AREA, STEEL, role="stirrup"))
    return longitudinal_rebars() + tuple(stirrups)


def lateral_loads(model) -> list[Load]:
    """Reference lateral pattern: H_REF distributed over the beam-level nodes, +X (DisplacementControl
    drives the magnitude, so only the shape matters)."""
    ids = select_nodes(model, (-COL, SPAN + COL, H - EPS, H + BEAM + EPS))
    return [Load(nid, (H_REF / len(ids), 0.0)) for nid in ids]


def control_base_nodes(model) -> tuple[int, list[int]]:
    """The pushover control node (top-left beam-column joint, x=0 y=H) and the base nodes (both
    fixed column bases), queried off a built model."""
    ctrl = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(model, (-COL, SPAN + COL, -EPS, EPS))
    return ctrl, base


def add_axial_mass(model, per_column_mass: float) -> None:
    """Lump `per_column_mass` (= P/g, the axial gravity treated as tributary seismic mass) onto EACH
    column's top (beam-joint) nodes, in place — the dynamic counterpart of the gravity BoxLoads."""
    half = COL / 2.0
    for cx in (0.0, SPAN):
        ids = select_nodes(model, (cx - half - EPS, cx + half + EPS, H - EPS, H + EPS))
        for nid in ids:
            mx, my = model.masses[nid]
            model.masses[nid] = (mx + per_column_mass / len(ids), my + per_column_mass / len(ids))
