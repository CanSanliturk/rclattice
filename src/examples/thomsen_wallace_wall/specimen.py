"""Shared specimen definition for RW2, the slender rectangular wall of Thomsen & Wallace.

TEST SOURCE:  Thomsen, J. H., and J. W. Wallace (1995), "Displacement based design of reinforced
concrete structural walls: an experimental investigation of walls with rectangular and T-shaped
cross-sections", Report CU/CEE-95/06, Clarkson University; and Thomsen & Wallace (2004),
"Displacement-based design of slender reinforced concrete structural walls — experimental
verification", J. Struct. Eng. 130(4): 618-630.  **NEITHER IS IN THE REPO.**

LATTICE SOURCE:  Aydin, B. B., K. Tuncay, and B. Binici (2019), "Simulation of Reinforced Concrete
Member Response Using Lattice Model", J. Struct. Eng. 145(9): 04019091 — Table 1 (properties),
Table 2 (lattice size), Table 4 (measured vs computed), Fig. 4(e) (their model) and Fig. 9
(specimen, response, damage). The PDF sits in `aydins_paper/`.

EVERY NUMBER HERE IS SECOND-HAND THROUGH THE 2019 PAPER, as for the Aldemir wall (user decision,
2026-09-20: build from the 2019 paper; the primary sources were not supplied). `testdata.py` records
which table or figure each value came from and what the 2019 paper does NOT print — the measured
steel strengths, the loading protocol, the cycle count, the failure mode, the strain-hardening ratio.
Where the physical test is described below from general knowledge of RW2 it is marked RECOLLECTION
and is not used by the model.

WHY THIS SPECIMEN. It is the slender counterpart of the Aldemir wall: aspect ratio 3.0 against
0.75, a CONSTANT AXIAL LOAD of 378 kN (0.07 A_g f_c) where Aldemir had none, confined boundary
elements with hoops at 76 mm, and a flexural mechanism (Fig. 9(c): horizontal flexural cracks over
the lower third). It is the specimen the 2019 paper predicts BEST on force — 1.04x at horizon 1.5 —
against 1.21x on Aldemir, so the pair brackets the method: if the lattice reproduces RW2's peak AND
its post-yield plateau with the same modelling choices that gave Aldemir 0.997x (D97), the method
transfers across aspect ratio and axial load; the 2019 paper itself notes "overestimation of response
in the postyielding region", which is where the D101/D102 hardening finding lands.

GEOMETRY (Fig. 9(a), dimensions in mm — the drawing rounds the test's inch values: 1220 = 48 in,
3660 = 144 in, 102 = 4 in, 191 = 7.5 in, 153 = 6 in, 51 = 2 in, 19 = 3/4 in, 76 = 3 in, 4.76 = 3/16 in).

Units: N, mm (stresses MPa, forces N).
"""

from __future__ import annotations

import math
from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, Rebar, RectangleDomain, SteelGrade

SPECIMEN = "Thomsen & Wallace (1995, 2004) wall RW2"

# --- geometry (2019 paper, Fig. 9a) ---------------------------------------------------------------
LW = 1220.0        # wall length (in-plane), Fig. 9(a) "1220"
HW = 3660.0        # wall height to the load point = shear span, Fig. 9(a) "3660"
TW = 102.0         # thickness, Fig. 9(c) "102 mm"
A_SHEAR = HW       # drift denominator: drift = u_top / HW (the paper's 72 mm max = 1.97%)
ASPECT = HW / LW   # 3.00 — "aspect ratio 1:3" in the paper's own words

EPS = 1e-6
OUT = Path(__file__).resolve().parent.parent / "output" / "thomsen_wallace_wall"

# --- reinforcement (Fig. 9a) ----------------------------------------------------------------------
#
# BOUNDARY ELEMENTS: "8-#3 bars" at each end, drawn as four bar positions "3@51" from a "19" cover —
# so x = 19, 70, 121, 172 and their mirror 1048, 1099, 1150, 1201 — two curtains (8 = 4 x 2). Hoops
# "2-db = 4.76 mm (@ 76 mm)" close the boundary element: in the plane model they are the horizontal
# legs between the outer boundary bars, at 76 mm up the height.
#
# WEB: "8-#2 bars" vertical — four positions x two curtains — and "2-#2 bars (@191 mm)" horizontal.
# The vertical web positions are NOT dimensioned individually. The drawing's chain reads
# 19 | 153 | 3@191 | 153 | 19, which sums to 917, not 1220: the remaining 303 mm are the two gaps
# between the last boundary bar and the first web bar, taken as EQUAL (151.5 each), i.e. the four
# web bars at 191 mm centred between the boundary elements (x = 323.5, 514.5, 705.5, 896.5). This
# matches a 7.5 in web spacing centred in the 33 in web of the test — INFERRED, and marked so.
#
# The HORIZONTAL web bars' start height is not dimensioned; the first is placed at half a spacing
# (95.5 mm) and the rest every 191 mm (19 bars, the last at 3533.5) — ASSUMED. The hoops start at
# 38 mm (half a spacing): 48 hoops, the last at 3610 — ASSUMED.
COVER_X = 19.0
BE_PITCH = 51.0                       # "3@51"
N_BE_POS = 4                          # bar positions per boundary element
WEB_PITCH = 191.0                     # "@191 mm" — both directions
HOOP_PITCH = 76.0                     # "@ 76 mm"

X_BE_LEFT = tuple(COVER_X + i * BE_PITCH for i in range(N_BE_POS))            # 19 .. 172
X_BE_RIGHT = tuple(LW - x for x in reversed(X_BE_LEFT))                       # 1048 .. 1201
_web_span = X_BE_RIGHT[0] - X_BE_LEFT[-1]                                     # 876
_gap = (_web_span - 3.0 * WEB_PITCH) / 2.0                                    # 151.5
X_WEB = tuple(X_BE_LEFT[-1] + _gap + i * WEB_PITCH for i in range(4))         # 323.5 .. 896.5
X_BARS = X_BE_LEFT + X_WEB + X_BE_RIGHT                                       # 12 vertical lines

Y_HOOPS = tuple(HOOP_PITCH / 2.0 + i * HOOP_PITCH for i in range(int(HW // HOOP_PITCH)))   # 48
Y_WEB = tuple(WEB_PITCH / 2.0 + i * WEB_PITCH for i in range(int(HW // WEB_PITCH)))        # 19

# Bar areas (US sizes; the 2019 paper names the sizes, not the areas):
#   #3: d_b = 9.53 mm, 71 mm^2      #2: d_b = 6.35 mm, 32 mm^2      3/16 in wire: 4.76 mm, 17.8 mm^2
# Two curtains collapse onto one in-plane line per position (the standard plane-model collapse, as
# for WSH3, VK3 and Aldemir), so every in-plane area below is TWICE the bar area.
A_NO3 = math.pi * 9.53 ** 2 / 4.0          # 71.3 mm^2
A_NO2 = math.pi * 6.35 ** 2 / 4.0          # 31.7 mm^2
A_HOOP = math.pi * 4.76 ** 2 / 4.0         # 17.8 mm^2
A_BE_LINE = 2.0 * A_NO3                    # 142.6 mm^2 per boundary-bar position
A_WEB_LINE = 2.0 * A_NO2                   # 63.3 mm^2 per web position (vertical AND horizontal)
A_HOOP_LINE = 2.0 * A_HOOP                 # 35.6 mm^2 — two legs of one hoop in the plane

# --- materials (2019 paper, Table 1, "Wall / Thomsen and Wallace (1995, 2004)" column) ------------
FC = 42.8                       # MPa
FT = 2.03                       # MPa — NO footnote in Table 1; equals 0.31*sqrt(f_c) to three figures
EC = 31030.0                    # MPa — NO footnote; 4700*sqrt(f_c) would be 30,749, so likely a
                                #       measured or otherwise-derived value, quoted as printed
GF = 0.075                      # N/mm  (the table's 75 N/m)
NU = 0.20                       # not given by the paper; the repo's convention for every wall
RHO_C = 2.4e-9                  # t/mm^3

FY = 414.0                      # MPa — Table 1. The nominal Grade 60 (60 ksi = 413.7); the test's
                                #       measured coupon strengths are NOT in the 2019 paper
ES = 200000.0                   # MPa
N_AXIAL = 378.0e3               # N — Table 1 "N (kN) 378" = 0.071 A_g f_c, held constant

# epsc0 DERIVED, never quoted (D56): Concrete02's initial compressive tangent is 2*fc/epsc0
# regardless of the grade's E, so any other value silently mismatches the calibrated modulus.
EPSC0 = 2.0 * FC / EC           # 0.002759

# Aydin's own trilinear tension-softening shape for this specimen, Table 1: a1 = 1.5, a2 = 80,
# a3 = 350 (fitted at HIS 19 mm grid). a3/a2 = 4.375 here, against 5.14 for Aldemir.
AYDIN_A1, AYDIN_A2, AYDIN_A3 = 1.5, 80.0, 350.0

GRADES = {
    "wall": ConcreteGrade(name="wall", E=EC, nu=NU, rho=RHO_C, fc=FC, ft=FT, Gf=GF,
                          epsc0=EPSC0, fcu=0.2 * FC, epsU=8.0 * EPSC0),
}
# ONE steel grade for every bar: Table 1 prints a single f_y. b = 0.01 is the repo convention
# (D101/D102: unprinted, and the largest single control on drift capacity).
STEEL = SteelGrade(name="Gr60", fy=FY, E0=ES, b=0.01)

# --- discretisation ---------------------------------------------------------------------------------
#
# TWO GRID MODES (user decision, 2026-09-20; D104):
#   "uniform" — the structured grid every earlier study used. `mesh` must divide 1220 and 3660
#               (10, 12.2, 15.25, 20, 30.5, 61); bar lines SNAP to the nearest grid line.
#   "rebar"   — a GRADED structured grid whose x-lines are the 12 vertical bar axes and whose
#               y-lines are the 48 hoop axes, each gap filled at ~`mesh`; bars sit on nodes exactly.
#               The horizontal web bars (191 pitch) snap to the nearest y-line (<= half a cell).
# The default is the graded grid at a 25 mm target — close to Aydin's own 19 mm without the cost.
MESH, HORIZON = 25.0, 1.5
GRID = "rebar"
UNIFORM_MESHES = (10.0, 12.2, 15.25, 20.0, 30.5, 61.0)

# --- drive settings (repo conventions, D62/D64; NOT measured for this specimen) -------------------
QUASI_STATIC_RATE = 7.6         # mm/s
DAMPING_RATIO = 0.5


def zone_of(_x: float, _y: float) -> str:
    """One concrete zone: the paper reports a single f_c and models a single rectangle. The
    confined boundary elements get NO separate grade — the hoops are in-plane rebar struts that
    restrain the compressed boundary themselves, so a Mander grade would double-count (D66)."""
    return "wall"


def wall_problem(*, length: float = LW, height: float = HW) -> Problem:
    """The wall as a `Problem`: fixed base, the 378 kN axial load spread over the top row.

    The axial load IS the Problem's load — the runners apply the model's own loads as the constant
    gravity case before any lateral drive (`opensees._gravity_loads`), which is how WSH3 and VK3
    carry theirs. The lateral pattern is supplied separately by `lateral_loads`.
    """
    return Problem(
        ndm=2, ndf=2,
        domain=RectangleDomain(length=length, height=height, thickness=TW),
        material=GRADES["wall"],
        supports=(BoxSupport(box=(-EPS, length + EPS, -EPS, EPS), fix=(1, 1)),),
        loads=(BoxLoad(box=(-EPS, length + EPS, height - EPS, height + EPS),
                       total=(0.0, -N_AXIAL)),),
    )


# --- grids ------------------------------------------------------------------------------------------

def snap(v: float, lines) -> float:
    """The grid line nearest `v`."""
    return min(lines, key=lambda g: abs(g - v))


def grid_lines(grid: str = GRID, mesh_size: float = MESH, *, length: float = LW,
               height: float = HW) -> tuple[list[float], list[float]]:
    """The x- and y-lines of the node grid for either mode."""
    from rclattice.mesh import graded_lines
    if grid == "uniform":
        check_mesh_alignment(mesh_size, length=length, height=height)
        nx, ny = int(round(length / mesh_size)), int(round(height / mesh_size))
        return ([length * i / nx for i in range(nx + 1)], [height * j / ny for j in range(ny + 1)])
    if grid == "rebar":
        return (graded_lines(length, mesh_size, X_BARS), graded_lines(height, mesh_size, Y_HOOPS))
    raise SystemExit(f"unknown grid {grid!r}; use 'uniform' or 'rebar'")


def bar_layout(grid: str = GRID, mesh_size: float = MESH) -> dict:
    """The bar lines AS MODELLED on this grid, and how far each was moved from the drawing."""
    xs, ys = grid_lines(grid, mesh_size)
    xv = [snap(x, xs) for x in X_BARS]
    yh = [snap(y, ys) for y in Y_HOOPS]
    yw = [snap(y, ys) for y in Y_WEB]
    return {
        "x_vertical": xv, "x_vertical_drawn": list(X_BARS),
        "x_offset_max": max(abs(a - b) for a, b in zip(xv, X_BARS)),
        "y_hoops": yh, "y_hoops_offset_max": max(abs(a - b) for a, b in zip(yh, Y_HOOPS)),
        "y_web": yw, "y_web_offset_max": max(abs(a - b) for a, b in zip(yw, Y_WEB)),
        "n_x_lines": len(xs), "n_y_lines": len(ys),
        "spacing_x": (min(b - a for a, b in zip(xs, xs[1:])), max(b - a for a, b in zip(xs, xs[1:]))),
        "spacing_y": (min(b - a for a, b in zip(ys, ys[1:])), max(b - a for a, b in zip(ys, ys[1:]))),
    }


def rebars(grid: str = GRID, mesh_size: float = MESH, *, length: float = LW, height: float = HW,
           full_height: bool = True) -> tuple[Rebar, ...]:
    """Every bar as an in-plane line on the chosen grid.

    Vertical bars are tagged "longitudinal"; hoops and horizontal web bars "stirrup" (the
    visualizer's split). `full_height` runs the vertical bars to the top face (D78/D84 — the
    default here, since the actuator loads the top row and a plain-concrete band under it tears).
    """
    lay = bar_layout(grid, mesh_size)
    xs, ys = grid_lines(grid, mesh_size, length=length, height=height)
    top = height if full_height else snap(height - 2 * mesh_size, ys)
    bars: list[Rebar] = []
    for x, x0 in zip(lay["x_vertical"], X_BARS):
        area = A_BE_LINE if (x0 in X_BE_LEFT or x0 in X_BE_RIGHT) else A_WEB_LINE
        bars.append(Rebar(path=[(x, 0.0), (x, top)], area=area, steel=STEEL, role="longitudinal"))
    xl, xr = lay["x_vertical"][0], lay["x_vertical"][-1]
    xl_in, xr_in = lay["x_vertical"][N_BE_POS - 1], lay["x_vertical"][-N_BE_POS]
    for y in lay["y_hoops"]:                     # hoop legs across each boundary element
        bars.append(Rebar(path=[(xl, y), (xl_in, y)], area=A_HOOP_LINE, steel=STEEL, role="stirrup"))
        bars.append(Rebar(path=[(xr_in, y), (xr, y)], area=A_HOOP_LINE, steel=STEEL, role="stirrup"))
    for y in lay["y_web"]:                       # horizontal web bars, outer bar to outer bar
        bars.append(Rebar(path=[(xl, y), (xr, y)], area=A_WEB_LINE, steel=STEEL, role="stirrup"))
    return tuple(bars)


def check_mesh_alignment(mesh_size: float, *, length: float = LW, height: float = HW) -> None:
    """The UNIFORM grid must divide the panel; the bar lines are allowed to snap."""
    bad = [(n, v) for n, v in (("length", length), ("height", height))
           if abs(v / mesh_size - round(v / mesh_size)) > 1e-9]
    if bad:
        raise SystemExit(
            f"mesh {mesh_size:g} does not divide the " + ", ".join(f"{n} ({v:g})" for n, v in bad)
            + f" of the {length:g} x {height:g} panel; the uniform grid accepts "
            + ", ".join(f"{m:g}" for m in UNIFORM_MESHES) + " — or use --grid rebar")


# --- node selectors (post-build queries) ------------------------------------------------------------
# Driven across the WHOLE top row, as every wall in the repo is (no loading beam is modelled).

def top_nodes(model, *, height: float = HW) -> list[int]:
    """The loaded row — CONCRETE nodes only (the duplicate steel nodes of a bonded model follow
    through their links; driving them would pin zero slip at the loaded row, D78/D85)."""
    return select_nodes(model, (-LW, 2.0 * LW, height - EPS, height + EPS), kind="concrete")


def lateral_loads(model, *, height: float = HW) -> list[Load]:
    ids = top_nodes(model, height=height)
    return [Load(nid, (1.0 / len(ids), 0.0)) for nid in ids]


def control_node(model, *, height: float = HW) -> int:
    """Drift is read at the top-left corner of the loaded edge."""
    return select_nodes(model, (-EPS, EPS, height - EPS, height + EPS), kind="concrete")[0]


def drive_nodes(model, *, height: float = HW) -> list[int]:
    return top_nodes(model, height=height)


def base_nodes(model) -> list[int]:
    """The reaction set — EVERY node on the fixed row, steel duplicates included (D85)."""
    return select_nodes(model, (-LW, 2.0 * LW, -EPS, EPS), kind="any")
