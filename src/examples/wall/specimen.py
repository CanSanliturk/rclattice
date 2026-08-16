"""Shared specimen definition for the flexure-controlled RC shear-wall study.

The specimen is SW-NC-FF from Sahinkaya, Binbir, Orakcal & Ilki (2025), "Behavior of
Nonconforming Flexure-Controlled RC Structural Walls Under Reversed Cyclic Lateral Loading",
Buildings 15, 4501 — a full-scale, relatively slender (hw/lw = 3.0) nonconforming wall with
low-strength concrete, plain longitudinal bars and NO boundary confinement. Switch to its sibling
SW-IC-FF (20% boundary length, boundary ties) via `BOUNDARY_RATIO`/`boundary_rebars`.

This is the backend-agnostic half: geometry (wall + pedestal), concrete/steel grades,
reinforcement layout, the `Problem` (base fixity + axial load), material zoning, and the lateral
load shape. Model building and the Aydin energy-balance calibration live in `build.py`.

Units: N, mm (so stresses are MPa, forces N). `OUT` mirrors this package's location under the
output tree: figures land in `examples/output/wall/`.

SCOPE NOTE — the test's dominant deformation mechanism is ROCKING driven by progressive debonding
of the PLAIN longitudinal bars (74-76% of tip displacement at 2% drift; paper Fig. 21). The
lattice ties rebar to concrete at shared nodes, i.e. PERFECT bond (D5/D13), so it cannot reproduce
that mechanism. This does not affect the present ELASTIC calibration — rocking engages only after
cracking — but it caps what a later nonlinear/cyclic study on this specimen may claim.
"""

from __future__ import annotations

import math
from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, CompoundRectangles, ConcreteGrade, Problem, Rebar, SteelGrade

# --- geometry (paper Sec. 2.1, Fig. 3b) ---------------------------------------------------------
LW, HW, TW = 1000.0, 3000.0, 200.0        # wall length (in-plane), height, thickness (out-of-plane)
PED_L, PED_H, PED_W = 1900.0, 600.0, 800.0  # pedestal length (in-plane), depth (height), width
A_SHEAR = 2200.0                          # lateral load height above the W-P interface (SSDR = 2.2)
UPPER_Y = 2000.0                          # top 1 m is the stiffer, higher-strength loading head

MESH, HORIZON = 50.0, 1.5                 # grid spacing; horizon (D9) — 8 neighbours per node
EPS = 1e-6
OUT = Path(__file__).resolve().parent.parent / "output" / "wall"

# --- loading (paper Eq. 2 and Sec. 2.2) ---------------------------------------------------------
N_AXIAL = 600.0e3                         # constant axial compression, N (= 0.20 * Ag * fc')
H_REF = 100.0e3                           # reference lateral load pattern magnitude, N

# --- concrete grades (paper Sec. 2.1: three casts, Ec per Eq. 1) --------------------------------
# The pedestal is 800 mm wide against the wall's 200 mm, but a 2D plane model carries ONE
# thickness. Its extra width is folded into an equivalent modulus (E * PED_W / TW) so the
# constant-thickness model reproduces the real axial/flexural rigidity of the block. The paper
# measured pedestal sliding and rigid-body rotation as negligible, which this reproduces.
NU = 0.20
RHO = 2.4e-9                              # 2400 kg/m^3 in N-mm-s (tonne/mm^3)

# Compression backbones follow the measured cylinder curves (paper Fig. 6): the test unit peaks
# near eps=0.0022 and has lost ~15% of its strength by eps=0.004, falling to a few MPa by 0.012.
TEST_UNIT = ConcreteGrade("test-unit", E=16100.0, nu=NU, rho=RHO,
                          fc=15.1, epsc0=0.0022, fcu=4.0, epsU=0.012, ft=1.5)
UPPER = ConcreteGrade("upper", E=26400.0, nu=NU, rho=RHO,
                      fc=43.4, epsc0=0.0025, fcu=10.0, epsU=0.008, ft=3.3)
PEDESTAL = ConcreteGrade("pedestal", E=23400.0 * (PED_W / TW), nu=NU, rho=RHO,
                         fc=33.7, epsc0=0.0022, fcu=8.0, epsU=0.010, ft=2.7)

GRADES = {"test-unit": TEST_UNIT, "upper": UPPER, "pedestal": PEDESTAL}

# Tensile fracture energy per CEB-FIP MC90, Gf = Gf0*(fcm/10)^0.7 with Gf0 = 0.030 N/mm
# (16 mm max aggregate). Low-strength test-unit concrete gets a correspondingly low Gf.
GF = {"test-unit": 0.053, "upper": 0.093, "pedestal": 0.080}   # N/mm

# --- steel grades (paper Table 2: plain bars, per diameter) -------------------------------------
# Hardening ratio b is MEASURED from Table 2 rather than assumed: these are plain bars with a long
# yield plateau (eps_h >> eps_y), so the secant from yield to (f_max, eps_max) is the honest
# bilinear slope -- b = (f_max - fy) / (Es * (eps_max - eps_y)). It lands near 0.002-0.004, well
# below the 0.01 rule of thumb, which would overstate post-yield strength.
S14 = SteelGrade("phi14", fy=294.0, E0=198600.0, b=0.0024)   # f_max 396 @ eps_max 21.4%
S12 = SteelGrade("phi12", fy=311.0, E0=230300.0, b=0.0021)   # f_max 414 @ eps_max 21.0%
S10 = SteelGrade("phi10", fy=340.0, E0=216600.0, b=0.0019)   # f_max 425 @ eps_max 21.3%
S8 = SteelGrade("phi8", fy=378.0, E0=224800.0, b=0.0038)     # f_max 522 @ eps_max 16.9%
# phi16 (loading head only) is not in Table 2; phi14's properties are the nearest measured proxy.
S16 = SteelGrade("phi16", fy=294.0, E0=198600.0, b=0.0024)


def _bar_area(diameter: float, count: int) -> float:
    return count * math.pi * diameter * diameter / 4.0


# Two reinforcement curtains (front/back) collapse onto one line in a 2D plane model, so each
# in-plane bar line carries BOTH curtains' area.
BOUNDARY_AREA = _bar_area(14.0, 4)   # 4phi14 per boundary region
WEB_V_AREA = _bar_area(12.0, 2)      # phi12@200 vertical web, two curtains   -> rho_v = 0.57%
WEB_H_AREA = _bar_area(10.0, 2)      # phi10@200 horizontal web, two curtains -> rho_t = 0.39%

TIE_AREA = _bar_area(8.0, 2)         # phi8 hoop: 2 in-plane legs (SW-IC-FF boundary only)
HEAD_AREA = _bar_area(16.0, 2)       # phi16@100 loading-head reinforcement, Fig. 3c

# --- specimen variant (paper Table 1) -----------------------------------------------------------
# SW-NC-FF: boundary length 10% of lw, rho_b 2.37%, rho_con = "NC" -> NO boundary confinement.
# SW-IC-FF: boundary length 20%, rho_b 1.54%, rho_con 0.63% -> phi8@200 boundary hoops.
# The "closed hoop / 90-degree hook" columns of Table 1 describe the HORIZONTAL WEB reinforcement's
# anchorage, not the boundary — SW-NC-FF's web hoops close on themselves with 90-degree hooks, which
# is what lets them open up under load and is a defining feature of the nonconforming detail.
SPECIMEN = "SW-NC-FF"
_VARIANTS = {
    "SW-NC-FF": {"boundary_ratio": 0.10, "confined": False, "tie_spacing": None},
    "SW-IC-FF": {"boundary_ratio": 0.20, "confined": True, "tie_spacing": 200.0},
}
BOUNDARY_RATIO = _VARIANTS[SPECIMEN]["boundary_ratio"]
CONFINED_BOUNDARY = _VARIANTS[SPECIMEN]["confined"]
TIE_SPACING = _VARIANTS[SPECIMEN]["tie_spacing"]

WEB_SPACING = 200.0
WEB_V_X = (-300.0, -100.0, 100.0, 300.0)   # vertical web bar lines, @200 inside the web

# Loading head (Fig. 3c): phi16@100 over the top 500 mm, with welded joints. Well above the
# actuator level, so it stiffens the head rather than participating in the flexural response.
HEAD_DEPTH, HEAD_SPACING = 500.0, 100.0

# --- anchorage into the pedestal (Fig. 3a/b) ----------------------------------------------------
# The longitudinal bars do NOT stop at the wall-pedestal interface: they run down through the
# pedestal and turn 90 degrees OUTWARD near the soffit, forming the L-shaped hooks visible in
# Fig. 3. Modelling the hook matters twice over — it stops every bar force being dumped into the
# single node row at y = 0 (the section where flexural damage localizes), and the horizontal leg
# gives those pedestal nodes lateral restraint a straight bar cannot.
ANCHORAGE_DEPTH = 500.0              # bar runs to y = -500, i.e. 100 mm above the soffit
HOOK_LENGTH = 300.0                  # horizontal leg, turned away from the wall centreline


def wall_problem() -> Problem:
    """Wall on its pedestal, fixed at the pedestal soffit, with the constant axial load on top.

    The pedestal is meshed as a second rectangle sharing the y = 0 node row with the wall
    (coincident nodes merged), so the W-P interface is continuous. Fixing the pedestal SOFFIT
    rather than the wall base lets the base region deform, which is what modelling the pedestal
    buys over a clamped wall.
    """
    domain = CompoundRectangles(
        rects=[
            (-PED_L / 2.0, -PED_H, PED_L, PED_H),   # pedestal
            (-LW / 2.0, 0.0, LW, HW),               # wall
        ],
        thickness=TW,
    )
    supports = [BoxSupport(box=(-PED_L, PED_L, -PED_H - EPS, -PED_H + EPS), fix=(1, 1))]
    loads = [BoxLoad(box=(-LW / 2.0, LW / 2.0, HW - EPS, HW + EPS), total=(0.0, -N_AXIAL))]
    return Problem(ndm=2, ndf=2, domain=domain, material=TEST_UNIT, supports=supports, loads=loads)


def zone_of(x: float, y: float) -> str:
    """Material zone at a point: the three separate concrete casts of the specimen."""
    if y < 0.0:
        return "pedestal"
    return "upper" if y >= UPPER_Y else "test-unit"


def boundary_x() -> float:
    """In-plane centroid of a boundary longitudinal group."""
    return LW / 2.0 - BOUNDARY_RATIO * LW / 2.0


def boundary_rebars() -> tuple[Rebar, ...]:
    """The 4phi14 boundary longitudinal groups, lumped at each boundary region's centroid.

    Lumping at the centroid preserves the group's first moment about the section, so the flexural
    lever arm — the quantity that matters for a flexure-controlled wall — is exact.

    The bars run from the hooked anchorage in the pedestal up to the wall top — see `_hooked_path`.
    """
    xc = boundary_x()
    return tuple(Rebar(_hooked_path(dx), BOUNDARY_AREA, S14) for dx in (-xc, xc))


def _hooked_path(x: float) -> list[tuple[float, float]]:
    """A vertical bar with the Fig. 3 L-hook: down the wall, into the pedestal, then turned out.

    The horizontal leg is turned AWAY from the wall centreline, matching the drawing. Every vertex
    is a whole number of 50 mm cells, so the polyline lands on lattice nodes.
    """
    toe = x + HOOK_LENGTH * (1.0 if x >= 0.0 else -1.0)
    return [(toe, -ANCHORAGE_DEPTH), (x, -ANCHORAGE_DEPTH), (x, HW)]


def boundary_ties() -> tuple[Rebar, ...]:
    """Boundary confinement hoops — SW-IC-FF ONLY. Empty for SW-NC-FF, which has none.

    Table 1 gives rho_con = 0.63% for SW-IC-FF and "NC" for SW-NC-FF: the nonconforming unit has NO
    confinement reinforcement in its boundary regions at all. That absence is the point of the
    specimen, and it is why its compression toe crushed. Adding hoops here would model the wrong
    wall and would quietly improve exactly the behaviour the test was designed to expose.

    Where they DO exist (SW-IC-FF), a hoop lies in a horizontal plane wrapping the 4phi14 group, so
    a 2D in-plane model sees only its two legs running along x — 2 x A_phi8 spanning the boundary
    region. The two legs running through the thickness are out of plane and cannot be represented.
    """
    if not CONFINED_BOUNDARY:
        return ()
    half = BOUNDARY_RATIO * LW / 2.0
    xc = boundary_x()
    x0, x1 = xc - half, xc + half
    nrows = int(round(HW / TIE_SPACING))
    ties: list[Rebar] = []
    for i in range(nrows + 1):
        y = i * TIE_SPACING
        ties.append(Rebar([(-x1, y), (-x0, y)], TIE_AREA, S8, role="stirrup"))
        ties.append(Rebar([(x0, y), (x1, y)], TIE_AREA, S8, role="stirrup"))
    return tuple(ties)


def head_rebars() -> tuple[Rebar, ...]:
    """phi16@100 horizontal reinforcement in the top 500 mm loading head (Fig. 3c)."""
    n = int(round(HEAD_DEPTH / HEAD_SPACING))
    y0 = HW - HEAD_DEPTH
    return tuple(
        Rebar([(-LW / 2.0, y0 + i * HEAD_SPACING), (LW / 2.0, y0 + i * HEAD_SPACING)],
              HEAD_AREA, S16, role="stirrup")
        for i in range(n + 1)
    )


def web_rebars() -> tuple[Rebar, ...]:
    """Distributed web reinforcement: phi12@200 verticals and phi10@200 horizontals.

    Verticals are anchored into the pedestal like the boundary bars. Horizontals carry
    `role="stirrup"` so the visualiser draws them distinctly (D32); they also give the lattice the
    transverse tie path that keeps it from forming a local mechanism (D23).
    """
    verticals = tuple(Rebar(_hooked_path(dx), WEB_V_AREA, S12) for dx in WEB_V_X)
    nrows = int(round(HW / WEB_SPACING))
    horizontals = tuple(
        Rebar([(-LW / 2.0, i * WEB_SPACING), (LW / 2.0, i * WEB_SPACING)], WEB_H_AREA, S10,
              role="stirrup")
        for i in range(nrows + 1)
    )
    return verticals + horizontals


def check_mesh_alignment(mesh_size: float) -> None:
    """Fail loudly if `mesh_size` cannot place a node on every reinforcement line.

    The lattice ties rebar to SHARED nodes, so every bar coordinate must fall on the grid. The
    binding constraint here is the boundary group centroid at x = +-450 mm together with the
    200 mm web spacing, which restricts the grid to divisors of 50 mm (50, 25, 10, ...). A mesh
    that misses a bar line would otherwise produce a wall quietly missing that reinforcement.
    """
    half = BOUNDARY_RATIO * LW / 2.0
    xc = boundary_x()
    required = {xc, xc - half, xc + half, *(abs(x) for x in WEB_V_X), LW / 2.0, WEB_SPACING,
                ANCHORAGE_DEPTH, HOOK_LENGTH, HEAD_DEPTH, HEAD_SPACING, PED_L / 2.0, PED_H, HW,
                A_SHEAR, xc + HOOK_LENGTH, max(abs(x) for x in WEB_V_X) + HOOK_LENGTH}
    if TIE_SPACING:
        required.add(TIE_SPACING)
    bad = sorted(v for v in required if abs(v / mesh_size - round(v / mesh_size)) > 1e-9)
    if bad:
        raise ValueError(f"mesh_size={mesh_size:g} mm does not place nodes on {bad} — every rebar "
                         f"line and boundary must be a whole number of cells. Use a divisor of "
                         f"{int(min(required))} mm (e.g. 50 or 25).")


def rebars() -> tuple[Rebar, ...]:
    """Full reinforcement cage: boundary longitudinals (hooked into the pedestal), the web
    curtains, the loading-head bars, and — for SW-IC-FF only — the boundary confinement hoops."""
    return boundary_rebars() + boundary_ties() + web_rebars() + head_rebars()


# --- cyclic loading protocol (paper Fig. 9, per ACI 374.2R-13) ----------------------------------
# Displacement amplitudes at the actuator level, mm. The first two are 0.5*delta_y and delta_y
# (delta_y = 4.4 mm); the rest are the displacement-controlled drift levels 0.3 ... 4.0%.
# Three cycles are applied up to the yield displacement and two at every level after it.
PROTOCOL_PEAKS = (2.2, 4.4, 6.6, 11.0, 16.5, 22.0, 33.0, 44.0, 55.0, 66.0, 88.0)
PROTOCOL_CYCLES = (3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2)


def protocol(max_drift: float | None = None) -> tuple[tuple[float, ...], tuple[int, ...]]:
    """The test protocol truncated at `max_drift` (a ratio, e.g. 0.02 for 2%); None keeps it all."""
    if max_drift is None:
        return PROTOCOL_PEAKS, PROTOCOL_CYCLES
    keep = [i for i, p in enumerate(PROTOCOL_PEAKS) if p <= max_drift * A_SHEAR + EPS]
    if not keep:
        raise ValueError(f"max_drift={max_drift} is below the first protocol amplitude "
                         f"({PROTOCOL_PEAKS[0] / A_SHEAR:.4%})")
    return (tuple(PROTOCOL_PEAKS[i] for i in keep), tuple(PROTOCOL_CYCLES[i] for i in keep))


def lateral_loads(model) -> list[Load]:
    """Reference lateral pattern: a horizontal line load at the actuator level, y = 2200 mm."""
    ids = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))
    return [Load(nid, (H_REF / len(ids), 0.0)) for nid in ids]


def control_node(model) -> int:
    """The drift control node: wall centreline at the load level (drift = u / A_SHEAR)."""
    return select_nodes(model, (-EPS, EPS, A_SHEAR - EPS, A_SHEAR + EPS))[0]


def base_nodes(model) -> list[int]:
    """The supported pedestal soffit — where base shear is summed."""
    return select_nodes(model, (-PED_L, PED_L, -PED_H - EPS, -PED_H + EPS))
