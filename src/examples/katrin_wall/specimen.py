"""Shared specimen definition for the ETH Zurich WSH wall series.

Source: Dazio, Beyer & Bachmann (2009), "Quasi-static cyclic tests and plastic hinge analysis of
RC structural walls", Engineering Structures 31:1556-1571. Six half-scale slender cantilever walls
(WSH1-WSH6) with identical geometry, differing in longitudinal reinforcement content, steel
ductility, axial load and confinement.

WHY THIS SERIES IS A GOOD TARGET FOR A PERFECT-BOND LATTICE. The Sahinkaya SW-NC-FF wall is
dominated by ROCKING from plain-bar debonding (74% of displacement at 2% drift), which a
shared-node lattice cannot represent — that capped what the model could claim. The WSH walls use
DEFORMED bars that stay bonded; their displacement decomposes into flexure and shear with only a
small fixed-end (strain-penetration) component (paper Fig. 9a). Perfect bond is therefore a far
better assumption here, and drift capacity becomes a fair comparison rather than an unfair one.

NOT a shear-controlled series: all six failed in flexure, with shear only 5-12% of the flexural
displacement (paper Fig. 9b) at a shear span ratio of 2.28. For genuinely shear-controlled walls
see Terzioglu, Orakcal & Massone (2018) on squat walls.

Units: N, mm (stresses MPa, forces N).
"""

from __future__ import annotations

import math
from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, CompoundRectangles, ConcreteGrade, Problem, Rebar, SteelGrade

SPECIMEN = "WSH3"          # the paper's "model" ductile wall; see WSH_UNITS below

# --- geometry, identical for all six units (paper Sec. 2.1, Fig. 2) -----------------------------
LW, TW = 2000.0, 150.0     # wall length (in-plane), thickness (out-of-plane)
H_WALL = 4030.0            # wall height above the foundation
L_V = 4560.0               # shear span to the actuator (WSH6: 4520) -> L_V/LW = 2.28
HEAD_L, HEAD_H, TAPER_H = 400.0, 700.0, 220.0   # tapered loading head above the wall
FND_L, FND_H, FND_W = 2800.0, 600.0, 700.0      # foundation block: long x deep x thick

MESH, HORIZON = 50.0, 1.5
EPS = 1e-6
OUT = Path(__file__).resolve().parent.parent / "output" / "katrin_wall"

# The 400-wide head sits on a 220 taper above the 2000-wide wall — geometry a rectangle mesher
# cannot express. It is replaced by a STIFF ELASTIC CAP of full wall width spanning from the wall
# top to the actuator level, which preserves the shear span exactly (the quantity that sets the
# moment-to-shear ratio at every section) and keeps the load application rigid, as the real head
# does. The cap sits far above the plastic zone, where the wall is elastic, so the idealization is
# second-order. Heights are rounded to whole mesh cells; the error is reported by
# `check_mesh_alignment`.
H_WALL_MODEL = 4050.0      # +0.5% on the real 4030
L_V_MODEL = 4550.0         # -0.2% on the real 4560
CAP_STIFFNESS = 20.0       # cap modulus multiplier — rigid relative to the wall

# --- in-plane bar positions, read from Fig. 1c (WSH3) -------------------------------------------
# The dimension string across the section is
#   30|100|100|125|125|125|125|125|145|145|125|125|125|125|125|100|100|30 = 2000,
# giving 17 in-plane positions x 2 curtains = 34 bars = 6phi12 + 22phi8 + 6phi12.
_SPACINGS = (30, 100, 100, 125, 125, 125, 125, 125, 145, 145, 125, 125, 125, 125, 125, 100, 100, 30)


def _positions() -> tuple[list[float], list[float]]:
    """(boundary, web) in-plane bar x-coordinates, centred on the wall (x = 0 at mid-length)."""
    xs, acc = [], 0.0
    for s in _SPACINGS[:-1]:
        acc += s
        xs.append(acc - LW / 2.0)
    if len(xs) != 17:
        raise RuntimeError(f"expected 17 bar positions, got {len(xs)}")
    return [*xs[:3], *xs[-3:]], xs[3:-3]


BOUNDARY_X, WEB_X = _positions()
H_SPACING = 150.0          # phi6@150 horizontal web reinforcement, both faces -> rho_h = 0.25%
PLASTIC_ZONE = 1700.0      # height of the close-tied plastic zone (Fig. 2)
ANCHORAGE = 450.0          # bars run into the 600-deep foundation, stopping clear of the soffit


def _area(d: float, n: int) -> float:
    return n * math.pi * d * d / 4.0


# --- per-unit properties (Tables 1-3); bar POSITIONS are read from Fig. 1c and are WSH3's -------
# Adding another unit requires reading its own panel of Fig. 1: the spacing strings differ.
WSH_UNITS = {
    "WSH3": {
        "N": 686.0e3,                       # axial load, N   (N/Ag*fc' = 0.058)
        "fc": 39.2, "Ec": 35200.0, "rho_c": 2.381e-9,
        "boundary": (12.0, 6), "web": (8.0, 22), "horizontal": (6.0, 2),
        "tie": (6.0, 2), "tie_spacing": 75.0,
        "V_max": 454.0e3, "drift_u": 0.0203, "mu": 5.7,
    },
}
UNIT = WSH_UNITS[SPECIMEN]

# --- steel (Table 2). These bars have NO pronounced yield plateau, so Rp02 is used as fy — which
# suits Steel02's smooth Giuffre-Menegotto-Pinto curve. Hardening b is measured, not assumed:
# b = (Rm - Rp02) / (Es * (Agt - Rp02/Es)). Unlike the Sahinkaya plain bars (b ~ 0.002-0.004),
# these deformed bars harden genuinely.
ES = 200000.0


def _steel(name: str, rp02: float, rm: float, agt: float) -> SteelGrade:
    b = (rm - rp02) / (ES * (agt - rp02 / ES))
    return SteelGrade(name, fy=rp02, E0=ES, b=b)


S12 = _steel("phi12", 601.0, 725.5, 0.0769)   # boundary longitudinal
S8 = _steel("phi8", 569.2, 700.2, 0.0734)     # web longitudinal
S6 = _steel("phi6", 489.0, 552.2, 0.0645)     # horizontal web + confinement hoops
S42 = _steel("phi4.2", 562.2, 615.0, 0.0306)  # crossties

BOUNDARY_AREA = _area(*UNIT["boundary"]) / 3.0   # 6phi12 spread over 3 in-plane positions
WEB_AREA = _area(*UNIT["web"]) / len(WEB_X)      # 22phi8 spread over the web positions
H_AREA = _area(*UNIT["horizontal"])              # phi6, two curtains, per horizontal line
TIE_AREA = _area(*UNIT["tie"])                   # boundary hoop: 2 in-plane legs of phi6

# --- concrete grades ----------------------------------------------------------------------------
NU = 0.20
FC, EC = UNIT["fc"], UNIT["Ec"]
FT = 0.30 * FC ** (2.0 / 3.0)          # CEB-FIP MC90 mean tensile strength
GF_C = 0.030 * (FC / 10.0) ** 0.7      # MC90 fracture energy, 16 mm aggregate

WALL = ConcreteGrade("wall", E=EC, nu=NU, rho=UNIT["rho_c"],
                     fc=FC, epsc0=0.0022, fcu=0.2 * FC, epsU=0.010, ft=FT)
# The foundation is 700 thick against the wall's 150; a plane model carries one thickness, so its
# extra width is folded into an equivalent modulus. It is a restraint device, kept elastic.
FOUNDATION = ConcreteGrade("foundation", E=EC * (FND_W / TW), nu=NU, rho=UNIT["rho_c"],
                           fc=FC, epsc0=0.0022, fcu=0.2 * FC, epsU=0.010, ft=FT)
CAP = ConcreteGrade("cap", E=EC * CAP_STIFFNESS, nu=NU, rho=UNIT["rho_c"],
                    fc=FC, epsc0=0.0022, fcu=0.2 * FC, epsU=0.010, ft=FT)
GRADES = {"wall": WALL, "foundation": FOUNDATION, "cap": CAP}
GF = {"wall": GF_C, "foundation": GF_C, "cap": GF_C}


def wall_problem() -> Problem:
    """Wall + stiff loading cap on the foundation block, clamped at the foundation soffit."""
    domain = CompoundRectangles(
        rects=[
            (-FND_L / 2.0, -FND_H, FND_L, FND_H),                    # foundation
            (-LW / 2.0, 0.0, LW, L_V_MODEL),                         # wall + cap (one rectangle)
        ],
        thickness=TW,
    )
    supports = [BoxSupport(box=(-FND_L, FND_L, -FND_H - EPS, -FND_H + EPS), fix=(1, 1))]
    loads = [BoxLoad(box=(-LW / 2.0, LW / 2.0, L_V_MODEL - EPS, L_V_MODEL + EPS),
                     total=(0.0, -UNIT["N"]))]
    return Problem(ndm=2, ndf=2, domain=domain, material=WALL, supports=supports, loads=loads)


def zone_of(x: float, y: float) -> str:
    if y < 0.0:
        return "foundation"
    return "cap" if y >= H_WALL_MODEL else "wall"


def snap(value: float, mesh_size: float) -> float:
    return round(value / mesh_size) * mesh_size


def rebars(mesh_size: float = MESH) -> tuple[Rebar, ...]:
    """Longitudinal bars (boundary + web) and the horizontal web curtains.

    Bar positions are SNAPPED to the grid: the true spacings (100/125/145 mm) share no common
    divisor coarse enough to mesh, so exact placement would need a 10 mm grid. `check_mesh_alignment`
    reports the resulting error in the boundary-group centroid, which is what sets the flexural
    lever arm.
    """
    top = H_WALL_MODEL
    bars = [Rebar([(snap(x, mesh_size), -ANCHORAGE), (snap(x, mesh_size), top)],
                  BOUNDARY_AREA, S12) for x in BOUNDARY_X]
    bars += [Rebar([(snap(x, mesh_size), -ANCHORAGE), (snap(x, mesh_size), top)],
                   WEB_AREA, S8) for x in WEB_X]
    n = int(round(H_WALL_MODEL / H_SPACING))
    bars += [Rebar([(-LW / 2.0, snap(i * H_SPACING, mesh_size)),
                    (LW / 2.0, snap(i * H_SPACING, mesh_size))], H_AREA, S6, role="stirrup")
             for i in range(n + 1)]
    bars += list(boundary_ties(mesh_size))
    return tuple(bars)


def tie_spacing(mesh_size: float) -> tuple[float, float]:
    """(spacing used, area scale) for the boundary hoops on this grid.

    The true 75 mm spacing is not a whole number of 50 mm cells, so it is snapped — but the AREA is
    then rescaled by the same factor, holding the transverse steel per unit height (and hence the
    confinement ratio) exactly at the specimen's value. Preserving the mechanical quantity matters
    more than preserving the geometric spacing. A 25 mm grid divides 75 exactly, so the scale is 1.
    """
    s = max(mesh_size, round(UNIT["tie_spacing"] / mesh_size) * mesh_size)
    return s, s / UNIT["tie_spacing"]


def boundary_ties(mesh_size: float = MESH) -> tuple[Rebar, ...]:
    """Boundary confinement hoops (phi6@75) over the plastic zone, projected onto the in-plane model.

    A hoop lies in a HORIZONTAL plane wrapping the boundary bars, so a 2D in-plane model sees only
    its two legs running along the wall length; the two legs through the thickness, and the phi4.2
    crossties (which run through the thickness), are out of plane and cannot be represented here.
    The legs span between the outermost and innermost boundary bar, which is what the hoop ties
    together.
    """
    spacing, scale = tie_spacing(mesh_size)
    x_out = snap(abs(BOUNDARY_X[0]), mesh_size)
    x_in = snap(abs(BOUNDARY_X[2]), mesh_size)
    n = int(round(PLASTIC_ZONE / spacing))
    ties: list[Rebar] = []
    for i in range(n + 1):
        y = snap(i * spacing, mesh_size)
        ties.append(Rebar([(-x_out, y), (-x_in, y)], TIE_AREA * scale, S6, role="stirrup"))
        ties.append(Rebar([(x_in, y), (x_out, y)], TIE_AREA * scale, S6, role="stirrup"))
    return tuple(ties)


def check_mesh_alignment(mesh_size: float) -> None:
    """Fail on a grid that cannot host the bars, and report the snapping error."""
    for name, value in (("wall+cap height", L_V_MODEL), ("wall height", H_WALL_MODEL),
                        ("foundation depth", FND_H), ("foundation length", FND_L),
                        ("wall length", LW), ("horizontal spacing", H_SPACING)):
        if abs(value / mesh_size - round(value / mesh_size)) > 1e-9:
            raise ValueError(f"mesh_size={mesh_size:g} does not divide {name} ({value:g} mm)")
    true_c = sum(abs(x) for x in BOUNDARY_X[:3]) / 3.0
    snap_c = sum(abs(snap(x, mesh_size)) for x in BOUNDARY_X[:3]) / 3.0
    print(f"bar snapping at mesh={mesh_size:g} mm: boundary centroid {true_c:.1f} -> {snap_c:.1f} mm "
          f"({abs(snap_c - true_c) / true_c * 100:.2f}% on the flexural lever arm)")


def lateral_loads(model) -> list[Load]:
    ids = select_nodes(model, (-LW, LW, L_V_MODEL - EPS, L_V_MODEL + EPS))
    return [Load(nid, (UNIT["V_max"] / len(ids), 0.0)) for nid in ids]


def control_node(model) -> int:
    return select_nodes(model, (-EPS, EPS, L_V_MODEL - EPS, L_V_MODEL + EPS))[0]


def base_nodes(model) -> list[int]:
    return select_nodes(model, (-FND_L, FND_L, -FND_H - EPS, -FND_H + EPS))
