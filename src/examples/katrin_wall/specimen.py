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

# The analysis scripts call the actuator height A_SHEAR, as the SW-NC-FF study does, so the two
# example packages read the same. Drift is u / A_SHEAR throughout — computed on the MODEL's shear
# span, which is 0.22% shorter than the paper's 4560, so a model drift is 0.22% (relative) higher
# than the paper would call the same displacement. Below the digitization error, but stated.
A_SHEAR = L_V_MODEL

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


# Table 2's A_gt, the strain at MAXIMUM force, kept per grade as the grades are built so the
# number that sets hardening and the number that breaks a bar cannot drift apart (D93).
#
# A_gt IS NOT THE FRACTURE STRAIN. It is uniform elongation, where necking starts; a bar fractures
# somewhat past it. Using A_gt as a rupture strain therefore removes the bar at its peak force
# rather than at its true fracture, which is CONSERVATIVE (early) — and it is measured, which is
# the point: on this specimen the rupture strain is data, not the assumption it has to be on the
# Aldemir wall (D91).
AGT: dict[str, float] = {}


def _steel(name: str, rp02: float, rm: float, agt: float) -> SteelGrade:
    b = (rm - rp02) / (ES * (agt - rp02 / ES))
    AGT[name] = agt
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
# TENSILE STRENGTH — evaluated at the CHARACTERISTIC strength, not the mean, and the paper is what
# settles it (D65/D66). It never prints an f_t, but Table 5's cracking moment pins the one it used:
# M_cr = (f_ctm + N/A_g)*t*l_w^2/6 with N/A_g = 2.287 MPa needs f_ctm = 2.97 MPa, which is EC2's
# 0.30*f_ck^(2/3) at f_ck = f'_c - 8 = 31.2. The same reading reproduces WSH1/WSH5/WSH6's M_cr to
# 1%, so it is the paper's convention rather than a coincidence.
#
# Applying the same formula to the MEAN f'_c = 39.2 instead gives 3.46 MPa, and that 17% error is
# not cosmetic: the regularized softening slope goes as Ets = f_t^2*L/(2*Gf), so it costs 36% in
# softening steepness. At 3.46 a strut fell from cracking to zero stress in 10.1x its cracking
# strain, against 16.2x for the SW-NC-FF wall that ran to completion, and the first cyclic attempt
# lost convergence at 0.065% drift (D66).
#
# NOTE the two MC90 quantities take DIFFERENT strengths on purpose: f_ctm is a characteristic-based
# formula, while G_f = G_f0*(f_cm/10)^0.7 is defined on the MEAN, which is f'_c here.
FCK = FC - 8.0                         # EC2 / MC90 characteristic from mean
FT = 0.30 * FCK ** (2.0 / 3.0)         # 2.97 MPa — matches the paper's own M_cr
GF_C = 0.030 * (FC / 10.0) ** 0.7      # MC90 fracture energy on the MEAN strength, 16 mm aggregate

# epsc0 is DERIVED, not quoted: Concrete02's initial compressive tangent is 2*fc/epsc0 whatever the
# grade's E says, so any other value silently makes the compressive modulus disagree with the
# measured Ec (D56). Here 2*39.2/35200 = 0.00223, within 1.2% of the 0.0022 a table would give —
# but deriving it costs nothing and removes the inconsistency.
EPSC0 = 2.0 * FC / EC

WALL = ConcreteGrade("wall", E=EC, nu=NU, rho=UNIT["rho_c"],
                     fc=FC, epsc0=EPSC0, fcu=0.2 * FC, epsU=0.010, ft=FT)
# The foundation is 700 thick against the wall's 150; a plane model carries one thickness, so its
# extra width is folded into an equivalent modulus. It is a restraint device, kept elastic.
FOUNDATION = ConcreteGrade("foundation", E=EC * (FND_W / TW), nu=NU, rho=UNIT["rho_c"],
                           fc=FC, epsc0=EPSC0, fcu=0.2 * FC, epsU=0.010, ft=FT)
CAP = ConcreteGrade("cap", E=EC * CAP_STIFFNESS, nu=NU, rho=UNIT["rho_c"],
                    fc=FC, epsc0=EPSC0, fcu=0.2 * FC, epsU=0.010, ft=FT)
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


# --- cyclic loading protocol (Fig. 6 + Sec. 2.3) ------------------------------------------------
# Two FORCE-controlled cycles to 0.75*F_y, then two DISPLACEMENT-controlled cycles at each
# ductility level mu = 2..6, amplitude mu * delta_y. Every cycle goes South (+) first.
#
# delta_y here is the 3/4-RULE value, 15.4 mm — deliberately, even though the paper later judged it
# the poorer of its three yield estimates (Sec. 2.3, Table 4). It is the one the actuator history
# was actually built on, so it is the one that reproduces the amplitudes the wall really saw. The
# better estimate, 16.2 mm, belongs in a ductility calculation, not in the drive.
DELTA_Y = 15.4
PROTOCOL_PEAKS = (11.4, 30.8, 46.2, 61.6, 77.0, 92.4)   # 0.75Fy cycle, then mu = 2, 3, 4, 5, 6
PROTOCOL_CYCLES = (2, 2, 2, 2, 2, 2)


def protocol(max_drift: float | None = None) -> tuple[tuple[float, ...], tuple[int, ...]]:
    """The test protocol truncated at `max_drift` (a ratio, e.g. 0.01 for 1%); None keeps it all."""
    if max_drift is None:
        return PROTOCOL_PEAKS, PROTOCOL_CYCLES
    keep = [i for i, p in enumerate(PROTOCOL_PEAKS) if p <= max_drift * A_SHEAR + EPS]
    if not keep:
        raise ValueError(f"max_drift={max_drift} is below the first protocol amplitude "
                         f"({PROTOCOL_PEAKS[0] / A_SHEAR:.4%})")
    return (tuple(PROTOCOL_PEAKS[i] for i in keep), tuple(PROTOCOL_CYCLES[i] for i in keep))


# --- how the protocol is DRIVEN (dynamic relaxation) --------------------------------------------
# Shared by pushover and cyclic so the two stay comparable, exactly as in the SW-NC-FF study.
#
# THE REAL TEST CANNOT BE REPRODUCED AT ITS OWN SPEED, and this is the sharpest analysis difference
# between the two walls. Fig. 6 gives the actuator at **1.2 to 3.6 mm per MINUTE** — 0.02-0.06 mm/s,
# against SW-NC-FF's published 7.6 mm/s. Driving this protocol at 0.06 mm/s would take ~43,000 s of
# simulated time at a ~1e-3 s step: some 4x10^7 steps, months of wall-clock. So the drive here is
# ~130x faster than the experiment and the licence for that is NOT the rate. It is the MEASURED
# residual: `--quasi-static` reports the inertia + damping the drive and base reactions fail to
# balance, and that is the number to read.
#
# 7.6 mm/s is kept — the same ABSOLUTE speed as the SW-NC-FF runs rather than the same drift rate.
# Two reasons. The contaminating forces (C.v, and the ringing released at each crack) scale with the
# speed, not with the drift; holding the speed holds them, while this wall's base shear is ~2x
# larger (454 vs 213 kN), so the same absolute contamination reads about HALF as much of the peak.
# And the total drive path is similar — 2555 mm for this 2%-drift protocol against 2816 mm for
# SW-NC-FF's 4% one — so run time stays in the same class, modulo this model being ~3x larger.
QUASI_STATIC_RATE = 7.6
#
# DAMPING is the knob that matters, not rate (D64): the contamination is crack-release ringing,
# which damping suppresses and a slower drive does not. 0.2 is both the measured optimum there and
# a defensible physical value (10-20% equivalent viscous damping for heavily cracked RC).
DAMPING_RATIO = 0.2


# --- vertical strain gauge: the test's own LVDT chain ---------------------------------------------
# Rows of nodes across the wall width, differenced pairwise: the same construction as the chain of
# vertical LVDTs on the tested wall's west face (Fig. 5b), from which the paper derives its strain
# profiles (Fig. 10), curvature profiles (Fig. 11c) and base curvature (Fig. 15d).
#
# LVDT_ROWS are the boundaries between successive devices, RECONSTRUCTED from Fig. 5b's dimension
# string — 50, then 3x100, then 4x300, then 2x700. The reconstruction is checkable: it totals
# 2950 mm, and Fig. 11's height axis runs to exactly 3.0 m. Every boundary is a whole number of
# 50 mm cells, so each row lands on nodes.
#
# The lowest device (0-50 mm) is DELIBERATELY ABSENT. The paper discards it: "due to the strain
# penetration into the foundation it is not possible to define a meaningful base length" (Sec. 4.2).
# Including it would let the model's own missing strain-penetration show up as a spurious profile
# point rather than as the known omission it is.
LVDT_ROWS = (50.0, 150.0, 250.0, 350.0, 650.0, 950.0, 1250.0, 1550.0, 2250.0, 2950.0)

# The BASE profile — the one quantity the SW-NC-FF study also has — is rows[0] to rows[3], i.e.
# y = 50 to 350 mm. That brackets the h = 60-360 mm window the paper's own Demec tensile strains
# are averaged over (Table 6 note), and starts clear of the contaminated base row.
GAUGE_Y0, GAUGE_H = LVDT_ROWS[0], LVDT_ROWS[3] - LVDT_ROWS[0]
BASE_SEGMENT = (0, 3)

# Yield curvature from the paper's section analysis, ~the same for all six units (Sec. 5.1). It is
# what defines the plastic-zone height L_pz — the height at which curvature falls to phi_y — and so
# which segments the base-curvature fit is allowed to use.
PHI_YIELD = 3.0e-6                      # 3 /km, in 1/mm


def gauge_nodes(model, *, rows: tuple[float, ...] = LVDT_ROWS) -> tuple[list[float], list[list[int]]]:
    """Aligned node columns at each row height: `(xs, ids_per_row)`.

    `ids_per_row[k][i]` sits at x = `xs[i]`, y = `rows[k]`, so differencing any two rows column by
    column gives the average vertical strain of that column over that segment.

    Restricted to the WALL width, and every row is checked for alignment — a row that landed off the
    grid would silently return a different node count and corrupt every profile after it.
    """
    half = LW / 2.0
    ids, xs0 = [], None
    for y in rows:
        got = select_nodes(model, (-half - EPS, half + EPS, y - EPS, y + EPS))
        x = [model.nodes[n].coords[0] for n in got]
        if xs0 is None:
            xs0 = x
        elif len(x) != len(xs0) or max(abs(a - b) for a, b in zip(x, xs0)) > 1e-6:
            raise ValueError(f"gauge row y={y:g} has {len(x)} nodes against {len(xs0)} at "
                             f"y={rows[0]:g}. Every row height must be a whole number of cells.")
        ids.append(got)
    return xs0, ids


def gauge_probe(ids_per_row: list[list[int]]) -> list[int]:
    """Flatten the rows into the single node list a runner's `node_history` takes."""
    return [n for row in ids_per_row for n in row]


def segment_strain(sample, n_x: int, lo: int, hi: int,
                   rows: tuple[float, ...] = LVDT_ROWS) -> list[float]:
    """Average vertical strain per column over the segment between rows `lo` and `hi`.

    `sample` is one `node_history` record: the rows CONCATENATED in `LVDT_ROWS` order, so each row
    is a fixed slice of it. Positive = tension (the column got longer).
    """
    a = sample[lo * n_x:(lo + 1) * n_x]
    b = sample[hi * n_x:(hi + 1) * n_x]
    length = rows[hi] - rows[lo]
    return [(t - s) / length for s, t in zip(a, b)]


def gauge_strains(node_history, xs, segment: tuple[int, int] = BASE_SEGMENT) -> list[list[float]]:
    """Strain profile across the width for one segment, for every sample in a `node_history`."""
    n = len(xs)
    return [segment_strain(v, n, *segment) for v in node_history["values"]]


def curvature_profile(node_history, xs, *, rows: tuple[float, ...] = LVDT_ROWS):
    """`(mid_heights, [[curvature per segment] per sample])` — the model's Fig. 11c.

    One curvature per LVDT segment, obtained the way the test obtains it: least-squares slope of the
    vertical-strain profile across the width, reported at the segment's mid-height. Units 1/mm.
    """
    n = len(xs)
    mids = [0.5 * (rows[k] + rows[k + 1]) for k in range(len(rows) - 1)]
    out = []
    for vals in node_history["values"]:
        out.append([neutral_axis(xs, segment_strain(vals, n, k, k + 1, rows))[0]
                    for k in range(len(rows) - 1)])
    return mids, out


def base_curvature(mids, curvature, *, phi_y: float = PHI_YIELD) -> tuple[float, float]:
    """`(phi_base, L_pz)` by the paper's own construction (Sec. 5.1) — its Fig. 15d quantity.

    Not a reading at the base, which strain penetration makes meaningless. Instead: a best-fit
    LINEAR curvature profile over the height of the plastic zone, extrapolated down to y = 0, where
    L_pz is the height at which the curvature has fallen to the yield curvature phi_y.

    Reproducing the construction is the point — comparing the model's base-segment curvature
    against the paper's extrapolated value would be comparing two different quantities.
    """
    used = [(h, c) for h, c in zip(mids, curvature) if abs(c) >= phi_y]
    if len(used) < 2:
        return float("nan"), 0.0
    lpz = max(h for h, _c in used)
    hs = [h for h, _c in used]
    cs = [c for _h, c in used]
    mh, mc = sum(hs) / len(hs), sum(cs) / len(cs)
    shh = sum((h - mh) ** 2 for h in hs)
    shc = sum((h - mh) * (c - mc) for h, c in zip(hs, cs))
    slope = shc / shh if shh else 0.0
    return mc - slope * mh, lpz


def neutral_axis(xs, strain) -> tuple[float, float, float]:
    """Least-squares line through one strain profile: `(curvature, strain_at_centre, x_zero)`.

    A plane section would make the profile exactly linear, so the fit is both a summary (curvature
    = the slope, in 1/mm) and a test of that assumption — which matters more here than it did for
    SW-NC-FF, because this paper measures curvature directly (Fig. 11c, Fig. 15d) and warns in
    Sec. 5.1 that inclined flexure-shear cracks break plane sections near the base.

    `x_zero` is where the fitted line crosses zero — the neutral-axis position — and is `nan` for a
    profile with no gradient (a purely axial state, before any cracking localizes).
    """
    n = len(xs)
    mx = sum(xs) / n
    me = sum(strain) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxe = sum((x - mx) * (e - me) for x, e in zip(xs, strain))
    slope = sxe / sxx if sxx else 0.0
    intercept = me - slope * mx
    x0 = -intercept / slope if abs(slope) > 1e-18 else float("nan")
    return slope, intercept + slope * 0.0, x0
