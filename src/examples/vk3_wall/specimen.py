"""Shared specimen definition for ETH Zurich wall-type bridge pier Test Unit VK3.

Source: Bimschas, M. (2010), "Displacement Based Seismic Assessment of Existing Bridges in Regions
of Moderate Seismicity", IBK Bericht Nr. 326, Institut fuer Baustatik und Konstruktion, ETH Zuerich;
vdf Hochschulverlag. doi:10.3929/ethz-a-006237119. Chapter 5, "Experimental Test Campaign".

Three large-scale (1:2) quasi-static cyclic tests on squat wall-type bridge piers, VK1/VK2/VK3,
differing
ONLY in their longitudinal reinforcement: VK1 28phi14 continuous, VK2 the same with a 600 mm
lap-splice at the base, VK3 42phi14 continuous.

WHY VK3 IS A DIFFERENT TARGET FROM THE OTHER TWO WALLS IN THIS REPO. SW-NC-FF was rocking-dominated
(plain bars, capped the scope) and WSH3 was flexure-dominated with shear at ~12% of displacement.
VK3 is the first SHEAR-relevant specimen: its aspect ratio is essentially WSH3's (2.20 vs 2.28) but
it carries THREE TIMES LESS transverse steel (rho_sw = 0.08% vs 0.25%), shear reaches 20-22% of the
top displacement (Fig. 5.19 right), and the test ended in a combined SHEAR + AXIAL-LOAD failure at
mu_prov = 5 (52.5 mm, 1.59% drift). A lattice represents diagonal action explicitly, so this is the
mechanism it is best placed to say something new about — and also where its calibration is weakest
(the Aydin energy balance matches C11 exactly but never matches shear stiffness independently; see
`build.report_calibration` and D53).

The paper's own reading of the failure (Sec. 5.4) is lattice-shaped: the shear failure was NOT
diagonal-tension or web crushing, but the loss of the compression zone that had been SUPPORTING the
diagonal strut, forcing a switch to a hoop-based mechanism the 0.08% steel could not deliver. That
is a load-path collapse, which `element_groups` can measure directly.

BOND IS FAIR HERE. VK3's longitudinal bars are continuous deformed bars (topar-S 500C, ductility
class C) with NO lap-splice, so perfect bond is a reasonable idealization — as it was for WSH3 and
was NOT for the plain-bar SW-NC-FF wall.

Units: N, mm (stresses MPa, forces N).
"""

from __future__ import annotations

import math
import sys
from dataclasses import replace
from datetime import datetime
from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, CompoundRectangles, ConcreteGrade, Problem, Rebar, SteelGrade

SPECIMEN = "VK3"

# --- geometry, identical for all three units (Tab. 5.1, Fig. 5.1) -------------------------------
LW, TW = 1500.0, 350.0     # section depth (in-plane), section width (out-of-plane)
L_V = 3300.0               # shear span: horizontal actuator axis above the pier base -> Lv/lw = 2.20
H_PIER = 3700.0            # full pier height above the foundation (400 mm continues above L_V)
FND_L, FND_H = 3000.0, 900.0
COVER_L = 26.0             # clear cover to the longitudinal bars (Tab. 5.1)

# The foundation block's OUT-OF-PLANE thickness is never given in the chapter. 1000 mm is assumed
# and folded into an equivalent modulus, exactly as WSH3's 700 mm block is. The block is a restraint
# device held elastic, and `elastic.py --foundation-sensitivity` prices the assumption: the global
# stiffness barely moves across 700-1200 mm.
FND_W = 1000.0

MESH, HORIZON = 50.0, 1.5
EPS = 1e-6
OUT = Path(__file__).resolve().parent.parent / "output" / "vk3_wall"
RUNS = OUT / "runs"


def run_dir(kind: str, *, tag: str | None = None) -> Path:
    """A fresh TIMESTAMPED directory for one run. Nothing is ever overwritten.

    The other two wall packages write every run to a FIXED stem, so a two-minute diagnostic
    silently destroys a twenty-hour result unless someone remembers to archive it first — which is
    exactly the trap `katrin_wall/runs/` had to be created after the fact to close. Here each run
    gets its own directory from the start, and the archiving discipline is structural rather than
    a rule someone has to follow.

    Also drops a `command.txt` next to the results, because a figure whose command line is lost is
    a figure that cannot be reproduced.
    """
    name = f"{datetime.now():%Y-%m-%d_%H%M%S}_{kind}" + (f"_{tag}" if tag else "")
    d = RUNS / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "command.txt").write_text(" ".join(sys.argv) + "\n")
    link = RUNS / f"latest_{kind}"
    try:                                   # convenience pointer for replot/compare
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(d.name)
    except OSError:
        pass                               # a filesystem without symlinks is not a reason to fail
    return d


def find_run(kind: str | None = None, *, name: str | None = None, complete: bool = True) -> Path:
    """Locate a run directory: an explicit `name`, else the NEWEST one of `kind`.

    `complete=True` skips directories with no `data.json`. That matters because a run that is still
    in flight — or one that was cancelled — has already created its directory and written its
    `command.txt`, so a bare "newest of this kind" would hand `replot.py` a run with nothing in it.
    Results are written only when the analysis finishes.
    """
    if name:
        d = RUNS / name
        if not d.is_dir():
            raise SystemExit(f"no run directory {d}")
        return d
    got = sorted((p for p in RUNS.glob("*")
                  if p.is_dir() and not p.name.startswith("latest_")
                  and (kind is None or kind in p.name)), key=lambda p: p.name)
    done = [p for p in got if (p / "data.json").exists()]
    if complete and done:
        return done[-1]
    if not got:
        raise SystemExit(f"no runs under {RUNS}"
                         + (f" matching {kind!r}" if kind else "")
                         + " — run an analysis first")
    if complete and not done:
        raise SystemExit(f"no COMPLETED {kind or 'run'} under {RUNS} — "
                         f"{len(got)} directory(ies) exist but none has a data.json "
                         f"(still running, or cancelled). Newest: {got[-1].name}")
    return got[-1]

# NO GEOMETRIC ROUNDING IS NEEDED, unlike WSH3. Every model dimension — 1500, 3300, 3700, 3000, 900
# — is a whole number of 50 mm cells, so the shear span is EXACT and drift needs no correction.
A_SHEAR = L_V             # drift denominator, the paper's own definition delta = Delta_top / L_v

# The top 1000 mm is the LOAD-INTRODUCTION zone: Fig. 5.1 details it with phi6@75 hoops against the
# phi6@200 used everywhere else, i.e. it is detailed to stay elastic while the base hinges. It is
# held elastic in the model for the same reason WSH3's loading head is, and it is where the
# horizontal actuator (y = 3300) and the vertical actuators (y = 3700) both attach. It sits well
# above the diagonal cracking, which in the test reached ~2000 mm (Fig. 5.13).
HEAD_Y0 = 2700.0          # tight-hoop zone starts at 2675 in Fig. 5.1; rounded to a cell boundary

# --- in-plane bar positions, read from Fig. 5.2 (right, VK3) -------------------------------------
# 2x17phi14 at s = 80 over a 1280 mm span, plus 2x4phi14 at the section ends. The four end bars are
# stacked THROUGH the 284 mm thickness (350 - 2*20 cover - 2*6 hoop - 14), so in a plane model they
# collapse onto a single in-plane position at the extreme fibre — the same collapse WSH3's two
# curtains undergo. The reading is confirmed by three of the paper's own printed ratios: rho_sl,
# rho_sw and the axial-load ratio all reproduce to <1% (see `summary.py`).
N_MID, S_MID = 17, 80.0
D_LONG, D_HOOP = 14.0, 6.0


def _area(d: float, n: int) -> float:
    return n * math.pi * d * d / 4.0


def _positions() -> tuple[list[float], float]:
    """(web bar x-positions, end bar |x|), centred on the section (x = 0 at mid-depth)."""
    span = (N_MID - 1) * S_MID                              # 1280 mm
    xs = [-span / 2.0 + i * S_MID for i in range(N_MID)]
    x_end = LW / 2.0 - (COVER_L + D_LONG / 2.0)             # 717.0 mm
    return xs, x_end


WEB_X, END_X = _positions()
MID_AREA = _area(D_LONG, 2)      # 2phi14 per in-plane position (one per curtain)
END_AREA = _area(D_LONG, 4)      # 4phi14 stacked through the thickness at each end
HOOP_AREA = _area(D_HOOP, 2)     # the hoop's TWO in-plane legs; see `hoops()`

# Hoop layout over the height (Fig. 5.1). TWO spacings, and the loose one is the one that matters:
# the plastic hinge and the whole diagonal-crack region sit in the s = 200 zone.
HOOP_BASE = (75.0, 2475.0, 200.0)      # (y_start, y_end, spacing) — 13phi6@200, rho_sw = 0.08%
HOOP_HEAD = (2675.0, 3630.0, 75.0)     # 14phi6@75, the load-introduction zone
ANCHORAGE = 750.0                      # bars run to 150 mm clear of the 900 mm foundation soffit

# --- materials -----------------------------------------------------------------------------------
# The chapter gives f_c,cyl = 34 MPa for VK3 and NOTHING ELSE about the concrete. E_c, f_t and G_f
# are therefore CONVENTIONS, not measurements, and are flagged as such in `summary.py`.
FC = 34.0                              # Tab. 5.4, 150x300 cylinders at the day of testing
NU = 0.20

# E_c: SIA 262's 10000*f_cm^(1/3). The specimen is Swiss and the chapter's whole frame of reference
# is SIA 262; EC2/MC90's 22000*(f_cm/10)^0.3 lands within 2% of it, so the choice is robust. (ACI's
# 4700*sqrt(f_c) gives 27.4 GPa, 15% lower — rejected as the wrong code frame, not as wrong.)
EC = 10000.0 * FC ** (1.0 / 3.0)       # 32,396 MPa

# f_t: EC2's 0.30*f_ck^(2/3) at the CHARACTERISTIC strength f_ck = f_cm - 8 = 26 MPa. This is the
# convention D67 settled on for WSH3, where the paper's own M_cr pinned it. HERE THERE IS NO SUCH
# PIN — the chapter never prints a cracking moment — so this is consistency with the repo, not
# evidence. It matters: crack-band regularization makes the softening slope go as f_t^2, and the
# resulting strut life eps_ult/eps_cr = 13.9 sits between SW-NC-FF's 19.8 (ran clean) and WSH3's
# 13.4 (needed --gf-factor 2). Budget tension stiffening from the start.
FCK = FC - 8.0
FT = 0.30 * FCK ** (2.0 / 3.0)         # 2.633 MPa
GF_C = 0.030 * (FC / 10.0) ** 0.7      # MC90 on the MEAN strength, 16 mm aggregate assumed

# epsc0 DERIVED, not quoted: Concrete02's initial compressive tangent is 2*fc/epsc0 whatever the
# grade's E says, so any other value silently mismatches the calibrated modulus (D56).
EPSC0 = 2.0 * FC / EC                  # 0.00210
RHO_C = 2.4e-9                         # t/mm^3; wall self-weight ~47 kN of the 70 kN in Sec. 5.2.4

# --- steel (Tab. 5.4), hardening b measured the same way WSH3's is ------------------------------
# NOTE a source inconsistency: Sec. 5.2.1c's text says f_s,l = 520 MPa while Tab. 5.4 says 515. The
# table is used (b changes by 0.0002 — immaterial), and `summary.py` records the discrepancy.
#
# The longitudinal steel has a PRONOUNCED YIELD PLATEAU to about eps_s = 0.025 (Sec. 5.2.3a), after
# which it hardens to 630 MPa at 12.6%. Steel02's smooth Giuffre-Menegotto-Pinto curve smears the
# plateau and the hardening into one slope; the paper itself neglected hardening entirely for its
# moment-curvature work, calling that justified.
ES = 200000.0


def _steel(name: str, fy: float, fu: float, agt: float) -> SteelGrade:
    b = (fu - fy) / (ES * (agt - fy / ES))
    return SteelGrade(name, fy=fy, E0=ES, b=b)


S_LONG = _steel("phi14", 515.0, 630.0, 0.126)    # longitudinal, topar-S 500C  -> b = 0.0047
S_HOOP = _steel("phi6", 518.0, 681.0, 0.084)     # hoops, cold-worked, no plateau -> b = 0.0100


def set_steel_R0(r0: float) -> None:
    """Override Steel02's yield-corner sharpness on both steel grades. Call BEFORE building.

    THE DISTINCTION THAT MAKES THIS LEGITIMATE. `fy` and `b` are MEASURED — b = 0.00466 follows
    from Tab. 5.4's f_y, f_u and eps_su — and must not be tuned. `R0` is different: it controls only
    how abruptly the Giuffre-Menegotto-Pinto curve turns the elastic-plastic corner, has no measured
    counterpart, and OpenSees itself recommends anywhere in 10-20. Lowering it spreads the
    transition over more strain WITHOUT moving the yield stress or the hardening slope.

    Why it is worth a try here: the VK3 model breaks up at 0.310-0.339% drift, and the chapter's own
    two estimates of first yield bracket that (0.306% predicted, 0.324% measured). With the concrete
    already cracked, the tension chord's stiffness drops by 215x at that corner — and VK3's post-
    yield modulus is E/215 against WSH3's E/119, because its real steel has a pronounced plateau to
    eps_s ~ 0.025 that Steel02 smears into one very flat slope. A sharp corner makes that drop
    happen inside a single step.
    """
    global S_LONG, S_HOOP
    S_LONG = replace(S_LONG, R0=r0)
    S_HOOP = replace(S_HOOP, R0=r0)


def set_steel_b(b: float) -> None:
    """DIAGNOSTIC ONLY — override the longitudinal steel's hardening ratio. Call BEFORE building.

    THIS OVERRIDES A MEASURED VALUE and must never be used for a production result. `b` = 0.00466
    follows directly from Tab. 5.4's f_y = 515, f_u = 630 and eps_su = 12.6%, and it faithfully
    reflects a real steel with a pronounced yield plateau to eps_s ~ 0.025 (Sec. 5.2.3a).

    It exists to test ONE hypothesis. A reinforcing bar in this lattice is a CHAIN of `mesh`-long
    rebar struts. With b this low the chain is nearly perfectly plastic, and a flat chain LOCALIZES:
    once one strut yields, every further increment of plastic strain goes into that strut alone
    rather than spreading along the bar. The model's post-yield branch needs the extreme bar to
    reach ~6*eps_y (1.55%) spread over a plastic hinge of a few hundred mm; localization instead
    drives it into a single 50 mm element.

    Raising `b` suppresses localization. If a run with an artificially hardening bar clears the
    0.310% drift barrier that every other lever failed to move, the mechanism is confirmed — and the
    legitimate fix is then a modelling decision (hinge-length regularization of the steel, or an
    explicit plateau law), NOT this switch.
    """
    global S_LONG
    S_LONG = replace(S_LONG, b=b)

# --- concrete grades ------------------------------------------------------------------------------
PIER = ConcreteGrade("pier", E=EC, nu=NU, rho=RHO_C,
                     fc=FC, epsc0=EPSC0, fcu=0.2 * FC, epsU=0.010, ft=FT)
# The foundation is FND_W thick against the pier's 350; a plane model carries one thickness, so the
# extra width is folded into an equivalent modulus. A restraint device, kept elastic.
FOUNDATION = ConcreteGrade("foundation", E=EC * (FND_W / TW), nu=NU, rho=RHO_C,
                           fc=FC, epsc0=EPSC0, fcu=0.2 * FC, epsU=0.010, ft=FT)
HEAD = ConcreteGrade("head", E=EC, nu=NU, rho=RHO_C,
                     fc=FC, epsc0=EPSC0, fcu=0.2 * FC, epsU=0.010, ft=FT)
GRADES = {"pier": PIER, "foundation": FOUNDATION, "head": HEAD}
GF = {"pier": GF_C, "foundation": GF_C, "head": GF_C}

# --- loading (Sec. 5.2.4) -------------------------------------------------------------------------
# Constant axial compression from two vertical actuators on a cross-beam at the pier TOP, tied to
# the strong floor by prestressing tendons. N_top = 1300 kN; the chapter's N_base = 1370 kN adds
# ~70 kN of pier, actuator and cross-beam weight, of which the model supplies the pier's own ~47 kN
# through its density. The tendons keep the force essentially vertical as the top displaces, which
# is what a constant nodal force does.
N_TOP = 1300.0e3
N_BASE = 1370.0e3
FY_PRIME = 644.0e3         # first-yield load, numerically computed (Tab. 5.6/5.8) — sets the
                           # amplitudes of the elastic force-controlled cycles
K0_MEASURED = 60.0e3       # N/mm, Tab. 5.10 — a SECANT to first yield, NOT an elastic stiffness


def pier_problem() -> Problem:
    """Pier on its foundation block, clamped at the foundation soffit."""
    domain = CompoundRectangles(
        rects=[
            (-FND_L / 2.0, -FND_H, FND_L, FND_H),        # foundation
            (-LW / 2.0, 0.0, LW, H_PIER),                # pier, full height including the head
        ],
        thickness=TW,
    )
    supports = [BoxSupport(box=(-FND_L, FND_L, -FND_H - EPS, -FND_H + EPS), fix=(1, 1))]
    loads = [BoxLoad(box=(-LW / 2.0, LW / 2.0, H_PIER - EPS, H_PIER + EPS), total=(0.0, -N_TOP))]
    return Problem(ndm=2, ndf=2, domain=domain, material=PIER, supports=supports, loads=loads)


def zone_of(x: float, y: float) -> str:
    if y < 0.0:
        return "foundation"
    return "head" if y >= HEAD_Y0 else "pier"


def snap(value: float, mesh_size: float) -> float:
    return round(value / mesh_size) * mesh_size


def rebars(mesh_size: float = MESH) -> tuple[Rebar, ...]:
    """Longitudinal bars (19 in-plane positions) and the hoops (as their in-plane legs).

    Bar positions are SNAPPED to the grid. At mesh 50 the true 80 mm spacing becomes an uneven
    50/100 pattern but all 19 positions stay DISTINCT (mesh 100 collapses four of them, which is
    why 50 is the working grid); `check_mesh_alignment` reports the resulting error in the
    end-bar-group centroid, which is what sets the flexural lever arm.
    """
    bars = [Rebar([(snap(x, mesh_size), -ANCHORAGE), (snap(x, mesh_size), H_PIER)],
                  MID_AREA, S_LONG) for x in WEB_X]
    bars += [Rebar([(snap(sgn * END_X, mesh_size), -ANCHORAGE),
                    (snap(sgn * END_X, mesh_size), H_PIER)], END_AREA, S_LONG) for sgn in (-1.0, 1.0)]
    bars += list(hoops(mesh_size))
    return tuple(bars)


def hoop_levels(mesh_size: float = MESH) -> tuple[list[float], list[float], float]:
    """`(base zone y-levels, head zone y-levels, head area scale)` on this grid.

    The BASE zone's 200 mm spacing is an exact 4 cells at mesh 50, so only its PHASE shifts (the
    first hoop moves 75 -> 100 mm); the spacing that carries the shear is preserved exactly. The
    head zone's 75 mm is not a whole number of cells, so it is snapped and its area rescaled by the
    same factor, holding the steel per unit height at the specimen's value.
    """
    y0, y1, s = HOOP_BASE
    base = [snap(y0, mesh_size) + i * s for i in range(int(round((y1 - y0) / s)) + 1)]
    y0h, y1h, sh = HOOP_HEAD
    s_head = max(mesh_size, snap(sh, mesh_size))
    head = [snap(y0h, mesh_size) + i * s_head
            for i in range(int(round((y1h - y0h) / s_head)) + 1)]
    return base, [y for y in head if y <= H_PIER], s_head / sh


def hoops(mesh_size: float = MESH) -> tuple[Rebar, ...]:
    """The phi6 hoops, projected onto the in-plane model.

    A hoop lies in a HORIZONTAL plane wrapping the 1500x350 section, so it has two LONG legs running
    along the section depth (in plane, and both projecting onto the same horizontal line) and two
    SHORT legs running through the 350 mm thickness (out of plane, and not representable here).
    The in-plane object is therefore ONE horizontal bar of area 2*A_phi6 per hoop level — which is
    exactly the steel the paper's rho_sw = 0.08% counts.

    These bars ARE the shear reinforcement, and they are the reason VK3 is worth modelling: their
    failure to take over the shear when the base compression zone was destroyed is what ended the
    test. NOTE the model's hoops can never lose anchorage, while the real ones were closed by simple
    90-degree hooks that "could open after spalling of the cover concrete" (Sec. 5.2.1b) — see the
    scope note in `summary.py`.
    """
    base, head, scale = hoop_levels(mesh_size)
    half = LW / 2.0
    out = [Rebar([(-half, y), (half, y)], HOOP_AREA, S_HOOP, role="stirrup") for y in base]
    out += [Rebar([(-half, y), (half, y)], HOOP_AREA * scale, S_HOOP, role="stirrup") for y in head]
    return tuple(out)


def check_mesh_alignment(mesh_size: float) -> None:
    """Fail on a grid that cannot host the model, and report the bar-snapping error."""
    for name, value in (("shear span L_v", L_V), ("pier height", H_PIER),
                        ("foundation depth", FND_H), ("foundation length", FND_L),
                        ("section depth l_w", LW), ("head zone start", HEAD_Y0),
                        ("hoop spacing (base)", HOOP_BASE[2])):
        if abs(value / mesh_size - round(value / mesh_size)) > 1e-9:
            raise ValueError(f"mesh_size={mesh_size:g} does not divide {name} ({value:g} mm)")
    xs = [snap(x, mesh_size) for x in WEB_X] + [snap(s * END_X, mesh_size) for s in (-1.0, 1.0)]
    if len(set(xs)) != len(xs):
        raise ValueError(f"mesh_size={mesh_size:g} collapses longitudinal bars onto shared "
                         f"positions ({len(set(xs))} distinct of {len(xs)}) — use 50 mm or finer")
    snapped = snap(END_X, mesh_size)
    print(f"bar snapping at mesh={mesh_size:g} mm: end bar {END_X:.1f} -> {snapped:.1f} mm "
          f"({abs(snapped - END_X) / END_X * 100:.2f}% on the flexural lever arm); "
          f"{len(set(xs))} distinct in-plane positions")


def lateral_loads(model) -> list[Load]:
    ids = select_nodes(model, (-LW, LW, L_V - EPS, L_V + EPS))
    return [Load(nid, (FY_PRIME / len(ids), 0.0)) for nid in ids]


def control_node(model) -> int:
    """The actuator axis at y = L_v — the point every displacement in the chapter refers to."""
    return select_nodes(model, (-EPS, EPS, L_V - EPS, L_V + EPS))[0]


def drive_nodes(model) -> list[int]:
    return select_nodes(model, (-LW, LW, L_V - EPS, L_V + EPS))


def base_nodes(model) -> list[int]:
    return select_nodes(model, (-FND_L, FND_L, -FND_H - EPS, -FND_H + EPS))


# --- cyclic loading protocol (Fig. 5.9 + Sec. 5.2.4) ---------------------------------------------
# Two symmetric cycles at every level, South (+) first.
#
#   ELASTIC phase: FORCE controlled at 0.25 / 0.50 / 0.75 / 1.00 * F_y', with F_y' = 644 kN taken
#   from the paper's own section analysis (Tab. 5.8), NOT from the experiment. The runners here are
#   displacement-driven, so these are converted through the MEASURED elastic stiffness k0 = 60 MN/m
#   (Tab. 5.10). The conversion validates itself: 1.00*F_y' comes out at 10.73 mm against the
#   measured first-yield displacement of 10.7 mm (Tab. 5.8).
#
#   INELASTIC phase: displacement controlled at mu_prov * 10.5 mm, mu_prov = 1, 1.5, 2, 3, 4, 5.
#   10.5 mm is the PROVISIONAL yield displacement the test was actually run on — deliberately shared
#   by all three units so their load steps matched. VK3's own nominal yield displacement is 14.1 mm
#   (Tab. 5.8), so its real ductilities are 0.74x the provisional ones (Tab. 5.9).
#
# THE "SMALL" INTERMEDIATE CYCLES ARE OMITTED. Fig. 5.9 inserts two small cycles between the two
# large ones from mu_prov = 1.5 to 5. For VK1 they were defined by force (+-0.50F_y', then +0.75 /
# -0.25 F_y'), but for VK2 and VK3 they were defined by "the top displacements MEASURED during the
# corresponding cycles of VK1" — numbers this chapter does not print. They are left out rather than
# invented; the consequence is that the model sees less small-amplitude cycling than the specimen
# did, which flatters cumulative damage slightly and is stated in `summary.py`.
DELTA_Y_PROV = 10.5                     # provisional nominal yield displacement, shared by VK1-VK3
MU_PROV = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
ELASTIC_FRACTIONS = (0.25, 0.50, 0.75, 1.00)
CYCLES_PER_LEVEL = 2

ELASTIC_PEAKS = tuple(f * FY_PRIME / K0_MEASURED for f in ELASTIC_FRACTIONS)
INELASTIC_PEAKS = tuple(m * DELTA_Y_PROV for m in MU_PROV)
PROTOCOL_PEAKS = ELASTIC_PEAKS + INELASTIC_PEAKS
PROTOCOL_CYCLES = (CYCLES_PER_LEVEL,) * len(PROTOCOL_PEAKS)


def protocol(max_drift: float | None = None) -> tuple[tuple[float, ...], tuple[int, ...]]:
    """The test protocol truncated at `max_drift` (a ratio, e.g. 0.01 for 1%); None keeps it all."""
    if max_drift is None:
        return PROTOCOL_PEAKS, PROTOCOL_CYCLES
    keep = [i for i, p in enumerate(PROTOCOL_PEAKS) if p <= max_drift * A_SHEAR + EPS]
    if not keep:
        raise ValueError(f"max_drift={max_drift} is below the first protocol amplitude "
                         f"({PROTOCOL_PEAKS[0] / A_SHEAR:.4%})")
    return tuple(PROTOCOL_PEAKS[i] for i in keep), tuple(PROTOCOL_CYCLES[i] for i in keep)


# --- how the protocol is DRIVEN (dynamic relaxation) ---------------------------------------------
# THE CHAPTER NEVER STATES AN ACTUATOR RATE. WSH3 at least printed 1.2-3.6 mm/min, which made the
# ~130x speed-up explicit; here there is not even a number to be embarrassed about. So the licence
# is entirely the MEASURED residual: `quasi_static=True` reports the inertia + damping the drive and
# base reactions fail to balance, and that is the number to read (D62/D64).
#
# 7.6 mm/s is kept — the same ABSOLUTE speed as both other wall studies, so the three are comparable
# and this pier's larger base shear (~890 vs 454 kN) makes the same contamination read about half as
# much of the peak.
QUASI_STATIC_RATE = 7.6
# DAMPING is the knob that controls contamination, not rate (D64). 0.2 is both the measured optimum
# there and a defensible physical value for heavily cracked RC.
DAMPING_RATIO = 0.2


# --- instrumentation: the test's own measurement chains -------------------------------------------
# Two instruments, both reproduced as the test built them, and both fed by ONE `node_history` probe
# that records u_x and u_y at every node of the grid below (D68).
#
#  1. VERTICAL CHAIN — ten LVDTs on each pier face measuring edge deformations, from which the
#     paper computes flexural and axial deformations (Sec. 5.2.2), backed by a 150x150 mm Demec grid
#     of 218 points to a height of 2.75 m. `LVDT_ROWS` sits on that Demec pitch, and 150 mm is
#     exactly 3 cells at mesh 50, so the model grid and the measurement grid coincide.
#
#  2. DIAGONAL CHAIN — four string potentiometers on the backside measuring diagonal elongations,
#     "in order to determine the shear deformations of the test units" (Sec. 5.2.2). That is the
#     measurement behind Fig. 5.19-right, and reproducing it is the whole reason this study exists.
#     See `panel_shear`.
LVDT_ROWS = (0.0, 150.0, 300.0, 450.0, 600.0, 750.0, 1050.0, 1350.0, 1650.0, 1950.0, 2250.0,
             2550.0, 2700.0)

# The BASE profile skips the lowest segment on purpose. Segment 0 (y = 0-150) is the BASE CRACK —
# the paper reports its contribution as a separate deformation component (Fig. 5.19 right), because
# it lumps strain penetration into the foundation together with reinforcement strain just above the
# joint. The model has no strain penetration at all, so including that segment in a "strain profile"
# would advertise a quantity the model does not possess. The base profile is rows 1..3, y = 150-450.
BASE_SEGMENT = (1, 3)
GAUGE_Y0, GAUGE_H = LVDT_ROWS[BASE_SEGMENT[0]], LVDT_ROWS[BASE_SEGMENT[1]] - LVDT_ROWS[BASE_SEGMENT[0]]

# Nominal yield curvature from the paper's own moment-curvature analysis (Tab. 5.5), used to define
# the plastic-zone height the base-curvature fit is allowed to use.
PHI_YIELD = 3.36e-6                     # 3.36 /km, in 1/mm
PHI_YIELD_FIRST = 2.54e-6               # first-yield curvature, Tab. 5.5
PHI_ULT = 14.5e-6                       # ultimate curvature at eps_cu = -0.005, Tab. 5.5


def gauge_nodes(model, *, rows: tuple[float, ...] = LVDT_ROWS) -> tuple[list[float], list[list[int]]]:
    """Aligned node columns at each row height: `(xs, ids_per_row)`.

    `ids_per_row[k][i]` sits at x = `xs[i]`, y = `rows[k]`, so differencing any two rows column by
    column gives the average vertical strain of that column over that segment, and the four corners
    of any two rows give a shear panel.
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


def split_sample(sample, n_rows: int, n_x: int) -> tuple[list[float], list[float]]:
    """One `node_history` record -> `(ux, uy)`, each `n_rows * n_x` long in `LVDT_ROWS` order.

    The probe is asked for dofs `(1, 2)`, and `_probe_sample` stores them dof-MAJOR, so the record
    is the whole u_x field followed by the whole u_y field (D68).
    """
    n = n_rows * n_x
    if len(sample) == n:                       # single-dof record (a pre-D68 uy-only probe)
        return [0.0] * n, list(sample)
    return list(sample[:n]), list(sample[n:2 * n])


def segment_strain(sample, n_x: int, lo: int, hi: int,
                   rows: tuple[float, ...] = LVDT_ROWS) -> list[float]:
    """Average vertical strain per column over the segment between rows `lo` and `hi`.

    Positive = tension (the column got longer).
    """
    _ux, uy = split_sample(sample, len(rows), n_x)
    a = uy[lo * n_x:(lo + 1) * n_x]
    b = uy[hi * n_x:(hi + 1) * n_x]
    length = rows[hi] - rows[lo]
    return [(t - s) / length for s, t in zip(a, b)]


def gauge_strains(node_history, xs, segment: tuple[int, int] = BASE_SEGMENT) -> list[list[float]]:
    """Strain profile across the width for one segment, for every sample in a `node_history`."""
    n = len(xs)
    return [segment_strain(v, n, *segment) for v in node_history["values"]]


def curvature_profile(node_history, xs, *, rows: tuple[float, ...] = LVDT_ROWS):
    """`(mid_heights, [[curvature per segment] per sample])`. Units 1/mm."""
    n = len(xs)
    mids = [0.5 * (rows[k] + rows[k + 1]) for k in range(len(rows) - 1)]
    out = []
    for vals in node_history["values"]:
        out.append([neutral_axis(xs, segment_strain(vals, n, k, k + 1, rows))[0]
                    for k in range(len(rows) - 1)])
    return mids, out


# --- the shear instrument (D68): the test's four string potentiometers ---------------------------
def panel_shear(sample, xs, *, rows: tuple[float, ...] = LVDT_ROWS) -> list[float]:
    """Average shear strain of each panel between successive rows, from its two DIAGONALS.

    This is the string-potentiometer measurement of Sec. 5.2.2, and it is why `node_history` had to
    learn to carry two dofs: a diagonal's length change depends on both u_x and u_y at each end, so
    neither component alone determines it.

    For a panel of width `b` and height `h` whose diagonals change length by `dd1` (the one rising
    to the right) and `dd2`,

        gamma = sqrt(b^2 + h^2) * (dd1 - dd2) / (2 * b * h)

    which is the standard reduction used with diagonal gauges on RC walls. The panel's shear
    displacement is then `gamma * h`, and summing those over the height gives the shear part of the
    top displacement — the quantity Fig. 5.19-right plots as a percentage.
    """
    n = len(xs)
    ux, uy = split_sample(sample, len(rows), n)
    b = xs[-1] - xs[0]
    out = []
    for k in range(len(rows) - 1):
        h = rows[k + 1] - rows[k]
        lo, hi = k * n, (k + 1) * n
        # corners: (left,bottom) (right,bottom) (left,top) (right,top)
        lb, rb = (ux[lo], uy[lo]), (ux[lo + n - 1], uy[lo + n - 1])
        lt, rt = (ux[hi], uy[hi]), (ux[hi + n - 1], uy[hi + n - 1])
        d0 = math.hypot(b, h)
        d1 = math.hypot(b + rt[0] - lb[0], h + rt[1] - lb[1])     # left-bottom -> right-top
        d2 = math.hypot(b + rb[0] - lt[0], -h + rb[1] - lt[1])    # left-top -> right-bottom
        out.append(d0 * (d1 - d2) / (2.0 * b * h))          # = d0*(dd1 - dd2)/(2bh)
    return out


def deformation_components(sample, xs, *, rows: tuple[float, ...] = LVDT_ROWS,
                           height: float = A_SHEAR) -> dict:
    """Split one sample's top displacement the way Fig. 5.19-right splits the measured one.

    Returns absolute contributions in mm: `shear`, `flexure_base_crack`, `flexure_above`, `total`
    (the four components the paper distinguishes, minus base sliding — see below).

    FLEXURE is integrated from the curvature profile the vertical chain gives: a segment of height
    `h_k` at mid-height `y_k` rotates by `phi_k * h_k` and moves the top by `phi_k * h_k *
    (height - y_k)`. Segment 0 is reported separately because the paper reports it separately: it is
    the BASE CRACK, whose opening it treats as its own deformation component.

    SHEAR comes from `panel_shear`, summed as `gamma_k * h_k` over the panels.

    BASE SLIDING is NOT returned. The paper measures it (two LVDTs across the construction joint)
    and finds it negligible, "hardly exceeding the range of the measuring tolerance". The model has
    no joint at all — the pier and foundation share nodes — so its sliding is identically zero by
    construction, and reporting a zero as though it were a result would be misleading.
    """
    n = len(xs)
    gam = panel_shear(sample, xs, rows=rows)
    shear = sum(g * (rows[k + 1] - rows[k]) for k, g in enumerate(gam))
    flex = []
    for k in range(len(rows) - 1):
        phi = neutral_axis(xs, segment_strain(sample, n, k, k + 1, rows))[0]
        h = rows[k + 1] - rows[k]
        # phi here is d(eps)/dx, which is MINUS the beam curvature d2u/dy2 (eps = -kappa*x), hence
        # the sign: a pier pushed to +x carries tension on its -x face, giving phi < 0.
        flex.append(-phi * h * (height - 0.5 * (rows[k] + rows[k + 1])))
    return {"shear": shear, "flexure_base_crack": flex[0], "flexure_above": sum(flex[1:]),
            "flexure": sum(flex), "modelled_total": shear + sum(flex)}


def neutral_axis(xs, strain) -> tuple[float, float, float]:
    """Least-squares line through one strain profile: `(curvature, strain_at_centre, x_zero)`.

    A plane section would make the profile exactly linear, so the fit is both a summary (curvature
    = the slope, in 1/mm) and a test of that assumption. That test matters more for VK3 than for
    either previous wall: this specimen's inclined flexure-shear cracks are large (up to 4 mm at
    mu_prov = 4) and Sec. 5.3.1 reports the shear cracks growing WIDER than the flexural ones, which
    is precisely the condition under which plane sections stop holding.

    `x_zero` is where the fitted line crosses zero — the neutral axis — and is `nan` for a profile
    with no gradient (a purely axial state, before cracking localizes).
    """
    n = len(xs)
    mx = sum(xs) / n
    me = sum(strain) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxe = sum((x - mx) * (e - me) for x, e in zip(xs, strain))
    slope = sxe / sxx if sxx else 0.0
    intercept = me - slope * mx
    x0 = -intercept / slope if abs(slope) > 1e-18 else float("nan")
    return slope, intercept, x0


# --- the MEASURED drive history (recovered by `digitize.py`) --------------------------------------
# `PROTOCOL_PEAKS` above is the protocol as DESIGNED. The protocol as RUN is recoverable from the
# digitized Fig. 5.13 loops, because those are vector paths in drawing order, and it differs from
# the design in two ways that matter:
#
#   1. THE ELASTIC AMPLITUDES ARE NOT WHAT k0 PREDICTS. Converting 0.25/0.50/0.75 F_y' through the
#      measured k0 = 60 kN/mm gives 2.68 / 5.37 / 8.05 mm; the pier actually reached 0.85 / 2.10 /
#      5.08 mm. k0 is a secant to FIRST YIELD on an already-cracked pier, so using it for pre-
#      cracking amplitudes overstates them by up to 3x. The measured values are the honest ones.
#   2. THE SMALL INTERMEDIATE CYCLES ARE RECOVERABLE. The chapter defines VK3's small cycles only as
#      "the top displacements measured during the corresponding cycles of VK1", numbers it never
#      prints — but they are plainly in the figure, including the deliberately ASYMMETRIC second one
#      of each pair. So they need not be invented or omitted.
#
# `measured_protocol` therefore returns the pier's own turning points as a ready-made `history` for
# the cyclic runners, and `protocol()` above stays as the fallback for when the digitized file has
# not been produced yet.
MEASURED_LOOPS = Path(__file__).resolve().parent / "data" / "vk3_fig513.npz"


def measured_protocol(max_drift: float | None = None, *, include_small: bool = True):
    """The test's OWN turning points as a reversal history, or None if not digitized yet.

    `include_small=False` drops the small intermediate cycles, which is the like-for-like comparison
    against `protocol()`; keeping them is the faithful drive.
    """
    if not MEASURED_LOOPS.exists():
        return None
    import numpy as np
    d = np.load(MEASURED_LOOPS, allow_pickle=False)
    tips = [float(v) for v in d["tip_disp_mm"]]
    if not include_small:
        flags = [bool(v) for v in d["tip_is_primary"]]
        tips = [u for u, keep in zip(tips, flags) if keep]
    if max_drift is not None:
        limit = max_drift * A_SHEAR + EPS
        out = []
        for u in tips:
            if abs(u) > limit:
                break
            out.append(float(u))
        tips = out
    return [float(u) for u in tips]
