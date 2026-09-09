"""A REPLICA of Aydin, Tuncay & Binici's (2019) own lattice model of the Aldemir wall.

Not our model with his parameters — his model, as closely as this codebase can express it. The
point is to separate two questions that six failed experiments have been conflating:

    Does OUR model collapse because of how WE built it, or because of something the lattice
    method does to this specimen?

If the replica reproduces his 1164 kN plateau, the difference is ours and is findable. If the
replica collapses too, then either the difference is in the one thing we cannot yet replicate
(bond, see below) or his published curve is not reproducible from his published parameters.

EVERY CHOICE BELOW IS HIS, and each is sourced. Where the paper is silent, that is said plainly
rather than filled in from our study.

Units: N, mm.
"""

from __future__ import annotations

import math
import os

from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, Rebar, RectangleDomain, SteelGrade

SPECIMEN = "Aldemir wall — Aydin (2019) replica"

# --- GEOMETRY, from his Table 2 rather than from any figure ---------------------------------------
# Table 2 reports 20,385 particles and 80,684 elements at d = 20 mm. For a horizon-1.5 structured
# grid those two integers have a UNIQUE solution (up to transposition): 150 x 134 cells.
#
#     nodes    = 151 * 135                                     = 20,385   exactly
#     elements = 150*135 + 151*134 + 2*150*134                 = 80,684   exactly
#
# So his analysed panel was 3000 x 2680 mm, matching his Fig. 10(a) drawing (3000 x 2250) no better
# than it matches anything else. We take the arithmetic over the drawing: it is the only statement
# about his MODEL geometry that is internally verifiable. It also proves the Table 2 element count
# is CONCRETE STRUTS ONLY — steel and bond elements are not in it.
#
# THE TRANSPOSITION IS UNRESOLVED, AND IT IS A REAL RISK TO EVERY NUMBER BELOW (D82). 150 x 134 and
# 134 x 150 reproduce both integers equally well, i.e. 3000 x 2680 (squat) and 2680 x 3000 (slender).
# The line taken here is the first. Evidence for the second: Fig. 4(f) dimensions 2680 along the
# HORIZONTAL, and that figure is a DETAIL VIEW of the bottom 1500 mm of the model (author,
# 2026-09-05) rather than a panel, so its width label would be the model's true full width.
# If the second reading is right, this replica has his wall on its side and every published ratio
# here — K/K_sim, F/F_sim, the drift denominator — moves. Needs the author before it is changed.
LW, HW = 3000.0, 2680.0
# NEVER STATED by the paper. 120 is what the Fig. 10(a) plan brackets; 210 is what his own
# K_sim = 943.16 kN/mm requires (stiffness is exactly linear in thickness) and is also 50+100+60
# from that plan's stack. Override for the comparison: ALDEMIR_TW=210.
TW = float(os.environ.get("ALDEMIR_TW", 120.0))
MESH = 20.0       # Table 1 and Table 2 (his Fig. 10b legend says 25 — unreconciled)
HORIZON = 1.5     # Table 4 row we are targeting; he also ran 3.01
A_SHEAR = HW
EPS = 1e-6

# --- MATERIALS, Table 1, "Wall / Aldemir et al. (2017)" column ------------------------------------
FC, FT, EC = 28.0, 1.85, 24870.0    # f_t is TS 500, E_t is ACI 4700*sqrt(f_c) — both code values
GF = 0.075                          # N/mm (his 75 N/m)
NU = 1.0 / 3.0                      # HIS value, from the Appendix (Hrennikoff); ours is 0.20
RHO_C = 2.4e-9
FY, ES = 360.0, 200000.0
EPSC0 = 2.0 * FC / EC

# His trilinear tension tail, taken as PUBLISHED rather than solved from Gf. At d = 20 the energy
# identity would give a2 = 33.7 / a3 = 173.3; he prints 70 / 360, about 2.1x longer. Using his
# numbers is the whole point of a replica, and it deliberately gives up mesh objectivity.
# a1 IS THE HORIZON. Table 1's footnote d ties them: a1 = 1.5 at delta = 1.5d and a1 = 3.01 at
# delta = 3.01d. So a horizon sweep must move a1 with it — holding a1 = 1.5 at horizon 3.01 would
# silently run his 1.5d softening law on a 3.01d lattice. `build.wall_lattice` takes a1 from the
# horizon it is actually given, not from this constant.
A1, A2, A3 = HORIZON, 70.0, 360.0
B1, B2 = 0.6, 0.2                   # NOT PUBLISHED — the paper never prints b1/b2. Repo defaults.

CONCRETE = ConcreteGrade(name="aydin-wall", E=EC, nu=NU, rho=RHO_C, fc=FC, ft=FT, Gf=GF,
                         epsc0=EPSC0, fcu=0.2 * FC, epsU=8.0 * EPSC0)

# ELASTIC–PERFECTLY PLASTIC, his Fig. 1(d): no hardening at all. Our own studies use b = 0.01.
STEEL = SteelGrade(name="aydin-steel", fy=FY, E0=ES, b=0.0)

# --- BOND, his Fig. 1(c) ---------------------------------------------------------------------
# `a` = 0.7 is the ONLY bond number the paper prints (his Sec. "Reinforced Concrete Lattice
# Modeling", carried over from Aydin 2017, where it was calibrated for DEFORMED bars). It is the
# residual plateau the link drops to after the brittle fall, and it is load-bearing in a second,
# structural sense: a steel node is held only by its own ring of links, so a = 0 leaves it free.
#
# EVERYTHING ELSE ABOUT HIS BOND IS UNPUBLISHED — link area, peak slip, peak force. Those are
# inferred in `build.py` and are the largest remaining free parameters in this replica.
BOND_RESIDUAL = 0.7

# --- REINFORCEMENT, Fig. 10(a): "3- Ø8 bars (@100 mm)" both directions ---------------------------
D_BAR, N_PER_POSITION, S_BAR = 8.0, 3, 100.0
BAR_AREA = N_PER_POSITION * math.pi * D_BAR ** 2 / 4.0
# Cover must land on his 20 mm grid. The Fig. 10(a) measurement gives 50 (vertical) and 75
# (horizontal); the nearest grid values are 40 and 80, both of which preserve the drawn bar COUNT.
COVER_X, COVER_Y = 40.0, 80.0

# --- what CANNOT be replicated ---------------------------------------------------------------
# A FUNCTION, not a constant: two of these entries are switchable now, and a run that prints the
# unconditional list would over-state its own departures. Every run prints the applicable set.
_ALWAYS = (
    "b1 and b2 of the trilinear tail, which the paper never prints.",
    "His model THICKNESS, which the paper never states; 120 mm is read off the Fig. 10(a) plan "
    "and 210 is what his own K_sim requires (D73). Set by ALDEMIR_TW.",
)

_NO_BOND = (
    "BOND ELEMENTS. His steel nodes are separate from the concrete nodes and tied to them by a ring "
    "of elastic-brittle links with residual a = 0.7. This run uses SHARED NODES instead. This is "
    "the single largest known departure and the first thing to suspect if the replica behaves "
    "differently — pass --bond to remove it.",
)

_WITH_BOND = (
    "BOND PARAMETERS (not the elements — those are now present). Only his residual ratio a = 0.7 "
    "is published; the link AREA, PEAK SLIP and PEAK FORCE are inferred here (D75), so agreement "
    "or disagreement on bond is a statement about our inference, not about his numbers.",
    "BOND TOPOLOGY. His links tie a steel node to the concrete; ours are a horizon RING, which "
    "cannot reduce to perfect bond as its stiffness grows (D72). Calibration bounds that artefact "
    "to ~1% of elastic stiffness rather than removing it; a coincident-node zeroLength spring "
    "would remove it and is unbuilt.",
)

_IMPLICIT = (
    "EXPLICIT TIME INTEGRATION with PID control (dt = 5e-8 s, Table 2). We solve implicitly with "
    "dynamic relaxation. An explicit scheme needs no tangent, so it walks through softening and "
    "local instability that an implicit Newton has to converge past.",
)

_EXPLICIT = (
    "His PID FORCE control (Table 2) — we prescribe displacement at the top row instead — and his "
    "dt = 5e-8 s, which is ~2 orders below the stable step our own element sizes allow.",
)


_CYCLIC = (
    "HIS CONCRETE LAW, in a cyclic run. His tension-only material is an ElasticMultiLinear, which "
    "is path-independent and therefore meaningless under reversal — a cracked strut would recover "
    "full stiffness on reload. A cyclic run substitutes `concrete_lattice_aydin_cyclic`, which "
    "keeps his trilinear tension SHAPE but replaces his tension-only compression with a "
    "Concrete02-like envelope. So a cyclic run is OUR lattice using HIS softening law, not a "
    "replica of his model, and its compressive behaviour is ours.",
    "THE LOADING PROTOCOL. The paper never prints the number, amplitude or order of cycles "
    "(`testdata.NOT_REPORTED`), and the Fig. 10(b) cloud does not reveal it either — the "
    "displacement histogram is flat, with no amplitude ladder. Any protocol is our assumption.",
)


def not_replicated(*, bond: bool = False, explicit: bool = False,
                   cyclic: bool = False) -> tuple[str, ...]:
    """The departures that still apply, given what this run switched on."""
    return ((_WITH_BOND if bond else _NO_BOND) + (_EXPLICIT if explicit else _IMPLICIT)
            + (_CYCLIC if cyclic else ()) + _ALWAYS)


def zone_of(_x: float, _y: float) -> str:
    return "wall"


GRADES = {"wall": CONCRETE}


def bar_lines(extent: float, cover: float, spacing: float = S_BAR) -> list[float]:
    n = int((extent - 2.0 * cover) // spacing) + 1
    return [cover + i * spacing for i in range(n)]


def rebars(mesh_size: float = MESH, *, full_height: bool = False) -> tuple[Rebar, ...]:
    """Two orthogonal bar layers. `full_height` runs the LONGITUDINAL bars to the very top.

    WHY THIS IS AN OPTION AT ALL. The bar lines are placed on a cover + spacing grid, so the top
    line lands at y = 2580 and leaves a **100 mm band of PLAIN CONCRETE** above it — five rows of
    mesh-20 cells with no longitudinal steel. That band is also where the drive is applied, so the
    whole base shear is introduced through unreinforced material.

    All three collapses measured so far tore inside it (y = 2630-2670 at horizon 3.01, 2570-2650 at
    horizon 1.5, and 2550-2610 with an elastic cap), each after a peak at which damage was still
    correctly distributed over the middle and base. `full_height=True` carries the longitudinal bars
    to y = HW so the loaded region is reinforced like the rest of the wall.

    It is an OPTION and not the default because the paper says nothing about bar termination — the
    grid placement above is an inference, and so is this. Which one is right is unknown; what is
    measurable is how much the answer depends on it.
    """
    xs = bar_lines(LW, COVER_X)
    ys = bar_lines(HW, COVER_Y)
    top = HW if full_height else ys[-1]
    bars = [Rebar(path=[(x, 0.0), (x, top)], area=BAR_AREA, steel=STEEL, role="longitudinal")
            for x in xs]
    bars += [Rebar(path=[(xs[0], y), (xs[-1], y)], area=BAR_AREA, steel=STEEL, role="stirrup")
             for y in ys]
    return tuple(bars)


def wall_problem() -> Problem:
    """Fixed base, lateral drive across the top row, NO axial load (Table 1 leaves N blank)."""
    return Problem(
        ndm=2, ndf=2,
        domain=RectangleDomain(length=LW, height=HW, thickness=TW),
        material=CONCRETE,
        supports=(BoxSupport(box=(-EPS, LW + EPS, -EPS, EPS), fix=(1, 1)),),
        loads=(BoxLoad(box=(-EPS, LW + EPS, HW - EPS, HW + EPS), total=(1.0, 0.0)),),
    )


# --- his published results, the targets -----------------------------------------------------------
# HIS TABLE 4 HAS TWO ROWS, ONE PER HORIZON, and they differ by 14% in K and 14% in F. A run must
# be compared against its OWN row: quoting the 1.5 targets at horizon 3.01 understates the model by
# 12%. Keyed by horizon for exactly that reason.
PAPER_BY_HORIZON = {
    1.5:  {"K_sim_kN_per_mm": 943.16,  "F_sim_kN": 1164.413, "elements": 80_684},
    3.01: {"K_sim_kN_per_mm": 1052.37, "F_sim_kN": 1325.675, "elements": 280_260},
}
MEASURED_PAPER = {"K_exp_kN_per_mm": 1038.44, "F_exp_kN": 963.592,
                  "nodes": 20_385, "max_disp_mm": 20.0}


def paper_for(horizon: float) -> dict:
    """His Table 4 row for this horizon, merged with the measured values (horizon-independent)."""
    key = min(PAPER_BY_HORIZON, key=lambda h: abs(h - horizon))
    if abs(key - horizon) > 1e-6:
        raise SystemExit(f"no published Table 4 row for horizon {horizon}; he ran "
                         f"{sorted(PAPER_BY_HORIZON)} only — comparing against another row would "
                         f"be meaningless")
    return {**MEASURED_PAPER, **PAPER_BY_HORIZON[key], "horizon": key}


PAPER = paper_for(HORIZON)      # backwards-compatible default (horizon 1.5)
