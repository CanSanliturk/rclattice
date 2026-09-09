"""Shared specimen definition for the squat RC wall of Aldemir, Binici & Canbay (2017).

TEST SOURCE:  Aldemir, A., B. Binici, and E. Canbay (2017), "Cyclic testing of reinforced concrete
double walls", ACI Structural Journal 114(2): 395-406, doi:10.14359/51689432.

LATTICE SOURCE:  Aydin, B. B., K. Tuncay, and B. Binici (2019), "Simulation of Reinforced Concrete
Member Response Using Lattice Model", J. Struct. Eng. 145(9): 04019091,
doi:10.1061/(ASCE)ST.1943-541X.0002381 — Table 1 (properties), Table 2 (lattice size), Table 4
(measured vs computed), Fig. 4(f) (their model) and Fig. 10 (specimen, response, damage).

The PDF of the 2019 paper sits next to this module; the 2017 test paper is NOT in the repo, so every
number here is traced to the 2019 paper and `testdata.py` records which figure or table it came from.

WHY THIS SPECIMEN. It is the squattest wall in the repo by a wide margin — 3000 long against a
2250 shear span, ratio 0.75, against SW-NC-FF's 3.00, WSH3's 2.28 and VK3's 2.20. It carries NO
axial load and was "designed to yield in shear" (2019 paper, "RC Wall Simulations"), reaching ~1%
drift with severe inclined cracking and NO strength degradation. Where VK3 was the first specimen
whose shear mattered, this is the first whose shear is the whole mechanism.

It is also the specimen the 2019 paper predicts WORST: +21% on peak force at horizon 1.5d and +38%
at 3.01d, against 4-11% for the other five. That over-prediction is the thing worth attacking, and
the two features this repo just gained — the published `field="equibiaxial"` calibration route and
optional bond elements — are both plausible contributors.

BOND. Unknown. The 2017 title says "double walls", a precast twin-shell system, and neither the
bar type (deformed vs plain) nor the shell/core interface is described in the 2019 paper. Perfect
bond is therefore an ASSUMPTION here, not a justified idealization the way it is for WSH3 and VK3.
This is the study to run WITH the new bond elements (`--bond`) precisely because the honest default
is unknown; see `build.wall_lattice`.

Units: N, mm (stresses MPa, forces N).
"""

from __future__ import annotations

import math
import os
import sys
from datetime import datetime
from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, Rebar, RectangleDomain, SteelGrade

SPECIMEN = "Aldemir et al. (2017) squat wall"

# --- geometry (2019 paper, Fig. 10a) -------------------------------------------------------------
#
# Fig. 10(a) is a raster panel with no dimension text beyond the four numbers below, so it was
# measured rather than read: the outer rectangle spans 1756 x 1317 px against the printed 3000 and
# 2250, giving a consistent 1.708 mm/px, and the drawn bar mesh comes to 30 vertical and ~22-23
# horizontal lines at 55.4 px = 94.7 mm spacing. That is Ø8 @ 100 mm at 50 mm cover, and it is the
# whole basis for `COVER` below.
LW = 3000.0        # wall length (in-plane), Fig. 10(a) elevation
HW = 2250.0        # wall height = shear span; the load is applied at the top
# WALL THICKNESS = 210 mm, NOT the 120 that the Fig. 10(a) plan brackets (D73, 2026-08-29).
#
# 120 is what the plan dimensions against ONE panel of a precast "double wall"; the structural
# section is the whole stack that carries load. Five independent checks put it at 210:
#
#   1. The Aydin REPLICA (replica/) reproduces his Table 2 counts EXACTLY — 20,385 nodes and
#      80,684 concrete struts — and then misses his published K_sim = 943.16 kN/mm by a factor of
#      1.743. Lattice stiffness is EXACTLY linear in thickness (A_t = C*d*w, so every strut area
#      is linear in w), which makes that factor directly readable: 120 * 1.743 = 209.1 mm.
#   2. Run at 210 the replica returns K = 923.7 kN/mm = 0.979 of his K_sim — a number it was
#      never fitted to.
#   3. 210 = 50 + 100 + 60 of the Fig. 10(a) plan's own 50|140|100|60 stack: two precast shells
#      plus a cast-in-place core, with 140 the clear gap between shell faces.
#   4. At 210 the UNCRACKED transformed section finally sits ABOVE the measured 1,038 kN/mm
#      (1.31x, as it must for a fixed-base cantilever). At 120 it sat BELOW, at 0.75x, which is
#      impossible and was the anomaly `summary.py` flagged for two days.
#   5. At 210 the distributed-steel flexural capacity lands on the measured force: V_flex/F_exp
#      = 1.00 (at 120 it was 0.92, and rho was overstated 1.257% against the true 0.718%).
#
# CONSEQUENCE FOR EVERY RUN BEFORE 2026-08-29: they modelled a wall 43% too thin, so every strut
# carried 1.75x too little area while the reinforcement stayed at full size. Peak base shear went
# 591 -> 903 kN (0.61 -> 0.94 of the measured) on the thickness change alone.
# `ALDEMIR_TW` overrides it for A/B work.
TW = float(os.environ.get("ALDEMIR_TW", 210.0))

# THE PANEL IS SETTLED: 3000 WIDE x 2250 HIGH (user, 2026-08-28), which is what Fig. 10(a) prints
# and what the drawing independently measures. The elevation is to scale — its two printed dimensions
# give 1.3529 and 1.3505 mm/px, agreeing to 0.17% — and every feature reconstructs from those pixels:
# 30 vertical bar lines at 100.0 mm with 52.1/46.7 mm cover rebuild the 3000, and 22 horizontal lines
# at 99.8 mm with 72.9/80.4 mm cover rebuild the 2250. Four independent quantities, one scale.
#
# FIG. 4(f) IS A DETAIL VIEW, NOT A PANEL (author, 2026-09-05; D82). It shows the BOTTOM 1500 mm of
# his model, which is why that number appears at all — so the "1500" is the height of a CROP and
# carries no model dimension. It was never a candidate geometry and is no longer offered as one.
#
# This DISSOLVES the anomaly recorded here until 2026-09-05 rather than explaining it: the earlier
# note called the figure "NOT to scale" because its drawn block's aspect (2.235) matched none of its
# own labels. A crop is exactly what produces that mismatch, so no drafting error need be assumed and
# the figure's status changes from unreliable to reliable about something else.
#
# THE "2680" DIMENSIONED ALONG THE HORIZONTAL THERE MAY STILL BE REAL, AND IT BEARS ON AN OPEN
# QUESTION. Inverting his Table 2 counts gives 150 x 134 cells, unique only UP TO TRANSPOSITION:
# 3000 x 2680 and 2680 x 3000 are equally consistent with both integers. `replica/specimen.py` takes
# the first. A FULL-WIDTH crop labelled 2680 implies the second, i.e. a model taller than wide.
# NOT RESOLVED — it needs the author, and it would move every published replica ratio. Both readings
# are exposed below so `elastic.py --geometry-sweep` can price the difference.
#
# ONE DISCREPANCY REMAINS, AND IT IS NOT ABOUT THE PANEL SIZE. The measured initial stiffness of
# 1,038 kN/mm exceeds what an UNCRACKED gross section of 3000 x 2250 x 120 can give (804 kN/mm
# transformed), which is physically impossible for a fixed-base cantilever. With the panel settled,
# the candidates are the THICKNESS (the 2017 title says "double walls"; two 120 mm shells would give
# 1,608), the TOP BOUNDARY (restrained against rotation gives 1,156), or how Table 4 defines "initial
# stiffness". `summary.py` prices all three. Expect K_lattice/K_measured ~ 0.71 and judge the
# calibration on K_lattice/K_transformed instead — the other three walls land at 0.93.
FIG4F_CROP_H = 1500.0                  # height of the Fig. 4(f) DETAIL VIEW — not a model dimension
PAPER_GRID = (3000.0, 2680.0)          # Table 2 inversion, as `replica/` resolves the transposition
PAPER_GRID_T = (2680.0, 3000.0)        # the other resolution, equally consistent with both integers

A_SHEAR = HW      # drift denominator: drift = u_top / HW

# 50 mm DIVIDES EVERYTHING. 3000/50 = 60, 2250/50 = 45, bar spacing 100/50 = 2, cover 50/50 = 1 —
# no geometric rounding anywhere, which is not true of the paper's own d = 20 (2250/20 = 112.5).
MESH, HORIZON = 50.0, 1.5
EPS = 1e-6

OUT = Path(__file__).resolve().parent.parent / "output" / "aydin_aldemir_wall"
RUNS = OUT / "runs"


def run_dir(kind: str, *, tag: str | None = None) -> Path:
    """A fresh TIMESTAMPED directory for one run (D70) — nothing is ever overwritten.

    Same contract as `vk3_wall.specimen.run_dir`: `command.txt` beside the results, and a
    `latest_<kind>` symlink for `replot`/`compare` to follow.
    """
    name = f"{datetime.now():%Y-%m-%d_%H%M%S}_{kind}" + (f"_{tag}" if tag else "")
    d = RUNS / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "command.txt").write_text(" ".join(sys.argv) + "\n")
    link = RUNS / f"latest_{kind}"
    try:
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(d.name)
    except OSError:
        pass
    return d


def find_run(kind: str | None = None, *, name: str | None = None, complete: bool = True) -> Path:
    """Locate a run directory: an explicit `name`, else the NEWEST COMPLETED one of `kind`."""
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
                         + (f" matching {kind!r}" if kind else "") + " — run an analysis first")
    if complete and not done:
        raise SystemExit(f"no COMPLETED {kind or 'run'} under {RUNS} — {len(got)} directory(ies) "
                         f"exist but none has a data.json. Newest: {got[-1].name}")
    return got[-1]


# --- reinforcement (2019 paper, Fig. 10a) --------------------------------------------------------
#
# "3- Ø8 bars (@100 mm)" is called out TWICE in Fig. 10(a), once against the vertical bars and once
# against the horizontal, so the mesh is the same both ways: an orthogonal grid at 100 mm with three
# Ø8 at every position. Three bars in a 120 mm wall is unusual, and it is what the drawing says; the
# ratio it produces is the reason to believe it — see `summary.py`, where rho = 1.26% reproduces the
# measured peak force to within a few percent and 2 bars (0.84%) falls ~25% short.
D_BAR = 8.0
N_PER_POSITION = 3
S_BAR = 100.0

# COVER IS DIFFERENT IN THE TWO DIRECTIONS, and both values are MEASURED off Fig. 10(a) rather than
# assumed. The elevation is drawn to scale — its two printed dimensions give 1.3529 and 1.3505 mm/px,
# agreeing to 0.17% — so pixel positions convert to millimetres reliably:
#
#   vertical bars    30 lines @ 100.0 mm, cover 52.1 / 46.7  -> 50 mm, cage 2900, panel 3000 ✓
#   horizontal bars  22 lines @  99.8 mm, cover 72.9 / 80.4  -> 75 mm, cage 2100, panel 2250 ✓
#
# Both reconstruct the printed panel exactly, which is the check that the reading is right.
COVER_X = 50.0        # cover to the VERTICAL bars (measured off the 3000 mm dimension)
COVER_Y = 75.0        # cover to the HORIZONTAL bars (measured off the 2250 mm dimension)
COVER = COVER_X       # backwards-compatible alias; prefer the explicit pair

# 75 mm is NOT a multiple of the 50 mm grid, so the horizontal bars snap. Snapping DOWN to 50 keeps
# all 22 lines (50 ... 2150); snapping up to 100 would drop one to 21 and change rho by 4.5%. The
# count is what sets the ratio, so the count is what is preserved — at the price of the whole
# horizontal cage sitting 25 mm low. `check_bar_layout` reports the offset on every build; `--mesh 25`
# removes it.

BAR_AREA = N_PER_POSITION * math.pi * D_BAR * D_BAR / 4.0     # 150.8 mm^2 per in-plane position

# The foundation (4Ø16 + Ø6@100 in the Fig. 10a plan) is NOT modelled. That follows the paper, which
# says of its column simulations that "the foundation of the specimens and the strong floor was not
# modeled, and fixed supports were placed at the bottom boundary", and whose Fig. 4(f) wall model is
# a bare rectangle. It also removes an unknown: the plan gives the foundation's length but never its
# depth or its out-of-plane thickness.

# --- materials (2019 paper, Table 1, "Wall / Aldemir et al. (2017)" column) -----------------------
FC = 28.0                       # MPa
FT = 1.85                       # MPa — footnote a: TS 500 (TSE 2000), not a measured value
EC = 24870.0                    # MPa — footnote b: ACI 318, Ec = 4700*sqrt(fc) = 24,869.7
GF = 0.075                      # N/mm  (the table's 75 N/m)
NU = 0.20                       # not given by the paper; the repo's convention for all four walls
RHO_C = 2.4e-9                  # t/mm^3

FY = 360.0                      # MPa
ES = 200000.0                   # MPa

# epsc0 DERIVED, never quoted (D56): Concrete02's initial compressive tangent is 2*fc/epsc0
# regardless of the grade's E, so any other value silently mismatches the calibrated modulus.
EPSC0 = 2.0 * FC / EC           # 0.002252

# Aydin's own trilinear tension-softening shape for this specimen, Table 1: a1 = 1.5, a2 = 70,
# a3 = 360. `materials.aydin_lattice_softening` solves a2/a3 from Gf instead of taking them as
# inputs, so these are the paper's values to CHECK against, not to feed in. Note a3/a2 = 5.14, which
# is what the repo's `a3_over_a2 = 5.0` default assumes.
AYDIN_A1, AYDIN_A2, AYDIN_A3 = 1.5, 70.0, 360.0

GRADES = {
    "wall": ConcreteGrade(name="wall", E=EC, nu=NU, rho=RHO_C, fc=FC, ft=FT, Gf=GF,
                          epsc0=EPSC0, fcu=0.2 * FC, epsU=8.0 * EPSC0),
}
STEEL = SteelGrade(name="S420-mesh", fy=FY, E0=ES, b=0.01)

# --- drive settings ------------------------------------------------------------------------------
# NEITHER IS MEASURED FOR THIS SPECIMEN. The 2019 paper gives a loading velocity only for the two
# sensitivity cases (Table 3: 0.25 m/s for the column, 0.025 m/s for the beam) and none for the
# walls; the 2017 test paper is not in the repo. These are the repo's cross-study constants (D62/D64)
# so this wall stays comparable with SW-NC-FF and WSH3 — damping is the knob that matters, not rate.
QUASI_STATIC_RATE = 7.6         # mm/s
DAMPING_RATIO = 0.2


def zone_of(_x: float, _y: float) -> str:
    """One concrete zone. The paper reports a single f_c and models a single rectangle."""
    return "wall"


def wall_problem(*, length: float = LW, height: float = HW) -> Problem:
    """The wall as a `Problem`: fixed base, lateral drive across the top row, no axial load.

    `N` is blank in Table 1 for this specimen — unlike the two columns (632 and 600 kN) and the
    Thomsen & Wallace wall (378 kN) — so there is no gravity case to run. The unit load below is a
    PATTERN for the displacement-controlled runners, not a force to be applied.
    """
    return Problem(
        ndm=2, ndf=2,
        domain=RectangleDomain(length=length, height=height, thickness=TW),
        material=GRADES["wall"],
        supports=(BoxSupport(box=(-EPS, length + EPS, -EPS, EPS), fix=(1, 1)),),
        loads=(BoxLoad(box=(-EPS, length + EPS, height - EPS, height + EPS), total=(1.0, 0.0)),),
    )


def bar_lines(extent: float, *, cover: float = COVER_X, spacing: float = S_BAR,
              mesh_size: float | None = None) -> list[float]:
    """Bar positions across `extent`: first at `cover`, then every `spacing` while it fits.

    With `mesh_size` given, the starting cover is snapped DOWN to a grid line if it is not already on
    one — down rather than to-nearest, because snapping up can drop the last bar off the far face and
    change the reinforcement ratio. The COUNT is the quantity that must survive.
    """
    n = int((extent - 2.0 * cover) // spacing) + 1
    if mesh_size is not None and abs(cover / mesh_size - round(cover / mesh_size)) > 1e-9:
        cover = (cover // mesh_size) * mesh_size
        n = int((extent - 2.0 * cover) // spacing) + 1
    return [cover + i * spacing for i in range(n)]


def rebars(mesh_size: float = MESH, *, length: float = LW, height: float = HW,
           full_height: bool = False) -> tuple[Rebar, ...]:
    """The orthogonal Ø8 @ 100 mesh as in-plane bar lines.

    Both curtains of a two-curtain wall (or all three layers of this one) collapse onto a single
    in-plane line per position — the standard plane-model collapse, as for WSH3 and VK3. Vertical
    bars are tagged "longitudinal" and horizontal ones "stirrup" so the visualizer separates them;
    in a wall this dense the distinction is what makes `draw.py` readable.

    `full_height` RUNS THE LONGITUDINAL BARS TO THE VERY TOP (D78, ported from the replica). Bar
    lines are placed on a cover + spacing grid, so the top line lands at y = 2150 of 2250 and leaves
    a **100 mm band of PLAIN CONCRETE** — and that band is where the drive is applied, so the whole
    base shear is introduced through unreinforced material. Every collapse measured on this wall
    tore inside it, after a peak at which damage was still correctly distributed over the base and
    middle. It is an OPTION and not the default because the paper says nothing about bar
    termination: the grid placement is an inference and so is this one, and what is measurable is
    how much the answer depends on which is chosen.
    """
    check_mesh_alignment(mesh_size, length=length, height=height)
    xs = bar_lines(length, cover=COVER_X, mesh_size=mesh_size)
    ys = bar_lines(height, cover=COVER_Y, mesh_size=mesh_size)
    top = height if full_height else ys[-1]
    bars: list[Rebar] = []
    for x in xs:
        bars.append(Rebar(path=[(x, 0.0), (x, top)], area=BAR_AREA, steel=STEEL,
                          role="longitudinal"))
    for y in ys:
        bars.append(Rebar(path=[(xs[0], y), (xs[-1], y)], area=BAR_AREA, steel=STEEL,
                          role="stirrup"))
    return tuple(bars)


def check_bar_layout(mesh_size: float = MESH, *, length: float = LW,
                     height: float = HW) -> dict[str, object]:
    """Compare the MODELLED bar layout against what Fig. 10(a) measures, and report any snap.

    Returned so `summary.py` can print it and a test can assert on it. The counts must match exactly
    — they set rho — while a position offset of one grid step is tolerated and reported.
    """
    xs = bar_lines(length, cover=COVER_X, mesh_size=mesh_size)
    ys = bar_lines(height, cover=COVER_Y, mesh_size=mesh_size)
    return {
        "n_vertical": len(xs), "n_vertical_drawn": 30,
        "n_horizontal": len(ys), "n_horizontal_drawn": 22,
        "cover_x_model": xs[0], "cover_x_drawn": COVER_X,
        "cover_y_model": ys[0], "cover_y_drawn": COVER_Y,
        "cover_y_offset": ys[0] - COVER_Y,
        "spacing_model": xs[1] - xs[0] if len(xs) > 1 else float("nan"),
        "spacing_drawn": S_BAR,
    }


def mesh_fits(mesh_size: float, *, length: float = LW, height: float = HW) -> bool:
    """Whether `mesh_size` divides this panel and its bar layout — the non-raising form."""
    try:
        check_mesh_alignment(mesh_size, length=length, height=height)
    except SystemExit:
        return False
    return True


def check_mesh_alignment(mesh_size: float, *, length: float = LW, height: float = HW) -> None:
    """Fail loudly if the grid cannot carry the geometry or the bar lines.

    A bar line that falls between grid lines silently loses its struts (`rebar_node_chain` raises)
    or, worse, snaps to the wrong row. Every dimension here is a whole number of 50 mm cells, so
    this should never fire at the default — it fires when someone passes `--mesh 30`.
    """
    dims = (("length", length), ("height", height), ("bar spacing", S_BAR), ("cover", COVER))
    bad = [(n, v) for n, v in dims if abs(v / mesh_size - round(v / mesh_size)) > 1e-9]
    if bad:
        coarsest = math.gcd(*(round(v) for _n, v in dims))
        nodes = int((length / coarsest + 1) * (height / coarsest + 1))
        raise SystemExit(
            f"mesh {mesh_size:g} does not divide the "
            + ", ".join(f"{n} ({v:g})" for n, v in bad)
            + f" of the {length:g} x {height:g} panel.\n"
              f"  The coarsest mesh that divides every dimension of THIS panel is {coarsest:g} mm "
              f"({nodes:,} nodes).\n"
              f"  The default {MESH:g} carries the Fig. 10(a) panel ({LW:g} x {HW:g}); the other "
              f"candidate panels do not\n  share its divisors, which is why `mesh_fits` guards the "
              f"geometry sweep rather than letting it fail here.")


# --- node selectors (post-build queries) ----------------------------------------------------------
# The wall is driven across its WHOLE top row, not at a point. There is no loading beam in the model
# — the paper's Fig. 4(f) shows a bare rectangle — so a point load would punch a local mechanism into
# a lattice this squat before the wall ever responded globally.


def top_nodes(model, *, height: float = HW) -> list[int]:
    """The loaded row — CONCRETE nodes only.

    Under the bond scheme a bar node is duplicated at a concrete node's coordinates, so an
    unqualified box query returns both (D80 item 4). The actuator loads the concrete and the bars
    follow through their bond links; driving the duplicate as well would pin zero slip at exactly
    the row whose load introduction is the question (D78). Without bond this is every top node, as
    before.
    """
    return select_nodes(model, (-LW, 2.0 * LW, height - EPS, height + EPS), kind="concrete")


def lateral_loads(model, *, height: float = HW) -> list[Load]:
    """Unit lateral load pattern spread evenly over the top row (the pushover reference pattern)."""
    ids = top_nodes(model, height=height)
    return [Load(nid, (1.0 / len(ids), 0.0)) for nid in ids]


def control_node(model, *, height: float = HW) -> int:
    """Displacement control point: the top-left corner, so drift is read at the loaded edge."""
    return select_nodes(model, (-EPS, EPS, height - EPS, height + EPS), kind="concrete")[0]


def drive_nodes(model, *, height: float = HW) -> list[int]:
    return top_nodes(model, height=height)


def base_nodes(model) -> list[int]:
    """The reaction set — EVERY node on the fixed row, steel duplicates included.

    Deliberately `kind="any"`, the opposite of `top_nodes`: since 2026-09-05 the builder anchors a
    steel node wherever its host concrete node is supported, so those nodes carry part of the base
    shear and dropping them would lose the bars' contribution to it.
    """
    return select_nodes(model, (-LW, 2.0 * LW, -EPS, EPS), kind="any")
