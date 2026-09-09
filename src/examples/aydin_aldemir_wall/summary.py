"""Trace every modelling number back to its source, and price the paper's own contradictions.

Same role as `katrin_wall/summary.py` and `vk3_wall/summary.py`: run it before any analysis, read
what is a measurement and what is a convention, and check the derived ratios against whatever the
source prints. Here it does one extra job — the 2019 paper contradicts itself on the panel size, and
this module is where that is ADJUDICATED rather than asserted, by computing the initial stiffness
and flexural capacity of each candidate geometry and comparing both against Table 4.

    uv run python examples/aydin_aldemir_wall/summary.py
"""

from __future__ import annotations

import math

import testdata
from build import calibrate
from specimen import (
    BAR_AREA, COVER_X, COVER_Y, EC, EPSC0, FC, FIG4F_CROP_H, FT, FY, GF, HORIZON, HW, LW, MESH, NU,
    PAPER_GRID, PAPER_GRID_T, S_BAR, TW, bar_lines, check_bar_layout,
)

K_EXP = 1038.44e3      # N/mm  (Table 4)
F_EXP = 963.592e3      # N     (Table 4)


def _rho() -> float:
    """Reinforcement ratio of the Ø8 @ 100 mesh, identical in both directions."""
    return BAR_AREA / (S_BAR * TW)


def _cantilever_stiffness(length: float, height: float,
                          thickness: float = TW) -> tuple[float, float, float]:
    """(flexural, shear, combined) tip stiffness of an uncracked gross-section cantilever, N/mm."""
    I = thickness * length ** 3 / 12.0
    G = EC / (2.0 * (1.0 + NU))
    k_f = 3.0 * EC * I / height ** 3
    k_v = G * (5.0 / 6.0) * thickness * length / height
    return k_f, k_v, 1.0 / (1.0 / k_f + 1.0 / k_v)


def _invert_paper_grid(nodes: int = 20_385, elements: int = 80_684) -> list[tuple[int, int]]:
    """Every structured grid whose node AND element counts match Table 2 exactly.

    The counts are for a horizon-1.5 lattice on an nx-by-ny cell grid: orthogonal struts along each
    row and column plus two diagonals per cell. Solving both simultaneously over the integers is a
    much stronger constraint than either alone — and it is the only quantitative statement about this
    specimen's geometry that does not depend on reading a raster figure.
    """
    out = []
    for nx in range(1, 1000):
        if nodes % (nx + 1):
            continue
        ny = nodes // (nx + 1) - 1
        if ny >= 1 and nx * (ny + 1) + (nx + 1) * ny + 2 * nx * ny == elements:
            out.append((nx, ny))
    return out


def _flexural_shear(length: float, height: float, thickness: float = TW,
                    rho: float | None = None) -> tuple[float, float]:
    """(M_n, V_flex) for a wall with uniformly distributed vertical steel and no axial load.

    Distributed-steel wall approximation: all bars beyond the neutral axis yield in tension, the
    compression block is Whitney's. Crude — no strain hardening, no compression steel, no
    confinement — so treat it as an order check on the geometry, not a capacity prediction.
    """
    rho = _rho() if rho is None else rho
    a_st = rho * thickness * length                # total vertical steel area
    t_cap = a_st * FY                              # full-yield tension resultant
    c = t_cap / (0.85 * FC * thickness * 0.85 + t_cap / length)
    m_n = 0.5 * t_cap * length * (1.0 - c / length)
    return m_n, m_n / height


def main() -> None:
    print(f"{testdata.__doc__.splitlines()[0]}\n")
    print("=" * 100)
    print("GEOMETRY — SETTLED: 3000 wide x 2250 high x 210 thick")
    print("=" * 100)
    print(f"  Adopted: {LW:.0f} x {HW:.0f} x {TW:.0f}   aspect ratio {HW / LW:.2f}")
    print("  Panel confirmed by the user 2026-08-28; THICKNESS corrected 120 -> 210 on 2026-08-29.")

    print("\n" + "-" * 100)
    print("  Why 210 and not the 120 the plan brackets (D73)")
    print("-" * 100)
    print("  120 dimensions ONE panel of a precast 'double wall'. The load-carrying section is the")
    print("  whole stack. Five independent checks, none of them a fit to the collapse:")
    print("    1. replica/ reproduces Aydin's Table 2 counts exactly (20,385 nodes / 80,684 struts)")
    print("       and then misses his K_sim = 943.16 kN/mm by 1.743. Stiffness is EXACTLY linear in")
    print("       thickness, so 120 x 1.743 = 209.1 mm.")
    print("    2. run at 210 the replica gives K = 923.7 = 0.979 of his K_sim, unfitted.")
    print("    3. 210 = 50 + 100 + 60 of the plan's own 50|140|100|60 stack (two shells + core;")
    print("       140 is the clear gap).")
    print("    4. the uncracked section finally sits ABOVE the measured stiffness, as it must.")
    print("    5. V_flex lands on the measured force.")

    print("\n" + "-" * 100)
    print("  The two checks that decide it (uncracked gross section)")
    print("-" * 100)
    print(f"  {'thickness':<26s}{'rho %':>8s}{'k_total':>10s}{'/K_exp':>9s}{'V_flex':>10s}{'/F_exp':>9s}")
    for label, t in (("120 (one panel)", 120.0), ("190", 190.0), ("210 (ADOPTED)", 210.0),
                     ("240 (two 120 shells)", 240.0), ("350 (whole plan stack)", 350.0)):
        rho = BAR_AREA / (S_BAR * t)
        _kf, _kv, k = _cantilever_stiffness(LW, HW, t)
        _m, v = _flexural_shear(LW, HW, t, rho)
        flag = "  <- impossible: below measured" if k < K_EXP else ""
        print(f"  {label:<26s}{rho * 100:>8.3f}{k / 1e3:>10,.0f}{k / K_EXP:>9.2f}"
              f"{v / 1e3:>10,.0f}{v / F_EXP:>9.2f}{flag}")
    print(f"\n  measured (Table 4): K = {K_EXP / 1e3:,.1f} kN/mm, F = {F_EXP / 1e3:,.1f} kN")
    print("  An uncracked gross section is an UPPER BOUND on a fixed-base cantilever's initial")
    print("  stiffness, so k_total must exceed the measurement. Only 210 and above do, and 240+")
    print("  overshoot the force.")

    print("\n" + "-" * 100)
    print("  Fig. 4(f) and the Table 2 inversion")
    print("-" * 100)
    print(f"  Fig. 4(f) is a DETAIL VIEW of the bottom {FIG4F_CROP_H:.0f} mm of his model (author,")
    print("  2026-09-05), so its '1500' is a CROP height and not a panel dimension. It was never a")
    print("  candidate geometry. This also dissolves the aspect anomaly recorded here until then:")
    print("  a drawn block matching none of its own labels is what a crop looks like.")
    print("  Inverting his Table 2 counts gives a unique cell grid:")
    for nx, ny in _invert_paper_grid():
        print(f"      {nx} x {ny} cells = {nx * 20:,.0f} x {ny * 20:,.0f} mm")
    print("  — matching neither printed geometry, and proving the element count is CONCRETE ONLY.")
    print("  That inversion is what made the thickness measurable, since it fixed his panel size.")
    print("\n  BUT IT IS UNIQUE ONLY UP TO TRANSPOSITION, and that is OPEN:")
    print(f"      {PAPER_GRID[0]:,.0f} x {PAPER_GRID[1]:,.0f} mm   <- what replica/ assumes "
          f"(squat, {PAPER_GRID[1] / PAPER_GRID[0]:.2f})")
    print(f"      {PAPER_GRID_T[0]:,.0f} x {PAPER_GRID_T[1]:,.0f} mm   <- equally consistent "
          f"(slender, {PAPER_GRID_T[1] / PAPER_GRID_T[0]:.2f})")
    print("  The 2680 dimensioned along the HORIZONTAL in Fig. 4(f) points at the second, since a")
    print("  full-width crop would show the model's true width. Not resolved: it needs the author,")
    print("  and it would move every published replica ratio. `elastic.py --geometry-sweep` prices it.")

    print("\n" + "=" * 100)
    print("REINFORCEMENT — '3- Ø8 bars (@100 mm)', called out in both directions in Fig. 10(a)")
    print("=" * 100)
    rho = _rho()
    print(f"  {len(bar_lines(LW, cover=COVER_X))} vertical bar lines @ {S_BAR:.0f} mm, cover "
          f"{COVER_X:.0f} mm   (Fig. 10a measures 30 lines @ 100.0, cover 52.1/46.7)")
    print(f"  {len(bar_lines(HW, cover=COVER_Y))} horizontal bar lines @ {S_BAR:.0f} mm, cover "
          f"{COVER_Y:.0f} mm   (Fig. 10a measures 22 lines @  99.8, cover 72.9/80.4)")
    print(f"  A_s = {BAR_AREA:.1f} mm^2 per position  ->  rho = {rho * 100:.3f}% each way")
    print("\n  '3-Ø8' in a 120 mm wall is unusual, so the ratio is the evidence for it. Against the")
    print("  measured peak force, with the adopted geometry:")
    for n, label in ((2, "2 bars (one curtain each face)"), (3, "3 bars (as drawn)")):
        r = n * math.pi * 8.0 ** 2 / 4.0 / (S_BAR * TW)
        a_st = r * TW * LW
        t_cap = a_st * FY
        c = t_cap / (0.85 * FC * TW * 0.85 + t_cap / LW)
        v = 0.5 * t_cap * LW * (1.0 - c / LW) / HW
        v_s = r * FY * TW * (0.8 * LW)                 # 45-degree truss, horizontal steel only
        print(f"    {label:<32s} rho = {r * 100:.3f}%   V_flex = {v / 1e3:>6,.0f} kN "
              f"({v / F_EXP:.2f} x)   V_s = {v_s / 1e3:>6,.0f} kN ({v_s / F_EXP:.2f} x)")
    print("  Both routes point the same way: 3 bars lands near the measurement, 2 falls well short.")

    print("\n  MODELLED LAYOUT vs the drawing (Fig. 10a is to scale: its two printed dimensions")
    print("  give 1.3529 and 1.3505 mm/px, agreeing to 0.17%, so pixel positions are millimetres):")
    lay = check_bar_layout(MESH)
    ok_v = lay["n_vertical"] == lay["n_vertical_drawn"]
    ok_h = lay["n_horizontal"] == lay["n_horizontal_drawn"]
    print(f"    vertical bars    model {lay['n_vertical']:2d} lines @ {lay['spacing_model']:.0f} mm, "
          f"cover {lay['cover_x_model']:.0f}    drawn {lay['n_vertical_drawn']} @ "
          f"{lay['spacing_drawn']:.0f}, cover {lay['cover_x_drawn']:.0f}   {'OK' if ok_v else 'MISMATCH'}")
    print(f"    horizontal bars  model {lay['n_horizontal']:2d} lines @ {lay['spacing_model']:.0f} mm, "
          f"cover {lay['cover_y_model']:.0f}    drawn {lay['n_horizontal_drawn']} @ "
          f"{lay['spacing_drawn']:.0f}, cover {lay['cover_y_drawn']:.0f}   {'OK' if ok_h else 'MISMATCH'}")
    if lay["cover_y_offset"]:
        print(f"    -> the horizontal cage sits {abs(lay['cover_y_offset']):.0f} mm LOW: 75 is not a "
              f"multiple of the {MESH:.0f} mm grid, and the cover")
        print(f"       is snapped DOWN so all {lay['n_horizontal']} lines survive (snapping up drops "
              f"one and moves rho by 4.5%). --mesh 25 removes the offset.")

    print("\n" + "=" * 100)
    print("MATERIALS — what is measured, what is a code value")
    print("=" * 100)
    print(f"  f_c   = {FC:>8.2f} MPa   MEASURED (Table 1)")
    print(f"  f_t   = {FT:>8.2f} MPa   CONVENTION — TS 500, Table 1 footnote a. There is no M_cr in")
    print(f"                          the 2019 paper to pin it against, unlike WSH3 (D67).")
    print(f"  E_c   = {EC:>8,.0f} MPa   CONVENTION — ACI 318 4700*sqrt(f_c), Table 1 footnote b")
    print(f"  G_f   = {GF:>8.3f} N/mm  CONVENTION — Table 1 gives 75 N/m with no source footnote")
    print(f"  nu    = {NU:>8.2f}       NOT GIVEN — the repo's constant for all four walls")
    print(f"  epsc0 = {EPSC0:>8.5f}       DERIVED as 2*f_c/E_c, never quoted (D56)")
    print(f"  f_y   = {FY:>8.0f} MPa   MEASURED (Table 1) — the lowest of the paper's six specimens")

    print("\n" + "=" * 100)
    print("STRUT LIFE — the D67 pre-flight check")
    print("=" * 100)
    eps_cr = FT / EC
    for label, L in (("orthogonal", MESH), ("diagonal", MESH * math.sqrt(2.0))):
        ets = FT * FT * L / (2.0 * GF)
        eps_ult = eps_cr + FT / ets
        print(f"  {label:<12s} L = {L:>6.1f} mm   E_ts = {ets:>9,.0f} MPa   "
              f"eps_ult/eps_cr = {eps_ult / eps_cr:>5.1f}")
    print("  For comparison: SW-NC-FF 19.8 (ran clean), VK3 14.2, WSH3 13.4 (needed --gf-factor 2).")

    print("\n" + "=" * 100)
    print("CALIBRATION — two routes, and this is the specimen that motivates asking")
    print("=" * 100)
    uni, equi = calibrate(field="uniaxial"), calibrate(field="equibiaxial")
    for cal in (uni, equi):
        print(f"  {cal.field:<12s} A_t = {cal.area:>8,.1f} mm^2   EA = {cal.EA:.4g} N   "
              f"C = {cal.EA / (EC * MESH * TW):.4f}")
    print(f"  ratio uniaxial/equibiaxial = {uni.area / equi.area:.3f}")
    print(f"  shear-stiffness error at nu = {NU}: {uni.isotropy_error * 100:.2f}% — "
          f"the one to watch at aspect ratio {HW / LW:.2f}")

    print("\n" + "=" * 100)
    print("NOT MODELLED")
    print("=" * 100)
    for item in ("the foundation (following the paper: fixed base, bare rectangle in Fig. 4f)",
                 "any axial load (Table 1 reports none)",
                 "bar buckling and fracture (Steel02 has neither) — but the test showed no "
                 "degradation at all, so this may not bind",
                 "whatever the 'double wall' section actually is: the 50|140|100|60 stack in the "
                 "Fig. 10(a) plan is modelled as one 120 mm panel",
                 "bond, UNLESS --bond is passed; the bar type is not reported either way"):
        print(f"  - {item}")


if __name__ == "__main__":
    main()
