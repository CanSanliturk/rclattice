"""Trace every modelling number back to its source, and price the paper's own anomaly.

Same role as `aydin_aldemir_wall/summary.py`: run it before any analysis, read what is a
measurement and what is a convention, and check the derived ratios against what the source prints.
The one thing to adjudicate here is Table 4's "initial stiffness" of 35.19 kN/mm, which exceeds
what the printed geometry and modulus can give.

    uv run python examples/thomsen_wallace_wall/summary.py [--grid rebar|uniform] [--mesh 25]
"""

from __future__ import annotations

import argparse
import math

import testdata
from build import calibrate, cantilever_stiffness, gross_inertia, strut_life, transformed_inertia
from specimen import (
    A_BE_LINE, A_HOOP_LINE, A_WEB_LINE, EC, EPSC0, FC, FT, FY, GF, GRID, HW, LW, MESH, N_AXIAL, NU,
    TW, X_BARS, X_BE_LEFT, X_WEB, Y_HOOPS, Y_WEB, bar_layout, grid_lines,
)

K_EXP = 35.19e3        # N/mm  (Table 4)
F_EXP = 163.284e3      # N     (Table 4)
K_SIM = 32.26e3        # N/mm  (Table 4, horizon 1.5)


def _flexural_shear() -> tuple[float, float]:
    """(M_n, V_flex): boundary bars yielding, web bars yielding, Whitney block, N included.

    An order check on the geometry and steel — no hardening, no confinement, no compression steel."""
    a_be = 4 * A_BE_LINE                      # one boundary element, both curtains
    a_web = 4 * A_WEB_LINE
    t_be, t_web = a_be * FY, a_web * FY
    # tension: far boundary element + web; compression block balances T + N
    comp = t_be + t_web + N_AXIAL
    a_blk = comp / (0.85 * FC * TW)
    c = a_blk / 0.85
    x_be = sum(X_BE_LEFT) / 4.0               # centroid of the boundary bars, from the edge
    x_web = sum(X_WEB) / 4.0
    m_n = (t_be * (LW - x_be - a_blk / 2.0) + t_web * (LW - x_web - a_blk / 2.0)
           + N_AXIAL * (LW / 2.0 - a_blk / 2.0))
    return m_n, m_n / HW


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", choices=("rebar", "uniform"), default=GRID)
    ap.add_argument("--mesh", type=float, default=MESH)
    a = ap.parse_args()

    print(f"{testdata.__doc__.splitlines()[0]}\n")
    print("=" * 100)
    print("GEOMETRY — Fig. 9(a): 1220 long x 3660 high x 102 thick, aspect ratio 3.00")
    print("=" * 100)
    print(f"  Every dimension is a rounded inch value (48 x 144 x 4 in). Drift = u_top / {HW:.0f}.")
    print(f"  Axial load N = {N_AXIAL / 1e3:.0f} kN = {N_AXIAL / (LW * TW * FC):.3f} A_g f_c, held constant.")

    print("\n" + "-" * 100)
    print("  The anomaly: Table 4's initial stiffness exceeds the uncracked section")
    print("-" * 100)
    i_g = gross_inertia()
    i_tr = transformed_inertia(a.grid, a.mesh)
    k_g, sh_g = cantilever_stiffness(shear_span=HW, inertia=i_g)
    k_tr, sh_tr = cantilever_stiffness(shear_span=HW, inertia=i_tr)
    print(f"  {'section':<28s}{'I (mm^4)':>14s}{'K (kN/mm)':>12s}{'/K_exp':>9s}{'shear share':>13s}")
    print(f"  {'gross concrete':<28s}{i_g:>14.4g}{k_g / 1e3:>12.2f}{k_g / K_EXP:>9.3f}{sh_g:>13.1%}")
    print(f"  {'transformed (n-1)As':<28s}{i_tr:>14.4g}{k_tr / 1e3:>12.2f}{k_tr / K_EXP:>9.3f}{sh_tr:>13.1%}")
    print(f"\n  measured K (Table 4)      {K_EXP / 1e3:>8.2f} kN/mm")
    print(f"  Aydin's K_sim (Table 4)   {K_SIM / 1e3:>8.2f} kN/mm = {K_SIM / k_tr:.3f} x the transformed section")
    print("  An uncracked fixed-base cantilever is an UPPER BOUND on a top-displacement secant, and")
    print("  the measurement sits 1.18x ABOVE it — so either Table 4's 'initial stiffness' is not a")
    print("  top-displacement secant of this geometry, or E_c / a dimension differs from what is")
    print("  printed. Aydin's OWN lattice also sits above the section (his equibiaxial closed form")
    print("  at nu = 1/3 gives 1.17x the uniaxial-field EA, D72), which is consistent with his")
    print("  K_sim/K_exp = 0.917 being two overshoots that partly cancel. The primary sources would")
    print("  settle it; they are not in the repo. JUDGE THE CALIBRATION ON K_lattice/K_continuum.")

    print("\n" + "=" * 100)
    print("REINFORCEMENT — Fig. 9(a)")
    print("=" * 100)
    print(f"  boundary elements: 8-#3 each end = 4 positions x 2 curtains at x = "
          + ", ".join(f"{x:g}" for x in X_BE_LEFT) + f" (mirrored); {A_BE_LINE:.1f} mm^2 per line")
    print(f"  web vertical:      8-#2 = 4 positions x 2 curtains at x = "
          + ", ".join(f"{x:g}" for x in X_WEB) + f" (INFERRED: 191 centred); {A_WEB_LINE:.1f} mm^2 per line")
    print(f"  web horizontal:    2-#2 @ 191 -> {len(Y_WEB)} lines from y = {Y_WEB[0]:g} (ASSUMED start); "
          f"{A_WEB_LINE:.1f} mm^2")
    print(f"  hoops:             2 x d_b 4.76 @ 76 -> {len(Y_HOOPS)} per boundary element from y = "
          f"{Y_HOOPS[0]:g} (ASSUMED start); {A_HOOP_LINE:.1f} mm^2 per hoop line")
    a_v = 8 * A_BE_LINE + 4 * A_WEB_LINE
    print(f"  total vertical steel {a_v:.0f} mm^2 -> rho_v = {a_v / (LW * TW):.3%}; "
          f"boundary element rho = {4 * A_BE_LINE / (191.0 * TW):.2%} over 191 mm")
    m_n, v_flex = _flexural_shear()
    print(f"  flexural check: M_n = {m_n / 1e6:.0f} kN.m -> V_flex = {v_flex / 1e3:.1f} kN = "
          f"{v_flex / F_EXP:.3f} x the measured peak (nominal f_y, no hardening, N included)")

    lay = bar_layout(a.grid, a.mesh)
    print(f"\n  MODELLED LAYOUT on the {a.grid} grid at {a.mesh:g} mm: {lay['n_x_lines']} x "
          f"{lay['n_y_lines']} lines, spacing x {lay['spacing_x'][0]:.1f}-{lay['spacing_x'][1]:.1f}, "
          f"y {lay['spacing_y'][0]:.1f}-{lay['spacing_y'][1]:.1f} mm")
    print(f"    bar-line offsets from the drawing: vertical {lay['x_offset_max']:.1f} mm, hoops "
          f"{lay['y_hoops_offset_max']:.1f} mm, web horizontals {lay['y_web_offset_max']:.1f} mm")

    print("\n" + "=" * 100)
    print("MATERIALS — what is measured, what is a convention")
    print("=" * 100)
    print(f"  f_c   = {FC:>8.2f} MPa   AS PRINTED (Table 1)")
    print(f"  f_t   = {FT:>8.2f} MPa   AS PRINTED (Table 1, NO footnote) — equals 0.31*sqrt(f_c)")
    print(f"  E_c   = {EC:>8,.0f} MPa   AS PRINTED (Table 1, NO footnote) — ACI would give 30,749")
    print(f"  G_f   = {GF:>8.3f} N/mm  CONVENTION — Table 1 gives 75 N/m with no source")
    print(f"  nu    = {NU:>8.2f}       NOT GIVEN — the repo's constant for every wall")
    print(f"  epsc0 = {EPSC0:>8.5f}       DERIVED as 2*f_c/E_c (D56)")
    print(f"  f_y   = {FY:>8.0f} MPa   NOMINAL Grade 60 (Table 1); measured coupons NOT in the 2019 paper")
    print(f"  b, eps_su          NOT GIVEN — the D101/D102 controls on capacity; b = 0.01 convention")

    print("\n" + "=" * 100)
    print("STRUT LIFE — the D67 pre-flight check")
    print("=" * 100)
    for label, diag in (("orthogonal", False), ("diagonal", True)):
        print(f"  {label:<12s} eps_ult/eps_cr = {strut_life(1.0, a.mesh, diagonal=diag, grid=a.grid):.1f}"
              f"   (longest strut of the family on this grid)")
    print("  For comparison: Aldemir 22.8 at mesh 50, SW-NC-FF 19.8, VK3 14.2, WSH3 13.4 (needed x2).")
    print("  A finer grid gives LONGER strut life (shorter struts soften more slowly per strain).")

    print("\n" + "=" * 100)
    print("CALIBRATION — two routes")
    print("=" * 100)
    uni, equi = calibrate(mesh_size=a.mesh, field="uniaxial"), calibrate(mesh_size=a.mesh, field="equibiaxial")
    for cal in (uni, equi):
        print(f"  {cal.field:<12s} A_t = {cal.area:>8,.1f} mm^2   EA = {cal.EA:.4g} N   "
              f"C = {cal.EA / (EC * a.mesh * TW):.4f}")
    print(f"  ratio equibiaxial/uniaxial = {equi.area / uni.area:.3f}; nu_effective = {uni.nu_effective:.3f}")
    print("  On a FLEXURAL wall the uniaxial field under-reads the continuum by (1-nu_eff^2)/(1-nu^2)")
    print(f"  = {(1 - uni.nu_effective ** 2) / (1 - NU ** 2):.3f} on the concrete alone (D53): bending is a"
          " uniaxial-STRESS state, and the")
    print("  balance pins the CONFINED modulus. Stage 0 measures both routes against the continuum.")

    print("\n" + "=" * 100)
    print("NOT MODELLED")
    print("=" * 100)
    for item in ("the foundation (fixed base, as the paper's own Fig. 4(e) model)",
                 "the loading beam (the top row is driven directly)",
                 "confinement as a separate concrete grade (the hoops are in-plane struts, D66)",
                 "bar buckling (Steel02 has none) — the likely end of the physical test",
                 "bond, unless the study's --bond axis is used; bar type not reported"):
        print(f"  - {item}")


if __name__ == "__main__":
    main()
