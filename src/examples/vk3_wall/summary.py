"""Model-vs-chapter verification report for the VK3 lattice model.

Source: Bimschas, M. (2010), IBK Bericht Nr. 326, ETH Zuerich, Chapter 5.

Reads what is ACTUALLY in the built model — areas off the `Rebar` objects, counts off the assembled
FE model — and tabulates it against Chapter 5 of the source thesis. The point is to catch a model
that quietly disagrees with its source, so every row shows the chapter's value, the model's value
and the error, and the parts that cannot be compared numerically are listed explicitly rather than
omitted.

Run from src/:  python examples/vk3_wall/summary.py [--mesh 50]
Output: examples/output/vk3_wall/model_summary.md (and the same table on stdout)
"""

from __future__ import annotations

import math

import build
import specimen as sp
import testdata as td

# Chapter values with their source, so a reader can check them without opening the PDF.
PAPER = {
    "lw": (1500.0, "Tab. 5.1"), "bw": (350.0, "Tab. 5.1"), "Lv": (3300.0, "Tab. 5.1"),
    "Lv/lw": (2.2, "Tab. 5.1"), "H_pier": (3700.0, "Fig. 5.1"),
    "fnd_l": (3000.0, "Fig. 5.1"), "fnd_h": (900.0, "Fig. 5.1"),
    "n_bars": (42.0, "Tab. 5.3"), "rho_sl": (1.23, "Tab. 5.3"), "rho_sw": (0.08, "Tab. 5.3"),
    "hoop_s": (200.0, "Tab. 5.3"),
    "fc": (34.0, "Tab. 5.4"), "fy_l": (515.0, "Tab. 5.4"), "fy_w": (518.0, "Tab. 5.4"),
    "N_top": (1300.0, "Tab. 5.7"), "N_base": (1370.0, "Tab. 5.7"), "nu_base": (0.077, "Tab. 5.7"),
    "Fn": (851.0, "Tab. 5.6"), "Fy_prime": (644.0, "Tab. 5.6"),
    "phi_y": (3.36, "Tab. 5.5"), "phi_u": (14.5, "Tab. 5.5"),
    "k0_meas": (60.0, "Tab. 5.10"), "drift_u": (1.59, "Tab. 5.9"), "mu": (3.7, "Tab. 5.9"),
}


def _row(label, paper, model, unit="", src="", tol=1.0):
    if paper is None:
        return f"| {label} | — | {model} | — | {src} |"
    err = abs(model - paper) / abs(paper) * 100 if paper else 0.0
    flag = "" if err <= tol else "  ⚠"
    return f"| {label} | {paper:g} {unit} | {model:.4g} {unit} | {err:.2f}%{flag} | {src} |"


def report(mesh_size: float = sp.MESH) -> str:
    bars = sp.rebars(mesh_size)
    vertical = [r for r in bars if abs(r.path[-1][0] - r.path[-2][0]) < 1e-9]
    n_bars = (sum(1 for r in vertical if abs(r.area - sp.MID_AREA) < 1e-6) * 2
              + sum(1 for r in vertical if abs(r.area - sp.END_AREA) < 1e-6) * 4)
    a_long = sum(r.area for r in vertical)
    rho_sl = a_long / (sp.LW * sp.TW) * 100
    base_lv, head_lv, scale = sp.hoop_levels(mesh_size)
    s_base = base_lv[1] - base_lv[0]
    rho_sw = sp.HOOP_AREA / (s_base * sp.TW) * 100

    cal = build.calibrate(mesh_size=mesh_size)
    model = build.pier_lattice(cal.area, mesh_size=mesh_size)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    i_g, i_tr = build.gross_inertia(), build.transformed_inertia(mesh_size)
    k_g, share = build.cantilever_stiffness(shear_span=sp.L_V, inertia=i_g)
    k_t, _ = build.cantilever_stiffness(shear_span=sp.L_V, inertia=i_tr)
    m_cr = (sp.FT + sp.N_BASE / (sp.LW * sp.TW)) * sp.TW * sp.LW ** 2 / 6.0
    snapped = sp.snap(sp.END_X, mesh_size)

    L = [f"# {sp.SPECIMEN} lattice model — verification against Chapter 5",
         "",
         "Bimschas, M. (2010), *Displacement Based Seismic Assessment of Existing Bridges in "
         "Regions of Moderate Seismicity*, IBK Bericht Nr. 326, ETH Zürich; vdf Hochschulverlag. "
         "[doi:10.3929/ethz-a-006237119](https://doi.org/10.3929/ethz-a-006237119), Chapter 5.",
         "",
         f"Wall-type bridge pier VK3. Mesh {mesh_size:g} mm, horizon {sp.HORIZON}.",
         "", "## Geometry", "",
         "| Quantity | Chapter | Model | Error | Source |", "|---|---|---|---|---|",
         _row("section depth l_w", *PAPER["lw"][:1], sp.LW, "mm", PAPER["lw"][1]),
         _row("section width b_w", PAPER["bw"][0], sp.TW, "mm", PAPER["bw"][1]),
         _row("shear span L_v", PAPER["Lv"][0], sp.L_V, "mm", PAPER["Lv"][1]),
         _row("aspect ratio L_v/l_w", PAPER["Lv/lw"][0], sp.L_V / sp.LW, "", PAPER["Lv/lw"][1]),
         _row("pier height", PAPER["H_pier"][0], sp.H_PIER, "mm", PAPER["H_pier"][1]),
         _row("foundation length", PAPER["fnd_l"][0], sp.FND_L, "mm", PAPER["fnd_l"][1]),
         _row("foundation depth", PAPER["fnd_h"][0], sp.FND_H, "mm", PAPER["fnd_h"][1]),
         "",
         "Every model dimension is a whole number of 50 mm cells, so — unlike the WSH3 model — the "
         "shear span needs **no rounding** and drift needs no correction.",
         "", "## Reinforcement", "",
         "| Quantity | Chapter | Model | Error | Source |", "|---|---|---|---|---|",
         _row("longitudinal bar count", PAPER["n_bars"][0], n_bars, "ϕ14", PAPER["n_bars"][1]),
         _row("ρ_sl", PAPER["rho_sl"][0], rho_sl, "%", PAPER["rho_sl"][1]),
         _row("hoop spacing (base zone)", PAPER["hoop_s"][0], s_base, "mm", PAPER["hoop_s"][1]),
         _row("ρ_sw", PAPER["rho_sw"][0], rho_sw, "%", PAPER["rho_sw"][1], tol=2.0),
         "",
         f"In-plane bar lines: {len(sp.WEB_X)} web positions of {sp.MID_AREA:.1f} mm² "
         f"(2ϕ14 each) + 2 end positions of {sp.END_AREA:.1f} mm² (4ϕ14 each, the four bars "
         f"stacked through the 284 mm thickness collapsing onto one in-plane point).",
         "", "## Materials and loading", "",
         "| Quantity | Chapter | Model | Error | Source |", "|---|---|---|---|---|",
         _row("f'_c", PAPER["fc"][0], sp.FC, "MPa", PAPER["fc"][1]),
         _row("f_y longitudinal", PAPER["fy_l"][0], sp.S_LONG.fy, "MPa", PAPER["fy_l"][1]),
         _row("f_y hoops", PAPER["fy_w"][0], sp.S_HOOP.fy, "MPa", PAPER["fy_w"][1]),
         _row("axial load N_top", PAPER["N_top"][0], sp.N_TOP / 1e3, "kN", PAPER["N_top"][1]),
         _row("axial ratio at base", PAPER["nu_base"][0],
              sp.N_BASE / (sp.LW * sp.TW * sp.FC), "", PAPER["nu_base"][1], tol=2.0),
         "",
         "**Three material quantities are CONVENTIONS, not measurements** — the chapter gives f'_c "
         "and nothing else about the concrete, and never prints a cracking moment to pin f_t "
         "against (WSH3's paper did, which is how D67 settled that convention):",
         "",
         "| Quantity | Model | Basis |", "|---|---|---|",
         f"| E_c | {sp.EC:,.0f} MPa | SIA 262, 10000·f_cm^(1/3). The specimen is Swiss and the "
         f"chapter's frame of reference is SIA 262; EC2/MC90 gives "
         f"{22000 * (sp.FC / 10) ** 0.3:,.0f} MPa, "
         f"{abs(22000 * (sp.FC / 10) ** 0.3 / sp.EC - 1) * 100:.1f}% away |",
         f"| f_t | {sp.FT:.3f} MPa | EC2 0.30·f_ck^(2/3) at f_ck = f_cm − 8 = {sp.FCK:.0f} MPa — "
         f"the D67 convention. Implies M_cr = {m_cr / 1e6:.0f} kN·m, V_cr = "
         f"{m_cr / sp.L_V / 1e3:.0f} kN, i.e. cracking at {m_cr / sp.L_V / k_t / sp.L_V:.3%} drift |",
         f"| G_f | {sp.GF_C:.4f} N/mm | MC90 on the MEAN strength, 16 mm aggregate assumed (the "
         f"chapter does not state aggregate size) |",
         f"| ε_c0 | {sp.EPSC0:.5f} | derived as 2f_c/E_c, because Concrete02's compressive tangent "
         f"is 2f_c/ε_c0 whatever E says (D56) |",
         "",
         "Steel (Tab. 5.4), hardening b measured as (f_u − f_y)/(E_s(ε_su − f_y/E_s)):", "",
         "| Bar | f_y | f_u | ε_su | b |", "|---|---|---|---|---|",
         f"| ϕ14 longitudinal | {sp.S_LONG.fy:.0f} MPa | 630 | 12.6% | {sp.S_LONG.b:.4f} |",
         f"| ϕ6 hoops | {sp.S_HOOP.fy:.0f} MPa | 681 | 8.4% | {sp.S_HOOP.b:.4f} |",
         "",
         "The longitudinal steel has a **pronounced yield plateau to ε_s ≈ 0.025** (Sec. 5.2.3a). "
         "`Steel02`'s smooth curve smears the plateau and the hardening into one slope; the chapter "
         "itself neglected hardening entirely for its moment-curvature work and called that "
         "justified.",
         "", "## Source inconsistencies (recorded, not silently resolved)", "",
         f"- **f_s,l**: Sec. 5.2.1c's text says {td.FS_L_TEXT_MPa:.0f} MPa, Tab. 5.4 says "
         f"{PAPER['fy_l'][0]:.0f}. The table is used; b changes by "
         f"{abs(sp.S_LONG.b - (630 - 520) / (200000 * (0.126 - 520 / 200000))):.4f}, which is "
         f"immaterial.",
         "- **f'_c for VK1/VK2**: the text gives 39/35 MPa, Tab. 5.4 gives 35/39 — swapped. VK3 is "
         "34 MPa in both, so this specimen is unaffected; it is a caution about the chapter.",
         "", "## Derived elastic response (a consistency check, not a chapter value)", "",
         f"- gross-section cantilever stiffness **{k_g / 1e3:.1f} kN/mm** "
         f"(shear = {share * 100:.1f}% of the flexibility at L_v/l_w = {sp.L_V / sp.LW:.2f})",
         f"- transformed-section stiffness **{k_t / 1e3:.1f} kN/mm** (I_tr/I_g = {i_tr / i_g:.3f}, "
         f"n = E_s/E_c = {sp.S_LONG.E0 / sp.EC:.1f})",
         f"- the chapter's measured k₀ = {PAPER['k0_meas'][0]:.0f} kN/mm is "
         f"**{PAPER['k0_meas'][0] * 1e3 / k_t:.2f}× the uncracked value** — it is a SECANT TO FIRST "
         f"YIELD on an already cracked pier, *not* an elastic stiffness, so it is not a valid "
         f"target for an elastic model",
         f"- the digitized loops supply one that is: **≈174 kN/mm** measured at the first elastic "
         f"level (±0.86 mm, ~151 kN), which is below the model's own V_cr = "
         f"{m_cr / sp.L_V / 1e3:.0f} kN and therefore genuinely uncracked",
         f"- strut life ε_ult/ε_cr = **{build.strut_life(1.0, mesh_size):.1f}** at G_f×1 "
         f"({build.strut_life(2.0, mesh_size):.1f} at ×2). SW-NC-FF ran clean at 19.8; WSH3 "
         f"diverged at 10.1 and needed ×2 to reach 13.4 (D67), so budget tension stiffening.",
         "", "## Discretization", "",
         f"- mesh {mesh_size:g} mm, horizon {sp.HORIZON} → **{len(model.nodes)} nodes, "
         f"{struts} concrete struts, {len(model.elements) - struts} bar struts**",
         f"- Aydin energy-balance strut area **A_t = {cal.area:,.1f} mm²** = "
         f"{cal.area / (sp.TW * mesh_size):.4f}·(t·mesh)",
         f"- shear-stiffness error at ν = {cal.nu:.2f}: **{cal.isotropy_error * 100:.2f}%**. "
         f"`nu_consistent` = {cal.nu_consistent:.3f} is where the normal and shear calibration "
         f"routes agree — **not** a Poisson ratio; the lattice's own is "
         f"{cal.nu_effective:.3f} with cubic anisotropy {cal.cubic_anisotropy:.2f} (D53). This "
         f"matters more here than in either previous study, because VK3 carries 20–22% of its "
         f"displacement in shear.",
         f"- end bar snapped {sp.END_X:.1f} → {snapped:.1f} mm "
         f"(**{abs(snapped - sp.END_X) / sp.END_X * 100:.2f}%** on the flexural lever arm); "
         f"all 19 in-plane positions stay distinct (mesh 100 collapses four of them)",
         f"- head-zone hoops snapped {sp.HOOP_HEAD[2]:.0f} → "
         f"{sp.hoop_levels(mesh_size)[1][1] - sp.hoop_levels(mesh_size)[1][0]:.0f} mm with the "
         f"area rescaled ×{scale:.3f}, holding the steel per unit height exactly. The BASE zone's "
         f"200 mm spacing is exact — and that is the zone that carries the shear.",
         "", "## Deliberately NOT modelled", "",
         "- **The hoops' 90° hooks can never open.** Sec. 5.2.1b: the hoops were closed by simple "
         "90° hooks of 100 mm \"which could open after spalling of the cover concrete\". A "
         "perfectly-bonded rebar strut never loses anchorage, so the model's transverse steel keeps "
         "working forever. Since the test's failure was precisely the moment the shear had to hand "
         "over to the hoops and they could not deliver, **expect the model to over-retain strength "
         "after ~1.27% drift.** This is the single most consequential omission.",
         "- **Bar buckling and cover spalling.** `Steel02` has neither, and both were observed from "
         "0.95% drift onward. They are what destroyed the base compression zone, which the chapter "
         "identifies as the trigger of the shear failure (Sec. 5.4).",
         "- **Sliding on an open diagonal crack.** The failure kinematic is four wedges sliding on "
         "crossed cracks; a shared-node truss lattice has no aggregate-interlock or dowel element, "
         "so it can show the diagonal load path shedding but not the sliding itself.",
         "- **The out-of-plane hoop legs.** A hoop's two legs through the 350 mm thickness are out "
         "of plane. Only the two in-plane legs are represented — which is exactly the steel ρ_sw "
         "counts, so the shear reinforcement is right, but the 3D confinement is absent.",
         "- **The base construction joint.** The chapter notes the base crack formed at 0.25F_y' "
         "from \"decompression and the absence of tension strength in the construction joint\". The "
         "model has no joint: pier and foundation share nodes, so its base crack forms only when "
         "f_t is reached, and its base sliding is identically zero (the test measured sliding as "
         "negligible, so the second half costs nothing).",
         "- **Foundation and head reinforcement.** Both blocks are elastic; the foundation's "
         f"assumed {sp.FND_W:.0f} mm out-of-plane thickness is folded into its modulus "
         f"(×{sp.FND_W / sp.TW:.2f}) and is priced by `elastic.py --foundation-sensitivity`.",
         "- **Strain penetration into the foundation.** Bars are anchored "
         f"{sp.ANCHORAGE:.0f} mm into an elastic block with perfect bond, so the fixed-end rotation "
         "the chapter separates out as \"flexure in the base crack\" (10–20% of the top "
         "displacement) is only partly representable.",
         "", "## What the digitized loops add", "",
         "VK3's peak base shear is **not tabulated anywhere in the chapter**, so without "
         "`digitize.py` there is barely a backbone to compare against. Because Fig. 5.13's panels "
         "are vector paths in drawing order rather than a raster, the digitization recovers:",
         "",
         "- the complete ordered hysteresis (3,184 points) → loop shape, energy dissipation and "
         "residual drift all become fair comparisons;",
         "- peak shear **+891 / −876 kN** at ~0.95% drift, against the chapter's predicted F_n = "
         f"{PAPER['Fn'][0]:.0f} kN;",
         "- the **measured elastic amplitudes** (0.85 / 2.10 / 5.08 mm), which the force-controlled "
         "elastic block never published and which k₀ mispredicts by up to 3×;",
         "- the **small intermediate cycles**, which the chapter defines for VK3 only as \"the top "
         "displacements measured during the corresponding cycles of VK1\" and never prints.",
         ""]
    return "\n".join(L)


def main(mesh_size: float = sp.MESH) -> None:
    sp.OUT.mkdir(parents=True, exist_ok=True)
    text = report(mesh_size)
    path = sp.OUT / "model_summary.md"
    path.write_text(text)
    print(text)
    print(f"\nsaved {path}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="VK3 model-vs-chapter verification report")
    p.add_argument("--mesh", type=float, default=sp.MESH)
    a = p.parse_args()
    main(a.mesh)
