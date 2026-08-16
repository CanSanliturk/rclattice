"""Model-vs-paper verification report for the WSH lattice model.

Reads what is ACTUALLY in the built model — areas off the `Rebar` objects, counts off the assembled
FE model — and tabulates it against Dazio, Beyer & Bachmann (2009). The point is to catch a model
that quietly disagrees with its source, so every row shows the paper's value, the model's value and
the error, and the sections that cannot be compared numerically are listed explicitly rather than
omitted.

Run from src/:  python examples/katrin_wall/summary.py [--mesh 50]
Output: examples/output/katrin_wall/model_summary.md (and the same table on stdout)
"""

from __future__ import annotations

import build
import specimen as sp

# Paper values, each with its source, so a reader can check them without opening the PDF.
PAPER = {
    "lw": (2000.0, "Sec. 2.1"), "tw": (150.0, "Sec. 2.1"),
    "Lv": (4560.0, "Sec. 2.1"), "Lv/lw": (2.28, "Sec. 2.1"),
    "H_wall": (4030.0, "Fig. 2"), "fnd_l": (2800.0, "Fig. 2"), "fnd_h": (600.0, "Fig. 2"),
    "rho_bound": (1.54, "Table 1"), "rho_web": (0.54, "Table 1"),
    "rho_tot": (0.82, "Table 1"), "rho_h": (0.25, "Table 1"),
    "A_bound": (678.58, "6phi12, Fig. 1c"), "A_web": (1105.84, "22phi8, Fig. 1c"),
    "N": (686.0, "Table 1"), "N_ratio": (0.058, "Table 1"),
    "fc": (39.2, "Table 3"), "Ec": (35.2, "Table 3"),
    "s_tie": (75.0, "Table 1"), "s_Dnom": (6.25, "Table 1"),
    "V_max": (454.0, "Table 5"), "drift_u": (2.03, "Table 4"), "mu": (5.7, "Table 4"),
}


def _row(label, paper, model, unit="", src="", tol=1.0):
    if paper is None:
        return f"| {label} | — | {model} | — | {src} |"
    err = abs(model - paper) / abs(paper) * 100 if paper else 0.0
    flag = "" if err <= tol else "  ⚠"
    return (f"| {label} | {paper:g} {unit} | {model:.4g} {unit} | {err:.2f}%{flag} | {src} |")


def report(mesh_size: float = sp.MESH) -> str:
    bars = sp.rebars(mesh_size)
    vertical = [r for r in bars if abs(r.path[-1][0] - r.path[-2][0]) < 1e-9]
    a_bound = sum(r.area for r in vertical if r.steel is sp.S12) / 2.0     # per boundary
    a_web = sum(r.area for r in vertical if r.steel is sp.S8)
    horiz = [r for r in bars if r.role == "stirrup" and abs(r.area - sp.H_AREA) < 1e-6]
    s_tie, scale = sp.tie_spacing(mesh_size)

    rho_tot = (2 * a_bound + a_web) / (sp.LW * sp.TW) * 100
    rho_h = sp.H_AREA / (sp.H_SPACING * sp.TW) * 100
    bound_len = PAPER["A_bound"][0] / (PAPER["rho_bound"][0] / 100 * sp.TW)
    rho_bound = a_bound / (bound_len * sp.TW) * 100
    rho_web = a_web / ((sp.LW - 2 * bound_len) * sp.TW) * 100

    cal = build.calibrate(mesh_size=mesh_size)
    model = build.wall_lattice(cal.area, mesh_size=mesh_size, nonlinear=False)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    k_g, share = build.cantilever_stiffness(shear_span=sp.L_V, inertia=build.gross_inertia())
    k_t, _ = build.cantilever_stiffness(shear_span=sp.L_V,
                                        inertia=build.transformed_inertia(mesh_size))
    true_c = sum(abs(x) for x in sp.BOUNDARY_X[:3]) / 3.0
    snap_c = sum(abs(sp.snap(x, mesh_size)) for x in sp.BOUNDARY_X[:3]) / 3.0

    L = [f"# {sp.SPECIMEN} lattice model — verification against Dazio, Beyer & Bachmann (2009)",
         "",
         f"*Engineering Structures* 31:1556–1571. Mesh {mesh_size:g} mm, horizon {sp.HORIZON}.",
         "", "## Geometry", "",
         "| Quantity | Paper | Model | Error | Source |", "|---|---|---|---|---|",
         _row("wall length l_w", *PAPER["lw"][:1], sp.LW, "mm", PAPER["lw"][1]),
         _row("wall thickness t_w", PAPER["tw"][0], sp.TW, "mm", PAPER["tw"][1]),
         _row("shear span L_v", PAPER["Lv"][0], sp.L_V_MODEL, "mm", PAPER["Lv"][1]),
         _row("shear span ratio", PAPER["Lv/lw"][0], sp.L_V_MODEL / sp.LW, "", PAPER["Lv/lw"][1]),
         _row("wall height", PAPER["H_wall"][0], sp.H_WALL_MODEL, "mm", PAPER["H_wall"][1]),
         _row("foundation length", PAPER["fnd_l"][0], sp.FND_L, "mm", PAPER["fnd_l"][1]),
         _row("foundation depth", PAPER["fnd_h"][0], sp.FND_H, "mm", PAPER["fnd_h"][1]),
         "", "## Reinforcement", "",
         "| Quantity | Paper | Model | Error | Source |", "|---|---|---|---|---|",
         _row("boundary steel / end", PAPER["A_bound"][0], a_bound, "mm²", PAPER["A_bound"][1]),
         _row("web steel total", PAPER["A_web"][0], a_web, "mm²", PAPER["A_web"][1]),
         _row("rho_bound", PAPER["rho_bound"][0], rho_bound, "%", PAPER["rho_bound"][1]),
         _row("rho_web", PAPER["rho_web"][0], rho_web, "%", PAPER["rho_web"][1], tol=4.0),
         _row("rho_tot", PAPER["rho_tot"][0], rho_tot, "%", PAPER["rho_tot"][1]),
         _row("rho_h", PAPER["rho_h"][0], rho_h, "%", PAPER["rho_h"][1]),
         _row("tie spacing", PAPER["s_tie"][0], s_tie, "mm", PAPER["s_tie"][1], tol=40.0),
         _row("tie steel / height", sp.TIE_AREA / PAPER["s_tie"][0],
              sp.TIE_AREA * scale / s_tie, "mm²/mm", "area rescaled to hold the ratio"),
         "", "## Materials and loading", "",
         "| Quantity | Paper | Model | Error | Source |", "|---|---|---|---|---|",
         _row("f'c", PAPER["fc"][0], sp.FC, "MPa", PAPER["fc"][1]),
         _row("E_c", PAPER["Ec"][0], sp.EC / 1000, "GPa", PAPER["Ec"][1]),
         _row("axial load N", PAPER["N"][0], sp.UNIT["N"] / 1e3, "kN", PAPER["N"][1]),
         _row("N / (A_g f'c)", PAPER["N_ratio"][0],
              sp.UNIT["N"] / (sp.LW * sp.TW * sp.FC), "", PAPER["N_ratio"][1], tol=3.0),
         "",
         "Steel (Table 2), hardening b measured as (R_m − R_p02)/(E_s(A_gt − R_p02/E_s)):", "",
         "| Bar | f_y = R_p02 | R_m | A_gt | b |", "|---|---|---|---|---|",
         f"| ϕ12 boundary | {sp.S12.fy:.1f} MPa | 725.5 | 7.69% | {sp.S12.b:.4f} |",
         f"| ϕ8 web | {sp.S8.fy:.1f} MPa | 700.2 | 7.34% | {sp.S8.b:.4f} |",
         f"| ϕ6 horiz/ties | {sp.S6.fy:.1f} MPa | 552.2 | 6.45% | {sp.S6.b:.4f} |",
         "", "## Derived elastic response (not a paper value — a consistency check)", "",
         f"- gross-section cantilever stiffness **{k_g / 1e3:.1f} kN/mm** (shear = {share * 100:.0f}% of flexibility)",
         f"- transformed-section stiffness **{k_t / 1e3:.1f} kN/mm** (I_tr/I_g = {build.transformed_inertia(mesh_size) / build.gross_inertia():.3f})",
         f"- elastic tip displacement at V_max = {PAPER['V_max'][0]:.0f} kN: "
         f"**{PAPER['V_max'][0] * 1e3 / k_t:.2f} mm** = {PAPER['V_max'][0] * 1e3 / k_t / sp.L_V * 100:.3f}% drift",
         f"- measured yield displacement 15.4 mm (0.34% drift) ⇒ effective stiffness ≈ "
         f"{PAPER['V_max'][0] * 1e3 / 15.4 / k_t:.2f} × the uncracked value, i.e. heavy cracking by yield",
         f"- measured ultimate drift {PAPER['drift_u'][0]}%, ductility {PAPER['mu'][0]}",
         "", "## Discretization", "",
         f"- mesh {mesh_size:g} mm, horizon {sp.HORIZON} → **{len(model.nodes)} nodes, {struts} concrete struts, "
         f"{len(model.elements) - struts} bar struts**",
         f"- Aydin energy-balance strut area **A_t = {cal.area:,.1f} mm²** = "
         f"{cal.area / (sp.TW * mesh_size):.4f}·(t·mesh); lattice Poisson ratio {cal.nu_consistent:.3f}",
         f"- boundary-group centroid snapped {true_c:.1f} → {snap_c:.1f} mm "
         f"(**{abs(snap_c - true_c) / true_c * 100:.2f}%** on the flexural lever arm)",
         f"- geometric rounding: wall height {PAPER['H_wall'][0]:.0f}→{sp.H_WALL_MODEL:.0f} mm, "
         f"shear span {PAPER['Lv'][0]:.0f}→{sp.L_V_MODEL:.0f} mm",
         "", "## Deliberately NOT modelled", "",
         "- **ϕ4.2@75 crossties** and the hoops' through-thickness legs — both run out of plane, so a",
         "  2D in-plane model cannot carry them. Only the hoops' two in-plane legs are represented.",
         "- **Tapered 400×700 loading head** — replaced by a full-width stiff elastic cap that",
         "  preserves the shear span exactly. It sits above the plastic zone, where the wall is elastic.",
         "- **Foundation and head reinforcement** — both blocks are modelled elastic (the foundation's",
         "  700 mm thickness folded into its modulus), so their bars carry no meaning.",
         "- **Bar buckling and bar fracture**, which ended the real test (D12 corner bar ruptured at",
         "  1.79% drift after buckling). Steel02 has neither, so the model cannot reproduce the failure",
         "  mode — only the response up to it.",
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

    p = argparse.ArgumentParser(description="WSH model-vs-paper verification report")
    p.add_argument("--mesh", type=float, default=sp.MESH)
    a = p.parse_args()
    main(a.mesh)
