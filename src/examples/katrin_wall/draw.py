"""Render the lattice discretization of the WSH wall (Dazio, Beyer & Bachmann 2009), for the report.

Builds the same Aydin-calibrated model the analyses use, but runs NO analysis. Two panels: the full
elevation and a zoom of the wall-foundation region where the plastic hinge forms.

Run from src/:  python examples/katrin_wall/draw.py [--mesh 50] [--horizon 1.5]
Output: examples/output/katrin_wall/katrin_wall_drawing.png
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from build import calibrate, wall_lattice  # noqa: E402
from specimen import (  # noqa: E402
    BOUNDARY_AREA, FND_H, FND_L, H_AREA, HORIZON, H_WALL_MODEL, L_V, L_V_MODEL, LW, MESH, OUT,
    PLASTIC_ZONE, SPECIMEN, UNIT, WEB_AREA,
)

C_CONC, C_FIX = "0.78", "#d62728"
FAMILIES = [
    ("bound", BOUNDARY_AREA, "#1f77b4", 2.0, "boundary 6ϕ12 (3 lines/end)"),
    ("web", WEB_AREA, "#9467bd", 1.1, "web vertical 22ϕ8 (11 lines)"),
    ("horiz", H_AREA, "#2ca02c", 0.9, "horizontal web 2ϕ6 @150"),
]
ZONES = [
    (-FND_H, 0.0, "#e8e0d0", f"foundation 2800×600 (elastic, E×{700/150:.2f} for thickness)"),
    (0.0, H_WALL_MODEL, "#dce8f0", f"wall (fc {UNIT['fc']:.1f} MPa, Ec {UNIT['Ec']/1000:.1f} GPa)"),
    (H_WALL_MODEL, L_V_MODEL, "#cdd0e0", "stiff cap ≙ tapered loading head"),
]


def classify(model):
    groups = {n: [] for n, *_ in FAMILIES}
    conc = []
    for el in model.elements:
        if len(el.nodes) != 2:
            continue
        seg = [model.nodes[el.nodes[0]].coords, model.nodes[el.nodes[1]].coords]
        if el.kind not in ("longitudinal", "stirrup"):
            conc.append(seg)
            continue
        for name, area, *_ in FAMILIES:
            if abs(el.args[0] - area) < 1e-6:
                groups[name].append(seg)
                break
    return conc, groups


def draw(ax, parts, base_xy, *, title, ylim, show_load=True):
    conc, groups = parts
    for y0, y1, fc, _lab in ZONES:
        x0 = -FND_L / 2.0 if y1 <= 0.0 else -LW / 2.0
        w = FND_L if y1 <= 0.0 else LW
        ax.add_patch(Rectangle((x0, y0), w, y1 - y0, facecolor=fc, edgecolor="none", zorder=0))
    ax.add_collection(LineCollection(conc, colors=C_CONC, linewidths=0.25, zorder=1))
    for z, (name, _a, colour, lw, _l) in enumerate(FAMILIES):
        if groups[name]:
            ax.add_collection(LineCollection(groups[name], colors=colour, linewidths=lw, zorder=2 + z))
    if len(base_xy):
        ax.scatter(base_xy[:, 0], base_xy[:, 1], marker="^", s=14, c=C_FIX, zorder=6)
    ax.axhline(PLASTIC_ZONE, color="#ff7f0e", lw=1.0, ls="-.", zorder=5)
    ax.text(-LW / 2.0 + 30, PLASTIC_ZONE + 40, f"confined plastic zone {PLASTIC_ZONE:.0f}",
            fontsize=7, color="#ff7f0e")
    if show_load:
        ax.annotate("", xy=(-LW / 2.0 - 40, L_V_MODEL), xytext=(-LW / 2.0 - 640, L_V_MODEL),
                    arrowprops=dict(arrowstyle="-|>", color=C_FIX, lw=2.0))
        ax.text(-LW / 2.0 - 660, L_V_MODEL, f"V\nLv={L_V:.0f}", color=C_FIX, fontsize=8,
                ha="right", va="center")
    ax.set_aspect("equal")
    ax.set_xlim(-FND_L / 2.0 - 750, FND_L / 2.0 + 80)
    ax.set_ylim(*ylim)
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_title(title, fontsize=9)


def main(*, mesh_size: float = MESH, horizon: float = HORIZON) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    model = wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon, nonlinear=False)
    parts = classify(model)
    base_xy = np.array([n.coords for n in model.nodes.values() if abs(n.coords[1] + FND_H) < 1e-6])
    struts, bars = len(parts[0]), sum(len(v) for v in parts[1].values())

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 8.6), gridspec_kw={"width_ratios": [1, 1]})
    draw(axes[0], parts, base_xy, ylim=(-FND_H - 120, L_V_MODEL + 120),
         title=f"Full elevation — {len(model.nodes)} nodes, {struts} concrete struts, {bars} bar struts")
    draw(axes[1], parts, base_xy, ylim=(-FND_H - 60, 2100), show_load=False,
         title="Plastic-hinge region: orthogonal + diagonal struts on the horizon lattice")

    handles = [plt.Line2D([], [], color=C_CONC, lw=1.4,
                          label=f"concrete struts (horizon {horizon}×{mesh_size:.0f} mm)")]
    handles += [plt.Line2D([], [], color=c, lw=max(lw, 1.4), label=lab)
                for n, _a, c, lw, lab in FAMILIES if parts[1][n]]
    handles += [plt.Line2D([], [], color=C_FIX, marker="^", ls="", label="clamped foundation soffit")]
    handles += [plt.Line2D([], [], color=fc, lw=7, label=lab) for _y0, _y1, fc, lab in ZONES]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle(f"{SPECIMEN} (Dazio, Beyer & Bachmann 2009) — lattice model; "
                 f"Lv/lw = {L_V/LW:.2f}, N = {UNIT['N']/1e3:.0f} kN", fontsize=11)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    path = OUT / "katrin_wall_drawing.png"
    fig.savefig(path, dpi=150)
    print(f"saved {path}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Draw the WSH wall lattice model (no analysis)")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon)
