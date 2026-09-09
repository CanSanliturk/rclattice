"""Render the lattice discretization of pier VK3 (Bimschas 2010, IBK Bericht 326), for `--draw`.

Builds the same Aydin-calibrated model the analyses use, but runs NO analysis. Two panels: the full
elevation and a zoom of the base region where the plastic hinge and the diagonal cracking develop.

Run from src/:  python examples/vk3_wall/draw.py [--mesh 50] [--horizon 1.5] [--nonlinear]
Output: examples/output/vk3_wall/vk3_model.png
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from build import calibrate, nonlinear_pier_lattice, pier_lattice  # noqa: E402
from specimen import (  # noqa: E402
    END_AREA, FC, FND_H, FND_L, FND_W, HEAD_Y0, HOOP_AREA, HORIZON, H_PIER, LW, L_V, MESH,
    MID_AREA, N_TOP, OUT, SPECIMEN, TW, hoop_levels,
)

C_CONC, C_FIX = "0.78", "#d62728"
ZONES = [
    (-FND_H, 0.0, "#e8e0d0", f"foundation {FND_L:.0f}x{FND_H:.0f} (elastic, "
                             f"E x {FND_W / TW:.2f} for the assumed {FND_W:.0f} mm thickness)"),
    (0.0, HEAD_Y0, "#dce8f0", f"pier (f'c {FC:.0f} MPa, phi6@200 hoops, rho_sw = 0.08%)"),
    (HEAD_Y0, H_PIER, "#cdd0e0", "load-introduction zone (phi6@75, held elastic)"),
]


def families(mesh_size: float):
    _b, _h, scale = hoop_levels(mesh_size)
    return [
        ("end", END_AREA, "#1f77b4", 2.2, "end bars 4phi14 (one line per section end)"),
        ("web", MID_AREA, "#9467bd", 1.1, f"web verticals 2phi14 @ {80:.0f} mm ({17} lines)"),
        ("hoop", HOOP_AREA, "#2ca02c", 1.0, "hoops 2phi6 @ 200 (in-plane legs) — THE shear steel"),
        ("hoop_head", HOOP_AREA * scale, "#8c8c8c", 0.8, "hoops @ 75 in the elastic head"),
    ]


def classify(model, mesh_size: float):
    fam = families(mesh_size)
    groups = {n: [] for n, *_ in fam}
    conc = []
    for el in model.elements:
        if len(el.nodes) != 2:
            continue
        seg = [model.nodes[el.nodes[0]].coords, model.nodes[el.nodes[1]].coords]
        if el.kind not in ("longitudinal", "stirrup"):
            conc.append(seg)
            continue
        for name, area, *_ in fam:
            if abs(el.args[0] - area) < 1e-6:
                groups[name].append(seg)
                break
    return conc, groups


def draw(ax, parts, base_xy, mesh_size, *, title, ylim, show_load=True):
    conc, groups = parts
    for y0, y1, fc, _lab in ZONES:
        x0 = -FND_L / 2.0 if y1 <= 0.0 else -LW / 2.0
        w = FND_L if y1 <= 0.0 else LW
        ax.add_patch(Rectangle((x0, y0), w, y1 - y0, facecolor=fc, edgecolor="none", zorder=0))
    ax.add_collection(LineCollection(conc, colors=C_CONC, linewidths=0.25, zorder=1))
    for z, (name, _a, colour, lw, _l) in enumerate(families(mesh_size)):
        if groups[name]:
            ax.add_collection(LineCollection(groups[name], colors=colour, linewidths=lw,
                                             zorder=2 + z))
    if len(base_xy):
        ax.scatter(base_xy[:, 0], base_xy[:, 1], marker="^", s=14, c=C_FIX, zorder=6)
    if show_load:
        ax.annotate("", xy=(-LW / 2.0 - 40, L_V), xytext=(-LW / 2.0 - 640, L_V),
                    arrowprops=dict(arrowstyle="-|>", color=C_FIX, lw=2.0))
        ax.text(-LW / 2.0 - 660, L_V, f"V\nLv={L_V:.0f}", color=C_FIX, fontsize=8,
                ha="right", va="center")
        # The two actuators act at DIFFERENT heights — unlike either previous wall, where the
        # lateral load and the axial load shared a level.
        ax.annotate("", xy=(0.0, H_PIER + 40), xytext=(0.0, H_PIER + 420),
                    arrowprops=dict(arrowstyle="-|>", color="#1a1a1a", lw=2.0))
        ax.text(30, H_PIER + 300, f"N = {N_TOP / 1e3:.0f} kN\nat y = {H_PIER:.0f}",
                fontsize=8, color="#1a1a1a", va="center")
    ax.set_aspect("equal")
    ax.set_xlim(-FND_L / 2.0 - 750, FND_L / 2.0 + 80)
    ax.set_ylim(*ylim)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title, fontsize=9)


def main(*, mesh_size: float = MESH, horizon: float = HORIZON, nonlinear: bool = False,
         savepath=None) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    build = nonlinear_pier_lattice if nonlinear else pier_lattice
    model = build(cal.area, mesh_size=mesh_size, horizon=horizon)
    parts = classify(model, mesh_size)
    base_xy = np.array([n.coords for n in model.nodes.values()
                        if abs(n.coords[1] + FND_H) < 1e-6])
    struts, bars = len(parts[0]), sum(len(v) for v in parts[1].values())

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 8.6))
    draw(axes[0], parts, base_xy, mesh_size, ylim=(-FND_H - 200, H_PIER + 500),
         title=f"Full elevation — {len(model.nodes)} nodes, {struts} concrete struts, "
               f"{bars} bar struts")
    draw(axes[1], parts, base_xy, mesh_size, ylim=(-FND_H - 60, 2000), show_load=False,
         title="Base region: where the plastic hinge and the diagonal cracks form")

    handles = [plt.Line2D([], [], color=C_CONC, lw=1.4,
                          label=f"concrete struts (horizon {horizon} x {mesh_size:.0f} mm)")]
    handles += [plt.Line2D([], [], color=c, lw=max(lw, 1.4), label=lab)
                for n, _a, c, lw, lab in families(mesh_size) if parts[1][n]]
    handles += [plt.Line2D([], [], color=C_FIX, marker="^", ls="",
                           label="clamped foundation soffit")]
    handles += [plt.Line2D([], [], color=fc, lw=7, label=lab) for _y0, _y1, fc, lab in ZONES]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle(f"{SPECIMEN} wall-type bridge pier — lattice analysis model "
                 f"({'nonlinear' if nonlinear else 'elastic'}); "
                 f"Lv/lw = {L_V / LW:.2f}, N = {N_TOP / 1e3:.0f} kN", fontsize=11)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    path = savepath or (OUT / "vk3_model.png")
    pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"saved model drawing to {path}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Draw the VK3 pier lattice model (no analysis)")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--nonlinear", action="store_true",
                   help="build the nonlinear model instead (same geometry; verifies it assembles)")
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon, nonlinear=a.nonlinear)
