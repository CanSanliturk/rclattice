"""Render the lattice discretization of the SW-NC-FF flexure-controlled wall, for the report.

Builds the SAME Aydin-calibrated model the analyses use (via `build.py`) and draws it — but runs
NO analysis, so it is fast. Two panels: the full elevation, and a zoom of the wall-pedestal
interface where the flexural damage localizes.

What the drawing shows that the generic `--draw` figure does not: the three concrete casts as
shaded zones (test unit / upper loading head / pedestal), reinforcement separated by BAR TYPE
rather than lumped as "longitudinal", the fixed pedestal soffit, and the actuator level at
y = 2200 mm that sets the 2.2 shear-span ratio.

Output: examples/output/wall/wall_drawing.png. Units: N, mm.
Run from src/:  python examples/wall/draw.py [--mesh 50] [--horizon 1.5]
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from build import calibrate, nonlinear_wall_lattice  # noqa: E402
from specimen import (  # noqa: E402
    A_SHEAR, ANCHORAGE_DEPTH, BOUNDARY_AREA, HEAD_AREA, HOOK_LENGTH, HORIZON, HW, LW, MESH, OUT,
    PED_H, PED_L, SPECIMEN, TIE_AREA, UPPER_Y, WEB_H_AREA, WEB_V_AREA,
)

# zone shading: (y0, y1, facecolor, label)
ZONES = [
    (-PED_H, 0.0, "#e8e0d0", "pedestal (fc 33.7 MPa, restrained)"),
    (0.0, UPPER_Y, "#dce8f0", "test unit (fc 15.1, Ec 16.1 GPa)"),
    (UPPER_Y, HW, "#cdd0e0", "upper head (fc 43.4, Ec 26.4 GPa)"),
]
C_CONC, C_BOUND, C_WEBV, C_WEBH, C_FIX = "0.78", "#1f77b4", "#9467bd", "#2ca02c", "#d62728"
C_HOOP = "#ff7f0e"


#  bar family -> (cross-section area, colour, linewidth, legend label)
FAMILIES = [
    ("bound", BOUNDARY_AREA, "#1f77b4", 2.2, "boundary 4ϕ14 (lumped at centroid, hooked)"),
    ("webv", WEB_V_AREA, "#9467bd", 1.3, "web vertical 2ϕ12 @200 (hooked)"),
    ("webh", WEB_H_AREA, "#2ca02c", 0.9, "web horizontal 2ϕ10 @200 (closed hoop, 90° hooks)"),
    ("hoop", TIE_AREA, "#ff7f0e", 2.0, "boundary hoop 2ϕ8 (SW-IC-FF only)"),
    ("head", HEAD_AREA, "#8c564b", 1.4, "loading head 2ϕ16 @100"),
]


def classify(model):
    """Split elements into concrete struts and the reinforcement families, keyed on BAR AREA.

    Area is the reliable discriminator here: the hooked bars turn horizontal in the pedestal, so
    orientation no longer identifies them, and three separate families all carry `role="stirrup"`.
    Every family has a distinct area, and the builder stores it as the element's first argument.
    """
    groups = {name: [] for name, *_ in FAMILIES}
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
        x0 = -PED_L / 2.0 if y1 <= 0.0 else -LW / 2.0
        w = PED_L if y1 <= 0.0 else LW
        ax.add_patch(Rectangle((x0, y0), w, y1 - y0, facecolor=fc, edgecolor="none", zorder=0))

    ax.add_collection(LineCollection(conc, colors=C_CONC, linewidths=0.30, zorder=1))
    for z, (name, _area, colour, lw, _lab) in enumerate(FAMILIES):
        if groups[name]:
            ax.add_collection(LineCollection(groups[name], colors=colour, linewidths=lw, zorder=2 + z))
    if len(base_xy):
        ax.scatter(base_xy[:, 0], base_xy[:, 1], marker="^", s=18, c=C_FIX, zorder=5)

    if show_load:
        ax.annotate("", xy=(-LW / 2.0 - 30, A_SHEAR), xytext=(-LW / 2.0 - 330, A_SHEAR),
                    arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=2.0))
        ax.text(-LW / 2.0 - 340, A_SHEAR, f"V\ny={A_SHEAR:.0f}", color="#d62728", fontsize=8,
                ha="right", va="center")
        ax.axhline(A_SHEAR, color="#d62728", lw=0.6, ls="--", zorder=1)

    ax.set_aspect("equal")
    ax.set_xlim(-PED_L / 2.0 - 420, PED_L / 2.0 + 60)
    ax.set_ylim(*ylim)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title, fontsize=9)


def main(*, mesh_size: float = MESH, horizon: float = HORIZON) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon)
    parts = classify(model)
    base_xy = np.array([n.coords for n in model.nodes.values() if abs(n.coords[1] + PED_H) < 1e-6])

    struts = len(parts[0])
    bars = sum(len(v) for v in parts[1].values())
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 8.0), gridspec_kw={"width_ratios": [1, 1]})
    draw(axes[0], parts, base_xy, ylim=(-PED_H - 80, HW + 80),
         title=f"Full elevation — {len(model.nodes)} nodes, {struts} concrete struts, {bars} bar struts")
    draw(axes[1], parts, base_xy, ylim=(-PED_H - 40, 700), show_load=False,
         title=f"Wall–pedestal interface: bars hooked {ANCHORAGE_DEPTH:.0f} mm down, {HOOK_LENGTH:.0f} mm out")

    handles = [plt.Line2D([], [], color=C_CONC, lw=1.4,
                          label=f"concrete struts (horizon {horizon}×{mesh_size:.0f} mm)")]
    handles += [plt.Line2D([], [], color=c, lw=max(lw, 1.4), label=lab)
                for name, _a, c, lw, lab in FAMILIES if parts[1][name]]
    handles += [plt.Line2D([], [], color=C_FIX, marker="^", ls="", label="fixed pedestal soffit")]
    handles += [plt.Line2D([], [], color=fc, lw=7, label=lab) for _y0, _y1, fc, lab in ZONES]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(f"{SPECIMEN} flexure-controlled wall — lattice model "
                 f"(h$_w$/l$_w$=3.0, shear span {A_SHEAR:.0f} mm, N=600 kN)", fontsize=11)
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    savepath = OUT / "wall_drawing.png"
    fig.savefig(savepath, dpi=150)
    print(f"saved {savepath}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Draw the SW-NC-FF wall lattice model (no analysis)")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON, help=f"horizon (default {HORIZON})")
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon)
