"""Render the RC lattice discretization of the cantilever column for the report.

Builds the SAME calibrated nonlinear RC lattice used by the pushover/dynamic studies (via the
column package's build.py) and draws it: concrete horizon struts, longitudinal rebar struts,
transverse stirrup struts, and the fixed base nodes. Two panels — full elevation and a zoom of
the base region. Output: examples/doc/figures/column_lattice_mesh.png. Units: kip, in.

Run from src/:  python examples/doc/make_model_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

COLUMN = Path(__file__).resolve().parent.parent / "column"
sys.path.insert(0, str(COLUMN))

from build import beamcolumn_reference, calibrate_area, rc_lattice  # noqa: E402
from specimen import GF, GFC, H, W  # noqa: E402

OUT = Path(__file__).resolve().parent / "figures" / "column_lattice_mesh.png"


def _classify(model):
    """Split struts into concrete / longitudinal-rebar / stirrup-rebar segment lists."""
    steel = {m.id for m in model.uniaxial_materials if m.mtype == "Steel02"}
    concrete, longit, stirrup = [], [], []
    for el in model.elements:
        (xi, yi) = model.nodes[el.nodes[0]].coords
        (xj, yj) = model.nodes[el.nodes[1]].coords
        seg = [(xi, yi), (xj, yj)]
        if el.args[1] in steel:
            (longit if abs(xj - xi) < 1e-6 else stirrup).append(seg)
        else:
            concrete.append(seg)
    return concrete, longit, stirrup


def _draw(ax, concrete, longit, stirrup, base_xy, *, title, ylim=None):
    ax.add_collection(LineCollection(concrete, colors="0.78", linewidths=0.35, zorder=1))
    ax.add_collection(LineCollection(stirrup, colors="#2ca02c", linewidths=0.9, zorder=2))
    ax.add_collection(LineCollection(longit, colors="#1f77b4", linewidths=1.8, zorder=3))
    if len(base_xy):
        ax.scatter(base_xy[:, 0], base_xy[:, 1], marker="^", s=26, c="#d62728", zorder=4)
    ax.set_aspect("equal")
    ax.set_xlim(-W / 2 - 2, W / 2 + 2)
    ax.set_ylim(*(ylim if ylim else (-3, H + 3)))
    ax.set_xlabel("x (in)")
    ax.set_ylabel("y (in)")
    ax.set_title(title, fontsize=9)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    bc = beamcolumn_reference()
    k = bc["shear"][1] / bc["disp"][1]
    area, _ctrl, base = calibrate_area(k)
    model = rc_lattice(True, GF, GFC, area)

    concrete, longit, stirrup = _classify(model)
    base_xy = np.array([model.nodes[n].coords for n in base])

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 7.4), gridspec_kw={"width_ratios": [1, 1.25]})
    _draw(axes[0], concrete, longit, stirrup, base_xy,
          title=f"Full elevation — {len(model.nodes)} nodes, {len(model.elements)} struts")
    _draw(axes[1], concrete, longit, stirrup, base_xy,
          title="Base region (0–18 in): orthogonal + diagonal struts", ylim=(-2, 18))

    handles = [
        plt.Line2D([], [], color="0.78", lw=1.2, label="concrete struts"),
        plt.Line2D([], [], color="#1f77b4", lw=2.0, label="longitudinal rebar"),
        plt.Line2D([], [], color="#2ca02c", lw=1.4, label="stirrup ties"),
        plt.Line2D([], [], color="#d62728", marker="^", ls="", label="fixed base nodes"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, frameon=False)
    fig.suptitle("RC lattice discretization of the cantilever column (kip, in)", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(OUT, dpi=140)
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
