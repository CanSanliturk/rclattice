"""Reduction, persistence and figures for VK3's two measurement chains (D63 + D68).

The pier carried two independent deformation instruments, and this module reproduces both from the
SAME `node_history` probe — which is why that probe had to learn to record more than one dof:

  1. VERTICAL CHAIN — ten LVDTs on each face plus a 150x150 mm Demec grid of 218 points, from which
     the chapter computes flexural and axial deformations. Two aligned node rows differenced column
     by column give the average vertical strain across the section: its slope is the curvature, its
     zero crossing is the neutral axis, and its departure from a straight line tests plane sections.

  2. DIAGONAL CHAIN — four string potentiometers measuring diagonal elongations "in order to
     determine the shear deformations" (Sec. 5.2.2). This is what Fig. 5.19-right is built from, and
     it needs u_x AND u_y at each panel corner, so a one-dof probe cannot express it.

THE SHEAR SPLIT IS THE POINT OF THIS STUDY. Neither previous wall package has an equivalent: the
SW-NC-FF paper reports a flexure/shear/rocking split the model could not fairly be compared against
(74% rocking), and the WSH3 paper reports curvature instead. VK3's chapter reports the split
directly, as a percentage of top displacement, and states the answer: shear is 20-22%, the highest
of the three units, and stays a roughly CONSTANT fraction through the inelastic range.

Everything computed here is also written to the run's JSON, so a figure can be redrawn without
re-running an analysis that costs hours.
"""

from __future__ import annotations

from rclattice import viz

from specimen import (
    A_SHEAR, BASE_SEGMENT, GAUGE_H, LVDT_ROWS, LW, PHI_YIELD, curvature_profile,
    deformation_components, gauge_strains, neutral_axis, panel_shear,
)
import testdata as td


def _drift(res, nh):
    return [res["disp"][i] / A_SHEAR * 100.0 for i in nh["index"]]


def reduce_history(res, gx) -> tuple[list[float], list[list[float]]]:
    """`(drift_pct, strain_profiles)` from a runner result carrying `node_history`."""
    nh = res.get("node_history")
    if not nh or not nh["values"]:
        return [], []
    return _drift(res, nh), gauge_strains(nh, gx)


def select_levels(drift_pct, targets_pct, *, both_directions: bool = True):
    """Indices of the samples closest to each target drift level, one per level and direction."""
    if not drift_pct:
        return []
    out = []
    for t in targets_pct:
        for sgn in ((1.0, -1.0) if both_directions else (1.0,)):
            goal = sgn * t
            i = min(range(len(drift_pct)), key=lambda k: abs(drift_pct[k] - goal))
            if i not in out:
                out.append(i)
    return out


def reduce_curvature(res, gx):
    """`(drift_pct, mid_heights, [[curvature per segment] per sample])`, units 1/mm."""
    nh = res.get("node_history")
    if not nh or not nh["values"]:
        return [], [], []
    mids, curv = curvature_profile(nh, gx)
    return _drift(res, nh), mids, curv


def reduce_components(res, gx):
    """`(drift_pct, top_disp_mm, [component dict per sample])` — the model's Fig. 5.19-right.

    `top_disp_mm` is the runner's OWN control displacement at each sample, kept alongside the
    decomposition so the two can be compared. That comparison is the instrument's self-check: the
    components are integrated from strains and diagonals over the gauged height only, so their sum
    should track the measured top displacement but need not equal it — the ungauged height above
    `LVDT_ROWS[-1]` is elastic and contributes a little more. A sum that DIVERGES from the control
    displacement means the reduction is wrong, not that the pier is interesting.
    """
    nh = res.get("node_history")
    if not nh or not nh["values"]:
        return [], [], []
    tops = [res["disp"][i] for i in nh["index"]]
    comps = [deformation_components(v, gx) for v in nh["values"]]
    return _drift(res, nh), tops, comps


def payload(res, gx, *, gauge: float = GAUGE_H) -> dict:
    """The gauge block to embed in a run's JSON: enough to redraw anything without re-running."""
    drift, strains = reduce_history(res, gx)
    _d, mids, curv = reduce_curvature(res, gx)
    _d2, tops, comps = reduce_components(res, gx)
    nh = res.get("node_history") or {}
    return {"x_mm": list(gx), "gauge_mm": gauge, "drift_pct": drift, "strain": strains,
            "rows_mm": list(LVDT_ROWS), "mid_mm": list(mids), "curvature": curv,
            "top_mm": tops, "components": comps,
            "shear_strain": [panel_shear(v, gx) for v in nh.get("values", [])]}


# --- vertical chain: strain profile across the section --------------------------------------------
def figure(drift_pct, strains, gx, *, indices, savepath, title, gauge: float = GAUGE_H,
           show_fit: bool = True) -> list[tuple[int, float, float, float]]:
    """Draw the selected profiles. Returns `(index, drift%, curvature, x_neutral)` per profile."""
    profiles, fits, rows = [], [], []
    for i in sorted(indices, key=lambda k: drift_pct[k]):
        fit = neutral_axis(gx, strains[i])
        profiles.append({"x": list(gx), "strain": strains[i],
                         "label": f"{drift_pct[i]:+.2f}% drift"})
        fits.append(fit if show_fit else None)
        rows.append((i, drift_pct[i], fit[0], fit[2]))
    viz.figure_strain_profile(profiles, savepath=str(savepath), gauge=gauge, width=LW,
                              title=title, fits=fits if show_fit else None)
    return rows


def print_table(rows) -> None:
    print(f"  {'drift':>8}  {'curvature (1/mm)':>18}  {'neutral axis x (mm)':>20}  "
          f"{'NA depth from compressed edge':>30}")
    for _i, d, curv, x0 in rows:
        if x0 != x0:
            na, depth = "        -- (no gradient)", "                            --"
        else:
            na = f"{x0:20.1f}"
            edge = LW / 2.0 if curv < 0.0 else -LW / 2.0
            depth = f"{abs(edge - x0):30.1f}"
        print(f"  {d:+7.2f}%  {curv:18.3e}  {na}  {depth}")


# --- curvature over height ------------------------------------------------------------------------
MODEL_C, INK, INK_2, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#d8d7d2"


def figure_curvature(drift_pct, mids, curvature, *, indices, savepath, title,
                     phi_y: float = PHI_YIELD) -> list[tuple[int, float, float]]:
    """Curvature profiles over height. Returns `(index, drift%, base-segment curvature)`."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = sorted(indices, key=lambda i: drift_pct[i])
    fig, ax = plt.subplots(figsize=(6.4, 7.2))
    cmap = plt.get_cmap("viridis")
    rows = []
    for k, i in enumerate(order):
        c = [v * 1e6 for v in curvature[i]]
        colour = cmap(k / max(1, len(order) - 1) * 0.85)
        ax.plot(c, [h / 1000.0 for h in mids], "-o", ms=3.5, lw=1.3, color=colour,
                label=f"{drift_pct[i]:+.2f}% drift")
        rows.append((i, drift_pct[i], curvature[i][BASE_SEGMENT[0]]))
    for sign in (1, -1):
        ax.axvline(sign * phi_y * 1e6, color=GRID, lw=0.9, ls=":")
    ax.axvline(0, color=GRID, lw=0.8)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.55)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("curvature [1/km]", fontsize=10, color=INK_2)
    ax.set_ylabel("height above the pier base [m]", fontsize=10, color=INK_2)
    ax.set_ylim(0.0, mids[-1] / 1000.0 * 1.05)
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    ax.text(0.02, 0.02, f"dotted: nominal $\\phi_y$ = {phi_y * 1e6:.2f} /km (Tab. 5.5); "
                        f"$\\phi_u$ = {td.PREDICTION['phi_u_per_km']:.1f} /km",
            transform=ax.transAxes, fontsize=8, color=INK_2)
    ax.legend(fontsize=8, loc="upper right", frameon=True, framealpha=0.9, edgecolor=GRID)
    fig.tight_layout()
    fig.savefig(str(savepath), dpi=180, facecolor="white")
    plt.close(fig)
    return rows


def print_curvature_table(rows) -> None:
    phi_y, phi_u = td.PREDICTION["phi_y_per_km"], td.PREDICTION["phi_u_per_km"]
    print(f"  {'drift':>8}  {'base phi (1/km)':>16}  {'phi/phi_y':>10}  {'phi/phi_u':>10}")
    for _i, d, phi in rows:
        p = phi * 1e6
        print(f"  {d:+7.2f}%  {p:16.2f}  {p / phi_y:10.2f}  {p / phi_u:10.2f}")


# --- diagonal chain: the shear / flexure split, i.e. Fig. 5.19-right ------------------------------
def figure_components(drift_pct, tops, comps, *, savepath, title, indices=None):
    """Cumulative relative displacement vs top displacement — the chapter's own presentation.

    Fig. 5.19-right stacks the components as PERCENTAGES of the measured top displacement and plots
    them against that displacement, so the reader sees at a glance that the shear share is roughly
    constant through the inelastic range. This draws the same axes, with the chapter's measured
    bands shaded behind the model's curves.

    Percentages are taken against the runner's own control displacement, not against the sum of the
    components, so a shortfall shows up as the curves failing to reach 100% rather than being
    normalized away.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keep = [i for i in (indices if indices is not None else range(len(tops)))
            if abs(tops[i]) > 0.2]
    keep.sort(key=lambda i: abs(tops[i]))
    u = [abs(tops[i]) for i in keep]
    sgn = [1.0 if tops[i] >= 0 else -1.0 for i in keep]
    sh = [comps[i]["shear"] * s / abs(tops[i]) * 100.0 for i, s in zip(keep, sgn)]
    bc = [comps[i]["flexure_base_crack"] * s / abs(tops[i]) * 100.0 for i, s in zip(keep, sgn)]
    fa = [comps[i]["flexure_above"] * s / abs(tops[i]) * 100.0 for i, s in zip(keep, sgn)]

    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    lo, hi = td.COMPONENTS["shear_pct_range"]
    ax.axhspan(lo, hi, color="#2ca02c", alpha=0.13, zorder=0)
    ax.text(0.985, (lo + hi) / 2.0, f" measured shear {lo:.0f}-{hi:.0f}%", transform=
            ax.get_yaxis_transform(), ha="right", va="center", fontsize=8, color="#2ca02c")
    lo_b, hi_b = td.COMPONENTS["base_crack_pct_range"]
    ax.axhspan(lo_b, hi_b, color="#ff7f0e", alpha=0.10, zorder=0)

    cum1 = sh
    cum2 = [a + b for a, b in zip(sh, bc)]
    cum3 = [a + b for a, b in zip(cum2, fa)]
    ax.fill_between(u, 0, cum1, color="#2ca02c", alpha=0.55, label="shear (diagonal gauges)")
    ax.fill_between(u, cum1, cum2, color="#ff7f0e", alpha=0.55, label="flexure in the base crack")
    ax.fill_between(u, cum2, cum3, color="#2a78d6", alpha=0.45, label="flexure above the base crack")
    ax.plot(u, cum3, color=INK, lw=1.2, label="sum of the modelled components")
    ax.axhline(100.0, color=GRID, lw=1.0, ls="--")

    ax.set_xlabel("top displacement [mm]", fontsize=10, color=INK_2)
    ax.set_ylabel("cumulative relative displacement [%]", fontsize=10, color=INK_2)
    ax.set_ylim(-5, 125)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.55)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    ax.text(0.01, -0.16, "base sliding is not drawn: the model has no construction joint, so its "
                         "sliding is zero by construction (the test measured it as negligible)",
            transform=ax.transAxes, fontsize=7.5, color=INK_2)
    ax.legend(fontsize=8, loc="lower right", frameon=True, framealpha=0.9, edgecolor=GRID)
    fig.tight_layout()
    fig.savefig(str(savepath), dpi=180, facecolor="white")
    plt.close(fig)
    return [(u[k], sh[k], bc[k], fa[k], cum3[k]) for k in range(len(u))]


def print_components_table(rows, *, every: int = 1) -> None:
    lo, hi = td.COMPONENTS["shear_pct_range"]
    print(f"  {'|u| (mm)':>9}  {'drift':>7}  {'shear %':>8}  {'base crack %':>13}  "
          f"{'flexure above %':>16}  {'sum %':>7}   test shear {lo:.0f}-{hi:.0f}%")
    for k in range(0, len(rows), max(1, every)):
        u, sh, bc, fa, tot = rows[k]
        mark = "  <-- in band" if lo <= sh <= hi else ""
        print(f"  {u:9.2f}  {u / A_SHEAR * 100:6.2f}%  {sh:8.1f}  {bc:13.1f}  {fa:16.1f}  "
              f"{tot:7.1f}{mark}")
