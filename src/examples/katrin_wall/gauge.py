"""Vertical strain gauge over the wall base — reduction, persistence and figure (D63).

The measurement: two aligned rows of nodes, `GAUGE_H` apart over the wall base, differenced column
by column and divided by the gauge length, giving the average vertical strain of each column across
the wall width. That is the same construction as the chain of vertical LVDTs on the tested wall's
face, from which the paper derives its curvature profiles (Fig. 11c) and base curvature (Fig. 15d),
so the model quantity and the measured one are the same thing rather than merely analogous.

The window is the PAPER'S, not the base — y = 50 to 350 mm, clear of the lowest LVDT the paper
itself discards to strain penetration, and bracketing the h = 60-360 mm window its Demec strains are
averaged over (see `specimen.GAUGE_Y0`). This is the one place the measurement differs from the
SW-NC-FF study's, which gauges from the wall-pedestal face upward.

Why this profile: a flexure-controlled wall does its work in the base region, and this is what
characterises it — its slope is the curvature, its zero crossing is the neutral axis, and its
departure from a straight line is the direct test of "plane sections remain plane", which a lattice
is under no obligation to satisfy and which the paper warns (Sec. 5.1) the real wall violates near
the base wherever inclined flexure-shear cracks form.

Everything here works off a runner's `node_history` result, and everything it computes is also
written to the run's JSON, so a figure can be redrawn without re-running the analysis — the whole
point, given the analysis costs hours.
"""

from __future__ import annotations

from rclattice import viz

from specimen import (
    A_SHEAR, GAUGE_H, LVDT_ROWS, LW, PHI_YIELD, base_curvature, curvature_profile, gauge_strains,
    neutral_axis,
)


def reduce_history(res, gx) -> tuple[list[float], list[list[float]]]:
    """`(drift_pct, strain_profiles)` from a runner result carrying `node_history`.

    Each sample's drift is looked up through the probe's `index`, which points into the run's own
    `disp` array — so a profile is always tied to the drift it was actually taken at, whatever the
    sampling stride was.
    """
    nh = res.get("node_history")
    if not nh or not nh["values"]:
        return [], []
    disp = res["disp"]
    drift = [disp[i] / A_SHEAR * 100.0 for i in nh["index"]]
    return drift, gauge_strains(nh, gx)


def select_levels(drift_pct, targets_pct, *, both_directions: bool = True):
    """Indices of the samples closest to each target drift level, one per level and direction.

    Picks by nearest drift rather than by protocol bookkeeping, so it works the same on a monotonic
    pushover, on a cyclic run, and on a run that stopped early — a level the run never reached
    simply resolves to its nearest neighbour, which `label_for` then reports honestly.
    """
    if not drift_pct:
        return []
    out = []
    for t in targets_pct:
        signs = (1.0, -1.0) if both_directions else (1.0,)
        for sgn in signs:
            goal = sgn * t
            i = min(range(len(drift_pct)), key=lambda k: abs(drift_pct[k] - goal))
            if i not in out:
                out.append(i)
    return out


def reduce_curvature(res, gx):
    """`(drift_pct, mid_heights, [[curvature per segment] per sample])` from a runner result.

    The model's counterpart to the paper's Fig. 11c: one curvature per LVDT segment, at every
    sample the probe recorded.
    """
    nh = res.get("node_history")
    if not nh or not nh["values"]:
        return [], [], []
    disp = res["disp"]
    drift = [disp[i] / A_SHEAR * 100.0 for i in nh["index"]]
    mids, curv = curvature_profile(nh, gx)
    return drift, mids, curv


def payload(res, gx, *, gauge: float = GAUGE_H) -> dict:
    """The gauge block to embed in a run's JSON: enough to redraw anything without re-running.

    Stores the two DERIVED products (base strain profile, curvature profile) rather than the raw
    nodal displacements. Ten rows of 41 columns would make the raw record ~8x larger, and nothing
    downstream needs anything the two products do not already carry.
    """
    drift, strains = reduce_history(res, gx)
    _d, mids, curv = reduce_curvature(res, gx)
    return {"x_mm": list(gx), "gauge_mm": gauge, "drift_pct": drift, "strain": strains,
            "rows_mm": list(LVDT_ROWS), "mid_mm": list(mids), "curvature": curv}


def figure(drift_pct, strains, gx, *, indices, savepath, title, gauge: float = GAUGE_H,
           show_fit: bool = True) -> list[tuple[int, float, float, float]]:
    """Draw the selected profiles. Returns `(index, drift%, curvature, x_neutral)` per profile."""
    profiles, fits, rows = [], [], []
    order = sorted(indices, key=lambda i: drift_pct[i])
    for i in order:
        fit = neutral_axis(gx, strains[i])
        profiles.append({"x": list(gx), "strain": strains[i],
                         "label": f"{drift_pct[i]:+.2f}% drift"})
        fits.append(fit if show_fit else None)
        rows.append((i, drift_pct[i], fit[0], fit[2]))
    viz.figure_strain_profile(profiles, savepath=str(savepath), gauge=gauge, width=LW,
                              title=title, fits=fits if show_fit else None)
    return rows


def print_table(rows) -> None:
    """Curvature and neutral-axis depth per profile — what the profile is FOR."""
    print(f"  {'drift':>8}  {'curvature (1/mm)':>18}  {'neutral axis x (mm)':>20}  "
          f"{'NA depth from compressed edge':>30}")
    for _i, d, curv, x0 in rows:
        if x0 != x0:                                   # nan — no gradient to cross zero
            na = "        -- (no gradient)"
            depth = "                            --"
        else:
            na = f"{x0:20.1f}"
            # the compressed edge is the one the fitted line makes negative
            edge = LW / 2.0 if curv < 0.0 else -LW / 2.0
            depth = f"{abs(edge - x0):30.1f}"
        print(f"  {d:+7.2f}%  {curv:18.3e}  {na}  {depth}")


# --- curvature over height: the model's Fig. 11c, and its Fig. 15d point -------------------------
# Drawn here rather than in `viz` because the layout is the paper's, not the library's: curvature on
# the horizontal axis and height on the vertical, so a model profile can be laid directly beside
# Fig. 11c. The base-curvature fit is overlaid as the dashed line it is in Fig. 15(b).

MODEL_C, INK, INK_2, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#d8d7d2"


def figure_curvature(drift_pct, mids, curvature, *, indices, savepath, title,
                     phi_y: float = PHI_YIELD) -> list[tuple[int, float, float, float]]:
    """Curvature profiles over height. Returns `(index, drift%, phi_base 1/mm, L_pz mm)` per profile.

    `phi_base` and `L_pz` come from `specimen.base_curvature`, which reproduces the paper's own
    construction: a linear fit over the plastic zone extrapolated to y = 0. That is the quantity
    Fig. 15(d) plots, so it is the quantity to compare — NOT the base segment's own curvature.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = sorted(indices, key=lambda i: drift_pct[i])
    fig, ax = plt.subplots(figsize=(6.4, 7.2))
    cmap = plt.get_cmap("viridis")
    rows = []
    for k, i in enumerate(order):
        c = [v * 1e6 for v in curvature[i]]                    # 1/mm -> 1/km
        colour = cmap(k / max(1, len(order) - 1) * 0.85)
        ax.plot(c, [h / 1000.0 for h in mids], "-o", ms=3.5, lw=1.3, color=colour,
                label=f"{drift_pct[i]:+.2f}% drift")
        phi_b, lpz = base_curvature(mids, curvature[i], phi_y=phi_y)
        if phi_b == phi_b and lpz > 0:                          # not nan
            ax.plot([phi_b * 1e6, 0.0], [0.0, lpz / 1000.0], ls="--", lw=0.9, color=colour,
                    alpha=0.7)
            ax.plot([phi_b * 1e6], [0.0], marker="*", ms=10, color=colour, mec="white", mew=0.6)
        rows.append((i, drift_pct[i], phi_b, lpz))

    for sign in (1, -1):
        ax.axvline(sign * phi_y * 1e6, color=GRID, lw=0.9, ls=":")
    ax.axvline(0, color=GRID, lw=0.8)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.55)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("curvature [1/km]", fontsize=10, color=INK_2)
    ax.set_ylabel("height above foundation [m]", fontsize=10, color=INK_2)
    ax.set_ylim(0.0, mids[-1] / 1000.0 * 1.05)
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    ax.text(0.02, 0.02, f"dotted: $\\phi_y$ = {phi_y * 1e6:.0f} /km   ★: $\\phi_{{base}}$, "
                        f"linear fit over $L_{{pz}}$ extrapolated to y = 0",
            transform=ax.transAxes, fontsize=8, color=INK_2)
    ax.legend(fontsize=8, loc="upper right", frameon=True, framealpha=0.9, edgecolor=GRID)
    fig.tight_layout()
    fig.savefig(str(savepath), dpi=180, facecolor="white")
    plt.close(fig)
    return rows


# Measured base curvature, digitized from Fig. 15(d) — see `testdata.BASE_CURVATURE`.
MEASURED_BASE_CURVATURE = ((31.0, 4.2), (46.5, 9.1), (62.0, 13.4), (77.9, 18.5), (91.2, 23.4))


def print_curvature_table(rows) -> None:
    """Base curvature and plastic-zone height per profile, against the test's own Fig. 15d."""
    print(f"  {'drift':>8}  {'phi_base (1/km)':>16}  {'L_pz (mm)':>10}  {'test phi_base':>14}")
    for _i, d, phi_b, lpz in rows:
        u = abs(d) / 100.0 * A_SHEAR
        near = min(MEASURED_BASE_CURVATURE, key=lambda t: abs(t[0] - u))
        test = f"{near[1]:.1f} @ {near[0]:.0f}mm" if abs(near[0] - u) < 8.0 else "--"
        got = "     -- (elastic)" if phi_b != phi_b else f"{phi_b * 1e6:16.2f}"
        print(f"  {d:+7.2f}%  {got}  {lpz:10.0f}  {test:>14}")
