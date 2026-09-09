"""Vertical strain gauge over the wall base — reduction, persistence and figure (D63).

The measurement: two aligned rows of nodes, one on the wall-pedestal face (y = 0) and one
`GAUGE_H` above it, differenced column by column and divided by the gauge length, giving the
average vertical strain of each 250 mm column across the wall width. That is the same construction
as the line of vertical LVDTs on the tested wall's face (the paper's Fig. 21 flexure/shear/rocking
split is derived from them), so the model quantity and the measured one are the same thing rather
than merely analogous.

Why the base and why 250 mm: a flexure-controlled wall does its work in the base region, and this
profile is what characterises it — its slope is the curvature, its zero crossing is the neutral
axis, and its departure from a straight line is the direct test of "plane sections remain plane",
which a lattice is under no obligation to satisfy.

Everything here works off a runner's `node_history` result, and everything it computes is also
written to the run's JSON, so a figure can be redrawn without re-running the analysis — the whole
point, given the analysis costs hours.
"""

from __future__ import annotations

from rclattice import viz

from specimen import A_SHEAR, GAUGE_H, LW, gauge_strains, neutral_axis


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


def payload(res, gx, *, gauge: float = GAUGE_H) -> dict:
    """The gauge block to embed in a run's JSON: enough to redraw anything without re-running."""
    drift, strains = reduce_history(res, gx)
    return {"x_mm": list(gx), "gauge_mm": gauge, "drift_pct": drift, "strain": strains}


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
