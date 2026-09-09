"""Re-draw the cyclic figures from a saved run, without repeating the analysis.

`cyclic.py` persists its raw drift/shear response next to the figures. These runs cost tens of
minutes to hours, so any change to how the results are PRESENTED — a new overlay, different
labels, a rescaled axis — should never require paying for the analysis again.

Run from src/:  python examples/wall/replot.py [--stem wall_cyclic_dynamic]
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from rclattice import viz

import gauge
from cyclic import PAPER_PEAK_KN, backbone, paper_backbone_curves
from specimen import A_SHEAR, OUT, PROTOCOL_PEAKS

# Digitized test LOOPS from Fig. 14b (see digitize.py). A point cloud, not an ordered path —
# overlapping loops cannot be re-sequenced from pixels — so it is drawn with markers, never lines.
DIGITIZED = Path(__file__).resolve().parent / "data" / "swncff_fig14b.npz"
TEST_STYLE = {"color": "#d62728", "ls": "", "marker": ".", "ms": 0.7, "alpha": 0.5}


def measured_loops():
    """(drift %, load kN) of the digitized test hysteresis, or None if not yet extracted."""
    if not DIGITIZED.exists():
        return None
    d = np.load(DIGITIZED)
    return d["drift_pct"], d["load_kN"]


def main(*, stem: str = "wall_cyclic_dynamic", levels_pct=None) -> None:
    data = json.loads((OUT / f"{stem}_data.json").read_text())
    d_pct, s_kn = data["drift_pct"], data["shear_kN"]
    compression, converged = data["compression"], data["converged"]
    print(f"{stem}: {len(d_pct)} points, compression={compression}, converged={converged}, "
          f"max drift +-{max(abs(d) for d in d_pct):.3f}%")

    exp_pos, exp_neg = paper_backbone_curves()
    exp_d = [d for d, _v in exp_neg][::-1] + [d for d, _v in exp_pos]
    exp_s = [v for _d, v in exp_neg][::-1] + [v for _d, v in exp_pos]
    exp_style = {"color": "C3", "ls": "--", "lw": 2, "marker": "o", "ms": 5}

    loops = measured_loops()
    model_series = {"disp": d_pct, "shear": s_kn,
                    "label": f"lattice, compression={compression}", "style": {"color": "C0", "lw": 0.9}}
    if loops is not None:
        td, tl = loops
        test_series = {"disp": list(td), "shear": list(tl),
                       "label": "SW-NC-FF test, digitized Fig. 14b", "style": TEST_STYLE}
        print(f"overlaying {len(td)} digitized test points")
    else:
        test_series = {"disp": exp_d, "shear": exp_s,
                       "label": "SW-NC-FF measured backbone", "style": exp_style}
        print("digitized loops not found — falling back to the reported backbone milestones")

    viz.figure_hysteresis(
        [model_series, test_series],
        savepath=str(OUT / f"{stem}.png"), drift_label="drift ratio (%)",
        shear_label="base shear (kN)",
        title=f"SW-NC-FF — Aydin-calibrated lattice vs test (Fig. 14b digitized); "
              f"test peak {PAPER_PEAK_KN:.0f} kN")

    # envelope of the loops, in model units, against the measured backbone
    pos, neg = backbone([d * A_SHEAR / 100.0 for d in d_pct], [s * 1e3 for s in s_kn])
    mdl_d = [u / A_SHEAR * 100.0 for u, _s in neg][::-1] + [u / A_SHEAR * 100.0 for u, _s in pos]
    mdl_s = [s / 1e3 for _u, s in neg][::-1] + [s / 1e3 for _u, s in pos]
    peak = max(abs(s) for s in mdl_s)
    viz.figure_pushover(
        [{"disp": mdl_d, "shear": mdl_s, "label": f"lattice envelope, peak {peak:.0f} kN",
          "style": {"color": "C0", "lw": 2}},
         {"disp": exp_d, "shear": exp_s,
          "label": f"SW-NC-FF measured backbone, peak {PAPER_PEAK_KN:.0f} kN", "style": exp_style}],
        savepath=str(OUT / f"{stem}_backbone.png"), xlabel="drift ratio (%)",
        ylabel="base shear (kN)",
        title="SW-NC-FF — lattice cyclic envelope vs the measured backbone\n"
              "(test points reconstructed from reported values, not digitized loops)")
    print(f"redrew {stem}.png and {stem}_backbone.png")

    # Base strain profile, redrawn from the saved gauge series (D63). This is the whole reason the
    # gauge is persisted: choosing which drift levels to show, or restyling the figure, must never
    # cost another run.
    g = data.get("gauge")
    if not g or not g["strain"]:
        print("no gauge series in this run — re-run to record one")
        return
    gx, gd, gs = g["x_mm"], g["drift_pct"], g["strain"]
    # Tolerance is RELATIVE: a run stops a hair short of its nominal amplitude (the drive lands
    # between steps), so an exact `<= reached` test drops the very level the run was aiming at.
    reached = max(abs(d) for d in gd) * 1.005 + 1e-6
    levels = levels_pct or [p / A_SHEAR * 100.0 for p in PROTOCOL_PEAKS
                            if p / A_SHEAR * 100.0 <= reached]
    idx = gauge.select_levels(gd, levels)
    sp = OUT / f"{stem}_strain_profile.png"
    rows = gauge.figure(gd, gs, gx, indices=idx, savepath=sp, gauge=g["gauge_mm"],
                        title="SW-NC-FF — vertical strain across the base section, "
                              "at each protocol level")
    print(f"\nbase strain gauge: {len(gd)} profiles saved, {len(idx)} drawn "
          f"(levels {', '.join(f'{lv:g}%' for lv in levels)})")
    gauge.print_table(rows)
    print(f"redrew {sp}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Re-draw cyclic figures from a saved run")
    p.add_argument("--stem", default="wall_cyclic_dynamic", help="output stem of the saved run")
    p.add_argument("--levels", type=float, nargs="*", default=None,
                   help="drift levels (%%) to draw strain profiles at; default = every protocol "
                        "level the run reached, e.g. --levels 0.3 1.0 2.0")
    a = p.parse_args()
    main(stem=a.stem, levels_pct=a.levels)
