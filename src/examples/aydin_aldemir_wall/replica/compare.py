"""Overlay every replica run against Aydin's published result and the measured values.

    uv run python examples/aydin_aldemir_wall/replica/compare.py

Re-reads saved runs and redraws — it never re-runs an analysis, so it is safe to call while a
pushover is still going (running runs are read from their live console log and drawn dashed).

REFERENCES COME FROM `../digitize.py` (Fig. 10b), falling back to Table 4 levels if it has not been
run. What each one is:

  * AYDIN'S OWN CURVE, horizon 1.5 — a clean coloured line in his figure, recovered directly. This
    is the true like-for-like for the replica, which is built to be his model.
  * THE MEASURED BACKBONE — a loop-tip envelope of his figure's experimental dot cloud, NOT a
    published curve. The cloud is drawn as discrete dots, so ORDER is gone and with it per-cycle
    energy, degradation and the load path; extremes survive unordered, so the envelope does not.
    Its peak reproduces Table 4's measured maximum to 1.006, which is a check, not a fit.
  * THE LEGEND in his figure occludes real data at positive displacement below about -500 kN, so
    the negative backbone is a lower bound there.

HEIGHT MISMATCH — the reason the x axis is DISPLACEMENT, not drift. The replica panel is 3000x2680,
solved from his Table 2 element counts; the specimen is 3000x2250. Equal displacement is therefore
NOT equal drift, and the test's 20 mm is 0.89% on the specimen against 0.75% on the replica. Drift
appears only as a secondary axis, labelled for the replica.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent / "output" / "aydin_aldemir_wall" / "replica"
HW_REPLICA, HW_SPECIMEN = 2680.0, 2250.0

# Table 4 of Aydin, Tuncay & Binici (2019). No curve — two numbers each.
K_SIM, F_SIM = 943.16, 1164.413        # his lattice, horizon 1.5d
K_EXP, F_EXP = 1038.44, 963.592        # measured, as reported by that paper
D_EXP = 20.0                           # mm, Table 2 — the displacement the test reached
FIG10B = HERE.parent / "data" / "fig10b.npz"


def reference_curves():
    """Digitized Fig. 10(b): (aydin_x, aydin_y, test_x, test_y) or None if not digitized yet."""
    if not FIG10B.exists():
        return None
    z = np.load(FIG10B)
    ax_, ay_ = z["aydin15_x"], z["aydin15_y"]
    k = ax_ >= 0                                   # monotonic push, to match our pushover
    # RAW upper envelope, not the monotone one. Accumulating from the origin starts the curve at a
    # flat ~430 kN, which is not a capacity: near zero displacement every loop crosses, so the
    # spread there is loop-crossing scatter. The raw envelope is ragged but it is what the cloud
    # actually shows, and it is an UPPER BOUND on the backbone — equal to it only at loop tips,
    # and above it in between, where the boundary is the unloading side of a larger cycle.
    bx, be = z["backbone_x"], z["backbone"]
    j = bx >= 0
    return ax_[k], ay_[k], bx[j], be[j, 0]

STEP = re.compile(r"drift\s+\+([\d.]+)%\s+shear\s+\+?(-?[\d.]+) kN")


def from_log(path: Path) -> tuple[list[float], list[float]]:
    """Sampled (mm, kN) from a progress log — coarse, but it is all an aborted run leaves."""
    d, s = [], []
    for line in path.read_text(errors="ignore").splitlines():
        m = STEP.search(line)
        if m:
            d.append(float(m.group(1)) / 100.0 * HW_REPLICA)
            s.append(float(m.group(2)))
    return d, s


def from_json(path: Path) -> tuple[list[float], list[float]]:
    j = json.loads(path.read_text())
    return list(j["disp"]), [v / 1e3 for v in j["shear"]]


def collect() -> list[dict]:
    """Every run we have, richest source first: data.json (full) then console.log (sampled)."""
    runs = []
    for d in sorted(OUT.glob("*")):
        if not d.is_dir() or "elastic" in d.name:
            continue
        dj, cl = d / "data.json", d / "console.log"
        if dj.exists():
            x, y = from_json(dj)
            runs.append(dict(dir=d, x=x, y=y, full=True, live=False))
        elif cl.exists():
            x, y = from_log(cl)
            if x:
                runs.append(dict(dir=d, x=x, y=y, full=False, live=False))
    for live in OUT.glob("_live_*.log"):                    # a run still going
        x, y = from_log(live)
        if x:
            runs.append(dict(dir=live, x=x, y=y, full=False, live=True))
    return runs


def label_for(name: str) -> tuple[str, str, str]:
    """(label, colour, linestyle) — keyed off the directory/log name."""
    if "bond0.004" in name or "bond004" in name:
        return "bond, ratio 0.004 (t=210)", "#c0392b", "-"
    if "bond0.01" in name or "_live_bond." in name:
        return "bond, ratio 0.01 (t=210)", "#e67e22", "-"
    if "2026-08-29" in name:
        return "no bond, t=120  [superseded]", "#95a5a6", "-"
    if "d0.003" in name:
        return "no bond, t=210 — to collapse", "#1b2a41", "-"
    return "no bond, t=210  (stopped at 0.15%)", "#5d7a9e", "-"


def draw(ax, runs, xmax, ymax, *, star: bool) -> None:
    """One panel. `ymax` is fixed BEFORE annotating so a reference outside it is not mislabelled."""
    for r in runs:
        lab, col, ls = label_for(r["dir"].name)
        if r["live"]:
            lab += "  [running]"
        ax.plot(r["x"], r["y"], ls, color=col, lw=1.8 if r["full"] else 1.4,
                marker="" if r["full"] else "o", ms=3.5,
                alpha=0.55 if "superseded" in lab else 1.0, label=lab, zorder=3)
        if r["full"]:      # mark where a completed run STOPPED, which may not be a capacity
            ax.plot(r["x"][-1], r["y"][-1], "s", color=col, ms=6, zorder=4)

    # --- references, digitized from Fig. 10(b) where available ---
    ref = reference_curves()
    if ref is not None:
        ax_, ay_, tx, ty = ref
        ax.plot(ax_, ay_, "-", color="#8e44ad", lw=2.0, alpha=.9, zorder=2,
                label="Aydin 2019, horizon 1.5  [Fig. 10b]")
        ax.plot(tx, ty, "-", color="#16a085", lw=2.4, alpha=.9, zorder=2,
                label="measured: cloud envelope  [upper bound on backbone]")
    else:                                   # not digitized yet: Table 4 levels only
        for k, f, col, who in ((K_SIM, F_SIM, "#8e44ad", "Aydin 2019 lattice"),
                               (K_EXP, F_EXP, "#16a085", "measured (Table 4)")):
            ax.plot([0, min(f / k, xmax)], [0, min(f, k * xmax)], ":", color=col, lw=1.6, zorder=2)
            if f <= ymax:
                ax.axhline(f, color=col, ls="--", lw=1.4, alpha=0.8, zorder=2)
                ax.text(xmax * 0.015, f, f"{who}: {f:,.0f} kN", color=col, fontsize=7.5,
                        va="bottom", ha="left")
    if star:
        ax.plot([D_EXP], [F_EXP], "*", color="#16a085", ms=15, zorder=5,
                label=f"Table 4 max: {F_EXP:,.0f} kN")

    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("top displacement (mm)")
    ax.grid(alpha=0.25, lw=0.5)
    sec = ax.secondary_xaxis("top", functions=(lambda v: v / HW_REPLICA * 100,
                                               lambda v: v * HW_REPLICA / 100))
    sec.set_xlabel("replica drift (%)", fontsize=8)
    sec.tick_params(labelsize=7)


def single(runs) -> None:
    """One summary panel: the test, Aydin's two models, and our two headline runs.

    Deliberately drawn over the FULL test range rather than the range our runs reached, because the
    truncation is itself a result — our pushover stops at 4 mm of the test's 20.
    """
    z = np.load(FIG10B) if FIG10B.exists() else None
    if z is None:
        raise SystemExit("run examples/aydin_aldemir_wall/digitize.py first")

    fig, ax = plt.subplots(figsize=(10, 7))

    # --- the measured test, positive quadrant ---
    cl = z["cloud"]
    k = (cl[:, 0] >= 0) & (cl[:, 1] >= 0)
    ax.plot(cl[k, 0], cl[k, 1], ".", ms=1.1, color="#b9c6c3", zorder=1,
            label="Aldemir test — measured hysteresis (Fig. 10b dots)")
    bx, be = z["backbone_x"], z["backbone"]
    j = (bx >= 0)
    x_env, y_env = bx[j], be[j, 0].copy()
    # Beyond the pinch the loop tips grow with amplitude, so the envelope is non-decreasing there;
    # a dip is a bin that happened to catch no tip, not a loss of strength. Take a running max from
    # the pinch outward. NOT applied from the origin — near zero every loop crosses, so accumulating
    # from there would invent a flat ~430 kN plateau that is loop-crossing scatter, not capacity.
    pinch = 1.5
    m = x_env >= pinch
    y_env[m] = np.maximum.accumulate(y_env[m])
    ax.plot(x_env, y_env, "-", lw=3.0, color="#16a085", zorder=4,
            label="Aldemir test — backbone (cloud envelope, upper bound)")
    mkx = float(z["marker"][0])          # where the paper plots its own "experiment" point
    ax.plot([mkx], [F_EXP], "*", ms=18, color="#0e6655", zorder=6,
            label=f"Aldemir test — Table 4 max, {F_EXP:,.0f} kN")

    # --- Aydin's own lattice, both horizons ---
    for key, hz, col, lw, ls in (("aydin15", "1.5", "#8e44ad", 2.6, "-"),
                                 ("aydin301", "3.01", "#c39bd3", 1.9, "--")):
        x, y = z[f"{key}_x"], z[f"{key}_y"]
        m = x >= 0
        ax.plot(x[m], y[m], ls, lw=lw, color=col, zorder=3,
                label=f"Aydin 2019 lattice — horizon {hz}")

    # --- our runs ---
    style = {"no bond, t=210 — to collapse": ("#1b2a41", 2.8, "-",
                                                "our lattice — no bond (collapses at 0.153%)"),
             "bond, ratio 0.004 (t=210)": ("#c0392b", 2.2, "-", "our lattice — with bond")}
    for r in runs:
        lab, _c, _ls = label_for(r["dir"].name)
        if lab not in style:
            continue
        col, lw, ls, nice = style[lab]
        ax.plot(r["x"], r["y"], ls, lw=lw, color=col, zorder=5,
                marker="" if r["full"] else "o", ms=4, label=nice)
        ax.plot(r["x"][-1], r["y"][-1], "s", ms=7, color=col, zorder=6)
    ax.annotate("tracks Aydin's own curve to 1.018 (mean, 0.5-4 mm),\n"
                "then loses its load path ENTIRELY in 0.006% of drift:\n"
                "868.6 kN at 4.0 mm -> 28.6 kN at 4.2 mm",
                xy=(4.15, 260), xytext=(5.7, 470), fontsize=8.5, color="#1b2a41",
                arrowprops=dict(arrowstyle="->", color="#1b2a41", lw=1.1))
    ax.annotate("bond turns over at 0.4 mm\nat BOTH ratios tried", xy=(0.9, 250),
                xytext=(2.2, 120), fontsize=8.5, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.1))

    ax.set_xlim(0, 17.5); ax.set_ylim(0, 1420)
    ax.set_xlabel("lateral displacement (mm)")
    ax.set_ylabel("lateral load / base shear (kN)")
    ax.grid(alpha=.25, lw=.5)
    ax.legend(fontsize=8.6, loc="lower right", framealpha=.96, borderpad=0.8)
    ax.set_title("Aldemir wall (Aldemir, Binici & Canbay 2017) — test, Aydin's lattice, and ours",
                 fontsize=12)
    fig.text(0.5, 0.014,
             "Test data and both Aydin curves digitized from Fig. 10(b) of Aydin, Tuncay & Binici "
             "(2019). Calibration uses gridline geometry only, then self-checks to 0.995 / 0.985\n"
             "on his Table 4 plateaus and 1.006 on the measured maximum — three numbers it never "
             "saw. The backbone is the cloud's outer ENVELOPE and an UPPER BOUND: the dots are\n"
             "unordered, so a boundary point may be a loop tip or the unloading side of a larger "
             "cycle. The square marks where our bond run was stopped; the no-bond run was stopped AFTER its collapse.",
             ha="center", va="bottom", fontsize=7.6, color="#555", linespacing=1.5)
    fig.tight_layout(rect=(0, 0.088, 1, 1))
    dst = OUT / "summary.png"
    fig.savefig(dst, dpi=175)
    print(f"saved {dst}")


def main() -> None:
    runs = collect()
    if not runs:
        raise SystemExit(f"no runs found under {OUT}")

    reach = max(max(r["x"]) for r in runs)
    bond = [r for r in runs if "bond" in r["dir"].name]
    bx = max(max(r["x"]) for r in bond) * 1.15 if bond else reach

    if "--single" in __import__("sys").argv:
        single(runs)
        return

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16.5, 5.6))

    draw(a1, runs, D_EXP * 1.05, max(F_SIM, F_EXP) * 1.18, star=True)
    a1.set_ylabel("base shear (kN)")
    a1.set_title("Full test range — how far the runs got", fontsize=10)
    a1.legend(fontsize=7.2, loc="center right", framealpha=0.95)

    draw(a2, runs, reach * 1.05, max(max(r["y"]) for r in runs) * 1.22, star=False)
    a2.set_title(f"Range reached (0–{reach * 1.05:.1f} mm)", fontsize=10)

    draw(a3, runs, bx, max(max(r["y"]) for r in bond) * 1.35 if bond else 100, star=False)
    a3.set_title(f"Bond detail (0–{bx:.2f} mm) — where both bond runs turn over", fontsize=10)

    fig.suptitle("Aldemir wall — Aydin (2019) replica: every run vs the published values",
                 fontsize=12.5, y=0.985)
    fig.text(0.5, 0.008,
             "References digitized from the paper's Fig. 10(b) — Aydin's own horizon-1.5 curve, and a "
             "the OUTER ENVELOPE of its experimental dot cloud (peak 1.006x Table 4's measured "
             "maximum).\nThat envelope is an UPPER BOUND on the backbone, not the backbone: the cloud "
             "is unordered, so a point on it may be a loop tip or the unloading side of a larger "
             "cycle. Per-cycle energy and degradation are NOT recoverable, and the legend occludes "
             "data at positive displacement below about -500 kN. Squares mark where a completed "
             "run stopped: the t=210 control was still ASCENDING at its target, so its 868.6 kN is a "
             "LOWER BOUND. Replica panel is 2680 mm tall against the specimen's 2250.",
             ha="center", fontsize=7.4, color="#444")
    fig.tight_layout(rect=(0, 0.062, 1, 0.955))

    dst = OUT / "comparison.png"
    fig.savefig(dst, dpi=170)
    print(f"saved {dst}")
    for r in runs:
        print(f"  {r['dir'].name:<58s} {len(r['x']):>7,} pts  "
              f"max {max(r['y']):>7.1f} kN at {r['x'][r['y'].index(max(r['y']))]:.3f} mm"
              f"{'  [sampled]' if not r['full'] else ''}{'  [running]' if r['live'] else ''}")


if __name__ == "__main__":
    main()
