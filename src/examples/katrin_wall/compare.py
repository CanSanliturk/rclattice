"""Overlay the lattice cyclic response on the DIGITIZED WSH3 test loops (paper Fig. 7c).

Reads two things and draws one figure — it never runs an analysis:

  * the saved lattice response, `examples/output/katrin_wall/<stem>_data.json` (from cyclic.py);
  * the digitized test hysteresis, `examples/katrin_wall/data/wsh3_fig7.npz` (from digitize.py).

Both are expressed as (drift %, base shear kN). The model's drift is computed on the MODEL's shear
span (4550 mm) and the test's on the paper's (4560 mm) — a 0.22% relative difference, well inside
the digitization error, not corrected for.

The test data is a point CLOUD, not an ordered path — overlapping loops cannot be re-sequenced from
pixels — so it is drawn with fine markers and never joined by lines. Its BACKBONE need not be
recovered from the cloud at all: `digitize.backbone` reads the tips at the protocol's own
amplitudes, which are known exactly, and lands within one pixel of Table 5. The slice-based
envelope below is kept as an independent cross-check of that.

READ BOTH PANELS. That is the difference from the sibling SW-NC-FF study, where loop SHAPE was a
known-unfair comparison because 74% of that test's displacement was rocking on debonded plain bars.
WSH3's deformed bars stay bonded and its displacement is flexure plus ~12% shear, so the loops are
fair too — up to 1.70% drift, where the real wall's boundary bars began to buckle and this model's
`Steel02` cannot follow.

Run from src/:  python examples/katrin_wall/compare.py [--stem wsh3_cyclic_dynamic]
Output: examples/output/katrin_wall/<stem>_vs_test.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cyclic import PAPER_PEAK_KN, PAPER_PEAK_NEG_KN, backbone
from specimen import A_SHEAR, OUT, PROTOCOL_PEAKS, SPECIMEN

DIGITIZED = Path(__file__).resolve().parent / "data" / "wsh3_fig7.npz"

# Validated categorical pair (dataviz slots 1 and 2); worst-pair CVD dE 24.7, normal-vision 33.6.
MODEL_C, TEST_C = "#2a78d6", "#eb6834"
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#d8d7d2"

# Drift levels of the test protocol (Fig. 6): the 0.75*F_y cycle, then mu = 2..6 at multiples of
# delta_y = 15.4 mm. Read off the specimen so the two can never drift apart.
PROTOCOL_LEVELS = tuple(round(p / A_SHEAR * 100.0, 4) for p in PROTOCOL_PEAKS)


def test_loops():
    """(drift %, load kN) point cloud of the digitized test hysteresis."""
    if not DIGITIZED.exists():
        raise FileNotFoundError(f"{DIGITIZED} not found — run `python examples/katrin_wall/digitize.py`")
    d = np.load(DIGITIZED)
    return d["drift_pct"], d["load_kN"]


def test_backbone(drift, load, levels=PROTOCOL_LEVELS, inner: float = 0.10, outer: float = 0.05):
    """Loop-tip envelope of the digitized cloud: the peak load of the cycle at each drift level.

    Only points near a level's own drift extreme belong to that cycle's tip, so the load is taken
    over the slice [(1-inner)*d, (1+outer)*d]. The slice must reach slightly INSIDE the extreme
    because a degrading wall peaks just before its peak displacement — at the tip proper the load
    has already come off by a few percent. `inner = 0.10` is the widest slice that still excludes
    the reloading branches of LARGER loops passing through: at 0.15 the 3.0% level jumps from 152
    to 174 kN, which is the 4% loop, not the 3% cycle.

    Levels the figure does not reach are skipped rather than guessed.
    """
    pos, neg = [(0.0, 0.0)], [(0.0, 0.0)]
    for d in levels:
        for out, sign in ((pos, +1), (neg, -1)):
            band = (np.abs(drift) >= d * (1.0 - inner)) & (np.abs(drift) <= d * (1.0 + outer)) \
                & (np.sign(drift) == sign)
            if band.sum() < 5:
                continue
            v = load[band]
            out.append((sign * d, v.max() if sign > 0 else v.min()))
    return pos, neg


def _style(ax, *, xlim, ylim):
    ax.axhline(0, color=GRID, lw=0.8, zorder=1)
    ax.axvline(0, color=GRID, lw=0.8, zorder=1)
    ax.grid(True, color=GRID, lw=0.5, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


def main(*, stem: str = "wsh3_cyclic_dynamic") -> None:
    data = json.loads((OUT / f"{stem}_data.json").read_text())
    md, ms = np.asarray(data["drift_pct"]), np.asarray(data["shear_kN"])
    td, tl = test_loops()

    m_pos, m_neg = backbone(md * A_SHEAR / 100.0, ms * 1e3)
    m_pos = [(u / A_SHEAR * 100.0, s / 1e3) for u, s in m_pos]
    m_neg = [(u / A_SHEAR * 100.0, s / 1e3) for u, s in m_neg]
    t_pos, t_neg = test_backbone(td, tl)

    m_peak, m_peak_n = max(s for _d, s in m_pos), min(s for _d, s in m_neg)
    t_peak, t_peak_n = max(s for _d, s in t_pos), min(s for _d, s in t_neg)
    print(f"{stem}: {len(md)} model points, compression={data['compression']}, "
          f"converged={data['converged']}, max drift +-{np.abs(md).max():.2f}%")
    print(f"digitized test: {len(td)} points, drift {td.min():+.2f}..{td.max():+.2f}%")
    print(f"  peak base shear   model {m_peak:+.1f} / {m_peak_n:+.1f} kN")
    print(f"                    test  {t_peak:+.1f} / {t_peak_n:+.1f} kN "
          f"(Table 5 V_max {PAPER_PEAK_KN:+.1f})")
    print(f"  model / test      {m_peak / t_peak:.3f} push, {m_peak_n / t_peak_n:.3f} pull")

    lim_d = max(np.abs(md).max(), np.abs(td).max()) * 1.06
    lim_s = max(np.abs(ms).max(), np.abs(tl).max()) * 1.12

    fig = plt.figure(figsize=(13.2, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.32, 1.0), wspace=0.20,
                          left=0.055, right=0.985, top=0.715, bottom=0.105)
    ax, bx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # --- left: the loops, on top of each other ------------------------------------------------
    ax.scatter(td, tl, s=0.5, c=TEST_C, alpha=0.55, linewidths=0, zorder=2, rasterized=True)
    ax.plot(md, ms, color=MODEL_C, lw=0.55, alpha=0.85, solid_joinstyle="round", zorder=3)
    _style(ax, xlim=(-lim_d, lim_d), ylim=(-lim_s, lim_s))
    ax.set_xlabel("drift ratio (%)", fontsize=10, color=INK_2)
    ax.set_ylabel("base shear (kN)", fontsize=10, color=INK_2)
    ax.set_title("Hysteresis loops, overlaid", fontsize=11, color=INK, loc="left", pad=34)

    top = ax.secondary_xaxis("top", functions=(lambda d: d / 100.0 * A_SHEAR,
                                               lambda u: u / A_SHEAR * 100.0))
    top.set_xlabel(f"lateral displacement at {A_SHEAR:.0f} mm (mm)", fontsize=9,
                   color=INK_2, labelpad=4)
    top.tick_params(colors=INK_2, labelsize=8.5)
    top.spines["top"].set_color(GRID)

    # Legend sits INSIDE the axes, in the upper-left quadrant: for drift <= -1.4% the response is
    # entirely below zero shear, so nothing is hidden. Labels are kept short to keep the box small.
    leg = ax.legend([plt.Line2D([], [], color=MODEL_C, lw=1.6),
                     plt.Line2D([], [], color=TEST_C, ls="", marker="o", ms=4)],
                    ["lattice model", f"{SPECIMEN} test (digitized)"],
                    fontsize=8, loc="upper left", frameon=True, framealpha=0.9,
                    edgecolor=GRID, borderpad=0.45, handlelength=1.5, handletextpad=0.6,
                    labelspacing=0.4)
    leg.get_frame().set_linewidth(0.5)
    for txt in leg.get_texts():
        txt.set_color(INK_2)

    # --- right: the envelopes, where strength is actually comparable ---------------------------
    for pos, neg, color, lw in ((m_pos, m_neg, MODEL_C, 2.0), (t_pos, t_neg, TEST_C, 2.0)):
        d = [x for x, _y in neg][::-1] + [x for x, _y in pos]
        s = [y for _x, y in neg][::-1] + [y for _x, y in pos]
        bx.plot(d, s, color=color, lw=lw, zorder=3, solid_capstyle="round")
    bx.plot([d for d, _s in t_pos] + [d for d, _s in t_neg],
            [s for _d, s in t_pos] + [s for _d, s in t_neg],
            ls="", marker="o", ms=4.5, mfc="white", mec=TEST_C, mew=1.4, zorder=4)
    # mark the two peaks; the VALUES go in the legend, where they cannot collide with the curves
    for drift, value, color in ((next(d for d, s in m_pos if s == m_peak), m_peak, MODEL_C),
                                (next(d for d, s in t_pos if s == t_peak), t_peak, TEST_C)):
        bx.plot([drift], [value], marker="*", ms=13, color=color, mec="white", mew=0.8, zorder=6)
    _style(bx, xlim=(-lim_d, lim_d), ylim=(-lim_s, lim_s))
    bx.set_xlabel("drift ratio (%)", fontsize=10, color=INK_2)
    bx.set_ylabel("base shear (kN)", fontsize=10, color=INK_2)
    bx.set_title("Loop envelopes", fontsize=11, color=INK, loc="left", pad=34)

    # peak values as direct labels in the empty upper-left quadrant — no box, nothing to collide with
    for y, text, color in ((0.94, f"model peak  {m_peak:.0f} kN", MODEL_C),
                           (0.87, f"test peak  {t_peak:.0f} kN", TEST_C)):
        bx.text(0.035, y, text, transform=bx.transAxes, fontsize=10, color=color, va="top")
    bx.text(0.035, 0.795, "(★ marked on the envelopes)", transform=bx.transAxes,
            fontsize=8.5, color=INK_2, va="top")

    fig.suptitle(f"{SPECIMEN} shear wall — Aydin-calibrated lattice vs the measured test hysteresis",
                 fontsize=14, color=INK, x=0.055, ha="left", y=0.972)
    fig.text(0.055, 0.905,
             f"Lattice model (Aydin-calibrated, compression = {data['compression']}) vs the "
             f"{SPECIMEN} test loops digitized from Fig. 7c — digitized peaks land within one "
             f"pixel of Table 5's {PAPER_PEAK_KN:.0f} kN.",
             fontsize=9.5, color=INK_2, ha="left")
    fig.text(0.055, 0.868,
             f"Peak base shear: model {m_peak:.0f} kN vs test {t_peak:.0f} kN — "
             f"{m_peak / t_peak:.2f}x push, {m_peak_n / t_peak_n:.2f}x pull.",
             fontsize=9.5, color=INK_2, ha="left")
    fig.text(0.055, 0.831,
             "Deformed bars, bonded: flexure plus ~12% shear, so loop SHAPE is comparable too — "
             "until 1.70% drift, where the test's boundary bars buckled and Steel02 cannot follow.",
             fontsize=9.5, color=INK_2, ha="left")

    savepath = OUT / f"{stem}_vs_test.png"
    fig.savefig(savepath, dpi=200, facecolor="white")
    print(f"\nsaved {savepath}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Overlay the lattice cyclic run on the digitized test")
    p.add_argument("--stem", default="wsh3_cyclic_dynamic", help="output stem of the saved run")
    a = p.parse_args()
    main(stem=a.stem)
