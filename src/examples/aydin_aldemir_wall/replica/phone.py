"""One tall, large-type chart of everything that matters — sized to be read on a phone.

    uv run python examples/aydin_aldemir_wall/replica/phone.py

Same data as `compare.py --single`, but portrait, with fewer series and bigger text. Deliberately
shows only the runs worth comparing: the two dead ends (bond, and the baseline that tore at its
unreinforced top band) are omitted except the baseline, which is kept precisely because its collapse
is the thing the current run is being read against.
"""
from __future__ import annotations

import json, re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent / "output" / "aydin_aldemir_wall" / "replica"
FIG = HERE.parent / "data" / "fig10b.npz"
HW = 2680.0
F_EXP = 963.592
S = re.compile(r"drift\s+\+([\d.]+)%\s+shear\s+\+?(-?[\d.]+) kN")


def from_log(p: Path):
    a = np.array([[float(m.group(1)), float(m.group(2))]
                  for L in p.read_text(errors="ignore").splitlines() if (m := S.search(L))])
    return a[:, 0] / 100 * HW, a[:, 1]


def from_json(p: Path):
    j = json.loads(p.read_text())
    return np.array(j["disp"]), np.array(j["shear"]) / 1e3


def main() -> None:
    z = np.load(FIG)
    fig, ax = plt.subplots(figsize=(7.6, 9.6))

    # --- the test ---
    bx, be = z["backbone_x"], z["backbone"]
    j = bx >= 0
    xe, ye = bx[j], be[j, 0].copy()
    m = xe >= 1.5
    ye[m] = np.maximum.accumulate(ye[m])
    ax.plot(xe, ye, color="#16a085", lw=3.4, zorder=3, label="TEST — measured envelope")
    ax.axhline(F_EXP, color="#16a085", ls=":", lw=2.0, alpha=.85, zorder=2)
    ax.text(16.6, F_EXP + 14, f"test max {F_EXP:,.0f} kN", color="#0e6655",
            fontsize=11, ha="right", fontweight="bold")

    # --- Aydin ---
    for key, hz, col, ls, lw in (("aydin15", "1.5", "#8e44ad", "-", 2.8),
                                 ("aydin301", "3.01", "#c39bd3", "--", 2.0)):
        x, y = z[f"{key}_x"], z[f"{key}_y"]
        k = x >= 0
        ax.plot(x[k], y[k], ls, color=col, lw=lw, zorder=3,
                label=f"AYDIN 2019 — horizon {hz}")

    # --- ours ---
    base = OUT / "2026-08-31_213636_pushover_d0.003_CentralDifference" / "console.log"
    if base.exists():
        x, y = from_log(base)
        ax.plot(x, y, color="#95a5a6", lw=2.4, zorder=4,
                label="OURS — bars stop 100 mm short\n(tears at the bare band, 0.153%)")
    # Both bars-to-top runs, now COMPLETE. Match each horizon EXPLICITLY: the h=3.01 directory
    # also ends in "rebartop_CentralDifference", so a loose glob silently picks whichever sorts
    # last and drops the other series.
    ours = [("*_pushover_d0.003_rebartop_CentralDifference/data.json", "h=1.5",
             "#c0392b", "-", 3.6, "peak"),
            ("*_h3.01_rebartop_CentralDifference/data.json", "h=3.01",
             "#e67e22", "--", 3.0, "end")]
    for pattern, hz, col, ls, lw, mark in ours:
        hit = sorted(OUT.glob(pattern))
        if not hit:
            continue
        x, y = from_json(hit[-1])
        i = int(np.argmax(y))
        note = ("peaked then softened" if hz == "h=1.5" else "no peak — monotonic")
        ax.plot(x, y, ls, color=col, lw=lw, zorder=6,
                label=f"OURS — bars to the top, {hz}\n({note})")
        ax.plot(x[i], y[i], "o", ms=10, color=col, zorder=7)
        ax.annotate(f"{hz}: {y[i]:,.0f} kN\n{y[i] / F_EXP:.2f} x test",
                    xy=(x[i], y[i]),
                    xytext=(x[i] + 2.2, y[i] - (250 if hz == "h=1.5" else 90)),
                    fontsize=11.5, color=col, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=col, lw=2))

    ax.set_xlim(0, 17); ax.set_ylim(0, 1400)
    ax.set_xlabel("top displacement (mm)", fontsize=14)
    ax.set_ylabel("base shear (kN)", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(alpha=.3, lw=.7)
    ax.legend(fontsize=11.5, loc="lower right", framealpha=.97, borderpad=0.9, labelspacing=0.9)
    ax.set_title("Aldemir wall — where we stand", fontsize=17, fontweight="bold", pad=14)
    fig.text(0.5, 0.012,
             "Test + Aydin curves digitized from his Fig. 10(b); the measured envelope is an UPPER\n"
             "BOUND on the backbone. Grey vs red differ ONLY in whether the longitudinal bars reach\n"
             "the top of the wall. Both horizons run to target once the bars are extended.",
             ha="center", fontsize=10, color="#555", linespacing=1.6)
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    dst = OUT / "phone.png"
    fig.savefig(dst, dpi=170)
    print(f"saved {dst}")


if __name__ == "__main__":
    main()
