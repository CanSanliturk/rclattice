"""Plot the load-displacement loops of a cyclic run — LIVE, from its rolling log.

    uv run python examples/aydin_aldemir_wall/replica/hyst.py [--run <dir>]

Reads `hysteresis.npz`, which `run.py` rewrites on every progress tick, so the loops can be watched
as they form and survive a kill. The monotonic backbone of the same model is overlaid where one
exists, since the cyclic envelope should track it.
"""
from __future__ import annotations

import argparse, json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent.parent.parent / "output" / "aydin_aldemir_wall" / "replica"
FIG10B = Path(__file__).resolve().parent.parent / "data" / "fig10b.npz"
HW = 2680.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, default=None)
    a = ap.parse_args()
    run = a.run or max((d for d in OUT.glob("*cyclic*") if (d / "hysteresis.npz").exists()),
                       key=lambda d: d.stat().st_mtime)
    z = np.load(run / "hysteresis.npz")
    u, s = z["disp"], z["shear"] / 1e3
    par = json.loads((run / "params.json").read_text())

    z10 = np.load(FIG10B) if FIG10B.exists() else None
    lv = [float(v) for v in str(par.get("cyclic", "")).split(",") if v]
    mono = sorted(OUT.glob("*_pushover_d0.003_rebartop_CentralDifference/data.json"))
    mx = my = None
    if mono:
        j = json.loads(mono[-1].read_text())
        mx, my = np.array(j["disp"]), np.array(j["shear"]) / 1e3

    # Two panels: the test's own range, then a zoom on ours. Our protocol tops out at 0.15% drift
    # (4.0 mm) against the test's ~19 mm, so on one axis our loops would be a smudge at the origin.
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.0, 11.6))

    def furnish(ax, zoom):
        if z10 is not None:
            cl = z10["cloud"]
            if zoom:
                k = (np.abs(cl[:, 0]) <= zoom * 1.08)
                cl = cl[k]
            ax.plot(cl[:, 0], cl[:, 1], ".", ms=1.6 if zoom else 1.0, color="#9fb8b3",
                    zorder=1, label="TEST — measured hysteresis (Fig. 10b)")
        if mx is not None:
            ax.plot(mx, my, color="#7f8c8d", lw=2.0, zorder=2, label="our monotonic pushover")
            ax.plot(-mx, -my, ":", color="#7f8c8d", lw=2.0, alpha=.75, zorder=2)
        ax.plot(u, s, "-", color="#c0392b", lw=1.9, zorder=4, label="our CYCLIC response")
        ax.plot(u[-1], s[-1], "o", ms=9, color="#c0392b", zorder=6)
        ax.axhline(0, color="0.65", lw=.8); ax.axvline(0, color="0.65", lw=.8)
        for d in lv:
            for sgn in (1, -1):
                ax.axvline(sgn * d * HW, color="#16a085", ls="--", lw=1.0, alpha=.5)
        ax.set_ylabel("base shear (kN)", fontsize=12)
        ax.grid(alpha=.28, lw=.6)
        ax.tick_params(labelsize=11)

    furnish(ax0, None)
    ax0.set_xlim(-20, 18); ax0.set_ylim(-1150, 1150)
    ax0.set_title("against the measured test — full range", fontsize=11)
    ax0.legend(fontsize=9, loc="lower right", framealpha=.95)

    span = max(max(abs(u)) * 1.35, (max(lv) * HW * 1.15 if lv else 1.0))
    furnish(ax1, span)
    ax1.set_xlim(-span, span)
    inwin = np.abs(z10["cloud"][:, 0]) <= span if z10 is not None else None
    ylim = max(abs(s).max() * 1.3,
               (abs(z10["cloud"][inwin, 1]).max() * 1.1 if z10 is not None and inwin.any() else 0))
    ax1.set_ylim(-ylim, ylim)
    ax1.set_xlabel("top displacement (mm)", fontsize=12)
    ax1.set_title(f"zoom on our protocol (levels {', '.join(f'{d:.3%}' for d in lv)})", fontsize=11)

    fig.suptitle(f"Aldemir wall — cyclic, horizon {par.get('horizon')}, bars to top\n"
                 f"{len(u):,} logged points, reached {max(abs(u)) / HW:.4%} drift "
                 f"of the test's ~0.71%", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    dst = run / "hysteresis.png"
    fig.savefig(dst, dpi=165)
    print(f"points {len(u):,}   drift range {u.min() / HW:+.4%} .. {u.max() / HW:+.4%}")
    print(f"shear  range {s.min():+.1f} .. {s.max():+.1f} kN")
    print(f"saved {dst}")


if __name__ == "__main__":
    main()
