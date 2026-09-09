"""Damage history of a cyclic run: one subfigure per captured instant, on a shared strain scale.

    uv run python examples/aydin_aldemir_wall/replica/damage_history.py [--run <dir>] [--cols 4]

Reads the archived `frames/*.npz` (a complete history if the archiver ran) and falls back to the
live `snapshots.npz` ring (the last 10 frames only). Each panel is the crack pattern at that
instant, captioned with where it sits in the load history, and a small load-displacement inset
marks the instant on the hysteresis so a panel can be located at a glance.

SHARED COLOUR SCALE across panels — otherwise each panel self-normalizes and a lightly damaged
early state looks identical to a heavily damaged late one. Struts past `ft/E` are drawn as cracked;
compression is NOT classified, because these struts have no compressive strength to exceed.
"""
from __future__ import annotations

import argparse, json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from build import calibrate, wall_lattice
from specimen import EC, FT

OUT = Path(__file__).resolve().parent.parent.parent / "output" / "aydin_aldemir_wall" / "replica"
HW = 2680.0
EPS_CRACK = FT / EC


def frames(run: Path):
    fd = run / "frames"
    got = []
    if fd.exists():
        for p in sorted(fd.glob("*.npz"), key=lambda q: int(q.stem)):
            z = np.load(p)
            got.append((z["ids"], z["u"], z["meta"]))
    if not got and (run / "snapshots.npz").exists():
        z = np.load(run / "snapshots.npz", allow_pickle=True)
        got = [(z["ids"], u, m) for u, m in zip(z["ring_u"], z["ring_meta"])]
    return got


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, default=None)
    ap.add_argument("--cols", type=int, default=2,
                    help="2 = portrait, readable on a phone; raise for a wide sheet")
    ap.add_argument("--max-panels", type=int, default=8)
    a = ap.parse_args()

    run = a.run or max((d for d in OUT.glob("*cyclic*") if d.is_dir()),
                       key=lambda d: d.stat().st_mtime)
    par = json.loads((run / "params.json").read_text())
    fr = frames(run)
    if not fr:
        raise SystemExit(f"{run.name}: no frames")
    if len(fr) > a.max_panels:                       # keep the ends, thin the middle
        idx = np.unique(np.linspace(0, len(fr) - 1, a.max_panels).astype(int))
        fr = [fr[i] for i in idx]

    cal = calibrate(mesh_size=par["mesh"], horizon=par["horizon"])
    model, _ = wall_lattice(cal.area, mesh_size=par["mesh"], horizon=par["horizon"],
                            nonlinear=True, full_height_rebar=bool(par.get("full_height_rebar")),
                            cyclic=bool(par.get("cyclic")),
                            cyclic_law=par.get("cyclic_law", "aydin"),
                            gf_factor=par.get("gf_factor", 1.0))
    from rclattice import viz
    ids0, idx, pts, _l, _q = viz._arrays(model)
    els = [e for e in model.elements if len(e.nodes) == 2]
    segs = np.array([[pts[idx[e.nodes[0]]], pts[idx[e.nodes[1]]]] for e in els])

    strains, metas = [], []
    for ids, u, m in fr:
        disp = {int(i): list(map(float, r)) for i, r in zip(ids, u)}
        _e, st = viz.strut_strains(model, disp)
        strains.append(st); metas.append(m)
    vmax = max(float(np.abs(s).max()) for s in strains)
    norm = matplotlib.colors.SymLogNorm(linthresh=EPS_CRACK, vmin=-vmax, vmax=vmax, base=10)

    n = len(fr); cols = min(a.cols, n); rows = int(np.ceil(n / cols))
    pw = 5.2 if cols <= 2 else 3.4
    fig, axes = plt.subplots(rows, cols, figsize=(pw * cols, pw * 0.95 * rows + 1.9),
                             squeeze=False)
    for k in range(rows * cols):
        ax = axes[k // cols][k % cols]
        if k >= n:
            ax.axis("off"); continue
        st, m = strains[k], metas[k]
        lc = LineCollection(list(segs), cmap="coolwarm", norm=norm, linewidths=.35)
        lc.set_array(st); ax.add_collection(lc)
        ncr = int((st > EPS_CRACK).sum())
        ax.set_xlim(pts[:, 0].min(), pts[:, 0].max())
        ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"step {int(m[0]):,}\n{m[1] / HW * 100:+.4f}% drift, {m[2] / 1e3:+,.0f} kN\n"
                     f"{ncr:,} cracked ({ncr / len(st):.1%})",
                     fontsize=12 if cols <= 2 else 8.5)

    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="coolwarm"),
                 ax=axes, fraction=.022, pad=.015,
                 label="axial strain (compression < 0 < tension)")
    fig.suptitle(f"Damage history — cyclic, horizon {par.get('horizon')}, bars to top\n"
                 f"shared symlog strain scale — panels ARE comparable\n"
                 f"tension-only struts: compression is not classified as damage",
                 fontsize=14 if cols <= 2 else 12)
    dst = run / "damage_history.png"
    fig.savefig(dst, dpi=145, bbox_inches="tight")
    print(f"{n} panels from {len(frames(run))} frames -> {dst}")


if __name__ == "__main__":
    main()
