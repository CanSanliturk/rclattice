"""Two extra damage views for the status page: panels AT THE REVERSALS, and damage evolution.

    uv run python examples/aydin_aldemir_wall/replica/damage_extra.py [--run <dir>]

`damage_history.py` samples frames evenly, which lands them at arbitrary points in the protocol.
These two are keyed to the load history instead:

  * `damage_reversals.png` — the crack pattern at each displacement REVERSAL, i.e. the loop tips,
    which are the instants the wall is at its worst on each half-cycle;
  * `damage_evolution.png` — cracked fraction against drift and against step, with the reversals
    marked, so accumulation can be read as a curve rather than inferred from panel captions.

Strains are computed vectorized here rather than through `viz.strut_strains`: that loops in Python
over every element, which is fine for one frame and far too slow for 130.
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, default=None)
    a = ap.parse_args()
    run = a.run or max((d for d in OUT.glob("*cyclic*") if (d / "frames").exists()),
                       key=lambda d: d.stat().st_mtime)
    par = json.loads((run / "params.json").read_text())

    frames = []
    for p in sorted((run / "frames").glob("*.npz"), key=lambda q: int(q.stem)):
        z = np.load(p)
        frames.append((z["ids"], z["u"], z["meta"]))
    if not frames:
        raise SystemExit("no frames")

    cal = calibrate(mesh_size=par["mesh"], horizon=par["horizon"])
    model, _ = wall_lattice(cal.area, mesh_size=par["mesh"], horizon=par["horizon"],
                            nonlinear=True, full_height_rebar=bool(par.get("full_height_rebar")),
                            cyclic=bool(par.get("cyclic")),
                            cyclic_law=par.get("cyclic_law", "aydin"),
                            gf_factor=par.get("gf_factor", 1.0))

    # vectorized strain: precompute topology once
    ids = frames[0][0]
    pos = {int(n): k for k, n in enumerate(ids)}
    els = [e for e in model.elements if len(e.nodes) == 2]
    ia = np.array([pos[e.nodes[0]] for e in els])
    ib = np.array([pos[e.nodes[1]] for e in els])
    pts = np.array([model.nodes[int(n)].coords for n in ids], float)
    vec = pts[ib] - pts[ia]
    lsq = np.einsum("ij,ij->i", vec, vec)
    segs = np.stack([pts[ia], pts[ib]], axis=1)

    def strain(u):
        return np.einsum("ij,ij->i", u[ib] - u[ia], vec) / lsq

    strains = [strain(u) for _i, u, _m in frames]
    metas = np.array([m for _i, _u, m in frames])
    frac = np.array([(st > EPS_CRACK).mean() for st in strains])
    drift = metas[:, 1] / HW * 100

    # --- reversals, from the hysteresis log (finer than the frame stride) ---
    z = np.load(run / "hysteresis.npz")
    hu, hs = z["disp"], z["shear"] / 1e3
    turns = np.nonzero(np.diff(np.sign(np.diff(hu))))[0] + 1
    rev_u = hu[turns]

    # --- FIGURE 1: the frame nearest each reversal ---
    pick = [int(np.argmin(np.abs(metas[:, 1] - ru))) for ru in rev_u]
    vmax = max(float(np.abs(strains[k]).max()) for k in pick) if pick else 1.0
    norm = matplotlib.colors.SymLogNorm(linthresh=EPS_CRACK, vmin=-vmax, vmax=vmax, base=10)
    n = len(pick); cols = 2; rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5.2 * cols, 5.0 * rows + 1.4), squeeze=False)
    for k in range(rows * cols):
        ax = axes[k // cols][k % cols]
        if k >= n:
            ax.axis("off"); continue
        i = pick[k]; st = strains[i]
        lc = LineCollection(list(segs), cmap="coolwarm", norm=norm, linewidths=.35)
        lc.set_array(st); ax.add_collection(lc)
        ax.set_xlim(pts[:, 0].min(), pts[:, 0].max()); ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"reversal {k + 1} — {metas[i, 1]:+.2f} mm, {metas[i, 2] / 1e3:+,.0f} kN\n"
                     f"{(st > EPS_CRACK).sum():,} cracked ({frac[i]:.1%})", fontsize=12)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap="coolwarm"), ax=axes,
                 fraction=.022, pad=.015, label="axial strain (compression < 0 < tension)")
    fig.suptitle("Damage at each loop tip — the worst instant of every half-cycle",
                 fontsize=14)
    fig.savefig(run / "damage_reversals.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # --- FIGURE 2: accumulation as a curve ---
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.4, 8.6))
    a1.plot(drift, frac * 100, "-", color="#c0392b", lw=1.8)
    a1.plot(drift, frac * 100, ".", color="#c0392b", ms=3)
    for ru in rev_u:
        a1.axvline(ru / HW * 100, color="#16a085", ls=":", lw=1.0, alpha=.7)
    a1.set_xlabel("drift (%)", fontsize=12); a1.set_ylabel("struts cracked (%)", fontsize=12)
    a1.set_title("damage vs drift — the loops trace back over themselves", fontsize=12)
    a1.grid(alpha=.3, lw=.6)

    a2.plot(metas[:, 0] / 1e6, frac * 100, "-", color="#c0392b", lw=2.0)
    for t_ in turns:
        a2.axvline(0, color="none")
    a2.set_xlabel("step (millions)", fontsize=12); a2.set_ylabel("struts cracked (%)", fontsize=12)
    a2.set_title("damage vs time — monotonic, and it never heals", fontsize=12)
    a2.grid(alpha=.3, lw=.6)
    fig.suptitle(f"Damage accumulation — {len(frames)} frames, "
                 f"{frac[0]:.1%} to {frac[-1]:.1%} cracked", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(run / "damage_evolution.png", dpi=140)
    print(f"{len(frames)} frames | cracked {frac[0]:.1%} -> {frac[-1]:.1%} | "
          f"{n} reversal panels")
    print(f"saved {run / 'damage_reversals.png'}")
    print(f"saved {run / 'damage_evolution.png'}")


if __name__ == "__main__":
    main()
