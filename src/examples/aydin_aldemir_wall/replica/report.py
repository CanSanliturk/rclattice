"""ONE portrait image with both the hysteresis and the damage panels — for reading on a phone.

    uv run python examples/aydin_aldemir_wall/replica/report.py [--run <dir>] [--panels 4]

`hyst.py` and `damage_history.py` still write their own figures; this composes the same content
into a single file so one image carries the whole picture.
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
FIG10B = Path(__file__).resolve().parent.parent / "data" / "fig10b.npz"
HW = 2680.0
F_EXP = 963.592
EPS_CRACK = FT / EC


def load_frames(run: Path):
    fd, got = run / "frames", []
    if fd.exists():
        for p in sorted(fd.glob("*.npz"), key=lambda q: int(q.stem)):
            z = np.load(p); got.append((z["ids"], z["u"], z["meta"]))
    if not got and (run / "snapshots.npz").exists():
        z = np.load(run / "snapshots.npz", allow_pickle=True)
        got = [(z["ids"], u, m) for u, m in zip(z["ring_u"], z["ring_meta"])]
    return got


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, default=None)
    ap.add_argument("--panels", type=int, default=4)
    a = ap.parse_args()
    run = a.run or max((d for d in OUT.glob("*cyclic*") if d.is_dir()),
                       key=lambda d: d.stat().st_mtime)
    par = json.loads((run / "params.json").read_text())
    z = np.load(run / "hysteresis.npz"); u, s = z["disp"], z["shear"] / 1e3
    z10 = np.load(FIG10B) if FIG10B.exists() else None
    lv = [float(v) for v in str(par.get("cyclic", "")).split(",") if v]

    def mono_curve(pattern):
        hit = sorted(OUT.glob(pattern))
        if not hit:
            return None
        j = json.loads(hit[-1].read_text())
        return np.array(j["disp"]), np.array(j["shear"]) / 1e3

    ours15 = mono_curve("*_pushover_d0.003_rebartop_CentralDifference/data.json")
    ours301 = mono_curve("*_h3.01_rebartop_CentralDifference/data.json")

    # every OTHER cyclic run, so a status figure carries the whole series rather than only the
    # live one — the HystereticSM run is the controlled comparison for the Concrete02 one
    prev = []
    for d in sorted(OUT.glob("*cyclic*")):
        if d == run or not (d / "hysteresis.npz").exists():
            continue
        zz = np.load(d / "hysteresis.npz")
        pj = json.loads((d / "params.json").read_text()) if (d / "params.json").exists() else {}
        prev.append((pj.get("cyclic_law", "aydin"), zz["disp"], zz["shear"] / 1e3))

    fr = load_frames(run)
    if len(fr) > a.panels:
        fr = [fr[i] for i in np.unique(np.linspace(0, len(fr) - 1, a.panels).astype(int))]

    cal = calibrate(mesh_size=par["mesh"], horizon=par["horizon"])
    model, _ = wall_lattice(cal.area, mesh_size=par["mesh"], horizon=par["horizon"],
                            nonlinear=True, full_height_rebar=bool(par.get("full_height_rebar")),
                            cyclic=bool(par.get("cyclic")),
                            cyclic_law=par.get("cyclic_law", "aydin"),
                            gf_factor=par.get("gf_factor", 1.0))
    from rclattice import viz
    _i, idx, pts, _l, _q = viz._arrays(model)
    els = [e for e in model.elements if len(e.nodes) == 2]
    segs = np.array([[pts[idx[e.nodes[0]]], pts[idx[e.nodes[1]]]] for e in els])
    strains, metas = [], []
    for ids, uu, m in fr:
        d = {int(i): list(map(float, r)) for i, r in zip(ids, uu)}
        _e, st = viz.strut_strains(model, d); strains.append(st); metas.append(m)
    vmax = max(float(np.abs(x).max()) for x in strains) if strains else 1.0
    norm = matplotlib.colors.SymLogNorm(linthresh=EPS_CRACK, vmin=-vmax, vmax=vmax, base=10)

    prows = int(np.ceil(len(fr) / 2))
    fig = plt.figure(figsize=(8.6, 10.4 + 4.4 * prows))
    gs = fig.add_gridspec(2 + prows, 2, height_ratios=[6.2, 5.0] + [4.4] * prows, hspace=.52,
                          wspace=.06)

    def hyst(ax, zoom, legend=False, only_test_and_live=False):
        # TEST: the cloud is context, the envelope is the readable line. Drawn faintly and
        # underneath everything so seven series stay separable.
        if z10 is not None:
            cl = z10["cloud"]
            if zoom:
                cl = cl[np.abs(cl[:, 0]) <= zoom * 1.08]
            ax.plot(cl[:, 0], cl[:, 1], ".", ms=1.8 if zoom else 1.0, color="#b9cbc7",
                    alpha=.55, zorder=1, label="test — measured cloud")
            # RAW envelope, unsmoothed. A running extremum applied only beyond the pinch leaves a
            # visible seam where it starts, and applied from the origin it invents a plateau out of
            # loop-crossing scatter. Ragged is the honest option.
            bx, be = z10["backbone_x"], z10["backbone"]
            for half, sgn in ((0, 1), (1, -1)):        # upper for +x, lower for -x
                k = (bx > 0) if sgn > 0 else (bx < 0)
                ax.plot(bx[k], be[k, half], "-", color="#0f8a78", lw=2.1, zorder=3,
                        label="test — envelope" if sgn > 0 else None)

        # AYDIN's own lattice, both horizons, full range straight from his figure
        if not only_test_and_live:
            for key, hz, ls, lw, al in (("aydin15", "1.5", "-", 1.8, .95),
                                        ("aydin301", "3.01", "--", 1.5, .8)):
                ax.plot(z10[f"{key}_x"], z10[f"{key}_y"], ls, color="#7d3c98", lw=lw, alpha=al,
                        zorder=2, label=f"Aydin — h {hz}")

        # OUR monotonic pushovers, mirrored into the pull quadrant for comparison
        if not only_test_and_live:
            for cur, hz, ls in ((ours15, "1.5", "-"), (ours301, "3.01", "--")):
                if cur is None:
                    continue
                xx, yy = cur
                ax.plot(xx, yy, ls, color="#b9770e", lw=1.6, zorder=2,
                        label=f"ours pushover — h {hz}")
                ax.plot(-xx, -yy, ls, color="#b9770e", lw=1.2, alpha=.5, zorder=2)

            for law, pu, ps in prev:
                ax.plot(pu, ps, "-", color="#d98880", lw=1.3, alpha=.9, zorder=4,
                        label=f"ours cyclic — {law} (stopped)")
        ax.plot(u, s, "-", color="#c0392b", lw=2.0, zorder=5,
                label=f"ours CYCLIC — {par.get('cyclic_law', 'aydin')}")
        ax.plot(u[-1], s[-1], "o", ms=8.5, color="#c0392b", zorder=6)
        ax.axhline(0, color="0.65", lw=.8); ax.axvline(0, color="0.65", lw=.8)
        for d_ in lv:
            for sg in (1, -1):
                ax.axvline(sg * d_ * HW, color="#16a085", ls=":", lw=.9, alpha=.45)
        ax.set_ylabel("base shear (kN)", fontsize=13); ax.grid(alpha=.25, lw=.6)
        ax.tick_params(labelsize=11)
        if legend == "below":
            ax.legend(fontsize=8.8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
                      frameon=False, columnspacing=1.4, handlelength=2.0)
        elif legend:
            ax.legend(fontsize=9, loc="upper left", framealpha=.96,
                      columnspacing=1.0, handlelength=1.9, borderpad=.6)

    ax0 = fig.add_subplot(gs[0, :]); hyst(ax0, None, legend="below")
    ax0.set_xlim(-20, 18); ax0.set_ylim(-1150, 1150)
    ax0.set_title("full test range — everything on one axis", fontsize=13)
    span = max(max(abs(u)) * 1.35, (max(lv) * HW * 1.15 if lv else 1.0))
    ax1 = fig.add_subplot(gs[1, :]); hyst(ax1, span, legend=True,
                                         only_test_and_live=True)
    ax1.set_xlim(-span, span)
    # y limit from EVERY series inside the window, not just the live run — early in a run its own
    # range is a few tens of kN and everything else falls off the panel.
    cand = [abs(s).max()]
    if z10 is not None:
        cl = z10["cloud"]; k = np.abs(cl[:, 0]) <= span
        if k.any():
            cand.append(np.abs(cl[k, 1]).max())
    ylim = max(cand) * 1.15
    ax1.set_ylim(-ylim, ylim)
    ax1.set_xlabel("top displacement (mm)", fontsize=13)
    ax1.set_title("zoom — measured test vs this run only", fontsize=13)

    for k, (st, m) in enumerate(zip(strains, metas)):
        ax = fig.add_subplot(gs[2 + k // 2, k % 2])
        lc = LineCollection(list(segs), cmap="coolwarm", norm=norm, linewidths=.35)
        lc.set_array(st); ax.add_collection(lc)
        ncr = int((st > EPS_CRACK).sum())
        ax.set_xlim(pts[:, 0].min(), pts[:, 0].max()); ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"step {int(m[0]):,} — {m[1] / HW * 100:+.4f}%, {m[2] / 1e3:+,.0f} kN\n"
                     f"{ncr:,} cracked ({ncr / len(st):.1%})", fontsize=11)

    fig.suptitle(f"Aldemir wall — cyclic, horizon {par.get('horizon')}, bars to top\n"
                 f"reached {max(abs(u)) / HW:.4%} drift of the test's ~0.71%"
                 + ("   |   damage panels share one symlog strain scale" if fr else
                    "   |   no damage frames yet"), fontsize=14, y=.995)
    dst = run / "report.png"
    fig.savefig(dst, dpi=140, bbox_inches="tight")
    print(f"saved {dst}  ({len(fr)} damage panels, {len(u):,} hysteresis points)")


if __name__ == "__main__":
    main()
