"""Overlay the VK3 lattice run against the digitized test — the study's headline comparison.

Reads a saved run's JSON and `data/vk3_fig513.npz`, and reports the three comparisons VK3 makes
possible, in decreasing order of how fair they are:

  1. STRENGTH AND STIFFNESS — peak base shear each way, the backbone level by level, and the
     initial stiffness. Fair: these follow from section equilibrium and the elastic calibration.
  2. THE SHEAR / FLEXURE SPLIT — the model's diagonal-gauge shear share against the chapter's
     measured 20-22% (Fig. 5.19 right). Fair, and unique to this specimen among the three walls.
  3. CYCLIC DEGRADATION — the 2nd-cycle / 1st-cycle strength ratio level by level. Fair up to about
     1.27% drift and NOT fair beyond it: the test degraded because bars buckled, the cover spalled
     and the hoops' 90-degree hooks could open, none of which the model has.

Run from src/:  python examples/vk3_wall/compare.py [--run vk3_cyclic_dynamic]
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from rclattice import viz

from specimen import A_SHEAR, MEASURED_LOOPS, OUT, SPECIMEN, find_run
import testdata as td


def envelope(disp, shear):
    """Loop-tip envelope of a run: `(push, pull)` as lists of (drift %, kN)."""
    pos, neg, hi, lo = [], [], 0.0, 0.0
    for u, s in zip(disp, shear):
        if u > hi:
            hi = u
            pos.append((u, s))
        elif u < lo:
            lo = u
            neg.append((u, s))
    return pos, neg


def main(name: str | None, kind: str) -> None:
    out = find_run(kind, name=name)
    path = out / "data.json"
    if not path.exists():
        raise SystemExit(f"{out.name} has no data.json — it is not a completed analysis run")
    if not MEASURED_LOOPS.exists():
        raise SystemExit("no digitized test data — run digitize.py first")
    d = json.loads(path.read_text())
    m = np.load(MEASURED_LOOPS)

    drift, shear = d["drift_pct"], d["shear_kN"]
    pos, neg = envelope(drift, shear)
    peak_p = max((s for _u, s in pos), default=0.0)
    peak_n = min((s for _u, s in neg), default=0.0)
    tp, tn = float(m["peak_push_kN"]), float(m["peak_pull_kN"])

    print(f"{SPECIMEN}: lattice ({out.name}) vs the digitized test\n")
    print("1. STRENGTH — fair (section equilibrium)")
    print(f"  peak base shear   model +{peak_p:7.1f} / {peak_n:7.1f} kN")
    print(f"                    test  +{tp:7.1f} / {tn:7.1f} kN")
    print(f"                    ratio  {peak_p / tp:7.3f} / {peak_n / tn:7.3f}")
    print(f"  chapter's predicted F_n = {td.PREDICTION['Fn_kN']:.0f} kN "
          f"-> test/F_n = {tp / td.PREDICTION['Fn_kN']:.3f}, "
          f"model/F_n = {peak_p / td.PREDICTION['Fn_kN']:.3f}")

    print("\n2. BACKBONE, level by level (1st cycles)")
    lvl, up, dn = m["bb_level_mm"], m["bb_push_kN"], m["bb_pull_kN"]
    print(f"  {'drift':>7}  {'test push':>10}  {'model push':>11}  {'ratio':>7}   "
          f"{'test pull':>10}  {'model pull':>11}  {'ratio':>7}")
    for a, b, c in zip(lvl, up, dn):
        dr = float(a) / A_SHEAR * 100.0
        mp = next((s for u, s in pos if u >= dr), None)
        mn = next((s for u, s in neg if u <= -dr), None)
        f1 = f"{mp / float(b):7.3f}" if (mp is not None and b == b) else f"{'--':>7}"
        f2 = f"{mn / float(c):7.3f}" if (mn is not None and c == c) else f"{'--':>7}"
        s1 = f"{mp:11.1f}" if mp is not None else f"{'--':>11}"
        s2 = f"{mn:11.1f}" if mn is not None else f"{'--':>11}"
        print(f"  {dr:6.2f}%  {float(b):10.1f}  {s1}  {f1}   {float(c):10.1f}  {s2}  {f2}")

    g = d.get("gauge") or {}
    if g.get("components"):
        lo, hi = td.COMPONENTS["shear_pct_range"]
        shares = []
        for top, c in zip(g["top_mm"], g["components"]):
            if abs(top) > 2.0:
                shares.append(abs(c["shear"] / top) * 100.0)
        if shares:
            print(f"\n3. SHEAR SHARE — fair, and unique to this specimen")
            print(f"  model  {min(shares):.1f}-{max(shares):.1f}% of the top displacement "
                  f"(median {sorted(shares)[len(shares) // 2]:.1f}%)")
            print(f"  test   {lo:.0f}-{hi:.0f}% (Fig. 5.19 right; the highest of the three units)")

    if "bb2_push_kN" in m:
        print("\n4. CYCLIC DEGRADATION, 2nd cycle / 1st cycle (test only above ~1.27% drift is "
              "NOT a fair model comparison — bar buckling and hook opening drive it)")
        print(f"  {'drift':>7}  {'test 2nd/1st':>13}")
        for a, b, c in zip(lvl, up, m["bb2_push_kN"]):
            if b == b and c == c:
                print(f"  {float(a) / A_SHEAR * 100:6.2f}%  {float(c) / float(b):13.3f}")

    exp_d = list(-lvl[::-1] / A_SHEAR * 100.0) + list(lvl / A_SHEAR * 100.0)
    exp_s = list(dn[::-1]) + list(up)
    viz.figure_pushover(
        [{"disp": [u for u, _s in neg][::-1] + [u for u, _s in pos],
          "shear": [s for _u, s in neg][::-1] + [s for _u, s in pos],
          "label": f"lattice envelope (peak {peak_p:.0f} kN)", "style": {"color": "C0", "lw": 2}},
         {"disp": exp_d, "shear": exp_s,
          "label": f"{SPECIMEN} measured backbone (peak {tp:.0f} kN)",
          "style": {"color": "C3", "ls": "--", "lw": 2, "marker": "o", "ms": 5}}],
        savepath=str(out / "vs_test.png"), xlabel="drift ratio (%)",
        ylabel="base shear (kN)",
        title=f"{SPECIMEN} — lattice vs test backbone")
    print(f"\nsaved {out / 'vs_test.png'}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Compare a VK3 run against the digitized test")
    p.add_argument("--run", default=None,
                   help="run directory name under examples/output/vk3_wall/runs "
                        "(default: the newest of --kind)")
    p.add_argument("--kind", default="cyclic", choices=("cyclic", "pushover"))
    a = p.parse_args()
    main(a.run, a.kind)
