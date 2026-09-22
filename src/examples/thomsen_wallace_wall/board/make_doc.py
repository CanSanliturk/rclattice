"""Generate the RW2 status doc widget: ONE figure per run.

Each panel = grey cyclic loops (background) + test backbone + that ONE run's force-deformation
curve (base shear vs top displacement). Finished and running runs get the same treatment; the
running ones simply have shorter curves that grow each refresh.

Emits the widget `code` string (a <div> of <claude.Visualize> tags, single-quoted so it drops into
the Docs JSON with no escaping) to stdout. Push it with the Claude Docs `update` tool.

    uv run python examples/thomsen_wallace_wall/board/make_doc.py > /tmp/doc_widget.js
"""
from __future__ import annotations
import glob
import re

import numpy as np

H = 3660.0
ROOT = "examples/output/thomsen_wallace_wall/study/"

# (glob, label, done-substring-or-None). Panels render in this order.
RUNS = [
    ("2026-09-20_174437*", "S1 uniform 30.5 - fy414, crushing"),
    ("2026-09-21_085207*", "S2 fy454 - crushing"),
    ("2026-09-21_085208*", "S2 fy454 + EPP - crushing"),
    ("2026-09-22_013850*", "S2 fy454 - comp=linear (Aydin)"),
    ("2026-09-20_174435*", "S1 graded 25 - fy414, crushing"),
]


def samples(pat: str):
    g = sorted(glob.glob(ROOT + pat))
    if not g:
        return None, False
    txt = open(g[-1] + "/console.log").read()
    pts = [(round(float(m.group(1)) / 100 * H, 1), round(float(m.group(2)), 1))
           for m in re.finditer(r"drift\s+\+([\d.]+)%\s+shear\s+\+?([-\d.]+) kN", txt)]
    done = "traced to" in txt
    return (pts or None), done


def peak_note(pts):
    if not pts:
        return "no samples yet"
    smax = max(s for _, s in pts)
    return f"reached {pts[-1][0] / H * 100:.2f} pct, max sample {smax:.0f} kN / test 163.3 = {smax / 163.3:.3f}x"


def spos(x, y, n=None):
    x = np.asarray(x); y = np.asarray(y)
    o = np.argsort(x); x, y = x[o], y[o]; m = x >= 0; x, y = x[m], y[m]
    if n:
        i = np.linspace(0, len(x) - 1, min(n, len(x))).round().astype(int)
        x, y = x[i], y[i]
    return list(zip(np.round(x, 1), np.round(y, 1)))


def js(points):
    return "[" + ",".join("{disp:%g,shear:%g}" % (a, b) for a, b in points) + "]"


def main():
    d = np.load("examples/thomsen_wallace_wall/data/fig9b.npz")
    both = d["cloud"]
    pos = both[both[:, 0] >= 0]    # positive-direction loops only (the pushover direction)
    st = max(1, len(pos) // 35)    # sparse: it is only background, and keeps each refresh light
    cloud = [(round(float(pos[i, 0]), 1), round(float(pos[i, 1]), 1)) for i in range(0, len(pos), st)]
    bb = spos(d["mono_x"], d["mono_y"], 26)
    cloud_js, bb_js = js(cloud), js(bb)

    tags = []
    for k, (pat, label) in enumerate(RUNS):
        pts, done = samples(pat)
        if not pts:
            continue
        run_js = js([(dd, vv) for dd, vv in pts])
        status = "DONE" if done else "RUNNING - live"
        title = label[:40]                      # chart title hard cap is 40 chars
        note = f"{status}. {peak_note(pts)}"
        tags.append(
            "<claude.Visualize data-claude-component='r%d' "
            "sources={{cloud:{kind:'data', data:%s}, bb:{kind:'data', data:%s}, run:{kind:'data', data:%s}}} "
            "definition={{kind:'chart', source:'run', title:'%s', note:'%s', "
            "axes:{x:[{column:'disp', type:'linear', title:'Top displacement (mm)', format:{kind:'number', suffix:' mm'}}], "
            "y:[{title:'Base shear (kN)', format:{kind:'number', spec:'~s'}}]}, "
            "marks:[{type:'point', y:'shear', source:'cloud', title:'Test cyclic loops', tone:{column:'shear', cutoffs:[-100000,100000]}}, "
            "{type:'line', y:'shear', source:'bb', title:'Test backbone'}, "
            "{type:'line', y:'shear', source:'run', title:'model'}]}} />"
            % (k, cloud_js, bb_js, run_js, title, note)
        )
    code = "export default () => <div>" + "".join(tags) + "</div>;"
    print(code)


if __name__ == "__main__":
    main()
