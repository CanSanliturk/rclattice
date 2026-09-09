"""The matrix-aware master report (PLAN.md §10).

Builds the matrix by READING `params.json`, never by parsing directory names, which is what lets a
new parameter or a new value extend the matrix without touching this file. Cells that have not been
run are shown as GAPS rather than omitted, so the matrix reads as a plan and a result at once.
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import params as P    # noqa: E402
import metrics

AXES = ("comp", "tail", "bond")


def runs(root: Path) -> list[dict]:
    """Every run directory that has both a parameter record and a result."""
    out = []
    for d in sorted(root.glob("*")):
        pj, dj = d / "params.json", d / "data.json"
        if not (pj.is_file() and dj.is_file()):
            continue
        params = json.loads(pj.read_text())
        # Runs made before a parameter existed stay readable: fill its default and say so.
        missing = [p.name for p in P.REGISTRY if p.name not in params]
        for name in missing:
            params[name] = P.BY_NAME[name].default
        P.normalize(params)      # pre-D94 comp/bond spellings, so old runs share a cell with new
        out.append({"dir": d, "params": params, "data": json.loads(dj.read_text()),
                    "filled_defaults": missing})
    return out


def cell_row(params, rec) -> list:
    if rec is None:
        return ["—", "—", "—", "—", "not run"]
    d = rec["data"]
    # The headline peak is the SMOOTHED one wherever the run has it (D92) — a raw sample maximum
    # reports the ringing of a released bar, not resistance. A `*` marks a cell still quoting a raw
    # maximum that `rescore.py` has not corrected.
    peak, peak_drift, smoothed = metrics.headline_peak(d)
    peak = peak or None
    mark = "" if smoothed or not peak else " \\*"
    return [
        f"{peak / 1e3:,.1f}{mark}" if peak else "—",
        f"{peak / 1e3 / 963.592:.3f}" if peak else "—",
        f"{peak_drift:.3%}" if peak else "—",
        f"{d.get('end_drift', 0):.3%}" if d.get("end_drift") else "—",
        f"[{rec['dir'].name}]({rec['dir'].name}/report.md)",
    ]


def cross_run_figures(root: Path, have: list[dict]) -> list[str]:
    """One figure per AXIS, overlaying the runs that differ only in that axis (PLAN §10).

    Comparability is decided by `params.comparable_key`, which excludes report-only parameters — so
    adding a figure or a probe never splits a comparison in two, and a new axis needs no edit here.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = []
    figs = root / "figures"
    for axis in AXES + ("panel", "tw", "mesh", "horizon", "gf"):
        groups: dict[tuple, list[dict]] = {}
        for r in have:
            if not r["data"].get("disp"):
                continue
            groups.setdefault(P.comparable_key(r["params"], axis=axis), []).append(r)
        families = [g for g in groups.values() if len({r["params"][axis] for r in g}) > 1]
        if not families:
            continue
        figs.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=170)
        for family in families:
            for r in sorted(family, key=lambda r: str(r["params"][axis])):
                d, h = r["data"], r["data"]["meta"]["drift_denominator_mm"]
                ax.plot([u / h * 100.0 for u in d["disp"]], [s / 1e3 for s in d["shear"]],
                        lw=1.8, label=f"{axis}={r['params'][axis]}")
        ax.axhline(963.592, color="C3", ls="--", lw=1.2, label="test peak (Table 4)")
        ax.set_xlabel("drift ratio (%)")
        ax.set_ylabel("base shear (kN)")
        ax.set_title(f"Varying {axis}, everything else held")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        name = f"by_{axis}.png"
        fig.savefig(figs / name)
        plt.close(fig)
        made.append(name)
    return made


def build(root: Path) -> str:
    have = runs(root)
    index = {}
    for r in have:
        index.setdefault(tuple(r["params"][a] for a in AXES), []).append(r)

    values = {a: P.BY_NAME[a].choices for a in AXES}
    lines = ["# Aldemir wall — parametric study, master report", "",
             f"{len(have)} run(s) found under `{root}`. Cells with no run are gaps, not omissions.",
             "", "## The matrix", "",
             "| comp | tail | bond | peak (kN) | /test | drift@peak | traced to | run |",
             "|---|---|---|---|---|---|---|---|"]
    for comp, tail, bond in product(*(values[a] for a in AXES)):
        got = index.get((comp, tail, bond))
        newest = max(got, key=lambda r: r["dir"].name) if got else None
        lines.append("| " + " | ".join([comp, tail, bond] + cell_row(None, newest)) + " |")

    lines += ["", "## Analyses present", ""]
    by_analysis: dict[str, int] = {}
    for r in have:
        by_analysis[r["params"]["analysis"]] = by_analysis.get(r["params"]["analysis"], 0) + 1
    for k, v in sorted(by_analysis.items()):
        lines.append(f"* {k}: {v}")

    stale = [r for r in have if r["filled_defaults"]]
    if stale:
        lines += ["", "## Runs predating a parameter", "",
                  "Their `params.json` is older than the registry; the missing values are shown "
                  "filled with today's defaults.", ""]
        for r in stale:
            lines.append(f"* `{r['dir'].name}` — filled: {', '.join(r['filled_defaults'])}")

    made = cross_run_figures(root, have)
    lines += ["", "## Cross-run figures", ""]
    lines += ([f"* `figures/{f}`" for f in made] if made else
              ["No axis yet has two runs that differ in it alone, so there is nothing to overlay."])

    lines += ["", "## Comparability", "",
              "Two runs answer the same question when every model- and analysis-affecting "
              "parameter matches except the axis under study; report-only parameters are excluded, "
              "so adding a figure never splits a comparison in two."]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=None)
    a = ap.parse_args()
    import os
    os.environ.setdefault("ALDEMIR_TW", "210.0")
    import specimen
    root = a.root or (specimen.OUT / "study")
    root.mkdir(parents=True, exist_ok=True)
    text = build(root)
    (root / "master_report.md").write_text(text)
    print(text)
    print(f"written to {root / 'master_report.md'}")


if __name__ == "__main__":
    main()
