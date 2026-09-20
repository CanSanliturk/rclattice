"""The matrix-aware master report (Aldemir PLAN.md §10; shared in D103).

Builds the matrix by READING `params.json`, never by parsing directory names, which is what lets a
new parameter or a new value extend the matrix without touching this file. Cells that have not been
run are shown as GAPS rather than omitted, so the matrix reads as a plan and a result at once.

**A matrix that reduces N runs to one cell must always say which run it picked and how it
differs** — `variant_note` — because the shadowing defect was found twice (D92, D98).
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path

from . import metrics
from .spec import StudySpec


def runs(spec: StudySpec, root: Path) -> list[dict]:
    """Every run directory that has both a parameter record and a result."""
    out = []
    for d in sorted(root.glob("*")):
        pj, dj = d / "params.json", d / "data.json"
        if not (pj.is_file() and dj.is_file()):
            continue
        params = json.loads(pj.read_text())
        # Runs made before a parameter existed stay readable: fill its default and say so.
        missing = spec.registry.fill_defaults(params)
        out.append({"dir": d, "params": params, "data": json.loads(dj.read_text()),
                    "filled_defaults": missing})
    return out


def variant_note(spec: StudySpec, params: dict) -> str:
    """`mesh 25, damping 0.05` — how the shown run differs from the registry defaults.

    Deliberately a SHORT list (`spec.variant_keys`): `drift` and the protocol already appear in the
    run name, and listing everything would bury the signal.
    """
    out = []
    for name in spec.variant_keys:
        p = spec.registry.by_name.get(name)
        if p is None or name not in params:
            continue
        v, d = params[name], p.default
        if v == d:
            continue
        if isinstance(v, bool):
            out.append(name.replace("_", " "))
        elif isinstance(v, float):
            out.append(f"{name} {v:g}")
        else:
            out.append(f"{name} {v}")
    return ", ".join(out)


def _reach(d: dict) -> str:
    """How far the run actually took the model.

    `end_drift` is the LAST drift, which is right for a monotonic push and wrong for a cyclic
    history: that ends back at the origin and reads as "-0.000%", as though the run went nowhere.
    For a cyclic run the honest reach is the largest amplitude it reached.
    """
    disp = d.get("disp")
    height = (d.get("meta") or {}).get("drift_denominator_mm")
    if disp and height:
        peak_amp = max(abs(u) for u in disp) / float(height)
        end = abs(d.get("end_drift") or 0.0)
        if peak_amp > end * 1.5:                      # a history that came back — quote the amplitude
            return f"±{peak_amp:.3%}"
    return f"{d.get('end_drift', 0):.3%}" if d.get("end_drift") else "—"


def cell_row(spec: StudySpec, rec) -> list:
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
        f"{peak / 1e3 / spec.measured_F_kN:.3f}" if peak else "—",
        f"{peak_drift:.3%}" if peak else "—",
        _reach(d),
        f"[{rec['dir'].name}]({rec['dir'].name}/report.md)",
    ]


def cross_run_figures(spec: StudySpec, root: Path, have: list[dict]) -> list[str]:
    """One figure per AXIS, overlaying the runs that differ only in that axis (PLAN §10).

    Comparability is decided by `registry.comparable_key`, which excludes report-only parameters
    — so adding a figure or a probe never splits a comparison in two, and a new axis needs no edit
    here.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = []
    figs = root / "figures"
    for axis in spec.axes + spec.cross_run_axes:
        groups: dict[tuple, list[dict]] = {}
        for r in have:
            if not r["data"].get("disp"):
                continue
            groups.setdefault(spec.registry.comparable_key(r["params"], axis=axis), []).append(r)
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
        ax.axhline(spec.measured_F_kN, color="C3", ls="--", lw=1.2,
                   label=f"test peak ({spec.references.f_source})")
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


ORDER = ["elastic", "static", "pushover", "cyclic"]


def build(spec: StudySpec, root: Path) -> str:
    have = runs(spec, root)
    index = {}
    for r in have:
        index.setdefault(spec.cell(r["params"]), []).append(r)

    axes = spec.axes
    values = {a: spec.registry.by_name[a].choices for a in axes}
    head = "| " + " | ".join(axes) + " | analysis | peak (kN) | /test | drift@peak | traced to | run |"
    sep = "|" + "---|" * (len(axes) + 6)
    lines = [f"# {spec.title} — parametric study, master report", "",
             f"{len(have)} run(s) found under `{root}`. Cells with no run are gaps, not omissions.",
             ""]
    if spec.matrix_note:
        lines += [spec.matrix_note, ""]
    lines += ["## The matrix", "", head, sep]
    for cell in product(*(values[a] for a in axes)):
        got = index.get(cell) or []
        by_kind: dict[str, list] = {}
        for r in got:
            by_kind.setdefault(r["params"]["analysis"], []).append(r)
        if not by_kind:
            lines.append("| " + " | ".join(list(cell) + ["—"] + cell_row(spec, None)) + " |")
            continue
        # A stable, meaningful order rather than alphabetical: how far each analysis takes the model.
        for kind in sorted(by_kind, key=lambda k: (ORDER.index(k) if k in ORDER else 99, k)):
            newest = max(by_kind[kind], key=lambda r: r["dir"].name)
            n = len(by_kind[kind])
            label = kind if n == 1 else f"{kind} ({n}, newest)"
            # The matrix axes are few, so ANY other parameter shadows: the newest run in a cell may
            # be a mesh or damping variant standing where the baseline used to be. Name what is
            # non-default about the run actually shown (D98).
            odd = variant_note(spec, newest["params"])
            if odd:
                label += f" — {odd}"
            lines.append("| " + " | ".join(list(cell) + [label] + cell_row(spec, newest)) + " |")

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

    made = cross_run_figures(spec, root, have)
    lines += ["", "## Cross-run figures", ""]
    lines += ([f"* `figures/{f}`" for f in made] if made else
              ["No axis yet has two runs that differ in it alone, so there is nothing to overlay."])

    lines += ["", "## Comparability", "",
              "Two runs answer the same question when every model- and analysis-affecting "
              "parameter matches except the axis under study; report-only parameters are excluded, "
              "so adding a figure never splits a comparison in two."]
    return "\n".join(lines) + "\n"


def main(spec: StudySpec, argv=None) -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=None)
    a = ap.parse_args(argv)
    root = a.root or spec.out_root
    root.mkdir(parents=True, exist_ok=True)
    text = build(spec, root)
    (root / "master_report.md").write_text(text)
    print(text)
    print(f"written to {root / 'master_report.md'}")
