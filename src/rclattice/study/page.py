"""The advisor-facing HTML page and the JSON payload it renders itself from (D103).

The page is an advisor-facing deliverable on a PUBLICLY SHARED artifact, and a public artifact
cannot be granted the `db` capability — so it cannot fetch anything at load time and instead
carries its own data. That makes rebuilding it a step in the study's own workflow rather than a
one-off piece of hand-written HTML:

    python master.py && python report_page.py
    # then republish the emitted file to the same artifact URL

NOTHING HERE IS HAND-MAINTAINED. Rows, ratios and curves are all derived from the run directories,
so a run that lands after this was written appears without anyone editing the page. The template
(`spec.page.template`) holds only prose and layout, with `__PAYLOAD__` and `__REFDATA__` slots.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import master, metrics
from .spec import StudySpec

ORDER = ["elastic", "static", "pushover", "cyclic"]


def curve(data: dict, n: int = 160) -> list:
    """The run's smoothed response, resampled for the page (drift %, kN)."""
    sh = np.asarray(data["shear"], float)
    dp = np.asarray(data["disp"], float)
    h = float(data["meta"]["drift_denominator_mm"])
    sm = metrics.moving_average(sh, int(data.get("peak_smoothing_samples") or 1))
    i = np.linspace(0, len(sm) - 1, min(n, len(sm))).astype(int)
    return [[round(float(dp[k] / h * 100), 4), round(float(sm[k] / 1e3), 1)] for k in i]


def loops(data: dict, n: int = 2600) -> list:
    """A cyclic run's RAW path — smoothing a hysteresis would round off its corners."""
    sh = np.asarray(data["shear"], float)
    dp = np.asarray(data["disp"], float)
    h = float(data["meta"]["drift_denominator_mm"])
    i = np.linspace(0, len(dp) - 1, min(n, len(dp))).astype(int)
    return [[round(float(dp[k] / h * 100), 4), round(float(sh[k] / 1e3), 1)] for k in i]


def payload(spec: StudySpec, root: Path) -> dict:
    have = master.runs(spec, root)
    index: dict[tuple, list] = {}
    for r in have:
        index.setdefault(spec.cell(r["params"]), []).append(r)

    axes = spec.axes
    pg = spec.page
    first = pg.comp_order if pg and pg.comp_order else spec.registry.by_name[axes[0]].choices
    rest = [spec.registry.by_name[a].choices for a in axes[1:]]
    from itertools import product

    rows, series, best = [], {}, None
    for cell in product(first, *rest):
        got = index.get(tuple(cell)) or []
        by_kind: dict[str, list] = {}
        for r in got:
            by_kind.setdefault(r["params"]["analysis"], []).append(r)
        for kind in sorted(by_kind, key=lambda k: (ORDER.index(k) if k in ORDER else 99, k)):
            rec = max(by_kind[kind], key=lambda r: r["dir"].name)
            d = rec["data"]
            peak, at_drift, smoothed = metrics.headline_peak(d)
            if not peak:
                continue
            key = rec["dir"].name
            series[key] = loops(d) if kind == "cyclic" else curve(d)
            ratio = peak / 1e3 / spec.measured_F_kN
            end = abs(d.get("end_drift") or 0.0)
            amp = max(abs(u) for u in d["disp"]) / float(d["meta"]["drift_denominator_mm"])
            cyclic = amp > end * 1.5
            row = {a: v for a, v in zip(axes, cell)}
            row.update({"analysis": kind,
                        "peak_kN": round(peak / 1e3, 1), "ratio": round(ratio, 3),
                        "at_drift": at_drift, "reach": amp if cyclic else end,
                        "cyclic": cyclic, "run": key, "series": key,
                        "xmax": round(max(amp, end) * 100, 3),
                        "smoothed": bool(smoothed), "n_of_kind": len(by_kind[kind]),
                        "variant": master.variant_note(spec, rec["params"]),
                        # A peak below 0.1% drift on a bonded model is the START-UP TRANSIENT,
                        # not sustained resistance — real force, but not strength.
                        "spike": bool(spec.is_bonded(rec["params"]) and at_drift < 0.001)})
            if row["spike"]:
                w = np.asarray(d["disp"], float) / float(d["meta"]["drift_denominator_mm"])
                sm = metrics.moving_average(np.asarray(d["shear"], float),
                                            int(d.get("peak_smoothing_samples") or 1))
                m = w >= 0.001
                row["sustained"] = [round(float(sm[m].max() / 1e3), 1),
                                    round(float(sm[-1] / 1e3), 1)]
            rows.append(row)
            if kind == "pushover" and (best is None or abs(ratio - 1) < abs(best - 1)):
                best = ratio
    for r in rows:
        r["best"] = (r["analysis"] == "pushover" and best is not None
                     and r["ratio"] == round(best, 3))

    counts: dict[str, int] = {}
    for r in have:
        counts[r["params"]["analysis"]] = counts.get(r["params"]["analysis"], 0) + 1
    return {"generated": __import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M"),
            "run_count": len(have), "counts": counts, "rows": rows, "series": series,
            "comp_note": (pg.comp_note if pg else {}), "measured_kN": spec.measured_F_kN}


def reference_data(spec: StudySpec, height: float) -> dict:
    """The digitized experiment — fixed data, so it is baked in and never regenerated."""
    d = spec.references.load()
    if d is None:
        return {"cloud": [], "env": []}
    cloud = np.asarray(d["cloud"], float)
    mx, my = np.asarray(d["mono_x"], float), np.asarray(d["mono_y"], float)
    order = np.argsort(mx)
    return {"cloud": [[round(float(a / height * 100), 4), round(float(b), 1)] for a, b in cloud[::5]],
            "env": [[round(float(mx[k] / height * 100), 4), round(float(my[k]), 1)] for k in order]}


def write_payload(spec: StudySpec) -> dict:
    pl = payload(spec, spec.out_root)
    out = spec.page.out_payload
    out.write_text(json.dumps(pl))
    print(f"{out}  {out.stat().st_size:,} bytes")
    print(f"  {pl['run_count']} runs, {len(pl['rows'])} rows with a peak, "
          f"{len(pl['series'])} curves")
    for r in pl["rows"]:
        print("    " + " ".join(f"{r[a]:9s}" for a in spec.axes)
              + f" {r['analysis']:9s} {r['peak_kN']:>8,.1f} kN  {r['ratio']:.3f} x"
              + ("  SPIKE" if r["spike"] else "") + ("  <-- best" if r["best"] else ""))
    return pl


def render(spec: StudySpec) -> Path:
    pg = spec.page
    pl = payload(spec, spec.out_root)
    height = spec.drift_denominator({p.name: p.default for p in spec.registry})
    html = (pg.template.read_text(encoding="utf8")
            .replace("__PAYLOAD__", json.dumps(pl))
            .replace("__REFDATA__", json.dumps(reference_data(spec, height))))
    pg.out_html.write_text(html, encoding="utf8")
    print(f"{pg.out_html}  {pg.out_html.stat().st_size:,} bytes")
    print(f"  {pl['run_count']} runs, {len(pl['rows'])} rows, "
          f"{len(pl['series'])} curves, generated {pl['generated']}")
    return pg.out_html
