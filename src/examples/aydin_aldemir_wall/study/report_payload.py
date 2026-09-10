"""Build the JSON the published master-report page renders itself from.

The page is deployed once with a baked-in copy of this payload as its fallback; the live copy
lives in the artifact's `db` under `report/latest`. Re-run this after `master.py` and push the
result, and a reader's refresh shows the new state without the page being republished:

    python master.py && python report_payload.py
    # then, from a Claude session:  Artifact write_db  collection=report doc_id=latest
    #                               file_path=.../report_payload.json

NOTHING HERE IS HAND-MAINTAINED. Rows, ratios and curves are all derived from the run directories,
so a run that lands after this was written appears without anyone editing the page.

The digitized experiment is deliberately NOT in the payload: it is fixed data that never changes,
so it stays baked into the page and costs the store nothing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import master                                                             # noqa: E402
import metrics                                                            # noqa: E402
import params as P                                                        # noqa: E402
import references                                                         # noqa: E402

OUT = HERE / "report_payload.json"
ORDER = ["elastic", "static", "pushover", "cyclic"]
COMP_NOTE = {"crushing": "Concrete02: yields, softens, crushes",
             "capped": "elastic-perfectly-plastic: yields at f<sub>c</sub>, cannot crush",
             "linear": "Aydin's own law: compression never yields"}
COMP_ORDER = ["crushing", "capped", "linear"]


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


def build(root: Path) -> dict:
    have = master.runs(root)
    index: dict[tuple, list] = {}
    for r in have:
        index.setdefault(tuple(r["params"][a] for a in master.AXES), []).append(r)

    rows, series, best = [], {}, None
    for comp in COMP_ORDER:
        for tail in P.BY_NAME["tail"].choices:
            for bond in P.BY_NAME["bond"].choices:
                got = index.get((comp, tail, bond)) or []
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
                    ratio = peak / 1e3 / master.MEASURED if hasattr(master, "MEASURED") else \
                        peak / 1e3 / 963.592
                    end = abs(d.get("end_drift") or 0.0)
                    amp = max(abs(u) for u in d["disp"]) / float(d["meta"]["drift_denominator_mm"])
                    cyclic = amp > end * 1.5
                    row = {"comp": comp, "tail": tail, "bond": bond, "analysis": kind,
                           "peak_kN": round(peak / 1e3, 1), "ratio": round(ratio, 3),
                           "at_drift": at_drift, "reach": amp if cyclic else end,
                           "cyclic": cyclic, "run": key, "series": key,
                           "xmax": round(max(amp, end) * 100, 3),
                           "smoothed": bool(smoothed), "n_of_kind": len(by_kind[kind]),
                           "variant": master.variant_note(rec["params"]),
                           # A peak below 0.1% drift on a bonded model is the START-UP TRANSIENT,
                           # not sustained resistance — real force, but not strength.
                           "spike": bool(bond != "perfect" and at_drift < 0.001)}
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
        r["best"] = (r["analysis"] == "pushover" and r["ratio"] == round(best, 3))

    counts: dict[str, int] = {}
    for r in have:
        counts[r["params"]["analysis"]] = counts.get(r["params"]["analysis"], 0) + 1
    return {"generated": __import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M"),
            "run_count": len(have), "counts": counts, "rows": rows, "series": series,
            "comp_note": COMP_NOTE, "measured_kN": 963.592}


def main() -> None:
    import os
    os.environ.setdefault("ALDEMIR_TW", "210.0")
    import specimen
    payload = build(specimen.OUT / "study")
    OUT.write_text(json.dumps(payload))
    print(f"{OUT}  {OUT.stat().st_size:,} bytes")
    print(f"  {payload['run_count']} runs, {len(payload['rows'])} rows with a peak, "
          f"{len(payload['series'])} curves")
    for r in payload["rows"]:
        print(f"    {r['comp']:9s} {r['tail']:7s} {r['bond']:8s} {r['analysis']:9s} "
              f"{r['peak_kN']:>8,.1f} kN  {r['ratio']:.3f} x"
              + ("  SPIKE" if r["spike"] else "") + ("  <-- best" if r["best"] else ""))


if __name__ == "__main__":
    main()
