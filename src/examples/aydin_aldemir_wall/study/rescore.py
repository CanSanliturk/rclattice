"""Rescore finished runs with the derived metrics of `metrics.py`, without re-running anything.

`peak_shear` was recorded as a maximum over samples, which reports the ringing of a released bar
rather than resistance on any run carrying failure switches (D92). Every path-following runner
already stores the whole `shear` series, so the correction is a post-processing pass over the
records that exist — 2 h of compute per cell, recovered for a few seconds of arithmetic.

    python rescore.py                 # every run under the study output root, skipping done ones
    python rescore.py --force         # recompute even where the metrics are already present
    python rescore.py --dry-run       # say what would change, write nothing
    python rescore.py <dir> [<dir>…]  # just these

`data.json` is rewritten through a temporary file in the same directory and an atomic rename, so an
interruption cannot cost a run record. Only keys are ADDED; `peak_shear` keeps its recorded value,
so nothing that reads the old key changes meaning.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import metrics                                                                # noqa: E402
import params as P                                                            # noqa: E402
import report as report_mod                                                   # noqa: E402
import specimen                                                               # noqa: E402

ROOT = specimen.OUT / "study"      # the same root run.py writes to; never a second definition


def rewrite_report(d: Path, data: dict) -> None:
    """Regenerate `report.md` so it quotes the corrected numbers, not the ones it was written with.

    Parameters missing from an older `params.json` are filled with their registry defaults, the
    same way `master.py` keeps pre-parameter runs readable.
    """
    pj = d / "params.json"
    if not pj.is_file():
        return
    params = json.loads(pj.read_text())
    for spec in P.REGISTRY:
        params.setdefault(spec.name, spec.default)
    P.normalize(params)          # pre-D94 comp/bond spellings
    report_mod.write(d, params, data)


def rescore(d: Path, *, force: bool = False, dry_run: bool = False) -> str:
    dj = d / "data.json"
    if not dj.is_file():
        return "no data.json"
    data = json.loads(dj.read_text())
    if data.get("peak_shear") is None:
        return "no response (elastic run)"
    if data.get("peak_shear_smooth") is not None and not force:
        return f"already rescored ({data['peak_ringing_ratio']:.3f}x)"
    if not data.get("shear") or not data.get("disp"):
        return "no series stored"

    height = (data.get("meta") or {}).get("drift_denominator_mm")
    if not height:
        return "no drift denominator in meta"
    # `dt` is recorded from D92 on; before that, reconstruct it the way the runner formed it.
    dt = data.get("dt")
    if not dt:
        spp = data.get("steps_per_period_used")
        t1 = data.get("T1")
        if not (spp and t1):
            return "no dt and no T1/steps_per_period to rebuild it"
        dt = t1 / spp
        data["dt"] = dt
        data["dt_note"] = "reconstructed as T1/steps_per_period_used; the runner recorded no dt"

    m = metrics.response_metrics(data["shear"], data["disp"], dt=dt, height=float(height))
    if not m:
        return "smoothing failed"
    data.update(m)
    raw, sm = m["peak_shear_raw"] / 1e3, m["peak_shear_smooth"] / 1e3
    msg = (f"{raw:,.1f} -> {sm:,.1f} kN  ({m['peak_ringing_ratio']:.3f}x, "
           f"{m['drift_at_peak_raw']:.4%} -> {m['drift_at_peak_smooth']:.4%} drift)"
           + ("   <-- RINGING" if m["peak_ringing_ratio"] - 1.0 > metrics.RINGING_FLAG else ""))
    if dry_run:
        return "would write: " + msg
    tmp = dj.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data))
    os.replace(tmp, dj)         # atomic: an interruption cannot cost the record
    rewrite_report(d, data)
    return msg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", type=Path, help="run directories (default: all of them)")
    ap.add_argument("--force", action="store_true", help="recompute where metrics already exist")
    ap.add_argument("--dry-run", action="store_true", help="report, write nothing")
    a = ap.parse_args(argv)

    dirs = a.dirs or sorted(p for p in ROOT.glob("*") if (p / "data.json").is_file())
    if not dirs:
        print(f"no runs under {ROOT}")
        return 1
    for d in dirs:
        print(f"{d.name}\n    {rescore(d, force=a.force, dry_run=a.dry_run)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
