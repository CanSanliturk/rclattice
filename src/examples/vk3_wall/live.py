"""Build the hysteresis from a RUNNING (or dead) VK3 run's log, without waiting for it to finish.

`cyclic.py` writes its JSON and figures only AFTER the analysis completes, so a run still in flight
— or one that crashed at hour 12 of 13 — leaves no machine-readable result and no picture. That is
the reason this exists, and it matters more here than for WSH3: the VK3 cyclic run to failure is a
~13-17 h overnight job driven through 60 measured reversals.

The progress lines carry everything a loop needs — `step / drift / shear / elapsed`, printed every
200 steps — so a coarse hysteresis can be reconstructed from them at any moment.

WHAT THE SAMPLING COSTS. 200 steps at dt = T1/30 = 0.605 ms and 7.6 mm/s is ~0.93 mm of drive per
sample, so:
  * loop CORNERS are rounded — a reversal is caught within a sample, not at it;
  * loop TIPS read up to ~0.9 mm inside the true peak, which at the 31.5 mm level is ~3% of the
    amplitude and a few kN of shear;
  * anything faster than ~0.9 mm of drive (a single crack releasing) is invisible.
So this is for SHAPE and TREND. For final numbers use the run's own `data.json`, which records
every step. Figures made here are labelled reconstructed for that reason.

ENERGY is integrated per closed loop. Unlike the WSH3 study, the VK3 test data CAN be integrated
too — `digitize.py` recovers ordered vector paths rather than a point cloud — so a measured
counterpart is available; see `compare.py`.

Run from src/:  python examples/vk3_wall/live.py [--log <path>] [--kind cyclic]
Output: examples/output/vk3_wall/runs/<run>/live_hysteresis.png  (+ live_data.json)
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

from rclattice import viz

from specimen import A_SHEAR, MEASURED_LOOPS, OUT, RUNS, SPECIMEN, find_run

PROGRESS = re.compile(r"\s+step\s+(\d+)/(\d+)\s+drift\s+([+-][\d.]+)%\s+shear\s+([+-][\d.]+) kN"
                      r"\s+\[\s*(\d+)s\]")


def newest_log(kind: str = "cyclic") -> Path:
    """The most recently touched console log — /tmp while a run is in flight, else the run dir.

    A run's log only lands in its directory when the wrapping shell moves it on exit, so a live run
    must be found in /tmp. Both are checked, newest wins.
    """
    cands = glob.glob("/tmp/vk3_*.log") + glob.glob(str(RUNS / "*" / "console.log"))
    if not cands:
        raise SystemExit("no VK3 console log found")
    return Path(max(cands, key=lambda p: Path(p).stat().st_mtime))


def parse(log: Path) -> dict:
    steps, drift, shear, elapsed = [], [], [], []
    total, compression, gf, run_dir = None, "crushing", 1.0, None
    for line in log.read_text(errors="ignore").splitlines():
        m = PROGRESS.match(line)
        if m:
            steps.append(int(m[1])); total = int(m[2])
            drift.append(float(m[3])); shear.append(float(m[4])); elapsed.append(int(m[5]))
        elif line.startswith("nonlinear lattice (compression="):
            compression = line.split("compression=")[1].split(")")[0]
        elif "strut life" in line and "Gf scaled x" in line:
            gf = float(line.split("Gf scaled x")[1].split(")")[0])
        elif line.startswith("output directory:"):
            run_dir = line.split(":", 1)[1].strip()
    if not steps:
        raise SystemExit(f"no progress lines in {log} — the run has not reached step 200 yet")
    return {"step": steps, "drift_pct": drift, "shear_kN": shear, "elapsed_s": elapsed,
            "total_steps": total, "compression": compression, "gf_factor": gf,
            "run_dir": run_dir}


def reversals(drift):
    """Indices where the direction of travel changes — the loop tips."""
    return [i - 1 for i in range(2, len(drift))
            if (drift[i] - drift[i - 1]) * (drift[i - 1] - drift[i - 2]) < 0]


def _work(drift, shear, a, b) -> float:
    return sum(0.5 * (shear[k] + shear[k + 1])
               * (drift[k + 1] - drift[k]) / 100.0 * A_SHEAR for k in range(a, b))


def cycles(drift, shear):
    """Closed loops and their dissipated energy: `[(i0, i1, amplitude %, +tip, -tip, E)]`.

    A cycle runs tip to tip to tip, so it starts and ends at nearly the same displacement and the
    integral over it is the enclosed area — the energy dissipated. Pairs straddling a protocol level
    change are rejected: the path never closes, so its "area" is not a dissipation.
    """
    tips = reversals(drift)
    pos = [i for i in tips if drift[i] > 0]
    out = []
    for a, b in zip(pos, pos[1:]):
        neg = [i for i in tips if a < i < b and drift[i] < 0]
        if not neg:
            continue
        j = min(neg, key=lambda k: drift[k])
        amps = (abs(drift[a]), abs(drift[j]), abs(drift[b]))
        if max(amps) > 1.10 * min(amps):
            continue
        out.append((a, b, sum(amps) / 3.0, shear[a], shear[j], abs(_work(drift, shear, a, b))))
    return out


def main(*, log: Path | None = None, kind: str = "cyclic") -> None:
    path = Path(log) if log else newest_log(kind)
    d = parse(path)
    n, tot = d["step"][-1], d["total_steps"]
    print(f"{path.name}: {len(d['step'])} samples over {n:,}/{tot:,} steps ({n / tot:.1%}), "
          f"{d['elapsed_s'][-1] / 3600:.2f} h elapsed")
    print(f"  compression={d['compression']}, Gf x{d['gf_factor']:g}, "
          f"drift {min(d['drift_pct']):+.3f}..{max(d['drift_pct']):+.3f}%, "
          f"peak shear {max(d['shear_kN']):+.1f} / {min(d['shear_kN']):+.1f} kN")
    if n < tot:
        rate = d["elapsed_s"][-1] / n
        print(f"  ETA {(tot - n) * rate / 3600:.1f} h at the rate so far "
              f"({rate * 1e3:.0f} ms/step)")

    dr, sh = d["drift_pct"], d["shear_kN"]
    loops = cycles(dr, sh)
    if loops:
        print(f"\n  {'cycle':>20} {'amplitude':>10} {'+tip':>9} {'-tip':>9} "
              f"{'dissipated':>13} {'equiv. zeta':>12}")
        for a, b, amp, tp, tn, e in loops:
            u = amp / 100.0 * A_SHEAR
            estrain = 0.25 * (abs(tp) + abs(tn)) * u
            print(f"  {d['step'][a]:8,}->{d['step'][b]:8,} {amp:9.3f}% {tp:+8.1f} {tn:+8.1f} "
                  f"{e / 1e3:10.2f} kN.m {e / (4.0 * 3.141592653589793 * estrain):11.1%}")
        print("  (dissipated = loop area; equivalent zeta = E/(4*pi*E_strain), the usual measure of\n"
              "   how fat a loop is. Unlike WSH3, the TEST's counterpart IS computable here — the\n"
              "   digitized VK3 loops are ordered paths, not a point cloud.)")

    # Write beside the run, not into the flat output dir. A live log sits in /tmp until the
    # wrapping shell moves it on exit, but every script prints its run directory on the first line,
    # so the destination is recoverable even mid-run.
    out = Path(d["run_dir"]) if d.get("run_dir") else path.parent
    if not out.is_dir():
        out = OUT
    out.mkdir(parents=True, exist_ok=True)
    (out / "live_data.json").write_text(json.dumps({
        "drift_pct": dr, "shear_kN": sh, "compression": d["compression"], "solver": "dynamic",
        "converged": f"PARTIAL — {n:,}/{tot:,} steps, reconstructed from the log at 200-step stride",
    }))

    series = [{"disp": dr, "shear": sh,
               "label": f"lattice, reconstructed ({n / tot:.0%} of the run)",
               "style": {"color": "C0", "lw": 0.9}}]
    if MEASURED_LOOPS.exists():
        import numpy as np
        z = np.load(MEASURED_LOOPS)
        series.append({"disp": list(z["drift_pct"]), "shear": list(z["load_kN"]),
                       "label": f"{SPECIMEN} measured (digitized Fig. 5.13)",
                       "style": {"color": "C3", "lw": 0.6, "alpha": 0.75}})
    fig = out / "live_hysteresis.png"
    viz.figure_hysteresis(series, savepath=str(fig), drift_label="drift ratio (%)",
                          shear_label="base shear (kN)",
                          title=f"{SPECIMEN} — LIVE reconstruction at {n / tot:.0%} "
                                f"({d['elapsed_s'][-1] / 3600:.1f} h), 200-step stride")
    print(f"\nsaved {fig}")
    print(f"saved {out / 'live_data.json'}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Hysteresis from a running VK3 run's log")
    p.add_argument("--log", default=None, help="console log (default: the newest VK3 one)")
    p.add_argument("--kind", default="cyclic")
    a = p.parse_args()
    main(log=a.log, kind=a.kind)
