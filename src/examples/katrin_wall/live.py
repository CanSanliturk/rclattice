"""Build the hysteresis from a RUNNING (or dead) cyclic run's log, without waiting for it to finish.

`cyclic.py` writes its JSON only after the analysis completes, so a run that is still going — or one
that crashed — leaves no machine-readable result. But its progress lines carry everything a loop
needs: `step / drift / shear / elapsed`, printed every 200 steps. This reconstructs
(drift %, base shear kN) from them, writes it in the shape `compare.py` expects, and hands over to
that module so the figure is the SAME figure, drawn by the same validated code.

WHAT THE SAMPLING COSTS. 200 steps at dt = T1/30 and 7.6 mm/s is ~1.0 mm of drive per sample, so:
  * loop CORNERS are rounded — the reversal is caught within a sample, not at it;
  * loop TIPS read up to ~1 mm inside the true peak, which at the 46.2 mm level is ~2% of the
    amplitude and a few kN of shear;
  * anything faster than ~1 mm of drive (a single crack releasing) is invisible.
So this is for SHAPE and TREND. For the final numbers use the run's own JSON, which records every
step. Every figure it makes is labelled as reconstructed for that reason.

ENERGY: the model's dissipated energy per half-cycle is integrated here, but there is no test
counterpart and there cannot be — the digitized test data is a point CLOUD, and overlapping loops
cannot be re-sequenced from pixels into the ordered path an integral needs (see `digitize.py`).

Run from src/:  python examples/katrin_wall/live.py [--log <path>] [--stem wsh3_cyclic_live]
Output: examples/output/katrin_wall/<stem>_vs_test.png  (+ <stem>_data.json)
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from specimen import A_SHEAR, OUT

PROGRESS = re.compile(r"\s+step\s+(\d+)/(\d+)\s+drift\s+([+-][\d.]+)%\s+shear\s+([+-][\d.]+) kN"
                      r"\s+\[\s*(\d+)s\]")
DEFAULT_LOG = OUT / "wsh3_cyclic_1p02_gf2.log"


def parse(log: Path) -> dict:
    """(step, drift %, shear kN, elapsed s) from every progress line, plus the run's own header."""
    steps, drift, shear, elapsed = [], [], [], []
    total, compression, gf = None, "crushing", 1.0
    for line in log.read_text(errors="ignore").splitlines():
        m = PROGRESS.match(line)
        if m:
            steps.append(int(m[1])); total = int(m[2])
            drift.append(float(m[3])); shear.append(float(m[4])); elapsed.append(int(m[5]))
        elif line.startswith("nonlinear lattice (compression="):
            compression = line.split("compression=")[1].split(")")[0]
        elif line.startswith("tension stiffening: Gf scaled x"):
            gf = float(line.split("x")[-1].split()[0])
    if not steps:
        raise RuntimeError(f"no progress lines in {log} — the run has not reached step 200 yet")
    return {"step": steps, "drift_pct": drift, "shear_kN": shear, "elapsed_s": elapsed,
            "total_steps": total, "compression": compression, "gf_factor": gf}


def reversals(drift):
    """Indices where the direction of travel changes — the loop tips."""
    return [i - 1 for i in range(2, len(drift))
            if (drift[i] - drift[i - 1]) * (drift[i - 1] - drift[i - 2]) < 0]


def _work(drift, shear, a, b) -> float:
    """Trapezoid integral of V du from sample a to b, in kN.mm."""
    return sum(0.5 * (shear[k] + shear[k + 1])
               * (drift[k + 1] - drift[k]) / 100.0 * A_SHEAR for k in range(a, b))


def cycles(drift, shear):
    """Closed loops and their dissipated energy: [(i0, i1, amplitude %, +tip kN, -tip kN, E)].

    A cycle runs tip to tip to tip — from one POSITIVE reversal to the next — so it starts and ends
    at nearly the same displacement and `∮ V du` over it is the area enclosed, i.e. the energy the
    wall dissipated in that cycle. Half-cycles are not used for this: the integral over an open
    path is not an area and is not dissipation, it is just work done getting there.
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
        # Reject a pair of tips that straddles a PROTOCOL LEVEL CHANGE. Those three points still
        # look like a loop, but the path never closes — it ends at a larger amplitude than it
        # started — so its enclosed "area" is not a dissipated energy. Requiring the three tips to
        # agree within 10% keeps only genuine same-amplitude cycles.
        if max(amps) > 1.10 * min(amps):
            continue
        out.append((a, b, sum(amps) / 3.0, shear[a], shear[j],
                    abs(_work(drift, shear, a, b))))
    return out
def main(*, log: Path = DEFAULT_LOG, stem: str = "wsh3_cyclic_live") -> None:
    d = parse(Path(log))
    n, tot = d["step"][-1], d["total_steps"]
    print(f"{Path(log).name}: {len(d['step'])} samples over {n:,}/{tot:,} steps "
          f"({n / tot:.1%}), {d['elapsed_s'][-1] / 3600:.2f} h elapsed")
    print(f"  compression={d['compression']}, Gf x{d['gf_factor']:g}, "
          f"drift {min(d['drift_pct']):+.3f}..{max(d['drift_pct']):+.3f}%")

    dr, sh = d["drift_pct"], d["shear_kN"]
    tips = reversals(dr)
    print(f"\n  loop tips ({len(tips)} reversals):")
    for i in tips:
        print(f"    step {d['step'][i]:8,}  drift {dr[i]:+7.3f}%  shear {sh[i]:+7.1f} kN")

    loops = cycles(dr, sh)
    if loops:
        print(f"\n  {'cycle':>20} {'amplitude':>10} {'+tip':>9} {'-tip':>9} "
              f"{'dissipated':>12} {'equiv. zeta':>12}")
        for a, b, amp, tp, tn, e in loops:
            # equivalent viscous damping of the loop: E / (2*pi*E_strain), the standard measure
            u = amp / 100.0 * A_SHEAR
            estrain = 0.5 * 0.5 * (abs(tp) + abs(tn)) * u
            print(f"  {d['step'][a]:8,}->{d['step'][b]:8,} {amp:9.3f}% {tp:+8.1f} {tn:+8.1f} "
                  f"{e / 1e3:9.2f} kN.m {e / (4.0 * 3.14159265 * estrain):11.1%}")
        print("  (dissipated = the loop area; equivalent zeta = E / (4*pi*E_strain), the usual\n"
              "   measure of how fat a loop is. The TEST's cannot be computed — its digitized data\n"
              "   is a point cloud and overlapping loops cannot be re-sequenced into a path.)")

    # Hand the reconstruction to compare.py in its own format, so the figure is the same figure.
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{stem}_data.json").write_text(json.dumps({
        "drift_pct": d["drift_pct"], "shear_kN": d["shear_kN"],
        "compression": d["compression"], "solver": "dynamic",
        "converged": f"partial — {n:,}/{tot:,} steps, reconstructed from the log at 200-step stride",
    }))
    import compare
    compare.main(stem=stem)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Hysteresis from a running cyclic run's log")
    p.add_argument("--log", default=str(DEFAULT_LOG))
    p.add_argument("--stem", default="wsh3_cyclic_live")
    a = p.parse_args()
    main(log=Path(a.log), stem=a.stem)
