"""TEMPORARY sibling of `cyclic.py`: the SAME 4% cyclic run at the FASTER damping.

Why this exists. The run in flight uses `specimen.DAMPING_RATIO = 0.2` (D64). At 4% drift that
leaves the transient tangent under-regularized, ~70% of steps fail their first Newton attempt and
are rescued by sub-stepping, and the run is heading for ~60 h instead of the ~7 h the committed
study took at zeta = 0.8. This script re-runs the identical analysis at zeta = 0.8 — the value
behind `examples/doc/reports/kutay_wall/report.tex` — so a complete, report-comparable result is
available without waiting for the slow one.

It is NOT a copy of `cyclic.py`. It imports and calls `cyclic.main`, so the model, calibration,
materials, protocol, solver, integrator, tolerances and gauge are the same code, not a
same-looking transcription. Exactly two things differ from the running command
(`--solver dynamic --drift 0.04 --draw`):

  1. damping 0.2 -> 0.8;
  2. every output is redirected into `examples/output/wall/temp_kutay_wall/` by rebinding
     `cyclic.OUT` before the call, so NOTHING this writes can collide with the run in flight.
     `--draw` is also off, since the model figure would be byte-identical and shares a filename.

SAFE TO RUN ALONGSIDE THE SLOW ANALYSIS. Separate process, separate OpenSees domain, separate
output directory; a running Python process does not re-read source files, so adding this file
cannot perturb it. The only shared resource is CPU, and OpenSees here is single-threaded on an
8-core machine at load ~2.5, so this takes an idle core.

INHERITED LIMITATION: like `cyclic.py`, results are written only when the run finishes — there is
no checkpointing, so a kill loses everything.

Run from src/:
    uv run python examples/wall/temp_kutay_wall.py
Delete this file once the damping question is settled; it is scaffolding, not part of the study.
"""

from __future__ import annotations

from pathlib import Path

import cyclic
from specimen import HORIZON, MESH, OUT, QUASI_STATIC_RATE

# The damping of the committed study (report.tex, 463,534 steps in 7.1 h), not the D64 default.
FAST_DAMPING = 0.8

# Redirect every output of `cyclic.main` into a subdirectory of its own. `cyclic` resolves `OUT` as
# a module global at call time, so rebinding it here is enough to isolate the whole run.
# The default name carries the damping, so runs at different zeta cannot overwrite each other —
# the mistake this script exists to avoid making twice.
def out_dir(damping: float) -> Path:
    return OUT / f"damping_{int(round(damping * 100)):03d}"


def main(*, drift: float = 0.04, damping: float = FAST_DAMPING, mesh_size: float = MESH,
         horizon: float = HORIZON, rate: float = QUASI_STATIC_RATE,
         gauge_every: int = 200, out: str | None = None) -> None:
    dest = (OUT / out) if out else out_dir(damping)
    dest.mkdir(parents=True, exist_ok=True)
    cyclic.OUT = dest                          # isolate before anything is written
    print(f"TEMP RUN — damping {damping:.0%} of critical, "
          f"rate {rate:g} mm/s, outputs -> {dest}")
    print("the concurrent zeta = 0.2 analysis is untouched: separate process, separate directory\n")
    cyclic.main(compression="crushing", drift=drift, mesh_size=mesh_size, horizon=horizon,
                solver="dynamic", damping=damping, rate=rate, quasi_static=True,
                gauge_every=gauge_every, draw_model=False)
    print(f"\nTEMP RUN complete. Compare against the report's zeta = 0.8 numbers "
          f"(peak +227.2 / -221.5 kN at +0.36 / -0.30% drift).")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Temporary 4% cyclic run at the study's zeta = 0.8 damping, isolated outputs")
    p.add_argument("--drift", type=float, default=0.04, help="highest protocol level (default 0.04)")
    p.add_argument("--damping", type=float, default=FAST_DAMPING,
                   help=f"damping ratio (default {FAST_DAMPING:g}, the committed study's value)")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON, help=f"strut horizon (default {HORIZON})")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"drive speed, mm/s (default {QUASI_STATIC_RATE:g})")
    p.add_argument("--gauge-every", type=int, default=200, help="strain-gauge stride (default 200)")
    p.add_argument("--out", default=None,
                   help="output subdirectory under examples/output/wall/ "
                        "(default: damping_<zeta*100>, e.g. damping_050)")
    a = p.parse_args()
    main(drift=a.drift, damping=a.damping, mesh_size=a.mesh, horizon=a.horizon, rate=a.rate,
         gauge_every=a.gauge_every, out=a.out)
