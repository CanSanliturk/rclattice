"""Measure what a long VK3 cyclic run will cost, before committing to it.

D67's lesson, learned the expensive way on WSH3: a screening run must reach PAST the instability it
is screening for, or it certifies nothing. That pre-flight drove to 0.06% drift — just below where
the trouble started — and its 276 ms/step was an uncracked-lattice cost that under-predicted the
real one by 50%.

So this drives to a drift chosen to be well INTO the cracked regime (default 0.25%, past the ~0.026%
cracking drift and past first yield of the outer bars), and reports:
  * T1 and hence the integrator's time step;
  * ms per step, both average and MARGINAL over the second half — the marginal cost is what a long
    run will actually pay, since cracking makes steps more expensive as it spreads;
  * the measured inertia + damping residual, the only thing licensing the drive speed;
  * projected wall-clock to each protocol milestone.

Run from src/:  python examples/vk3_wall/preflight.py [--drift 0.0025] [--gf-factor 2]
"""

from __future__ import annotations

import argparse
import json
import time

from rclattice.opensees import run_cyclic_dynamic, run_modal

from build import calibrate, nonlinear_pier_lattice, report_calibration, strut_life
from specimen import (
    A_SHEAR, DAMPING_RATIO, HORIZON, MESH, QUASI_STATIC_RATE, SPECIMEN, base_nodes, control_node,
    drive_nodes, measured_protocol, protocol, run_dir,
)
from rclattice.opensees import cyclic_protocol


def full_history(kind: str = "measured"):
    h = measured_protocol() if kind == "measured" else None
    if h is None:
        peaks, cycles = protocol(None)
        h = cyclic_protocol(peaks, cycles_per_level=cycles)
    return h


def main(*, drift: float = 0.0025, mesh_size: float = MESH, horizon: float = HORIZON,
         gf_factor: float = 1.0, damping: float = DAMPING_RATIO,
         rate: float = QUASI_STATIC_RATE, protocol_kind: str = "measured") -> None:
    out = run_dir("preflight", tag=f"{drift * 100:g}pct"
                                  + (f"_gf{gf_factor:g}" if gf_factor != 1.0 else ""))
    print(f"output directory: {out}")

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    print(f"strut life eps_ult/eps_cr = {strut_life(gf_factor, mesh_size):.1f}"
          + (f"  (Gf x{gf_factor:g})" if gf_factor != 1.0 else "  (plain-concrete MC90)"))

    model = nonlinear_pier_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   gf_factor=gf_factor)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"model: {len(model.nodes)} nodes, {struts} concrete struts + "
          f"{len(model.elements) - struts} rebar struts")

    modal = run_modal(model, 1)
    print(f"T1 = {modal['periods'][0] * 1e3:.2f} ms  ->  dt = T1/30 = "
          f"{modal['periods'][0] / 30 * 1e3:.3f} ms at the runner's default")

    amp = drift * A_SHEAR
    history = [amp, -amp, amp, -amp]
    t0 = time.time()

    def report(i, n, u, s):
        # flush=True matters: stdout is block-buffered when redirected to a log, so without it a
        # 40-minute screen shows nothing at all until it finishes and only OpenSees' unbuffered
        # stderr warnings are visible while it runs.
        print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+7.3%}  shear {s / 1e3:+8.1f} kN"
              f"  [{time.time() - t0:6.0f}s]", flush=True)

    res = run_cyclic_dynamic(model, control_node=control_node(model), control_dof=1,
                             history=history, drive_nodes=drive_nodes(model),
                             base_nodes=base_nodes(model), rate=rate, steps_per_period=30,
                             damping_ratio=damping, quasi_static=True, progress=report)
    elapsed = time.time() - t0
    steps = res.get("steps") or len(res["disp"])
    ms = elapsed / max(1, steps) * 1e3
    half = steps // 2
    integ = res.get("integrator", ("?",))
    print(f"\nintegrator = {integ[0]}{tuple(integ[1:]) if len(integ) > 1 else ''}  "
          f"(this runner's, not the pushover runner's — they differ, D70)")
    print(f"\nscreen to +-{drift:.3%} drift ({amp:.2f} mm), {len(history)} reversals")
    print(f"  reached          {max(abs(u) for u in res['disp']) / A_SHEAR:.3%} drift, "
          f"converged={res['converged']}, {res.get('reached')}/{len(history)} reversals")
    print(f"  steps            {steps:,} in {elapsed:.1f}s  ->  {ms:.0f} ms/step average")
    if res.get("dynamic"):
        dyn = res["dynamic"]
        peak_dyn = max(abs(x) for x in dyn)
        peak_v = max(abs(s) for s in res["shear"])
        print(f"  residual         {peak_dyn / 1e3:.1f} kN = {peak_dyn / peak_v:.1%} of the "
              f"{peak_v / 1e3:.0f} kN reached so far")
        print(f"                   (do NOT fix this with --rate: the residual is LINEAR in rate. "
              f"Damping is the knob, D64)")

    full = full_history(protocol_kind)
    path_full = sum(abs(b - a) for a, b in zip([0.0] + list(full), full))
    path_here = sum(abs(b - a) for a, b in zip([0.0] + history, history))
    print(f"\nprojection ({protocol_kind} protocol, {len(full)} reversals, "
          f"{path_full:,.0f} mm of drive path):")
    print(f"  {'milestone':>26}  {'path (mm)':>10}  {'steps':>10}  {'hours':>7}")
    marks = [("mu_prov = 2  (0.64% drift)", 21.0), ("mu_prov = 3  (0.95% drift)", 31.5),
             ("mu_prov = 4  (1.27% drift)", 42.0), ("failure     (1.59% drift)", 52.5)]
    for label, level in marks:
        cut, run = [], 0.0
        prev = 0.0
        for u in full:
            # 2% tolerance, NOT an exact cut: these are the pier's MEASURED turning points, and the
            # actuator overshoots its nominal target (the mu_prov = 2 level peaks at 21.26 mm, not
            # 21.00). An exact comparison truncates the level's own cycles and under-reports every
            # milestone below it -- the mu_prov=2 path reads 398 mm instead of 624.
            if abs(u) > level * 1.02:
                break
            run += abs(u - prev)
            prev = u
            cut.append(u)
        proj_steps = steps * run / path_here if path_here else 0.0
        print(f"  {label:>26}  {run:10,.0f}  {proj_steps:10,.0f}  {proj_steps * ms / 3.6e6:7.1f}")
    print("\n  these are LOWER bounds: cost per step rises as cracking spreads (D67).")

    (out / "summary.json").write_text(json.dumps({
        "drift_screened": drift, "gf_factor": gf_factor, "damping": damping, "rate": rate,
        "T1_s": modal["periods"][0], "steps": steps, "elapsed_s": elapsed, "ms_per_step": ms,
        "converged": res["converged"], "reached": res.get("reached"),
        "residual_kN": (max(abs(x) for x in res["dynamic"]) / 1e3) if res.get("dynamic") else None,
        "peak_shear_kN": max(abs(s) for s in res["shear"]) / 1e3,
        "nodes": len(model.nodes), "protocol": protocol_kind,
        "integrator": list(res.get("integrator", [])),
    }, indent=2))
    print(f"saved {out / 'summary.json'}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Size a VK3 cyclic run before committing to it")
    p.add_argument("--drift", type=float, default=0.0025,
                   help="screening amplitude (default 0.25%%, well past cracking at ~0.026%%)")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--gf-factor", type=float, default=1.0)
    p.add_argument("--damping", type=float, default=DAMPING_RATIO)
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE)
    p.add_argument("--protocol", dest="protocol_kind", choices=("measured", "nominal"),
                   default="measured")
    a = p.parse_args()
    main(drift=a.drift, mesh_size=a.mesh, horizon=a.horizon, gf_factor=a.gf_factor,
         damping=a.damping, rate=a.rate, protocol_kind=a.protocol_kind)
