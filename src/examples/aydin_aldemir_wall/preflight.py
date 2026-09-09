"""Measure what a quasi-static run of Aldemir's wall will COST, before committing to it.

Same role as `katrin_wall/preflight.py` and `vk3_wall/preflight.py`. A dynamic-relaxation pushover
on this wall is a multi-hour commitment and its length is set by two numbers that can only be
measured, not guessed:

  * T1, the fundamental period, which fixes the time step (dt = T1 / steps_per_period). This wall is
    STIFF and LIGHT — 739 kN/mm carrying 1.94 t — so T1 is small and the step count is large. That is
    the opposite of the usual intuition that a squat wall is cheap.
  * the wall-clock cost per step, which depends on the element count and how hard the cracked
    tangent is to solve.

It runs a SHORT dynamic segment (default 400 steps) and extrapolates. The extrapolation is a LOWER
BOUND: those first steps are elastic or barely cracked, and a cracking lattice iterates harder.

Run from src/:
    uv run python examples/aydin_aldemir_wall/preflight.py [--drift 0.01] [--steps 400]
                                                           [--gf-factor 1] [--rate 7.6]
"""

from __future__ import annotations

import argparse
import time

from rclattice.opensees import run_modal, run_pushover_dynamic

from build import calibrate, describe, strut_life, wall_lattice
from specimen import (
    A_SHEAR, DAMPING_RATIO, HORIZON, MESH, QUASI_STATIC_RATE, SPECIMEN, base_nodes, control_node,
    drive_nodes,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drift", type=float, default=0.01, help="drift the real run would target")
    ap.add_argument("--steps", type=int, default=400, help="probe steps to time")
    ap.add_argument("--mesh", type=float, default=MESH)
    ap.add_argument("--horizon", type=float, default=HORIZON)
    ap.add_argument("--gf-factor", type=float, default=1.0)
    ap.add_argument("--rate", type=float, default=QUASI_STATIC_RATE)
    ap.add_argument("--damping", type=float, default=DAMPING_RATIO)
    ap.add_argument("--steps-per-period", type=int, default=40)
    a = ap.parse_args()

    cal = calibrate(mesh_size=a.mesh, horizon=a.horizon)
    model, _e = wall_lattice(cal.area, mesh_size=a.mesh, horizon=a.horizon, nonlinear=True,
                             gf_factor=a.gf_factor)
    print(f"{SPECIMEN}  mesh {a.mesh:g}, horizon {a.horizon:g}, Gf x{a.gf_factor:g}")
    print(f"  {describe(model)}")
    print(f"  strut life {strut_life(a.gf_factor, a.mesh):.1f} orthogonal / "
          f"{strut_life(a.gf_factor, a.mesh, diagonal=True):.1f} diagonal")

    modal = run_modal(model, 3)
    t1 = modal["periods"][0]
    total_mass = sum(m[0] for m in model.masses.values())
    print(f"\n  T1 = {t1 * 1e3:.4f} ms   (modes: "
          + ", ".join(f"{p * 1e3:.3f}" for p in modal["periods"]) + " ms)")
    print(f"  total mass {total_mass:.3f} t — a stiff, light wall, so T1 is SHORT and the step")
    print("  count is large; squat does not mean cheap here.")

    dt = t1 / a.steps_per_period
    target = a.drift * A_SHEAR
    sim_time = target / a.rate
    n_steps = sim_time / dt
    print(f"\n  dt = T1/{a.steps_per_period} = {dt * 1e6:.3f} us")
    print(f"  to {a.drift:.2%} drift = {target:.2f} mm at {a.rate:g} mm/s = {sim_time:.3f} s "
          f"of simulated time")
    print(f"  -> {n_steps:,.0f} steps")

    probe_target = min(target, a.rate * dt * a.steps)
    t0 = time.time()
    res = run_pushover_dynamic(model, control_node=control_node(model), control_dof=1,
                               target=probe_target, drive_nodes=drive_nodes(model),
                               base_nodes=base_nodes(model), rate=a.rate,
                               steps_per_period=a.steps_per_period, damping_ratio=a.damping,
                               quasi_static=True)
    elapsed = time.time() - t0
    done = len(res["disp"])
    per_step = elapsed / max(done, 1)
    print(f"\n  timed {done} steps in {elapsed:.1f}s  ->  {per_step * 1e3:.1f} ms/step")
    print(f"  reached {res['disp'][-1] / A_SHEAR:.5%} drift, shear "
          f"{res['shear'][-1] / 1e3:.1f} kN, converged={res['converged']}")

    hours = n_steps * per_step / 3600.0
    print(f"\n  ESTIMATE for {a.drift:.2%} drift: {n_steps:,.0f} steps x {per_step * 1e3:.1f} ms "
          f"= {hours:.1f} h")
    print("  This is a LOWER BOUND: the timed steps are elastic or barely cracked, and a cracking")
    print("  lattice iterates harder. WSH3's equivalent estimate was 11 h and the run took 21.7 h.")
    for d in (0.0025, 0.005, 0.01):
        h = (d * A_SHEAR / a.rate / dt) * per_step / 3600.0
        print(f"    {d:.2%} drift -> {d * A_SHEAR / a.rate / dt:>10,.0f} steps, {h:6.1f} h")


if __name__ == "__main__":
    main()
