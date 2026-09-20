"""Measure what a quasi-static EXPLICIT run of RW2 will COST, before committing to it.

Same role as `aydin_aldemir_wall/preflight.py`, sized the D74 way: the explicit step comes from
`critical_time_step` (the stiffest, lightest element), never from T1. It runs a short explicit
segment, times it, and extrapolates to the requested drift — a LOWER BOUND, since the timed steps
are elastic and a cracking lattice is not cheaper.

Run from src/:
    uv run python examples/thomsen_wallace_wall/preflight.py [--drift 0.025] [--grid rebar]
        [--mesh 25] [--steps 3000] [--rate 7.6] [--damping 0.5]
"""

from __future__ import annotations

import argparse
import time

from rclattice.builders import critical_time_step
from rclattice.opensees import run_modal, run_pushover_dynamic

from build import calibrate, describe, strut_life, wall_lattice
from specimen import (
    A_SHEAR, DAMPING_RATIO, EC, GRID, HORIZON, MESH, QUASI_STATIC_RATE, SPECIMEN, STEEL,
    base_nodes, control_node, drive_nodes,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drift", type=float, default=0.025, help="drift the real run would target")
    ap.add_argument("--steps", type=int, default=3000, help="probe steps to time")
    ap.add_argument("--grid", choices=("rebar", "uniform"), default=GRID)
    ap.add_argument("--mesh", type=float, default=MESH)
    ap.add_argument("--horizon", type=float, default=HORIZON)
    ap.add_argument("--gf-factor", type=float, default=1.0)
    ap.add_argument("--rate", type=float, default=QUASI_STATIC_RATE)
    ap.add_argument("--damping", type=float, default=DAMPING_RATIO)
    a = ap.parse_args()

    cal = calibrate(mesh_size=a.mesh, horizon=a.horizon)
    model, _e = wall_lattice(cal.area, grid=a.grid, mesh_size=a.mesh, horizon=a.horizon,
                             nonlinear=True, gf_factor=a.gf_factor)
    print(f"{SPECIMEN}  {a.grid} grid at {a.mesh:g}, horizon {a.horizon:g}, Gf x{a.gf_factor:g}")
    print(f"  {describe(model)}")
    print(f"  strut life {strut_life(a.gf_factor, a.mesh, grid=a.grid):.1f} orthogonal / "
          f"{strut_life(a.gf_factor, a.mesh, diagonal=True, grid=a.grid):.1f} diagonal")

    modulus = lambda e: STEEL.E0 if e.kind in ("longitudinal", "stirrup") else EC   # noqa: E731
    dt_crit, w_max, worst = critical_time_step(model, modulus)
    t1 = run_modal(model, 1)["periods"][0]
    spp = int(t1 / (0.8 * dt_crit)) + 1
    dt = t1 / spp
    mass = sum(m[0] for m in model.masses.values())
    print(f"\n  T1 = {t1 * 1e3:.3f} ms, total mass {mass:.3f} t")
    print(f"  dt_crit = {dt_crit * 1e6:.3f} us (2/w_max, element {worst})  ->  dt = {dt * 1e6:.3f} us "
          f"= T1/{spp:,}")
    target = a.drift * A_SHEAR
    sim_time = target / a.rate
    n_steps = sim_time / dt
    print(f"  to {a.drift:.2%} drift = {target:.1f} mm at {a.rate:g} mm/s = {sim_time:.2f} s "
          f"-> {n_steps:,.0f} steps")

    probe_target = min(target, a.rate * dt * a.steps)
    t0 = time.time()
    res = run_pushover_dynamic(model, control_node=control_node(model), control_dof=1,
                               target=probe_target, drive_nodes=drive_nodes(model),
                               base_nodes=base_nodes(model), rate=a.rate,
                               steps_per_period=spp, damping_ratio=a.damping,
                               quasi_static=True, integrator=("CentralDifference",))
    elapsed = time.time() - t0
    done = len(res["disp"])
    per_step = elapsed / max(done, 1)
    print(f"\n  timed {done:,} steps in {elapsed:.1f}s  ->  {per_step * 1e3:.2f} ms/step "
          f"({1 / per_step:.0f} steps/s)")
    print(f"  reached {res['disp'][-1] / A_SHEAR:.5%} drift, shear {res['shear'][-1] / 1e3:.2f} kN, "
          f"converged={res['converged']}")
    hours = n_steps * per_step / 3600.0
    print(f"\n  ESTIMATE for {a.drift:.2%} drift: {n_steps:,.0f} steps x {per_step * 1e3:.2f} ms "
          f"= {hours:.1f} h  ({hours / target:.4f} h/mm of drive) — a LOWER BOUND")
    for d in (0.005, 0.01, 0.015, 0.02, 0.025):
        print(f"    {d:.1%} drift -> {d * A_SHEAR / a.rate / dt:>12,.0f} steps, "
              f"{(d * A_SHEAR / a.rate / dt) * per_step / 3600:6.1f} h")


if __name__ == "__main__":
    main()
