"""Size the WSH3 cyclic run before paying for it: measure T1, cost per step and contamination.

A cyclic run on this model is a multi-hour commitment, and three of the numbers that decide how
long it takes — the fundamental period, the seconds per step, and how much of the recorded base
shear is the solver rather than the wall — can only be MEASURED, not estimated. This runs the
smallest thing that produces all three: a short dynamic pushover at the production drive settings,
far enough to be past first cracking (~0.025% drift) so the timing reflects a CRACKING lattice
rather than an elastic one, which is a different cost per step entirely.

It then extrapolates to each protocol scope, so the choice of how far to run is made on measured
numbers. The extrapolation is deliberately crude — cost per step rises as more struts crack — so
treat the projections as lower bounds and read the multiplier, not the hour count.

Run from src/:  python examples/katrin_wall/preflight.py [--drift 0.0006] [--rate 7.6]
Nothing is saved; this only prints.
"""

from __future__ import annotations

import time

from rclattice.builders import select_nodes
from rclattice.opensees import run_pushover_dynamic

from build import calibrate, nonlinear_wall_lattice, report_calibration
from specimen import (
    A_SHEAR, DAMPING_RATIO, EPS, HORIZON, LVDT_ROWS, LW, MESH, PROTOCOL_CYCLES, PROTOCOL_PEAKS,
    QUASI_STATIC_RATE, base_nodes, control_node, gauge_nodes, gauge_probe,
)

# Drive path of each candidate scope, mm — 4 * amplitude per full cycle.
SCOPES = {
    "0.68% (mu=2)": 2,
    "1.02% (mu=3)": 3,
    "1.35% (mu=4)": 4,
    "2.03% (full)": 6,
}


def main(*, drift: float = 0.0006, mesh_size: float = MESH, horizon: float = HORIZON,
         rate: float = QUASI_STATIC_RATE, damping: float = DAMPING_RATIO,
         compression: str = "crushing", gf_factor: float = 1.0) -> None:
    t_build = time.time()
    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor)
    t_build = time.time() - t_build
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    from specimen import EC, FT, GF_C
    ets = min(FT * FT * mesh_size / (2.0 * GF_C * gf_factor), 0.5 * EC)
    print(f"model: {len(model.nodes)} nodes, {struts} concrete struts + "
          f"{len(model.elements) - struts} rebar struts  [built in {t_build:.1f}s]")
    print(f"tension: f_t = {FT:.2f} MPa, Gf = {GF_C * gf_factor:.4f} N/mm (x{gf_factor:g}) -> "
          f"Ets = {ets:.0f} MPa = {ets / EC:.3f}*E, strut fails at "
          f"{1.0 + EC / ets:.1f}x its cracking strain")
    print(f"         (the SW-NC-FF wall that ran to completion sits at 0.066*E and 16.2x; the first "
          f"WSH3 attempt at 0.110*E and 10.1x lost convergence at 0.065% drift)")

    ctrl, base = control_node(model), base_nodes(model)
    gx, g_rows = gauge_nodes(model)
    probe = (gauge_probe(g_rows), 2)
    drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))
    print(f"LVDT chain: {len(gx)} columns x {len(g_rows)} rows to y = {LVDT_ROWS[-1]:.0f} mm "
          f"({len(probe[0])} probed DOFs)")

    target = drift * A_SHEAR
    print(f"\ndriving {target:.2f} mm ({drift:.3%} drift) at {rate:g} mm/s, damping {damping:.0%} "
          f"— past the ~0.025% cracking drift, so the timing is a CRACKING lattice's")
    t0 = time.time()
    res = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=target,
                               drive_nodes=drive, base_nodes=base, rate=rate,
                               steps_per_period=30, damping_ratio=damping,
                               quasi_static=True, node_history=probe, node_history_every=20)
    elapsed = time.time() - t0

    steps = len(res["disp"]) - 1
    if steps <= 0:
        raise RuntimeError("pre-flight produced no steps — the gravity stage failed")
    t1 = res["T1"]
    dt = t1 / 30.0
    per_step = elapsed / steps
    peak = max(abs(s) for s in res["shear"])

    print(f"\nMEASURED")
    print(f"  T1                     {t1:.5f} s   ->  dt = T1/30 = {dt * 1e3:.3f} ms")
    print(f"  drive                  {res['rate']:g} mm/s = {res['rate'] * t1 * 1e3:.3f} um per "
          f"fundamental period")
    print(f"  steps                  {steps:,} in {elapsed:.1f}s  ->  "
          f"**{per_step * 1e3:.0f} ms/step**")
    print(f"  peak shear so far      {peak / 1e3:.1f} kN at {res['disp'][-1] / A_SHEAR:.3%} drift "
          f"(converged={res['converged']})")
    if res.get("dynamic"):
        dyn = res["dynamic"]
        n_start = min(len(dyn) - 1, 10 * 30)
        start = max(abs(d) for d in dyn[:n_start + 1])
        tail = dyn[n_start:] or dyn
        steady = max(abs(d) for d in tail)
        print(f"  inertia + damping      {steady / 1e3:.1f} kN steady = {steady / peak:.1%} of the "
              f"shear so far; {start / 1e3:.1f} kN start-up transient")
        print(f"                         (read it against the SW-NC-FF wall's 4.0% at its peak — "
              f"early-drift shares run high because the shear is still small)")

    print(f"\nPROJECTED  (path / rate / dt * ms-per-step; a LOWER bound — cost per step rises as "
          f"struts crack)")
    print(f"  {'scope':<16} {'levels':>7} {'path mm':>9} {'sim s':>8} {'steps':>10} {'hours':>8}")
    for label, n in SCOPES.items():
        path = 4.0 * sum(c * p for p, c in zip(PROTOCOL_PEAKS[:n], PROTOCOL_CYCLES[:n]))
        sim = path / rate
        n_steps = sim / dt
        print(f"  {label:<16} {n:>7} {path:>9,.0f} {sim:>8.1f} {n_steps:>10,.0f} "
              f"{n_steps * per_step / 3600.0:>8.1f}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Measure T1, ms/step and contamination before "
                                            "committing to a WSH3 cyclic run")
    p.add_argument("--drift", type=float, default=0.0006,
                   help="how far to drive the sizing pushover (default 0.06%%, just past cracking)")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE)
    p.add_argument("--damping", type=float, default=DAMPING_RATIO)
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing")
    p.add_argument("--gf-factor", type=float, default=1.0)
    a = p.parse_args()
    main(drift=a.drift, mesh_size=a.mesh, horizon=a.horizon, rate=a.rate, damping=a.damping,
         compression=a.compression, gf_factor=a.gf_factor)
