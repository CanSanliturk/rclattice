"""Measure what a replica pushover costs before committing to it.

    uv run python examples/aydin_aldemir_wall/replica/preflight.py [--bond] [--explicit]

The replica is 20,385 nodes / 88,324 elements — 1.85x the largest model this study has run, and
`--bond` roughly doubles it again. Its cost is set by the time step (which under an explicit
integrator comes from the stiffest, lightest ELEMENT, not from T1) and the per-step solve, and
neither can be guessed from element count alone: the parent study's mesh-25 run was projected at
3.6x mesh-50's cost and came in at 13x.

PASS THE SAME FLAGS YOU INTEND TO RUN WITH. An implicit probe does not price an explicit run — the
stable step is ~20x smaller, so the step COUNT is ~20x larger while each step is ~45x cheaper.
"""
from __future__ import annotations
import argparse, sys, time
sys.path.insert(0, ".")
from rclattice.builders import select_nodes
from rclattice.opensees import run_modal, run_pushover_dynamic
from build import calibrate, describe, explicit_steps_per_period, report_bond, wall_lattice
from specimen import A_SHEAR, EPS, HORIZON, HW, LW, MESH

RATE, SPP, DAMP, PROBE = 7.6, 40, 0.5, 300

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--bond", action="store_true", help="price HIS bond scheme, not shared nodes")
ap.add_argument("--bond-area-ratio", type=float, default=0.01)
ap.add_argument("--explicit", nargs="?", const="CentralDifference", default=None,
                help="price an EXPLICIT run (dt from critical_time_step, D74)")
ap.add_argument("--horizon", type=float, default=HORIZON, help="1.5 or his 3.01")
ap.add_argument("--steps", type=int, default=PROBE, help="timed steps")
a = ap.parse_args()

t0 = time.time()
cal = calibrate(horizon=a.horizon)
bond_area = a.bond_area_ratio * cal.area if a.bond else None
if a.bond:
    report_bond(cal.area, bond_area, mesh_size=MESH)
model, _ = wall_lattice(cal.area, horizon=a.horizon, bond=a.bond, bond_area=bond_area)
print(f"\n{describe(model)}   [built {time.time() - t0:.0f}s]", flush=True)

if a.explicit:
    spp, dt_crit, t1 = explicit_steps_per_period(model)
    dt = t1 / spp
    print(f"T1 = {t1 * 1e3:.4f} ms   dt_crit = {dt_crit * 1e6:.3f} us   "
          f"dt = {dt * 1e6:.3f} us (spp {spp:,})   "
          f"mass {sum(m[0] for m in model.masses.values()):.3f} t", flush=True)
else:
    spp = SPP
    t1 = run_modal(model, 1)["periods"][0]
    dt = t1 / spp
    print(f"T1 = {t1 * 1e3:.4f} ms   dt = {dt * 1e6:.2f} us   "
          f"mass {sum(m[0] for m in model.masses.values()):.3f} t", flush=True)

top = select_nodes(model, (-1., LW + 1., HW - EPS, HW + EPS))
ctrl = select_nodes(model, (-EPS, EPS, HW - EPS, HW + EPS))[0]
base = select_nodes(model, (-1., LW + 1., -EPS, EPS))
t0 = time.time()
r = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=RATE * dt * a.steps,
                         drive_nodes=top, base_nodes=base, rate=RATE, steps_per_period=spp,
                         damping_ratio=DAMP, quasi_static=True,
                         integrator=(a.explicit,) if a.explicit else ("Newmark", 0.5, 0.25))
el = time.time() - t0; ps = el / max(len(r["disp"]), 1)
print(f"\ntimed {len(r['disp'])} steps in {el:.1f}s -> {ps * 1e3:.0f} ms/step  "
      f"(reached {r['disp'][-1] / A_SHEAR:.5%}, shear {r['shear'][-1] / 1e3:.1f} kN, "
      f"converged={r['converged']})", flush=True)
print(f"\nESTIMATES ({'explicit ' + a.explicit if a.explicit else 'implicit Newmark'}, "
      f"{'BOND' if a.bond else 'shared nodes'}) — LOWER BOUNDS, these steps are barely cracked:")
for d in (0.0012, 0.0015, 0.0025):
    n = d * A_SHEAR / RATE / dt
    print(f"  to {d:.2%} drift: {n:,.0f} steps x {ps * 1e3:.0f} ms = {n * ps / 3600:.1f} h")
print("\n  collapse points seen so far: 0.097% / 0.110% / 0.122% / 0.199% drift (no bond).")
print("  In the PARENT study bond moved collapse 0.2030% -> 0.2374%, a 1.21x gain, so budget")
print("  past the no-bond number rather than to it.")
