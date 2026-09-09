"""Aldemir's squat wall — monotonic pushover, STATIC solver (stage 1 of 3).

The agreed sequence is static pushover -> quasi-static pushover -> quasi-static cyclic. This is
stage 1, and its job is as much DIAGNOSTIC as predictive: it establishes how far a path-following
static solve gets on this wall before the cracking lattice defeats it, which is what decides how the
next two stages are run.

EXPECT IT TO STALL EARLY. With no axial load and f_t = 1.85 MPa, this wall's own cracking shear is
about 148 kN, reached at roughly 0.008% drift — an ORDER OF MAGNITUDE earlier than WSH3 (0.025%),
whose static solver then gave up at 0.067% and forced the whole D66 study onto dynamic relaxation.
SW-NC-FF, which cracks later still, managed 0.3%. So a static stall here is the expected outcome and
not a failure of the model; what matters is WHERE it stalls and what the load path looks like when it
does.

THE LOAD-PATH SPLIT IS THE POINT OF THIS SPECIMEN. At an aspect ratio of 0.75 with no axial load the
wall must carry its shear largely on an inclined compression field, and a lattice represents that
explicitly. `--groups` decomposes the base shear into vertical / horizontal / diagonal struts and
rebar at every converged step — the same instrument D55 used to show the compression cube shedding
its inclined path. For a wall "designed to yield in shear" this is the measurement worth having.

Output: examples/output/aydin_aldemir_wall/runs/<stamp>_pushover_static.../  Units: N, mm.
Run from src/:
    uv run python examples/aydin_aldemir_wall/pushover.py [--drift 0.01] [--gf-factor 2]
        [--mesh 50] [--horizon 1.5] [--groups] [--capture] [--steps 400]
"""

from __future__ import annotations

import json
import time

from rclattice import viz
from rclattice.opensees import run_pushover, run_pushover_dynamic

from build import (
    calibrate, cracking_shear, describe, report_calibration, strut_groups, strut_life, wall_lattice,
)
from specimen import (
    A_SHEAR, DAMPING_RATIO, EC, EPSC0, FT, HORIZON, HW, LW, MESH, QUASI_STATIC_RATE, SPECIMEN,
    base_nodes, control_node, drive_nodes, lateral_loads, run_dir,
)

F_MEASURED = 963.592e3      # N, Table 4
K_MEASURED = 1038.44e3      # N/mm, Table 4


def main(*, drift: float = 0.01, mesh_size: float = MESH, horizon: float = HORIZON,
         gf_factor: float = 1.0, steps: int = 400, groups: bool = False,
         capture: bool = False, field: str = "uniaxial", solver: str = "static",
         rate: float = QUASI_STATIC_RATE, damping: float = DAMPING_RATIO,
         steps_per_period: int = 40, progress_every: int = 2000,
         material: str = "concrete02", reinforced: bool = True,
         integrator: tuple = ("Newmark", 0.5, 0.25), dt_override: float | None = None,
         bond: bool = False, bond_area_ratio: float = 0.01,
         rebar_to_top: bool = False) -> None:
    # Damping and rate go in the run TAG: this study's open question is which drive settings make a
    # cracking squat wall traceable, and that is a comparison ACROSS runs (D70). A fixed stem would
    # overwrite the very thing being compared.
    tag = solver + (f"_gf{gf_factor:g}" if gf_factor != 1.0 else "")
    if material != "concrete02":
        tag += f"_{material}"
    if not reinforced:
        tag += "_norebar"
    if rebar_to_top:
        tag += "_rebartop"
    if bond:
        tag += f"_bond{bond_area_ratio:g}"
    if integrator[0] != "Newmark":
        tag += f"_{integrator[0]}"
    if solver == "dynamic":
        tag += f"_z{damping:g}_r{rate:g}"
    out = run_dir("pushover", tag=tag)
    print(f"output directory: {out}")

    cal = calibrate(mesh_size=mesh_size, horizon=horizon, field=field)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)

    v_cr, d_cr = cracking_shear()
    print(f"\nstrut life eps_ult/eps_cr: orthogonal {strut_life(gf_factor, mesh_size):.1f}, "
          f"diagonal {strut_life(gf_factor, mesh_size, diagonal=True):.1f}"
          + (f"  (Gf x{gf_factor:g})" if gf_factor != 1.0 else "  (as given, Gf = 75 N/m)"))
    print(f"first cracking: V_cr = {v_cr / 1e3:.0f} kN at {d_cr:.4%} drift — "
          f"{drift / d_cr:.0f}x below the {drift:.2%} requested")
    print("  for scale: WSH3 cracked at 0.025% and its STATIC solver stalled at 0.067%;")
    print("  SW-NC-FF cracked later and reached 0.3%. A stall here is expected, not a failure.")

    model, _edges = wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon, nonlinear=True,
                                 gf_factor=gf_factor, material=material,
                                 reinforced=reinforced, bond=bond,
                                 bond_area=bond_area_ratio * cal.area if bond else None,
                                 full_height_rebar=rebar_to_top)
    print(f"\n{SPECIMEN}  {LW:.0f} x {HW:.0f}, mesh {mesh_size:g}, horizon {horizon:g}")
    print(f"  {describe(model)}")

    target = drift * A_SHEAR
    eg = None
    if groups:
        shear_g = strut_groups(model, quantity="shear", mesh_size=mesh_size)
        mom_g = strut_groups(model, quantity="moment", mesh_size=mesh_size)
        eg = {f"V:{k}": v for k, v in shear_g.items()}
        eg.update({f"M:{k}": v for k, v in mom_g.items()})
        print("  load-path probe across y = %.0f mm: " % (0.5 * mesh_size)
              + ", ".join(f"{len(v):,} {k}" for k, v in shear_g.items()))

    t0 = time.time()

    def report(i, n, u, s_):
        print(f"    step {i:7d}/{n}  drift {u / A_SHEAR:+8.4%}  shear {s_ / 1e3:+8.1f} kN"
              f"  [{time.time() - t0:7.0f}s]", flush=True)

    if solver == "dynamic":
        print(f"\n  QUASI-STATIC (dynamic relaxation): drive {rate:g} mm/s across the top row, "
              f"damping {damping:.0%}")
        print("  NEITHER CONSTANT IS MEASURED FOR THIS SPECIMEN — the 2019 paper gives loading rates")
        print("  only for its two sensitivity cases and the 2017 test paper is not in the repo. These")
        print("  are the repo's cross-study values (D62/D64), so the licence is the residual below.")
        res = run_pushover_dynamic(model, control_node=control_node(model), control_dof=1,
                                   target=target, drive_nodes=drive_nodes(model),
                                   base_nodes=base_nodes(model), rate=rate,
                                   steps_per_period=steps_per_period, damping_ratio=damping,
                                   quasi_static=True, element_groups=eg, capture=capture,
                                   progress=report, progress_every=progress_every,
                                   integrator=integrator)
    else:
        res = run_pushover(model, lateral_loads=lateral_loads(model),
                           control_node=control_node(model), control_dof=1, dU=target / steps,
                           target=target, base_nodes=base_nodes(model),
                           algorithm=("ModifiedNewton", "-initial"), element_groups=eg,
                           capture=capture)
    elapsed = time.time() - t0

    if not res["disp"]:
        raise RuntimeError("pushover produced no steps at all")
    end_drift = res["disp"][-1] / A_SHEAR
    peak = max(res["shear"])
    at_peak = res["disp"][res["shear"].index(peak)] / A_SHEAR
    k0 = res["shear"][0] / res["disp"][0] if res["disp"][0] else float("nan")

    print(f"\ntraced to {end_drift:.4%} drift of the {drift:.2%} requested   "
          f"(converged={res['converged']}, {len(res['disp'])} steps, {elapsed:.0f}s)")
    print(f"  initial secant         {k0 / 1e3:8.1f} kN/mm   ({k0 / K_MEASURED:.3f} x measured)")
    print(f"  peak base shear        {peak / 1e3:8.1f} kN at {at_peak:.4%} drift")
    print(f"  test peak (Table 4)    {F_MEASURED / 1e3:8.1f} kN   -> model/test = "
          f"{peak / F_MEASURED:.3f}")
    if not res["converged"]:
        print(f"\n  STALLED at {end_drift:.4%} drift. That is {end_drift / d_cr:.0f}x the cracking")
        print("  drift, so the lattice did crack and redistribute before giving up. Whether the")
        print("  peak above is a CAPACITY or just where the solver stopped is decided by stage 2:")
        print("  if the quasi-static run passes this shear while still ascending, it was a stall.")

    if res.get("dynamic"):
        dyn = res["dynamic"]
        n_start = min(len(dyn) - 1, 10 * steps_per_period)
        start = max(abs(d) for d in dyn[:n_start + 1])
        tail = dyn[n_start:] or dyn
        steady = max(abs(d) for d in tail)
        share = steady / abs(peak) if peak else float("inf")
        print(f"\n  inertia + damping      {steady / 1e3:8.1f} kN = {share:.1%} of peak shear"
              + ("" if share < 0.06 else "   [damping is the knob, not rate: D64]"))
        print(f"  start-up transient     {start / 1e3:8.1f} kN once, at near-zero drift")
        print(f"  drive                  {res['rate']:8.4g} mm/s, damping {damping:.0%}, "
              f"T1 = {res['T1'] * 1e3:.3f} ms")

    if eg and res.get("groups"):
        g = res["groups"]
        v_tot = sum(s_[-1] for k, s_ in g.items() if k.startswith("V:"))
        m_tot = sum(s_[-1] for k, s_ in g.items() if k.startswith("M:"))
        print("\n  LOAD PATH at the last converged step, across a cut at y = "
              f"{0.5 * mesh_size:.0f} mm")
        print(f"    {'':<12s}{'shear (kN)':>13s}{'%':>8s}{'moment (kN.m)':>16s}{'%':>8s}")
        for label in ("vertical", "horizontal", "diagonal", "rebar"):
            vs = g.get(f"V:{label}")
            ms = g.get(f"M:{label}")
            if vs is None and ms is None:
                continue
            v = vs[-1] if vs else 0.0
            m = ms[-1] if ms else 0.0
            print(f"    {label:<12s}{v / 1e3:>13.1f}{v / v_tot * 100 if v_tot else 0:>8.1f}"
                  f"{m / 1e6:>16.1f}{m / m_tot * 100 if m_tot else 0:>8.1f}")
        print(f"    {'sum':<12s}{v_tot / 1e3:>13.1f}{'':>8s}{m_tot / 1e6:>16.1f}")
        print(f"    reconciles to base shear {res['shear'][-1] / 1e3:.1f} kN "
              f"({v_tot / res['shear'][-1] if res['shear'][-1] else float('nan'):+.4f}) and to "
              f"V*h = {res['shear'][-1] * A_SHEAR / 1e6:.1f} kN.m "
              f"({m_tot / (res['shear'][-1] * A_SHEAR) if res['shear'][-1] else float('nan'):+.4f})")
        print("    the shear row is a SELF-CHECK, not a mechanism: a vertical strut cut")
        print("    horizontally has no horizontal component, so diagonals take 100% by geometry.")

    drift_pct = [u / A_SHEAR * 100.0 for u in res["disp"]]
    span = [0.0, max(drift_pct)]
    curves = [
        {"disp": drift_pct, "shear": [s / 1e3 for s in res["shear"]],
         "label": f"lattice, static pushover (peak {peak / 1e3:.0f} kN)",
         "style": {"color": "C0", "lw": 2}},
        {"disp": span, "shear": [F_MEASURED / 1e3] * 2,
         "label": f"test peak, Table 4 ({F_MEASURED / 1e3:.0f} kN)",
         "style": {"color": "C3", "ls": "--", "lw": 1.5}},
    ]
    viz.figure_pushover(curves, savepath=str(out / "backbone.png"), xlabel="drift ratio (%)",
                        ylabel="base shear (kN)",
                        title=f"{SPECIMEN} — static monotonic pushover"
                              + (f", Gf x{gf_factor:g}" if gf_factor != 1.0 else ""))

    if capture and res.get("disps_peak") is not None:
        from rclattice.viz import figure_damage
        panels = [(f"at peak shear ({peak / 1e3:.0f} kN)", model, res["disps_peak"]),
                  (f"at the last step ({end_drift:.4%} drift)", model, res["disps_final"])]
        # eps_crush is a POSITIVE magnitude — figure_damage negates it internally. Passing
        # -EPSC0 inverts the test to `eps < +epsc0` and flags ~99% of struts as crushed.
        figure_damage(panels, eps_crack=FT / EC, eps_crush=EPSC0,
                      savepath=str(out / "damage.png"))

    (out / "data.json").write_text(json.dumps({
        "specimen": SPECIMEN, "solver": solver, "panel": [LW, HW], "mesh": mesh_size,
        "horizon": horizon, "gf_factor": gf_factor, "calibration_field": field,
        "material": material, "reinforced": reinforced, "integrator": list(integrator),
        "bond": bond, "bond_area_ratio": bond_area_ratio if bond else None,
        "rebar_to_top": rebar_to_top,
        "area": cal.area, "converged": res["converged"], "steps": len(res["disp"]),
        "end_drift": end_drift, "peak_shear": peak, "drift_at_peak": at_peak,
        "k_initial": k0, "F_measured": F_MEASURED, "K_measured": K_MEASURED,
        "V_cr": v_cr, "drift_cr": d_cr, "elapsed_s": elapsed,
        "disp": res["disp"], "shear": res["shear"],
        "groups": res.get("groups"), "rate": res.get("rate"), "T1": res.get("T1"),
        "damping": damping if solver == "dynamic" else None,
        # The residual SERIES, not just its max: a single number cannot distinguish a steady
        # contamination from a one-off spike at a collapse, and that distinction decides whether the
        # recorded shear is usable at all.
        "dynamic_residual": res.get("dynamic"),
        # Captured displacement fields, so damage figures can be redrawn without re-running.
        "disps_peak": res.get("disps_peak"), "disps_final": res.get("disps_final"),
    }, indent=2))
    print(f"\nsaved to {out}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Aldemir squat wall — static monotonic pushover")
    p.add_argument("--drift", type=float, default=0.01, help="target drift (default 1%%, the test's)")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale Gf for tension stiffening; WSH3 needed 2 (D67)")
    p.add_argument("--steps", type=int, default=400, help="displacement-control steps to the target")
    p.add_argument("--groups", action="store_true", help="decompose base shear by strut orientation")
    p.add_argument("--capture", action="store_true", help="also save the damage/crack figure")
    p.add_argument("--calibration", choices=("uniaxial", "equibiaxial"), default="uniaxial")
    p.add_argument("--no-rebar", dest="reinforced", action="store_false",
                   help="strip ALL reinforcement — the falsification test for the "
                        "bar-chain hypothesis (D72)")
    p.add_argument("--rebar-to-top", action="store_true",
                   help="run the longitudinal bars to the top face instead of stopping one "
                        "spacing short of it (D78): removes the plain-concrete band the drive "
                        "is applied through")
    p.add_argument("--bond", action="store_true",
                   help="Aydin bond elements (steel on its own nodes, tied by a horizon ring). "
                        "The ring stiffness is set by --bond-area-ratio, calibrated so the ELASTIC "
                        "stiffness matches perfect bond (D75).")
    p.add_argument("--bond-area-ratio", type=float, default=0.01,
                   help="bond link area as a fraction of the calibrated strut area A_t. 0.01 gives "
                        "K/K_perfect = 1.006; 1.0 (the D72 default) gives 2.37 and is wrong; below "
                        "~0.003 the steel nodes lose restraint and the tangent goes singular.")
    p.add_argument("--material", choices=("concrete02", "aydin"), default="concrete02",
                   help="aydin = his trilinear tension backbone, tension-only struts "
                        "(2019 Fig. 1c). Monotonic only — the law is path-independent.")
    p.add_argument("--solver", choices=("static", "dynamic"), default="static",
                   help="dynamic = quasi-static relaxation (stage 2); static stalls at ~0.03%% drift")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic: drive speed mm/s (default {QUASI_STATIC_RATE:g}, NOT measured "
                        "for this specimen)")
    p.add_argument("--damping", type=float, default=DAMPING_RATIO,
                   help=f"dynamic: Rayleigh damping ratio (default {DAMPING_RATIO:g}; D64 found "
                        "this is the knob that matters, not rate)")
    p.add_argument("--progress-every", type=int, default=2000,
                   help="dynamic: print progress every N steps (drop it for a diagnostic run)")
    p.add_argument("--explicit", nargs="?", const="CentralDifference", default=None,
                   help="use an EXPLICIT integrator (default CentralDifference). "
                        "Conditionally stable: --steps-per-period must put dt below "
                        "2/w_max, NOT T1/40. Damping becomes mass-proportional only.")
    p.add_argument("--steps-per-period", type=int, default=40,
                   help="dynamic: integration steps per fundamental period, i.e. dt = T1/this")
    a = p.parse_args()
    main(drift=a.drift, mesh_size=a.mesh, horizon=a.horizon, gf_factor=a.gf_factor, steps=a.steps,
         groups=a.groups, capture=a.capture, field=a.calibration, solver=a.solver, rate=a.rate,
         damping=a.damping, steps_per_period=a.steps_per_period,
         progress_every=a.progress_every, material=a.material,
         reinforced=a.reinforced, bond=a.bond, bond_area_ratio=a.bond_area_ratio,
         rebar_to_top=a.rebar_to_top,
         integrator=((a.explicit,) if a.explicit else ("Newmark", 0.5, 0.25)))
