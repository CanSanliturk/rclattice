"""VK3 bridge pier — REVERSED CYCLIC analysis of the Aydin-calibrated lattice.

Drives the nonlinear lattice through the pier's own loading history to the combined shear + axial
failure at 1.59% drift. This is the deliverable of the VK3 study: it is the only one of the three
wall packages whose specimen is shear-relevant, and the only one whose chapter reports a
shear/flexure split to compare against.

TWO PROTOCOLS ARE AVAILABLE, and the default is the measured one:

  `--protocol measured` (default) drives the pier's OWN turning points, recovered by `digitize.py`
  from the vector geometry of Fig. 5.13. This is possible only because those loops are ordered
  paths rather than a raster point cloud, and it fixes two things a nominal protocol gets wrong:
  the elastic block was FORCE controlled, so its displacements were 0.85/2.10/5.08 mm rather than
  the 2.68/5.37/8.05 that the measured k0 predicts (k0 is a yield secant, not an elastic stiffness);
  and the "small intermediate cycles" of Fig. 5.9, which the chapter never tabulates for VK3, are
  recovered directly — including the deliberately ASYMMETRIC second one of each pair.

  `--protocol nominal` is the designed protocol from `specimen.protocol()`, kept as the fallback for
  when the digitized file has not been produced and as the like-for-like control.

WHAT CAN AND CANNOT BE JUDGED. VK3's longitudinal bars are continuous and deformed, so perfect bond
is a fair idealization — unlike SW-NC-FF, whose plain-bar rocking capped that study. Fair here:
peak strength, stiffness, loop shape, energy dissipation, residual displacement, where damage
localizes, and above all the SHEAR / FLEXURE SPLIT. Not fair: anything at or past the failure. The
pier failed by sliding on crossed diagonal cracks once bar buckling and core spalling had destroyed
the base compression zone, and the model has no aggregate interlock, no bar buckling in `Steel02`,
and hoops that can never lose the 90-degree hook anchorage the real ones could (Sec. 5.2.1b). Expect
the model to OVER-retain strength after ~1.27% drift.

HOW THE PROTOCOL IS DRIVEN. A static path-follower stalls early on a cracking lattice, so the
history is imposed as a transient solve. THE CHAPTER STATES NO ACTUATOR RATE ANYWHERE — WSH3 at
least printed 1.2-3.6 mm/min — so unlike the other two studies there is not even a number to be
embarrassed about, and the licence for the drive is entirely the MEASURED residual that
`--quasi-static` reports. Damping, not rate, is the knob that controls it (D64).

Output: examples/output/vk3_wall/vk3_cyclic[_elastic][_dynamic].png. Units: N, mm.
Run from src/: python examples/vk3_wall/cyclic.py --solver dynamic [--drift 0.0095] [--draw]
"""

from __future__ import annotations

import json
import time

from rclattice import viz
from rclattice.opensees import cyclic_protocol, run_cyclic, run_cyclic_dynamic

import gauge
import response as resp
from build import calibrate, nonlinear_pier_lattice, report_calibration, strut_life
from specimen import (
    A_SHEAR, DAMPING_RATIO, HORIZON, LVDT_ROWS, MEASURED_LOOPS, MESH, OUT, QUASI_STATIC_RATE,
    SPECIMEN, base_nodes, control_node, drive_nodes, gauge_nodes, gauge_probe, lateral_loads,
    measured_protocol, protocol, run_dir,
)
import testdata as td

# Hand section analysis of the AS-BUILT model section (bar areas/positions + axial load), which
# reproduces the chapter's own nominal moment capacity to 0.4% (2798 vs 2808 kN.m). A valid
# pushover CANNOT exceed this, so it is the physical bound `response.summarize` checks against.
SECTION_CAPACITY_KN = 848.0



def measured_backbone():
    """`(push, pull, label)` loop-tip backbone of the test, from the digitized Fig. 5.13.

    Falls back to the chapter's printed milestones, which are much coarser: VK3's peak base shear is
    not tabulated anywhere, so without the digitization there is barely a backbone to compare to.
    """
    if MEASURED_LOOPS.exists():
        import numpy as np
        d = np.load(MEASURED_LOOPS)
        lv = float(d["Lv_mm"])
        lvl, up, dn = d["bb_level_mm"], d["bb_push_kN"], d["bb_pull_kN"]
        pos = [(0.0, 0.0)] + [(float(a) / lv * 100.0, float(b))
                              for a, b in zip(lvl, up) if b == b]
        neg = [(0.0, 0.0)] + [(-float(a) / lv * 100.0, float(b))
                              for a, b in zip(lvl, dn) if b == b]
        return pos, neg, "digitized Fig. 5.13, 1st cycles"
    fn = td.PREDICTION["Fn_kN"]
    pos = [(0.0, 0.0), (td.drift(td.PREDICTION["delta_y_mm"]), td.PREDICTION["Fy_prime_kN"]),
           (td.MEASURED["drift_u_pct"], fn)]
    return pos, [(-d, -v) for d, v in pos], "chapter milestones (no digitized loops)"


def backbone(disp, shear):
    """Peak base shear reached at each new displacement extreme, both directions."""
    pos, neg, hi, lo = [], [], 0.0, 0.0
    for u, s in zip(disp, shear):
        if u > hi:
            hi = u
            pos.append((u, s))
        elif u < lo:
            lo = u
            neg.append((u, s))
    return pos, neg


def main(*, compression: str = "crushing", drift: float | None = None, mesh_size: float = MESH,
         horizon: float = HORIZON, steps_per_mm: float = 4.0, gf_factor: float = 1.0,
         solver: str = "dynamic", periods: float | None = None, damping: float = DAMPING_RATIO,
         rate: float = QUASI_STATIC_RATE, quasi_static: bool = True, gauge_every: int = 200,
         draw_model: bool = False, protocol_kind: str = "measured",
         small_cycles: bool = True, capture: bool = False) -> None:
    tag = f"{solver}" + (f"_gf{gf_factor:g}" if gf_factor != 1.0 else "")
    if drift is not None:
        tag += f"_{drift * 100:g}pct"
    out = run_dir("cyclic", tag=tag + ("_elastic" if compression == "elastic" else ""))
    print(f"output directory: {out}")

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    print(f"strut life eps_ult/eps_cr = {strut_life(gf_factor, mesh_size):.1f}"
          + (f"  (Gf scaled x{gf_factor:g})" if gf_factor != 1.0 else "  (plain-concrete MC90)"))

    model = nonlinear_pier_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    if draw_model:
        # Drawn BEFORE the run: on a multi-hour analysis the model figure is the thing to check
        # first, not after.
        import draw as draw_mod
        draw_mod.main(mesh_size=mesh_size, horizon=horizon, nonlinear=True,
                      savepath=out / "model.png")

    history = measured_protocol(drift, include_small=small_cycles) \
        if protocol_kind == "measured" else None
    if history:
        source = f"MEASURED turning points ({'with' if small_cycles else 'without'} small cycles)"
    else:
        if protocol_kind == "measured":
            print("  (no digitized loops found — falling back to the nominal protocol)")
        peaks, cycles = protocol(drift)
        history = cyclic_protocol(peaks, cycles_per_level=cycles)
        source = "NOMINAL protocol (specimen.protocol)"
    path = sum(abs(b - a) for a, b in zip([0.0] + list(history), history))
    print(f"protocol: {source} — {len(history)} reversals to "
          f"{max(abs(u) for u in history) / A_SHEAR:.2%} drift, {path:,.0f} mm of drive path")

    ctrl, base = control_node(model), base_nodes(model)
    gx, g_rows = gauge_nodes(model)
    probe = (gauge_probe(g_rows), (1, 2))          # u_x AND u_y: the diagonal chain needs both
    print(f"gauge chains: {len(gx)} columns x {len(g_rows)} rows to y = {LVDT_ROWS[-1]:.0f} mm, "
          f"every {gauge_every} steps (+ every reversal)")

    t0 = time.time()
    if solver == "dynamic":
        def report(i, n, u, s):
            print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+7.3%}  shear {s / 1e3:+8.1f} kN"
                  f"  [{time.time() - t0:6.0f}s]", flush=True)

        res = run_cyclic_dynamic(model, control_node=ctrl, control_dof=1, history=history,
                                 drive_nodes=drive_nodes(model), base_nodes=base,
                                 periods_to_peak=periods if periods is not None else 48.0,
                                 rate=rate if rate > 0.0 else None,
                                 steps_per_period=30, damping_ratio=damping,
                                 quasi_static=quasi_static, node_history=probe,
                                 node_history_every=gauge_every, progress=report,
                                 capture=capture)
        integ = res.get("integrator", ("?",))
        print(f"  integrator = {integ[0]}{tuple(integ[1:]) if len(integ) > 1 else ''}, "
              f"dt = T1/30")
        print(f"  T1 = {res['T1']:.4f} s, drive rate = {res['rate']:.4g} mm/s, "
              f"damping {damping:.0%} of critical, {res['steps']} steps")
    else:
        res = run_cyclic(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                         control_dof=1, history=history, dU=1.0 / steps_per_mm, base_nodes=base,
                         node_history=probe, node_history_every=gauge_every)
    elapsed = time.time() - t0

    if not res["disp"]:
        raise RuntimeError("cyclic run produced no steps — the gravity stage failed")
    reached = res.get("reached")
    if reached is None:
        print(f"\ncompleted {res['steps']} steps in {elapsed:.1f}s "
              f"(converged={res['converged']})")
    else:
        print(f"\ncompleted {reached}/{len(history)} reversals in {elapsed:.1f}s "
              f"(converged={res['converged']})")

    turned = max(res["disp"]) > 0 and min(res["disp"]) < 0
    if not turned:
        print("\n*** THE RUN NEVER REVERSED — there is no hysteresis in this result. ***")
        print(f"    It stopped at {max(abs(u) for u in res['disp']) / A_SHEAR:.3%} drift, short of "
              f"the first reversal at {history[0] / A_SHEAR:.3%}. Every figure below is a MONOTONIC "
              f"ramp; fix the solve before reading anything else.")

    pos, neg = backbone(res["disp"], res["shear"])
    d_pct_all = [u / A_SHEAR * 100.0 for u in res["disp"]]
    s_kn_all = [s_ / 1e3 for s_ in res["shear"]]
    summary = resp.summarize(d_pct_all, s_kn_all, capacity_kN=SECTION_CAPACITY_KN, cyclic=True)
    summary_neg = resp.summarize([-x for x in d_pct_all], [-x for x in s_kn_all],
                                 capacity_kN=SECTION_CAPACITY_KN, cyclic=True)
    peak_p = summary.reliable.kN * 1e3        # RELIABLE, never raw max(shear) (D70)
    peak_n = -summary_neg.reliable.kN * 1e3
    max_drift = max(abs(u) for u in res["disp"]) / A_SHEAR
    print(f"  drift reached          +-{max_drift:.3%}")
    print("  --- push (South) ---")
    print(summary.report())
    print("  --- pull (North) ---")
    print(summary_neg.report())
    if MEASURED_LOOPS.exists():
        import numpy as np
        d = np.load(MEASURED_LOOPS)
        tp, tn = float(d["peak_push_kN"]), float(d["peak_pull_kN"])
        print(f"  test VK3 (digitized)   +{tp:.1f} / {tn:.1f} kN"
              f"   -> model/test = {peak_p / 1e3 / tp:.3f} / {peak_n / 1e3 / tn:.3f}")

    contam = None
    if res.get("dynamic"):
        dyn = res["dynamic"]
        peak_dyn = max(abs(x) for x in dyn)
        rms = (sum(x * x for x in dyn) / len(dyn)) ** 0.5
        ref = max(abs(peak_p), abs(peak_n))
        contam = peak_dyn / ref if ref else float("inf")
        head = dyn[:max(1, len(dyn) // 5)]
        tail = dyn[-max(1, len(dyn) // 5):]
        ratio = (max(abs(x) for x in tail) / max(abs(x) for x in head)) if head else float("nan")
        verdict = ("steady — the HHT dissipation term, not a diverging solve" if ratio < 1.5
                   else "GROWING with the response — treat the shear as suspect")
        print(f"  inertia + damping      {peak_dyn / 1e3:.1f} kN peak, {rms / 1e3:.1f} kN rms "
              f"= {contam:.1%} of peak shear")
        print(f"                         last fifth / first fifth = {ratio:.2f}x: {verdict}")
        print(f"                         (HHT inflates this ~60x against Newmark without biasing "
              f"the shear, D64 — compare runs of THIS solver, never absolutely)")

    print(f"  residual displacement  {res['disp'][-1]:+.2f} mm at the end of the last unload")

    d_pct = [u / A_SHEAR * 100.0 for u in res["disp"]]
    s_kn = [s / 1e3 for s in res["shear"]]

    datapath = out / "data.json"
    datapath.write_text(json.dumps({
        "drift_pct": d_pct, "shear_kN": s_kn, "compression": compression, "solver": solver,
        "drift_target": drift, "mesh": mesh_size, "horizon": horizon, "gf_factor": gf_factor,
        "protocol": protocol_kind, "small_cycles": small_cycles,
        "rate": res.get("rate", rate), "damping": damping, "T1": res.get("T1"),
        "integrator": list(res.get("integrator", [])),
        "strut_element": "Truss",
        "dynamic_kN_decimated": [x / 1e3 for x in res.get("dynamic", [])[::max(1, gauge_every)]],
        "dynamic_every": gauge_every,
        "dynamic_peak_kN": (max(abs(x) for x in res["dynamic"]) / 1e3)
                           if res.get("dynamic") else None,
        "dynamic_share": contam,
        "gauge": gauge.payload(res, gx),
        "converged": res["converged"], "steps": res.get("steps"), "elapsed_s": elapsed,
    }))
    print(f"saved raw response to {datapath}")

    exp_pos, exp_neg, source_bb = measured_backbone()
    exp_d = [d for d, _v in exp_neg][::-1] + [d for d, _v in exp_pos]
    exp_s = [v for _d, v in exp_neg][::-1] + [v for _d, v in exp_pos]
    exp_style = {"color": "C3", "ls": "--", "lw": 2, "marker": "o", "ms": 5}

    viz.figure_hysteresis(
        [{"disp": d_pct, "shear": s_kn, "label": f"lattice, compression={compression}",
          "style": {"color": "C0", "lw": 0.9}},
         {"disp": exp_d, "shear": exp_s, "label": f"{SPECIMEN} backbone ({source_bb})",
          "style": exp_style}],
        savepath=str(out / "hysteresis.png"), drift_label="drift ratio (%)",
        shear_label="base shear (kN)",
        title=f"{SPECIMEN} — reversed cyclic, Aydin-calibrated lattice "
              f"(compression={compression})")
    print(f"\nsaved hysteresis to {out / 'hysteresis.png'}")

    viz.figure_pushover(
        [{"disp": [u / A_SHEAR * 100.0 for u, _s in neg][::-1]
                  + [u / A_SHEAR * 100.0 for u, _s in pos],
          "shear": [s / 1e3 for _u, s in neg][::-1] + [s / 1e3 for _u, s in pos],
          "label": f"lattice backbone (loop envelope), peak {peak_p / 1e3:.0f} kN",
          "style": {"color": "C0", "lw": 2}},
         {"disp": exp_d, "shear": exp_s, "label": f"{SPECIMEN} measured backbone ({source_bb})",
          "style": exp_style}],
        savepath=str(out / "backbone.png"), xlabel="drift ratio (%)",
        ylabel="base shear (kN)",
        title=f"{SPECIMEN} — lattice cyclic envelope vs the measured backbone")
    print(f"saved backbone to {out / 'backbone.png'}")

    g_drift, g_strain = gauge.reduce_history(res, gx)
    if g_drift:
        levels = sorted({abs(u) / A_SHEAR * 100.0 for u in history})
        idx = gauge.select_levels(g_drift, levels)
        if (len(g_drift) - 1) not in idx:
            idx.append(len(g_drift) - 1)
        sp_ = out / "strain_profile.png"
        rows = gauge.figure(g_drift, g_strain, gx, indices=idx, savepath=sp_,
                            title=f"{SPECIMEN} — vertical strain across the base section")
        print(f"\nbase strain gauge: {len(g_drift)} profiles recorded, {len(idx)} drawn")
        gauge.print_table(rows)
        print(f"saved strain profile to {sp_}")

        c_drift, mids, curv = gauge.reduce_curvature(res, gx)
        cp = out / "curvature.png"
        crows = gauge.figure_curvature(c_drift, mids, curv, indices=idx, savepath=cp,
                                       title=f"{SPECIMEN} — curvature over height")
        print()
        gauge.print_curvature_table(crows)
        print(f"saved curvature profile to {cp}")

        k_drift, tops, comps = gauge.reduce_components(res, gx)
        comp_path = out / "components.png"
        crow = gauge.figure_components(k_drift, tops, comps, savepath=comp_path,
                                       title=f"{SPECIMEN} — deformation components "
                                             f"(cf. chapter Fig. 5.19 right)")
        print(f"\ndeformation components (test: shear "
              f"{td.COMPONENTS['shear_pct_range'][0]:.0f}-"
              f"{td.COMPONENTS['shear_pct_range'][1]:.0f}%, base crack "
              f"{td.COMPONENTS['base_crack_pct_range'][0]:.0f}-"
              f"{td.COMPONENTS['base_crack_pct_range'][1]:.0f}%):")
        gauge.print_components_table(crow, every=max(1, len(crow) // 14))
        print(f"saved components to {comp_path}")

    if capture and res.get("disps_peak") is not None:
        from rclattice.viz import figure_damage, strut_strains
        from specimen import EC, EPSC0, FT
        eps_cr = FT / EC
        life = strut_life(gf_factor, mesh_size)
        panels = [(f"at peak shear ({peak_p / 1e3:.0f} kN)", model, res["disps_peak"]),
                  (f"at the end ({res['disp'][-1] / A_SHEAR:+.3%} drift)", model,
                   res["disps_final"])]
        dmg = out / "damage.png"
        figure_damage(panels, eps_crack=eps_cr, eps_crush=-EPSC0, savepath=str(dmg),
                      suptitle=f"{SPECIMEN} — strut damage (crack at {eps_cr:.2e}, zero tensile "
                               f"stress at {life:.0f}x that)")
        print(f"\nsaved damage map to {dmg}")
        for label, _m, d in panels:
            els, eps = strut_strains(model, d)
            conc = [(e, v) for e, v in zip(els, eps)
                    if e.kind not in ("longitudinal", "stirrup")]
            n = len(conc)
            print(f"  {label}: {sum(1 for _e, v in conc if v > eps_cr):,}/{n:,} cracked, "
                  f"{sum(1 for _e, v in conc if v > life * eps_cr):,} with zero tensile stiffness, "
                  f"{sum(1 for _e, v in conc if v < -EPSC0):,} past epsc0 in compression")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="VK3 reversed-cyclic analysis, Aydin-calibrated lattice")
    p.add_argument("--drift", type=float, default=None,
                   help="highest drift to include (default: the whole history, 1.59%%). "
                        "e.g. 0.0095 stops after the mu_prov = 3 level")
    p.add_argument("--protocol", dest="protocol_kind", choices=("measured", "nominal"),
                   default="measured",
                   help="'measured' (default) drives the pier's own digitized turning points; "
                        "'nominal' uses the designed protocol")
    p.add_argument("--no-small-cycles", dest="small_cycles", action="store_false",
                   help="drop the small intermediate cycles from the measured history")
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON})")
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale the tensile fracture energy (tension stiffening); VK3's strut life "
                        "is 14.2, close to WSH3's 13.4, which needed 2.0 (D67)")
    p.add_argument("--steps-per-mm", type=float, default=4.0)
    p.add_argument("--solver", choices=("static", "dynamic"), default="dynamic")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic solver: drive speed mm/s (default {QUASI_STATIC_RATE:g}). The "
                        "chapter states NO actuator rate, so read the measured residual")
    p.add_argument("--periods", type=float, default=None)
    p.add_argument("--damping", type=float, default=DAMPING_RATIO)
    p.add_argument("--no-quasi-static", dest="quasi_static", action="store_false")
    p.add_argument("--capture", action="store_true",
                   help="also snapshot the displacement field at peak shear and at the end, and "
                        "draw the strut-level damage map")
    p.add_argument("--draw", dest="draw_model", action="store_true",
                   help="also save the analysis-model lattice figure before the run starts")
    p.add_argument("--gauge-every", type=int, default=200)
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         steps_per_mm=a.steps_per_mm, gf_factor=a.gf_factor, solver=a.solver, periods=a.periods,
         damping=a.damping, rate=a.rate, quasi_static=a.quasi_static, gauge_every=a.gauge_every,
         draw_model=a.draw_model, protocol_kind=a.protocol_kind,
         small_cycles=a.small_cycles, capture=a.capture)
