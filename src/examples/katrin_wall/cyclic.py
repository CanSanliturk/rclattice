"""WSH3 shear wall — REVERSED CYCLIC analysis of the Aydin-calibrated lattice.

Drives the Aydin-calibrated nonlinear lattice through the test's own loading protocol (Fig. 6, the
standard Park protocol): two force-controlled cycles to 0.75*F_y, then two displacement-controlled
cycles at each ductility level mu = 2..6, i.e. drift levels 0.25 / 0.68 / 1.02 / 1.35 / 1.69 /
2.03%. South (+) first in every cycle, as the actuator did.

Strut areas come from Aydin's elastic energy balance (D47) unchanged. Tension softening is
regularized by strut length — his crack-band principle. `--compression` selects the real crushing
backbone (default) or his literal assumption that concrete stays elastic in compression.

WHAT THIS MODEL CAN AND CANNOT BE JUDGED ON — and the list is far more favourable than SW-NC-FF's.
WSH3 is reinforced with DEFORMED bars that stay bonded; the paper decomposes its displacement into
flexure and shear with only a small fixed-end (strain-penetration) part, and shear is just ~12% of
the flexural displacement (Fig. 9b). A perfect-bond lattice is therefore a fair idealization of the
mechanism, which it was NOT for the plain-bar SW-NC-FF wall (74% rocking at 2% drift). So:

  * fair: peak strength, stiffness, the curvature and neutral-axis depth at the base, where damage
    localizes, hysteresis loop SHAPE and energy dissipation, and drift capacity up to the point the
    real wall's bars buckled;
  * still unfair: everything after 1.70% drift. Buckling of the boundary bars began there and one
    D12 corner bar ruptured at 1.79%, ending the test. `Steel02` has neither mechanism, so the model
    will keep carrying load where the wall began to lose it. That is the model behaving as built.

A residual bond effect survives even so: the model has no strain penetration into the foundation
(bars are anchored into an ELASTIC block with perfect bond), so the fixed-end rotation the test
measured separately is missing and the model should read STIFFER than the test at any given load.

HOW THE PROTOCOL IS DRIVEN (`--solver dynamic`). A static path-follower stalls at 0.067% drift on
this wall — earlier than SW-NC-FF's 0.3%, because WSH3 cracks at ~0.025% drift — so the history is
imposed as a transient solve. The drive runs at `specimen.QUASI_STATIC_RATE`, and the important
thing to understand is that THIS IS NOT THE TEST'S RATE and cannot be: the actuator moved at 1.2-3.6
mm per MINUTE (Fig. 6), some 130x slower than anything a transient solve can afford. What licences
the drive is the MEASURED residual — `--quasi-static` reports the inertia + damping the drive and
base reactions fail to balance — together with damping set to the physical 0.2 rather than the
near-critical value dynamic relaxation is usually run at (D64).

Output: examples/output/katrin_wall/wsh3_cyclic[_elastic][_dynamic].png. Units: N, mm.
Run as `python examples/katrin_wall/cyclic.py --solver dynamic [--drift 0.0102]`.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import cyclic_protocol, run_cyclic, run_cyclic_dynamic

import gauge
from build import calibrate, nonlinear_wall_lattice, report_calibration
from specimen import (
    A_SHEAR, DAMPING_RATIO, EPS, GAUGE_H, GAUGE_Y0, HORIZON, LVDT_ROWS, L_V, LW, MESH, OUT,
    QUASI_STATIC_RATE, SPECIMEN, base_nodes, control_node, gauge_nodes, gauge_probe, lateral_loads, protocol,
)

# Measured response of WSH3, for annotation only — never a calibration target.
PAPER_PEAK_KN = 454.0
PAPER_PEAK_NEG_KN = -454.0       # Table 5 reports one V_max; Fig. 7c is near-symmetric
PAPER_PEAK_DRIFT = 0.0135
PAPER_YIELD_DRIFT = 0.0025

# The measured BACKBONE is DIGITIZED, not reconstructed — the single biggest data difference from
# the SW-NC-FF study, where only a handful of reported milestones were available and the backbone
# had to be assembled from them. Here `digitize.py` reads the loop tips straight off Fig. 7c at
# known protocol amplitudes, and the result is checkable: its peak lands within one pixel of Table
# 5's 454 kN. `PAPER_BACKBONE` below is only the fallback for when that file has not been made yet.
DIGITIZED = Path(__file__).resolve().parent / "data" / "wsh3_fig7.npz"

# Fallback milestones, each from a printed number: (drift %, load kN, label, derived?)
PAPER_BACKBONE = [
    (0.0, 0.0, None, False),
    (0.025, 115.6, "cracking (M_cr, Table 5)", True),     # V_cr = M_cr / L_v; drift is derived
    (0.25, 371.0, "first yield, outer bar", False),       # Table 4 delta'_y = 11.3 mm
    (1.35, 454.0, "peak (Table 5)", False),
    (2.03, 454.0, "ultimate (Table 4; 20%-drop never met)", False),
]


def measured_backbone():
    """(push, pull) loop-tip backbone of the test, digitized from Fig. 7c if available.

    Returns lists of (drift %, load kN). Falls back to `PAPER_BACKBONE`'s printed milestones, which
    are coarser but need no digitization.
    """
    if DIGITIZED.exists():
        import numpy as np
        d = np.load(DIGITIZED)
        lv = float(d["Lv_mm"])
        pos = [(0.0, 0.0)] + [(u / lv * 100.0, v) for u, v in zip(d["bb_disp_mm"], d["bb_push_kN"])]
        neg = [(0.0, 0.0)] + [(-u / lv * 100.0, v) for u, v in zip(d["bb_disp_mm"], d["bb_pull_kN"])]
        return pos, neg, "digitized Fig. 7c"
    pos = [(d, v) for d, v, _l, _x in PAPER_BACKBONE]
    neg = [(-d, -v) for d, v in pos]
    return pos, neg, "reported milestones"


def backbone(disp, shear):
    """Peak base shear reached at each new displacement extreme, both directions — the loop envelope."""
    pos, neg = [], []
    hi, lo = 0.0, 0.0
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
         steel_rupture=None, concrete_residual: float = 0.2,
         solver: str = "dynamic", periods: float | None = None,
         damping: float = DAMPING_RATIO, rate: float = QUASI_STATIC_RATE,
         quasi_static: bool = True, gauge_every: int = 200,
         draw_model: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    if gf_factor != 1.0:
        print(f"tension stiffening: Gf scaled x{gf_factor:g} (softening gentler than plain concrete)")

    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor,
                                   steel_rupture=steel_rupture,
                                   concrete_residual=concrete_residual)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    if draw_model:
        # Drawn BEFORE the run: on a multi-hour analysis the model figure is the thing you want to
        # check first, not after.
        import draw as draw_mod
        draw_mod.main(mesh_size=mesh_size, horizon=horizon)

    peaks, cycles = protocol(drift)
    history = cyclic_protocol(peaks, cycles_per_level=cycles)
    path = 4.0 * sum(c * p for p, c in zip(peaks, cycles))
    print(f"protocol to {max(peaks) / A_SHEAR:.2%} drift: {len(peaks)} levels "
          f"({', '.join(f'{p / A_SHEAR:.2%}' for p in peaks)}), {sum(cycles)} cycles, "
          f"{len(history)} reversals, {path:,.0f} mm of drive path")

    ctrl, base = control_node(model), base_nodes(model)

    gx, g_rows = gauge_nodes(model)
    probe = (gauge_probe(g_rows), 2)
    print(f"LVDT chain: {len(gx)} columns across x = {gx[0]:.0f}..{gx[-1]:.0f} mm x "
          f"{len(g_rows)} rows to y = {LVDT_ROWS[-1]:.0f} mm; base profile over "
          f"y = {GAUGE_Y0:.0f}-{GAUGE_Y0 + GAUGE_H:.0f} mm, sampled every {gauge_every} steps "
          f"(+ every reversal)")

    t0 = time.time()
    if solver == "dynamic":
        drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))

        def report(i, n, u, s):
            print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+7.3%}  shear {s / 1e3:+8.1f} kN"
                  f"  [{time.time() - t0:6.0f}s]", flush=True)

        res = run_cyclic_dynamic(model, control_node=ctrl, control_dof=1, history=history,
                                 drive_nodes=drive, base_nodes=base,
                                 periods_to_peak=periods if periods is not None else 48.0,
                                 rate=rate if rate > 0.0 else None,
                                 steps_per_period=30, damping_ratio=damping,
                                 quasi_static=quasi_static, node_history=probe,
                                 node_history_every=gauge_every, progress=report)
        print(f"  T1 = {res['T1']:.4f} s, drive rate = {res['rate']:.4g} mm/s "
              f"= {res['rate'] * res['T1']:.4g} mm per fundamental period, "
              f"damping {damping:.0%} of critical, {res['steps']} steps")
        print(f"  the TEST ran at 0.02-0.06 mm/s (Fig. 6) — this drive is "
              f"~{res['rate'] / 0.06:.0f}x faster, which is why the residual below is the licence "
              f"and the rate is not")
    else:
        res = run_cyclic(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                         control_dof=1, history=history, dU=1.0 / steps_per_mm, base_nodes=base,
                         node_history=probe, node_history_every=gauge_every)
    elapsed = time.time() - t0

    if not res["disp"]:
        raise RuntimeError("cyclic run produced no steps — the gravity stage failed")
    reached = res.get("reached")
    if reached is None:
        print(f"\ncompleted {res['steps']} steps in {elapsed:.1f}s (converged={res['converged']})")
    else:
        print(f"\ncompleted {reached}/{len(history)} reversals in {elapsed:.1f}s "
              f"(converged={res['converged']})")
        if reached < len(history):
            print(f"  stopped inside the "
                  f"{abs(history[min(reached, len(history) - 1)]) / A_SHEAR:.2%} drift level — "
                  f"reported below is what was actually traced")

    # A run that never reversed has no hysteresis to report, and every figure below would be a
    # monotonic ramp drawn as if it were a loop. Say so loudly and once, at the top of the results.
    reversals_done = sum(1 for a, b in zip(res["disp"], res["disp"][1:])
                         if (a - res["disp"][0]) * (b - res["disp"][0]) < 0) if res["disp"] else 0
    turned = max(res["disp"]) > 0 and min(res["disp"]) < 0
    if not turned:
        print("\n*** THE RUN NEVER REVERSED — there is no hysteresis in this result. ***")
        print(f"    It stopped at {max(abs(u) for u in res['disp']) / A_SHEAR:.3%} drift, short of "
              f"the first reversal at {history[0] / A_SHEAR:.3%}. Every figure below is a MONOTONIC "
              f"ramp;\n    the 'backbone' is one branch of a pushover and the loop comparison is "
              f"meaningless. Fix the solve before reading anything else.")

    pos, neg = backbone(res["disp"], res["shear"])
    peak_p = max((s for _u, s in pos), default=0.0)
    peak_n = min((s for _u, s in neg), default=0.0)
    max_drift = max(abs(u) for u in res["disp"]) / A_SHEAR
    print(f"  drift reached          +-{max_drift:.3%}")
    print(f"  peak base shear        +{peak_p / 1e3:.1f} / {peak_n / 1e3:.1f} kN")
    print(f"  test WSH3              +{PAPER_PEAK_KN:.1f} / {PAPER_PEAK_NEG_KN:.1f} kN at "
          f"{PAPER_PEAK_DRIFT:.2%}   -> model/test = {peak_p / 1e3 / PAPER_PEAK_KN:.3f}")

    # The inertia + damping the drive and base reactions fail to balance (D62). Reported without a
    # pass/fail threshold on purpose — this runner's HHT(0.7) leaves a large numerical-dissipation
    # term in the residual that does NOT bias the recorded shear. Compare runs of this solver.
    contam = None
    if res.get("dynamic"):
        dyn = res["dynamic"]
        peak_dyn = max(abs(d) for d in dyn)
        rms = (sum(d * d for d in dyn) / len(dyn)) ** 0.5
        ref = max(abs(peak_p), abs(peak_n))
        contam = peak_dyn / ref if ref else float("inf")
        # A near-CONSTANT residual is the HHT(0.7) numerical-dissipation term, not a diverging
        # solve: it rides on the drive speed, so it does not grow with the shear (D64 measures it
        # at ~60x what Newmark reports on the same model). A residual that GROWS with the shear is
        # the one to worry about, so report the trend, not just the peak.
        head = dyn[:max(1, len(dyn) // 5)]
        tail = dyn[-max(1, len(dyn) // 5):]
        drift_ratio = (max(abs(d) for d in tail) / max(abs(d) for d in head)) if head else float("nan")
        verdict = ("steady — the HHT dissipation term, not a diverging solve"
                   if drift_ratio < 1.5 else "GROWING with the response — treat the shear as suspect")
        print(f"  inertia + damping      {peak_dyn / 1e3:.1f} kN peak, {rms / 1e3:.1f} kN rms "
              f"= {contam:.1%} of peak shear")
        print(f"                         last fifth / first fifth = {drift_ratio:.2f}x: {verdict}")
        print(f"                         (HHT inflates this ~60x against Newmark without biasing "
              f"the shear, D64 — compare runs of THIS solver, never absolutely)")

    # Residual displacement — a FAIR comparison here, unlike SW-NC-FF, where the test's
    # self-centering came from plain-bar rocking this model cannot reproduce. WSH3's bars stay
    # bonded, so its residual is a flexural quantity the model is entitled to be judged on.
    print(f"  residual displacement  {res['disp'][-1]:+.2f} mm at the end of the last unload")

    d_pct = [u / A_SHEAR * 100.0 for u in res["disp"]]
    s_kn = [s / 1e3 for s in res["shear"]]

    # Persist the raw response BEFORE plotting: these runs cost hours, and without this any change
    # to the figure means paying for the whole analysis again.
    stem = ("wsh3_cyclic" + ("_elastic" if compression == "elastic" else "")
            + ("_dynamic" if solver == "dynamic" else ""))
    datapath = OUT / f"{stem}_data.json"
    datapath.write_text(json.dumps({
        "drift_pct": d_pct, "shear_kN": s_kn, "compression": compression, "solver": solver,
        "drift_target": drift, "mesh": mesh_size, "horizon": horizon,
        # PROVENANCE: a record that does not carry the switches it was run with cannot
        # be compared to one that was (D70/D80). Gf is not a neutral knob (x2 = +14.2%
        # base shear, D87) and the failure switches decide whether a capacity exists.
        "gf_factor": gf_factor, "steel_rupture": steel_rupture,
        "concrete_residual": concrete_residual,
        "rate": res.get("rate", rate), "damping": damping, "T1": res.get("T1"),
        # The residual is a DIAGNOSTIC, so it is stored decimated at the gauge stride.
        "dynamic_kN_decimated": [d / 1e3 for d in res.get("dynamic", [])[::max(1, gauge_every)]],
        "dynamic_every": gauge_every,
        "dynamic_peak_kN": (max(abs(d) for d in res["dynamic"]) / 1e3) if res.get("dynamic") else None,
        "dynamic_share": contam,
        "gauge": gauge.payload(res, gx),
        "converged": res["converged"], "steps": res.get("steps"), "elapsed_s": elapsed,
    }))
    print(f"saved raw response to {datapath}")

    exp_pos, exp_neg, source = measured_backbone()
    exp_d = [d for d, _v in exp_neg][::-1] + [d for d, _v in exp_pos]
    exp_s = [v for _d, v in exp_neg][::-1] + [v for _d, v in exp_pos]
    exp_style = {"color": "C3", "ls": "--", "lw": 2, "marker": "o", "ms": 5}

    viz.figure_hysteresis(
        [{"disp": d_pct, "shear": s_kn, "label": f"lattice, compression={compression}",
          "style": {"color": "C0", "lw": 0.9}},
         {"disp": exp_d, "shear": exp_s, "label": f"{SPECIMEN} backbone ({source})",
          "style": exp_style}],
        savepath=str(OUT / f"{stem}.png"), drift_label="drift ratio (%)",
        shear_label="base shear (kN)",
        title=f"{SPECIMEN} — reversed cyclic, Aydin-calibrated lattice "
              f"(compression={compression}); test peak {PAPER_PEAK_KN:.0f} kN")
    print(f"\nsaved hysteresis to {OUT / f'{stem}.png'}")

    viz.figure_pushover(
        [{"disp": [u / A_SHEAR * 100.0 for u, _s in neg][::-1]
                  + [u / A_SHEAR * 100.0 for u, _s in pos],
          "shear": [s / 1e3 for _u, s in neg][::-1] + [s / 1e3 for _u, s in pos],
          "label": f"lattice backbone (loop envelope), peak {peak_p / 1e3:.0f} kN",
          "style": {"color": "C0", "lw": 2}},
         {"disp": exp_d, "shear": exp_s,
          "label": f"{SPECIMEN} measured backbone ({source}), peak {PAPER_PEAK_KN:.0f} kN",
          "style": exp_style}],
        savepath=str(OUT / f"{stem}_backbone.png"), xlabel="drift ratio (%)",
        ylabel="base shear (kN)",
        title=f"{SPECIMEN} — lattice cyclic envelope vs the measured backbone ({source})")
    print(f"saved backbone to {OUT / f'{stem}_backbone.png'}")

    # Vertical strain profile across the base section, one curve per protocol drift level.
    g_drift, g_strain = gauge.reduce_history(res, gx)
    if g_drift:
        levels = [p / A_SHEAR * 100.0 for p in peaks]
        idx = gauge.select_levels(g_drift, levels)
        if (len(g_drift) - 1) not in idx:
            idx.append(len(g_drift) - 1)
        sp = OUT / f"{stem}_strain_profile.png"
        rows = gauge.figure(
            g_drift, g_strain, gx, indices=idx, savepath=sp,
            title=f"{SPECIMEN} — vertical strain across the base section, at each protocol level")
        print(f"\nbase strain gauge: {len(g_drift)} profiles recorded, {len(idx)} drawn "
              f"(protocol levels + the latest state, {g_drift[-1]:+.3f}% drift)")
        gauge.print_table(rows)
        print(f"saved strain profile to {sp}")

        # Curvature over height — the model's counterpart to Fig. 11c, with the base-curvature fit
        # of Fig. 15d marked on each profile and tabulated against the digitized measurement. This
        # is the comparison the SW-NC-FF study has no equivalent of: that paper reports a
        # flexure/shear/rocking split, this one reports curvature.
        c_drift, mids, curv = gauge.reduce_curvature(res, gx)
        cp = OUT / f"{stem}_curvature.png"
        crows = gauge.figure_curvature(
            c_drift, mids, curv, indices=idx, savepath=cp,
            title=f"{SPECIMEN} — curvature over height at each protocol level (cf. Fig. 11c)")
        print()
        gauge.print_curvature_table(crows)
        print(f"saved curvature profile to {cp}")


def _rupture_strain(v: str):
    """`--steel-rupture` takes the literal 'measured' or a strain."""
    return v if v == "measured" else float(v)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="WSH3 reversed-cyclic analysis, Aydin-calibrated lattice")
    p.add_argument("--drift", type=float, default=None,
                   help="highest protocol drift level to include (default: the whole protocol, "
                        "2.03%%). e.g. 0.0102 stops after the mu = 3 level")
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing",
                   help="'crushing' (default) or 'elastic' (Aydin's literal assumption)")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON}). 3.01 gives ~28 "
                        "neighbours per node — the redundant bracing that keeps cracked "
                        "tension-side nodes from going singular")
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale the tensile fracture energy (default 1.0 = plain-concrete MC90). "
                        ">1 represents tension stiffening from embedded reinforcement")
    p.add_argument("--steel-rupture", type=_rupture_strain, default=None,
                   help="bar rupture strain: 'measured' uses each bar's own published A_gt (phi12 7.69%%, phi8 7.34%%, phi6 6.45%%), a number applies one strain to all, omit to leave bars unbreakable. The test's corner bar ruptured at 1.79%% drift, which is what this can be validated against (D93)")
    p.add_argument("--concrete-residual", type=float, default=0.2,
                   help="floor the crushing strength at this fraction of fc (default 0.20, D22's plateau). 0.0 lets struts actually crush — available under the dynamic solver, which has no tangent to go singular")
    p.add_argument("--steps-per-mm", type=float, default=4.0,
                   help="static solver only: displacement-control resolution (default 4 steps/mm)")
    p.add_argument("--solver", choices=("static", "dynamic"), default="dynamic",
                   help="'dynamic' (default: dynamic relaxation) or 'static' — the static "
                        "path-follower stalls at 0.067%% drift on this wall, so it is a diagnostic "
                        "only")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic solver: drive speed in mm/s at the actuator (default "
                        f"{QUASI_STATIC_RATE:g}). NOT the test's rate, which was 0.02-0.06 mm/s "
                        "and is unreachable; read the measured residual instead")
    p.add_argument("--periods", type=float, default=None,
                   help="dynamic solver: LEGACY alternative to --rate. Only used with --rate 0")
    p.add_argument("--damping", type=float, default=DAMPING_RATIO,
                   help=f"dynamic solver: damping ratio (default {DAMPING_RATIO:g}). This is the "
                        "knob that controls the contamination, not --rate (D64)")
    p.add_argument("--no-quasi-static", dest="quasi_static", action="store_false",
                   help="skip the inertia+damping measurement")
    p.add_argument("--draw", dest="draw_model", action="store_true",
                   help="also save the analysis-model lattice figure before the run starts")
    p.add_argument("--gauge-every", type=int, default=200,
                   help="base strain gauge: sample every N steps (default 200). Reversals are "
                        "ALWAYS sampled on top of this, so the loop tips are exact")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         steps_per_mm=a.steps_per_mm, gf_factor=a.gf_factor, steel_rupture=a.steel_rupture,
         concrete_residual=a.concrete_residual, solver=a.solver, periods=a.periods,
         damping=a.damping, rate=a.rate, quasi_static=a.quasi_static, gauge_every=a.gauge_every,
         draw_model=a.draw_model)
