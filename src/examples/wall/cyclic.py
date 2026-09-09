"""SW-NC-FF shear wall — REVERSED CYCLIC analysis of the Aydin-calibrated lattice.

Drives the Aydin-calibrated nonlinear lattice through the test's own loading protocol (paper
Fig. 9, per ACI 374.2R-13): displacement-controlled cycles at drift ratios 0.3 ... 4.0%, three
cycles up to the yield displacement and two at each level after it.

Strut areas come from Aydin's elastic energy balance (D47) unchanged. Tension softening is
regularized by strut length — his crack-band principle. `--compression` selects the real crushing
backbone (default) or his literal assumption that concrete stays elastic in compression.

WHAT THIS MODEL CAN AND CANNOT BE JUDGED ON. The test's displacement at 2% drift is 74% ROCKING,
driven by progressive debonding of the plain longitudinal bars; only 20% is flexure and 4% shear
(paper Fig. 21b). This lattice ties rebar to concrete at shared nodes — PERFECT bond — so it has
no access to that mechanism. Therefore:

  * fair comparisons: PEAK STRENGTH (section equilibrium: bar areas, fy, fc, axial load, lever
    arm — perfect bond barely moves it), initial stiffness, neutral-axis depth, and where damage
    localizes;
  * unfair comparisons: drift capacity, loop pinching, self-centering, energy dissipation and
    residual displacement — all four are consequences of the debonding this model omits. They are
    reported for honesty, NOT as targets to tune toward.

A predictable, diagnostic direction of error follows: with perfect bond the bars strain more per
unit drift, so this model should yield EARLIER than the measured 0.30% drift and run STIFFER after
cracking. Seeing that is the model behaving correctly given its assumptions, not a defect.

HOW THE PROTOCOL IS DRIVEN (`--solver dynamic`). A static path-follower stalls near 0.3% drift on
a cracking lattice, so the history is imposed as a transient solve. Two settings then decide whether
the recorded base shear is the wall's or the solver's, and both are physical rather than tuned
(`specimen.QUASI_STATIC_RATE` / `DAMPING_RATIO`, D62): the drive runs at 2.0 mm/s at the actuator,
and damping is 5% of critical — the value for cracked RC, not the near-critical setting dynamic
relaxation is usually run at. The leftover is MEASURED, not assumed: `--quasi-static` (on by
default) reports the inertia + damping the drive and base reactions fail to balance as a percentage
of the peak shear, which is the number that licences the two settings above.

Output: examples/output/wall/wall_cyclic[_elastic].png. Units: N, mm.
Run as `python examples/wall/cyclic.py [--drift 0.01] [--compression crushing|elastic]`.
"""

from __future__ import annotations

import json
import time

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import cyclic_protocol, run_cyclic, run_cyclic_dynamic

import gauge
from build import calibrate, nonlinear_wall_lattice, report_calibration
from specimen import (
    A_SHEAR, DAMPING_RATIO, EPS, GAUGE_H, HORIZON, LW, MESH, OUT, QUASI_STATIC_RATE, base_nodes,
    control_node, gauge_nodes, gauge_strains, lateral_loads, neutral_axis, protocol,
)

# Measured response of SW-NC-FF, for annotation only — never a calibration target.
PAPER_PEAK_KN = 212.6
PAPER_PEAK_NEG_KN = -209.4
PAPER_PEAK_DRIFT = 0.0100
PAPER_YIELD_DRIFT = 0.0030

# Measured BACKBONE milestones, reconstructed from the paper's REPORTED VALUES — NOT digitized
# from the Fig. 14b loops, which would need the raw test data. Each point cites its source, and
# `derived` flags the one coordinate that is computed rather than read off:
#   (drift %, load kN, label, derived?)
PAPER_BACKBONE = [
    (0.0, 0.0, None, False),
    # 45.4 kN cracking load is reported; its DRIFT is derived from the reported uncracked rigidity
    # of 0.67 Ec*Ig (Table 4) applied to the gross-section cantilever stiffness.
    (0.047, 45.4, "first cracking", True),
    (1.00, 212.6, "peak", False),                   # Table 3
    # delta_u = 1.5% is reported (Table 3); the LOAD is 0.8*Vu, the conventional 20%-loss
    # definition of ultimate drift, so the ordinate here is definitional rather than measured.
    (1.50, 170.1, "ultimate (0.8 Vu)", True),
    (4.00, 138.2, "4% drift", False),               # 35% degradation reported, Sec. 5
]


def paper_backbone_curves():
    """The measured backbone as two curves (push and pull) for overlay on the model envelope.

    The pull branch mirrors the push branch scaled by the reported negative/positive peak ratio —
    the paper reports both peaks but only one degradation figure, so the shape is shared.
    """
    ratio = abs(PAPER_PEAK_NEG_KN) / PAPER_PEAK_KN
    pos = [(d, v) for d, v, _l, _x in PAPER_BACKBONE]
    neg = [(-d, -v * ratio) for d, v in pos]
    return pos, neg


def backbone(disp, shear):
    """Peak base shear reached at each new displacement extreme, both directions.

    Walks the history keeping the running max |displacement| in each direction and the shear at
    that point, which is what the paper's Fig. 14c backbone is: the envelope of the loops.
    """
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


def main(*, compression: str = "crushing", drift: float = 0.01, mesh_size: float = MESH,
         horizon: float = HORIZON, steps_per_mm: float = 4.0, gf_factor: float = 1.0,
         solver: str = "static", periods: float | None = None,
         damping: float = DAMPING_RATIO, rate: float = QUASI_STATIC_RATE,
         quasi_static: bool = True, gauge_every: int = 200,
         draw_model: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    if gf_factor != 1.0:
        print(f"tension stiffening: Gf scaled x{gf_factor:g} (softening gentler than plain concrete)")

    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    if draw_model:
        # Drawn BEFORE the run: on a multi-hour analysis the model figure is the thing you want to
        # check first, not after. `draw.main` rebuilds the same model from the same calibration.
        import draw as draw_mod
        draw_mod.main(mesh_size=mesh_size, horizon=horizon)

    peaks, cycles = protocol(drift)
    history = cyclic_protocol(peaks, cycles_per_level=cycles)
    print(f"protocol to {drift:.2%} drift: {len(peaks)} levels "
          f"({', '.join(f'{p / A_SHEAR:.2%}' for p in peaks)}), {sum(cycles)} cycles, "
          f"{len(history)} reversals")

    ctrl, base = control_node(model), base_nodes(model)

    # Vertical strain gauge over the wall base (D63): two aligned node rows GAUGE_H apart, asked
    # for as ONE concatenated list (bottom then top) so a single probe covers both. dof 2 = uy.
    gx, g_bot, g_top = gauge_nodes(model)
    probe = (list(g_bot) + list(g_top), 2)
    print(f"base strain gauge: {len(gx)} columns across x = {gx[0]:.0f}..{gx[-1]:.0f} mm, "
          f"{GAUGE_H:.0f} mm gauge, sampled every {gauge_every} steps (+ every reversal)")

    t0 = time.time()
    if solver == "dynamic":
        # Dynamic relaxation through the whole reversing history (D49). A static path-follower
        # cannot cross the softening instability of a cracking lattice; mass and damping can.
        drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))

        def report(i, n, u, s):
            print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+7.3%}  shear {s / 1e3:+8.1f} kN"
                  f"  [{time.time() - t0:6.0f}s]", flush=True)

        res = run_cyclic_dynamic(model, control_node=ctrl, control_dof=1, history=history,
                                 drive_nodes=drive, base_nodes=base,
                                 periods_to_peak=periods if periods is not None else 48.0,
                                 # --rate 0 hands control back to the legacy --periods route
                                 rate=rate if rate > 0.0 else None,
                                 steps_per_period=30, damping_ratio=damping,
                                 quasi_static=quasi_static, node_history=probe,
                                 node_history_every=gauge_every, progress=report)
        print(f"  T1 = {res['T1']:.4f} s, drive rate = {res['rate']:.4g} mm/s "
              f"= {res['rate'] * res['T1']:.4g} mm per fundamental period, "
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
        if reached < len(history):
            print(f"  stopped inside the "
                  f"{abs(history[min(reached, len(history) - 1)]) / A_SHEAR:.2%} drift level — "
                  f"reported below is what was actually traced")

    pos, neg = backbone(res["disp"], res["shear"])
    peak_p = max((s for _u, s in pos), default=0.0)
    peak_n = min((s for _u, s in neg), default=0.0)
    max_drift = max(abs(u) for u in res["disp"]) / A_SHEAR
    print(f"  drift reached          +-{max_drift:.3%}")
    print(f"  peak base shear        +{peak_p / 1e3:.1f} / {peak_n / 1e3:.1f} kN")
    print(f"  paper SW-NC-FF         +{PAPER_PEAK_KN:.1f} / -209.4 kN at {PAPER_PEAK_DRIFT:.2%}"
          f"   -> model/test = {peak_p / 1e3 / PAPER_PEAK_KN:.3f}")

    # How much of that shear is the wall and how much is the solver: the inertia + damping the drive
    # and base reactions fail to balance (D62). Reported without a pass/fail threshold on purpose —
    # this runner's HHT(0.7) leaves a large numerical-dissipation term in the residual that does NOT
    # bias the recorded shear (on an elastic wall it reads ~2x the base shear while reproducing the
    # static answer to 0.24%). The number is for comparing runs of THIS solver, where it is linear
    # in --rate and falls with --damping once cracking starts.
    contam = None
    if res.get("dynamic"):
        dyn = res["dynamic"]
        peak_dyn = max(abs(d) for d in dyn)
        rms = (sum(d * d for d in dyn) / len(dyn)) ** 0.5
        ref = max(abs(peak_p), abs(peak_n))
        contam = peak_dyn / ref if ref else float("inf")
        print(f"  inertia + damping      {peak_dyn / 1e3:.1f} kN peak, {rms / 1e3:.1f} kN rms "
              f"= {contam:.1%} of peak shear  (compare across runs of this solver, not absolutely "
              f"— HHT inflates it)")

    # Residual displacement at the end of the last completed unload — a consequence of the
    # debonding mechanism this model omits, so reported as a known-unfair comparison.
    print(f"  residual displacement  {res['disp'][-1]:+.2f} mm  (test: negligible — its "
          f"self-centering comes from plain-bar rocking, which perfect bond cannot reproduce)")

    d_pct = [u / A_SHEAR * 100.0 for u in res["disp"]]
    s_kn = [s / 1e3 for s in res["shear"]]

    # Persist the raw response BEFORE plotting: these runs cost tens of minutes to hours, and
    # without this any change to the figure means paying for the whole analysis again.
    datapath = OUT / (("wall_cyclic" + ("_elastic" if compression == "elastic" else "")
                       + ("_dynamic" if solver == "dynamic" else "")) + "_data.json")
    datapath.write_text(json.dumps({
        "drift_pct": d_pct, "shear_kN": s_kn, "compression": compression, "solver": solver,
        "drift_target": drift, "mesh": mesh_size, "horizon": horizon,
        "rate": res.get("rate", rate), "damping": damping, "T1": res.get("T1"),
        # The residual is a DIAGNOSTIC, so it is stored decimated at the gauge stride: a 4%
        # protocol is ~1.8M steps, and three full-length arrays would make this file ~100 MB.
        # Its peak/rms summary is kept exactly, which is what the reporting above actually uses.
        "dynamic_kN_decimated": [d / 1e3 for d in res.get("dynamic", [])[::max(1, gauge_every)]],
        "dynamic_every": gauge_every,
        "dynamic_peak_kN": (max(abs(d) for d in res["dynamic"]) / 1e3) if res.get("dynamic") else None,
        "dynamic_share": contam,
        "gauge": gauge.payload(res, gx),
        "converged": res["converged"], "steps": res.get("steps"), "elapsed_s": elapsed,
    }))
    print(f"saved raw response to {datapath}")
    exp_pos, exp_neg = paper_backbone_curves()
    series = [
        {"disp": d_pct, "shear": s_kn, "label": f"lattice, compression={compression}",
         "style": {"color": "C0", "lw": 0.9}},
        {"disp": [d for d, _v in exp_neg][::-1] + [d for d, _v in exp_pos],
         "shear": [v for _d, v in exp_neg][::-1] + [v for _d, v in exp_pos],
         "label": "SW-NC-FF measured backbone",
         "style": {"color": "C3", "ls": "--", "lw": 2, "marker": "o", "ms": 4}},
    ]
    stem = ("wall_cyclic" + ("_elastic" if compression == "elastic" else "")
            + ("_dynamic" if solver == "dynamic" else ""))
    savepath = OUT / f"{stem}.png"
    viz.figure_hysteresis(series, savepath=str(savepath), drift_label="drift ratio (%)",
                          shear_label="base shear (kN)",
                          title=f"SW-NC-FF — reversed cyclic, Aydin-calibrated lattice "
                                f"(compression={compression}); test peak {PAPER_PEAK_KN:.0f} kN")
    print(f"\nsaved hysteresis to {savepath}")

    bb = OUT / f"{stem}_backbone.png"
    exp_pos, exp_neg = paper_backbone_curves()
    curves = [
        {"disp": [u / A_SHEAR * 100.0 for u, _s in neg][::-1] + [u / A_SHEAR * 100.0 for u, _s in pos],
         "shear": [s / 1e3 for _u, s in neg][::-1] + [s / 1e3 for _u, s in pos],
         "label": f"lattice backbone (loop envelope), peak {peak_p / 1e3:.0f} kN",
         "style": {"color": "C0", "lw": 2}},
        {"disp": [d for d, _v in exp_neg][::-1] + [d for d, _v in exp_pos],
         "shear": [v for _d, v in exp_neg][::-1] + [v for _d, v in exp_pos],
         "label": f"SW-NC-FF measured backbone (reported milestones), peak {PAPER_PEAK_KN:.0f} kN",
         "style": {"color": "C3", "ls": "--", "lw": 2, "marker": "o", "ms": 5}},
    ]
    viz.figure_pushover(
        curves, savepath=str(bb), xlabel="drift ratio (%)", ylabel="base shear (kN)",
        title="SW-NC-FF — lattice cyclic envelope vs the measured backbone\n"
              "(test points reconstructed from reported values, not digitized loops)")
    print(f"saved backbone to {bb}")

    # Vertical strain profile across the base section, one curve per protocol drift level (D63).
    g_drift, g_strain = gauge.reduce_history(res, gx)
    if g_drift:
        levels = [p / A_SHEAR * 100.0 for p in peaks]
        idx = gauge.select_levels(g_drift, levels)
        # The LATEST state, always drawn whatever the protocol levels resolved to. On a run that
        # stopped early this is the only profile that reflects where the analysis actually got to.
        if (len(g_drift) - 1) not in idx:
            idx.append(len(g_drift) - 1)
        sp = OUT / f"{stem}_strain_profile.png"
        rows = gauge.figure(
            g_drift, g_strain, gx, indices=idx, savepath=sp,
            title="SW-NC-FF — vertical strain across the base section, at each protocol level")
        print(f"\nbase strain gauge: {len(g_drift)} profiles recorded, {len(idx)} drawn "
              f"(protocol levels + the latest state, {g_drift[-1]:+.3f}% drift)")
        gauge.print_table(rows)
        print(f"saved strain profile to {sp}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="SW-NC-FF reversed-cyclic analysis, Aydin-calibrated lattice")
    p.add_argument("--drift", type=float, default=0.01,
                   help="highest protocol drift level to include (default 0.01 = 1%%, where the "
                        "test specimen peaked)")
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing",
                   help="'crushing' (default) or 'elastic' (Aydin's literal assumption)")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON}). 3.01 gives ~28 "
                        "neighbours per node — Aydin's recommendation for RC, and the redundant "
                        "bracing that keeps cracked tension-side nodes from going singular")
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale the tensile fracture energy (default 1.0 = plain-concrete MC90). "
                        ">1 represents tension stiffening from embedded reinforcement")
    p.add_argument("--steps-per-mm", type=float, default=4.0,
                   help="displacement-control resolution (default 4 steps/mm)")
    p.add_argument("--solver", choices=("static", "dynamic"), default="static",
                   help="'static' (DisplacementControl; stalls near 0.3%% drift) or 'dynamic' "
                        "(dynamic relaxation — rides through the softening instability, D49)")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic solver: drive speed in mm/s at the actuator level (default "
                        f"{QUASI_STATIC_RATE:g}). The inertial and damping forces riding in the "
                        "recorded base shear are LINEAR in this, and so is the run time; 7.6 was "
                        "the value once checked against the static cyclic to 1%%")
    p.add_argument("--periods", type=float, default=None,
                   help="dynamic solver: LEGACY alternative to --rate — derive the speed as "
                        "max|history| / (periods * T1). Discouraged: it scales with the protocol "
                        "amplitude, so the same setting drives the 4%% protocol four times faster "
                        "than the 1%% one. Only used when --rate is given as 0")
    p.add_argument("--damping", type=float, default=DAMPING_RATIO,
                   help=f"dynamic solver: damping ratio (default {DAMPING_RATIO:g} = 5%% of "
                        "critical, the physical value for cracked RC). Raise it only if the solve "
                        "will not stay together — it is a viscous drag that inflates the shear")
    p.add_argument("--no-quasi-static", dest="quasi_static", action="store_false",
                   help="skip the inertia+damping measurement (one extra reaction sum per step)")
    p.add_argument("--draw", dest="draw_model", action="store_true",
                   help="also save the analysis-model lattice figure (wall_drawing.png) before "
                        "the run starts")
    p.add_argument("--gauge-every", type=int, default=200,
                   help="base strain gauge: sample every N steps (default 200). Reversals are "
                        "ALWAYS sampled on top of this, so the loop tips are exact whatever the "
                        "stride. Keep it coarse — a 4%% protocol is ~1.8M steps")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         steps_per_mm=a.steps_per_mm, gf_factor=a.gf_factor, solver=a.solver, periods=a.periods,
         damping=a.damping, rate=a.rate, quasi_static=a.quasi_static, gauge_every=a.gauge_every,
         draw_model=a.draw_model)
