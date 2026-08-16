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

Output: examples/output/wall/wall_cyclic[_elastic].png. Units: N, mm.
Run as `python examples/wall/cyclic.py [--drift 0.01] [--compression crushing|elastic]`.
"""

from __future__ import annotations

import json
import time

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import cyclic_protocol, run_cyclic, run_cyclic_dynamic

from build import calibrate, nonlinear_wall_lattice
from specimen import (
    A_SHEAR, EPS, HORIZON, LW, MESH, OUT, base_nodes, control_node, lateral_loads, protocol,
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
         solver: str = "static", periods: float = 48.0, damping: float = 0.8,
         rate: float | None = None) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    print(f"Aydin energy balance: A_t = {cal.area:,.1f} mm^2 "
          f"= {cal.area / (mesh_size * 200.0):.4f} * (thickness * mesh)   "
          f"[horizon {horizon}, lattice nu = {cal.nu_consistent:.3f}]")
    if gf_factor != 1.0:
        print(f"tension stiffening: Gf scaled x{gf_factor:g} (softening gentler than plain concrete)")

    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    peaks, cycles = protocol(drift)
    history = cyclic_protocol(peaks, cycles_per_level=cycles)
    print(f"protocol to {drift:.2%} drift: {len(peaks)} levels "
          f"({', '.join(f'{p / A_SHEAR:.2%}' for p in peaks)}), {sum(cycles)} cycles, "
          f"{len(history)} reversals")

    ctrl, base = control_node(model), base_nodes(model)
    t0 = time.time()
    if solver == "dynamic":
        # Dynamic relaxation through the whole reversing history (D49). A static path-follower
        # cannot cross the softening instability of a cracking lattice; mass and damping can.
        drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))

        def report(i, n, u, s):
            print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+7.3%}  shear {s / 1e3:+8.1f} kN"
                  f"  [{time.time() - t0:6.0f}s]", flush=True)

        res = run_cyclic_dynamic(model, control_node=ctrl, control_dof=1, history=history,
                                 drive_nodes=drive, base_nodes=base, periods_to_peak=periods,
                                 rate=rate, steps_per_period=30, damping_ratio=damping,
                                 progress=report)
        print(f"  T1 = {res['T1']:.4f} s, drive rate = {res['rate']:.4g} mm/s, "
              f"{res['steps']} steps")
    else:
        res = run_cyclic(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                         control_dof=1, history=history, dU=1.0 / steps_per_mm, base_nodes=base)
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
        "drift_target": drift, "mesh": mesh_size, "horizon": horizon, "rate": rate,
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
    p.add_argument("--periods", type=float, default=48.0,
                   help="dynamic solver: fundamental periods to traverse the largest amplitude. "
                        "Higher = slower = less inertial contamination of the recorded base shear")
    p.add_argument("--damping", type=float, default=0.8, help="dynamic solver: damping ratio")
    p.add_argument("--rate", type=float, default=None,
                   help="dynamic solver: drive speed in mm/s, set DIRECTLY (preferred over "
                        "--periods, which scales with protocol amplitude and so silently drives a "
                        "larger protocol faster). 7.6 mm/s matched the static cyclic within 1%%")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         steps_per_mm=a.steps_per_mm, gf_factor=a.gf_factor, solver=a.solver, periods=a.periods,
         damping=a.damping, rate=a.rate)
