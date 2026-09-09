"""SW-NC-FF shear wall — MONOTONIC pushover of the Aydin-calibrated nonlinear lattice (Stage 0).

The diagnostic that comes before any cyclic run: it converges much further than a reversed-cyclic
analysis, produces the backbone that is directly comparable to the paper's Fig. 14c, and pins down
PEAK LATERAL STRENGTH — the one quantity a perfect-bond lattice can fairly be judged on, because
it follows from section equilibrium (bar areas, fy, fc, axial load, lever arm) rather than from the
bond behaviour this model does not represent.

Strut areas come from Aydin's elastic energy balance (D47), unchanged: that calibration fixes the
initial tangent, which cracking does not alter. Tension softening is regularized by strut length,
Aydin's crack-band principle. `--compression` selects whether concrete may crush (the real
behaviour, and how this specimen failed) or stays elastic (Aydin's literal assumption, his Sec.
2.2) — running both prices that assumption on a specimen governed by crushing.

Paper reference points for SW-NC-FF: peak +212.6 / -209.4 kN at +-1.00% drift; first bar yield at
0.30% drift; ACI flexural prediction V@Mn = 199.8 kN; ultimate drift 1.5%.

`--solver dynamic` imposes the ramp as a transient solve instead, at the same physical drive
settings the cyclic run uses (2.0 mm/s, 5% of critical damping — `specimen.QUASI_STATIC_RATE` /
`DAMPING_RATIO`, D62) so the two are comparable. It reports the inertia + damping riding in the
recorded shear as a percentage of the peak, which is what makes "quasi-static" a measurement here
rather than a claim. This is the cheap place to check a drive setting before paying for a cyclic run.

Output: examples/output/wall/wall_pushover[_elastic].png. Units: N, mm.
Run as `python examples/wall/pushover.py [--compression crushing|elastic] [--drift 0.02]`.
"""

from __future__ import annotations

import json

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import run_pushover, run_pushover_dynamic

import gauge
from build import calibrate, nonlinear_wall_lattice, report_calibration
from specimen import (
    A_SHEAR, DAMPING_RATIO, EPS, GAUGE_H, HORIZON, LW, MESH, OUT, QUASI_STATIC_RATE, base_nodes,
    control_node, gauge_nodes, lateral_loads,
)

# Measured response of SW-NC-FF, for annotation only — never a calibration target.
PAPER_PEAK_KN = 212.6
PAPER_PEAK_DRIFT = 0.0100
PAPER_YIELD_DRIFT = 0.0030
PAPER_ACI_KN = 199.8


def main(*, compression: str = "crushing", drift: float = 0.02, mesh_size: float = MESH,
         horizon: float = HORIZON, gf_factor: float = 1.0, solver: str = "static",
         periods: float | None = None, damping: float = DAMPING_RATIO,
         rate: float = QUASI_STATIC_RATE, quasi_static: bool = True,
         gauge_every: int = 20) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    if gf_factor != 1.0:
        print(f"tension stiffening: Gf scaled x{gf_factor:g}")

    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    target = drift * A_SHEAR
    ctrl, base = control_node(model), base_nodes(model)

    # Vertical strain gauge over the wall base (D63): two aligned node rows GAUGE_H apart, asked
    # for as ONE concatenated list (bottom then top) so a single probe covers both. dof 2 = uy.
    gx, g_bot, g_top = gauge_nodes(model)
    probe = (list(g_bot) + list(g_top), 2)
    print(f"base strain gauge: {len(gx)} columns across x = {gx[0]:.0f}..{gx[-1]:.0f} mm, "
          f"{GAUGE_H:.0f} mm gauge, sampled every {gauge_every} steps")
    if solver == "dynamic":
        # Dynamic relaxation (D22/D46): the whole load-level row is driven by an imposed ramp, slow
        # enough that inertia is negligible, so the recorded curve is quasi-static — but the mass
        # and damping regularize the softening instability that stops the static solver dead.
        drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))
        res = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=target,
                                   drive_nodes=drive, base_nodes=base,
                                   periods_to_target=periods if periods is not None else 12.0,
                                   # --rate 0 hands control back to the legacy --periods route
                                   rate=rate if rate > 0.0 else None,
                                   steps_per_period=30, damping_ratio=damping,
                                   quasi_static=quasi_static, node_history=probe,
                                   node_history_every=gauge_every)
    else:
        res = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                           control_dof=1, dU=target / 400.0, target=target, base_nodes=base,
                           algorithm=("ModifiedNewton", "-initial"), node_history=probe,
                           node_history_every=gauge_every)

    if not res["disp"]:
        raise RuntimeError("pushover produced no steps — the gravity stage failed")
    peak = max(res["shear"])
    at = res["disp"][res["shear"].index(peak)]
    end_drift = res["disp"][-1] / A_SHEAR
    print(f"\ntraced to {end_drift:.3%} drift of {drift:.2%} requested  (converged={res['converged']})")
    print(f"  peak base shear        {peak / 1e3:7.1f} kN at {at / A_SHEAR:.3%} drift")
    print(f"  paper SW-NC-FF         {PAPER_PEAK_KN:7.1f} kN at {PAPER_PEAK_DRIFT:.2%} drift"
          f"   -> model/test = {peak / 1e3 / PAPER_PEAK_KN:.3f}")
    print(f"  ACI flexural V@Mn      {PAPER_ACI_KN:7.1f} kN"
          f"                     -> model/ACI  = {peak / 1e3 / PAPER_ACI_KN:.3f}")
    if res.get("dynamic"):
        # How much of that shear is the wall and how much is the solver: the inertia + damping the
        # drive and base reactions fail to balance (D62). This runner marches on Newmark, where the
        # residual DOES track the error in the recorded shear (measured: a 27% residual share went
        # with a 2.6% overshoot of the static answer, 7% with 0.5%) — so unlike the cyclic script's
        # HHT reading it is worth a rough threshold. It is what --rate and --damping are set by:
        # at the defaults this wall reads ~4.0% at 0.3% drift, against 10.3% at the old zeta = 0.8.
        # Split the START-UP transient from the rest. The ramp begins from rest with a velocity
        # step, which rings once at near-zero drift where the shear is small; it is an artifact of
        # starting, not a contamination of the traced curve, and it scales with --rate. What
        # contaminates the curve is the steady part, so that is the headline. (The cyclic run
        # cannot make this split — every one of its reversals is another velocity step.)
        dyn = res["dynamic"]
        n_start = min(len(dyn) - 1, 10 * 30)          # ~10 fundamental periods at 30 steps each
        start = max(abs(d) for d in dyn[:n_start + 1])
        tail = dyn[n_start:] or dyn
        steady = max(abs(d) for d in tail)
        rms = (sum(d * d for d in tail) / len(tail)) ** 0.5
        share = steady / abs(peak) if peak else float("inf")
        flag = "" if share < 0.06 else "   [above the 4.0% baseline — check --damping before --rate]"
        print(f"  inertia + damping      {steady / 1e3:7.1f} kN, {rms / 1e3:.1f} kN rms "
              f"= {share:.1%} of peak shear{flag}")
        print(f"  start-up transient     {start / 1e3:7.1f} kN once, at near-zero drift "
              f"(scales with --rate; not part of the traced curve)")
        print(f"  drive                  {res['rate']:7.4g} mm/s, damping {damping:.0%} of "
              f"critical, T1 = {res['T1']:.4f} s")
    if res["shear"][-1] < peak:
        print(f"  post-peak degradation  {(1 - res['shear'][-1] / peak) * 100:5.1f}% by "
              f"{end_drift:.2%} drift")

    drift_pct = [u / A_SHEAR * 100.0 for u in res["disp"]]
    span = [0.0, max(drift_pct)]
    curves = [
        {"disp": drift_pct, "shear": [s / 1e3 for s in res["shear"]],
         "label": f"lattice backbone, compression={compression} (peak {peak / 1e3:.0f} kN)",
         "style": {"color": "C0", "lw": 2}},
        {"disp": span, "shear": [PAPER_PEAK_KN] * 2,
         "label": f"test peak, SW-NC-FF ({PAPER_PEAK_KN:.0f} kN)",
         "style": {"color": "C3", "ls": "--", "lw": 1.5}},
        {"disp": span, "shear": [PAPER_ACI_KN] * 2,
         "label": f"ACI V@Mn ({PAPER_ACI_KN:.0f} kN)",
         "style": {"color": "0.55", "ls": ":", "lw": 1.5}},
    ]
    stem = ("wall_pushover" + ("_elastic" if compression == "elastic" else "")
            + ("_dynamic" if solver == "dynamic" else ""))
    savepath = OUT / f"{stem}.png"
    viz.figure_pushover(
        curves, savepath=str(savepath), xlabel="drift ratio (%)", ylabel="base shear (kN)",
        title=f"SW-NC-FF — monotonic backbone, Aydin-calibrated lattice (compression={compression})",
    )
    print(f"\nsaved backbone to {savepath}")

    # Vertical strain profile across the base section at a spread of drift levels (D63), plus the
    # raw gauge series as JSON so the figure can be redrawn without paying for the run again.
    g_drift, g_strain = gauge.reduce_history(res, gx)
    if g_drift:
        reached = max(g_drift)
        levels = [lv for lv in (0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0) if lv <= reached]
        idx = gauge.select_levels(g_drift, levels or [reached], both_directions=False)
        sp = OUT / f"{stem}_strain_profile.png"
        rows = gauge.figure(
            g_drift, g_strain, gx, indices=idx, savepath=sp,
            title="SW-NC-FF — vertical strain across the base section")
        print(f"\nbase strain gauge: {len(g_drift)} profiles recorded, {len(idx)} drawn")
        gauge.print_table(rows)
        print(f"saved strain profile to {sp}")
        datapath = OUT / f"{stem}_data.json"
        datapath.write_text(json.dumps({
            "drift_pct": [u / A_SHEAR * 100.0 for u in res["disp"]],
            "shear_kN": [v / 1e3 for v in res["shear"]],
            "compression": compression, "solver": solver, "mesh": mesh_size, "horizon": horizon,
            "gauge": gauge.payload(res, gx),
        }))
        print(f"saved raw response to {datapath}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="SW-NC-FF monotonic pushover, Aydin-calibrated lattice")
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing",
                   help="'crushing' (default, real Concrete02 backbone) or 'elastic' "
                        "(Aydin's literal assumption: concrete never crushes)")
    p.add_argument("--drift", type=float, default=0.02, help="target drift ratio (default 0.02)")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON}; 3.01 = ~28 neighbours)")
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale the tensile fracture energy (tension stiffening); default 1.0")
    p.add_argument("--solver", choices=("static", "dynamic"), default="static",
                   help="'static' (DisplacementControl) or 'dynamic' (dynamic relaxation — rides "
                        "through the softening instability that stops the static solver, D22/D46)")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic solver: drive speed in mm/s at the actuator level (default "
                        f"{QUASI_STATIC_RATE:g}). The inertial and damping forces riding in the "
                        "recorded base shear are LINEAR in this, and so is the run time")
    p.add_argument("--periods", type=float, default=None,
                   help="dynamic solver: LEGACY alternative to --rate — derive the speed as "
                        "|target| / (periods * T1). Discouraged: it scales with the target, so the "
                        "same setting drives a larger pushover faster. Used only with --rate 0")
    p.add_argument("--damping", type=float, default=DAMPING_RATIO,
                   help=f"dynamic solver: damping ratio (default {DAMPING_RATIO:g} = 5%% of "
                        "critical, the physical value for cracked RC). Raise it only if the solve "
                        "will not stay together — it is a viscous drag that inflates the shear")
    p.add_argument("--no-quasi-static", dest="quasi_static", action="store_false",
                   help="skip the inertia+damping measurement (one extra reaction sum per step)")
    p.add_argument("--gauge-every", type=int, default=20,
                   help="base strain gauge: sample every N steps (default 20)")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         gf_factor=a.gf_factor, solver=a.solver, periods=a.periods, damping=a.damping,
         rate=a.rate, quasi_static=a.quasi_static, gauge_every=a.gauge_every)
