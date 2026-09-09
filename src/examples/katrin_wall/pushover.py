"""WSH3 shear wall — MONOTONIC pushover of the Aydin-calibrated nonlinear lattice.

The diagnostic that comes before any cyclic run: it converges much further than a reversed-cyclic
analysis, produces the backbone directly comparable to the loop-tip backbone digitized from the
paper's Fig. 7c, and pins down PEAK LATERAL STRENGTH — which follows from section equilibrium (bar
areas, f_y, f'_c, axial load, lever arm), the part of the response a lattice is best placed to get
right.

Strut areas come from Aydin's elastic energy balance (D47), unchanged: that calibration fixes the
initial tangent, which cracking does not alter. Tension softening is regularized by strut length,
Aydin's crack-band principle. `--compression` selects whether concrete may crush (the real
behaviour) or stays elastic (Aydin's literal assumption, his Sec. 2.2).

WHAT IS FAIR TO JUDGE, and it is a LONGER list than for SW-NC-FF. That wall's displacement was 74%
rocking on debonded PLAIN bars, which capped the comparison at strength and stiffness. WSH3 uses
DEFORMED bars that stay bonded — its displacement decomposes into flexure and shear with only a
small fixed-end component (Fig. 9a) — so peak strength, stiffness, curvature, neutral-axis depth
AND drift capacity are all fair here. What remains out of reach is the FAILURE MODE: `Steel02` has
neither bar buckling nor fracture, and those are what ended the test at 1.79% drift.

Paper reference points (Dazio, Beyer & Bachmann 2009): V_max 454 kN (Table 5); backbone plateau
reached near 1.35% drift; first yield of the outer bar at 0.25% drift, delta_y (3/4-rule) 0.34%;
M_cr 527 kN.m -> V_cr 116 kN (Table 5); ultimate drift 2.03%, but the 20%-drop failure criterion was
never met, so that is a LOWER bound.

`--solver dynamic` imposes the ramp as a transient solve at the same drive settings the cyclic run
uses (`specimen.QUASI_STATIC_RATE` / `DAMPING_RATIO`), so the two are comparable, and reports the
inertia + damping riding in the recorded shear as a percentage of the peak — which is what makes
"quasi-static" a measurement rather than a claim. THAT MEASUREMENT MATTERS MORE HERE than it did for
SW-NC-FF: this test ran at 1.2-3.6 mm per MINUTE, so the drive is ~130x the real speed and the rate
itself licences nothing.

Output: examples/output/katrin_wall/wsh3_pushover[_elastic][_dynamic].png. Units: N, mm.
Run as `python examples/katrin_wall/pushover.py [--compression crushing|elastic] [--drift 0.025]`.
"""

from __future__ import annotations

import json

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import run_pushover, run_pushover_dynamic

import gauge
from build import calibrate, nonlinear_wall_lattice, report_calibration
from specimen import (
    A_SHEAR, DAMPING_RATIO, EPS, GAUGE_H, GAUGE_Y0, HORIZON, LVDT_ROWS, LW, MESH, OUT, PROTOCOL_PEAKS,
    QUASI_STATIC_RATE, SPECIMEN, base_nodes, control_node, gauge_nodes, gauge_probe, lateral_loads,
)

# Measured response of WSH3, for annotation only — never a calibration target.
PAPER_PEAK_KN = 454.0            # Table 5
PAPER_PEAK_DRIFT = 0.0135        # backbone plateau, digitized Fig. 7c (mu = 4)
PAPER_YIELD_DRIFT = 0.0025       # first yield of the outer longitudinal bar, Table 4
PAPER_VCR_KN = 115.6             # M_cr 527 kN.m / L_v 4560, Table 5
PAPER_ULT_DRIFT = 0.0203         # Table 4 — a LOWER bound, the 20%-drop criterion was never met


def main(*, compression: str = "crushing", drift: float = 0.025, mesh_size: float = MESH,
         horizon: float = HORIZON, gf_factor: float = 1.0,
         steel_rupture=None, concrete_residual: float = 0.2, solver: str = "static",
         periods: float | None = None, damping: float = DAMPING_RATIO,
         rate: float = QUASI_STATIC_RATE, quasi_static: bool = True,
         gauge_every: int = 20) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    if gf_factor != 1.0:
        print(f"tension stiffening: Gf scaled x{gf_factor:g}")

    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor,
                                   steel_rupture=steel_rupture,
                                   concrete_residual=concrete_residual)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    target = drift * A_SHEAR
    ctrl, base = control_node(model), base_nodes(model)

    # Vertical strain gauge over the wall base: two aligned node rows GAUGE_H apart, asked for as
    # ONE concatenated list (bottom then top) so a single probe covers both. dof 2 = uy.
    gx, g_rows = gauge_nodes(model)
    probe = (gauge_probe(g_rows), 2)
    print(f"LVDT chain: {len(gx)} columns across x = {gx[0]:.0f}..{gx[-1]:.0f} mm x "
          f"{len(g_rows)} rows to y = {LVDT_ROWS[-1]:.0f} mm; base profile over "
          f"y = {GAUGE_Y0:.0f}-{GAUGE_Y0 + GAUGE_H:.0f} mm, sampled every {gauge_every} steps")

    if solver == "dynamic":
        drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))
        res = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=target,
                                   drive_nodes=drive, base_nodes=base,
                                   periods_to_target=periods if periods is not None else 12.0,
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
    print(f"  test WSH3              {PAPER_PEAK_KN:7.1f} kN at {PAPER_PEAK_DRIFT:.2%} drift"
          f"   -> model/test = {peak / 1e3 / PAPER_PEAK_KN:.3f}")
    print(f"  cracking V_cr (EC2)    {PAPER_VCR_KN:7.1f} kN"
          f"                     -> model/V_cr = {peak / 1e3 / PAPER_VCR_KN:.3f}")
    if res.get("dynamic"):
        # How much of that shear is the wall and how much is the solver (D62). Split the START-UP
        # transient from the rest: the ramp begins from rest with a velocity step that rings once at
        # near-zero drift, an artifact of starting rather than a contamination of the traced curve.
        dyn = res["dynamic"]
        n_start = min(len(dyn) - 1, 10 * 30)          # ~10 fundamental periods at 30 steps each
        start = max(abs(d) for d in dyn[:n_start + 1])
        tail = dyn[n_start:] or dyn
        steady = max(abs(d) for d in tail)
        rms = (sum(d * d for d in tail) / len(tail)) ** 0.5
        share = steady / abs(peak) if peak else float("inf")
        flag = "" if share < 0.06 else "   [check --damping before --rate: D64]"
        print(f"  inertia + damping      {steady / 1e3:7.1f} kN, {rms / 1e3:.1f} kN rms "
              f"= {share:.1%} of peak shear{flag}")
        print(f"  start-up transient     {start / 1e3:7.1f} kN once, at near-zero drift "
              f"(scales with --rate; not part of the traced curve)")
        print(f"  drive                  {res['rate']:7.4g} mm/s, damping {damping:.0%} of "
              f"critical, T1 = {res['T1']:.4f} s")
        print(f"                         the TEST ran at 0.02-0.06 mm/s (1.2-3.6 mm/min, Fig. 6), "
              f"so this is ~{res['rate'] / 0.06:.0f}x the real speed")
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
         "label": f"test peak, {SPECIMEN} ({PAPER_PEAK_KN:.0f} kN)",
         "style": {"color": "C3", "ls": "--", "lw": 1.5}},
        {"disp": span, "shear": [PAPER_VCR_KN] * 2,
         "label": f"EC2 cracking V$_{{cr}}$ ({PAPER_VCR_KN:.0f} kN)",
         "style": {"color": "0.55", "ls": ":", "lw": 1.5}},
    ]
    stem = ("wsh3_pushover" + ("_elastic" if compression == "elastic" else "")
            + ("_dynamic" if solver == "dynamic" else ""))
    savepath = OUT / f"{stem}.png"
    viz.figure_pushover(
        curves, savepath=str(savepath), xlabel="drift ratio (%)", ylabel="base shear (kN)",
        title=f"{SPECIMEN} — monotonic backbone, Aydin-calibrated lattice "
              f"(compression={compression})",
    )
    print(f"\nsaved backbone to {savepath}")

    # Vertical strain profile across the base section, plus the raw gauge series as JSON so the
    # figure can be redrawn without paying for the run again.
    g_drift, g_strain = gauge.reduce_history(res, gx)
    if g_drift:
        reached = max(g_drift)
        # The PROTOCOL levels, so a pushover profile and a cyclic profile are drawn at the same
        # drifts and can be laid side by side.
        levels = [p / A_SHEAR * 100.0 for p in PROTOCOL_PEAKS if p / A_SHEAR * 100.0 <= reached]
        idx = gauge.select_levels(g_drift, levels or [reached], both_directions=False)
        sp = OUT / f"{stem}_strain_profile.png"
        rows = gauge.figure(
            g_drift, g_strain, gx, indices=idx, savepath=sp,
            title=f"{SPECIMEN} — vertical strain across the base section")
        print(f"\nbase strain gauge: {len(g_drift)} profiles recorded, {len(idx)} drawn")
        gauge.print_table(rows)
        print(f"saved strain profile to {sp}")

        # Curvature over height — the model's counterpart to the paper's Fig. 11c, with the
        # base-curvature fit of Fig. 15d marked on each profile.
        c_drift, mids, curv = gauge.reduce_curvature(res, gx)
        cp = OUT / f"{stem}_curvature.png"
        crows = gauge.figure_curvature(
            c_drift, mids, curv, indices=idx, savepath=cp,
            title=f"{SPECIMEN} — curvature over height (cf. paper Fig. 11c)")
        gauge.print_curvature_table(crows)
        print(f"saved curvature profile to {cp}")

        datapath = OUT / f"{stem}_data.json"
        datapath.write_text(json.dumps({
            "drift_pct": drift_pct,
            "shear_kN": [v / 1e3 for v in res["shear"]],
            "compression": compression, "solver": solver, "mesh": mesh_size, "horizon": horizon,
            # PROVENANCE: a record that does not carry the switches it was run with cannot
            # be compared to one that was (D70/D80). Gf is not a neutral knob (x2 = +14.2%
            # base shear, D87) and the failure switches decide whether a capacity exists.
            "gf_factor": gf_factor, "steel_rupture": steel_rupture,
            "concrete_residual": concrete_residual,
            "converged": res["converged"],
            "gauge": gauge.payload(res, gx),
        }))
        print(f"saved raw response to {datapath}")


def _rupture_strain(v: str):
    """`--steel-rupture` takes the literal 'measured' or a strain."""
    return v if v == "measured" else float(v)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="WSH3 monotonic pushover, Aydin-calibrated lattice")
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing",
                   help="'crushing' (default, real Concrete02 backbone) or 'elastic' "
                        "(Aydin's literal assumption: concrete never crushes)")
    p.add_argument("--drift", type=float, default=0.025,
                   help="target drift ratio (default 0.025 — just past the test's 2.03%%)")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON}; 3.01 = ~28 neighbours)")
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale the tensile fracture energy (tension stiffening); default 1.0")
    p.add_argument("--steel-rupture", type=_rupture_strain, default=None,
                   help="bar rupture strain: 'measured' uses each bar's own published A_gt (phi12 7.69%%, phi8 7.34%%, phi6 6.45%%), a number applies one strain to all, omit to leave bars unbreakable. The test's corner bar ruptured at 1.79%% drift, which is what this can be validated against (D93)")
    p.add_argument("--concrete-residual", type=float, default=0.2,
                   help="floor the crushing strength at this fraction of fc (default 0.20, D22's plateau). 0.0 lets struts actually crush — available under the dynamic solver, which has no tangent to go singular")
    p.add_argument("--solver", choices=("static", "dynamic"), default="static",
                   help="'static' (DisplacementControl) or 'dynamic' (dynamic relaxation — rides "
                        "through the softening instability that stops the static solver, D22/D46)")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic solver: drive speed in mm/s at the actuator level (default "
                        f"{QUASI_STATIC_RATE:g}). The TEST ran at 0.02-0.06 mm/s, which is not "
                        "reachable; the licence is the measured residual, not this number")
    p.add_argument("--periods", type=float, default=None,
                   help="dynamic solver: LEGACY alternative to --rate — derive the speed as "
                        "|target| / (periods * T1). Used only with --rate 0")
    p.add_argument("--damping", type=float, default=DAMPING_RATIO,
                   help=f"dynamic solver: damping ratio (default {DAMPING_RATIO:g}). This is the "
                        "knob that controls the contamination, not --rate (D64)")
    p.add_argument("--no-quasi-static", dest="quasi_static", action="store_false",
                   help="skip the inertia+damping measurement (one extra reaction sum per step)")
    p.add_argument("--gauge-every", type=int, default=20,
                   help="base strain gauge: sample every N steps (default 20)")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         gf_factor=a.gf_factor, steel_rupture=a.steel_rupture,
         concrete_residual=a.concrete_residual, solver=a.solver, periods=a.periods, damping=a.damping,
         rate=a.rate, quasi_static=a.quasi_static, gauge_every=a.gauge_every)
