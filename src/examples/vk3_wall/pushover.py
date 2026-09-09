"""VK3 bridge pier — MONOTONIC pushover of the Aydin-calibrated nonlinear lattice.

The diagnostic that comes before any cyclic run: it converges much further than a reversed-cyclic
analysis, produces a backbone directly comparable to the loop tips digitized from Fig. 5.13, and
pins down PEAK LATERAL STRENGTH — which follows from section equilibrium (bar areas, f_y, f'_c,
axial load, lever arm), the part of the response a lattice is best placed to get right.

VK3's peak strength is NOT tabulated anywhere in the chapter; the digitized loops supply it
(+891 / -876 kN, peaking at ~0.95% drift), against the chapter's own predicted nominal capacity of
851 kN. So the comparison target here comes from `digitize.py`, not from a printed table.

WHAT IS FAIR TO JUDGE. VK3's bars are continuous and deformed, so perfect bond is a reasonable
idealization and strength, stiffness, the shear/flexure split and where damage localizes are all
fair. What is NOT reachable is the failure itself: the pier failed by sliding on crossed diagonal
cracks after bar buckling and core spalling destroyed the base compression zone, and a shared-node
truss lattice has no element for aggregate interlock, no bar buckling in `Steel02`, and no hoop that
can lose its 90-degree hook anchorage the way the real ones could.

`--element-groups` decomposes the base shear into the orthogonal and diagonal strut families. That
is the direct test of the chapter's own explanation of the failure (Sec. 5.4): the diagonal strut
lost the compression zone that supported it. `compression_cube` measured exactly this kind of load-
path collapse (D55, diagonals 28% -> 0%), so the instrument already exists.

Output: examples/output/vk3_wall/vk3_pushover[_elastic][_dynamic].png. Units: N, mm.
Run from src/: python examples/vk3_wall/pushover.py [--drift 0.02] [--gf-factor 2] [--draw]
"""

from __future__ import annotations

import json
import time

from rclattice import viz
from rclattice.opensees import run_pushover, run_pushover_dynamic

import gauge
import response as resp
from build import calibrate, nonlinear_pier_lattice, report_calibration, strut_life
from specimen import (
    A_SHEAR, DAMPING_RATIO, EC, EPSC0, FT, HORIZON, LVDT_ROWS, MEASURED_LOOPS, MESH, OUT,
    PROTOCOL_PEAKS,
    QUASI_STATIC_RATE, SPECIMEN, base_nodes, control_node, drive_nodes, gauge_nodes, gauge_probe,
    lateral_loads, run_dir, set_steel_R0, set_steel_b,
)
import testdata as td

PAPER_FN_KN = td.PREDICTION["Fn_kN"]

# Hand section analysis of the AS-BUILT model section (bar areas/positions + axial load), which
# reproduces the chapter's own nominal moment capacity to 0.4% (2798 vs 2808 kN.m). A valid
# pushover CANNOT exceed this, so it is the physical bound `response.summarize` checks against.
SECTION_CAPACITY_KN = 848.0
              # 851 kN, predicted nominal capacity


def measured_peak():
    """`(push_kN, pull_kN, drift_at_peak)` from the digitized loops, or None."""
    if not MEASURED_LOOPS.exists():
        return None
    import numpy as np
    d = np.load(MEASURED_LOOPS)
    u, v = d["disp_mm"], d["load_kN"]
    i = int(np.argmax(v))
    return float(v.max()), float(v.min()), float(u[i]) / A_SHEAR


def main(*, compression: str = "crushing", drift: float = 0.02, mesh_size: float = MESH,
         horizon: float = HORIZON, gf_factor: float = 1.0, solver: str = "static",
         periods: float | None = None, damping: float = DAMPING_RATIO,
         rate: float = QUASI_STATIC_RATE, quasi_static: bool = True, gauge_every: int = 20,
         draw_model: bool = False, steel_r0: float | None = None,
         steel_b: float | None = None, capture: bool = False,
         steps_per_period: int = 30, integrator: str = "newmark",
         corot: bool = False) -> None:
    if steel_r0 is not None:
        set_steel_R0(steel_r0)
    if steel_b is not None:
        set_steel_b(steel_b)
    tag = f"{solver}" + (f"_gf{gf_factor:g}" if gf_factor != 1.0 else "")
    if corot:
        tag += "_corot"
    if steel_r0 is not None:
        tag += f"_R0-{steel_r0:g}"
    if steel_b is not None:
        tag += f"_DIAG-b{steel_b:g}"
    out = run_dir("pushover", tag=tag + ("_elastic" if compression == "elastic" else ""))
    print(f"output directory: {out}")

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    integ = ("HHT", 0.7) if integrator == "hht" else ("Newmark", 0.5, 0.25)
    if solver == "dynamic":
        print(f"integration: {integ[0]}{integ[1:]}, dt = T1/{steps_per_period}"
              + ("" if steps_per_period == 30 else "   *** REFINED from the default 30 ***"))
    print(f"strut life eps_ult/eps_cr = {strut_life(gf_factor, mesh_size):.1f}"
          + (f"  (Gf scaled x{gf_factor:g})" if gf_factor != 1.0 else "  (plain-concrete MC90)"))
    from specimen import S_LONG as _SL
    print(f"steel yield corner: Steel02 R0 = {_SL.R0:g} "
          f"({'default' if steel_r0 is None else 'OVERRIDDEN — smoother transition'}), "
          f"b = {_SL.b:.5f}"
          + ("  (measured, untouched)" if steel_b is None else
             "  *** DIAGNOSTIC OVERRIDE — measured value is 0.00466; this result is NOT "
             "a production number ***"))

    model = nonlinear_pier_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor,
                                   strut_element="corotTruss" if corot else "Truss")
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"strut element: {'corotTruss (geometrically consistent)' if corot else 'Truss (small displacement)'}")
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    if draw_model:
        import draw as draw_mod
        draw_mod.main(mesh_size=mesh_size, horizon=horizon, nonlinear=True)

    target = drift * A_SHEAR
    ctrl, base = control_node(model), base_nodes(model)

    # ONE probe feeds both instruments: dofs (1, 2) = (u_x, u_y) at every gauge node, recorded
    # dof-major (D68). The vertical chain reads the u_y block; the diagonal chain needs both.
    gx, g_rows = gauge_nodes(model)
    probe = (gauge_probe(g_rows), (1, 2))
    print(f"gauge chains: {len(gx)} columns across x = {gx[0]:.0f}..{gx[-1]:.0f} mm x "
          f"{len(g_rows)} rows to y = {LVDT_ROWS[-1]:.0f} mm, sampled every {gauge_every} steps")

    t0 = time.time()

    def report(i, n, u, s):
        print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+7.3%}  shear {s / 1e3:+8.1f} kN"
              f"  [{time.time() - t0:6.0f}s]", flush=True)

    if solver == "dynamic":
        res = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=target,
                                   drive_nodes=drive_nodes(model), base_nodes=base,
                                   periods_to_target=periods if periods is not None else 12.0,
                                   rate=rate if rate > 0.0 else None,
                                   steps_per_period=steps_per_period, damping_ratio=damping,
                                   quasi_static=quasi_static, node_history=probe,
                                   node_history_every=gauge_every, progress=report,
                                   capture=capture, integrator=integ)
    else:
        res = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                           control_dof=1, dU=target / 400.0, target=target, base_nodes=base,
                           algorithm=("ModifiedNewton", "-initial"), node_history=probe,
                           node_history_every=gauge_every)

    if not res["disp"]:
        raise RuntimeError("pushover produced no steps — the gravity stage failed")
    end_drift = res["disp"][-1] / A_SHEAR
    summary = resp.summarize([u / A_SHEAR * 100.0 for u in res["disp"]],
                             [s_ / 1e3 for s_ in res["shear"]],
                             capacity_kN=SECTION_CAPACITY_KN, cyclic=False)
    peak = summary.reliable.kN * 1e3          # the RELIABLE peak, never raw max(shear) (D70)
    at = summary.reliable.drift_pct / 100.0 * A_SHEAR
    print(f"\ntraced to {end_drift:.3%} drift of {drift:.2%} requested  "
          f"(converged={res['converged']})")
    print(summary.report())
    meas = measured_peak()
    if meas:
        print(f"  test VK3 (digitized)   {meas[0]:7.1f} kN at {meas[2]:.2%} drift"
              f"   -> model/test = {peak / 1e3 / meas[0]:.3f}")
    print(f"  predicted F_n (Tab 5.6){PAPER_FN_KN:7.1f} kN"
          f"                     -> model/F_n = {peak / 1e3 / PAPER_FN_KN:.3f}")

    if res.get("dynamic"):
        dyn = res["dynamic"]
        n_start = min(len(dyn) - 1, 10 * 30)
        start = max(abs(d) for d in dyn[:n_start + 1])
        tail = dyn[n_start:] or dyn
        steady = max(abs(d) for d in tail)
        share = steady / abs(peak) if peak else float("inf")
        print(f"  inertia + damping      {steady / 1e3:7.1f} kN = {share:.1%} of peak shear"
              + ("" if share < 0.06 else "   [check --damping before --rate: D64]"))
        print(f"  start-up transient     {start / 1e3:7.1f} kN once, at near-zero drift")
        print(f"  drive                  {res['rate']:7.4g} mm/s, damping {damping:.0%}, "
              f"T1 = {res['T1']:.4f} s")
        print(f"                         the chapter states NO actuator rate, so the licence is "
              f"this residual and nothing else")
    if res["shear"][-1] < peak:
        print(f"  post-peak degradation  {(1 - res['shear'][-1] / peak) * 100:5.1f}% by "
              f"{end_drift:.2%} drift")

    drift_pct = [u / A_SHEAR * 100.0 for u in res["disp"]]
    span = [0.0, max(drift_pct)]
    curves = [
        {"disp": drift_pct, "shear": [s / 1e3 for s in res["shear"]],
         "label": f"lattice backbone, compression={compression} (peak {peak / 1e3:.0f} kN)",
         "style": {"color": "C0", "lw": 2}},
        {"disp": span, "shear": [PAPER_FN_KN] * 2,
         "label": f"predicted F$_n$ (Tab. 5.6, {PAPER_FN_KN:.0f} kN)",
         "style": {"color": "0.55", "ls": ":", "lw": 1.5}},
    ]
    if meas:
        curves.append({"disp": span, "shear": [meas[0]] * 2,
                       "label": f"test peak, digitized ({meas[0]:.0f} kN)",
                       "style": {"color": "C3", "ls": "--", "lw": 1.5}})
    savepath = out / "backbone.png"
    viz.figure_pushover(curves, savepath=str(savepath), xlabel="drift ratio (%)",
                        ylabel="base shear (kN)",
                        title=f"{SPECIMEN} — monotonic backbone, Aydin-calibrated lattice "
                              f"(compression={compression})")
    print(f"\nsaved backbone to {savepath}")

    if capture and res.get("disps_peak") is not None:
        # The direct observation, not another inference: axial strain of every strut at peak shear
        # and at the last converged step. eps_crack = ft/E is where a strut first cracks; the
        # regularized law takes it to zero stress at `strut_life` times that, so anything beyond
        # the second contour has NO tensile stiffness left at all.
        from rclattice.viz import figure_damage
        eps_cr = FT / EC
        life = strut_life(gf_factor, mesh_size)
        panels = [(f"at peak shear ({peak / 1e3:.0f} kN)", model, res["disps_peak"]),
                  (f"at the last step ({end_drift:.3%} drift)", model, res["disps_final"])]
        dmg = out / "damage.png"
        figure_damage(panels, eps_crack=eps_cr, eps_crush=-EPSC0, savepath=str(dmg),
                      suptitle=f"{SPECIMEN} — strut damage (crack at {eps_cr:.2e}, "
                               f"zero tensile stress at {life:.0f}x that = {life * eps_cr:.2e})")
        print(f"saved damage map to {dmg}")
        for label, _m, d in panels:
            from rclattice.viz import strut_strains
            els, eps = strut_strains(model, d)
            conc = [(e, v) for e, v in zip(els, eps)
                    if e.kind not in ("longitudinal", "stirrup")]
            n = len(conc)
            cracked = sum(1 for _e, v in conc if v > eps_cr)
            dead = sum(1 for _e, v in conc if v > life * eps_cr)
            crushed = sum(1 for _e, v in conc if v < -EPSC0)
            print(f"  {label}: {cracked:,}/{n:,} struts cracked ({cracked / n:.1%}), "
                  f"{dead:,} with ZERO tensile stiffness ({dead / n:.1%}), "
                  f"{crushed:,} past epsc0 in compression ({crushed / n:.1%})")

    g_drift, g_strain = gauge.reduce_history(res, gx)
    if g_drift:
        reached = max(g_drift)
        levels = [p / A_SHEAR * 100.0 for p in PROTOCOL_PEAKS if p / A_SHEAR * 100.0 <= reached]
        idx = gauge.select_levels(g_drift, levels or [reached], both_directions=False)
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
        gauge.print_curvature_table(crows)
        print(f"saved curvature profile to {cp}")

        # The shear / flexure split — Fig. 5.19-right, and the reason this study exists.
        k_drift, tops, comps = gauge.reduce_components(res, gx)
        comp_path = out / "components.png"
        crow = gauge.figure_components(k_drift, tops, comps, savepath=comp_path,
                                       title=f"{SPECIMEN} — deformation components "
                                             f"(cf. chapter Fig. 5.19 right)")
        print(f"\ndeformation components (test: shear "
              f"{td.COMPONENTS['shear_pct_range'][0]:.0f}-"
              f"{td.COMPONENTS['shear_pct_range'][1]:.0f}% of the top displacement):")
        gauge.print_components_table(crow, every=max(1, len(crow) // 12))
        print(f"saved components to {comp_path}")

        datapath = out / "data.json"
        datapath.write_text(json.dumps({
            "drift_pct": drift_pct, "shear_kN": [v / 1e3 for v in res["shear"]],
            "compression": compression, "solver": solver, "mesh": mesh_size, "horizon": horizon,
            "gf_factor": gf_factor, "converged": res["converged"],
            "integrator": list(res.get("integrator", [])),
            "strut_element": "corotTruss" if corot else "Truss",
            "steps_per_period": steps_per_period, "damping": damping,
            "response_summary": summary.to_dict(),
            "gauge": gauge.payload(res, gx),
        }))
        print(f"saved raw response to {datapath}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="VK3 monotonic pushover, Aydin-calibrated lattice")
    p.add_argument("--compression", choices=("crushing", "elastic"), default="crushing")
    p.add_argument("--drift", type=float, default=0.02,
                   help="target drift ratio (default 0.02 — past the test's 1.59%% failure)")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--gf-factor", type=float, default=1.0,
                   help="scale the tensile fracture energy (tension stiffening). VK3's strut life "
                        "is 14.2, close to WSH3's 13.4, which needed 2.0 to run (D67)")
    p.add_argument("--solver", choices=("static", "dynamic"), default="static")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"dynamic solver: drive speed mm/s (default {QUASI_STATIC_RATE:g}). The "
                        "chapter states no actuator rate at all")
    p.add_argument("--periods", type=float, default=None)
    p.add_argument("--damping", type=float, default=DAMPING_RATIO)
    p.add_argument("--no-quasi-static", dest="quasi_static", action="store_false")
    p.add_argument("--gauge-every", type=int, default=20)
    p.add_argument("--steel-r0", type=float, default=None,
                   help="Steel02 yield-corner sharpness R0 (default 18). LOWER is smoother; "
                        "OpenSees recommends 10-20. Numerical only — f_y and b are measured and "
                        "are never changed by this")
    p.add_argument("--corot", action="store_true",
                   help="use corotTruss instead of the small-displacement Truss. VK3 is a squat "
                        "pier carrying 1300 kN to ~1%% drift, so P-delta and strut-axis rotation "
                        "are first-order (D22/D60)")
    p.add_argument("--integrator", choices=("newmark", "hht"), default="newmark",
                   help="'newmark' (0.5, 0.25) is average-acceleration: energy-conserving, so it "
                        "sustains spurious high-frequency modes. 'hht' (alpha 0.7) adds numerical "
                        "damping that annihilates them — and is what run_cyclic_dynamic ALREADY "
                        "uses, so a screen with newmark certifies nothing about the cyclic run")
    p.add_argument("--steps-per-period", type=int, default=30,
                   help="dynamic solver: integration steps per fundamental period, i.e. dt = "
                        "T1/N (default 30). T1 is a GLOBAL period; a strut that cracks creates a "
                        "much faster LOCAL mode, and once that mode is finer than dt the "
                        "integrator stops resolving it. Raising this is the direct test")
    p.add_argument("--capture", action="store_true",
                   help="snapshot the full nodal displacement field at peak shear and at the last "
                        "step, and draw the strut-level damage map (D53). The direct way to SEE "
                        "which struts have failed and where, rather than infer it")
    p.add_argument("--steel-b", type=float, default=None,
                   help="DIAGNOSTIC ONLY: override the measured hardening ratio b (0.00466). "
                        "Tests whether the bar chain is localizing. Never a production setting")
    p.add_argument("--draw", dest="draw_model", action="store_true",
                   help="also save the analysis-model lattice figure before the run starts")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         gf_factor=a.gf_factor, solver=a.solver, periods=a.periods, damping=a.damping,
         rate=a.rate, quasi_static=a.quasi_static, gauge_every=a.gauge_every,
         draw_model=a.draw_model, steel_r0=a.steel_r0, steel_b=a.steel_b,
         capture=a.capture, steps_per_period=a.steps_per_period,
         integrator=a.integrator, corot=a.corot)
