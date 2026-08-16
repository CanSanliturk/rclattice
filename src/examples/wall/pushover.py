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

Output: examples/output/wall/wall_pushover[_elastic].png. Units: N, mm.
Run as `python examples/wall/pushover.py [--compression crushing|elastic] [--drift 0.02]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import run_pushover, run_pushover_dynamic

from build import calibrate, nonlinear_wall_lattice
from specimen import A_SHEAR, EPS, HORIZON, LW, MESH, OUT, base_nodes, control_node, lateral_loads

# Measured response of SW-NC-FF, for annotation only — never a calibration target.
PAPER_PEAK_KN = 212.6
PAPER_PEAK_DRIFT = 0.0100
PAPER_YIELD_DRIFT = 0.0030
PAPER_ACI_KN = 199.8


def main(*, compression: str = "crushing", drift: float = 0.02, mesh_size: float = MESH,
         horizon: float = HORIZON, gf_factor: float = 1.0, solver: str = "static",
         periods: float = 12.0, damping: float = 0.8, rate: float | None = None) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    print(f"Aydin energy balance: A_t = {cal.area:,.1f} mm^2 = {cal.area / (mesh_size * 200.0):.4f}"
          f" * (thickness * mesh)   [horizon {horizon}, lattice nu = {cal.nu_consistent:.3f}]")
    if gf_factor != 1.0:
        print(f"tension stiffening: Gf scaled x{gf_factor:g}")

    model = nonlinear_wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon,
                                   compression=compression, gf_factor=gf_factor)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    print(f"nonlinear lattice (compression={compression}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {len(model.elements) - struts} rebar struts")

    target = drift * A_SHEAR
    ctrl, base = control_node(model), base_nodes(model)
    if solver == "dynamic":
        # Dynamic relaxation (D22/D46): the whole load-level row is driven by an imposed ramp, slow
        # enough that inertia is negligible, so the recorded curve is quasi-static — but the mass
        # and damping regularize the softening instability that stops the static solver dead.
        drive = select_nodes(model, (-LW, LW, A_SHEAR - EPS, A_SHEAR + EPS))
        res = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=target,
                                   drive_nodes=drive, base_nodes=base,
                                   periods_to_target=periods, rate=rate,
                                   steps_per_period=30, damping_ratio=damping)
    else:
        res = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                           control_dof=1, dU=target / 400.0, target=target, base_nodes=base,
                           algorithm=("ModifiedNewton", "-initial"))

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
    p.add_argument("--periods", type=float, default=12.0,
                   help="dynamic solver: fundamental periods taken to reach the target. Higher = "
                        "slower ramp = less inertial contamination of the recorded base shear")
    p.add_argument("--damping", type=float, default=0.8, help="dynamic solver: damping ratio")
    p.add_argument("--rate", type=float, default=None,
                   help="dynamic solver: drive speed in mm/s, set DIRECTLY (preferred over "
                        "--periods, which scales with the target). 7.6 mm/s matched the static "
                        "solver within 1%%")
    a = p.parse_args()
    main(compression=a.compression, drift=a.drift, mesh_size=a.mesh, horizon=a.horizon,
         gf_factor=a.gf_factor, solver=a.solver, periods=a.periods, damping=a.damping,
         rate=a.rate)
