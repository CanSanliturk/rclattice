"""VK3 wall-type bridge pier as a lattice — ELASTIC, calibrated by Aydin's energy balance.

Builds the pier-on-foundation lattice and calibrates the one uniform strut area by Aydin's (2017)
overlapping-lattice energy balance: equate the continuum's stored elastic energy under an affine
strain field to the lattice's with EA = 1, and divide. No reference model, no FE solve, no optimiser.

WHAT THE ELASTIC TARGET IS, AND WHAT IT IS NOT. The chapter's k0 = 60 kN/mm (Tab. 5.10) is a SECANT
TO FIRST YIELD of an already heavily cracked pier — about 0.24x the uncracked transformed-section
cantilever — so comparing an uncracked elastic lattice against it would be a category error. Two
better targets exist:

  * the closed-form transformed-section cantilever, the like-for-like check both other wall studies
    use (they land at K_lattice/K_transformed = 0.93);
  * the MEASURED initial stiffness, ~174 kN/mm, recovered by `digitize.py` from the first elastic
    load level (+-0.86 mm at ~151 kN). That level sits below the pier's own cracking shear of
    ~208 kN, so its secant genuinely is an elastic stiffness — a real measured target that the
    chapter itself never prints. This is the payoff from digitizing the loops.

THE SHEAR SHARE IS A RESULT HERE, not a footnote. At Lv/lw = 2.20 shear is ~14% of the ELASTIC
flexibility before any cracking, and the chapter measures 20-22% of the top displacement through the
inelastic range (Fig. 5.19). Since the energy balance matches the lattice's normal stiffness exactly
but never its shear stiffness independently, the elastic shear share is the first check on whether
the calibration is fit for a shear-relevant specimen.

Output: examples/output/vk3_wall/vk3_elastic.png. Units: N, mm.
Run from src/: python examples/vk3_wall/elastic.py [--mesh 50] [--horizon 1.5] [--draw]
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import (
    calibrate, cantilever_stiffness, gross_inertia, pier_lattice, report_calibration,
    transformed_inertia,
)
from specimen import (
    A_SHEAR, EC, FT, HORIZON, LW, MEASURED_LOOPS, MESH, N_BASE, OUT, S_LONG, SPECIMEN, TW,
    base_nodes, control_node, lateral_loads, run_dir,
)
import testdata as td

TARGET_DRIFT = 0.0002   # 0.02% — below the ~0.026% drift at which the pier's own f_t is reached


def measured_initial_stiffness():
    """Uncracked stiffness measured at the first elastic load level, or None if not digitized."""
    if not MEASURED_LOOPS.exists():
        return None
    import numpy as np
    d = np.load(MEASURED_LOOPS)
    k = float(d["k_initial_kN_per_mm"]) if "k_initial_kN_per_mm" in d else None
    return k * 1e3 if k else None


def cracking_shear() -> float:
    """V_cr from the model's own f_t and the axial load — the drift beyond which elastic stops."""
    m_cr = (FT + N_BASE / (LW * TW)) * TW * LW ** 2 / 6.0
    return m_cr / A_SHEAR


def report_energy_balance(cal, mesh_size: float, horizon: float) -> None:
    nominal = TW * mesh_size
    print("Aydin (2017) OLM energy balance  [Eqs 2.1-2.3, plane stress]")
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    print(f"  lattice patch                {cal.n_nodes} nodes, {cal.n_struts} struts")
    print(f"  physical tributary area      {nominal:,.1f} mm^2  ->  lattice would be "
          f"{nominal / cal.area:.3f}x too stiff")
    print(f"  grid anisotropy |Ax-Ay|/Ax   {cal.anisotropy * 100:.2f}%")
    print(f"  nu_consistent                {cal.nu_consistent:.4f}   (where the normal and shear "
          f"balances agree — NOT a Poisson ratio)")
    print(f"  lattice's own Poisson ratio  {cal.nu_effective:.4f}   (cubic anisotropy "
          f"{cal.cubic_anisotropy:.3f}; 1.0 would be isotropic)")
    print(f"  using nu = {cal.nu:.2f}  ->  shear-stiffness error {cal.isotropy_error * 100:.2f}%")


def foundation_sensitivity(mesh_size: float, horizon: float) -> None:
    """Price the assumed foundation thickness by rebuilding at 700 / 1000 / 1200 mm.

    `specimen.FND_W` is an assumption — the chapter never gives the block's out-of-plane thickness —
    so the honest thing is to measure what it costs rather than to defend the number.
    """
    import importlib
    import specimen as sp
    original = sp.FND_W
    print("\nfoundation-thickness sensitivity (the block's out-of-plane width is ASSUMED):")
    for w in (700.0, 1000.0, 1200.0):
        sp.FND_W = w
        sp.FOUNDATION = sp.ConcreteGrade("foundation", E=sp.EC * (w / sp.TW), nu=sp.NU,
                                         rho=sp.RHO_C, fc=sp.FC, epsc0=sp.EPSC0,
                                         fcu=0.2 * sp.FC, epsU=0.010, ft=sp.FT)
        sp.GRADES["foundation"] = sp.FOUNDATION
        import build as bd
        importlib.reload(bd)
        cal = bd.calibrate(mesh_size=mesh_size, horizon=horizon)
        model = bd.pier_lattice(cal.area, mesh_size=mesh_size, horizon=horizon)
        res = run_pushover(model, lateral_loads=sp.lateral_loads(model),
                           control_node=sp.control_node(model), control_dof=1,
                           dU=TARGET_DRIFT * A_SHEAR / 20.0, target=TARGET_DRIFT * A_SHEAR,
                           base_nodes=sp.base_nodes(model))
        k = res["shear"][-1] / res["disp"][-1]
        print(f"  FND_W = {w:6.0f} mm  ->  K_lattice = {k / 1e3:7.2f} kN/mm")
    sp.FND_W = original


def main(*, mesh_size: float = MESH, horizon: float = HORIZON, draw: bool = False,
         sensitivity: bool = False) -> None:
    out = run_dir("elastic")
    print(f"output directory: {out}")

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_energy_balance(cal, mesh_size, horizon)

    model = pier_lattice(cal.area, mesh_size=mesh_size, horizon=horizon)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    bars = len(model.elements) - struts
    print(f"\n{SPECIMEN} lattice (mesh={mesh_size:.0f} mm, horizon={horizon}): "
          f"{len(model.nodes)} nodes, {struts} concrete struts + {bars} rebar struts")

    if draw:
        import draw as draw_mod
        draw_mod.main(mesh_size=mesh_size, horizon=horizon, savepath=out / "model.png")

    ctrl, base = control_node(model), base_nodes(model)
    target = TARGET_DRIFT * A_SHEAR
    res = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=target / 40.0, target=target, base_nodes=base)
    if not res["converged"]:
        raise RuntimeError("elastic pushover did not reach the target drift")

    k_lat = res["shear"][-1] / res["disp"][-1]
    i_gross, i_tr = gross_inertia(), transformed_inertia(mesh_size)
    k_gross, _ = cantilever_stiffness(shear_span=A_SHEAR, inertia=i_gross)
    k_tr, shear_share = cantilever_stiffness(shear_span=A_SHEAR, inertia=i_tr)
    secants = [s / u for u, s in zip(res["disp"][1:], res["shear"][1:])]
    print(f"\nelastic response to {TARGET_DRIFT * 100:.2f}% drift ({target:.2f} mm at "
          f"y = {A_SHEAR:.0f} mm)")
    print(f"  lattice secant stiffness     {k_lat / 1e3:7.2f} kN/mm")
    print(f"  cantilever, plain concrete   {k_gross / 1e3:7.2f} kN/mm  (gross section)")
    print(f"  cantilever, transformed      {k_tr / 1e3:7.2f} kN/mm  (+ vertical bars at n = "
          f"Es/Ec = {S_LONG.E0 / EC:.1f}; I_tr/I_g = {i_tr / i_gross:.3f})")
    print(f"  ratio  K_lattice / K_transf. {k_lat / k_tr:7.4f}   "
          f"(< 1 expected: the foundation adds base flexibility the fixed-base hand calc omits)")
    print(f"  shear share of flexibility   {shear_share * 100:7.1f}%   at Lv/lw = "
          f"{A_SHEAR / LW:.2f}; the test measured {td.COMPONENTS['shear_pct_range'][0]:.0f}-"
          f"{td.COMPONENTS['shear_pct_range'][1]:.0f}% of DISPLACEMENT once cracked (Fig. 5.19)")
    print(f"  secant spread over the pull  {(max(secants) - min(secants)) / k_lat * 100:7.4f}%  "
          f"({'LINEAR' if (max(secants) - min(secants)) / k_lat < 1e-3 else 'nonlinear'})")

    k_meas = measured_initial_stiffness()
    if k_meas:
        print(f"\n  MEASURED initial stiffness   {k_meas / 1e3:7.2f} kN/mm  (digitized Fig. 5.13, "
              f"first elastic level)")
        print(f"  ratio  K_lattice / K_measured{k_lat / k_meas:7.4f}")
    v_cr = cracking_shear()
    print(f"\n  for context — the model's own f_t = {FT:.2f} MPa with N = {N_BASE / 1e3:.0f} kN "
          f"gives V_cr = {v_cr / 1e3:.0f} kN,")
    print(f"  reached at {v_cr / k_tr:.2f} mm = {v_cr / k_tr / A_SHEAR:.3%} drift;")
    print(f"  the chapter's k0 = {td.MEASURED['k0_MN_per_m']:.0f} kN/mm is a secant to FIRST YIELD "
          f"({td.MEASURED['k0_MN_per_m'] * 1e3 / k_tr:.2f}x the uncracked value), i.e. the pier is "
          f"heavily cracked by yield")

    drift = [u / A_SHEAR * 100.0 for u in res["disp"]]
    curves = [
        {"disp": drift, "shear": [s / 1e3 for s in res["shear"]],
         "label": f"lattice, Aydin-calibrated A$_t$={cal.area:,.0f} mm² ({k_lat / 1e3:.1f} kN/mm)",
         "style": {"color": "C0", "lw": 2, "marker": ".", "markevery": 5}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, k_tr * res["disp"][-1] / 1e3],
         "label": f"transformed-section cantilever ({k_tr / 1e3:.1f} kN/mm)",
         "style": {"color": "C3", "ls": "--", "lw": 2}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, k_gross * res["disp"][-1] / 1e3],
         "label": f"plain-concrete cantilever ({k_gross / 1e3:.1f} kN/mm)",
         "style": {"color": "0.6", "ls": ":", "lw": 1.5}},
    ]
    if k_meas:
        curves.append({"disp": [0.0, drift[-1]], "shear": [0.0, k_meas * res["disp"][-1] / 1e3],
                       "label": f"measured initial stiffness ({k_meas / 1e3:.1f} kN/mm)",
                       "style": {"color": "C2", "ls": "-.", "lw": 2}})
    savepath = out / "elastic.png"
    viz.figure_pushover(curves, savepath=str(savepath), xlabel="drift ratio (%)",
                        ylabel="lateral load (kN)",
                        title=f"{SPECIMEN} bridge pier — elastic lattice, Aydin energy-balance "
                              f"calibration")
    print(f"\nsaved elastic response to {savepath}")

    if sensitivity:
        foundation_sensitivity(mesh_size, horizon)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="VK3 elastic lattice, Aydin OLM energy balance")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--draw", action="store_true", help="also save the analysis-model drawing")
    p.add_argument("--foundation-sensitivity", dest="sensitivity", action="store_true",
                   help="rebuild at FND_W = 700/1000/1200 mm to price the assumed block thickness")
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon, draw=a.draw, sensitivity=a.sensitivity)
