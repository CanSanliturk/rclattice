"""WSH3 shear wall as a lattice — ELASTIC, calibrated by Aydin's energy balance.

Builds the wall-on-foundation lattice for Test Unit WSH3 of Dazio, Beyer & Bachmann (2009) and
calibrates the one uniform strut area by Aydin's (2017) overlapping-lattice energy balance (OLM
Eqs 2.1-2.3): equate the continuum's stored elastic energy under an affine strain field to the
lattice's with EA = 1, and divide. No reference model, no FE solve, no optimiser — the affine field
makes every strut elongation exact, so the calibration is a closed-form sum over the strut list.

WHY A CALIBRATED AREA IS NEEDED: a horizon lattice puts EVERY strut through a node — orthogonal
plus the diagonals that brace it — in parallel, so a physically-sized tributary strut
(area = MESH * TW) over-counts the material and makes the assembly ~1.5x too stiff. The energy
balance fixes the scale as a MATERIAL property: A_t / (thickness * mesh) is a constant of the
(horizon, nu) lattice, so the same calibration holds at any mesh size and across all three zones.

THE ELASTIC TARGET IS SET BY CRACKING, NOT BY YIELD, and this is where the WSH3 run differs from
the SW-NC-FF one. That wall was pushed to 0.10% drift, comfortably under its 0.30% first yield.
WSH3 yields at 0.25% drift but CRACKS at about 0.025%: Table 5's M_cr = 527 kN.m is only 25% of
M_max, so V_cr ~ 116 kN, which the transformed section reaches in ~1.1 mm. Pushing an elastic model
past that would still be exactly linear — the materials are linear — but the number would no longer
be comparable to anything the wall did. So the target is 0.02% drift, and the check below is against
the closed-form cantilever, not against a measured stiffness.

Output: examples/output/katrin_wall/wsh3_elastic.png. Units: N, mm.
Run as `python examples/katrin_wall/elastic.py [--mesh 50] [--horizon 1.5] [--draw]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import (
    calibrate, cantilever_stiffness, gross_inertia, report_calibration, transformed_inertia,
    wall_lattice,
)
from specimen import (
    A_SHEAR, EC, HORIZON, L_V, MESH, OUT, S12, SPECIMEN, TW, UNIT, base_nodes, control_node,
    lateral_loads,
)

TARGET_DRIFT = 0.0002   # 0.02% — below the ~0.025% drift at which Table 5's M_cr is reached

# Measured reference points, for context only — never a calibration target.
PAPER_V_MAX = UNIT["V_max"] / 1e3          # 454 kN, Table 5
PAPER_M_CR = 527.0                         # kN.m, Table 5
PAPER_DELTA_Y = 15.4                       # mm, Table 4 (3/4-rule)


def report_energy_balance(cal, mesh_size: float, horizon: float) -> None:
    """The full calibration report: `build.report_calibration` (area + EA per zone) plus the
    isotropy diagnostics that only this deep-dive script needs."""
    nominal = TW * mesh_size
    print("Aydin (2017) OLM energy balance  [Eqs 2.1-2.3, plane stress]")
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)
    print(f"  lattice patch                {cal.n_nodes} nodes, {cal.n_struts} struts")
    print(f"  physical tributary area      {nominal:,.1f} mm^2  ->  lattice would be "
          f"{nominal / cal.area:.3f}x too stiff")
    print(f"  grid anisotropy |Ax-Ay|/Ax   {cal.anisotropy * 100:.2f}%   (eps_x vs eps_y balance)")
    # nu_consistent is NOT the lattice's Poisson ratio — it is where the two CALIBRATION ROUTES
    # agree. The lattice's actual lateral response is nu_effective, and it is cubic-symmetric rather
    # than isotropic at any nu (D53), which is why both numbers are printed.
    print(f"  nu_consistent                {cal.nu_consistent:.4f}   (where the normal and shear "
          f"balances agree — not a Poisson ratio)")
    print(f"  lattice's own Poisson ratio  {cal.nu_effective:.4f}   (cubic anisotropy "
          f"{cal.cubic_anisotropy:.3f}; 1.0 would be isotropic)")
    print(f"  using nu = {cal.nu:.2f}  ->  shear-stiffness error {cal.isotropy_error * 100:.2f}% "
          f"(the cost of the pinned nu)")


def main(*, mesh_size: float = MESH, horizon: float = HORIZON, draw: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_energy_balance(cal, mesh_size, horizon)

    model = wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    bars = len(model.elements) - struts
    print(f"\n{SPECIMEN} lattice (mesh={mesh_size:.0f} mm, horizon={horizon}): "
          f"{len(model.nodes)} nodes, {struts} concrete struts + {bars} rebar struts")

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
    print(f"  cantilever, plain concrete   {k_gross / 1e3:7.2f} kN/mm  (gross section — NOT the "
          f"like-for-like check)")
    print(f"  cantilever, transformed      {k_tr / 1e3:7.2f} kN/mm  (+ the vertical bars at n = "
          f"Es/Ec = {S12.E0 / EC:.1f}; I_tr/I_g = {i_tr / i_gross:.3f})")
    print(f"  ratio  K_lattice / K_transf. {k_lat / k_tr:7.4f}   "
          f"(< 1 expected: the foundation adds base flexibility the fixed-base hand calc omits)")
    print(f"  shear share of flexibility   {shear_share * 100:7.0f}%   (flexure-dominated: the "
          f"test measured shear at ~12% of the flexural displacement, Fig. 9b)")
    print(f"  secant spread over the pull  {(max(secants) - min(secants)) / k_lat * 100:7.4f}%  "
          f"({'LINEAR' if (max(secants) - min(secants)) / k_lat < 1e-3 else 'nonlinear'})")

    # Context, not a calibration target. Two milestones bracket where an UNCRACKED elastic model
    # stops meaning anything, and both are the paper's own numbers.
    v_cr = PAPER_M_CR * 1e6 / L_V
    print(f"\n  for context — Table 5's cracking moment {PAPER_M_CR:.0f} kN.m is V_cr = "
          f"{v_cr / 1e3:.0f} kN,")
    print(f"  reached at {v_cr / k_tr:.2f} mm = {v_cr / k_tr / A_SHEAR:.3%} drift on this stiffness;")
    print(f"  the measured yield displacement is {PAPER_DELTA_Y:.1f} mm "
          f"({PAPER_DELTA_Y / L_V:.2%} drift), i.e. an effective stiffness of "
          f"{PAPER_V_MAX * 1e3 / PAPER_DELTA_Y / k_tr:.2f} x the uncracked value —")
    print(f"  the wall is heavily cracked well before yield, so an elastic model is a CALIBRATION "
          f"check, not a response prediction.")

    drift = [u / A_SHEAR * 100.0 for u in res["disp"]]
    shear = [s / 1e3 for s in res["shear"]]
    curves = [
        {"disp": drift, "shear": shear,
         "label": f"lattice, Aydin-calibrated A$_t$={cal.area:,.0f} mm² ({k_lat / 1e3:.1f} kN/mm)",
         "style": {"color": "C0", "lw": 2, "marker": ".", "markevery": 5}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, k_tr * res["disp"][-1] / 1e3],
         "label": f"transformed-section cantilever ({k_tr / 1e3:.1f} kN/mm)",
         "style": {"color": "C3", "ls": "--", "lw": 2}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, k_gross * res["disp"][-1] / 1e3],
         "label": f"plain-concrete cantilever ({k_gross / 1e3:.1f} kN/mm)",
         "style": {"color": "0.6", "ls": ":", "lw": 1.5}},
    ]
    savepath = OUT / "wsh3_elastic.png"
    viz.figure_pushover(curves, savepath=str(savepath), xlabel="drift ratio (%)",
                        ylabel="lateral load (kN)",
                        title=f"{SPECIMEN} shear wall — elastic lattice, Aydin energy-balance "
                              f"calibration")
    print(f"\nsaved elastic response to {savepath}")

    if draw:
        drawpath = OUT / "wsh3_model.png"
        viz.figure_model([(f"{SPECIMEN} lattice", model)], savepath=str(drawpath),
                         suptitle=f"{SPECIMEN} shear wall on foundation — lattice analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="WSH3 elastic lattice, calibrated by Aydin's OLM "
                                            "energy balance")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing in mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon as a multiple of the grid spacing (default {HORIZON})")
    p.add_argument("--draw", action="store_true", help="also save a drawing of the lattice model")
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon, draw=a.draw)
