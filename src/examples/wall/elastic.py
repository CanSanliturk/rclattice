"""Flexure-controlled RC shear wall (SW-NC-FF) as a lattice — ELASTIC, calibrated by Aydin's
energy balance.

Builds the wall-on-pedestal lattice for the Sahinkaya et al. (2025) specimen and calibrates the
one uniform strut area by Aydin's (2017) overlapping-lattice energy balance (OLM Eqs 2.1-2.3):
equate the continuum's stored elastic energy under an affine strain field to the lattice's with
EA = 1, and divide. No reference model, no FE solve, no optimiser — the affine field makes every
strut elongation exact, so the calibration is a closed-form sum over the strut list.

WHY A CALIBRATED AREA IS NEEDED: a horizon lattice puts EVERY strut through a node — orthogonal
plus the diagonals that brace it — in parallel, so a physically-sized tributary strut
(area = MESH * TW) over-counts the material and makes the assembly ~1.5x too stiff. The energy
balance is what fixes the scale, and it fixes it as a MATERIAL property: A_t / (thickness * mesh)
is a constant of the (horizon, nu) lattice, so the same calibration holds at any mesh size and
across all three concrete zones of the specimen.

The run then pushes the calibrated lattice to 0.10% drift — below the 0.30% drift at which the
test first yielded, so the comparison stays in the uncracked elastic range the calibration targets
— and reports its secant stiffness against the gross-section cantilever hand calculation.

Output: examples/output/wall/wall_elastic.png. Units: N, mm.
Run as `python examples/wall/elastic.py [--mesh 50] [--horizon 1.5] [--draw]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import calibrate, cantilever_stiffness, gross_inertia, transformed_inertia, wall_lattice
from specimen import (
    A_SHEAR, HORIZON, MESH, OUT, S14, TEST_UNIT, TW, base_nodes, control_node, lateral_loads,
)

TARGET_DRIFT = 0.0010   # 0.10% — the test first yielded at 0.30%, so this is safely pre-cracking


def report_calibration(cal, mesh_size: float) -> None:
    """Print the energy balance and the isotropy diagnostics it implies."""
    nominal = TW * mesh_size
    print("Aydin (2017) OLM energy balance  [Eqs 2.1-2.3, plane stress]")
    print(f"  lattice patch                {cal.n_nodes} nodes, {cal.n_struts} struts")
    print(f"  calibrated strut area A_t    {cal.area:,.1f} mm^2  =  {cal.area / nominal:.4f} * (thickness * mesh)")
    print(f"  physical tributary area      {nominal:,.1f} mm^2  ->  lattice would be "
          f"{nominal / cal.area:.3f}x too stiff")
    print(f"  grid anisotropy |Ax-Ay|/Ax   {cal.anisotropy * 100:.2f}%   (eps_x vs eps_y balance)")
    print(f"  lattice's own Poisson ratio  {cal.nu_consistent:.4f}   (where the normal and shear "
          f"balances agree)")
    print(f"  using nu = {cal.nu:.2f}  ->  shear-stiffness error {cal.isotropy_error * 100:.2f}% "
          f"(the cost of the pinned nu)")


def main(*, mesh_size: float = MESH, horizon: float = HORIZON, draw: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    report_calibration(cal, mesh_size)

    model = wall_lattice(cal.area, mesh_size=mesh_size, horizon=horizon)
    struts = sum(1 for e in model.elements if e.kind not in ("longitudinal", "stirrup"))
    bars = len(model.elements) - struts
    print(f"\nwall lattice (mesh={mesh_size:.0f} mm, horizon={horizon}): {len(model.nodes)} nodes, "
          f"{struts} concrete struts + {bars} rebar struts")

    ctrl, base = control_node(model), base_nodes(model)
    target = TARGET_DRIFT * A_SHEAR
    res = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=target / 40.0, target=target, base_nodes=base)
    if not res["converged"]:
        raise RuntimeError("elastic pushover did not reach the target drift")

    k_lat = res["shear"][-1] / res["disp"][-1]
    i_gross, i_tr = gross_inertia(), transformed_inertia()
    k_gross, _ = cantilever_stiffness(shear_span=A_SHEAR, inertia=i_gross)
    k_tr, shear_share = cantilever_stiffness(shear_span=A_SHEAR, inertia=i_tr)
    secants = [s / u for u, s in zip(res["disp"][1:], res["shear"][1:])]
    print(f"\nelastic response to {TARGET_DRIFT * 100:.2f}% drift ({target:.2f} mm at y = {A_SHEAR:.0f} mm)")
    print(f"  lattice secant stiffness     {k_lat / 1e3:7.2f} kN/mm")
    print(f"  cantilever, plain concrete   {k_gross / 1e3:7.2f} kN/mm  (gross section — NOT the "
          f"like-for-like check)")
    print(f"  cantilever, transformed      {k_tr / 1e3:7.2f} kN/mm  (+ the vertical bars at n = "
          f"Es/Ec ~ {S14.E0 / TEST_UNIT.E:.1f}; I_tr/I_g = {i_tr / i_gross:.3f})")
    print(f"  ratio  K_lattice / K_transf. {k_lat / k_tr:7.4f}   "
          f"(< 1 expected: the pedestal adds base flexibility the fixed-base hand calc omits)")
    print(f"  shear share of flexibility   {shear_share * 100:7.0f}%   (flexure-dominated, as designed)")
    print(f"  secant spread over the pull  {(max(secants) - min(secants)) / k_lat * 100:7.4f}%  "
          f"({'LINEAR' if (max(secants) - min(secants)) / k_lat < 1e-3 else 'nonlinear'})")

    # Context, not a calibration target: the paper's measured uncracked rigidity is 0.62-0.67
    # Ec*Ig, i.e. already reduced by micro-cracking and plain-bar slip. An elastic gross-section
    # model is expected to sit ABOVE it.
    print(f"\n  for context — paper Table 4 measured 0.62-0.67 Ec*Ig uncracked, i.e. an effective "
          f"{0.645 * k_gross / 1e3:.1f} kN/mm;\n  the gap to an uncracked elastic model is damage "
          f"and plain-bar slip, not a modelling error.")

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
    savepath = OUT / "wall_elastic.png"
    viz.figure_pushover(curves, savepath=str(savepath), xlabel="drift ratio (%)",
                        ylabel="lateral load (kN)",
                        title="SW-NC-FF shear wall — elastic lattice, Aydin energy-balance calibration")
    print(f"\nsaved elastic response to {savepath}")

    if draw:
        drawpath = OUT / "wall_model.png"
        viz.figure_model([("SW-NC-FF lattice", model)], savepath=str(drawpath),
                         suptitle="SW-NC-FF shear wall on pedestal — lattice analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="RC shear wall (SW-NC-FF) elastic lattice, calibrated "
                                            "by Aydin's OLM energy balance")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing in mm (default {MESH:.0f})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon as a multiple of the grid spacing (default {HORIZON})")
    p.add_argument("--draw", action="store_true", help="also save a drawing of the lattice model")
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon, draw=a.draw)
