"""Concrete cube, uniaxial tension: lattice vs a single displacement-based beam-column.

The simplest lattice verification: model ONE 0.1 m concrete cube two ways and compare their
uniaxial-tension responses —
  - reference: a single displacement-based fiber `dispBeamColumn` (constant axial strain -> the
    clean 1D material coupon, Concrete02);
  - lattice:   the same specimen as a Concrete02 strut lattice (gmsh nodes + horizon struts).

The lattice strut area is calibrated so its initial AXIAL stiffness equals the beam-column's
(E*A/L), so both start from the same elastic slope. Fully fixed base, top edge pulled up (+Y).
Output is the coupon stress-strain curve (axial stress = base reaction / cube-face area, strain =
top disp / edge). Units: N, mm.

FOLLOWING THE REFERENCE (two things are needed past the elastic branch, D39):
  1. Concrete law (`--plain` toggles): DEFAULT = fracture-energy length-regularized Concrete02
     struts (crack-band, D20) matched to the reference's implied Gf. Without this the crack
     localizes into ONE strut row (~10 mm) vs the reference's smeared 100 mm, so the lattice drops
     ~10x too steeply and does NOT follow. Regularization makes both dissipate the same Gf, so the
     descending branches span the same strain range. `--plain` recovers the D37 K0-only baseline.
     The PEAK still runs higher than ft — the diagonal struts add a parallel tensile path, a known
     lattice overstrength the K0 calibration does not remove. `--peak-correction` optionally scales it
     onto the reference: `none` (default, leave it), `area` (reduce strut area — matches peak but
     lowers K0 so the elastic branch no longer overlaps, D40), or `strength` (reduce concrete ft —
     matches peak while keeping the area, so K0 stays matched, D41).
  2. Solver (`--solver`): how the path is traced past the peak.
     - `static`    DisplacementControl. Dies at the crack (the softened band -> 0 stiffness -> the
                   block above is a kinematic mechanism).
     - `arclength` cylindrical arc-length (`run_pushover_arclength`). The proper snapback solver, but
                   it CANNOT cross this crack either: the tangent is SINGULAR (a mechanism), not a
                   snapback, and no static path-follower passes a singular tangent (stalls at peak).
     - `dynamic`   dynamic-relaxation (`run_pushover_dynamic`, D22, DEFAULT): a quasi-static transient
                   whose inertia keeps the system well-posed THROUGH the mechanism, so it traces the
                   full pre- and post-peak branch. The only option that follows the branch to zero.
     See D38.

The reference (single dispBeamColumn) always traces its full backbone: one element has a CONSTANT
axial strain field, so it softens uniformly (a limit point, no localization) and DisplacementControl
follows it to zero stress. Output: examples/output/cube/pushover/cube_tension_<solver>[_plain].png.

Run as `python examples/cube/pushover.py [--solver static|arclength|dynamic] [--plain]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.opensees import run_pushover, run_pushover_arclength, run_pushover_dynamic

from build import PEAK_RATIO, REF_GF, beamcolumn_reference, calibrate_area, cube_lattice
from specimen import A_FACE, CONCRETE, DU, EPS, L, OUT, TARGET, axial_loads


def _lattice_pushover(solver: str, model, ctrl: int, base: list[int]) -> dict:
    """Run the calibrated lattice tension pushover with the chosen path-following solver."""
    if solver == "static":
        return run_pushover(model, lateral_loads=axial_loads(model), control_node=ctrl,
                            control_dof=2, dU=DU, target=TARGET, base_nodes=base)
    if solver == "arclength":
        return run_pushover_arclength(model, lateral_loads=axial_loads(model), control_node=ctrl,
                                      control_dof=2, dU=DU, target=TARGET, base_nodes=base)
    if solver == "dynamic":
        top = select_nodes(model, (-EPS, L + EPS, L - EPS, L + EPS))   # drive the whole top edge
        return run_pushover_dynamic(model, control_node=ctrl, control_dof=2, target=TARGET,
                                    drive_nodes=top, base_nodes=base,
                                    periods_to_target=12.0, steps_per_period=30, damping_ratio=0.8)
    raise ValueError(f"unknown solver {solver!r} (expected static/arclength/dynamic)")


def main(*, solver: str = "dynamic", regularize: bool = True, peak_correction: str = "none",
         draw: bool = False) -> None:
    outdir = OUT / "pushover"
    outdir.mkdir(parents=True, exist_ok=True)

    ref = beamcolumn_reference()
    k_ref = ref["shear"][1] / ref["disp"][1]           # initial (pre-cracking) axial stiffness, N/mm
    print(f"dispBeamColumn: K0={k_ref:.1f} N/mm | peak N={max(ref['shear']):.1f} N "
          f"({max(ref['shear']) / A_FACE:.2f} MPa) | u->{ref['disp'][-1]:.4f} mm (conv={ref['converged']})")

    area, ctrl, base = calibrate_area(k_ref)
    print(f"scalar K0 calibration: strut area A = {area:.4f} mm^2 (matches dispBeamColumn K0)")
    if peak_correction == "area":
        print(f"peak correction [area]: strut area /= PEAK_RATIO({PEAK_RATIO}) -> {area / PEAK_RATIO:.4f} mm^2 "
              f"(matches peak; lowers K0 by the same factor)")
    elif peak_correction == "strength":
        print(f"peak correction [strength]: concrete ft /= PEAK_RATIO({PEAK_RATIO}) -> "
              f"{CONCRETE.ft / PEAK_RATIO:.3f} MPa (matches peak; area / K0 unchanged)")
    else:
        print(f"peak correction: none (K0 matched; lattice peak overshoots the reference by ~{PEAK_RATIO})")

    law = f"length-regularized (Gf={REF_GF:.3f} N/mm, matched to reference)" if regularize else "plain Concrete02"
    model = cube_lattice(area, regularize=regularize, peak_correction=peak_correction)
    print(f"cube lattice: {len(model.nodes)} nodes, {len(model.elements)} struts | {law} | solver={solver}")
    lat = _lattice_pushover(solver, model, ctrl, base)
    print(f"lattice: peak N={max(lat['shear']):.1f} N ({max(lat['shear']) / A_FACE:.2f} MPa) | "
          f"u->{lat['disp'][-1]:.4f} mm (conv={lat['converged']})")

    # coupon stress-strain: axial stress (MPa) = base reaction / cube-face area, strain = u / edge
    tag = "regularized" if regularize else "plain"
    pk = "" if peak_correction == "none" else f", peak-fix:{peak_correction}"
    curves = [
        {"disp": [u / L for u in ref["disp"]], "shear": [s / A_FACE for s in ref["shear"]],
         "label": "dispBeamColumn (single element)", "style": {"color": "C3", "lw": 2}},
        {"disp": [u / L for u in lat["disp"]], "shear": [s / A_FACE for s in lat["shear"]],
         "label": f"Concrete02 lattice ({tag}, {solver}{pk})",
         "style": {"color": "C0", "ls": "--", "lw": 2, "marker": "."}},
    ]
    stem = (f"cube_tension_{solver}" + ("" if regularize else "_plain")
            + ("" if peak_correction == "none" else f"_{peak_correction}"))
    savepath = outdir / f"{stem}.png"
    viz.figure_pushover(curves, savepath=str(savepath),
                        xlabel="axial strain (-)", ylabel="axial stress (MPa)",
                        title=f"Concrete cube, uniaxial tension: lattice ({tag}, {solver}{pk}) vs single dispBeamColumn")
    print(f"saved stress-strain curve to {savepath}")

    if draw:
        drawpath = outdir / "cube_tension_model.png"
        viz.figure_model([("cube lattice", model)], savepath=str(drawpath),
                         suptitle="Concrete cube — lattice analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Concrete cube uniaxial tension: lattice vs single dispBeamColumn")
    parser.add_argument("--solver", choices=("static", "arclength", "dynamic"), default="dynamic",
                        help="path-following solver for the lattice pushover: DisplacementControl "
                             "(static) / cylindrical arc-length / dynamic relaxation (default — the "
                             "only one that traces the full descending branch; static & arclength die "
                             "at the crack mechanism, D38)")
    parser.add_argument("--plain", action="store_true",
                        help="plain Concrete02 struts (D37 K0-only baseline) instead of the default "
                             "fracture-energy-regularized struts (regularization makes the lattice "
                             "descending branch FOLLOW the reference; D39)")
    parser.add_argument("--peak-correction", choices=("none", "area", "strength"), default="none",
                        help="optional strength calibration to scale the lattice tensile PEAK onto the "
                             "reference (it overshoots by ~PEAK_RATIO from diagonal overstrength): "
                             "'none' (default, K0 matched / peak overshoots), 'area' (reduce strut "
                             "area — matches peak, lowers K0; D40), 'strength' (reduce concrete ft — "
                             "matches peak, keeps K0; D41)")
    parser.add_argument("--draw", action="store_true", help="also save a drawing of the lattice model")
    args = parser.parse_args()
    main(solver=args.solver, regularize=not args.plain, peak_correction=args.peak_correction, draw=args.draw)
