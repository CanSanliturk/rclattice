"""Concrete cube, uniaxial tension — PLASTIC pushover: nonlinear Concrete02 lattice vs dispBeamColumn.

Building on the elastic baseline (elastic.py / D42), this traces the SAME calibrated-area cube lattice
PAST the cracking strain into the nonlinear regime: Concrete02 struts (tension cracking + softening),
DisplacementControl tension pull, compared against the Concrete02 dispBeamColumn backbone.

NUMERICS — the primary solution algorithm is ModifiedNewton with the INITIAL (elastic) stiffness
(`ModifiedNewton -initial`) BY DEFAULT: iterating on the constant elastic tangent never re-forms the
singular/negative tangent a cracked truss band produces (the D38 mechanism), so it is the robust
choice for a softening pass. `--full-newton` switches to plain Newton for contrast (it dies earlier).

CRACK-BAND REGULARIZATION (`--regularize`, D20/D45) — plain Concrete02 struts localize the crack into
ONE ~10 mm strut row, dissipating ~10x less fracture energy than the reference's smeared 100 mm element,
so the plain lattice descending branch drops ~10x too steeply. Regularization sets each strut's tension
softening slope from its own length h so the dissipated energy per crack area stays = Gf regardless of
h: `Ets = ft^2 * h / (2*Gf)` (crack opening w = 2*Gf/ft, band-width-independent). With Gf matched to the
reference's IMPLIED fracture energy `REF_GF` (the reference is plain Concrete02, Ets=0.1*E, smeared over
L), the two global descending branches line up. Regularization fixes the localization but NOT the
mechanism: at full crack the band still softens to zero (singular tangent), so a STATIC solver still
cannot reach zero stress — the complete branch needs regularization + dynamic relaxation (D38/D39).

SOLVER (`--solver`, D46) — 'static' DisplacementControl (default) dies at the crack mechanism; 'dynamic'
dynamic relaxation (`run_pushover_dynamic`, D22) rides through the mechanism via inertia and is the only
one that traces the full descending branch to zero. Regularized + dynamic is what makes the lattice
descending branch FOLLOW the reference (D39).

The strut area is the elastic K0 calibration (elastic.py) — matched to the beam-column's E*A/L and NOT
re-tuned for the nonlinear range. Fully fixed base, top edge pulled +Y.
Output: examples/output/cube/plastic/cube_tension_plastic[_reg][_dynamic|_newton].png. Units: N, mm.

Run as `python examples/cube/plastic.py [--regularize] [--gf G] [--solver static|dynamic] [--full-newton] [--draw]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import build_lattice_rc, select_nodes
from rclattice.materials import concrete_uniaxial_nonlinear, concrete_uniaxial_regularized
from rclattice.opensees import run_pushover, run_pushover_dynamic

from elastic import beamcolumn_reference, calibrate_area
from specimen import A_FACE, CONCRETE, DU, EPS, HORIZON, L, MESH, OUT, TARGET, axial_loads, cube_problem

OUT_PLASTIC = OUT / "plastic"

EPS_CR = CONCRETE.ft / CONCRETE.E   # cracking strain ft/E (elastic-branch end); where "plastic" begins

# Reference's IMPLIED fracture energy: the single dispBeamColumn uses plain Concrete02 (default
# Ets = 0.1*E) smeared over its whole length L, so Gf_ref = ft^2 * L / (2*Ets_ref) (= 0.15 N/mm here).
# Regularizing the lattice struts to this SAME Gf makes the localized lattice crack dissipate the same
# energy per crack area as the smeared reference -> aligned global descending branches (D45).
REF_ETS = 0.1 * CONCRETE.E
REF_GF = CONCRETE.ft ** 2 * L / (2.0 * REF_ETS)


def cube_lattice(area: float, *, regularize: bool, Gf: float = REF_GF, horizon: float = HORIZON):
    """The calibrated Concrete02 strut lattice (uniform `area`, no rebar). `regularize=True` gives each
    strut a crack-band length-regularized softening slope `Ets = ft^2*L/(2*Gf)` (D20) at fracture energy
    `Gf`; `regularize=False` is plain Concrete02 (same slope on every strut, so the crack localizes)."""
    def material_for(zone: str, length: float):
        if regularize:
            return concrete_uniaxial_regularized(CONCRETE, 0, length, Gf=Gf, Gfc=250.0 * Gf,
                                                 residual_ratio=0.0)
        return concrete_uniaxial_nonlinear(CONCRETE, 0)

    model, _ = build_lattice_rc(
        cube_problem(), MESH, material_for=material_for,
        zone_of=lambda x, y: "concrete", rebars=(), strut_area=area,
        horizon=horizon, strut_element="Truss",
    )
    return model


def _lattice_pushover(solver: str, model, ctrl: int, base: list, algorithm: tuple) -> dict:
    """Trace the tension pushover with the chosen path-following solver:
      - `static`  DisplacementControl (dies at the crack mechanism; `algorithm` primary, default
                  ModifiedNewton -initial);
      - `dynamic` dynamic relaxation (`run_pushover_dynamic`, D22) — inertia keeps the system
                  well-posed THROUGH the cracked-band mechanism, so it traces the full descending
                  branch to zero. The top edge is driven (imposed ramp) so the pull stays uniform."""
    if solver == "static":
        return run_pushover(model, lateral_loads=axial_loads(model), control_node=ctrl,
                            control_dof=2, dU=DU, target=TARGET, base_nodes=base, algorithm=algorithm)
    if solver == "dynamic":
        top = select_nodes(model, (-EPS, L + EPS, L - EPS, L + EPS))   # drive the whole top edge
        return run_pushover_dynamic(model, control_node=ctrl, control_dof=2, target=TARGET,
                                    drive_nodes=top, base_nodes=base,
                                    periods_to_target=12.0, steps_per_period=30, damping_ratio=0.8)
    raise ValueError(f"unknown solver {solver!r} (expected static/dynamic)")


def main(*, regularize: bool = False, Gf: float = REF_GF, solver: str = "static",
         algorithm: tuple = ("ModifiedNewton", "-initial"), draw: bool = False) -> None:
    outdir = OUT_PLASTIC
    outdir.mkdir(parents=True, exist_ok=True)

    # reference: Concrete02 dispBeamColumn — one element, constant axial strain, so it traces the
    # whole backbone (rise to ft, then softening) to zero with plain DisplacementControl (D37).
    ref = beamcolumn_reference("concrete02", TARGET)
    k_ref = ref["shear"][1] / ref["disp"][1]                       # E*A/L (initial tangent), N/mm
    print(f"dispBeamColumn (Concrete02): peak N = {max(ref['shear']):.1f} N "
          f"({max(ref['shear']) / A_FACE:.2f} MPa = ft) | u->{ref['disp'][-1]:.4f} mm (conv={ref['converged']})")

    area, ctrl, base = calibrate_area(k_ref)
    law = f"crack-band regularized (Gf={Gf:.3f} N/mm)" if regularize else "plain Concrete02"
    print(f"scalar K0 calibration (from the elastic baseline): strut area A = {area:.2f} mm^2 | {law}")

    model = cube_lattice(area, regularize=regularize, Gf=Gf)
    algo_lbl = "dynamic relaxation" if solver == "dynamic" else " ".join(algorithm)
    print(f"cube lattice: {len(model.nodes)} nodes, {len(model.elements)} Concrete02 struts | "
          f"solver = {solver} ({algo_lbl}) | cracking strain eps_cr = {EPS_CR:.2e}")
    lat = _lattice_pushover(solver, model, ctrl, base, algorithm)
    print(f"lattice: peak N = {max(lat['shear']):.1f} N ({max(lat['shear']) / A_FACE:.2f} MPa) | "
          f"traced to u = {lat['disp'][-1]:.4f} mm (eps = {lat['disp'][-1] / L:.2e}, "
          f"{lat['disp'][-1] / L / EPS_CR:.1f}x eps_cr) | conv = {lat['converged']}")

    # coupon stress-strain: axial stress (MPa) = base reaction / cube-face area, strain = u / edge
    tag = "regularized" if regularize else "plain"
    curves = [
        {"disp": [u / L for u in ref["disp"]], "shear": [s / A_FACE for s in ref["shear"]],
         "label": "dispBeamColumn (Concrete02, full backbone)", "style": {"color": "C3", "lw": 2}},
        {"disp": [u / L for u in lat["disp"]], "shear": [s / A_FACE for s in lat["shear"]],
         "label": f"Concrete02 lattice ({tag}, {algo_lbl})",
         "style": {"color": "C0", "ls": "--", "lw": 2, "marker": ".", "markevery": 10}},
    ]
    solver_tag = "_dynamic" if solver == "dynamic" else ("_newton" if algorithm[0] == "Newton" else "")
    stem = "cube_tension_plastic" + ("_reg" if regularize else "") + solver_tag
    savepath = outdir / f"{stem}.png"
    viz.figure_pushover(curves, savepath=str(savepath),
                        xlabel="axial strain (-)", ylabel="axial stress (MPa)",
                        title=f"Concrete cube, uniaxial tension — plastic pushover ({tag}, {algo_lbl}): "
                              f"Concrete02 lattice vs single dispBeamColumn")
    print(f"saved stress-strain curve to {savepath}")

    if draw:
        drawpath = outdir / f"{stem}_model.png"
        viz.figure_model([("cube lattice (Concrete02)", model)], savepath=str(drawpath),
                         suptitle="Concrete cube — Concrete02 lattice analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Concrete cube uniaxial tension, PLASTIC pushover: "
                                                 "Concrete02 lattice vs single dispBeamColumn")
    parser.add_argument("--regularize", action="store_true",
                        help="crack-band (fracture-energy) regularize each strut's softening "
                             "(Ets = ft^2*L/(2*Gf), D20/D45) so the localized lattice crack dissipates "
                             "the reference's Gf and the descending branches align; default off (plain "
                             "Concrete02, which localizes and drops ~10x too steeply)")
    parser.add_argument("--gf", type=float, default=REF_GF,
                        help=f"tension fracture energy Gf (N/mm) for --regularize; default {REF_GF:.3f} "
                             f"= the reference dispBeamColumn's implied Gf (ft^2*L/(2*0.1*E))")
    parser.add_argument("--solver", choices=("static", "dynamic"), default="static",
                        help="path-following solver: 'static' DisplacementControl (default; dies at the "
                             "crack mechanism) or 'dynamic' dynamic relaxation (rides through the "
                             "mechanism via inertia — the only one that traces the full descending "
                             "branch to zero; regularized + dynamic makes the lattice FOLLOW the "
                             "reference, D38/D39)")
    parser.add_argument("--full-newton", action="store_true",
                        help="(static solver only) use plain Newton instead of the default "
                             "ModifiedNewton -initial — for contrast; it re-forms the singular cracked "
                             "tangent and dies earlier")
    parser.add_argument("--draw", action="store_true", help="also save a drawing of the lattice model")
    args = parser.parse_args()
    algo = ("Newton",) if args.full_newton else ("ModifiedNewton", "-initial")
    main(regularize=args.regularize, Gf=args.gf, solver=args.solver, algorithm=algo, draw=args.draw)
