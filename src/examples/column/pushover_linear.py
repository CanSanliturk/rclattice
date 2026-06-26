"""Single RC cantilever column, LINEAR-ELASTIC: lattice vs a selectable reference (linear sibling
of pushover.py).

Same specimen, geometry, reinforcement topology, lattice graph (HORIZON) and stiffness calibration
as pushover.py — only the constitutive law changes: every material is linear `Elastic`.

  - reference beamcolumn: the same fiber `forceBeamColumn` cantilever, but its core/cover/steel
    fibers are Elastic (concrete E, steel E0) — the EXACT same section, just linear;
  - reference continuum: the 2D plane-stress continuum with ElasticIsotropic+PlaneStress quads
    and elastic rebar — the like-for-like elastic linear analog of the nonlinear continuum (D29);
  - lattice: the SAME RC lattice topology (concrete struts + longitudinal bars + stirrups) as the
    nonlinear case, but with Elastic concrete and Elastic rebar.

With linear materials the pushover of each model is a straight line whose slope is its lateral
stiffness, so this script is a stiffness-calibration check: the lattice concrete strut area is
calibrated so the FULL elastic lattice has the same K0 as the selected reference.
Units: kip, in. Output: examples/output/column/pushover_linear/{reference}/.
Run as `python examples/column/pushover_linear.py [--reference {beamcolumn,continuum}]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import (
    calibrate_area_linear, lattice_k0, make_reference_linear,
    modal_calibration_figure, rc_lattice_linear,
)
from specimen import DU, OUT, TARGET, lateral_loads

REFERENCE_LABEL = {"beamcolumn": "elastic fiber beam-column", "continuum": "elastic 2D continuum"}


def main(*, reference: str = "beamcolumn", draw: bool = False) -> None:
    outdir = OUT / "pushover_linear" / reference
    outdir.mkdir(parents=True, exist_ok=True)
    label = REFERENCE_LABEL[reference]

    ref = make_reference_linear(reference)
    k_ref = ref["shear"][1] / ref["disp"][1]
    print(f"{label}: K0={k_ref:.2f} kip/in | V@{ref['disp'][-1]:.1f}in={ref['shear'][-1]:.1f} kip "
          f"(conv={ref['converged']})")

    area, ctrl, base = calibrate_area_linear(k_ref)
    model = rc_lattice_linear(area)
    k_lat = lattice_k0(area, ctrl, base)
    print(f"calibrated concrete strut area A = {area:.3f} in^2 -> lattice K0={k_lat:.2f} kip/in "
          f"({(k_lat / k_ref - 1) * 100:+.2f}% vs {label})")
    print(f"elastic RC lattice: {len(model.nodes)} nodes, {len(model.elements)} struts")

    lat = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)
    print(f"lattice: V@{lat['disp'][-1]:.1f}in={lat['shear'][-1]:.1f} kip (conv={lat['converged']})")

    viz.figure_pushover(
        [
            {"disp": ref["disp"], "shear": ref["shear"], "label": label,
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat["disp"], "shear": lat["shear"], "label": "elastic RC lattice",
             "style": {"color": "C0", "ls": "--", "lw": 2, "marker": "."}},
        ],
        savepath=str(outdir / f"column_pushover_linear_{reference}.png"),
        xlabel="tip displacement (in)", ylabel="base shear (kip)",
        title=f"LINEAR RC cantilever column pushover: lattice vs {label}",
    )
    print(f"saved pushover curve to {outdir / f'column_pushover_linear_{reference}.png'}")

    # calibration output (D35): first N mode shapes (selected reference vs lattice) + periods
    caption = (f"linear scalar calibration — strut area A={area:.3f} in^2 -> K0={k_lat:.2f} kip/in "
               f"({(k_lat / k_ref - 1) * 100:+.2f}% vs {label})")
    modalpath = outdir / f"column_modal_linear_{reference}.png"
    t_ref, t_lat = modal_calibration_figure(reference=reference, lattice_model=model,
                                            label=label, caption=caption,
                                            savepath=str(modalpath), linear=True)
    print(f"modal calibration: {label} T={[f'{t:.4f}' for t in t_ref]} s | "
          f"lattice T={[f'{t:.4f}' for t in t_lat]} s -> saved {modalpath}")

    if draw:
        panels = [("RC lattice (linear)", model)]
        if reference == "continuum":
            from build import _continuum_model_linear
            panels.insert(0, ("elastic 2D continuum", _continuum_model_linear()[0]))
        drawpath = outdir / f"column_pushover_linear_{reference}_model.png"
        viz.figure_model(panels, savepath=str(drawpath),
                         suptitle="LINEAR RC cantilever column — analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="LINEAR RC cantilever column pushover")
    parser.add_argument("--reference", choices=("beamcolumn", "continuum"), default="beamcolumn",
                        help="verification reference: elastic fiber beam-column (fast, default) "
                             "or linear 2D plane-stress continuum (ElasticIsotropic quads)")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model (rebar highlighted)")
    args = parser.parse_args()
    main(reference=args.reference, draw=args.draw)
