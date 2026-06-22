"""Single RC cantilever column, LINEAR-ELASTIC: lattice vs fiber beam-column (linear sibling of
pushover.py).

Same specimen, geometry, reinforcement topology, lattice graph (HORIZON) and stiffness calibration
as pushover.py — only the constitutive law changes: every material is linear `Elastic`.

  - reference: the same fiber `forceBeamColumn` cantilever, but its core/cover/steel fibers are
    Elastic (concrete E, steel E0) — the EXACT same section, just linear;
  - lattice: the SAME RC lattice topology (concrete struts + longitudinal bars + stirrups) as the
    nonlinear case, but with Elastic concrete and Elastic rebar.

With linear materials the pushover of each model is a straight line whose slope is its lateral
stiffness, so this script is really a stiffness-calibration check: the lattice concrete strut area
is calibrated so the FULL elastic lattice (concrete + rebar, which adds ~flexural stiffness) has the
same initial lateral stiffness K0 as the elastic fiber column — after which the two lines overlap.
That same matched K0 (hence matched period) is what makes the linear DYNAMIC comparison
(dynamic_linear.py) agree. Units: kip, in. Output: examples/output/column/column_pushover_linear.png.
Run as `python examples/column/pushover_linear.py`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import beamcolumn_reference_linear, calibrate_area_linear, lattice_k0, rc_lattice_linear
from specimen import DU, OUT, TARGET, lateral_loads


def main(*, draw: bool = False) -> None:
    outdir = OUT / "pushover_linear" / "beamcolumn"
    outdir.mkdir(parents=True, exist_ok=True)

    bc = beamcolumn_reference_linear()
    k_bc = bc["shear"][1] / bc["disp"][1]
    print(f"elastic fiber beam-column: K0={k_bc:.2f} kip/in | V@{bc['disp'][-1]:.1f}in={bc['shear'][-1]:.1f} kip "
          f"(conv={bc['converged']})")

    area, ctrl, base = calibrate_area_linear(k_bc)
    model = rc_lattice_linear(area)
    k_lat = lattice_k0(area, ctrl, base)
    print(f"calibrated concrete strut area A = {area:.3f} in^2 -> lattice K0={k_lat:.2f} kip/in "
          f"({(k_lat / k_bc - 1) * 100:+.2f}% vs beam-column)")
    print(f"elastic RC lattice: {len(model.nodes)} nodes, {len(model.elements)} struts")

    lat = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)
    print(f"lattice: V@{lat['disp'][-1]:.1f}in={lat['shear'][-1]:.1f} kip (conv={lat['converged']})")

    viz.figure_pushover(
        [
            {"disp": bc["disp"], "shear": bc["shear"], "label": "elastic fiber beam-column",
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat["disp"], "shear": lat["shear"], "label": "elastic RC lattice",
             "style": {"color": "C0", "ls": "--", "lw": 2, "marker": "."}},
        ],
        savepath=str(outdir / "column_pushover_linear.png"),
        xlabel="tip displacement (in)", ylabel="base shear (kip)",
        title="LINEAR RC cantilever column pushover: lattice vs fiber beam-column",
    )
    print(f"saved pushover curve to {outdir / 'column_pushover_linear.png'}")

    if draw:
        drawpath = outdir / "column_pushover_linear_model.png"
        viz.figure_model([("RC lattice (linear)", model)], savepath=str(drawpath),
                         suptitle="LINEAR RC cantilever column — analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="LINEAR RC cantilever column pushover")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model (rebar highlighted)")
    args = parser.parse_args()
    main(draw=args.draw)
