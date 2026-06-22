"""RC frame pushover — Stage 1 (ELASTIC pipeline de-risk), D18/D19.

Reproduces the geometry/loading of the OpenSees RCFrameGravity -> RCFramePushOver benchmark
(1 bay x 1 storey, units kip-in) with the lattice builder, but ELASTICALLY: the goal of this
stage is to prove the analysis machinery (gravity LoadControl, DisplacementControl pushover,
control node, base-shear-from-reactions, the curve plumbing), not the nonlinear physics.
Stage 2 (pushover_rc.py) adds the actual RC behaviour.

An elastic pushover is a straight line; its slope is the frame's elastic lateral stiffness. We
overlay the lattice against the continuum (same node grid) as a sanity check, after a one-shot
strut-area calibration so their initial stiffnesses agree (K ~ strut area).

Specimen imported from this package (specimen.py). Output: examples/output/frame/frame_pushover.png.
Run as `python examples/frame/pushover.py`. Units: kip, in.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import build_continuum, build_lattice, select_nodes
from rclattice.opensees import run_pushover

from specimen import COL, DU, EPS, FRAME_ELASTIC, H, OUT, SPAN, TARGET, frame_problem, lateral_loads

# Stage 1 is elastic with no rebar, so the mesh need not hit the bar lines (x=+-10.5); a coarse
# mesh keeps the demo fast. Stage 2 requires mesh = 1.5 so nodes land on the rebar.
MESH = 6.0


def main(*, draw: bool = False) -> None:
    outdir = OUT / "pushover" / "continuum"
    outdir.mkdir(parents=True, exist_ok=True)
    problem = frame_problem(FRAME_ELASTIC)

    # lattice and continuum share the identical node grid -> same ids for control/base/loads.
    lattice, _ = build_lattice(problem, MESH, strut_area=1.0)
    continuum, _ = build_continuum(problem, MESH)
    control = select_nodes(lattice, (-EPS, EPS, H - EPS, H + EPS))[0]  # top-left joint (~node 3)
    base = select_nodes(lattice, (-COL, SPAN + COL, -EPS, EPS))
    lat_loads = lateral_loads(lattice)
    print(f"nodes = {len(lattice.nodes)} | control node = {control} | base nodes = {len(base)}")

    kw = dict(lateral_loads=lat_loads, control_node=control, control_dof=1,
              dU=DU, target=TARGET, base_nodes=base)

    # continuum reference pushover, and a trial unit-area lattice -> calibrate strut area by
    # matching the elastic lateral stiffness (base shear / roof drift at the final step).
    cont_res = run_pushover(continuum, **kw)
    trial_res = run_pushover(lattice, **kw)
    k_cont = cont_res["shear"][-1] / cont_res["disp"][-1]
    k_trial = trial_res["shear"][-1] / trial_res["disp"][-1]
    area = k_cont / k_trial
    print(f"K_continuum = {k_cont:.3f} kip/in | calibrated strut area A = {area:.4e} in^2")

    lattice_cal, _ = build_lattice(problem, MESH, strut_area=area)
    lat_res = run_pushover(lattice_cal, **kw)
    k_lat = lat_res["shear"][-1] / lat_res["disp"][-1]
    print(f"K_lattice (calibrated) = {k_lat:.3f} kip/in  (target {k_cont:.3f})")
    print(f"converged: lattice={lat_res['converged']} continuum={cont_res['converged']}")

    viz.figure_pushover(
        [
            {"disp": cont_res["disp"], "shear": cont_res["shear"], "label": "continuum",
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat_res["disp"], "shear": lat_res["shear"], "label": "lattice (calibrated)",
             "style": {"color": "C0", "ls": "--", "lw": 2}},
        ],
        savepath=str(outdir / "frame_pushover.png"),
        xlabel="roof displacement (in)", ylabel="base shear (kip)",
        title="Elastic pushover (Stage 1): lattice vs continuum",
    )
    print(f"saved pushover curve to {outdir / 'frame_pushover.png'}")

    if draw:
        drawpath = outdir / "frame_pushover_model.png"
        viz.figure_model([("lattice", lattice_cal), ("continuum", continuum)],
                         savepath=str(drawpath), suptitle="Elastic frame (Stage 1) — analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Elastic frame pushover (Stage 1): lattice vs continuum")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis models (lattice + continuum)")
    args = parser.parse_args()
    main(draw=args.draw)
