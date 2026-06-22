"""RC frame pushover — Stage 2 (NONLINEAR RC), D18/D19.

Builds the OpenSees RCFrameGravity -> RCFramePushOver benchmark frame (1 bay x 1 storey, kip-in)
as a reinforced-concrete LATTICE and pushes it over, then overlays the result on the original
fiber `forceBeamColumn` benchmark (run_benchmark_rc_frame) as the verification reference.

Lattice idealization (per the resolved forks, D19):
  - Columns are 24-in-deep (X) x 144-tall (Y) 2D regions, out-of-plane width 15 in. Their struts
    use Concrete02: an unconfined COVER law (fc=5) on the outer 1.5-in ring, a confined CORE law
    (fc=6) inside — keeping the benchmark's two concrete grades (fork 3, "fc=6/5").
  - The beam is an ELASTIC concrete region (fork 2).
  - Longitudinal reinforcement: vertical Steel02 struts at x=+-10.5 (1.8 in^2 each) and x=0
    (1.2 in^2) in each column, full height, on shared lattice nodes (perfect bond, fork 4 / D13).
  - Mesh = 1.5 in so nodes land exactly on x = +-10.5 and the 1.5-in cover ring is one cell.

The concrete strut area is calibrated elastically (lattice <-> continuum lateral stiffness, the
Stage-1 machinery) so the initial stiffness is physical before the nonlinear run.

Specimen imported from this package (specimen.py). Exact replication is NOT expected (axial lattice
+ Concrete02 vs a fiber beam-column); the aim is the right shape (elastic -> yield -> plateau) and
capacity ballpark. Output: examples/output/frame/frame_pushover_rc.png. Units: kip, in. NOTE: mesh
1.5 -> several-thousand-node model; the nonlinear run takes a few minutes.
Run as `python examples/frame/pushover_rc.py`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import build_continuum, build_lattice, build_lattice_rc, select_nodes
from rclattice.materials import concrete_uniaxial_elastic, concrete_uniaxial_regularized
from rclattice.opensees import run_benchmark_rc_frame, run_pushover

from specimen import (
    BEAM_C, COL, CORE, COVER_C, DU, EPS, GF, GFC, H, OUT, P_COL, SPAN, TARGET,
    frame_problem, lateral_loads, rebars, zone_of,
)

MESH = 1.5


def calibrate_area() -> float:
    """Elastic strut area matching the continuum lateral stiffness (Stage-1 machinery)."""
    problem = frame_problem(CORE)
    lat0, _ = build_lattice(problem, MESH, strut_area=1.0)
    cont, _ = build_continuum(problem, MESH)
    ctrl = select_nodes(lat0, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(lat0, (-COL, SPAN + COL, -EPS, EPS))
    kw = dict(lateral_loads=lateral_loads(lat0), control_node=ctrl, control_dof=1,
              dU=DU, target=DU, base_nodes=base)  # one elastic step is enough for the slope
    k_cont = (lambda r: r["shear"][-1] / r["disp"][-1])(run_pushover(cont, **kw))
    k_unit = (lambda r: r["shear"][-1] / r["disp"][-1])(run_pushover(lat0, **kw))
    return k_cont / k_unit


def main(*, draw: bool = False) -> None:
    outdir = OUT / "pushover_rc" / "beamcolumn"
    outdir.mkdir(parents=True, exist_ok=True)

    area = calibrate_area()
    print(f"calibrated concrete strut area A = {area:.4e} in^2")

    problem = frame_problem(CORE)

    def material_for(zone: str, length: float):
        if zone == "beam":
            return concrete_uniaxial_elastic(BEAM_C, 0)
        grade = CORE if zone == "core" else COVER_C
        return concrete_uniaxial_regularized(grade, 0, length, Gf=GF, Gfc=GFC)

    model, _ = build_lattice_rc(problem, MESH, material_for=material_for, zone_of=zone_of,
                               rebars=rebars(), strut_area=area)
    ctrl = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(model, (-COL, SPAN + COL, -EPS, EPS))
    print(f"RC lattice: {len(model.nodes)} nodes, {len(model.elements)} struts | control {ctrl}")

    if draw:
        drawpath = outdir / "frame_pushover_rc_model.png"
        viz.figure_model([("RC lattice", model)], savepath=str(drawpath),
                         suptitle="RC frame (Stage 2) — analysis model")
        print(f"saved model drawing to {drawpath}")

    lat = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)
    print(f"lattice pushover: converged={lat['converged']} | drift {lat['disp'][-1]:.2f} in | "
          f"peak V {max(lat['shear']):.1f} kip")

    bench = run_benchmark_rc_frame(dU=DU, target=TARGET, gravity_P=P_COL)
    print(f"benchmark pushover: converged={bench['converged']} | drift {bench['disp'][-1]:.2f} in | "
          f"peak V {max(bench['shear']):.1f} kip")

    viz.figure_pushover(
        [
            {"disp": bench["disp"], "shear": bench["shear"], "label": "benchmark (fiber beam-column)",
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat["disp"], "shear": lat["shear"], "label": "RC lattice",
             "style": {"color": "C0", "ls": "--", "lw": 2}},
        ],
        savepath=str(outdir / "frame_pushover_rc.png"),
        xlabel="roof displacement (in)", ylabel="base shear (kip)",
        title="RC frame pushover (Stage 2): lattice vs benchmark",
    )
    print(f"saved pushover curve to {outdir / 'frame_pushover_rc.png'}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RC frame pushover (Stage 2): lattice vs benchmark")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model (rebar highlighted)")
    args = parser.parse_args()
    main(draw=args.draw)
