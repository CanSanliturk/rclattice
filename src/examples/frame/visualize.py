"""Visualize lattice vs continuum on a one-bay, one-storey portal frame (D-frame).

Builds the SAME portal-frame Problem two ways (lattice + continuum), calibrates the lattice
strut area to the continuum static response, then renders:
  - static deflected shapes side-by-side (should match well),
  - first N mode shapes side-by-side (with periods),
  - an animated GIF of the mode shapes.

This is a standalone SI study (different geometry/units from the kip-in benchmark), so it defines
its own frame; only the package output dir `OUT` is shared with the benchmark scripts. Outputs to
examples/output/frame/. Run as `python examples/frame/visualize.py`. Units: SI.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import build_continuum, build_lattice
from rclattice.calibration import static_response
from rclattice.opensees import run_modal, run_static
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, portal_frame

from specimen import OUT

MESH = 0.1
PLANE = "PlaneStress"
N_MODES = 3
H, SPAN, DC, DB, T = 3.0, 4.0, 0.3, 0.3, 0.3
FX = 50e3  # N, lateral (sway) load at beam level


def frame_problem() -> Problem:
    domain = portal_frame(height=H, span=SPAN, col_depth=DC, beam_depth=DB, thickness=T)
    grade = ConcreteGrade(name="C30", E=30e9, nu=0.2, rho=2400.0)
    supports = [BoxSupport(box=(-1.0, SPAN + 1.0, -1e-3, 1e-3), fix=(1, 1))]  # fix both bases
    loads = [BoxLoad(box=(-1.0, SPAN + 1.0, H - 1e-3, H + DB + 1e-3), total=(FX, 0.0))]  # sway
    return Problem(ndm=2, ndf=2, domain=domain, material=grade, supports=supports, loads=loads)


def main(*, draw: bool = False) -> None:
    outdir = OUT / "visualize" / "continuum"
    outdir.mkdir(parents=True, exist_ok=True)
    problem = frame_problem()

    # static-scalar area calibration (lateral dof = 0): deflection ~ 1/area
    cont0, _ = build_continuum(problem, MESH, plane=PLANE)
    delta_c = static_response(cont0, 0)
    lat0, _ = build_lattice(problem, MESH, strut_area=1.0)
    area = static_response(lat0, 0) / delta_c
    print(f"nodes = {len(cont0.nodes)} | calibrated strut area A = {area:.4e} m^2")

    lattice, _ = build_lattice(problem, MESH, strut_area=area)
    continuum, _ = build_continuum(problem, MESH, plane=PLANE)

    if draw:
        viz.figure_model([("lattice", lattice), ("continuum", continuum)],
                         savepath=str(outdir / "frame_model.png"), suptitle="Portal frame — analysis model")
        print(f"saved model drawing to {outdir / 'frame_model.png'}")

    lat_static = run_static(lattice)["disps"]
    cont_static = run_static(continuum)["disps"]
    viz.figure_static(lattice, lat_static, continuum, cont_static, savepath=str(outdir / "frame_static.png"))

    lat_modal = run_modal(lattice, N_MODES)
    cont_modal = run_modal(continuum, N_MODES)
    panels = [
        {
            "lattice": (lattice, lat_modal["shapes"][i]),
            "continuum": (continuum, cont_modal["shapes"][i]),
            "T_lat": lat_modal["periods"][i],
            "T_cont": cont_modal["periods"][i],
        }
        for i in range(N_MODES)
    ]
    viz.figure_modes(panels, savepath=str(outdir / "frame_modes.png"))
    viz.animate_modes(panels, savepath=str(outdir / "frame_modes.gif"))

    print("continuum periods (s):", [f"{t:.4f}" for t in cont_modal["periods"]])
    print("lattice   periods (s):", [f"{t:.4f}" for t in lat_modal["periods"]])
    print(f"saved figures + gif to {outdir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Portal frame lattice-vs-continuum visualization")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis models (undeformed)")
    args = parser.parse_args()
    main(draw=args.draw)
