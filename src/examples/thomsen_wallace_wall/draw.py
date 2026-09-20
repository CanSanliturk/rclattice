"""Render the lattice discretization of RW2 on either grid, for `--draw`.

Run from src/:
    uv run python examples/thomsen_wallace_wall/draw.py [--grid rebar|uniform] [--mesh 25]
                                                         [--horizon 1.5] [--nonlinear]
Output: examples/output/thomsen_wallace_wall/rw2_model_<grid><mesh>.png
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from rclattice.viz import draw_model_kinds  # noqa: E402

from build import calibrate, describe, wall_lattice  # noqa: E402
from specimen import GRID, HORIZON, HW, LW, MESH, OUT, SPECIMEN, TW  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", choices=("rebar", "uniform"), default=GRID)
    ap.add_argument("--mesh", type=float, default=MESH)
    ap.add_argument("--horizon", type=float, default=HORIZON)
    ap.add_argument("--nonlinear", action="store_true")
    args = ap.parse_args()

    cal = calibrate(mesh_size=args.mesh, horizon=args.horizon)
    model, _edges = wall_lattice(cal.area, grid=args.grid, mesh_size=args.mesh,
                                 horizon=args.horizon, nonlinear=args.nonlinear)
    print(f"{SPECIMEN}  {LW:.0f} x {HW:.0f} x {TW:.0f}, {args.grid} grid at {args.mesh:g} mm")
    print(f"  {describe(model)}")

    zooms = [("base, left boundary element", (0.0, 400.0, 0.0, 400.0), "all"),
             ("reinforcement only, base", (0.0, 1220.0, 0.0, 600.0), "rebar")]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 7.2), gridspec_kw={"width_ratios": [0.6, 1.0, 1.4]})
    draw_model_kinds(axes[0], model, which="all",
                     title=f"{LW:.0f} x {HW:.0f} x {TW:.0f}")
    for ax, (label, (x0, x1, y0, y1), which) in zip(axes[1:], zooms):
        draw_model_kinds(ax, model, which=which, title=label, lim=((x0, x1), (y0, y1)))
    fig.suptitle(f"{SPECIMEN} — lattice model  |  {args.grid} grid, target {args.mesh:g} mm, "
                 f"horizon {args.horizon:g}", fontsize=10)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"rw2_model_{args.grid}{args.mesh:g}.png"
    fig.savefig(path, dpi=170)
    print(f"  -> {path}")


if __name__ == "__main__":
    main()
