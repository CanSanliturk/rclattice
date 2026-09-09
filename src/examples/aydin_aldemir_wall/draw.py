"""Render the lattice discretization of the Aldemir et al. (2017) squat wall, for `--draw`.

Builds the same calibrated model the analyses use and runs NO analysis. Three panels: the full
elevation, a zoom on the base corner, and — when `--bond` is given — a zoom tight enough to actually
see the bond rings, which is the whole point of looking at this model rather than reading a number.

Run from src/:
    uv run python examples/aydin_aldemir_wall/draw.py [--mesh 50] [--horizon 1.5]
                                                      [--bond] [--panel table2] [--nonlinear]
Output: examples/output/aydin_aldemir_wall/aldemir_model.png
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from rclattice.viz import draw_model_kinds  # noqa: E402

from build import calibrate, describe, wall_lattice  # noqa: E402
from specimen import (COVER, HORIZON, HW, LW, MESH, OUT, PAPER_GRID, PAPER_GRID_T,  # noqa: E402
                      S_BAR, SPECIMEN, TW)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mesh", type=float, default=MESH)
    ap.add_argument("--horizon", type=float, default=HORIZON)
    ap.add_argument("--calibration", choices=("uniaxial", "equibiaxial"), default="uniaxial")
    ap.add_argument("--panel", choices=("fig10a", "table2", "table2t"), default="fig10a",
                    help="fig10a = the drawn specimen (default); table2 = the grid his Table 2 "
                         "counts imply, table2t = its transposition, which is equally consistent "
                         "with both integers and is unresolved (D82). Fig. 4(f) is NOT a panel — it "
                         "is a detail view of the bottom 1500 mm of his model")
    ap.add_argument("--nonlinear", action="store_true", help="draw the nonlinear model instead")
    args = ap.parse_args()

    length, height = {"fig10a": (LW, HW), "table2": PAPER_GRID,
                      "table2t": PAPER_GRID_T}[args.panel]
    cal = calibrate(mesh_size=args.mesh, horizon=args.horizon, field=args.calibration)
    model, _edges = wall_lattice(cal.area, mesh_size=args.mesh, horizon=args.horizon,
                                 nonlinear=args.nonlinear, bond=False,
                                 length=length, height=height)
    print(f"{SPECIMEN} [{args.panel}]  {length:.0f} x {height:.0f} x {TW:.0f}")
    print(f"  {describe(model)}")
    print(f"  A_t = {cal.area:,.1f} mm^2 [{cal.field}], mesh {args.mesh:g}, horizon {args.horizon:g}"
          + (f", bond ON" if False else ", perfect bond"))

    # Zoom windows: the base corner (where the diagonal strut lands) and, with bond on, a few cells.
    zooms = [("base corner", (0.0, 6 * args.mesh, 0.0, 6 * args.mesh), "all")]
    if False:
        x0 = COVER - 1.2 * args.mesh
        zooms.append(("bond detail — steel + bond only\n(each link lies along a concrete strut)",
                      (x0, x0 + 3.2 * args.mesh, 0.0, 3.2 * args.mesh), "rebar"))

    fig, axes = plt.subplots(1, 1 + len(zooms),
                             figsize=(5.5 + 4.2 * len(zooms), 6.4),
                             gridspec_kw={"width_ratios": [1.35] + [1.0] * len(zooms)})
    draw_model_kinds(axes[0], model, which="all",
                     title=f"{length:.0f} x {height:.0f} x {TW:.0f}  "
                           f"(Ø8 @ {S_BAR:.0f} both ways, cover {COVER:.0f})")
    for ax, (label, (x0, x1, y0, y1), which) in zip(axes[1:], zooms):
        draw_model_kinds(ax, model, which=which, title=label, lim=((x0, x1), (y0, y1)))

    fig.suptitle(f"{SPECIMEN} — lattice model  |  mesh {args.mesh:g} mm, horizon {args.horizon:g}, "
                 f"{args.calibration} calibration, "
                 f"{'bond elements' if False else 'perfect bond'}", fontsize=10)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "aldemir_model.png"
    fig.savefig(path, dpi=170)
    print(f"  -> {path}")


if __name__ == "__main__":
    main()
