"""Draw the strain field and damage pattern of a replica run — including runs that FAILED.

    uv run python examples/aydin_aldemir_wall/replica/damage.py            # newest run
    uv run python examples/aydin_aldemir_wall/replica/damage.py --run <dir> [--stage peak|latest]

Reads `snapshots.npz` (written on every progress tick while the run is going) or `capture.npz`
(written only if the run returned), rebuilds the model geometry, and renders:

  * the axial-STRAIN field, on a symmetric-log scale — once a crack localizes its strut strain is
    two to three orders above the surrounding elastic field, so a linear scale shows one red strut
    and nothing else;
  * the discrete DAMAGE pattern, struts past the cracking strain drawn as cracked.

WHY SNAPSHOTS AND NOT JUST `capture`. `capture=True` returns the field when a run RETURNS. Every
informative run in this study was killed mid-way — after a collapse, once a hypothesis was settled,
after a stall — so `capture` produced nothing for exactly the runs worth looking at. `run.py` now
dumps the field on every progress tick and keeps the peak alongside the latest, so an interrupted
run still renders. A collapsed run's LATEST field is the most interesting picture in the study.

TENSION ONLY, BY CONSTRUCTION. The replica's struts are Aydin's tension-only law: linear elastic in
compression forever, no crushing at any strain. So compression is NOT classified as damage here —
`eps_crush` is left None and the compressive side of the field is shown but never labelled
"crushed", which would assert a mechanism this material does not have.

CAUTION carried from D75: `viz.figure_damage`'s `eps_crush` wants a POSITIVE magnitude. Passing a
negative one flags ~99% of struts as crushed; the tell is a `< --2.3e-03` legend.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from rclattice import viz

from build import calibrate, wall_lattice
from specimen import A_SHEAR, EC, FT, HW, LW

OUT = Path(__file__).resolve().parent.parent.parent / "output" / "aydin_aldemir_wall" / "replica"
EPS_CRACK = FT / EC          # strut cracking strain, 7.44e-5


def newest_run() -> Path:
    runs = [d for d in OUT.glob("*") if d.is_dir()
            and ((d / "snapshots.npz").exists() or (d / "capture.npz").exists())]
    if not runs:
        raise SystemExit(f"no run under {OUT} has a snapshots.npz or capture.npz")
    return max(runs, key=lambda d: d.stat().st_mtime)


def load(run: Path) -> tuple[dict, list]:
    """Return (params, [(label, ids, u, meta), ...]) — richest source first."""
    par = json.loads((run / "params.json").read_text()) if (run / "params.json").exists() else {}
    stages = []
    if (run / "snapshots.npz").exists():
        z = np.load(run / "snapshots.npz", allow_pickle=True)
        ids = z["ids"]
        if "ring_u" in z.files and z["ring_u"].ndim == 3:
            for u, meta in zip(z["ring_u"], z["ring_meta"]):
                stages.append((f"t{int(meta[0]):,}", ids, u, meta))
        if z["u_peak"].ndim == 2:
            stages.append(("peak", ids, z["u_peak"], z["meta_peak"]))
        stages.append(("latest", ids, z["u_latest"], z["meta_latest"]))
    elif (run / "capture.npz").exists():
        z = np.load(run / "capture.npz")
        ids = z["ids"]
        stages.append(("peak", ids, z["u_peak"], None))
        stages.append(("final", ids, z["u_final"], None))
    else:
        raise SystemExit(f"{run.name}: no snapshots.npz or capture.npz")
    return par, stages


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, default=None, help="run directory (default: newest)")
    ap.add_argument("--stage", choices=("peak", "latest", "final", "both", "ring"), default="both",
                    help="'ring' renders the last N ticks in order — use it to watch a collapse "
                         "happen instead of only seeing its aftermath")
    ap.add_argument("--ring-last", type=int, default=4, help="how many ring frames to draw")
    ap.add_argument("--scale", type=float, default=0.0, help="deformed-shape amplification")
    a = ap.parse_args()

    run = a.run or newest_run()
    par, stages = load(run)
    ringed = [s for s in stages if s[0].startswith("t")]
    if a.stage == "ring":
        stages = ringed[-a.ring_last:] or stages[-1:]
    else:
        stages = [s for s in stages if not s[0].startswith("t")]
        if a.stage != "both":
            stages = [s for s in stages if s[0] == a.stage] or stages[-1:]
    # Early in a run the peak IS the latest tick; drawing both gives two identical rows.
    if (len(stages) == 2 and stages[0][3] is not None and stages[1][3] is not None
            and np.allclose(stages[0][3], stages[1][3])):
        stages = stages[-1:]
        print("  (peak == latest at this point in the run — drawing one stage)")

    mesh = par.get("mesh", 20.0)
    horizon = par.get("horizon", 1.5)
    print(f"run     : {run.name}")
    print(f"model   : mesh {mesh:g}, horizon {horizon:g}, t = {par.get('thickness', '?')}, "
          f"bond {par.get('bond', False)}")
    cal = calibrate(mesh_size=mesh, horizon=horizon)
    model, _e = wall_lattice(cal.area, mesh_size=mesh, horizon=horizon,
                             nonlinear=not par.get("elastic", False),
                             cap_rows=int(par.get("cap_rows", 0)),
                             full_height_rebar=bool(par.get("full_height_rebar", False)),
                             bond=bool(par.get("bond")),
                             bond_area=(par["bond_area_ratio"] * cal.area
                                        if par.get("bond_area_ratio") else None))
    print(f"rebuilt : {len(model.nodes):,} nodes, {len(model.elements):,} elements")

    panels = []
    for label, ids, u, meta in stages:
        disp = {int(i): list(map(float, row)) for i, row in zip(ids, u)}
        title = label
        if meta is not None and np.ndim(meta) == 1 and len(meta) == 3:
            step, dsp, shear = meta
            title = (f"{label}\nstep {int(step):,}\n{dsp / A_SHEAR:+.4%} drift\n"
                     f"{shear / 1e3:+,.0f} kN")
        _els, st = viz.strut_strains(model, disp)
        cracked = int((st > EPS_CRACK).sum())
        print(f"  {label:7s}: max eps {st.max():+.3e}, min {st.min():+.3e}, "
              f"{cracked:,} of {len(st):,} struts past eps_cr ({cracked / len(st):.2%})")
        panels.append((title, model, disp))

    dst = run / "damage.png"
    # eps_crush deliberately None: these struts have NO compressive strength to exceed (D60).
    viz.figure_damage(panels, eps_crack=EPS_CRACK, eps_crush=None, savepath=str(dst),
                      scale=a.scale, crack_label="cracked (eps > ft/E)",
                      suptitle=f"{run.name}\nAydin replica — mesh {mesh:g}, horizon {horizon:g}, "
                               f"tension-only struts (no crushing by construction)")
    print(f"\nsaved {dst}")


if __name__ == "__main__":
    main()
