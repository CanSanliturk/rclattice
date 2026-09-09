"""Redraw VK3 figures from a saved run — never re-run the analysis.

A cyclic run costs many hours, so every analysis script writes its full response (and both gauge
chains) to its own timestamped run directory. This reads one back and regenerates the hysteresis,
the backbone, the strain profile, the curvature profile and the deformation-component split, so a
change to a figure costs seconds instead of a night.

Run from src/:
    python examples/vk3_wall/replot.py [--run vk3_cyclic_dynamic] [--levels 0.5 0.95 1.27]
"""

from __future__ import annotations

import argparse
import json
import pathlib

from rclattice import viz

import gauge
from specimen import A_SHEAR, LW, MEASURED_LOOPS, OUT, SPECIMEN, find_run, neutral_axis


def load(name: str | None, kind: str) -> tuple[dict, "pathlib.Path"]:
    """The run's saved response plus the directory it came from (figures go back beside it)."""
    d = find_run(kind, name=name)
    path = d / "data.json"
    if not path.exists():
        raise SystemExit(f"{d.name} has no data.json — it is not a completed analysis run")
    return json.loads(path.read_text()), d


def measured():
    if not MEASURED_LOOPS.exists():
        return None
    import numpy as np
    return np.load(MEASURED_LOOPS)


def main(name: str | None, kind: str, levels: list[float] | None) -> None:
    d, out = load(name, kind)
    g = d.get("gauge") or {}
    gx = g.get("x_mm", [])
    drift, shear = d["drift_pct"], d["shear_kN"]
    print(f"{out.name}: {len(drift)} steps, solver={d.get('solver')}, "
          f"compression={d.get('compression')}, converged={d.get('converged')}, "
          f"peak {max(shear):.1f} / {min(shear):.1f} kN")

    m = measured()
    series = [{"disp": drift, "shear": shear, "label": "lattice",
               "style": {"color": "C0", "lw": 0.9}}]
    if m is not None:
        series.append({"disp": list(m["drift_pct"]), "shear": list(m["load_kN"]),
                       "label": f"{SPECIMEN} measured (digitized Fig. 5.13)",
                       "style": {"color": "C3", "lw": 0.6, "alpha": 0.75}})
    viz.figure_hysteresis(series, savepath=str(out / "replot.png"),
                          drift_label="drift ratio (%)", shear_label="base shear (kN)",
                          title=f"{SPECIMEN} — {out.name} (redrawn)")
    print(f"saved {out / 'replot.png'}")

    if not gx or not g.get("drift_pct"):
        print("no gauge data in this run — nothing else to redraw")
        return

    g_drift = g["drift_pct"]
    targets = levels or sorted({round(abs(v), 2) for v in
                                (0.25, 0.48, 0.64, 0.95, 1.27, 1.59)})
    idx = gauge.select_levels(g_drift, targets)
    rows = gauge.figure(g_drift, g["strain"], gx, indices=idx,
                        savepath=out / "replot_strain.png",
                        gauge=g.get("gauge_mm"),
                        title=f"{SPECIMEN} — vertical strain across the base section")
    gauge.print_table(rows)
    print(f"saved {out / 'replot_strain.png'}")

    if g.get("curvature"):
        crows = gauge.figure_curvature(g_drift, g["mid_mm"], g["curvature"], indices=idx,
                                       savepath=out / "replot_curvature.png",
                                       title=f"{SPECIMEN} — curvature over height")
        gauge.print_curvature_table(crows)
        print(f"saved {out / 'replot_curvature.png'}")

    if g.get("components"):
        crow = gauge.figure_components(g_drift, g["top_mm"], g["components"],
                                       savepath=out / "replot_components.png",
                                       title=f"{SPECIMEN} — deformation components "
                                             f"(cf. Fig. 5.19 right)")
        gauge.print_components_table(crow, every=max(1, len(crow) // 12))
        print(f"saved {out / 'replot_components.png'}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Redraw VK3 figures from a saved run")
    p.add_argument("--run", default=None,
                   help="run directory name under examples/output/vk3_wall/runs "
                        "(default: the newest of --kind)")
    p.add_argument("--kind", default="cyclic", choices=("cyclic", "pushover", "elastic"),
                   help="which kind of run to pick when --run is not given")
    p.add_argument("--levels", type=float, nargs="*", default=None,
                   help="drift levels (%%) to draw profiles at")
    a = p.parse_args()
    main(a.run, a.kind, a.levels)
