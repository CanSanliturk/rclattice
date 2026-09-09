"""Digitize the measured VK3 hysteresis from Figs. 5.10-5.13 of Bimschas (2010).

Bimschas, M. (2010), "Displacement Based Seismic Assessment of Existing Bridges in Regions of
Moderate Seismicity", IBK Bericht Nr. 326, ETH Zuerich. doi:10.3929/ethz-a-006237119.

UNLIKE the WSH3 digitizer, this reads VECTOR geometry, not pixels. The hysteresis panels in this
chapter are not rasters: they are PDF path objects with explicit stroke colours. Three consequences,
all of them improvements:

  * no digitization error at all — the recovered coordinates are the ones the plot was drawn with,
    to within the calibration of the axes;
  * no text/curve separation problem — the measured curve is the only GREEN geometry on the panel,
    so the axis rules, tick stubs, load-step numbers and legend never have to be filtered out (the
    three-stage erosion filter the WSH3 digitizer needs exists solely to solve that);
  * the loops stay ORDERED. A raster hysteresis is a point cloud that cannot be re-sequenced, which
    is why the WSH3 data carries the warning that it must never be integrated for hysteretic energy.
    Here the drawing order IS the loading order, so energy dissipation, residual displacement and
    per-cycle degradation are all recoverable.

HOW THE CURVE IS RECONSTRUCTED. Each panel's curve is drawn as ~11 separate green polylines, one
per load-step block, and they are stored in REVERSE order. They are re-chained by matching each
series' start point to another's end point (`_chain`), which is order-independent and self-checking:
every junction must close to well under a plotted line width, or the chain is rejected.

AXES ARE CALIBRATED ON THE TICK MARKS, NOT ON THE FRAME, and that is not fussiness. The frame of
these panels spans +-110 mm x +-950 kN while the printed tick labels stop at +-100 and +-800, so
assuming the frame equals the axis limits would inflate every load by 19%. Major ticks are twice the
length of minor ones, which is what separates them.

The calibration then CHECKS ITSELF against a completely independent axis: each panel carries a
secondary base-moment axis, and 200 kN of shear at L_v = 3.3 m must be 0.66 MN.m of moment. The two
axes are drawn from different data, so agreement to a fraction of a percent means the load scale is
right. A second check compares the recovered peak displacement against the protocol amplitude the
figure is captioned with.

Requires `pdftocairo` (poppler). Run from src/:
    python examples/vk3_wall/digitize.py [--unit VK3] [--figure 5.13] [--pdf PATH] [--plot]
Output: examples/vk3_wall/data/<unit>_fig<n>.npz
        examples/output/vk3_wall/<unit>_fig<n>_digitized.png   with --plot
"""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

DEFAULT_PDF = ("/Users/cansanliturk/Desktop/Master/Literature/Shearwall/"
               "Displacement_based_seismic_assessment_of_existing_bridges_in_regions_of_"
               "moderate_seismicity_CHAPTER_5.pdf")
DATA = Path(__file__).resolve().parent / "data"
OUT = Path(__file__).resolve().parent.parent / "output" / "vk3_wall"

# The chapter is paginated from 301, so PDF page = printed page - 300.
PAGE_OFFSET = 300
# figure -> (printed page, panel units left-to-right, captioned peak displacement mm)
FIGURES = {
    "5.10": (324, ("VK1", "VK2", "VK3"), 21.0),
    "5.11": (325, ("VK1", "VK2", "VK3"), 31.5),
    "5.12": (326, ("VK1", "VK2", "VK3"), 42.0),
    "5.13": (327, ("VK1", "VK2", "VK3"), 52.5),      # VK3's LAST level — the complete history
}
CURVE_RGB = "0%, 100%, 0%"          # the measured hysteresis; blue/red are the peak markers
L_V = 3300.0                        # mm, the chapter's shear span and drift denominator

# Axis metadata, identical on every panel (read off the printed tick labels).
DISP_PER_TICK = 20.0                # mm between major ticks on the displacement axis
LOAD_PER_TICK = 200.0               # kN between major ticks on the shear-force axis
MOMENT_PER_TICK = 1.0               # MN.m between major ticks on the secondary moment axis

# The protocol's STRUCTURE (Fig. 5.9), which is documented even where its numbers are not:
#   * an ELASTIC block of four FORCE-controlled levels (0.25/0.50/0.75/1.00 F_y'), two cycles each
#     and no small cycles -> the first 16 turning points. Their DISPLACEMENTS are whatever the pier
#     gave, which is exactly why they are worth recovering here.
#   * an INELASTIC block at mu_prov * 10.5 mm, two cycles each, with two small intermediate cycles
#     inserted between them from mu_prov = 1.5 up.
ELASTIC_TIPS = 16
INELASTIC_LEVELS = (10.5, 15.75, 21.0, 31.5, 42.0, 52.5)
NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


# --- SVG parsing ----------------------------------------------------------------------------------
def _paths(svg_text: str):
    """`[(stroke_rgb_or_None, [(x, y), ...])]` in page points.

    pdftocairo emits every path with its own `transform="matrix(a,0,0,-a,tx,ty)"`, so each path's
    coordinates are local and must be mapped individually. Beziers are reduced to their endpoints:
    these are polyline plots, so any `C` present is a rendering artefact of a straight segment.
    """
    body = svg_text.split("</defs>", 1)[1]
    out = []
    for p in re.findall(r'<path\b.*?/>', body, re.S):
        d = re.search(r'\sd="([^"]*)"', p)
        if not d:
            continue
        t = re.search(r'transform="matrix\(([^)]*)\)"', p)
        a, b, c, dd, tx, ty = ([float(v) for v in t.group(1).split(',')] if t
                               else [1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        col = re.search(r'stroke="rgb\(([^)]*)\)"', p)
        pts = []
        for op, nums in re.findall(rf'([MLC])\s*((?:{NUM}[\s,]*)+)', d.group(1)):
            v = [float(x) for x in re.findall(NUM, nums)]
            if op == 'C':
                v = v[-2:]
            for i in range(0, len(v) - 1, 2):
                pts.append((a * v[i] + c * v[i + 1] + tx, b * v[i] + dd * v[i + 1] + ty))
        if pts:
            out.append((col.group(1) if col else None, pts))
    return out


def _page_svg(pdf: Path, page: int) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        dst = Path(tmp) / "page.svg"
        subprocess.run(["pdftocairo", "-svg", "-f", str(page), "-l", str(page), str(pdf), str(dst)],
                       check=True, capture_output=True)
        return dst.read_text()


# --- panel geometry and axis calibration -----------------------------------------------------------
def _frames(paths):
    """The plot frames on the page, left to right: `[(x0, x1, y0, y1)]`."""
    got = set()
    for col, p in paths:
        if col != "0%, 0%, 0%" or len(p) > 8:
            continue
        xs = [x for x, _y in p]
        ys = [y for _x, y in p]
        if max(xs) - min(xs) > 80 and max(ys) - min(ys) > 50:
            got.add((round(min(xs), 2), round(max(xs), 2), round(min(ys), 2), round(max(ys), 2)))
    return sorted(got)


def _major_ticks(paths, frame, axis: str):
    """Positions of the MAJOR ticks on one axis of one panel.

    Major ticks are drawn at twice the length of minor ticks, which is the only thing distinguishing
    them; everything longer than 1.5x the shortest tick found is taken as major.
    """
    x0, x1, y0, y1 = frame
    found = []
    for col, p in paths:
        if col != "0%, 0%, 0%" or len(p) > 3:
            continue
        xs = [x for x, _y in p]
        ys = [y for _x, y in p]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if axis == "x" and h < 6 and w < 0.6 and abs(max(ys) - y1) < 4:
            # constrain to THIS panel: all panels share y0/y1, so without an x window the ticks of
            # all three are pooled and the spacing test sees the gaps between panels
            if x0 - 2.0 <= max(xs) <= x1 + 2.0:
                found.append((max(xs), h))
        elif axis in ("y", "y2") and w < 6 and h < 0.6:
            anchor = x0 if axis == "y" else x1
            if abs((min(xs) if axis == "y" else max(xs)) - anchor) < 4:
                found.append((max(ys), w))
    if not found:
        raise RuntimeError(f"no {axis}-axis ticks found for frame {frame}")
    short = min(L for _p, L in found if L > 0)
    major = sorted({round(pos, 3) for pos, L in found if L > 1.5 * short})
    if len(major) < 3:
        raise RuntimeError(f"only {len(major)} major ticks on {axis} for frame {frame}")
    return major


def _scale(major, per_tick):
    """`(origin_in_points, units_per_point)` from evenly spaced major ticks."""
    diffs = [b - a for a, b in zip(major, major[1:])]
    spacing = sorted(diffs)[len(diffs) // 2]
    spread = (max(diffs) - min(diffs)) / spacing
    if spread > 0.02:
        raise RuntimeError(f"major ticks are not evenly spaced (spread {spread:.1%})")
    return 0.5 * (major[0] + major[-1]), per_tick / spacing


def calibrate(paths, frame):
    """`(x0, mm_per_pt, y0, kN_per_pt, moment_check)` for one panel.

    `moment_check` is the ratio between the load scale implied by the secondary base-moment axis and
    the one implied by the shear axis. The two are drawn from independent tick sets, so a value of
    1.000 means the load calibration is right for a reason unrelated to how it was obtained.
    """
    x0, mm_per_pt = _scale(_major_ticks(paths, frame, "x"), DISP_PER_TICK)
    y0, kn_per_pt = _scale(_major_ticks(paths, frame, "y"), LOAD_PER_TICK)
    try:
        _y2, mnm_per_pt = _scale(_major_ticks(paths, frame, "y2"), MOMENT_PER_TICK)
        implied_kn_per_pt = mnm_per_pt * 1e3 / (L_V / 1e3)      # MN.m/pt -> kN/pt at lever L_v
        check = implied_kn_per_pt / kn_per_pt
    except RuntimeError:
        check = float("nan")
    return x0, mm_per_pt, y0, kn_per_pt, check


# --- curve reconstruction ---------------------------------------------------------------------------
def _chain(series):
    """Re-order disjoint polylines into one continuous path by matching endpoints.

    Starts from the series whose first point is nearest the origin (the plot starts at zero load and
    zero displacement) and repeatedly appends whichever unused series begins closest to the current
    end. Returns `(points, worst_junction_gap)`; the caller rejects a large gap.
    """
    left = list(series)
    start = min(range(len(left)), key=lambda i: abs(left[i][0][0]) + abs(left[i][0][1]))
    chain = list(left.pop(start))
    worst = 0.0
    while left:
        end = chain[-1]
        j = min(range(len(left)),
                key=lambda i: (left[i][0][0] - end[0]) ** 2 + (left[i][0][1] - end[1]) ** 2)
        nxt = left.pop(j)
        gap = ((nxt[0][0] - end[0]) ** 2 + (nxt[0][1] - end[1]) ** 2) ** 0.5
        worst = max(worst, gap)
        chain.extend(nxt[1:])                 # drop the duplicated junction point
    return chain, worst


def backbone(disp, load):
    """Peak load at each new displacement extreme, both directions — the loop-tip envelope."""
    pos, neg = [], []
    hi, lo = 0.0, 0.0
    for u, s in zip(disp, load):
        if u > hi:
            hi, _ = u, pos.append((u, s))
        elif u < lo:
            lo, _ = u, neg.append((u, s))
    return pos, neg


def excursions(disp, load, deadband: float = 0.5):
    """Every turning point of the ordered history: `[(peak_disp_mm, load_there_kN)]`.

    Only possible because these are VECTOR paths in DRAWING order — a rastered hysteresis is an
    unordered cloud that cannot be cut into excursions at all, which is why the WSH3 data carries a
    warning against integrating it.

    A turning point is confirmed only once the displacement has retreated by `deadband` from the
    running extreme, so the actuator's small hunting movements inside a branch do not register as
    tips. 0.5 mm is ~2% of the smallest genuine amplitude and ~1% of the largest.

    WHAT THIS RECOVERS. The full VK3 history comes out as 60 turning points, and they contain the
    "small intermediate cycles" of Fig. 5.9 that the chapter does NOT tabulate for VK3 — it defines
    them only as "the top displacements measured during the corresponding cycles of VK1". They are
    plainly visible here, including the deliberately ASYMMETRIC second one (e.g. +11.8 / -2.0 mm at
    the 15.75 mm level, the displacement counterpart of VK1's +0.75F_y' / -0.25F_y'). So the model
    can be driven through the history the pier actually saw, rather than through an idealization of
    it — see `specimen.measured_protocol`.
    """
    out, ext, ei, dirn = [], disp[0], 0, 1
    for i, v in enumerate(disp):
        if dirn > 0:
            if v >= ext:
                ext, ei = v, i
            elif ext - v > deadband:
                out.append(ei); dirn = -1; ext, ei = v, i
        else:
            if v <= ext:
                ext, ei = v, i
            elif v - ext > deadband:
                out.append(ei); dirn = 1; ext, ei = v, i
    out.append(ei)
    return [(float(disp[i]), float(load[i])) for i in out]


def classify_tips(tips, rel_tol: float = 0.12):
    """Split turning points into PRIMARY (a protocol level) and INTERMEDIATE (a small cycle).

    Uses the protocol's documented STRUCTURE rather than trying to infer levels from the amplitudes,
    which does not work: the small cycles of several different levels happen to land near one
    another, so any clustering of amplitudes alone invents levels that were never targeted.

    The elastic block is taken as the first `ELASTIC_TIPS` turning points, grouped four to a level
    (two cycles x two directions) in time order, and each level is LABELLED WITH ITS OWN MEASURED
    amplitude — the whole point, since those displacements are not published. The inelastic block is
    matched against the known displacement targets; anything more than `rel_tol` away from all of
    them is a small intermediate cycle. The separation is clean rather than marginal: at the 15.75
    mm level the four primary tips land at 15.7-16.0 mm while the largest intermediate is 11.8 mm.

    Returns `(primary, intermediate, flags)`. `flags` is `(u, s, level, is_primary)` IN TIP ORDER —
    which is what gets stored, because a primary/intermediate split written in split order rather
    than tip order silently mis-aligns with `tip_disp_mm` and makes `--no-small-cycles` a no-op.
    """
    flags = []
    head, tail = tips[:ELASTIC_TIPS], tips[ELASTIC_TIPS:]
    for k in range(0, len(head), 4):
        block = head[k:k + 4]
        lv = round(sum(abs(u) for u, _s in block) / len(block), 2)
        flags += [(u, s, lv, True) for u, s in block]
    for u, s in tail:
        lv = min(INELASTIC_LEVELS, key=lambda v: abs(abs(u) - v))
        flags.append((u, s, lv, abs(abs(u) - lv) <= rel_tol * lv))
    primary = [(u, s, lv) for u, s, lv, ok in flags if ok]
    intermediate = [(u, s, lv) for u, s, lv, ok in flags if not ok]
    return primary, intermediate, flags


def cycle_backbones(primary):
    """The chapter's own presentation: 1st- and 2nd-cycle backbones in each direction.

    Fig. 5.19-left plots exactly four measured curves — 1st/2nd cycle x south/north — because the
    difference between them IS the cyclic degradation. Returns `{(cycle, 'S'|'N'): [(level, kN)]}`,
    cycles counted in TIME order within each level and direction.
    """
    out, seen = {}, {}
    for u, s, lv in primary:
        d = "S" if u > 0 else "N"
        n = seen.get((lv, d), 0) + 1
        seen[(lv, d)] = n
        if n <= 2:
            out.setdefault((n, d), []).append((lv, s))
    return {k: sorted(v) for k, v in out.items()}


def extract(pdf: Path, figure: str, unit: str):
    """Digitize one panel. Returns a dict of arrays plus the self-check numbers."""
    printed, units, captioned = FIGURES[figure]
    if unit not in units:
        raise SystemExit(f"{unit} is not on Fig. {figure} (panels: {', '.join(units)})")
    paths = _paths(_page_svg(pdf, printed - PAGE_OFFSET))
    frames = _frames(paths)
    if len(frames) != len(units):
        raise RuntimeError(f"found {len(frames)} panel frames on page {printed}, expected "
                           f"{len(units)}")
    frame = frames[units.index(unit)]
    x0, mm_per_pt, y0, kn_per_pt, moment_check = calibrate(paths, frame)

    fx0, fx1 = frame[0] - 2.0, frame[1] + 2.0
    green = [p for col, p in paths
             if col == CURVE_RGB and min(x for x, _y in p) >= fx0 and max(x for x, _y in p) <= fx1]
    if not green:
        raise RuntimeError(f"no curve geometry in the {unit} panel of Fig. {figure}")
    pts = [[((x - x0) * mm_per_pt, -(y - y0) * kn_per_pt) for x, y in p] for p in green]
    chained, gap = _chain(pts)
    disp = np.array([u for u, _s in chained])
    load = np.array([s for _u, s in chained])
    return {
        "disp_mm": disp, "load_kN": load, "drift_pct": disp / L_V * 100.0,
        "n_series": len(green), "junction_gap_mm": gap, "moment_check": moment_check,
        "captioned_mm": captioned, "frame": frame,
    }


def main(pdf: Path, figure: str, unit: str, plot: bool) -> None:
    d = extract(pdf, figure, unit)
    disp, load = d["disp_mm"], d["load_kN"]
    peak_p, peak_n = float(load.max()), float(load.min())
    reach_p, reach_n = float(disp.max()), float(disp.min())

    print(f"{unit}, Fig. {figure}: {d['n_series']} polylines chained into {len(disp):,} ordered "
          f"points")
    print(f"  peak base shear        +{peak_p:.1f} / {peak_n:.1f} kN")
    print(f"  displacement reached   +{reach_p:.2f} / {reach_n:.2f} mm  "
          f"(= {reach_p / L_V:+.2%} / {reach_n / L_V:.2%} drift)")
    print("  --- self-checks ---")
    print(f"  moment axis vs load axis   {d['moment_check']:.4f}   (independent axes; 1.0000 means "
          f"the load scale is right)")
    print(f"  captioned peak {d['captioned_mm']:.1f} mm    recovered "
          f"{max(reach_p, -reach_n):.2f} mm  -> {max(reach_p, -reach_n) / d['captioned_mm']:.4f}")
    print(f"  worst chain junction       {d['junction_gap_mm']:.4f} mm  (a plotted line is ~1 mm "
          f"wide at this scale)")
    if d["junction_gap_mm"] > 2.0:
        raise SystemExit("chain junctions are too wide — the series were not re-sequenced correctly")

    tips = excursions(list(disp), list(load))
    primary, intermediate, flags = classify_tips(tips)
    cyc = cycle_backbones(primary)
    tip_levels = sorted({lv for _u, _s, lv in primary})

    # The FIRST elastic level is an uncracked-stiffness measurement, and a valuable one: the
    # chapter's k0 = 60 kN/mm is a SECANT TO FIRST YIELD on a heavily cracked pier, so it is not
    # comparable to an elastic model at all. At the first level the pier carries ~150 kN, well under
    # the ~208 kN its own cracking moment implies, so the secant there IS an elastic stiffness.
    first = [(u, v) for u, v, lv in primary if lv == min(l for _a, _b, l in primary)]
    k_init = sum(abs(v / u) for u, v in first) / len(first) if first else float("nan")
    print(f"\n  initial (uncracked) secant  {k_init:.1f} kN/mm over the first elastic level "
          f"(+-{abs(first[0][0]):.2f} mm, ~{abs(first[0][1]):.0f} kN)")
    print(f"    the chapter's k0 = 60 kN/mm is a secant to FIRST YIELD, {60.0 / k_init:.2f}x this "
          f"— not an elastic stiffness")

    bb_pos, bb_neg = backbone(list(disp), list(load))

    DATA.mkdir(parents=True, exist_ok=True)
    path = DATA / f"{unit.lower()}_fig{figure.replace('.', '')}.npz"
    np.savez(
        path,
        disp_mm=disp, load_kN=load, drift_pct=d["drift_pct"],
        tip_disp_mm=np.array([u for u, _s in tips]),
        tip_load_kN=np.array([s for _u, s in tips]),
        tip_level_mm=np.array([lv for _u, _s, lv, _ok in flags]),
        tip_is_primary=np.array([ok for _u, _s, _lv, ok in flags]),
        bb_level_mm=np.array(tip_levels),
        bb_push_kN=np.array([dict(cyc.get((1, "S"), [])).get(lv, np.nan) for lv in tip_levels]),
        bb_pull_kN=np.array([dict(cyc.get((1, "N"), [])).get(lv, np.nan) for lv in tip_levels]),
        bb2_push_kN=np.array([dict(cyc.get((2, "S"), [])).get(lv, np.nan) for lv in tip_levels]),
        bb2_pull_kN=np.array([dict(cyc.get((2, "N"), [])).get(lv, np.nan) for lv in tip_levels]),
        env_push_disp_mm=np.array([u for u, _s in bb_pos]),
        env_push_kN=np.array([s for _u, s in bb_pos]),
        env_pull_disp_mm=np.array([u for u, _s in bb_neg]),
        env_pull_kN=np.array([s for _u, s in bb_neg]),
        **{f"cyc{c}{d}_{what}": np.array(vals) for (c, d), rows in cyc.items()
           for what, vals in (("mm", [r[0] for r in rows]), ("kN", [r[1] for r in rows]))},
        k_initial_kN_per_mm=k_init,
        Lv_mm=L_V, unit=unit, figure=figure,
        peak_push_kN=peak_p, peak_pull_kN=peak_n,
        moment_check=d["moment_check"], junction_gap_mm=d["junction_gap_mm"],
        ordered=True,
    )
    print(f"\n  turning points: {len(tips)} total = {len(primary)} at protocol levels + "
          f"{len(intermediate)} small intermediate cycles")
    print("\n  per-cycle backbones (Fig. 5.19-left's own four curves):")
    print(f"  {'level':>8}  {'drift':>7}  {'1st S':>8}  {'2nd S':>8}  {'1st N':>8}  {'2nd N':>8}"
          f"  {'2nd/1st S':>10}")
    for lv in tip_levels:
        def at(c, d, lv=lv):
            return dict(cyc.get((c, d), [])).get(lv)
        a, b, c_, e = at(1, "S"), at(2, "S"), at(1, "N"), at(2, "N")
        deg = f"{b / a:10.3f}" if (a and b) else f"{'--':>10}"
        cells = "  ".join(f"{v:8.1f}" if v is not None else f"{'--':>8}" for v in (a, b, c_, e))
        print(f"  {lv:8.2f}  {lv / L_V:6.2%}  {cells}  {deg}")
    print(f"\nsaved {path}")

    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        OUT.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8.2, 6.0))
        ax.plot(disp, load, lw=0.7, color="#2a78d6", label=f"{unit} measured (ordered)")
        for (c, dirn), rows in sorted(cyc.items()):
            sgn = 1.0 if dirn == "S" else -1.0
            ax.plot([sgn * lv for lv, _v in rows], [v for _lv, v in rows],
                    "o--" if c == 1 else "s:", ms=4, lw=1.3,
                    color="#d62728" if c == 1 else "#ff7f0e",
                    label=f"{c}{'st' if c == 1 else 'nd'} cycle, {dirn}")
        if intermediate:
            ax.plot([u for u, _s, _l in intermediate], [v for _u, v, _l in intermediate],
                    "x", color="#7a7a7a", ms=6, label="small intermediate cycles")
        ax.axhline(0, color="0.6", lw=0.8)
        ax.axvline(0, color="0.6", lw=0.8)
        ax.set_xlabel("top displacement (mm)")
        ax.set_ylabel("base shear (kN)")
        ax.grid(True, color="#d8d7d2", lw=0.5)
        ax.set_axisbelow(True)
        ax.set_title(f"{unit} — measured hysteresis digitized from Fig. {figure} "
                     f"(vector paths, exact)", fontsize=11, loc="left")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        p = OUT / f"{unit.lower()}_fig{figure.replace('.', '')}_digitized.png"
        fig.savefig(p, dpi=170, facecolor="white")
        plt.close(fig)
        print(f"saved {p}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Digitize VK hysteresis loops (vector, not raster)")
    ap.add_argument("--pdf", type=Path, default=Path(DEFAULT_PDF))
    ap.add_argument("--figure", default="5.13", choices=sorted(FIGURES),
                    help="5.13 (default) is VK3's last level and contains its whole history")
    ap.add_argument("--unit", default="VK3", choices=("VK1", "VK2", "VK3"))
    ap.add_argument("--plot", action="store_true", help="also save a check figure")
    a = ap.parse_args()
    if not a.pdf.exists():
        raise SystemExit(f"source PDF not found: {a.pdf}")
    main(a.pdf, a.figure, a.unit, a.plot)
