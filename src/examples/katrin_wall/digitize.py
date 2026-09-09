"""Digitize the measured WSH hysteresis loops from Figure 7 of Dazio, Beyer & Bachmann (2009).

Extracts the actual LOOPS, so a cyclic lattice run can be overlaid loop-for-loop against the test.
Figure 7 is a 3x2 grid of panels (a..f = WSH1..WSH6) embedded in the PDF as ONE 1418x774 grayscale
raster at 200 dpi; this reads that image at its native resolution rather than re-rendering the page.

It is a point CLOUD, not an ordered path. Overlapping loops cannot be re-sequenced from pixels, so
the result is drawn as a scatter and reduced to an envelope, never integrated for hysteretic energy.

THE HARD PART IS TEXT, not the curve. The figure is monochrome, so colour cannot separate the
annotation (unit name, load-step numbers, ductility labels, legend) from the data. Three filters
run in sequence, each targeting a different way a glyph presents itself:

  1. a 2x2 EROSION finds glyph-like seeds — compact blobs that survive erosion, unlike the 1-2 px
     curve stroke — and each seed's padded bounding box is cut out. This is what kills glyphs that
     TOUCH the curve and would otherwise be inseparable from it.
  2. leftover glyph fragments — and the ductility tick stubs on the zero line — are dropped as
     small COMPACT components (bbox <= 20 px square), once the axis rules that held them onto the
     curve have been cut away. A loop arc severed by step 1 is long and thin, so it survives this
     test where a digit or a tick stub does not.
  3. the cloud is clipped to |V| <= 1.03 * V_max from Table 5. V_max is the largest shear MEASURED
     during the test, so no genuine pixel can lie above it; the 3% allows for the ~4 kN half-width
     of the printed line. This removes the legend row and any glyph that survived 1-2 above a tip.

Calibration is self-checking, and the check is the paper's own numbers: the recovered peak load is
printed against Table 5's V_max and the recovered peak displacement against Table 4's delta_u. At
200 dpi one pixel is 0.59 mm of displacement and 4.3 kN of load, which bounds the accuracy.

A note on the two checks: the LOAD check should always pass to within a pixel. The DISPLACEMENT
check passes for WSH3 and WSH6 but reads high for WSH1/2/4/5, and correctly so — delta_u is the
displacement at FAILURE (a 20% strength drop), and the figure keeps plotting the post-failure
branch that the test continued onto. It is a sanity bound, not an equality.

Requires `pdfimages` (poppler). Run from src/:
    python examples/katrin_wall/digitize.py [--unit WSH3] [--pdf PATH] [--plot]
Output: examples/katrin_wall/data/<unit>_fig7.npz  (drift_pct, disp_mm, load_kN + backbone)
        examples/output/katrin_wall/<unit>_fig7_digitized.png  with --plot
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

DEFAULT_PDF = ("/Users/cansanliturk/Desktop/Master/Literature/Shearwall/"
               "Quasi_static_cyclic_tests_and_plastic_hinge_analysis_of_RC_structural_walls.pdf")
FIG7_PAGE = 6                       # journal page 1561; Fig. 7 is the only image on it
DATA = Path(__file__).resolve().parent / "data"

# Panel axes, read from the figure's own tick labels (identical for all six panels).
MARGIN = 10                         # px stripped inside the frame: the major ticks point inward
DISP_TICKS = 25.0                   # mm between major x ticks (-100..100)
LOAD_TICKS = 200.0                  # kN between major y ticks (-600..600)

# (panel row, panel column) in the 2x3 grid, plus the paper's own values for the calibration check.
# V_max: Table 5. delta_u, L_v: Table 4 / Sec. 2.1 (WSH6's shear span is 4520, the rest 4560).
# dy is the 3/4-rule yield displacement the LOADING HISTORY was built on (Table 4a) — the cycle
# amplitudes are its integer multiples, which is what anchors the backbone windows below.
UNITS = {
    "WSH1": dict(cell=(0, 0), V_max=336.0, du=47.5, Lv=4560.0, dy=10.5),
    "WSH2": dict(cell=(0, 1), V_max=359.0, du=63.0, Lv=4560.0, dy=10.5),
    "WSH3": dict(cell=(0, 2), V_max=454.0, du=92.4, Lv=4560.0, dy=15.4),
    "WSH4": dict(cell=(1, 0), V_max=443.0, du=61.6, Lv=4560.0, dy=15.4),
    "WSH5": dict(cell=(1, 1), V_max=439.0, du=62.0, Lv=4560.0, dy=6.2),
    "WSH6": dict(cell=(1, 2), V_max=597.0, du=93.7, Lv=4520.0, dy=12.8),
}

# Annotation that survives all three generic filters, cut out by hand. Each entry is a data-space
# box (disp_lo, disp_hi, load_lo, load_hi) and each is verified to clear the curve: the load band
# starts ABOVE the highest genuine pixel in that column range, with ~2 px to spare.
ANNOTATION_BOXES = {
    "WSH3": [(43.0, 50.5, 452.0, 500.0)],   # the "26" load-step label, drawn over its leader line
}


def extract_figure(pdf: str, page: int) -> np.ndarray:
    """Pull Fig. 7's embedded raster out of the PDF at its native 200 dpi, as grayscale."""
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdfimages", "-f", str(page), "-l", str(page), "-png", pdf, f"{td}/f"],
                       check=True)
        pngs = sorted(Path(td).glob("f*.png"))
        if len(pngs) != 1:
            raise RuntimeError(f"expected one image on page {page}, found {len(pngs)}")
        return np.asarray(Image.open(pngs[0]).convert("L")).astype(int)


def _rules(index: list[int]) -> list[float]:
    """Collapse runs of adjacent indices into one centre each — a drawn rule is 2-3 px wide."""
    out, run = [], [index[0]]
    for i in index[1:]:
        if i - run[-1] <= 2:
            run.append(i)
        else:
            out.append(sum(run) / len(run))
            run = [i]
    out.append(sum(run) / len(run))
    return out


def _longest_run(flags: np.ndarray) -> int:
    best = cur = 0
    for f in flags:
        cur = cur + 1 if f else 0
        best = max(best, cur)
    return best


def panel_frames(fig: np.ndarray) -> tuple[list[float], list[float]]:
    """(horizontal rule rows, vertical rule columns) of the whole 2x3 grid.

    A frame or zero rule is one UNBROKEN dark line spanning most of a panel. Selecting on the
    longest contiguous run rather than on the total dark count is what separates the rules from the
    axis-label text, whose columns are equally dark but only ever a few pixels tall.
    """
    dark = fig < 200
    h, w = dark.shape
    cols = _rules([j for j in range(w) if _longest_run(dark[:, j]) >= 0.25 * h])
    rows = _rules([i for i in range(h) if _longest_run(dark[i, :]) >= 0.22 * w])
    if len(rows) != 6 or len(cols) != 9:
        raise RuntimeError(f"expected 6 rows and 9 columns of rules, got {len(rows)}, {len(cols)}")
    return rows, cols


def _ticks(dark: np.ndarray, along: str, base: int, step: int, lo: int, hi: int, want: int) -> list[float]:
    """Major tick centres on one axis, found as short stubs growing inward from the frame line.

    `base`/`step` set where the stub starts and which way it grows INWARD from that frame edge.
    Length-bounded on BOTH sides: a tick is 4-9 px, so neither a frame/zero rule nor a curve branch
    that happens to run into the frame (WSH6 does) can be mistaken for one. The bound also excludes
    the ZERO tick, which is drawn as part of the zero rule — hence 8 x ticks and 6 y ticks, not
    9 and 7; the origin comes from the zero rules themselves, not from the ticks.
    """
    hits = []
    for k in range(lo, hi + 1):
        line = dark[:, k] if along == "x" else dark[k, :]
        n = 0
        while 0 <= base + step * n < len(line) and line[base + step * n]:
            n += 1
        if 4 <= n <= 9:
            hits.append(k)
    got = _rules(hits)
    if len(got) > want:
        got = _on_grid(got, want)
    if len(got) != want:
        raise RuntimeError(f"expected {want} {along} ticks, got {len(got)}")
    return got


def _on_grid(candidates: list[float], want: int) -> list[float]:
    """Keep the `want` candidates that lie on a uniform grid; drop the rest.

    Major ticks are evenly spaced by construction. A curve branch that grazes the frame with a
    stub of tick-like length (WSH6 has two) lands off that grid, so the grid is the discriminator.
    """
    step = np.median(np.diff(candidates))
    best: list[float] = []
    for anchor in candidates:
        on = [c for c in candidates if abs((c - anchor) / step - round((c - anchor) / step)) < 0.15]
        if len(on) > len(best):
            best = on
    return best


def calibrate(fig: np.ndarray, panel: tuple[float, float, float, float, float, float]) -> tuple[float, float]:
    """(px_per_kN, px_per_mm) from the panel's own major tick marks."""
    top, bot, _zr, left, right, _zc = (int(round(v)) for v in panel)
    dark = fig < 200
    xs = _ticks(dark, "x", bot - 2, -1, left, right, 8)     # -100..-25, +25..+100 mm
    ys = _ticks(dark, "y", left + 2, +1, top, bot, 6)       # +600..+200, -200..-600 kN
    return (ys[-1] - ys[0]) / (6 * LOAD_TICKS), (xs[-1] - xs[0]) / (8 * DISP_TICKS)


def _curve_pixels(panel: np.ndarray, zero_row: float, zero_col: float) -> np.ndarray:
    """Boolean mask of data pixels: everything dark that is neither annotation nor axis furniture."""
    dark = panel < 210
    conn = np.ones((3, 3), bool)

    # (1) glyph seeds -> cut out their padded bounding boxes
    seeds, _ = ndimage.label(ndimage.binary_erosion(dark, structure=np.ones((2, 2))), structure=conn)
    text = np.zeros_like(dark)
    for i, sl in enumerate(ndimage.find_objects(seeds)):
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        size = int((seeds[sl] == i + 1).sum())
        if h <= 15 and w <= 15 and size >= 6 and size / (h * w) >= 0.30:
            text[max(0, sl[0].start - 2):sl[0].stop + 2, max(0, sl[1].start - 2):sl[1].stop + 2] = True
    keep = dark & ~text

    # (2) axis furniture, removed BEFORE the component pass and not after — the two zero rules carry
    # the ductility tick stubs, which sit at exactly the loop-tip displacements. Detaching them
    # first is what lets step (3) see them as the small free-standing marks they are.
    zr, zc = int(round(zero_row)), int(round(zero_col))
    keep[zr - 2:zr + 3, :] = False
    keep[:, zc - 2:zc + 3] = False
    keep[:MARGIN, :] = keep[-MARGIN:, :] = False
    keep[:, :MARGIN] = keep[:, -MARGIN:] = False

    # (3) leftover fragments: small AND compact. A loop arc severed above is small but NOT compact.
    comps, _ = ndimage.label(keep, structure=conn)
    for i, sl in enumerate(ndimage.find_objects(comps)):
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        size = int((comps[sl] == i + 1).sum())
        if size < 8 or (h <= 20 and w <= 20 and size <= 220):
            keep &= comps != i + 1

    return keep


def digitize(unit: str = "WSH3", pdf: str = DEFAULT_PDF, page: int = FIG7_PAGE) -> dict:
    spec = UNITS[unit]
    fig = extract_figure(pdf, page)
    rows, cols = panel_frames(fig)
    r, c = spec["cell"]
    geom = (rows[3 * r], rows[3 * r + 2], rows[3 * r + 1], cols[3 * c], cols[3 * c + 2], cols[3 * c + 1])
    px_kN, px_mm = calibrate(fig, geom)
    top, bot, zero_row, left, right, zero_col = (int(round(v)) for v in geom)

    panel = fig[top:bot + 1, left:right + 1]
    keep = _curve_pixels(panel, zero_row - top, zero_col - left)
    ys, xs = np.nonzero(keep)
    load = (zero_row - top - ys) / px_kN
    disp = (xs - (zero_col - left)) / px_mm

    # (3) physical clip: no pixel can sit above the largest shear the test actually recorded.
    inside = np.abs(load) <= 1.03 * spec["V_max"]
    for d_lo, d_hi, v_lo, v_hi in ANNOTATION_BOXES.get(unit, []):
        inside &= ~((disp >= d_lo) & (disp <= d_hi) & (load >= v_lo) & (load <= v_hi))
    load, disp = load[inside], disp[inside]

    print(f"{unit}: panel rows {top}-{bot}, cols {left}-{right}; "
          f"{1 / px_mm:.2f} mm/px, {1 / px_kN:.2f} kN/px")
    print(f"  {len(load)} curve pixels ({(~inside).sum()} removed by the V_max clip and "
          f"{len(ANNOTATION_BOXES.get(unit, []))} hand-cut annotation boxes)")
    print(f"  CALIBRATION CHECK  peak load  +{load.max():.0f} / {load.min():.0f} kN "
          f"vs Table 5 V_max {spec['V_max']:.0f} kN")
    print(f"  CALIBRATION CHECK  peak disp  +{disp.max():.1f} / {disp.min():.1f} mm "
          f"vs Table 4 delta_u {spec['du']:.1f} mm")
    return {"disp_mm": disp, "load_kN": load, "drift_pct": disp / spec["Lv"] * 100.0}


def backbone(disp: np.ndarray, load: np.ndarray, dy: float, window: float = 3.0) -> dict:
    """Peak load at each cycle amplitude — the loop-TIP locus, which is the comparable curve.

    Anchored on the protocol rather than found from the cloud: the amplitudes are known exactly
    (integer multiples of the 3/4-rule yield displacement, Sec. 2.3), so each tip is read from a
    narrow window around a known displacement instead of from a running maximum that any stray
    pixel can lift. `window` covers the actuator's overshoot, which ran to about 2 mm.
    """
    amps, hi, lo, mus = [], [], [], []
    for mu in range(1, 9):
        a = mu * dy
        if a > np.abs(disp).max():
            break
        pos = load[np.abs(disp - a) <= window]
        neg = load[np.abs(disp + a) <= window]
        if len(pos) < 3 or len(neg) < 3:
            continue
        mus.append(mu)
        amps.append(a)
        hi.append(pos.max())
        lo.append(neg.min())
    return {"bb_mu": np.array(mus), "bb_disp_mm": np.array(amps),
            "bb_push_kN": np.array(hi), "bb_pull_kN": np.array(lo)}


def plot(unit: str, d: dict) -> Path:
    """Digitized cloud + envelope against the paper's own tabulated peaks — the visual check."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spec = UNITS[unit]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.scatter(d["disp_mm"], d["load_kN"], s=0.6, c="#1f4e79", alpha=0.5, lw=0, label="digitized")
    ax.plot(d["bb_disp_mm"], d["bb_push_kN"], "o-", c="#c0392b", lw=1.4, ms=4, label="backbone (loop tips)")
    ax.plot(-d["bb_disp_mm"], d["bb_pull_kN"], "o-", c="#c0392b", lw=1.4, ms=4)
    for sign in (+1, -1):
        ax.axhline(sign * spec["V_max"], color="0.45", ls="--", lw=0.8)
        ax.axvline(sign * spec["du"], color="0.45", ls=":", lw=0.8)
    ax.plot([], [], color="0.45", ls="--", lw=0.8, label=f"Table 5  V_max {spec['V_max']:.0f} kN")
    ax.plot([], [], color="0.45", ls=":", lw=0.8, label=f"Table 4  $\\Delta_u$ {spec['du']:.1f} mm")
    ax.set_xlabel("top displacement [mm]")
    ax.set_ylabel("actuator force [kN]")
    ax.set_title(f"{unit} — Fig. 7 digitized (Dazio, Beyer & Bachmann 2009)")
    top = ax.secondary_xaxis("top", functions=(lambda x: x / spec["Lv"] * 100,
                                               lambda x: x * spec["Lv"] / 100))
    top.set_xlabel("drift [%]")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    out = Path(__file__).resolve().parent.parent / "output" / "katrin_wall"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{unit.lower()}_fig7_digitized.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def main(unit: str = "WSH3", pdf: str = DEFAULT_PDF, page: int = FIG7_PAGE,
         make_plot: bool = False) -> None:
    d = digitize(unit, pdf, page)
    d |= backbone(d["disp_mm"], d["load_kN"], UNITS[unit]["dy"])
    d["Lv_mm"] = np.array(UNITS[unit]["Lv"])
    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA / f"{unit.lower()}_fig7.npz"
    np.savez_compressed(out, **d)
    print(f"saved {out} ({out.stat().st_size / 1024:.0f} kB)")
    if make_plot:
        print(f"saved {plot(unit, d)}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Digitize a WSH hysteresis from Fig. 7")
    p.add_argument("--unit", default="WSH3", choices=sorted(UNITS))
    p.add_argument("--pdf", default=DEFAULT_PDF)
    p.add_argument("--page", type=int, default=FIG7_PAGE)
    p.add_argument("--plot", action="store_true", help="also save a check figure")
    a = p.parse_args()
    main(a.unit, a.pdf, a.page, a.plot)
