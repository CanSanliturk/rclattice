"""Digitize Fig. 10(b) of Aydin, Tuncay & Binici (2019) — the Aldemir wall load-deflection panel.

    uv run python examples/aydin_aldemir_wall/digitize.py [--check]

Writes `data/fig10b.npz`: the measured hysteresis point cloud, a backbone reduced from it, and
BOTH of the paper's own lattice curves.

WHY THIS EXISTS, AND WHAT AN EARLIER NOTE GOT WRONG. The repo previously recorded this figure as
unrecoverable ("an unordered raster point cloud, both standard routes return noise"). That is half
right and was over-generalized. What is genuinely lost is ORDER: the experimental data is drawn as
discrete dot markers, not a polyline, so which dot follows which is gone and with it per-cycle
energy, degradation and the loading path — the same limit as SW-NC-FF's Fig. 14b. What is NOT lost
is the cloud's OUTLINE. Loop tips are extreme points, extremes survive unordered, and a backbone is
an envelope. The figure also carries both of Aydin's OWN curves as clean coloured lines, which are
straightforwardly recoverable and are the direct comparison for our replica pushover.

CALIBRATION IS SELF-CHECKING. Axes are fixed from detected gridlines alone — nine horizontal rules
for the load labels (+-1400 by 350) and eight vertical for the displacement labels (-19..16 by 5) —
with no reference to any published value. The resulting scale then reproduces Table 4's two
simulated peaks to 0.999 (horizon 1.5) and 0.993 (horizon 3.01), and the plotted "experiment"
marker to within 1% of Table 4's measured maximum. Three independent numbers the calibration never
saw, so agreement is a test rather than a fit.

TWO LOSSES, both structural and both reported by `--check`:
  * THE LEGEND BOX occludes the paper's own data in the lower right (positive displacement,
    load below about -500 kN). Nothing can recover what the legend was printed over.
  * TICK LABELS sit INSIDE the axes. Glyphs are separated from dots morphologically (a dot is a
    ~6 px blob, a glyph stroke is 2-3 px, so an opening keeps one and drops the other) and then by
    position. Where dots merge INTO a label the component is sacrificed, which costs interior dots
    in the dense pinch near the origin. The envelope is unaffected, which is what the backbone uses.

Units: mm and kN, as printed on the figure.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
PDF = HERE / "Simulation_of_reinforced_concrete_member_response_using_lattice_model.pdf"
DATA = HERE / "data"
PAGE, PANEL = 12, 2          # Fig. 10(b) is the 3rd image object on page 12

# Printed axis labels. The ONLY external input to the calibration.
LOADS = np.arange(1400, -1401, -350.0)          # 9 horizontal rules, top to bottom
DISPS = np.arange(-19, 17, 5.0)                 # 8 vertical rules, left to right

# Table 4, for validation only — never used to set the scale.
F_SIM_15, F_SIM_301, F_EXP = 1164.413, 1325.675, 963.592


def panel() -> np.ndarray:
    """Extract the Fig. 10(b) bitmap from the PDF (cached under data/)."""
    DATA.mkdir(exist_ok=True)
    # NOT IN THE REPOSITORY, deliberately: this is a verbatim crop of Fig. 10(b) of the ASCE
    # paper, and the repo is public. The digitized OUTPUT (`fig10b.npz`) is committed, so nothing
    # downstream needs this file. To re-run the digitizing, crop the figure from the PDF yourself
    # and save it here at the same name.
    png = DATA / "fig10b_panel.png"
    if not png.exists():
        if not PDF.exists():
            sys.exit(f"missing {PDF}")
        stem = DATA / "_p"
        subprocess.run(["pdfimages", "-f", str(PAGE), "-l", str(PAGE), "-png", str(PDF), str(stem)],
                       check=True)
        src = sorted(DATA.glob("_p-*.png"))[PANEL]
        Image.open(src).save(png)
        for f in DATA.glob("_p-*.png"):
            f.unlink()
    return np.array(Image.open(png).convert("RGB")).astype(int)


def calibrate(im: np.ndarray) -> tuple:
    """Find the gridlines and return (x_of_px, y_of_px, geometry dict). No published value used."""
    gray = im.mean(2)
    H, W = gray.shape
    light = (gray > 150) & (gray < 235)
    dark = gray < 120

    def rules(mask_l, mask_d, axis, frac):
        """Gridline indices along `axis`, merging runs of adjacent lines to their centre."""
        hits = sorted(set(np.nonzero(mask_l.sum(axis) > frac * mask_l.shape[axis])[0])
                      | set(np.nonzero(mask_d.sum(axis) > 0.5 * mask_d.shape[axis])[0]))
        out, run = [], [hits[0]]
        for v in hits[1:]:
            (run.append(v) if v - run[-1] <= 2 else (out.append(np.mean(run)), run.clear(),
                                                     run.append(v)))
        out.append(np.mean(run))
        return np.array(out, float)

    rows = rules(light, dark, 1, 0.40)
    cols = rules(light, dark, 0, 0.40)
    # Keep the evenly spaced family: the true gridlines share one spacing; strays (legend edge,
    # frame) do not. Take the largest subset consistent with the modal spacing.
    def regular(v, n):
        best = None
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                s = (v[j] - v[i]) / (j - i)
                if s <= 0:
                    continue
                pred = v[i] + s * (np.arange(n) - i)
                if pred[0] < -5 or pred[-1] > max(v) + 5:
                    continue
                err = sum(min(abs(v - p)) for p in pred)
                if best is None or err < best[0]:
                    best = (err, pred, s)
        return best[1]

    yr, xr = regular(rows, len(LOADS)), regular(cols, len(DISPS))
    py2y = lambda p: np.interp(p, yr, LOADS)                                    # noqa: E731
    px2x = lambda p: np.interp(p, xr, DISPS)                                    # noqa: E731
    geo = dict(rows=yr, cols=xr,
               kN_per_px=(LOADS[0] - LOADS[-1]) / (yr[-1] - yr[0]),
               mm_per_px=(DISPS[-1] - DISPS[0]) / (xr[-1] - xr[0]),
               plot=(xr[0], xr[-1], yr[0], yr[-1]))
    return px2x, py2y, geo


def masks(im: np.ndarray, geo: dict) -> dict:
    """Colour classes, with the legend box removed (it occludes data, irrecoverably)."""
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]
    H, W = R.shape
    inside = np.zeros((H, W), bool)
    y0, y1 = int(geo["rows"][0]) + 3, int(geo["rows"][-1]) - 2
    x0, x1 = int(geo["cols"][0]) + 2, W - 5
    inside[y0:y1, x0:x1] = True
    legend = np.zeros((H, W), bool)
    legend[306:450, 358:683] = True
    ok = inside & ~legend
    return dict(
        red=((R > 120) & (R - G > 50) & (R - B > 50)) & ok,
        green=((G > 110) & (G - B > 40) & (R - B > 20) & (R < G + 40)) & ok,
        black=((R < 110) & (G < 110) & (B < 110)) & ok,
        legend=legend,
    )


def curve(mask: np.ndarray, px2x, py2y) -> tuple[np.ndarray, np.ndarray]:
    """A single-valued coloured line -> (x, y), one point per pixel column it occupies."""
    xs, ys = [], []
    for c in np.unique(np.nonzero(mask)[1]):
        rows = np.nonzero(mask[:, c])[0]
        xs.append(px2x(c))
        ys.append(py2y(rows.mean()))
    o = np.argsort(xs)
    return np.asarray(xs)[o], np.asarray(ys)[o]


def cloud(black: np.ndarray, geo: dict, px2x, py2y):
    """Experimental dots, separated from tick-label glyphs. Returns (pts, marker, dropped)."""
    disk = lambda r: (np.add.outer(*[(np.arange(-r, r + 1)) ** 2] * 2) <= r * r)   # noqa: E731
    # A dot is a ~6 px blob and a glyph stroke is 2-3 px, so an opening keeps dots and drops text.
    dots = ndi.binary_propagation(ndi.binary_erosion(black, disk(2)), mask=black)

    lab, n = ndi.label(dots, structure=np.ones((3, 3)))
    objs = ndi.find_objects(lab)
    area = ndi.sum(dots, lab, range(1, n + 1))

    # The plotted "experiment" marker is a filled disc, larger and denser than a dot. It cannot be
    # found by size alone — a chain of merged dots near the origin outweighs it — so it is taken as
    # the densest large blob in the upper-right quadrant, where the legend's third entry is drawn.
    #
    # TREAT THE VALUE IT YIELDS AS A COARSE CHECK, NOT A PRECISE ONE, for two reasons: its location
    # is assumed rather than derived, and the disc touches neighbouring data dots, so its centroid
    # is pulled off the true marker centre. The two coloured-curve plateaus are the strong checks.
    best, big = None, None
    for i, ob in enumerate(objs, 1):
        h, w = ob[0].stop - ob[0].start, ob[1].stop - ob[1].start
        fill = area[i - 1] / (h * w)
        cy, cx = (ob[0].start + ob[0].stop) / 2, (ob[1].start + ob[1].stop) / 2
        if area[i - 1] >= 300 and fill >= 0.45 and px2x(cx) > 10.0 and py2y(cy) > 500.0:
            if best is None or fill * area[i - 1] > best:
                best, big = fill * area[i - 1], i
    if big is None:
        big = int(np.argmax(area)) + 1
    o = objs[big - 1]
    marker = (px2x((o[1].start + o[1].stop) / 2), py2y((o[0].start + o[0].stop) / 2))

    # Tick labels live in two strips: left of the load axis, and just under the zero-load rule.
    ystrip = (geo["cols"][0] + 0.30 * (geo["cols"][-1] - geo["cols"][0]),
              geo["cols"][int(len(DISPS) / 2)] - 3)
    xstrip = (geo["rows"][len(LOADS) // 2] + 28, geo["rows"][len(LOADS) // 2] + 62)
    drop = {big}
    for i, ob in enumerate(objs, 1):
        h, w = ob[0].stop - ob[0].start, ob[1].stop - ob[1].start
        cy, cx = (ob[0].start + ob[0].stop) / 2, (ob[1].start + ob[1].stop) / 2
        in_y = ystrip[0] <= cx <= ystrip[1] and min(abs(geo["rows"] - cy)) <= 16
        in_x = xstrip[0] <= cy <= xstrip[1]
        if (in_y or in_x) and h >= 9:
            drop.add(i)
    clean = dots & ~np.isin(lab, list(drop))

    ys, xs = np.nonzero(clean)
    return np.column_stack([px2x(xs), py2y(ys)]), marker, len(drop) - 1


def backbone(pts: np.ndarray, step: float = 0.75) -> tuple[np.ndarray, np.ndarray]:
    """Loop-tip envelope: per displacement bin, the extreme load of each sign.

    This is what survives an unordered cloud. It is NOT a measured backbone curve — the test never
    published one — and it cannot give energy or degradation, only the outline the loops reach.
    """
    x, y = pts[:, 0], pts[:, 1]
    edges = np.arange(np.floor(x.min()), np.ceil(x.max()) + step, step)
    xs, up, lo = [], [], []
    for a, b in zip(edges, edges[1:]):
        k = (x >= a) & (x < b)
        if k.sum() < 3:
            continue
        xs.append(0.5 * (a + b))
        up.append(np.percentile(y[k], 99))
        lo.append(np.percentile(y[k], 1))
    return np.array(xs), np.column_stack([up, lo])


def monotone_backbone(bx: np.ndarray, be: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The half of each envelope that is physically a backbone, made monotone.

    Only two of the four envelope branches mean anything. Pushing in +x, the loop TIPS are the
    upper envelope; the lower envelope over the same range is the unloading side of a loop, not a
    capacity. Mirrored for -x. So the backbone is `upper for x>0` and `lower for x<0`.

    The protocol drove increasing amplitudes, so tips grow outward from the origin — a running
    extremum away from zero removes bin-to-bin dropouts where a bin happens to catch no tip,
    without inventing strength the cloud does not show.
    """
    pos, neg = bx > 0, bx < 0
    xp, yp = bx[pos], np.maximum.accumulate(be[pos, 0])
    xn = bx[neg][::-1]
    yn = np.minimum.accumulate(be[neg, 1][::-1])
    return np.concatenate([xn[::-1], xp]), np.concatenate([yn[::-1], yp])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="print validation and write an overlay")
    a = ap.parse_args()

    im = panel()
    px2x, py2y, geo = calibrate(im)
    mk = masks(im, geo)
    rx, ry = curve(mk["red"], px2x, py2y)
    gx, gy = curve(mk["green"], px2x, py2y)
    pts, marker, ndrop = cloud(mk["black"], geo, px2x, py2y)
    bx, be = backbone(pts)
    mx, my = monotone_backbone(bx, be)

    plateau = lambda x, y: y[x >= x.max() - 0.08 * (x.max() - x.min())].mean()    # noqa: E731
    p15, p301 = plateau(rx, ry), plateau(gx, gy)

    print(f"calibration: {geo['mm_per_px']:.4f} mm/px, {geo['kN_per_px']:.4f} kN/px")
    print(f"  gridlines: {len(geo['rows'])} load rules, {len(geo['cols'])} displacement rules")
    print("\nVALIDATION — three published numbers the calibration never saw:")
    print(f"  horizon 1.5  plateau {p15:9.1f} kN  vs Table 4 {F_SIM_15:9.3f}  -> {p15 / F_SIM_15:.4f}")
    print(f"  horizon 3.01 plateau {p301:9.1f} kN  vs Table 4 {F_SIM_301:9.3f}  -> {p301 / F_SIM_301:.4f}")
    print(f"  'experiment' marker  {marker[1]:9.1f} kN  vs Table 4 {F_EXP:9.3f}  -> {marker[1] / F_EXP:.4f}"
          f"   (at {marker[0]:+.2f} mm)")
    print("    ^ COARSE: its position is assumed and the disc touches nearby dots. The two curve")
    print("      plateaus above are the real checks — neither uses any published number.")
    print(f"\ncloud: {len(pts):,} px in the dot mask, {ndrop} glyph components removed")
    print(f"backbone: {len(bx)} bins, peak +{be[:, 0].max():.1f} / {be[:, 1].min():.1f} kN"
          f"   (vs Table 4 measured {F_EXP:,.1f} -> {be[:, 0].max() / F_EXP:.4f})")
    print("  NOTE this envelope is NOT a published backbone and carries no order, so per-cycle")
    print("  energy and degradation are NOT recoverable. The legend also occludes the lower right.")

    DATA.mkdir(exist_ok=True)
    out = DATA / "fig10b.npz"
    np.savez(out, cloud=pts, backbone_x=bx, backbone=be, mono_x=mx, mono_y=my,
             aydin15_x=rx, aydin15_y=ry, aydin301_x=gx, aydin301_y=gy,
             marker=np.array(marker), grid_rows=geo["rows"], grid_cols=geo["cols"])
    print(f"\nsaved {out}")

    if a.check:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(pts[:, 0], pts[:, 1], ".", ms=1.2, color="#333", label=f"cloud ({len(pts):,} px)")
        ax.plot(bx, be[:, 0], "-", lw=1, color="#9ad", alpha=.9, label="raw envelope (both halves)")
        ax.plot(bx, be[:, 1], "-", lw=1, color="#9ad", alpha=.9)
        ax.plot(mx, my, "-", lw=2.4, color="#16a085", label="backbone (monotone, meaningful half)")
        ax.plot(rx, ry, "-", lw=2, color="#c0392b", label=f"Aydin horizon 1.5 ({p15:,.0f} kN)")
        ax.plot(gx, gy, "-", lw=2, color="#7f9c3a", label=f"Aydin horizon 3.01 ({p301:,.0f} kN)")
        ax.plot(*marker, "o", ms=11, color="k", label=f"experiment marker ({marker[1]:,.0f} kN)")
        ax.set_xlabel("lateral displacement (mm)"); ax.set_ylabel("lateral load (kN)")
        ax.grid(alpha=.3); ax.legend(fontsize=8)
        ax.set_title("Fig. 10(b) digitized — check against the original", fontsize=10)
        fig.tight_layout(); fig.savefig(DATA / "fig10b_check.png", dpi=160)
        print(f"wrote {DATA / 'fig10b_check.png'}")


if __name__ == "__main__":
    main()
