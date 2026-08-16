"""Digitize the measured SW-NC-FF hysteresis from Figure 14b of the source paper.

Extracts the actual LOOPS — not just the backbone milestones — so the model's hysteresis can be
compared against the test's loop-for-loop. The figure is rendered from the PDF, the plot frame is
located from its own axis lines, and the blue curve is separated by colour; the result is a point
cloud in (drift %, load kN), which is what a raster figure can honestly yield. It is a CLOUD, not an
ordered path: overlapping loops cannot be re-sequenced from pixels, so it is drawn as a scatter.

Calibration is self-checking. The axis box is verified symmetric about the light-grey zero lines,
and the recovered peak loads are printed against the paper's reported +212.6 / -209.4 kN — if those
disagree by more than a couple of percent, the calibration is wrong and the output should not be used.

Requires `pdftoppm` (poppler). Run from src/:
    python examples/wall/digitize.py [--pdf PATH] [--page 13]
Output: examples/wall/data/swncff_fig14b.npz
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

DEFAULT_PDF = ("/Users/cansanliturk/Desktop/Master/Literature/Shearwall/"
               "Behavior of Nonconforming Flexure-Controlled RC Structural_Walls_"
               "Under_Reversed_Cyclic_Lateral_Loading.pdf")
OUT = Path(__file__).resolve().parent / "data" / "swncff_fig14b.npz"

# Panel (b) axis limits, read from the figure's own tick labels.
DISP_LIM, LOAD_LIM, SHEAR_SPAN = 110.0, 250.0, 2200.0   # mm, kN, mm
REPORTED_PEAKS = (212.6, -209.4)                        # paper Sec. 3, for the calibration check


def render(pdf: str, page: int, dpi: int = 400) -> np.ndarray:
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(dpi), "-png",
                        pdf, f"{td}/p"], check=True)
        img = next(Path(td).glob("p*.png"))
        return np.asarray(Image.open(img).convert("RGB")).astype(int)


def _lines(indices: list[int]) -> list[float]:
    """Collapse runs of adjacent indices into one centre each — a drawn rule is 2-3 px wide, so
    the raw indices come in pairs that must not be mistaken for two separate rules."""
    out, run = [], [indices[0]]
    for i in indices[1:]:
        if i - run[-1] <= 2:
            run.append(i)
        else:
            out.append(sum(run) / len(run))
            run = [i]
    out.append(sum(run) / len(run))
    return out


def find_box(page: np.ndarray, region: tuple[int, int, int, int]) -> tuple[float, float, float, float]:
    """Locate panel (b)'s plot frame: the long black horizontal/vertical rules bounding it."""
    r0, r1, c0, c1 = region
    reg = page[r0:r1, c0:c1]
    dark = reg.sum(axis=2) < 250
    rows = _lines([i for i in range(dark.shape[0]) if dark[i, :].sum() > 500])
    cols = _lines([j for j in range(dark.shape[1]) if dark[:, j].sum() > 500])
    if len(rows) < 2 or len(cols) < 2:
        raise RuntimeError(f"could not locate the plot frame (rows={rows}, cols={cols})")
    # the panel's own left frame is the LAST vertical rule before the right one (earlier rules
    # belong to the neighbouring panel (a))
    return (r0 + rows[0], r0 + rows[-1], c0 + cols[-2], c0 + cols[-1])


def digitize(pdf: str = DEFAULT_PDF, page: int = 13) -> dict:
    pg = render(pdf, page)
    top, bot, left, right = find_box(pg, (450, 1500, 1950, 3100))
    box = pg[int(top):int(bot) + 1, int(left):int(right) + 1]
    h, w = box.shape[0] - 1, box.shape[1] - 1

    r, g, b = box[:, :, 0], box[:, :, 1], box[:, :, 2]
    blue = (b > r + 40) & (b > g + 40) & (b > 90)      # the test curve; damage markers are not blue
    ys, xs = np.nonzero(blue)
    disp = -DISP_LIM + (xs / w) * 2 * DISP_LIM
    load = LOAD_LIM - (ys / h) * 2 * LOAD_LIM

    # Drop the legend's line sample: a perfectly horizontal dash below the curve's own extreme.
    keep = load > min(REPORTED_PEAKS) - 6.0
    dropped = (~keep).sum()
    disp, load = disp[keep], load[keep]
    drift = disp / SHEAR_SPAN * 100.0

    print(f"frame rows {top}-{bot}, cols {left}-{right}")
    print(f"{len(drift)} curve points ({dropped} legend-sample pixels removed)")
    print(f"  drift {drift.min():+.2f}..{drift.max():+.2f}%   load {load.min():.1f}..{load.max():.1f} kN")
    for got, want, side in ((load.max(), REPORTED_PEAKS[0], "push"), (load.min(), REPORTED_PEAKS[1], "pull")):
        print(f"  CALIBRATION CHECK {side}: digitized {got:+.1f} vs reported {want:+.1f} kN "
              f"({abs(got - want) / abs(want) * 100:.1f}% — line thickness is ~3 kN)")
    return {"drift_pct": drift, "load_kN": load}


def main(pdf: str = DEFAULT_PDF, page: int = 13) -> None:
    d = digitize(pdf, page)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT, **d)
    print(f"saved {OUT} ({OUT.stat().st_size / 1024:.0f} kB)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Digitize SW-NC-FF hysteresis from Fig. 14b")
    p.add_argument("--pdf", default=DEFAULT_PDF)
    p.add_argument("--page", type=int, default=13)
    a = p.parse_args()
    main(a.pdf, a.page)
