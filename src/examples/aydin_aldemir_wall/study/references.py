"""The digitized Fig. 10(b) and Table 4, and how a run is compared against them (PLAN.md §10).

**COMPARE AT MATCHED DISPLACEMENT, never peak-to-peak** (D77's method note): the model and the test
reach their peaks at different drifts, so a peak/peak ratio silently compares two different states.

TWO LIMITS OF THE FIGURE CARRY INTO EVERY COMPARISON, and are repeated wherever a number from here
is quoted:

* the frame **CLIPS the data at +16 mm** while the test reached about 20 mm — 582 dots pile up at
  the edge — so the envelope is an UPPER BOUND on the backbone past ~15 mm, not the backbone;
* the legend occludes the negative quadrant at positive displacement.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data" / "fig10b.npz"

MEASURED_F_KN = 963.592      # Table 4
MEASURED_K_KNMM = 1038.44    # Table 4
CLIP_MM = 16.0               # where the figure frame cuts the record
AYDIN = {1.5: (943.16, 1164.413), 3.01: (1052.37, 1325.675)}   # horizon -> (K kN/mm, F kN)


def load() -> dict | None:
    if not DATA.is_file():
        return None
    d = np.load(DATA)
    return {k: d[k] for k in d.files}


def envelope_at(x, y, x0: float, *, tol: float = 0.4) -> float | None:
    """The largest |force| within `tol` mm of `x0`, signed by the branch.

    A digitized hysteresis is a LOOP, not a function of displacement, so it cannot be interpolated
    directly; reading its envelope in a window is the honest reduction. `tol` is a few pixels wide.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.abs(x - x0) <= tol
    if not m.any():
        return None
    return float(y[m].max() if x0 >= 0 else y[m].min())


def compare(disp_mm: float, shear_kn: float, *, horizon: float) -> dict:
    """Model against every reference, AT THIS DISPLACEMENT."""
    out = {"at_mm": disp_mm, "model_kN": shear_kn, "clipped": abs(disp_mm) > CLIP_MM}
    ref = load()
    if ref is not None:
        test = envelope_at(ref["mono_x"], ref["mono_y"], disp_mm)
        aydin = envelope_at(ref["aydin15_x"], ref["aydin15_y"], disp_mm) if horizon == 1.5 else \
            envelope_at(ref["aydin301_x"], ref["aydin301_y"], disp_mm)
        out["test_envelope_kN"] = test
        out["aydin_kN"] = aydin
        if test:
            out["model_over_test"] = shear_kn / test
        if aydin:
            out["model_over_aydin"] = shear_kn / aydin
    out["table4_peak_kN"] = MEASURED_F_KN
    out["model_over_table4_peak"] = shear_kn / MEASURED_F_KN
    k, f = AYDIN.get(horizon, (None, None))
    out["aydin_table4_kN"] = f
    return out


def comparison_points(disp, shear, *, horizon: float, n: int = 5) -> list[dict]:
    """A ladder of matched-displacement comparisons up to the model's own maximum.

    **A CYCLIC RUN REVISITS EVERY DISPLACEMENT, so "the shear at u" is not a number.** This wall's
    1.0% ladder passes |u| = 4.5 mm 29,607 times carrying anywhere from -817 to +847 kN, and taking
    the nearest index — which is what this did before 2026-09-06 — returns whichever crossing came
    first, typically a point on an unloading branch. The result looked like the model carrying 39 kN
    where the test carried 969.

    So the model side is reduced the SAME way the test side is: the largest force within a window of
    the target displacement, i.e. its envelope. On a monotonic push each displacement occurs once
    and the window returns that point, so the two cases agree by construction.
    """
    dmax = max(abs(d) for d in disp) if disp else 0.0
    if dmax <= 0:
        return []
    win = max(dmax * 0.01, 0.05)          # a window, not a point
    cyclic = min(disp) < -1e-9
    rows = []
    for i in range(n):
        t = dmax * (i + 1) / n
        near = [abs(shear[j]) for j in range(len(disp)) if abs(abs(disp[j]) - t) <= win]
        if not near:
            continue
        rows.append(compare(t, max(near) / 1e3, horizon=horizon))
        rows[-1]["envelope_of"] = len(near) if cyclic else None
    return rows
