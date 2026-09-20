"""The digitized record and how a run is compared against it (Aldemir PLAN.md §10, D77; D103).

**COMPARE AT MATCHED DISPLACEMENT, never peak-to-peak**: the model and the test reach their peaks
at different drifts, so a peak/peak ratio silently compares two different states.

The `.npz` CONTRACT a specimen's `digitize.py` writes:

  * `cloud`     (N, 2) — every digitized test point, mm and kN, unordered
  * `mono_x`, `mono_y` — the test ENVELOPE (running maximum of the cloud per branch), mm / kN
  * `backbone`  (N, 2) — the same envelope as one array, for the run sheets
  * the author's own simulated curves, one `<key>_x` / `<key>_y` pair per horizon, named by
    `author_curves` below (e.g. `{1.5: "aydin15", 3.01: "aydin301"}`)

The limits of the figure (a clipping frame, an occluding legend) are carried by `clip_mm` and
`note`, and repeated wherever a number from here is quoted.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class ReferenceSet:
    path: Path                                  # the .npz
    measured_F_kN: float                        # the published peak
    measured_K_kNmm: float                      # the published stiffness
    author_by_horizon: dict[float, tuple[float, float]]   # horizon -> (K kN/mm, F kN) published
    author_curves: dict[float, str] = field(default_factory=dict)  # horizon -> npz key prefix
    clip_mm: float = float("inf")               # where the figure frame cuts the record
    record_extent_mm: tuple[float, float] = (0.0, 0.0)   # (min, max) displacement the record covers
    note: str = ""                              # the figure's limits, in one paragraph
    figure: str = ""                            # e.g. "Fig. 10(b)"
    f_source: str = "Table 4"                   # where the measured peak is printed
    k_source: str = "Table 4"                   # where the measured stiffness is printed
    test_label: str = ""                        # e.g. "Aldemir et al. 2017"
    author_label: str = "Aydin"                 # who made the published simulation

    def load(self) -> dict | None:
        if not self.path.is_file():
            return None
        d = np.load(self.path)
        return {k: d[k] for k in d.files}

    def author_key(self, horizon: float) -> str | None:
        return self.author_curves.get(float(horizon))

    @staticmethod
    def envelope_at(x, y, x0: float, *, tol: float = 0.4) -> float | None:
        """The largest |force| within `tol` mm of `x0`, signed by the branch.

        A digitized hysteresis is a LOOP, not a function of displacement, so it cannot be
        interpolated directly; reading its envelope in a window is the honest reduction.
        """
        x, y = np.asarray(x, float), np.asarray(y, float)
        m = np.abs(x - x0) <= tol
        if not m.any():
            return None
        return float(y[m].max() if x0 >= 0 else y[m].min())

    def compare(self, disp_mm: float, shear_kn: float, *, horizon: float) -> dict:
        """Model against every reference, AT THIS DISPLACEMENT."""
        out = {"at_mm": disp_mm, "model_kN": shear_kn, "clipped": abs(disp_mm) > self.clip_mm}
        ref = self.load()
        if ref is not None:
            test = self.envelope_at(ref["mono_x"], ref["mono_y"], disp_mm)
            key = self.author_key(horizon)
            aydin = (self.envelope_at(ref[key + "_x"], ref[key + "_y"], disp_mm)
                     if key and key + "_x" in ref else None)
            out["test_envelope_kN"] = test
            out["aydin_kN"] = aydin
            if test:
                out["model_over_test"] = shear_kn / test
            if aydin:
                out["model_over_aydin"] = shear_kn / aydin
        out["table4_peak_kN"] = self.measured_F_kN
        out["model_over_table4_peak"] = shear_kn / self.measured_F_kN
        _k, f = self.author_by_horizon.get(float(horizon), (None, None))
        out["aydin_table4_kN"] = f
        return out

    def comparison_points(self, disp, shear, *, horizon: float, n: int = 5) -> list[dict]:
        """A ladder of matched-displacement comparisons up to the model's own maximum.

        **A CYCLIC RUN REVISITS EVERY DISPLACEMENT, so "the shear at u" is not a number.** The
        model side is reduced the SAME way the test side is: the largest force within a window of
        the target displacement, i.e. its envelope. On a monotonic push each displacement occurs
        once and the window returns that point, so the two cases agree by construction.
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
            rows.append(self.compare(t, max(near) / 1e3, horizon=horizon))
            rows[-1]["envelope_of"] = len(near) if cyclic else None
        return rows
