"""Derived response metrics — the ones a raw series maximum gets wrong (D92).

`peak_shear` is a MAXIMUM OVER SAMPLES. Under an explicit march that samples every step (1.1 M
samples at dt = 8 us here) it is not a measure of resistance once a run carries failure switches:
when a bar ruptures or a strut releases, the energy it was holding rings in a handful of elements
at a LOCAL period — ~1.8 ms in this wall, against T1 = 5.8 ms — and the maximum lands on a crest of
that ringing rather than on the resistance the wall actually offered.

The fix is to smooth over a physical time window before taking the maximum. The window has to be
long against the local mode and short against the loading, which at 7.6 mm/s traverses the whole
ascending branch in ~1.5 s: 1 ms is three decades clear of both. The result is insensitive to the
choice — on the eps_su = 0.025 cell the smoothed peak moves 0.5% across a 10x range of window
widths (966.8 -> 956.2 kN over 0.5 -> 5 ms) while the raw maximum reads 996.6.

Nothing here re-runs anything: every path-following runner already stores `shear` per step, so an
existing run is rescored from its own `data.json` (`rescore.py`).
"""
from __future__ import annotations

import numpy as np

WINDOW_S = 1.0e-3          # the default smoothing window; see the module docstring
RINGING_FLAG = 0.02        # raw/smoothed above this and the raw peak is reporting ringing
DROP = 0.20                # drift capacity = where the smoothed response falls to 80% of peak
PLATEAU = 0.005            # the peak is "flat" over the drift range holding within 0.5% of it


def moving_average(y, w: int):
    """Centred moving average of `y` over `w` samples, correct at the ends.

    Each output is the mean of the samples that ACTUALLY fall in the window, so the first and last
    half-windows are not biased towards zero the way a zero-padded convolution makes them.
    """
    y = np.asarray(y, float)
    n = len(y)
    w = int(max(1, min(w, n)))
    if w == 1:
        return y.copy()
    c = np.concatenate(([0.0], np.cumsum(y)))
    i = np.arange(n)
    lo = np.maximum(0, i - w // 2)
    hi = np.minimum(n, lo + w)
    lo = np.maximum(0, hi - w)
    return (c[hi] - c[lo]) / (hi - lo)


def response_metrics(shear, disp, *, dt: float, height: float, window_s: float = WINDOW_S,
                     drop: float = DROP, plateau: float = PLATEAU) -> dict:
    """Peak base shear over a `window_s` moving average, with the raw maximum for comparison.

    Takes `max`, not `max(abs)`, to mirror what `run.py` records — on a monotonic pushover the two
    agree, and on a cyclic run both refer to the positive side.

    `ringing_ratio` is raw/smoothed: 1.00 means the maximum sits on a smooth curve and the raw
    number is fine; anything above `RINGING_FLAG` means the raw peak is a crest of a local release
    mode and should not be quoted as resistance.
    """
    shear = np.asarray(shear, float)
    disp = np.asarray(disp, float)
    if shear.size == 0 or not np.isfinite(dt) or dt <= 0:
        return {}
    w = int(round(window_s / dt))
    sm = moving_average(shear, w)
    drift = disp / height
    i_sm = int(np.argmax(sm))
    i_raw = int(np.argmax(shear))
    peak_sm, peak_raw = float(sm[i_sm]), float(shear[i_raw])
    out = {"peak_shear_smooth": peak_sm,
           "drift_at_peak_smooth": float(drift[i_sm]),
           "peak_smoothing_window_s": window_s,
           "peak_smoothing_samples": w,
           "peak_ringing_ratio": (peak_raw / peak_sm) if peak_sm else float("nan"),
           "peak_shear_raw": peak_raw,
           "drift_at_peak_raw": float(drift[i_raw])}

    # HOW SHARP IS THE PEAK? On these walls the top is often a PLATEAU — the no-rupture cell holds
    # within 0.5% of its maximum from 0.475% to 0.703% drift — and then "drift at peak" is an
    # argmax over a flat region, not a measurement. A cell truncated by bar rupture instead peaks
    # inside 0.0001% drift. Reporting the plateau alongside stops the first being read like the
    # second.
    on_top = np.flatnonzero(sm >= (1.0 - plateau) * peak_sm)
    if on_top.size:
        out["peak_plateau"] = [float(drift[on_top[0]]), float(drift[on_top[-1]])]
        out["peak_plateau_fraction"] = plateau

    # DRIFT CAPACITY: where the smoothed response first falls to (1-`drop`) of its peak AFTER the
    # peak — the same 20%-drop convention the test literature uses (and which WSH3 never met, so
    # its 2.04% is a lower bound). Measured on the smoothed curve because the raw one crosses the
    # threshold repeatedly while ringing.
    below = np.flatnonzero(sm[i_sm:] < (1.0 - drop) * peak_sm)
    out["capacity_drop"] = drop
    out["drift_capacity"] = float(drift[i_sm + below[0]]) if below.size else None
    if not below.size:
        out["drift_capacity_note"] = (f"never fell to {1 - drop:.0%} of peak — the run ended at "
                                      f"{float(drift[-1]):.4%} drift still above it, so this is a "
                                      f"LOWER BOUND, not a capacity")
    return out


def smoothed_peak(*a, **k) -> dict:
    """Backwards-compatible alias for `response_metrics`."""
    return response_metrics(*a, **k)


def headline_peak(data: dict) -> tuple[float, float, bool]:
    """(peak, drift_at_peak, smoothed?) — the number to QUOTE from a run record.

    Prefers the smoothed peak when the run has one, so a rescored old run and a new run read the
    same way, and falls back to the raw maximum for a record written before D92.
    """
    if data.get("peak_shear_smooth") is not None:
        return float(data["peak_shear_smooth"]), float(data.get("drift_at_peak_smooth", 0.0)), True
    return float(data.get("peak_shear") or 0.0), float(data.get("drift_at_peak") or 0.0), False


def ringing_note(data: dict) -> str | None:
    """One line naming the correction, or None when the raw maximum was already clean."""
    r = data.get("peak_ringing_ratio")
    if r is None or not np.isfinite(r):
        return None
    raw = data.get("peak_shear_raw", data.get("peak_shear"))
    sm = data.get("peak_shear_smooth")
    if r - 1.0 <= RINGING_FLAG:
        return (f"raw maximum {raw / 1e3:,.1f} kN is {r:.3f}x the smoothed peak — no correction "
                f"needed, the maximum sits on a smooth curve")
    return (f"RAW MAXIMUM {raw / 1e3:,.1f} kN IS RINGING, {r:.3f}x the smoothed {sm / 1e3:,.1f} kN. "
            f"A released bar or strut rings locally far below T1 and the sample maximum lands on a "
            f"crest of it; quote the smoothed value as resistance (D92)")
