"""Spike-safe, structured summary of a base-shear response (D70).

WHY THIS EXISTS. Every runner reported `peak = max(shear)` over the whole history and every script
printed it as though it were a capacity. On a run that collapses or diverges, the largest recorded
shear is a SPIKE inside the wreckage, not a strength — and it was quoted as a result three separate
times in one investigation:

  * Newmark pushover: "peak 967.9 kN" — 1.14x the model's own 848 kN section capacity, recorded
    while the base shear was oscillating +-1300 kN;
  * corotTruss pushover: "peak 774.3 kN at 0.884% drift, model/test = 0.869" — ONE sample out of
    424 in that drift window, whose median was 24.9 kN and whose minimum was -1010 kN. The pier had
    collapsed 0.55% of drift earlier;
  * the damage map captured at `disps_peak`, which for the same reason snapshotted the divergence
    spike and drew the aftermath rather than the cause.

None of those were subtle once looked at. They were missed because the number arrived pre-computed
and labelled "peak", so it was read rather than checked. The fix is to make the summary itself carry
the checks, and to make the RELIABLE number the one that gets printed first.

THREE INDEPENDENT TESTS, all cheap:
  1. COLLAPSE — the first point where the running maximum is lost by more than `drop` and never
     recovered. Everything after it is post-collapse and is reported separately, never as strength.
  2. SPIKE — is the peak an isolated sample? Compared against the MEDIAN of its own drift
     neighbourhood, a genuine peak sits on a smooth branch (ratio ~1) and a spike towers over it.
  3. CAPACITY — an optional physical bound (here the hand section analysis, 848 kN). A pushover
     cannot exceed its own section capacity; if it does, the number is arithmetic, not mechanics.

A CYCLIC response is summarised on its ENVELOPE (peak shear at each new drift extreme), because the
raw signal oscillates by design and every reversal would read as a collapse.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass
class Peak:
    kN: float
    drift_pct: float
    index: int
    is_spike: bool = False
    prominence: float = 1.0          # peak / median of its drift neighbourhood; ~1 = smooth branch


@dataclass
class Summary:
    reported: Peak                   # naive max(shear) — what the old code printed
    reliable: Peak                   # peak of the pre-collapse branch: the number to quote
    collapse_drift_pct: float | None
    collapse_drop: float
    capacity_kN: float | None
    exceeds_capacity: bool
    n_points: int
    cyclic: bool
    warnings: list = field(default_factory=list)

    @property
    def trustworthy(self) -> bool:
        return not self.warnings

    def to_dict(self) -> dict:
        d = asdict(self)
        d["trustworthy"] = self.trustworthy
        return d

    def report(self) -> str:
        L = []
        r = self.reliable
        L.append(f"  RELIABLE peak         {r.kN:8.1f} kN at {r.drift_pct:.3f}% drift"
                 f"   <-- quote this one")
        if self.reported.index != r.index:
            p = self.reported
            L.append(f"  raw max(shear)        {p.kN:8.1f} kN at {p.drift_pct:.3f}% drift"
                     f"   [NOT a capacity: {'isolated spike, ' if p.is_spike else ''}"
                     f"prominence {p.prominence:.1f}x its neighbourhood median]")
        if self.collapse_drift_pct is not None:
            L.append(f"  COLLAPSE at           {self.collapse_drift_pct:8.3f}% drift"
                     f"   ({self.collapse_drop:.0%} of the running maximum lost and not recovered)")
            L.append(f"                        everything beyond this is post-collapse and is NOT "
                     f"a strength result")
        if self.capacity_kN:
            L.append(f"  section capacity      {self.capacity_kN:8.1f} kN"
                     f"   -> reliable/capacity = {r.kN / self.capacity_kN:.3f}")
        for w in self.warnings:
            L.append(f"  ** {w}")
        return "\n".join(L)


def envelope(drift, shear):
    """`(drift, shear)` at each new POSITIVE drift extreme — the cyclic counterpart of a backbone."""
    out, hi = [], -1e30
    for i, (u, s) in enumerate(zip(drift, shear)):
        if u > hi:
            hi = u
            out.append((i, u, s))
    return out


def _spike(drift, shear, i, halfwidth=0.03, frac=2.0):
    """Is sample `i` an isolated spike? Returns `(is_spike, prominence)`.

    Compared against the MEDIAN of its own drift neighbourhood: a peak on a smooth branch has
    prominence ~1, a spike towers over a neighbourhood that is doing something else entirely.
    """
    u0 = drift[i]
    win = [s for u, s in zip(drift, shear) if abs(u - u0) <= halfwidth]
    if len(win) < 5:
        return False, 1.0
    med = sorted(win)[len(win) // 2]
    if abs(med) < 1e-9:
        return True, float("inf")
    prom = abs(shear[i] / med)
    return prom >= frac, prom


def _collapse(drift, shear, drop=0.40, recover=0.90, sustain=5):
    """First index where the running max is lost by `drop` and never SUSTAINABLY recovers.

    "Sustainably" is the whole point. An earlier version tested recovery with `max(shear[i:])`, and a
    single post-collapse spike satisfied it — so the detector walked straight past the real collapse
    at 0.338% drift and reported the spike at 0.884% as the reliable peak, which is the exact bug
    this module exists to prevent. Recovery now needs `sustain` samples above the threshold, which a
    lone spike cannot provide.

    The running maximum is also taken over a SPIKE-RESISTANT signal (`_median3`), so a single sample
    cannot raise the bar it is later judged against.
    """
    med = _median3(shear)
    run_max = -1e30
    for i, s in enumerate(med):
        if s > run_max:
            run_max = s
        elif run_max > 0 and s < (1.0 - drop) * run_max:
            tail = med[i:]
            if sum(1 for t in tail if t >= recover * run_max) < sustain:
                return i, 1.0 - s / run_max
    return None, 0.0


def _median3(x, w=2):
    """Running median over +-w samples — removes isolated spikes without shifting a smooth branch."""
    n = len(x)
    return [sorted(x[max(0, i - w):min(n, i + w + 1)])[
                (min(n, i + w + 1) - max(0, i - w)) // 2] for i in range(n)]


def summarize(drift, shear, *, capacity_kN: float | None = None, cyclic: bool = False,
              drop: float = 0.40) -> Summary:
    """Structured summary with the spike / collapse / capacity checks applied."""
    if cyclic:
        env = envelope(drift, shear)
        idx = [i for i, _u, _s in env]
        d = [u for _i, u, _s in env]
        v = [s for _i, _u, s in env]
    else:
        idx, d, v = list(range(len(drift))), list(drift), list(shear)
    if not v:
        raise ValueError("empty response")

    k = max(range(len(v)), key=lambda i: v[i])
    is_spike, prom = _spike(d, v, k)
    reported = Peak(v[k], d[k], idx[k], is_spike, prom)

    ci, cdrop = _collapse(d, v, drop=drop)
    if ci is not None:
        pre = v[:ci] or v[:1]
        collapse_drift = d[ci]
    else:
        pre, collapse_drift = v, None
    # search the pre-collapse branch on the SPIKE-RESISTANT signal, so an isolated sample inside an
    # otherwise valid branch cannot become the quoted capacity either
    pm = _median3(pre)
    j = max(range(len(pm)), key=lambda i: pm[i])
    reliable = Peak(pre[j], d[j], idx[j], *_spike(d[:len(pre)], pre, j))

    warn = []
    # only complain if the raw peak is MATERIALLY above the reliable one; a peak landing one sample
    # either side of the collapse boundary is the same peak, not a misreading
    if reported.kN > 1.02 * reliable.kN:
        warn.append(f"raw max(shear) = {reported.kN:.1f} kN is POST-COLLAPSE — do not quote it")
    if reported.is_spike:
        warn.append(f"raw max(shear) is an isolated spike ({reported.prominence:.1f}x its "
                    f"neighbourhood median)")
    over = bool(capacity_kN and reliable.kN > capacity_kN)
    if over:
        warn.append(f"reliable peak {reliable.kN:.1f} kN EXCEEDS the section capacity "
                    f"{capacity_kN:.1f} kN — that is arithmetic, not mechanics")
    return Summary(reported, reliable, collapse_drift, cdrop, capacity_kN, over,
                   len(v), cyclic, warn)
