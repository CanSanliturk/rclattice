"""Dynamic base excitation + intensity tuning for the frame seismic studies (shared by dynamic.py
and dynamic_linear.py). Identical to the column package's excitation module — each example package
keeps its own self-contained copy.

Provides the two excitation sources — the recorded El Centro 1940 NS (`load_elcentro`) and a harmonic
record (`sine_excitation`), unified by `make_excitation` — plus `tune_intensity`, which scales the
chosen record to a target peak roof drift on the cheap fiber-frame reference. Both records are
returned in g, so the same `scale = G * intensity` machinery in the runners applies to either.
Units: kip, in.
"""

from __future__ import annotations

import math
from pathlib import Path

G = 386.4                                   # gravity (in/s^2) — kip-in units; converts record g -> in/s^2
RHO = 2.25e-7                               # concrete mass density (kip-s^2/in^4), = column self-mass source
RECORD = Path(__file__).resolve().parent.parent / "data" / "elcentro.at2"
N_CYCLES = 20                               # sine excitation: number of cycles (duration = N*T1)
SAMPLES_PER_CYCLE = 100                     # sine excitation: record sampling resolution


def load_elcentro() -> tuple[list[float], float]:
    """Parse the PEER .at2 record -> (acceleration values in g, dt). Header/NPTS lines are skipped
    (they fail float parsing); dt is read from the 'DT= .02000 SEC' line."""
    accel: list[float] = []
    dt = 0.02
    for line in RECORD.read_text().splitlines():
        if "DT=" in line:
            dt = float(line.split("DT=")[1].split()[0])
        try:
            row = [float(t) for t in line.split()]   # whole row first: a partial row must not leak
        except ValueError:
            continue                                  # header / NPTS line
        accel.extend(row)
    return accel, dt


def sine_excitation(period: float, *, n_cycles: int = N_CYCLES, amplitude: float = 1.0,
                    samples_per_cycle: int = SAMPLES_PER_CYCLE) -> tuple[list[float], float]:
    """Harmonic base-acceleration record a_g(t) = amplitude * sin(2*pi*t/period), in g.

    Returned in g (the same unit as load_elcentro), so the SAME `scale = G * intensity` machinery
    converts it to in/s^2 and applies the tuned intensity. `n_cycles` sets the duration
    (n_cycles*period); `samples_per_cycle` the resolution (dt = period/samples_per_cycle). Starting
    at sin(0)=0 there is no initial acceleration jump, so no start-up window is needed. Returns
    (accel_g, dt)."""
    dt = period / samples_per_cycle
    npts = n_cycles * samples_per_cycle + 1
    accel = [amplitude * math.sin(2.0 * math.pi * k / samples_per_cycle) for k in range(npts)]
    return accel, dt


def make_excitation(kind: str, period: float, n_cycles: int = N_CYCLES) -> tuple[list[float], float, str]:
    """Build the dynamic base input. `kind="sine"` -> a harmonic acceleration RESONANT with `period`
    (the structure's fundamental T1), `n_cycles` long; `kind="elcentro"` -> the recorded El Centro
    1940 NS. Both are returned as (accel in g, dt_record, label) for the UniformExcitation."""
    if kind == "sine":
        accel, dt = sine_excitation(period, n_cycles=n_cycles)
        return accel, dt, f"sine resonance (T={period:.2f}s, {n_cycles} cyc)"
    accel, dt = load_elcentro()
    return accel, dt, "El Centro 1940 NS"


def tune_intensity(run_bc, target_drift: float, *, linear: bool, max_iter: int = 10,
                   cap_factor: float = 5.0) -> tuple[float, dict]:
    """Find the excitation intensity whose peak roof drift ~ `target_drift`, on the cheap fiber
    column (`run_bc(intensity) -> result dict with 'peak_disp', 'converged', 'periods'`).

    `linear=True`: the response is proportional to intensity, so a single trial at 1.0x plus one
    exact correction lands the target.

    `linear=False`: under a RESONANT sine a softening column ratchets to collapse past a threshold,
    so peak drift is non-monotonic / unbounded in intensity and the naive proportional update
    diverges (collapsing the scale toward 0 = elastic). Use a divergence-guarded bisection instead:
    treat a run that does not converge, or whose drift exceeds `cap_factor*target_drift` or simply
    exceeds the target, as "too strong" (upper bracket); keep the largest intensity that stays
    converged and at/below target as the answer. If the target sits past the largest stable drift,
    this returns the closest stable (still nonlinear) response rather than chasing collapse.
    Returns (intensity, bc_result_at_that_intensity)."""
    bc = run_bc(1.0)
    peak = abs(bc["peak_disp"])
    print(f"  [1.000x] BC peak drift = {peak:.3f} in  (T1_bc={bc['periods'][0]:.3f}s)")
    if linear:
        used = target_drift / peak if peak > 1e-9 else 1.0
        return used, run_bc(used)

    cap = cap_factor * target_drift
    lo_i, hi_i = 0.0, None          # bracket: lo = largest OK (<=target) intensity, hi = smallest too-strong
    best_i, best, cur_i = 0.0, None, 1.0
    for _ in range(max_iter):
        too_strong = (not bc["converged"]) or peak > cap or peak > target_drift
        if too_strong:
            hi_i = cur_i if hi_i is None else min(hi_i, cur_i)
        elif cur_i > lo_i:
            lo_i, best_i, best = cur_i, cur_i, bc
        if best is not None and abs(abs(best["peak_disp"]) - target_drift) < 0.1 * target_drift:
            break
        if hi_i is None:            # never reached target -> grow
            cur_i *= 2.0
        elif lo_i == 0.0:           # only too-strong samples so far -> shrink below hi
            cur_i = hi_i * 0.5
        else:                       # bracketed -> bisect
            cur_i = 0.5 * (lo_i + hi_i)
        bc = run_bc(cur_i)
        peak = abs(bc["peak_disp"])
        print(f"  [{cur_i:.3f}x] BC peak drift = {peak:.3f} in  (conv={bc['converged']})")
    if best is None:                # nothing stayed bounded -> return the last (smallest) trial
        return cur_i, bc
    return best_i, best
