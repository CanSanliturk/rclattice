# TODO — RC shear wall (SW-NC-FF) study

Open items as of 2026-08-16. The elastic calibration (D47) and the cyclic machinery are working;
what remains is documentation, tests, and one unresolved question about the specimen itself.

## Blocking a validation claim

- [ ] **Resolve the ρ_b discrepancy.** Table 1 reports ρ_b = 2.37 % for SW-NC-FF, but 4Ø14 over a
      10 % boundary gives 615.75 / (0.10·1000·200) = **3.08 %**. The identical calculation is exact
      for SW-IC-FF (1.539 % vs 1.54 %), so it is not a definitional difference. Working backwards,
      2.37 % implies ~474 mm² (≈3Ø14) or a ~130 mm boundary length. Check Section 2-2 for the bar
      count. This is ±10 % on boundary steel, feeding straight into peak strength — the number the
      whole validation rests on.

## Analyses

- [ ] **Full 4 % cyclic protocol** — running (`/tmp/cyc4.log`, ~2.9 h, 463 534 steps).
- [ ] **`--compression elastic` at 1 %** — prices Aydin's uncalibrated-compression assumption at
      drifts where crushing matters. Earlier runs stopped short of crushing, so the ~1 % cost
      measured then understates it.
- [ ] Re-run the 1 % case if full-resolution loops are wanted; the existing data is log-sampled
      (434 of 86 913 points), from before the persistence fix.

## Documentation (DECISIONS.md)

- [ ] **D48 — `run_cyclic` + `cyclic_protocol`.** Reversed-cyclic displacement control; stalls near
      0.3 % drift on a cracking lattice, which is why D49 exists.
- [ ] **D49 — `run_cyclic_dynamic` and two solver bugs.** The most transferable lesson of the day:
      - `ops.wipeAnalysis()` is REQUIRED before setting a transient integrator after a gravity
        stage. Without it OpenSees refuses the integrator and silently falls back to a default —
        no error, just a different analysis than the one written.
      - Rayleigh damping must use `betaKinit` (3rd arg), never `betaK` (2nd). The `betaK` term
        rides the committed tangent, so when struts crack the damping force explodes; symptom is
        `Norm deltaR ~1e7` appearing exactly as cracking spreads.
      - D33's answer for the seismic runner (modal damping) does NOT transfer here: dynamic
        relaxation needs near-critical damping, and `modalDamping(0.8)` fails to converge from the
        first step. Tried and rejected — record it so it is not retried.
      - Same two bugs were also present in `run_pushover_dynamic` and are now fixed there.
      - Parameterize the drive by `rate`, not `periods_to_*`: the latter scales speed with
        amplitude, so the same setting drives a larger protocol proportionally faster and
        contaminates the recorded base shear (reactions include inertial and C·v terms).
        7.6 mm/s reproduces the static solver within 1 % for this wall.
- [ ] **D50 — specimen corrections read off the paper figures.** Each changed the model:
      - SW-NC-FF has **no boundary confinement** (Table 1, ρ_con = "NC"). The Ø8@200 hoops
        initially added belong to SW-IC-FF; including them models the wrong wall and improves
        exactly the behaviour the test was designed to expose.
      - Table 1's "closed hoop / 90° hook" columns describe the HORIZONTAL WEB reinforcement's
        anchorage, not the boundary.
      - Longitudinal bars are **hooked into the pedestal** (Fig. 3a/b), not stopped at the interface.
      - Steel hardening b is **measurable** from Table 2 (ε_h, f_max, ε_max): 0.0019–0.0038, not
        the assumed 0.01, which overstated post-yield strength ~4×.
      - Ø16@100 loading head over the top 500 mm (Fig. 3c).
- [ ] **D51 — digitizing the measured hysteresis** (`digitize.py`): frame-based axis calibration,
      colour separation of the test curve, self-check against the reported peaks (±1 %). Note the
      limits: a point cloud, not an ordered path, so no per-cycle energy; accuracy capped at line
      thickness (~3 kN).
- [ ] **D52 — overlaying the model on the digitized loops** (`compare.py`, builds on D51). Both
      series share the 2200 mm shear span, so (drift %, kN) overlays with no rescaling. The test
      BACKBONE is recovered from the unordered cloud as a **loop-tip envelope**: at each protocol
      drift level, the extreme load over the slice [(1−0.10)·d, (1+0.05)·d]. The slice must reach
      slightly INSIDE the extreme because a degrading wall peaks just before its peak displacement;
      0.10 is the widest slice that still excludes LARGER loops' reloading branches (at 0.15 the
      3.0 % level jumps 152 → 174 kN, which is the 4 % loop). Recovers +209.9/−211.8 kN against the
      reported +212.6/−209.4. Result: model/test peak = **1.08 push, 1.05 pull**, but the model
      holds ~190 kN at 4 % drift where the test has degraded to ~142 kN — post-peak degradation is
      the visible gap, and it is downstream of the omitted debonding, so it is not a tuning target.
- [ ] Also record the earlier wrong turns so they are not repeated: horizon 3.01 made convergence
      WORSE (0.160 % vs 0.375 %) because length regularization makes longer struts more brittle;
      boundary ties at 200 mm did not help because they brace only every 4th node row at a 50 mm mesh.

## Tests

- [ ] `run_cyclic` — protocol expansion, reversal handling, partial-result reporting on failure.
- [ ] `run_cyclic_dynamic` — that the integrator is actually honoured (i.e. `wipeAnalysis` is
      present), and that a small-amplitude dynamic cycle reproduces the static one within a few
      percent. That equivalence is the whole basis for trusting the dynamic result.
- [ ] `transformed_inertia` with hooked bar paths — this regressed silently once (hooked paths made
      every longitudinal bar look horizontal, so I_tr collapsed to I_gross).

## Docs

- [ ] `CLAUDE.md` — status section, the new `examples/wall/` files (`pushover.py`, `cyclic.py`,
      `draw.py`, `digitize.py`, `replot.py`, `data/`), and the two new runners in `opensees.py`.

## Deferred / considered and not done

- **Aydin's own tension calibration (his §2.4).** His multilinear law needs a₁,a₂,a₃ fitted for OUR
  concrete via a direct-tension coupon; his Table 1 values are per-specimen. `HystereticSM` is the
  only available material that carries his exact 4-point envelope, and is the right vehicle for a
  MONOTONIC run (his SLA is monotonic anyway). Not for cyclic: it puts up to 50×F_cr of spurious
  compression into a cracked strut at closure; `-beta 2` cuts that by two orders of magnitude but
  still leaves ~40 % F_cr, where Concrete02 is exactly zero.
- **Bond-slip / rocking.** The specimen's dominant mechanism (74 % of displacement at 2 % drift) is
  rocking from plain-bar debonding. A perfect-bond lattice cannot represent it, so drift capacity,
  pinching, self-centering and energy dissipation are permanently out of reach without a bond
  interface. Strength and stiffness are not.
