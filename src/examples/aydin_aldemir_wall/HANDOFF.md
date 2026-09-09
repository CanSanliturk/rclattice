# Aldemir wall — handoff

Written 2026-08-29 to start a fresh session with the state of this study. Read
**DECISIONS.md D72–D75** for the full reasoning; this is the map.

## The specimen

Aldemir, Binici & Canbay (2017), *Cyclic testing of reinforced concrete double walls*,
ACI Struct. J. 114(2) 395–406 — **not in the repo**. Everything is second-hand via
Aydin, Tuncay & Binici (2019), *Simulation of RC Member Response Using Lattice Model*,
J. Struct. Eng. 145(9) 04019091, whose PDF sits in this directory. That paper's Table 1/2/4
and Fig. 4(f)/10 are the only sources.

**Geometry: 3000 wide × 2250 high × 210 thick.** Squattest wall in the repo (ratio 0.75),
**no axial load**, Ø8 @ 100 both ways ×3 bars (ρ = 0.718%), f_c 28 / f_t 1.85 (TS 500) /
E_c 24,870 (ACI) / f_y 360. Mesh 50 divides everything.

Measured: **K = 1,038.44 kN/mm, F = 963.592 kN at ~1% drift, no strength degradation.**
Aydin's own model: K_sim 943.16, F_sim 1,164.4 (**1.21× the test** — the worst of his six specimens).

## Two results, keep them apart

**SETTLED — the thickness was 120 and should be 210 (D73).** Five independent checks, none of
them a fit to the collapse. The `replica/` package reproduces his Table 2 counts *exactly*
(20,385 nodes / 80,684 struts) and then misses his K_sim by 1.743; stiffness is exactly linear
in thickness, so that factor *is* the thickness. 210 = 50+100+60 of his own plan stack. At 210
the replica returns 0.979 of his K_sim unfitted, the uncracked section finally sits **above**
the measured stiffness (1.31×, as it must — at 120 it was 0.75×, impossible), and V_flex/F_exp
= 1.00. Every run before 2026-08-29 was 43% too thin.

**OPEN — a premature collapse.** The wall loses its load path at **0.237% drift** against the
test's ~1%. Ten experiments have failed to remove it.

## What moved it, and what didn't

| | peak | drift capacity |
|---|---|---|
| thickness 120 → 210 (D73) | ×1.53 | ×1.9 |
| bond elements (D75) | ×1.14 | ×1.21 |
| **the other eight** | — | **timing only** |

Damping (0.2/0.5), fracture energy (×1/×2), crack-release ringing, lattice topology, mesh
(50/25), constitutive law (Concrete02 vs his trilinear), reinforcement presence, and
**integration scheme** all changed *when* it failed, never *that* it did.

Best model to date — mesh 50, t = 210, bond, explicit:
**peak 1,026.6 kN = 1.065× measured**, collapse 0.2374%, ascending-branch residual 1.51%.
Closer to the test than Aydin's own 1.21× overshoot.

## Capabilities built here (library, not study)

- **`field="equibiaxial"`** on the energy balance (D72) — his published Appendix route;
  `aydin_closed_form_C` reproduces his C = 0.621 / 0.102. **The repo default `"uniaxial"` is
  measurably better here** (0.92 vs 0.79 against a continuum), so it stays.
- **Optional bond** (D72/D75) — `--bond --bond-area-ratio 0.01`. Force–slip law.
- **EXPLICIT integration** (D74) — `--explicit [CentralDifference]`. **45× cheaper per step,
  2.7× faster end to end.** `builders.critical_time_step()` sizes it.

## Traps that cost real time

1. **`betaKinit` damping is fatal under an explicit integrator** — shrinks the stable step ~90×.
   Correct implicitly (D49), wrong explicitly. Fails as an 18-hour run, not an error.
2. **Explicit `dt` is set by `2/w_max` (stiffest, lightest element), NOT by T1.** `steps_per_period`
   sized off T1 is 15–20× too large.
3. **Bond: area is the wrong single knob.** A truss couples `k = EA/L` and `F = ft·A`. Calibrating
   area for stiffness destroys strength (collapse at 0.026%, *earlier* than no bond at all).
4. **A bond link must not inherit concrete's cracking strain** — that is 0.004 mm of slip. Real
   bond fails at 0.1–1 mm.
5. **`figure_damage`'s `eps_crush` wants a POSITIVE magnitude.** Passing `-EPSC0` flags ~99% of
   struts as crushed; the tell is a `< --2.3e-03` legend.
6. **Gf is not a neutral knob** — ×2 changes base shear 14.2%. That caveat also lands on D67's WSH3.
7. **His trilinear law is unusable with Newton** — path-independent `ElasticMultiLinear` makes
   struts flip branches; 1,510 solver failures vs 19. He runs it explicitly for this reason.
8. **PROCESS.** Three conclusions were reported from samples taken *before* the event being watched
   for, on a collapse already measured at ~30 steps against 250-step sampling. **Report what a run
   has shown, not what it implies about a drift it has not reached.**

## Fig. 10(b) — RECOVERED (D77, supersedes "not recoverable")

The earlier note said the test backbone could not be digitized. That was over-general and is
corrected: **`digitize.py` recovers the figure**, writing `data/fig10b.npz`.

What is genuinely lost is **order** — the experiment is drawn as discrete dots, so per-cycle energy,
degradation and the load path are gone. What is not lost is the **outline**: extremes survive
unordered. Also recovered are **both of Aydin's own curves**, which are clean coloured lines.

Calibration uses gridline geometry only, then reproduces four numbers it never saw:

| | digitized | published | ratio |
|---|---|---|---|
| Aydin horizon 1.5 plateau | 1158.0 kN | 1164.413 | 0.9945 |
| Aydin horizon 3.01 plateau | 1305.3 kN | 1325.675 | 0.9846 |
| cloud envelope peak | 969.4 kN | 963.592 | 1.0060 |
| "experiment" marker | 924.2 kN | 963.592 | 0.959 (coarse) |

**The envelope is an UPPER BOUND on the backbone, not the backbone** — a point on it may be a loop
tip or the unloading side of a larger cycle. The figure's legend also occludes real data at positive
displacement below about -500 kN.

`replica/compare.py` plots every run against these.

## What is left

- The drift gap is 3× and the remaining unknowns are all **bond parameters he never prints**
  (area, slip, strength — only the residual ratio a = 0.7 is published), plus b1/b2, his PID
  drive, and dt = 5e-8 s.
- `replica/` is a calibrated instrument (his counts exactly, his K_sim to 2%) and runs under
  `--explicit`. Without bond it reaches **0.42 of his published peak** and collapses at 0.098%.
  **Bond is now wired into it (D76) and NOT YET RUN** — `run.py --bond`, `preflight.py --bond`,
  same force-slip law as the parent. 20,385 nodes / 88,324 elements → **28,081 / 149,802**; his
  Table 2 concrete counts are untouched. Run in this order: `preflight.py --bond --explicit` to
  price it, `--elastic --bond` against `--elastic` to size the D72 artefact and check the area
  ratio (0.01 is the parent's, calibrated at mesh **50**; at mesh 20 the failure slip is 0.149 mm,
  on the FLOOR of the Model Code band, and 0.003 went singular at mesh 50 — the window is narrow
  and unmeasured here), then the pushover.
- A coincident-node `zeroLength` bond spring would fix the topology properly; unbuilt.
- Cyclic stage never started. `ElasticMultiLinear` is unsuitable for cyclic regardless (D60).

## Files

`specimen.py` (ALDEMIR_TW overrides thickness) · `testdata.py` (every published number + what is
absent) · `build.py` · `summary.py` (run first — adjudicates the geometry) · `elastic.py` ·
`pushover.py` (`--solver`, `--explicit`, `--bond`, `--material`, `--no-rebar`, `--groups`) ·
`preflight.py` · `draw.py` · `replica/`.
Runs under `examples/output/aydin_aldemir_wall/runs/`, each with `command.txt`; aborted ones carry
an `ABORTED.md` saying what they established.
