# Aldemir wall — parametric study plan

Written 2026-09-05. The programme is the author's (Aydin, first author of the 2019 paper and the
supervisor of this work), set on 2026-09-05 and refined by five follow-up answers the same day:

> * Pushover with the old method. Linear compression, tension softening calibration, run.
> * Compression-capable lattice pushover. Compression calibration after tension softening
>   calibration, run.
> * Cyclic response.
> * Don't forget bond, and regularization.

This document is the plan, not the implementation. Nothing here is built yet. Read
**DECISIONS.md D72–D81** for how the study reached this point; D78 (the collapse was a
load-introduction artefact), D80 (seven corrections to the record) and D81 (this programme) are the
load-bearing ones.

**Amendments of 2026-09-05, from a code inspection of the claims below**, are marked in place in
§2, §3, §7 and §8. Four of them tighten wording; the §3 one changes what Stages 5–6 can actually run.

---

## 1. What the study delivers

Two concrete models of the same wall, each with and without bond, each pushed and cycled, all
driven from one parameterized entry point, all reported against the same references. The
deliverable is a comparison the paper itself does not contain: what the linear-compression
assumption costs, what bond buys, and whether either explains the +21% over-prediction that makes
this the worst-predicted specimen of the paper's six.

**Fair claims from this study:** peak strength, elastic stiffness, damage location, the load-path
split between plane-sections and truss action, and the *controlled effect* of each parameter, since
every cell differs from its neighbour in exactly one thing.

**Not fair, and to be stated in every report:** drift capacity, because the loading protocol is
invented (§6); cyclic energy and per-cycle degradation, because the digitized test data is an
unordered point cloud; anything resting on the section, because the 210 mm thickness is inferred
(D73) and the 2017 test paper is not available.

---

## 2. Fixed inputs

The specimen and everything the matrix does not vary. Sources are the 2019 paper's Table 1, Table 4
and Fig. 10; the 2017 test paper is **not** in the repo, so every measured value is second-hand.

| Quantity | Value | Source |
|---|---|---|
| panel | 3000 × 2250 × **210** mm, aspect 0.75, no axial load | Fig. 10(a); thickness inferred, D73 |
| reinforcement | Ø8 @ 100 both ways, 3 bars/position, 150.8 mm²/line, ρ = 0.718% | Fig. 10(a) |
| cover | 50 vertical, 75 horizontal (snapped to 50 on the grid) | measured off Fig. 10(a), D72 |
| concrete | f_c 28, f_t 1.85 (TS 500), E_c 24,870 (ACI), G_f 0.075 N/mm | Table 1 |
| ε_c0 | 0.002252, derived as 2f_c/E_c | D56 |
| steel | f_y 360, E_s 200,000, elastic–perfectly plastic | Table 1, Fig. 1(d) |
| **measured** | K = 1038.44 kN/mm, F = 963.592 kN at ~1% drift, no degradation | Table 4 |
| **Aydin 1.5d** | K = 943.16 kN/mm, F = 1164.413 kN | Table 4 |
| **Aydin 3.01d** | K = 1052.37 kN/mm, F = 1325.675 kN | Table 4 |

Fixed model settings, all justified elsewhere and none of them a matrix axis:

| Setting | Value | Why |
|---|---|---|
| mesh | 50 mm | divides 3000, 2250, 100 and 50 exactly; the paper's own 20 does not divide 2250 |
| horizon | 1.5 | his own recommendation for RC; 3.01 is an extension slot |
| elastic calibration | uniaxial energy balance at ν = 0.20 → **A_t = 6324.7 mm², EA = 1.573e8 N** | D47; measurably better than the published equibiaxial route on this specimen (0.92 against 0.79) |
| longitudinal bars | run to the top of the wall | D78 — without this the wall tears at the driven row |
| strut element | `corotTruss` | large-displacement consistency |
| integrator | explicit `CentralDifference`, dt from `critical_time_step` | D74 — 45× cheaper per step, and the only scheme his path-independent law tolerates |
| drive | 7.6 mm/s, damping 0.5 (mass-proportional under explicit) | the repo's cross-study values; the licence is the measured residual, reported per run |
| model size | 2,806 nodes / 13,471 elements; **5,424 / 34,325 with bond** | disputed, see below |

**Amendment — SETTLED by a build, 2026-09-05: the table above is right and D72 was wrong.** At
mesh 50, horizon 1.5, t = 210 the model builds as **2,806 nodes / 13,471 elements** (10,905
concrete + 1,290 longitudinal + 1,276 stirrup) and **5,424 / 34,325 with bond** (+20,854 bond
links). D72 and CLAUDE.md gave 5,454 / 13,501 / 34,595 — 30 / 30 / 270 higher, on a configuration the
current code does not build, and **the provenance of those numbers is unexplained**. "Bars run to
the top" was the obvious candidate and is refuted by measurement: that option costs **+60**
elements, not +30, because the top bar line sits at y = 2150, i.e. TWO mesh-50 cells below the top
face, so each of the 30 bars gains two segments (and with bond, +60 steel nodes and +450 elements).
A +30 delta would need the bars to stop ONE cell short, i.e. a different cover than the one now
modelled. CLAUDE.md corrected to the measured numbers. Counts are independent of `--tw`.

### Geometry is a parameter, and the drawn panel is its default

The panel is `--panel`, not a constant, so a dimension question can be answered by a controlled run
instead of a rebuild. **The default is and stays `fig10a`.**

| `--panel` | Dimensions | What it is | Legal mesh | Drift denominator |
|---|---|---|---|---|
| **`fig10a`** (default) | 3000 × 2250 | the **specimen**, from a dimensioned drawing that is to scale and reconstructs from its own pixels to 0.17% | **50**, 25, 10, 5 | 2250 |
| `table2` | 3000 × 2680 | **his analysed grid**, from inverting Table 2's node and element counts | 10, 5 | 2680 |
| `table2t` | 2680 × 3000 | the **transposition** of the same inversion, which reproduces his 20,385 nodes and 80,684 struts just as exactly (D82) | 10, 5 | 3000 |

**Fig. 4(f) is not a panel option, and never was.** Its `1500` is not a model height: the figure
shows the **bottom 1500 mm** of the model, which is why that number appears at all (author,
2026-09-05). A cropped view carries no panel dimension, so it cannot be a candidate geometry.
This also explains the anomaly D72 recorded — the drawn block's aspect matching none of its own
labels — without needing the figure to be badly drafted.

Thickness is a separate number, `--tw`, defaulting to **210** (D73). Pass 120 for the one-panel
A/B, or 240 and 350 to price the other readings of the plan stack.

Three things follow, and each belongs in the report of any run that leaves the default:

* **Panel and mesh are coupled, and it is worse than it looks, so the non-default panels are for
  ELASTIC comparison only.** The mesh has to land on the 100 mm bar spacing and the 50 mm cover as
  well as on the panel, and `gcd(3000, 2680, 100, 50) = 10`. So `fig10a` runs at mesh 50 while both
  Table 2 panels force **mesh 10**. Elements go as `1/m²` and steps as `1/m`, so total cost goes as
  `1/m³`: mesh 25 is 8× mesh 50 and **mesh 10 is 125×**. A nonlinear run on his panel is therefore
  not affordable here, and it does not need to be — `replica/` already runs it at mesh 20, because
  its covers of 40 and 80 admit that grid. Use these panels for the elastic and geometry questions
  and leave the nonlinear work on his panel to the replica.
* **The drift denominator moves with the height**, so the test's 20 mm is 0.89% on `fig10a` and
  0.75% on `table2`. Every run prints its denominator, and cross-panel comparisons are made in
  millimetres, never in drift.
* **`table2` here is not a replacement for `replica/`.** The replica reproduces his model — his
  mesh, his calibration field, his ν, his material laws — and remains the instrument for comparing
  against his published curves. Running `--panel table2` in this study varies geometry alone while
  everything else stays ours, which is a different and narrower experiment.

**References every run is compared against:** Table 4's two numbers, and the digitized Fig. 10(b)
(`data/fig10b.npz`, D77) which supplies the measured cloud, its outer envelope, and both of Aydin's
own curves. Two limits of that figure carry into every comparison: the frame **clips the data at
+16 mm** while the test reached about 20 mm (582 dots pile up at the edge), and the legend occludes
the negative quadrant at positive displacement. The envelope is an **upper bound** on the backbone,
not the backbone.

---

## 3. The two concrete models

Named for their compression branch, as the author asked.

### `lincomp` — his published law

Tension linear to `f_t` at `ε_cr = f_t/E = 7.44e-5`, then trilinear through `(a1·ε_cr, b1·f_t)` and
`(a2·ε_cr, b2·f_t)` to zero at `a3·ε_cr`. Compression **linear at E forever**, no cap, no crushing;
compressive failure is meant to emerge as indirect tensile splitting followed by lattice
instability. Shape constants `a1 = 1.5` (tied to the horizon by Table 1's footnote), `b1 = 0.6`,
`b2 = 0.2`, all three **published** in the Fig. 2 flowchart box.

Implemented today by `materials.concrete_lattice_aydin` (`ElasticMultiLinear`). Monotonic only:
the law is path-independent, so a cracked strut recovers full stiffness on reload.

### `eppcomp` — compression-capable, new

Tension identical. Compression linear at E to `f_c,strut = f_c × fc_scale`, then **perfectly
plastic**. Consequences worth stating before anything is run:

* there is no compression softening, so `G_fc` and the crushing strain cease to be parameters and
  **regularization applies to tension only**;
* the damage figure's "crushed" class becomes "yielded in compression", which is a different
  statement and must be relabelled;
* nothing in the model can now lose compressive load-carrying capacity, so a compression-driven
  collapse cannot occur by material failure. That is a deliberate property of an EPP branch, not an
  oversight.

**No existing material does this.** It is the one new constitutive object the study needs.

### Cyclic carriers

Both laws are path-independent as written, so a cyclic run needs a hysteretic carrier with the same
envelope. `HystereticSM` is the available vehicle and it has a measured defect (D79): a cycled strut
holds **4.38 MPa at zero strain**, so the wall does not pinch, and the pinching parameters are inert
against it (0.3/0.1 → 0.10/0.03 moved it to 4.43). Aydin's own unloading is **origin-oriented**
(Fig. 1(c)), which is what `HystereticSM`'s `beta` is for, and both September runs used `beta = 0`.

**Prerequisite task, one second of solver time:** cycle a single strut under each candidate carrier
and check two things — the stress at zero strain after cracking, and whether the unloading path
points at the origin. No cyclic wall run is worth launching until that returns a carrier that
pinches. `examples/compression_cube` has the pattern for a single-material probe.

**Amendment — that check has already been run, and `HystereticSM` failed it (D79).** This section
was written as though the carrier were an open question; it is not. The single-strut probe was run
and then confirmed structurally on a 22 h replica cyclic run (t = 210, mesh 20, horizon 1.5, bars to
the top, explicit, three loops to ±4.02 mm):

| carrier | strut stress at zero strain | wall force retained at zero displacement | \|pull\|/push |
|---|---|---|---|
| `HystereticSM` (his trilinear tail) | **4.38 MPa** | 0.56 / 0.36 / 0.72 | 1.103 |
| `Concrete02` (regularized) | **0.29 MPa** | 0.23 / 0.05 / 0.40 | 0.915 |

Three things follow, and the third is the one that moves the plan:

* the pinching parameters are inert against this (0.3/0.1 → 0.10/0.03 moved the residual to 4.43):
  they shape the APPROACH to zero, not the value AT zero, and the value at zero is what decides
  whether a wall pinches;
* the push/pull asymmetry reverses with the carrier — 1.103 under `HystereticSM`, 0.915 under
  `Concrete02`, the ordinary direction for a wall already damaged by the push — so it was the
  material, not the wall;
* **`Concrete02` is the only carrier in the repo demonstrated to pinch**, which promotes `c02` from
  an extension slot on the compression axis (§7) to the default carrier of Stages 5–6. Its cost is
  that it is not Aydin's envelope, so a cyclic cell must report which envelope it actually ran.

**Both laws need a carrier, not just `lincomp`.** `eppcomp` is `ElasticMultiLinear` too, so it is
path-independent by construction and inherits this whole problem unchanged.

**What is still genuinely open** is the one item D79 left as a one-strut check: Aydin's unloading is
**origin-oriented**, `HystereticSM`'s `beta` is what that parameter is for, and both September runs
used `beta = 0`. Probe `beta = 1` (the secant through the origin) against `Concrete02` on a single
strut before Stage 5. If `beta = 1` does not bring the residual down, the Aydin envelope cannot be
cycled at all and Stages 5–6 run on `Concrete02` alone.

---

## 4. Bond

**The law, per the author (2026-09-05): a bond link is a concrete strut that drops to ~60% of its
capacity and stays there.** Elastic to `F_cr`, brittle fall to `0.6·F_cr`, flat plateau, symmetric.
That is exactly `materials.bond_elastic_brittle(grade, tag, residual=0.6)` at its **default full
concrete-strut area** — verified: elastic to 1.85 MPa at ε = 7.44e-5, then a flat 1.11 MPa. Nothing
new is needed.

This supersedes the force-slip formulation of D75, which existed only because the link area had
been scaled down 100× to fix stiffness, and which is a departure from the author's description
rather than a fix.

**Topology.** Each bar gets its own steel nodes, duplicated at the concrete node positions on its
path; bond links join each steel node to the concrete nodes within `horizon × mesh`, excluding the
coincident one — a ring of 8 at horizon 1.5.

**The artefact, and how to report it.** A full-area ring is a parallel load path beside an intact
concrete lattice, so a bonded model comes out *stiffer* than perfect bond: measured **2.19×** the
transformed section at mesh 50. That is not an implementation error to calibrate away; it is what
this topology does at this bar density, and the density is the variable:

| | bar lines | grid | nodes carrying a steel node |
|---|---|---|---|
| this study, mesh 50 | 30 vertical + 22 horizontal | 61 × 46 | **73.5%** |
| Aydin's own model, mesh 20 | 30 + 26 | 151 × 135 | **35.3%** |

His model should therefore show a much milder artefact than ours, consistent with his K_sim being
0.91× the measured stiffness rather than twice it. **Every bonded run reports `K/K_perfect`** from
its own elastic twin, so the artefact is measured per configuration rather than assumed.

The properly-conditioned alternative — a coincident-node `zeroLength` interface spring, which *can*
reduce to perfect bond in the limit — remains unbuilt and out of scope here.

---

## 5. Calibration, in the order the author specified

1. **Elastic** — the energy balance, already validated on this wall: `K_lattice / K_continuum` =
   **0.9985** at mesh 50 and t = 210 against a plane-stress continuum on the same grid, same rebar,
   same boundary conditions. The remaining gap to a cantilever formula (0.914) is beam theory's own
   error at aspect ratio 0.75, where shear carries 57% of the flexibility. Not a matrix axis.

2. **Tension softening** — two modes, both offered as *different tensile laws* per the author:

   | mode | a2, a3 at mesh 50 | energy per strut |
   |---|---|---|
   | `solved` | 13.6 / 70 orthogonal, 9.7 / 50 diagonal | G_f by construction, mesh-objective |
   | `paper` | 70 / 360 as printed | **5.2 × G_f** |

   The published pair is fitted at his 20 mm grid, where it dissipates 2.08 × G_f by the same
   identity; the tail's energy scales with strut length, so transplanting it to a 50 mm grid is not
   neutral. **Every run prints the implied G_f multiple** so this stays visible rather than buried.

   *Open question for the author:* does the 15% energy check in his Fig. 2 compare the coupon's
   total dissipation to G_f per cracking strut, or to G_f over a gauge spanning several? That would
   explain the factor and settle which mode is faithful.

3. **Compression** (`eppcomp` only) — a scalar `fc_scale` on the strut strength, exposed as a plain
   argument alongside `fc`, `ft` and `gf` so it can be swept rather than switched.

   **It is direction-dependent, and that is a real limitation of a single scalar.** Across a
   horizontal cut the vertical struts present 0.612 of the gross face; an inclined family along 45°
   presents 0.852. A squat wall carries a diagonal compression field, so:

   | target | fc_scale |
   |---|---|
   | material strength, uncorrected | 1.000 |
   | assembly reaches f_c under **axial** compression | **1.633** |
   | assembly reaches f_c along a **45° strut** | **1.174** |

   Whichever is chosen, the other direction is off by 1.39×. Recorded before any run rather than
   discovered after one.

---

## 6. The cyclic protocol is invented, and must be labelled as such

**The 2019 paper never prints the loading history** — not the number of cycles, the amplitudes, or
the order — and the 2017 test paper is unavailable. Unlike WSH3, where the model was driven through
the specimen's own published protocol, any cyclic run here drives a history we made up.

Fig. 10(b) does not give it back, but it constrains it:

| evidence from the digitized cloud | reading |
|---|---|
| 582 dots piled at +16 mm, 32 at the left edge | the frame **clips** the record; the test reached ~20 mm |
| dot density shows no clean tip spikes | amplitudes are not recoverable from actuator dwell |
| 2–4 distinct dot bands at sections from 6 to 14 mm | on the order of **4–6 full loops**, not 16 |

So **one cycle per level** is the defensible shape; an eight-level ladder at two cycles each would
put sixteen loops on a figure that does not have them. This is inference from band counts on an
unordered cloud and belongs in the report as an assumption, not a recovery.

Cost is set by the total drive **path**, which the largest amplitude dominates: about
**0.032 h/mm without bond and 0.110 h/mm with**, from measured step costs.

**Amendment — the bonded figure is ~0.178 h/mm, not 0.110 (2026-09-05, D87).** Bond costs **5.7x**
nobond, not 3.4x, and the correction is structural rather than a mis-measurement: bond links sit on
light steel nodes, so `dt_crit` falls from 10.0 to 4.5 us and the same drive path takes **248,088
steps instead of 111,017** (2.23x), on **2.57x** the elements. The step count is exact; the per-step
factor is the element ratio until a bonded run completes. **Every bonded row in the table below is
therefore ~1.6x its printed cost** — the eight-level ladder is ~49 h with bond, not 30.

| candidate protocol (1 cycle/level) | path | no bond | with bond |
|---|---|---|---|
| 8 levels to 1.0% drift | 274 mm | 8.8 h | 30 h |
| same at 2 cycles/level | 549 mm | 18 h | 60 h |
| 5 levels to 0.5% | 99 mm | 3.2 h | 11 h |
| 5 levels to 0.3% | 72 mm | 2.3 h | 8 h |

**Deferred by agreement.** The protocol is a parameter with named presets, and the values are chosen
after the pushovers (Stage 4) show where the model actually turns over.

---

## 7. The parameter matrix

Eight model configurations: compression × tail × bond. Materials, mesh, horizon and protocol ride on
top as arguments with defaults.

| axis | code | values now | extension slots |
|---|---|---|---|
| compression law | `comp` | `lincomp`, `eppcomp` | `c02` (the repo's own control law) |
| tension tail | `tail` | `solved`, `paper` | a real Fig. 2 coupon loop |
| bond | `bond` | `nobond`, `bond-a06` | other residuals, `zeroLength` interface |
| analysis | — | `elastic`, `static`, `pushover`, `cyclic` | |
| f_c, f_t, G_f, fc_scale | `fc` `ft` `gf` `fcx` | 28, 1.85, 0.075, 1.0 | any |
| panel, thickness | `panel` `tw` | `fig10a`, 210 | `table2`, `table2t`, any custom L × H; 120, 240, 350 |
| mesh, horizon | `m` `h` | 50, 1.5 | 25, 10; 3.01 |
| protocol | `proto` | named presets | any level list |

Elastic runs depend only on bond, so there are two of those and not eight.

**Amendment — the `comp` axis is not free in Stages 5–6.** `Concrete02` is not only a carrier; it
brings its own capped, softening compression branch. Adopting it per the §3 amendment therefore
partly COLLAPSES this axis for the cyclic cells: a cyclic cell labelled `lincomp` would not have
linear compression, so the label would lie. Two honest readings, to be settled when Stage 4 names
the winning cells — run the cyclic stage as one carrier and drop `comp` from the cyclic matrix, or
keep `comp` and accept that its cyclic cells answer a different question from its pushover cells.
The pushover stages are unaffected: `lincomp` and `eppcomp` are monotonic and mean exactly what
they say there.

---

## 8. Staging, with gates

Each stage answers something and gates the next. Nothing past Stage 0 is launched until Stage 0 is
merged, because three of its items silently corrupt bonded and cyclic results.

**Stage 0 — prerequisites. Code only, no solver.** **COMPLETE (2026-09-05, D83/D84/D85/D86).**
Marked per item below.

1. **Bond base anchorage.** Supports resolve on the concrete coordinate array, so no steel node is
   ever fixed: on this wall **30 steel nodes sit at y = 0 unanchored**, held only by their own ring,
   which fails at 0.0037 mm of slip. Meanwhile `select_nodes` queries the built model and *does*
   pick up coincident steel nodes, so with bars to the top the drive would be imposed on steel at
   the top while the base stays free. Every bonded result to date carries this. (D80 item 4)
   **The fix has to cover BOTH ends.** With bars to the top (item 4) the same defect becomes a
   DRIVE-side defect as well as a support-side one, which is exactly the combination Stage 6 runs.
   Verify it with an assertion in the build — every prescribed and every fixed DOF accounted for —
   not by inspecting a figure.
   **DONE (D85).** Supports propagate to coincident steel nodes (base steel unanchored 30 → 0), and
   `select_nodes(kind=)` splits the two sets in opposite directions: the drive is concrete-only, the
   base reaction set is everything. Verified with bond + bars-to-top; no change without bond.
2. **`strut_groups` classifies bond links as concrete** by orientation, so the bonded run's
   "rebar 0.7% of the overturning moment" is a probe artefact. (D80 item 5)
   **DONE (D85):** bond links get their own group.
3. **`run_cyclic_dynamic` reports `("HHT", 0.7)` unconditionally** even when marched explicitly, and
   keeps the implicit sub-step rescue under an explicit integrator, against D74's rule. (D80 item 6)
   **Worse than recorded, and the fix already exists next door.** The rescue ladder ends with an
   unconditional `ops.algorithm("Newton")`, so after ONE failed explicit step the run continues as
   Newton on a `Diagonal` system for every remaining step and nothing in `data.json` says so — a
   million-step cyclic run can be half-explicit and read as neither. `run_pushover_dynamic` already
   gates the ladder on `explicit` and reports `tuple(integrator)`; the cyclic runner is the stale
   copy. While in there: its `capture` block tests the PREVIOUS step's shear against the CURRENT
   displacement field, so `disps_peak` is one step stale — the same field D78 read to find the
   tear.
   **DONE (D83):** all three fixed, transcribed from `run_pushover_dynamic`, which was already
   correct.
4. Parent specimen gains a bars-to-top option (the replica already has one).
   **DONE (D84):** `rebars(full_height=)`, `wall_lattice(full_height_rebar=)`,
   `pushover.py --rebar-to-top`, tagged `_rebartop`. Costs +60 elements (two cells, not one), and
   with bond +60 steel nodes / +450 elements. `elastic.py` and the continuum twin deliberately
   still call `rebars` without it.
5. The new `eppcomp` material — a two-point edit to the `neg_eps`/`neg_sig` pair of
   `materials.concrete_lattice_aydin`, not a new object — and the carrier, which per the §3
   amendment is `Concrete02` unless the `beta = 1` probe rescues the Aydin envelope.
   **DONE (D86).** `eppcomp` is `concrete_lattice_aydin(fc_cap=)`, exposed as `--comp eppcomp`
   with `--fcx`. The probe (`carrier_probe.py`) was run: `beta = 1` moves the stress at zero strain
   only 7.70 -> 6.59 MPa, so it does NOT rescue his envelope, and **Stages 5-6 run on
   `Concrete02`** (0.000 MPa at zero, and it unloads at 1.000 x the origin secant — his own
   Fig. 1(c) rule, better satisfied than by `HystereticSM`'s beta).

**Stage 1 — confirm the D78 diagnosis on this wall.** One explicit pushover of the existing model
with bars to the top. **11 minutes.** *Gate:* no tear at the driven row, and drift capacity beyond
the 0.237% the bonded run reached. If the tear persists the diagnosis is wrong and the plan pauses.
**Set the target well past the gate** — D78's replica runs used 0.30%. A run launched at
`--drift 0.0025` ends exactly at the number it is being asked to beat, which is D78's own PROCESS
trap (its 2026-08-30 control ended "still ascending, converged" ~3,000 steps from the cliff), and it
would license nothing.

**Stage 2 — elastic gates.** Two runs, seconds. *Gate:* `K/K_continuum` ≈ 1.00 without bond;
`K/K_perfect` reported and understood with bond.

**Stage 3 — static pushovers.** Eight runs, minutes. Expected to stall near 0.03% drift on this
early-cracking wall; kept as a diagnostic of where and how, not as a capacity.
**DONE for the nobond half (D87):** all four stall while still ascending at **0.0098-0.0120%**,
three times earlier than predicted here. `lincomp` and `eppcomp` stall at the identical point, since
nothing has reached compressive yield that early.

**Stage 4 — quasi-static pushovers. The decisive stage.** Eight runs, **3.2 h total**. Compared on
peak against 963.6 kN, drift at peak against ~1%, load-path split, damage location, and each run's
own residual. *Gate:* this is what narrows tail and compression before any cyclic time is spent.
**DONE for the nobond half (D87), and it narrows the matrix.** All four converged to the 0.30%
target with ascending-branch residuals of 0.6-1.3%:

| comp | tail | peak (kN) | /test | at drift | Gf multiple |
|---|---|---|---|---|---|
| `c02` (control) | solved | 943.8 | 0.979 | 0.266% | 1.00 |
| `lincomp` | solved | 944.2 | 0.980 | 0.292% | 1.05 |
| `eppcomp` | solved | 968.3 | **1.005** | 0.2997% — see caution | 1.05 |
| `lincomp` | paper | 1,158.0 | 1.202 | 0.284% | **5.26** |
| `eppcomp` | paper | 1,158.0 | 1.202 | 0.221% | **5.26** |

* **the compression law is nearly irrelevant to peak strength** — three different branches span
  2.6% at the solved tail, and 0.003% at the paper tail — because only **2 to 19 of 10,905**
  concrete struts ever pass compressive yield. What "linear forever" costs is visible in the strain
  field instead: `lincomp` runs a strut to 7.69x the cap (~215 MPa), which `eppcomp` holds at 28;
* **the tail dominates and is a MESH ARTEFACT**: +22.6% for both compression laws, bought with a
  measured 5.26 x Gf. So `paper` is not a physical alternative at mesh 50, and bonded paper-tail
  cells are not worth their ~2.4 h;
* **CAUTION on `eppcomp`/solved**: it peaked at 0.2997% of a 0.3000% target, i.e. still ascending,
  so 1.005 is a LOWER BOUND and not a capacity. Being re-run further before it is quoted.

**The gate is peak force, load path and damage location — NOT drift at peak.** §1 lists drift
capacity as not a fair claim from this study, and the ~1% it would be gated against comes from a
record the frame clips at +16 mm (§2). Report drift at peak in every cell; do not let it decide a
cell.

**Stage 5 — cyclic without bond.** Four cells at a protocol chosen from Stage 4.

**Stage 6 — cyclic with bond.** Two cells, the best compression × tail combination.

Running all eight cyclic cells at the eight-level protocol is about seven days of solver time.
Staging as above, with a 0.5% protocol, is under a day and a half.

---

## 9. Layout

Scripts in a new subdirectory of the Aldemir package; the existing package and `replica/` are
untouched.

```
examples/aydin_aldemir_wall/study/
    PLAN.md          this document
    params.py        the parameter registry — one record per parameter
    models.py        resolved parameters -> Problem, materials, model
    protocols.py     named cyclic protocols
    references.py    Table 4 + the digitized Fig. 10(b), and matched-displacement comparison
    run.py           single entry point for every analysis
    report.py        per-run report
    master.py        matrix-aware master report

examples/output/aydin_aldemir_wall/study/
    <stamp>_<analysis>_<comp>-<tail>_<bond>[_extras]/
        command.txt      how it was invoked
        params.json      EVERY parameter, including defaults, plus a schema version
        console.log      teed, so OpenSees warnings land in it
        data.json        the response, groups, fields, residual series
        figures/         backbone, damage, load path, comparison
        report.md        generated, self-contained
    master_report.md     the matrix, every cell, cross-run figures
```

Example directory names, extras appended only when non-default:

```
2026-09-06_1030_pushover_lincomp-solved_nobond
2026-09-06_1412_cyclic_eppcomp-paper_bond-a06
2026-09-07_0900_pushover_eppcomp-solved_nobond_fcx1.63_h3.01
```

**Names are for humans; `params.json` is for the machine.** The master report builds the matrix by
reading parameter files, never by parsing directory names, which is what lets a new parameter or a
new value extend the matrix without touching the reporting.

### Designed for the parameters that do not exist yet

The author's note that parameters and their values will grow is a design constraint, not a caveat:

* **`params.py` is a registry**, one record per parameter carrying its name, short code, default,
  allowed values or type, help text, and what it *affects* (model, analysis, or reporting only).
  Adding a parameter is one record; the CLI, the directory name, `params.json` and the master report
  all derive from it.
* The `affects` field lets the master report know which runs are comparable, so a new
  reporting-only parameter never fragments the matrix.
* **Named presets** (protocols, and later material sets) live in their own module and grow without
  touching `run.py`.
* **Report sections are a list of builders**, so a new figure is one entry rather than an edit to a
  monolithic function.
* `params.json` carries a **schema version**, so runs made before a parameter existed stay readable
  and the master report can fill defaults for them.

---

## 10. Reporting

**Per-run report**, generated into the run directory, self-contained and readable alone: the full
resolved parameter set with each value's source (measured, code convention, inferred, assumed); the
model as built, with element counts by kind; calibration output including A_t, EA, the implied G_f
multiple and, for bonded runs, `K/K_perfect`; the response against Table 4 and the digitized
references at *matched displacement*; the measured residual, separated into start-up and steady;
damage and load-path figures; and an explicit statement of what the run does and does not license.

**Master report**, reading every `params.json` and `data.json`: the matrix as a table with one row
per cell, cross-run figures grouped by the axis being varied, and the cells not yet run shown as
gaps rather than omitted, so the matrix reads as a plan and a result at once.

---

## 11. Decisions still open

| | |
|---|---|
| cyclic protocol, levels and cycles | deferred to after Stage 4, by agreement |
| `fc_scale` target — 1.0, 1.633 axial, or 1.174 diagonal | a swept argument; may not need settling |
| does his 15% energy check use G_f per strut or over a gauge | a question for the author; decides which tail is faithful |
| mesh 25 as a convergence check on the winning cell | optional, monotonic only |

---

## 12. What would change the plan

The 2017 test paper, if it can be obtained, would supply the protocol, the section, the bar type and
possibly the ordered hysteresis — removing four assumptions at once and turning drift capacity and
cyclic energy into fair comparisons. It is the single highest-value input to this study and the
author's co-author is on both papers.
