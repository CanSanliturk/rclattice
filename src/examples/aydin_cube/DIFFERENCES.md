# `aydin_cube` vs `compression_cube` — every major difference

Two cubes, both plain concrete, both 2D horizon lattices, both squashed by dynamic relaxation,
both calibrated by the same elastic energy balance. They disagree about almost everything else,
and the disagreements are not incidental — each one follows from a single root choice.

**The root choice: is `fc` an input or an output?**

* `compression_cube` puts a compressive strength on every strut (Concrete02) and measures what
  fraction of it the assembly delivers. Answer: **0.634 fc**, every time, for geometric reasons.
* `aydin_cube` puts **no compressive strength anywhere**. Compression is linear elastic forever.
  The cube's strength emerges as the stability loss of a load path that transverse tension has
  cracked free, and it is calibrated to `fc` through a *geometric* parameter, `Rmax/d`.

Everything below is downstream of that.

**Since D71 both examples have a strength-calibration knob, and comparing the two knobs is the
sharpest statement of the difference.** `compression_cube --peak-correction strength` scales the
strut grade's `fc` by a factor derived from the section (1/0.6341 = 1.577 at horizon 1.5), so its
strut strength becomes a *calibrated* input rather than the material value — the same status
`Rmax/d` has in `aydin_cube`. The knobs are not equivalent, and how they differ is the point:

| | `compression_cube` (D71) | `aydin_cube` |
|---|---|---|
| what is calibrated | the struts' **constitutive** strength `fc` | the grid's **geometry**, `Rmax/d` |
| how the value is obtained | **derived** in closed form from the strut list, before any run | **searched** — repeated FE runs until the mean lands within 10% of `fc` |
| what it fixes | the peak **stress** only | the peak, the softening branch, and the failure mode together |
| what it costs | peak **strain** moves by the same factor; splitting still starts at the same absolute strain, so relatively *earlier* | run-to-run scatter (CoV 7-11%), 5 realizations per answer |
| the mechanism afterwards | unchanged — still strut crushing after splitting sheds the diagonals | is itself the answer — stability loss of split columns |

One corrects an *outcome*; the other changes what fails and why.

> **A third configuration now exists.** `examples/compression_cube/aydin_approach.py` runs the
> *right-hand column of every table below* on the *left-hand column's specimen* — our 200 mm cube,
> our fc = 15.1 MPa grade, our smooth platens, unchanged — so that the approach is the only variable.
> `specimen.py`, `build.py` and `pushover.py` are untouched by it. See §10.

---

## 1. Source

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| paper | Aydin (2017), METU MSc thesis — *Overlapping lattice modeling…* | Aydin, Binici & Tuncay (2021), *Mag. Concr. Res.* 73(8):394–409 |
| what it takes from it | the elastic energy balance (Eqs 2.1–2.3) only | the energy balance **and** the entire compression method |
| the thesis on compression | explicitly out of scope (Sec. 2.2, p. 30) | this paper *is* that missing study |

`compression_cube` exists to price the thesis's exclusion. `aydin_cube` exists to reproduce the
paper that lifted it. They are not competing models of the same thing.

## 2. Material model — the difference that generates all the others

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| OpenSees material | `Concrete02` | `ElasticMultiLinear` |
| compression | **capped at `fc`**, parabola to `epsc0`, softening to a `0.2 fc` residual plateau | **no cap, ever.** Linear `E`; with RSM, `0.4 E` beyond `epsc0/3` |
| tension | bilinear: linear to `ft`, one softening slope `Ets` to zero | trilinear: `ft → 0.6 ft → 0.2 ft → 0` at `a1, a2, a3 × eps_cr` |
| tension regularization | `Ets = ft²L/(2Gf)` | `a2, a3` solved per strut so `L·∫σ dε = Gf` |
| compression regularization | `epsU = epsc0 + 2Gfc/((fc+fcu)L)`, `Gfc = 250 Gf` | **none — there is nothing to regularize** |
| unload path | plastic, with `lambda = 0.1` unloading ratio | path-independent (nonlinear elastic) |
| free parameters beyond the grade | `Gfc/Gf`, `residual_ratio` | `Rmax/d`, and the fixed shape `a1, b1, b2` |

Consequences:

* `aydin_cube` can never crush a strut. The damage figure legends say so, and the run prints
  "NO strut ever crushed" as a fact about the model, not an observation about the result.
* `compression_cube` cannot produce a strength above `0.634 ×` **the struts' own `fc`**, because
  that is the vertical struts' projected share of the face and they cap at their `fc`. **Grid
  perturbation cannot fix this** — perturbation only ever *lowers* strength, so there is no
  `Rmax/d` that reaches `fc`. The two ingredients are not separable: Aydin's calibration
  presupposes his constitutive law.
  Since D71 the struts' `fc` is itself a knob (`--peak-correction strength`), so the assembly *can*
  be made to peak at the material `fc` — by scaling the strut strength up 1.577×, not by touching
  the grid. That is a correction to the outcome, not a repair of the mechanism: read the row for it
  in §8 and §10 with the caveat attached.
* The path-independence in `aydin_cube` is a deliberate trade, documented in
  `materials.concrete_lattice_aydin`. `HystereticSM` reproduces the same envelope but unloads on
  the initial stiffness, putting a fully cracked strut at **−10 MPa (RSM) / −44 MPa (no RSM)** when
  its strain returns to zero — it would shove the split columns apart and destroy the mechanism.
  The price paid instead is no damage memory: a softened strut recovers on reload, where Aydin's
  sequentially linear analysis holds the reduced secant. Fine under monotonic compression,
  **unusable for cyclic work**.

## 3. Grid

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| node layout | structured, uniform | perturbed: `R ∈ [0, Rmax]` at `θ ∈ [0, 2π)` per node |
| `Rmax/d` | — | 0.075 (LR), 0.060 (NR) — the paper's fitted values for this specimen |
| strut topology | horizon on the grid | horizon on the **unperturbed** grid, then nodes moved |
| determinism | deterministic, 1 run is the answer | random; 5 realizations, report the mean |

Why topology is frozen: re-running `connect_horizon` on moved nodes drops the struts that happen
to cross `horizon·mesh_size` — 5 of 420 at `Rmax/d = 0.08` — which would make connectivity a second
random variable on top of geometry. Aydin perturbs a grid and keeps its cells (his Fig. 2(b) still
shows both diagonals everywhere).

Boundary handling is a **documented deviation**: the paper does not say what it does with boundary
nodes. Here edge nodes slide only *along* their own face and corners are pinned, so the faces stay
flat, the specimen's dimensions stay deterministic, and box selection of supports and driven nodes
stays exact.

## 4. Element formulation

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| strut element | `Truss` (small displacement) | `corotTruss` (corotational) |

Not a stylistic default. Aydin's failure mode is a *geometric* event — the loss of stability of a
zigzag column. D53 already measured what happens without it: unbounded compression on straight
`Truss` struts climbs monotonically to **6·fc with no failure at all**, while the same run on
corotational struts peaks and collapses. Aydin's explicit scheme integrates the nodes' actual
positions, so he gets this for free and never mentions it.

## 5. Calibration

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| elastic | Aydin energy balance → one uniform `EA` | same, but **re-run per realization** (perturbation changes every strut length) |
| strength | none — `fc` is an input | `Rmax/d` searched until mean strength is within 10% of `fc` (his Fig. 4 loop, `--calibrate`) |
| verification target | closed-form: measured tangent vs `E·(1−ν_lattice²)/(1−ν²)` | the same balance applies; the reported check is the strength itself |

`compression_cube`'s finding that the balance pins the **confined** modulus `E/(1−ν²)` — so an
unconfined cube reads `0.88 E` — applies here unchanged. It is a property of Eq. 2.1, not of either
specimen.

## 6. Specimen and boundary conditions

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| size | 200 × 200 × 200 mm | 100 × 100 × 100 mm (the paper's) |
| grid spacing | 20 mm | 10 mm (the paper's) |
| `fc` | 15.1 MPa (SW-NC-FF test unit, shared with the wall) | 20 MPa (the paper's) |
| `ft`, `E` | measured: 1.5 MPa, 16 100 MPa | derived: `0.35√fc`, `4700√fc` |
| `Gf` | 0.053 N/mm (CEB-FIP MC90) | 0.050 N/mm (chosen, as in the paper) |
| end conditions | one: rollers + a single pin (smooth platens ≈ NR) | two: **LR** (pinned base + laterally held top) and **NR**, `--boundary` |
| platen study | `--rigid-platen`, found inert under a prescribed ramp | LR vs NR, which the paper shows changes both strength and crack pattern |

`compression_cube` argues for smooth platens as the low-friction bound and never models the other
one. `aydin_cube` models both, because the paper's whole Fig. 3 is a two-column comparison and
because LR is the case that develops the diagonal cracking.

## 7. Solution

| | `compression_cube` | `aydin_cube` |
|---|---|---|
| solver | dynamic relaxation (`--solver static` also available, stalls pre-peak) | dynamic relaxation only |
| drive rate | 0.5 mm/s on a 200 mm cube | 1.0 mm/s on a 100 mm cube |
| runtime | ~52 s | ~25 s per realization, ~2 min for the set |

Both stand in for Aydin's "dynamic explicit integration at a sufficiently slow speed". A static
path-follower is not offered in `aydin_cube` on purpose: the response it must trace is a stability
loss, which is exactly the mechanism D38/D44/D46/D57 established no static solver crosses.

## 8. What each one actually predicts

| | `compression_cube` | `aydin_cube` (5 realizations) |
|---|---|---|
| peak | **0.634 fc**, at `0.998 epsc0` | **LR 0.924 fc** at `1.00 epsc0`; **NR 1.031 fc** at `1.53 epsc0` |
| peak, strength-corrected | **1.00 fc** by construction, at `1.577 epsc0` — D71, *predicted, not yet run* | — (its calibration targets `fc` by search instead) |
| vs the paper's criterion | — | −7.6% and +3.1%, both inside his 10% (Fig. 4), **at his own fitted `Rmax/d`, no refitting** |
| after peak | a residual plateau — the material's, not the structure's | genuine softening to ~0.32 of peak by 0.4% strain |
| uniform-grid control | *is* the uniform-grid case | **1.506 fc** at `--rmax 0` — the locking he reports as ~2 fc |
| failure mode | strut crushing, once splitting has shed the diagonals | stability loss of split columns; **no strut ever fails in compression** |
| scatter | none | CoV 6.7% (LR) / 10.7% (NR), as in his Table 2 |

The two land on opposite sides of `fc` from the *same* lattice topology, and the reason is only the
constitutive law:

* uncapped compression on a uniform grid → **locks**, over-strength (1.51 fc measured here;
  6 fc without corotational struts);
* capped compression on a uniform grid → **0.634 fc**, because the verticals cap at `fc` and
  occupy 0.634 of the face;
* uncapped compression on a *perturbed* grid → **0.92-1.03 fc**, the paper's result.

`compression_cube`'s conclusion that "the section is short of a load path, not of material" is
correct arithmetic for its own model. This example is the other half of that sentence: Aydin's
answer is that being short of a load path is the entire point, and the residual error was putting
a compressive strength on the struts at all.

## 9. Deviations from the paper in this example (things Aydin does that this does not)

1. **RSM is applied to every strut**, not only to laterally cracked ones. That switch is
   state-dependent and cannot be a static material assignment. Under uniaxial compression
   essentially every compression-carrying strut does crack laterally before peak, so the two
   coincide over the part of the response that matters.
2. **No damage memory** (§2) — `ElasticMultiLinear` is nonlinear elastic where his SLA holds a
   reduced secant.
3. **Boundary nodes slide along their face rather than freely** (§3).
4. **`a2`, `a3` are solved from `Gf`** rather than fitted to the Cornelissen et al. (1986) curve.
   The algebra reproduces his Table 1 to ~12% (Jansen & Shah: 52.9/264 against his 60/300) and is
   what keeps the model mesh-objective at grid spacings he never ran.
5. **One extra `ux` restraint under NR.** His explicit scheme tolerates a free-floating rigid-body
   mode; the implicit Newmark solve here would see a singular stiffness.
6. **`ν = 0.20` is assumed** for the energy balance; the paper does not state one.
7. The paper is **self-inconsistent on `α` and `β`** (notation list and Fig. 1 caption vs the
   p. 397 text). The self-consistent reading — `α = 1/3` strain, `β = 0.4` modulus — is used.

## 10. The controlled comparison: our cube, both approaches

`examples/compression_cube/aydin_approach.py` holds the specimen fixed and swaps the six choices of
§2-§4 together. Everything is our cube: 200 mm, `fc = 15.1` MPa with its measured `E = 16 100` and
`ft = 1.5`, 20 mm grid, smooth platens (rollers plus one pin), the same energy-balance strut area.

| on OUR 200 mm cube | Concrete02 (`pushover.py`) | Aydin approach (`aydin_approach.py`) |
|---|---|---|
| initial tangent / expected | 1.008 | **0.997** (and 0.879 of grade `E`) |
| uniform grid, `Rmax/d = 0` | 0.634 fc — a *capacity* | **≥1.86 fc and still rising** — it LOCKS |
| calibrated, `Rmax/d = 0.075` | not applicable | **0.943 fc** (5 realizations, CoV 11.3%) |
| calibrated, strut `fc ×1.577` | **1.00 fc** by construction (D71, *unrun*) | not applicable |
| strain at peak | 0.998 epsc0 (1.577 epsc0 corrected) | 1.36 epsc0 |
| pre-peak shape | parabolic (Concrete02's own curvature) | near-linear with the RSM knee |
| post-peak | gradual softening to a 0.2 fc residual plateau | **total collapse to ~0** |
| inclined path at peak | 28% → 0% | 0.6% (uniform) / 3% (perturbed) |

`Rmax/d` was calibrated for this cube rather than borrowed: the Fig. 4 loop ran
0.050 → 1.235 fc, 0.062 → 1.171 fc, 0.093 → 0.765 fc, **0.075 → 0.919 fc**, converged. That it lands
on Aydin's own LR figure is a coincidence of two specimens differing in size, grade, spacing and end
condition.

Two things this controlled comparison establishes that neither example could alone:

1. **The elastic calibration is genuinely shared.** Both approaches hit the same closed-form
   modulus target on the same cube (1.008 and 0.998), so every difference in strength below is a
   difference of *approach*, not of two calibrations talking past each other.
2. **The load-path collapse is topological, not constitutive.** D54/D55 attributed the inclined
   struts' 28% → 0% collapse to splitting on a uniform horizon lattice. Swapping the entire
   constitutive law reproduces it (0.6% at peak on the same grid), which confirms it is a property
   of the lattice geometry — exactly as D55 claimed, now tested rather than argued. The Aydin run
   also reproduces the D53 confined-modulus result independently, reading **0.879 E** unconfined.

**Where it gets harder on our cube than on Aydin's.** Our 20 mm grid at `Gf = 53` N/m gives each
strut a much shorter softening tail than his 10 mm grid does — `a3 = 120` against his `269`, so our
struts fail at 1.1% strain where his last to 2.0%. Brittler struts make each failure a more violent
event for the transient solver, and some realizations diverge outright: one recorded a **139 MPa**
spike on a lattice with no compressive strength at all, reported `converged = True`, and dragged a
13 MPa mean to 56 MPa. `first_divergence` now truncates any record at the first physically
impossible upward jump (a quasi-static run's largest genuine per-step rise is ~0.0008 MPa against a
0.755 MPa threshold) and discards realizations that never turned over before it. A diverged solve is
a failed measurement, not a sample.
