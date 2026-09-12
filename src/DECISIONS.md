
### D47 — 2026-08-15 — RC shear-wall example (`examples/wall/`) calibrated by Aydin's OLM energy balance; a second, homogenization-based calibration route alongside D16
- **Context (user):** build a flexure-dominated shear-wall example that meshes the structure as a
  lattice, and calibrate its ELASTIC behaviour with Beyazit Bestami Aydin's method (METU MSc, 2017 —
  "Overlapping lattice modeling for concrete fracture simulations using sequentially linear analysis",
  Sec. 2.2, Eqs 2.1–2.3). Decided in discussion: specimen **SW-NC-FF** from Sahinkaya, Binbir, Orakcal
  & Ilki (2025), *Buildings* 15:4501; **model the pedestal** (not a clamped wall base); follow Aydin's
  methodology **only** — no continuum/beam-column reference model, and no D16 cross-comparison.
- **The method, and why it is a different animal from D16.** D16 (`calibrate_lattice`) is a STRUCTURAL
  match: run static + modal analyses of the whole structure and fit orthogonal/diagonal areas with
  `least_squares` against a reference model. Aydin's is a HOMOGENIZATION: impose an affine strain field
  `u = F.x`, equate the continuum's stored energy (Eq 2.1) to the lattice's with `EA = 1` (Eq 2.2), and
  divide (Eq 2.3). Because the field is affine, every strut elongation is exact — `dL = (F.v).v / L` —
  so it is a **closed-form sum over the strut list**: no FE solve, no optimiser, no reference model.
  New backend: `calibration.energy_balance_area` / `energy_balance_rectangle` → `EnergyBalanceResult`.
  D16 is untouched and remains the right tool for matching a specific BVP; the two answer different
  questions and are kept side by side (supersede-not-rewrite).
- **Two properties that shape how it is used.** (1) `EA_t` scales with `E`, so the calibrated AREA
  `A_t = EA_t/E` is **E-independent** — ONE area serves all three concrete zones of this specimen
  (test unit / upper head / pedestal), each zone's modulus riding on its material. (2) `A_t/(t*d)` is a
  constant of the (horizon, nu) lattice, i.e. **mesh-objective**: 0.5985 / 0.6043 / 0.6072 at mesh
  100 / 50 / 25 mm, monotonically approaching the interior closed form
  `[1/(2(1-nu^2))] / (1/2 + sqrt(2)/4) = 0.6102` from below (the deficit is the free-edge boundary
  layer). Verified in `tests/test_energy_balance.py` (10 tests).
- **Correction to Aydin's stated nu — the substantive finding.** Aydin takes `nu = 1/3` in Eq 2.1
  (from Hrennikoff/Silling). That constant belongs to a lattice whose orthogonal and diagonal areas
  DIFFER. With a single `EA_t` on the horizon=1.5 (8-neighbour) grid, matching E under normal strain
  and matching G under pure shear agree only at **nu = 0.1796** — a closed form falls out of the two
  balances, `nu = 1 - 2*W_shear/W_normal`. Assuming 1/3 instead overshoots the shear stiffness by
  **18.7%**. Fortunately concrete's nu = 0.20 sits almost on the lattice's own value (2.5% shear-
  stiffness error), and this wall is only 15% shear-compliant, so the residual effect is negligible
  HERE. `EnergyBalanceResult` reports `nu_consistent` and `isotropy_error` so the assumption stays
  visible rather than buried. Horizon 3.01 (Aydin's RC recommendation, 28 neighbours) sits at
  nu = 0.36 — worse for elastic fidelity on this specimen; his reason for 3.01 was nonlinear BOND
  behaviour, not elastic accuracy, so **horizon 1.5 is the default here** and 3.01 is left for a
  later nonlinear pass.
- **Specimen modelling.** Wall 3000x1000x200 on a 1900x600 pedestal as `CompoundRectangles` with the
  y=0 node row merged; fixed at the pedestal SOFFIT so the base region deforms. Three concrete casts
  as zones (Ec = 16.1 / 26.4 / 23.4 GPa). The pedestal is 800 mm wide against the wall's 200 mm and a
  plane model carries one thickness, so its extra width is folded into an equivalent modulus
  `E * PED_W/TW` — reproducing the block's real rigidity, consistent with the paper measuring its
  sliding and rotation as negligible. Reinforcement (two curtains collapsed onto one in-plane line):
  4phi14 boundary groups lumped at each region's centroid (x = +-450, preserving the first moment),
  phi12@200 web verticals, phi10@200 horizontals. Axial 600 kN on top, lateral line load at y = 2200.
- **Result (mesh 50, horizon 1.5; 1767 nodes, 6736 concrete + 680 rebar struts).** `A_t = 6043 mm^2`
  = 0.6043 * (thickness*mesh); a physically-sized tributary strut would be **1.66x too stiff**.
  Pushed to 0.10% drift (the test first yielded at 0.30%, so this is safely pre-cracking) the lattice
  gives **73.05 kN/mm**, against a transformed-section cantilever hand check of 77.29 kN/mm →
  **ratio 0.945**. Clamping the wall base instead of standing it on the pedestal gives 78.12 kN/mm,
  **ratio 1.011** — so the Aydin calibration reproduces the hand calculation to ~1%, and the 5.5% gap
  in the full model is the pedestal's base flexibility (6.5%), exactly what modelling it buys.
- **Two traps found and fixed during the build, both worth recording.**
  1. *The hand check must be TRANSFORMED, not gross.* First comparison put the lattice 11% ABOVE a
     plain-concrete cantilever and the code asserted the opposite sign. With Ec = 16.1 GPa the modular
     ratio is n ~ 12.3, so the vertical bars add 21% to I (`I_tr/I_g = 1.206`). Comparing a reinforced
     lattice to a plain section understates the target by that much. `build.transformed_inertia()` now
     computes it from the actual `Rebar` list (vertical bars only — horizontals run parallel to the
     cut and contribute nothing to flexural I).
  2. *A silent rebar drop-out (backend bug).* At mesh 100 the boundary bars at x = +-450 land on no
     WALL node, but the pedestal's node row happens to supply one at y=0; `rebar_node_chain` returned
     a 1-node chain, `zip(chain, chain[1:])` was empty, and the model was built with **no boundary
     steel and no error**. `reinforcement.rebar_node_chain` now raises when a path matches fewer than
     2 nodes, and `specimen.check_mesh_alignment` (called from `build.wall_lattice`) rejects a mesh
     that cannot place a node on every bar line — here, divisors of 50 mm.
- **Scope flag (stated in the specimen docstring).** The test's dominant mechanism is ROCKING from
  progressive debonding of the PLAIN longitudinal bars — 74–76% of tip displacement at 2% drift. The
  lattice ties rebar to concrete at shared nodes, i.e. PERFECT bond (D5/D13), so it cannot reproduce
  that. Irrelevant to this elastic calibration (rocking engages after cracking) but it caps what a
  later nonlinear/cyclic study on this specimen may claim; Aydin's own bond treatment (elastic-
  perfectly-plastic interface struts at 70% of concrete ft, his Sec. 3.4/3.5) is the hook if pursued.
- **Status:** accepted. New backend fns in `calibration.py` + a guard in `reinforcement.py`; new
  `examples/wall/{specimen,build,elastic}.py`; new `tests/test_energy_balance.py`. Suite: 35 pass, 1
  pre-existing known failure (`test_rc.py::test_nonlinear_pushover_runs_and_yields`, D34). Output →
  `examples/output/wall/wall_elastic.png` (+ `wall_model.png` with `--draw`). Not exported from
  `rclattice/__init__.py` by user preference — examples import from submodules, as they already do.

### D53 — 2026-08-17 — 200 mm compression cube (`examples/compression_cube/`); Aydin's balance pins the CONFINED modulus, and `nu_consistent` is not the lattice's Poisson ratio
- **Context (user):** a new example — a cube modelled in 2D, strut area from Aydin's energy balance,
  compression applied to the upper face, traced by quasi-static dynamic relaxation as in the wall
  study, with a damage figure. Settled in discussion: **200 mm** edge (a standard compression-test
  specimen; "200 m" as first written would have made self-weight 4.7 MPa ~ 16-31% of fc and needed
  body-force support the library does not have), **smooth platens**, **no reference model**,
  **Concrete02** struts, **mesh 20 mm**, **horizon 1.5** by default with 3.01 via `--horizon`.
- **Why this specimen is worth having.** Aydin excludes compression explicitly (Sec. 2.2 p. 30:
  "Compression performance of OLM is not studied within the scope of this thesis"; future work p. 72
  asks for compressive fracture models). This prices that exclusion on the simplest possible
  specimen. `--compression crushing|elastic` selects the real Concrete02 backbone or his literal
  assumption. NOTE this `elastic` mode differs from the wall's, which holds a plateau at fc: a
  plateau is elastic-perfectly-plastic, and on a specimen whose subject IS compression the faithful
  reading is unbounded, imposed by scaling fc and epsc0 together (the initial tangent 2fc/epsc0 is
  invariant under that scaling, so E is preserved).
- **The main finding: Eq. 2.1 balances a CONFINED energy.** Aydin's continuum energy is written for
  a field with the transverse strain restrained (eps_x = e, eps_y = 0), so what the balance pins is
  E/(1-nu^2), not E. Verified numerically: a calibrated elastic lattice reproduces the confined
  modulus to **1.000**. But an unconfined test — which is what a cube between smooth platens is —
  measures `E_confined * (1 - nu_lattice^2)`, using the LATTICE's Poisson ratio. So this cube reads
  **0.87*E**, mesh-independently (0.888 / 0.876 / 0.870 / 0.868 at mesh 20 / 10 / 5 / 4). It is a
  property of the calibration, not a solver artifact, and it narrows to ~0.94*E at horizon 3.01.
- **`nu_consistent` is a different quantity from the lattice's Poisson ratio, and was being misread.**
  `nu_consistent = 1 - 2*C66/C11` answers "at what nu do the normal-strain and shear CALIBRATION
  ROUTES return the same area" (~0.18 at horizon 1.5). The lattice's actual Poisson ratio is
  `C12/C22` = **0.41**. They coincide only for an isotropic lattice, and a uniform-EA horizon lattice
  is **cubic-symmetric, not isotropic**: C66/((C11-C12)/2) = 1.38 at horizon 1.5. There is therefore
  no nu at which it is isotropic at all. Confirmed three ways — analytic unit-cell algebra, the
  affine energy tensor, and a direct numerical solve (nu measured 0.4108 vs 0.4082 predicted; the
  displacement field is affine to 0.19%, so the affine assumption itself is exact here).
  `EnergyBalanceResult` gains `nu_effective` and `cubic_anisotropy` so the distinction is reported
  rather than inferred. **This contradicts wording in CLAUDE.md and `examples/wall/build.py`, which
  describe the horizon=1.5 lattice as "isotropic at nu ~ 0.18, close to concrete's 0.20". The wall's
  own results are not invalidated — flexure is governed by the C11 behaviour the balance matches
  exactly — but the justification as written is wrong. Superseding text is left to a later entry
  rather than rewritten here.**
- **Aydin's horizon 3.01 recommendation gets a quantitative reason he did not give:** it is markedly
  closer to isotropic (cubic anisotropy 1.38 -> 0.92, nu_effective 0.41 -> 0.31). But at mesh 20 the
  horizon reaches 30% of the specimen, so homogenization is invalid there and the tangent check
  catches it (ratio 0.905, improving to 0.947 at mesh 10). The check earns its place.
- **Results (mesh 20, horizon 1.5, rate 8.25 mm/s, ~4 s per run).** `crushing`: peak **19.0 MPa =
  0.63*fc** at 0.995*epsc0, 80% degradation by 1% strain, 16.7% of struts cracked but only 2.6%
  crushed — i.e. the lattice's compressive strength is governed by transverse tensile SPLITTING, not
  by fc. That is the mechanism Aydin points at ("compression failure in concrete is also related
  with cracking in different directions") without modelling. `elastic`: **no failure at all** —
  stress rises monotonically to 180 MPa (6*fc) at 1% strain. Splitting occurs (the same 16.7%) but a
  small-displacement `Truss` cannot BUCKLE, so the split columns keep carrying. Under `--corot` the
  same run peaks at 152 MPa and then collapses 99.9%, confirming the missing mechanism is geometric.
  **Conclusion: Aydin's elastic-compression assumption leaves the model with no material compressive
  failure mechanism whatsoever.**
- **Alternatives rejected:** extending `examples/cube/` (that is the 100 mm TENSION study with a
  beam-column reference baked into its specimen and a staged rebuild, D42-D46 — different specimen,
  different calibration, no reference). Rough/bonded platens as the default (they form confinement
  cones and read above fc, conflating platen friction with material strength; smooth is the
  low-friction bound of the van Vliet & van Mier study Aydin cites in Fig. 1.10).
- **Backend changes.** `run_pushover_dynamic` accepted only POSITIVE targets — the drive rate was
  `abs(target)/...`, `ops.sp` imposed a positive ramp, and the stop test `cdisp() >= target` fired
  after one step for a negative target; it now carries the sign on the constraint and compares
  magnitudes, matching `run_pushover`. It also gains `capture=True`, returning the full nodal
  displacement field at the last step and at peak |base shear| (off by default), which is what the
  damage figure consumes. New `viz.strut_strains` + `viz.figure_damage` (pure geometry, no `ops.*`):
  per-stage rows, strain field beside a discrete crack/crush map. The strain field uses a SYMLOG
  colour scale with its linear threshold at the cracking strain — a localized crack sits 2-3 decades
  above the elastic field, so a linear scale renders everything but the worst strut flat grey.
- **Status:** accepted. New `examples/compression_cube/{specimen,build,pushover}.py`; additive fields
  on `EnergyBalanceResult`; `run_pushover_dynamic` sign + capture; two new `viz` functions. Suite: 35
  pass, 1 pre-existing known failure (`test_rc.py::test_nonlinear_pushover_runs_and_yields`, D34).

### D54 — 2026-08-17 — Why the compression cube never reaches fc: the calibrated area fixes STRENGTH by section accounting, and transverse cracking sheds the diagonals
- **Context (user):** "why is capacity not achieved?" — the D53 cube peaks at 19.0 MPa against a
  grade fc of 30 MPa. Plus: provide a force-deformation curve alongside the stress-strain one.
- **Answer, verified by force decomposition rather than argued.** Cut the cube horizontally between
  two node rows. Exactly one vertical strut crosses per node COLUMN, so `n = L/mesh + 1 = 11` struts
  of the calibrated area present themselves against the gross face: `11 * 2305.7 / 40000` = **0.634**
  of the section. Measured peak / fc = 19.02/30 = **0.634**. The peak occurs at 0.995*epsc0, i.e.
  exactly when the vertical struts reach fc. `peak / (verticals alone at fc)` = **1.000**.
- **The diagonals should help, and do — until cracking.** Under an affine field a diagonal's strain is
  `-eps*(1-nu)/2`, so at nu = 0.39 it strains at only ~30% of the vertical struts'. Elastically the
  split across the cut measures **72% vertical / 28% diagonal** (hand calculation from the unit cell
  predicted 72/28 independently). But the transverse struts crack at an axial strain of only
  `eps_cr/nu = 2.6e-4` = **12.8% of epsc0**, and as the block splits its apparent Poisson ratio rises
  toward and past 1, driving the diagonal strain to zero and then into tension. By peak the diagonals
  carry **0.0%**. Same split in both compression modes (72/28 and 73/27 pre-crack), so it is a
  property of the lattice topology, not of the compression law.
- **The general point, and it is the same defect the tension study found with the opposite sign.**
  Aydin's balance sizes the strut area to match ENERGY (stiffness). The same area then silently fixes
  STRENGTH through the section accounting above, and nothing in the method constrains that to fc.
  Tension overshoots by 1.50x because the diagonals add a parallel path that stays intact (D40/D41);
  compression undershoots to 0.63x because the diagonals shed. **A strength calibration is a separate
  knob the method does not have** — D41's `peak_correction="strength"` (scale ft, leaving area and
  hence K0 alone) is the pattern if one is wanted here.
- **Numerical honesty note.** The first version of the decomposition took `abs()` of each group before
  summing, which makes the cut-vs-reaction reconciliation meaningless and reported 1017% error. Groups
  are now summed SIGNED and the reconciliation is checked only where the load exceeds 5% of peak
  (near zero load the ratio is start-up transient, not disagreement), with an explicit
  ok / CHECK FAILED verdict printed rather than left for a reader to evaluate. Likewise the pre-crack
  sample index is chosen by STRAIN (last step below `eps_cr/nu`), not as a fraction of peak load — a
  load-based threshold sits past cracking in `elastic` mode, where the peak is 6x larger, and
  misreported that split as 94/6.
- **Drive rate lowered to 2.0 mm/s (from 8.25), on user request, and what that settled.** A sweep at
  mesh 20 / horizon 1.5 shows the RESPONSE is already rate-independent: peak stress 19.02 MPa, strain
  at peak 1.99e-3 and the 28% -> 0% diagonal split are unchanged at 8.25 / 4.0 / 2.0 / 0.825 mm/s.
  What a slower drive buys is a cleaner measurement, not a different answer — the reconciliation error
  falls 1.33% -> 0.44% -> 0.16% -> 0.05%, confirming that residual IS the inertial/damping part of the
  dynamic reaction and nothing else. 2.0 chosen as default: tangent converged to 0.02% of its
  slow-rate limit, contamination negligible, 18 s per run, and heavier configurations still tractable.
  The sweep table lives in `specimen.QUASI_STATIC_RATE`. Note the rate does NOT move the horizon-3.01
  tangent ratio (0.906 at every rate) — that is the homogenization-validity problem of D53, cleanly
  separated from rate: rate governs reconciliation, mesh/horizon governs the tangent.
- **Section-accounting bug found by rerunning at horizon 3.01.** `vertical_strut_capacity` assumed one
  vertical strut per node column (`L/mesh + 1`), true only for horizon < 2. At horizon 3.01 vertical
  struts span one, two AND three rows, so **66** cross the cut, not 11, and the old formula reported a
  19.02 MPa reference for a run whose real reference is 20.01 MPa. It now counts crossing struts from
  the MODEL and reads each element's own area, so it is correct for any horizon and for a non-uniform
  `strut_area`. This also overturned the narrative for horizon 3.01: there the diagonals still carry
  **22% at peak** and the peak arrives at 0.48*epsc0, BEFORE the verticals reach fc, so "the capacity
  is the verticals' section share" is a horizon-1.5 result, not a general one. The printed conclusion
  is now branched on the measured diagonal share rather than asserted.
- **`--compression elastic --corot` fails the reconciliation check (25%), correctly.** That run
  collapses 100%, and a full collapse with corotational geometry is genuinely dynamic — inertia is no
  longer negligible and the mid-height cut is no longer horizontal in the deformed configuration. The
  check flagging it is the check working; that configuration's post-peak branch is not a quasi-static
  result and should not be read as one.
- **Backend change.** `run_pushover_dynamic` gains `element_groups`, the same force-decomposition
  probe `run_pushover` has carried since D19 (label -> [(element_id, dof, coef)], recorded per step).
  Its docstring already advertised "vertical vs diagonal struts" as the use case; this is that.
- **New outputs** per run: `*_force_deformation.png` (axial force kN vs shortening mm, with the gross
  section at fc and the verticals-alone-at-fc reference lines) and `*_loadpath.png` (vertical vs
  diagonal share vs strain, with the transverse-cracking onset marked). The saved JSON now carries the
  per-step group forces, so both redraw without re-running.
- **Status:** accepted. Suite: 35 pass, 1 pre-existing known failure (D34).

### D55 — 2026-08-17 — Why the compression cube's capacity is not `fc * A_face`: the lattice has MORE material than the continuum, and loses the load path rather than the material
- **Context (user):** a clear explanation of why the carried axial force is not `fc * A_face`.
- **Supersedes the framing of D54, not its numbers.** D54 answered with the section share — 11 vertical
  struts x the calibrated area = 0.634 of the gross face, and `peak / (verticals at fc)` = 1.000. That
  is correct arithmetic for the FINAL state but it describes the outcome, not the cause, and invites
  the wrong conclusion that the lattice is short of material. It is not. New `build.capacity_accounting`
  computes four capacities over one horizontal cut, and the measured peak is the last of them:
  ```
    continuum: fc x gross face          1200 kN   1.00x    material fills the cut plane
    every crossing strut at fc          1739 kN   1.45x    strut area is 1.79x the face
    compatibility-limited at epsc0      1266 kN   1.05x    diagonals strain (1-nu)/2 of verticals
    verticals alone at fc                761 kN   0.63x    after splitting sheds the inclined path
    MEASURED peak                        761 kN   0.63x
  ```
- **The mechanics, in the order they matter.** A continuum resists `fc*A_face` because material fills
  the cut plane. A lattice resists only what its struts transmit, each contributing
  `sigma*A*cos(theta)`, so the governing section is the STRUTS' area projected on the load direction.
  Three things then separate the lattice from the continuum:
  1. *Area: the lattice has MORE, not less.* Crossing struts total 1.79x the face area, and at fc they
     would carry 1.45x the continuum. Matching STIFFNESS across many directions demands more total
     area than a continuum has, because each strut is only partly aligned with any one load. Same
     over-provisioning the tension study measured as PEAK_RATIO = 1.50 (D40/D41).
  2. *Compatibility, not strength, sets the sharing.* A 45-degree diagonal strains at only
     `(1-nu)/2` ~ 0.30 of the verticals, so when the verticals reach epsc0 the diagonals sit near
     15 MPa, not 30. This alone brings the lattice to **1.05x the continuum** — i.e. an intact lattice
     is very nearly correctly proportioned, which is a non-trivial endorsement of the calibration.
  3. *Splitting removes the inclined load path — the dominant term.* Horizontal struts are pulled into
     tension at `+nu*eps` and crack at an axial strain of only `eps_cr/nu` = 12.8% of epsc0. The block
     splits into columns, apparent nu rises past 1, diagonal strain `-eps(1-nu)/2` passes through zero
     and reverses into tension, and the diagonals crack too. Measured share 28% -> 0%.
     **Cost: 505 kN, 42% of the continuum capacity.**
  So the one-line answer: **the section is not short of material, it is short of a load path.** Remove
  splitting and the lattice lands within 5% of `fc*A_face`; with splitting only the vertical struts
  remain, and they occupy 0.634 of the face.
- **Why this is a real modelling limitation and not merely bookkeeping.** In real concrete the inclined
  and transverse capacity is not lost at 13% of the peak strain — aggregate interlock, the 3rd
  dimension's confinement, and lateral tensile capacity distributed over an area rather than lumped in
  a few bar-like struts all keep the diagonal path alive far longer. A 2D uniform-EA horizon lattice
  with `ft` on every strut sheds it almost immediately. Directions if the strength matters: raise the
  transverse struts' cracking strength specifically (the D41 `ft`-scaling knob, applied per
  orientation), go to horizon 3.01 where the redundancy is greater, or accept the model for stiffness
  and damage LOCATION and not for compressive strength.
- **Horizon dependence.** At horizon 3.01 the accounting differs and the D54 story does not carry: 66
  vertical struts cross the cut, the diagonals still hold 22% at peak, and the peak arrives at
  0.48*epsc0 before the verticals reach fc — failure is progressive rather than governed by one strut
  family. The printed conclusion branches on the measured diagonal share for exactly this reason.
- **Status:** accepted. `build.capacity_accounting` + the printed table in `pushover.py`; no analysis
  change (the runs are unchanged, only their explanation). Suite: 35 pass, 1 known failure (D34).

### D56 — 2026-08-18 — Compression cube simplified to one material (fc = 15.1 MPa) and one load path; `--rigid-platen` added and found to change nothing
- **Context (user):** "simplify the model as much as possible" — standard truss, Concrete02 at
  fc = 15.1 MPa, Aydin's calibration, very slow dynamic loading as the pushover, rollers under the
  cube with one pinned for stability. Plus a CLI flag (default False) constraining the upper face to
  move down by the same amount. "Let's see how it behaves."
- **What was removed.** `--compression elastic` (the Aydin-literal unbounded-compression variant),
  `--corot` (corotational struts), and the `--gf-factor` / `--residual-ratio` / `--periods` knobs.
  One grade, one element type, one drive. `build.py` loses `effective_epsc0` and
  `ELASTIC_COMPRESSION_FACTOR`; `secant_factor` drops its epsc0 argument. Supports were already
  rollers-plus-one-pin (D53), so that part needed no change.
- **The material, and one trap in it.** fc = 15.1 MPa with E = 16.1 GPa and ft = 1.5 MPa — the
  SW-NC-FF test-unit grade, so this example and the wall now share a concrete. **epsc0 is DERIVED as
  `2*fc/E` = 1.876e-3, not quoted from the paper's 0.0022.** Concrete02's initial compressive tangent
  is `2*fc/epsc0` regardless of what the grade's `E` field says, so quoting 0.0022 would have given
  the struts a 13.7 GPa tangent while Aydin's balance was handed 16.1 GPa — a silent 15% error in
  precisely the quantity the calibration exists to reproduce. With epsc0 derived, the measured
  tangent lands at 1.008 of the closed-form expectation.
- **Everything scales with fc exactly, confirming the D55 result is geometric.** Against the previous
  fc = 30 MPa run: peak 9.57 MPa = **0.634 fc** (was 19.02 = 0.634 fc), strain at peak 0.998*epsc0,
  pre-crack load split 71.9% / 28.1% vertical/diagonal going to 100% / 0% at peak, transverse
  cracking at 12.7% of epsc0, capacity accounting 1.00 / 1.45 / 1.05 / 0.63 x continuum. Identical
  ratios throughout — the 0.634 is a property of the lattice geometry and the calibrated area, not of
  the strength.
- **Rate 0.5 mm/s ("very slow").** Reconciliation between the cut force and the base reaction falls to
  **0.05%**, from 1.33% at the original 8.25 mm/s — i.e. the inertial and damping content of the
  recorded reaction is now negligible. ~52 s per run at mesh 20.
- **`--rigid-platen` changes NOTHING, and that is the honest result.** The flag ties every top-face
  node's vertical DOF to the control node via `ops.equalDOF` and drives only that node, so the face
  is forced to stay flat. Measured against the default:
  peak 9.57 MPa both, strain at peak 1.87e-3 both, degradation 73.5% both, damage 55/420 cracked and
  11 crushed both, load split identical, initial tangent 14,315 vs 14,321 MPa (0.04% apart).
  The reason is that the DEFAULT already prescribes the same ramp to every top node individually, so
  the top face was flat all along; the constraint imposes as an MPC what was already imposed as N
  identical single-point constraints. The only differences are numerical: T1 drops 0.827 -> 0.720 ms
  (equalDOF removes DOFs, stiffening the fundamental mode) and the run takes 83 s instead of 52 s,
  since `dt = T1/30` shrinks with it.
  A `top_face_uy_spread` diagnostic is now printed and saved so this equivalence is measured rather
  than assumed, instead of the flag quietly appearing to do something. It lands slightly the wrong
  way round, which is worth knowing: the DEFAULT is flat to **0.000e+00 mm** (identical prescribed
  motions are exact), while `--rigid-platen` leaves **1.2e-05 mm** of residual spread — the
  `Transformation` handler enforces the tie to solver tolerance, not exactly. Both are flat to any
  meaningful precision; neither is 'more rigid' than the other in practice.
- **What the flag would need to be interesting.** A rigid platen only differs from a free one when the
  free case lets the face WARP — i.e. when the top is loaded by a uniform traction rather than by a
  prescribed displacement. That is a different runner contract (`run_pushover_dynamic` imposes
  displacement, not load), and force control cannot trace the post-peak branch, which is the whole
  point of the pushover. Left as-is with the equivalence documented.
- **Backend change.** `run_pushover_dynamic` gains `equal_dof` — a list of
  `(retained, constrained, dof)` emitted as `ops.equalDOF` after `wipeAnalysis` and BEFORE the eigen
  solve, so tied DOFs are reflected in T1 and the Rayleigh damping. Requires the `Transformation`
  handler, which the runner already sets.
- **Status:** accepted. Suite: 35 pass, 1 pre-existing known failure (D34).

### D57 — 2026-08-18 — Compression cube gains a STATIC pushover alongside the dynamic one; it stalls before the peak, and it makes `--rigid-platen` matter
- **Context (user):** add a nonlinear static pushover as an alternative to the dynamic relaxation.
- **Implementation.** `--solver dynamic|static`, default dynamic. Static uses `run_pushover`:
  DisplacementControl on the control node scaling a uniform downward reference traction over the top
  face (`specimen.axial_loads`), with `ModifiedNewton -initial` as the primary algorithm — iterating
  on the constant elastic tangent never re-forms the singular tangent a cracked strut band produces
  (D44). `run_pushover` gains `capture` and `equal_dof` so both solvers feed the SAME diagnostics
  (damage figure, load-path decomposition, platen tie) and are interchangeable for a caller.
- **Result: the static solver stalls well before the peak, still ascending.** Free platen dies at
  0.152% strain (81% of epsc0) having reached 9.65 MPa; rigid platen at 0.174% (92% of epsc0) at
  9.85 MPa. The dynamic run traces the whole curve to 1% strain, peaking at 9.57 MPa at 0.998*epsc0
  and falling to a residual plateau. This is the D38/D44/D45/D46 finding again on a new specimen: a
  cracking lattice's softening is a MECHANISM, not a snap-back, and a static path-follower cannot
  cross it however good its algorithm. Reporting now labels a stalled maximum sitting at the last
  step as "max before stall" with an explicit "STILL ASCENDING ... this is a lower bound on the
  capacity, not the capacity" warning, because 9.65 > the dynamic 9.57 and would otherwise read as a
  higher strength when it is simply a curve that never turned over.
- **The unexpected payoff: static makes `--rigid-platen` a real variable.** D56 found the flag inert,
  because the dynamic solver PRESCRIBES the same downward ramp at every top node, so the face is flat
  with or without a tie. The static solver applies a uniform TRACTION instead, which lets the face
  warp — measured `top_face_uy_spread` = 3.5e-02 mm, **11.5% of the imposed shortening**, against
  0.000e+00 with the tie. That is exactly the traction boundary D56 said the flag would need, and it
  arrived free with the static solver rather than needing a new runner contract.
  With the face held flat the block is slightly stronger and lasts longer before the solver dies:
  9.85 vs 9.65 MPa (+2%), stall at 0.174% vs 0.152% strain, 104 vs 92 struts cracked, and the
  diagonals retain 3.6% vs 4.9% of the load at the stall point. Physically sensible — a rigid platen
  redistributes load off the softening columns onto the stiffer ones, delaying localization.
- **Which solver to believe.** Dynamic, for anything at or past the peak: it is the only one that
  reaches it. Static is worth having as the cheap check — 0.3 s against 52 s, no rate to calibrate —
  and the two agree on the elastic branch (initial tangent 1.013 / 1.006 of the closed-form
  expectation, against the dynamic's 1.008), which is the part a static solver can be trusted for.
- **Status:** accepted. `run_pushover` gains `capture` + `equal_dof`; `specimen.axial_loads` added;
  `--solver`/`--steps` in the CLI; static outputs land under a `_static` stem so they never overwrite
  the dynamic ones. Suite: 35 pass, 1 pre-existing known failure (D34).

### D58 — 2026-08-18 — Solution algorithm exposed on BOTH cube solvers; benchmarked, and it barely matters — the static stall is a mechanism, not an algorithm failure
- **Context (user):** what are the solution-algorithm options for the cube's analysis types?
- **Before:** the static pass hard-coded `ModifiedNewton -initial` and `run_pushover_dynamic`
  hard-coded plain `Newton` with no way to change it. Now both take an `algorithm` tuple
  (`run_pushover` already did, D44), surfaced as `--algorithm`.
- **All 19 candidates are accepted by this `openseespymac` build** — Linear, Newton (+`-initial`,
  `-initialThenCurrent`), ModifiedNewton (+`-initial`), KrylovNewton (+`-iterate initial`,
  `-increment initial`), SecantNewton, BFGS, Broyden, NewtonLineSearch (+`-type` Bisection / Secant /
  RegulaFalsi / InitialInterpolated), ExpressNewton, PeriodicNewton. None raised.
- **Static benchmark (mesh 20, 500 steps to 1% strain) — the substantive result.** Fifteen of the
  nineteen stall at **exactly 0.1520% strain** at 9.64-9.65 MPa: Newton, Newton -initial,
  Newton -initialThenCurrent, ModifiedNewton -initial, all three KrylovNewton variants, SecantNewton,
  BFGS, Broyden, all four NewtonLineSearch variants, PeriodicNewton. Identical stall point across
  every iteration strategy is the evidence that **the barrier is a MECHANISM — a genuinely singular
  tangent — and not a convergence failure any algorithm can out-iterate** (D38/D46 again). Only the
  outliers differ, and each for a legible reason: `ModifiedNewton` on the CURRENT tangent dies far
  earlier (0.0285%, 3.63 MPa) because it re-forms and then reuses the cracked tangent;
  `ExpressNewton` similarly (0.0300%); `Linear` runs slightly further (0.1660%, 9.71 MPa) only
  because it never iterates, so its numbers are not equilibrium solutions and must not be read as
  strength. Cost is flat (0.2-0.3 s) except BFGS/Broyden at 0.5-0.6 s.
- **Dynamic benchmark (to 0.3% strain):** every algorithm completes — Newton 16.3 s, KrylovNewton
  **12.2 s**, ModifiedNewton -initial 14.5 s, NewtonLineSearch 16.4 s, all converged with the same
  response. Mass and damping keep the effective tangent well-conditioned, so the choice is purely a
  speed question here; KrylovNewton is ~25% faster than the Newton default and is the one to reach
  for on a long run.
- **Per-solver defaults, since the two want different things.** `Newton` for dynamic (no need for the
  initial tangent when the effective tangent is well-conditioned); `ModifiedNewton -initial` for
  static (never re-forms the singular cracked tangent). The retry ladder — primary -> KrylovNewton ->
  NewtonLineSearch — is unchanged in both and is printed alongside the chosen primary, since it does
  much of the work and would otherwise make the primary look more decisive than it is.
- **CLI trap worth recording.** `--algorithm` takes ONE QUOTED STRING, split in the example, not
  `nargs="+"`. OpenSees algorithm options begin with a dash (`-initial`, `-type Bisection`) and
  argparse claims those as flags of its own, so `--algorithm ModifiedNewton -initial` fails outright.
  The first implementation used `nargs="+"` and that exact invocation silently produced no run inside
  a benchmark loop; caught because the comparison table came back with a missing row.
- **Status:** accepted. Suite: 35 pass, 1 pre-existing known failure (D34).

### D59 — 2026-08-18 — `--max-iter` added; it exposes that the static stall is a property of the BOUNDARY CONDITION, not of the lattice — partially superseding D57/D58
- **Context (user):** run the same analysis with max_iter = 1000 per step. Both runners took a
  `max_iter` (100 static, 50 dynamic) that the example never surfaced; now `--max-iter`.
- **The result overturns a claim I made in D57/D58.** Those entries said the static solver's stall was
  a MECHANISM that no algorithm could cross, generalising from the free-platen case. A max_iter sweep
  shows that is only half true, and which half depends entirely on the upper boundary condition:
  ```
    FREE platen (uniform traction)      RIGID platen (--rigid-platen, tied face)
      100    -> stalls 0.1520%            100  -> stalls 0.1740%
      1000   -> stalls 0.1520%            200  -> stalls 0.3077%
      5000   -> stalls 0.1520%            500  -> reaches 1.0000%, CONVERGED
      20000  -> stalls 0.1520%            1000 -> reaches 1.0000%, CONVERGED
  ```
  A **200x** iteration budget moves the free-platen stall by nothing at all — same strain, same
  9.65 MPa, same 92 cracked struts, only 65x the runtime (0.2 s -> 13.3 s burning iterations before
  giving up). That is a singular tangent and the D38/D46 reading holds there. But the tied face
  simply needs a bigger budget: at >= 500 iterations the static solver traces the FULL curve,
  including the descending branch, in 0.3-5 s.
- **Why the boundary condition decides it.** Under a uniform traction the top face is free to warp
  (measured spread 11.5% of the shortening, D57), so when one column starts softening nothing stops
  it shortening unopposed — a local mechanism, and the tangent goes singular. Tying the face forces
  the load to redistribute onto neighbouring columns instead, which keeps the tangent non-singular;
  what remains is merely slow convergence near the peak.
- **The static full trace agrees with the dynamic one**, which is the check that it is real and not a
  numerical artefact: peak 9.85 vs 9.57 MPa (+3%), post-peak degradation 74.3% vs 73.5%,
  cut-vs-reaction reconciliation 0.01%, load path 71.4/28.6 pre-crack going to 96.4/3.6 at peak.
- **A counterintuitive cost note.** max_iter 500 takes 5.0 s while max_iter 1000 takes 0.3 s on the
  same problem. More headroom is FASTER because steps that would hit a 500 limit instead converge
  directly, avoiding the retry ladder's sub-stepping (5/20/50 divisions across three algorithms),
  which is far more expensive than the extra iterations. Raising max_iter is close to free when it
  works and merely slow when it does not, so it is a cheap first thing to try.
- **Practical guidance now printed with every run** (algorithm + max_iter on one line): for a static
  pass use `--rigid-platen --max-iter 1000`; without the tie no iteration budget will finish, and the
  dynamic solver is the only option.
- **Status:** accepted. Supersedes the over-general "no algorithm can cross it" wording of D57/D58 —
  that statement is correct for a traction boundary and wrong for a tied one. Suite: 35 pass, 1
  pre-existing known failure (D34).
### D60 — 2026-08-18 — `examples/aydin_cube/`: Aydin, Binici & Tuncay (2021) reproduced faithfully — a TENSION-ONLY lattice that predicts compressive strength, and the full difference list against `compression_cube`
- **Context (user):** create an `aydin_cube` example, stay faithful to the paper, and list every
  major difference between it and `examples/compression_cube/`.
- **The paper is a different source from D47's.** D47-D59 calibrate against Aydin's 2017 METU MSc
  thesis (OLM), which states compression is out of scope (Sec. 2.2 p. 30) — the gap `compression_cube`
  exists to price. This is Aydin BB, Binici B & Tuncay K (2021), *Lattice simulation of concrete
  compressive behaviour as indirect tension failure*, Mag. Concr. Res. 73(8): 394-409, which closes
  that gap. Same elastic energy balance (its p. 396); everything about compression is new.
- **His mechanism, in the order it matters.** (1) Struts have **no compressive strength at all** —
  linear elastic forever, optionally with a reduced-stiffness (RSM) knee at `epsc0/3` dropping to
  `0.4*Ec`. (2) Transverse struts crack in indirect tension and the block splits into columns.
  (3) The split columns lose STABILITY, and that event is the compressive strength. (4) It only
  works on a **perturbed grid**: each node moved by `R <= Rmax` at a random angle, with `Rmax/d` the
  single calibration parameter, fitted by his Fig. 4 loop until the mean of 5 runs is within 10% of
  the target `fc`.
- **His own control case is `compression_cube`, and D53 had already reproduced it.** His Fig. 3(a),
  uniform grid: NR climbs to 45 MPa against `fc = 20` with *no strength decrease* (that is merely
  where he stopped), LR fails at 43 MPa ~ 2.15*fc by diagonal shear. He concludes "uniform grid
  lattice simulations are incapable of simulating the compression response" and names it **vertical
  locking**. D53's `--compression elastic` run is that NR case independently: monotone to 6*fc with
  no failure. Measured here at `--rmax 0`: **1.51 x fc** (corotational; the peak is real, at
  0.36% strain), against his ~2 x fc.
- **Results — the reproduction holds at the paper's own fitted `Rmax/d`, with no refitting.**
  100 mm cube, d = 10 mm, `fc = 20` MPa and `Gf = 50` N/m chosen, `ft = 0.35*sqrt(fc)` and
  `Ec = 4700*sqrt(fc)` derived, 5 realizations, ~30 s each:
  ```
    LR, Rmax/d = 0.075 (his value):  mean 18.48 MPa = 0.924 x fc   CoV  6.7%   (16.98 - 20.22)
    NR, Rmax/d = 0.060 (his value):  mean 20.62 MPa = 1.031 x fc   CoV 10.7%   (17.03 - 22.91)
    uniform grid, Rmax/d = 0:             30.12 MPa = 1.506 x fc   — locks, as he reports
  ```
  His Fig. 4 acceptance criterion is 10%; **both boundary conditions pass at his own fitted
  `Rmax/d`, with no refitting** (-7.6% and +3.1%). LR peaks at exactly **1.00 x epsc0**, which is
  what the RSM knee exists to deliver (his Fig. 3(b): without RSM the strength is right and the
  failure displacement is not); NR peaks later, at 1.53 x epsc0. LR softens to 0.32 of peak by
  0.4% strain — a genuine post-peak branch, not a material residual plateau. The scatter is
  his too (Table 2 spans 42.9-66.2 MPa on one parameter set), which is why 5 realizations is the
  default and a single run is labelled as one sample.
- **His Fig. 4 calibration loop works, and finds his answer independently.** `--calibrate` searches
  `Rmax/d` until the mean of `n` realizations is within 10% of `fc`. Started deliberately far off at
  0.020 (NR, 3 realizations): **0.020 -> 1.322 fc, 0.026 -> 1.319 fc, 0.053 -> 1.044 fc, converged**
  — landing on `Rmax/d = 0.053` against the value he fitted for this specimen, **0.060**. The search
  needs care that a plain secant does not survive: the strength estimate is the mean of a few random
  meshes, so when two rounds differ by less than the scatter the apparent slope is noise and the
  secant divides by it (the first version extrapolated to `Rmax/d = 0.30`, which folds the mesh).
  It now falls back to a proportional step below a noise-scaled slope threshold, damps to at most a
  factor of two per round, and clamps to [0.005, 0.15] — above ~0.15 the perturbation approaches the
  node spacing and neighbours start swapping places.
- **Three implementation findings that the paper does not state and that cost real time.**
  1. **`ElasticMultiLinear`, not `HystereticSM`.** Both reproduce his Fig. 1 envelope exactly, but
     `HystereticSM` unloads on the initial stiffness, so a fully cracked strut reads **-10 MPa
     (RSM) / -44 MPa (no RSM)** when its strain returns to zero — it would shove the split columns
     apart and destroy the mechanism. `ElasticMultiLinear` is path-independent and returns exactly
     0. The price is the opposite error, no damage memory: a softened strut recovers on reload
     where his sequentially linear analysis holds the reduced secant. Acceptable under monotonic
     compression, **unusable for cyclic work**.
  2. **`corotTruss` is mandatory, and he never says so.** The failure mode is geometric. His
     explicit scheme integrates the nodes' actual positions and gets it free; a small-displacement
     `Truss` cannot express it (D53: 6*fc, no failure). `--truss` keeps the negative control.
  3. **Strut topology must be frozen on the UNPERTURBED grid.** Re-running `connect_horizon` on
     moved nodes drops the struts that happen to cross `horizon*mesh_size` — 5 of 420 at
     `Rmax/d = 0.08` — making connectivity a second random variable on top of geometry. His
     Fig. 2(b) still shows both diagonals in every cell, so keeping the cells is the faithful
     reading. `build_lattice_rc` gains `grid=` and `pairs=` overrides for this.
- **`a2`/`a3` are SOLVED from `Gf`, not tabulated.** His Table 1 fits them to the Cornelissen et al.
  (1986) curve per specimen. Requiring `L * integral(sigma d eps) == Gf` instead gives them in
  closed form per strut (`materials.aydin_lattice_softening`), which is what keeps the model
  mesh-objective at spacings he never ran — and it reproduces his tabulated values to **12%**
  (Jansen & Shah: 52.9/264 against his 60/300). Orthogonal and diagonal struts get different tails,
  exactly as `concrete_uniaxial_regularized` steepens `Ets` with length (D20).
- **The governing difference from `compression_cube`, which decides what can be combined.**
  `compression_cube` makes `fc` an INPUT (Concrete02 caps every strut) and measures **0.634 fc**;
  Aydin makes it an OUTPUT. His Figs. 5(a)-(b) are a monotone map — more perturbation, less
  strength — descending onto `fc` from an unbounded ceiling. `compression_cube` starts BELOW the
  target at 0.634, so **no `Rmax/d` reaches `fc`**: bolting perturbation onto Concrete02 would
  double-count compression, once as a strut strength and once as an emergent instability. The two
  ingredients are not separable. This does not invalidate D54/D55 — their arithmetic stands — but it
  answers the question they left open: being short of a load path IS the point, and the residual
  error was putting a compressive strength on the struts at all.
- **Documented deviations from the paper** (also listed in `examples/aydin_cube/DIFFERENCES.md`):
  RSM applied to every strut rather than only laterally cracked ones (state-dependent, not
  expressible as a static assignment; under uniaxial compression nearly every compression strut
  does crack laterally before peak); no damage memory (above); boundary nodes slide only ALONG
  their own face and corners are pinned, so the faces stay flat and box selection stays exact (the
  paper does not say what it does); one extra `ux` restraint under NR, because his explicit scheme
  tolerates a rigid-body mode the implicit Newmark solve cannot; `nu = 0.20` assumed for the
  balance, which he does not state.
- **The paper contradicts itself on `alpha`/`beta`.** Its notation list and Fig. 1's caption make
  `alpha` the strain multiplier (0.33) and `beta` the modulus multiplier (0.4); the p. 397 text says
  "alpha and beta were taken as 0.4 and 1/3" and then describes the opposite in the same sentence.
  The self-consistent reading — notation list, caption, and the text's own "40% stiffness reduction
  at a strain of 1/3 of epsc0" — is what is implemented.
- **Status:** accepted. New `examples/aydin_cube/{specimen,build,compression}.py` + `DIFFERENCES.md`;
  new `mesh.perturb_nodes`, `materials.aydin_lattice_softening`, `materials.concrete_lattice_aydin`
  (exported); `build_lattice_rc` gains `grid=`/`pairs=`. New `tests/test_aydin_cube.py` (6 tests).
  Suite: 41 pass, 1 pre-existing known failure (D34).

### D61 — 2026-08-19 — `examples/compression_cube/aydin_approach.py`: OUR cube, Aydin's approach — a controlled swap that isolates the constitutive law from the specimen
- **Context (user):** make the compression cube's modelling and analysis approach similar to Aydin's,
  keeping our cube intact.
- **What "intact" was taken to mean, literally.** `specimen.py`, `build.py` and `pushover.py` are
  **not modified**. The new script imports `specimen.py` unchanged and reuses `build.expected_modulus`,
  so the 200 mm cube, the fc = 15.1 MPa grade with its measured `E` and `ft`, the 20 mm grid, the
  smooth platens (rollers + one ux pin) and the energy-balance strut area are all the same objects
  the D53-D59 runs used. Every D53-D59 number stays reproducible by the command that produced it.
  The only concession: `perturb_nodes` holds the two mid-width nodes on the bottom and top faces
  fixed, because `specimen.py` finds them by an exact box at `x = L/2` (the ux restraint and the
  control node). Two of 121 nodes, and every top node is driven by the same ramp regardless.
- **The six choices swapped, together.** Concrete02 -> tension-only `ElasticMultiLinear` (no
  compressive strength, RSM knee at `epsc0/3`); bilinear -> trilinear tension tail solved from `Gf`;
  uniform -> perturbed grid; `Truss` -> `corotTruss`; deterministic -> 5 random realizations;
  `fc` as an INPUT -> `fc` as an OUTPUT with `Rmax/d` calibrated to it. D60 established they are not
  separable, which is why the script swaps all six rather than offering them as flags.
- **`Rmax/d` calibrated for THIS cube, not borrowed.** His fitted values (0.060 NR / 0.075 LR) are
  for a 100 mm cube at fc = 20, Gf = 50, d = 10 and were not assumed. The Fig. 4 loop ran
  **0.050 -> 1.235 fc, 0.062 -> 1.171 fc, 0.093 -> 0.765 fc, 0.075 -> 0.919 fc, converged**
  (-8.1%, inside his 10% criterion). The production set of 5 at the settled 0.075 gives
  **14.23 MPa = 0.943*fc, -5.7%, CoV 11.3%** (peaks 11.77 / 14.89 / 16.14 / 13.81 / 14.54). That the
  search lands on top of his LR figure is a coincidence of two specimens differing in size, grade,
  grid spacing and end condition.
- **Results on our cube, and the two things the controlled swap proves.**
  ```
                                   Concrete02 (pushover.py)   Aydin approach (this script)
    initial tangent / expected     1.008                      0.997  (0.879 of grade E)
    uniform grid (Rmax/d = 0)      0.634 fc — a CAPACITY      >=1.86 fc, STILL RISING — it LOCKS
    calibrated (Rmax/d = 0.075)    n/a                        0.943 fc  (5 runs, CoV 11.3%)
    strain at peak                 0.998*epsc0                1.36*epsc0
    pre-peak shape                 parabolic (Concrete02's)   near-linear + the RSM knee
    post-peak                      softening to a 0.2 fc      total collapse to ~0
                                     residual plateau
    inclined path at peak          28% -> 0%                  0.4% (perturbed) / 0.6% (uniform)
  ```
  1. **The elastic calibration is genuinely shared.** Both approaches hit the same closed-form
     modulus target on the same cube (1.008 and 0.996), and the Aydin run independently reproduces
     the D53 result that Eq. 2.1 pins the CONFINED modulus, reading **0.879*E** unconfined. So every
     strength difference is a difference of APPROACH, not of two calibrations talking past each other.
  2. **The load-path collapse is topological, not constitutive.** D54/D55 attributed the inclined
     struts' 28% -> 0% collapse to splitting on a uniform horizon lattice. Replacing the entire
     constitutive law reproduces it (0.4-0.6% at peak), which **tests** what D55 argued.
- **A numerical finding that had to be fixed before any of the above was trustworthy.** One
  realization diverged locally, recorded a **139 MPa** peak on a lattice whose struts have no
  compressive strength at all, still reported `converged = True`, and turned a 13 MPa mean into
  56 MPa — sending the `Rmax` search to its upper bound. New `first_divergence` truncates a record
  at the first physically impossible upward jump: a quasi-static run's largest genuine per-step rise
  is ~0.0008 MPa against a 0.755 MPa (= 0.05*fc) threshold, so the test is not close. Realizations
  that never turned over before truncation are DISCARDED rather than averaged — a diverged solve is
  a failed measurement, not a sample. **Locking and divergence are reported separately**: at
  `Rmax/d = 0` every realization legitimately rises to the end of the run, and that is the result,
  not a failure.
- **Why our cube is harder than his.** Our 20 mm grid at Gf = 53 N/m gives each strut a far shorter
  softening tail than his 10 mm grid — `a3 = 120` against `269`, so our struts fail at 1.1% strain
  where his last to 2.0%. Brittler struts make each failure a more violent event for the transient
  solver, which is why divergence appears here and not in `examples/aydin_cube/`. The remedies, in
  order of preference: refine the grid, lower `--rate`, raise `--damping`.
- **Status:** accepted. New `examples/compression_cube/aydin_approach.py` only; no existing example
  file touched. `tests/test_aydin_cube.py` gains the divergence-filter test (7 tests).
  `examples/aydin_cube/DIFFERENCES.md` gains section 10, the controlled comparison.

### D62 — 2026-08-19 — Wall example: the dynamic drive made measurable, and both knobs lowered — damping 0.8 -> 0.2 and rate 7.6 -> 2.0 mm/s; EA reported per zone
- **Context (user):** on `examples/wall/`, "lower damping to something meaningful", "lower the
  loading rate", and "report EA calculated for concrete lattices".
- **A measurement first, because neither knob could be set honestly without one.** The recorded base
  shear is a sum of reactions, so in a transient solve it carries the inertia and damping of
  everything above the base — but nothing in the repo said how much. Global equilibrium in the drive
  direction supplies it exactly: element resisting forces cancel over the whole domain and the
  horizontal external load is zero (gravity is vertical and held), so `S_base + S_drive =
  -sum_free (M.a + C.v)` where each `S` is the `-sum(reaction)` over that node set. Whatever the two
  CONSTRAINED sets fail to balance IS the dynamic contamination. New opt-in `quasi_static=True` on
  `run_pushover_dynamic` and `run_cyclic_dynamic` returns it as `"dynamic"`; it costs one extra
  reaction sum over the ~21 drive nodes per step, reusing the `ops.reactions()` the base shear
  already triggered. It assumes base + drive cover every constrained DOF in that direction, which
  holds whenever the supports are the base row.
- **Validated before it was trusted**, on an elastic lattice driven slowly, where the answer is
  known: the residual falls linearly with rate and goes to zero (Newmark: 21.9 kN at 7.6 mm/s ->
  2.7 kN at 1.0), and the recorded shear converges on the static reference of 80.49 kN
  (K = 73.17 kN/mm, matching D47's 73.05).
- **Result 1 — the old settings were never as quasi-static as assumed.** At the previous 0.8 / 7.6,
  the 0.3%-drift pushover carries **23.0 kN of solver in a 223.1 kN peak — 10.3%**.
- **Result 2 — the contamination is CRACK-RELEASE RINGING, not viscous drag, and that inverts the
  obvious move.** The steady `C.v` drag is only ~2.5 kN by hand at zeta = 0.8; the measured 23 kN is
  the ringing released every time a strut cracks. Damping is what suppresses that, so the two knobs
  do NOT trade off freely. Dropping zeta to 0.05 alone made it WORSE (12.4%), and dropping the rate
  as well made it worse again (**20.5%**, peak inflated 5% to 234.7 kN) — each crack rings for ~20
  periods into a record that now has 3.8x more steps.
  ```
    zeta   rate   peak kN   start-up ring   steady resid   share   steps
    0.8    7.6     223.1        22.2 kN        23.0 kN     10.3%    1086   <- was
    0.05   7.6     229.3        22.2           28.4        12.4%    1086
    0.05   2.0     234.7         5.2           48.1        20.5%    4128
    0.2    2.0     221.7         5.3           10.2         4.6%    4128   <- is
    0.05   1.0     226.8         2.4           10.6         4.7%    8257
  ```
- **Settled: `DAMPING_RATIO = 0.2`, `QUASI_STATIC_RATE = 2.0` (`specimen.py`, shared by pushover and
  cyclic so the two stay comparable).** Damping is 4x lower and now has a physical reading — 10-20%
  equivalent viscous damping is what a heavily cracked RC wall exhibits near capacity, which is the
  state this specimen spends the protocol in — rather than being a near-critical relaxation fudge.
  The rate is 3.8x lower, quartering the start-up and reversal transients (22.2 -> 5.2 kN).
  Contamination is **less than half** the old value. The two heavily damped runs agree on the peak
  to 0.7% (221.7 vs 223.1) where the zeta = 0.05 runs scatter over 8 kN, which is the signature of a
  converged answer. `zeta = 0.05, rate 1.0` reaches the same 4.6% but costs TWICE the run time.
  Price of the slower drive: the full 4% protocol goes from ~0.5M to ~1.8M steps, ~7 h to ~28 h.
- **A trap removed on the way.** `--rate` was DEFAULTED TO NONE, deriving the speed from
  `--periods 48`, which scales with protocol amplitude: the same setting drove the 1% protocol at
  ~20 mm/s and the 4% one at ~79 mm/s, an order of magnitude above the 7.6 that was ever checked.
  The published runs passed `--rate 7.6` explicitly, so they are unaffected, but the default path
  was not the verified one. `--rate` now defaults to the constant and `--periods` is the discouraged
  fallback, reached only via `--rate 0`.
- **The residual is RELATIVE, not absolute — and the reason is worth recording.** It is the
  instantaneous unbalance at the COMMITTED state, and `run_cyclic_dynamic`'s HHT(0.7) leaves a large
  numerical-dissipation term there that does NOT bias the recorded shear. On the elastic wall that
  runner reads a residual of ~2x the base shear while reproducing the static answer to **0.24%**,
  where Newmark reads 27% for a **2.6%** error — a 60x difference in the residual-to-error ratio
  between integrators. So the number compares runs of the SAME runner (it is linear in rate, and
  falls with damping once cracking starts); it is not a portable error estimate. The initial
  implementation printed an absolute "NOT quasi-static" verdict on that basis and was wrong — it
  flagged 277% on a cyclic run that is accurate to a fraction of a percent. Verdict removed from
  `cyclic.py`; `pushover.py` (Newmark, where the residual does track the error) keeps a threshold
  referenced to the measured 4.6% baseline.
- **EA reported (`build.report_calibration`, used by `elastic`/`pushover`/`cyclic`).** `A_t` alone
  is not what the struts are built with — a truss's stiffness is `EA/L`, and Aydin's balance returns
  EA and divides by E only to make the result transferable (D47). At mesh 50 / horizon 1.5,
  `A_t = 6042.9 mm^2` gives `EA` = **9.729e7 / 1.595e8 / 5.656e8 N** for test-unit / upper /
  pedestal. One area serves all three casts, but the test unit's struts are 1.6x softer than the
  head's and 5.8x softer than the pedestal's, whose E is an equivalent modulus (`E * PED_W/TW`) and
  so not a material property at all; diagonals are sqrt(2) softer than orthogonals from the same EA.
- **Also corrected while in the file:** `build.py` and `elastic.py` still described the lattice as
  "isotropic at nu = 0.18", the wording D53 flagged as wrong. `nu_consistent` is where the two
  CALIBRATION ROUTES agree, not a Poisson ratio; `elastic.py` now prints `nu_effective` (0.41) and
  `cubic_anisotropy` alongside it.
- **Status:** accepted. Backend: `quasi_static` on both dynamic runners (opt-in, existing results
  unchanged). Wall: new `specimen.QUASI_STATIC_RATE`/`DAMPING_RATIO`, `build.report_calibration`,
  reworked `--rate`/`--periods`/`--damping`/`--no-quasi-static` on pushover + cyclic. Suite: 42 pass,
  1 pre-existing known failure (`test_rc.py::test_nonlinear_pushover_runs_and_yields`, D34).
  NOTE the published cyclic results in CLAUDE.md predate this and were run at 0.8 / 7.6; re-running
  at the new defaults is a separate job and the peak should move by ~1%.

### D63 — 2026-08-19 — Wall example: vertical strain gauge over the base, and the strain profile across the section
- **Context (user):** measure vertical displacement at the wall–pedestal face and 250 mm above it,
  difference the aligned nodes, divide by 250 mm, and plot the resulting strain distribution across
  the wall width — to be used AFTER a long analysis run.
- **What it is.** Two aligned node rows, one on the W–P face (`y = 0`) and one at `y = GAUGE_H`,
  differenced column by column: `eps(x) = (uy(x, 250) - uy(x, 0)) / 250`, positive = tension. This
  is the same construction as the line of vertical LVDTs on the tested wall's face — the paper's
  Fig. 21 flexure/shear/rocking split is derived from them — so the model quantity and the measured
  one are the same thing, not merely analogous. 21 columns at 50 mm across the full 1000 mm width.
- **New backend probe: `node_history=(node_ids, dof)` + `node_history_every`,** on all four
  path-following runners (`run_pushover`, `run_pushover_dynamic`, `run_cyclic`,
  `run_cyclic_dynamic`). It follows a handful of nodes through the WHOLE run, which is the
  complement of the existing `capture` (the whole field at two instants). Returns
  `{"nodes", "dof", "index", "values"}` where `index` points into `disp`/`shear`, so every profile
  is tied to the drift it was actually taken at whatever the stride was. `run_cyclic`'s loop was
  refactored onto a `record()` closure in the process (it appended in four places); behaviour is
  unchanged.
- **Sampling: stride PLUS every reversal.** A 4% protocol at the D62 rate is ~1.8M steps, so
  `every=1` would hold tens of millions of floats; the default stride is 200. But a bare stride
  misses the loop tips by up to a full stride of drive — 0.32 mm, which is 15% of the FIRST
  protocol amplitude (2.2 mm) — so the cyclic runners also sample whenever the control
  displacement turns around. The tips are therefore exact at any stride, and they are the instants
  a profile is wanted at.
- **`gauge.py` (new, `examples/wall/`)** reduces a `node_history` to `(drift %, profiles)`, selects
  the samples nearest requested drift levels, builds the JSON payload, and draws. Shared by
  `pushover.py`, `cyclic.py` and `replot.py` so there is one implementation. `viz.figure_strain_profile`
  draws in MICROSTRAIN (a wall base runs to a few thousand; reading 2.4e-3 off an axis is worse
  than reading 2400).
- **The fit is included because it is what the profile is FOR.** `specimen.neutral_axis` puts a
  least-squares line through each profile: its slope is the curvature, its zero crossing is the
  neutral-axis position (which the paper reports), and its residual is a DIRECT test of "plane
  sections remain plane" — an assumption a lattice is under no obligation to satisfy. Drawn as a
  faint dotted line with the crossing ticked, and tabulated as curvature / NA position / NA depth.
- **Verified, and the verification is worth keeping.** The static `run_pushover` and the dynamic
  `run_cyclic_dynamic` — different solvers, different load paths — produce the SAME profile at the
  same drift: at +0.10%, curvature -1.052e-06 vs -1.052e-06 1/mm and neutral axis at
  x = -188.0 vs -188.1 mm. The physics also reads correctly: tension on the -x side for a +x push,
  the neutral axis migrating from x = -188 mm at 0.10% drift to +30 mm at 0.30% as the tension side
  cracks (compression zone 69% -> 47% of the section), and every cyclic profile rotating about a
  common point at x ~ 0, eps ~ -200e-6 — the uniform axial strain from the 600 kN axial load.
  The tension side departs visibly from the straight-line fit while the compression side stays
  linear, which is the discrete cracking localizing into individual strut columns.
- **Persistence, since the run costs hours.** The gauge series (x, gauge, drift, strains) goes into
  the run's `_data.json`, and `replot.py` redraws profiles at any levels (`--levels 0.3 1.0 2.0`)
  without touching the analysis. `pushover.py` now writes a `_data.json` too, which it did not
  before. The D62 `dynamic` residual is stored DECIMATED at the gauge stride (plus its exact
  peak/rms): three full-length arrays would have made a 4% run's JSON ~100 MB.
- **Status:** accepted. Backend: `node_history` on four runners + `_node_probe`;
  `viz.figure_strain_profile`. Wall: `specimen.GAUGE_H`/`gauge_nodes`/`gauge_strains`/
  `neutral_axis`, new `examples/wall/gauge.py`, wiring + `--gauge-every` on pushover and cyclic,
  `--levels` on replot. Suite: 42 pass, 1 pre-existing known failure (D34).
  Outputs: `wall_{pushover,cyclic_dynamic}_strain_profile.png` + the `gauge` block in each `_data.json`.

### D64 — 2026-08-19 — Drive rate returned to the published 7.6 mm/s; the missing measurement shows DAMPING, not rate, controls the contamination — partially superseding D62
- **Context (user):** "Make default loading rate to published, 7.6." D62 had lowered it to 2.0 at a
  cost of ~27 h against ~7 h for the full 4% protocol.
- **The combination D62 never measured.** Its sweep varied damping at rate 7.6 and rate at
  zeta = 0.05, but never ran **zeta = 0.2 at 7.6** — the pair the new defaults actually imply. It
  takes 39 s, and it inverts the conclusion:
  ```
    zeta   rate   peak kN   reversal ring   steady resid   share
    0.2    7.6     219.1        22.0 kN         8.8 kN      4.0%   <- best measured
    0.2    3.8     231.9        10.5           29.1        12.5%
    0.2    2.0     221.7         5.3           10.2         4.6%
  ```
- **Result: the published rate is not a compromise, it is the best point measured.** At zeta = 0.2
  the residual over rates 7.6 / 3.8 / 2.0 reads 4.0 / 12.5 / 4.6% — scatter between discrete
  cracking events, not a trend — so **slowing the drive buys nothing once the damping is right**.
  D62 drew its "2.0 more than halves the contamination" conclusion from single samples at
  zeta = 0.05, where lowering the rate genuinely did make things worse; it over-read that as
  evidence that rate mattered generally. It does not. The contamination is crack-release ringing,
  and DAMPING is what suppresses it: at 7.6 mm/s, 10.3% at zeta = 0.8, 12.4% at 0.05, **4.0% at
  0.2**. The whole D62 gain survives, and it was the damping change all along.
- **What the rate DOES control, cleanly and monotonically, is the start-up and reversal transient**
  — 22.0 / 10.5 / 5.3 kN at 7.6 / 3.8 / 2.0. Each of the 49 reversals in the 4% protocol carries
  that. So a study of LOOP SHAPE near the reversals (pinching, unloading stiffness) is a real reason
  to lower the rate; peak strength, envelope and damage location are not.
- **A caveat D62 should have stated.** These are single runs, and the peak shear itself scatters:
  219.1 / 231.9 / 221.7 kN over the three rates, ~±3% about the mean, driven by which struts crack
  in which step. Any conclusion drawn from a difference smaller than that needs repeats, which the
  D62 entry's 0.7%-agreement argument did not have.
- **Settled: `QUASI_STATIC_RATE = 7.6` (back to published), `DAMPING_RATIO = 0.2` (unchanged from
  D62).** The 4% protocol returns to ~464k steps / ~7 h, reproduces the committed runs' drive, and
  carries less solver contamination than those runs did. `--rate 2.0` remains one flag away for a
  reversal-shape study.
- **Status:** accepted; supersedes D62's rate choice, keeps its damping choice, its measurement
  (`quasi_static`) and its integrator caveat. `specimen.py` constants + rationale rewritten to the
  measured numbers; `pushover.py`'s advisory threshold re-referenced to the 4.0% baseline and now
  points at --damping before --rate. No backend change.

### D65 — 2026-08-24 — Katrin's wall verified against the source paper, and the measured WSH3 results captured as comparison data (`digitize.py` + `testdata.py`)
- **Context (user):** "Inspect the paper. See if katrins wall example matches with the model here.
  Also, digitize the hysteresis loop and all meaningful results so that when we run analysis in next
  steps, we have data to compare."
- **Verification — the model MATCHES the paper.** Every number in `specimen.py` was traced back to
  Dazio, Beyer & Bachmann (2009) and re-derived independently: geometry (Sec. 2.1 / Fig. 2), the
  Fig. 1c dimension string 30|100|100|125x5|145x2|125x5|100|100|30 = 2000 giving 17 positions x 2
  curtains = 6phi12 + 22phi8 + 6phi12, all four steel grades (Table 2), concrete (Table 3), axial
  load and every reinforcement ratio (Table 1, rho_tot 0.821 vs 0.82), V_max/M_max (Table 5) and the
  displacement set (Table 4). `summary.py` already reported all of this; nothing was found to be
  wrong. Three things were ADDED to that report rather than changed in the model:
  - **`f_t` is 17% high, and it is checkable.** The paper never prints a tensile strength, but
    Table 5's M_cr does: M_cr = (f_ctm + N/A_g)*t*l_w^2/6 with N/A_g = 2.287 needs **f_ctm = 2.97
    MPa**, which is EC2's 0.30*f_ck^(2/3) at the CHARACTERISTIC f_ck = f'_c - 8 = 31.2, not at the
    mean f'_c = 39.2. The same reading reproduces WSH1/5/6's M_cr to 1%, so it is the paper's
    convention, not a coincidence. `specimen.FT = 0.30*FC^(2/3)` = 3.46 MPa gives M_cr = 575 vs the
    paper's 527 kN.m (+9%). Left as-is and FLAGGED as a warned row, because it is a constitutive
    choice the user should make knowingly: everything the lattice does before yield is set by when
    its struts crack.
  - **Two omissions added to "Deliberately NOT modelled":** the SECOND concrete cast (Table 3
    characterises only phase 1, foundation + wall to 1.5 m; the concrete above is never reported,
    and the model applies the phase-1 grade throughout — defensible because the 1700 mm plastic zone
    is inside phase 1), and STRAIN PENETRATION into the foundation (bars anchored 450 mm into an
    ELASTIC block with perfect bond cannot produce the fixed-end rotation the test measured
    separately in Fig. 9a — which is why base CURVATURE is the fairer local comparison).
  - **`nu_consistent` mislabel fixed.** `summary.py` still printed it as "lattice Poisson ratio",
    the wording D53 corrected in `examples/wall/build.py` but not here. Now reports `nu_consistent`
    = 0.178 (the nu at which the normal-strain and shear routes agree, NOT a Poisson ratio),
    `nu_effective` = 0.412 and `cubic_anisotropy` = 1.38.
  - Known and already-documented approximations stand: boundary centroid snapped 870 -> 850 mm at
    mesh 50 (**2.3% on the flexural lever arm**, the largest single geometric error — mesh 25 gives
    875, i.e. 0.6%), ties 75 -> 100 mm with the area rescaled to hold the steel-per-height exactly,
    and the head replaced by a stiff full-width cap holding the shear span to 0.2%.
- **`digitize.py` — the hysteresis.** Fig. 7 is a 3x2 grid of panels (WSH1..WSH6) embedded as ONE
  1418x774 GRAYSCALE raster at 200 dpi, so unlike the SW-NC-FF figure (D51) colour cannot separate
  the annotation from the data. Three filters in sequence do it instead: (1) a 2x2 EROSION finds
  compact glyph seeds and cuts out their padded bounding boxes — this is what kills digits that
  TOUCH the curve; (2) the axis rules are removed BEFORE the component pass, not after, which is
  what detaches the ductility TICK STUBS that sit at exactly the loop-tip displacements, so they
  then fall to a small-and-compact component test that a severed loop arc (long, thin) survives;
  (3) the cloud is clipped to |V| <= 1.03*V_max from Table 5 — a physical bound, since no pixel can
  lie above the largest shear the test recorded. One hand-cut box remains (WSH3's "26" label), and
  it is documented with the clearance that makes it safe. Frames and ticks are found from the
  figure's own rules via LONGEST CONTIGUOUS RUN (total darkness fails: label text is just as dark).
  - **Result, WSH3: peak +453 / -449 kN against Table 5's 454, and +-93.2 mm against Table 4's 92.4**
    — inside one pixel (4.3 kN, 0.59 mm). All six units run; peak loads land within 2-3% of Table 5.
  - It is a point CLOUD, not an ordered path — overlapping loops cannot be re-sequenced from pixels,
    so **no per-cycle hysteretic energy**, same caveat as D51.
  - The reduced curve is a **backbone anchored on the PROTOCOL, not an envelope of the cloud**: the
    cycle amplitudes are known exactly (integer multiples of delta_y = 15.4 mm), so each loop tip is
    read from a +-3 mm window at a known displacement instead of from a running maximum any stray
    pixel can lift. WSH3: 371 / 423 / 449 / 453 / 453 / 453 kN at mu = 1..6 — a clean plateau at
    V_max.
- **`testdata.py` — everything the paper states as a number**, transcribed with its source and no
  digitization error: strength (Table 5), the three different yield displacements (Table 4 — the
  3/4-rule 15.4 mm is the one the LOADING HISTORY was built on, and so the one the protocol must
  use), the exact protocol (Fig. 6: two force-controlled cycles to 0.75F_y then two cycles at each
  mu = 2..6, South first, with the three actuator speeds), the load-step -> drift map (Fig. 11c
  legend) so a model state can be compared at the same instant the paper reports a profile or a
  photograph, the damage milestones (Sec. 3.1), Table 6's limit-state strains kept in their two
  separate flavours (plane-section vs Demec — the paper's own point in Sec. 5.1 is that they
  disagree by ~25%, and a lattice computes something closer to Demec), and shear/flexure = 0.12.
  Only Fig. 15(d) base curvature and 15(e) L_ph are digitized (marker centroids; no tabulated
  counterpart) and are marked as such.
- **Two findings that matter for the NEXT step, the cyclic run:**
  - **The test is ~130x slower than SW-NC-FF.** Fig. 6 gives 1.2-3.6 mm/MINUTE, i.e. 0.02-0.06 mm/s
    against that wall's 7.6 mm/s. A dynamic-relaxation run cannot be driven at the real rate, so
    D64's conclusion applies with more force here: DAMPING, not rate, is the knob.
  - **WSH3 never met the paper's own failure criterion** (Sec. 3.1: the drop was slightly under 20%
    when the test was stopped at LS63), so its 2.04% drift capacity is a LOWER bound.
- **Scope, and why this wall is the better target (already argued in `specimen.py`):** WSH3 uses
  DEFORMED bars that stay bonded, so perfect bond is a fair assumption and drift capacity becomes a
  fair comparison — unlike SW-NC-FF, where 74% of the drift at 2% was rocking on debonded plain bars
  (D52). What is still out of reach is the failure MODE: `Steel02` has no bar buckling and no
  fracture, and those are what ended the test at 1.79%.
- **Status:** accepted. Example-layer only — `specimen.py`, `build.py` and the backend are untouched;
  `summary.py` gained the f_t/M_cr rows, two omissions and the nu fix. New: `digitize.py`,
  `testdata.py`, `data/wsh{1..6}_fig7.npz`. Output -> `.../katrin_wall/<unit>_fig7_digitized.png`.

### D66 — 2026-08-24 — Katrin's wall given the Kutay-wall pipeline (elastic / pushover / cyclic / gauge / replot / compare), plus a curvature chain the SW-NC-FF study has no equivalent of; cyclic run SIZED but not launched
- **Context (user):** "Follow modelling and analysis approaches of kutay wall for katrin wall, and
  discuss what are the differences from modelling and analysis perspectives." Then, after the
  discussion: cyclic not pushover; staged to 1.02% drift; 7.6 mm/s; pre-flight first; add the
  height-wise curvature profile before the run; hand the long run over.
- **Ported file-for-file** from `examples/wall/`: `build.py` split into `wall_lattice` (elastic) and
  `nonlinear_wall_lattice` (`compression` + `gf_factor`) with `report_calibration` printing EA per
  zone; `elastic.py`, `pushover.py`, `cyclic.py`, `gauge.py`, `replot.py`, `compare.py`; and the
  `specimen.py` half — protocol, drive settings, gauge. `summary.py`'s builder call updated.
- **THE BOND ASYMMETRY RUNS THE OTHER WAY, and it is the headline modelling difference.** SW-NC-FF's
  PLAIN bars debonded — 74% of its displacement at 2% drift was rocking (D52) — which capped that
  study at strength, stiffness and damage location. WSH3 uses DEFORMED bars that stay bonded; the
  paper decomposes its displacement into flexure plus ~12% shear with a small fixed-end part
  (Fig. 9a/b). So **loop SHAPE, energy dissipation, residual drift and drift capacity all become
  fair comparisons here**, which is exactly why the cyclic run is the deliverable and the pushover
  is only a diagnostic. What is still out of reach is the FAILURE MODE: bar buckling from 1.70% and
  the corner-bar rupture at 1.79% that ended the test — `Steel02` has neither.
- **Other modelling differences, all consequences of the specimen:** n = Es/Ec is 5.7 here against
  ~12-14 there, so the bars add 5.5% to I rather than ~20%; axial load ratio 0.058 vs 0.20; measured
  steel hardening b = 0.005-0.009 (real hardening) vs 0.002-0.004 (plain-bar plateau); WSH3 is
  CONFINED where SW-NC-FF is not.
- **Confinement deliberately gets NO separate grade**, unlike the column study (D23). The hoops are
  already in the model as in-plane rebar struts between the boundary bars, and in a plane lattice
  those struts DO restrain the in-plane transverse expansion of the compressed boundary — the
  mechanism a confined law represents. A Mander grade on top would count the same restraint twice.
  Absent either way: the hoops' through-thickness legs and the phi4.2 crossties, both out of plane.
- **`epsc0` derived, not quoted** — 2*fc/Ec = 0.00223, per D56, since Concrete02's compressive
  tangent is 2*fc/epsc0 whatever the grade's E says. The quoted 0.0022 was only 1.2% off, but
  deriving it costs nothing.
- **Analysis differences.**
  - Calibration transfers UNCHANGED: same lattice constant A_t = 0.6058*(t*mesh). Elastic result
    `K_lattice / K_transformed` = **0.930**, the same as SW-NC-FF's, the gap being base-block
    flexibility the fixed-base hand calc omits. An unforced hit: the model puts **13%** of the
    flexibility in shear against the test's measured **12%** (Fig. 9b).
  - **The static solver is not merely worse here, it is unusable: it stalls at 0.067% drift**
    against SW-NC-FF's 0.3%, because WSH3 cracks at ~0.025% (M_cr is only 25% of M_max, Table 5).
    Dynamic relaxation is the only route, not a preference. `cyclic.py` therefore defaults to
    `--solver dynamic`, where the SW-NC-FF one defaults to static.
  - **The test's own rate is unreachable and this is the sharpest analysis difference.** Fig. 6
    gives 1.2-3.6 mm per MINUTE = 0.02-0.06 mm/s, against SW-NC-FF's published 7.6 mm/s. Any usable
    drive is ~130-380x the experiment, so the licence is the MEASURED residual, never the rate.
    7.6 mm/s is kept — the same ABSOLUTE speed, not the same drift rate: the contaminating forces
    scale with speed while this wall's shear is ~2x larger, so the same absolute contamination reads
    about half as much of the peak. Damping stays 0.2 (D64: damping is the knob, not rate).
  - **50 mm is the COARSEST LEGAL MESH** — gcd(4550, 4050, 600, 2800, 2000, 150) = 50 — so the usual
    lever for a too-expensive run is not available without retuning the model heights.
  - Gauge window is the PAPER'S, not the base: y = 50-350 mm, clear of the lowest LVDT the paper
    itself discards to strain penetration (Sec. 4.2) and bracketing its h = 60-360 Demec window.
- **NEW, and beyond SW-NC-FF parity: the LVDT chain and curvature over height.** That paper reports a
  flexure/shear/rocking split; this one reports CURVATURE (Fig. 11c) and base curvature vs
  displacement (Fig. 15d), which a single base gauge cannot match. `specimen.LVDT_ROWS` reconstructs
  the test's own chain from Fig. 5b — 50, 3x100, 4x300, 2x700 — and the reconstruction is checkable:
  it totals 2950 mm and Fig. 11's height axis runs to exactly 3.0 m. The lowest device (0-50 mm) is
  deliberately absent, as in the paper. `curvature_profile` gives one curvature per segment at its
  mid-height; `base_curvature` reproduces the paper's OWN construction — a linear fit over the
  plastic zone extrapolated to y = 0, with L_pz the height at which curvature falls to phi_y = 3/km
  — so the model's Fig. 15d quantity is the same quantity, not an analogue. `gauge.figure_curvature`
  draws it in the paper's layout with the fit and its extrapolation marked.
- **`preflight.py` (new): size the run by MEASUREMENT, not estimate.** A short dynamic pushover past
  first cracking, reporting T1, ms/step and the contamination, then projecting each protocol scope.
  **Measured at mesh 50: T1 = 19.5 ms, 276 ms/step, steady inertia+damping 16.4 kN** (9.4% of the
  173 kN shear at 0.040% drift — high because the shear is still small; at the 454 kN peak the same
  absolute residual is ~3.6%, comparable to SW-NC-FF's 4.0%). Projections, LOWER bounds since cost
  per step rises as struts crack: 0.68% -> 68k steps/5.2 h; **1.02% -> 143k steps/11.0 h**;
  1.35% -> 243k/18.7 h; full 2.03% -> 518k/39.7 h. Against SW-NC-FF's 4% protocol at 464k steps/7.1 h
  — same step count, ~5x the cost per step, because this model is 4513 nodes / 20,678 elements.
- **Post-processing validated against a SYNTHETIC response before committing to the long run** —
  gauge reduction, both figures, the JSON payload and `replot.py` all exercised on a plane-section
  field with a known neutral axis and base curvature, both of which come back exactly. An 11-hour
  run must not die in its plotting code. JSON size projects to ~1 MB at `--gauge-every 200`.
- **Status:** accepted. Example-layer only; backend untouched. NO cyclic run launched — the user runs
  it: `python examples/katrin_wall/cyclic.py --drift 0.0102` (defaults: dynamic, 7.6 mm/s, zeta 0.2,
  crushing, mesh 50, horizon 1.5). Staged deliberately: 1.02% covers cracking, first yield (0.25%),
  spalling onset (1.02%) and 449 of the 454 kN peak; whether the plateau and post-peak earn another
  ~29 h is a decision to make on the first run's evidence.

### D67 — 2026-08-24 — First WSH3 cyclic attempt DIVERGED before the first reversal; cause traced to strut brittleness, fixed by correcting f_t to the paper's own value (+ tension stiffening) — partially superseding D66
- **Context (user):** ran the D66 command; "Check the analysis output".
- **THE RUN FAILED AND PRODUCED NO HYSTERESIS.** 1,865 steps of 143,443 (1.3%) in 814 s,
  `converged=False`, stopping at **+0.202% drift** — short of the first reversal at +0.2505%. The
  tell is in the summary line itself: `peak base shear +355.3 / -0.0 kN`. Every figure was a
  monotonic ramp drawn as a loop, and `compare.py` dutifully reported "0.78x push, **0.00x pull**".
- **A red herring first: the reported "inertia + damping = 764% of peak shear" was NOT the
  failure.** It is ~2676 kN CONSTANT from step 200 while the shear grows 93 -> 310 kN — a
  rate-proportional offset, i.e. the documented HHT(0.7) numerical-dissipation term (D64: ~60x what
  Newmark reports; the Newmark pre-flight on the same model read 16.4 kN, and 16.4 x 60 is the same
  order). A residual that is flat while the response grows is not a diverging solve. `cyclic.py` now
  reports the residual's TREND (last fifth / first fifth) with a verdict, instead of a bare
  percentage that reads like a catastrophe.
- **The actual failure:** Newton starts failing at step ~600 = **0.065% drift**, cost per step
  climbs 0.20 -> 0.55 s as it falls back through KrylovNewton and NewtonLineSearch, and it gives up
  at 1,865. 0.065% is exactly where the STATIC solver stalled (0.067%, D66) — dynamic relaxation
  postponed the same mechanism by 3x rather than curing it.
- **CAUSE: this wall's struts are far more brittle than SW-NC-FF's.** With crack-band
  regularization `Ets = f_t^2*L/(2*Gf)`, a strut's tensile life is `eps_ult/eps_cr = 1 + E/Ets`:
  | | f_t | Gf | Ets/E | eps_ult/eps_cr |
  | SW-NC-FF test-unit (ran to completion) | 1.50 | 0.053 | 0.066 | 16.2 |
  | WSH3 as first run | 3.46 | 0.077 | 0.110 | **10.1** |
  60% more brittle, on top of cracking at 0.025% drift instead of 0.047%, a 2000 mm section at an
  axial load ratio of only 0.058 (wide, barely pre-compressed tension zone), and horizon 1.5.
- **`f_t` CORRECTED PERMANENTLY, and it was the flagged one.** D65 raised that `specimen.FT =
  0.30*f'_c^(2/3)` = 3.46 MPa is 17% above the 2.97 the paper's own M_cr implies, and left it
  flagged rather than changed. That was the wrong call: `Ets` goes as `f_t^2`, so a 17% error in f_t
  is 36% in softening steepness. `FT` is now `0.30*(f'_c - 8)^(2/3)` = **2.973 MPa**, and the
  model's implied cracking moment is **526 kN.m against Table 5's 527 — 0.2%**, where it was +9%.
  `summary.py`'s warned row becomes a passing check. NOTE the two MC90 quantities deliberately take
  DIFFERENT strengths: `f_ctm` is characteristic-based, `G_f = G_f0*(f_cm/10)^0.7` is defined on the
  mean, which stays `f'_c`.
- **SCREEN (`preflight.py --drift 0.0025 --gf-factor 2`) PASSES**, at f_t = 2.97 with tension
  stiffening x2 (Ets/E = 0.040, 25.9x — gentler than the SW-NC-FF wall that completed):
  reaches the full 0.250% drift with `converged=True`, and **28 Newton failures in 2,307 steps
  (1.2%) against the failed run's 68 in 1,865 (3.6%)**. Physics checks too: 344 kN at 0.250% drift
  against the digitized test's ~315 kN there — 1.09x, the stiffer direction `cyclic.py` predicts
  from perfect bond plus the missing strain penetration.
- **COST IS HIGHER THAN D66 PROJECTED, and D66's number was measured in the wrong place.** That
  pre-flight drove only to 0.06% drift — just BELOW where the trouble starts — so its 276 ms/step
  was an uncracked-lattice cost and it saw a healthy run. Measured properly to 0.25%: **413 ms/step
  average, 439 ms/step marginal** over 0.04-0.25%, still rising. Revised lower bounds: 0.68% ->
  7.8 h, **1.02% -> 16.4 h** (was 11.0), full 2.03% -> 59.4 h. Lesson: a screening run must reach
  PAST the instability it is screening for, or it certifies nothing.
- **New guard in `cyclic.py`:** a run whose displacement never changes sign now prints
  "*** THE RUN NEVER REVERSED — there is no hysteresis in this result ***" above the results, naming
  the drift it stopped at and the reversal it missed. The failed run silently produced a full set of
  loop figures; that must not happen twice.
- **Residual to watch, not a blocker:** the Newmark pre-flight reads 36.2 kN = 10.5% of the shear so
  far, tracking proportionally (16.4 kN at 173 kN earlier), so ~8% at the eventual peak against the
  SW-NC-FF wall's 4.0%. Do NOT address it with `--rate`: the residual is LINEAR in rate, so halving
  the runtime at 15.2 mm/s would double it to ~20%.
- **Status:** accepted; supersedes D66's `f_t` decision (flagged -> corrected) and its cost
  projection (276 -> 413 ms/step), keeps everything else. No run launched from here — the user runs
  them (`cyclic.py --drift 0.0102 --gf-factor 2`). Example-layer only; backend untouched.

### D68 — 2026-08-25 — `examples/vk3_wall/`: ETH bridge pier VK3, the first SHEAR-relevant specimen; loops digitized from VECTOR paths; `node_history` given multiple dofs; and a node-merge bug found in the shared mesher
- **Context (user):** "By following the modelling and analysis techniques kutays wall and katrins wall,
  I would love to model the wall in the chapter I shared, VK3." Discussed first; the user chose VK3
  alone (not the VK1/VK3 pair), asked for the backend extension enabling the shear measurement, and
  chose to leave the base construction joint unmodelled and flagged. Plus: add a `--draw` argument
  and save the digitized hysteresis for comparison.
- **Source:** Bimschas, M. (2010), "Displacement Based Seismic Assessment of Existing Bridges in
  Regions of Moderate Seismicity", IBK Bericht Nr. 326, Institut fuer Baustatik und Konstruktion,
  ETH Zuerich; vdf Hochschulverlag. doi:10.3929/ethz-a-006237119. Chapter 5.
- **The specimen and why it earns a third wall package.** Test Unit VK3 of a three-unit ETH Zurich
  campaign on squat wall-type bridge piers (1500x350, L_v = 3300, L_v/l_w = **2.20**, 42phi14
  continuous deformed bars, phi6@200 hoops, N = 1300 kN). Its aspect ratio is essentially WSH3's
  (2.28) but it carries **three times less transverse steel** (rho_sw 0.08% vs 0.25%), shear reaches
  **20-22%** of the top displacement against WSH3's ~12%, and the test ended in a combined SHEAR +
  AXIAL-LOAD failure at 1.59% drift. It is the first specimen in the repo where shear is an actor.
  The chapter's own reading of the failure is lattice-shaped: not diagonal tension, not web crushing,
  but the loss of the base compression zone that had been SUPPORTING the diagonal strut — a load-path
  collapse, which `element_groups` already measures (cf. D55's 28% -> 0%).
- **Reading of the drawings verified by the chapter's own printed ratios**, before any modelling:
  2x17phi14@80 + 2x4phi14 at the ends reproduces rho_sl = 1.232% (printed 1.23), 2phi6@200 gives
  rho_sw = 0.081% (printed 0.08), and N_base/(A_g f_c) = 0.0768 (printed 0.077).
- **Mesh 50 mm, and no geometric rounding at all.** Every model dimension (1500 / 3300 / 3700 /
  3000 / 900) is a whole number of cells, so the shear span is EXACT and drift needs no correction —
  unlike WSH3, which had to round its height. 50 mm is also the coarsest legal grid: it keeps all 19
  in-plane bar positions distinct, where 100 mm collapses four of them. The base hoop zone's 200 mm
  spacing is exact (only its phase shifts, 75 -> 100 mm); only the head zone's 75 mm needs the
  area-rescaling snap, and the head is elastic anyway.
- **THREE MATERIAL QUANTITIES ARE CONVENTIONS, NOT MEASUREMENTS, and this is the biggest data gap
  against WSH3.** The chapter gives f'_c = 34 MPa and nothing else about the concrete — no E_c, no
  f_t, and **no cracking moment to pin f_t against** (WSH3's M_cr is exactly what settled that
  argument in D67). Chosen and flagged in `summary.py`: E_c = SIA 262's 10000*f_cm^(1/3) = 32.4 GPa
  (the specimen is Swiss and the chapter's frame is SIA 262; EC2/MC90 agrees within 2%, ACI is 15%
  low and rejected as the wrong code frame); f_t = EC2 at f_ck = f_cm - 8 = 2.633 MPa, the D67
  convention; G_f = MC90 on the mean. Strut life eps_ult/eps_cr = **14.2**, between SW-NC-FF's 19.8
  (ran clean) and WSH3's 13.4 (needed `--gf-factor 2`), so budget tension stiffening.
- **k0 = 60 kN/mm IS NOT AN ELASTIC TARGET** and chasing it would be a category error: it is a secant
  to FIRST YIELD on an already cracked pier, 0.24x the uncracked transformed cantilever
  (246.6 kN/mm). The digitization supplies a real one instead — see below.
- **THE LOOPS ARE VECTOR, NOT RASTER, and that changes what the data can be used for.** Figs.
  5.10-5.13 are PDF path objects with explicit stroke colours, so `digitize.py` reads geometry rather
  than pixels via `pdftocairo -svg` (no new dependency; poppler was already required). Three
  consequences, each an improvement on the WSH3 digitizer: no digitization error; no text/curve
  separation problem (the curve is the only GREEN geometry, so the three-stage erosion filter WSH3
  needs is unnecessary); and **the loops stay ORDERED**, so energy dissipation, residual drift and
  per-cycle degradation all become recoverable, where the WSH3 cloud carries a warning against
  integrating it. The curve arrives as 11 disjoint polylines stored in reverse; `_chain` re-sequences
  them by endpoint matching, worst junction **0.007 mm**.
- **CALIBRATE ON TICKS, NOT ON THE FRAME.** These panels' frames span +-110 mm x +-950 kN while the
  printed tick labels stop at +-100 and +-800; assuming frame = axis limits would inflate every load
  by 19%. Major ticks are twice the length of minor ones, which is what separates them. The result
  then CHECKS ITSELF against a wholly independent axis: each panel carries a secondary base-moment
  axis, and 200 kN at L_v = 3.3 m must be 0.66 MN.m — recovered ratio **1.0000**.
- **WHAT THE DIGITIZATION BOUGHT, beyond a curve to plot against.** VK3's peak base shear is **not
  tabulated anywhere in the chapter**, so without this there is barely a backbone: it comes out at
  **+891 / -876 kN at ~0.95% drift**, against the chapter's predicted F_n = 851 kN. Two further
  recoveries neither planned nor expected:
  * **the measured ELASTIC amplitudes** (0.85 / 2.10 / 5.08 mm). The elastic block was FORCE
    controlled, so its displacements were never published; converting through k0 gives 2.68 / 5.37 /
    8.05, wrong by up to 3x for exactly the reason above. The first level (+-0.86 mm at ~151 kN) sits
    below the pier's own V_cr = 209 kN, so its secant — **~174 kN/mm** — is a genuine UNCRACKED
    stiffness measurement, and the elastic target the chapter never prints.
  * **the small intermediate cycles**, which the chapter defines for VK3 only as "the top
    displacements measured during the corresponding cycles of VK1" and never tabulates. All 13 are
    in the figure, including the deliberately ASYMMETRIC second one of each pair. `cyclic.py`
    therefore defaults to `--protocol measured`, driving the pier's own 60 turning points, with the
    designed protocol kept as `--protocol nominal`.
- **BACKEND (D68 proper): `node_history` now accepts several dofs.** `dof` may be one 1-based
  component or a sequence; samples are stored dof-MAJOR, so a one-dof probe records exactly what it
  did before and every existing reader (the D63 vertical gauge in `wall/` and `katrin_wall/`) is
  untouched. Needed because a diagonal chord's length change is not a function of either
  displacement component alone. All four path-following runners carry it.
- **The shear instrument, which is the point of the study.** The pier carried four string
  potentiometers measuring diagonal elongations "in order to determine the shear deformations"
  (Sec. 5.2.2) — the measurement behind Fig. 5.19-right, which reports the shear / flexure / base-
  crack split as percentages of top displacement. `specimen.panel_shear` reduces each panel's two
  diagonals as `gamma = sqrt(b^2+h^2)*(dd1-dd2)/(2bh)`, and `deformation_components` integrates
  flexure from the same probe's vertical strains. VALIDATED ON SYNTHETIC FIELDS FIRST (D66's lesson):
  pure shear recovers gamma to 1.000000, pure flexure gives shear 1e-13 mm and flexure to 1.000000,
  and a combined field separates both exactly. Neither other wall package has an equivalent — the
  SW-NC-FF split was 74% rocking the model could not fairly be compared against, and the WSH3 paper
  reports curvature instead.
- **BUG FOUND IN THE SHARED MESHER, pre-existing and affecting `katrin_wall`.**
  `mesh_compound_rectangles` merged coincident nodes on `round(coords, 9)`. gmsh returns a shared
  edge's node as -750.0 from one rectangle and -749.999999996 from the other — 4e-9 apart, a
  RELATIVE difference of 1e-12, i.e. ordinary double-precision noise — which rounds to two different
  9-decimal keys, so the nodes survived as a coincident pair and `connect_horizon` then joined them
  with a strut of length ~4e-9 mm. Its EA/L is some ten orders of magnitude above every real strut's:
  a near-rigid link that does kinematically what the merge should have done while wrecking the
  conditioning of the assembly. Counts: **katrin_wall 109, vk3_wall 50, wall (SW-NC-FF) 0** — the
  SW-NC-FF pedestal happened to land on identical floats. Fixed by keying the merge to the GRID
  (`coords / (mesh_size * 1e-6)`, i.e. 5e-5 mm at mesh 50 — absorbing the noise while separating
  real nodes by ~1e6), plus a backstop in `connect_horizon` that refuses pairs closer than
  `1e-4 * mesh_size`. After the fix all three models report zero, shortest strut = 50.000 mm;
  katrin_wall drops 4513 -> 4472 nodes, vk3_wall 3484 -> 3453. Test suite: 42 pass, 1 fail — the
  known D34 failure, unchanged.
  **CAUTION for D66/D67:** WSH3's published numbers were produced WITH these 109 links present, and
  D67 attributed its divergence at 0.065% drift entirely to strut brittleness. That diagnosis may be
  incomplete — ill-conditioning from near-rigid links is a plausible contributing cause, and the
  D67 run should be repeated on the fixed mesher before its cost projections are trusted.
- **Scope, stated as sharply as the other two studies state theirs.** FAIR: peak strength, stiffness,
  the shear / flexure split, loop shape, energy dissipation, residual drift, damage location. NOT
  FAIR: anything at or past failure. The hoops were closed by simple 90-degree hooks "which could
  open after spalling of the cover concrete" (Sec. 5.2.1b) and a perfectly-bonded rebar strut never
  loses anchorage, so the model's transverse steel works forever — **expect over-retention of
  strength beyond ~1.27% drift.** Also absent: bar buckling and cover spalling (from 0.95% drift,
  and the drivers of the compression-zone loss that triggered the failure), and aggregate interlock,
  without which the four-wedge sliding kinematic cannot be expressed at all.
- **Also decided:** no stiff cap — VK3's top 400 mm is ordinary pier, so the model runs the full
  3700 mm height with the horizontal actuator at y = 3300 and the vertical actuators at y = 3700,
  the first time in this repo the two loads act at different levels. The load-introduction zone
  (y >= 2700, phi6@75, detailed to stay elastic) is held elastic; it sits above the ~2000 mm the
  diagonal cracking reached. The foundation's out-of-plane thickness is never given, so 1000 mm is
  ASSUMED and folded into its modulus, priced by `elastic.py --foundation-sensitivity`.
- **MEASURED AFTERWARDS (2026-08-25):** the WSH3 1.02% run completed on the PRE-FIX mesher —
  143,443 steps, 21.7 h, converged, peak base shear **+431.4 / -434.0 kN against the test's 454
  (0.950)** and base curvature at 1.02% drift **9.00 /km against the measured 9.1 (0.99)**. Archived
  under `examples/output/katrin_wall/runs/2026-08-25_cyclic_1p02_gf2/`; `runs/README.md` states the
  rule (archive first, analyse second) because the analysis scripts write to a fixed stem and would
  otherwise overwrite a 21.7 h result with a two-minute diagnostic.
  An ELASTIC A/B of the two meshers on that model gives **K_fixed / K_prefix = 0.9861**, i.e. the 109
  near-rigid links added **1.4% stiffness**; nodes 4513 -> 4472. User's call: accept 1.4% and move on.
  **THE NONLINEAR EFFECT IS NOT MEASURED, and the elastic number is a weak proxy for it.** Those
  struts are not rigid once cracked — crack-band regularization sets `Ets = f_t^2*L/(2*Gf)`, so at
  L ~ 4e-9 mm `Ets ~ 0` and the tension branch is a FLAT PLASTIC PLATEAU that never fails, capping
  at `f_t*A_t` ~ 13.5 kN each. In the long run they were therefore ~109 spurious elastic-plastic ties
  across the WALL-FOUNDATION INTERFACE — where a flexural wall's tension chord anchors — with ~1.5 MN
  of combined capacity, able to yield and let coincident nodes separate. Whether that biased the
  21.7 h result is unquantified; the test that would settle it is a short nonlinear pushover past
  cracking on both meshers at `--gf-factor 2`, ~20 min, not run.
- **Status:** accepted. Example package complete (specimen / testdata / build / draw / digitize /
  summary / elastic / pushover / cyclic / gauge / replot / compare / preflight); `--draw` on every
  analysis script. Backend touched twice: `opensees._node_probe` (multi-dof) and `mesh.py` (the
  merge fix). NO FE RUN LAUNCHED — the user runs them. `summary.py` reproduces every checkable
  chapter value to <1%.

### D69 — 2026-08-25 — VK3 will not run at horizon 1.5: the dynamic solve breaks up at 0.31% drift. Horizon is NOT the fix, and the reason is a coupling between horizon and crack-band brittleness
- **Context (user):** after preflight sized the cyclic run, ran a dynamic pushover as the cheap check
  before committing overnight. It went unstable. User then asked why VK3 differs from the two walls
  that completed, and asked for a short screen before adopting any fix.
- **THE FAILURE.** `pushover.py --solver dynamic --gf-factor 2` rises smoothly to **699.5 kN at
  0.282% drift** at ~175 ms/step with 0.12% convergence failures, then between 0.282% and 0.310%
  breaks into violent oscillation — base shear swinging **+1355 / -1329 kN while drift increases
  monotonically**, 48.3% convergence failures, and a 27x slowdown to ~4,800 ms/step as the
  sub-stepping fallbacks fire every step. A pushover cannot go negative on a monotonic ramp, so this
  is numerical, not softening. Run preserved as
  `runs/2026-08-25_164857_pushover_dynamic_gf2_UNSTABLE/`.
- **MY PREFLIGHT CERTIFIED NOTHING, and it repeated D67's own mistake verbatim.** It screened to
  **0.25%**; the instability starts at **~0.31%**. Its clean result (0.12% failures, 112 ms/step,
  converged) described only the pre-instability regime, and the 13-17 h cost projection built on it
  is **withdrawn** — in this configuration the run does not complete at all. D67 wrote "a screening
  run must reach PAST the instability it is screening for, or it certifies nothing"; that sentence
  is quoted in `preflight.py`'s own docstring, above a default drift chosen below the trouble.
  `preflight.py --drift` now defaults unchanged but every screen since is run to 0.6%.
- **WHY VK3 AND NOT THE OTHER TWO** — three compounding differences, all measurable without running
  anything:
  | | thickness | A_t | f_t | **crack force per strut** | rho_transverse | L_v/l_w |
  | SW-NC-FF | 200 | 6,043 | 1.50 | **9.1 kN** | 0.39% | 3.0 |
  | WSH3 | 150 | 4,544 | 2.97 | **13.5 kN** | 0.25% | 2.27 |
  | **VK3** | **350** | **10,595** | 2.63 | **27.9 kN** | **0.08%** | 2.20 |
  A_t scales with THICKNESS (`0.605*t*mesh`), and VK3 is a 350 mm pier against 150-200 mm walls, so
  every cracking strut dumps 2.1x WSH3's and 3.1x SW-NC-FF's force into a mass-damped system. Its
  transverse steel is 3-5x lower and sits on a coarser pitch (phi6@200 = every 4th node row against
  WSH3's phi6@150 = every 3rd), so there is far less tie steel to restrain a cracked node. And it is
  the squattest, with the chapter reporting shear cracks growing WIDER than flexural ones by
  mu_prov = 4 — i.e. the cracking lands on the inclined struts, not the vertical ones. The ordering
  of transverse steel (0.39 / 0.25 / 0.08%) tracks the outcomes exactly (4% drift / 1.02% / fails at
  0.31%).
- **HORIZON WAS THE WRONG FIX, and the screen proved it.** `cyclic.py --horizon`'s help recommends
  3.01 as "the redundant bracing that keeps cracked tension-side nodes from going singular" — advice
  inherited from the RC column study (D31), where shear accuracy did not matter. TWO reasons not to
  take it here:
  1. **The shear-stiffness error is a SAWTOOTH in horizon**, because each threshold admits a
     different strut ring and an axis-aligned ring without its diagonal partners wrecks the
     directional balance: 1.5 -> **2.67%**, 2.05 -> 108%, 2.3 -> 25.9%, 2.9 -> **7.66%**, 3.01 ->
     **25.2%**, 3.2 -> 38.0%, 3.7 -> **5.32%**. Horizon 1.5 is the BEST of the whole range on the
     one quantity VK3 exists to measure, and the recommended 3.01 sits in a trough.
  2. **Horizon fights itself under crack-band regularization.** `Ets = f_t^2*L/(2*Gf)`, so strut
     life `1 + E/Ets` FALLS as struts lengthen — and a bigger horizon is exactly what creates long
     struts. Longest-strut life: h=1.5 -> **18.6**, h=2.9 -> **10.1**, h=3.7 -> **8.1**, against
     D67's markers of 19.8 (SW-NC-FF, clean) and 10.1 (WSH3, DIVERGED). So h=2.9 drags its longest
     struts to exactly the value at which WSH3 diverged and h=3.7 below it.
     Measured static reach confirms the non-monotonicity: **0.151% (h=1.5) -> 0.200% (h=2.9) ->
     0.100% (h=3.7)**. Restoring h=1.5's brittleness at h=2.9 would need `--gf-factor ~3.9`, twice
     what is already applied and beyond what `build.py` calls modest.
- **TWO DISTINCT FAILURES WERE BEING CONFLATED.** The STATIC stall at 0.151% is the ordinary
  softening-branch limit of displacement control — a lower bound, not a capacity, the same
  distinction D59 drew on `compression_cube`. The DYNAMIC break-up at 0.31% is separate: the
  dynamic solve passes straight through 0.151% cleanly. If a kinematic mechanism existed at 0.15%
  the dynamic run would ring from there; it does not. So bracing addressed neither.
- **REVISED DIAGNOSIS: crack-release RINGING, not a mechanism.** D64 established that the
  contamination in these runs is crack-release ringing and that DAMPING, not rate and not topology,
  is the knob. Two facts fit: VK3's per-strut crack release is 2.1x WSH3's, so its ringing should be
  about twice as violent; and preflight measured the residual at **143% of base shear at 0.25%
  drift**, before anything went wrong. `damping_ratio` is 0.2, tuned by D64 on SW-NC-FF — a wall
  with a third of VK3's crack-release energy. D62 originally used 0.8. Raising damping costs nothing
  in model size, adds no brittleness, and preserves horizon 1.5's 2.67% shear error.
- **RESOLVED — IT WAS THE INTEGRATOR, AND THE SCREEN WAS RUN ON THE WRONG RUNNER.**
  `run_cyclic_dynamic` (the deliverable, and `preflight.py`) uses `ops.integrator("HHT", 0.7)`.
  `run_pushover_dynamic` used `ops.integrator("Newmark", 0.5, 0.25)`. **Every one of the seven
  failed diagnostics was run on the pushover runner.** Newmark average-acceleration has gamma = 0.5,
  which produces EXACTLY ZERO numerical damping and conserves energy, so it sustains a spurious
  high-frequency mode indefinitely; HHT sets gamma = 1.5 - alpha = 0.80 and beta = (2-alpha)^2/4 =
  0.4225, dissipating high-frequency response while leaving the structural response untouched.
  Changing NOTHING but the integrator, on the same model, mesh, materials, damping and dt:
  | | Newmark | HHT(0.7) |
  | at 0.310% drift | 0.1 kN (collapse) | **727.1 kN** |
  | peak | 967.9 kN — ABOVE the 848 kN section capacity, unphysical | **730.6 kN at 0.314%** |
  | convergence failures | 6,568 | **30** |
  | runtime (0.6% screen) | 6.3 h at dt/4 | **8 min** |
  | past peak | +-600 kN chaos | smooth 46.7% degradation |
  Peak 730.6 kN = **0.820x the digitized test** and 0.858x the chapter's F_n — and 0.858 sits on the
  elastic `K_lattice/K_transformed` = 0.855, so the lattice appears uniformly ~14% softer/weaker
  than beam theory in both regimes.
- **WHY RAYLEIGH DAMPING COULD NOT SAVE IT, which is the transferable part.** Both runners use
  `ops.rayleigh(zeta*w1, 0, zeta/w1, 0)` — mass + betaKINIT, never betaK (D49). That gives the
  Rayleigh U-curve `zeta_eff(w) = zeta/2*(w1/w + w/w1)`, so for VK3 (EA/L = 6.86e6 N/mm against a
  2.1e-3 t nodal mass, w_strut/w1 ~ **165**) the stiffest strut modes already carry **1652% of
  critical damping** at zeta = 0.2 — and it still blew up. Rayleigh damps the RESPONSE while an
  energy-conserving integrator RE-INJECTS energy at every sudden stiffness change. Algorithmic
  damping removes it at source. This also explains the two partial successes: zeta 0.2 -> 0.5 and
  dt/4 each bought ~0.03% drift because both attack the same band, far less efficiently.
- **THE PROCESS FAILURE, recorded because it is the lesson.** Seven diagnostics — horizon x2,
  damping, compression, R0, b, dt — every one a MODEL or DRIVE parameter, every one testing a
  hypothesis about the structure, every one eliminated. The difference in integrators was noticed
  and written down BEFORE five of them were run ("every failure above was on the pushover runner...
  worth one screen") and then not acted on. **Screen with the integrator your deliverable actually
  uses**; a screen on a different solver certifies nothing. Corollary: the damage map (9.6% of
  struts cracked at 0.290% drift — a near-intact model) refuted all five constitutive/topological
  hypotheses in ONE short run and should have been the first diagnostic, not the seventh.
- **`run_pushover_dynamic` now takes `integrator=`** (default Newmark for back-compat) and
  `progress=`/`progress_every=`, mirroring `run_cyclic_dynamic`; `pushover.py` exposes
  `--integrator {newmark,hht}`, `--steps-per-period`, `--capture`, and the two DIAGNOSTIC-ONLY
  overrides `--steel-r0` / `--steel-b` (the latter overrides a MEASURED value and is labelled as
  never a production setting). Tests: 42 pass, 1 known D34 failure, unchanged throughout.
- **PREFLIGHT RE-RUN PAST THE BARRIER (D67's rule applied properly), on HHT:** `--drift 0.006
  --gf-factor 2` reaches 0.600% drift, `converged=True`, **29,812 steps in 3,284.7 s = 110 ms/step**
  with **123 failures (0.41%)**, residual 130.8% of the 731 kN reached. The cost does NOT blow up
  past cracking — 110 ms/step against the pre-barrier 112 — so **D68's withdrawn projection is
  reinstated and now MEASURED**: mu_prov = 2 -> 4.1 h, mu_prov = 3 (peak) -> 6.2 h, mu_prov = 4 ->
  9.1 h, **failure at 1.59% -> 12.7 h** (lower bounds).
- **Also fixed:** `preflight.py`'s milestone table cut at the NOMINAL drift level, but these are the
  pier's MEASURED turning points and the actuator overshoots (mu_prov = 2 peaks at 21.26 mm, not
  21.00), so every milestone was truncated — the mu_prov = 2 path read 398 mm against the true 624.
  Now cut at `level * 1.02`.
- **Status:** resolved; the barrier was an artefact of the screening tool, not a property of the
  model. VK3 cyclic run sized at 12.7 h and NOT yet launched.


### D70 — 2026-08-27 — `examples/compression_cube/` given VK3's timestamped run directories; `console.log` teed from `run_dir` itself rather than by a shell wrapper
- **Context (user):** "Have a look at compression cube example. Make run outputs as structured as
  vk3 wall." The cube wrote every run to a FIXED stem under `examples/output/compression_cube/`
  (`compression_cube[_h3][_rigid][_static]_*.png` + `_data.json`, `aydin_approach_rmax*.png`), the
  pattern D68 flagged as the trap `katrin_wall` nearly lost a 21.7 h result to.
- **Why it matters HERE, and it is not the D68 reason.** The cube's runs take ~52 s, so nothing
  irreplaceable is at risk. The problem is different and worse for this specimen: on the cube almost
  every question is a comparison ACROSS runs — free vs `--rigid-platen` (D53/D59), `--solver static`
  vs `dynamic` (D57/D59), `--algorithm` (D58), horizon 1.5 vs 3.01, Concrete02 vs `aydin_approach.py`
  (D61), one `--rmax` against the next. The stem encoded only SOME of those knobs, so an algorithm
  sweep or an `--rmax` step overwrote its own predecessor and a sweep became a memory rather than a
  listing.
- **Implementation (example-layer only; no backend, no results change):** `specimen.run_dir(kind,
  tag=)` and `specimen.find_run(kind, name=, complete=)` ported from `vk3_wall/specimen.py`.
  `pushover.py` → `<stamp>_pushover_<solver>_<free|rigid>[_h..][_m..]/`, `aydin_approach.py` →
  `<stamp>_aydin_[calib_]rmax<r>[_norsm][_Truss]/`, each holding `command.txt`, `console.log`,
  `data.json` and generically-named figures (`stress_strain.png`, `force_deformation.png`,
  `loadpath.png`, `damage.png`, `model.png`) — so the directory name carries the configuration and
  the file names carry the content, instead of both being crammed into one stem.
- **`--calibrate` is the one case a tag cannot name.** The converged `Rmax/d` is an OUTPUT of the
  search, so the directory records the value the search STARTED from (`calib_rmax0.0500`) and the
  converged value lives in `data.json`. Naming the directory after the answer would require creating
  it after the run, which is exactly what loses the `command.txt` of a run that dies.
- **`aydin_approach.py`'s Concrete02 baseline is now looked up, not hard-coded.** It read
  `OUT / "compression_cube_data.json"` — a path that only ever existed because of the fixed stem.
  It now takes `find_run("pushover")`, the newest run WITH a `data.json`, and PRINTS which run it
  overlaid; `data.json` doubles as the completion marker, so a cancelled or in-flight run is skipped
  rather than silently overlaid as an empty curve.
- **`console.log` is teed from `run_dir`, which is a departure from VK3 and a better one.** VK3 gets
  its log from a shell wrapper (log to /tmp, move in on exit) because the directory is a timestamp
  nobody can redirect into ahead of time. Here `run_dir` swaps `sys.stdout`/`sys.stderr` for a `_Tee`
  writing to both, line-buffered, guarded against double-installation. VERIFIED: openseespy routes
  its own output through `sys.stdout`, so the OpenSees convergence WARNINGs — the thing a stalled
  static run (D59) has to be read from — land in `console.log` too, as do tracebacks from a transient
  that blows up. No wrapper, no redirection, and `tail -f` still follows a run in flight.
- **Rationale for teeing at all:** on this specimen the printed diagnostics ARE the result. The
  modulus check against `expected_modulus`, the four-line capacity accounting (D55) and the load-path
  split are console-only — none of them is in a figure — so a run whose numbers live in a terminal
  scrollback cannot be quoted next to its own figures.
- **Status:** accepted. Scripts compile and `--help` clean; `run_dir`/`find_run`/the tee exercised
  against a redirected `RUNS` (fresh directory, `command.txt`, `latest_*` symlink, incomplete-run
  rejection, stdout+stderr+openseespy capture). No FE run launched — the first real run creates the
  first directory. `examples/output/compression_cube/runs/README.md` documents the layout. The flat
  output directory was empty, so nothing needed migrating. `aydin_cube/` still writes a flat stem.

### D71 — 2026-08-27 — Compression cube: a DERIVED `fc` scaling that lands the lattice peak on the material fc, the compression twin of D41 — and the peak STRAIN it necessarily moves
- **Context (user):** given the ranked options for the 0.634*fc shortfall, "let's scale fc to match".
  This is option 2 of that list: correct the strength, leave the stiffness calibration alone.
- **The knob.** `--peak-correction strength` on `pushover.py`; `build.cube_lattice(fc_scale=...)` and
  `build.scaled_grade`. It mirrors `examples/cube/build.py`'s `peak_correction="strength"` (D41),
  which fixes the TENSION cube's 1.50x overshoot by scaling `ft` and leaving the area alone — the
  same separation of knobs (area -> K0, strength -> peak), with the sign reversed.
- **The factor is DERIVED, not fitted.** `build.strength_correction_factor` reads the strut list:
  at horizon 1.5 / mesh 20, 11 vertical struts cross the mid-height cut and carry 0.6341 of the
  gross face, so the factor is 1/0.6341 = **1.5771**, computed before the analysis runs and with no
  reference to a measured curve. It is legitimate because D54 measured `peak/(verticals at fc)` =
  1.000 there: the inclined path is fully shed by peak, so the deficit IS the verticals' share.
- **`epsc0` must scale with `fc`, and that is forced rather than chosen.** Concrete02's initial
  compressive tangent is `2*fc/epsc0` whatever the grade's `E` says (D56), so holding `epsc0` while
  raising `fc` would multiply every strut's stiffness by 1.577 and destroy the energy-balance
  calibration this correction exists to preserve. Deriving `epsc0 = 2*fc'/E` holds the tangent at
  exactly `E` — **verified to 1.000000** at both horizons.
- **THE PRICE, and it is not avoidable.** Each strut now peaks at `1.577*epsc0`, so the lattice's
  peak STRESS is corrected onto fc and its peak STRAIN moves by the same factor (1.88e-3 -> 2.96e-3).
  Concrete02's parabola ties (E, peak stress, peak strain) together — three quantities, one free
  parameter — so no single scaling fixes the stress and holds the strain. The run prints this as an
  explicit `PRICE:` block rather than leaving it to be discovered in the curve.
- **The mechanism is NOT corrected, and gets slightly worse in relative terms.** `ft` is deliberately
  untouched (it is measured, and raising it is the *other* knob — per-orientation tension, option 3
  of the same list), so splitting still begins at the same absolute strain `eps_cr/nu` = 2.38e-4.
  Against the enlarged strut `epsc0` that is **8.0% instead of 12.7%**, i.e. relatively EARLIER. The
  correction buys the right peak, not the right failure process; the printed load-path split and the
  damage map are unchanged in character.
- **Self-check built in.** `capacity_accounting` now takes `fc_strut` separately from the continuum
  target, which always uses the MATERIAL fc. Under correction `verticals/continuum` reads **1.0000**
  by construction, so the accounting table verifies its own premise on every run. `fcu` scales too
  (the residual keeps its ratio to fc); `Gf`/`Gfc` do not, so the post-peak branch changes shape and
  is not claimed to be corrected.
- **A bias found while wiring it.** `secant_factor` removes Concrete02's parabolic curvature from the
  least-squares modulus fit using `epsc0` — the STRUT's, which the correction moves. Reading the
  material value there biases the modulus verdict by ~1.5% at the default fit window (0.9699 vs
  0.9809), and that verdict is quoted to four decimals. It now takes `epsc0=`.
- **Horizon 3.01 is the documented limit.** The factor comes out 1.4996 there, but D54/D55 measured
  the diagonals still holding 22% at peak with the peak arriving at 0.48*epsc0, BEFORE the verticals
  reach fc — so `verticals` is the wrong denominator and the correction over-corrects. `pushover.py`
  already prints the measured diagonal share at peak, so the assumption is checked per run rather
  than assumed; the docstring says so.
- **Verification, no FE run launched (D-standing user preference).** Algebraic: factor, tangent
  preservation, `verticals/continuum` = 1.0000, and the actual Concrete02 arg tuples for orthogonal
  and diagonal struts at both horizons (ft and Ets confirmed UNCHANGED at 1.50 / 424.5 / 600.4).
  Path: `main()` driven end-to-end with a SYNTHETIC response through both branches, exercising every
  print, the JSON and all four figures — which is how the title overflow was caught (the correction
  note now occupies its own two short lines; matplotlib does not wrap titles).
- **Status:** implemented, unrun. `data.json` gains `peak_correction`, `fc_scale`, `fc_material`,
  `fc_strut`, `epsc0_material`, `epsc0_strut`, `peak_over_fc`; runs tagged `_fcscaled` (D70). The
  prediction to test is `peak_over_fc` = 1.00 at horizon 1.5; anything else means the shedding
  assumption behind the factor does not hold as measured.

### D72 — 2026-08-28 — The 2019 JSE paper read in full: its PUBLISHED energy balance added as a second calibration route, its BOND elements added as an option, and `examples/aydin_aldemir_wall/` opened on the squat wall it simulates worst

**Source.** Aydin, B. B., K. Tuncay, and B. Binici (2019), "Simulation of Reinforced Concrete Member
Response Using Lattice Model", *J. Struct. Eng.* 145(9): 04019091,
doi:10.1061/(ASCE)ST.1943-541X.0002381. This is the **missing middle** of the Aydin trilogy the repo
already leans on: the 2017 MSc thesis supplies D47's energy balance, the 2021 paper supplies D60's
tension-only compression cube, and this 2019 paper is the only one of the three that models
*reinforced* concrete — and the only one where the closed-form `Et*At` is actually published.

**It independently confirms D53, which no other source had.** The paper states the horizon-1.5d
lattice's Poisson ratio "changes between 0.26 and 0.42 depending on the rotation of the loading axis"
and is "about 0.33 irrespective of direction" at 3.01d. Rotating the loading axis through the repo's
own `nu_effective` diagnostic reproduces this: **0.266 (45 deg) to 0.408 (0/90 deg)** at horizon 1.5,
and **0.314-0.353, mean ~0.33**, at 3.01. So D53's correction — that `nu_consistent` ~ 0.18 is NOT
the lattice's Poisson ratio — now has outside corroboration, and the paper's own range is a direct
measurement of the cubic anisotropy the repo reports as 1.38.

**1. A SECOND CALIBRATION ROUTE (`field=`), because the thesis and the paper do not use the same
affine field.** `calibration.energy_balance_area`/`energy_balance_rectangle` take
`field="uniaxial"` (unchanged default, D47: `eps_x` with the transverse strain restrained, continuum
side `E e^2 / (2(1-nu^2))`) or `field="equibiaxial"` (the 2019 Appendix Eqs 3-6: equal stresses give
`eps_x = eps_y = sigma(1-nu)/E`, continuum side `E e^2 / (1-nu)`). Under the equibiaxial field every
strut at a node sees the same axial strain whatever its direction, which is what collapses the
balance to their closed form `Et At = 4 Et A w / ((1-nu) sum L_i) = C Et d w`.
`calibration.aydin_closed_form_C` reproduces the published **C = 0.621 (horizon 1.5d) and C = 0.102
(3.01d)** to 0.05%.
- **They disagree, and by a structural amount.** At matched nu and horizon 1.5 the uniaxial route
  returns **1.056x** the EA at nu = 1/3 and **1.173x** at nu = 0.20. Since every truss stiffness
  scales linearly with EA, that is the factor by which each wall study's global stiffness would move.
- **The repo's existing numbers survive by a near-cancellation, not by agreement.** SW-NC-FF's
  published `EA = 9.729e7 N` (uniaxial at nu = 0.20) is **0.973x** the paper's closed form evaluated
  at ITS nu = 1/3 (9.998e7) — within 3%. The field difference and the nu difference nearly cancel.
  That is worth knowing precisely because it is a coincidence: `EnergyBalanceResult` now carries
  `field` and `area_equibiaxial`, so both numbers are on every result and a study has to say which
  one it used. **No existing study's calibration changes**; the default is untouched.

**2. BOND ELEMENTS, OPTIONAL (`build_lattice_rc(bond_material=...)`).** Perfect bond has been the
only option since D5/D13, and CLAUDE.md names it as SW-NC-FF's dominant error (1.32x over-strength at
0.3% drift). The paper's scheme is implementable against the existing `Rebar` machinery and needs no
new backend element:
- each `Rebar` gets its OWN steel nodes, duplicated at the coordinates of the lattice nodes on its
  path, so steel and concrete no longer share a DOF;
- **bond struts** join every steel node to the concrete nodes within `bond_horizon * mesh_size`,
  excluding the coincident one (a zero-length truss is not an element) — the ring of 8 at horizon 1.5;
- `materials.bond_elastic_brittle` is his Fig. 1(c) law: elastic to the concrete's own `ft` at
  `eps_cr`, a near-vertical drop to `a*ft`, then a flat plateau, symmetric in slip direction. His
  `a = 0.7` (from Aydin 2017, for DEFORMED bars) is the default.
- **The residual is structural, not just constitutive.** A steel node is held ONLY by its bond ring;
  fully brittle links would leave it free and the tangent singular once the ring cracks. `residual=0`
  is rejected with that message.
- **Bars that cross get separate steel nodes** and are coupled only through the concrete — which is
  what the physical bars do, sitting in different layers through the thickness.
- **Mass had to be handled explicitly.** Steel nodes have no tributary area, and a zero-mass DOF makes
  the dynamic runners singular. `bond_mass_share` (default 0.1) moves that fraction of a concrete
  node's lumped mass to each steel node sitting on it, conserving the TOTAL exactly; a share that
  would strip a host node raises rather than silently producing a massless DOF.
- **Cost.** On the Aldemir wall at mesh 50 the model goes 2,806 -> 5,454 nodes and 13,501 -> 34,595
  elements, i.e. bond struts outnumber concrete struts about 2:1. Budget for it before a cyclic run.
- Drawing: `viz` gains a `bond` kind, drawn dashed and UNDER the skeleton — bond links lie exactly
  along concrete struts by construction, so on top they hide the entire lattice.

**3. `examples/aydin_aldemir_wall/` — the target specimen.** Aldemir, A., B. Binici and E. Canbay
(2017), "Cyclic testing of reinforced concrete double walls", *ACI Struct. J.* 114(2): 395-406, as
reported by the 2019 paper (the 2017 paper is NOT in the repo, so `testdata.py` labels every
"measured" value as second-hand and lists what is simply absent). The squattest wall in the repo by
far — 3000 long on a 2250 shear span, **ratio 0.75** against SW-NC-FF's 3.00, WSH3's 2.28 and VK3's
2.20 — with NO axial load, "designed to yield in shear", reaching ~1% drift with severe inclined
cracking and **no strength degradation at all**. It is also the specimen the 2019 paper predicts
worst: **+20.8% on peak force at horizon 1.5d and +37.6% at 3.01d**, against 4-11% for the other five.
That over-prediction is the thing worth attacking, and both features above are candidate causes.
- **50 mm divides every dimension** (3000, 2250, bar spacing 100, cover 50), which the paper's own
  d = 20 does not (2250/20 = 112.5).
- **Reinforcement Ø8 @ 100 both ways, 3 bars per position (rho = 1.257%), cover 50.** Fig. 10(a) is a
  raster with no dimension text, so it was measured: the outer rectangle gives a consistent
  1.708 mm/px and the drawn mesh comes to 30 vertical and 22-23 horizontal bar lines at 94.7 mm. The
  "3-Ø8" call-out is unusual for a 120 mm wall, so the RATIO is the evidence: at 3 bars the
  distributed-steel flexural capacity is 0.92x the measured peak, at 2 bars 0.65x.

**GEOMETRY — SETTLED at 3000 wide x 2250 high x 120 thick (user, 2026-08-28).** The paper prints two
incompatible panels (Fig. 10(a) 3000 x 2250; Fig. 4(f) 2680 x 1500) and its body text dimensions
nothing. Measuring the figures settles which one can be trusted:
- **Fig. 10(a) is drawn to scale and reconstructs itself.** Its two printed dimensions give 1.3529
  and 1.3505 mm/px — agreeing to **0.17%** — and on that one scale the drawn reinforcement rebuilds
  both: **30 vertical bar lines at 100.0 mm** with 52.1/46.7 mm cover give the 3000, and **22
  horizontal lines at 99.8 mm** with 72.9/80.4 mm cover give the 2250. Four independent quantities,
  one scale.
- **Fig. 4(f) is not to scale.** Its drawn block has an aspect of **2.235**, matching none of its own
  candidate labels (2680/1500 = 1.787 is nearest, 25% off), so it cannot arbitrate geometry.
  `PANEL_FIG4F` is kept only to reproduce the paper's Table 4 on its own stated geometry.
- **A COVER CORRECTION the measurement forced.** Cover is NOT 50 mm both ways as first assumed:
  the vertical bars sit at 50 but the horizontal bars at **75**. `COVER_X`/`COVER_Y` are now separate
  measured constants. 75 is not a multiple of the 50 mm grid, so `bar_lines` snaps the cover DOWN —
  preserving all 22 lines (snapping up drops one and moves rho by 4.5%) at the price of the
  horizontal cage sitting 25 mm low; `--mesh 25` removes the offset. `specimen.check_bar_layout`
  reports modelled-vs-drawn counts, spacing and cover on demand and `summary.py` prints it.

**WHAT REMAINS OPEN IS NOT THE PANEL.** The measured initial stiffness of 1,038 kN/mm EXCEEDS the
804 kN/mm an uncracked transformed section of the settled panel can give, which is impossible for a
fixed-base cantilever. With the height fixed at 2250 the only candidate that lands on both checks is
a **TOP RESTRAINED AGAINST ROTATION** (1.11x stiffness, force unchanged at 0.92x) — a squat wall
loaded through a stiff head. Doubling the thickness, which was the best fit while the height was
still open, now OVERSHOOTS both (1.49x and 1.31x), so the "double wall" of the 2017 title is not the
explanation. The remaining possibility is that Table 4's "initial stiffness" is not a fixed-base
cantilever secant at all. Judge the calibration on `K_lattice/K_transformed`, not on
`K/K_measured`.

**A separate finding about the PAPER'S MODEL (not the specimen).** Inverting Table 2's own integers —
20,385 particles and 80,684 elements at d = 20 — over a horizon-1.5 grid gives a UNIQUE solution (up
to transposition): **134 x 150 cells = 2680 x 3000 mm**, matching neither printed geometry. Two
things follow: '2680' in Fig. 4(f) is a real dimension while '1500' is not the height, and **Table 2's
element count is concrete struts only** — no steel, no bond — which matters when comparing model
sizes against `build.describe`. `summary._invert_paper_grid` reproduces it.

**BOND IS AN ASSUMPTION HERE, NOT AN IDEALIZATION.** WSH3 and VK3 have documented continuous deformed
bars; SW-NC-FF has documented plain bars, which capped its scope. This specimen's bar type is simply
not reported, and the 2017 title says "double walls" — a precast twin-shell system whose interface is
never described. That is exactly why it is the right specimen to run BOTH ways.

**Verification, no FE run launched (standing user preference).** `aydin_closed_form_C` against the
published 0.621/0.102; the two routes' ratio; the nu-vs-loading-axis sweep against the paper's
0.26-0.42 and ~0.33; bond ring counts, lengths (exactly `d` and `d*sqrt2`), mass conservation and law
shape as `tests/test_bond.py` (8 tests) and four new cases in `tests/test_energy_balance.py`; the
36 non-OpenSees tests pass. `draw.py` rendered for all three panel/calibration/bond variants.
**STAGE 1 — STATIC PUSHOVER (run 2026-08-28). It stalls, as predicted, and `--gf-factor` does not
help.** `pushover.py`, mesh 50, horizon 1.5, no axial load, target 1% drift:

| run | stalled at | peak shear | model/test | steps |
|---|---|---|---|---|
| Gf as given | **0.0263% drift** | 356.4 kN | 0.370 | 11 |
| `--gf-factor 2` | **0.0175% drift** | 271.4 kN | 0.282 | — |

- **This wall cracks at 0.0082% drift** (V_cr = 148 kN from its own f_t with NO axial load to offset
  it) — an order of magnitude earlier than WSH3's 0.025%, which is the run whose static solver died
  at 0.067% and forced D66 onto dynamic relaxation. SW-NC-FF, cracking later still, reached 0.3%.
  The static solve here gets 3x past cracking and no further.
- **TENSION STIFFENING IS NOT THE REMEDY HERE, unlike WSH3 (D67).** Doubling Gf made the stall
  EARLIER (0.0175% vs 0.0263%), not later. That is consistent with strut life not being the binding
  constraint in the first place: 22.8 orthogonal / 16.4 diagonal, already above WSH3's 13.4. The
  barrier is the mechanism — a squat, axially unloaded wall forming a diagonal crack field — not
  strut brittleness, so stage 2 must be the dynamic solver rather than a re-tuned static one.
- **Whether 356 kN is a capacity or just where the solver stopped is undecided** and is exactly what
  stage 2 settles: if the quasi-static run passes it while still ascending, it was a stall.

**LOAD-PATH PROBE (`--groups`), and a bug worth recording.** `build.strut_groups` decomposes a
free-body cut by element orientation. First version summed `eleForce` over EVERY element and
reported a "diagonal share" of -16,000 kN against a 356 kN base shear — internal forces cancel in
equal and opposite pairs, so only elements CROSSING the cut may be counted. Second version probed
`dof` 1/3 believing them 1-based; they index `eleForce` directly (0-based), so it was reading the
VERTICAL components, which on an axially unloaded wall sum to ~0 and look like a dead probe. Both
now reconcile to **+1.0000** against base shear and against V*h. The shear split is geometrically
trivial (a vertical strut cut horizontally has no horizontal component, so diagonals take 100% by
construction) and is kept only as that self-check; the informative split is the OVERTURNING MOMENT,
which at 0.0263% drift is **40.8% vertical struts / 37.9% rebar / 21.2% diagonals** — i.e. ~79%
plane-sections couple against ~21% truss action at this early stage.

**Status:** library features implemented and tested. Wall package: `specimen`/`testdata`/`build`/
`draw`/`summary`/`elastic`/`pushover` (static); quasi-static pushover, cyclic, gauge, replot, compare
and preflight are not yet ported.

**Geometry decision (user, 2026-08-28): Fig. 10(a) as drawn** — 3000 x 2250 x 120, rho = 1.257%, on
the ground that it is the only reading taken from a dimensioned drawing of the specimen rather than
inferred. The consequence is accepted and documented rather than hidden: it is the weakest-scoring
candidate, and elastic results will read ~0.72x the measured stiffness.

**`elastic.py` (first analysis script) and what it is built to measure.** SHEAR carries **58.2% of
the uncracked flexibility** at this aspect ratio, against ~14% for VK3 and a few percent for
SW-NC-FF. Every earlier study could treat the calibration's shear weakness as a footnote; here the
energy balance's weakest property is the majority of the answer, so the script prints the shear share
and `isotropy_error` alongside the stiffness. Closed forms for the adopted panel: K_gross = 775.7,
**K_transformed = 804.3 kN/mm** (I_tr/I_g = 1.088), against the measured **1,038.44**. The
MEASUREMENT EXCEEDS THE GEOMETRY'S OWN UNCRACKED CEILING (0.775x), so the lattice cannot reach it
either — the script says so explicitly and directs the reader to `K_lattice/K_transformed`, the
like-for-like ratio the other three walls land at 0.93. `--geometry-sweep` re-solves on every
candidate panel to turn the open question into a measurement; panels the mesh cannot carry are
SKIPPED rather than rounded (rounding the height is precisely what the drift argument turns on), so
the sweep is single-row at mesh 50 and complete at `--mesh 10`. `--calibration` and `--bond` are
wired through, giving both new library features their first exercise on a real structure.
**Validated on a SYNTHETIC response** (no FE solve, standing user preference): every print, the
figure and `data.json` driven end to end, which is how the dead `--draw` branch and the sweep's
divisibility failure were both caught. `specimen.mesh_fits` is the non-raising form of the alignment
check that the sweep needs.

**ELASTIC RESULTS (run 2026-08-28).** Four solves at mesh 50, horizon 1.5, on the settled panel:

| run | K (kN/mm) | K/K_transformed | K/K_measured | secant spread |
|---|---|---|---|---|
| uniaxial (default) | 738.9 | **0.9187** | 0.7116 | 0.01% |
| equibiaxial | 635.2 | 0.7897 | 0.6117 | 0.00% |
| uniaxial + BOND | 1758.1 | 2.1858 | 1.6930 | **12.44%** |

- **The default calibration survives its hardest test.** 0.9187 against the 0.93 the other three
  walls land at — on the specimen where **shear carries 58.2% of the elastic flexibility** (VK3 ~14%,
  SW-NC-FF a few percent). The energy balance matches C11 exactly and never matches shear stiffness
  independently, so this was the run most likely to expose that; it did not.
- **The PUBLISHED route is measurably worse: 0.7897.** The 2019 paper's own equibiaxial balance makes
  the fit 13% worse on the 2019 paper's own specimen. `field="uniaxial"` stays the default, and that
  is now a measurement rather than an argument.
- **ELASTIC BEHAVIOUR IS ACCURATE, and the earlier 0.919 verdict was measured against the wrong
  reference.** Adding a plane-stress CONTINUUM on the same node grid, same rebar, same BCs (D12/D14's
  original purpose) gives **K_lattice/K_continuum = 1.0018 at mesh 50 and 0.9910 at mesh 25** — the
  lattice reproduces the continuum to a few tenths of a percent. Meanwhile
  **K_continuum/K_transformed = 0.9171**, i.e. essentially ALL of the 0.919 gap against the closed
  form is BEAM THEORY's error, not the lattice's: at aspect ratio 0.75 a cantilever formula assumes
  plane sections and a uniform shear distribution and this wall has neither. The other three walls'
  0.93 was measured against a reference that was fair for slender walls and is not fair here.
  `build.wall_continuum` and `elastic.py`'s continuum block exist so this is never conflated again.
- **The calibration is mesh-objective, as the balance claims.** `A_t` goes 3,614 -> 1,819 mm^2 from
  mesh 50 to 25 (ratio 0.5033 against the predicted 0.5 — the excess is the finite patch's edge
  nodes), and under that 2x refinement K moves **-0.70% for the lattice and +0.37% for the
  continuum**. Both are converged; neither drifts alone. `elastic.py --mesh-convergence`. Note the
  legal mesh set is small: it must divide gcd(3000, 2250, 100, 50) = 50.
- **THE BOND RESULT IS WRONG, AND THE DEFECT IS IN THIS IMPLEMENTATION, NOT THE PAPER.** 2.19x the
  transformed section, with a 12.4% secant spread in a fully LINEAR model. Bond can only ever ADD
  flexibility relative to perfect bond. Diagnosis: each bond link is given the concrete strut's area
  and modulus, so the 8-link star around every steel node forms a stiff parallel load path BETWEEN
  the surrounding concrete nodes — on top of the concrete lattice, which is still fully present. It
  fails the limit test that defines the feature: as bond stiffness goes to infinity the model must
  reduce to perfect bond, and instead the star over-constrains. **The horizon-star topology cannot
  satisfy that limit at all** while the concrete lattice remains intact; the topology that can is the
  coincident-node interface spring (`zeroLength`, bond law along the bar axis + a stiff transverse
  material) that this repo's CLAUDE.md anticipated from the start. `tests/test_bond.py` passes
  because it tests assembly — counts, lengths, mass, law shape — and never the stiffness the assembly
  produces; a limit test belongs there. **BOND IS DISABLED IN THIS STUDY BY USER INSTRUCTION
  (2026-08-28) and stays disabled until explicitly re-enabled**: `--bond` is gone from the scripts and
  `build.wall_lattice(bond=True)` raises.

### D73 — 2026-08-29 — Aldemir's wall is 210 mm thick, not 120: a REPLICA of Aydin's own model measured the error, and eight experiments show the premature collapse is a separate, still-open problem

**How the error was found — by building his model, not ours.** After six parameter experiments all
failed to move a premature collapse (below), `examples/aydin_aldemir_wall/replica/` was written to
answer a different question: does OUR model collapse because of how WE built it, or because of what
the lattice method does to this specimen? The replica is his model as closely as this codebase can
express it — his panel, his mesh, his calibration route, his material laws — with every departure
listed in `specimen.NOT_REPLICATED` and printed on every run.

It reproduced **his Table 2 counts exactly** (20,385 nodes, 80,684 concrete struts, from the unique
150x134 cell grid those two integers imply) and then missed his published `K_sim` = 943.16 kN/mm by a
factor of **1.743**. Lattice stiffness is EXACTLY linear in thickness — the energy balance gives
`A_t = C*d*w`, so every strut area is linear in `w` — which makes that factor directly readable as a
thickness: `120 * 1.743 = 209.1 mm`.

**FIVE INDEPENDENT CHECKS, none of them a fit to the collapse:**
1. the replica's stiffness miss, above;
2. re-run at 210 the replica returns **K = 923.7 kN/mm = 0.979 of his K_sim**, a number it was never
   fitted to;
3. **210 = 50 + 100 + 60** of the Fig. 10(a) plan's own `50|140|100|60` stack — two precast shells
   plus a cast-in-place core, with 140 the clear gap. The paper's title calls the specimen a
   "double wall"; 120 dimensions ONE panel of it, not the load-carrying section;
4. at 210 the UNCRACKED transformed section sits **1.31x ABOVE** the measured 1,038 kN/mm, as it must
   for a fixed-base cantilever. At 120 it sat at **0.75x, i.e. BELOW** — physically impossible, and
   the anomaly `summary.py` had been flagging for two days and D72 wrongly attributed to the top
   boundary condition;
5. at 210 the distributed-steel flexural capacity lands on the measured force, **V_flex/F_exp = 1.00**
   (0.92 at 120), because the bar COUNT is drawn while the RATIO it implies falls from an overstated
   1.257% to the true **0.718%**.

**WHAT IT COST.** Every run before this date modelled a wall **43% too thin**: struts carried 1.75x
too little area while the reinforcement stayed full size, so the lattice's concrete was badly
outmatched by its own steel. On the thickness change alone, at mesh 50 with everything else identical:

| | peak | at drift | collapsed | vs measured 963.6 kN |
|---|---|---|---|---|
| t = 120 | 591.4 kN | 0.109% | 0.110% (30 steps) | 0.61 |
| **t = 210** | **903.4 kN** | **0.195%** | 0.199% (84 steps) | **0.94** |

and for the first time a genuine DESCENDING BRANCH (880 -> 868 -> 645 kN) rather than a vertical drop.

**THE RESIDUAL QUESTION IS NOW ANSWERED BY MEASUREMENT.** D72 reported 74% of peak and could only
infer it was the collapse transient, because only the maximum was saved. `pushover.py` now persists
the residual SERIES, and on the t = 210 run: ascending branch **p95 = 1.6% of peak** (max 10.7%),
collapse window 49.4%. The peak is only lightly contaminated; the headline number is the transient,
as inferred but not previously demonstrated.

**THE COLLAPSE IS A SEPARATE, STILL-OPEN PROBLEM — eight experiments.** damping (0.2/0.5), fracture
energy (x1/x2), crack-release ringing, lattice topology, mesh (50/25), constitutive law (Concrete02
bilinear vs Aydin trilinear), reinforcement presence, thickness (120/210). Thickness moved it
furthest and did not remove it: the wall still loses its load path at **0.20% drift against the test's
~1%**. What remains that Aydin has and this study does not is BOND ELEMENTS and EXPLICIT INTEGRATION.
The replica can now test the second directly, since it matches his elastic stiffness to 2%.

**Notable negatives worth not repeating:**
- **Gf is not a neutral knob** — x2 raised base shear 14.2% at 0.097% drift, so it changes strength,
  not just stability. That caveat also lands on D67's WSH3 result, which used it.
- **Aydin's trilinear backbone made everything worse** (peak 451 kN, collapse at 0.097%, and 1,510
  solver failures against 19): its long tail is offset by an immediate drop to `b1*ft` at
  `1.5*eps_cr`, and `ElasticMultiLinear` is path-independent, so struts flip branches between
  iterations and Newton never settles. He runs it EXPLICITLY, with no tangent at all.
- **A plain-concrete wall (no rebar) collapses at 0.025% drift with 0.7% of struts cracked**, in a
  base-corner wedge, with ZERO solver failures — the load-path loss demonstrated in isolation. The
  reinforced wall's damage at peak is a DIFFERENT mechanism (28.2% cracked, diagonal field), so the
  two are not the same event; an earlier claim that they were was too strong.
- **`figure_damage`'s `eps_crush` takes a POSITIVE magnitude** and negates internally. Passing
  `-EPSC0` inverts the test and flags ~99% of struts as crushed; the giveaway is a `< --2.3e-03`
  legend. Fixed in `pushover.py`.

**PROCESS NOTE, recorded because it cost real time.** Three conclusions were reported prematurely,
each from a sample taken BEFORE the event being watched for — mesh 25 "fixed the collapse" (it
collapsed 250 steps later) and "the collapse is gone" at t = 210 (it came at 0.209%). The collapse
was already measured as occupying ~30 steps while progress was sampled every 250. **Report what a
run has shown, not what it implies about a drift it has not reached.**

**Status:** `specimen.TW` = 210 (override `ALDEMIR_TW` for A/B). `summary.py` rewritten to adjudicate
the thickness. All pre-2026-08-29 runs in `examples/output/aydin_aldemir_wall/runs/` were computed at
120 and are superseded on strength and stiffness; their NEGATIVE findings about damping, Gf, mesh and
constitutive law still stand, since all of them were internally consistent comparisons.

### D74 — 2026-08-29 — EXPLICIT time integration added to the dynamic runner: 45x cheaper per step, 2.7x faster end to end, and it cross-validates the Aldemir collapse as REAL rather than a Newton failure

- **Context (user):** after the implicit replica preflight failed to converge at 0.0062% drift with
  every fallback algorithm exhausted, and D73 closed with "explicit integration is not in this
  codebase", the user asked whether an explicit solver could be used. It can: this `openseespymac`
  build carries `CentralDifference`, `CentralDifferenceAlternative`, `CentralDifferenceNoDamping`,
  `NewmarkExplicit`, `HHTExplicit`, `AlphaOS`, `AlphaOSGeneralized` and `KRAlphaExplicit`.
  (`Explicitdifference` is NOT in this build — it reports "unknown integrator type" and silently
  falls back to the default, which is worth knowing because the fallback still returns 0.)

- **THREE CHANGES, and the second one decides whether the idea works at all.**
  `opensees.EXPLICIT_INTEGRATORS` / `is_explicit()` gate all of them inside `run_pushover_dynamic`,
  so a caller only passes `integrator=("CentralDifference",)`:
  1. **`system("Diagonal")` + `algorithm("Linear")`.** A lumped-mass diagonal system is the whole
     point — no factorization — and there is nothing to iterate, since the step IS the update.
  2. **THE `betaKinit` TERM MUST GO.** D49 established stiffness-proportional damping on the INITIAL
     stiffness as the right choice for these lattices, and it is — implicitly. Under an explicit
     integrator it adds `xi = betaK*w_max/2` at the highest mode and multiplies the stable step by
     `(sqrt(1+xi^2) - xi)`. For this wall `betaK = 4.6e-4` and `w_max = 2e5 rad/s` give `xi = 46`, so
     the stable step would shrink **~90x** — turning a 12-minute run into an 18-hour one. Explicit
     runs therefore use MASS-PROPORTIONAL damping only, which costs nothing in stability. Silent and
     fatal if missed.
  3. **No sub-step rescue.** A failed explicit step is INSTABILITY, not non-convergence; retrying
     with KrylovNewton or a tenth of the step cannot help, and doing so would hide the signal. The
     runner breaks with `converged=False` and leaves the fix (a smaller `dt`) to the caller.

- **SIZING IT: `builders.critical_time_step(model, modulus_of)`.** Explicit is conditionally stable
  at `dt < 2/w_max`, and `w_max` is set by the stiffest, lightest ELEMENT — not by T1. Sizing
  `steps_per_period` off T1, the implicit habit, is 15-20x too large here and diverges. Estimated
  element-wise as `w_e = 2*sqrt(EA/L / m_min)`. For the Aldemir wall at mesh 50: **dt_crit = 10.0 us
  against the 146 us the implicit run used**, so `--steps-per-period 750` rather than 40, and 95,310
  steps instead of 5,084.

- **THE ECONOMICS ARE THE OPPOSITE OF THE INTUITION.** 18x the steps and still faster, because a step
  costs a force recovery instead of a factorized Newton solve:

  | | steps | ms/step | wall clock |
  |---|---|---|---|
  | implicit Newmark | 5,084 | 347 | 29.4 min |
  | explicit CentralDifference | 95,310 | **7.8** | **10.8 min** |

  45x cheaper per step, **2.7x faster end to end**.

- **THE CROSS-VALIDATION IS THE REAL RESULT.** Two schemes sharing almost nothing — factorized
  tangent vs none, mass+stiffness damping vs mass-only, 5k vs 95k steps — on the same model:

  | | peak | at drift | collapse | residual, ascending p95 |
  |---|---|---|---|---|
  | implicit | 903.4 kN | 0.1950% | 0.1991% | 1.56% |
  | explicit | 904.6 kN | 0.1847% | 0.2030% | 0.97% |

  **Peak agreement 1.0013**, collapse drift within 0.004%. They also tracked each other to within
  ~1% across fifteen progress points of the ascending branch.
  **So the Aldemir collapse is NOT a Newton convergence failure** — that was the last remaining way
  it could have been numerical, and it is now closed. Nine experiments (damping, Gf, ringing,
  topology, mesh, constitutive law, reinforcement, thickness, integration scheme) have failed to
  remove it; it is a property of the model.
  Note also the residual is LOWER under explicit (0.97% vs 1.56%) despite dropping a damping term,
  which is the strongest evidence yet that the ~904 kN peak is genuine resistance and not damping.

- **WHY IT MATTERS BEYOND SPEED.** `ElasticMultiLinear` is PATH-INDEPENDENT, so a strut sitting on a
  segment boundary flips branches between Newton iterations and the tangent never settles — 1,510
  solver failures against Concrete02's 19 in the D73 trilinear run, and total non-convergence in the
  replica. An explicit march never forms a tangent, so that pathology cannot arise. This is why
  Aydin runs his own model explicitly at dt = 5e-8 s: it is not an arbitrary choice, it is the only
  scheme his constitutive law tolerates at scale. It also makes `replica/` runnable again, so
  "does HIS model collapse, or does his curve depend on his solver" becomes answerable.

- **Status:** `--explicit [INTEGRATOR]` on `examples/aydin_aldemir_wall/pushover.py` (run directories
  tagged with the integrator name). Verified by the cross-validation above. NOT yet used on
  `replica/`, and not yet wired into `run_cyclic_dynamic` — the cyclic runner still has only the
  implicit path, and `ElasticMultiLinear` remains unsuitable for cyclic work regardless (D60).

### D75 — 2026-08-29 — Bond re-enabled and made to work: the fix is a FORCE-SLIP law, because a truss couples stiffness and strength through the area. Peak 1,027 kN (1.07x measured) — and the collapse still arrives, at 0.237% drift

- **Context (user):** bond was disabled 2026-08-28 (D72) and re-enabled today after the replica (D73)
  reached only 0.42 of Aydin's published peak WITHOUT bond, making it the leading suspect.

- **THREE ATTEMPTS, THREE DIFFERENT BUGS, each found by a run rather than by reading code.**
  1. **Full-area links (D72).** Bond given the concrete strut's area and modulus: the ring around each
     steel node became a stiff parallel path and the model came out **2.37x stiffer than perfect
     bond** — impossible, since bond can only add flexibility.
  2. **Area calibrated for stiffness.** Sweeping the bond area against the perfect-bond elastic
     stiffness gives **0.01*A_t -> K/K_perfect = 1.0058**, a narrow window (at 0.003 the tangent goes
     singular, `U(i,i)=0` — the steel nodes lose restraint, exactly the fragility predicted in D72).
     But scaling the area down 100x scaled the STRENGTH down 100x too: `F_cr` fell from the strut's
     11.70 kN to **0.117 kN**. The ring cracked at once, the reinforcement disconnected, and the wall
     collapsed at **0.026% drift — 7.6x EARLIER than with no bond at all.**
  3. **Force-slip law — the fix.** In a truss `k = EA/L` and `F = ft*A`, so **AREA IS THE WRONG SINGLE
     KNOB**: it cannot satisfy stiffness and strength together. `materials.bond_elastic_brittle` now
     takes `peak_force` / `area` / `length` / `slip_at_peak` and sets them independently.

- **THE REAL DEFECT WAS UPSTREAM OF ALL THREE.** The original law inherited the CONCRETE grade's
  cracking strain, `ft/E` = 7.4e-5 — **0.0037 mm of slip over a 50 mm link**. A deformed-bar
  interface fails at `s1` ~ 0.1-1 mm (Model Code), two orders of magnitude larger. Sized properly, at
  `A_bond = 0.01*A_t` the link needs `ft_bond` = 185 MPa to carry the strut's `F_cr`, which puts its
  failure slip at **0.372 mm** — inside the Model Code range. `build._bond_material` derives
  `slip_at_peak` so the link's modulus comes out at the concrete's `E`, which is what the area
  calibration assumed, so both constraints hold at once. Elastic check: **K/K_perfect = 1.0129**.

- **RESULT (mesh 50, t = 210, explicit CentralDifference, zeta = 0.5):**

  | | peak | at drift | collapse | vs measured 963.6 | residual asc p95 |
  |---|---|---|---|---|---|
  | no bond | 904.6 kN | 0.1847% | 0.2030% | 0.939 | 0.97% |
  | **with bond** | **1,026.6 kN** | 0.2325% | **0.2374%** | **1.065** | 1.51% |
  | | **x1.135** | | **+21%** | | |

  The residual confirms the peak is real resistance, not inertia. The model now OVERSHOOTS the
  measured force by 6.5%, having undershot by 6% — and is closer than Aydin's own 1,164 kN (+21%).

- **BOND IS NOT THE MISSING PIECE EITHER.** He reports no collapse to 20 mm (0.75% drift on his
  panel); we reach 0.237%. Bond is the TENTH experiment and lands where thickness did — a real,
  substantial improvement that does not remove the failure.

- **WHAT HAS EVER MOVED THIS MODEL, after ten experiments:**
  | | peak | drift capacity |
  |---|---|---|
  | thickness 120 -> 210 (D73) | x1.53 | x1.9 |
  | bond (this entry) | x1.14 | x1.21 |
  | the other eight | — | timing only |
  Damping, Gf, crack-release ringing, lattice topology, mesh size, constitutive law, reinforcement
  presence and integration scheme all changed WHEN the wall failed, never THAT it did.

- **STILL UNPINNED, and all of it is bond parameters the paper never prints:** he gives the residual
  ratio `a` = 0.7 and nothing else — our bond area, slip and strength are inferred. Also unpinned:
  `b1`/`b2` of the trilinear tail, his PID force drive against our prescribed displacement, and his
  dt = 5e-8 s against our 8e-6.

- **CAVEAT THAT STANDS.** The horizon-star topology still cannot satisfy "rigid bond == perfect bond"
  in the limit; calibration bounds the artefact (0.6-1.3% on elastic stiffness) rather than removing
  it. A coincident-node `zeroLength` interface spring would, and remains unbuilt.

- **Status:** `--bond --bond-area-ratio 0.01` on `examples/aydin_aldemir_wall/pushover.py`;
  `tests/test_bond.py` still carries the xfail limit test, which is still correct to xfail.

### D76 — 2026-08-30 — Bond wired into the Aydin replica: the last big departure removed, the model doubles, and the area ratio does NOT transfer from mesh 50 to mesh 20

- **Context (user):** "add bond to replica". D75 named this the obvious next run and left it undone.
  `replica/specimen.NOT_REPLICATED` had called shared nodes "the single largest known departure and
  the first thing to suspect if the replica behaves differently"; without bond the replica reaches
  only **0.42 of his published peak** and collapses at 0.098% drift.

- **What was added.** `replica/build.py` gains `_bond_material` (the D75 force-slip law, his
  residual `a` = 0.7 from `specimen.BOND_RESIDUAL`), `report_bond`, and `wall_lattice(bond=,
  bond_area=)`; `run.py` and `preflight.py` gain `--bond` / `--bond-area-ratio`. Same scheme as the
  parent, unchanged: separate steel nodes, a horizon ring of elastic-brittle links,
  `i_accept_the_known_bond_defect=True`.

- **VERIFIED BY BUILD, NOT BY SOLVE** (no FE run launched, per standing instruction):

  | | shared nodes | bond |
  |---|---|---|
  | nodes | 20,385 | 28,081 |
  | elements | 88,324 | 149,802 |
  | of which bond | — | 61,478 |

  **His Table 2 counts are untouched** — 20,385 / 80,684 concrete struts, still exact — confirming
  again that his published element count is concrete only. Links come out at **7.99 per steel node**
  against the ring of 8 the horizon predicts (the shortfall is the panel edges), nodal mass is
  conserved to 9 digits by `bond_mass_share`, and the minimum nodal mass stays positive, so no DOF
  is left massless for the dynamic runners.

- **THE AREA RATIO DOES NOT TRANSFER ACROSS MESH SIZE.** The failure slip of a link is
  `ft/(ratio*E) * length`, i.e. **proportional to the link length**, so the parent's calibrated
  `0.01*A_t` gives **0.372 mm at mesh 50 and 0.149 mm at the replica's mesh 20**. Both are inside
  the Model Code `s1` band for a deformed bar (0.1–1 mm), but the replica sits on its FLOOR, and
  lowering the ratio to raise the slip walks toward the singular window (0.003 gave `U(i,i)=0` at
  mesh 50). The window is narrow and is **not measured at this mesh** — `--elastic --bond` against
  `--elastic` is the check, and `report_bond` prints the slip and flags it if it leaves the band.
  The other bond numbers are ratio-only and so mesh- and thickness-independent: `ft_bond` = 185 MPa,
  `E_bond` = `E_c` by construction.

- **Two smaller corrections that came with it.**
  - `NOT_REPLICATED` became `not_replicated(bond=, explicit=)`. A static list makes a run overstate
    its own departures: with `--bond --explicit` two of the four entries no longer apply, and the
    bond entry has to be *replaced* rather than deleted — the elements are now present, but his link
    area, peak slip and peak force are still unpublished, and the ring topology still cannot reduce
    to perfect bond in the limit (D72). Those are different claims and the printout now says so.
  - `explicit_steps_per_period` moved into `build.py` from a local block in `run.py`, so
    `preflight.py` can price an explicit run the same way the run itself sizes one. Bond links take
    the CONCRETE modulus in `element_modulus`, which is exact rather than approximate: the slip
    derivation sets `E_bond == grade.E`. **`preflight.py` must be given the same flags as `run.py`**
    — an implicit probe misprices an explicit run by ~20x in step count, and the no-bond model does
    not price the bond one.

- **Not run.** The staged sequence is `preflight.py --bond --explicit`, then `run.py --elastic
  --bond` against `run.py --elastic`, then `run.py --drift 0.0025 --bond --explicit`. Budget past
  the no-bond collapse of 0.098%, not to it: in the parent, bond moved collapse drift x1.21.

- **Status:** wired and build-verified; the D72 library guard and its xfail limit test stand
  unchanged.

### D77 — 2026-08-30 — Fig. 10(b) digitized after all: "not recoverable" was over-general. The replica reproduces AYDIN'S OWN CURVE to 1.018 at matched displacement, and the 0.42 / 0.746 shortfalls were two separate measurement errors

- **Context (user):** told plainly that the paper provides hysteretic data Aydin compares against.
  D72 had recorded Fig. 10(b) as "an unordered raster point cloud, both standard routes return
  noise", and that note had been propagated into `testdata.NOT_REPORTED` and the handoff. Checking
  the figure instead of the note showed the note was half right and had been generalized too far.

- **WHAT IS ACTUALLY LOST IS *ORDER*, NOT THE DATA.** The experiment is drawn as discrete dot
  markers, so which dot follows which is gone — per-cycle energy, degradation and the load path are
  unrecoverable, exactly as for SW-NC-FF's Fig. 14b. But a backbone is an ENVELOPE, extremes survive
  unordered, and the figure additionally carries **both of Aydin's own lattice curves as clean
  coloured lines**. The panel is a raster (the page's 257 SVG paths are glyph outlines), so this is
  a pixel job — but unlike WSH3's grayscale panels it has real COLOUR separation.

- **`examples/aydin_aldemir_wall/digitize.py`.** Axes are fixed from detected gridlines alone — nine
  load rules, eight displacement rules — with no published value entering the calibration. Dots are
  separated from tick-label glyphs morphologically (a dot is a ~6 px blob, a glyph stroke 2-3 px, so
  an opening keeps one and drops the other), then by position.

  | | digitized | published | ratio |
  |---|---|---|---|
  | Aydin horizon 1.5 plateau | 1158.0 kN | 1164.413 | **0.9945** |
  | Aydin horizon 3.01 plateau | 1305.3 kN | 1325.675 | **0.9846** |
  | cloud envelope peak | 969.4 kN | 963.592 | **1.0060** |
  | "experiment" marker | 924.2 kN | 963.592 | 0.959 (coarse) |

  Four numbers the calibration never saw, so these are tests rather than fits. The marker is the
  weak one — its position is assumed and the disc touches neighbouring dots.

- **THE ENVELOPE IS AN UPPER BOUND ON THE BACKBONE, NOT THE BACKBONE.** From an unordered cloud a
  boundary point may be a loop tip or the unloading side of a larger cycle. A first attempt made it
  monotone by accumulating outward from the origin; that was DISCARDED because it starts at a flat
  ~430 kN at zero displacement, which is not a capacity — near the origin every loop crosses, so the
  spread there is loop-crossing scatter. The raw envelope is ragged and honest. Two further losses:
  the figure's legend occludes real data at positive displacement below about -500 kN, and dots that
  merge INTO a tick label are sacrificed, costing interior dots in the dense pinch (not the envelope).

- **THE COMPARISON THIS ENABLES CORRECTS TWO ERRORS AT ONCE.**
  1. **t = 120 vs 210 (already caught).** The aborted no-bond run's `ABORTED.md` claimed t = 210
     while its own console showed `A_t = 1,486.8` = t = 120. Its "0.42 of his peak" was measured
     43% too thin.
  2. **A LIKE-FOR-LIKE ERROR THAT SURVIVED THE FIRST CORRECTION.** Replacing 0.42 with "0.746 of his
     peak" compared our 4 mm ENDPOINT against his 16 mm PLATEAU — different displacements. Against
     his actual digitized curve at MATCHED displacement the t = 210 control gives **mean 1.018 over
     0.5-4.0 mm (min 0.936, max 1.126)**. The replica is not undershooting him at all; it simply has
     not been pushed as far. Both shortfall figures were artefacts of measurement, not of the model.

  Against the measured envelope the control sits at 0.85-0.91 over the same range — and since that
  envelope bounds the backbone from ABOVE, the true gap is smaller.

- **METHOD NOTE, the general form of error 2.** Comparing a truncated run's endpoint against a
  reference's plateau is not a comparison. Where a reference CURVE exists, compare at matched
  abscissa; where only a peak exists, say explicitly that the run must reach the same displacement
  before the ratio means anything. Both of this study's headline shortfalls came from ignoring that.

- **Status:** `digitize.py` -> `data/fig10b.npz`; `replica/compare.py` plots every run against the
  digitized references and falls back to Table 4 levels if the npz is absent. `testdata.py` and the
  handoff corrected. The bond-ratio experiment (D76) is still running and is reported separately.

### D78 — 2026-09-05 — The Aldemir collapse is a LOAD-INTRODUCTION artefact: the bars stop one spacing short of the driven row. Bars to the top remove it — supersedes "still open after ten experiments" (D73/D74/D75)

- **Context:** six replica runs of 2026-08-31 to 09-03, plus a re-reading TODAY of the parent
  study's own saved displacement fields (no run: `data.json` already stores `disps_peak` and
  `disps_final`).

- **THE MECHANISM, from the snapshot ring.** Three replica configurations were instrumented
  through their collapse. All three peak with damage correctly distributed — worst struts at
  y = 40-450 mm, the flexural base of a wall in flexure — and then, within two progress ticks,
  tear along the row the displacement is prescribed to, after which everything below unloads:

  | run (t = 210, mesh 20) | peak | tore at | where the worst struts sit after |
  |---|---|---|---|
  | horizon 3.01, no cap | 785.9 kN @ 0.1001% | 0.1075% | y = 2630-2670 of 2680 |
  | horizon 1.5, elastic cap of 3 rows | 846.1 kN @ 0.1414% | 0.1441% | y = 2550-2610 (just under the cap) |
  | horizon 1.5, BASELINE | 866.6 kN @ 0.1504% | 0.1531% | y = 2570-2650 |

- **THE CAUSE IS THE BAR GRID, and it is present in BOTH models.** Bar lines are placed on a
  cover + spacing grid, so the top line lands one spacing short of the top face and leaves a band
  of PLAIN CONCRETE that the whole base shear must leave through: **replica y = 2580 of 2680,
  parent y = 2150 of 2250 — 100 mm in both.** The elastic cap did NOT fix it because it removed
  the failure criterion of the top rows without changing their strut areas, so the demand simply
  relocated to the first failable row beneath.

- **`--rebar-to-top` REMOVES THE COLLAPSE.** Same model otherwise:

  | replica, bars to the top | result |
  |---|---|
  | horizon 1.5 | peak 897.1 kN at 4.34 mm, traced to the full 0.30% target, converged |
  | horizon 3.01 | 861.2 kN at 8.0 mm, still ascending at target, converged |

  Against the digitized references at MATCHED displacement (D77's method note): 743 kN at 2.68 mm
  against his 760 and the test envelope's 847; 866 against 905 and 969 at 4.0 mm.

- **THE PARENT STUDY HAS THE SAME DEFECT, measured rather than inferred.** Its three t = 210 runs
  store the field at peak and at the last step. Share of concrete struts past `eps_ult`, by height:

  | run | at PEAK: base <750 / middle / top 150 | at the LAST STEP | worst 1% then sit at |
  |---|---|---|---|
  | implicit Newmark, no bond | 69.5% / 27.5% / 0.0% | 54.5% / 0.7% / 44.8% | y = 2175-2225 |
  | explicit, no bond | 71.4% / 28.6% / 0.0% | 55.9% / 1.0% / 43.1% | y = 2175-2225 |
  | explicit, bond 0.01 | 67.9% / 23.9% / 1.2% | 0% / 0% / **100%** | y = 2225-2225 |

  The wall is behaving correctly at peak in every case and tearing at the driven band afterwards.
  The saved `damage.png` of both explicit runs shows the same single red row with the wall below
  unloaded, and had been read as "the wall lost its load path".

- **WHAT THIS SUPERSEDES.** D73/D74/D75 report the premature collapse as a property of the model
  surviving ten experiments. Nine of those — damping, Gf, ringing, topology, mesh, constitutive
  law, reinforcement presence, integration scheme, and the elastic cap — cannot affect a boundary
  detail, which is why none of them moved it. What thickness (x1.53 peak) and bond (x1.14) bought
  is unaffected and stands.

- **NOT MEASURED: the parent wall with bars to the top.** `specimen.rebars` has no such option
  (the replica's is `full_height_rebar`). The parent explicit pushover is 95,310 steps at
  6.8 ms/step = **11 minutes**, so this is the cheapest open question in the study.

- **PROCESS.** The 2026-08-30 control was given `--drift 0.0015`, ended at exactly 0.1500%
  "still ascending, converged, no turnover", and was ~3,000 steps from the cliff. That is the
  D68/D73 trap one level up: a run that ends at its REQUESTED target says nothing about behaviour
  past it.

- **Status:** accepted. Example-layer only; no backend change. The replica carries `--rebar-to-top`
  and `--cap`; the parent carries neither.

### D79 — 2026-09-05 — Cyclic on the replica: `HystereticSM` cannot pinch and `Concrete02` can; the push/pull asymmetry was the MATERIAL, not the wall

- Two reversed-cyclic runs, same wall (t = 210, mesh 20, horizon 1.5, bars to the top, explicit
  CentralDifference), same protocol (0.05 / 0.10 / 0.15% drift, one cycle each), differing only in
  the cyclic concrete carrier. `concrete_lattice_aydin` itself is unusable here: it is
  `ElasticMultiLinear`, path-independent, so a cracked strut recovers full stiffness on reload.

  | | force retained at zero displacement, as a fraction of that cycle's peak | \|pull\|/push at level 1 |
  |---|---|---|
  | `HystereticSM` (his trilinear tail) | 0.56 / 0.36 / 0.72 | 1.103 |
  | `Concrete02` (regularized) | 0.23 / 0.05 / 0.40 | 0.915 |

  Predicted first at the single-strut level — a cycled strut holds **4.38 MPa at zero strain**
  under `HystereticSM` against **0.29 MPa** under Concrete02 — and then confirmed structurally.
  **The pinch parameters are inert for this**: 0.3/0.1 -> 0.10/0.03 moved the residual 4.38 -> 4.43
  MPa. They shape the APPROACH to zero, not the value at zero, and the value at zero is what
  decides whether a wall pinches. An earlier claim that loop shape "is governed by the pinching
  parameters" was wrong.

- **The asymmetry reverses with the material**, which settles what it was. `HystereticSM` made the
  pull half 10% STRONGER, attributed at the time to crack closure carrying compression on
  reversal; Concrete02 on the same wall, protocol and push-first ordering gives 0.915, the ordinary
  direction for a wall already damaged by the push.

- **The completed run:** 1,340,553 steps, 22.0 h at 59 ms/step, converged, three full loops to
  +-4.02 mm, peak **+744.5 / -671.7 kN = 0.77 of the measured 963.6**. Damage rises to 17.6% of
  struts cracked at the largest tip and falls back to ~5% at each zero crossing — that metric is
  the fraction currently past `eps_cr`, i.e. a STATE, not an accumulated damage count.

- **His unloading is ORIGIN-ORIENTED** ("The unloading rules are origin-oriented for the tension
  model and elastic for the reinforcement", Fig. 1(c)/(d)). `HystereticSM`'s `beta` degrades the
  unloading stiffness with ductility and beta = 1 is the secant through the origin; both runs used
  **beta = 0**. Not tested, and it is a one-strut check.

- **Status:** accepted. New replica tooling: rolling `hysteresis.npz`, a snapshot ring plus
  archived `frames/`, `hyst.py`, `damage_history.py`, `damage_extra.py`, `report.py`, `phone.py`,
  `status_page.py`. No backend change.

### D80 — 2026-09-05 — Seven corrections to this study's record of the 2019 paper and of its own bond implementation

1. **b1 AND b2 ARE PUBLISHED.** The Fig. 2 flowchart box prints `b1 = 0.6, b2 = 0.2, a = 0.7`.
   `replica/specimen.py`, D75, D76 and HANDOFF.md all list them as never printed. The published
   values equal the repo's defaults, so no result changes — but "unpublished bond and tail
   parameters" is a shorter list than claimed.

2. **BOND IS A CONCRETE STRUT WITH A RESIDUAL PLATEAU (author, 2026-09-05).** Same area and same
   `F_cr` as a concrete strut; once it reaches capacity the strength drops to ~60% and is
   perfectly plastic. This is what Fig. 1(c) draws — Bond and Concrete rising on one line to the
   same `F_cr` — and `materials.bond_elastic_brittle(grade, tag, residual=0.6)` at the DEFAULT
   full area is exactly it (verified: elastic to 1.85 MPa at `eps_cr` = 7.44e-5, drop to 1.11 MPa,
   flat, symmetric). **The D75 force-slip law was needed only because the area had been scaled
   100x to fix stiffness; against the author's description it is a departure, not a fix.**

3. **Consequently the D72 "2.19x too stiff" is not an implementation error to calibrate away.**
   It is what a full-area ring does at THIS bar density, and the density is the variable:
   **73.5% of the parent's nodes at mesh 50 carry a steel node, against 35.3% of the replica's at
   mesh 20** (Ø8@100 against a 50 mm and a 20 mm grid). The artefact should be markedly milder on
   his own model than on ours, which is consistent with his `K_sim` = 0.91x the measured stiffness
   rather than 2x it. NOT MEASURED — the elastic A/B at mesh 20 has not been run.

4. **BASE ANCHORAGE DEFECT, new and measured.** `build_lattice_rc` resolves supports on the
   CONCRETE coordinate array, so a steel node can never be fixed: on the bonded parent wall,
   **30 steel nodes sit at y = 0 and none of them is fixed**, each held only by its own ring of
   links, which fail at `ft/E * mesh` = **0.0037 mm** of slip. Meanwhile `select_nodes` queries the
   BUILT model, so it does pick up coincident steel nodes — with bars run to the top the drive
   would be imposed on steel at the top while the base stays free. Every bonded result to date
   carries this, and it is the most likely reason bond made the replica worse (peak 303 kN at
   0.4 mm at both area ratios) while the elastic gate passed.

5. **`build.strut_groups` classifies bond links by orientation**, so they are counted as concrete
   vertical/diagonal. The bond run's "rebar 0.7% of the overturning moment" against 21% without
   bond is therefore a probe artefact, not a finding about the bars.

6. **`run_cyclic_dynamic` reports `("HHT", 0.7)` unconditionally** even when marched explicitly,
   and keeps the implicit sub-step rescue ladder under an explicit integrator — against D74's rule
   that a failed explicit step is instability and must not be retried. Both cyclic runs above are
   labelled HHT in their `data.json` and were CentralDifference.

7. **`replica/specimen.TW` still defaults to 120 mm** via `ALDEMIR_TW`; every run since
   2026-08-29 set 210 by hand, and one ABORTED.md already recorded a 0.42 ratio measured at the
   wrong thickness.

- **Status:** recorded, none of them fixed in code yet. Items 4, 5 and 6 must be fixed before any
  further bonded or cyclic run is worth launching.

### D81 — 2026-09-05 — The study programme agreed with the author: two concrete models x bond x three analyses, driven from one parameterized entry point

- **Context:** the author (Aydin, the paper's first author and the user's supervisor) set the
  programme directly: (i) pushover with the OLD method — linear compression, tension-softening
  calibration; (ii) pushover with a COMPRESSION-CAPABLE lattice, compression calibration after the
  tension calibration; (iii) cyclic response with both. Plus: do not forget bond, and do not forget
  regularization. Answers to follow-up questions settled: no access to the 2017 test paper (use
  Aldemir's results as reported by the 2019 paper, i.e. Table 4 plus the digitized Fig. 10(b));
  compression is **elastic-perfectly plastic** and must be named in the model names; tension keeps
  his softening law in both models; bond is the concrete-strut law of D80 item 2.

- **Two concrete models, named for their compression branch:**
  - `lincomp` — his published model: tension-only trilinear softening, LINEAR compression forever.
  - `eppcomp` — same tension, compression capped at a calibrated strength and PERFECTLY PLASTIC.
    New: no existing material does this. It removes compression softening entirely, so `Gfc` and
    the crushing strain stop being parameters and regularization applies to tension only.

- **The compression calibration is direction-dependent and that is a real open question.** On a
  horizon-1.5 lattice the axial capacity across a horizontal cut is the vertical struts' share,
  **0.612 of the gross face at mesh 50** (factor 1.633), while an inclined family along 45 degrees
  presents **0.852**. A single scalar fitted on an axially-loaded cube therefore over-corrects the
  diagonal compression field a squat wall actually carries, by roughly 1.4x. Recorded before any
  run rather than discovered after one.

- **The tension calibration is ALSO not settled, and the paper's own numbers say so.** His Fig. 2
  loop fits `a1, a2, a3` until a simulated coupon dissipates `Gf` within 15%; the repo instead
  SOLVES `a2, a3` in closed form from `L * integral(sigma d eps) = Gf` per strut. The two disagree:
  his published `a2 = 70, a3 = 360` at `d = 20` dissipate **2.08x Gf** by that identity, and the
  identity at this grade gives `a2 = 33.7 / a3 = 173` at L = 20 and `a2 = 13.6 / a3 = 70` at
  L = 50. Either his loop measures energy over a gauge spanning several cracking struts, or his
  tail is deliberately ductile. Both routes are implemented and the choice is a matrix axis.

- **Cost, from MEASURED explicit runs, decides the shape of the matrix.** Parent mesh 50 explicit:
  6.8 ms/step, 23.4 with bond; replica mesh 20: 59-67 ms/step. An eight-level cyclic ladder to 1.0%
  drift at 7.6 mm/s is ~4.9M steps: **9 h at mesh 50 without bond, ~28 h with, and 210-250 h at
  mesh 20**. So cyclic work lives at mesh 50 and mesh 20 stays monotonic.

- **Layout:** one parameterized entry point under `examples/aydin_aldemir_wall/study/`, run
  directories named from the resolved parameters, a `params.json` per run carrying EVERY parameter
  including defaults, a per-run report, and a master report that discovers the matrix by reading
  those files rather than by parsing names — so new parameters and new values extend the matrix
  without touching the reporting.

- **Status:** agreed, not implemented. Blocked on the D80 fixes (bond anchorage, probe
  classification, cyclic runner provenance) and on the parent gaining `--rebar-to-top`.

### D82 — 2026-09-05 — Fig. 4(f) shows the BOTTOM 1500 mm of the model, so its "1500" is a crop height and not a panel dimension — corrects D72

- **Source:** the paper's first author, 2026-09-05, told plainly when the study plan offered
  `2680 x 1500` as a candidate geometry.

- **What D72 recorded.** That Fig. 4(f) "labels the analysed model 2680 x 1500", that it is "NOT to
  scale" because its drawn block has an aspect of 2.235 "matching none of its own candidate labels
  (2680/1500 = 1.787 is the closest, 25% off)", and that it therefore "cannot arbitrate geometry".
  `specimen.PANEL_FIG4F = (2680.0, 1500.0)` exists for that reading and `summary.py` prints it as a
  candidate panel; `elastic.py --geometry-sweep` offers it as a third row.

- **What it actually is.** A DETAIL VIEW of the lower part of the model, in the same style as the
  zoom insets beside it. The 1500 is the height of the crop, not of the wall. A cropped view carries
  no panel dimension at all, so it was never a candidate and the sweep should not offer it.

- **This also DISSOLVES the anomaly rather than explaining it.** D72 inferred "not to scale" from
  the aspect mismatch and used that to disqualify the figure. The mismatch is what a crop looks
  like, so no drafting error needs to be assumed, and the figure's status changes from "unreliable"
  to "reliable about something else".

- **ONE CONSEQUENCE IS NOT SETTLED AND IS FLAGGED RATHER THAN ACTED ON.** Inverting Table 2's counts
  gives 150 x 134 cells, i.e. 3000 x 2680 mm, and the inversion is UNIQUE ONLY UP TO TRANSPOSITION
  (D72 states this and reports the pair in both orders in different places). `replica/specimen.py`
  resolves it as `LW, HW = 3000, 2680` — wider than tall. If the 2680 dimensioned along the
  HORIZONTAL in Fig. 4(f) is the model's full width, as a full-width crop implies, the resolution
  runs the other way: **2680 wide x 3000 tall**, which is 134 x 150 cells and equally consistent
  with both Table 2 integers. Every replica result to date assumes the first. Not changed here: it
  needs the author's confirmation, and it would move the drift denominator, the aspect ratio and
  every published replica ratio.

- **Status:** recorded. The plan drops `fig4f` as a panel option (only `fig10a` and `table2` remain);
  `specimen.PANEL_FIG4F`, `summary.py`'s candidate list and `elastic.py --geometry-sweep` still carry
  the old reading and should be corrected when that file is next touched.

### D83 — 2026-09-05 — The study plan inspected against the code: D80 item 6 FIXED (three defects, not one), the bonded element counts settled by a build, and the cyclic carrier question closed against `HystereticSM`

- **Context:** `examples/aydin_aldemir_wall/study/PLAN.md` (D81) read against the source rather than
  against the record. Every defect it lists is real and was reproduced; four of its statements are
  amended in place, and one backend fix follows. No solver run except one model build.

- **`run_cyclic_dynamic` FIXED, and it was worse than D80 item 6 recorded.** Three changes:

  1. the sub-step rescue ladder is now gated on `explicit`. It was not merely running under an
     explicit march — it ends with an unconditional `ops.algorithm("Newton")`, so ONE failed
     explicit step converted the remainder of the march to an implicit solve on a `Diagonal`
     system, silently and for every step after it;
  2. the output reports `tuple(integrator)` instead of the literal `("HHT", 0.7)`;
  3. `capture` tested `shear[-1]` BEFORE this step's shear was appended, so `disps_peak` held the
     displacement field of step i against the shear of step i-1 — one step stale, in the field
     D78 used to diagnose the Aldemir tear.

  `run_pushover_dynamic` already did all three correctly; the cyclic runner was the stale copy, so
  the fix is a transcription. **D80's "items 4, 5 and 6 must be fixed before any further bonded or
  cyclic run" now reads 4 and 5.** Docstring corrected too: it claimed a Newmark march (the default
  is HHT(0.7)) and listed TWO automatic changes under an explicit integrator, now three.

- **ELEMENT COUNTS SETTLED BY A BUILD, and D72 was wrong.** At mesh 50, horizon 1.5 the parent wall
  builds as **2,806 nodes / 13,471 elements** (10,905 concrete + 1,290 longitudinal + 1,276 stirrup)
  and **5,424 / 34,325 with bond** (+20,854 links). D72 and CLAUDE.md gave 5,454 / 13,501 / 34,595.
  The old numbers are 30 / 30 / 270 higher and THEIR PROVENANCE IS UNEXPLAINED. "Bars run to the
  top" was the obvious candidate — 30 is the number of vertical bar lines — and it is REFUTED by
  building it (below): that option costs +60 elements, because the top bar line sits at y = 2150,
  two mesh-50 cells below the top face, so each bar gains TWO segments. A +30 delta needs bars
  stopping ONE cell short, i.e. a cover other than the one now modelled; the directory is untracked,
  so the history is not recoverable. CLAUDE.md corrected to the measured numbers. Counts do not
  depend on `--tw`. (Recorded because the arithmetic fit was reported as a decomposition before it
  was built — the D68/D73/D75 trap, one level down.)

- **The plan's cyclic-carrier section was written as though D79 had not happened.** It names
  `HystereticSM` "the available vehicle" and schedules a single-strut probe as a prerequisite; D79
  ran that probe AND a 22 h structural confirmation, and the answer is that `HystereticSM` does not
  pinch (4.38 MPa at zero strain) while `Concrete02` does (0.29). So Stages 5-6 have no validated
  material as written, and `c02` moves from an extension slot to the default cyclic carrier.
  **A consequence the plan now records:** `Concrete02` brings its own capped, softening compression
  branch, so adopting it partly COLLAPSES the `comp` axis for the cyclic cells — a cyclic cell
  labelled `lincomp` would not have linear compression. The pushover stages are unaffected. What
  stays open is D79's own one-strut item: whether `beta = 1` (origin-oriented, as Fig. 1(c)) rescues
  the Aydin envelope.

- **Three smaller amendments to PLAN.md.** Stage 0 item 1's fix must cover the DRIVE end as well as
  the support end (with bars to the top the anchorage defect becomes both) and be asserted in the
  build; Stage 1's gate must be run with a target well past the 0.237% it has to beat, or it ends at
  the number it is being asked to exceed (D78's own PROCESS trap); Stage 4 gates on peak force, load
  path and damage location, NOT on drift at peak, which §1 already lists as an unfair claim and
  which would be measured against a record the frame clips at +16 mm. Also recorded: `eppcomp` is a
  two-point edit to the `neg_eps`/`neg_sig` pair of `materials.concrete_lattice_aydin`, not a new
  constitutive object.

- **Status:** accepted. Backend change in `opensees.py` (`run_cyclic_dynamic` only); documentation
  changes in CLAUDE.md and PLAN.md. Nothing in the study is built yet, and no analysis has been run.

### D84 — 2026-09-05 — Bars-to-top ported to the parent Aldemir wall (D78's cheapest open question), and what it costs

- `specimen.rebars(..., full_height=True)`, `build.wall_lattice(..., full_height_rebar=True)` and
  `pushover.py --rebar-to-top` (tagged `_rebartop` in the run directory), mirroring the replica's
  `full_height_rebar`. Default unchanged: the paper says nothing about bar termination, so both the
  grid placement and this are inferences, and what is measurable is how much the answer depends on
  which is chosen.

- **COST, built and counted (no solver), mesh 50 / horizon 1.5:**

  | | nodes | elements | longitudinal |
  |---|---|---|---|
  | default | 2,806 | 13,471 | 1,290 |
  | `--rebar-to-top` | 2,806 | 13,531 | 1,350 |
  | `--rebar-to-top --bond` | 5,484 | 34,775 | 1,350 |

  **+60 elements, not +30**: the top bar line is at y = 2150, TWO mesh-50 cells below the top face
  (D78 measures that band as 100 mm), so each of the 30 bars gains two segments. With bond it is
  +60 steel nodes and +450 elements — the 60 new steel nodes carry 390 ring links between them
  rather than 60 x 8, because a ring truncates at the top face.

- **Stage 1 of the study plan is now unblocked** and is the only Stage-0 item it needed: it is a
  no-bond run, so the anchorage defect (D80 item 4, steel nodes exist only under bond) and the
  `strut_groups` mislabel (item 5) cannot reach it, and item 6 is fixed (D83). Not run here.

- **Status:** accepted. Example layer only. `elastic.py` and the continuum twin still call `rebars`
  without the option — deliberate, since the elastic gate is insensitive to bar termination, and
  worth revisiting if Stage 2 is ever compared across the two terminations.

### D85 — 2026-09-05 — Bond anchorage FIXED at both ends, and bond links stop being counted as concrete: D80 items 4 and 5 closed

- **Item 4, the support end (library).** `_apply_supports_loads` now propagates a support to the
  steel nodes duplicated at that concrete node (`steel_of_concrete`), so the Aldemir wall's 30 base
  bar nodes are anchored instead of being held only by a ring of bond links that fails at 0.0037 mm
  of slip. Measured: supports 61 -> 91 with bond, **0 unanchored base steel nodes** where there were
  30. Fixing the duplicate is also the physical reading — a bar crossing the support row is anchored
  into the foundation.

- **Item 4, the DRIVE end (the half D80 did not state).** With bars run to the top (D84) the same
  defect appears at the loaded row: a coordinate box returns both the concrete node and its steel
  duplicate, so the drive would have been imposed on the bar as well. `select_nodes` gains
  `kind="any" | "concrete" | "steel"` (default `"any"`, so every existing caller is unchanged), and
  the Aldemir specimen now picks the two sets DELIBERATELY and in opposite directions:

  | set | kind | why |
  |---|---|---|
  | `top_nodes` / `drive_nodes` / `control_node` | `concrete` | the actuator loads concrete; the bars follow through bond. Driving the duplicate pins ZERO SLIP at exactly the row whose load introduction is the question (D78) |
  | `base_nodes` | `any` | the anchored steel nodes now carry part of the base shear; dropping them loses the bars' contribution to the reaction sum |

  Verified with bond + bars-to-top: 30 steel nodes exist at the top row and **none is driven**;
  61 drive nodes, all concrete; the control node is concrete.

- **Item 5 (example layer).** `build.strut_groups` gets a `"bond"` group. Bond links sit at every
  orientation, so falling through to the geometric branches counted them as vertical/horizontal/
  diagonal struts — the reason a bonded run reported "rebar 0.7% of the overturning moment" against
  21% without bond. That number was a probe artefact and should not be quoted.

- **No behaviour change without bond**, checked rather than assumed: 61 supports / 61 base / 61
  drive nodes and no `bond` group, exactly as before. `Model` gains `steel_nodes: set[int]`, empty
  unless the bond scheme is used. Full suite 55 passed, 1 failed — the known D34 failure.

- **D80's "items 4, 5 and 6 must be fixed before any further bonded or cyclic run" is now closed**
  (6 in D83, 4 and 5 here). Every bonded result before today carries the anchorage defect and none
  of them should be quoted; the "2.19x stiffer than transformed section" of D72 in particular was
  measured with the base bars unanchored and needs remeasuring.

- **Status:** accepted. Library: `model.py` (one field), `builders.py` (support propagation +
  `select_nodes(kind=)`). Example: `aydin_aldemir_wall/specimen.py`, `build.py`. No analysis run.

### D86 — 2026-09-05 — Stage 0 CLOSED: `eppcomp` built, the cyclic carrier settled against `HystereticSM` at beta = 1 too, and `run_cyclic` was the one runner that never wiped its domain

- **`eppcomp` is a keyword, not a new object.** `materials.concrete_lattice_aydin(..., fc_cap=)`
  replaces the RSM knee with elastic-to-`fc_cap`-then-perfectly-plastic (it raises if `rsm` is also
  set — they are two different compression branches). Verified on the backbone: elastic at exactly
  E = 24,870 to 45.72 MPa = 28 x 1.633, then flat, tension trilinear untouched. Threaded through
  `build.wall_lattice(fc_cap=)` and the study's `--comp eppcomp --fcx`.
  The damage figure is relabelled per PLAN §3 — `figure_damage(crush_label=)` already existed, so
  `eppcomp` draws "yielded in compression" and `lincomp` passes `eps_crush=None`, since compression
  is linear forever there and no compressive limit exists to draw.

- **THE CARRIER QUESTION IS CLOSED, and beta = 1 does not rescue his envelope.** D79 left one
  one-strut check open: his unloading is origin-oriented and `HystereticSM`'s `beta` is that
  parameter, but both of its runs used beta = 0. `study/carrier_probe.py` cycles ONE strut
  (L = 50, A = 6,324.7) to 12 x eps_cr and back, and reads what it holds at zero strain:

  | carrier | stress at zero strain |
  |---|---|
  | `ElasticMultiLinear` (the control) | **-0.000 MPa** |
  | `HystereticSM`, beta = 0 | **-7.70 MPa** |
  | `HystereticSM`, beta = 1 (origin-oriented) | **-6.59 MPa** |
  | `Concrete02`, regularized | **-0.000 MPa** |

  beta = 1 moves it 14%, from 7.70 to 6.59 — the same order, and nowhere near pinching. The control
  reads exactly zero, which is what validates the probe: a path-independent law MUST return to zero.
  (These are larger than D79's 4.38 MPa because the excursion is deeper; the ordering, not the
  magnitude, is the result.) `Concrete02` also unloads at 1.000 x the origin secant, i.e. it
  satisfies his Fig. 1(c) unloading rule better than `HystereticSM`'s own beta does.
  **So Stages 5-6 run on `Concrete02`, and they report that they are not running his envelope.**

- **`run_cyclic` NEVER WIPED ITS DOMAIN — the only model-building runner that did not.** `build()`
  emits into the current domain and does not clear it (its own docstring says "after ops.wipe"),
  and every sibling — `run_static`, `run_modal`, `run_gravity`, `run_pushover`,
  `run_pushover_arclength`, `run_pushover_dynamic`, `run_cyclic_dynamic`, `run_dynamic` — calls
  `ops.wipe()` first. `run_cyclic` did not, so a second call in one process appended to the previous
  model and OpenSees rejected the repeated node tags. Found by the carrier probe, which is the first
  code here to call it twice.

- **The wipe is OPTIONAL, by user instruction, and it skips the build and the gravity stage with
  it.** The `run_gravity` -> `run_cyclic` pairing was ALREADY broken before the fix (double build,
  colliding tags) and the wipe alone would have made it fail differently, by silently discarding
  gravity. `wipe=False` therefore means the CALLER owns the domain and its loads: no wipe, no
  build, and no internal gravity — `_gravity_loads(model, None)` falls back to the model's own
  loads, so applying it would have re-applied gravity and collided on the time-series tag. What the
  branch keeps is the `loadConst`, so displacement control starts from a held gravity state. An
  empty domain under `wipe=False` is refused with a message rather than analysed. All three paths
  verified: gravity-first now works (1,601 steps, converged), back-to-back default runs work, and
  the empty-domain case raises.

- **Stage 2's gate was being reported wrong by the study harness.** It printed K/K_measured; PLAN §8
  gates on **K_lattice/K_continuum**, a plane-stress continuum on the same grid with the same rebar
  and BCs. Now measured per run: **0.9988** (PLAN §5 predicted 0.9985). K/K_measured = 1.2203 is the
  separate geometry question of D73 and is labelled as such.

- **Status:** accepted. Library: `materials.concrete_lattice_aydin(fc_cap=)`,
  `opensees.run_cyclic(wipe=)` + its missing wipe. Example: `build.wall_lattice(fc_cap=)`, the
  study's `models.py`/`run.py`/`report.py`, new `study/carrier_probe.py`. Tests 55 pass, 1 known
  D34 failure. **Stage 0 of the plan is now complete.**

### D87 — 2026-09-05 — Stages 1-4 (nobond half): the COMPRESSION LAW barely matters and the TENSION TAIL dominates. Bond wired per PLAN §4, and its full-area ring measures 2.39x too stiff

- **Stage 3, static diagnostics (4 cells).** All four stall while still ascending, at
  **0.0098-0.0120% drift** — three times earlier than the plan's predicted ~0.03%, and earlier than
  SW-NC-FF (0.3%) or WSH3 (0.067%). `lincomp` and `eppcomp` stall at the IDENTICAL point (268.9 kN,
  0.0098%), which is the expected result: the two laws differ only in compression and nothing has
  reached compressive yield at that drift. Diagnostics, not capacities.

- **Stage 4, quasi-static, nobond half (4 cells + the Stage 1 control). All converged to the full
  0.30% target with ascending-branch residuals of 0.6-1.3%, so these peaks are real resistance.**

  | comp | tail | peak (kN) | /test | at drift | Gf multiple |
  |---|---|---|---|---|---|
  | `c02` (control) | solved | 943.8 | 0.979 | 0.266% | 1.00 |
  | `lincomp` | solved | 944.2 | 0.980 | 0.292% | 1.05 |
  | `eppcomp` | solved | 968.3 | 1.005 | **0.2997%** | 1.05 |
  | `lincomp` | paper | 1,158.0 | 1.202 | 0.284% | 5.26 |
  | `eppcomp` | paper | 1,158.0 | 1.202 | 0.221% | 5.26 |

- **THE COMPRESSION LAW IS NEARLY IRRELEVANT TO PEAK STRENGTH HERE, and the reason is measurable.**
  Three fundamentally different compression branches — Concrete02's full softening, "linear at E
  forever", and an EPP cap at f_c — span **2.6%** at the solved tail; at the paper tail `lincomp`
  and `eppcomp` differ by **36 N, 0.003%**. The cause: at peak only **2 to 19 of 10,905 concrete
  struts** are past the compressive yield strain. The peak is set by TENSION cracking and the load
  path it leaves, exactly as D55 argued on the compression cube and D78 on this wall.
  What "linear forever" actually costs is visible in the strain field rather than the strength:
  `lincomp` runs a strut to **7.69x the cap (-8.65e-3, i.e. ~215 MPa of compression)**, which
  `eppcomp` caps at 28 — physically absurd, numerically inconsequential.

- **THE TENSION TAIL DOMINATES, and it is a MESH ARTEFACT, not a physical difference.** `paper`
  gives +22.6% base shear over `solved` for BOTH compression laws (944 -> 1,158, 0.98 -> 1.20 of
  the measured). It also dissipates a measured **5.26 x Gf** per cracking strut — the plan predicted
  5.2 — because Table 1's a2/a3 = 70/360 were fitted at his 20 mm grid and the tail's energy scales
  with strut length. So the tail axis moves the answer ~8x more than the compression axis, and it
  buys that with fracture energy the specimen does not have. This is the quantitative form of D75's
  "Gf is not a neutral knob".

- **CAUTION on `eppcomp`/solved: its peak is a LOWER BOUND.** It peaked at 0.2997% against a target
  of 0.3000%, i.e. at the last step, still ascending. Comparing it to cells that turned over is not
  like-for-like, and quoting 1.005 as "the model matches the test" would be the D78 trap. Re-run it
  further before using that number.

- **Bond wired per PLAN §4** (`build.wall_lattice(bond_law="aydin", bond_residual=)`): a bond link
  is a concrete strut with a residual plateau — elastic to f_t, brittle fall to 0.6 f_t, flat,
  symmetric, at the FULL concrete-strut area. The parent's own `--bond` keeps D75's force-slip form,
  so §9's "existing package untouched" holds and D75/D76 stay reproducible. `bond-aNN` parses the
  residual, so `bond-a07` (his 0.7 for deformed bars) is available.

- **The bond artefact, measured for the first time with the anchorage correct (D85):**
  **K_bond/K_perfect = 2.394** (3,033.3 vs 1,267.3 kN/mm), against a no-bond gate of 0.9988. Divided
  by the transformed section instead it is ~2.2, next to D72's 2.19 — but D72's was measured with
  the base steel UNANCHORED, so the agreement is reassurance, not confirmation. **This confounds
  what bond "buys":** at 2.4x elastic stiffness a bonded cell's stiffness and drift-at-peak are
  measuring the ring, not bond behaviour. Strength may survive better, since the law only bites past
  f_t. The plan's escape route — bar density, 73.5% of nodes here carrying a steel node against
  35.3% in his mesh-20 model — predicts a much milder artefact on the replica and is NOT yet
  measured.

- **Every cell still peaks at 0.22-0.30% drift against the test's ~1%**, the perfect-bond signature
  of D75. No cell in the nobond half addresses that, by construction.

- **Status:** accepted. 4 static + 4 quasi-static + 2 elastic runs under
  `examples/output/aydin_aldemir_wall/study/`, matrix and cross-run figures in `master_report.md`.
  One accidental kill mid-batch, its directory marked ABORTED and re-run.

### D88 — 2026-09-06 — Pushed to 2% drift, the Aldemir lattice NEVER TURNS OVER. Every "peak" in D87 was the drift target, and the early-peak finding is withdrawn

- **The run:** `lincomp/solved/nobond` to 2.0% drift, explicit, 740,116 steps, 79 min, converged,
  ascending-branch residual **0.4%** of peak. Prompted by the user asking whether the pushovers
  "ended early, drifted that much, or collapsed" — all five D87 cells had simply stopped at the
  0.30% target they were given, `converged=True`, `end/peak` 0.979-0.998.

- **THERE IS NO CAPACITY IN THIS RANGE.** The maximum, 1,054.6 kN, sits at **1.9994% drift — the
  last step**. Shape, as slope per 0.1% of drift:

  | drift | shear (kN) | /test | slope |
  |---|---|---|---|
  | 0.3% | 941.7 | 0.977 | — |
  | 0.5% | 983.4 | 1.021 | +20.9 |
  | 0.9% | 1,005.4 | 1.043 | +3.2 |
  | 1.1% | 1,008.9 | 1.047 | **+1.8** |
  | 1.5% | 1,027.8 | 1.067 | +5.8 |
  | 2.0% | 1,051.6 | 1.091 | +2.7 |

  A steep climb to ~0.5%, a SHOULDER at 0.9-1.1%, then **re-stiffening**. Not a plateau.

- **WHAT IS WITHDRAWN.** D87 reported that every cell "peaks at 0.22-0.40% drift against the test's
  ~1%" and read that as the perfect-bond signature. It was an artefact of the target: the same model
  passes 944 kN (its 0.30% "peak") and keeps going. The cross-specimen claim built on it — that the
  early peak is systematic across the repo — is half-retracted: SW-NC-FF's 0.36%-vs-1.00% stands
  (measured on a run driven to 4% drift, and its cause is plain-bar debonding), but this wall shows
  nothing of the kind once pushed far enough. **A mid-run reading of "it plateaus" was also wrong**
  and is corrected here; it was taken at 68% and landed exactly on the shoulder.

- **WHAT THE RISING BRANCH ACTUALLY IS.** Nothing in this model can shed load: perfect bond cannot
  slip, `Steel02` has no bar buckling or fracture, and the concrete is cracked through by ~1%. What
  remains is the reinforcement's b = 0.01 hardening. So **1.094 at 2% drift is a statement about the
  steel law, not a strength prediction**, and the only defensible comparison is at a drift the test
  actually reached.

- **THE FAIR COMPARISON, at matched displacement (D77's method):** at 9.00 mm the model reads 975.6
  kN against the test envelope's 969.4 — **1.006**. Past 16 mm the figure's frame clips and there is
  nothing to compare against, which is why the run's own report prints `nan` there rather than a
  ratio. At the test's own peak drift (~0.89%, 20 mm) the model is **1.043**.

- **PROCESS, and this is the third time.** D68 recorded a run that ended at its requested target
  saying nothing about behaviour past it; D78 recorded the same trap one level up; D87 then quoted
  five target-limited peaks as capacities, and this run's own maximum is target-limited too. The
  harness now writes the drift target into every run directory name (`_d2pct`), default or not, and
  the user's standing instruction is that the target is asked for before any run is launched.

- **Status:** accepted. Supersedes D87's peak-drift finding; D87's compression-law and tension-tail
  results are unaffected, since those compared cells at equal targets.

### D89 — 2026-09-06 — The first CYCLIC run: cycling produces the capacity monotonic pushing never does — 918.8 kN, 0.95 of the measured — and it exposed a matched-displacement comparison that is meaningless on cyclic data

- **The run:** `cyclic / c02 / solved / nobond`, the eight-level ladder to 1.0% drift, one cycle per
  level, 4,514,700 steps, 7.8 h, converged, ascending-branch residual **0.7%** of peak. `Concrete02`
  is the carrier because D86's probe showed neither Aydin law can pinch. **The protocol is
  invented** (PLAN §6), so drift capacity and cyclic energy are not fair claims against the test.

- **IT TURNS OVER, AND NOTHING ELSE IN THIS STUDY DOES.** Loop tips, from the reversal records:

  | level | push (kN) | pull (kN) | loop area (kN·m) |
  |---|---|---|---|
  | 0.30% | 867.2 | -817.2 | 2.78 |
  | 0.50% | **918.7** | -841.2 | 5.95 |
  | 0.75% | **918.8** | **-863.0** | 9.81 |
  | 1.00% | 913.4 | -860.2 | 11.56 |

  Both directions stop gaining and edge down across a doubling of drift, while loop area keeps
  growing (34.7 kN·m total) — the wall dissipates more while it stops getting stronger. End
  residual drift **-0.0000%**: no ratcheting. Against the test: **0.953 push / 0.896 pull**, and
  **0.89-0.91 at matched displacement** over 4.5-13.5 mm.

- **THE CYCLIC ENVELOPE RUNS 7-13% BELOW THE MONOTONIC BACKBONE** at equal drift (0.868 / 0.927 /
  0.897 at 0.3 / 0.5 / 0.7%). Prior cycles damage struts a single push never revisits, which is
  exactly the mechanism D88 found missing: monotonic runs here cannot lose load (perfect bond cannot
  slip, `Steel02` has no buckling or fracture), so they rise indefinitely. **This is the first
  number in the study that is a capacity rather than a stopping point.**

- **BUG FOUND AND FIXED — `references.comparison_points` was meaningless on cyclic data.** It took
  the shear at the index nearest a target displacement. A cyclic history revisits every
  displacement: this run passes |u| = 4.5 mm **29,607 times** carrying anywhere from -817 to
  +847 kN, so the "matched displacement" table reported the model at **39.4 kN where the test read
  969.4** — a point on an unloading branch. It now reduces the model the SAME way the test side is
  reduced, as an envelope in a window around the target, which on a monotonic push returns the same
  single point it always did. The affected run's `data.json` and `report.md` were regenerated;
  no earlier conclusion rested on it, since every previous run was monotonic.

- **Against Aydin.** His monotonic simulation overshoots this specimen by **1.208x** — his worst of
  six. Our cyclic envelope is **0.95x**. His model has no cyclic-damage mechanism either (it is a
  monotonic prediction of a cyclic test), and the 7-13% we measure between cyclic and monotonic
  envelopes is roughly half his overshoot. That makes missing cyclic degradation a plausible PARTIAL
  explanation of his error, not a demonstrated one: his tail differs too, and our own `paper`-tail
  cells reproduced 1.20 monotonically (D87).

- **Cross-specimen.** Pull/push settles at 0.93-0.94, the ordinary direction for a wall pushed
  first, matching D79's single-strut prediction for this carrier (0.915). SW-NC-FF's cyclic model
  gave 1.08/1.05 and WSH3 0.950; this wall at 0.953/0.896 is the third in that band, and the first
  here measured the same way those two were.

- **Status:** accepted. Library-adjacent fix in `study/references.py`. Board updated with the exact
  loops (1,817 points kept from 4.5 M, striding plus every reversal).

### D90 — 2026-09-06 — The bond artefact IS the ring's density: halving the mesh cuts the excess stiffness by 56%. D80 item 3's "NOT MEASURED" is closed, and the calibration is confirmed mesh-objective under bond

- **The measurement.** `elastic / lincomp / solved / bond-a06` at mesh 25 against the same cell at
  mesh 50, each run building its own bonded and perfect-bond elastic twins:

  | mesh | nodes carrying a steel node | K_bond/K_perfect | K_perfect | K_continuum |
  |---|---|---|---|---|
  | 50 | 73.5% | **2.394** | 1,267.3 | 1,268.8 |
  | 25 | 43.0% | **1.614** | 1,257.0 | 1,273.0 |

  The EXCESS falls 1.394 -> 0.614, i.e. **56%**, for a 41% fall in density. PLAN §4 predicted this
  ("the density is the variable") and D80 item 3 recorded it as NOT MEASURED. It is now measured.

- **THE CALIBRATION IS MESH-OBJECTIVE, including under bond.** `K_perfect` moves 0.8% and
  `K_continuum` 0.3% across a mesh halving, so the artefact is the ring and nothing else. This is a
  stronger check than D72's, which only ever compared unbonded models.

- **Aydin's own density is 35.5%, below both points**, so his artefact should be milder still —
  consistent with his `K_sim` = 0.91x the measured rather than 2x. **NOT extrapolated onto a number:
  two points do not fix a law, and mesh 20 is ILLEGAL on the Fig. 10(a) panel (2250/20 = 112.5), so
  his density is not reachable here at all.** That comparison belongs to `replica/`, which already
  runs his mesh — the elastic A/B there (`run.py --elastic --bond` against `--elastic`) is the
  measurement that would settle it, and is still unrun (D76).

- **What it changes.** The horizon-ring topology is not fundamentally broken, it is under-resolved
  at mesh 50: a bonded run there is contaminated 139%, at mesh 25 61%. Neither is small enough to
  read strength from directly, so **the perfect-bond results remain the study's usable ones**, but
  the bond axis is no longer a dead end. The unbuilt alternative — a coincident-node `zeroLength`
  interface, which CAN reduce to perfect bond in the limit — is still the principled fix.

- **A bug in the diagnostic that reported this.** `run.py` printed a hardcoded "73.5% of nodes here
  carry a steel node" on the mesh-25 run, stating the mesh-50 value for exactly the variable the
  sentence calls the variable. Now computed from the model's own bond elements. The densities in the
  table above come from geometry, not from that print, so no number here was affected.

- **Status:** accepted. One elastic run (~20 min at mesh 25: 11,011 concrete nodes, three elastic
  solves plus a continuum). Nothing in the nonlinear record changes.

### D91 — 2026-09-07 — Cycling to 2% degrades the wall by 2.9% and the pull side not at all; adding bar fracture + real crushing gives 0.998x the measured peak and a 1.03% DRIFT CAPACITY. The failure mechanisms, not the compression law, are what was missing

- **Two runs, and the second is the control for the first.**

  | run | law | to | peak | retained at 2.0% |
  |---|---|---|---|---|
  | `cyclic c02-solved nobond d2pct ext2p0` | stock | 2.0% | +918.8 / -884.0 kN | **97.1% push, pull still rising** |
  | `pushover c02-solved eps0.05 res0 d2pct` | + fracture + crushing | 2.0% | **961.7 kN** | **46.5%** |

  27.1 h and 2.4 h respectively, both converged, residuals 0.7% and 0.8%.

- **STOCK CYCLING BARELY DEGRADES THIS WALL.** Exact reversal tips: push peaks **918.7 kN at 0.50%**
  drift then 913.3 / 911.7 / 906.5 / 893.7 / 888.9 / **891.9** through 2.00% — a **2.9%** decline,
  flat within scatter. The PULL side never turns over at all: -843.3 -> -883.9 kN, rising
  monotonically to the end, and the push/pull ratio closes from 0.918 to 0.991 as damage saturates.
  D89's "it turns over" reading stands but is now quantified and is much weaker than it looked from
  four levels: **the model cannot fail cyclically either.**

- **WITH FAILURE ENABLED IT FAILS PROPERLY, AND AT THE RIGHT PLACE.** `MinMax` on `Steel02` at
  eps_su = 0.05 plus the crushing floor lowered to zero: peak **961.7 kN = 0.998x the measured
  963.6**, at 0.508% drift; **-20% at 1.033% drift** (the usual capacity criterion), -50% at 1.831%,
  ending at 0.465 of peak. At matched displacement 9.00 mm reads **0.974** of the test envelope.
  The test reached ~1% drift and reported NO degradation, so a capacity just past that is
  self-consistent rather than contradicted.

- **ATTRIBUTION, and it needed the control.** The failure run was first compared against a
  `lincomp` push, whose compression is linear at E forever — the wrong baseline, chosen by
  extrapolating a 2.6% law difference measured at 0.3% drift (D87) to 2%, where it was never
  established. The cyclic run above IS the right control: same law, same drift, switches off.
  **Cycling is the harsher loading, yet it retains 97.1% where the monotonic run with switches
  retains 46.5%** — so the collapse is the failure mechanisms, not c02's compression branch.

- **WHAT THE FAILURE RUN CANNOT SAY.** At its peak, **zero** bars had ruptured and only 7 struts of
  10,905 were on the compression descending branch (61 and 9 by the end). Failure localises, so
  small counts can still halve a wall — but the two mechanisms were enabled together and cannot be
  separated by this run. Two further 80-minute pushes would do it.

- **eps_su = 0.05 IS AN ASSUMPTION.** Table 1 gives only f_y = 360 and elastic-perfectly-plastic.
  The rupture strain, not the wall, sets the 1.03% capacity: welded mesh at 0.025 would fail near
  1%, hot-rolled bar at 0.075 well past 2%. Any drift-capacity claim must carry it.

- **New library capability.** `materials.steel_uniaxial_ruptured` (Steel02 wrapped in MinMax) and a
  `build_lattice_rc` that accepts a rebar factory returning SEVERAL materials, since a wrapper must
  reference its base material's tag; single-material factories are unchanged and
  `build_continuum_rc` raises rather than mishandling one. `build.wall_lattice(concrete_residual=)`
  also had to lower `fcu` ON THE GRADE — `residual_ratio` alone does nothing here, because the
  material takes `max(grade.fcu, ratio*fc)` and this specimen's grade already carries 0.2*fc.
  Caught before the run; otherwise it would have "tested crushing" while changing nothing.

- **This is the mechanism CLAUDE.md records as missing for SW-NC-FF** (no post-peak degradation,
  1.34x at 4% drift, cause logged as `Steel02` having no bar buckling). The same two switches should
  apply there, and WSH3 — which ended by bar fracture at 1.79% — is where the MinMax route is most
  directly justified. **Drift capacity moves from "outside the model by construction" to a
  measurable quantity conditioned on one assumed material property.**

- **Status:** accepted. Supersedes D88's "the model has no capacity" for the failure-enabled
  configuration; D88 and D89 remain correct for the stock one. Tests 55 pass, 1 known D34 failure.

### D92 — 2026-09-07 — `peak_shear` was a SAMPLE MAXIMUM, and on a run with failure switches it reports the ringing of a released bar rather than resistance. Smoothed over 1 ms the sweep is monotone and every failure cell reads 0.997–0.998x the measured peak

**The trigger.** The eps_su = 0.025 cell of the failure sweep (D91) reported 996.6 kN = **1.034x**
the measured 963.6, breaking a pattern in which every other failure cell read 0.998x. A cell with
LESS bar ductility returning MORE strength has no mechanism behind it, so the number was suspect
before it was explained.

**It is an artefact of the metric, not of the model.** Three things establish that, all from the
stored series — no re-run:

- The eps_su = 0.025 and no-rupture curves are **identical to four decimals up to 0.45% drift** and
  diverge at 0.46–0.47%, exactly where the first bar breaks. Nothing upstream of rupture differs.
- At the reported peak the shear **oscillates +-56 kN with a ~1.83 ms period**, against T1 = 5.82 ms.
  That is a LOCAL mode — the energy the ruptured bar was holding, ringing in a few elements — not a
  structural response.
- A **1 ms moving average** gives **961.0 kN at 0.4692% drift = 0.997x** the measured, back in family.
  The smoothed value is insensitive to the window: 966.8 / 961.0 / 958.1 / 956.2 kN over 0.5 / 1 / 2 /
  5 ms, a 1.1% spread across a tenfold range, while the raw maximum reads 996.6.

`shear` is sampled EVERY step (1.1 M samples at dt = 8 us), so `max()` over it is a maximum over
samples, and an explicit march resolves modes far above the one being loaded. The 1 ms window is
three decades clear of both ends: long against the 1.83 ms local mode, short against the ~1.5 s the
7.6 mm/s drive takes to traverse the ascending branch.

**Rescoring every finished run found a SECOND ringing cell that had never been reported:**
`2026-09-05_223007_pushover_lincomp-solved_bond-a06` reads 1,373.4 kN raw against **1,286.0 kN**
smoothed — **1.068x**, the largest correction in the study. The bonded peak was 6.8% ringing.
All other cells move under 0.7%.

**What the correction does NOT change: the drift capacities of D91 stand exactly.** Measured as the
80%-of-peak drop on the smoothed curve, they move 0.5911 -> 0.5919%, 1.1549 -> 1.1550%, 1.3739 ->
1.3739%. The peak VALUE was wrong; the capacity ordering, and D91's conclusion that capacity is
eps_su-dominated below 0.05, are untouched.

**A second reporting defect fell out, and it cuts the other way.** On the cells NOT truncated by
rupture the peak is a **PLATEAU**: the no-rupture curve holds within 0.5% of its maximum from 0.475%
to 0.703% drift. "Peak at 0.508%" and "peak at 0.647%" are the same measurement — an argmax over a
flat top — and the sweep table was quoting them at four significant figures. The eps_su = 0.025 cell
by contrast peaks inside 0.0002 drift-percentage-points, because rupture truncates the plateau, so
there the drift IS resolved. `peak_plateau` now reports the range and the report says which case a
cell is in.

**Built.** `study/metrics.py` — an edge-correct centred moving average (each output divides by the
samples actually in the window, so the ends are not biased towards zero), `response_metrics`
returning the smoothed peak, the raw maximum, `peak_ringing_ratio`, the plateau and the 80%-drop
`drift_capacity`, and `headline_peak`/`ringing_note` so every consumer quotes the same number.
`run.py`, `report.py` and `master.py` now lead with the smoothed peak and print the raw one beside
it whenever the two differ by more than 2%; `master.py` marks with `*` any cell still quoting an
uncorrected maximum. Both dynamic runners now return **`dt`**, the step they actually marched — the
series were sampled per step and there was no recorded way to turn a sample index into a time.

**`study/rescore.py` recovers this from records that already exist.** Every path-following runner
stores the whole `shear` series, so 2 h of compute per cell is re-scored in seconds: `data.json` is
rewritten through a temp file and an atomic rename (an interrupted rescore cannot cost a run record),
only keys are ADDED so `peak_shear` keeps its recorded value, and `report.md` is regenerated so a
report never outlives the numbers it quotes. Applied to all 17 runs carrying a response.

**The general lesson, and where it lands next.** A maximum over a finely-sampled series is a
measurement of the sampling as much as of the structure, and it gets worse exactly where the physics
gets interesting — every future run with `--steel-rupture` or `--concrete-residual` is exposed, and
so is the WSH3 validation these switches were built for. Cell 4 of the sweep (eps_su = 0.075) was
already running when this was found and will need `rescore.py`; its process holds the pre-D92
`run.py`.

**Status:** accepted. Corrects the reported peak of D91's eps_su = 0.025 cell (996.6 -> 961.0 kN,
1.034 -> 0.997x) and of the D75 bonded cell (1,373.4 -> 1,286.0 kN); supersedes neither entry's
conclusions, both of which rest on capacity and on ratios the correction leaves standing.

### D93 — 2026-09-07 — The failure switches are wired into WSH3, where the rupture strain is MEASURED rather than assumed — the one specimen in the repo that can VALIDATE the D91 mechanism instead of being fitted by it

**Why WSH3 and not another wall.** D91 gave the Aldemir wall a drift capacity by adding bar rupture
and real crushing, and the capacity it predicts is set by `eps_su`, which that source never prints:
the sweep runs 0.591% at 0.025 to 1.374% at 0.075 (D91/D92). A parameter that free, deciding the
headline answer, is a fit. WSH3 breaks that circle from both ends:

- **The rupture strain is data.** Dazio, Beyer & Bachmann's Table 2 gives A_gt per bar type, and
  `specimen.py` already read those numbers — they set the hardening slope `b`. They are now recorded
  in `specimen.AGT` AS THE GRADES ARE BUILT, so the number that hardens a bar and the number that
  breaks it are literally the same one and cannot drift apart: phi12 **7.69%**, phi8 7.34%, phi6
  6.45%, phi4.2 3.06%. `--steel-rupture measured` gives each bar its own.
- **The failure drift is published.** The corner bar ruptured at **1.79% drift**, ending the test.
  So the model can be checked against a failure it was not shown, which is what the Aldemir wall
  cannot offer.

**A_gt IS NOT THE FRACTURE STRAIN**, and the difference is recorded rather than smoothed over. A_gt
is uniform elongation, where necking begins; a bar fractures past it. Using A_gt removes the bar at
its peak force, which is CONSERVATIVE — early.

**The prediction, stated before the run.** Expect the model to rupture LATE. The test's corner bar
broke after buckling from 1.70% drift, and a bar that has buckled and straightened fractures at a
far lower nominal tensile strain than its coupon A_gt. Nothing in a plane lattice of truss struts
buckles. So a predicted capacity at or past 1.79% is the honest expectation and would confirm the
mechanism without confirming the timing; an EARLY rupture would be the surprising result and would
mean something other than bar ductility is ending the run. Writing this down first is the D75
process lesson.

**Wired.** `build.nonlinear_wall_lattice(steel_rupture=, concrete_residual=)`, with
`--steel-rupture {measured|<strain>}` and `--concrete-residual` on both `pushover.py` and
`cyclic.py`. Two things carried over from D91 rather than rediscovered: the crushing floor must be
lowered ON THE GRADE (`max(grade.fcu, ratio*fc)` defeats the argument alone, and these grades carry
`fcu = 0.2*fc`), and `compression="elastic"` now REJECTS a lowered residual instead of silently
ignoring it, since that mode holds the plateau at fc and has no crushing to floor. Verified on a
build: `fpcu` 7.84 -> 0.00 MPa, three `MinMax` wrappers at 0.0645 / 0.0734 / 0.0769, and a 30-second
dynamic pushover completes with the switches live.

**Provenance gap closed at the same time.** Neither script recorded `gf_factor` in its `data.json`,
though D87 measured Gf x2 at **+14.2%** base shear — so two WSH3 runs differing by the single most
consequential knob in the study produced records that could not be told apart. All three of
`gf_factor`, `steel_rupture` and `concrete_residual` are now saved by both scripts.

**PROCESS FAILURE, recorded because the rule it broke is written down.** The 30-second verification
above was launched against `examples/output/katrin_wall/`'s FIXED STEM, whose own `runs/README.md`
says: *"Never launch a run — not even a two-minute diagnostic — while an unarchived result is
sitting in the parent directory."* It also skipped the standing requirement to agree a drift target
before launching. Damage assessment: the four `wsh3_pushover_dynamic*` files were overwritten;
nothing else was touched, the 19-hour archived cyclic run in `runs/2026-08-25_cyclic_1p02_gf2/` is
intact, and the WSH3 report cites only cyclic figures, so no result any document depends on was
lost. **The structural fix is the one D68/D70 already name** — `katrin_wall` and `wall` should give
every run its own timestamped directory the way `vk3_wall` and `compression_cube` do, so the rule
stops being something to remember. Until that lands the fixed stem is a live hazard.

**Status:** accepted, code only — NO WSH3 validation run has been launched. Staged:
`cyclic.py --drift 0.0179 --gf-factor 2 --steel-rupture measured --concrete-residual 0.0`, whose
cost must be priced first (D66 measured 276 ms/step, and 2.03% drift was estimated at ~40 h before
these switches).

### D95 — 2026-09-08 — The bond law was ELASTIC: a broken bond healed on every reversal. Aydin's own word is "brittle", which means irreversible, so the faithful reading is a DAMAGING law — same envelope to four decimals, and it does not heal

**Found by stopping a run.** The 93-hour bonded cyclic run of 2026-09-08 was aborted at 2.4% (step
360,000 of 15,133,338) once it was noticed that its 21,244 bond links — **61% of all elements** —
use `ElasticMultiLinear`, which is PATH-INDEPENDENT. The record is kept at
`runs/2026-09-08_083416_cyclic_..._protoladder/` with an `ABORTED.md` saying why.

**The same objection had already retired both concrete laws from cyclic work** (D60 on `aydin_cube`:
"NOT for cyclic work"; D83/D86 choosing `Concrete02` over `HystereticSM` as the carrier). It was
never carried across to the bond law, and `carrier_probe.py`'s own docstring states it for the
concrete laws while the bond law sat in the same model untouched.

**The source supports the change rather than merely permitting it.** Aydin, Tuncay & Binici (2019)
say: *"For the bond, an elastic brittle response with residual strength was used as shown in
Fig. 1(c). The residual bond strength parameter (a) was chosen as 0.7 ... to ensure that the
deformed bar residual bond strength was reflected accurately."* **Brittle means irreversible.** A
law that retraces its own drop on unloading is not brittle; it is nonlinear-elastic. So the
damaging form is the more faithful reading of the paper, not a departure from it.

**Built as a PARALLEL pair** (`materials.bond_elastic_brittle_damaging`), which gives the envelope
exactly instead of approximately:

  * `Elastic` at `(1-a)E` wrapped in `MinMax` at +/- eps_cr — the part of the rise that DISAPPEARS
    at cracking, and `MinMax` never comes back;
  * `ElasticPP` at `a*E` yielding at eps_cr — the part that SURVIVES, holding `a*ft` forever.

Below eps_cr the two sum to `E` exactly and reach `ft` exactly at eps_cr; above it the first is dead
and the second holds the plateau. The `ElasticMultiLinear` form cannot fall vertically and fakes the
drop over a `drop=1e-3` strain band; **this one has an exact corner and needs no ramp**.

**MEASURED, on one link, with `carrier_probe`'s own harness (0 -> +3 eps_cr -> 0 -> +3 eps_cr,
sampled at +0.5 eps_cr):**

| law | virgin | reload after damage | ratio |
|---|---|---|---|
| `ElasticMultiLinear` (what ran) | 0.930 MPa | **0.930 MPa** | **1.00 — fully healed** |
| `Parallel` damaging (new) | 0.930 MPa | **-0.552 MPa** | -0.59 — cannot heal |

**And the monotonic envelopes are IDENTICAL** — 0.4625 / 0.9287 / 1.3875 / 1.8315 / 1.1100 /
1.1100 / 1.1100 at 0.25 / 0.5 / 0.75 / 0.99 / 1.01 / 2.0 / 3.9 eps_cr, `diff = 0.0000` at every
point including across the drop. So every monotonic bonded result already obtained STANDS, and the
two laws are directly comparable: they differ only where the elastic one was never entitled to an
answer.

**On reversal** the residual branch unloads elastically at `a*E` and yields at `-a*ft`, i.e. the
surviving bond is elastic-perfectly-plastic both ways. That is the right shape for what the residual
physically IS — rib-on-concrete friction after adhesion has gone — and it **dissipates**, where the
elastic form dissipates nothing. A cyclic bonded run on the old law would have reported hysteretic
energy that its own bond elements could not produce.

**Library change beyond the material:** the builder's bond path assumed ONE `UniaxialMaterial`
(`bmat.id = tag`). It now accepts a sequence and rebases the stack's internal tag references,
elements using the LAST — the same convention `build_lattice_rc` already applied to a wrapped rebar
material (D91). Single-material bond factories are unchanged.

**Status:** accepted, material and builder only. NOT yet exposed on the study's `bond` axis, and no
run has used it — the axis currently offers `perfect` / `bond60` / `bond70` (renamed from
`nobond` / `bond-a06` / `bond-a07` in D94) and the naming of a damaging variant is an open choice.
Tests 55 pass, 1 known D34 failure.

### D96 — 2026-09-08 — The bond residual `a` has TWO values, both from Aydin, and they disagree: the published 0.7 and the author's own 0.6 given to this project. Default switched to 0.7 by user instruction; the difference is now something to measure, not to settle by choosing

**The conflict**, surfaced while preparing to re-run bonded work on the damaging law (D95):

| `a` | source | date |
|---|---|---|
| **0.7** | PUBLISHED — Aydin, Tuncay & Binici (2019), "chosen based on the preliminary simulation results of Aydin (2017) to ensure that the deformed bar residual bond strength was reflected accurately" | 2019 |
| **0.6** | THE AUTHOR, direct to this project (PLAN §4): "a bond link is a concrete strut that drops to ~60% of its capacity and stays there" | 2026-09-05 |

The second is LATER and SPECIFIC to this replication; the first is the one a reader of the paper
can check. Neither is our invention, and there is no basis in the repo for calling either wrong.

**Changed:** `build.wall_lattice(bond_residual=)` default 0.6 -> **0.7**, on user instruction
2026-09-08. `materials.bond_elastic_brittle` and `..._damaging` already defaulted to 0.7.

**NOT changed:** the study's `bond` axis still offers `bond60` and `bond70`, so this is a
one-flag comparison. **Every bonded run before 2026-09-08 used 0.6** — including D75's replica work,
D90's ring-density measurements and the two bonded pushovers of 2026-09-08 (peak 1,452 / 1,321 kN,
sustained 780 -> 680 kN) — so a `bond70` result is NOT directly comparable to any of them until the
axis is measured. The cost of that measurement is one bonded pushover.

**Process note.** The switch was requested on a summary of mine that called 0.7 "the paper's value,
not the 0.6 every bonded run here has used", which implied 0.6 was a choice of ours. It was not;
`models.py` and PLAN §4 both attribute it to the author. The provenance of both values is now
recorded at the point of use in `build.py` so the fork cannot be lost again by a one-line summary.

**Status:** accepted. The sensitivity of the wall to `a` between 0.6 and 0.7 is UNMEASURED.

### D97 — 2026-09-09 — Four measurements close the bond and capacity questions: the damaging law reproduces the elastic one to 0.49% on 21,244 elements, `a` is a near-proportional lever (1.166x), damping is immaterial to strength (0.36% over a 10x change), and DRIFT CAPACITY DOES NOT EXIST INSIDE AYDIN'S CONSTITUTIVE FRAMEWORK

Four runs, launched 2026-09-08, all converged, all landed by 02:11 on 2026-09-09.

**1. The damaging bond law (D95) is verified at assembly scale.** Two bonded pushovers to 1.5%,
identical but for the bond law, 10.2 h each. Across 120 points from 0.10% to 1.49% drift the
damaging law differs from the reversible one by **mean +0.49%, rms 0.51%, worst +0.84%** — under 1%
everywhere, which is what identical envelopes require. D95 proved the envelopes match on ONE link;
this shows the `Parallel(MinMax(Elastic), ElasticPP)` stack marches correctly under an explicit
integrator on 21,244 of them. The one place they part is the start-up spike, **1,775 vs 1,633 kN**,
where the damaging law's exact corner replaces a `drop = 1e-3` strain ramp. **The damaging law is
now licensed for cyclic work**, which is what the aborted 93-hour run lacked.

**2. `a` is a near-proportional lever, and it revises D96's framing.** At matched drift over
0.10-1.49%, **a = 0.7 carries 1.166 x the shear of a = 0.6** (min 1.159, max 1.176) — mean 836 vs
717 kN. A 17% change in residual bond buys a 17% change in wall strength, with no threshold. So
"bond makes this wall weaker" was really "a = 0.6 makes this wall weaker": against the measured
963.6 kN the bonded model averages **0.744 x at a = 0.6 and 0.868 x at his published 0.7**, and
reaches **0.956 x at its sustained peak**. Bond still costs accuracy against perfect bond's 0.997,
but 9 points rather than 22.

**3. Damping is immaterial to strength; our 50% is safe.** Three identical cells (crushing,
eps_su 0.05, residual 0, 1.5%) at zeta = 0.05 / 0.20 / 0.50 give peaks of **964.5 / 961.0 / 961.0
kN = 1.001 / 0.997 / 0.997** — a **0.36% spread across a tenfold change**, so the gap between our
50% and Aydin's published 5% does not matter where it counts. The residual falls monotonically with
damping (25.9 / 15.3 / 7.7 kN = 2.7% / 1.6% / 0.8% of peak), reproducing D64's finding that damping
suppresses crack-release ringing. **CAPACITY IS THE SOFTER NUMBER**: 1.113 / 1.152 / 1.033%, an 11%
spread and non-monotonic, so **every drift capacity in this study carries about +/-6% scatter** from
which struts crack in which step. The sweep's 0.591 / 1.033 / 1.155 / 1.374 ordering survives, but
its precision does not.

**4. DRIFT CAPACITY IS NOT AVAILABLE INSIDE AYDIN'S FRAMEWORK — the headline.** His compression law
is elastic and never crushes (2019: *"Concrete in compression is assumed to be elastic"*). Run with
that law PLUS bar rupture at eps_su = 0.05, the wall traces to 3.0% drift and **never falls to 80%
of peak**. Stronger: **the rupture switch never engages.** Against the same model with rupture
disabled the two curves are identical to **0.000% — mean, rms AND worst — over 150 points to 2%
drift**. Without crushing the compression zone never degrades, the neutral axis does not migrate,
and no bar reaches 5% strain.

So D92/D93's "crushing sets the capacity ceiling at 1.155%" understates it: **crushing is the only
reason this model has a capacity at all.** Every drift-capacity number in this study belongs to OUR
model, not to a replication of his, and must be reported that way. Peak strength is unaffected —
that is a tension-cracking result (D87) and it survives every compression law tried.

**CAVEAT on run 4, flagged rather than buried:** the smoothed curve peaks at 1,079.3 kN at 2.989%
and reads 958 kN at 3.000% — an 11% fall in the last ~4,000 steps. It never met the 80% criterion so
the result stands, but whether it was about to turn over is exactly the D75 process trap and only a
longer run settles it.

**Status:** accepted. Supersedes nothing; extends D92/D93 (capacity), D95 (bond law) and D96 (`a`).
Tests 55 pass, 1 known D34 failure.

---

### D98 — 2026-09-10 — Halving the mesh moves drift capacity 7.3%, inside the +/-6% damping scatter: the capacity ceiling is a property of the MODEL, not of the discretization. Peak strength does carry a small real mesh dependence (+5.6%). And the master matrix hid a run again, one axis down

**Context.** D97 left every drift-capacity number in this study resting on one untested support: they
were all measured at mesh 50. Crack-band regularization (D20) guarantees that the *dissipation* of a
single softening strut is mesh-independent; it guarantees nothing about a **load-path collapse**,
which is what sets capacity here (D55/D78/D87) and which depends on how many struts exist to
redistribute among. So the question could not be answered by argument.

**The run.** The `crushing / solved / perfect` cell with both failure switches (`eps_su = 0.05`,
`concrete_residual = 0`) to 1.5% drift, rebuilt at **mesh 25** — 13,531 -> **48,662 elements**,
2,806 -> **11,011 nodes** (3.60x the elements). Everything else identical, including zeta = 0.50 and
the 7.6 mm/s drive. Explicit CentralDifference sized from `critical_time_step` as always: the stable
step halves with the strut length, so 728 -> **1,464 steps/period** and 555,087 -> **1,113,594
steps**. Cost **10.7 h against 2.3 h** (4.70x, i.e. 3.60x the elements x 2.01x the steps, less
super-linear overhead than feared).

| | mesh 50 | mesh 25 | ratio |
|---|---|---|---|
| peak (1 ms smoothed) | 961.0 kN | 1,015.1 kN | **1.056** |
| /measured 963.6 kN | 0.997 | **1.053** | |
| drift at peak | 0.646% | 0.770% | |
| **80%-drop drift capacity** | **1.0328%** | **1.1084%** | **1.073** |
| ascending-branch residual p95 | 7.73 kN (0.8%) | 6.47 kN (0.6%) | |

**1. Capacity is mesh-objective at the precision this study can claim.** The 7.3% shift is the same
size as the **+/-6% scatter D97 measured from damping alone** (1.113 / 1.152 / 1.033% at zeta =
0.05 / 0.20 / 0.50). Two independent perturbations that touch nothing constitutive move capacity by
about the same amount, which is the signature of a quantity set by *which struts crack in which
step* rather than by the discretization. **The sweep's 0.591 / 1.033 / 1.155 / 1.374% ordering
therefore stands** — those separations are 1.8x to 2.3x apart and survive a 7% wobble — but no
capacity in this study should be quoted to better than about +/-7%.

**2. Peak strength DOES have a small real mesh dependence, and it is not scatter.** +5.6%, against
the ~+/-3% run-to-run scatter established for peak shear. At matched drift over 0.30-1.49% the
mesh-25 curve sits at **mean 1.0515** of the mesh-50 curve (min 0.916, max 1.149 — the extremes are
reversal/crack-release excursions, the mean is the signal). The finer grid is slightly stronger:
more struts crossing the same section, and a finer crack band releasing in smaller increments. Worth
one sentence in the write-up; it changes no conclusion, since **0.997 and 1.053 are both far closer
to the measurement than Aydin's own 1.208x** on this specimen.

**3. It also cuts the other way on his replication.** Aydin ran at **20 mm** — finer than either of
ours — where our own trend says the model gets *stronger*. That weakens "mesh" as an explanation for
his +21% overshoot and strengthens the constitutive reading already indicated by D87 and D97: an
elastic compression law with no crushing, plus his bond parameters, is where his direction of error
comes from.

**4. Reporting defect, second instance.** Regenerating `master_report.md` put the mesh-25 run in the
`crushing/solved/perfect` cell as the newest, so the matrix advertised **1.053x where the cell's
headline is 0.997x**. Same class as D92's shadowed cyclic runs, one axis down: the matrix keys on
comp x tail x bond x analysis, so *any* other parameter — mesh, eps_su, residual, damping — silently
substitutes. Fixed by `master.variant_note()`: a row now names what is non-default about the run it
is showing (`pushover (12, newest) — mesh 25, steel_rupture 0.05, concrete_residual 0`), rendered as
an amber chip on the advisor page. **The general lesson, now twice learned: a matrix that reduces N
runs to one cell must always say which run it picked and how it differs.**

**Status:** accepted. Extends D92/D93/D97 (capacity) and D20 (crack-band regularization: dissipation
is regularized, collapse drift is not — it is merely insensitive). Tests 55 pass, 1 known D34
failure.

---

### D99 — 2026-09-10 — The monotonic collapse does NOT survive cycling: at 1.5% drift the cyclic wall carries 1.50x what the push carries. The 80%-drop capacity metric is meaningless on a cyclic trace and reported a number below the drift at peak. And the reason nothing degrades is that `epsU` is not a failure strain

**The run.** `crushing / solved / perfect`, both failure switches (`eps_su = 0.05`,
`concrete_residual = 0`), ladder to 1.5%, 8 levels x 1 cycle, 7.6 mm/s, zeta = 0.5. Converged,
6,772,049 steps, 11.9 h. Built to be the CONTROLLED TWIN of the finished monotonic push
(`2026-09-08_215009`): same 2,806 nodes / 13,531 elements, same dt_crit 10.006 us, same 728
steps/period, same target. The only variable is monotonic against cyclic.

**1. THE PUSH DEGRADES; THE CYCLING DOES NOT.** Every loop tip, 1 ms smoothed, against the
monotonic twin at the same drift:

| level | tip + | tip - | monotonic | cyc/mono | +tip/test |
|---|---|---|---|---|---|
| +/-0.075% | 794.8 | 781.0 | 794.8 | 1.00 | 0.825 |
| +/-0.150% | 821.5 | 827.6 | 839.0 | 0.98 | 0.853 |
| +/-0.225% | 880.1 | 845.3 | 907.0 | 0.97 | 0.913 |
| +/-0.300% | 868.0 | 845.0 | 930.0 | 0.93 | 0.901 |
| +/-0.450% | 924.2 | 837.3 | 951.1 | 0.97 | 0.959 |
| +/-0.750% | **930.4** | 864.0 | 933.1 | 1.00 | **0.966** |
| +/-1.125% | 925.8 | 882.1 | 740.6 | **1.25** | 0.961 |
| +/-1.500% | 915.4 | 893.6 | 609.2 | **1.50** | 0.950 |

(kN.) The push falls **36% below its peak by 1.5%**; the cyclic envelope is FLAT from 0.45% on —
924 / 930 / 926 / 915, 1.6% off its maximum at the end. Peak **931.1 kN = 0.966x** the measured
963.6 against the push's 0.997x, so cycling costs ~3% of STRENGTH and, measurably, **nothing of
capacity**. Ascending-branch residual 0.6% of peak, so the peak is real resistance.

So **the monotonic collapse is DIRECTIONAL** — a load path that must accumulate in one direction to
fail, which every reversal re-shuffles — not strength the wall has lost. The pull side never
degrades at all (781, 828, 845, 845, 837, 864, 882, 894, still rising at the last level): the first
excursion of each level is positive, so the negative half always follows fresh damage and never gets
to be the side that accumulates. **This lands on D97: drift capacity is not just "ours not his", it
is a property of the LOADING PATH, and the 0.591/1.033/1.155/1.374% sweep is a MONOTONIC sweep.**

**2. `drift_capacity` IS INVALID ON A CYCLIC RUN** and the report printed 0.6978%. The 80%-drop
metric was written for a push; a cyclic trace crosses zero shear at every reversal, so it fires on
an unloading branch. The tell is that it landed BELOW the drift at peak (0.744%), which is
impossible for a capacity. On the tip envelope 80% of 930.4 is 744 kN and no later tip approaches
it: capacity here is **> 1.5%, i.e. not reached**. Same class as D92 — a metric reported outside the
domain it was written for. `metrics.py` must reduce a cyclic series to its TIP ENVELOPE first.

**3. WHY NOTHING DEGRADES — `epsU` is not a failure strain.** Concrete02's compression branch is
parabolic to (`fpc`, `epsc0`), LINEAR to (`fpcu`, `epsU`), then **constant at `fpcu` forever**.
Verified through the real code path for this grade (fc 28, epsc0 0.002252, Gf 0.075, Gfc = 250 Gf):

| `concrete_residual` | strut | `fpcu` | `epsU` / `epsc0` |
|---|---|---|---|
| 0.2 (default) | L = 50 | -5.600 | 10.91 |
| 0.2 | L = 70.71 | -5.600 | **8.01** |
| 0.0 | L = 50 | **-0.000** | 12.90 |
| 0.0 | L = 70.71 | **-0.000** | 9.41 |

Three compounding reasons: (a) at the default residual every "crushed" strut carries 5.6 MPa
FOREVER; (b) `epsU` is 9-13x `epsc0`, i.e. 2.1-2.9% compressive strain; (c) the descent is linear
across that whole range, so a strut at 2x `epsc0` still holds 91.6% of fc and at 5x still holds 66%.
D87's count — only 2-19 of 10,905 struts ever pass `epsc0` — completes it. Removing the residual
also LENGTHENS the branch, since `f_cu` sits in the denominator of `2*Gfc/((fc+fcu)*L)`.

**4. A LATENT TRAP found while verifying (3): the compression clamp binds first on DIAGONALS.**
`epsU = max(epsc0 + 2Gfc/((fc+fcu)L), grade.epsU)`, and at the DEFAULT residual the diagonal comes
out at **8.01x `epsc0` against the grade floor of 8.0** — 0.1% of margin. So in every
default-residual run the diagonals' compression crack-band regularization was riding the clamp
rather than being set by strut length, and at mesh 100 it binds outright (4.5x against the 8.0
floor). It affects no current result — the failure runs are all at residual 0, where diagonals sit
at 9.41x — and tension `Ets` is untouched either way. Record it before a coarser mesh is tried.

**5. NEW PARAMETER `steel_b` (schema v5).** Hardening was hardcoded at `b = 0.01`, an UNPRINTED
convention exactly like `eps_su`: the paper gives f_y = 360 and nothing else. It is not negligible —
a bar at `eps_su = 0.05` carries **1.27x f_y** — and it is one of the few things that could hold a
flat cyclic envelope up. `--steel-b 0.0` gives elastic-perfectly-plastic ties. Default reproduces
every earlier run byte-identically (`models.steel_b` returns None when it matches the specimen, so
grade names and material caching are unchanged; verified: same 13,531 elements, Steel02 args
identical at 0.01). A `b = 0` cyclic twin launched 23:20.

**Status:** accepted. Extends D92 (metrics reported outside their domain), D97 (capacity is ours,
and now also path-dependent) and D22/D91 (the residual floor). Tests 55 pass, 1 known D34 failure.

---

### D100 — 2026-09-10 — Read from the 2019 paper: Aydin calibrates exactly TWO things, elastic stiffness in closed form and tension against Cornelissen. THERE IS NO COMPRESSION CALIBRATION ANYWHERE. He also states he never calibrated bond per horizon — and his "Gf is the least important parameter" does not contradict D75, because he RE-FITS the tension tail whenever Gf moves

**Why this was checked.** The study keeps hitting parameters the paper never prints — `eps_su`,
`b`, `Gfc/Gf`, `residual_ratio`, every bond number. Before inventing more, the question was whether
he calibrates anything we have been guessing at. Text extracted from the PDF in
`examples/aydin_aldemir_wall/`, quotations verbatim.

**What he calibrates — two things, and the abstract says so in one line:**
> *"The force-deformation response of each element is calibrated from direct tension tests."*

1. **Elastic, closed form, no fitting.** The energy balance gives `Et*At = 0.621*Et*d*w` at horizon
   1.5d and `0.102` at 3.01d (his Appendix). `calibration.aydin_closed_form_C` reproduces both to
   0.05% (D72).
2. **Tension, an iterative fit done ONCE.** His Fig. 2 is a workflow: `a1, a2, a3` are adjusted
   until the lattice's uniaxial stress-average-displacement response matches the **Cornelissen
   et al. (1986)** stress-displacement model, iterating *"until the energy error is reduced below a
   threshold value"*. Then: *"simulations of the test specimens were conducted as blind predictions
   with the calibrated parameters."*

**COMPRESSION: there is nothing to calibrate, by construction.**
> *"Concrete in compression is assumed to be elastic. However, because of the mesoscale nature of
> modeling, splitting cracks and cover spalling due to compression were simulated as indirect
> tensile failures followed by lattice instability (Kendall 1978; Bazant et al. 1993; Bazant and
> Xiang 1997)."*

The 2017 thesis is the same and more explicit (Sec. 2.2, p. 30: out of scope). **A TRAP:** the paper
does say *"beyond about 30% of its compressive strength"* — that sentence is about **Poisson's ratio**
becoming variable after cracking (Kupfer et al. 1969), used to argue elastic-nu matching is not a
prerequisite. It is NOT a compression constitutive knee and it reads like one.

The only compression calibration in the trilogy is the **2021 cube paper**, and it is GEOMETRIC, not
constitutive: `Rmax/d`, searched by repeated FE runs until the mean strength lands within 10% of fc,
5 realizations per answer (D60). D61 already measured why it cannot be transplanted — perturbation
only ever LOWERS strength, so no `Rmax/d` reaches fc on a Concrete02 lattice. His calibration
presupposes his constitutive law.

**BOND: he says outright that he did not calibrate it per horizon.**
> *"This was remedied by reducing bond strength for the 3.01d horizon; however, no attempt was made
> to calibrate the bond parameter separately for different horizons."*

So the bond numbers D75/D80 could not find were never pinned by him either. Stop looking.

**THE Gf CONTRADICTION IS RESOLVED, and it is a real trap.** His Fig. 12 sensitivity concludes Gf is
the **LEAST important** parameter, while D75 measured `--gf-factor 2` at **+14.2% base shear**. Both
are right, because of one sentence of his:
> *"It should be noted that, upon changing Gf , the softening parameters of the tension model
> [Fig. 1(c)] were calibrated separately for each simulation."*

He RE-FITS `a1, a2, a3` every time Gf moves, so the change is absorbed by the refit and the tail
keeps matching Cornelissen. We change Gf and let the tail move with it. **Different experiments — his
"Gf doesn't matter" cannot be cited against D75, and D67's WSH3 `--gf-factor 2` is not licensed by
it either.** His separate +/-10% CoV study on `a1, a2, a3, b1, b2` moves demands <= 2%, with
`a1, a3, b1` dominating; that is about the FITTED parameters, not about Gf.

**Status:** accepted. Constrains D75/D80 (stop hunting his bond parameters), reframes D67/D75 (the
Gf sensitivity comparison was never like-for-like) and confirms D97 (no capacity exists inside his
framework because no compression law does).

---

### D101 — 2026-09-11 — Strain hardening is the largest single control on drift capacity in this study, and it is a number the paper never prints: `b = 0` gives 0.50% where `b = 0.01` gives > 1.5%. Peak strength barely moves. And the cyclic capacity metric was wrong a SECOND time, so it is now read off the tip envelope

**The experiment.** D99 left the flat cyclic envelope with two candidate causes — strain hardening
in the ties, or the directional load-path argument. `--steel-b 0.0` (D99 §5) isolates the first:
the exact twin of the 1.5% ladder, elastic-perfectly-plastic ties, nothing else changed. Converged,
11.4 h, ascending-branch residual 0.6% of peak.

**1. THE WALL IS DESTROYED.** Loop tips, 1 ms smoothed, against the `b = 0.01` twin:

| level | b = 0 | b = 0.01 | ratio | b=0 / test |
|---|---|---|---|---|
| +/-0.075% | 807.3 / 776.6 | 794.8 / 781.0 | 1.02 / 0.99 | 0.838 |
| +/-0.150% | 808.4 / 781.2 | 821.5 / 827.6 | 0.98 / 0.94 | 0.839 |
| +/-0.225% | **852.7** / 779.2 | 880.1 / 845.3 | 0.97 / 0.92 | **0.885** |
| +/-0.300% | 843.6 / 805.4 | 868.0 / 845.0 | 0.97 / 0.95 | 0.875 |
| +/-0.450% | 790.9 / 759.3 | 924.2 / 837.3 | 0.86 / 0.91 | 0.821 |
| +/-0.750% | **327.2** / 321.9 | 930.4 / 864.0 | **0.35 / 0.37** | 0.340 |
| +/-1.125% | 120.2 / 48.8 | 925.8 / 882.1 | 0.13 / 0.06 | 0.125 |
| +/-1.500% | 48.8 / 55.6 | 915.4 / 893.6 | 0.05 / 0.06 | 0.051 |

(kN.) Peak **853.8 kN = 0.886x** the measured 963.6 against the twin's 0.966, and it arrives at
**0.210% drift instead of 0.744%** — a sharp peak, not a plateau. The envelope then falls off a
cliff in ONE level, 790.9 -> 327.2, and crosses 80% of peak at **0.503% drift** (bracket
0.450-0.750%, since the envelope is only sampled at the protocol's amplitudes). The twin never
crossed 80% at all through 1.5%. At matched displacement: 0.870x the test envelope at 6.75 mm,
**0.491x at 13.5 mm** — half the wall already gone where the test was still carrying.

**So one unprinted parameter is worth a factor of three on drift capacity, and 9% on peak.** The
split is itself the result: peak strength is a TENSION-CRACKING quantity (D87) and barely notices;
capacity is a localization quantity and is dominated by it.

**2. WHY IT IS FAR BIGGER THAN THE ARITHMETIC SAID — the useful part.** The prediction made before
the run was +4.6% at 1% bar strain, from `E_h*eps`. That is the wrong quantity. `b` does not merely
add force, it adds **post-yield stiffness**, and post-yield stiffness is what prevents strain
LOCALIZING. Two ways to see it: (a) in a chain of EPP links the post-yield strain distribution is
INDETERMINATE — every distribution satisfies equilibrium at the same force — so the solver's answer
is set by whatever asymmetry it meets first; any `b > 0` makes force strictly increasing in strain
and restores uniqueness. (b) Mechanically, a hardened bar pushes its NEIGHBOUR into yield instead of
straining further itself, spreading yield over many rows; an EPP bar does not, so one row runs away
to `eps_su = 0.05` and ruptures. This is exactly why codes impose minimum hardening and uniform
elongation, and the model reproduced it unprompted. **A stability problem was estimated as an
arithmetic one.**

The same reasoning predicts `b = 1e-4` behaves like `b = 0`, not like 0.01: `E_h` = 20 MPa supplies
0.4 MPa at 2% bar strain = 0.11% of f_y, ~100x too weak against the disturbances an explicit march
on cracking concrete generates. Uniqueness that weak is uniqueness the solver cannot find. NOTE the
usual reason for `b = 1e-4` — avoiding a singular tangent in implicit Newton — does not apply here
at all, since CentralDifference never assembles a tangent (D74). A sweep at
**b = 0 / 1e-4 / 1e-3 / 1e-2** as parallel pushovers was launched 11:08 to measure this instead of
arguing it.

**3. THE TEST BRACKETS THE TWO.** Fig. 10(b) reaches -18.9 / +16.0 mm and the specimen got to about
20 mm, i.e. **~0.89% drift**. `b = 0` is finished by 0.50%; `b = 0.01` passes 1.5% undamaged. The
real wall lies between them — and this specimen is reinforced with **O8 welded mesh**, which
genuinely has low hardening and low uniform elongation, so 0.01 may be generous for it. `b` is
therefore not just an uncertainty, it is a FITTABLE parameter with a measured target, which is
exactly the status `eps_su` has (D91) and the reason both must be declared whenever a capacity is
quoted.

**4. THE CYCLIC CAPACITY METRIC FAILED AGAIN, AND IS NOW FIXED.** `drift_capacity` applied an
80%-drop rule to the trace itself; a cyclic trace crosses zero shear at EVERY reversal, so it fires
on the first unloading branch. It printed **0.6978%** on the `b = 0.01` ladder and **0.1934%** on
this one, both BELOW their own drift at peak, which is impossible for a capacity. Reported twice
before it was caught.

`metrics.py` now detects reversals (`turning_points`, a peak-valley filter with a 5%-of-amplitude
hysteresis — a monotonic push yields <= 1 candidate, the 8-level ladder yields 17), pairs the +/-
tips into one envelope point per AMPLITUDE LEVEL, and reads the 20%-drop crossing off that envelope
with linear interpolation, reporting `capacity_basis`, `peak_tip_shear` and
`drift_capacity_bracket`. Rescoring the finished cyclic runs: **all four earlier ones go from a
bogus 0.0035 / 0.453 / 0.453 / 0.698% to `None` — they never fell to 80%** — and this run reads
0.503%. Same class as D92: a metric quoted outside the domain it was written for. The basis is now
NAMED in every record so the domain is visible.

**5. WHAT THIS DOES TO EARLIER ENTRIES.** D97's capacity sweep (0.591 / 1.033 / 1.155 / 1.374%) and
D98's mesh-objectivity check were all run at `b = 0.01`; their ORDERING stands but every value
carries a hardening assumption larger than the +/-6% damping scatter and the 7.3% mesh shift
combined. D99's directional load-path reading is not refuted — both effects can hold — but
hardening is plainly the larger term. **SW-NC-FF and WSH3 carry the same `b = 0.01` convention**, and
SW-NC-FF's unexplained "no post-peak degradation, 1.34x at 4% drift" now has a second suspect that
costs one flag to test.

**Status:** accepted. Extends D91 (unprinted failure parameters), D92/D99 (metrics outside their
domain) and D97/D98 (capacity). Tests 55 pass, 1 known D34 failure.

---

### D102 — 2026-09-12 — The hardening sweep and the `b` x `eps_su` grid: `b = 1e-4` is `b = 0` to three decimals, the transition is a climb that accelerates, the two unprinted parameters INTERACT so the uncertainty stays two-dimensional, and peak strength refuses to move across a grid whose capacity spans 3.3x

**Why.** D101 left `b` as the largest known control on drift capacity and a number the paper never prints.
Two questions followed: where the transition sits (and whether `b = 1e-4`, the value people use to keep
an implicit solver out of a singular tangent, behaves like zero), and whether `b` and `eps_su` — both
unprinted, both controlling capacity — are one lever or two. Seven pushovers, mesh 50, `crushing/solved/
perfect`, `concrete_residual = 0`, zeta = 0.5, to 1.5% drift. **Parallel is free on this machine: three
at once ran at 71 steps/s against 68 solo**, so a sweep costs one run's wall time, not N.

**1. THE SWEEP (eps_su = 0.05).**

| b | peak kN | /test | drift at peak | capacity |
|---|---|---|---|---|
| 0 | 929.0 | 0.964 | 0.327% | **0.476%** |
| 1e-4 | 929.9 | 0.965 | 0.328% | **0.476%** |
| 1e-3 | 936.8 | 0.972 | 0.442% | 0.575% |
| 0.01 | 961.0 | 0.997 | 0.646% | **1.033%** |

**`b = 1e-4` IS `b = 0`** — capacity identical to three decimals, peak within 0.1%, drift at peak within
0.001%, curves superposed at every matched drift to 0.45% (ratios 1.0000 / 1.0013 / 1.0004 / 1.0007).
D101 predicted exactly this: any `b > 0` restores UNIQUENESS of the post-yield strain distribution in
principle, but `E_h` = 20 MPa supplies 0.4 MPa at 2% bar strain = 0.11% of f_y, ~100x too weak to
compete with the disturbances an explicit march on cracking concrete generates. **Uniqueness that weak
is uniqueness the solver cannot find.** NOTE the usual reason for `b = 1e-4` does not apply here at all:
CentralDifference never assembles a tangent (D74).

**THE TRANSITION IS NOT A THRESHOLD.** Capacity goes 0.476 -> 0.476 -> 0.575 -> 1.033% across four
decades: nothing from 0 to 1e-4, +21% to 1e-3, then **+80% in the last decade alone**. Almost everything
hardening buys this wall is bought between 0.001 and 0.01 — exactly where real reinforcement sits, which
is why the convention was never innocuous.

**2. A RESULT THAT TIES D99 TO D101.** At `b = 0` the monotonic and cyclic capacities AGREE — 0.476%
pushed against 0.503% cycled — while at `b = 0.01` they disagree completely (1.033% pushed, never
reached cycled). **So D99's cyclic enhancement exists ONLY when the ties harden.** Reversals can
re-shuffle a load path only if hardening gives them somewhere to shuffle it into; with EPP ties the wall
collapses the same way whichever path it is driven along. D99 must be read with this attached.

**3. THE GRID — the two levers INTERACT.**

| b | eps_su = 0.025 | eps_su = 0.05 | ratio |
|---|---|---|---|
| 0 | 0.308% | 0.476% | 1.546 |
| 1e-4 | 0.311% | 0.476% | 1.529 |
| 1e-3 | 0.365% | 0.575% | 1.576 |
| 0.0036 | **0.480%** | — | — |
| 0.01 | 0.592% | 1.033% | 1.745 |

Independent levers would give a CONSTANT ratio column. It does not stay constant — it reads
**1.546 / 1.529 / 1.576 / 1.745** — and the `b`-dependence is correspondingly flatter at the lower
ductility (1.92x across the row at eps_su = 0.025 against 2.17x at 0.05).

**CORRECTION, same day, when the missing `b = 0` cell landed (see 6).** This paragraph first read that
the ratio "climbs 1.529 -> 1.576 -> 1.745", written from three points. With the fourth the shape is
different: **1.546 / 1.529 / 1.576 are flat within scatter** (a 3% spread against the +/-6% run-to-run
scatter D97 measured), and only `b = 0.01` at **1.745** stands clear. So the interaction is NOT a gradual
drift with `b` — it is **concentrated in the last decade**, the same decade that carries +80% of the
hardening effect itself. The conclusion is unchanged (the levers interact, the uncertainty stays
two-dimensional); the mechanism reading sharpens, because whatever couples them only switches on where
hardening starts forcing rows to share. A trend was read into three points that four do not support. **THE UNCERTAINTY THEREFORE STAYS TWO-DIMENSIONAL — there is no single composite
parameter to quote a capacity against, and both must be declared every time one is.** Two obvious
candidates were tested and both fail: capacity is not a function of `b*eps_su`, and not a function of the
hardening stress gain `b*E*(eps_su - eps_y)` either — that one runs BACKWARDS (Deltasigma = 0.96 MPa
gives 0.476% while 4.64 MPa gives 0.365%).

**The interaction has the D101 mechanism behind it.** Every ratio is BELOW 2.0, where doubling eps_su
would double capacity exactly if the same number of bar rows participated. They do not: with more
ductility available the wall LOCALIZES FURTHER before failing, so fewer rows share the demand. Hardening
is what forces sharing, which is why the ratio climbs toward 2.0 as `b` rises. **Capacity ~ eps_su x
(rows participating), and `b` sets the second factor.**

**4. A BLIND PREDICTION, AND IT HELD.** Code-minimum O8 cold-worked welded mesh (EN 1992 Class A:
f_t/f_y ~ 1.05 at A_gt ~ 2.5%) gives `E_h` = 0.05*360/0.025 = 720 MPa, i.e. **b = 0.0036, eps_su =
0.025** — derived from the reinforcement class, not tuned. ~0.45% was put on the record before the run;
it returned **0.480%**, within 7%, and undershoots the test's ~0.89% by **1.85x**, also as predicted.
Note `b = 0.0036` lands inside the 0.003-0.005 window the sweep independently pointed at.
**THE FORK THIS OPENS CANNOT BE CLOSED WITH THIS SPECIMEN:** either the real mesh comfortably exceeded
code minimums (normal, and unverifiable from anything in the repo — the 2017 test paper is absent), or
the model carries a systematic capacity deficit. Say so wherever a capacity is quoted.

**5. PEAK STRENGTH REFUSES TO MOVE.** 926.6 to 961.0 kN — **0.962 to 0.997x** the measured 963.6 —
across a grid whose capacity spans **3.3x**. This is now the best-evidenced claim in the study and the
sharpest form of D87's split: peak is a TENSION-CRACKING quantity that no reinforcement parameter
reaches, capacity is a LOCALIZATION quantity that they jointly determine. It also strengthens the
reading of Aydin's 1.208x as constitutive rather than parametric — our peak will not go there for any
`b` or `eps_su`.

**6. A DESIGN SLIP, recorded — and its answer.** The eps_su = 0.025 row was framed as ALSO testing
whether the "`1e-4` = `0`" floor survives at lower ductility, but it was launched with `b = 1e-4` and no
`b = 0` partner, so that question went unanswered by construction. The general point is that a sweep
designed to answer two questions has to carry the control for both.

The missing cell returned **capacity 0.308% against 1e-4's 0.311% (1.0% apart) and peak 926.4 against
926.6 kN (0.02%)**. Against the +/-6% run-to-run scatter on capacity those are indistinguishable, so
**the floor is a property of `b` alone and holds at both ductilities** — 0.476/0.476 at eps_su = 0.05 and
0.308/0.311 at 0.025. The indeterminacy argument of D101 needs no qualification.

**Status:** accepted. Extends D101 (hardening), D99 (cyclic enhancement — now conditional on hardening),
D97 (the eps_su sweep was run at b = 0.01, the high corner of this grid) and D87 (peak vs capacity).
SW-NC-FF and WSH3 both carry `b = 0.01`. Tests 55 pass, 1 known D34 failure.
