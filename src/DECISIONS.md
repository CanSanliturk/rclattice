
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
