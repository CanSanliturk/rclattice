"""Uniaxial compression of Aydin, Binici & Tuncay's (2021) 100 mm cube — tension-only lattice (D60).

Reproduces the paper's central claim: a lattice whose struts have NO compressive strength predicts
a concrete cube's compressive strength, because compression failure is an indirect tension failure
followed by loss of stability. Nothing in the model is told what fc is; the strength that comes out
is compared with the fc that was used only to pick ft, Ec and Gf.

Two things drive the result, and both are switchable so the paper's own control cases can be run:

  * `--rmax` — grid perturbation `Rmax/d`, the model's single calibration parameter. At `--rmax 0`
    the grid is uniform and the specimen LOCKS: his Fig. 3(a) climbs past 2*fc without failing.
  * `--rsm` / `--no-rsm` — the reduced-stiffness compression knee. Without it the paper still gets
    the strength but misses the strain at peak (his Fig. 3(b) plots both).

`--calibrate` runs his Fig. 4 outer loop: search `Rmax/d` until the mean strength of `n`
realizations lands within 10% of fc.

Because the mesh is random, one run is one sample. The default is 5 realizations, as in the paper
(p. 399), reported as a Table-2-style list with the mean.

Outputs: examples/output/aydin_cube/*.png + *_data.json. Units: N, mm.
"""

from __future__ import annotations

import json
import time

import numpy as np

from rclattice import viz
from rclattice.opensees import run_pushover_dynamic

from build import calibrate, cube_lattice, perturbed_grid, softening_report
from specimen import (
    A_FACE, CONCRETE, FC, GF, HORIZON, L, MESH, OUT, QUASI_STATIC_RATE, REALIZATIONS,
    RMAX_RATIO, RSM_ALPHA, TARGET_STRAIN, base_nodes, boundary_note, control_node, top_nodes,
)

# Search bounds for `--calibrate`. Above ~0.15 the perturbation approaches the node spacing and
# neighbouring nodes start to swap places; the paper's fitted values span 0.02-0.12 (his Fig. 5).
RMAX_BOUNDS = (0.005, 0.15)

EPS_CRACK = CONCRETE.ft / CONCRETE.E        # 7.447e-5 — a strut is cracked past this
EPS_KNEE = RSM_ALPHA * CONCRETE.epsc0       # 6.343e-4 — RSM stiffness break, NOT a crushing strain


def run_one(seed: int, *, rmax_ratio: float, boundary: str, rsm: bool, mesh_size: float,
            horizon: float, rate: float, target_strain: float, damping: float,
            strut_element: str, Gf: float, capture: bool = False) -> dict:
    """One realization: perturb, calibrate, assemble, squash. Returns the recorded response."""
    coords, quads, pairs = perturbed_grid(mesh_size, rmax_ratio=rmax_ratio, seed=seed,
                                          horizon=horizon, boundary=boundary)
    cal = calibrate(coords, pairs)
    model = cube_lattice(coords, quads, pairs, cal.area, mesh_size=mesh_size, horizon=horizon,
                         boundary=boundary, rsm=rsm, Gf=Gf, strut_element=strut_element)

    t0 = time.time()
    res = run_pushover_dynamic(
        model, control_node=control_node(model), control_dof=2,
        target=-target_strain * L, drive_nodes=top_nodes(model), base_nodes=base_nodes(model),
        rate=rate, steps_per_period=30, damping_ratio=damping, capture=capture,
        algorithm=("Newton",), max_iter=100,
    )
    elapsed = time.time() - t0
    if not res["disp"]:
        raise RuntimeError(f"seed {seed}: the compression run produced no steps")

    # Compression is a negative control displacement and a negative reaction sum; report magnitudes.
    shortening = [abs(u) for u in res["disp"]]
    stress = [abs(s) / A_FACE for s in res["shear"]]
    peak = max(stress)
    ipk = int(np.argmax(stress))
    # A maximum sitting at (or next to) the last step is where the run stopped, not a capacity.
    # It matters here more than usual: the uniform-grid control case never turns over at all, and
    # its "strength" would otherwise be read as a number rather than as the absence of one.
    still_rising = ipk >= len(stress) - 2
    return {
        "still_rising": still_rising,
        "seed": seed, "model": model, "cal": cal, "res": res, "elapsed": elapsed,
        "shortening": shortening, "strain": [u / L for u in shortening], "stress": stress,
        "peak": peak, "disp_at_peak": shortening[ipk], "strain_at_peak": shortening[ipk] / L,
        "converged": res["converged"],
    }


def thin(values, keep: int = 3000) -> list[float]:
    """Every k-th sample, plus the last, so a saved curve stays plottable without being enormous.

    A quasi-static dynamic-relaxation run records ~80,000 steps; at five realizations that is a
    16 MB JSON of points no plot can resolve. The peak is preserved to within one stride, which is
    ~0.0002% strain — far finer than the realization-to-realization scatter the file exists to hold.
    """
    v = list(values)
    if len(v) <= keep:
        return [float(x) for x in v]
    step = len(v) // keep + 1
    out = v[::step]
    if out[-1] != v[-1]:
        out.append(v[-1])
    return [float(x) for x in out]


def mean_strength(runs) -> float:
    return float(np.mean([r["peak"] for r in runs]))


def main(*, rmax_ratio: float | None = None, boundary: str = "LR", rsm: bool = True,
         realizations: int = REALIZATIONS, mesh_size: float = MESH, horizon: float = HORIZON,
         rate: float = QUASI_STATIC_RATE, target_strain: float = TARGET_STRAIN,
         damping: float = 0.8, seed0: int = 0, strut_element: str = "corotTruss",
         Gf: float = GF, do_calibrate: bool = False, tol: float = 0.10,
         max_rounds: int = 6) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if rmax_ratio is None:
        rmax_ratio = RMAX_RATIO[boundary]

    print(f"Aydin, Binici & Tuncay (2021) — {L:.0f} mm cube, uniaxial compression")
    print(f"  {boundary_note(boundary)}")
    print(f"  fc = {FC:g} MPa and Gf = {Gf * 1000:g} N/m CHOSEN; ft = 0.35*sqrt(fc) = "
          f"{CONCRETE.ft:.4f} MPa, Ec = 4700*sqrt(fc) = {CONCRETE.E:,.0f} MPa (his p. 396)")
    print(f"  grid d = {mesh_size:g} mm, horizon {horizon:g}, struts = {strut_element}, "
          f"compression = linear elastic{' + RSM knee at ' + f'{EPS_KNEE:.3e}' if rsm else ''} "
          f"(NO compressive strength)")
    print(f"  tension tail solved from Gf: {softening_report(mesh_size, Gf=Gf)}")
    print(f"  drive {rate:g} mm/s to {target_strain:.2%} strain, {realizations} realizations")

    if rmax_ratio == 0.0 and realizations > 1:
        print(f"  NOTE: Rmax/d = 0 is the uniform grid — it has no randomness, so all "
              f"{realizations} realizations would be the same model. Running 1.")
        realizations = 1

    rounds: list[tuple[float, float, list]] = []
    for it in range(max_rounds):
        print(f"\n=== Rmax/d = {rmax_ratio:.4f}  (Rmax = {rmax_ratio * mesh_size:.3f} mm) ===")
        runs = []
        for k in range(realizations):
            r = run_one(seed0 + 1000 * it + k, rmax_ratio=rmax_ratio, boundary=boundary, rsm=rsm,
                        mesh_size=mesh_size, horizon=horizon, rate=rate,
                        target_strain=target_strain, damping=damping,
                        strut_element=strut_element, Gf=Gf, capture=(k == 0))
            runs.append(r)
            flag = "" if r["converged"] else "   [solver stopped early]"
            flag += "   [STILL ASCENDING - no peak reached]" if r["still_rising"] else ""
            print(f"  trial {k + 1}: peak {r['peak']:6.2f} MPa at {r['disp_at_peak']:.4f} mm "
                  f"({r['strain_at_peak']:.3%} strain)   {r['elapsed']:5.1f}s{flag}")
        fmean = mean_strength(runs)
        err = (fmean - FC) / FC
        print(f"  mean of {realizations}: {fmean:.2f} MPa vs fc = {FC:g} MPa   -> "
              f"{fmean / FC:.3f} x fc  (error {err:+.1%})")
        # Only the newest round needs its full records (the report is written from it, and the
        # search reads the scatter of the current round). Older rounds keep their peaks alone —
        # six rounds of five realizations would otherwise pin ~100 MB of step histories in memory
        # for the sake of two numbers.
        if rounds:
            r_old, f_old, runs_old = rounds[-1]
            rounds[-1] = (r_old, f_old, [{"peak": x["peak"]} for x in runs_old])
        rounds.append((rmax_ratio, fmean, runs))

        if not do_calibrate or abs(err) < tol:
            if do_calibrate:
                print(f"  CONVERGED: |error| < {tol:.0%}, Aydin's Fig. 4 acceptance criterion")
            break
        rmax_ratio = _next_rmax(rounds, mesh_size)
        if rmax_ratio is None:
            print("  search exhausted — reporting the last round")
            break
    else:
        print(f"  reached the {max_rounds}-round limit without meeting the {tol:.0%} criterion")

    rmax_ratio, fmean, runs = rounds[-1]
    _report(runs, rmax_ratio=rmax_ratio, fmean=fmean, boundary=boundary, rsm=rsm,
            mesh_size=mesh_size, horizon=horizon, rate=rate, target_strain=target_strain,
            strut_element=strut_element, Gf=Gf, rounds=rounds)


def _next_rmax(rounds, mesh_size: float) -> float | None:
    """Next `Rmax/d` for the Fig. 4 loop: a BOUNDED, DAMPED secant on (Rmax/d -> strength).

    Strength falls monotonically with perturbation (his Figs. 5(a)-(b)), so two rounds give a
    secant. Three things stop it running away, all of which a plain secant did:

      * the strength estimate is NOISY — it is the mean of a handful of random meshes — so when two
        rounds differ by less than the scatter, the apparent slope is noise and the secant divides
        by it. Below `MIN_SLOPE` the step falls back to the proportional rule.
      * the step is damped to at most a factor of two per round in either direction, so one bad
        slope cannot throw the search across the whole range.
      * the result is clamped to `RMAX_BOUNDS`. Above roughly 0.15 the perturbation approaches the
        node spacing itself and neighbouring nodes begin to swap places, which is no longer the
        model being calibrated; the paper's own fitted values span 0.02-0.12 (his Fig. 5).
    """
    r1, f1, runs1 = rounds[-1]
    lo, hi = RMAX_BOUNDS

    step = None
    if len(rounds) >= 2:
        r0, f0, _ = rounds[-2]
        slope = (f1 - f0) / (r1 - r0) if abs(r1 - r0) > 1e-12 else 0.0
        # Scale the "is this slope real?" threshold by the scatter of the current round.
        peaks = [r["peak"] for r in runs1]
        noise = float(np.std(peaks)) if len(peaks) > 1 else 0.05 * FC
        if slope < 0.0 and abs(slope) * max(abs(r1 - r0), 1e-12) > max(noise, 1e-9):
            step = (FC - f1) / slope
    if step is None:
        step = r1 * (f1 / FC) - r1 if r1 > 1e-9 else 0.02

    # Damp: never move further than the current value (a factor of two either way).
    step = float(np.clip(step, -0.5 * max(r1, lo), max(r1, lo)))
    nxt = float(np.clip(r1 + step, lo, hi))
    if abs(nxt - r1) < 1e-4:
        nxt = float(np.clip(r1 * 1.5 if f1 > FC else r1 * 0.7, lo, hi))
    if abs(nxt - r1) < 1e-6:
        return None
    return nxt


def _report(runs, *, rmax_ratio, fmean, boundary, rsm, mesh_size, horizon, rate, target_strain,
            strut_element, Gf, rounds) -> None:
    peaks = np.array([r["peak"] for r in runs])
    strains = np.array([r["strain_at_peak"] for r in runs])
    lead = runs[0]
    cal = lead["cal"]

    print(f"\n--- results ({len(runs)} realizations, Rmax/d = {rmax_ratio:.4f}) ---")
    print(f"  strength   mean {peaks.mean():6.2f} MPa   min {peaks.min():6.2f}   "
          f"max {peaks.max():6.2f}   CoV {peaks.std(ddof=1) / peaks.mean():.1%}"
          if len(peaks) > 1 else f"  strength   {peaks[0]:.2f} MPa")
    print(f"  mean strength / fc = {fmean / FC:.3f}")
    print(f"  strain at peak      mean {strains.mean():.3e} = {strains.mean() / CONCRETE.epsc0:.2f}"
          f" x epsc0 (= 2fc/Ec)")
    print(f"  calibrated strut area {cal.area:,.1f} mm^2; lattice Poisson ratio "
          f"{cal.nu_effective:.3f}, cubic anisotropy {cal.cubic_anisotropy:.3f}")

    els, eps_end = viz.strut_strains(lead["model"], lead["res"]["disps_final"])
    ncr = int((eps_end >= EPS_CRACK).sum())
    print(f"  realization 1 damage: {ncr}/{len(els)} struts cracked ({ncr / len(els):.1%}); "
          f"NO strut ever crushed — the material has no compressive strength")

    if any(r["still_rising"] for r in runs):
        n = sum(r["still_rising"] for r in runs)
        print(f"  ** {n}/{len(runs)} realizations were STILL ASCENDING at the end of the run: for "
              f"those the number above is a lower bound on the strength, not the strength. **")
    if rmax_ratio == 0.0:
        print(f"\n  UNIFORM GRID CONTROL CASE (his Fig. 3(a)). With no perturbation the vertical "
              f"struts run straight from platen to platen, so once splitting has cracked away the "
              f"transverse restraint there is still nothing to destabilize them: the specimen "
              f"LOCKS. Measured {fmean / FC:.2f} x fc against his 'about twice the compressive "
              f"strength of the selected concrete grade'. This is the control case the grid "
              f"perturbation exists to fix, and it is what a structured-grid lattice does.")

    stem = f"aydin_cube_{boundary}{'_rsm' if rsm else '_norsm'}_rmax{rmax_ratio:.4f}"
    (OUT / f"{stem}_data.json").write_text(json.dumps({
        "paper": "Aydin, Binici & Tuncay (2021) MCR 73(8):394-409",
        "boundary": boundary, "rsm": rsm, "rmax_ratio": rmax_ratio, "mesh": mesh_size,
        "horizon": horizon, "strut_element": strut_element, "rate": rate, "Gf": Gf,
        "fc_target": FC, "ft": CONCRETE.ft, "Ec": CONCRETE.E, "epsc0": CONCRETE.epsc0,
        "target_strain": target_strain, "strut_area": cal.area,
        "nu_effective": cal.nu_effective, "cubic_anisotropy": cal.cubic_anisotropy,
        "peaks_MPa": peaks.tolist(), "mean_MPa": float(peaks.mean()),
        "mean_over_fc": float(peaks.mean() / FC),
        "strain_at_peak": strains.tolist(),
        "calibration_rounds": [{"rmax_ratio": r, "mean_MPa": f} for r, f, _ in rounds],
        "curves": [{"seed": r["seed"], "shortening_mm": thin(r["shortening"]),
                    "stress_MPa": thin(r["stress"]),
                    "converged": r["converged"], "still_rising": r["still_rising"]}
                   for r in runs],
    }))
    print(f"\nsaved raw response to {OUT / f'{stem}_data.json'}")

    # --- Fig. 3 style: stress vs shortening, every realization on shared axes -------------------
    curves = [{"disp": r["shortening"], "shear": r["stress"],
               "label": f"trial {i + 1} (peak {r['peak']:.1f} MPa)",
               "style": {"lw": 1.6, "alpha": 0.9}} for i, r in enumerate(runs)]
    span = [0.0, max(max(r["shortening"]) for r in runs)]
    curves.append({"disp": span, "shear": [FC, FC], "label": f"fc = {FC:g} MPa (target)",
                   "style": {"color": "k", "ls": "--", "lw": 1.8}})
    viz.figure_pushover(
        curves, savepath=str(OUT / f"{stem}.png"),
        xlabel="axial shortening (mm)", ylabel="axial stress (MPa)",
        title=f"Aydin et al. (2021) {L:.0f} mm cube — struts have NO compressive strength\n"
              f"{boundary}{' + RSM' if rsm else ', no RSM'},  Rmax/d = {rmax_ratio:.3f},  "
              f"d = {mesh_size:g} mm   |   mean {fmean:.1f} MPa = {fmean / FC:.2f} x fc")
    print(f"saved stress-shortening curves to {OUT / f'{stem}.png'}")

    # --- Fig. 3 insets: the crack pattern that IS the compression failure -----------------------
    scale = viz.autoscale_factor(lead["model"], lead["res"]["disps_final"], frac=0.10)
    viz.figure_damage(
        [(f"at peak load ({lead['peak']:.1f} MPa)", lead["model"], lead["res"]["disps_peak"]),
         (f"end of run ({lead['strain'][-1]:.2%} strain)", lead["model"],
          lead["res"]["disps_final"])],
        eps_crack=EPS_CRACK, eps_crush=EPS_KNEE, crush_label="past the RSM knee",
        scale=scale, savepath=str(OUT / f"{stem}_damage.png"),
        suptitle=f"Compression failure as indirect TENSION failure — {boundary}, "
                 f"Rmax/d = {rmax_ratio:.3f}, deformation x{scale:.0f}\n"
                 f"every red strut is a tensile crack; no strut can fail in compression")
    print(f"saved crack pattern to {OUT / f'{stem}_damage.png'}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Aydin, Binici & Tuncay (2021) 100 mm cube in uniaxial compression")
    p.add_argument("--boundary", choices=("LR", "NR"), default="LR",
                   help="end condition of his Fig. 2(a): LR = pinned base + laterally held top "
                        "platen (default), NR = rollers under the base and a free top")
    p.add_argument("--rmax", type=float, default=None, metavar="RATIO",
                   help="grid perturbation Rmax/d, the model's one calibration parameter. Default "
                        f"is the paper's fitted value for this specimen ({RMAX_RATIO['LR']} for LR, "
                        f"{RMAX_RATIO['NR']} for NR). Pass 0 for his uniform-grid control case")
    p.add_argument("--no-rsm", dest="rsm", action="store_false",
                   help="drop the reduced-stiffness compression knee (his Fig. 1(b)), leaving "
                        "compression purely linear elastic. Strength should survive; the strain at "
                        "peak should not")
    p.add_argument("--calibrate", action="store_true",
                   help="run his Fig. 4 outer loop: search Rmax/d until the mean strength of the "
                        "realizations is within 10%% of fc")
    p.add_argument("--realizations", type=int, default=REALIZATIONS,
                   help=f"simulations per Rmax/d (default {REALIZATIONS}, as in the paper). The "
                        "mesh is random, so a single run is a single sample")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing d, mm (default {MESH:g})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON:g} = 8 neighbours)")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"drive speed, mm/s (default {QUASI_STATIC_RATE:g})")
    p.add_argument("--strain", type=float, default=TARGET_STRAIN,
                   help=f"target axial strain (default {TARGET_STRAIN:g})")
    p.add_argument("--damping", type=float, default=0.8, help="damping ratio (default 0.8)")
    p.add_argument("--seed", type=int, default=0, help="base RNG seed (default 0)")
    p.add_argument("--truss", dest="strut_element", action="store_const", const="Truss",
                   default="corotTruss",
                   help="use small-displacement Truss struts instead of corotTruss. Expect the "
                        "specimen to lock: without geometric nonlinearity the split columns cannot "
                        "lose stability, which is the failure mode being modelled")
    p.add_argument("--gf", type=float, default=GF,
                   help=f"tensile fracture energy, N/mm (default {GF:g} = {GF * 1000:g} N/m)")
    a = p.parse_args()
    main(rmax_ratio=a.rmax, boundary=a.boundary, rsm=a.rsm, realizations=a.realizations,
         mesh_size=a.mesh, horizon=a.horizon, rate=a.rate, target_strain=a.strain,
         damping=a.damping, seed0=a.seed, strut_element=a.strut_element, Gf=a.gf,
         do_calibrate=a.calibrate)
