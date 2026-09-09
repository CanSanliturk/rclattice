"""OUR 200 mm cube, modelled and analysed Aydin's way (D61).

The specimen is untouched: `specimen.py` supplies the same 200 mm cube, the same fc = 15.1 MPa
grade with its measured `E` and `ft`, the same 20 mm grid, the same smooth platens (rollers plus one
pin) and the same Aydin energy-balance strut area. Nothing about the cube changes. What changes is
the MODELLING AND ANALYSIS APPROACH, swapped wholesale for the one in Aydin, Binici & Tuncay (2021):

    pushover.py (D53-D59)                    this script (D61)
    ---------------------------------------  ---------------------------------------
    Concrete02: struts capped at fc          NO compressive strength anywhere; linear
                                             elastic + the RSM knee at epsc0/3
    bilinear tension softening (Ets)         trilinear tension tail, solved from Gf
    uniform structured grid                  grid perturbed by Rmax/d
    small-displacement Truss                 corotTruss (the failure is geometric)
    deterministic, one run is the answer     random mesh; realizations, reported as a mean
    fc is an INPUT; the answer is 0.634 fc   fc is an OUTPUT; Rmax/d is calibrated to it

Read `examples/aydin_cube/DIFFERENCES.md` for why those six move together rather than separately.
The point of running them on the SAME cube is that every remaining difference in the result is
attributable to the approach, since the specimen is held fixed.

`--calibrate` runs Aydin's Fig. 4 loop to find `Rmax/d` for THIS cube: his fitted values (0.06-0.075)
are for his fc = 20 MPa, Gf = 50 N/m, d = 10 mm specimen and do not transfer to ours, which has a
lower target strength (15.1) but a HIGHER tensile strength (measured 1.5 against his 0.35*sqrt(fc) =
1.36) — both of which push the required perturbation up.

Each run writes its own timestamped directory under examples/output/compression_cube/runs/,
named `<stamp>_aydin_rmax<...>`, alongside the `..._pushover_...` directories of the Concrete02
runs; the stress-strain figure overlays the NEWEST completed pushover.py run when there is one.
Units: N, mm.
"""

from __future__ import annotations

import json
import time

import numpy as np

from rclattice import viz
from rclattice.calibration import energy_balance_area
from rclattice.builders import build_lattice_rc
from rclattice.materials import aydin_lattice_softening, concrete_lattice_aydin
from rclattice.mesh import connect_horizon, mesh_rectangle_grid, perturb_nodes
from rclattice.opensees import run_pushover_dynamic

from build import expected_modulus
from specimen import (
    A_FACE, CONCRETE, EPS, GF, HORIZON, L, MESH, NU, QUASI_STATIC_RATE, THK,
    base_nodes, check_mesh_alignment, control_node, cube_problem, find_run, run_dir, top_nodes,
)

# --- Aydin's model parameters (2021, Fig. 1 caption and Table 1) --------------------------------
A1, B1, B2 = 1.5, 0.6, 0.2       # trilinear tension shape; a2/a3 are solved per strut from Gf
A3_OVER_A2 = 5.0                 # Jansen & Shah column of his Table 1 (300/60)
RSM_ALPHA, RSM_BETA = 1.0 / 3.0, 0.4   # reduced stiffness: 0.4*E beyond epsc0/3

# Calibrated for THIS cube by `--calibrate` (D61): the search ran 0.050 -> 0.062 -> 0.093 -> 0.075
# and converged at 0.919*fc, inside Aydin's 10% criterion. His own values (0.060 NR, 0.075 LR) are
# for a different specimen and were not assumed — that the answer lands on top of his LR figure is a
# coincidence of two specimens that differ in size, grade, grid spacing and end condition.
RMAX_RATIO = 0.075
RMAX_BOUNDS = (0.005, 0.20)
REALIZATIONS = 5

# The Aydin model fails and sheds load far earlier than Concrete02's residual plateau, so the run
# does not need pushover.py's 1% strain. 0.5% is ~2.7*epsc0 — well past peak and well into softening.
TARGET_STRAIN = 0.005

EPS_CRACK = CONCRETE.ft / CONCRETE.E          # 9.317e-5 — a strut is cracked past this
EPS_KNEE = RSM_ALPHA * CONCRETE.epsc0         # 6.253e-4 — RSM break, NOT a crushing strain
STEEP_DEG = 30.0                              # "near-vertical" cut-off for the load-path split

# Upper strain of the initial-tangent fit. It has to clear TWO events, not one: the RSM knee at
# `EPS_KNEE` and, earlier, the first transverse cracking at `eps_cr/nu ~ 2.4e-4`. `pushover.py` uses
# the same 2.0e-4 for the second reason alone, so the two approaches measure their modulus over the
# same window and the check below compares approaches rather than fit ranges.
LINEAR_FIT_STRAIN = min(2.0e-4, 0.5 * EPS_KNEE)


def perturbed_grid(mesh_size: float = MESH, *, rmax_ratio: float, seed: int | None = None,
                   horizon: float = HORIZON):
    """`(coords, quads, pairs)` for one realization of OUR cube with Aydin's perturbed mesh.

    Two things are pinned so the specimen stays the one `specimen.py` describes:

      * strut TOPOLOGY comes from the unperturbed grid, so perturbation randomizes geometry only
        (re-running `connect_horizon` on moved nodes drops struts that cross `horizon*mesh_size`);
      * the two mid-width nodes on the bottom and top faces do not move. `specimen.py` finds both by
        an exact box at `x = L/2`: the bottom one carries the single ux restraint that keeps the
        roller-supported cube non-singular, and the top one is the displacement-control node. Their
        neighbours sit a full grid spacing away and move by at most `Rmax`, so the three support
        boxes still partition the bottom edge exactly as they do on the regular grid, and
        `specimen.py` needs no change at all — which is what "the cube is intact" has to mean.
        Holding two of 121 nodes fixed is not a modelling choice with consequences; every top node
        is driven by the same ramp regardless, so the control node is only a place to read.
    """
    check_mesh_alignment(mesh_size)
    coords0, quads = mesh_rectangle_grid(L, L, mesh_size)
    pairs = connect_horizon(coords0, mesh_size, horizon)

    fix = set()
    for y in (0.0, L):
        face = np.flatnonzero(np.abs(coords0[:, 1] - y) < EPS)
        fix.add(int(face[np.abs(coords0[face, 0] - L / 2.0).argmin()]))
    coords = perturb_nodes(coords0, mesh_size, rmax_ratio, seed=seed, fix_nodes=fix)
    return coords, quads, pairs


def calibrate_area(coords, pairs):
    """Aydin's elastic energy balance on this realization's own geometry.

    Identical method to `build.calibrate` — the same Eq. 2.1, the same one uniform `EA` — but
    evaluated on the perturbed node positions, because perturbation changes every strut length.
    The elastic calibration is the one part of the approach that does NOT change between
    `pushover.py` and this script, which is why the initial-tangent check is shared.
    """
    return energy_balance_area(coords, pairs, E=CONCRETE.E, nu=NU, thickness=THK,
                               area_inplane=L * L)


def aydin_lattice(coords, quads, pairs, area: float, *, mesh_size: float = MESH,
                  horizon: float = HORIZON, rsm: bool = True, Gf: float = GF,
                  strut_element: str = "corotTruss"):
    """Our cube assembled with Aydin's tension-only struts on the perturbed grid."""
    def material_for(_zone: str, length: float):
        return concrete_lattice_aydin(CONCRETE, 0, length, Gf=Gf, a1=A1, b1=B1, b2=B2,
                                      a3_over_a2=A3_OVER_A2, rsm=rsm,
                                      alpha=RSM_ALPHA, beta=RSM_BETA)

    return build_lattice_rc(
        cube_problem(mesh_size), mesh_size,
        material_for=material_for, zone_of=lambda _x, _y: "concrete", rebars=(),
        strut_area=area, horizon=horizon, strut_element=strut_element,
        grid=(coords, quads), pairs=pairs,
    )[0]


def cut_element_groups(model, y_cut: float, *, steep_deg: float = STEEP_DEG) -> dict:
    """Struts crossing `y_cut`, split into near-vertical and inclined, as `element_groups`.

    `pushover.py` splits on `abs(dx) < tol`, which is exact on a structured grid and useless here:
    perturbation leaves NO strut exactly vertical. The split is therefore by ANGLE from vertical,
    with `steep_deg` chosen so the struts that were vertical before perturbation stay in the
    "near-vertical" group (a 2 mm move on a 20 mm strut tilts it by at most ~11 degrees) and the
    45-degree diagonals stay out of it.

    As in `pushover.py`, each entry records the VERTICAL force the strut exerts on its LOWER node —
    `eleForce` returns [Fxi, Fyi, Fxj, Fyj], so index 1 or 3 depending on which end sits below.
    """
    cos_lim = np.cos(np.deg2rad(steep_deg))
    groups: dict[str, list] = {"near-vertical struts": [], "inclined struts": []}
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        if (ya - y_cut) * (yb - y_cut) >= 0.0:
            continue
        dx, dy = xb - xa, yb - ya
        steep = abs(dy) / (dx * dx + dy * dy) ** 0.5 >= cos_lim
        dof = 1 if ya < yb else 3
        groups["near-vertical struts" if steep else "inclined struts"].append((e.id, dof, 1.0))
    return groups


# A quasi-static record can FALL fast — that is what an instability looks like — but it cannot RISE
# fast: the stress increment per step is bounded by the stiffness times one step of imposed
# displacement, and the run takes tens of thousands of steps to cross the whole response. A jump of
# this size in a single step is therefore a diverged transient, not mechanics.
MAX_RISE_PER_STEP = 0.05      # fraction of fc

# How far below its maximum the curve must fall for that maximum to count as a strength.
DESCENT_FRAC = 0.05


def first_divergence(stress, fc: float) -> int:
    """Index of the first spurious upward jump, or `len(stress)` if the record is clean.

    Needed because a blown-up step is not always a failed one: on our cube a realization diverged
    locally, recorded a 139 MPa spike — 9x fc, on a lattice whose struts have no compressive
    strength to give it — recovered, and still reported `converged = True`. Averaged in with the
    others it turned a 13 MPa mean into 56 MPa and sent the Rmax search to its upper bound. A
    diverged solve is a failed MEASUREMENT, not a sample of the model, and has to be discarded
    rather than averaged.
    """
    s = np.asarray(stress)
    if len(s) < 2:
        return len(s)
    bad = np.flatnonzero(np.diff(s) > MAX_RISE_PER_STEP * fc)
    return int(bad[0]) + 1 if len(bad) else len(s)


def initial_tangent(strain, stress, upper: float) -> float:
    """Least-squares secant through the origin over the uncracked window."""
    x, y = np.asarray(strain), np.asarray(stress)
    m = (x > 0.0) & (x <= upper)
    if m.sum() < 2:
        return float("nan")
    return float(np.dot(x[m], y[m]) / np.dot(x[m], x[m]))


def run_one(seed: int, *, rmax_ratio: float, rsm: bool, mesh_size: float, horizon: float,
            rate: float, target_strain: float, damping: float, strut_element: str, Gf: float,
            capture: bool = False) -> dict:
    """One realization of our cube under Aydin's approach."""
    coords, quads, pairs = perturbed_grid(mesh_size, rmax_ratio=rmax_ratio, seed=seed,
                                          horizon=horizon)
    cal = calibrate_area(coords, pairs)
    model = aydin_lattice(coords, quads, pairs, cal.area, mesh_size=mesh_size, horizon=horizon,
                          rsm=rsm, Gf=Gf, strut_element=strut_element)

    y_cut = L / 2.0 + mesh_size / 2.0
    t0 = time.time()
    res = run_pushover_dynamic(
        model, control_node=control_node(model), control_dof=2,
        target=-target_strain * L, drive_nodes=top_nodes(model), base_nodes=base_nodes(model),
        rate=rate, steps_per_period=30, damping_ratio=damping, capture=capture,
        element_groups=cut_element_groups(model, y_cut), algorithm=("Newton",), max_iter=100,
    )
    elapsed = time.time() - t0
    if not res["disp"]:
        raise RuntimeError(f"seed {seed}: the compression run produced no steps")

    strain = [abs(u) / L for u in res["disp"]]
    stress = [abs(s) / A_FACE for s in res["shear"]]

    # Truncate at the first diverged step, then judge the record that survives.
    icut = first_divergence(stress, CONCRETE.fc)
    diverged = icut < len(stress)
    strain, stress = strain[:icut], stress[:icut]
    ipk = int(np.argmax(stress))
    # "Turned over" has to mean a measurable DESCENT, not merely a maximum that is not the last
    # sample. A record holds tens of thousands of steps, so a peak three samples from the end
    # satisfies the naive test while telling us nothing; requiring the curve to fall `DESCENT_FRAC`
    # below its maximum afterwards is what distinguishes a strength from a record that stopped.
    descended = float(np.min(stress[ipk:])) < (1.0 - DESCENT_FRAC) * stress[ipk]
    still_rising = not descended
    return {
        "seed": seed, "model": model, "cal": cal, "res": res, "elapsed": elapsed, "y_cut": y_cut,
        "strain": strain, "stress": stress, "peak": stress[ipk], "strain_at_peak": strain[ipk],
        "converged": res["converged"], "still_rising": still_rising,
        "diverged": diverged,
        # A realization counts only if it turned over BEFORE anything went wrong numerically.
        "valid": not still_rising,
        "E0": initial_tangent(strain, stress, LINEAR_FIT_STRAIN),
    }


def thin(values, keep: int = 3000) -> list[float]:
    """Every k-th sample plus the last — a quasi-static run records tens of thousands of steps."""
    v = list(values)
    if len(v) <= keep:
        return [float(x) for x in v]
    out = v[::len(v) // keep + 1]
    if out[-1] != v[-1]:
        out.append(v[-1])
    return [float(x) for x in out]


def next_rmax(rounds) -> float | None:
    """Bounded, damped secant for Aydin's Fig. 4 loop (same guards as `examples/aydin_cube`).

    Strength falls monotonically with perturbation, but the estimate is the mean of a few random
    meshes, so a slope smaller than the scatter is noise and a plain secant divides by it. Below
    that threshold the step falls back to a proportional rule; either way it is damped to a factor
    of two per round and clamped, since beyond `RMAX_BOUNDS[1]` the perturbation approaches the node
    spacing and neighbours begin to swap places.
    """
    r1, f1, runs1 = rounds[-1]
    lo, hi = RMAX_BOUNDS
    target = CONCRETE.fc

    step = None
    if len(rounds) >= 2:
        r0, f0, _ = rounds[-2]
        slope = (f1 - f0) / (r1 - r0) if abs(r1 - r0) > 1e-12 else 0.0
        peaks = [r["peak"] for r in runs1]
        noise = float(np.std(peaks)) if len(peaks) > 1 else 0.05 * target
        if slope < 0.0 and abs(slope) * max(abs(r1 - r0), 1e-12) > max(noise, 1e-9):
            step = (target - f1) / slope
    if step is None:
        step = r1 * (f1 / target) - r1 if r1 > 1e-9 else 0.02

    step = float(np.clip(step, -0.5 * max(r1, lo), max(r1, lo)))
    nxt = float(np.clip(r1 + step, lo, hi))
    if abs(nxt - r1) < 1e-4:
        nxt = float(np.clip(r1 * 1.5 if f1 > target else r1 * 0.7, lo, hi))
    return None if abs(nxt - r1) < 1e-6 else nxt


def main(*, rmax_ratio: float = RMAX_RATIO, rsm: bool = True, realizations: int = REALIZATIONS,
         mesh_size: float = MESH, horizon: float = HORIZON, rate: float = QUASI_STATIC_RATE,
         target_strain: float = TARGET_STRAIN, damping: float = 0.8, seed0: int = 0,
         strut_element: str = "corotTruss", Gf: float = GF, do_calibrate: bool = False,
         tol: float = 0.10, max_rounds: int = 6) -> None:
    # The tag records the REQUESTED configuration. Under --calibrate the final Rmax/d is an OUTPUT
    # of the search, so it cannot name the directory; it is printed, and saved into data.json.
    tag = ("calib_" if do_calibrate else "") + f"rmax{rmax_ratio:.4f}"
    if not rsm:
        tag += "_norsm"
    if strut_element != "corotTruss":
        tag += f"_{strut_element}"
    out = run_dir("aydin", tag=tag)
    print(f"output directory: {out}\n")
    fc = CONCRETE.fc

    print(f"OUR {L:.0f} mm cube, modelled Aydin's way (2021) — specimen unchanged, approach swapped")
    print(f"  cube: {L:.0f} mm, d = {mesh_size:g} mm, horizon {horizon:g}, smooth platens "
          f"(rollers + one pin) — exactly as in specimen.py")
    print(f"  grade: fc = {fc:g} MPa, E = {CONCRETE.E:,.0f} MPa, ft = {CONCRETE.ft:g} MPa "
          f"(measured, not 0.35*sqrt(fc) = {0.35 * fc ** 0.5:.2f}), Gf = {Gf * 1000:g} N/m")
    print(f"  APPROACH: tension-only struts (NO compressive strength)"
          f"{f' + RSM knee at {EPS_KNEE:.3e}' if rsm else ', no RSM'}, "
          f"perturbed grid, {strut_element}")
    for label, length in (("orthogonal", mesh_size), ("diagonal", mesh_size * 2 ** 0.5)):
        a2, a3 = aydin_lattice_softening(CONCRETE.ft, CONCRETE.E, length, Gf=Gf,
                                         a1=A1, b1=B1, b2=B2, a3_over_a2=A3_OVER_A2)
        print(f"    tension tail, {label} L = {length:5.2f} mm -> a2 = {a2:5.1f}, a3 = {a3:6.1f}")
    print(f"  drive {rate:g} mm/s to {target_strain:.2%} strain, {realizations} realizations")

    if rmax_ratio == 0.0 and realizations > 1:
        print(f"  NOTE: Rmax/d = 0 is the uniform grid and has no randomness — running 1.")
        realizations = 1

    rounds: list[tuple[float, float, list]] = []
    for it in range(max_rounds):
        print(f"\n=== Rmax/d = {rmax_ratio:.4f}  (Rmax = {rmax_ratio * mesh_size:.3f} mm) ===")
        runs = []
        for k in range(realizations):
            r = run_one(seed0 + 1000 * it + k, rmax_ratio=rmax_ratio, rsm=rsm,
                        mesh_size=mesh_size, horizon=horizon, rate=rate,
                        target_strain=target_strain, damping=damping,
                        strut_element=strut_element, Gf=Gf, capture=(k == 0))
            runs.append(r)
            flag = "" if r["converged"] else "   [solver stopped early]"
            flag += f"   [diverged at {r['strain'][-1]:.3%} strain — record truncated there]" \
                if r["diverged"] else ""
            flag += "   [NO PEAK -> DISCARDED]" if not r["valid"] else ""
            print(f"  trial {k + 1}: peak {r['peak']:6.2f} MPa at {r['strain_at_peak']:.3%} strain "
                  f"({r['strain_at_peak'] / CONCRETE.epsc0:.2f} epsc0)   {r['elapsed']:5.1f}s{flag}")

        # Two very different reasons a realization can lack a peak, and they must not be conflated:
        #   LOCKED    — the curve rose monotonically to the end of the run with nothing numerically
        #               wrong. That is a physical result, and at Rmax/d = 0 it is THE result.
        #   DIVERGED  — the transient blew up before the specimen turned over. A failed measurement.
        good = [r for r in runs if r["valid"]]
        locked = [r for r in runs if r["still_rising"] and not r["diverged"]]
        lower_bound = False
        if good:
            if len(good) < len(runs):
                print(f"  {len(runs) - len(good)}/{len(runs)} realizations discarded (no peak "
                      f"before divergence); the mean below uses the remaining {len(good)}")
            runs = good
        elif len(locked) == len(runs):
            lower_bound = True
            print(f"  every realization LOCKED — the curve was still rising at {target_strain:.2%} "
                  f"strain, so the values below are LOWER BOUNDS, not strengths.")
        else:
            raise RuntimeError(
                f"every realization at Rmax/d = {rmax_ratio:.4f} diverged before reaching a peak. "
                f"Lower --rate, raise --damping, or refine --mesh: our 20 mm grid gives struts a "
                f"much shorter softening tail than Aydin's 10 mm one, so each strut failure is a "
                f"more violent event for the transient solver to absorb.")
        realizations_used = len(runs)
        fmean = float(np.mean([r["peak"] for r in runs]))
        err = (fmean - fc) / fc
        print(f"  mean of {realizations_used}: {fmean:.2f} MPa vs fc = {fc:g} MPa   -> "
              f"{fmean / fc:.3f} x fc  (error {err:+.1%})")

        if rounds:   # superseded rounds keep their peaks only; the step histories are large
            r_old, f_old, runs_old = rounds[-1]
            rounds[-1] = (r_old, f_old, [{"peak": x["peak"]} for x in runs_old])
        rounds.append((rmax_ratio, fmean, runs))

        if not do_calibrate or abs(err) < tol:
            if do_calibrate:
                print(f"  CONVERGED: |error| < {tol:.0%}, Aydin's Fig. 4 acceptance criterion")
            break
        nxt = next_rmax(rounds)
        if nxt is None:
            print("  search exhausted — reporting the last round")
            break
        rmax_ratio = nxt
    else:
        print(f"  reached the {max_rounds}-round limit without meeting the {tol:.0%} criterion")

    rmax_ratio, fmean, runs = rounds[-1]
    _report(runs, out=out, rmax_ratio=rmax_ratio, fmean=fmean, lower_bound=lower_bound, rsm=rsm,
            mesh_size=mesh_size,
            horizon=horizon, rate=rate, target_strain=target_strain,
            strut_element=strut_element, Gf=Gf, rounds=rounds)


def _report(runs, *, out, rmax_ratio, fmean, lower_bound, rsm, mesh_size, horizon, rate,
            target_strain, strut_element, Gf, rounds) -> None:
    fc = CONCRETE.fc
    peaks = np.array([r["peak"] for r in runs])
    strains = np.array([r["strain_at_peak"] for r in runs])
    tangents = np.array([r["E0"] for r in runs])
    lead, cal = runs[0], runs[0]["cal"]
    e_expected = expected_modulus(cal)

    print(f"\n--- results ({len(runs)} realizations, Rmax/d = {rmax_ratio:.4f}) ---")
    if len(peaks) > 1:
        print(f"  strength   mean {peaks.mean():6.2f} MPa   min {peaks.min():6.2f}   "
              f"max {peaks.max():6.2f}   CoV {peaks.std(ddof=1) / peaks.mean():.1%}")
    else:
        print(f"  strength   {peaks[0]:.2f} MPa")
    print(f"  mean strength / fc = {fmean / fc:.3f}"
          f"{'   (LOWER BOUND — never turned over)' if lower_bound else ''}")
    if lower_bound and rmax_ratio == 0.0:
        print(f"\n  UNIFORM-GRID CONTROL. With no perturbation the near-vertical struts run "
              f"straight from platen to platen, and since this approach gives them no compressive "
              f"strength there is nothing to stop them: the specimen LOCKS and simply keeps taking "
              f"load. This is Aydin's Fig. 3(a) on our cube, and it is what makes the grid "
              f"perturbation necessary rather than decorative.")
    print(f"  strain at peak     mean {strains.mean():.3e} = {strains.mean() / CONCRETE.epsc0:.2f}"
          f" x epsc0")

    # The elastic calibration is unchanged from pushover.py, so its check carries over verbatim —
    # and NO secant correction is needed here, because Aydin's compression is linear below the RSM
    # knee where Concrete02's is parabolic. That the same target is met by both approaches is what
    # makes the strength comparison below a comparison of the approaches and not of two calibrations.
    print(f"  initial tangent    mean {np.nanmean(tangents):9,.0f} MPa   vs expected "
          f"{e_expected:,.0f} MPa  -> ratio {np.nanmean(tangents) / e_expected:.4f}")
    print(f"                                          vs grade E {CONCRETE.E:,.0f} MPa"
          f"  -> ratio {np.nanmean(tangents) / CONCRETE.E:.4f}")
    print(f"  calibrated strut area {cal.area:,.1f} mm^2; lattice nu {cal.nu_effective:.3f}")

    els, eps_end = viz.strut_strains(lead["model"], lead["res"]["disps_final"])
    ncr = int((eps_end >= EPS_CRACK).sum())
    print(f"  realization 1 damage: {ncr}/{len(els)} struts cracked ({ncr / len(els):.1%}); "
          f"NO strut ever crushed — under this approach the material cannot fail in compression")
    if any(r["still_rising"] for r in runs):
        n = sum(r["still_rising"] for r in runs)
        print(f"  ** {n}/{len(runs)} realizations were STILL ASCENDING at the end — for those the "
              f"peak above is a lower bound, not a strength **")

    # --- who carries the load, and does the inclined path still collapse? -----------------------
    gv = np.abs(np.asarray(lead["res"]["groups"]["near-vertical struts"]))
    gd = np.abs(np.asarray(lead["res"]["groups"]["inclined struts"]))
    tot = np.where(gv + gd > 1e-9, gv + gd, 1e-9)
    ipk = int(np.argmax(lead["stress"]))
    eps_split = EPS_CRACK / max(cal.nu_effective, 1e-9)
    pre = np.flatnonzero(np.asarray(lead["strain"]) < eps_split)
    iel = int(pre[-1]) if len(pre) else 1
    print(f"\n  load path across y = {lead['y_cut']:.0f} mm (realization 1, split at "
          f"{STEEP_DEG:.0f} deg from vertical)")
    print(f"    transverse cracking begins at eps_cr/nu = {eps_split:.2e} "
          f"= {eps_split / CONCRETE.epsc0:.1%} of epsc0")
    for label, i in (("pre-crack", iel), ("at peak", ipk), ("end of run", len(tot) - 1)):
        print(f"    {label:<11} near-vertical {gv[i] / 1e3:8.1f} kN ({gv[i] / tot[i]:6.1%})   "
              f"inclined {gd[i] / 1e3:8.1f} kN ({gd[i] / tot[i]:6.1%})")
    print(f"    -> Concrete02 on the uniform grid sheds the inclined path completely "
          f"(28% -> 0%, D54/D55); here it holds {gd[ipk] / tot[ipk]:.0%} at peak.")

    # --- outputs --------------------------------------------------------------------------------
    (out / "data.json").write_text(json.dumps({
        "approach": "Aydin, Binici & Tuncay (2021) tension-only lattice on the compression_cube "
                    "specimen",
        "rmax_ratio": rmax_ratio, "rsm": rsm, "mesh": mesh_size, "horizon": horizon,
        "strut_element": strut_element, "rate": rate, "Gf": Gf, "target_strain": target_strain,
        "fc": fc, "ft": CONCRETE.ft, "E": CONCRETE.E, "epsc0": CONCRETE.epsc0,
        "strut_area": cal.area, "nu_effective": cal.nu_effective,
        "peaks_MPa": peaks.tolist(), "mean_MPa": float(peaks.mean()),
        "mean_over_fc": float(peaks.mean() / fc), "strain_at_peak": strains.tolist(),
        "E0_measured": tangents.tolist(), "E_expected": e_expected,
        "calibration_rounds": [{"rmax_ratio": r, "mean_MPa": f} for r, f, _ in rounds],
        "curves": [{"seed": r["seed"], "strain": thin(r["strain"]), "stress_MPa": thin(r["stress"]),
                    "converged": r["converged"]} for r in runs],
    }))
    print(f"\nsaved raw response to {out / 'data.json'}")

    # Overlay the Concrete02 baseline on the SAME cube: the NEWEST completed pushover.py run, if
    # there is one. Named in the caption, because which baseline was overlaid is part of the result.
    curves = [{"disp": r["strain"], "shear": r["stress"],
               "label": f"Aydin approach, trial {i + 1} ({r['peak']:.1f} MPa)",
               "style": {"lw": 1.4, "alpha": 0.85}} for i, r in enumerate(runs)]
    note = ""
    try:
        baseline = find_run("pushover") / "data.json"
    except SystemExit:
        baseline = None
        print("  (no completed pushover.py run to overlay — run it for the Concrete02 baseline)")
    if baseline is not None:
        b = json.loads(baseline.read_text())
        curves.append({"disp": b["strain"], "shear": b["stress_MPa"],
                       "label": f"Concrete02, uniform grid ({b['peak_MPa']:.1f} MPa "
                                f"= {b['peak_MPa'] / fc:.2f} fc)",
                       "style": {"color": "k", "lw": 2.4}})
        note = f"   vs   Concrete02 {b['peak_MPa'] / fc:.2f} fc"
        print(f"  Concrete02 baseline overlaid from {baseline.parent.name}")
    span = [0.0, max(max(r["strain"]) for r in runs)]
    curves.append({"disp": span, "shear": [fc, fc], "label": f"fc = {fc:g} MPa",
                   "style": {"color": "C3", "ls": "--", "lw": 1.6}})
    viz.figure_pushover(
        curves, savepath=str(out / "stress_strain.png"),
        xlabel="axial strain (-)", ylabel="axial stress (MPa)",
        title=f"Our {L:.0f} mm cube — same specimen, two modelling approaches\n"
              f"Aydin tension-only (Rmax/d = {rmax_ratio:.3f}): {fmean / fc:.2f} fc{note}")
    print(f"saved stress-strain comparison to {out / 'stress_strain.png'}")

    scale = viz.autoscale_factor(lead["model"], lead["res"]["disps_final"], frac=0.10)
    viz.figure_damage(
        [(f"at peak load ({lead['peak']:.1f} MPa)", lead["model"], lead["res"]["disps_peak"]),
         (f"end of run ({lead['strain'][-1]:.2%} strain)", lead["model"],
          lead["res"]["disps_final"])],
        eps_crack=EPS_CRACK, eps_crush=EPS_KNEE, crush_label="past the RSM knee",
        scale=scale, savepath=str(out / "damage.png"),
        suptitle=f"Our {L:.0f} mm cube, Aydin's approach — Rmax/d = {rmax_ratio:.3f}, "
                 f"deformation x{scale:.0f}\nevery red strut is a tensile crack; "
                 f"no strut can fail in compression")
    print(f"saved crack pattern to {out / 'damage.png'}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Our 200 mm compression cube, modelled and analysed Aydin's (2021) way")
    p.add_argument("--rmax", type=float, default=RMAX_RATIO, metavar="RATIO",
                   help=f"grid perturbation Rmax/d (default {RMAX_RATIO}, calibrated for THIS "
                        "cube). Pass 0 for the uniform-grid locking control")
    p.add_argument("--calibrate", action="store_true",
                   help="run Aydin's Fig. 4 loop: search Rmax/d until the mean strength is within "
                        "10%% of fc")
    p.add_argument("--no-rsm", dest="rsm", action="store_false",
                   help="drop the reduced-stiffness compression knee, leaving compression purely "
                        "linear elastic")
    p.add_argument("--realizations", type=int, default=REALIZATIONS,
                   help=f"simulations per Rmax/d (default {REALIZATIONS}); the mesh is random")
    p.add_argument("--mesh", type=float, default=MESH, help=f"grid spacing, mm (default {MESH:g})")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON:g})")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"drive speed, mm/s (default {QUASI_STATIC_RATE:g}, the cube's own)")
    p.add_argument("--strain", type=float, default=TARGET_STRAIN,
                   help=f"target axial strain (default {TARGET_STRAIN:g})")
    p.add_argument("--damping", type=float, default=0.8, help="damping ratio (default 0.8)")
    p.add_argument("--seed", type=int, default=0, help="base RNG seed (default 0)")
    p.add_argument("--truss", dest="strut_element", action="store_const", const="Truss",
                   default="corotTruss",
                   help="small-displacement Truss instead of corotTruss; expect locking, since "
                        "without geometric nonlinearity the split columns cannot lose stability")
    p.add_argument("--gf", type=float, default=GF, help=f"fracture energy, N/mm (default {GF:g})")
    a = p.parse_args()
    main(rmax_ratio=a.rmax, rsm=a.rsm, realizations=a.realizations, mesh_size=a.mesh,
         horizon=a.horizon, rate=a.rate, target_strain=a.strain, damping=a.damping,
         seed0=a.seed, strut_element=a.strut_element, Gf=a.gf, do_calibrate=a.calibrate)
