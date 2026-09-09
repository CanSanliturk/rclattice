"""200 mm concrete cube, uniaxial COMPRESSION — simplified Aydin-calibrated lattice (D56).

The simplest model that still exercises Aydin's calibration in compression: one Concrete02 grade at
fc = 15.1 MPa, standard `Truss` struts, rollers under the base with a single pinned node, and a very
slow dynamic-relaxation drive so the recorded curve is quasi-static while mass and damping regularize
the softening instability a static path-follower cannot cross.

VERIFICATION WITHOUT A REFERENCE MODEL. The check is closed-form: the measured initial tangent must
match `build.expected_modulus`. That target is deliberately NOT the grade's E — Aydin's Eq. 2.1
balances energies with the transverse strain restrained, so it pins the CONFINED modulus E/(1-nu^2),
while a cube between smooth platens expands freely and reads `E_confined*(1 - nu_lattice^2)`. With
nu_lattice ~ 0.41 at horizon 1.5 that is about 0.87*E (D53).

WHY THE PEAK IS NOT fc*A_face. Printed as a four-line accounting: the lattice has MORE material than
the continuum (crossing struts ~1.8x the face area), compatibility keeps the diagonals far from fc
when the verticals peak, and transverse cracking then sheds the inclined path entirely. The section
is short of a LOAD PATH, not of material (D55).

`--rigid-platen` ties the whole top face to one vertical DOF (`ops.equalDOF`) so it must stay flat,
against the default where each top node is driven independently.

Every run writes its OWN timestamped directory under examples/output/compression_cube/runs/
(`command.txt`, `data.json`, `stress_strain.png`, `force_deformation.png`, `loadpath.png`,
`damage.png`), so a sweep across platen / solver / horizon keeps every variant instead of only the
last one. Units: N, mm.
"""

from __future__ import annotations

import json
import time

import numpy as np

from rclattice import viz
from rclattice.opensees import run_pushover, run_pushover_dynamic

from build import (
    calibrate, capacity_accounting, cube_lattice, expected_modulus, scaled_grade, secant_factor,
    strength_correction_factor, vertical_strut_capacity,
)
from specimen import (
    A_FACE, CONCRETE, EPS_CRACK, GF, GFC_FACTOR, HORIZON, L, MESH,
    QUASI_STATIC_RATE, RESIDUAL_RATIO, TARGET_STRAIN, axial_loads, base_nodes, control_node,
    run_dir, top_nodes,
)

# Axial strain below which the response is still uncracked, used to measure the initial tangent.
# The first transverse struts crack at an axial strain of about eps_crack/nu, so a window capped
# well under that stays on the elastic branch.
LINEAR_FIT_STRAIN = 2.0e-4


def cut_element_groups(model, y_cut: float, tol: float = 1e-6) -> dict:
    """Struts crossing the horizontal cut at `y_cut`, grouped by orientation, as `element_groups`.

    For every strut spanning the cut, record the VERTICAL force it exerts on its LOWER node.
    `eleForce` returns a 2-node truss's global nodal forces as [Fxi, Fyi, Fxj, Fyj], so index 1 is
    Fy at the element's first node and 3 is Fy at its second — pick whichever sits below.

    Summed over every crossing strut this must reconcile with the axial reaction, which is what
    makes the split trustworthy; `main` prints that reconciliation with an explicit verdict.
    Horizontal struts never cross a horizontal cut, so the split is vertical vs diagonal — the
    horizontals matter indirectly, by cracking and releasing the restraint the diagonals need.
    """
    groups: dict[str, list] = {"vertical struts": [], "diagonal struts": []}
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        if (ya - y_cut) * (yb - y_cut) >= 0.0:
            continue
        dof = 1 if ya < yb else 3
        label = "vertical struts" if abs(xb - xa) < tol else "diagonal struts"
        groups[label].append((e.id, dof, 1.0))
    return groups


def initial_tangent(strain, stress) -> tuple[float, float]:
    """Least-squares secant through the origin over the uncracked window, and its upper strain."""
    x, y = np.asarray(strain), np.asarray(stress)
    m = (x > 0.0) & (x <= LINEAR_FIT_STRAIN)
    if m.sum() < 2:
        m = (x > 0.0) & (x <= max(x.max() * 0.05, 1e-12))
    if m.sum() < 2:
        return float("nan"), 0.0
    return float(np.dot(x[m], y[m]) / np.dot(x[m], x[m])), float(x[m].max())


def main(*, mesh_size: float = MESH, horizon: float = HORIZON,
         target_strain: float = TARGET_STRAIN, rate: float = QUASI_STATIC_RATE,
         rigid_platen: bool = False, damping: float = 0.8, solver: str = "dynamic",
         steps: int = 500, algorithm: "tuple[str, ...]" | None = None,
         max_iter: int = 100, draw: bool = False, peak_correction: str = "none") -> None:
    if peak_correction not in ("none", "strength"):
        raise ValueError(f"peak_correction must be none/strength, got {peak_correction!r}")
    tag = solver + ("_rigid" if rigid_platen else "_free")
    if peak_correction != "none":
        tag += "_fcscaled"
    if horizon != HORIZON:
        tag += f"_h{horizon:g}"
    if mesh_size != MESH:
        tag += f"_m{mesh_size:g}"
    out = run_dir("pushover", tag=tag)
    print(f"output directory: {out}")

    cal = calibrate(mesh_size=mesh_size, horizon=horizon)
    e_expected = expected_modulus(cal)
    print(f"Aydin energy balance: A_t = {cal.area:,.1f} mm^2 "
          f"= {cal.area / (mesh_size * 200.0):.4f} * (thickness * mesh)")
    print(f"  horizon {horizon}: {cal.n_nodes} nodes, {cal.n_struts} struts")
    print(f"  lattice Poisson ratio  nu_eff = {cal.nu_effective:.4f}   (concrete: {CONCRETE.nu:.2f});"
          f"  cubic anisotropy {cal.cubic_anisotropy:.3f} "
          f"({'isotropic' if abs(cal.cubic_anisotropy - 1) < 0.05 else 'NOT isotropic'})")
    print(f"  -> the balance pins the CONFINED modulus E/(1-nu^2); an unconfined cube should "
          f"measure {e_expected:,.0f} MPa = {e_expected / CONCRETE.E:.3f}*E")

    # Build once at the physical grade to MEASURE the section, derive the correction from it, then
    # rebuild. The correction changes neither the area nor the geometry, so the second lattice is
    # the same structure with a different strut strength — the factor is read off the strut list,
    # not fitted to a curve.
    y_cut = L / 2.0 + mesh_size / 2.0
    model = cube_lattice(cal.area, mesh_size=mesh_size, horizon=horizon)
    fc_scale = strength_correction_factor(model, y_cut) if peak_correction == "strength" else 1.0
    if fc_scale != 1.0:
        model = cube_lattice(cal.area, mesh_size=mesh_size, horizon=horizon, fc_scale=fc_scale)
    grade = scaled_grade(fc_scale)
    eps_crush = grade.epsc0                 # the STRUT's peak-compression strain, scaled with fc

    platen = "RIGID (top face tied to one vertical DOF)" if rigid_platen else "free (each top node driven)"
    print(f"\nlattice: {len(model.nodes)} nodes, {len(model.elements)} Concrete02 Truss struts")
    print(f"  fc = {CONCRETE.fc} MPa, E = {CONCRETE.E:,.0f} MPa, epsc0 = {CONCRETE.epsc0:.4e} "
          f"(= 2fc/E), ft = {CONCRETE.ft} MPa")
    print(f"  Gf = {GF:g} N/mm, Gfc = {GFC_FACTOR * GF:g} N/mm, residual_ratio = {RESIDUAL_RATIO:g}"
          f"   [these set the post-peak branch — results, not just settings]")
    print(f"  upper platen: {platen}")

    n_vert, cap_force, cap_stress = vertical_strut_capacity(model, y_cut, fc_strut=grade.fc)
    print(f"  section at y = {y_cut:.0f} mm: {n_vert} vertical struts = "
          f"{n_vert * cal.area / A_FACE:.3f} of the gross face -> verticals alone at "
          f"{'fc_strut' if fc_scale != 1.0 else 'fc'} give {cap_stress:.2f} MPa")
    if fc_scale != 1.0:
        print(f"\n  PEAK CORRECTION (D71, --peak-correction strength)")
        print(f"    strut fc     {CONCRETE.fc:6.2f} -> {grade.fc:6.2f} MPa   "
              f"x{fc_scale:.4f} = 1 / {1.0 / fc_scale:.4f}, the verticals' share of the gross face")
        print(f"    strut epsc0  {CONCRETE.epsc0:.4e} -> {grade.epsc0:.4e}   FORCED: Concrete02's "
              f"tangent is 2fc/epsc0, so epsc0 must scale with fc to leave E at {CONCRETE.E:,.0f} MPa")
        print(f"    strut fcu    {CONCRETE.fcu:6.2f} -> {grade.fcu:6.2f} MPa   (residual keeps its "
              f"ratio to fc);  ft UNCHANGED at {CONCRETE.ft:g} MPa — a different knob")
        print(f"    area, and therefore K0, is UNTOUCHED — that is the point of correcting strength "
              f"here rather than area")
        print(f"    PRICE: each strut now peaks at {fc_scale:.2f}*epsc0, so the lattice peak STRESS "
              f"is corrected and its peak STRAIN moves by the same factor.")
        print(f"           Concrete02's parabola ties (E, peak stress, peak strain) together; one "
              f"scaling cannot fix the stress and hold the strain.")

    target = -target_strain * L
    ctrl, base, top = control_node(model), base_nodes(model), top_nodes(model)
    # RIGID platen: tie every other top node's vertical DOF to the control node and drive only that
    # node, so the face is forced to stay flat. Default: drive each top node independently.
    if rigid_platen:
        tie = [(ctrl, nid, 2) for nid in top if nid != ctrl]
        drive = [ctrl]
    else:
        tie, drive = None, top

    # Per-solver default: the static pass needs the constant ELASTIC tangent to avoid re-forming a
    # cracked (singular/indefinite) one; the dynamic pass does not, because mass and damping keep the
    # effective tangent well-conditioned. `--algorithm` overrides either.
    if algorithm is None:
        algorithm = ("Newton",) if solver == "dynamic" else ("ModifiedNewton", "-initial")
    print(f"  algorithm: {' '.join(algorithm)}  (retry ladder: -> KrylovNewton -> NewtonLineSearch)"
          f"  max_iter = {max_iter} per step")

    groups_spec = cut_element_groups(model, y_cut)
    t0 = time.time()
    if solver == "dynamic":
        res = run_pushover_dynamic(model, control_node=ctrl, control_dof=2, target=target,
                                   drive_nodes=drive, base_nodes=base, rate=rate,
                                   steps_per_period=30, damping_ratio=damping, capture=True,
                                   element_groups=groups_spec, equal_dof=tie,
                                   algorithm=algorithm, max_iter=max_iter)
    else:
        # Static path-follower: a uniform downward reference traction on the top face, scaled by
        # DisplacementControl on the control node. `ModifiedNewton -initial` iterates on the constant
        # ELASTIC tangent, so it never re-forms the singular/negative tangent a cracked strut band
        # produces — the robust choice for a softening pass (D44). It is still expected to stall at
        # the softening instability, which is exactly what the dynamic solver exists to get past.
        res = run_pushover(model, lateral_loads=axial_loads(model), control_node=ctrl,
                           control_dof=2, dU=target / steps, target=target, base_nodes=base,
                           algorithm=algorithm, element_groups=groups_spec,
                           capture=True, equal_dof=tie, max_iter=max_iter)
        res.setdefault("rate", 0.0)
        res.setdefault("T1", float("nan"))
    elapsed = time.time() - t0
    if not res["disp"]:
        raise RuntimeError("compression pushover produced no steps")

    # Compression is a negative control displacement and a negative reaction sum; report magnitudes.
    strain = [abs(u) / L for u in res["disp"]]
    stress = [abs(s) / A_FACE for s in res["shear"]]

    e0, fit_max = initial_tangent(strain, stress)
    sf = secant_factor(fit_max, epsc0=grade.epsc0)
    e0_corr = e0 / sf
    peak = max(stress)
    at = strain[stress.index(peak)]
    how = (f"T1 = {res['T1']:.4g} s, rate = {res['rate']:g} mm/s" if solver == "dynamic"
           else f"DisplacementControl, dU = {abs(target / steps):.4g} mm")
    print(f"\ntraced to {strain[-1]:.4%} axial strain of {target_strain:.2%} requested in "
          f"{elapsed:.1f}s  [{solver}: {how}]  converged={res['converged']}")
    if not res["converged"]:
        print(f"  ** the {solver} solver STALLED before the target — everything below describes only "
              f"the portion actually traced **")
    print(f"  initial tangent        {e0_corr:9,.0f} MPa   vs expected {e_expected:,.0f} MPa"
          f"     -> ratio {e0_corr / e_expected:.4f}   [quasi-static + calibration check]")
    print(f"    (raw fit {e0:,.0f} MPa over eps <= {fit_max:.1e}, divided by the Concrete02 secant "
          f"factor {sf:.4f})")
    print(f"                                    vs grade E  {CONCRETE.E:,.0f} MPa"
          f"     -> ratio {e0_corr / CONCRETE.E:.4f}   [the unconfined-vs-confined gap]")
    # A stalled run whose maximum sits at (or next to) its last step never turned over — the value
    # is where the solver died, not a capacity, and must not be read as one.
    still_rising = (not res["converged"]) and stress.index(peak) >= len(stress) - 2
    peak_lbl = "max before stall" if still_rising else "peak stress"
    print(f"  {peak_lbl:<21} {peak:9.2f} MPa   vs fc = {CONCRETE.fc:.1f} MPa"
          f"          -> ratio {peak / CONCRETE.fc:.3f}")
    if still_rising:
        print(f"    ** the curve was STILL ASCENDING when the solver stalled — this is a lower "
              f"bound on the capacity, not the capacity **")
    print(f"  strain there           {at:9.2e}       vs strut epsc0 = {grade.epsc0:.2e}"
          f" -> ratio {at / grade.epsc0:.3f}")
    if fc_scale != 1.0:
        print(f"                                    vs material epsc0 {CONCRETE.epsc0:.2e}"
              f" -> ratio {at / CONCRETE.epsc0:.3f}   *** the corrected peak arrives this much "
              f"LATE — the price named above ***")
    if stress[-1] < peak:
        print(f"  post-peak degradation  {(1 - stress[-1] / peak) * 100:8.1f}% by "
              f"{strain[-1]:.3%} strain")

    # Did the platen constraint do anything? Measure how far from flat the top face ends up. With
    # --rigid-platen the tie forces this to zero by construction; without it each top node is given
    # the SAME prescribed ramp, so it should also be flat — the number says whether the two
    # idealizations actually differ on this specimen or merely differ in how they are imposed.
    uy_top = [res["disps_final"][n][1] for n in top]
    spread = max(uy_top) - min(uy_top)
    print(f"  top-face flatness      spread of uy = {spread:.3e} mm over {len(top)} nodes "
          f"({abs(spread / res['disp'][-1]):.2e} of the imposed shortening)")

    els, eps_final = viz.strut_strains(model, res["disps_final"])
    ncr = int((eps_final >= EPS_CRACK).sum())
    ncu = int((eps_final <= -eps_crush).sum())
    print(f"  final damage           {ncr}/{len(els)} struts cracked ({ncr / len(els):.1%}), "
          f"{ncu} crushed ({ncu / len(els):.1%})")

    # --- who carries the load -------------------------------------------------------------------
    gv = np.asarray(res["groups"]["vertical struts"])
    gd = np.asarray(res["groups"]["diagonal struts"])
    total_cut = gv + gd
    axial = np.abs(np.asarray(res["shear"]))
    ipk = int(np.argmax(axial))
    mask = axial > 0.05 * axial[ipk]
    recon = float(np.max(np.abs(np.abs(total_cut[mask]) - axial[mask]) / axial[mask]))
    eps_split = EPS_CRACK / max(cal.nu_effective, 1e-9)
    pre = np.flatnonzero(np.asarray(strain) < eps_split)
    iel = int(pre[-1]) if len(pre) else 1
    print(f"\n  load path across y = {y_cut:.0f} mm   (cut vs reaction: max error {recon:.2%} — "
          f"{'ok' if recon < 0.02 else 'CHECK FAILED'})")
    print(f"    transverse struts start cracking at axial strain eps_cr/nu = {eps_split:.2e} "
          f"= {eps_split / grade.epsc0:.1%} of the strut epsc0"
          + (f" (was {eps_split / CONCRETE.epsc0:.1%} before the correction — scaling fc makes "
             f"splitting RELATIVELY earlier, not later)" if fc_scale != 1.0 else ""))
    for label, i in (("pre-crack", iel), ("at peak", ipk), ("end of run", len(axial) - 1)):
        tot = total_cut[i] if abs(total_cut[i]) > 1e-9 else 1e-9
        print(f"    {label:<11} vertical {abs(gv[i]) / 1e3:8.1f} kN ({gv[i] / tot:6.1%})   "
              f"diagonal {abs(gd[i]) / 1e3:8.1f} kN ({gd[i] / tot:6.1%})")
    f_el, f_pk = gd[iel] / total_cut[iel], gd[ipk] / total_cut[ipk]

    # --- why the capacity is not fc*A_face ------------------------------------------------------
    acc = capacity_accounting(model, y_cut, cal.nu_effective, fc_strut=grade.fc)
    fc_lbl = "fc_strut" if fc_scale != 1.0 else "fc"
    measured = peak * A_FACE
    print(f"\n  axial capacity accounting (why it is not fc*A_face = {acc['continuum'] / 1e3:,.0f} kN)")
    rows = (
        ("continuum: fc x gross face", acc["continuum"], "material fills the cut plane"),
        (f"every crossing strut at {fc_lbl}", acc["all_at_fc"],
         f"crossing struts total {acc['area_crossing'] / A_FACE:.2f}x the face area"),
        ("compatibility-limited at epsc0", acc["compatible"],
         "diagonals strain (1-nu)/2 of the verticals, so are far from fc when they peak"),
        (f"verticals alone at {fc_lbl}", acc["verticals"],
         f"after splitting sheds the inclined path ({acc['area_vertical'] / A_FACE:.3f} of the face)"),
        ("MEASURED peak", measured, "the analysis result"),
    )
    for label, val, note in rows:
        print(f"      {label:<32} {val / 1e3:8.0f} kN  ({val / acc['continuum']:5.2f} x continuum)"
              f"   {note}")
    lost = acc["compatible"] - acc["verticals"]
    print(f"    -> shedding the inclined path costs {lost / 1e3:,.0f} kN "
          f"({lost / acc['continuum']:.0%} of the continuum capacity); the diagonals carry "
          f"{f_el:.0%} before cracking and {f_pk:.0%} at peak. The section is NOT short of "
          f"material — it is short of a LOAD PATH.")
    if measured < 0.95 * acc["verticals"]:
        print(f"       The measured peak is {(acc['verticals'] - measured) / 1e3:,.0f} kN below even "
              f"the verticals-alone bound: it arrives at {at / grade.epsc0:.2f}*(strut epsc0), so "
              f"the verticals never reach their fc either and failure is progressive rather than "
              f"sectional.")

    # --- outputs --------------------------------------------------------------------------------
    datapath = out / "data.json"
    datapath.write_text(json.dumps({
        "strain": strain, "stress_MPa": stress, "mesh": mesh_size, "horizon": horizon,
        "rigid_platen": rigid_platen, "solver": solver, "algorithm": list(algorithm),
        "max_iter": max_iter, "rate": res["rate"], "T1": res["T1"],
        "area": cal.area,
        "nu_effective": cal.nu_effective, "E_measured": e0_corr, "E_expected": e_expected,
        "peak_MPa": peak, "strain_at_peak": at, "still_rising_at_stall": still_rising,
        "cracked": ncr, "crushed": ncu,
        "n_struts": len(els), "top_face_uy_spread": spread, "cut_y": y_cut, "cut_vertical_N": list(gv),
        "cut_diagonal_N": list(gd), "cut_reconciliation": recon,
        "capacity": {k: v for k, v in acc.items()},
        "peak_correction": peak_correction, "fc_scale": fc_scale,
        "fc_material": CONCRETE.fc, "fc_strut": grade.fc,
        "epsc0_material": CONCRETE.epsc0, "epsc0_strut": grade.epsc0,
        "peak_over_fc": peak / CONCRETE.fc,
        "converged": res["converged"], "elapsed_s": elapsed,
    }))
    print(f"\nsaved raw response to {datapath}")

    plat_lbl = ("rigid platen" if rigid_platen else "free platen") + f", {solver}"
    corr_lbl = f", fc x{fc_scale:.3f}" if fc_scale != 1.0 else ""
    # Two short lines rather than one long one: matplotlib does not wrap a title, and a single
    # sentence here runs past the axes on both of the figures that carry it.
    corr_title = (f"\nstrut fc scaled x{fc_scale:.3f} (D71)"
                  f"\npeak stress corrected onto fc; peak strain moves with it"
                  if fc_scale != 1.0 else "")
    span = [0.0, max(strain)]
    viz.figure_pushover(
        [{"disp": strain, "shear": stress,
          "label": f"lattice, {plat_lbl}{corr_lbl} (peak {peak:.1f} MPa)",
          "style": {"color": "C0", "lw": 2}},
         {"disp": span, "shear": [CONCRETE.fc] * 2,
          "label": (f"material fc — TARGET ({CONCRETE.fc:.1f} MPa)" if fc_scale != 1.0
                    else f"strut grade fc ({CONCRETE.fc:.1f} MPa)"),
          "style": {"color": "C3", "ls": "--", "lw": 1.5}},
         *([{"disp": span, "shear": [grade.fc] * 2,
             "label": f"strut fc, scaled x{fc_scale:.3f} ({grade.fc:.1f} MPa)",
             "style": {"color": "C4", "ls": "--", "lw": 1.2}}] if fc_scale != 1.0 else []),
         {"disp": [grade.epsc0] * 2, "shear": [0.0, peak],
          "label": f"strut epsc0 ({grade.epsc0:.2e})",
          "style": {"color": "0.55", "ls": ":", "lw": 1.5}}],
        savepath=str(out / "stress_strain.png"), xlabel="axial strain (-)", ylabel="axial stress (MPa)",
        title=f"200 mm cube, uniaxial compression, smooth platens\n"
              f"Aydin-calibrated lattice, fc = {CONCRETE.fc} MPa ({plat_lbl}, horizon={horizon:g})"
              + corr_title)
    print(f"saved stress-strain curve to {out / 'stress_strain.png'}")

    short_mm = [abs(u) for u in res["disp"]]
    force_kn = [abs(s) / 1e3 for s in res["shear"]]
    fd_span = [0.0, max(short_mm)]
    viz.figure_pushover(
        [{"disp": short_mm, "shear": force_kn,
          "label": f"lattice, {plat_lbl}{corr_lbl} (peak {max(force_kn):,.0f} kN)",
          "style": {"color": "C0", "lw": 2}},
         {"disp": fd_span, "shear": [acc["continuum"] / 1e3] * 2,
          "label": f"gross section at fc ({acc['continuum'] / 1e3:,.0f} kN)",
          "style": {"color": "C3", "ls": "--", "lw": 1.5}},
         {"disp": fd_span, "shear": [cap_force / 1e3] * 2,
          "label": (f"vertical struts alone at strut fc ({cap_force / 1e3:,.0f} kN)" if fc_scale != 1.0
                    else f"vertical struts alone at fc ({cap_force / 1e3:,.0f} kN)"),
          "style": {"color": "C2", "ls": "-.", "lw": 1.5}}],
        savepath=str(out / "force_deformation.png"),
        xlabel="axial shortening (mm)", ylabel="axial force (kN)",
        title=f"200 mm cube, uniaxial compression — force-deformation\n"
              f"({plat_lbl}, fc = {CONCRETE.fc} MPa, horizon={horizon:g}, mesh={mesh_size:g} mm)"
              + corr_title)
    print(f"saved force-deformation curve to {out / 'force_deformation.png'}")

    viz.figure_pushover(
        [{"disp": strain, "shear": list(np.abs(total_cut) / 1e3), "label": "total across the cut",
          "style": {"color": "0.35", "lw": 2}},
         {"disp": strain, "shear": list(np.abs(gv) / 1e3), "label": "vertical struts",
          "style": {"color": "C0", "lw": 2}},
         {"disp": strain, "shear": list(np.abs(gd) / 1e3), "label": "diagonal struts",
          "style": {"color": "C1", "lw": 2}},
         {"disp": [eps_split] * 2, "shear": [0.0, float(np.max(np.abs(total_cut))) / 1e3],
          "label": f"transverse cracking begins (eps_cr/nu = {eps_split:.1e})",
          "style": {"color": "C3", "ls": ":", "lw": 1.5}}],
        savepath=str(out / "loadpath.png"), xlabel="axial strain (-)",
        ylabel="axial force at the cut (kN)",
        title=f"Who carries the load, and when they stop (cut at y = {y_cut:.0f} mm, {plat_lbl})\n"
              f"diagonals shed their share once transverse cracking frees the lateral restraint")
    print(f"saved load-path decomposition to {out / 'loadpath.png'}")

    scale = viz.autoscale_factor(model, res["disps_final"], frac=0.10)
    viz.figure_damage(
        [(f"at peak load ({peak:.1f} MPa)", model, res["disps_peak"]),
         (f"end of run ({strain[-1]:.2%} strain)", model, res["disps_final"])],
        eps_crack=EPS_CRACK, eps_crush=eps_crush, scale=scale,
        savepath=str(out / "damage.png"),
        suptitle=f"200 mm cube in compression (fc = {CONCRETE.fc} MPa, {plat_lbl}{corr_lbl}, "
                 f"horizon={horizon:g}) — deformation x{scale:.0f}\ntransverse tension from Poisson "
                 f"expansion is what cracks: the lattice fails by splitting")
    print(f"saved damage pattern to {out / 'damage.png'}")

    if draw:
        viz.figure_model([("compression cube lattice", model)],
                         savepath=str(out / "model.png"),
                         suptitle="200 mm compression cube — analysis model")
        print(f"saved model drawing to {out / 'model.png'}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="200 mm cube uniaxial compression, Aydin-calibrated lattice")
    p.add_argument("--peak-correction", choices=("none", "strength"), default="none",
                   help="'strength' (D71) scales the strut grade's fc so the LATTICE peaks at the "
                        "material fc, leaving the calibrated area — and hence K0 — alone, the "
                        "compression twin of the tension cube's peak_correction. The factor is "
                        "DERIVED from the section (1 / the verticals' share of the gross face, "
                        "~1.58 at horizon 1.5), not fitted. It corrects the peak STRESS only: "
                        "epsc0 must scale with fc to hold the tangent, so the peak strain moves by "
                        "the same factor, and the splitting mechanism is unchanged")
    p.add_argument("--rigid-platen", action="store_true",
                   help="tie every top-face node to one vertical DOF (ops.equalDOF) so the upper "
                        "face must stay flat and move down by the same amount, and drive only the "
                        "control node. Default off: each top node is driven independently")
    p.add_argument("--mesh", type=float, default=MESH,
                   help=f"grid spacing, mm (default {MESH:.0f}). Must divide {L / 2:.0f} so a node "
                        f"lands at the face centre")
    p.add_argument("--horizon", type=float, default=HORIZON,
                   help=f"strut horizon in grid spacings (default {HORIZON}, 8 neighbours; 3.01 "
                        "gives ~28 — Aydin's recommendation)")
    p.add_argument("--strain", type=float, default=TARGET_STRAIN,
                   help=f"target axial strain (default {TARGET_STRAIN})")
    p.add_argument("--rate", type=float, default=QUASI_STATIC_RATE,
                   help=f"drive speed in mm/s (default {QUASI_STATIC_RATE}, very slow). The response "
                        "is rate-independent; lowering it shrinks the inertial/damping share of the "
                        "recorded reaction, which the reconciliation check reports")
    p.add_argument("--solver", choices=("dynamic", "static"), default="dynamic",
                   help="'dynamic' (default) = dynamic relaxation, which rides through the softening "
                        "instability; 'static' = DisplacementControl with ModifiedNewton -initial on "
                        "a uniform top-face traction. The static one is cheaper and needs no rate "
                        "calibration, but is expected to stall at the peak (D38/D44/D46)")
    # ONE quoted string, split here — not nargs="+". OpenSees algorithm options start with a dash
    # ("-initial", "-type Bisection"), and argparse reads those as flags of its own, so
    # `--algorithm ModifiedNewton -initial` fails outright instead of doing what it plainly means.
    p.add_argument("--algorithm", type=str, default=None, metavar="SPEC",
                   help='solution algorithm as one quoted OpenSees spec, e.g. --algorithm '
                        '"KrylovNewton" or --algorithm "ModifiedNewton -initial" or --algorithm '
                        '"NewtonLineSearch -type Bisection". Default: Newton for --solver dynamic, '
                        '"ModifiedNewton -initial" for static. Whatever is chosen is the PRIMARY; '
                        'both runners still fall back to KrylovNewton then NewtonLineSearch on a '
                        'failed step')
    p.add_argument("--max-iter", type=int, default=100,
                   help="iterations allowed per step by the NormDispIncr test (default 100). Raising "
                        "it helps a step that is merely converging slowly; it cannot help one whose "
                        "tangent has gone singular, which is what ends this analysis (D58)")
    p.add_argument("--steps", type=int, default=500,
                   help="static solver: displacement-control steps to the target (default 500)")
    p.add_argument("--damping", type=float, default=0.8,
                   help="dynamic solver: damping ratio (default 0.8)")
    p.add_argument("--draw", action="store_true", help="also save a drawing of the analysis model")
    a = p.parse_args()
    main(peak_correction=a.peak_correction, mesh_size=a.mesh, horizon=a.horizon,
         target_strain=a.strain, rate=a.rate,
         rigid_platen=a.rigid_platen, damping=a.damping, solver=a.solver, steps=a.steps,
         algorithm=tuple(a.algorithm.split()) if a.algorithm else None, max_iter=a.max_iter,
         draw=a.draw)
