"""Single RC cantilever column: nonlinear SEISMIC time-history, lattice vs fiber beam-column.

Same specimen, materials, stiffness calibration, and lattice graph (HORIZON) as pushover.py — the
setup is imported from this package, so the ONLY new ingredients are dynamic. Both models are driven
by the SAME scaled UniformExcitation (seismic base input) in X:

  - reference: the fiber `forceBeamColumn` column, SUBDIVIDED to the lattice's row heights so its
    self-mass distributes identically (run_beamcolumn_dynamic);
  - lattice: the calibrated RC lattice (run_dynamic).

The base excitation is selectable via `--excitation` (default `sine`): a harmonic acceleration
RESONANT with the structure's fundamental period T1 (the simplest input that reliably drives the
column into the nonlinear range), or the recorded `elcentro` 1940 NS. Either way it is applied as a
UniformExcitation. The sine record is built AFTER the modal step so its period equals the computed
T1.

Consistency: lateral stiffness is calibrated (lattice strut area -> fiber K0); both carry the same
distributed self-mass PLUS a top lumped mass P/g (the axial load treated as tributary seismic
weight); both use 5% Rayleigh damping at modes 1 & 2 and a constant gravity P. The excitation scale
factor is tuned on the (cheap) fiber column to a target peak roof drift, then the SAME scaled record
drives the lattice. A period check gates the comparison. Output:
examples/output/column/column_dynamic_<excitation>.png (roof-drift & base-shear histories + hysteresis).

Run as `python examples/column/dynamic.py`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.materials import concrete_uniaxial_nonlinear, steel_uniaxial
from rclattice.opensees import run_beamcolumn_dynamic, run_dynamic, run_modal

from build import (
    _continuum_model, beamcolumn_reference, calibrate_area, continuum_dynamic, continuum_k0, rc_lattice,
)
from excitation import G, N_CYCLES, RHO, make_excitation, tune_intensity
from specimen import CORE, COVER_C, EPS, GF, GFC, H, HORIZON, MESH, OUT, P, STEEL, THK, W

TARGET_DRIFT = 2.0                          # in (~1.5% of H): moderate, well past yield
DT = 0.005                            # integration step: <= the sine record spacing (~T1/100)
#                                             so the resonant input is resolved (coarser dt under-
#                                             resolves it and triggers the hard steps that spike V)
REFERENCE_LABEL = {"beamcolumn": "fiber beam-column", "continuum": "2D continuum"}
CONTINUUM_CYCLES = 8                         # cap the sine cycles for the heavy continuum reference
#                                             (~1.5 s/step over 1536 ASDConcrete3D quads, D30)


def main(*, reference: str = "beamcolumn", horizon: float = HORIZON, regularize: bool = True,
         Gf: float = GF, Gfc: float = GFC, target_drift: float = TARGET_DRIFT, dt: float = DT,
         excitation: str = "sine", n_cycles: int = N_CYCLES, draw: bool = False) -> None:
    outdir = OUT / "dynamic" / reference
    outdir.mkdir(parents=True, exist_ok=True)
    ref_label = REFERENCE_LABEL[reference]

    self_mass = RHO * W * H * THK           # total column self-mass (rho * volume)
    top_mass = P / G                        # axial load P treated as tributary seismic mass
    nelem = int(round(H / MESH))            # subdivide the fiber column to the lattice rows
    print(f"mass: self={self_mass:.3e}, top P/g={top_mass:.3e} kip-s^2/in  (top/self = {top_mass/self_mass:.0f}x)")

    bc_materials = (concrete_uniaxial_nonlinear(CORE, 1),
                    concrete_uniaxial_nonlinear(COVER_C, 2),
                    steel_uniaxial(STEEL, 3))

    # stiffness calibration: strut area -> the SELECTED reference's K0 (mirrors the pushover). For the
    # continuum the lattice is matched to the continuum K0 so both share mass+stiffness+T1, isolating
    # the constitutive (Concrete02 vs ASDConcrete3D) hysteresis difference. The continuum is heavy, so
    # cap its sine cycles.
    if reference == "continuum":
        k_ref = continuum_k0(Gf=Gf, Gfc=Gfc)
        n_cyc = min(n_cycles, CONTINUUM_CYCLES)
    else:
        k_ref = (lambda r: r["shear"][1] / r["disp"][1])(beamcolumn_reference())
        n_cyc = n_cycles
    area, ctrl, base = calibrate_area(k_ref, horizon=horizon)
    print(f"calibrated strut area A = {area:.3f} in^2 (K0 match to {ref_label}) | horizon={horizon:g}")

    # lattice + top seismic mass baked into the model (so modal AND dynamic both see it)
    model = rc_lattice(regularize, Gf, Gfc, area, horizon=horizon)
    top_nodes = select_nodes(model, (-W, W, H - EPS, H + EPS))
    for nid in top_nodes:
        mx, my = model.masses[nid]
        model.masses[nid] = (mx + top_mass / len(top_nodes), my + top_mass / len(top_nodes))

    # period gate (cheap eigen on the lattice) — also the resonant period for the sine input
    lat_T = run_modal(model, 2)["periods"]
    print(f"lattice periods: T1={lat_T[0]:.3f}s  T2={lat_T[1]:.3f}s")

    # draw the analysis model(s) up front (before the expensive time-history)
    if draw:
        panels = [("RC lattice", model)]
        if reference == "continuum":
            panels.insert(0, ("2D continuum", _continuum_model(Gf, Gfc)[0]))
        hz = "" if abs(horizon - HORIZON) < 1e-9 else f"_h{horizon:g}"
        drawpath = outdir / f"column_dynamic_{reference}_{excitation}{hz}_model.png"
        viz.figure_model(panels, savepath=str(drawpath),
                         suptitle="RC cantilever column — analysis model")
        print(f"saved model drawing to {drawpath}")

    # base excitation (resonant sine at T1 by default; built AFTER the modal step)
    accel, dt_rec, exc_label = make_excitation(excitation, lat_T[0], n_cyc)
    pga = max(abs(a) for a in accel)
    print(f"excitation: {exc_label} | {len(accel)} pts @ {dt_rec:.4f}s, peak={pga:.3f}g, "
          f"duration={(len(accel) - 1) * dt_rec:.1f}s")

    # tune the excitation scale on the (cheap) fiber column to ~target_drift peak roof drift. The fiber
    # is the only cheap proxy (the lattice/continuum are both ~1.4 s/step); for the continuum reference
    # its K0 differs from the matched lattice/continuum, so the achieved drift is approximate — the
    # lattice-vs-continuum comparison is unaffected (both get the identical scaled record).
    def run_bc(intensity: float) -> dict:
        return run_beamcolumn_dynamic(height=H, P=P, materials=bc_materials, nelem=nelem,
                                      self_mass=self_mass, top_mass=top_mass, accel=accel,
                                      dt_record=dt_rec, scale=G * intensity, dt=dt)

    used, bc = tune_intensity(run_bc, target_drift, linear=False)
    print(f"=> scale {used:.3f}x (BC peak drift {abs(bc['peak_disp']):.3f} in) | "
          f"period match: T1_lat={lat_T[0]:.3f}s vs T1_bc={bc['periods'][0]:.3f}s")

    # lattice seismic run at the SAME scaled record (expensive)
    lat = run_dynamic(model, accel=accel, dt_record=dt_rec, scale=G * used,
                      control_node=ctrl, control_dof=1, base_nodes=base, dt=dt)

    # reference seismic run at the SAME scaled record: reuse the fiber tuning run, or run the continuum
    if reference == "continuum":
        print("running continuum seismic reference (heavy: ~1.5 s/step)...", flush=True)
        ref = continuum_dynamic(accel=accel, dt_record=dt_rec, scale=G * used, top_mass=top_mass,
                                dt=dt, Gf=Gf, Gfc=Gfc)
    else:
        ref = bc

    for name, r in ((ref_label, ref), ("RC lattice", lat)):
        print(f"{name:18s}: peak drift={abs(r['peak_disp']):.3f} in @ t={r['peak_time']:.2f}s | "
              f"residual={r['residual_disp']:+.3f} in | peakV={abs(r['peak_shear']):.1f} kip | conv={r['converged']}")

    series = [
        {"t": ref["t"], "disp": ref["disp"], "shear": ref["shear"],
         "label": ref_label, "style": {"color": "C3", "lw": 1.2}},
        {"t": lat["t"], "disp": lat["disp"], "shear": lat["shear"],
         "label": "RC lattice", "style": {"color": "C0", "lw": 1.0, "alpha": 0.85}},
    ]
    base_title = f"{exc_label} {used:.2f}x (T1~{lat_T[0]:.2f}s): RC lattice vs {ref_label}"
    hz = "" if abs(horizon - HORIZON) < 1e-9 else f"_h{horizon:g}"
    stem = f"column_dynamic_{reference}_{excitation}{hz}"
    viz.figure_timehistory(series, include_hysteresis=False,
                           savepath=str(outdir / f"{stem}.png"), title=base_title)
    # separate hysteresis loops (lattice, reference, overlaid) in addition to the overlaid histories
    viz.figure_hysteresis(list(reversed(series)),
                          savepath=str(outdir / f"{stem}_hyst.png"),
                          title=f"Hysteresis — {base_title}")
    print(f"saved to {outdir / f'{stem}.png'} (+ _hyst.png)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RC cantilever column nonlinear seismic time-history")
    parser.add_argument("--reference", choices=("beamcolumn", "continuum"), default="beamcolumn",
                        help="seismic reference: fiber beam-column (fast, default) or 2D continuum "
                             "(D30, heavy ~1.5 s/step; El Centro = full ~31 s record is ~78 min)")
    parser.add_argument("--horizon", type=float, default=HORIZON,
                        help="strut connectivity horizon * mesh_size (default 1.5; larger = more "
                             "redundant bracing against the post-peak mechanism, D31)")
    parser.add_argument("--no-regularize", dest="regularize", action="store_false",
                        help="use plain Concrete02 instead of the crack-band length-regularized law")
    parser.add_argument("--gf", type=float, default=GF, help="tensile fracture energy (kip, in)")
    parser.add_argument("--gfc", type=float, default=GFC, help="compressive fracture energy (kip, in)")
    parser.add_argument("--target-drift", type=float, default=TARGET_DRIFT, help="peak roof drift target (in)")
    parser.add_argument("--dt", type=float, default=DT, help="Newmark integration step (s)")
    parser.add_argument("--excitation", choices=("sine", "elcentro"), default="sine",
                        help="base excitation: resonant sine at T1 (default) or the El Centro record")
    parser.add_argument("--n-cycles", type=int, default=N_CYCLES, help="number of cycles for the sine excitation")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model(s) (rebar highlighted)")
    args = parser.parse_args()
    main(reference=args.reference, horizon=args.horizon, regularize=args.regularize, Gf=args.gf,
         Gfc=args.gfc, target_drift=args.target_drift, dt=args.dt, excitation=args.excitation,
         n_cycles=args.n_cycles, draw=args.draw)
