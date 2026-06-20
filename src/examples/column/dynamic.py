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

from build import beamcolumn_reference, calibrate_area, rc_lattice
from excitation import G, N_CYCLES, RHO, make_excitation, tune_intensity
from specimen import CORE, COVER_C, EPS, GF, GFC, H, MESH, OUT, P, STEEL, THK, W

TARGET_DRIFT = 2.0                          # in (~1.5% of H): moderate, well past yield
DT = 0.01 #0.005                            # integration step: <= the sine record spacing (~T1/100)
#                                             so the resonant input is resolved (coarser dt under-
#                                             resolves it and triggers the hard steps that spike V)


def main(*, regularize: bool = True, Gf: float = GF, Gfc: float = GFC,
         target_drift: float = TARGET_DRIFT, dt: float = DT,
         excitation: str = "sine", n_cycles: int = N_CYCLES) -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    self_mass = RHO * W * H * THK           # total column self-mass (rho * volume)
    top_mass = P / G                        # axial load P treated as tributary seismic mass
    nelem = int(round(H / MESH))            # subdivide the fiber column to the lattice rows
    print(f"mass: self={self_mass:.3e}, top P/g={top_mass:.3e} kip-s^2/in  (top/self = {top_mass/self_mass:.0f}x)")

    bc_materials = (concrete_uniaxial_nonlinear(CORE, 1),
                    concrete_uniaxial_nonlinear(COVER_C, 2),
                    steel_uniaxial(STEEL, 3))

    # stiffness calibration (identical to the pushover): strut area -> fiber column K0
    k_bc = (lambda r: r["shear"][1] / r["disp"][1])(beamcolumn_reference())
    area, ctrl, base = calibrate_area(k_bc)
    print(f"calibrated strut area A = {area:.3f} in^2")

    # lattice + top seismic mass baked into the model (so modal AND dynamic both see it)
    model = rc_lattice(regularize, Gf, Gfc, area)
    top_nodes = select_nodes(model, (-W, W, H - EPS, H + EPS))
    for nid in top_nodes:
        mx, my = model.masses[nid]
        model.masses[nid] = (mx + top_mass / len(top_nodes), my + top_mass / len(top_nodes))

    # period gate (cheap eigen on the lattice) — also the resonant period for the sine input
    lat_T = run_modal(model, 2)["periods"]
    print(f"lattice periods: T1={lat_T[0]:.3f}s  T2={lat_T[1]:.3f}s")

    # base excitation (resonant sine at T1 by default; built AFTER the modal step)
    accel, dt_rec, exc_label = make_excitation(excitation, lat_T[0], n_cycles)
    pga = max(abs(a) for a in accel)
    print(f"excitation: {exc_label} | {len(accel)} pts @ {dt_rec:.4f}s, peak={pga:.3f}g, "
          f"duration={(len(accel) - 1) * dt_rec:.1f}s")

    # tune the excitation scale on the (cheap) fiber column to ~target_drift peak roof drift
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

    for name, r in (("BC ", bc), ("lat", lat)):
        print(f"{name}: peak drift={abs(r['peak_disp']):.3f} in @ t={r['peak_time']:.2f}s | "
              f"residual={r['residual_disp']:+.3f} in | peakV={abs(r['peak_shear']):.1f} kip | conv={r['converged']}")

    series = [
        {"t": bc["t"], "disp": bc["disp"], "shear": bc["shear"],
         "label": "fiber beam-column", "style": {"color": "C3", "lw": 1.2}},
        {"t": lat["t"], "disp": lat["disp"], "shear": lat["shear"],
         "label": "RC lattice", "style": {"color": "C0", "lw": 1.0, "alpha": 0.85}},
    ]
    base_title = f"{exc_label} {used:.2f}x (T1~{lat_T[0]:.2f}s): RC lattice vs fiber beam-column"
    viz.figure_timehistory(series, include_hysteresis=False,
                           savepath=str(OUT / f"column_dynamic_{excitation}.png"), title=base_title)
    # separate hysteresis loops (lattice, fiber, overlaid) in addition to the overlaid histories
    viz.figure_hysteresis(list(reversed(series)),
                          savepath=str(OUT / f"column_dynamic_{excitation}_hyst.png"),
                          title=f"Hysteresis — {base_title}")
    print(f"saved to {OUT / f'column_dynamic_{excitation}.png'} (+ _hyst.png)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RC cantilever column nonlinear seismic time-history")
    parser.add_argument("--no-regularize", dest="regularize", action="store_false",
                        help="use plain Concrete02 instead of the crack-band length-regularized law")
    parser.add_argument("--gf", type=float, default=GF, help="tensile fracture energy (kip, in)")
    parser.add_argument("--gfc", type=float, default=GFC, help="compressive fracture energy (kip, in)")
    parser.add_argument("--target-drift", type=float, default=TARGET_DRIFT, help="peak roof drift target (in)")
    parser.add_argument("--dt", type=float, default=DT, help="Newmark integration step (s)")
    parser.add_argument("--excitation", choices=("sine", "elcentro"), default="sine",
                        help="base excitation: resonant sine at T1 (default) or the El Centro record")
    parser.add_argument("--n-cycles", type=int, default=N_CYCLES, help="number of cycles for the sine excitation")
    args = parser.parse_args()
    main(regularize=args.regularize, Gf=args.gf, Gfc=args.gfc, target_drift=args.target_drift,
         dt=args.dt, excitation=args.excitation, n_cycles=args.n_cycles)
