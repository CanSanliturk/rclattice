"""RC portal frame, LINEAR-ELASTIC seismic time-history: lattice vs fiber beam-column frame (linear
sibling of dynamic.py, the frame analog of examples/column/dynamic_linear.py).

The nonlinear dynamic comparison (dynamic.py) matches worse than the pushover because the hysteresis /
energy dissipation / residual drift depend on constitutive details that differ between axial struts
and a fiber section. This script takes a step back: every material is linear `Elastic`, so the two
models can agree only through their LINEAR dynamic properties — mass, lateral stiffness (hence
period), and damping. Those are matched by construction here, so the responses should track closely.

Same specimen / mass / damping / excitation plumbing as dynamic.py; the elastic setup is imported
from this package (build.py):
  - reference: the fiber `forceBeamColumn` portal frame (Elastic fibers), each member subdivided to
    the lattice row/column spacing so its self-mass distributes identically (run_beamcolumn_frame_dynamic);
  - lattice: the SAME elastic RC lattice (build_lattice_rc, Elastic concrete + rebar), strut area
    calibrated so its K0 (hence T1) matches the elastic fiber frame.

Both carry the same distributed self-mass + a per-column top lumped mass P/g, 5% modal damping at
modes 1 & 2, constant gravity P, and the SAME scaled UniformExcitation in X. The base input is
selectable via `--excitation` (default `sine`): a harmonic acceleration RESONANT with T1, or the
`elcentro` record. Because the system is linear, the response scales linearly with intensity, so the
scale is set in a single correction to hit a target peak roof drift (no iteration). A period check
gates the comparison. Output:
examples/output/frame/dynamic_linear/beamcolumn/frame_dynamic_linear_<excitation>.png.
Run as `python examples/frame/dynamic_linear.py`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.materials import concrete_uniaxial_elastic, steel_uniaxial_elastic
from rclattice.opensees import run_beamcolumn_frame_dynamic, run_dynamic, run_modal

from build import (
    beamcolumn_reference_linear, calibrate_area_linear, modal_calibration_figure, rc_lattice_linear,
)
from excitation import G, N_CYCLES, RHO, make_excitation, tune_intensity
from specimen import BEAM, COL, CORE, COVER_C, H, MESH, OUT, P, SPAN, STEEL, THK, add_axial_mass

TARGET_DRIFT = 1.0     # in — keep well inside the elastic range (no yield to chase); just sets amplitude
DT = 0.01              # Newmark step (record is at 0.02 s)


def main(*, target_drift: float = TARGET_DRIFT, dt: float = DT,
         excitation: str = "sine", n_cycles: int = N_CYCLES, draw: bool = False) -> None:
    outdir = OUT / "dynamic_linear" / "beamcolumn"
    outdir.mkdir(parents=True, exist_ok=True)

    self_mass_col = RHO * THK * COL * H               # per-column self-mass (rho * volume)
    self_mass_beam = RHO * THK * (SPAN + COL) * BEAM  # beam self-mass
    top_mass = P / G                                  # axial load P per column as tributary seismic mass
    ncol = int(round(H / MESH))                       # subdivide the fiber columns to the lattice rows
    nbeam = int(round(SPAN / MESH))                   # subdivide the fiber beam to the lattice columns
    print(f"mass: self(col)={self_mass_col:.3e}, self(beam)={self_mass_beam:.3e}, "
          f"top P/g={top_mass:.3e} kip-s^2/in/col")

    bc_materials = (concrete_uniaxial_elastic(CORE, 1),
                    concrete_uniaxial_elastic(COVER_C, 2),
                    steel_uniaxial_elastic(STEEL, 3))

    # stiffness calibration to the LINEAR elastic fiber-frame K0 (full RC-topology lattice)
    k_bc = (lambda r: r["shear"][1] / r["disp"][1])(beamcolumn_reference_linear())
    area, ctrl, base = calibrate_area_linear(k_bc)
    print(f"calibrated concrete strut area A = {area:.3f} in^2 (K0 match to {k_bc:.2f} kip/in)")

    # lattice + per-column top seismic mass baked into the model (so modal AND dynamic both see it)
    model = rc_lattice_linear(area)

    # calibration output (D35): mode shapes (elastic fiber frame vs lattice) + periods table, drawn up
    # front on the self-mass calibration basis (BEFORE the seismic top mass is added below)
    cal_caption = f"linear scalar calibration — strut area A={area:.3f} in^2 (K0 match to {k_bc:.2f} kip/in)"
    modalpath = outdir / "frame_modal_linear.png"
    t_ref, t_lat = modal_calibration_figure(reference="beamcolumn", lattice_model=model,
                                            label="elastic fiber frame", caption=cal_caption,
                                            savepath=str(modalpath), linear=True)
    print(f"modal calibration: fiber frame T={[f'{t:.4f}' for t in t_ref]} s | "
          f"lattice T={[f'{t:.4f}' for t in t_lat]} s -> saved {modalpath}")

    add_axial_mass(model, top_mass)

    lat_T = run_modal(model, 2)["periods"]
    print(f"lattice periods: T1={lat_T[0]:.3f}s  T2={lat_T[1]:.3f}s")

    if draw:
        drawpath = outdir / f"frame_dynamic_linear_{excitation}_model.png"
        viz.figure_model([("RC lattice (linear)", model)], savepath=str(drawpath),
                         suptitle="LINEAR RC portal frame — analysis model")
        print(f"saved model drawing to {drawpath}")

    # base excitation (resonant sine at T1 by default; built AFTER the modal step)
    accel, dt_rec, exc_label = make_excitation(excitation, lat_T[0], n_cycles)
    pga = max(abs(a) for a in accel)
    print(f"excitation: {exc_label} | {len(accel)} pts @ {dt_rec:.4f}s, peak={pga:.3f}g, "
          f"duration={(len(accel) - 1) * dt_rec:.1f}s")

    def run_bc(intensity: float) -> dict:
        return run_beamcolumn_frame_dynamic(height=H, span=SPAN, beam_depth=BEAM, beam_width=THK, P=P,
                                            materials=bc_materials, ncol=ncol, nbeam=nbeam,
                                            self_mass_col=self_mass_col, self_mass_beam=self_mass_beam,
                                            top_mass=top_mass, accel=accel, dt_record=dt_rec,
                                            scale=G * intensity, dt=dt)

    # linear -> response is proportional to intensity, so one trial + one exact correction hits target
    used, bc = tune_intensity(run_bc, target_drift, linear=True)
    print(f"=> scale {used:.3f}x | period match: T1_lat={lat_T[0]:.3f}s vs T1_bc={bc['periods'][0]:.3f}s")

    # lattice seismic run at the SAME scaled record
    lat = run_dynamic(model, accel=accel, dt_record=dt_rec, scale=G * used,
                      control_node=ctrl, control_dof=1, base_nodes=base, dt=dt)

    for name, r in (("fiber frame", bc), ("RC lattice", lat)):
        print(f"{name:12s}: peak drift={abs(r['peak_disp']):.3f} in @ t={r['peak_time']:.2f}s | "
              f"residual={r['residual_disp']:+.4f} in | peakV={abs(r['peak_shear']):.1f} kip | conv={r['converged']}")

    series = [
        {"t": bc["t"], "disp": bc["disp"], "shear": bc["shear"],
         "label": "elastic fiber frame", "style": {"color": "C3", "lw": 1.2}},
        {"t": lat["t"], "disp": lat["disp"], "shear": lat["shear"],
         "label": "elastic RC lattice", "style": {"color": "C0", "lw": 1.0, "alpha": 0.85}},
    ]
    base_title = f"LINEAR {exc_label} {used:.2f}x (T1~{lat_T[0]:.2f}s): RC lattice vs fiber frame"
    viz.figure_timehistory(series, include_hysteresis=False,
                           savepath=str(outdir / f"frame_dynamic_linear_{excitation}.png"), title=base_title)
    viz.figure_hysteresis(list(reversed(series)),
                          savepath=str(outdir / f"frame_dynamic_linear_{excitation}_hyst.png"),
                          title=f"Hysteresis — {base_title}")
    print(f"saved to {outdir / f'frame_dynamic_linear_{excitation}.png'} (+ _hyst.png)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="LINEAR RC portal frame seismic time-history")
    parser.add_argument("--target-drift", type=float, default=TARGET_DRIFT, help="peak roof drift target (in)")
    parser.add_argument("--dt", type=float, default=DT, help="Newmark integration step (s)")
    parser.add_argument("--excitation", choices=("sine", "elcentro"), default="sine",
                        help="base excitation: resonant sine at T1 (default) or the El Centro record")
    parser.add_argument("--n-cycles", type=int, default=N_CYCLES, help="number of cycles for the sine excitation")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model (rebar highlighted)")
    args = parser.parse_args()
    main(target_drift=args.target_drift, dt=args.dt, excitation=args.excitation, n_cycles=args.n_cycles,
         draw=args.draw)
