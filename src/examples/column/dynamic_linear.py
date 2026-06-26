"""Single RC cantilever column, LINEAR-ELASTIC seismic time-history: lattice vs a selectable
reference (linear sibling of dynamic.py).

The nonlinear dynamic comparison (dynamic.py) matched far worse than the pushover because the
hysteresis / energy dissipation / residual drift depend on constitutive details that differ between
axial struts and a fiber section. This script takes a step back: every material is linear `Elastic`,
so the two models can agree only through their LINEAR dynamic properties — mass, lateral stiffness
(hence period), and damping. Those are matched by construction here, so the responses should track
closely.

  - reference beamcolumn: a SINGLE-element fiber `forceBeamColumn` (Elastic fibers) carrying its
    self-mass as distributed element mass;
  - reference continuum: the 2D plane-stress continuum with ElasticIsotropic+PlaneStress quads and
    elastic rebar — the like-for-like elastic linear analog of the nonlinear continuum (D29/D30);
  - lattice: the SAME elastic RC lattice (build_lattice_rc, Elastic concrete + rebar), strut area
    calibrated so its K0 (hence T1) matches the selected reference.

Both carry the same distributed self-mass + a top lumped mass P/g, 5% modal damping at modes 1 & 2,
constant gravity P, and the SAME scaled UniformExcitation in X. The base input is selectable via
`--excitation` (default `sine`): a harmonic acceleration RESONANT with T1, or the `elcentro` record.
Because the system is linear, the response scales linearly with intensity, so the scale is set in a
single correction to hit a target peak roof drift (no iteration). A period check gates the
comparison. Output: examples/output/column/dynamic_linear/{reference}/.
Run as `python examples/column/dynamic_linear.py [--reference {beamcolumn,continuum}]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.materials import concrete_uniaxial_elastic, steel_uniaxial_elastic
from rclattice.opensees import run_beamcolumn_dynamic, run_dynamic, run_modal

from build import (
    beamcolumn_reference_linear, calibrate_area_linear, continuum_dynamic_linear,
    continuum_k0_linear, modal_calibration_figure, rc_lattice_linear,
)
from excitation import G, N_CYCLES, RHO, make_excitation, tune_intensity
from specimen import CORE, COVER_C, EPS, H, OUT, P, STEEL, THK, W

TARGET_DRIFT = 1.0     # in — keep well inside the elastic range (no yield to chase); just sets amplitude
DT = 0.01              # Newmark step (record is at 0.02 s)
REFERENCE_LABEL = {"beamcolumn": "elastic fiber beam-column", "continuum": "elastic 2D continuum"}


def main(*, reference: str = "beamcolumn", target_drift: float = TARGET_DRIFT, dt: float = DT,
         excitation: str = "sine", n_cycles: int = N_CYCLES, draw: bool = False) -> None:
    outdir = OUT / "dynamic_linear" / reference
    outdir.mkdir(parents=True, exist_ok=True)
    ref_label = REFERENCE_LABEL[reference]

    self_mass = RHO * W * H * THK            # total column self-mass (rho * volume)
    top_mass = P / G                         # axial load P treated as tributary seismic mass
    print(f"mass: self={self_mass:.3e}, top P/g={top_mass:.3e} kip-s^2/in (top/self={top_mass/self_mass:.0f}x)")

    bc_materials = (concrete_uniaxial_elastic(CORE, 1),
                    concrete_uniaxial_elastic(COVER_C, 2),
                    steel_uniaxial_elastic(STEEL, 3))

    # stiffness calibration to the SELECTED reference's K0 (mirrors the pushover)
    if reference == "continuum":
        k_ref = continuum_k0_linear()
    else:
        k_ref = (lambda r: r["shear"][1] / r["disp"][1])(beamcolumn_reference_linear())
    area, ctrl, base = calibrate_area_linear(k_ref)
    print(f"calibrated concrete strut area A = {area:.3f} in^2 (K0 match to {k_ref:.2f} kip/in, {ref_label})")

    # lattice + top seismic mass baked into the model (so modal AND dynamic both see it)
    model = rc_lattice_linear(area)

    # calibration output (D35): mode shapes (selected reference vs lattice) + periods table,
    # drawn up front on the self-mass calibration basis (BEFORE the seismic top mass is added)
    cal_caption = (f"linear scalar calibration — strut area A={area:.3f} in^2 "
                   f"(K0 match to {k_ref:.2f} kip/in, {ref_label})")
    modalpath = outdir / f"column_modal_linear_{reference}.png"
    t_ref_modal, t_lat_modal = modal_calibration_figure(
        reference=reference, lattice_model=model, label=ref_label,
        caption=cal_caption, savepath=str(modalpath), linear=True,
    )
    print(f"modal calibration: {ref_label} T={[f'{t:.4f}' for t in t_ref_modal]} s | "
          f"lattice T={[f'{t:.4f}' for t in t_lat_modal]} s -> saved {modalpath}")

    top_nodes = select_nodes(model, (-W, W, H - EPS, H + EPS))
    for nid in top_nodes:
        mx, my = model.masses[nid]
        model.masses[nid] = (mx + top_mass / len(top_nodes), my + top_mass / len(top_nodes))

    lat_T = run_modal(model, 2)["periods"]
    print(f"lattice periods: T1={lat_T[0]:.3f}s  T2={lat_T[1]:.3f}s")

    if draw:
        panels = [("RC lattice (linear)", model)]
        if reference == "continuum":
            from build import _continuum_model_linear
            panels.insert(0, ("elastic 2D continuum", _continuum_model_linear()[0]))
        drawpath = outdir / f"column_dynamic_linear_{reference}_{excitation}_model.png"
        viz.figure_model(panels, savepath=str(drawpath),
                         suptitle="LINEAR RC cantilever column — analysis model")
        print(f"saved model drawing to {drawpath}")

    # base excitation (resonant sine at T1 by default; built AFTER the modal step)
    accel, dt_rec, exc_label = make_excitation(excitation, lat_T[0], n_cycles)
    pga = max(abs(a) for a in accel)
    print(f"excitation: {exc_label} | {len(accel)} pts @ {dt_rec:.4f}s, peak={pga:.3f}g, "
          f"duration={(len(accel) - 1) * dt_rec:.1f}s")

    def run_bc(intensity: float) -> dict:
        return run_beamcolumn_dynamic(height=H, P=P, materials=bc_materials,
                                      self_mass=self_mass, top_mass=top_mass, accel=accel,
                                      dt_record=dt_rec, scale=G * intensity, dt=dt)

    # linear -> response is proportional to intensity; one trial + one correction hits the target
    used, bc = tune_intensity(run_bc, target_drift, linear=True)
    print(f"=> scale {used:.3f}x | period match: T1_lat={lat_T[0]:.3f}s vs T1_bc={bc['periods'][0]:.3f}s")

    # lattice seismic run at the SAME scaled record
    lat = run_dynamic(model, accel=accel, dt_record=dt_rec, scale=G * used,
                      control_node=ctrl, control_dof=1, base_nodes=base, dt=dt)

    # reference seismic run: reuse the fiber run for beamcolumn; run the linear continuum for continuum
    if reference == "continuum":
        print("running linear elastic continuum seismic reference...", flush=True)
        ref = continuum_dynamic_linear(accel=accel, dt_record=dt_rec, scale=G * used,
                                       top_mass=top_mass, dt=dt)
    else:
        ref = bc

    for name, r in ((ref_label, ref), ("RC lattice", lat)):
        print(f"{name:26s}: peak drift={abs(r['peak_disp']):.3f} in @ t={r['peak_time']:.2f}s | "
              f"residual={r['residual_disp']:+.4f} in | peakV={abs(r['peak_shear']):.1f} kip | conv={r['converged']}")

    series = [
        {"t": ref["t"], "disp": ref["disp"], "shear": ref["shear"],
         "label": ref_label, "style": {"color": "C3", "lw": 1.2}},
        {"t": lat["t"], "disp": lat["disp"], "shear": lat["shear"],
         "label": "elastic RC lattice", "style": {"color": "C0", "lw": 1.0, "alpha": 0.85}},
    ]
    base_title = f"LINEAR {exc_label} {used:.2f}x (T1~{lat_T[0]:.2f}s): RC lattice vs {ref_label}"
    stem = f"column_dynamic_linear_{reference}_{excitation}"
    viz.figure_timehistory(series, include_hysteresis=False,
                           savepath=str(outdir / f"{stem}.png"), title=base_title)
    # separate hysteresis loops (lattice, reference, overlaid) in addition to the overlaid histories
    viz.figure_hysteresis(list(reversed(series)),
                          savepath=str(outdir / f"{stem}_hyst.png"),
                          title=f"Hysteresis — {base_title}")
    print(f"saved to {outdir / f'{stem}.png'} (+ _hyst.png)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="LINEAR RC cantilever column seismic time-history")
    parser.add_argument("--reference", choices=("beamcolumn", "continuum"), default="beamcolumn",
                        help="verification reference: elastic fiber beam-column (fast, default) "
                             "or linear 2D plane-stress continuum (ElasticIsotropic quads)")
    parser.add_argument("--target-drift", type=float, default=TARGET_DRIFT, help="peak roof drift target (in)")
    parser.add_argument("--dt", type=float, default=DT, help="Newmark integration step (s)")
    parser.add_argument("--excitation", choices=("sine", "elcentro"), default="sine",
                        help="base excitation: resonant sine at T1 (default) or the El Centro record")
    parser.add_argument("--n-cycles", type=int, default=N_CYCLES, help="number of cycles for the sine excitation")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model (rebar highlighted)")
    args = parser.parse_args()
    main(reference=args.reference, target_drift=args.target_drift, dt=args.dt,
         excitation=args.excitation, n_cycles=args.n_cycles, draw=args.draw)
