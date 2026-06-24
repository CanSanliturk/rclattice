"""Single RC cantilever column: RC lattice vs a SINGLE-element force-based fiber beam-column —
nonlinear SEISMIC time-history OR nonlinear pushover, selectable via `--analysis`.

Companion to `dynamic.py` / `pushover.py`, with ONE difference: the fiber-beam-column reference is
idealized as a SINGLE force-based `forceBeamColumn` element with 5 Gauss-Lobatto integration points
instead of `dynamic.py`'s column subdivided into ~96 elements. Everything else mirrors the siblings
so the studies are directly comparable. Two analyses share one specimen, one fiber section
(15x24, 3+2+3 bars), the SAME Concrete02 core/cover + Steel02 grades, and one scalar K0 calibration
(lattice strut area matched to the single-element beam-column's initial lateral stiffness):

  - `--analysis dynamic` (default): nonlinear seismic time-history. Both models carry the column
    self-mass PLUS a top lumped mass P/g (axial load as tributary seismic weight), constant gravity
    P, and 5% modal damping at modes 1 & 2. The single column's self-mass is carried as DISTRIBUTED
    element mass (forceBeamColumn -mass massDens); the lattice gets the builder's lumped tributary
    self-mass. The base input is selectable via `--excitation` (default `sine`): a harmonic
    acceleration RESONANT with the lattice's T1, or the recorded El Centro 1940 NS — the SAME scaled
    UniformExcitation drives both. The intensity scale is tuned on the cheap single beam-column to a
    target peak roof drift, then drives the lattice unchanged. A period gate (T1_lat vs T1_bc)
    confirms the comparison is apples-to-apples. Output:
    single_beamcolumn_<excitation>.png (+ _hyst, model drawings).
  - `--analysis pushover`: nonlinear pushover, IDENTICAL machinery to `pushover.py --reference
    beamcolumn` (whose beam-column reference IS this single-element fiber column). Fixed base,
    constant axial P held, then a DisplacementControl lateral pushover of the top; base shear vs tip
    displacement for the single beam-column reference and the calibrated lattice. Output:
    single_beamcolumn_pushover.png (+ model drawings).

Run as `python examples/column/single_beamcolumn.py [--analysis dynamic|pushover]
[--excitation sine|elcentro]`.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openseespy.opensees as ops

from rclattice import viz
from rclattice.builders import select_nodes
from rclattice.materials import concrete_uniaxial_nonlinear, steel_uniaxial
from rclattice.opensees import run_dynamic, run_modal, run_pushover

from build import beamcolumn_reference, calibrate_area, rc_lattice
from excitation import G, N_CYCLES, RHO, make_excitation, tune_intensity
from specimen import CORE, COVER_C, DU, EPS, GF, GFC, H, HORIZON, OUT, P, STEEL, TARGET, THK, W, lateral_loads

TARGET_DRIFT = 2.0                     # in (~1.5% of H): moderate, well past yield
DT = 0.005                            # Newmark step (<= the resonant sine spacing ~T1/100)
ZETA = 0.05                            # modal damping ratio
MODES = (1, 2)                         # modes that receive modal damping
CONTROL, BASE = 2, 1                   # single-column: top (free) node / base (fixed) node

# The single beam-column's fiber-section materials, sourced from the SAME physical grades the lattice
# uses (so the two models are provably material-matched): Concrete02 core/cover + Steel02.
CORE_MAT = concrete_uniaxial_nonlinear(CORE, 1)
COVER_MAT = concrete_uniaxial_nonlinear(COVER_C, 2)
STEEL_MAT = steel_uniaxial(STEEL, 3)

SELF_MASS = RHO * W * H * THK          # total column self-mass (rho * volume)
TOP_MASS = P / G                       # axial load P treated as tributary seismic mass


# --- single-element fiber beam-column (self-contained, direct openseespy) ------------------------

def _fiber_section(tag: int) -> None:
    """The benchmark RC column fiber section: 15 wide (z) x 24 deep (y), 1.5 cover, core (mat 1) /
    cover (mat 2) / steel (mat 3), 3+2+3 bars @ 0.60 in^2 — identical geometry to the lattice column
    section. Emits the Concrete02 core/cover + Steel02 materials (the SAME grade mapping as the
    lattice) so the fiber section is material-matched to the struts."""
    ops.uniaxialMaterial(CORE_MAT.mtype, 1, *CORE_MAT.args)    # confined core
    ops.uniaxialMaterial(COVER_MAT.mtype, 2, *COVER_MAT.args)  # unconfined cover
    ops.uniaxialMaterial(STEEL_MAT.mtype, 3, *STEEL_MAT.args)
    y1, z1, cover, As = 12.0, 7.5, 1.5, 0.60
    ops.section("Fiber", tag)
    ops.patch("quad", 1, 1, 10, -y1 + cover, z1 - cover, -y1 + cover, -z1 + cover,
              y1 - cover, -z1 + cover, y1 - cover, z1 - cover)               # core
    ops.patch("quad", 2, 1, 1, -y1, z1, -y1, -z1, -y1 + cover, -z1 + cover, -y1 + cover, z1 - cover)
    ops.patch("quad", 2, 1, 1, y1 - cover, z1 - cover, y1 - cover, -z1 + cover, y1, -z1, y1, z1)
    ops.patch("quad", 2, 1, 1, -y1 + cover, z1 - cover, -y1 + cover, z1, y1 - cover, z1, y1 - cover, z1 - cover)
    ops.patch("quad", 2, 1, 1, -y1 + cover, -z1, -y1 + cover, -z1 + cover, y1 - cover, -z1 + cover, y1 - cover, -z1)
    ops.layer("straight", 3, 3, As, -y1 + cover, z1 - cover, -y1 + cover, -z1 + cover)  # 3 bars y=-10.5
    ops.layer("straight", 3, 2, As, 0.0, z1 - cover, 0.0, -z1 + cover)                  # 2 bars y=0
    ops.layer("straight", 3, 3, As, y1 - cover, z1 - cover, y1 - cover, -z1 + cover)    # 3 bars y=+10.5


def _build_single_column() -> None:
    """Wipe + build the single-element column domain: fixed base (node 1), free top (node 2), one
    force-based fiber `forceBeamColumn` with 5 Lobatto points and a P-Delta transform. The column
    self-mass is carried as DISTRIBUTED element mass (massDens = self_mass/H) and the axial-load
    tributary seismic mass P/g is lumped at the free top node (translational only)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, H)
    ops.fix(1, 1, 1, 1)
    ops.mass(2, TOP_MASS, TOP_MASS, 0.0)                 # axial-load tributary seismic mass (top)

    _fiber_section(1)
    ops.geomTransf("PDelta", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)              # 5 integration points, Gauss-Lobatto
    ops.element("forceBeamColumn", 1, 1, 2, 1, 1, "-mass", SELF_MASS / H)


def _hold_gravity(tol: float = 1e-8) -> bool:
    """Apply the constant axial load P at the top in 10 LoadControl steps, then hold it (loadConst).
    Returns True on convergence."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -P, 0.0)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", tol, 100)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if ops.analyze(10) != 0:
        return False
    ops.loadConst("-time", 0.0)
    return True


def run_single(intensity: float, accel: list[float], dt_record: float, dt: float = DT,
               tol: float = 1e-6, max_iter: int = 50) -> dict:
    """Full seismic run of the single beam-column at excitation `scale = G * intensity`: build ->
    gravity (held) -> 5% modal damping at modes 1 & 2 -> UniformExcitation HHT-alpha transient in X.
    Records relative roof drift (control node X) and base shear (-base reaction X) each step. Returns
    the same dict shape as run_dynamic (t/disp/shear histories + peaks + modal "periods" + converged).
    This is the cheap tuning proxy AND the comparison reference."""
    _build_single_column()
    if not _hold_gravity():
        return {"converged": False, "stage": "gravity", "t": [], "disp": [], "shear": [],
                "peak_disp": 0.0, "peak_time": 0.0, "residual_disp": 0.0, "peak_shear": 0.0,
                "periods": [math.inf, math.inf]}

    scale = G * intensity
    ops.timeSeries("Path", 90, "-dt", dt_record, "-values", *accel, "-factor", scale)
    ops.pattern("UniformExcitation", 90, 1, "-accel", 90)

    # Drop the gravity stage's Static analysis so the transient integrator (HHT) is honored, then
    # re-declare the solver. Modal damping is formed on the gravity-held tangent and stored on the
    # domain so the transient integrator created next uses it (same recipe as opensees.run_dynamic).
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    nmode = max(MODES)
    try:
        lams = ops.eigen(nmode)
    except Exception:
        lams = ops.eigen("-fullGenLapack", nmode)
    ops.modalDamping(ZETA)
    periods = [2.0 * math.pi / math.sqrt(l) if l > 0 else math.inf for l in lams]

    ops.test("NormDispIncr", tol, max_iter)
    ops.algorithm("Newton")
    ops.integrator("HHT", 0.7)
    ops.analysis("Transient")

    def cdisp() -> float:
        return ops.nodeDisp(CONTROL, 1)

    def base_shear() -> float:
        ops.reactions()
        return -ops.nodeReaction(BASE)[0]

    total_time = (len(accel) - 1) * dt_record
    nsteps = int(round(total_time / dt))
    t = [0.0]; disp = [cdisp()]; shear = [base_shear()]
    converged = True
    algos = (("Newton", ()), ("KrylovNewton", ()),
             ("NewtonLineSearch", ("-type", "Bisection")), ("ModifiedNewton", ()))
    for k in range(nsteps):
        if ops.analyze(1, dt) != 0:
            ok = False
            ops.test("NormDispIncr", tol * 10, max_iter * 4)
            for algo, a in algos:
                ops.algorithm(algo, *a)
                for nsub in (10, 50, 200):
                    if ops.analyze(nsub, dt / nsub) == 0:
                        ok = True
                        break
                if ok:
                    break
            ops.test("NormDispIncr", tol, max_iter)
            ops.algorithm("Newton")
            if not ok:
                converged = False
                break
        t.append((k + 1) * dt); disp.append(cdisp()); shear.append(base_shear())
    ipk = max(range(len(disp)), key=lambda m: abs(disp[m]))
    return {"t": t, "disp": disp, "shear": shear, "converged": converged, "periods": periods,
            "peak_disp": disp[ipk], "peak_time": t[ipk], "residual_disp": disp[-1],
            "peak_shear": max(shear, key=abs), "control_node": CONTROL, "base_nodes": [BASE]}


def draw_single_column(savepath: Path, *, top_label: str = "top mass P/g") -> None:
    """Schematic of the single beam-column model: column elevation (5 Lobatto points, fixed base,
    top node) + the fiber section (core/cover patches and longitudinal rebar). `top_label` labels
    the top-node marker ("top mass P/g" for the dynamic study, "axial load P" for the pushover)."""
    fig, (axE, axS) = plt.subplots(1, 2, figsize=(8, 6),
                                   gridspec_kw={"width_ratios": [1, 1.2]})
    axE.plot([0, 0], [0, H], color="0.4", lw=3, zorder=1)
    lob = [0.5 - 0.5 * math.cos(math.pi * i / 4) for i in range(5)]  # Lobatto points (endpoints incl.)
    axE.scatter([0] * 5, [p * H for p in lob], s=40, color="C1", zorder=3, label="5 Lobatto pts")
    axE.scatter([0], [H], s=160, marker="s", color="C3", zorder=4, label=top_label)
    axE.scatter([0], [0], s=160, marker="^", color="k", zorder=4, label="fixed base")
    axE.set_xlim(-W, W); axE.set_ylim(-8, H + 8); axE.set_aspect("equal")
    axE.set_title("forceBeamColumn (single element)"); axE.legend(loc="center right", fontsize=8)
    axE.set_xlabel("x (in)"); axE.set_ylabel("y (in)")

    y1, z1, c = 12.0, 7.5, 1.5
    axS.add_patch(plt.Rectangle((-z1, -y1), 2 * z1, 2 * y1, fc="0.85", ec="0.5", label="cover"))
    axS.add_patch(plt.Rectangle((-z1 + c, -y1 + c), 2 * (z1 - c), 2 * (y1 - c),
                                fc="0.65", ec="0.4", label="core"))
    bars = [(zz, yy) for yy in (-y1 + c, 0.0, y1 - c) for zz in (-z1 + c, z1 - c)]
    axS.scatter([b[0] for b in bars], [b[1] for b in bars], s=80, color="C0", zorder=3,
                label="longitudinal rebar")
    axS.set_xlim(-z1 - 2, z1 + 2); axS.set_ylim(-y1 - 2, y1 + 2); axS.set_aspect("equal")
    axS.set_title("fiber section (15 x 24, 3+2+3 @0.6 in$^2$)")
    axS.legend(loc="upper right", fontsize=8); axS.set_xlabel("z (in)"); axS.set_ylabel("y (in)")

    fig.suptitle("Single-element RC cantilever column — analysis model")
    fig.tight_layout()
    fig.savefig(savepath, dpi=130)
    plt.close(fig)


# --- comparison drivers -------------------------------------------------------------------------

def main_dynamic(*, regularize: bool = True, Gf: float = GF, Gfc: float = GFC, horizon: float = HORIZON,
                 target_drift: float = TARGET_DRIFT, dt: float = DT, excitation: str = "sine",
                 n_cycles: int = N_CYCLES, draw: bool = True) -> None:
    """Nonlinear seismic time-history: RC lattice vs the single-element fiber beam-column."""
    outdir = OUT / "single_beamcolumn"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"single_beamcolumn_{excitation}"

    print(f"mass: self={SELF_MASS:.3e}, top P/g={TOP_MASS:.3e} kip-s^2/in  "
          f"(top/self = {TOP_MASS / SELF_MASS:.0f}x)")

    # stiffness calibration: lattice strut area -> the single beam-column's K0 (mirrors dynamic.py's
    # beamcolumn path; beamcolumn_reference IS the single-element fiber column, same model as run_single)
    k_ref = (lambda r: r["shear"][1] / r["disp"][1])(beamcolumn_reference())
    area, ctrl, base = calibrate_area(k_ref, horizon=horizon)
    print(f"calibrated strut area A = {area:.3f} in^2 (K0 match to single beam-column) | horizon={horizon:g}")

    # calibrated lattice + top seismic mass baked in (so modal AND dynamic both see it)
    model = rc_lattice(regularize, Gf, Gfc, area, horizon=horizon)
    top_nodes = select_nodes(model, (-W, W, H - EPS, H + EPS))
    for nid in top_nodes:
        mx, my = model.masses[nid]
        model.masses[nid] = (mx + TOP_MASS / len(top_nodes), my + TOP_MASS / len(top_nodes))

    # period gate (cheap eigen on the lattice) — also the resonant period for the sine input
    lat_T = run_modal(model, 2)["periods"]
    print(f"lattice periods: T1={lat_T[0]:.3f}s  T2={lat_T[1]:.3f}s")

    # draw both analysis models up front (before the expensive lattice time-history)
    if draw:
        viz.figure_model([("RC lattice", model)], savepath=str(outdir / f"{stem}_lattice_model.png"),
                         suptitle="RC cantilever column — lattice (vs single beam-column)")
        draw_single_column(outdir / f"{stem}_bc_model.png")
        print(f"saved model drawings to {outdir / f'{stem}_lattice_model.png'} and _bc_model.png")

    # base excitation (resonant sine at the lattice T1 by default; built AFTER the modal step)
    accel, dt_rec, exc_label = make_excitation(excitation, lat_T[0], n_cycles)
    pga = max(abs(a) for a in accel)
    print(f"excitation: {exc_label} | {len(accel)} pts @ {dt_rec:.4f}s, peak={pga:.3f}g, "
          f"duration={(len(accel) - 1) * dt_rec:.1f}s")

    # tune the excitation scale on the cheap single beam-column to ~target_drift peak roof drift,
    # then drive the lattice with the SAME scaled record (exactly as dynamic.py tunes on its fiber column)
    used, bc = tune_intensity(lambda i: run_single(i, accel, dt_rec, dt), target_drift, linear=False)
    print(f"=> scale {used:.3f}x (single-BC peak drift {abs(bc['peak_disp']):.3f} in) | "
          f"period match: T1_lat={lat_T[0]:.3f}s vs T1_bc={bc['periods'][0]:.3f}s")

    # lattice seismic run at the SAME scaled record (expensive)
    lat = run_dynamic(model, accel=accel, dt_record=dt_rec, scale=G * used,
                      control_node=ctrl, control_dof=1, base_nodes=base, dt=dt)

    ref = bc
    for name, r in (("single beam-column", ref), ("RC lattice", lat)):
        print(f"{name:18s}: peak drift={abs(r['peak_disp']):.3f} in @ t={r['peak_time']:.2f}s | "
              f"residual={r['residual_disp']:+.3f} in | peakV={abs(r['peak_shear']):.1f} kip | conv={r['converged']}")

    series = [
        {"t": ref["t"], "disp": ref["disp"], "shear": ref["shear"],
         "label": "single beam-column", "style": {"color": "C3", "lw": 1.2}},
        {"t": lat["t"], "disp": lat["disp"], "shear": lat["shear"],
         "label": "RC lattice", "style": {"color": "C0", "lw": 1.0, "alpha": 0.85}},
    ]
    base_title = f"{exc_label} {used:.2f}x (T1~{lat_T[0]:.2f}s): RC lattice vs single beam-column"
    viz.figure_timehistory(series, include_hysteresis=False,
                           savepath=str(outdir / f"{stem}.png"), title=base_title)
    viz.figure_hysteresis(list(reversed(series)),
                          savepath=str(outdir / f"{stem}_hyst.png"),
                          title=f"Hysteresis — {base_title}")
    print(f"saved to {outdir / f'{stem}.png'} (+ _hyst.png)")


def main_pushover(*, regularize: bool = True, Gf: float = GF, Gfc: float = GFC,
                  horizon: float = HORIZON, draw: bool = True) -> None:
    """Nonlinear pushover: RC lattice vs the single-element fiber beam-column.

    IDENTICAL machinery to `pushover.py --reference beamcolumn` — whose beam-column reference IS this
    single-element force-based fiber column (`beamcolumn_reference` -> `run_beamcolumn_cantilever`,
    one `forceBeamColumn`, 5 Lobatto points). Both: fixed base, constant axial P held, then a
    DisplacementControl lateral pushover of the top. The lattice strut area is scalar-calibrated so
    its initial lateral stiffness K0 equals the single beam-column's (the SAME calibration the
    dynamic path uses), then both are pushed to TARGET and their base-shear-vs-tip-drift curves
    compared."""
    outdir = OUT / "single_beamcolumn"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = "single_beamcolumn_pushover"

    # reference: the single-element fiber beam-column pushover curve (== pushover.py's beamcolumn ref)
    ref = beamcolumn_reference()
    k_ref = ref["shear"][1] / ref["disp"][1]   # initial (small-drift) lateral stiffness
    print(f"single beam-column: K0={k_ref:.2f} kip/in | peakV={max(ref['shear']):.2f} kip | "
          f"drift->{ref['disp'][-1]:.2f} in (conv={ref['converged']})")

    # scalar calibration: lattice strut area -> single beam-column K0 (same as the dynamic sibling)
    area, ctrl, base = calibrate_area(k_ref, horizon=horizon)
    mode = f"length-regularized (Gf={Gf:g}, Gfc={Gfc:g})" if regularize else "plain (no regularization)"
    print(f"calibrated strut area A = {area:.3f} in^2 (K0 match to single beam-column) | "
          f"concrete law: Concrete02 {mode} | horizon={horizon:g}")

    model = rc_lattice(regularize, Gf, Gfc, area, horizon=horizon)
    print(f"RC lattice column: {len(model.nodes)} nodes, {len(model.elements)} struts")

    # draw both analysis models up front (before the expensive lattice pushover)
    if draw:
        viz.figure_model([("RC lattice", model)], savepath=str(outdir / f"{stem}_lattice_model.png"),
                         suptitle="RC cantilever column — lattice (vs single beam-column)")
        draw_single_column(outdir / f"{stem}_bc_model.png", top_label="axial load P")
        print(f"saved model drawings to {outdir / f'{stem}_lattice_model.png'} and _bc_model.png")

    lat = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)
    print(f"lattice: peakV={max(lat['shear']):.2f} kip | drift->{lat['disp'][-1]:.2f} in "
          f"(conv={lat['converged']})")

    viz.figure_pushover(
        [
            {"disp": ref["disp"], "shear": ref["shear"], "label": "single beam-column",
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat["disp"], "shear": lat["shear"], "label": "RC lattice (scalar calib.)",
             "style": {"color": "C0", "ls": "--", "lw": 2, "marker": "."}},
        ],
        savepath=str(outdir / f"{stem}.png"),
        xlabel="tip displacement (in)", ylabel="base shear (kip)",
        title="RC cantilever column pushover: lattice vs single beam-column",
    )
    print(f"saved pushover curve to {outdir / f'{stem}.png'}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="RC lattice vs single-element fiber beam-column — nonlinear seismic time-history "
                    "or nonlinear pushover")
    parser.add_argument("--analysis", choices=("dynamic", "pushover"), default="dynamic",
                        help="analysis type: seismic time-history (default) or nonlinear pushover "
                             "(same machinery as pushover.py --reference beamcolumn)")
    parser.add_argument("--excitation", choices=("sine", "elcentro"), default="sine",
                        help="[dynamic] base excitation: resonant sine at the lattice T1 (default) "
                             "or El Centro")
    parser.add_argument("--no-regularize", dest="regularize", action="store_false",
                        help="use plain Concrete02 instead of the crack-band length-regularized law")
    parser.add_argument("--gf", type=float, default=GF, help="tensile fracture energy (kip, in)")
    parser.add_argument("--gfc", type=float, default=GFC, help="compressive fracture energy (kip, in)")
    parser.add_argument("--horizon", type=float, default=HORIZON,
                        help="strut connectivity horizon * mesh_size (default 1.5)")
    parser.add_argument("--target-drift", type=float, default=TARGET_DRIFT,
                        help="[dynamic] peak roof drift target (in)")
    parser.add_argument("--dt", type=float, default=DT, help="[dynamic] Newmark integration step (s)")
    parser.add_argument("--n-cycles", type=int, default=N_CYCLES,
                        help="[dynamic] number of cycles for the sine excitation")
    parser.add_argument("--no-draw", dest="draw", action="store_false",
                        help="skip the analysis-model drawings")
    args = parser.parse_args()
    if args.analysis == "pushover":
        main_pushover(regularize=args.regularize, Gf=args.gf, Gfc=args.gfc, horizon=args.horizon,
                      draw=args.draw)
    else:
        main_dynamic(regularize=args.regularize, Gf=args.gf, Gfc=args.gfc, horizon=args.horizon,
                     target_drift=args.target_drift, dt=args.dt, excitation=args.excitation,
                     n_cycles=args.n_cycles, draw=args.draw)
