"""Aldemir's squat wall as an ELASTIC lattice — the first real check on this specimen.

Two things make this run worth doing before anything expensive.

1. SHEAR IS THE ELASTIC RESPONSE HERE. At an aspect ratio of 0.75, shear carries ~58% of the
   uncracked flexibility, against ~14% for VK3 (2.20) and ~5% for a slender wall. The Aydin energy
   balance matches the lattice's NORMAL stiffness (C11) exactly and never matches its shear stiffness
   independently — `nu_consistent` is not a Poisson ratio and the lattice is cubic-symmetric at any
   nu (D53). Every previous wall study could treat that as a footnote. This one cannot: the
   calibration's weakest property is now the majority of the answer.

2. THE GEOMETRY IS UNRESOLVED and this run can price it. The 2019 paper gives two incompatible
   panels and its own Table 2 counts imply a third (see `summary.py`). `--geometry-sweep` rebuilds
   and re-solves across all of them and prints K against the measured 1,038 kN/mm, which turns a
   hand calculation into a measurement.

WHAT THE TARGET IS. The measured initial stiffness of 1,038.44 kN/mm (Table 4) is the only elastic
number the paper reports, and it is NOT a safe target for the adopted Fig. 10(a) geometry: an
uncracked transformed-section cantilever of 3000 x 2250 x 120 gives only ~804 kN/mm, i.e. the
measurement exceeds what that geometry can produce. Expect K_lattice/K_measured ~ 0.7 and read it as
the open geometry question, NOT as a calibration failure. `K_lattice/K_transformed` is the
like-for-like ratio to judge the calibration by; the other three walls land at 0.93.

Output: examples/output/aydin_aldemir_wall/runs/<stamp>_elastic/. Units: N, mm.
Run from src/:
    uv run python examples/aydin_aldemir_wall/elastic.py [--mesh 50] [--horizon 1.5]
        [--calibration uniaxial|equibiaxial] [--bond] [--geometry-sweep]
    (the analysis-model drawing is its own script: draw.py)
"""

from __future__ import annotations

import json

from rclattice import viz
from rclattice.opensees import run_pushover

from build import (
    calibrate, cantilever_stiffness, describe, gross_inertia, report_calibration,
    transformed_inertia, wall_continuum, wall_lattice,
)
from specimen import (
    A_SHEAR, EC, ES, HORIZON, HW, LW, MESH, PAPER_GRID, PAPER_GRID_T, SPECIMEN, TW, base_nodes,
    control_node, lateral_loads, mesh_fits, run_dir,
)

TARGET_DRIFT = 0.0002        # 0.02% — an order below the ~0.07% drift at which this wall cracks
K_MEASURED = 1038.44e3       # N/mm, Table 4
F_MEASURED = 963.592e3       # N,    Table 4


def solve(area, *, mesh_size, horizon, bond, length, height):
    """Build and push the elastic lattice to `TARGET_DRIFT`; return (model, result, K)."""
    model, _edges = wall_lattice(area, mesh_size=mesh_size, horizon=horizon, nonlinear=False,
                                 bond=bond, length=length, height=height)
    target = TARGET_DRIFT * height
    res = run_pushover(model, lateral_loads=lateral_loads(model, height=height),
                       control_node=control_node(model, height=height), control_dof=1,
                       dU=target / 40.0, target=target, base_nodes=base_nodes(model))
    if not res["converged"]:
        raise RuntimeError(f"elastic pushover on {length:.0f}x{height:.0f} did not reach target")
    return model, res, res["shear"][-1] / res["disp"][-1]


def push(model, *, height: float = HW):
    """Push a built model to `TARGET_DRIFT` and return (result, secant stiffness)."""
    target = TARGET_DRIFT * height
    res = run_pushover(model, lateral_loads=lateral_loads(model, height=height),
                       control_node=control_node(model, height=height), control_dof=1,
                       dU=target / 40.0, target=target, base_nodes=base_nodes(model))
    if not res["converged"]:
        raise RuntimeError("elastic pushover did not reach the target drift")
    return res, res["shear"][-1] / res["disp"][-1]


def mesh_convergence(field: str, horizon: float, meshes) -> list[dict]:
    """K for lattice and continuum across meshes — the test of the balance's mesh objectivity.

    The energy balance claims the calibrated AREA is a material property of (mesh, horizon, nu,
    thickness), which implies the assembled stiffness should be mesh-INDEPENDENT. That is a claim
    about the lattice, so it is checked against a continuum refined on the same grid: if both drift
    together the mesh is under-resolved, and if only the lattice drifts the calibration is not
    objective.

    Legal meshes divide gcd(3000, 2250, 100, 50) = 50, so the usable set is small: 50, 25, 10, ...
    """
    from rclattice.calibration import energy_balance_rectangle
    rows = []
    for m in meshes:
        if not mesh_fits(m):
            continue
        cal = energy_balance_rectangle(LW, HW, m, E=EC, nu=0.20, thickness=TW, horizon=horizon,
                                       field=field)
        lat, _e = wall_lattice(cal.area, mesh_size=m, horizon=horizon, nonlinear=False)
        _r, k_lat = push(lat)
        con, _e2 = wall_continuum(mesh_size=m)
        _r2, k_con = push(con)
        rows.append({"mesh": m, "area": cal.area, "K_lattice": k_lat, "K_continuum": k_con,
                     "nodes": len(lat.nodes)})
    return rows


def geometry_sweep(area_of, *, mesh_size, horizon, bond) -> list[dict]:
    """Re-solve on every candidate panel, so the open geometry question gets a MEASURED answer.

    `area_of(length, height)` re-runs the energy balance per panel — the calibrated area is a
    property of the grid and thickness, not of the rectangle, so it should barely move; printing it
    per row is the check that it does not.

    Returns `(rows, skipped)`; a panel whose dimensions the mesh cannot carry is skipped rather than
    silently rounded, because rounding the height is exactly what the drift comparison turns on.
    """
    # Fig. 4(f) is NOT a candidate and never was: it is a DETAIL VIEW of the bottom 1500 mm of his
    # model, so that number is a crop height (D82). What replaces it is the real open question —
    # inverting his Table 2 counts is unique only UP TO TRANSPOSITION, and the two resolutions are a
    # squat panel and a slender one. `replica/` assumes the first; the sweep prices both.
    cands = [("Fig. 10(a)  3000x2250", LW, HW),
             ("Table 2     3000x2680", *PAPER_GRID),
             ("Table 2 T   2680x3000", *PAPER_GRID_T)]
    rows, skipped = [], []
    for label, L, H in list(cands):
        if not mesh_fits(mesh_size, length=L, height=H):
            skipped.append((label, L, H))
            cands.remove((label, L, H))
    for label, L, H in cands:
        cal = area_of(L, H)
        _m, _r, k = solve(cal.area, mesh_size=mesh_size, horizon=horizon, bond=bond,
                          length=L, height=H)
        k_tr, share = cantilever_stiffness(
            shear_span=H, inertia=transformed_inertia(mesh_size, length=L, height=H), length=L)
        rows.append({"label": label, "length": L, "height": H, "area": cal.area,
                     "K": k, "K_transformed": k_tr, "shear_share": share})
    return rows, skipped


def main(*, mesh_size: float = MESH, horizon: float = HORIZON, field: str = "uniaxial",
         bond: bool = False, sweep: bool = False, convergence: bool = False) -> None:
    out = run_dir("elastic", tag=f"{field}{'_bond' if bond else ''}")
    print(f"output directory: {out}")

    def area_of(length, height):
        from rclattice.calibration import energy_balance_rectangle
        return energy_balance_rectangle(length, height, mesh_size, E=EC, nu=0.20, thickness=TW,
                                        horizon=horizon, field=field)

    cal = calibrate(mesh_size=mesh_size, horizon=horizon, field=field)
    report_calibration(cal, mesh_size=mesh_size, horizon=horizon)

    model, res, k_lat = solve(cal.area, mesh_size=mesh_size, horizon=horizon, bond=bond,
                              length=LW, height=HW)
    print(f"\n{SPECIMEN}  {LW:.0f} x {HW:.0f} x {TW:.0f}, mesh {mesh_size:g}, horizon {horizon:g}, "
          f"{'BOND elements' if bond else 'perfect bond'}")
    print(f"  {describe(model)}")

    i_g = gross_inertia()
    i_tr = transformed_inertia(mesh_size)
    k_gross, _ = cantilever_stiffness(shear_span=A_SHEAR, inertia=i_g)
    k_tr, share = cantilever_stiffness(shear_span=A_SHEAR, inertia=i_tr)
    secants = [s / u for u, s in zip(res["disp"][1:], res["shear"][1:])]
    spread = (max(secants) - min(secants)) / k_lat

    print(f"\nelastic response to {TARGET_DRIFT * 100:.2f}% drift "
          f"({TARGET_DRIFT * A_SHEAR:.3f} mm at y = {A_SHEAR:.0f} mm)")
    print(f"  lattice secant stiffness      {k_lat / 1e3:8.2f} kN/mm")
    print(f"  cantilever, plain concrete    {k_gross / 1e3:8.2f} kN/mm  (gross section)")
    print(f"  cantilever, transformed       {k_tr / 1e3:8.2f} kN/mm  (n = Es/Ec = {ES / EC:.1f}; "
          f"I_tr/I_g = {i_tr / i_g:.3f})")
    print(f"  ratio  K_lattice / K_transf.  {k_lat / k_tr:8.4f}   (the other three walls land at 0.93 "
          f"— but see the continuum below)")
    print(f"  secant spread over the pull   {spread * 100:8.4f}%  "
          f"({'LINEAR' if spread < 1e-3 else 'NONLINEAR — investigate'})")

    con_model, _ce = wall_continuum(mesh_size=mesh_size)
    _cres, k_con = push(con_model)
    print(f"\n  CONTINUUM reference, same grid / rebar / BCs ({describe(con_model)})")
    print(f"    plane-stress continuum      {k_con / 1e3:8.2f} kN/mm")
    print(f"    ratio  K_lattice / K_cont.  {k_lat / k_con:8.4f}   <- JUDGE THE CALIBRATION ON THIS")
    print(f"    ratio  K_cont. / K_transf.  {k_con / k_tr:8.4f}   <- how wrong BEAM THEORY is here; "
          f"at aspect {HW / LW:.2f} a")
    print("      cantilever formula assumes plane sections and uniform shear, and this wall has")
    print("      neither. Any gap here is the reference's error, not the lattice's.")

    print(f"\n  SHEAR SHARE of the elastic flexibility  {share * 100:5.1f}%  at Lv/lw = "
          f"{A_SHEAR / LW:.2f}")
    print("  For scale: VK3 (2.20) is ~14% and SW-NC-FF (3.00) is a few percent. The energy balance")
    print(f"  never matches shear stiffness independently — its error at nu = {cal.nu:.2f} is "
          f"{cal.isotropy_error * 100:.2f}% —")
    print("  and on this specimen that error is applied to the majority of the response.")

    print(f"\n  MEASURED initial stiffness    {K_MEASURED / 1e3:8.2f} kN/mm  (Table 4)")
    print(f"  ratio  K_lattice / K_measured {k_lat / K_MEASURED:8.4f}")
    print(f"  ratio  K_transf. / K_measured {k_tr / K_MEASURED:8.4f}   <- the geometry's OWN ceiling")
    if k_tr < K_MEASURED:
        print("  The measurement EXCEEDS what an uncracked section of the adopted geometry can give,")
        print("  so the lattice cannot reach it either. That is the open geometry question (see")
        print("  summary.py), not a calibration failure. Read the K_transformed ratio instead.")

    rows = []
    if sweep:
        print("\n" + "-" * 96)
        print("  GEOMETRY SWEEP — the same calibration and solver on every candidate panel")
        print("-" * 96)
        print(f"  {'panel':<24s}{'A_t':>10s}{'K_lat':>10s}{'K_transf':>10s}{'K/K_tr':>9s}"
              f"{'K/K_meas':>10s}{'shear %':>9s}")
        rows, skipped = geometry_sweep(area_of, mesh_size=mesh_size, horizon=horizon, bond=bond)
        for r in rows:
            print(f"  {r['label']:<24s}{r['area']:>10,.0f}{r['K'] / 1e3:>10,.1f}"
                  f"{r['K_transformed'] / 1e3:>10,.1f}{r['K'] / r['K_transformed']:>9.3f}"
                  f"{r['K'] / K_MEASURED:>10.3f}{r['shear_share'] * 100:>9.1f}")
        for label, L, H in skipped:
            print(f"  {label:<24s}  SKIPPED — mesh {mesh_size:g} does not divide {L:g} x {H:g}; "
                  f"use --mesh 10 to include it")
        print("\n  A panel whose K/K_meas lands near 1.0 with K/K_tr near 0.93 is the one the")
        print("  measurement supports. K/K_tr should be ~constant across panels — it is the")
        print("  calibration; K/K_meas is the geometry.")

    conv = []
    if convergence:
        print("\n" + "-" * 96)
        print("  MESH CONVERGENCE — is the calibrated area really a material property?")
        print("-" * 96)
        print(f"  {'mesh':>6s}{'nodes':>9s}{'A_t':>10s}{'K_lattice':>12s}{'K_continuum':>13s}"
              f"{'lat/cont':>10s}")
        conv = mesh_convergence(field, horizon, (50.0, 25.0))
        for r in conv:
            print(f"  {r['mesh']:>6.0f}{r['nodes']:>9,d}{r['area']:>10,.0f}"
                  f"{r['K_lattice'] / 1e3:>12,.1f}{r['K_continuum'] / 1e3:>13,.1f}"
                  f"{r['K_lattice'] / r['K_continuum']:>10.4f}")
        if len(conv) > 1:
            dl = conv[-1]["K_lattice"] / conv[0]["K_lattice"] - 1.0
            dc = conv[-1]["K_continuum"] / conv[0]["K_continuum"] - 1.0
            print(f"\n  refining 50 -> 25 moves the lattice {dl * 100:+.2f}% and the continuum "
                  f"{dc * 100:+.2f}%.")
            print("  Both drifting together is mesh resolution; only the lattice drifting would mean")
            print("  the calibrated area is not mesh-objective.")

    drift = [u / A_SHEAR * 100.0 for u in res["disp"]]
    curves = [
        {"disp": drift, "shear": [s / 1e3 for s in res["shear"]],
         "label": f"lattice [{field}], A$_t$={cal.area:,.0f} mm² ({k_lat / 1e3:.0f} kN/mm)",
         "style": {"color": "C0", "lw": 2, "marker": ".", "markevery": 5}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, k_tr * res["disp"][-1] / 1e3],
         "label": f"transformed-section cantilever ({k_tr / 1e3:.0f} kN/mm)",
         "style": {"color": "C3", "ls": "--", "lw": 2}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, k_gross * res["disp"][-1] / 1e3],
         "label": f"plain-concrete cantilever ({k_gross / 1e3:.0f} kN/mm)",
         "style": {"color": "0.6", "ls": ":", "lw": 1.5}},
        {"disp": [0.0, drift[-1]], "shear": [0.0, K_MEASURED * res["disp"][-1] / 1e3],
         "label": f"MEASURED initial stiffness ({K_MEASURED / 1e3:.0f} kN/mm, Table 4)",
         "style": {"color": "C2", "ls": "-.", "lw": 2}},
    ]
    viz.figure_pushover(curves, savepath=str(out / "elastic.png"), xlabel="drift ratio (%)",
                        ylabel="lateral load (kN)",
                        title=f"{SPECIMEN} — elastic lattice, {field} energy balance, "
                              f"{'bond elements' if bond else 'perfect bond'}")

    (out / "data.json").write_text(json.dumps({
        "specimen": SPECIMEN, "panel": [LW, HW, TW], "mesh": mesh_size, "horizon": horizon,
        "calibration_field": field, "bond": bond, "area": cal.area, "EA": cal.EA,
        "K_lattice": k_lat, "K_gross": k_gross, "K_transformed": k_tr, "K_continuum": k_con,
        "K_measured": K_MEASURED, "shear_share": share,
        "isotropy_error": cal.isotropy_error, "nu_effective": cal.nu_effective,
        "disp": res["disp"], "shear": res["shear"],
        "geometry_sweep": rows, "mesh_convergence": conv,
    }, indent=2))
    print(f"\nsaved to {out}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Aldemir squat wall — elastic lattice")
    p.add_argument("--mesh", type=float, default=MESH)
    p.add_argument("--horizon", type=float, default=HORIZON)
    p.add_argument("--calibration", choices=("uniaxial", "equibiaxial"), default="uniaxial",
                   help="which energy-balance affine field sets the strut area (D72)")
    p.add_argument("--geometry-sweep", dest="sweep", action="store_true",
                   help="re-solve on every candidate panel and score them against the measured K")
    p.add_argument("--mesh-convergence", dest="convergence", action="store_true",
                   help="re-solve lattice AND continuum at mesh 50 and 25 to test mesh objectivity")
    a = p.parse_args()
    main(mesh_size=a.mesh, horizon=a.horizon, field=a.calibration, bond=False, sweep=a.sweep,
         convergence=a.convergence)
