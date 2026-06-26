"""Single RC cantilever column: lattice vs force-based fiber beam-column (D12 single-member check).

The cleanest verification of the RC lattice: model ONE column two ways and compare their plastic
pushover curves —
  - reference: a force-based fiber `forceBeamColumn` cantilever (run_beamcolumn_cantilever),
  - lattice: the same 24-in-deep x 144-tall column as a regularized RC lattice (build_lattice_rc).

Both: fixed base, constant axial P at the free top, then a DisplacementControl lateral pushover
of the top. The lattice concrete strut area is calibrated so its initial lateral stiffness equals
the beam-column's (so the two start from the same elastic slope). The lattice concrete softening
law is CONFIGURABLE: the crack-band fracture-energy length-regularized Concrete02 (D20, default)
or plain Concrete02 — toggle via `main(regularize=...)` or the `--no-regularize` CLI flag; `Gf`/`Gfc`
are likewise parametrized (`main(Gf=, Gfc=)` / `--gf` / `--gfc`). Units: kip, in. Output:
examples/output/column/column_pushover.png.

The specimen + builders are imported from this package (specimen.py / build.py), shared with the
linear and dynamic siblings. Run as `python examples/column/pushover.py`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import (
    _continuum_model, beamcolumn_reference, calibrate_area, calibrate_groups, make_reference,
    modal_calibration_figure, rc_lattice,
)
from specimen import DU, GF, GFC, H, HORIZON, MESH, OUT, TARGET, lateral_loads

REFERENCE_LABEL = {"beamcolumn": "fiber beam-column", "continuum": "2D continuum"}


def main(*, reference: str = "beamcolumn", calibration: str = "scalar", horizon: float = HORIZON,
         regularize: bool = True, Gf: float = GF, Gfc: float = GFC, draw: bool = False) -> None:
    outdir = OUT / "pushover" / reference
    outdir.mkdir(parents=True, exist_ok=True)

    ref = make_reference(reference)
    label = REFERENCE_LABEL[reference]
    k_ref = ref["shear"][1] / ref["disp"][1]   # initial (small-drift) lateral stiffness
    print(f"{label}: K0={k_ref:.2f} kip/in | peakV={max(ref['shear']):.2f} kip | "
          f"drift->{ref['disp'][-1]:.2f} in (conv={ref['converged']})")

    if calibration == "groups":
        area, ctrl, base, cal = calibrate_groups(horizon=horizon)   # strong 2-group fit (D16)
        o, d = cal.areas["orthogonal"], cal.areas["diagonal"]
        caption = (f"groups calibration — A_ortho={o:.3f}, A_diag={d:.3f} in^2 (d/o {d/o:.2f}), "
                   f"RMS {cal.rms:.3f}, success={cal.success}")
        print(f"strong 2-group calibration (continuum static+modal): orthogonal A={o:.3f}, "
              f"diagonal A={d:.3f} in^2 (d/o ratio {d/o:.2f}, RMS {cal.rms:.3f}, success={cal.success})")
    else:
        area, ctrl, base = calibrate_area(k_ref, horizon=horizon)
        caption = f"scalar calibration — strut area A={area:.3f} in^2 (K0 match to {label})"
        print(f"scalar calibration: strut area A = {area:.3f} in^2 (K0 match to {label})")

    mode = f"length-regularized (Gf={Gf:g}, Gfc={Gfc:g})" if regularize else "plain (no regularization)"
    print(f"concrete law: Concrete02 {mode} | horizon={horizon:g}")

    model = rc_lattice(regularize, Gf, Gfc, area, horizon=horizon)
    print(f"RC lattice column: {len(model.nodes)} nodes, {len(model.elements)} struts")
    lat = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)
    print(f"lattice: peakV={max(lat['shear']):.2f} kip | drift->{lat['disp'][-1]:.2f} in "
          f"(conv={lat['converged']})")

    suffix = "" if calibration == "scalar" else f"_{calibration}"
    hz = "" if abs(horizon - HORIZON) < 1e-9 else f"_h{horizon:g}"
    savepath = outdir / f"column_pushover_{reference}{suffix}{hz}.png"
    viz.figure_pushover(
        [
            {"disp": ref["disp"], "shear": ref["shear"], "label": label,
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat["disp"], "shear": lat["shear"], "label": f"RC lattice ({calibration} calib.)",
             "style": {"color": "C0", "ls": "--", "lw": 2, "marker": "."}},
        ],
        savepath=str(savepath),
        xlabel="tip displacement (in)", ylabel="base shear (kip)",
        title=f"RC cantilever column pushover: lattice ({calibration} calib.) vs {label}",
    )
    print(f"saved pushover curve to {savepath}")

    # calibration output (D35): first N mode shapes (selected reference vs lattice) + periods table
    modalpath = outdir / f"column_modal_{reference}{suffix}{hz}.png"
    t_ref, t_lat = modal_calibration_figure(reference=reference, lattice_model=model, label=label,
                                            caption=caption, savepath=str(modalpath), Gf=Gf, Gfc=Gfc)
    print(f"modal calibration: {label} T={[f'{t:.4f}' for t in t_ref]} s | "
          f"lattice T={[f'{t:.4f}' for t in t_lat]} s -> saved {modalpath}")

    if draw:
        panels = [("RC lattice", model)]
        if reference == "continuum":
            panels.insert(0, ("2D continuum", _continuum_model(Gf, Gfc)[0]))
        drawpath = outdir / f"column_pushover_{reference}{suffix}{hz}_model.png"
        viz.figure_model(panels, savepath=str(drawpath),
                         suptitle="RC cantilever column — analysis model")
        print(f"saved model drawing to {drawpath}")


def _base_cut_groups(model, y_cut: float):
    """Force-decomposition probe (diagnosis) for `run_pushover`'s `element_groups`.

    For every strut crossing the horizontal cut at `y_cut` (one node below, one above), takes the
    GLOBAL force the element exerts on its LOWER node (`eleForce` index fx/fy of that node) — so the
    corotational/P-Δ geometry is included exactly. Free-body equilibrium of everything below the cut
    then gives, summed over the crossing struts:
      - "V:<cat>"  fx_lower            -> base shear (balances the base horizontal reactions),
      - "M:<cat>"  x_lower * fy_lower  -> overturning moment about the base center (≈ V*H).
    Category <cat> is concrete-vert / concrete-diag / rebar (longitudinal), separating the flexural
    couple (vertical struts + rebar) from the diagonal truss action the fiber section has no analog
    for. `eleForce` for a 2D truss returns [fx_i, fy_i, fx_j, fy_j].
    """
    steel = {m.id for m in model.uniaxial_materials if m.mtype == "Steel02"}
    groups: dict[str, list[tuple[int, int, float]]] = {}
    for el in model.elements:
        ni, nj = el.nodes
        (xi, yi), (xj, yj) = model.nodes[ni].coords, model.nodes[nj].coords
        if not (min(yi, yj) < y_cut < max(yi, yj)):
            continue
        fx, fy, x_lo = (0, 1, xi) if yi < yj else (2, 3, xj)   # eleForce indices at the lower node
        orient = "concrete-vert" if abs(xj - xi) < 1e-6 else "concrete-diag"
        cat = "rebar" if el.args[1] in steel else orient
        groups.setdefault(f"V:{cat}", []).append((el.id, fx, 1.0))
        groups.setdefault(f"M:{cat}", []).append((el.id, fy, x_lo))
    return groups


def diagnose(*, regularize: bool = True, Gf: float = GF, Gfc: float = GFC) -> None:
    """Decompose the lattice's base-cut resistance by element category to locate the ~17%
    overstrength vs the fiber beam-column, WITHOUT changing any mechanical parameter."""
    bc = beamcolumn_reference()
    k_bc = bc["shear"][1] / bc["disp"][1]
    area, ctrl, base = calibrate_area(k_bc)
    model = rc_lattice(regularize, Gf, Gfc, area)

    groups = _base_cut_groups(model, MESH / 2.0)   # cut between the base row and the first row up
    res = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base, element_groups=groups)

    sh, gp = res["shear"], res["groups"]
    ip = max(range(len(sh)), key=lambda i: sh[i])          # lattice peak step
    v = {k[2:]: gp[k][ip] for k in gp if k.startswith("V:")}
    m = {k[2:]: gp[k][ip] for k in gp if k.startswith("M:")}
    vtot, mtot = sum(v.values()), sum(m.values())

    print(f"\n=== base-cut decomposition at lattice peak: drift={res['disp'][ip]:.2f} in, "
          f"V={sh[ip]:.2f} kip ===")
    print("base shear (horizontal force across cut):")
    for k in sorted(v):
        print(f"  {k:16s} {v[k]:8.2f} kip  ({(v[k]/vtot*100 if vtot else 0):5.1f}%)")
    print(f"  {'cut total':16s} {vtot:8.2f} kip  (recorded base shear {sh[ip]:.2f})")
    print("overturning moment about base center (the flexural quantity, V = M/H):")
    for k in sorted(m):
        print(f"  {k:16s} {m[k]:10.0f} k-in ({(m[k]/mtot*100 if mtot else 0):5.1f}%)")
    print(f"  {'cut total':16s} {mtot:10.0f} k-in (V*H check = {sh[ip]*H:.0f})")
    print(f"fiber beam-column peak: V={max(bc['shear']):.2f} kip, M=V*H={max(bc['shear'])*H:.0f} k-in")
    diag = m.get("concrete-diag", 0.0)
    print(f"\n=> diagonal-strut truss action carries {abs(diag/mtot*100) if mtot else 0:.1f}% of the "
          f"base moment — the path the fiber section lacks (likely overstrength source).")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RC cantilever column pushover: lattice vs a selectable reference")
    parser.add_argument("--reference", choices=("beamcolumn", "continuum"), default="beamcolumn",
                        help="verification reference: fiber beam-column (fast, default) or 2D continuum (D29)")
    parser.add_argument("--calibration", choices=("scalar", "groups"), default="scalar",
                        help="lattice calibration: single K0 scalar (default) or the strong 2-group "
                             "orthogonal/diagonal fit to the continuum static+modal (D16)")
    parser.add_argument("--horizon", type=float, default=HORIZON,
                        help="strut connectivity horizon * mesh_size (default 1.5; larger = more "
                             "redundant bracing against the post-peak mechanism, D31)")
    parser.add_argument("--no-regularize", dest="regularize", action="store_false",
                        help="use plain Concrete02 instead of the crack-band length-regularized law (D20)")
    parser.add_argument("--gf", type=float, default=GF, help="tensile fracture energy (kip, in); regularized law only")
    parser.add_argument("--gfc", type=float, default=GFC, help="compressive fracture energy (kip, in); regularized law only")
    parser.add_argument("--diagnose", action="store_true",
                        help="decompose the lattice base-cut resistance by element category (overstrength diagnosis)")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model(s) (rebar highlighted)")
    args = parser.parse_args()
    if args.diagnose:
        diagnose(regularize=args.regularize, Gf=args.gf, Gfc=args.gfc)
    else:
        main(reference=args.reference, calibration=args.calibration, horizon=args.horizon,
             regularize=args.regularize, Gf=args.gf, Gfc=args.gfc, draw=args.draw)
