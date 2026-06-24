"""RC portal frame: lattice vs a selectable reference — nonlinear pushover (the frame analog of
examples/column/pushover.py).

The frame is two RC cantilever columns (24-in deep x 144-tall, confined CORE + unconfined COVER
Concrete02, 3+2+3 Steel02 bars + stirrups) connected by a thinner 18-in RC beam (same grades, top +
bottom bars + stirrups). It is modelled as a regularized RC lattice and pushed over, compared against
a SELECTABLE reference — `--reference {beamcolumn,continuum}`:
  - beamcolumn: a fiber `forceBeamColumn` portal frame (1D), the SAME Concrete02/Steel02 sections;
  - continuum: a 2D plane-stress continuum frame (D29), material-matched at the grade level.

Both: fixed bases, constant axial P at each column top, then a DisplacementControl lateral pushover
of the top-left joint. The lattice concrete strut area is calibrated so its initial lateral stiffness
equals the chosen reference's (same elastic slope). The lattice concrete softening law is CONFIGURABLE
(crack-band length-regularized Concrete02, D20, default; or plain) and so are Gf/Gfc/horizon. Units:
kip, in. Output: examples/output/frame/pushover/<reference>/frame_pushover_<reference>.png.

The specimen + builders are imported from this package (specimen.py / build.py), shared with the
linear and dynamic siblings. Run as `python examples/frame/pushover.py`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.opensees import run_pushover

from build import (
    _continuum_model, beamcolumn_reference, calibrate_area, calibrate_groups, make_reference, rc_lattice,
)
from specimen import DU, GF, GFC, HORIZON, MESH, OUT, TARGET, lateral_loads

REFERENCE_LABEL = {"beamcolumn": "fiber beam-column frame", "continuum": "2D continuum frame"}


def main(*, reference: str = "beamcolumn", calibration: str = "scalar", horizon: float = HORIZON,
         regularize: bool = True, Gf: float = GF, Gfc: float = GFC, beam_nonlinear: bool = False,
         draw: bool = False) -> None:
    outdir = OUT / "pushover" / reference
    outdir.mkdir(parents=True, exist_ok=True)

    ref = make_reference(reference, beam_nonlinear)
    label = REFERENCE_LABEL[reference]
    beam_mode = "nonlinear" if beam_nonlinear else "elastic"
    k_ref = ref["shear"][1] / ref["disp"][1]   # initial (small-drift) lateral stiffness
    print(f"beam concrete: {beam_mode} (default elastic — the softening lattice beam forms a local "
          f"mechanism the static pushover can't trace)")
    print(f"{label}: K0={k_ref:.2f} kip/in | peakV={max(ref['shear']):.2f} kip | "
          f"drift->{ref['disp'][-1]:.2f} in (conv={ref['converged']})")

    if calibration == "groups":
        area, ctrl, base, cal = calibrate_groups(horizon=horizon)   # strong 2-group fit (D16)
        o, d = cal.areas["orthogonal"], cal.areas["diagonal"]
        print(f"strong 2-group calibration (continuum static+modal): orthogonal A={o:.3f}, "
              f"diagonal A={d:.3f} in^2 (d/o ratio {d/o:.2f}, RMS {cal.rms:.3f}, success={cal.success})")
    else:
        area, ctrl, base = calibrate_area(k_ref, horizon=horizon)
        print(f"scalar calibration: strut area A = {area:.3f} in^2 (K0 match to {label})")

    mode = f"length-regularized (Gf={Gf:g}, Gfc={Gfc:g})" if regularize else "plain (no regularization)"
    print(f"concrete law: Concrete02 {mode} | horizon={horizon:g}")

    model = rc_lattice(regularize, Gf, Gfc, area, beam_nonlinear=beam_nonlinear, horizon=horizon)
    print(f"RC lattice frame: {len(model.nodes)} nodes, {len(model.elements)} struts")
    lat = run_pushover(model, lateral_loads=lateral_loads(model), control_node=ctrl,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)
    print(f"lattice: peakV={max(lat['shear']):.2f} kip | drift->{lat['disp'][-1]:.2f} in "
          f"(conv={lat['converged']})")

    suffix = "" if calibration == "scalar" else f"_{calibration}"
    nb = "_nlbeam" if beam_nonlinear else ""
    hz = "" if abs(horizon - HORIZON) < 1e-9 else f"_h{horizon:g}"
    savepath = outdir / f"frame_pushover_{reference}{suffix}{nb}{hz}.png"
    viz.figure_pushover(
        [
            {"disp": ref["disp"], "shear": ref["shear"], "label": label,
             "style": {"color": "C3", "lw": 2}},
            {"disp": lat["disp"], "shear": lat["shear"], "label": f"RC lattice ({calibration} calib.)",
             "style": {"color": "C0", "ls": "--", "lw": 2, "marker": "."}},
        ],
        savepath=str(savepath),
        xlabel="roof displacement (in)", ylabel="base shear (kip)",
        title=f"RC portal frame pushover ({beam_mode} beam): lattice ({calibration} calib.) vs {label}",
    )
    print(f"saved pushover curve to {savepath}")

    if draw:
        panels = [("RC lattice", model)]
        if reference == "continuum":
            panels.insert(0, ("2D continuum", _continuum_model(Gf, Gfc)[0]))
        drawpath = outdir / f"frame_pushover_{reference}{suffix}{nb}{hz}_model.png"
        viz.figure_model(panels, savepath=str(drawpath),
                         suptitle="RC portal frame — analysis model")
        print(f"saved model drawing to {drawpath}")


def _base_cut_groups(model, y_cut: float):
    """Force-decomposition probe (diagnosis) for `run_pushover`'s `element_groups`.

    For every strut crossing the horizontal cut at `y_cut` (one node below, one above) — across BOTH
    columns — takes the GLOBAL force the element exerts on its LOWER node (`eleForce` fx/fy of that
    node, so the corotational/P-Delta geometry is exact). Free-body equilibrium below the cut then
    gives, summed over the crossing struts: "V:<cat>" base shear and "M:<cat>" overturning moment
    about x=0. Category <cat> is concrete-vert / concrete-diag / rebar (longitudinal)."""
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


def diagnose(*, regularize: bool = True, Gf: float = GF, Gfc: float = GFC,
             beam_nonlinear: bool = False) -> None:
    """Decompose the lattice frame's base-cut resistance by element category to locate the
    diagonal-strut overstrength vs the fiber-frame reference, WITHOUT changing any parameter."""
    bc = beamcolumn_reference(beam_nonlinear)
    k_bc = bc["shear"][1] / bc["disp"][1]
    area, ctrl, base = calibrate_area(k_bc)
    model = rc_lattice(regularize, Gf, Gfc, area, beam_nonlinear=beam_nonlinear)

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
    print("overturning moment about x=0:")
    for k in sorted(m):
        print(f"  {k:16s} {m[k]:10.0f} k-in ({(m[k]/mtot*100 if mtot else 0):5.1f}%)")
    print(f"fiber-frame peak: V={max(bc['shear']):.2f} kip")
    diag = m.get("concrete-diag", 0.0)
    print(f"\n=> diagonal-strut truss action carries {abs(diag/mtot*100) if mtot else 0:.1f}% of the "
          f"base moment — the path the fiber frame lacks (likely overstrength source).")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RC portal frame pushover: lattice vs a selectable reference")
    parser.add_argument("--reference", choices=("beamcolumn", "continuum"), default="beamcolumn",
                        help="verification reference: fiber beam-column frame (fast, default) or 2D continuum (D29)")
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
    parser.add_argument("--nonlinear-beam", dest="beam_nonlinear", action="store_true",
                        help="model the thin beam concrete with the SAME nonlinear Concrete02 as the "
                             "columns (default: elastic beam — the softening lattice beam forms a local "
                             "mechanism the static pushover can't trace; the dynamic runs are stable)")
    parser.add_argument("--diagnose", action="store_true",
                        help="decompose the lattice base-cut resistance by element category (overstrength diagnosis)")
    parser.add_argument("--draw", action="store_true",
                        help="also save a drawing of the analysis model(s) (rebar highlighted)")
    args = parser.parse_args()
    if args.diagnose:
        diagnose(regularize=args.regularize, Gf=args.gf, Gfc=args.gfc, beam_nonlinear=args.beam_nonlinear)
    else:
        main(reference=args.reference, calibration=args.calibration, horizon=args.horizon,
             regularize=args.regularize, Gf=args.gf, Gfc=args.gfc, beam_nonlinear=args.beam_nonlinear,
             draw=args.draw)
