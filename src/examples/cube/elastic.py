"""Concrete cube, uniaxial tension — LINEAR-ELASTIC baseline: lattice vs a single dispBeamColumn.

The foundational verification, before any cracking/softening: model ONE 0.1 m concrete cube two
ways and confirm the lattice reproduces the beam-column's uniaxial stress-strain line (slope =
Young's modulus E) in the elastic range. No regularization, no dynamic solver — just the elastic
modulus, matched.

Two materials are selectable (`--material`), and BOTH stay in the elastic branch:
  - `elastic`    (default) — a linear `Elastic` material at E. The pure-elastic D42 baseline; the
                 stress-strain line is straight forever, traced over the full range.
  - `concrete02` — the real `Concrete02` law, whose TENSION branch is linear-elastic up to the
                 tensile strength ft at the cracking strain eps_cr = ft/E, then softens. This mode
                 pulls ONLY to eps_cr, so the comparison stays on the elastic branch and shows that
                 Concrete02's initial branch reproduces the SAME E-line as the beam-column, up to ft
                 (the point where cracking — the next modelling step — would begin).

Either way:
  - reference: a single displacement-based fiber `dispBeamColumn` -> the clean 1D coupon (sigma = E*eps);
  - lattice:   the same specimen as struts (gmsh nodes + horizon struts), with ONE uniform strut area
               calibrated so the lattice's initial axial stiffness equals the beam-column's E*A/L.

WHY A CALIBRATED AREA IS NEEDED: a horizon-strut lattice sums the axial stiffness of EVERY strut
(orthogonal + the diagonal struts that brace it — an orthogonal-only lattice is a singular mechanism,
the D4 caveat) in parallel, so a physically-sized strut (area = tributary MESH*THK = 1000 mm^2) makes
the assembly ~1.58x too stiff. The one knob this baseline uses is the single uniform strut area, tuned
so the assembly's K0 = the beam-column's E*A/L. With that, both curves lie exactly on sigma = E*eps.
Fully fixed base, top edge pulled +Y.

Output: examples/output/cube/elastic/cube_tension_elastic[_concrete02].png. Units: N, mm.
Run as `python examples/cube/elastic.py [--material elastic|concrete02] [--draw]`.
"""

from __future__ import annotations

from rclattice import viz
from rclattice.builders import build_lattice, build_lattice_rc, select_nodes
from rclattice.materials import concrete_uniaxial_elastic, concrete_uniaxial_nonlinear
from rclattice.opensees import run_dispbeamcolumn_tension, run_pushover

from specimen import (
    A_FACE, CONCRETE, DU, EPS, HORIZON, L, MESH, OUT, THK, TARGET, axial_loads, cube_problem,
)

OUT_ELASTIC = OUT / "elastic"

# Concrete02 tension is linear-elastic up to ft at the cracking strain eps_cr = ft/E, then softens.
# The concrete02 mode pulls only to (just shy of) this elastic limit, so the whole comparison stays
# on the elastic branch. (For this grade: eps_cr = 3/30000 = 1e-4, u_cr = 0.01 mm.)
EPS_CR = CONCRETE.ft / CONCRETE.E
U_CR = EPS_CR * L


def _concrete_material(material: str, tag: int):
    """The uniaxial concrete material for the chosen mode: linear `Elastic` (at E) or `Concrete02`
    (whose tension branch is elastic to ft, exercised here only up to eps_cr)."""
    if material == "concrete02":
        return concrete_uniaxial_nonlinear(CONCRETE, tag)
    return concrete_uniaxial_elastic(CONCRETE, tag)


def beamcolumn_reference(material: str, target: float) -> dict:
    """Single displacement-based fiber beam-column in uniaxial tension — the 1D reference. The fiber
    section is the cube face (in-plane L x out-of-plane THK), so its axial stiffness is exactly E*A/L.
    Returns the tension pushover curve {"disp","shear","converged"} (shear = axial base reaction, N)."""
    return run_dispbeamcolumn_tension(length=L, width=THK, depth=L,
                                      material=_concrete_material(material, 1), dU=DU, target=target)


def calibrate_area(k_ref: float, *, horizon: float = HORIZON):
    """SCALAR calibration: the single uniform strut area giving the ELASTIC lattice the same initial
    axial stiffness as `k_ref`. K is linear in strut area, so one tiny tension solve fixes the scale.
    The elastic builder is used regardless of the study material — Concrete02's INITIAL tangent is the
    same E, so the calibrated area matches K0 in either mode. Returns (area, control_node, base_nodes)."""
    lat0, _ = build_lattice(cube_problem(), MESH, strut_area=1.0, horizon=horizon)
    ctrl = select_nodes(lat0, (L / 2.0 - EPS, L / 2.0 + EPS, L - EPS, L + EPS))[0]
    base = select_nodes(lat0, (-EPS, L + EPS, -EPS, EPS))
    r = run_pushover(lat0, lateral_loads=axial_loads(lat0), control_node=ctrl,
                     control_dof=2, dU=1e-3, target=1e-3, base_nodes=base)
    return k_ref / (r["shear"][-1] / r["disp"][-1]), ctrl, base


def cube_lattice(area: float, material: str, *, horizon: float = HORIZON):
    """The calibrated cube lattice (uniform `area`, no rebar). `material` = 'elastic' (Elastic struts,
    via `build_lattice`) or 'concrete02' (Concrete02 `Truss` struts, one plain-concrete zone, via
    `build_lattice_rc`)."""
    if material == "elastic":
        model, _ = build_lattice(cube_problem(), MESH, strut_area=area, horizon=horizon)
        return model
    model, _ = build_lattice_rc(
        cube_problem(), MESH,
        material_for=lambda zone, length: concrete_uniaxial_nonlinear(CONCRETE, 0),
        zone_of=lambda x, y: "concrete", rebars=(), strut_area=area,
        horizon=horizon, strut_element="Truss",
    )
    return model


def main(*, material: str = "elastic", draw: bool = False) -> None:
    outdir = OUT_ELASTIC
    outdir.mkdir(parents=True, exist_ok=True)
    # concrete02 pulls just shy of the cracking strain (stay strictly in the elastic branch);
    # elastic can trace the full range (the line is straight forever).
    target = 0.98 * U_CR if material == "concrete02" else TARGET

    ref = beamcolumn_reference(material, target)
    k_ref = ref["shear"][1] / ref["disp"][1]                       # E*A/L, N/mm (initial tangent)
    e_ref = (ref["shear"][1] / A_FACE) / (ref["disp"][1] / L)      # stress/strain slope, MPa
    print(f"dispBeamColumn ({material}): K0 = {k_ref:.4e} N/mm | effective E = {e_ref:.1f} MPa "
          f"(grade E = {CONCRETE.E:.0f}) | traced to eps={target / L:.2e} | conv={ref['converged']}")

    area, ctrl, base = calibrate_area(k_ref)
    a_phys = MESH * THK
    print(f"scalar K0 calibration: strut area A = {area:.2f} mm^2  "
          f"(physical tributary MESH*THK = {a_phys:.0f} mm^2 -> lattice would be {a_phys / area:.3f}x too stiff)")

    model = cube_lattice(area, material)
    lat = run_pushover(model, lateral_loads=axial_loads(model), control_node=ctrl,
                       control_dof=2, dU=DU, target=target, base_nodes=base)
    e_lat = (lat["shear"][1] / A_FACE) / (lat["disp"][1] / L)
    sec = [(s / A_FACE) / (u / L) for u, s in zip(lat["disp"][1:], lat["shear"][1:])]
    strut_kind = "Concrete02" if material == "concrete02" else "Elastic"
    print(f"cube lattice ({material}): {len(model.nodes)} nodes, {len(model.elements)} {strut_kind} struts | "
          f"effective E = {e_lat:.1f} MPa | secant-E spread over the elastic pull = "
          f"{max(sec) - min(sec):.3f} MPa ({'LINEAR' if max(sec) - min(sec) < 1.0 else 'nonlinear'}) | conv={lat['converged']}")
    print(f"MATCH: E_lattice / E_ref = {e_lat / e_ref:.5f}  (|error| = {abs(e_lat - e_ref) / e_ref * 100:.3f}%)")

    # coupon stress-strain: axial stress (MPa) = base reaction / cube-face area, strain = u / edge
    ref_lbl = f"dispBeamColumn ({strut_kind}, E={e_ref:.0f} MPa)"
    lat_lbl = f"{strut_kind} lattice (calibrated A={area:.0f} mm², E={e_lat:.0f} MPa)"
    curves = [
        {"disp": [u / L for u in ref["disp"]], "shear": [s / A_FACE for s in ref["shear"]],
         "label": ref_lbl, "style": {"color": "C3", "lw": 2}},
        {"disp": [u / L for u in lat["disp"]], "shear": [s / A_FACE for s in lat["shear"]],
         "label": lat_lbl, "style": {"color": "C0", "ls": "--", "lw": 2, "marker": ".", "markevery": 8}},
    ]
    branch = "elastic branch to ft" if material == "concrete02" else "linear elastic"
    stem = "cube_tension_elastic" + ("_concrete02" if material == "concrete02" else "")
    savepath = outdir / f"{stem}.png"
    viz.figure_pushover(curves, savepath=str(savepath),
                        xlabel="axial strain (-)", ylabel="axial stress (MPa)",
                        title=f"Concrete cube, uniaxial tension — {branch}: lattice vs single dispBeamColumn")
    print(f"saved stress-strain curve to {savepath}")

    if draw:
        drawpath = outdir / f"{stem}_model.png"
        viz.figure_model([(f"cube lattice ({material})", model)], savepath=str(drawpath),
                         suptitle=f"Concrete cube — {material} lattice analysis model")
        print(f"saved model drawing to {drawpath}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Concrete cube uniaxial tension, LINEAR-ELASTIC baseline: "
                                                 "lattice vs single dispBeamColumn")
    parser.add_argument("--material", choices=("elastic", "concrete02"), default="elastic",
                        help="strut/fiber material: linear 'elastic' (default, the D42 pure-elastic "
                             "baseline) or 'concrete02' (the real law, pulled only to the cracking "
                             "strain eps_cr=ft/E so only its elastic branch is exercised)")
    parser.add_argument("--draw", action="store_true", help="also save a drawing of the lattice model")
    args = parser.parse_args()
    main(material=args.material, draw=args.draw)
