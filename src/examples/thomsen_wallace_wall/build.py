"""Model builders + energy-balance calibration for Thomsen & Wallace's RW2.

Same recipe as the other wall studies — strut areas from Aydin's overlapping-lattice energy balance
(D47, `field="uniaxial"`), tension softening regularized by strut length (D20) — with one addition
this specimen is the reason for: a GRADED grid whose lines are the bar axes (`grid="rebar"`, D104).

HOW THE GRADED GRID IS CALIBRATED. The balance is struck on a UNIFORM patch at the target mesh,
exactly as before, and gives one `A_t` for that spacing. On the graded grid every strut then takes
`A_t` scaled by its own tributary width (`mesh.tributary_area_scale`: an orthogonal strut by the
mean perpendicular spacing beside it, a diagonal by the isotropic scale of its cell), which is what
Aydin's closed form `E_t A_t = C E_t d w` says a strut of local spacing `d` should carry. Under a
uniform rescaling of the grid this reduces exactly to the uniform rule; on this grid the spacings
run 19.0-25.5 mm against the 25 mm target, so the scaling stays within 0.76-1.02. The Stage 0 elastic
gate — lattice against a plane-stress continuum on the SAME graded grid — is what verifies it.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

from rclattice.builders import build_continuum_rc, build_lattice_rc
from rclattice.calibration import EnergyBalanceResult, aydin_closed_form_C, energy_balance_rectangle
from rclattice.materials import (
    bond_elastic_brittle,
    bond_elastic_brittle_damaging,
    concrete_lattice_aydin,
    concrete_nd_elastic_planestress,
    concrete_uniaxial_elastic,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
    steel_uniaxial_ruptured,
)
from rclattice.mesh import (connect_index_horizon, mesh_rectangle_lines, tributary_area_scale)

from specimen import (
    EC, FT, GF, GRADES, GRID, HORIZON, HW, LW, MESH, N_AXIAL, NU, TW, X_BARS, Y_HOOPS,
    check_mesh_alignment, grid_lines, rebars, wall_problem, zone_of,
)


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON,
              field: str = "uniaxial") -> EnergyBalanceResult:
    """The energy balance over a wall-sized patch of a UNIFORM grid at `mesh_size`.

    The result is a property of (mesh_size, horizon, nu, thickness, field), not of this rectangle;
    on the graded grid it is scaled per strut (module docstring).
    """
    return energy_balance_rectangle(LW, HW, mesh_size, E=EC, nu=NU, thickness=TW,
                                    horizon=horizon, field=field)


def report_calibration(cal: EnergyBalanceResult, *, mesh_size: float = MESH,
                       horizon: float = HORIZON, grid: str = GRID) -> None:
    nominal = mesh_size * TW
    published = aydin_closed_form_C(horizon)
    print(f"Energy balance [{cal.field}] (horizon {horizon}, nu = {cal.nu:.2f}): "
          f"A_t = {cal.area:,.1f} mm^2 = {cal.area / nominal:.4f} * (thickness * mesh)")
    print(f"  EA = {cal.EA:.4g} N   ->  EA/L = {cal.EA / mesh_size:.4g} N/mm (orthogonal), "
          f"{cal.EA / (mesh_size * 2 ** 0.5):.4g} N/mm (diagonal) at the target spacing")
    if grid == "rebar":
        print(f"  GRADED GRID: each strut takes A_t x (local spacing / {mesh_size:g}) — "
              f"EA/L is held per family, not A (D104)")
    print(f"  paper's published closed form: C = {published:.4f} -> A_t = "
          f"{published * nominal:,.1f} mm^2 at ITS nu = 1/3 "
          f"({published * nominal / cal.area:.3f} x this run)")
    print(f"  nu_consistent = {cal.nu_consistent:.3f} (NOT a Poisson ratio — nu_effective = "
          f"{cal.nu_effective:.3f}, cubic anisotropy {cal.cubic_anisotropy:.2f})")
    print(f"  shear-stiffness error at nu = {cal.nu:.2f}: {cal.isotropy_error * 100:.2f}%  "
          f"— shear is only ~7% of this wall's elastic flexibility (aspect 3.0)")


def grid_for(grid: str = GRID, mesh_size: float = MESH, *, length: float = LW,
             height: float = HW, horizon: float = HORIZON):
    """`(coords, quads, pairs)` for either grid mode — None for the uniform default, which the
    builders generate themselves; the graded triple otherwise."""
    if grid == "uniform":
        check_mesh_alignment(mesh_size, length=length, height=height)
        return None
    if grid != "rebar":
        raise SystemExit(f"unknown grid {grid!r}; use 'uniform' or 'rebar'")
    coords, quads = mesh_rectangle_lines(length, height, mesh_size, x_lines=X_BARS, y_lines=Y_HOOPS)
    pairs = connect_index_horizon(coords, horizon)
    return coords, quads, pairs


def _bars(grid, mesh_size, length, height, full_height, reinforced, steel_b):
    if not reinforced:
        return ()
    bars = rebars(grid, mesh_size, length=length, height=height, full_height=full_height)
    if steel_b is None:
        return bars
    return tuple(dataclasses.replace(
        r, steel=dataclasses.replace(r.steel, name=f"{r.steel.name}-b{steel_b:g}", b=float(steel_b)))
        for r in bars)


def wall_lattice(
    area: float,
    *,
    grid: str = GRID,
    mesh_size: float = MESH,
    horizon: float = HORIZON,
    nonlinear: bool = True,
    gf_factor: float = 1.0,
    strut_element: str = "corotTruss",
    length: float = LW,
    height: float = HW,
    material: str = "concrete02",
    reinforced: bool = True,
    full_height_rebar: bool = True,
    tail_a2a3: "tuple[float, float] | None" = None,
    fc_cap: float | None = None,
    bond: bool = False,
    bond_residual: float = 0.7,
    bond_damage: bool = False,
    steel_rupture: float | None = None,
    concrete_residual: float = 0.2,
    steel_b: float | None = None,
):
    """The calibrated RC lattice on either grid. `nonlinear=False` gives the elastic twin.

    The material switches are the Aldemir study's (D86/D87/D91/D95): `material` is "concrete02"
    (crushing) or "aydin" (his trilinear tension tail, compression linear forever unless `fc_cap`
    caps it); `steel_rupture` wraps Steel02 in MinMax; `concrete_residual` lowers the crushing floor
    ON THE GRADE (D91); `bond` gives each bar its own nodes tied by a full-area ring with Aydin's
    elastic-brittle law at residual `bond_residual`, irreversible when `bond_damage`.
    """
    if nonlinear and material == "aydin":
        cap = {} if fc_cap is None else {"fc_cap": fc_cap, "rsm": False}
        material_for = lambda zone, L: concrete_lattice_aydin(            # noqa: E731
            GRADES[zone], 0, L, Gf=GF * gf_factor, a3_over_a2=350.0 / 80.0, **cap,
            **({} if tail_a2a3 is None else {"a2": tail_a2a3[0], "a3": tail_a2a3[1]}))
        rebar_material = steel_uniaxial
    elif nonlinear:
        crush_grade = {z: dataclasses.replace(g, fcu=concrete_residual * g.fc)
                       for z, g in GRADES.items()} if concrete_residual < 0.2 else GRADES
        material_for = lambda zone, L: concrete_uniaxial_regularized(     # noqa: E731
            crush_grade[zone], 0, L, Gf=GF * gf_factor, residual_ratio=concrete_residual)
        rebar_material = steel_uniaxial
    else:
        material_for = lambda zone, _L: concrete_uniaxial_elastic(GRADES[zone], 0)  # noqa: E731
        rebar_material = steel_uniaxial_elastic

    if nonlinear and steel_rupture is not None:
        rebar_material = lambda grade, tag: steel_uniaxial_ruptured(      # noqa: E731
            grade, tag, eps_rupture=steel_rupture)

    extra = {}
    if bond:
        _law = bond_elastic_brittle_damaging if bond_damage else bond_elastic_brittle
        extra = dict(bond_material=lambda zone, _L: _law(GRADES[zone], 0, residual=bond_residual),
                     bond_horizon=horizon, i_accept_the_known_bond_defect=True)

    g = grid_for(grid, mesh_size, length=length, height=height, horizon=horizon)
    if g is not None:
        coords, quads, pairs = g
        scale = tributary_area_scale(coords, pairs, mesh_size)
        by_pair = {p: float(area * s) for p, s in zip(pairs, scale)}
        extra.update(grid=(coords, quads), pairs=pairs,
                     strut_area_of_pair=lambda i, j, _L: by_pair[(i, j)])

    model, edges = build_lattice_rc(
        wall_problem(length=length, height=height), mesh_size,
        material_for=material_for, zone_of=zone_of,
        rebars=_bars(grid, mesh_size, length, height, full_height_rebar, reinforced, steel_b),
        strut_area=area, horizon=horizon, strut_element=strut_element,
        rebar_material=rebar_material, **extra)
    return model, edges


def wall_continuum(*, grid: str = GRID, mesh_size: float = MESH, length: float = LW,
                   height: float = HW):
    """The SAME wall as an elastic plane-stress CONTINUUM on the SAME grid — the Stage 0 reference.

    `K_lattice / K_continuum` isolates one thing: whether the (scaled) energy-balance areas make the
    lattice an equivalent elastic continuum on this grid, which is the only claim the calibration
    makes. Beam theory is a separate, weaker check (`cantilever_stiffness`).
    """
    g = grid_for(grid, mesh_size, length=length, height=height)
    model, edges = build_continuum_rc(
        wall_problem(length=length, height=height), mesh_size,
        nd_material_for=lambda zone: concrete_nd_elastic_planestress(GRADES[zone], 0),
        zone_of=zone_of,
        rebars=rebars(grid, mesh_size, length=length, height=height),
        rebar_material=steel_uniaxial_elastic,
        grid=None if g is None else (g[0], g[1]),
    )
    return model, edges


def describe(model) -> str:
    kinds: dict[str, int] = {}
    for e in model.elements:
        kinds[e.kind] = kinds.get(e.kind, 0) + 1
    parts = ", ".join(f"{n:,} {k}" for k, n in sorted(kinds.items()))
    return f"{len(model.nodes):,} nodes, {len(model.elements):,} elements ({parts})"


# --- closed-form elastic references -----------------------------------------------------------------

def gross_inertia(*, length: float = LW, thickness: float = TW) -> float:
    return thickness * length ** 3 / 12.0


def transformed_inertia(grid: str = GRID, mesh_size: float = MESH, *, length: float = LW,
                        height: float = HW, thickness: float = TW) -> float:
    """Gross concrete I plus the vertical bars' (n-1)*A*x^2 about the section centroid."""
    inertia = gross_inertia(length=length, thickness=thickness)
    centroid = length / 2.0
    for rb in rebars(grid, mesh_size, length=length, height=height):
        (xa, _ya), (xb, _yb) = rb.path[-2], rb.path[-1]
        if abs(xb - xa) > 1e-9:
            continue
        inertia += (rb.steel.E0 / EC - 1.0) * rb.area * (xb - centroid) ** 2
    return inertia


def cantilever_stiffness(*, shear_span: float, inertia: float, length: float = LW,
                         thickness: float = TW, kappa: float = 1.2) -> tuple[float, float]:
    """Fixed-base elastic tip stiffness: `(K in N/mm, shear share of the flexibility)`."""
    g = EC / (2.0 * (1.0 + NU))
    flex = shear_span ** 3 / (3.0 * EC * inertia)
    shear = kappa * shear_span / (g * thickness * length)
    return 1.0 / (flex + shear), shear / (flex + shear)


def strut_life(gf_factor: float = 1.0, mesh_size: float = MESH, *, diagonal: bool = False,
               grid: str = GRID) -> float:
    """`eps_ult / eps_cr` for a regularized strut — how far past cracking it carries tension.

    On the graded grid the LONGEST strut of each family is the brittle one (crack-band: the
    softening slope grows with L), so the life is quoted at the largest spacing the grid has.
    """
    if grid == "rebar":
        xs, ys = grid_lines(grid, mesh_size)
        dmax = max(max(np.diff(xs)), max(np.diff(ys)))
        L = math.hypot(max(np.diff(xs)), max(np.diff(ys))) if diagonal else dmax
    else:
        L = mesh_size * (2.0 ** 0.5 if diagonal else 1.0)
    ets = FT * FT * L / (2.0 * GF * gf_factor)
    return (FT / EC + FT / ets) / (FT / EC)


def cracking_shear(*, length: float = LW, height: float = HW, grid: str = GRID,
                   mesh_size: float = MESH) -> tuple[float, float]:
    """(V_cr, drift at V_cr) with the axial load OFFSETTING the tensile stress: M_cr = (f_t + N/A) S."""
    i_g = gross_inertia(length=length)
    s_mod = i_g / (length / 2.0)
    m_cr = (FT + N_AXIAL / (TW * length)) * s_mod
    v_cr = m_cr / height
    k_tr, _share = cantilever_stiffness(
        shear_span=height, inertia=transformed_inertia(grid, mesh_size, length=length, height=height),
        length=length)
    return v_cr, (v_cr / k_tr) / height


def strut_groups(model, *, quantity: str = "shear", cut_y: float | None = None,
                 mesh_size: float = MESH, tol: float = 1e-6) -> dict:
    """Free-body probe across a horizontal CUT, split by element orientation (Aldemir `build.py`).

    `dof` indexes `eleForce` 0-based; `coef` carries -1 for shear because `eleForce` is the force
    the element exerts ON its node. Bond links get their own group (D85).
    """
    if quantity not in ("shear", "moment"):
        raise ValueError(f"quantity must be 'shear' or 'moment', got {quantity!r}")
    cut = 0.5 * mesh_size if cut_y is None else cut_y
    centroid = LW / 2.0
    groups: dict[str, list] = {"vertical": [], "horizontal": [], "diagonal": [], "rebar": [],
                               "bond": []}
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        if (ya - cut) * (yb - cut) >= 0.0:
            continue
        below_first = ya < cut
        x = xa if below_first else xb
        if quantity == "shear":
            dof, coef = (0 if below_first else 2), -1.0
        else:
            dof, coef = (1 if below_first else 3), (x - centroid)
        if e.kind in ("longitudinal", "stirrup"):
            label = "rebar"
        elif e.kind == "bond":
            label = "bond"
        elif abs(xb - xa) < tol:
            label = "vertical"
        elif abs(yb - ya) < tol:
            label = "horizontal"
        else:
            label = "diagonal"
        groups[label].append((e.id, dof, coef))
    return {k: v for k, v in groups.items() if v}
