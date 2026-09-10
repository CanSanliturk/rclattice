"""Model builders + energy-balance calibration for the Aldemir et al. (2017) squat wall.

Same two-part recipe as the other three wall studies — strut areas from Aydin's overlapping-lattice
energy balance, tension softening regularized by strut length — with two switches this specimen is
the reason for.

WHICH ENERGY BALANCE (`--calibration`). The repo has used the thesis route (uniaxial strain with the
transverse strain restrained, D47) for SW-NC-FF, WSH3 and VK3. The 2019 paper that simulates THIS
wall publishes a different one: an equal-stress equibiaxial field, whose closed form is the
`Et*At = 0.621*Et*d*w` its Appendix prints for horizon 1.5d. They disagree by ~18% in EA at matched
nu, so on the specimen the 2019 paper itself over-predicts by 21% it is worth running both. Note the
two errors nearly cancel at the repo's working nu: `--calibration uniaxial` at nu = 0.20 lands within
3% of the paper's published closed form at its own nu = 1/3.

BOND (`--bond`) — **BROKEN, kept only to reproduce the defect.** The elastic run gives 2.19x the
transformed-section stiffness with a 12.4% secant spread in a LINEAR model; bond can only ever add
flexibility. The horizon-star links carry the concrete strut's EA and so form a parallel load path
beside a concrete lattice that is still fully present, and the topology cannot satisfy
"rigid bond == perfect bond" at any stiffness (D72). The fix is a coincident-node zeroLength
interface spring. Until then `--bond` produces numbers that are not usable.

Why it matters for THIS specimen: every other wall in the repo has a documented reason for perfect
bond (WSH3 and VK3 have continuous deformed bars) or a documented reason it fails (SW-NC-FF's plain
bars, which capped that study's scope). Here the bar type is simply not reported, so perfect bond is
an assumption with nothing behind it — which is why a working bond option is worth having.
"""

from __future__ import annotations

from rclattice.builders import build_lattice_rc
from rclattice.calibration import EnergyBalanceResult, aydin_closed_form_C, energy_balance_rectangle
from rclattice.builders import build_continuum_rc
from rclattice.materials import (
    bond_elastic_brittle,
    bond_elastic_brittle_damaging,
    steel_uniaxial_ruptured,
    concrete_lattice_aydin,
    concrete_nd_elastic_planestress,
    concrete_uniaxial_elastic,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
)

from specimen import (
    EC, FT, GF, GRADES, HORIZON, HW, LW, MESH, NU, TW, check_mesh_alignment, rebars, wall_problem,
    zone_of,
)

BOND_RESIDUAL = 0.7      # Aydin's `a`, calibrated in Aydin (2017) for DEFORMED bars


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON,
              field: str = "uniaxial") -> EnergyBalanceResult:
    """The energy balance over a wall-sized patch of the working grid.

    `field` is "uniaxial" (the repo's D47 route) or "equibiaxial" (the 2019 paper's published one).
    The result is a property of (mesh_size, horizon, nu, thickness, field), not of this rectangle.
    """
    return energy_balance_rectangle(LW, HW, mesh_size, E=EC, nu=NU, thickness=TW,
                                    horizon=horizon, field=field)


def report_calibration(cal: EnergyBalanceResult, *, mesh_size: float = MESH,
                       horizon: float = HORIZON) -> None:
    """Print the calibrated area, the EA it gives the struts, and BOTH routes side by side."""
    nominal = mesh_size * TW
    published = aydin_closed_form_C(horizon)
    print(f"Energy balance [{cal.field}] (horizon {horizon}, nu = {cal.nu:.2f}): "
          f"A_t = {cal.area:,.1f} mm^2 = {cal.area / nominal:.4f} * (thickness * mesh)")
    print(f"  EA = {cal.EA:.4g} N   ->  EA/L = {cal.EA / mesh_size:.4g} N/mm (orthogonal), "
          f"{cal.EA / (mesh_size * 2 ** 0.5):.4g} N/mm (diagonal)")
    print(f"  the OTHER route would give A_t = "
          f"{cal.area_equibiaxial if cal.field == 'uniaxial' else cal.area_x:,.1f} mm^2 "
          f"({(cal.area_equibiaxial if cal.field == 'uniaxial' else cal.area_x) / cal.area:.3f} x)")
    print(f"  paper's published closed form: C = {published:.4f} -> A_t = "
          f"{published * nominal:,.1f} mm^2 at ITS nu = 1/3 "
          f"({published * nominal / cal.area:.3f} x this run)")
    print(f"  nu_consistent = {cal.nu_consistent:.3f} (NOT a Poisson ratio — nu_effective = "
          f"{cal.nu_effective:.3f}, cubic anisotropy {cal.cubic_anisotropy:.2f})")
    print(f"  shear-stiffness error at nu = {cal.nu:.2f}: {cal.isotropy_error * 100:.2f}%  "
          f"— READ THIS ONE: this wall's aspect ratio is 0.75 and it was designed to yield in shear")


def _bond_material(bond_area: float, strut_area: float):
    """Bond law factory matching `material_for`'s (zone, length) signature (D75).

    Sets the link's STRENGTH and STIFFNESS independently, which a single area cannot do:

      * strength — `peak_force` is the concrete strut's own cracking force `ft * A_t`, because his
        Fig. 1(c) draws Bond and Concrete rising to the SAME `F_cr`;
      * stiffness — `slip_at_peak` is chosen so the link's modulus comes out at the concrete's `E`,
        which is exactly what the area calibration assumed when it was fitted for elastic
        equivalence. So the two constraints are satisfied together rather than traded off.

    The slip that falls out, ~0.37 mm over a 50 mm link, is where a deformed-bar interface actually
    fails (Model Code s1 ~ 0.1-1 mm). Inheriting the concrete grade's own cracking strain instead put
    it at 0.004 mm and disconnected the reinforcement at 0.026% drift.
    """
    def make(zone, length):
        grade = GRADES[zone]
        f_cr = grade.ft * strut_area
        slip = (f_cr / bond_area) / grade.E * length      # -> E_bond == grade.E
        return bond_elastic_brittle(grade, 0, residual=BOND_RESIDUAL, peak_force=f_cr,
                                    area=bond_area, length=length, slip_at_peak=slip)
    return make


def _bars(mesh_size, length, height, full_height, reinforced, steel_b):
    """The bar set, with an optional override of the steel hardening ratio `b`.

    `b` is NOT a number the paper gives — Table 1 prints f_y = 360 and nothing about hardening, so
    the repo's 0.01 is a convention, exactly like eps_su (D91). `steel_b=0.0` makes the ties
    elastic-perfectly-plastic, which is the controlled way to ask how much of a flat cyclic
    envelope is strain hardening rather than an absent failure mechanism.
    """
    import dataclasses

    if not reinforced:
        return ()
    bars = rebars(mesh_size, length=length, height=height, full_height=full_height)
    if steel_b is None:
        return bars
    return tuple(dataclasses.replace(
        r, steel=dataclasses.replace(r.steel, name=f"{r.steel.name}-b{steel_b:g}", b=float(steel_b)))
        for r in bars)


def wall_lattice(
    area: float,
    *,
    mesh_size: float = MESH,
    horizon: float = HORIZON,
    nonlinear: bool = True,
    gf_factor: float = 1.0,
    strut_element: str = "corotTruss",
    bond: bool = False,
    bond_area: float | None = None,
    length: float = LW,
    height: float = HW,
    material: str = "concrete02",
    reinforced: bool = True,
    full_height_rebar: bool = False,
    tail_a2a3: "tuple[float, float] | None" = None,
    fc_cap: float | None = None,
    bond_law: str = "force_slip",
    bond_residual: float = 0.7,
    bond_damage: bool = False,
    steel_rupture: float | None = None,
    concrete_residual: float = 0.2,
    steel_b: float | None = None,
):
    """The calibrated RC lattice. `nonlinear=False` gives the linear-material twin for `elastic.py`.

    `gf_factor` scales the fracture energy used for crack-band regularization. It is the knob D67
    needed for WSH3 (`--gf-factor 2`): the regularized softening slope goes as f_t^2/(Gf) per strut,
    so a low-Gf/high-f_t combination produces struts that fail almost immediately after cracking and
    a lattice that cannot redistribute. This specimen's strut life is printed by `summary.py`.
    """
    check_mesh_alignment(mesh_size, length=length, height=height)
    if nonlinear and material == "aydin":
        # Aydin's OWN trilinear tension backbone (2019 Fig. 1c / 2021 Fig. 1a) instead of
        # Concrete02's bilinear Ets tail. Two differences that matter here:
        #   * the TAIL IS ~3x LONGER at the same Gf and strut length (a3 = 70.1 vs a bilinear
        #     eps_ult/eps_cr of 22.8 at L = 50), because the trilinear puts most of the fracture
        #     energy into a long low-stress tail rather than a straight ramp;
        #   * the strut has NO COMPRESSIVE STRENGTH — compression is linear elastic forever, which
        #     is exactly what the 2019 paper assumes ("Concrete in compression is assumed to be
        #     elastic"), with crushing emerging as indirect tensile splitting instead.
        # `a3_over_a2` = 5.14 is HIS ratio for this specimen (Table 1: a2 = 70, a3 = 360).
        # PATH-INDEPENDENT (ElasticMultiLinear): fine for a monotonic pushover, NOT for cyclic.
        # `tail_a2a3` overrides the Gf-solved tail with an explicit (a2, a3) — the study's `paper`
        # tail, i.e. Table 1's printed 70 / 360. Those were fitted at HIS 20 mm grid and the tail's
        # energy scales with strut length, so at mesh 50 they dissipate ~5.2x Gf: not neutral, and
        # the study prints the implied multiple for exactly that reason (PLAN §5).
        # `fc_cap` switches the compression branch from his linear-forever assumption to
        # elastic-perfectly-plastic (the study's `eppcomp`); it replaces the RSM knee, so rsm=False
        # travels with it. Tension is identical either way.
        cap = {} if fc_cap is None else {"fc_cap": fc_cap, "rsm": False}
        material_for = lambda zone, L: concrete_lattice_aydin(            # noqa: E731
            GRADES[zone], 0, L, Gf=GF * gf_factor, a3_over_a2=5.14, **cap,
            **({} if tail_a2a3 is None else {"a2": tail_a2a3[0], "a3": tail_a2a3[1]}))
        rebar_material = steel_uniaxial
    elif nonlinear:
        # `concrete_residual` floors the crushing strength at that fraction of fc. D22 set 0.2 to
        # remove the zero-tangent local mechanism that terminated a NEWTON pushover just past
        # yield; under the explicit march this study uses (D74) there is no tangent to go singular,
        # so a lower floor is available again — and it is what lets a strut actually crush.
        # `residual_ratio` alone CANNOT lower the floor here: the material takes
        # `max(grade.fcu, residual_ratio*fc)`, and this specimen's grade already carries
        # `fcu = 0.2*fc`, so the grade's own value wins and the argument does nothing. The floor has
        # to be lowered on the GRADE as well, which is what this does.
        import dataclasses
        crush_grade = {z: dataclasses.replace(g, fcu=concrete_residual * g.fc)
                       for z, g in GRADES.items()} if concrete_residual < 0.2 else GRADES
        material_for = lambda zone, L: concrete_uniaxial_regularized(     # noqa: E731
            crush_grade[zone], 0, L, Gf=GF * gf_factor, residual_ratio=concrete_residual)
        rebar_material = steel_uniaxial
    else:
        material_for = lambda zone, _L: concrete_uniaxial_elastic(GRADES[zone], 0)  # noqa: E731
        rebar_material = steel_uniaxial_elastic

    extra = {}
    if bond:
        # RE-ENABLED 2026-08-29 by user instruction, with the D72 defect addressed by CALIBRATION
        # rather than by a topology change. The star still cannot satisfy "rigid bond == perfect
        # bond" in the limit, but at a bond area small enough to leave the ELASTIC stiffness intact
        # the artefact is bounded and measurable — and elastic equivalence is the right constraint,
        # since real bond and perfect bond differ only after the bond itself starts to fail.
        if bond_law == "aydin":
            # PLAN §4, from the author (2026-09-05) and D80 item 2: a bond link IS a concrete strut
            # that drops to a residual and stays there — same area, same F_cr, elastic to ft, a
            # brittle fall to `residual * ft`, then a flat symmetric plateau. `bond_area` is left
            # OFF so the links take the concrete strut area, which is the whole point of this law.
            #
            # THIS SUPERSEDES the force-slip form below FOR THE STUDY ONLY. D75's version exists
            # because the link area had been scaled 100x down to fix stiffness, and against the
            # author's description that is a departure rather than a fix. The parent's own
            # `pushover.py --bond` keeps `force_slip`, so D75/D76's runs stay reproducible.
            #
            # KNOWN AND DELIBERATE: a full-area ring is a parallel load path beside an intact
            # lattice, so a bonded model comes out STIFFER than perfect bond (2.19x the transformed
            # section at mesh 50, D72). That is what this topology does at this bar density, not an
            # error to calibrate away — which is why every bonded run reports K/K_perfect.
            # `bond_residual` = Aydin's `a`. TWO VALUES ARE BOTH HIS, and they disagree (D96):
            #   0.7  PUBLISHED — 2019 paper, "chosen based on the preliminary simulation results of
            #        Aydin (2017) to ensure that the deformed bar residual bond strength was
            #        reflected accurately". The default here since 2026-09-08, by user instruction.
            #   0.6  THE AUTHOR, direct to this project 2026-09-05 (PLAN §4): "a bond link is a
            #        concrete strut that drops to ~60% of its capacity and stays there" — later
            #        than the paper and specific to this replication.
            # Every bonded run before 2026-09-08 used 0.6. The study axis carries both (`bond60`,
            # `bond70`), so this is a one-flag comparison, not a commitment.
            # `bond_damage` swaps the reversible ElasticMultiLinear form for the irreversible
            # Parallel(MinMax(Elastic), ElasticPP) one (D95). SAME ENVELOPE to four decimals, so a
            # monotonic run is a correctness check: any difference is a bug, not a result.
            _bond_law = bond_elastic_brittle_damaging if bond_damage else bond_elastic_brittle
            extra = dict(bond_material=lambda zone, _L: _bond_law(
                             GRADES[zone], 0, residual=bond_residual),
                         bond_horizon=horizon,
                         i_accept_the_known_bond_defect=True)
        elif bond_law == "force_slip":
            extra = dict(bond_material=_bond_material(bond_area or area, area),
                         bond_horizon=horizon,
                         bond_area=bond_area, i_accept_the_known_bond_defect=True)
        else:
            raise ValueError(f"bond_law must be 'aydin' or 'force_slip', got {bond_law!r}")

    # A rupture strain turns the bars into something that can break (MinMax over Steel02).
    if nonlinear and steel_rupture is not None:
        rebar_material = lambda grade, tag: steel_uniaxial_ruptured(      # noqa: E731
            grade, tag, eps_rupture=steel_rupture)

    model, edges = build_lattice_rc(
        wall_problem(length=length, height=height), mesh_size,
        material_for=material_for,
        zone_of=zone_of,
        rebars=_bars(mesh_size, length, height, full_height_rebar, reinforced, steel_b),
        strut_area=area,
        horizon=horizon,
        strut_element=strut_element,
        rebar_material=rebar_material,
        **extra,
    )
    return model, edges


def describe(model) -> str:
    """One line of model size, by element kind — the thing to check before a long run."""
    kinds: dict[str, int] = {}
    for e in model.elements:
        kinds[e.kind] = kinds.get(e.kind, 0) + 1
    parts = ", ".join(f"{n:,} {k}" for k, n in sorted(kinds.items()))
    return f"{len(model.nodes):,} nodes, {len(model.elements):,} elements ({parts})"


# --- closed-form elastic references ---------------------------------------------------------------
#
# The like-for-like check both other wall studies use. NOTE what it means HERE: at an aspect ratio of
# 0.75 this wall is squat enough that SHEAR carries most of the elastic flexibility, where for
# SW-NC-FF (ratio 3.00) it carried almost none. The energy balance matches the lattice's normal
# stiffness exactly and never its shear stiffness independently (D53), so on this specimen the
# closed-form comparison is testing the weakest part of the calibration — which is the point.


def gross_inertia(*, length: float = LW, thickness: float = TW) -> float:
    return thickness * length ** 3 / 12.0


def transformed_inertia(mesh_size: float = MESH, *, length: float = LW, height: float = HW,
                        thickness: float = TW) -> float:
    """Gross concrete I plus the vertical bars' (n-1)*A*x^2 about the section centroid.

    Only VERTICAL bars count — a horizontal cut severs them. Unlike VK3's section this one is not
    centred on x = 0, so the lever arm is measured from `length/2`.
    """
    inertia = gross_inertia(length=length, thickness=thickness)
    centroid = length / 2.0
    for rb in rebars(mesh_size, length=length, height=height):
        (xa, _ya), (xb, _yb) = rb.path[-2], rb.path[-1]
        if abs(xb - xa) > 1e-9:                     # horizontal bar: not cut by a horizontal plane
            continue
        arm = xb - centroid
        inertia += (rb.steel.E0 / EC - 1.0) * rb.area * arm * arm
    return inertia


def cantilever_stiffness(*, shear_span: float, inertia: float, length: float = LW,
                         thickness: float = TW, kappa: float = 1.2) -> tuple[float, float]:
    """Fixed-base elastic tip stiffness: `(K in N/mm, shear share of the flexibility)`."""
    g = EC / (2.0 * (1.0 + NU))
    flex = shear_span ** 3 / (3.0 * EC * inertia)
    shear = kappa * shear_span / (g * thickness * length)
    return 1.0 / (flex + shear), shear / (flex + shear)


def wall_continuum(*, mesh_size: float = MESH, length: float = LW, height: float = HW):
    """The SAME wall as an elastic plane-stress CONTINUUM — the like-for-like elastic reference.

    Why this and not only the closed-form cantilever (D12/D14, the builders' original purpose): at an
    aspect ratio of 0.75 a cantilever formula is a poor reference in its own right. `3EI/H^3` plus a
    kappa-corrected shear term assumes plane sections and a uniform shear distribution, and neither
    holds for a wall two-thirds as tall as it is long — shear lag, the non-uniform base restraint and
    the loaded-edge boundary all matter. Judging the lattice against it conflates the lattice's
    discretization error with the beam theory's own.

    The continuum shares the node grid, the rebar, the grades and the boundary conditions exactly, so
    `K_lattice / K_continuum` isolates one thing: whether the energy-balance area makes the lattice
    an equivalent elastic continuum, which is the only claim that calibration makes.
    """
    check_mesh_alignment(mesh_size, length=length, height=height)
    model, edges = build_continuum_rc(
        wall_problem(length=length, height=height), mesh_size,
        nd_material_for=lambda zone: concrete_nd_elastic_planestress(GRADES[zone], 0),
        zone_of=zone_of,
        rebars=rebars(mesh_size, length=length, height=height),
        rebar_material=steel_uniaxial_elastic,
    )
    return model, edges


def strut_life(gf_factor: float = 1.0, mesh_size: float = MESH, *, diagonal: bool = False) -> float:
    """`eps_ult / eps_cr` for a regularized strut — how far past cracking it carries tension.

    Crack-band regularization sets the softening slope from the strut LENGTH, so the diagonal struts
    (longer by sqrt2) are the brittle ones and are what limit a run. D67's rule of thumb from WSH3:
    below ~14 a static solve breaks up almost immediately and `--gf-factor 2` is needed.
    """
    L = mesh_size * (2.0 ** 0.5 if diagonal else 1.0)
    ets = FT * FT * L / (2.0 * GF * gf_factor)
    return (FT / EC + FT / ets) / (FT / EC)


def cracking_shear(*, length: float = LW, height: float = HW) -> tuple[float, float]:
    """(V_cr, drift at V_cr) from the model's own f_t — where a static solve starts to struggle.

    No axial load on this specimen, so the cracking moment is f_t * I / c with nothing to offset it.
    THIS WALL CRACKS EXTREMELY EARLY: the number below is an order of magnitude below WSH3's 0.025%,
    which is the run whose static solver stalled at 0.067%. Expect the static pushover here to be a
    diagnostic, not the deliverable.
    """
    i_g = gross_inertia(length=length)
    m_cr = FT * i_g / (length / 2.0)
    v_cr = m_cr / height
    k_tr, _share = cantilever_stiffness(shear_span=height,
                                        inertia=transformed_inertia(length=length, height=height),
                                        length=length)
    return v_cr, (v_cr / k_tr) / height


def strut_groups(model, *, quantity: str = "shear", cut_y: float | None = None,
                 mesh_size: float = MESH, tol: float = 1e-6) -> dict:
    """Free-body probe across a horizontal CUT, split by element orientation.

    ONLY ELEMENTS CROSSING THE CUT COUNT. Summing `eleForce` over every element in the model is
    meaningless — internal forces cancel in equal and opposite pairs, and doing that here produced a
    "diagonal share" of -16,000 kN against a 356 kN base shear. What reconciles is a free-body cut:
    for each element with one node below `cut_y` and one above, take the component at the node BELOW.

    `quantity` picks what is decomposed:

      "shear"  — horizontal force transmitted across the cut; reconciles to the base shear.
                 On a lattice this split is GEOMETRICALLY TRIVIAL and is kept only as the probe's
                 self-check: a vertical strut cut horizontally has no horizontal component at all, so
                 the diagonals carry 100% of it by construction, not by mechanism.

      "moment" — overturning moment about the section centroid, from the VERTICAL force at each
                 crossing node times its lever arm. This is the informative one: it separates the
                 FLEXURAL couple (vertical struts and vertical rebar, a plane-sections mechanism)
                 from the TRUSS action (diagonals), and watching the diagonals' share change as the
                 wall cracks is the load-path measurement D55 made on the compression cube.

    `dof` INDEXES `eleForce` DIRECTLY, i.e. 0-based: 0/1 are (Fx, Fy) at the element's first node and
    2/3 at its second. Using 1/3 for "shear" silently probes the VERTICAL components, which on an
    axially unloaded wall sum to about zero and look like a dead probe.

    `coef` carries a -1 because `eleForce` returns the force the element exerts ON its node, while
    base shear is the reaction; with it the sums reconcile to +1.0000.
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
        if (ya - cut) * (yb - cut) >= 0.0:          # does not cross the cut
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
            # BOND LINKS ARE NOT CONCRETE (D80 item 5). They sit at every orientation, so falling
            # through to the geometric branches below counted them as vertical/horizontal/diagonal
            # struts — which is what made a bonded run report "rebar 0.7% of the overturning moment"
            # against 21% without bond: a probe artefact, not a finding about the bars.
            label = "bond"
        elif abs(xb - xa) < tol:
            label = "vertical"
        elif abs(yb - ya) < tol:
            label = "horizontal"
        else:
            label = "diagonal"
        groups[label].append((e.id, dof, coef))
    return {k: v for k, v in groups.items() if v}
