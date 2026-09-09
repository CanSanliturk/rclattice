"""Build Aydin's model of the Aldemir wall, using HIS calibration route and HIS material laws."""

from __future__ import annotations

from rclattice.builders import build_lattice_rc
from rclattice.calibration import aydin_closed_form_C, energy_balance_rectangle
from rclattice.materials import (
    bond_elastic_brittle,
    concrete_lattice_aydin,
    concrete_lattice_aydin_cyclic,
    concrete_uniaxial_regularized,
    concrete_uniaxial_elastic,
    steel_uniaxial,
    steel_uniaxial_elastic,
)

from specimen import (
    A1, A2, A3, B1, B2, BOND_RESIDUAL, CONCRETE, EC, ES, GF, GRADES, HORIZON, HW, LW, MESH, NU, TW,
    rebars, wall_problem, zone_of,
)


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON):
    """HIS energy balance, not ours: the equibiaxial equal-stress field of his Appendix, at nu = 1/3.

    Our own studies use the thesis route (`field="uniaxial"`) at nu = 0.20, which returns ~18% more
    EA at matched nu. Replicating his result means replicating the calibration that produced it, so
    this is deliberately not the repo default.
    """
    return energy_balance_rectangle(LW, HW, mesh_size, E=EC, nu=NU, thickness=TW,
                                    horizon=horizon, field="equibiaxial")


def report_calibration(cal, *, mesh_size: float = MESH, horizon: float = HORIZON) -> None:
    closed = aydin_closed_form_C(horizon, nu=NU)
    print(f"Aydin's Appendix balance (equibiaxial, nu = 1/3, horizon {horizon}):")
    print(f"  A_t = {cal.area:,.1f} mm^2   C = {cal.EA / (EC * mesh_size * TW):.4f}   "
          f"his closed form C = {closed:.4f} (paper prints 0.621)")
    print(f"  EA = {cal.EA:.4g} N   ->  EA/L = {cal.EA / mesh_size:.4g} N/mm orthogonal")


def _bond_material(bond_area: float, strut_area: float):
    """His Fig. 1(c) bond law, written as a FORCE-SLIP law rather than a stress-strain one (D75).

    A truss couples stiffness and strength through one number — `k = EA/L` and `F = ft*A` — so the
    link AREA cannot set both. Sizing it for stiffness (0.01*A_t, the ratio that leaves the elastic
    stiffness intact) divides the strength by 100 as a side effect, the ring cracks at once, and the
    reinforcement disconnects EARLIER than with no bond at all. So the two are set independently:

      * strength — `peak_force` is the concrete strut's own cracking force `ft * A_t`, because his
        Fig. 1(c) draws Bond and Concrete rising together to the SAME `F_cr`;
      * stiffness — `slip_at_peak` is chosen so the link's modulus comes out at the concrete's `E`.

    NOTE THE MESH DEPENDENCE, which matters more here than in the parent study. The failure slip is
    `ft/(ratio*E) * length`, i.e. proportional to the LINK LENGTH, so his 20 mm mesh gives 2.5x less
    slip than the parent study's 50 mm at the same area ratio: 0.149 mm against 0.372. Both sit
    inside the Model Code `s1` range for a deformed bar (0.1-1 mm), but this one sits at its floor.
    Lowering the ratio raises the slip and drops the stiffness — and 0.003 went singular at mesh 50
    (D75), so the usable window is narrow and is NOT known for this mesh. `report_bond` prints it.
    """
    def make(zone, length):
        grade = GRADES[zone]
        f_cr = grade.ft * strut_area
        slip = (f_cr / bond_area) / grade.E * length      # -> E_bond == grade.E
        return bond_elastic_brittle(grade, 0, residual=BOND_RESIDUAL, peak_force=f_cr,
                                    area=bond_area, length=length, slip_at_peak=slip)
    return make


def report_bond(area: float, bond_area: float, *, mesh_size: float = MESH) -> None:
    """Print the inferred bond numbers — none of which he publishes — so a run can be judged."""
    grade = GRADES["wall"]
    f_cr = grade.ft * area
    ft_bond = f_cr / bond_area
    slip = ft_bond / grade.E * mesh_size
    print(f"\nBOND (his Fig. 1c; only a = {BOND_RESIDUAL} is published, the rest is inferred):")
    print(f"  A_bond = {bond_area:,.1f} mm^2 = {bond_area / area:.4f} * A_t   "
          f"F_cr = {f_cr / 1e3:.2f} kN (the strut's own)")
    print(f"  ft_bond = {ft_bond:,.0f} MPa   E_bond = {grade.E:,.0f} MPa (= concrete's, by "
          f"construction)")
    print(f"  failure slip = {slip:.3f} mm over a {mesh_size:g} mm link "
          f"({slip / mesh_size * 1e2:.2f}% strain) — Model Code s1 for a deformed bar is 0.1-1 mm")
    if not 0.1 <= slip <= 1.0:
        print(f"  ** OUTSIDE the Model Code range — a link failing at {slip:.3f} mm is not bond **")
    print("  the AREA RATIO was calibrated against perfect bond at mesh 50 in the parent study, "
          "NOT here.")
    print("  Check it before trusting a pushover:  run.py --elastic --bond  vs  run.py --elastic")


def element_modulus(e) -> float:
    """E for `critical_time_step`, which sees only material TAGS and so has to be told (D74).

    Bond links get the CONCRETE modulus because `_bond_material` derives their slip so that
    `E_bond == grade.E` exactly — if that derivation changes, this must change with it.
    """
    return ES if e.kind in ("longitudinal", "stirrup") else EC


def explicit_steps_per_period(model, *, safety: float = 0.8) -> tuple[int, float, float]:
    """Size an EXPLICIT run: `(steps_per_period, dt_crit, T1)`.

    The stable step comes from the stiffest/lightest ELEMENT (`2/w_max`), never from T1 — sizing
    off T1 is the implicit habit and is 15-20x too large here (D74). T1 is fetched only because
    `run_pushover_dynamic` takes its step as a fraction of a period.
    """
    from rclattice.builders import critical_time_step
    from rclattice.opensees import run_modal

    dt_crit, _w, _worst = critical_time_step(model, element_modulus)
    t1 = run_modal(model, 1)["periods"][0]
    return int(t1 / (safety * dt_crit)) + 1, dt_crit, t1


def wall_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON,
                 nonlinear: bool = True, strut_element: str = "corotTruss",
                 bond: bool = False, bond_area: float | None = None, cap_rows: int = 0,
                 full_height_rebar: bool = False, cyclic: bool = False,
                 cyclic_law: str = "concrete02", gf_factor: float = 1.0):
    """His lattice: tension-only concrete on his published tail, elastic-perfectly-plastic steel.

    Three departures from this repo's own wall models, all his:
      * concrete has NO COMPRESSIVE STRENGTH (his Sec. "Reinforced Concrete Lattice Modeling":
        "Concrete in compression is assumed to be elastic"), with crushing emerging as indirect
        tensile splitting followed by lattice instability;
      * the tension tail is his PUBLISHED a2 = 70 / a3 = 360, pinned rather than solved from Gf;
      * steel is elastic-perfectly plastic (b = 0), his Fig. 1(d).

    `bond=False` uses SHARED NODES, the one departure `specimen.not_replicated` calls the largest.
    `bond=True` switches to his actual scheme: separate steel nodes tied to the concrete by a ring
    of elastic-brittle links (D72/D75). It roughly doubles the model — measure the cost with
    `preflight.py --bond` before committing to a pushover.
    """
    # `a1 = horizon` is not a convenience: Table 1's footnote d defines it that way (a1 = 1.5 at
    # delta = 1.5d, 3.01 at 3.01d), so the softening law travels with the connectivity.
    # LOAD-INTRODUCTION CAP (`cap_rows`). The drive is a prescribed displacement on ONE row of top
    # nodes, so the whole base shear leaves the wall through the struts immediately below that row.
    # The h=3.01 run tore exactly there — 93% of post-collapse strain in the top row at 25% strain,
    # base and middle fully unloaded. Making the top `cap_rows` rows ELASTIC (they cannot crack at
    # any strain) is the stiff loading beam the specimen actually had, and turns "is the collapse a
    # load-introduction artefact?" into a controlled experiment: same wall, same everything else,
    # only the top few rows made unfailable.
    #
    # NOTE this is the RIGHT form of that test and `equalDOF` across the driven row is NOT: the
    # drive already imposes an identical `ops.sp` on every top node, so they translate together
    # regardless, and tying them changes nothing (D53 measured exactly this on the compression
    # cube). What tears is the interface BELOW the driven row, which only a cap can address.
    zone_fn = zone_of
    if cap_rows:
        cap_y = HW - cap_rows * mesh_size - 1e-6
        zone_fn = lambda _x, y: "cap" if y >= cap_y else "wall"          # noqa: E731

    if nonlinear:
        def material_for(zone, L):
            if zone == "cap":
                return concrete_uniaxial_elastic(CONCRETE, 0)     # unfailable: the loading beam
            if cyclic:
                # `concrete_lattice_aydin` is ElasticMultiLinear — PATH-INDEPENDENT, so a cracked
                # strut recovers full stiffness on reload and carries no damage memory. Fine for a
                # monotonic push, MEANINGLESS cyclically (D60). Two substitutes, and the choice
                # decides whether the wall can PINCH:
                #
                #   "aydin"      — his trilinear tension shape on a HystereticSM. Keeps the softening
                #                  shape but a cycled strut still holds ~4.4 MPa of compression AT
                #                  ZERO STRAIN, so cracked struts never release and the loops stay
                #                  fat. MEASURED at the strut level; `pinch_x`/`pinch_y` do NOT fix
                #                  it (0.3/0.1 -> 0.10/0.03 moved 4.38 -> 4.43 MPa, i.e. nothing).
                #   "concrete02" — returns to 0.29 MPa at zero strain, 15x less, which IS the
                #                  pinching mechanism: a cracked strut carries nothing until its
                #                  faces re-contact. Gives up his trilinear tail, which on VK3 was
                #                  worth 3x in strain capacity — so expect earlier failure and
                #                  reach for `gf_factor` if it collapses early.
                if cyclic_law == "aydin":
                    return concrete_lattice_aydin_cyclic(GRADES[zone], 0, L, Gf=GF * gf_factor,
                                                         a1=horizon, b1=B1, b2=B2,
                                                         a3_over_a2=A3 / A2)
                return concrete_uniaxial_regularized(GRADES[zone], 0, L, Gf=GF * gf_factor)
            return concrete_lattice_aydin(GRADES[zone], 0, L, Gf=GF, a1=horizon,
                                          b1=B1, b2=B2, a2=A2, a3=A3, rsm=False)
        rebar_material = steel_uniaxial
    else:
        material_for = lambda _zone, _L: concrete_uniaxial_elastic(CONCRETE, 0)   # noqa: E731
        rebar_material = steel_uniaxial_elastic

    extra = {}
    if bond:
        extra = dict(bond_material=_bond_material(bond_area or area, area),
                     bond_horizon=horizon,
                     bond_area=bond_area,
                     i_accept_the_known_bond_defect=True)

    model, edges = build_lattice_rc(
        wall_problem(), mesh_size,
        material_for=material_for, zone_of=zone_fn,
        rebars=rebars(mesh_size, full_height=full_height_rebar),
        strut_area=area, horizon=horizon, strut_element=strut_element,
        rebar_material=rebar_material,
        **extra,
    )
    return model, edges


def describe(model) -> str:
    kinds: dict[str, int] = {}
    for e in model.elements:
        kinds[e.kind] = kinds.get(e.kind, 0) + 1
    return (f"{len(model.nodes):,} nodes, {len(model.elements):,} elements ("
            + ", ".join(f"{n:,} {k}" for k, n in sorted(kinds.items())) + ")")
