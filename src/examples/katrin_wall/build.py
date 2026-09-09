"""Model builders + Aydin energy-balance calibration for the WSH wall series.

Same two-part recipe as the SW-NC-FF study: strut areas from Aydin's (2017) overlapping-lattice
energy balance (D47) — a homogenization needing no reference model and no FE solve — and tension
softening regularized by strut length, his crack-band principle.

The calibrated AREA is independent of E, so one value serves the wall, the stiff loading cap and
the foundation alike; each zone's modulus rides on its own material. What differs per zone is EA,
which is what a truss element is actually built with, so that is what `report_calibration` prints.

`nu_consistent` (~0.18 for the horizon = 1.5 grid) is NOT the lattice's Poisson ratio — it is the nu
at which the normal-strain and shear calibration routes agree. The lattice's actual lateral response
is `nu_effective` ~ 0.41 and it is cubic-symmetric rather than isotropic at any nu (D53).
"""

from __future__ import annotations

from dataclasses import replace

from rclattice.builders import build_lattice_rc
from rclattice.calibration import EnergyBalanceResult, energy_balance_rectangle
from rclattice.materials import (
    concrete_uniaxial_elastic,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
    steel_uniaxial_ruptured,
)

from specimen import (
    AGT, EC, FND_W, GF, GRADES, HORIZON, H_WALL_MODEL, LW, MESH, NU, TW, check_mesh_alignment, rebars,
    wall_problem, zone_of,
)


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON) -> EnergyBalanceResult:
    """Aydin's energy balance on a wall-sized patch of the working grid.

    Run over the WALL rectangle alone: the foundation and cap are stiffness devices, not the
    material being homogenized, and their re-entrant corners would bias the strut sum. The result
    is a property of (mesh_size, horizon, nu, thickness), so it transfers to the compound model.
    """
    return energy_balance_rectangle(LW, H_WALL_MODEL, mesh_size, E=EC, nu=NU, thickness=TW,
                                    horizon=horizon)


def report_calibration(cal: EnergyBalanceResult, *, mesh_size: float = MESH,
                       horizon: float = HORIZON) -> None:
    """Print the calibration outcome: the strut AREA, and the axial rigidity EA it gives each zone.

    `A_t` alone is not the property the struts are built with — a truss element's stiffness is
    `EA/L` — so EA is the number that actually enters the model, and it differs per zone even
    though the area does not. Aydin's balance returns EA and divides by E only to make the result
    transferable (D47/D62), so reporting EA is reporting the calibration's own output.

    Both non-wall zones carry a deliberately non-physical modulus: the foundation's stands in for a
    block 700 mm thick in a 150 mm plane model, the cap's is 20x the wall's to make it rigid.
    """
    nominal = mesh_size * TW
    print(f"Aydin energy balance (horizon {horizon}, nu = {cal.nu:.2f}): "
          f"A_t = {cal.area:,.1f} mm^2 = {cal.area / nominal:.4f} * (thickness * mesh); "
          f"nu_consistent = {cal.nu_consistent:.3f} (not a Poisson ratio — nu_effective = "
          f"{cal.nu_effective:.3f}, cubic anisotropy {cal.cubic_anisotropy:.2f})")
    print(f"  concrete strut axial rigidity  EA = E * A_t   "
          f"[orthogonal L = {mesh_size:.0f} mm, diagonal L = {mesh_size * 2 ** 0.5:.1f} mm]")
    notes = {"foundation": f"  (equivalent modulus, E * {FND_W:.0f}/{TW:.0f})",
             "cap": "  (stiffened loading cap, not a material)"}
    for zone, grade in GRADES.items():
        ea = grade.E * cal.area
        print(f"    {zone:<11s} E = {grade.E:>9,.0f} MPa   EA = {ea:.4g} N   "
              f"EA/L = {ea / mesh_size:.4g} / {ea / (mesh_size * 2 ** 0.5):.4g} N/mm"
              f"{notes.get(zone, '')}")


def wall_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON):
    """The calibrated ELASTIC RC lattice: per-zone linear concrete struts + linear steel rebar.

    Every concrete strut gets the uniform calibrated `area` and its zone's `Elastic` material;
    rebar struts keep their physical bar areas and a linear steel material, so the whole model is
    linear and its tangent IS its secant.
    """
    check_mesh_alignment(mesh_size)
    model, _edges = build_lattice_rc(
        wall_problem(), mesh_size,
        material_for=lambda zone, _length: concrete_uniaxial_elastic(GRADES[zone], 0),
        zone_of=zone_of,
        rebars=rebars(mesh_size),
        strut_area=area,
        horizon=horizon,
        rebar_material=steel_uniaxial_elastic,
    )
    return model


def nonlinear_wall_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON,
                           compression: str = "crushing", gf_factor: float = 1.0,
                           steel_rupture: float | str | None = None,
                           concrete_residual: float = 0.2):
    """The calibrated NONLINEAR RC lattice: Concrete02 struts + Steel02 rebar.

    Strut areas still come from Aydin's ELASTIC energy balance — that calibration fixes the initial
    tangent and is unchanged by what happens after cracking.

    TENSION follows Aydin's principle exactly: a softening branch whose strain axis is scaled by the
    strut LENGTH, so the stress-crack-opening curve, and hence the fracture energy per unit crack
    area, is the same for every strut (`concrete_uniaxial_regularized`, D20).

    COMPRESSION departs from Aydin (his Sec. 2.2 assumes it stays elastic), and the choice is left
    open so the cost of that assumption is measurable rather than argued about:
      * `"crushing"` (default) — the real Concrete02 backbone, regularized the same crack-band way.
        Needed here: WSH3's boundary bars buckled and its cover spalled, so the compression zone
        matters even though the wall never crushed outright.
      * `"elastic"` — Aydin's literal assumption, imposed by pushing the crushing strain out of range.

    CONFINEMENT IS NOT GIVEN A SEPARATE GRADE, unlike the RC column study (D23), and that is a
    decision rather than an omission. WSH3's boundary hoops are already in the model as rebar struts
    spanning between the boundary bars (`specimen.boundary_ties`), and in a plane lattice those
    struts DO restrain the in-plane transverse expansion of the compressed boundary — the mechanism
    a confined concrete law represents. Layering a Mander-type grade on top would count the same
    restraint twice. What the model still cannot see is the through-thickness leg of each hoop and
    the phi4.2 crossties, both out of plane; that part of the confinement is simply absent.

    The foundation and the loading cap stay LINEAR in both modes: both are restraint/loading devices
    whose moduli already encode geometry they do not physically have, so a nonlinear law on their
    struts would be meaningless.

    `steel_rupture` and `concrete_residual` are the two FAILURE SWITCHES (D91), and this specimen
    is where they can be validated rather than merely applied (D93). Without them the bars cannot
    break and a floored crushing strength cannot crush, so the model has no drift capacity by
    construction — which is exactly what CLAUDE.md records as missing for SW-NC-FF.

      * `steel_rupture="measured"` wraps each bar in `MinMax` at ITS OWN published A_gt
        (`specimen.AGT`, Table 2): 7.69% for the phi12 boundary bars, 7.34% phi8 web, 6.45% phi6,
        3.06% phi4.2. On the Aldemir wall this number is an assumption that sets the answer; here
        it is measured, and Dazio, Beyer & Bachmann report the corner bar rupturing at 1.79% drift,
        so the model can be checked against a published failure rather than fitted to it. A float
        applies one strain to every bar instead; `None` (default) leaves the bars unbreakable.
      * `concrete_residual` floors crushing at that fraction of fc. It must be lowered ON THE GRADE
        as well as passed on: the material takes `max(grade.fcu, residual_ratio*fc)` and these
        grades already carry `fcu = 0.2*fc`, so passing the argument alone changes NOTHING (D91).

    EXPECT THE MODEL TO RUPTURE LATE. The test's corner bar broke after buckling from 1.70% drift,
    and a bar that has buckled and straightened fractures at a far lower nominal tensile strain
    than A_gt. Nothing here buckles, so a prediction near or past 1.79% is the honest expectation
    and an EARLY rupture would be the surprising result.

    `gf_factor` scales the tensile fracture energy of every nonlinear zone. It exists because the
    plain-concrete MC90 value understates what concrete EMBEDDED in reinforcement does — between
    cracks the bars carry load across the crack, the tension-stiffening effect smeared-crack RC
    models represent the same way. It also delays the numerical mechanism that ends these runs, so
    it is BOTH a physical correction and a stabiliser: keep it modest and report it.
    """
    if compression not in ("crushing", "elastic"):
        raise ValueError(f"compression must be 'crushing' or 'elastic', got {compression!r}")
    if compression == "elastic" and concrete_residual < 0.2:
        raise ValueError("concrete_residual is meaningless with compression='elastic' — that mode "
                         "holds the plateau at fc so there is no crushing to floor")
    check_mesh_alignment(mesh_size)

    def material_for(zone: str, length: float):
        grade = GRADES[zone]
        if zone in ("foundation", "cap"):
            return concrete_uniaxial_elastic(grade, 0)
        if compression == "elastic":
            # Aydin: compression never crushes. Push epsU far out and hold the plateau at fc.
            grade = replace(grade, fcu=grade.fc, epsU=1.0)
        elif concrete_residual < 0.2:
            # The floor has to come down on the GRADE; `residual_ratio` alone is defeated by
            # `max(grade.fcu, ratio*fc)` (D91).
            grade = replace(grade, fcu=concrete_residual * grade.fc)
        gf = GF[zone] * gf_factor
        return concrete_uniaxial_regularized(grade, 0, length, Gf=gf, Gfc=250.0 * gf,
                                             residual_ratio=concrete_residual)

    rebar_material = steel_uniaxial
    if steel_rupture is not None:
        # Per-grade when "measured", so each bar breaks at the ductility ITS OWN coupon showed.
        def rebar_material(grade, tag):                            # noqa: F811
            eps = AGT[grade.name] if steel_rupture == "measured" else float(steel_rupture)
            return steel_uniaxial_ruptured(grade, tag, eps_rupture=eps)

    model, _edges = build_lattice_rc(
        wall_problem(), mesh_size,
        material_for=material_for,
        zone_of=zone_of,
        rebars=rebars(mesh_size),
        strut_area=area,
        horizon=horizon,
        rebar_material=rebar_material,
    )
    return model


def gross_inertia() -> float:
    return TW * LW ** 3 / 12.0


def transformed_inertia(mesh_size: float = MESH) -> float:
    """Gross concrete I plus the vertical bars' (n-1)*A*x^2, uncracked transformed section.

    Only vertical bars count — a horizontal cut severs them. The LAST path segment is tested, not
    first-to-last, so a bar whose path bends would still be classified correctly.

    The effect is far smaller here than for SW-NC-FF: WSH3's concrete is strong (Ec = 35.2 GPa), so
    n = Es/Ec is only 5.7 against that wall's ~12-14, and the bars add ~5% to I rather than ~20%.
    """
    inertia = gross_inertia()
    for rb in rebars(mesh_size):
        (xa, _ya), (xb, _yb) = rb.path[-2], rb.path[-1]
        if abs(xb - xa) > 1e-9:
            continue
        inertia += (rb.steel.E0 / EC - 1.0) * rb.area * xb * xb
    return inertia


def cantilever_stiffness(*, shear_span: float, inertia: float, kappa: float = 1.2) -> tuple[float, float]:
    """Fixed-base elastic tip stiffness: (K in N/mm, shear share of the flexibility)."""
    g = EC / (2.0 * (1.0 + NU))
    flex = shear_span ** 3 / (3.0 * EC * inertia)
    shear = kappa * shear_span / (g * TW * LW)
    return 1.0 / (flex + shear), shear / (flex + shear)
