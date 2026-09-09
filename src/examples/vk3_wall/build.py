"""Model builders + Aydin energy-balance calibration for wall-type bridge pier VK3.

Same two-part recipe as the other two wall studies: strut areas from Aydin's (2017)
overlapping-lattice energy balance (D47) — a homogenization needing no reference model and no FE
solve — and tension softening regularized by strut length, his crack-band principle.

The calibrated AREA is independent of E, so one value serves the pier, the elastic head and the
foundation alike; each zone's modulus rides on its own material. What differs per zone is EA, which
is what a truss element is actually built with, so that is what `report_calibration` prints.

A CAVEAT THAT MATTERS MORE HERE THAN IN EITHER PREVIOUS STUDY. The balance matches the lattice's
NORMAL stiffness (C11) exactly; it does not independently match its SHEAR stiffness. `nu_consistent`
(~0.18 at horizon 1.5) is merely the nu at which the normal-strain and shear routes agree — it is
NOT the lattice's Poisson ratio, which is `nu_effective` ~ 0.41, and the lattice is cubic-symmetric
rather than isotropic at any nu (D53). For SW-NC-FF and WSH3, both flexure-dominated, that was
second-order. VK3 carries 20-22% of its displacement in shear and FAILS in shear, so
`cal.isotropy_error` — the shear-stiffness error the pinned nu costs — stops being a footnote and
becomes something `elastic.py` reports against the measured shear share.
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
)

from specimen import (
    EC, FND_W, GF, GRADES, HORIZON, H_PIER, LW, MESH, NU, TW, check_mesh_alignment, pier_problem,
    rebars, zone_of,
)

ELASTIC_ZONES = ("foundation", "head")


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON) -> EnergyBalanceResult:
    """Aydin's energy balance on a pier-sized patch of the working grid.

    Run over the PIER rectangle alone: the foundation and the head are stiffness/loading devices,
    not the material being homogenized, and their re-entrant corners would bias the strut sum. The
    result is a property of (mesh_size, horizon, nu, thickness), so it transfers to the compound
    model.
    """
    return energy_balance_rectangle(LW, H_PIER, mesh_size, E=EC, nu=NU, thickness=TW,
                                    horizon=horizon)


def report_calibration(cal: EnergyBalanceResult, *, mesh_size: float = MESH,
                       horizon: float = HORIZON) -> None:
    """Print the calibration outcome: the strut AREA, and the axial rigidity EA it gives each zone.

    `A_t` alone is not the property the struts are built with — a truss element's stiffness is
    `EA/L` — so EA is the number that actually enters the model, and it differs per zone even though
    the area does not (D47/D62).
    """
    nominal = mesh_size * TW
    print(f"Aydin energy balance (horizon {horizon}, nu = {cal.nu:.2f}): "
          f"A_t = {cal.area:,.1f} mm^2 = {cal.area / nominal:.4f} * (thickness * mesh); "
          f"nu_consistent = {cal.nu_consistent:.3f} (not a Poisson ratio — nu_effective = "
          f"{cal.nu_effective:.3f}, cubic anisotropy {cal.cubic_anisotropy:.2f})")
    print(f"  shear-stiffness error at nu = {cal.nu:.2f}: {cal.isotropy_error * 100:.2f}%  "
          f"— READ THIS ONE: VK3 carries 20-22% of its displacement in shear (Fig. 5.19)")
    print(f"  concrete strut axial rigidity  EA = E * A_t   "
          f"[orthogonal L = {mesh_size:.0f} mm, diagonal L = {mesh_size * 2 ** 0.5:.1f} mm]")
    notes = {"foundation": f"  (equivalent modulus, E * {FND_W:.0f}/{TW:.0f} for the block's "
                           f"out-of-plane thickness — ASSUMED, see specimen.FND_W)",
             "head": "  (load-introduction zone, held elastic)"}
    for zone, grade in GRADES.items():
        ea = grade.E * cal.area
        print(f"    {zone:<11s} E = {grade.E:>9,.0f} MPa   EA = {ea:.4g} N   "
              f"EA/L = {ea / mesh_size:.4g} / {ea / (mesh_size * 2 ** 0.5):.4g} N/mm"
              f"{notes.get(zone, '')}")


def pier_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON):
    """The calibrated ELASTIC RC lattice: per-zone linear concrete struts + linear steel rebar.

    Every concrete strut gets the uniform calibrated `area` and its zone's `Elastic` material; rebar
    struts keep their physical bar areas and a linear steel material, so the whole model is linear
    and its tangent IS its secant.
    """
    check_mesh_alignment(mesh_size)
    model, _edges = build_lattice_rc(
        pier_problem(), mesh_size,
        material_for=lambda zone, _length: concrete_uniaxial_elastic(GRADES[zone], 0),
        zone_of=zone_of,
        rebars=rebars(mesh_size),
        strut_area=area,
        horizon=horizon,
        rebar_material=steel_uniaxial_elastic,
    )
    return model


def nonlinear_pier_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON,
                           compression: str = "crushing", gf_factor: float = 1.0,
                           strut_element: str = "Truss"):
    """The calibrated NONLINEAR RC lattice: Concrete02 struts + Steel02 rebar.

    Strut areas still come from Aydin's ELASTIC energy balance — that calibration fixes the initial
    tangent and is unchanged by what happens after cracking.

    TENSION follows Aydin's principle exactly: a softening branch whose strain axis is scaled by the
    strut LENGTH, so the stress-crack-opening curve, and hence the fracture energy per unit crack
    area, is the same for every strut (`concrete_uniaxial_regularized`, D20).

    COMPRESSION departs from Aydin (his Sec. 2.2 assumes it stays elastic), and the choice is left
    open so the cost of that assumption is measurable rather than argued about:
      * `"crushing"` (default) — the real Concrete02 backbone, regularized the same crack-band way.
        NOT optional here in practice: the chapter's own reading of the failure is that the base
        COMPRESSION ZONE was destroyed and took the diagonal strut's support with it (Sec. 5.4), so
        a model whose compression never fails cannot express the mechanism under study.
      * `"elastic"` — Aydin's literal assumption, imposed by pushing the crushing strain out of
        range. Kept as the control that prices the choice.

    CONFINEMENT GETS NO SEPARATE GRADE, and here that is barely even a choice: VK3's hoops are
    phi6@200 with 90-degree hooks, rho_sw = 0.08%. The chapter says "no special measures were taken
    concerning the detailing of the hoops for ductile behavior" and that they could open once the
    cover spalled. There is no confined core to model. The hoops are present as in-plane rebar
    struts, which is the whole of their in-plane contribution.

    The foundation and the head stay LINEAR in both modes: both are restraint/loading devices whose
    moduli or detailing already encode something they do not physically have, so a nonlinear law on
    their struts would be meaningless.

    `strut_element` is "Truss" (small-displacement, the default everywhere in this repo) or
    "corotTruss". The distinction is NOT cosmetic for VK3: it is a squat pier carrying 1300 kN of
    axial load to ~1% drift, so the P-delta term and the rotation of every strut's axis are
    first-order effects that a small-displacement element simply omits. D22 introduced corotTruss
    for exactly this regime, and D60 found the Aydin cube's failure mode is unreachable without it.

    `gf_factor` scales the tensile fracture energy of the pier zone. It exists because the plain-
    concrete MC90 value understates what concrete EMBEDDED in reinforcement does — between cracks
    the bars carry load across the crack, the tension-stiffening effect smeared-crack RC models
    represent the same way. It also delays the numerical mechanism that ends these runs, so it is
    BOTH a physical correction and a stabiliser: keep it modest and report it. VK3's strut life
    eps_ult/eps_cr is 13.9, close to WSH3's 13.4, which needed `--gf-factor 2` to run (D67).
    """
    if compression not in ("crushing", "elastic"):
        raise ValueError(f"compression must be 'crushing' or 'elastic', got {compression!r}")
    check_mesh_alignment(mesh_size)

    def material_for(zone: str, length: float):
        grade = GRADES[zone]
        if zone in ELASTIC_ZONES:
            return concrete_uniaxial_elastic(grade, 0)
        if compression == "elastic":
            grade = replace(grade, fcu=grade.fc, epsU=1.0)
        gf = GF[zone] * gf_factor
        return concrete_uniaxial_regularized(grade, 0, length, Gf=gf, Gfc=250.0 * gf,
                                             residual_ratio=0.2)

    model, _edges = build_lattice_rc(
        pier_problem(), mesh_size,
        material_for=material_for,
        zone_of=zone_of,
        rebars=rebars(mesh_size),
        strut_area=area,
        horizon=horizon,
        strut_element=strut_element,
        rebar_material=steel_uniaxial,
    )
    return model


def strut_life(gf_factor: float = 1.0, length: float = MESH) -> float:
    """`eps_ult / eps_cr` for a pier strut — the brittleness number D67 traced a lost run to.

    With crack-band regularization the softening modulus is `Ets = f_t^2 * L / (2*Gf)`, so a strut
    falls from cracking to zero stress in `1 + E/Ets` times its cracking strain. Bigger is gentler:
    SW-NC-FF ran clean at 19.8, WSH3 diverged at 10.1 and needed tension stiffening to reach 13.4.
    """
    from specimen import FT, GF_C
    ets = FT ** 2 * length / (2.0 * GF_C * gf_factor)
    return 1.0 + EC / ets


def gross_inertia() -> float:
    return TW * LW ** 3 / 12.0


def transformed_inertia(mesh_size: float = MESH) -> float:
    """Gross concrete I plus the vertical bars' (n-1)*A*x^2, uncracked transformed section.

    Only vertical bars count — a horizontal cut severs them. The LAST path segment is tested, not
    first-to-last, so a bar whose path bends would still be classified correctly.
    """
    inertia = gross_inertia()
    for rb in rebars(mesh_size):
        (xa, _ya), (xb, _yb) = rb.path[-2], rb.path[-1]
        if abs(xb - xa) > 1e-9:
            continue
        inertia += (rb.steel.E0 / EC - 1.0) * rb.area * xb * xb
    return inertia


def cantilever_stiffness(*, shear_span: float, inertia: float,
                         kappa: float = 1.2) -> tuple[float, float]:
    """Fixed-base elastic tip stiffness: `(K in N/mm, shear share of the flexibility)`.

    The shear share is worth watching here: at Lv/lw = 2.20 this pier is squat enough that shear is
    ~14% of the ELASTIC flexibility before any cracking, against ~5% for a slender wall.
    """
    g = EC / (2.0 * (1.0 + NU))
    flex = shear_span ** 3 / (3.0 * EC * inertia)
    shear = kappa * shear_span / (g * TW * LW)
    return 1.0 / (flex + shear), shear / (flex + shear)
