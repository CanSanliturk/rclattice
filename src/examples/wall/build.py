"""Model builders + Aydin energy-balance calibration for the RC shear-wall study.

Calibration follows Aydin (2017), OLM Sec. 2.2 (Eqs 2.1-2.3): equate the elastic energy a
continuum stores under an affine strain field to the energy the same-geometry lattice stores with
EA = 1, and solve for the one uniform EA. It needs NO reference model and NO FE solve — the
affine field makes every strut elongation exact, so the balance is a closed-form sum over the
strut list (`rclattice.calibration.energy_balance_area`, D47).

Two consequences shape how it is used here:

  * A_t = EA_t / E is INDEPENDENT of E. One calibrated area therefore serves all three concrete
    zones of this specimen (test unit / upper head / pedestal) — each zone's own modulus rides on
    its material, not on its strut area.
  * A horizon-strut lattice is an isotropic solid at only ONE Poisson ratio. For the horizon=1.5
    (8-neighbour) grid that value is ~0.18, close enough to concrete's 0.20 that the residual
    error is small; `calibrate()` reports it so the assumption stays visible.
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
    GF, GRADES, HORIZON, HW, LW, MESH, NU, TEST_UNIT, TW, check_mesh_alignment, rebars,
    wall_problem, zone_of,
)


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON) -> EnergyBalanceResult:
    """Aydin's energy balance on a wall-sized patch of the working grid.

    The balance is run over the WALL rectangle (the member being homogenised) rather than the
    full wall+pedestal domain: the pedestal is a stiffness device, not part of the material being
    calibrated, and its re-entrant corner would bias the strut sum. Because the result is a
    property of (mesh_size, horizon, nu, thickness) alone, it transfers to the compound model.
    """
    return energy_balance_rectangle(LW, HW, mesh_size, E=TEST_UNIT.E, nu=NU, thickness=TW,
                                    horizon=horizon)


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
        rebars=rebars(),
        strut_area=area,
        horizon=horizon,
        rebar_material=steel_uniaxial_elastic,
    )
    return model


def nonlinear_wall_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON,
                           compression: str = "crushing", gf_factor: float = 1.0):
    """The calibrated NONLINEAR RC lattice for the cyclic study: Concrete02 struts + Steel02 rebar.

    Strut areas still come from Aydin's ELASTIC energy balance — that calibration fixes the initial
    tangent and is unchanged by what happens after cracking.

    TENSION follows Aydin's principle exactly: a softening branch whose strain axis is scaled by the
    strut LENGTH so the stress-crack-opening curve, and hence the fracture energy per unit crack
    area, is the same for every strut. Aydin does this by multiplying his multilinear strain
    multipliers a1,a2,a3 by d/L (his Sec. 2.2); `concrete_uniaxial_regularized` (D20) does it by
    sizing the softening slope as `Ets = ft^2*L/(2*Gf)`. Same crack-band idea; his backbone is
    trilinear where Concrete02's is bilinear.

    COMPRESSION is where this DEPARTS from Aydin, deliberately:
      * `"crushing"` (default) — the real Concrete02 compression backbone, regularized in the same
        crack-band way. Needed because this specimen fails by concrete core crushing, which is the
        whole point of the cyclic test.
      * `"elastic"` — Aydin's literal assumption (his Sec. 2.2: compression is elastic and
        uncalibrated). Kept selectable so the cost of that assumption is measurable rather than
        argued about; it is imposed by pushing the crushing strain far out of range.

    The pedestal stays LINEAR in both modes: it is a stiffness device standing in for a block the
    test measured as effectively rigid, and its equivalent modulus already encodes a thickness it
    does not physically have — a nonlinear law on those struts would be meaningless.

    `gf_factor` scales the tensile fracture energy of every concrete zone. It exists because the
    plain-concrete MC90 values in `specimen.GF` understate what concrete EMBEDDED in reinforcement
    does: between cracks the bars carry load across the crack, so the effective softening of the
    concrete is far gentler than a plain coupon's — the tension-stiffening effect that smeared-crack
    RC models represent the same way. Raising it also delays the numerical mechanism that ends these
    runs (cracked tension-side nodes losing their lateral bracing), so it is BOTH a physical
    correction and a stabiliser; keep it modest and report it, since it is the one knob here that
    can quietly buy convergence at the price of strength.
    """
    if compression not in ("crushing", "elastic"):
        raise ValueError(f"compression must be 'crushing' or 'elastic', got {compression!r}")
    check_mesh_alignment(mesh_size)

    def material_for(zone: str, length: float):
        if zone == "pedestal":
            return concrete_uniaxial_elastic(GRADES[zone], 0)
        grade = GRADES[zone]
        if compression == "elastic":
            # Aydin: compression never crushes. Push epsU far out and hold the plateau at fc.
            grade = replace(grade, fcu=grade.fc, epsU=1.0)
        gf = GF[zone] * gf_factor
        return concrete_uniaxial_regularized(grade, 0, length, Gf=gf, Gfc=250.0 * gf,
                                             residual_ratio=0.2)

    model, _edges = build_lattice_rc(
        wall_problem(), mesh_size,
        material_for=material_for,
        zone_of=zone_of,
        rebars=rebars(),
        strut_area=area,
        horizon=horizon,
        rebar_material=steel_uniaxial,
    )
    return model


def gross_inertia() -> float:
    """Second moment of the plain concrete section about its in-plane bending axis."""
    return TW * LW ** 3 / 12.0


def transformed_inertia() -> float:
    """Gross concrete I plus the reinforcement's (n-1)*A*x^2, uncracked transformed section.

    Only the VERTICAL bars count: a horizontal cut through the wall severs them, so they carry
    flexural stress. The horizontal web bars run parallel to that cut and contribute nothing to I
    (they do stiffen the wall in shear, which this hand check does not attempt to transform).

    This matters far more than usual here — the specimen's low-strength test-unit concrete gives a
    modular ratio n = Es/Ec of about 12-14, so the boundary groups alone add roughly a fifth of the
    gross stiffness. Comparing a REINFORCED lattice against a plain-concrete section would
    understate the target by that much.
    """
    inertia = gross_inertia()
    for rb in rebars():
        # Test the LAST segment, not first-to-last: a hooked bar starts at its hook toe in the
        # pedestal, so comparing path endpoints would misread every longitudinal bar as horizontal.
        (xa, _ya), (xb, _yb) = rb.path[-2], rb.path[-1]
        if abs(xb - xa) > 1e-9:          # horizontal bar — parallel to the cut, no contribution
            continue
        n = rb.steel.E0 / TEST_UNIT.E
        inertia += (n - 1.0) * rb.area * xb * xb
    return inertia


def cantilever_stiffness(*, shear_span: float, inertia: float, kappa: float = 1.2) -> tuple[float, float]:
    """Fixed-base elastic cantilever tip stiffness: (K in N/mm, shear share of the flexibility).

    The independent hand check for the calibration — Euler-Bernoulli flexure on `inertia` plus a
    Timoshenko shear term on the gross web area.
    """
    E = TEST_UNIT.E
    G = E / (2.0 * (1.0 + NU))
    flex = shear_span ** 3 / (3.0 * E * inertia)
    shear = kappa * shear_span / (G * TW * LW)
    return 1.0 / (flex + shear), shear / (flex + shear)
