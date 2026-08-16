"""Model builders + Aydin energy-balance calibration for the WSH wall series.

Same two-part recipe as the SW-NC-FF study: strut areas from Aydin's (2017) overlapping-lattice
energy balance (D47) — a homogenization needing no reference model and no FE solve — and tension
softening regularized by strut length, his crack-band principle.

The calibrated AREA is independent of E, so one value serves the wall, the stiff loading cap and
the foundation alike; each zone's modulus rides on its own material.
"""

from __future__ import annotations

from rclattice.builders import build_lattice_rc
from rclattice.calibration import EnergyBalanceResult, energy_balance_rectangle
from rclattice.materials import (
    concrete_uniaxial_elastic,
    concrete_uniaxial_regularized,
    steel_uniaxial,
    steel_uniaxial_elastic,
)

from specimen import (
    EC, FND_W, GF, GRADES, HORIZON, H_WALL_MODEL, LW, MESH, NU, TW, check_mesh_alignment, rebars,
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


def wall_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON,
                 nonlinear: bool = True, gf_factor: float = 1.0):
    """The calibrated RC lattice. `nonlinear` selects Concrete02+Steel02 or a fully linear model.

    The foundation and the loading cap stay LINEAR in both modes: both are restraint/loading
    devices whose moduli already encode geometry they do not physically have (the foundation's
    700 mm thickness, the cap's rigidity), so a nonlinear law on their struts would be meaningless.
    """
    check_mesh_alignment(mesh_size)

    def material_for(zone: str, length: float):
        grade = GRADES[zone]
        if zone in ("foundation", "cap") or not nonlinear:
            return concrete_uniaxial_elastic(grade, 0)
        gf = GF[zone] * gf_factor
        return concrete_uniaxial_regularized(grade, 0, length, Gf=gf, Gfc=250.0 * gf,
                                             residual_ratio=0.2)

    model, _edges = build_lattice_rc(
        wall_problem(), mesh_size,
        material_for=material_for,
        zone_of=zone_of,
        rebars=rebars(mesh_size),
        strut_area=area,
        horizon=horizon,
        rebar_material=steel_uniaxial if nonlinear else steel_uniaxial_elastic,
    )
    return model


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


def cantilever_stiffness(*, shear_span: float, inertia: float, kappa: float = 1.2) -> tuple[float, float]:
    """Fixed-base elastic tip stiffness: (K in N/mm, shear share of the flexibility)."""
    g = EC / (2.0 * (1.0 + NU))
    flex = shear_span ** 3 / (3.0 * EC * inertia)
    shear = kappa * shear_span / (g * TW * LW)
    return 1.0 / (flex + shear), shear / (flex + shear)
