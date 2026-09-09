"""Model builder + Aydin energy-balance calibration for the compression cube (D53-D56).

Calibration follows Aydin (2017), OLM Sec. 2.2 (Eqs 2.1-2.3): equate the elastic energy a continuum
stores under an affine strain field to the energy the same-geometry lattice stores with EA = 1, and
solve for the one uniform EA. No FE solve, no optimiser, no reference model — a closed-form sum over
the strut list (`rclattice.calibration.energy_balance_rectangle`, D47).

A uniform-EA horizon lattice is NOT an isotropic solid at any Poisson ratio — it is cubic-symmetric,
with C66 about 1.38x the isotropic (C11-C12)/2 at horizon 1.5. Its actual Poisson ratio is ~0.41,
and the ~0.18 that `EnergyBalanceResult.nu_consistent` reports is a different quantity: the nu at
which the normal-strain and shear CALIBRATION ROUTES return the same area. Both are printed by
`pushover.py` so the distinction stays visible — see `expected_modulus` for its consequence.
"""

from __future__ import annotations

from dataclasses import replace

from rclattice.builders import build_lattice_rc
from rclattice.calibration import EnergyBalanceResult, energy_balance_rectangle
from rclattice.materials import concrete_uniaxial_regularized
from rclattice.problem import ConcreteGrade

from specimen import (
    A_FACE, CONCRETE, GF, GFC_FACTOR, HORIZON, L, MESH, NU, RESIDUAL_RATIO, THK,
    check_mesh_alignment, cube_problem,
)


def scaled_grade(fc_scale: float) -> ConcreteGrade:
    """The strut grade with its compressive strength scaled by `fc_scale` (D71).

    `epsc0` scales WITH `fc`, and that is forced rather than chosen. Concrete02's initial
    compressive tangent is `2*fc/epsc0` whatever the grade's `E` says (D56), so holding `epsc0`
    while raising `fc` would raise every strut's stiffness by the same factor and destroy the
    energy-balance calibration the correction exists to preserve. Deriving `epsc0 = 2*fc'/E` keeps
    the tangent at exactly `E`.

    The consequence is unavoidable and is the price of this knob: each strut now peaks at
    `fc_scale*epsc0` instead of `epsc0`. Concrete02's parabola ties (E, peak stress, peak strain)
    together, so no single scaling can fix the peak STRESS and leave the peak STRAIN alone.

    `fcu` scales too, so the residual plateau keeps its ratio to `fc`. `ft` deliberately does NOT:
    it is a measured property, and it governs the transverse splitting that causes the shortfall in
    the first place — raising it is a different knob with different consequences (per-orientation
    tension), not part of this one.
    """
    if fc_scale == 1.0:
        return CONCRETE
    fc = CONCRETE.fc * fc_scale
    return replace(CONCRETE, name=f"{CONCRETE.name}_fcx{fc_scale:.3f}", fc=fc,
                   epsc0=2.0 * fc / CONCRETE.E, fcu=CONCRETE.fcu * fc_scale)


def calibrate(*, mesh_size: float = MESH, horizon: float = HORIZON) -> EnergyBalanceResult:
    """Aydin's energy balance over the cube face itself — the member being homogenised."""
    return energy_balance_rectangle(L, L, mesh_size, E=CONCRETE.E, nu=NU, thickness=THK,
                                    horizon=horizon)


def cube_lattice(area: float, *, mesh_size: float = MESH, horizon: float = HORIZON,
                 Gf: float = GF, gfc_factor: float = GFC_FACTOR,
                 residual_ratio: float = RESIDUAL_RATIO, fc_scale: float = 1.0):
    """The Aydin-calibrated Concrete02 cube lattice — plain concrete, standard `Truss` struts.

    Strut areas come from the ELASTIC energy balance unchanged: that calibration fixes the initial
    tangent, and nothing after cracking alters it.

    Both softening branches are length-regularized, Aydin's crack-band principle: the tension slope
    from `Ets = ft^2*L/(2*Gf)` and the crushing strain from `epsU = epsc0 + 2*Gfc/((fc+fcu)*L)`, so
    the energy dissipated per unit crack area is the same for every strut regardless of grid spacing.
    (Aydin scales his multilinear strain multipliers by d/L for the same reason; his backbone is
    trilinear where Concrete02's is bilinear.)

    `residual_ratio` and `gfc_factor` set the post-peak plateau and the softening slope directly, so
    they are RESULTS here rather than the numerical stabilisers they were in the column and wall
    (D22); `pushover.py` prints them with every run.

    `fc_scale` applies the strength correction of D71: the strut grade's `fc` is multiplied by it
    while the calibrated `area` — and therefore K0 — is left alone, so stiffness and strength are
    set by independent knobs. Use `strength_correction_factor` to derive it; 1.0 is off.
    """
    check_mesh_alignment(mesh_size)
    grade = scaled_grade(fc_scale)

    def material_for(_zone: str, length: float):
        return concrete_uniaxial_regularized(grade, 0, length, Gf=Gf, Gfc=gfc_factor * Gf,
                                             residual_ratio=residual_ratio)

    model, _edges = build_lattice_rc(
        cube_problem(mesh_size), mesh_size,
        material_for=material_for,
        zone_of=lambda _x, _y: "concrete",
        rebars=(),
        strut_area=area,
        horizon=horizon,
        strut_element="Truss",
    )
    return model


def expected_modulus(cal: EnergyBalanceResult) -> float:
    """The unconfined axial modulus this calibrated lattice should ACTUALLY show (MPa) — the
    verification target for the measured initial tangent (D53).

    It is NOT the grade's E, and why not is the main elastic result of this example.

    Aydin's Eq. 2.1 equates energies under a field with the transverse strain RESTRAINED
    (eps_x = e, eps_y = 0), so what it pins is the CONFINED modulus E/(1-nu^2) — numerically exact,
    verified to 0.1%. But a cube between smooth platens expands freely, and an unconfined test
    measures `E_confined * (1 - nu_lattice^2)`, using the LATTICE's Poisson ratio rather than the
    concrete value fed into Eq. 2.1:

        E_free = E * (1 - nu_lattice^2) / (1 - nu_concrete^2)

    At horizon 1.5 nu_lattice ~ 0.41, not the ~0.18 `nu_consistent` reports, so this cube should
    measure about 0.87*E. The shortfall is a property of the calibration, not a solver artifact: it
    is mesh-independent, and narrows to ~0.94*E at horizon 3.01 where the lattice is nearer isotropic.
    """
    return CONCRETE.E * (1.0 - cal.nu_effective ** 2) / (1.0 - NU ** 2)


def secant_factor(fit_max_strain: float, *, epsc0: float | None = None) -> float:
    """How far below the true tangent a straight-line fit over [0, `fit_max_strain`] lands.

    Concrete02's pre-peak compression branch is parabolic, sigma = fc*(2r - r^2) with r = eps/epsc0,
    so its SECANT at strain eps is E*(1 - eps/(2*epsc0)), not E. A least-squares-through-origin fit
    over roughly uniformly spaced points in [0, a] returns E*(1 - 3a/(8*epsc0)). Dividing it out
    leaves a check on the calibration and the drive rate alone, with the known material curvature
    removed rather than silently absorbed into the verdict.

    `epsc0` is the STRUT's, which is the grade's only when the D71 correction is off — the
    correction scales epsc0 with fc, flattening the parabola over the fit window. Reading the
    material value there would bias the modulus verdict by ~1.5% at the default fit window, and
    that verdict is quoted to four decimals.
    """
    e0 = CONCRETE.epsc0 if epsc0 is None else epsc0
    return 1.0 - 3.0 * fit_max_strain / (8.0 * e0)


def vertical_strut_capacity(model, y_cut: float, *, fc_strut: float | None = None,
                            tol: float = 1e-6) -> tuple[int, float, float]:
    """What the VERTICAL struts crossing `y_cut` can carry at fc: (count, force N, stress MPa).

    Counted from the MODEL rather than assumed: at horizon 1.5 exactly one vertical strut crosses
    per node column (`L/mesh + 1`), but at horizon 3.01 vertical struts span one, two AND three
    rows, so several cross the same cut and an `L/mesh + 1` formula understates the section badly.
    Areas are read off each element's own args, so a non-uniform `strut_area` also works.
    """
    total_area, n = 0.0, 0
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        if (ya - y_cut) * (yb - y_cut) >= 0.0 or abs(xb - xa) > tol:
            continue
        total_area += float(e.args[0])
        n += 1
    force = total_area * (CONCRETE.fc if fc_strut is None else fc_strut)
    return n, force, force / A_FACE


def strength_correction_factor(model, y_cut: float, *, tol: float = 1e-6) -> float:
    """The factor to multiply the strut grade's `fc` by so the LATTICE peaks at the material fc (D71).

    DERIVED, not fitted. At horizon 1.5 the peak arrives exactly when the vertical struts reach fc,
    with the inclined path already fully shed (measured `peak/(verticals at fc)` = 1.000, D54), so
    the deficit is precisely the verticals' share of the gross face — 0.634 at mesh 20 — and the
    factor is its reciprocal, ~1.58. Nothing is tuned to a measured curve: it is read off the strut
    list and the calibrated area, before the analysis runs.

    Consistency argument for scaling `fc` at all: this method already accepts that a strut's area is
    not a physical area (crossing struts total 1.79x the cut face, because matching stiffness across
    many directions demands more area than a continuum has). A strut is a homogenization device, not
    a piece of concrete, so requiring its STRENGTH to be the material strength while its AREA is
    already not the material area is the inconsistency — not this correction.

    LIMIT, and it is a real one: `verticals` is the right denominator only where the diagonals shed
    completely before peak. At horizon 3.01 they still hold 22% at peak and the peak arrives at
    0.48*epsc0, before the verticals reach fc, so this factor over-corrects there. `pushover.py`
    prints the measured diagonal share at peak so the assumption is checked on every run rather
    than assumed.
    """
    _n, _force, cap_stress = vertical_strut_capacity(model, y_cut, tol=tol)
    return CONCRETE.fc / cap_stress


def capacity_accounting(model, y_cut: float, nu: float, *, fc_strut: float | None = None,
                        tol: float = 1e-6) -> dict:
    """Why the axial capacity is not `fc * A_face`, as four numbers over one horizontal cut (D55).

    A continuum resists `fc * A_face` because material fills the whole cut plane. A lattice resists
    only what its struts transmit, and a strut at angle theta from vertical contributes
    `sigma * A * cos(theta)`. So the capacity is `sum(sigma_i * A_i * cos_i)` over crossing struts —
    and the three ways that sum can differ from the continuum's are separated here:

      `continuum`   fc * A_face — the reference.
      `all_at_fc`   every crossing strut simultaneously at fc, a lower-bound / strut-and-tie style
                    idealization. It comes out ABOVE the continuum: matching STIFFNESS across many
                    directions requires more total strut area than a continuum has, since each strut
                    is only partly aligned with any one load (the D40/D41 over-provisioning).
      `compatible`  strains from the affine field at eps_y = -epsc0, eps_x = +nu*epsc0, each strut's
                    stress read off Concrete02's pre-peak parabola. Compatibility, not strength,
                    decides the sharing: a diagonal strains at only `(1-nu)/2` of the verticals.
                    Struts the field puts in TENSION are excluded — they crack rather than carry.
      `verticals`   the vertical struts alone at fc — what is left once transverse cracking has shed
                    the inclined load path entirely.

    The gap between `compatible` and `verticals` is the cost of splitting, and it dominates: without
    it the lattice would land close to the continuum.

    `fc_strut` is the struts' compressive strength, which is the grade's `fc` only when the D71
    correction is off. `continuum` always uses the MATERIAL fc, because it is the target being
    matched — so under correction `verticals/continuum` should read 1.00, which is the correction's
    own self-check.
    """
    fcs = CONCRETE.fc if fc_strut is None else fc_strut
    out = {"continuum": CONCRETE.fc * A_FACE, "all_at_fc": 0.0, "compatible": 0.0, "verticals": 0.0,
           "area_vertical": 0.0, "area_crossing": 0.0}
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        if (ya - y_cut) * (yb - y_cut) >= 0.0:
            continue
        dx, dy = xb - xa, yb - ya
        lsq = dx * dx + dy * dy
        area = float(e.args[0])
        cos = abs(dy) / lsq ** 0.5
        out["all_at_fc"] += area * cos * fcs
        out["area_crossing"] += area

        f = (nu * dx * dx - dy * dy) / lsq       # affine strain / epsc0; negative = compression
        if f < 0.0:
            r = min(-f, 1.0)
            out["compatible"] += area * cos * fcs * (2.0 * r - r * r)
        if abs(dx) < tol:
            out["verticals"] += area * fcs
            out["area_vertical"] += area
    return out
