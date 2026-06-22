"""Material mapping layer (D15): physical grade -> per-builder OpenSees material.

This is also where lattice calibration will eventually live (e.g. fracture-energy
regularization of strut softening). For now only the elastic mappings are implemented;
the nonlinear bundle (Concrete02 / ASDConcrete3D / Steel02) will be added here.
"""

from __future__ import annotations

from .model import NDMaterial, UniaxialMaterial
from .problem import ConcreteGrade, SteelGrade


def concrete_uniaxial_elastic(grade: ConcreteGrade, tag: int) -> UniaxialMaterial:
    """Uniaxial Elastic material (lattice struts, beam-column fibers)."""
    return UniaxialMaterial(tag, "Elastic", (grade.E,))


def concrete_nd_elastic(grade: ConcreteGrade, tag: int) -> NDMaterial:
    """nD ElasticIsotropic material (continuum solids / plane elements)."""
    return NDMaterial(tag, "ElasticIsotropic", (grade.E, grade.nu))


def concrete_nd_nonlinear(
    grade: ConcreteGrade,
    base_tag: int,
    wrapper_tag: int,
    *,
    lch: float,
    Gf: float | None = None,
    Gfc: float | None = None,
    tension_residual: float = 0.02,
    max_damage: float = 0.95,
    plastic_frac: float = 0.0,
) -> tuple[NDMaterial, NDMaterial]:
    """Nonlinear nD concrete for plane-stress quads (D29): ASDConcrete3D + a PlaneStress wrapper.

    The continuum analog of the lattice's length-regularized uniaxial Concrete02 (D20): the SAME
    physical grade (E, nu, ft, fc, epsc0, fcu, epsU) maps here to a 3D damage law whose tension and
    compression softening branches are crack-band regularized by the element characteristic length
    `lch` (the quad size), so dissipation is mesh-objective and matches the lattice's per-strut
    regularization. ASDConcrete3D is driven by uniaxial curves given as POSITIVE magnitudes:
      - tension `-Te/-Ts/-Td`: elastic to `ft` at `eps_cr=ft/E`, then linear softening to a small
        residual `tension_residual*ft` at `eps_tu = eps_cr + 2*Gf/(ft*lch)`;
      - compression `-Ce/-Cs/-Cd`: elastic to ~0.4*fc, peak `fc` at `epsc0`, softening to `fcu` at
        `eps_cu = max(epsc0 + 2*Gfc/((fc+fcu)*lch), epsU)`.
    `Gf` defaults to DEFAULT_GF, `Gfc` to DEFAULT_GFC_FACTOR*Gf (same as the strut law). Returns the
    (ASDConcrete3D base, PlaneStress wrapper) pair; the quad element uses the wrapper tag.

    HYSTERESIS / `plastic_frac` (D30) — controls the damage↔plasticity split, hence the cyclic
    unloading and the match to the lattice's Concrete02 under DYNAMIC loading (irrelevant to a
    monotonic pushover, which follows the backbone Ts/Cs regardless). ASDConcrete3D back-computes the
    plastic strain from the supplied damage as `eps_p = eps - sigma/((1-d)*E)`, so:
      - `plastic_frac=0` → `d = clip(1 - sigma/(E*eps), 0, max_damage)`: PURE isotropic damage,
        unloading toward the origin, NO plastic strain. A cyclic single-quad coupon shows this is the
        closest available match to Concrete02 (dissipated energy ~108%, comparable residual strain);
      - `plastic_frac>0` scales the damage DOWN (`d *= 1-plastic_frac`), introducing residual strain /
        plasticity — biases toward residual-drift matching but OVERSHOOTS Concrete02's dissipation
        (~134% at 0.3). Pure damage (the default) is the recommended dynamic-hysteresis configuration.
    """
    if grade.fc is None or grade.epsc0 is None or grade.fcu is None or grade.epsU is None:
        raise ValueError(f"grade {grade.name!r} lacks nonlinear params (fc/epsc0/fcu/epsU)")
    E, nu, fc, ft = grade.E, grade.nu, grade.fc, grade.ft if grade.ft is not None else 0.1 * grade.fc
    epsc0, fcu, epsU = grade.epsc0, grade.fcu, grade.epsU
    gf = Gf if Gf is not None else (grade.Gf if grade.Gf is not None else DEFAULT_GF)
    gfc = Gfc if Gfc is not None else DEFAULT_GFC_FACTOR * gf

    def damage(eps: list[float], sig: list[float]) -> list[float]:
        return [0.0 if e <= 0.0 else (1.0 - plastic_frac) * max(0.0, min(max_damage, 1.0 - s / (E * e)))
                for e, s in zip(eps, sig)]

    eps_cr = ft / E
    eps_tu = eps_cr + 2.0 * gf / (ft * lch)
    Te = [0.0, eps_cr, eps_tu]
    Ts = [0.0, ft, tension_residual * ft]
    Td = damage(Te, Ts)

    eps_a = 0.4 * fc / E   # strain where the elastic line (slope E) reaches 0.4*fc (< epsc0 always,
    #                        since the initial tangent E exceeds the secant fc/epsc0 for concrete);
    #                        using 0.4*epsc0 would overshoot fc because E*0.4*epsc0 > fc.
    eps_cu = max(epsc0 + 2.0 * gfc / ((fc + fcu) * lch), epsU)
    Ce = [0.0, eps_a, epsc0, eps_cu]
    Cs = [0.0, 0.4 * fc, fc, fcu]
    Cd = damage(Ce, Cs)

    base = NDMaterial(base_tag, "ASDConcrete3D",
                      (E, nu, "-Te", *Te, "-Ts", *Ts, "-Td", *Td, "-Ce", *Ce, "-Cs", *Cs, "-Cd", *Cd))
    wrapper = NDMaterial(wrapper_tag, "PlaneStress", (base_tag,))
    return base, wrapper


def concrete_uniaxial_nonlinear(grade: ConcreteGrade, tag: int) -> UniaxialMaterial:
    """Uniaxial Concrete02 for lattice struts (D19, fork: tension + softening).

    OpenSees Concrete02 args: (fpc, epsc0, fpcu, epsU, lambda, ft, Ets), with compression
    negative. Tensile strength `ft` and softening slope `Ets` default to ~0.1*fc and ~0.1*E
    when the grade leaves them None — a modest tension branch that keeps axial truss struts
    stable (a compression-only law forms mechanisms; D4 caveat).
    """
    if grade.fc is None or grade.epsc0 is None or grade.fcu is None or grade.epsU is None:
        raise ValueError(f"grade {grade.name!r} lacks Concrete02 params (fc/epsc0/fcu/epsU)")
    ft = grade.ft if grade.ft is not None else 0.1 * grade.fc
    ets = grade.Ets if grade.Ets is not None else 0.1 * grade.E
    args = (-grade.fc, -grade.epsc0, -grade.fcu, -grade.epsU, grade.lam, ft, ets)
    return UniaxialMaterial(tag, "Concrete02", args)


def steel_uniaxial(grade: SteelGrade, tag: int) -> UniaxialMaterial:
    """Uniaxial Steel02 for rebar struts (D19). Args: (Fy, E0, b, R0, cR1, cR2)."""
    return UniaxialMaterial(tag, "Steel02", (grade.fy, grade.E0, grade.b, grade.R0, grade.cR1, grade.cR2))


def steel_uniaxial_elastic(grade: SteelGrade, tag: int) -> UniaxialMaterial:
    """Uniaxial Elastic material for rebar struts / steel fibers (linear-elastic studies).

    Drop-in replacement for `steel_uniaxial` (same `(grade, tag)` signature) that emits a linear
    `Elastic` material at the steel modulus `E0` — used by the linear-material verification where
    both the lattice and the fiber beam-column are kept fully elastic to isolate the elastic
    dynamic equivalence from any constitutive difference."""
    return UniaxialMaterial(tag, "Elastic", (grade.E0,))


# default concrete fracture energies (kip, in); Gf ~ 0.1 N/mm, compression Gfc ~ 250*Gf (D20)
DEFAULT_GF = 6.0e-4
DEFAULT_GFC_FACTOR = 250.0


def concrete_uniaxial_regularized(
    grade: ConcreteGrade,
    tag: int,
    length: float,
    *,
    Gf: float | None = None,
    Gfc: float | None = None,
    max_ets_ratio: float = 0.5,
    residual_ratio: float = 0.2,
) -> UniaxialMaterial:
    """Length-regularized Concrete02 for a lattice strut (crack-band / fracture energy, D20).

    Both softening branches are regularized by the strut `length` L so dissipation is
    mesh-objective and, for small struts, gentle enough to stay convergent past yield:
      - tension softening slope `Ets = ft^2 * L / (2*Gf)`, capped at `max_ets_ratio*E` (no snap-back);
      - compression crushing strain `epsU = epsc0 + 2*Gfc/((fc+fcu)*L)`, never below the grade's
        nominal epsU.
    `residual_ratio` floors the crushing strength at `residual_ratio*fc` (D22): beyond epsU,
    Concrete02 holds fpcu as a FLAT residual, so a crushed strut keeps positive stiffness instead
    of dropping to zero — this removes the zero-tangent local mechanism that otherwise terminates
    the lattice pushover just past yield (set 0.0 to recover the raw grade crushing strength).
    `Gf` (tension) defaults to DEFAULT_GF; `Gfc` (compression) to DEFAULT_GFC_FACTOR*Gf.
    """
    if grade.fc is None or grade.epsc0 is None or grade.fcu is None or grade.epsU is None:
        raise ValueError(f"grade {grade.name!r} lacks Concrete02 params (fc/epsc0/fcu/epsU)")
    gf = Gf if Gf is not None else (grade.Gf if grade.Gf is not None else DEFAULT_GF)
    gfc = Gfc if Gfc is not None else DEFAULT_GFC_FACTOR * gf
    ft = grade.ft if grade.ft is not None else 0.1 * grade.fc
    fcu = max(grade.fcu, residual_ratio * grade.fc)  # residual compression plateau (D22)

    ets = min(ft * ft * length / (2.0 * gf), max_ets_ratio * grade.E)
    epsU = max(grade.epsc0 + 2.0 * gfc / ((grade.fc + fcu) * length), grade.epsU)
    args = (-grade.fc, -grade.epsc0, -fcu, -epsU, grade.lam, ft, ets)
    return UniaxialMaterial(tag, "Concrete02", args)
