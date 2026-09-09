"""Material mapping layer (D15): physical grade -> per-builder OpenSees material.

This is also where lattice calibration lives (fracture-energy / crack-band regularization of
strut softening by strut length, D20). Both the elastic mappings and the nonlinear bundle are
implemented here: Concrete02 for uniaxial struts/fibers (plain and length-regularized),
ASDConcrete3D + PlaneStress for the nD continuum (D29/D30), and Steel02 for rebar (plus Elastic
variants for the linear-material studies).
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


def concrete_nd_elastic_planestress(grade: ConcreteGrade, tag: int) -> tuple[NDMaterial, NDMaterial]:
    """Elastic nD material pair for plane-stress quads (linear analog of ``concrete_nd_nonlinear``):
    ElasticIsotropic base + a PlaneStress wrapper. Tags and wrapper args are assigned by the
    builder (same contract as ``concrete_nd_nonlinear``)."""
    base = NDMaterial(tag, "ElasticIsotropic", (grade.E, grade.nu))
    wrapper = NDMaterial(tag + 1, "PlaneStress", ())   # args (base tag) filled in by the builder
    return base, wrapper


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


def steel_uniaxial_ruptured(grade: SteelGrade, tag: int, *, eps_rupture: float) -> list:
    """`Steel02` that BREAKS: the same bar law wrapped in `MinMax` at +/- `eps_rupture`.

    Plain `Steel02` has no rupture strain and no buckling, so a lattice built from it cannot lose
    reinforcement at any drift — the reason the monotonic pushes in the Aldemir study rise
    indefinitely and the cyclic envelope sags a few percent and then holds (D88/D89). Past the
    rupture strain `MinMax` returns zero stress and zero stiffness, so the bar strut drops out of
    the assembly and its force redistributes.

    RETURNS TWO MATERIALS — the base `Steel02` at `tag` and the `MinMax` wrapper at `tag + 1` that
    references it. Elements must use the LAST one; `build_lattice_rc` handles that.

    `eps_rupture` IS AN ASSUMPTION unless the source gives it. Welded wire mesh is typically
    low-ductility (Class A, ~2.5%); hot-rolled deformed bar reaches 7.5% and beyond. The value
    chosen decides the predicted drift capacity, so every run records it.
    """
    if eps_rupture <= grade.fy / grade.E0:
        raise ValueError(f"eps_rupture {eps_rupture} is at or below yield "
                         f"{grade.fy / grade.E0:.5f} — the bar would break before it yields")
    base = steel_uniaxial(grade, tag)
    return [base, UniaxialMaterial(tag + 1, "MinMax",
                                   (tag, "-min", -eps_rupture, "-max", eps_rupture))]


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


def aydin_lattice_softening(
    ft: float,
    E: float,
    length: float,
    *,
    Gf: float,
    a1: float = 1.5,
    b1: float = 0.6,
    b2: float = 0.2,
    a3_over_a2: float = 5.0,
) -> tuple[float, float]:
    """Solve Aydin's trilinear tension tail for `(a2, a3)` at the fracture energy `Gf` (D60).

    Aydin's tensile backbone (2021, Fig. 1(a)) is linear to `(eps_cr, ft)` and then trilinear
    through `(a1*eps_cr, b1*ft)` and `(a2*eps_cr, b2*ft)` to zero at `a3*eps_cr`. The SHAPE
    constants `a1, b1, b2` are fixed across every specimen in his Table 1; what varies is how far
    the tail runs, and he fits it so the simulated tension test dissipates the target `Gf` (his
    Fig. 4 inner loop, "Energy error < Tol").

    Doing that as algebra rather than as a fitted table is what makes the model mesh-objective. The
    energy a strut of length `L` dissipates per unit CRACK AREA is `L * integral(sigma d eps)` over
    the post-peak branch, so requiring it to equal `Gf` gives a2 directly:

        Phi_req = Gf / (L * ft * eps_cr)          (the required area under the normalized tail)
        Phi(a2) = C0 + C1*a2                       with a3 = a3_over_a2 * a2
        C0 = (1+b1)(a1-1)/2 - (b1+b2)*a1/2
        C1 = (b1+b2)/2 + b2*(a3_over_a2-1)/2

    Longer struts get a shorter tail, exactly as `concrete_uniaxial_regularized` steepens
    Concrete02's `Ets` with length (D20), and as Aydin scales his multipliers by `d/L`.

    Sanity check against his Table 1 (Jansen & Shah: ft = 2.55, E = 34 217, d = 10, Gf = 0.08):
    this returns a2 = 52.9, a3 = 265 against his tabulated 60 and 300 — 12% apart, which is the
    gap between an exact energy identity and his fit to the Cornelissen et al. (1986) curve.

    `a3_over_a2` defaults to 5.0, the Jansen & Shah ratio; his other specimens run 5.3-6.8.
    """
    eps_cr = ft / E
    phi_req = Gf / (length * ft * eps_cr)
    c0 = 0.5 * (1.0 + b1) * (a1 - 1.0) - 0.5 * (b1 + b2) * a1
    c1 = 0.5 * (b1 + b2) + 0.5 * b2 * (a3_over_a2 - 1.0)
    a2 = (phi_req - c0) / c1
    if a2 <= a1:
        raise ValueError(
            f"strut length {length:g} is too long for Gf = {Gf:g} and ft = {ft:g}: the trilinear "
            f"tail collapses (a2 = {a2:.3g} <= a1 = {a1:g}). Refine the grid or raise Gf — a strut "
            f"cannot dissipate less than its own peak-to-first-knee energy."
        )
    return a2, a3_over_a2 * a2


def concrete_lattice_aydin(
    grade: ConcreteGrade,
    tag: int,
    length: float,
    *,
    Gf: float,
    a1: float = 1.5,
    b1: float = 0.6,
    b2: float = 0.2,
    a3_over_a2: float = 5.0,
    a2: float | None = None,
    a3: float | None = None,
    rsm: bool = True,
    alpha: float = 1.0 / 3.0,
    beta: float = 0.4,
    eps_max: float = 0.2,
    fc_cap: float | None = None,
) -> UniaxialMaterial:
    """Aydin's TENSION-ONLY lattice strut material (2021, Fig. 1) as an `ElasticMultiLinear` (D60).

    This is the paper's central modelling idea, so read the sign convention carefully: the strut
    has **no compressive strength at all**. Compression is linear elastic at `E`, optionally with
    the reduced-stiffness model (RSM) knee, and it never peaks, never softens and never crushes.
    The specimen's compressive strength is not a material property here — it EMERGES, as the
    stability loss of a load path whose lateral restraint has been cracked away by transverse
    (indirect) tension. That is why the calibration parameter is a geometric one (`Rmax/d`,
    `mesh.perturb_nodes`) rather than a strength.

    Backbone, with `eps_cr = ft/E`:
      * tension      `(eps_cr, ft) -> (a1*eps_cr, b1*ft) -> (a2*eps_cr, b2*ft) -> (a3*eps_cr, 0)`,
                     then flat zero. `a2`/`a3` come from `aydin_lattice_softening` at `Gf`.
      * compression   slope `E`, and with `rsm=True` slope `beta*E` beyond `alpha*epsc0`, forever.

    RSM (his Fig. 1(b)) is a lateral-cracking stiffness reduction in the spirit of Vecchio &
    Collins: 40% of the modulus (`beta = 0.4`) past a third of the strain at peak stress
    (`alpha = 1/3`), where `epsc0 = 2*fc/E`. Without it the paper still gets the STRENGTH right
    but not the strain at peak or the post-cracking stiffness (his Fig. 3(b)); both curves are
    worth reproducing, hence the switch.
    NOTE the paper is self-inconsistent on these two symbols: its notation list and Fig. 1's caption
    make `alpha` the strain multiplier (0.33) and `beta` the modulus multiplier (0.4), while the
    p. 397 text says "alpha and beta were taken as 0.4 and 1/3" and then describes the opposite in
    the same sentence. The self-consistent reading — notation list, caption, and the text's own
    "40% stiffness reduction at a strain of 1/3 of epsc0" — is the one used here.
    NOTE also that RSM is applied to EVERY strut, whereas Aydin applies it only at laterally cracked
    nodes. That state-dependent switch cannot be expressed as a static material assignment; in a
    uniaxial compression test essentially every compression-carrying strut does crack laterally
    before peak, so the two coincide over the part of the response that matters.

    WHY `ElasticMultiLinear`. It is path-INDEPENDENT: stress is a pure function of strain, so a
    cracked strut returns exactly zero at zero strain. `HystereticSM` reproduces the same envelope
    but unloads on the initial stiffness, which puts a fully cracked strut at -10 MPa (RSM) or
    -44 MPa (no RSM) when its strain comes back to zero — it would push the split columns apart and
    destroy the very mechanism being modelled. The cost of the path-independent choice is the
    opposite error, with no damage memory: a strut that softens and then unloads recovers its full
    stiffness on reload, where Aydin's sequentially linear analysis holds the reduced secant. Under
    monotonic compression each strut's strain is essentially monotonic, so the two agree; this
    material is NOT suitable for cyclic work.
    """
    if grade.ft is None:
        raise ValueError(f"grade {grade.name!r} needs ft for Aydin's lattice material")
    ft, E = grade.ft, grade.E
    eps_cr = ft / E
    # `a2`/`a3` are normally SOLVED from Gf so the law is mesh-objective (a longer strut gets a
    # shorter tail). Passing them explicitly pins the published tail instead, which is what a
    # replica of a specific paper needs — and it deliberately breaks mesh objectivity, so the
    # caller owns the strut length it is valid at.
    if (a2 is None) != (a3 is None):
        raise ValueError("pass a2 and a3 together, or neither")
    if a2 is None:
        a2, a3 = aydin_lattice_softening(ft, E, length, Gf=Gf, a1=a1, b1=b1, b2=b2,
                                         a3_over_a2=a3_over_a2)

    # Compression: capped (EPP), one knee at alpha*epsc0 (RSM), or none at all. `eps_max` only has
    # to sit beyond anything the analysis reaches — ElasticMultiLinear holds the last segment's
    # stress outside its range.
    eps_max = max(eps_max, 2.0 * a3 * eps_cr)
    if fc_cap is not None:
        # `eppcomp`: linear at E to `fc_cap`, then PERFECTLY PLASTIC. Three consequences, all
        # deliberate and all reported by the study that uses it:
        #   * there is no compression softening, so Gfc and a crushing strain cease to exist as
        #     parameters and regularization applies to TENSION ONLY;
        #   * a damage figure's "crushed" class becomes "yielded in compression" — a different
        #     statement, and `figure_damage(crush_label=...)` is how it gets relabelled;
        #   * nothing can lose compressive load-carrying capacity, so a compression-driven collapse
        #     cannot occur by material failure. That is a property of an EPP branch, not an
        #     oversight — a load-path collapse can still happen, by cracking away the restraint.
        if rsm:
            raise ValueError("fc_cap and rsm are two different compression branches; pass "
                             "rsm=False with fc_cap (the cap replaces the RSM knee)")
        if fc_cap <= 0.0:
            raise ValueError(f"fc_cap must be positive, got {fc_cap}")
        eps_y = fc_cap / E
        eps_max = max(eps_max, 2.0 * eps_y)
        neg_eps = [-eps_max, -eps_y]
        neg_sig = [-fc_cap, -fc_cap]
    elif rsm:
        if grade.epsc0 is None:
            raise ValueError(f"grade {grade.name!r} needs epsc0 (= 2*fc/E) for the RSM knee")
        eps_knee = alpha * grade.epsc0
        sig_knee = E * eps_knee
        neg_eps = [-eps_max, -eps_knee]
        neg_sig = [-(sig_knee + beta * E * (eps_max - eps_knee)), -sig_knee]
    else:
        neg_eps = [-eps_max]
        neg_sig = [-E * eps_max]

    strains = neg_eps + [0.0, eps_cr, a1 * eps_cr, a2 * eps_cr, a3 * eps_cr, eps_max]
    stresses = neg_sig + [0.0, ft, b1 * ft, b2 * ft, 0.0, 0.0]
    args = ("-strain", *strains, "-stress", *stresses)
    return UniaxialMaterial(tag, "ElasticMultiLinear", args)


def concrete_lattice_aydin_cyclic(
    grade: ConcreteGrade,
    tag: int,
    length: float,
    *,
    Gf: float,
    Gfc: float | None = None,
    a1: float = 1.5,
    b1: float = 0.6,
    b2: float = 0.2,
    a3_over_a2: float = 5.0,
    residual_ratio: float = 0.2,
    tail_ratio: float = 0.01,
    pinch_x: float = 0.3,
    pinch_y: float = 0.1,
    damage1: float = 0.0,
    damage2: float = 0.0,
    beta: float = 0.0,
) -> UniaxialMaterial:
    """Aydin's TRILINEAR tension softening on a `HystereticSM`, with a REAL compressive envelope.

    This is the cyclic counterpart of `concrete_lattice_aydin` (D60), and it exists because that
    material cannot be used for cyclic work: `ElasticMultiLinear` is path-independent, so a softened
    strut recovers its full stiffness on reload and carries no damage memory at all.

    WHAT IS BORROWED FROM AYDIN, AND WHAT IS NOT (D71). Only the tension SOFTENING SHAPE is taken:
    `(eps_cr, ft) -> (a1*eps_cr, b1*ft) -> (a2*eps_cr, b2*ft) -> (a3*eps_cr, ~0)`, with `a2`/`a3`
    solved from `Gf` by `aydin_lattice_softening` and calibrated by Aydin against Cornelissen et al.
    (1986). His TENSION-ONLY compression is NOT taken: that is inseparable from the perturbed grid
    and corotTruss of his compressive mechanism (D60/D61), and this material keeps a Concrete02-like
    compressive envelope instead.

    WHY THE SHAPE MATTERS SO MUCH. At the SAME fracture energy, a trilinear tail keeps a strut
    carrying to roughly 3x the strain of Concrete02's linear softening, because the same area sits
    under a long low-stress tail rather than a short ramp. For the VK3 pier this is decisive: the
    extreme fibre needs 0.269% strain at first yield, Concrete02 at plain-concrete Gf dies at
    0.115%, and this material at the SAME Gf survives to 0.346%. The `gf_factor` inflation that was
    otherwise needed is compensating for the wrong softening shape, not for a shortfall in energy.

    PINCHING IS THE PHYSICS. `pinch_x`/`pinch_y` make a cracked strut reload along a pinched path
    rather than on the initial stiffness, which is what a crack does: it closes carrying almost
    nothing until the faces make contact. This is also what bounds the artefact that ruled
    `HystereticSM` out for the tension-only cube, where a fully cracked strut returned to -44 MPa at
    zero strain; here the compressive envelope caps it and the pinching softens the approach.
    They have no measured counterpart for any specimen in this repository and are therefore a
    calibration -- but one that shapes HYSTERESIS rather than inflating strength.
    """
    if grade.fc is None or grade.epsc0 is None or grade.fcu is None or grade.epsU is None:
        raise ValueError(f"grade {grade.name!r} lacks Concrete02 params (fc/epsc0/fcu/epsU)")
    ft = grade.ft if grade.ft is not None else 0.1 * grade.fc
    gfc = Gfc if Gfc is not None else DEFAULT_GFC_FACTOR * Gf
    eps_cr = ft / grade.E
    a2, a3 = aydin_lattice_softening(ft=ft, E=grade.E, length=length, Gf=Gf,
                                     a1=a1, b1=b1, b2=b2, a3_over_a2=a3_over_a2)
    # compression: same crack-band regularization as concrete_uniaxial_regularized, so the two
    # materials differ ONLY in the tension branch and a comparison isolates that
    fcu = max(grade.fcu, residual_ratio * grade.fc)
    epsU = max(grade.epsc0 + 2.0 * gfc / ((grade.fc + fcu) * length), grade.epsU)
    pos = [ft, eps_cr, b1 * ft, a1 * eps_cr, b2 * ft, a2 * eps_cr, tail_ratio * ft, a3 * eps_cr]
    neg = [-grade.fc, -grade.epsc0, -fcu, -epsU, -residual_ratio * grade.fc, -10.0 * epsU]
    return UniaxialMaterial(
        tag, "HystereticSM",
        ("-posEnv", *pos, "-negEnv", *neg,
         "-pinch", pinch_x, pinch_y, "-damage", damage1, damage2, "-beta", beta),
    )


def bond_elastic_brittle_damaging(
    grade: ConcreteGrade,
    tag: int,
    *,
    residual: float = 0.7,
) -> list:
    """Aydin's bond law with the BRITTLE DROP MADE IRREVERSIBLE — the cyclic-safe form (D95).

    Same backbone as `bond_elastic_brittle`, and deliberately so: elastic on the concrete strut's
    own line to `ft`, a brittle fall to `residual * ft`, then a flat plateau. What changes is that
    the fall is DAMAGE rather than a shape traced on the way out and back:

        `ElasticMultiLinear` is path-INDEPENDENT. A bond link driven past `eps_cr` and then
        unloaded recovers its full elastic branch, so a broken bond heals on every reversal. That
        is invisible in a monotonic push and fatal in a cyclic one, where bond links are 61% of the
        elements. The same objection retired both concrete laws from cyclic work (D60, D83/D86);
        it was never carried across to the bond law.

    The paper says "elastic brittle response with residual strength" (2019, Fig. 1c), and BRITTLE
    means irreversible — so this is the more faithful reading of the source, not a departure from
    it. The envelope is identical, which is what makes the two directly comparable.

    Built as a PARALLEL pair, which gives the envelope exactly rather than approximately:

      * `Elastic` at `(1 - a) * E`, wrapped in `MinMax` at +/- eps_cr — carries the part of the
        rise that DISAPPEARS at cracking, and `MinMax` never comes back;
      * `ElasticPP` at `a * E` yielding at eps_cr — carries the part that SURVIVES, and holds
        `a * ft` forever after.

    Below eps_cr the two sum to `E` exactly and reach `ft` exactly at eps_cr; above it the first is
    dead and the second holds the plateau. No `drop` ramp is needed (the `ElasticMultiLinear` form
    cannot fall vertically and fakes it over a small strain band), so the corner is exact.

    ON REVERSAL the residual branch unloads elastically at `a * E` and yields at `-a * ft`, i.e.
    the surviving bond is elastic-perfectly-plastic in both directions. That is the right shape for
    what the residual physically IS — rib-on-concrete friction once the chemical and mechanical
    adhesion has gone — and it dissipates, where the elastic form dissipates nothing.

    RETURNS FOUR MATERIALS at `tag .. tag+3`; elements must use the LAST (the `Parallel`), which
    `build_lattice_rc` handles. `residual` is Aydin's `a`, his 0.7 the default; it must stay
    non-zero for the same structural reason as in `bond_elastic_brittle` — a steel node is held
    only by its ring.
    """
    if grade.ft is None:
        raise ValueError(f"grade {grade.name!r} needs ft for the bond material")
    if not 0.0 < residual <= 1.0:
        raise ValueError(f"residual must be in (0, 1]; got {residual} "
                         "(0 leaves steel nodes unrestrained once the bond ring cracks)")
    ft, E, a = grade.ft, grade.E, float(residual)
    eps_cr = ft / E
    return [
        UniaxialMaterial(tag, "Elastic", ((1.0 - a) * E,)),
        UniaxialMaterial(tag + 1, "MinMax", (tag, "-min", -eps_cr, "-max", eps_cr)),
        UniaxialMaterial(tag + 2, "ElasticPP", (a * E, eps_cr)),
        UniaxialMaterial(tag + 3, "Parallel", (tag + 1, tag + 2)),
    ]


def bond_elastic_brittle(
    grade: ConcreteGrade,
    tag: int,
    *,
    residual: float = 0.7,
    drop: float = 1e-3,
    eps_max: float = 0.2,
    peak_force: float | None = None,
    area: float | None = None,
    length: float | None = None,
    slip_at_peak: float | None = None,
) -> UniaxialMaterial:
    """Aydin's BOND element law: elastic-brittle with a residual plateau (2019, Fig. 1c), D72.

    A bond element is a truss joining a steel node to a nearby concrete node, so its "slip" is the
    axial stretch of that link and its law is written in the same stress-strain terms as a concrete
    strut. It shares the concrete strut's cracking point — the paper draws Bond and Concrete rising
    on the same line to the same `F_cr` — and then, instead of softening to zero, DROPS to a
    residual `a * F_cr` and stays there:

        0 -> (eps_cr, ft)  -> (eps_cr*(1+drop), residual*ft) -> flat at residual*ft

    `residual` is Aydin's `a`; his 0.7 is the default, chosen (his Sec. "Reinforced Concrete Lattice
    Modeling", from Aydin 2017) to reproduce the residual bond strength of a DEFORMED bar. Lower it
    for plain bars — that residual plateau is the only thing representing bar-rib interlock after
    the surrounding concrete has cracked.

    The law is SYMMETRIC: bond resists relative slip equally either way along the link.

    NOTE the residual is load-bearing in a second, structural sense. A steel node is held in place
    ONLY by its ring of bond elements; if they were fully brittle, cracking the ring would leave the
    steel node free and the tangent singular. A non-zero `residual` is what keeps the assembly
    well-posed, so do not set it to 0.

    `drop` sets the strain increment over which the brittle fall happens (`ElasticMultiLinear` needs
    strictly increasing strain points, so the fall cannot be exactly vertical); it is relative to
    `eps_cr`. `eps_max` only has to sit beyond anything the analysis reaches — the last segment's
    stress is held outside the range.
    """
    if grade.ft is None:
        raise ValueError(f"grade {grade.name!r} needs ft for the bond material")
    if not 0.0 < residual <= 1.0:
        raise ValueError(f"residual must be in (0, 1]; got {residual} "
                         "(0 leaves steel nodes unrestrained once the bond ring cracks)")
    ft, E = grade.ft, grade.E
    eps_cr = ft / E

    # FORCE-SLIP FORM (D75). A truss couples stiffness and strength through the area, so a bond link
    # sized for the right STIFFNESS gets the wrong STRENGTH and vice versa — sizing the area alone
    # cannot satisfy both. Given the link's `area` and `length`, plus the `peak_force` it must carry
    # and the `slip_at_peak` at which it should reach it, the two are set independently:
    #
    #     ft_bond = peak_force / area          (strength)      eps_peak = slip_at_peak / length
    #     E_bond  = ft_bond / eps_peak         (stiffness)
    #
    # THIS IS THE CORRECTION THAT MATTERS. Inheriting the concrete grade's eps_cr = ft/E gives a bond
    # link that fails at ~0.004 mm of slip over a 50 mm link. A deformed-bar interface fails at
    # s1 ~ 0.1-1 mm (Model Code), two orders of magnitude larger, so a concrete-strain bond law
    # disconnects the reinforcement almost immediately — measured as a collapse at 0.026% drift
    # against 0.199% without bond at all.
    if peak_force is not None:
        if area is None or length is None:
            raise ValueError("peak_force needs `area` and `length` to become a stress and a strain")
        ft = peak_force / float(area)
        if slip_at_peak is not None:
            eps_cr = float(slip_at_peak) / float(length)

    eps_drop = eps_cr * (1.0 + drop)
    res = residual * ft
    eps_max = max(eps_max, 2.0 * eps_drop)

    strains = [-eps_max, -eps_drop, -eps_cr, 0.0, eps_cr, eps_drop, eps_max]
    stresses = [-res, -res, -ft, 0.0, ft, res, res]
    return UniaxialMaterial(tag, "ElasticMultiLinear", ("-strain", *strains, "-stress", *stresses))
