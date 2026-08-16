"""Model builders + stiffness calibration for the plain-concrete cube (shared by the entry script).

Two models of the SAME `specimen`, compared under uniaxial tension:
  - `beamcolumn_reference` — one displacement-based fiber beam-column (the 1D coupon, Concrete02);
  - `cube_lattice`         — gmsh nodes + horizon struts, plain Concrete02 struts (no rebar).

`calibrate_area` finds the single strut area giving the elastic lattice the same initial axial
stiffness (K0) as the beam-column, so the two curves start from the same elastic slope (the only
calibration the study asks for — the post-peak softening is left as plain Concrete02 on both sides).
"""

from __future__ import annotations

from dataclasses import replace

from rclattice.builders import build_lattice, build_lattice_rc, select_nodes
from rclattice.materials import concrete_uniaxial_nonlinear, concrete_uniaxial_regularized
from rclattice.opensees import run_dispbeamcolumn_tension, run_pushover

from specimen import (
    CONCRETE, DU, EPS, HORIZON, L, MESH, TARGET, THK, axial_loads, cube_problem,
)

# The single-element reference uses plain Concrete02 (default Ets = 0.1*E) smeared over its WHOLE
# length L, so it dissipates a fracture energy Gf_ref = ft^2 * L / (2*Ets_ref). Regularizing the
# lattice struts to this SAME Gf (crack-band, D20) makes the localized lattice crack dissipate the
# same energy per unit area, so the two GLOBAL stress-strain descending branches span the same strain
# range (reach zero together) — i.e. the lattice FOLLOWS the reference. The PEAK still differs: the
# diagonal struts add a parallel tensile path (a known lattice overstrength the K0 calibration does
# not remove) — accepted here (matching strength would need a separate strut-strength calibration).
REF_ETS = 0.1 * CONCRETE.E
REF_GF = CONCRETE.ft ** 2 * L / (2.0 * REF_ETS)

# Peak (tensile-strength) mismatch of the K0-calibrated regularized lattice vs the reference:
# lattice_peak / reference_peak, measured ≈ 4.50 / 3.00 MPa on the default regularized+dynamic run.
# The diagonal struts add a parallel tensile load path, so matching K0 (stiffness) overshoots the
# tensile STRENGTH by this factor. It is NOT applied by default; `cube_lattice(peak_correction=...)`
# opts in to scale the lattice peak onto the reference two different ways (see there).
PEAK_RATIO = 1.50


def beamcolumn_reference() -> dict:
    """Single displacement-based fiber beam-column in uniaxial tension (Concrete02, plain grade) —
    the 1D reference. The fiber section is the cube face (in-plane L x out-of-plane THK), so its
    initial axial stiffness is E*A/L. Returns the tension pushover curve {"disp","shear","converged"}
    (shear = axial base reaction, N)."""
    mat = concrete_uniaxial_nonlinear(CONCRETE, 1)
    return run_dispbeamcolumn_tension(length=L, width=THK, depth=L, material=mat, dU=DU, target=TARGET)


def calibrate_area(k_ref: float, *, horizon: float = HORIZON):
    """SCALAR calibration: the single uniform strut area giving the ELASTIC lattice the same initial
    AXIAL stiffness as `k_ref`. The elastic concrete strut (Elastic material) makes K linear in area,
    so one tiny tension solve fixes the scale. Returns (area, control_node, base_nodes). The control
    node is the top-edge centre; the base nodes are the fixed bottom edge (the reaction sum = the
    axial force)."""
    lat0, _ = build_lattice(cube_problem(), MESH, strut_area=1.0, horizon=horizon)
    ctrl = select_nodes(lat0, (L / 2.0 - EPS, L / 2.0 + EPS, L - EPS, L + EPS))[0]
    base = select_nodes(lat0, (-EPS, L + EPS, -EPS, EPS))
    r = run_pushover(lat0, lateral_loads=axial_loads(lat0), control_node=ctrl,
                     control_dof=2, dU=1e-3, target=1e-3, base_nodes=base)
    return k_ref / (r["shear"][-1] / r["disp"][-1]), ctrl, base


def cube_lattice(area, *, regularize: bool = True, Gf: float = REF_GF, peak_correction: str = "none"):
    """The calibrated cube lattice (no rebar). One plain-concrete zone; struts of either length
    (orthogonal ~MESH, diagonal ~MESH*sqrt2). Small-displacement `Truss` struts (the tension coupon
    has no large-drift / P-Delta demand).

    `peak_correction` selects an OPTIONAL strength calibration that scales the lattice tensile PEAK
    (which overshoots the reference by ~PEAK_RATIO from the diagonal-strut overstrength) onto the
    reference. Default "none" — no correction: K0 stays matched, the peak overshoots, and (with
    regularization) the descending branches still follow (D39). The two opt-ins:
      - "area": reduce the calibrated strut `area` by PEAK_RATIO. Peak matches, but area scales
        stiffness AND strength together, so this also lowers K0 by PEAK_RATIO — the elastic branch no
        longer overlaps the reference (strength match traded for the stiffness match). (D40)
      - "strength": reduce the concrete grade's tensile strength `ft` by PEAK_RATIO for the struts.
        The area (hence K0) is UNCHANGED, so the elastic branch still overlaps AND the peak matches —
        the two are calibrated by independent knobs (area -> K0, ft -> strength). (D41)

    `regularize=True` (default) length-regularizes each strut's Concrete02 softening (crack-band, D20)
    to the fracture energy `Gf` (defaults to the reference's implied `REF_GF`), so the localized
    lattice crack dissipates the same energy as the smeared single-element reference and the two
    descending branches align — the lattice FOLLOWS the reference. `regularize=False` -> plain
    Concrete02 (the D37 K0-only baseline): the crack localizes into one strut row and drops ~10x too
    steeply, so the lattice does NOT follow the reference past the peak."""
    if peak_correction not in ("none", "area", "strength"):
        raise ValueError(f"peak_correction must be none/area/strength, got {peak_correction!r}")
    if peak_correction == "area":
        area = area / PEAK_RATIO
    grade = replace(CONCRETE, ft=CONCRETE.ft / PEAK_RATIO) if peak_correction == "strength" else CONCRETE

    def material_for(zone: str, length: float):
        if regularize:
            return concrete_uniaxial_regularized(grade, 0, length, Gf=Gf, Gfc=250.0 * Gf,
                                                 residual_ratio=0.0)
        return concrete_uniaxial_nonlinear(grade, 0)

    model, _ = build_lattice_rc(
        cube_problem(), MESH, material_for=material_for,
        zone_of=lambda x, y: "concrete", rebars=(), strut_area=area,
        horizon=HORIZON, strut_element="Truss",
    )
    return model
