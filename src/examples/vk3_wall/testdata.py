"""Measured VK3 results, transcribed EXACTLY from Chapter 5 of Bimschas (2010).

Bimschas, M. (2010), "Displacement Based Seismic Assessment of Existing Bridges in Regions of
Moderate Seismicity", IBK Bericht Nr. 326, Institut fuer Baustatik und Konstruktion, ETH Zuerich;
vdf Hochschulverlag. doi:10.3929/ethz-a-006237119.

Everything here is a number the chapter prints — in a table, in a figure legend, or in the
narrative of Sec. 5.3.1 — so it carries no digitization error. The measured HYSTERESIS is not here:
it is a dense curve, so it lives in `digitize.py` -> `data/vk3_fig513.npz`.

WHAT EACH GROUP IS FOR:
  GEOMETRY / REINFORCEMENT / MATERIALS  the specimen, for `summary.py` to check the model against
  PREDICTION    the chapter's OWN plastic-hinge prediction. NOT a measurement — but it is what the
                loading protocol was built on (F_y' sets the elastic cycle amplitudes), and it is
                the natural yardstick for a section-equilibrium quantity like peak strength
  MEASURED      what the pier actually did: yield displacement, stiffness, ductility at failure
  LOAD_STEPS    load step -> displacement, drift and REAL ductility, so a model state can be
                compared at the same instant the chapter reports damage
  DAMAGE        drift at each observed damage event, with what the model can and cannot be asked
                to reproduce marked explicitly
  COMPONENTS    the shear / flexure / base-crack split of Fig. 5.19-right — the comparison this
                study exists for, and the one neither other wall package has

TWO INTERNAL INCONSISTENCIES IN THE SOURCE, both recorded rather than silently resolved:
  * f_c: Sec. 5.2.1c's text gives VK1/VK2 as 39/35 MPa, Tab. 5.4 gives 35/39 — swapped. VK3 is
    34 MPa in BOTH, so this specimen is unaffected; it is a caution about the chapter, not a
    problem here.
  * f_s,l: the text says 520 MPa, Tab. 5.4 says 515. The table is used.

Units: mm, kN, kN.m, MPa; curvature in 1/km. Displacements positive to SOUTH, the direction the
actuator pushed first in every cycle.
"""

from __future__ import annotations

UNIT = "VK3"
L_V = 3300.0                 # shear span, mm — the chapter's own drift denominator


def drift(displacement_mm: float) -> float:
    """Top displacement -> drift in %, the chapter's definition delta = Delta_top / L_v."""
    return displacement_mm / L_V * 100.0


# --- GEOMETRY (Tab. 5.1, Fig. 5.1) ---------------------------------------------------------------
GEOMETRY = {
    "L_v_mm": 3300.0, "l_w_mm": 1500.0, "b_w_mm": 350.0,
    "aspect_ratio": 2.2,
    "cover_transverse_mm": 20.0, "cover_longitudinal_mm": 26.0,
    "pier_height_mm": 3700.0,          # Fig. 5.1: 4600 overall - 900 foundation
    "foundation_l_mm": 3000.0, "foundation_h_mm": 900.0,
}

# --- REINFORCEMENT (Tab. 5.3, Fig. 5.2 right) -----------------------------------------------------
REINFORCEMENT = {
    "n_longitudinal": 42, "d_longitudinal_mm": 14.0, "rho_sl_pct": 1.23,
    "hoop_d_mm": 6.0, "hoop_spacing_mm": 200.0, "hoop_legs": 2, "rho_sw_pct": 0.08,
    "lap_splice_mm": None,             # VK3 has continuous bars; only VK2 is spliced
    # Fig. 5.2 (right): 2x17phi14 at s = 80 over 1280 mm, plus 2x4phi14 at the section ends.
    "web_layout": "2x17phi14 @ 80 mm over 1280 mm",
    "end_layout": "2x4phi14 at each end, stacked through the 284 mm thickness",
}

# --- MATERIALS (Tab. 5.4) --------------------------------------------------------------------------
# The chapter gives NO E_c, NO f_t and NO cracking moment, so `specimen` supplies all three by
# convention. That is the single biggest data gap against WSH3, whose paper printed E_c and whose
# M_cr pinned f_t.
MATERIALS = {
    "fc_cyl_MPa": 34.0,                # 150x300 cylinders at the day of testing
    "fs_l_MPa": 515.0, "ft_l_MPa": 630.0, "eps_su_l": 0.126,
    "fs_w_MPa": 518.0, "ft_w_MPa": 681.0, "eps_su_w": 0.084,
    "yield_plateau_to": 0.025,         # Sec. 5.2.3a — the longitudinal steel's plateau
    "Ec_MPa": None, "ft_MPa": None, "M_cr_kNm": None,     # never stated
}
FS_L_TEXT_MPa = 520.0                  # the Sec. 5.2.1c text value; Tab. 5.4's 515 is used

# --- AXIAL LOAD (Tab. 5.7) -------------------------------------------------------------------------
AXIAL = {"N_top_kN": 1300.0, "N_base_kN": 1370.0, "nu_base": 0.077}

# --- PREDICTION: the chapter's own plastic-hinge analysis (Tab. 5.5, Tab. 5.6) ---------------------
# Numerically computed, NOT measured. F_y' is what the elastic cycles of the protocol were scaled
# on, so it enters `specimen.FY_PRIME` directly.
PREDICTION = {
    "My_prime_kNm": 2124.0, "Mn_kNm": 2808.0,
    "phi_y_prime_per_km": 2.54, "phi_y_per_km": 3.36, "phi_u_per_km": 14.5,
    "Fy_prime_kN": 644.0, "Fn_kN": 851.0,
    "delta_y_prime_mm": 10.1, "delta_y_mm": 13.4, "delta_u_mm": 28.8,
    "k0_MN_per_m": 64.0,
    "eps_cu_for_phi_u": -0.005,
}

# --- MEASURED (Tab. 5.8, Tab. 5.9, Tab. 5.10) ------------------------------------------------------
# NOTE k0 is a SECANT to first yield, not an elastic stiffness: 60 MN/m is about 0.26x the uncracked
# transformed-section cantilever. Comparing an uncracked elastic lattice against it would be a
# category error — `elastic.py` checks against the closed-form cantilever instead.
MEASURED = {
    "delta_y_prime_south_mm": 10.4, "delta_y_prime_north_mm": 10.9, "delta_y_prime_mm": 10.7,
    "delta_y_mm": 14.1,                # nominal yield displacement from the 1st-cycle average
    "k0_MN_per_m": 60.0,               # secant to first yield — the LOWEST of the three units,
                                       # despite VK3 having the HIGHEST flexural strength (Sec. 5.3.2a)
    "delta_u_mm": 52.5,                # failure level; 2nd cycle to mu_prov = 5
    "drift_u_pct": 1.59,
    "mu_delta": 3.7,                   # real displacement ductility at failure, Tab. 5.9
    "mu_delta_80pct": 3.0,             # ductility at the 80%-of-peak criterion (42 mm, 2nd cycle)
    "delta_y_provisional_mm": 10.5,    # shared by VK1-VK3 to keep their load steps comparable
}

# Peak base shear is NOT tabulated for VK3 anywhere in the chapter — Tab. 5.5's equivalent does not
# exist here. It comes from the digitized Fig. 5.13 loops instead (`digitize.py`), which is why that
# digitization is the deliverable rather than a convenience.
PEAK_SHEAR_SOURCE = "digitized from Fig. 5.13 (right panel); not printed in the chapter"

# --- LOAD STEPS (Tab. 5.9) -------------------------------------------------------------------------
# (top displacement mm, mu_prov, real mu for VK3). The provisional ductilities drove the test; the
# real ones are 0.74x because VK3's own nominal yield displacement is 14.1 mm, not the 10.5 used.
LOAD_STEPS = (
    (10.5, 1.0, 0.74), (15.75, 1.5, 1.1), (21.0, 2.0, 1.5),
    (31.5, 3.0, 2.2), (42.0, 4.0, 3.0), (52.5, 5.0, 3.7),
)

# --- DAMAGE (Sec. 5.3.1) ---------------------------------------------------------------------------
# (drift %, what was observed, whether a perfect-bond Concrete02/Steel02 lattice can show it)
DAMAGE = (
    (0.08, "uncracked; an almost invisible base crack from decompression and the absence of "
           "tension strength in the construction joint", "partly — the model has no joint, so its "
           "base crack forms only when f_t is reached"),
    (0.16, "first flexural cracks", "yes"),
    (0.24, "flexural cracks begin to incline; first shear cracks; widths <= 0.3 mm", "yes"),
    (0.48, "complete crack pattern developed; widths to 0.6 mm, one shear crack 0.9 mm", "yes"),
    (0.64, "minor longitudinal cracks in the compression zones — first sign of critical compressive "
           "strain; base crack 1-1.5 mm", "partly — as strut crushing, not as a splitting crack"),
    (0.95, "pronounced spalling of cover concrete plus onset of longitudinal bar buckling; some "
           "strength degradation in the 2nd cycle; shear cracks to 1.8 mm now exceed the flexural "
           "cracks", "NO for buckling and spalling; the crack pattern yes"),
    (1.27, "severe damage to the compression zones including core spalling; pronounced bar "
           "buckling; two diagonal shear cracks open -> onset of shear-critical behaviour; 80%%-of-"
           "peak criterion reached in the 2nd cycle", "NO — bar buckling is the driver"),
    (1.59, "complete SHEAR failure in the 2nd cycle with coinciding complete loss of AXIAL load "
           "capacity; pier separates into four wedges sliding on crossed diagonal cracks; the "
           "diagonal crack opens 5 mm while flexural cracks close to <= 0.5 mm",
     "NO — sliding on an open crack needs aggregate interlock the lattice has no element for"),
)

# --- DEFORMATION COMPONENTS (Fig. 5.19 right) -------------------------------------------------------
# The chapter's own split of the top displacement, measured from the Demec grid. This is the
# comparison `gauge.py --components` targets, and the reason `node_history` learned two dofs (D68).
#
# The chapter states these as ranges in Sec. 5.3.2a rather than tabulating them, so they are stored
# as the ranges they are.
COMPONENTS = {
    "shear_pct_range": (20.0, 22.0),        # highest of the three units, matching VK3's shear demand
    "flexure_pct_range": (80.0, 90.0),
    "base_crack_pct_range": (10.0, 20.0),   # part of flexure; stable through the inelastic range
    "base_sliding": "negligible — hardly exceeded the measuring tolerance",
    "note": "shear stays a roughly CONSTANT fraction through the inelastic range, i.e. the shear "
            "mechanism softens at about the same rate as the flexural one",
}

# --- FAILURE (Sec. 5.3.1e, Sec. 5.4) ----------------------------------------------------------------
FAILURE = {
    "mode": "combined shear and axial-load failure",
    "drift_pct": 1.59,
    "cycle": "2nd cycle to mu_prov = 5",
    "mechanism": "sliding on two crossed diagonal cracks separating the pier into four wedges; the "
                 "axial load pushed the upper wedge down and drove the side wedges out",
    # The chapter's own interpretation, and the one a lattice is best placed to test:
    "trigger": "NOT diagonal tension and NOT web crushing. The diagonal compression strut lost the "
               "base compression zone that supported it, and the 0.08% transverse steel was far "
               "too little to carry the shear across the diagonal crack by an alternative "
               "strut-and-tie path (Sec. 5.4).",
}


def load_step(displacement_mm: float):
    """The LOAD_STEPS entry nearest a displacement — for labelling a model state honestly."""
    return min(LOAD_STEPS, key=lambda r: abs(r[0] - displacement_mm))
