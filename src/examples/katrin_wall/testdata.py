"""Measured WSH3 results, transcribed EXACTLY from Dazio, Beyer & Bachmann (2009).

*Engineering Structures* 31:1556-1571. Everything here is a number the paper prints — in a table,
in a figure legend, or in the narrative of Sec. 3.1 — so it carries no digitization error. The two
exceptions are marked `# digitized`: Fig. 15(d) and 15(e) are marker plots with no tabulated
counterpart, and their values were read from marker centroids in the embedded raster (+-0.3 mm on
displacement, +-0.1 km^-1 on curvature).

The measured HYSTERESIS is not here: it is a dense curve, so it lives in `digitize.py` ->
`data/wsh3_fig7.npz`. This module is for everything that can be stated as a number.

WHAT EACH GROUP IS FOR, when the lattice run comes:
  PROTOCOL      the drive history — exactly what the actuator did, so the model is driven the same
  STRENGTH      peak shear and moment; the headline strength comparison
  DISPLACEMENT  yield and ultimate displacements; the headline stiffness / ductility comparison
  LOAD_STEPS    load step -> drift, so a model state can be compared at the SAME instant the paper
                reports a strain profile, a curvature profile or a photograph
  DAMAGE        drift at each observed damage event; what the model can and cannot be asked to hit
  CURVATURE     base curvature vs top displacement — the local check, comparable to `gauge.py`
  COMPONENTS    shear/flexure split; bounds how much of the response is NOT flexure

Units: mm, kN, kN.m, and strain in absolute (not per mille). Displacements positive to South,
which is the direction the actuator pushed first in every cycle.
"""

from __future__ import annotations

UNIT = "WSH3"
L_V = 4560.0                 # shear span, mm (Sec. 2.1) — the paper's own drift denominator


def drift(displacement_mm: float) -> float:
    """Top displacement -> drift in %, the paper's definition delta = Delta / L_v (Sec. 2.3)."""
    return displacement_mm / L_V * 100.0


# --- STRENGTH (Table 5) -------------------------------------------------------------------------
# V_max / M_max are the largest values MEASURED during the test; V_Rd,c and M_cr are theoretical
# EC2 values the paper computes for comparison, not measurements.
STRENGTH = {
    "V_max_kN": 454.0,       # max base shear, either direction
    "M_max_kNm": 2072.0,     # = V_max * L_v
    "M_cr_kNm": 527.0,       # EC2 cracking moment (theoretical)
    "V_Rdc_kN": 272.0,       # EC2 shear resistance without shear reinforcement (theoretical)
}

# The paper's M_cr pins the tensile strength it assumed: M_cr = (f_ctm + N/A_g) * t*l_w^2/6 with
# N/A_g = 2.287 MPa gives f_ctm = 2.97 MPa — which is EC2 0.30*f_ck^(2/3) evaluated at the
# CHARACTERISTIC strength f_ck = f'_c - 8 = 31.2, not at the mean f'_c = 39.2. Reproduced to 1% on
# WSH1/5/6 as well. Worth knowing because the lattice's cracking response scales directly with f_t.
F_CTM_IMPLIED = 2.97         # MPa, back-figured from Table 5's M_cr

# --- DISPLACEMENT (Table 4) ---------------------------------------------------------------------
# Three different yield displacements, because the paper found the first one unsatisfactory and
# re-evaluated afterwards. DELTA_Y_34 is the one the LOADING HISTORY was actually built on, so it
# is the one the protocol below must use; DELTA_Y_NOMINAL is the better estimate for judging
# ductility (Sec. 2.3).
DISPLACEMENT = {
    "delta_y_34_mm": 15.4,        # 3/4-rule; the loading history is built on this
    "delta_y_first_mm": 11.3,     # displacement at first yield of the outer longitudinal bar
    "delta_y_nominal_mm": 16.2,   # nominal yield displacement from strain limits
    "delta_u_mm": 92.4,           # ultimate (test ended at LS63)
    "mu_delta": 5.7,              # delta_u / delta_y(3/4)
}

# --- PROTOCOL (Fig. 6 + Sec. 2.3) ---------------------------------------------------------------
# Two force-controlled cycles to 0.75*F_y, then two displacement-controlled cycles at each
# ductility level mu = 2..6, each level's amplitude being mu * delta_y(3/4). Every cycle goes
# SOUTH (+) first, then NORTH (-). The actuator speed was stepped up three times.
CYCLES_PER_LEVEL = 2
DUCTILITY_LEVELS = (2, 3, 4, 5, 6)
FIRST_BLOCK_AMPLITUDE_MM = 11.4   # the 0.75*F_y cycles reached delta = 0.25% (Fig. 11c legend)

# Actuator speed in mm/min (Fig. 6). NOTE how slow this is: 0.02-0.06 mm/s, some 130x slower than
# the 7.6 mm/s of the SW-NC-FF wall. A dynamic-relaxation run cannot be driven at the real rate;
# it has to run faster and rely on damping to keep inertia out of the base shear.
RATE_MM_PER_MIN = {"0.75Fy": 1.2, 2: 1.2, 3: 2.4, 4: 2.4, 5: 3.6, 6: 3.6}


def protocol() -> tuple[float, ...]:
    """Peak displacements in order, as the actuator applied them: (+, -) per cycle.

    Returns the TURNING POINTS only. Feed straight into a reversed-cyclic driver; the paths between
    them were linear ramps at the speeds in RATE_MM_PER_MIN.
    """
    peaks: list[float] = []
    for _ in range(CYCLES_PER_LEVEL):
        peaks += [FIRST_BLOCK_AMPLITUDE_MM, -FIRST_BLOCK_AMPLITUDE_MM]
    for mu in DUCTILITY_LEVELS:
        amp = mu * DISPLACEMENT["delta_y_34_mm"]
        for _ in range(CYCLES_PER_LEVEL):
            peaks += [amp, -amp]
    return tuple(peaks)


# --- LOAD_STEPS (Fig. 6 numbering; drifts printed in the Fig. 11c legend) ------------------------
# The FIRST cycle of each level, both directions. These are the instants at which the paper reports
# strain profiles (Fig. 10c), curvature profiles (Fig. 11c) and crack photographs (Fig. 8c) — so
# they are the instants at which a model state is directly comparable.
LOAD_STEPS = {
    # LS: (direction, drift %, top displacement mm)
    5:  ("South", 0.25, +11.4), 9:  ("North", 0.25, -11.4),
    16: ("South", 0.68, +31.0), 19: ("North", 0.68, -31.0),
    26: ("South", 1.02, +46.5), 29: ("North", 1.02, -46.5),
    36: ("South", 1.36, +62.0), 39: ("North", 1.36, -62.0),
    46: ("South", 1.74, +79.3), 49: ("North", 1.70, -77.5),
    56: ("South", 2.04, +93.0), 59: ("North", 2.04, -93.0),
    63: ("North", 1.79, -81.6),   # test ended here: corner bar rupture (Sec. 3.1)
}

# --- DAMAGE (Sec. 3.1 narrative + Table 6) -------------------------------------------------------
# The drift at which each damage state was OBSERVED. A perfect-bond Steel02 lattice can be asked
# about the first two and about strength; it has no mechanism for the last two.
DAMAGE = {
    "spalling_onset": dict(drift_pct=1.02, ls=26, note="cover spalling begins; bars not yet visible"),
    "bars_visible": dict(drift_pct=1.70, ls=51, note="longitudinal bars exposed, first signs of buckling"),
    "bar_buckling": dict(drift_pct=1.70, ls=51, note="onset; worsens over the following cycles"),
    "bar_fracture": dict(drift_pct=1.79, ls=63,
                         note="one D12 corner bar ruptures ~140 mm above the base, where buckling "
                              "curvature had been largest; ends the test"),
    "max_drift": dict(drift_pct=2.04, ls=56, note="largest drift reached"),
}
# Force capacity had dropped slightly LESS than 20% when the test was stopped, so WSH3 never met
# the paper's own failure definition (Sec. 3.1) — its drift capacity is a lower bound.
FAILURE_CRITERION_MET = False

# --- Table 6: strains at limit states, WSH3 row --------------------------------------------------
# Two DIFFERENT kinds of number, and they must not be mixed:
#   "section": strains from a plane-section moment-curvature analysis, evaluated at the curvature
#              the test actually showed. These are what plastic-hinge analysis uses.
#   "demec":   strains measured directly on the corner bars over a 150 mm gauge length.
# The paper's own point (Sec. 5.1) is that these disagree systematically — measured compressive
# strains run ~25% higher — because inclined flexure-shear cracks break the plane-section
# assumption. A lattice model computes something closer to the DEMEC column.
LIMIT_STATE_STRAINS = {
    "spalling":        dict(top_mm=46.0, section=dict(eps_c=-3.0e-3),
                            demec=dict(eps_c=-3.7e-3)),
    "buckling":        dict(top_mm=78.0, section=dict(eps_c=-4.8e-3, eps_s=31.0e-3,
                                                      excursion=35.9e-3),
                            demec=dict(eps_c=-5.6e-3, eps_s=22.7e-3, excursion=28.3e-3)),
    "boundary_bar_fracture": dict(top_mm=92.0, section=dict(eps_s=40.0e-3, excursion=46.2e-3)),
}
PHI_YIELD_PER_KM = 3.0       # yield curvature, ~the same for all six units (Sec. 5.1)

# --- CURVATURE: Fig. 15(d), base curvature vs top displacement ----------------------------------
# Best-fit linear curvature profile over the plastic zone, extrapolated to the wall base — NOT a
# direct measurement at the base, which strain penetration into the foundation makes meaningless.
# The comparable model quantity is the base curvature from a vertical strain gauge line
# (`gauge.py`), fitted the same way.
BASE_CURVATURE = (       # digitized (Fig. 15d): (top displacement mm, phi_base 1/km)
    (31.0, 4.2),
    (46.5, 9.1),
    (62.0, 13.4),
    (77.9, 18.5),
    (91.2, 23.4),
)

# --- Fig. 15(e), plastic hinge length from Eq. (2) -----------------------------------------------
# Derived, not measured: L_ph = Delta_pf / (phi_p * L_v). Approximate above 70 mm, where WSH3's
# markers overlap WSH6's in the figure.
PLASTIC_HINGE_LENGTH_M = (   # digitized (Fig. 15e): (top displacement mm, L_ph m)
    (31.0, 1.20),
    (46.5, 0.77),
    (61.9, 0.74),
    (78.0, 0.70),
    (90.9, 0.66),
)

# --- COMPONENTS (Fig. 9b) ------------------------------------------------------------------------
# WSH3 has the LARGEST shear share of the six units, and it is still only ~12% of the flexural
# displacement — this is a flexure-controlled wall at a shear span ratio of 2.28. The ratio stayed
# approximately constant over the whole loading history.
SHEAR_TO_FLEXURE = 0.12


def summary() -> str:
    lines = [f"{UNIT} — measured results (Dazio, Beyer & Bachmann 2009)", ""]
    lines += [f"  V_max            {STRENGTH['V_max_kN']:.0f} kN   "
              f"M_max {STRENGTH['M_max_kNm']:.0f} kN.m   M_cr {STRENGTH['M_cr_kNm']:.0f} kN.m",
              f"  delta_y (3/4)    {DISPLACEMENT['delta_y_34_mm']:.1f} mm "
              f"({drift(DISPLACEMENT['delta_y_34_mm']):.2f}% drift)",
              f"  delta_u          {DISPLACEMENT['delta_u_mm']:.1f} mm "
              f"({drift(DISPLACEMENT['delta_u_mm']):.2f}% drift), mu = {DISPLACEMENT['mu_delta']}",
              f"  shear / flexure  {SHEAR_TO_FLEXURE:.2f}",
              "", f"  protocol: {len(protocol())} turning points, "
              f"+-{max(protocol()):.1f} mm max, {CYCLES_PER_LEVEL} cycles per level", ""]
    lines.append("  damage:")
    for name, d in DAMAGE.items():
        lines.append(f"    {name:<16} {d['drift_pct']:.2f}% drift (LS{d['ls']}) — {d['note']}")
    lines += ["", "  base curvature (Fig. 15d):"]
    for u, k in BASE_CURVATURE:
        lines.append(f"    {u:5.1f} mm ({drift(u):.2f}%)  ->  {k:.1f} /km")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
