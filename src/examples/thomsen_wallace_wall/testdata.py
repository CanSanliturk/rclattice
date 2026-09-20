"""Every number the 2019 paper prints for Thomsen & Wallace's RW2, with its source.

Companion to `specimen.py`, which holds what the MODEL needs; this holds what the PAPER SAYS, so a
comparison can never quietly drift away from the published values. Same role as
`aydin_aldemir_wall/testdata.py`, and with the same limitation: the primary sources (Thomsen &
Wallace 1995 report, 2004 JSE paper) are NOT in the repo, so every "measured" value below is as
REPORTED BY the 2019 lattice paper. What the 2019 paper does not print is listed under
`NOT_REPORTED` rather than guessed.

Run this module to print the table.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Value:
    name: str
    value: float | str
    unit: str
    source: str
    note: str = ""


# --- measured (as reported in the 2019 paper) ----------------------------------------------------
MEASURED = (
    Value("initial stiffness", 35.19, "kN/mm", "Table 4",
          "EXCEEDS the uncracked gross-section cantilever (~27 kN/mm with the printed E and "
          "dimensions), so it cannot be a top-displacement secant of this geometry as printed — "
          "see summary.py; the same kind of anomaly that exposed Aldemir's thickness (D73)"),
    Value("maximum force", 163.284, "kN", "Table 4",
          "the digitized Fig. 9(b) envelope peaks at +163.8 / -162.2 kN (1.003x), at ~55 mm"),
    Value("maximum displacement", 72.0, "mm", "Table 2",
          "the value the simulation was run to; the plotted 'experiment' marker sits at +70.9 mm "
          "and ~144 kN, i.e. the test's LAST point, not its peak"),
    Value("axial load N", 378.0, "kN", "Table 1", "0.071 A_g f_c with the printed f_c"),
)

# --- specimen geometry and reinforcement (Fig. 9a and 9c) ----------------------------------------
GEOMETRY = (
    Value("length", 1220.0, "mm", "Fig. 9(a)", "48 in"),
    Value("height to load", 3660.0, "mm", "Fig. 9(a)", "144 in; aspect ratio 3.0 ('1:3')"),
    Value("thickness", 102.0, "mm", "Fig. 9(c)", "4 in"),
    Value("boundary bars", "8-#3 per end, 3@51 from a 19 cover", "", "Fig. 9(a)",
          "four positions x two curtains at x = 19, 70, 121, 172 (mirrored)"),
    Value("boundary hoops", "2 x d_b 4.76 @ 76", "mm", "Fig. 9(a)", "3/16 in wire at 3 in"),
    Value("web vertical bars", "8-#2", "", "Fig. 9(a)",
          "four positions x two curtains; spacing INFERRED as 191 centred between the boundary "
          "elements (the drawn chain 19|153|3@191|153|19 sums to 917 of 1220)"),
    Value("web horizontal bars", "2-#2 @ 191", "mm", "Fig. 9(a)", "start height not dimensioned"),
)

# --- materials (Table 1) ------------------------------------------------------------------------
MATERIALS = (
    Value("f_t", 2.03, "MPa", "Table 1", "NO footnote (the other specimens' f_t carry 'a: TS 500'); "
                                          "equals 0.31*sqrt(f_c)"),
    Value("f_c", 42.8, "MPa", "Table 1", ""),
    Value("E_t", 31.03, "GPa", "Table 1", "NO footnote; ACI 4700*sqrt(f_c) would give 30.75"),
    Value("G_f", 75.0, "N/m", "Table 1", "no source footnote"),
    Value("f_y", 414.0, "MPa", "Table 1", "= 60 ksi nominal, for every bar size"),
    Value("E_s", 200.0, "GPa", "Table 1", ""),
)

# --- Aydin's own lattice (Tables 1, 2, 4) ---------------------------------------------------------
AYDIN_MODEL = (
    Value("grid d", 19.0, "mm", "Table 1/2", "1220/19 = 64.2 and 3660/19 = 192.6 — not exact"),
    Value("horizons", "1.5d and 3.01d", "", "Table 2", ""),
    Value("particles", 12998, "", "Table 2", "horizon 1.5"),
    Value("elements (1.5d)", 51211, "", "Table 2", "concrete only, cf. D82 for Aldemir"),
    Value("elements (3.01d)", 162903, "", "Table 2", ""),
    Value("time step", 1.0e-8, "s", "Table 2", ""),
    Value("PID Kp, Ki, Kd", "1e9, 1e8, 1e4", "", "Table 1", ""),
    Value("a1, a2, a3", "1.5, 80, 350", "", "Table 1", "tension-softening knees, fitted at d = 19"),
    Value("K_sim (1.5d)", 32.26, "kN/mm", "Table 4", "0.917 of the measured 35.19"),
    Value("F_sim (1.5d)", 169.834, "kN", "Table 4", "1.040 of the measured 163.284"),
    Value("F_sim (3.01d)", 174.544, "kN", "Table 4", "1.069"),
    Value("K_sim (3.01d)", 38.38, "kN/mm", "Table 4", "1.091"),
    Value("loading", "static monotonic", "", "text, 'RC Wall Simulations'",
          "the test was cyclic; the paper simulated a monotonic push and compares to the envelope"),
)

# --- what the 2019 paper does NOT print ----------------------------------------------------------
NOT_REPORTED = (
    "the loading protocol (drift levels, cycles per level) — the test was 'increasing cyclic "
    "lateral displacement excursions'; the Fig. 9(b) cloud shows ~8 loop amplitudes",
    "measured coupon strengths of the #3, #2 and wire — only the nominal 414 MPa",
    "the steel hardening ratio or ultimate strain (the D101/D102 controls on capacity)",
    "the failure mode and the drift at which the test ended — Fig. 9(c) shows a flexural crack "
    "pattern; the record ends at ~72 mm (1.97%) still carrying ~144 kN",
    "how 'initial stiffness' in Table 4 is defined (it exceeds the uncracked cantilever)",
    "whether the loading rate or any axial-load variation matters",
    "the foundation block (the model is fixed at the wall base, as the paper's own is)",
)

# RECOLLECTION, NOT USED BY THE MODEL: RW2 in the primary sources is reported to have been cycled
# at increasing drift levels to about 2.5%, with boundary-bar buckling and crushing at the end, and
# the coupon strengths were above nominal. None of that is in the repo and none of it enters
# `specimen.py`; it is written here only so a reader knows what the missing sources would settle.


def main() -> None:
    print(__doc__.splitlines()[0])
    for title, rows in (("MEASURED (as reported)", MEASURED), ("GEOMETRY", GEOMETRY),
                        ("MATERIALS", MATERIALS), ("AYDIN'S OWN MODEL", AYDIN_MODEL)):
        print(f"\n{title}")
        for v in rows:
            val = f"{v.value:g}" if isinstance(v.value, float) else str(v.value)
            print(f"  {v.name:<24s} {val:>22s} {v.unit:<6s} {v.source:<28s} {v.note}")
    print("\nNOT REPORTED by the 2019 paper")
    for s in NOT_REPORTED:
        print(f"  - {s}")


if __name__ == "__main__":
    main()
