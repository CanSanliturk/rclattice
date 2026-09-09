"""Every number the 2019 paper prints for the Aldemir et al. (2017) wall, with its source.

Companion to `specimen.py`, which holds what the MODEL needs; this holds what the PAPER SAYS, so a
comparison can never quietly drift away from the published values. Same role as
`katrin_wall/testdata.py` and `vk3_wall/testdata.py`.

One difference from those two, and it is the important one: for WSH3 and VK3 the primary source was
in hand and the numbers are the TEST's own. Here the primary source (Aldemir, Binici & Canbay 2017,
ACI Struct. J. 114(2): 395-406) is NOT in the repo, so every "measured" value below is as REPORTED
BY the 2019 lattice paper. Anything the 2019 paper does not print — the loading protocol, the axial
load history, the bar type, the number of cycles, the concrete cast sequence, the foundation
dimensions — is simply absent, and is listed under `NOT_REPORTED` rather than guessed.

Run this module to print the table.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Value:
    """One published quantity: what it is, its value, and exactly where it came from."""

    name: str
    value: float | str
    unit: str
    source: str
    note: str = ""


# --- measured (as reported in the 2019 paper) -----------------------------------------------------
MEASURED = (
    Value("initial stiffness", 1038.44, "kN/mm", "Table 4",
          "by far the stiffest of the six specimens; the next is the Foster & Gilbert deep beam "
          "at 553.9"),
    Value("maximum force", 963.592, "kN", "Table 4", ""),
    Value("maximum displacement", 20.0, "mm", "Table 2",
          "the value the simulation was run to, quoted as the specimen's"),
    Value("drift reached", "~1%", "-", "text, 'RC Wall Simulations'",
          "'tested up to an approximate 1% drift ratio' — 20/2250 = 0.89%, 20/1500 = 1.33%, which "
          "is the evidence for the Fig. 10(a) height of 2250 over the 2680 the Table 2 counts imply"),
    Value("strength degradation", "none", "-", "text, 'RC Wall Simulations'",
          "'did not sustain any strength degradation despite severe inclined cracking' — so the "
          "test is a LOWER BOUND on drift capacity, as WSH3's was"),
    Value("failure mode", "flexural-shear", "-", "text, 'RC Wall Simulations'",
          "'designed to yield in shear'"),
)

# --- lattice inputs the paper used (Table 1) ------------------------------------------------------
INPUTS = (
    Value("f_t", 1.85, "MPa", "Table 1, footnote a", "TS 500 (TSE 2000) — a code value, not a test"),
    Value("f_c", 28.0, "MPa", "Table 1", ""),
    Value("E_t", 24.87e3, "MPa", "Table 1, footnote b",
          "ACI 318's 4700*sqrt(28) = 24,869.7, so the table's 24.87 GPa is self-consistent"),
    Value("G_f", 75.0, "N/m", "Table 1", ""),
    Value("f_y", 360.0, "MPa", "Table 1", "the lowest of the six specimens (others 414-460)"),
    Value("E_s", 200e3, "MPa", "Table 1", ""),
    Value("N (axial load)", "none", "kN", "Table 1", "blank, against 632/600/378 kN elsewhere"),
    Value("grid size d", 20.0, "mm", "Tables 1 and 2",
          "Fig. 10(b)'s legend says 25 mm for both curves — an unreconciled third value"),
    Value("a1", 1.5, "-", "Table 1", "footnote d ties it to the horizon: a1 = 3.01 at delta = 3.01d"),
    Value("a2", 70.0, "-", "Table 1", ""),
    Value("a3", 360.0, "-", "Table 1", "a3/a2 = 5.14; the repo's aydin_lattice_softening default is 5.0"),
    Value("K_p", 1.0e9, "-", "Table 1", "PID gains; this repo drives displacement directly instead"),
    Value("K_i", 1.0e8, "-", "Table 1", ""),
    Value("K_d", 1.0e4, "-", "Table 1", ""),
)

# --- what their lattice produced (Tables 2 and 4) -------------------------------------------------
@dataclass(frozen=True)
class LatticeRun:
    horizon: float
    stiffness: float      # kN/mm
    force: float          # kN
    elements: int
    seconds_per_10k: float


PAPER_RUNS = (
    LatticeRun(1.5, 943.16, 1164.413, 80_684, 66.88),
    LatticeRun(3.01, 1052.37, 1325.675, 280_260, 208.09),
)
PAPER_PARTICLES = 20_385           # Table 2, both horizons share the node set
PAPER_TIMESTEP = 5.0e-8            # s, Table 2

# --- absent from the 2019 paper -------------------------------------------------------------------
NOT_REPORTED = (
    "loading protocol (number, amplitude and order of cycles)",
    "loading rate for the wall specimens (Table 3 gives rates only for the two sensitivity cases)",
    "bar type — deformed or plain — and therefore whether perfect bond is defensible",
    "what 'double wall' means for this specimen's section: the Fig. 10(a) plan shows a 50|140|100|60 "
    "= 350 mm stack whose parts are never named, with the wall itself dimensioned 120",
    "foundation depth and out-of-plane thickness (the plan gives only its 2950 mm length)",
    "concrete cast sequence / whether the wall and foundation share a grade",
    "measured hysteresis as an ORDERED path — Fig. 10(b) draws the experiment as discrete dots, so "
    "which dot follows which is gone and per-cycle energy and degradation are not recoverable (the "
    "same limit as SW-NC-FF's Fig. 14b). NOTE what IS recoverable, and is, by `digitize.py` (D77): "
    "the cloud's OUTER ENVELOPE — extremes survive unordered — and both of the paper's own lattice "
    "curves, which are clean coloured lines. The envelope is an UPPER BOUND on the backbone, since "
    "a point on it may be a loop tip or the unloading side of a larger cycle",
    "the loading protocol behind that cloud, so the envelope cannot be anchored to known amplitudes "
    "the way WSH3's backbone was (D65)",
    "axial load: reported as absent, but no load cell record is shown",
)

# --- discrepancies inside the 2019 paper itself ----------------------------------------------------
DISCREPANCIES = (
    ("panel size", "Fig. 10(a) draws 3000 x 2250; inverting Table 2's counts gives 150 x 134 cells",
     "SETTLED at 3000 x 2250 (user, 2026-08-28): Fig. 10(a) is to scale and reconstructs from its "
     "own pixels to 0.17%. Fig. 4(f) is NOT a third geometry — it is a DETAIL VIEW of the bottom "
     "1500 mm of his model (author, 2026-09-05), so its '1500' is a crop height"),
    ("his analysed panel", "the Table 2 inversion is unique only UP TO TRANSPOSITION: 3000 x 2680 "
     "and 2680 x 3000 both reproduce his 20,385 nodes and 80,684 struts exactly",
     "OPEN. `replica/` takes 3000 x 2680 (squat). The 2680 dimensioned along the HORIZONTAL in "
     "Fig. 4(f) points at 2680 x 3000 (slender), since a full-width crop shows the true width. "
     "Resolving it would move every published replica ratio; it needs the author"),
    ("initial stiffness", "the measured 1,038 kN/mm exceeds the 804 kN/mm an UNCRACKED 3000x2250x120 "
     "transformed section can give — impossible for a fixed-base cantilever",
     "OPEN, and no longer a question about panel size: candidates are the thickness (two 120 mm "
     "shells -> 1,608), the top boundary (rotation-restrained -> 1,156), or Table 4's definition"),
    ("grid size", "Tables 1 and 2 say d = 20 mm; the Fig. 10(b) legend says 25 mm",
     "neither is used: 50 mm divides every dimension of this specimen and 20 does not (2250/20 = 112.5)"),
    ("thickness", "120 mm is dimensioned in the plan, but the model thickness is never stated",
     "120 adopted; 240 or 350 would put the computed initial stiffness 2.5-3.6x above the measured"),
)


def _rows(values):
    for v in values:
        val = f"{v.value:,.4g}" if isinstance(v.value, (int, float)) else v.value
        yield f"  {v.name:<22s} {val:>12s} {v.unit:<6s}  {v.source:<28s} {v.note}"


def main() -> None:
    print(f"Aldemir et al. (2017) squat wall — as reported by Aydin, Tuncay & Binici (2019)\n")
    print("MEASURED (via the 2019 paper; the 2017 ACI paper is not in the repo)")
    print("\n".join(_rows(MEASURED)))
    print("\nLATTICE INPUTS the paper used")
    print("\n".join(_rows(INPUTS)))
    print(f"\nWHAT THEIR LATTICE PRODUCED  ({PAPER_PARTICLES:,} nodes, dt = {PAPER_TIMESTEP:g} s)")
    print(f"  {'horizon':<10s}{'K (kN/mm)':>12s}{'F (kN)':>12s}{'K/K_exp':>10s}{'F/F_exp':>10s}"
          f"{'elements':>12s}{'s/10k steps':>13s}")
    k_exp, f_exp = MEASURED[0].value, MEASURED[1].value
    for r in PAPER_RUNS:
        print(f"  {r.horizon:<10.2f}{r.stiffness:>12,.2f}{r.force:>12,.2f}"
              f"{r.stiffness / k_exp:>10.3f}{r.force / f_exp:>10.3f}"
              f"{r.elements:>12,d}{r.seconds_per_10k:>13.2f}")
    print("\n  Their 1.5d run is the WORST of the six specimens on peak force (+20.8%); the six "
          "range 0.93-1.38.")
    print("\nINTERNAL DISCREPANCIES IN THE 2019 PAPER")
    for what, says, taken in DISCREPANCIES:
        print(f"  {what}:\n      {says}\n      -> {taken}")
    print("\nNOT REPORTED (do not invent these)")
    for item in NOT_REPORTED:
        print(f"  - {item}")


if __name__ == "__main__":
    main()
