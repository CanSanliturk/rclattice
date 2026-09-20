"""The per-run report (PLAN.md §10): generated into the run directory, readable on its own.

SECTIONS ARE A LIST OF BUILDERS, so a new figure or a new comparison is one entry rather than an
edit inside a monolithic function. Each builder takes `(params, data)` and returns markdown or None
(None = nothing to say for this run, and the section is dropped).

The report states what a run does and does NOT license. That is not decoration: this study's fair
claims and unfair ones are separated by the specimen's own gaps (an invented protocol, an inferred
thickness, a second-hand test record), and a number quoted without them travels further than it
should.
"""
from __future__ import annotations

from pathlib import Path

# LIFTED (D103): the section builders live in `rclattice.study.report`; this module keeps the
# Aldemir-specific text — where each parameter's value comes from, and what a run does and does
# not license — and re-exports `write` for `rescore.py`-era callers.

# Where each parameter's value comes from. "measured" = from the test or the paper's tables;
# "convention" = a code/standard formula; "inferred" = reconstructed from a figure or a count;
# "assumed" = a modelling choice this study made.
SOURCES = {
    "panel": ("inferred", "Fig. 10(a), a dimensioned drawing that reconstructs to 0.17% (D72)"),
    "tw": ("inferred", "210 = 50+100+60 of the Fig. 10(a) plan stack; 120 dimensions ONE panel of "
                       "a precast double wall (D73)"),
    "mesh": ("assumed", "50 divides 3000, 2250, 100 and 50 exactly; the paper's own 20 does not "
                        "divide 2250"),
    "horizon": ("measured", "his own recommendation for RC (Table 4 runs 1.5 and 3.01)"),
    "gf": ("convention", "MC90-style value; NOT a neutral knob — x2 moved base shear +14.2% (D75)"),
    "fcx": ("assumed", "strut-strength scale; direction-dependent, see PLAN §5"),
    "comp": ("assumed", "which compression law is being tested"),
    "tail": ("measured", "'paper' = Table 1's printed a2/a3; 'solved' = the same Gf, re-solved at "
                         "this mesh"),
    "bond": ("assumed", "PLAN §4"),
    "rebar_top": ("assumed", "the paper says nothing about bar termination (D78/D84)"),
    "rate": ("assumed", "the repo's cross-study drive speed; the licence is the residual (D64)"),
    "damping": ("assumed", "damping is the knob that matters, not rate (D64)"),
    "proto": ("assumed", "THE PROTOCOL IS INVENTED — the paper never prints one (PLAN §6)"),
    "cycles": ("assumed", "PLAN §6 argues ONE cycle per level: the digitized cloud shows "
                          "2-4 dot bands between 6 and 14 mm, i.e. ~4-6 loops, not 16"),
    "fc": ("measured", "Table 1"),
    "ft": ("convention", "TS 500, Table 1; epsc0 is re-derived as 2fc/E whenever fc moves (D56)"),
    "analysis": ("convention", "which runner is used"),
    "rebar": ("measured", "Ø8 @ 100 both ways, 3 bars/position, rho = 0.718% (Fig. 10(a))"),
    "drift": ("assumed", "how far this run was asked to go — NOT a capacity, and a run that ends "
                         "at its own target licenses nothing about behaviour past it (D78)"),
    "integrator": ("convention", "explicit CentralDifference: 45x cheaper per step, and implicit "
                                 "and explicit agree on peak to 1.0013 here (D74/D75)"),
    "steps_per_period": ("convention", "0 = sized from critical_time_step; sizing off T1 is "
                                       "15-20x too large and diverges (D74)"),
    "steps": ("convention", "displacement-control steps, static solver only"),
    "groups": ("convention", "reporting only"),
    "capture": ("convention", "reporting only"),
    "progress_every": ("convention", "reporting only"),
}

MEASURED_F_KN = 963.592      # Table 4
MEASURED_K_KNMM = 1038.44    # Table 4
AYDIN = {1.5: (943.16, 1164.413), 3.01: (1052.37, 1325.675)}   # horizon -> (K kN/mm, F kN)


def licence(params: dict) -> tuple[list, list]:
    """(fair, unfair) for this run — the specimen's own gaps, stated on every report."""
    fair = ["peak strength", "elastic stiffness", "damage location",
            "the load-path split between plane-sections and truss action",
            "the controlled effect of this cell's own axis, since neighbouring cells differ in "
            "exactly one parameter"]
    unfair = [("drift capacity", "the loading protocol is invented — the 2019 paper never prints "
                                 "one and the 2017 test paper is not in the repo (PLAN §6)"),
              ("cyclic energy and per-cycle degradation", "the digitized test record is an "
                                                          "unordered point cloud"),
              ("anything resting on the section", "the 210 mm thickness is inferred (D73)")]
    if params["comp"] == "capped":
        unfair.append(("a compression-driven collapse",
                       "the EPP branch cannot lose compressive load-carrying capacity by material "
                       "failure — that is a property of the law, not a finding. A load-path "
                       "collapse can still occur, by cracking away the restraint (D55)"))
        fair.append("the compressive demand, since struts now yield at a stated fc rather than "
                    "carrying stress without limit")
    if params["comp"] == "linear":
        unfair.append(("anything about compressive capacity",
                       "compression is LINEAR AT E FOREVER in his published law: a strut carries "
                       "unbounded compressive stress, and the specimen's strength is meant to "
                       "emerge as indirect tensile splitting instead"))
    if params["analysis"] in ("pushover", "static"):
        unfair.append(("loop shape, pinching, residual drift",
                       "this is a monotonic run, and both concrete laws are path-independent "
                       "anyway (PLAN §3)"))
    return fair, unfair




def write(out: Path, params: dict, data: dict) -> Path:
    from rclattice.study import report as shared
    from study_spec import SPEC
    return shared.write(SPEC, out, params, data)
