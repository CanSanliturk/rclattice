"""The parameter registry — one record per parameter (PLAN.md §9).

WHY A REGISTRY AND NOT AN ARGPARSE BLOCK. The author's note that parameters and their values will
grow is a design constraint, not a caveat: adding a parameter here is ONE record, and the CLI, the
run-directory name, `params.json` and the master report all derive from it. Nothing downstream
parses directory names — names are for humans, `params.json` is for the machine — so a new axis or a
new value extends the matrix without touching the reporting.

`affects` is what lets the master report know which runs are comparable: two runs differing only in
a `report` parameter answer the same question, and a new reporting-only parameter therefore never
fragments the matrix.
"""
from __future__ import annotations

from typing import Any

from rclattice.study.registry import Param, Registry

# LIFTED (D103): `Param` and every derived function now live in `rclattice.study.registry`; this
# module keeps the Aldemir records VERBATIM and re-exports the module-level API the rest of the
# study (and `doc/reports/aydin_aldemir_runs/generate.py`) was written against.

# Bumped whenever a parameter is added, removed or changes meaning. Runs made before a parameter
# existed stay readable: the master report fills its default and says so.
# v2 (2026-09-05): added `fc`, `ft` and `cycles`.
# v3 (2026-09-07): added `steel_rupture` and `concrete_residual`. MISSED AT THE TIME — two runs went
# out stamped v2 while carrying v3 parameters, so a reader could not tell from the stamp whether a
# run predated the failure switches or simply left them off. Bump the version IN THE SAME EDIT that
# adds a parameter.
# Runs written under an older version stay readable: the master report fills a missing parameter
# with today's default and says which it filled.
# v5 (2026-09-10): added `steel_b`. Hardening was hardcoded at 0.01, an unprinted convention that
# raises a yielded bar 27% by eps_su = 0.05 and is therefore a candidate for the flat cyclic
# envelope of D99. Default reproduces every earlier run exactly.
SCHEMA_VERSION = 5    # 4: comp/bond values renamed (D94); older records are read through LEGACY_VALUES, not migrated


REGISTRY: tuple[Param, ...] = (
    # --- the stem: what the run IS -----------------------------------------------------------
    Param("analysis", "", "pushover", "what to run", "analysis", stem=True,
          choices=("elastic", "static", "pushover", "cyclic")),
    Param("comp", "", "linear", "what compression does: 'linear' never yields (the paper's own "
          "law), 'capped' yields at fc and holds (elastic-perfectly-plastic), 'crushing' yields "
          "then softens and crushes (Concrete02). PLAN §3", "model", stem=True,
          choices=("linear", "capped", "crushing")),
    Param("tail", "", "solved", "tension tail: Gf-solved a2/a3, or the paper's printed pair",
          "model", stem=True, choices=("solved", "paper")),
    Param("bond", "", "perfect",
          "bond model (PLAN §4): 'perfect' is perfect bond on shared nodes — the bar and the "
          "concrete cannot slip; 'bondNN' gives each bar its OWN nodes tied to the concrete by a "
          "full-area ring, whose residual plateau after the brittle drop is NN% of f_t "
          "(the paper's own value is 70)",
          "model", stem=True, choices=("perfect", "bond60", "bond70")),
    Param("bond_damage", "bonddmg", False,
          "make the bond law's brittle drop IRREVERSIBLE (D95). The published law is 'elastic "
          "brittle', and brittle means a broken bond stays broken; the default ElasticMultiLinear "
          "form heals on every reversal, which is invisible in a push and wrong in a cycle. "
          "Monotonic envelopes are identical to four decimals, so this only changes cyclic work",
          "model", flag=True),      # an EXTRA, not a stem: the stem is analysis/comp/tail/bond

    # --- geometry and discretization ---------------------------------------------------------
    Param("panel", "panel", "fig10a", "which panel geometry (PLAN §2)", "model",
          choices=("fig10a", "table2", "table2t")),
    Param("tw", "tw", 210.0, "wall thickness in mm (D73)", "model", type=float),
    Param("mesh", "m", 50.0, "lattice spacing in mm", "model", type=float),
    Param("horizon", "h", 1.5, "strut horizon in multiples of the mesh", "model", type=float),

    # --- materials ---------------------------------------------------------------------------
    Param("fc", "fc", 28.0, "concrete compressive strength in MPa (Table 1)", "model", type=float),
    Param("ft", "ft", 1.85, "concrete tensile strength in MPa (Table 1, TS 500)", "model",
          type=float),
    Param("gf", "gf", 0.075, "fracture energy in N/mm (NOT a neutral knob: D75)", "model",
          type=float),
    Param("fcx", "fcx", 1.0, "strut compressive-strength scale, `eppcomp` only (PLAN §5)",
          "model", type=float),

    # --- failure mechanisms (2026-09-07): without these NOTHING in the model can fail, so a
    # monotonic push rises forever and a cyclic envelope sags a few percent and holds (D88/D89).
    Param("steel_rupture", "eps", 0.0,
          "bar rupture strain; 0 = bars never break. NOT given by the paper — an assumption that "
          "decides the predicted drift capacity (mesh ~0.025, hot-rolled bar ~0.075)",
          "model", type=float),
    Param("concrete_residual", "res", 0.2,
          "crushing strength floor as a fraction of fc; D22's 0.2 was for Newton solvers, 0 lets a "
          "strut crush to nothing under the explicit march", "model", type=float),
    Param("steel_b", "sb", 0.01,
          "steel strain-hardening ratio; 0 = elastic-perfectly-plastic ties. NOT given by the paper "
          "(Table 1 prints f_y only), so 0.01 is a convention like eps_su, and it props up a flat "
          "cyclic envelope: a bar at eps_su = 0.05 carries 1.27x f_y", "model", type=float),

    # --- reinforcement -----------------------------------------------------------------------
    Param("rebar_top", "rebartop", True, "run the longitudinal bars to the top face (D78/D84)",
          "model", flag=True),
    Param("rebar", "", True, "include reinforcement at all", "model", flag=True),

    # --- what the analysis does --------------------------------------------------------------
    # ALWAYS IN THE NAME, default or not (user instruction, 2026-09-06). A pushover stops where its
    # target says, so a run whose name omits the target is unreadable later: five cells ending at
    # 6.75 mm looked like a property of the wall when it was only the number they were given.
    Param("drift", "d", 0.003, "target drift for a pushover, e.g. 0.012 for 1.2%", "analysis",
          type=float, always=True),
    Param("proto", "proto", "none",
          "cyclic protocol: a named preset, 'ladder' (the eight-level shape scaled to --drift), "
          "or an explicit comma-separated list of drift amplitudes", "analysis"),
    Param("cycles", "c", 1, "cycles per level; PLAN §6 argues ONE is the defensible shape, since "
                            "the digitized cloud shows ~4-6 loops and not 16", "analysis",
          type=int),
    Param("rate", "r", 7.6, "drive speed in mm/s", "analysis", type=float),
    Param("damping", "z", 0.5, "damping ratio for dynamic relaxation", "analysis", type=float),
    Param("integrator", "", "explicit", "explicit CentralDifference or implicit Newmark",
          "analysis", choices=("explicit", "newmark")),
    Param("steps_per_period", "", 0,
          "explicit sub-steps per fundamental period; 0 = SIZE IT FROM critical_time_step "
          "(D74 — an explicit step sized off T1 is 15-20x too large and diverges)",
          "analysis", type=int),
    Param("steps", "", 400, "displacement-control steps (static solver only)", "analysis",
          type=int),

    # --- reporting only: never fragments the matrix -------------------------------------------
    Param("groups", "", True, "decompose the base shear by strut group", "report", flag=True),
    Param("capture", "", True, "save displacement fields for the damage figure", "report",
          flag=True),
    Param("progress_every", "", 2000, "progress print interval in steps", "report", type=int),
)

# The names these axes USED to carry. `c02` was named after an OpenSees material while its
# siblings were named after the physics, and `nobond` said the opposite of what it meant — it is
# perfect bond, the stiffest possible, not the absence of bond. Renamed 2026-09-08 (D94).
#
# Every run recorded before that carries the old value in its `params.json`, and the shell scripts
# in this directory still pass the old flags, so both are accepted forever: `normalize` maps a
# record on the way in, and the CLI maps a flag on the way in. Nothing has to be migrated, and the
# 26 archived runs keep grouping with the new ones.
LEGACY_VALUES = {
    "comp": {"lincomp": "linear", "eppcomp": "capped", "c02": "crushing"},
    "bond": {"nobond": "perfect", "bond-a06": "bond60", "bond-a07": "bond70"},
}



REG = Registry(REGISTRY, SCHEMA_VERSION, LEGACY_VALUES,
               stem_format="{analysis}_{comp}-{tail}_{bond}")
BY_NAME = REG.by_name
STEM = REG.stem

# --- the module-level API, delegated -------------------------------------------------------------
normalize = REG.normalize
add_arguments = REG.add_arguments
resolve = REG.resolve
run_name = REG.run_name
comparable_key = REG.comparable_key
describe = REG.describe
