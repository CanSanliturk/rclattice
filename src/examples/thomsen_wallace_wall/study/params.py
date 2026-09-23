"""The RW2 parameter registry — composed from the shared factories plus this specimen's own axes.

The shared records (`rclattice.study.registry`) carry the names and codes every study uses, so the
master matrix, run sheets and advisor page read RW2 runs the way they read Aldemir's. What is
specimen-specific: the DEFAULTS (mesh 25 on the graded grid, 2.5% drift target, damping 0.5) and
the axes below.

v1 (2026-09-20): first registry — `comp` / `tail` / `bond` as for Aldemir (D94 names from the
start, no legacy spellings), plus `grid` (D104), `field` (the calibration route, measured in
Stage 0) and `nu`.
v2 (2026-09-22): `top_band` / `top_band_pitch` — the Fig. 9(a) dense top load-introduction band
(D108), OFF by default so no v1 run changes.
v3 (2026-09-23): `steel_iso` — Steel02 isotropic hardening (D109), OFF by default so the bar stays
purely kinematic and no v1/v2 run changes. Only bites under reversals, i.e. cyclic runs.
"""
from __future__ import annotations

from rclattice.study.registry import (Param, Registry, analysis_param, analysis_params,
                                      concrete_params, discretisation_params, failure_params,
                                      rebar_params, report_params)

SCHEMA_VERSION = 3

REGISTRY: tuple[Param, ...] = (
    analysis_param(),
    Param("comp", "", "crushing", "what compression does: 'linear' never yields (Aydin's own law), "
          "'capped' yields at fc and holds (elastic-perfectly-plastic), 'crushing' yields then "
          "softens and crushes (Concrete02). Default is the repo's own law, which on Aldemir gave "
          "0.997x (D97)", "model", stem=True, choices=("linear", "capped", "crushing")),
    Param("tail", "", "solved", "tension tail: Gf-solved a2/a3, or the paper's printed pair "
          "(80/350, fitted at his 19 mm grid — NOT neutral at another mesh, D87)",
          "model", stem=True, choices=("solved", "paper")),
    Param("bond", "", "perfect",
          "bond model: 'perfect' shares nodes; 'bondNN' gives each bar its own nodes tied by a "
          "full-area ring whose residual after the brittle drop is NN% of f_t (D72/D96)",
          "model", stem=True, choices=("perfect", "bond60", "bond70")),
    Param("bond_damage", "bonddmg", False, "make the bond law's brittle drop IRREVERSIBLE (D95)",
          "model", flag=True),

    # --- discretisation ---------------------------------------------------------------------
    Param("grid", "g", "rebar", "node grid: 'rebar' = graded lines on the bar axes (D104), "
          "'uniform' = the structured grid (mesh must divide 1220 and 3660)", "model",
          choices=("rebar", "uniform")),
    *discretisation_params(mesh=25.0, horizon=1.5),

    # --- calibration --------------------------------------------------------------------------
    Param("field", "", "uniaxial", "energy-balance affine field: the thesis 'uniaxial' route "
          "(D47) or the 2019 paper's 'equibiaxial' (D72); Stage 0 measured 0.951 vs 0.821 "
          "against the same-grid continuum", "model", choices=("uniaxial", "equibiaxial")),
    Param("nu", "nu", 0.20, "Poisson ratio in the balance (a repo convention; not printed)",
          "model", type=float),

    # --- materials ------------------------------------------------------------------------------
    *concrete_params(fc=42.8, ft=2.03, gf=0.075, fc_source="Table 1",
                     ft_source="Table 1, no footnote; = 0.31 sqrt(fc)"),
    Param("fcx", "fcx", 1.0, "strut compressive-strength scale, `capped` only", "model", type=float),
    *failure_params(steel_b=0.01),
    Param("fy", "fy", 414.0, "steel yield in MPa (Table 1: the NOMINAL 60 ksi; coupons unprinted)",
          "model", type=float),
    Param("steel_iso", "iso", 0.0, "Steel02 ISOTROPIC hardening (D109): growth of BOTH yield "
          "envelopes as a fraction of fy per unit accumulated plastic strain (sets a1 = a3, with "
          "a2 = a4 = 1). 0 = off, the default, leaving the bar purely KINEMATIC. Only observable "
          "under REVERSALS — it changes nothing in a monotonic pushover. Typical RC practice is "
          "0.01-0.03; UNPRINTED for this specimen, so it is an assumption like steel_b",
          "model", type=float),

    # --- reinforcement ------------------------------------------------------------------------
    *rebar_params(),
    Param("top_band", "tb", 0.0, "dense top load-introduction band: height in mm of denser "
          "full-width horizontal web bars at the top (Fig. 9a, D108); 0 = off, the default. "
          "The drawing suggests ~1000 mm — INFERRED, not dimensioned", "model", type=float),
    Param("top_band_pitch", "tbp", 76.0, "spacing in mm of the dense top-band bars (default = the "
          "76 mm hoop pitch); only meaningful when --top-band > 0", "model", type=float),

    # --- the analysis -------------------------------------------------------------------------
    *analysis_params(drift=0.025, rate=7.6, damping=0.5, cycles=1),

    # --- reporting only -----------------------------------------------------------------------
    *report_params(),
)

REG = Registry(REGISTRY, SCHEMA_VERSION, legacy_values={},
               stem_format="{analysis}_{comp}-{tail}_{bond}")
BY_NAME = REG.by_name
STEM = REG.stem
normalize = REG.normalize
add_arguments = REG.add_arguments
resolve = REG.resolve
run_name = REG.run_name
comparable_key = REG.comparable_key
describe = REG.describe
