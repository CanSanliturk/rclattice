"""The parameter registry — one record per parameter (Aldemir PLAN.md §9, D83, lifted in D103).

WHY A REGISTRY AND NOT AN ARGPARSE BLOCK. Parameters and their values grow: adding a parameter is
ONE record, and the CLI, the run-directory name, `params.json` and the master report all derive
from it. Nothing downstream parses directory names — names are for humans, `params.json` is for the
machine — so a new axis or a new value extends the matrix without touching the reporting.

`affects` is what lets the master report know which runs are comparable: two runs differing only in
a `report` parameter answer the same question, and a new reporting-only parameter therefore never
fragments the matrix.

Rules learned the hard way, and enforced by convention here:

* **Bump `schema_version` in the same edit that adds a parameter.** Two Aldemir runs once went out
  stamped v2 while carrying v3 parameters, so a reader could not tell from the stamp whether a run
  predated the failure switches or simply left them off.
* Runs predating a parameter stay readable: the master report fills today's default and SAYS SO
  (`filled_defaults`).
* Renamed values are mapped through `legacy_values` + `normalize()`, never migrated on disk (D94).
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Callable


@dataclass(frozen=True)
class Param:
    name: str                      # long name; the CLI flag is --<name with _ as ->
    code: str                      # short code used in run-directory names ("" = never in a name)
    default: Any
    help: str
    affects: str                   # "model" | "analysis" | "report"
    type: Callable = str
    choices: tuple | None = None
    stem: bool = False             # part of the fixed <analysis>_<axis>-<axis>_<axis> stem
    flag: bool = False             # a boolean switch rather than a value
    always: bool = False           # appears in the run name even at its default value

    @property
    def cli(self) -> str:
        return "--" + self.name.replace("_", "-")


@dataclass
class Registry:
    """A specimen's parameter set, with everything derived from it.

    `stem_format` renders the fixed part of a run name from the stem parameters, e.g. for the
    Aldemir study `"{analysis}_{comp}-{tail}_{bond}"`. It is a format string over the stem names so
    the specimen decides how its axes read; the default joins them with underscores.
    """

    params: tuple[Param, ...]
    schema_version: int
    legacy_values: dict[str, dict[str, str]] = field(default_factory=dict)
    stem_format: str | None = None

    def __post_init__(self) -> None:
        self.by_name = {p.name: p for p in self.params}
        self.stem = tuple(p.name for p in self.params if p.stem)
        if self.stem_format is None:
            self.stem_format = "_".join("{" + n + "}" for n in self.stem)

    def __iter__(self):
        return iter(self.params)

    def __getitem__(self, name: str) -> Param:
        return self.by_name[name]

    def normalize(self, params: dict) -> dict:
        """Map any legacy axis values in `params` to their current names, in place."""
        for name, mapping in self.legacy_values.items():
            v = params.get(name)
            if v in mapping:
                params[name] = mapping[v]
        return params

    def _value(self, name: str):
        """argparse `type` that accepts a legacy value and hands back the current one."""
        mapping = self.legacy_values.get(name, {})
        return lambda v: mapping.get(v, v)

    def add_arguments(self, parser: argparse.ArgumentParser) -> None:
        """Derive the whole CLI from the registry.

        Help text is %-ESCAPED on the way into argparse, which expands `%(default)s`-style tokens
        and dies on a bare `%` — `--drift`'s "e.g. 0.012 for 1.2%" was enough to make `--help`
        raise while every actual run worked. The registry keeps the unescaped text, since the
        reports render it.
        """
        def h(text: str) -> str:
            return text.replace("%", "%%")

        for p in self.params:
            if p.flag:
                # A boolean gets both switches, so a True default is still overridable.
                parser.add_argument(p.cli, dest=p.name, action="store_true", default=None,
                                    help=h(p.help) + f" (default {p.default})")
                parser.add_argument("--no-" + p.name.replace("_", "-"), dest=p.name,
                                    action="store_false", help=argparse.SUPPRESS)
            else:
                parser.add_argument(p.cli, dest=p.name,
                                    type=self._value(p.name) if p.name in self.legacy_values
                                    else p.type,
                                    choices=p.choices, default=None,
                                    help=h(p.help) + f" (default {p.default})")

    def resolve(self, args: argparse.Namespace) -> dict:
        """Every parameter, defaults included — the record `params.json` carries."""
        out = {"schema_version": self.schema_version}
        for p in self.params:
            v = getattr(args, p.name, None)
            out[p.name] = p.default if v is None else v
        return self.normalize(out)

    def fill_defaults(self, params: dict) -> list[str]:
        """Fill parameters that postdate a record with today's defaults; return what was filled.

        A filled value is an ASSUMPTION about the run rather than a record of it, which is why the
        names come back: every display of such a run says which were filled.
        """
        missing = [p.name for p in self.params if p.name not in params]
        for name in missing:
            params[name] = self.by_name[name].default
        self.normalize(params)
        return missing

    def run_name(self, params: dict, *, stamp: str | None = None) -> str:
        """`<stamp>_<stem>[_extras]`, extras only when non-default (or `always`).

        For humans. Nothing reads it back — the master report builds the matrix from `params.json`.
        """
        stamp = stamp or f"{datetime.now():%Y-%m-%d_%H%M%S}"
        name = f"{stamp}_" + self.stem_format.format(**{n: params[n] for n in self.stem})
        for p in self.params:
            if p.stem or not p.code:
                continue
            v = params[p.name]
            if v == p.default and not p.always:
                continue
            if p.name == "drift":
                # As a PERCENT, which is how the target is asked for and read: _d0.30pct.
                name += f"_{p.code}{float(v) * 100:g}pct"
            elif p.flag:
                name += f"_{'' if v else 'no'}{p.code}"
            else:
                name += f"_{p.code}{v:g}" if isinstance(v, float) else f"_{p.code}{v}"
        return name

    def comparable_key(self, params: dict, *, axis: str) -> tuple:
        """What must match for two runs to differ ONLY in `axis`.

        Report-only parameters are excluded, so adding a figure never splits a comparison in two.
        """
        return tuple((p.name, params.get(p.name, p.default))
                     for p in self.params
                     if p.affects != "report" and p.name != axis)

    def describe(self, params: dict) -> list[tuple[str, Any, Any, str]]:
        """(name, value, default, affects) rows, for the per-run report's parameter table."""
        return [(p.name, params.get(p.name, p.default), p.default, p.affects)
                for p in self.params]


# --- parameters every specimen shares ------------------------------------------------------------
#
# Factories rather than constants, so a specimen can set the default that suits it (a drift target,
# a mesh, a damping ratio) while the NAME, CODE and MEANING stay identical across studies — which
# is what lets the shared reporting read every study the same way. The Aldemir registry predates
# these and keeps its own literal records; a new study composes its registry from these plus its
# own axes.


def analysis_param(choices=("elastic", "static", "pushover", "cyclic")) -> Param:
    return Param("analysis", "", "pushover", "what to run", "analysis", stem=True, choices=choices)


def discretisation_params(*, mesh: float, horizon: float = 1.5) -> tuple[Param, ...]:
    return (
        Param("mesh", "m", float(mesh), "lattice spacing in mm", "model", type=float),
        Param("horizon", "h", float(horizon), "strut horizon in multiples of the mesh", "model",
              type=float),
    )


def concrete_params(*, fc: float, ft: float, gf: float, fc_source: str = "",
                    ft_source: str = "") -> tuple[Param, ...]:
    return (
        Param("fc", "fc", float(fc), "concrete compressive strength in MPa"
              + (f" ({fc_source})" if fc_source else ""), "model", type=float),
        Param("ft", "ft", float(ft), "concrete tensile strength in MPa"
              + (f" ({ft_source})" if ft_source else ""), "model", type=float),
        Param("gf", "gf", float(gf), "fracture energy in N/mm (NOT a neutral knob: D75)", "model",
              type=float),
    )


def failure_params(*, steel_b: float = 0.01) -> tuple[Param, ...]:
    """The failure switches: without these NOTHING in the model can fail (D88/D89/D91)."""
    return (
        Param("steel_rupture", "eps", 0.0,
              "bar rupture strain; 0 = bars never break. Decides the predicted drift capacity "
              "(mesh ~0.025, hot-rolled bar ~0.075) — D97/D102", "model", type=float),
        Param("concrete_residual", "res", 0.2,
              "crushing strength floor as a fraction of fc; D22's 0.2 was for Newton solvers, 0 "
              "lets a strut crush to nothing under the explicit march", "model", type=float),
        Param("steel_b", "sb", float(steel_b),
              "steel strain-hardening ratio; 0 = elastic-perfectly-plastic ties. The largest "
              "single control on drift capacity (D101/D102): a bar at eps_su = 0.05 carries "
              "1.27x f_y at b = 0.01", "model", type=float),
    )


def rebar_params() -> tuple[Param, ...]:
    return (
        Param("rebar_top", "rebartop", True, "run the longitudinal bars to the top face (D78/D84)",
              "model", flag=True),
        Param("rebar", "", True, "include reinforcement at all", "model", flag=True),
    )


def analysis_params(*, drift: float, rate: float = 7.6, damping: float = 0.5,
                    cycles: int = 1) -> tuple[Param, ...]:
    return (
        # ALWAYS IN THE NAME, default or not (user instruction, 2026-09-06). A pushover stops where
        # its target says, so a run whose name omits the target is unreadable later.
        Param("drift", "d", float(drift), "target drift for a pushover, e.g. 0.012 for 1.2%",
              "analysis", type=float, always=True),
        Param("proto", "proto", "none",
              "cyclic protocol: a named preset, 'ladder' (the eight-level shape scaled to "
              "--drift), or an explicit comma-separated list of drift amplitudes", "analysis"),
        Param("cycles", "c", int(cycles), "cycles per level", "analysis", type=int),
        Param("rate", "r", float(rate), "drive speed in mm/s", "analysis", type=float),
        Param("damping", "z", float(damping), "damping ratio for dynamic relaxation", "analysis",
              type=float),
        Param("integrator", "", "explicit", "explicit CentralDifference or implicit Newmark",
              "analysis", choices=("explicit", "newmark")),
        Param("steps_per_period", "", 0,
              "explicit sub-steps per fundamental period; 0 = SIZE IT FROM critical_time_step "
              "(D74 — an explicit step sized off T1 is 15-20x too large and diverges)",
              "analysis", type=int),
        Param("steps", "", 400, "displacement-control steps (static solver only)", "analysis",
              type=int),
    )


def report_params() -> tuple[Param, ...]:
    """Reporting only: never fragments the matrix."""
    return (
        Param("groups", "", True, "decompose the base shear by strut group", "report", flag=True),
        Param("capture", "", True, "save displacement fields for the damage figure", "report",
              flag=True),
        Param("progress_every", "", 2000, "progress print interval in steps", "report", type=int),
    )
