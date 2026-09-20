"""`StudySpec` — the contract a specimen fills in to run through the shared harness (D103).

Everything the harness needs to know about ONE specimen, in one object. The harness never imports
a specimen module; the specimen builds a `StudySpec` and hands it to `runner.main`, `master.build`,
`page.render`, `runsheet.main` and friends.

The split follows what was specimen-bound in the Aldemir study when it was lifted:

  * the REGISTRY (its axes, its defaults, its schema version, its legacy spellings);
  * `build(params)` — the parameter dict to a built model, plus calibration and a `meta` record;
  * the NODE SELECTORS (drive row, control node, base reaction set, the unit lateral pattern);
  * the MEASURED references (peak, stiffness, the digitized record, the author's own simulation);
  * the PROTOCOLS and their measured cost;
  * the REPORTING TEXT (sources of every parameter, what a run does and does not license, the
    calibration-field note, the matrix prose, the run-sheet vocabulary).

Anything not on this object is the harness's own and identical across studies.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from .protocols import ProtocolSet
from .references import ReferenceSet
from .registry import Registry


@dataclass
class Selectors:
    """Post-build node queries, all taking the built model (and the drift denominator where needed)."""

    lateral_loads: Callable[..., list]      # (model, height=) -> [Load] unit pattern on the drive row
    control_node: Callable[..., int]        # (model, height=) -> node id the drift is read at
    drive_nodes: Callable[..., list[int]]   # (model, height=) -> the driven row
    base_nodes: Callable[..., list[int]]    # (model) -> the reaction set


@dataclass
class RunsheetConfig:
    """The vocabulary of the LaTeX run sheets — everything specimen-flavoured in `runsheet.py`."""

    nu: float                                    # Poisson ratio used in the balance
    steel: tuple[float, float, float]            # (f_y, E_s, b) of the specimen's own grade
    axis_desc: dict[str, dict[str, str]]         # axis -> value -> LaTeX description (factsheet)
    axis_rows: tuple[tuple[str, str, str], ...]  # (row label, axis, "material" | "bond") in order
    title_fmt: dict[str, Callable[[Any], str]]   # non-default parameter -> title fragment
    record_label: str                            # legend: "test, measured record (…, Fig. 10b)"
    author_label: str                            # legend: "Aydin et al. (2019) simulation, …"
    figure_caption_push: str
    figure_caption_cyclic: str
    section_title: Callable[[tuple], str]        # cell tuple -> section heading
    sort_key: Callable[[dict], tuple] | None = None   # extra ordering within a cell
    doc_root: Path | None = None                 # the report directory (…/doc/reports/<study>_runs)
    date_glob: str = "2026-*"                    # run-directory prefix to discover


@dataclass
class PageConfig:
    """The advisor page: its template and the notes it renders per compression law."""

    template: Path
    comp_note: dict[str, str]
    comp_order: tuple[str, ...]
    out_html: Path
    out_payload: Path


@dataclass
class StudySpec:
    key: str                                   # output directory name, e.g. "aydin_aldemir_wall"
    title: str                                 # "Aldemir wall"
    specimen: str                              # one line naming the specimen, printed by the runner
    citation: str                              # the report header's second line
    registry: Registry
    out_root: Path                             # …/examples/output/<key>/study
    build: Callable[[dict], tuple]             # params -> (model, cal, meta)
    selectors: Selectors
    references: ReferenceSet
    protocols: ProtocolSet
    drift_denominator: Callable[[dict], float]
    panel: Callable[[dict], tuple[float, float]]
    axes: tuple[str, ...]                      # the matrix axes, in order
    # --- optional hooks, each with a sensible default --------------------------------------------
    prepare: Callable[[dict], None] | None = None     # runs BEFORE build (e.g. env for a module)
    continuum: Callable[[dict], object] | None = None # params -> elastic continuum twin (Stage 2)
    report_calibration: Callable[[Any, dict], None] | None = None   # (cal, params) -> prints
    load_path_groups: Callable[[object, dict], dict | None] | None = None
    element_modulus: Callable[[Any], float] | None = None   # for critical_time_step (D74)
    damage_thresholds: Callable[[dict, dict], tuple] | None = None  # (params, meta) ->
    #                                            (eps_crack, eps_crush | None, crush_label)
    bond_axis: str | None = "bond"             # the axis whose non-`perfect_bond` values are bonded
    perfect_bond: str = "perfect"
    elastic_target_drift: float = 0.0002       # an order below where the wall cracks
    cell_label: Callable[[dict], str] | None = None   # "linear/solved, perfect" for legends
    variant_keys: tuple[str, ...] = ("mesh", "damping", "steel_rupture", "concrete_residual",
                                     "bond_damage", "gf")
    cross_run_axes: tuple[str, ...] = ("mesh", "horizon", "gf")   # beyond `axes`
    sources: dict[str, tuple[str, str]] = field(default_factory=dict)   # report: where a value is from
    licence: Callable[[dict], tuple[list, list]] | None = None   # params -> (fair, unfair)
    calibration_note: str = ""                 # the "field" row of the calibration table
    matrix_note: str = ""                      # one paragraph under the matrix heading
    author_note: str = ""                      # one line about the author's own overshoot
    field_of: Callable[[dict], str] | None = None   # meta -> energy-balance field name
    runsheet: RunsheetConfig | None = None
    page: PageConfig | None = None

    # --- derived conveniences --------------------------------------------------------------------
    def is_bonded(self, params: dict) -> bool:
        if not self.bond_axis or self.bond_axis not in params:
            return False
        return params[self.bond_axis] != self.perfect_bond

    def label(self, params: dict) -> str:
        if self.cell_label:
            return self.cell_label(params)
        return "/".join(str(params[a]) for a in self.axes)

    def cell(self, params: dict) -> tuple:
        return tuple(params[a] for a in self.axes)

    @property
    def measured_F_kN(self) -> float:
        return self.references.measured_F_kN

    @property
    def measured_K_kNmm(self) -> float:
        return self.references.measured_K_kNmm
