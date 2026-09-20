"""The Aldemir wall as a `StudySpec` for the shared harness (`rclattice.study`, D103).

Everything specimen-flavoured that used to be spread over `run.py`, `report.py`, `master.py`,
`references.py`, `protocols.py`, `report_payload.py` and the run-sheet generator is gathered here;
the harness supplies the rest. The registry stays in `params.py` (its records are verbatim), the
parameter-to-model mapping stays in `models.py`.

`specimen.py` reads `ALDEMIR_TW` at IMPORT time, so every path into the parent package is lazy and
`prepare()` sets the environment from the run's own `tw` before anything is imported.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
os.environ.setdefault("ALDEMIR_TW", "210.0")

from rclattice.study.protocols import ProtocolSet                              # noqa: E402
from rclattice.study.references import ReferenceSet                            # noqa: E402
from rclattice.study.spec import PageConfig, RunsheetConfig, Selectors, StudySpec  # noqa: E402

import params as P                                                             # noqa: E402

OUT_ROOT = HERE.parent.parent / "output" / "aydin_aldemir_wall" / "study"
MEASURED_F_KN = 963.592      # Table 4
MEASURED_K_KNMM = 1038.44    # Table 4
AYDIN = {1.5: (943.16, 1164.413), 3.01: (1052.37, 1325.675)}   # horizon -> (K kN/mm, F kN)


# --- lazy bridges into the parent package ---------------------------------------------------------

def prepare(params: dict) -> None:
    os.environ["ALDEMIR_TW"] = repr(float(params["tw"]))


def _models():
    import models
    return models


def _specimen():
    import specimen
    return specimen


def build(params: dict):
    return _models().build(params)


def continuum(params: dict):
    from build import wall_continuum
    length, height = _models().panel(params)
    con, _e = wall_continuum(mesh_size=float(params["mesh"]), length=length, height=height)
    return con


def report_calibration(cal, params: dict) -> None:
    from build import report_calibration as rc
    rc(cal, mesh_size=float(params["mesh"]), horizon=float(params["horizon"]))


def load_path_groups(model, params: dict):
    return _models().load_path_groups(model, params)


def element_modulus(e) -> float:
    return _models().element_modulus(e)


def damage_thresholds(params: dict, meta: dict) -> tuple:
    """Under `capped` there IS no crushing: the branch is perfectly plastic, so the class is
    compressive YIELDING at fc_cap/E and the legend has to say so (PLAN §3). Under `linear`
    compression is linear forever and there is no compressive limit at all."""
    sp = _specimen()
    if params["comp"] == "capped":
        return sp.FT / sp.EC, meta["fc_cap_MPa"] / sp.EC, "yielded in compression"
    if params["comp"] == "linear":
        return sp.FT / sp.EC, None, "crushed"
    return sp.FT / sp.EC, sp.EPSC0, "crushed"


SELECTORS = Selectors(
    lateral_loads=lambda model, **kw: _specimen().lateral_loads(model, **kw),
    control_node=lambda model, **kw: _specimen().control_node(model, **kw),
    drive_nodes=lambda model, **kw: _specimen().drive_nodes(model, **kw),
    base_nodes=lambda model: _specimen().base_nodes(model),
)


# --- the measured record ---------------------------------------------------------------------------
REFERENCES = ReferenceSet(
    path=HERE.parent / "data" / "fig10b.npz",
    measured_F_kN=MEASURED_F_KN, measured_K_kNmm=MEASURED_K_KNMM,
    author_by_horizon=AYDIN, author_curves={1.5: "aydin15", 3.01: "aydin301"},
    clip_mm=16.0, record_extent_mm=(-18.91, 16.00),
    figure="Fig. 10(b)", test_label="Aldemir et al. 2017", author_label="Aydin",
    f_source="Table 4", k_source="Table 4",
    note=("Two limits of Fig. 10(b) carry into every row — the frame **clips at +16 mm** while "
          "the test reached ~20 mm, so the envelope is an UPPER BOUND past ~15 mm rather than the "
          "backbone; and the legend occludes the negative quadrant at positive displacement."),
)

# --- protocols (PLAN §6): THE PROTOCOL IS INVENTED ------------------------------------------------
import protocols as _pro                                                       # noqa: E402

PROTOCOLS = ProtocolSet(
    presets=_pro.PRESETS, h_per_mm=_pro.H_PER_MM, ladder=_pro.LADDER, invented=True,
    invented_note=("THE PROTOCOL IS INVENTED: the 2019 paper never prints a loading history and the\n"
                   "2017 test paper is not in the repo, so drift capacity and cyclic energy are NOT\n"
                   "fair claims from this run (PLAN §6)."))

# --- reporting text ------------------------------------------------------------------------------
import report as _rep                                                          # noqa: E402

MATRIX_NOTE = ("One row per (cell, ANALYSIS). Keying the matrix on the three model axes alone hid "
               "every run that was not the newest in its cell: all three cyclic runs sat in "
               "`crushing/solved/perfect` behind a later pushover, and a cell whose newest run was "
               "elastic showed a row of dashes as though it had produced nothing.")

# --- the run sheets (doc/reports/aydin_aldemir_runs) ----------------------------------------------
COMP = {"crushing": "Concrete02 in compression: parabola to $f_c$, linear softening to $f_{cu}$, "
                    "flat residual beyond",
        "capped": "elastic--perfectly-plastic in compression, capped at $f_c$, no crushing",
        "linear": "linear elastic in compression, never yields"}
TAIL = {"solved": "$a_2,a_3$ re-solved at this mesh so the tail dissipates $G_f$ per strut",
        "paper": "the printed $a_2,a_3$ of Table 1, transplanted from the 20\\,mm grid"}
BOND_SHORT = {"perfect": "perfect bond", "bond60": "bond links, $0.6f_t$ residual",
              "bond70": "bond links, $0.7f_t$ residual"}
BOND = {"perfect": "perfect bond: rebar shares the concrete nodes",
        "bond60": "bond links, residual plateau $0.6\\,f_t$",
        "bond70": "bond links, residual plateau $0.7\\,f_t$"}

# How a non-default parameter is written into a title. Anything not named here falls back to
# "<name> <value>", so a parameter added to the registry tomorrow still appears.
TITLE_FMT = {
    "steel_rupture": lambda v: f"$\\varepsilon_{{su}}={float(v):g}$",
    "concrete_residual": lambda v: f"residual {float(v):g}",
    "steel_b": lambda v: f"$b={float(v):g}$",
    "damping": lambda v: f"$\\zeta={float(v):g}$",
    "rate": lambda v: f"{float(v):g}\\,mm/s",
    "gf": lambda v: f"$G_f={float(v):g}$",
    "fc": lambda v: f"$f_c={float(v):g}$",
    "ft": lambda v: f"$f_t={float(v):g}$",
    "horizon": lambda v: f"horizon {float(v):g}",
    "tw": lambda v: f"$t={float(v):g}$\\,mm",
    "cycles": lambda v: f"{int(v)} cycles/level",
    "bond_damage": lambda v: "bond damage",
    "rebar": lambda v: "no reinforcement" if not v else "",
    "rebar_top": lambda v: "" if v else "bars stop below the top",
    "mesh": lambda v: f"mesh {float(v):g}\\,mm",
    "drift": lambda v: "",      # written last, by the caller, as a target or an amplitude
    "proto": lambda v: "",      # ditto
}


def cell_title(cell: tuple) -> str:
    comp, tail, bond = cell
    return f"{comp} compression --- {tail} tension tail --- {BOND_SHORT.get(bond, bond)}"


def _runsheet() -> RunsheetConfig:
    sp = _specimen()
    return RunsheetConfig(
        nu=sp.NU, steel=(sp.FY, sp.ES, sp.STEEL.b),
        axis_desc={"comp": COMP, "tail": TAIL, "bond": BOND},
        axis_rows=(("compression branch", "comp", "material"),
                   ("tension tail", "tail", "material"),
                   ("bond", "bond", "bond")),
        title_fmt=TITLE_FMT,
        record_label="test, measured record (Aldemir et al. 2017, Fig. 10b)",
        author_label="Aydin et al. (2019) simulation, horizon 1.5d and 3.01d",
        figure_caption_push=("Pushover over the measured cyclic record, its envelope, and the two "
                             "simulated curves of Aydin et al. (2019); full range above, the "
                             "record's own extent below."),
        figure_caption_cyclic=("Hysteresis and its backbone over the measured record and its "
                               "envelope; full range above, the record's own extent below."),
        section_title=cell_title,
        sort_key=lambda prm: (float(prm.get("steel_rupture") or 0),
                              float(prm.get("concrete_residual", 0.2)),
                              float(prm["steel_b"] if prm.get("steel_b") is not None
                                    else sp.STEEL.b)),
        doc_root=HERE.parents[2] / "doc" / "reports" / "aydin_aldemir_runs",
    )


PAGE = PageConfig(
    template=HERE / "report_page.template.html",
    comp_note={"crushing": "Concrete02: yields, softens, crushes",
               "capped": "elastic-perfectly-plastic: yields at f<sub>c</sub>, cannot crush",
               "linear": "Aydin's own law: compression never yields"},
    comp_order=("crushing", "capped", "linear"),
    out_html=HERE / "report_page.html",
    out_payload=HERE / "report_payload.json",
)

SPEC = StudySpec(
    key="aydin_aldemir_wall",
    title="Aldemir wall",
    specimen="Aldemir et al. (2017) squat wall",
    citation=("Aldemir, Binici & Canbay (2017) squat wall, as simulated by Aydin, Tuncay & Binici "
              "(2019)."),
    registry=P.REG,
    out_root=OUT_ROOT,
    build=build,
    selectors=SELECTORS,
    references=REFERENCES,
    protocols=PROTOCOLS,
    drift_denominator=lambda params: _models().drift_denominator(params),
    panel=lambda params: _models().panel(params),
    axes=("comp", "tail", "bond"),
    prepare=prepare,
    continuum=continuum,
    report_calibration=report_calibration,
    load_path_groups=load_path_groups,
    element_modulus=element_modulus,
    damage_thresholds=damage_thresholds,
    bond_axis="bond", perfect_bond="perfect",
    cell_label=lambda p: f"{p['comp']}/{p['tail']}, {p['bond']}",
    variant_keys=("mesh", "damping", "steel_rupture", "concrete_residual", "bond_damage", "gf", "tw"),
    cross_run_axes=("panel", "tw", "mesh", "horizon", "gf"),
    sources=_rep.SOURCES,
    licence=_rep.licence,
    calibration_note=("uniaxial (thesis default; measurably better than the published "
                      "equibiaxial route on this specimen, 0.92 vs 0.79)"),
    matrix_note=MATRIX_NOTE,
    author_note="the worst of his six specimens, which is why this wall is the target",
    runsheet=_runsheet(),
    page=PAGE,
)
