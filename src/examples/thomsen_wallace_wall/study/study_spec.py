"""RW2 as a `StudySpec` for the shared harness (`rclattice.study`, D103/D105)."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from rclattice.study.protocols import LADDER, ProtocolSet                     # noqa: E402
from rclattice.study.references import ReferenceSet                            # noqa: E402
from rclattice.study.spec import PageConfig, RunsheetConfig, Selectors, StudySpec  # noqa: E402

import models                                                                  # noqa: E402
import params as P                                                             # noqa: E402
import specimen                                                                # noqa: E402

OUT_ROOT = specimen.OUT / "study"
MEASURED_F_KN = 163.284      # Table 4
MEASURED_K_KNMM = 35.19      # Table 4 — NOT a safe target: exceeds the uncracked section (summary.py)
AYDIN = {1.5: (32.26, 169.834), 3.01: (38.38, 174.544)}   # horizon -> (K kN/mm, F kN), Table 4


def continuum(params: dict):
    from build import wall_continuum
    con, _e = wall_continuum(grid=params["grid"], mesh_size=float(params["mesh"]))
    return con


def report_calibration(cal, params: dict) -> None:
    from build import report_calibration as rc
    rc(cal, mesh_size=float(params["mesh"]), horizon=float(params["horizon"]), grid=params["grid"])


def damage_thresholds(params: dict, meta: dict) -> tuple:
    g = specimen.GRADES["wall"]
    if params["comp"] == "capped":
        return g.ft / g.E, meta["fc_cap_MPa"] / g.E, "yielded in compression"
    if params["comp"] == "linear":
        return g.ft / g.E, None, "crushed"
    return g.ft / g.E, g.epsc0, "crushed"


SELECTORS = Selectors(lateral_loads=specimen.lateral_loads, control_node=specimen.control_node,
                      drive_nodes=specimen.drive_nodes, base_nodes=specimen.base_nodes)

REFERENCES = ReferenceSet(
    path=HERE.parent / "data" / "fig9b.npz",
    measured_F_kN=MEASURED_F_KN, measured_K_kNmm=MEASURED_K_KNMM,
    author_by_horizon=AYDIN, author_curves={1.5: "aydin15", 3.01: "aydin301"},
    clip_mm=float("inf"), record_extent_mm=(-72.0, 72.0),
    figure="Fig. 9(b)", test_label="Thomsen & Wallace 1995/2004", author_label="Aydin",
    f_source="Table 4", k_source="Table 4",
    note=("The Fig. 9(b) record is NOT clipped by its frame (loops reach ~+-72 mm inside +-80), "
          "so the envelope is the whole outline; the legend occludes the lower right (positive "
          "displacement below about -100 kN), which is the unloading side of the positive loops, "
          "not their tips. The cloud carries no order, so per-cycle energy is not recoverable."),
)

# The test was cyclic at increasing amplitudes (2019 paper, "increasing cyclic lateral displacement
# excursions"), but the LEVELS are not printed; the Fig. 9(b) cloud shows about eight loop
# amplitudes out to ~72 mm. So the ladder is invented in SHAPE and constrained in reach.
PRESETS = {
    "none": (),
    "probe": (0.001, 0.0025, 0.005),
    "to1p0": (0.001, 0.0025, 0.005, 0.0075, 0.010),
    "to2p0": (0.001, 0.0025, 0.005, 0.0075, 0.010, 0.015, 0.020),
}
# hours per mm of drive path: MEASURED by preflight.py on the graded 25 mm grid (see D105); the
# bonded value is scaled by the Aldemir bond/perfect ratio 0.227/0.0284 = 8.0 until measured.
H_PER_MM = {"perfect": 0.115, "bond": 0.92}
PROTOCOLS = ProtocolSet(
    presets=PRESETS, h_per_mm=H_PER_MM, ladder=LADDER, invented=True,
    invented_note=("THE PROTOCOL LEVELS ARE INVENTED: the 2019 paper says the test was cycled at\n"
                   "increasing amplitudes but prints no levels; the cloud shows ~8 amplitudes to\n"
                   "~72 mm. Drift capacity from a cyclic run is therefore a model property, not a\n"
                   "comparison, until the primary source is in hand."))

SOURCES = {
    "grid": ("assumed", "'rebar' puts every bar axis on a node exactly (D104); 'uniform' snaps "
                        "them by up to half a cell"),
    "mesh": ("assumed", "25 mm target on the graded grid (spacings 19.0-25.5); Aydin used 19"),
    "horizon": ("measured", "his own recommendation for RC (Table 4 runs 1.5 and 3.01)"),
    "field": ("assumed", "Stage 0: uniaxial 0.951, equibiaxial 0.821 against the same-grid continuum"),
    "nu": ("assumed", "not printed; the repo's convention"),
    "gf": ("convention", "Table 1 gives 75 N/m with no source; NOT a neutral knob (D75)"),
    "fcx": ("assumed", "strut-strength scale under `capped` only"),
    "comp": ("assumed", "which compression law is being tested"),
    "tail": ("measured", "'paper' = Table 1's printed a2/a3 (80/350); 'solved' = the same Gf, "
                         "re-solved at this mesh"),
    "bond": ("assumed", "bar type and bond are not reported; perfect bond is the default"),
    "bond_damage": ("assumed", "D95"),
    "rebar_top": ("assumed", "bars run to the driven top row (D78/D84)"),
    "rebar": ("measured", "Fig. 9(a): 8-#3 per boundary element, 8-#2 web, #2 @ 191 horizontal, "
                          "4.76 mm hoops @ 76"),
    "fc": ("measured", "Table 1"),
    "ft": ("measured", "Table 1, no footnote (= 0.31 sqrt(fc))"),
    "fy": ("measured", "Table 1 — the NOMINAL Grade 60; measured coupons are not printed"),
    "steel_rupture": ("assumed", "not printed; decides capacity (D97/D102)"),
    "concrete_residual": ("assumed", "D91: the model has no capacity until this is lowered"),
    "steel_b": ("assumed", "not printed; the largest control on capacity (D101/D102)"),
    "rate": ("assumed", "the repo's cross-study drive speed; the licence is the residual (D64)"),
    "damping": ("assumed", "damping is the knob that matters, not rate (D64)"),
    "proto": ("assumed", "levels invented; the test's are not printed"),
    "cycles": ("assumed", "one per level"),
    "drift": ("assumed", "how far this run was asked to go — NOT a capacity"),
    "analysis": ("convention", "which runner is used"),
    "integrator": ("convention", "explicit CentralDifference (D74)"),
    "steps_per_period": ("convention", "0 = sized from critical_time_step (D74)"),
    "steps": ("convention", "displacement-control steps, static solver only"),
    "groups": ("convention", "reporting only"),
    "capture": ("convention", "reporting only"),
    "progress_every": ("convention", "reporting only"),
}


def licence(params: dict) -> tuple[list, list]:
    fair = ["peak strength", "elastic stiffness against the same-grid continuum",
            "damage location and the flexural crack pattern (Fig. 9(c))",
            "the load-path split between plane-sections and truss action",
            "the controlled effect of this cell's own axis"]
    unfair = [("drift capacity", "the protocol levels are invented and the measured coupon "
                                 "strengths, hardening and rupture strain are not printed"),
              ("cyclic energy and per-cycle degradation", "the digitized record is an unordered "
                                                          "point cloud"),
              ("the measured initial stiffness", "Table 4's 35.19 kN/mm exceeds the uncracked "
                                                 "section; judge stiffness on K/K_continuum"),
              ("bar buckling", "Steel02 has none, and it is the likely end of the physical test")]
    if params["comp"] == "linear":
        unfair.append(("anything about compressive capacity",
                       "compression is linear forever in Aydin's published law"))
    if params["comp"] == "capped":
        unfair.append(("a compression-driven collapse", "the EPP branch cannot crush"))
    if params["analysis"] in ("pushover", "static"):
        unfair.append(("loop shape, pinching, residual drift", "a monotonic run"))
    return fair, unfair


COMP = {"crushing": "Concrete02 in compression: parabola to $f_c$, linear softening to $f_{cu}$, "
                    "flat residual beyond",
        "capped": "elastic--perfectly-plastic in compression, capped at $f_c$, no crushing",
        "linear": "linear elastic in compression, never yields"}
TAIL = {"solved": "$a_2,a_3$ re-solved at this mesh so the tail dissipates $G_f$ per strut",
        "paper": "the printed $a_2,a_3$ of Table 1 (80/350), transplanted from the 19\\,mm grid"}
BOND_SHORT = {"perfect": "perfect bond", "bond60": "bond links, $0.6f_t$ residual",
              "bond70": "bond links, $0.7f_t$ residual"}
BOND = {"perfect": "perfect bond: rebar shares the concrete nodes",
        "bond60": "bond links, residual plateau $0.6\\,f_t$",
        "bond70": "bond links, residual plateau $0.7\\,f_t$"}
TITLE_FMT = {
    "steel_rupture": lambda v: f"$\\varepsilon_{{su}}={float(v):g}$",
    "concrete_residual": lambda v: f"residual {float(v):g}",
    "steel_b": lambda v: f"$b={float(v):g}$",
    "damping": lambda v: f"$\\zeta={float(v):g}$",
    "rate": lambda v: f"{float(v):g}\\,mm/s",
    "gf": lambda v: f"$G_f={float(v):g}$",
    "fc": lambda v: f"$f_c={float(v):g}$",
    "ft": lambda v: f"$f_t={float(v):g}$",
    "fy": lambda v: f"$f_y={float(v):g}$",
    "nu": lambda v: f"$\\nu={float(v):g}$",
    "horizon": lambda v: f"horizon {float(v):g}",
    "grid": lambda v: f"{v} grid",
    "field": lambda v: f"{v} field",
    "cycles": lambda v: f"{int(v)} cycles/level",
    "bond_damage": lambda v: "bond damage",
    "rebar": lambda v: "no reinforcement" if not v else "",
    "rebar_top": lambda v: "" if v else "bars stop below the top",
    "mesh": lambda v: f"mesh {float(v):g}\\,mm",
    "drift": lambda v: "",
    "proto": lambda v: "",
}


def cell_title(cell: tuple) -> str:
    comp, tail, bond = cell
    return f"{comp} compression --- {tail} tension tail --- {BOND_SHORT.get(bond, bond)}"


RUNSHEET = RunsheetConfig(
    nu=specimen.NU, steel=(specimen.FY, specimen.ES, specimen.STEEL.b),
    axis_desc={"comp": COMP, "tail": TAIL, "bond": BOND},
    axis_rows=(("compression branch", "comp", "material"), ("tension tail", "tail", "material"),
               ("bond", "bond", "bond")),
    title_fmt=TITLE_FMT,
    record_label="test, measured record (Thomsen & Wallace 1995/2004, via Aydin 2019 Fig. 9b)",
    author_label="Aydin et al. (2019) simulation, horizon 1.5d and 3.01d",
    figure_caption_push=("Pushover over the measured cyclic record, its envelope, and the two "
                         "simulated curves of Aydin et al. (2019); full range above, the record's "
                         "own extent below."),
    figure_caption_cyclic=("Hysteresis and its backbone over the measured record and its envelope; "
                           "full range above, the record's own extent below."),
    section_title=cell_title,
    sort_key=lambda prm: (float(prm.get("steel_rupture") or 0),
                          float(prm.get("concrete_residual", 0.2)),
                          float(prm["steel_b"] if prm.get("steel_b") is not None else 0.01)),
    doc_root=HERE.parents[2] / "doc" / "reports" / "thomsen_wallace_runs",
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
    key="thomsen_wallace_wall",
    title="Thomsen & Wallace RW2",
    specimen=specimen.SPECIMEN,
    citation=("Thomsen & Wallace (1995, 2004) wall RW2, as simulated by Aydin, Tuncay & Binici "
              "(2019); every number second-hand through the 2019 paper."),
    registry=P.REG,
    out_root=OUT_ROOT,
    build=models.build,
    selectors=SELECTORS,
    references=REFERENCES,
    protocols=PROTOCOLS,
    drift_denominator=models.drift_denominator,
    panel=models.panel,
    axes=("comp", "tail", "bond"),
    continuum=continuum,
    report_calibration=report_calibration,
    load_path_groups=models.load_path_groups,
    element_modulus=models.element_modulus,
    damage_thresholds=damage_thresholds,
    bond_axis="bond", perfect_bond="perfect",
    cell_label=lambda p: f"{p['comp']}/{p['tail']}, {p['bond']}",
    variant_keys=("grid", "mesh", "field", "damping", "steel_rupture", "concrete_residual",
                  "steel_b", "bond_damage", "gf", "fy"),
    cross_run_axes=("grid", "mesh", "field", "horizon", "gf", "steel_b", "steel_rupture"),
    sources=SOURCES,
    licence=licence,
    calibration_note=("uniaxial (D47 default; Stage 0 measured it at 0.951 of the same-grid "
                      "continuum against 0.821 for the published equibiaxial route — the D53 "
                      "flexural under-read, not a fit)"),
    matrix_note=("One row per (cell, ANALYSIS); the run shown in a cell is the newest of its kind, "
                 "and everything non-default about it is named (D92/D98)."),
    author_note="his BEST specimen on force (1.04x), against 1.21x on Aldemir",
    runsheet=RUNSHEET,
    page=PAGE,
)
