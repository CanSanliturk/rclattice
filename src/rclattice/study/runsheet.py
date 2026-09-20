"""LaTeX run sheets — one per pushover or cyclic run (Aldemir doc/reports, lifted in D103).

    runsheet.main(SPEC, HERE)              # from a report directory's own generate.py
      --only 2026-09-10                    just the runs whose directory name contains this
      --list                               show what would be generated, write nothing

WHAT A RUN SHEET IS. A record of the model and the analysis, not a discussion of either: geometry
and discretisation, the calibration, the OpenSees materials with their actual argument lists, the
solver configuration, and the numbers the run produced. No interpretation.

WHERE THE FACTS COME FROM. `params.json` holds every parameter of a run and `data.json` holds the
solver configuration, the counts and the response series. The one thing neither records is the
OpenSees material ARGUMENTS, so this rebuilds each model from its own parameters and reads the
materials off the assembled model. The build is deterministic, so a rebuild reproduces exactly what
the run used — but it IS a rebuild, and that is stated on every sheet.

Contents entries derive from the registry so no two runs can read alike while differing in
something real; where they still would, the run's own timestamp is appended.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from . import metrics
from .spec import StudySpec

PCT = "\\%"                    # the percent unit, as LaTeX; `num` does not escape its unit


# ------------------------------------------------------------------------------------------------
# LaTeX helpers
# ------------------------------------------------------------------------------------------------
def tex(s) -> str:
    """Escape a value for LaTeX text. Applied to every string that reaches the document."""
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"),
                 ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def num(x, d=3, unit="") -> str:
    """A number as typeset. `unit` is LaTeX, not text — it is NOT escaped, so pass `\\%` for a
    percent sign.

    Trailing zeros are stripped only when there is a decimal point to strip them from: doing it
    unconditionally turned 3000 into "3," and 360 into "36".
    """
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "---"
    if isinstance(x, float):
        if abs(x) < 0.5 * 10.0 ** -d:
            x = 0.0                      # anything that rounds to zero prints as 0, not -0
        s = f"{x:,.{d}f}"
        if "." in s:
            s = s.rstrip("0").rstrip(".")
    else:
        s = f"{x:,}"
    return f"\\val{{{s}{unit}}}"


# OpenSees argument names, in the order the command takes them. A material not listed here falls
# back to positional names, so an unfamiliar type still reads rather than failing.
ARGNAMES = {
    "Concrete02": ["fpc", "epsc0", "fpcu", "epsU", "lambda", "ft", "Ets"],
    "Concrete01": ["fpc", "epsc0", "fpcu", "epsU"],
    "Steel02": ["Fy", "E0", "b", "R0", "cR1", "cR2", "a1", "a2", "a3", "a4", "sigInit"],
    "Steel01": ["Fy", "E0", "b", "a1", "a2", "a3", "a4"],
    "Elastic": ["E", "eta", "Eneg"],
    "ElasticPP": ["E", "epsyP", "epsyN", "eps0"],
    "MinMax": ["otherTag"],
    "Parallel": ["matTag"],          # repeated: matTag1, matTag2, ...
    "Series": ["matTag"],
    "ElasticMultiLinear": [],        # entirely flag-driven
}
REPEATING = {"Parallel", "Series"}


def fmt(a) -> str:
    if isinstance(a, float):
        return f"{a + 0.0:.6g}" if a != 0.0 else "0"
    return str(a)


def sci(x, d: int = 3) -> str:
    """A number in scientific notation — EA runs from 1e5 to 1e9 across the families."""
    if x is None or not math.isfinite(x):
        return "---"
    m, e = f"{x:.{d}e}".split("e")
    return f"\\val{{{m}e{int(e)}}}"


def arglist(mtype: str, args) -> str:
    """An OpenSees argument tuple as `name:value` pairs, in the order the command takes them.

    POSITIONAL arguments take their name from `ARGNAMES`. FLAG arguments — a string beginning with
    `-` — name the values that follow them, so `-min, -0.05` prints as `-min:-0.05` and a flag
    carrying a list prints as `-strain:[...]`; that is the only way `ElasticMultiLinear`, whose
    whole record is flags, reads at all.
    """
    names = ARGNAMES.get(mtype, [])
    out, i, flag, held = [], 0, None, []

    def flush():
        if flag is None:
            return
        v = held[0] if len(held) == 1 else "[" + ", ".join(held) + "]"
        out.append(f"{flag}:{v}" if held else flag)

    for a in args:
        if isinstance(a, str) and a.startswith("-"):
            flush()
            flag, held = a, []
            continue
        if flag is not None:
            held.append(fmt(a))
            continue
        if mtype in REPEATING:
            name = f"{names[0]}{i + 1}"
        elif i < len(names):
            name = names[i]
        else:
            name = f"arg{i + 1}"
        out.append(f"{name}:{fmt(a)}")
        i += 1
    flush()
    return "\\val{" + tex(", ".join(out)) + "}"


# ------------------------------------------------------------------------------------------------
# the runs
# ------------------------------------------------------------------------------------------------
def discover(spec: StudySpec, only: list[str] | None) -> list[Path]:
    out = []
    for d in sorted(spec.out_root.glob(spec.runsheet.date_glob)):
        pj, dj = d / "params.json", d / "data.json"
        if not (pj.is_file() and dj.is_file()):
            continue
        if json.loads(pj.read_text()).get("analysis") not in ("pushover", "cyclic"):
            continue
        if only and not any(o in d.name for o in only):
            continue
        out.append(d)
    return out


def slug(d: Path) -> str:
    return d.name.replace(".", "p")


def load_params(spec: StudySpec, d: Path) -> tuple[dict, list[str]]:
    """A run's parameters, with any that postdate it filled from today's defaults.

    The registry grows, so a run made before a parameter existed has no record of it. Filling the
    default is what makes such a run rebuildable at all; the names filled are returned and printed
    on the sheet, because a filled value is an assumption about the run rather than a record of it.
    """
    prm = json.loads((d / "params.json").read_text())
    missing = spec.registry.fill_defaults(prm)
    return prm, missing


def initial_E(mats: dict, tag: int) -> float | None:
    """The initial axial modulus a material record implies, in MPa.

    Concrete02 does not carry E as an argument — its initial tangent is `2*fpc/epsc0`, which is the
    number that actually stiffens the strut and the one the calibration was struck against. The
    wrappers are resolved to whatever they wrap.
    """
    m = mats.get(tag)
    if m is None:
        return None
    a = m.args
    try:
        if m.mtype in ("Concrete02", "Concrete01"):
            return abs(2.0 * float(a[0]) / float(a[1]))
        if m.mtype in ("Steel02", "Steel01"):
            return float(a[1])
        if m.mtype in ("Elastic", "ElasticPP"):
            return float(a[0])
        if m.mtype == "MinMax":
            return initial_E(mats, int(a[0]))
        if m.mtype == "Parallel":
            parts = [initial_E(mats, int(t)) for t in a if isinstance(t, int)]
            return sum(v for v in parts if v is not None) or None
        if m.mtype == "ElasticMultiLinear":
            st = list(a)
            i, j = st.index("-strain"), st.index("-stress")
            eps = [float(v) for v in st[i + 1:j]]
            sig = [float(v) for v in st[j + 1:]]
            pos = [k for k, e in enumerate(eps) if e > 0]
            k = pos[0]
            return sig[k] / eps[k]
    except (ValueError, IndexError, TypeError, ZeroDivisionError):
        return None
    return None


def inspect_model(spec: StudySpec, prm: dict) -> tuple[list[dict], list[dict]]:
    """Rebuild the model and read its uniaxial materials, with what each one is applied to.

    Elements carry their material tag as the last argument, so the element list gives both the
    usage counts and — for concrete struts, whose law is regularized per strut length — the lengths
    each material serves. Materials referenced by OTHER materials (Steel02 inside MinMax, the bond
    stack inside Parallel) are pulled in transitively and marked as such, because they never appear
    on an element.
    """
    model, _cal, _meta = spec.build(prm)
    mats = {m.id: m for m in model.uniaxial_materials}
    used: dict[int, dict] = {}
    for e in model.elements:
        tag = e.args[-1]
        if not isinstance(tag, int) or tag not in mats:
            continue
        rec = used.setdefault(tag, {"count": 0, "kinds": set(), "lengths": set()})
        rec["count"] += 1
        rec["kinds"].add(e.kind or "?")
        if len(e.nodes) == 2:
            (xa, ya), (xb, yb) = (model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords)
            rec["lengths"].add(round(math.hypot(xb - xa, yb - ya), 2))

    referenced: dict[int, set] = {}
    seen, frontier = set(used), list(used)
    while frontier:
        t = frontier.pop()
        for a in mats[t].args:
            if isinstance(a, int) and a in mats and a != t:
                referenced.setdefault(a, set()).add(t)
                if a not in seen:
                    seen.add(a)
                    frontier.append(a)

    rows = []
    for tag in sorted(set(used) | set(referenced)):
        m = mats[tag]
        if tag in used:
            r = used[tag]
            kinds = "/".join(sorted(r["kinds"]))
            where = f"{r['count']:,} {kinds} elements"
            if r["lengths"] and len(r["lengths"]) <= 3:
                where += r", L = " + ", ".join(f"{v:g}" for v in sorted(r["lengths"])) + r" mm"
        else:
            where = "referenced by tag " + ", ".join(str(t) for t in sorted(referenced[tag]))
        rows.append({"tag": tag, "type": m.mtype, "args": m.args, "where": where})

    # --- the strut families: EA and EA/L per (kind, length) --------------------------------------
    # A truss element carries its area as its first argument and its material tag as its last, so
    # the axial stiffness of every family in the model is recoverable without assuming anything.
    fam: dict[tuple, dict] = {}
    mesh = float(prm.get("mesh", spec.registry.by_name["mesh"].default))
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = (model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords)
        L = round(math.hypot(xb - xa, yb - ya), 2)
        A = float(e.args[0]) if isinstance(e.args[0], (int, float)) else float("nan")
        tag = e.args[-1]
        key = (e.kind or "?", L, round(A, 4), tag)
        f = fam.setdefault(key, {"kind": e.kind or "?", "L": L, "A": A, "tag": tag,
                                 "etype": e.etype, "count": 0})
        f["count"] += 1
    fams = []
    for f in sorted(fam.values(), key=lambda v: (v["kind"], v["L"])):
        E = initial_E(mats, f["tag"])
        f["E"] = E
        f["EA"] = (E * f["A"]) if E else None
        f["EA_L"] = (E * f["A"] / f["L"]) if E and f["L"] else None
        # Name the CONCRETE families by orientation — that split is what the horizon rule creates,
        # and it is why they carry different Ets and different EA/L.
        o = ""
        if f["kind"] == "concrete":
            if abs(f["L"] - mesh) < 0.01 * mesh:
                o = ", orthogonal"
            elif abs(f["L"] - mesh * math.sqrt(2.0)) < 0.01 * mesh:
                o = ", diagonal"
        f["name"] = f["kind"] + o
        fams.append(f)
    return rows, fams


# ------------------------------------------------------------------------------------------------
# figures
# ------------------------------------------------------------------------------------------------
def reference_layers(spec: StudySpec, ax, *, aydin: bool, mm_per_pct: float):
    """The measured record, its envelope and — for a pushover — the author's own simulated curves."""
    import matplotlib.pyplot as plt  # noqa: F401  (backend chosen by the caller)
    refs, rs = spec.references, spec.runsheet
    ref = refs.load()
    if ref is None:
        return []
    handles = []
    c = ref["cloud"]
    ax.scatter(c[:, 0] / mm_per_pct, c[:, 1], s=0.6, c="0.62", alpha=0.30,
               linewidths=0, zorder=1, label="test, measured record")
    handles.append(rs.record_label)
    bb = ref["backbone"]
    o = np.argsort(bb[:, 0])
    ax.plot(bb[o, 0] / mm_per_pct, bb[o, 1], color="#B03A2E", lw=1.5, ls="--", zorder=3,
            label="test backbone")
    handles.append("test backbone (running maximum of the record)")
    if aydin:
        styles = [("#2E6DA4", 1.2, "-."), ("#7A5195", 1.0, ":")]
        for (h, key), (col, lw, ls) in zip(sorted(refs.author_curves.items()), styles):
            if key + "_x" not in ref:
                continue
            ax.plot(ref[key + "_x"] / mm_per_pct, ref[key + "_y"], color=col, lw=lw, ls=ls,
                    zorder=3, label=f"{refs.author_label}, horizon {h:g}d")
        handles += [rs.author_label]
    return handles


def make_figure(spec: StudySpec, d: Path, prm: dict, dat: dict, out: Path) -> None:
    """Two stacked panels: the run's full range, then the range the measured record covers."""
    import matplotlib.pyplot as plt

    height = spec.drift_denominator(prm)
    mm_per_pct = height / 100.0
    lo_mm, hi_mm = spec.references.record_extent_mm
    disp = np.asarray(dat["disp"], float)
    shear = np.asarray(dat["shear"], float) / 1e3
    drift = disp / height * 100.0
    cyclic = prm["analysis"] == "cyclic"

    turns = metrics.turning_points(disp) if cyclic else []
    sm = metrics.moving_average(shear, int(dat.get("peak_smoothing_samples") or 125))

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 7.4))
    for ax, focused in zip(axes, (False, True)):
        reference_layers(spec, ax, aydin=not cyclic, mm_per_pct=mm_per_pct)
        if cyclic:
            ax.plot(drift, shear, color="0.20", lw=0.45, alpha=0.85, zorder=4,
                    label="this run, hysteresis")
            if len(turns) >= 2:
                td, tv = drift[turns], sm[turns]
                for sgn in (+1, -1):
                    m = (td * sgn) > 0
                    o = np.argsort(np.abs(td[m]))
                    ax.plot(td[m][o], tv[m][o], color="#1F3864", lw=1.8, marker="o", ms=3.0,
                            zorder=5, label="this run, backbone" if sgn > 0 else None)
        else:
            ax.plot(drift, shear, color="0.10", lw=1.6, zorder=5, label="this run, pushover")

        if focused:
            lo, hi = lo_mm / mm_per_pct * 1.06, hi_mm / mm_per_pct * 1.06
            ax.set_xlim(lo, hi)
            inside = (drift >= lo) & (drift <= hi)
            ym = max(np.abs(shear[inside]).max() if inside.any() else 0.0,
                     spec.measured_F_kN) * 1.10
            ax.set_ylim(-ym, ym)
        else:
            ext = max(abs(lo_mm), abs(hi_mm))
            xm = max(np.abs(drift).max(), ext / mm_per_pct) * 1.03
            ax.set_xlim(-xm if cyclic else -ext / mm_per_pct * 1.03, xm)
        ax.axhline(0, color="0.7", lw=0.6, zorder=0)
        ax.axvline(0, color="0.7", lw=0.6, zorder=0)
        ax.set_ylabel("base shear (kN)")
        ax.grid(True, color="0.90", lw=0.5)
        ax.set_title("focused on the measured record's extent" if focused else "full range",
                     fontsize=8, color="0.30", loc="left")
    axes[0].legend(fontsize=7, loc="upper left", framealpha=0.92)
    axes[0].set_title(d.name, fontsize=6.5, color="0.45", loc="right")
    axes[1].set_xlabel(f"drift (%)   [top displacement / {height:g} mm]")
    fig.tight_layout()
    fig.savefig(out, format="pdf")
    plt.close(fig)


# ------------------------------------------------------------------------------------------------
# the sheet
# ------------------------------------------------------------------------------------------------
# In the heading but never in the contents entry: the part and section already carry them.
def run_titles(spec: StudySpec, prm: dict, dat: dict) -> tuple[str, str]:
    """(heading, contents entry).

    The entry names the mesh, the calibration, and EVERY parameter of this run that differs from
    the registry default — derived from the registry rather than from a hand-picked list, so two
    runs can never come out with the same entry while differing in something real.
    """
    rs = spec.runsheet
    in_section = {"analysis", *spec.axes}
    bits = [f"mesh {float(prm['mesh']):g}\\,mm",
            f"energy balance, {dat.get('meta', {}).get('field') or 'uniaxial'} field"]
    for q in spec.registry:
        if q.name in in_section or q.affects == "report" or q.name == "mesh":
            continue
        v = prm.get(q.name)
        if v == q.default or v is None:
            continue
        f = rs.title_fmt.get(q.name)
        piece = f(v) if f else f"{tex(q.name.replace('_', ' '))} {tex(v)}"
        if piece:
            bits.append(piece)
    if prm["analysis"] == "cyclic":
        bits.append(f"{prm.get('proto')} to $\\pm${float(prm['drift'])*100:g}\\%")
    else:
        bits.append(f"to {float(prm['drift'])*100:g}\\% drift")
    short = " --- ".join(bits)
    head = f"{prm['analysis'].capitalize()} --- {'/'.join(spec.cell(prm))} --- " + short
    return head, short


def sheet(spec: StudySpec, d: Path, *, first_in_section: bool = False, stamp: str = "") -> str:
    rs = spec.runsheet
    prm, filled = load_params(spec, d)
    dat = json.loads((d / "data.json").read_text())
    meta = dat.get("meta", {})
    L = []
    w = L.append

    head, short = run_titles(spec, prm, dat)
    if stamp:                      # only when two runs in a section would otherwise read alike
        short += f" --- {stamp}"
    if not first_in_section:
        w(r"\clearpage")
    w(r"\subsection[" + short + "]{" + head + "}")
    w(r"\runid{" + tex(d.name) + "}")

    # --- specimen and discretisation -----------------------------------------------------------
    panel = meta.get("panel_mm") or list(spec.panel(prm))
    w(r"\begin{factsheet}{Specimen and discretisation}")
    w(rf"panel (L $\times$ H $\times$ t) & {num(panel[0],0)} $\times$ {num(panel[1],0)} $\times$ "
      rf"{num(meta.get('thickness_mm'),0,' mm')} \\")
    w(rf"drift denominator & {num(meta.get('drift_denominator_mm'),0,' mm')} \\")
    w(rf"mesh size & {num(float(prm['mesh']),1,' mm')} \\")
    w(rf"strut horizon & {num(float(prm['horizon']),2)} $\times$ mesh \\")
    w(rf"nodes / elements & {num(dat.get('nodes'))} / {num(dat.get('elements'))} \\")
    w(rf"composition & {tex(meta.get('describe',''))} \\")
    w(rf"reinforcement & {tex('included' if prm.get('rebar', True) else 'omitted')}"
      rf"{tex(', bars run to the top face' if prm.get('rebar_top') else '')} \\")
    w(rf"strut life $\varepsilon_u/\varepsilon_{{cr}}$ & orthogonal "
      rf"{num(meta.get('strut_life_orthogonal'),1)}, diagonal "
      rf"{num(meta.get('strut_life_diagonal'),1)} \\")
    w(rf"first cracking & {num((meta.get('V_cr_N') or 0)/1e3,0,' kN')} \\")
    w(r"\end{factsheet}")

    # --- calibration ---------------------------------------------------------------------------
    # The formula is repeated on every sheet on purpose: a run sheet is meant to be readable on its
    # own, and only two of its terms (the affine field, and nu) vary between runs.
    field = meta.get("field", "uniaxial")
    w(r"\begin{factsheet}{Calibration}")
    w(r"strategy & Aydin's overlapping-lattice energy balance: the stored elastic energy of the "
      r"continuum and of the lattice are equated under an affine strain field and solved for one "
      r"uniform $E_tA_t$. A closed-form sum over the strut list --- no FE solve, no reference "
      r"model, no optimiser. \\")
    w(rf"affine field & {tex(field)} \\")
    w(rf"Poisson ratio $\nu$ in the balance & {num(rs.nu,2)} \\")
    w(rf"strut area $A_t$ & {num(meta.get('area_mm2'),1)}\,mm$^2$ \\")
    w(rf"$E_tA_t$ per strut & {num((meta.get('EA_N') or 0)/1e6,1,'e6 N')} \\")
    w(rf"$G_f$ multiple of the tail & {num(dat.get('gf_multiple'),2)} "
      rf"({tex(dat.get('gf_multiple_note',''))}) \\")
    w(r"\end{factsheet}")

    w(r"Under an affine displacement field $\mathbf{u}=\mathbf{F}\mathbf{x}$ a strut spanning "
      r"$\mathbf{v}_s=\mathbf{x}_b-\mathbf{x}_a$, of length $L_s=\lVert\mathbf{v}_s\rVert$, "
      r"carries the axial strain")
    w(r"\begin{equation*}"
      r"\varepsilon_s=\frac{(\mathbf{F}\mathbf{v}_s)\cdot\mathbf{v}_s}{L_s^{2}},"
      r"\qquad "
      r"W_\text{lattice}=\frac{E_tA_t}{2}\sum_s \varepsilon_s^{2}L_s ."
      r"\end{equation*}")
    if field == "equibiaxial":
        w(r"The continuum energy is taken under equal stresses $\sigma_x=\sigma_y$ "
          r"(2019 Appendix, Eqs 3--4), giving $\varepsilon_x=\varepsilon_y=\varepsilon$ and")
        w(r"\begin{equation*}"
          r"W_\text{continuum}=\frac{E\,\varepsilon^{2}V}{1-\nu},\qquad V=t\,A ."
          r"\end{equation*}")
    else:
        w(r"The continuum energy is taken under uniaxial strain with the transverse strain "
          r"restrained ($\varepsilon_x=\varepsilon$, $\varepsilon_y=0$), plane stress:")
        w(r"\begin{equation*}"
          r"W_\text{continuum}=\frac{E\,\varepsilon^{2}V}{2\left(1-\nu^{2}\right)},"
          r"\qquad V=t\,A ."
          r"\end{equation*}")
    w(r"Equating the two and solving for the one unknown,")
    w(r"\begin{equation*}"
      r"E_tA_t=\frac{W_\text{continuum}}"
      r"{\tfrac{1}{2}\sum_s \varepsilon_s^{2}L_s},\qquad A_t=\frac{E_tA_t}{E} ."
      r"\end{equation*}")
    w(r"Both sides are quadratic in $\varepsilon$, so the strain magnitude cancels exactly and "
      r"the result is a property of the grid, the horizon, $\nu$ and the thickness alone --- it "
      r"transfers to any model built on the same grid. For an interior node the equibiaxial "
      r"balance collapses to the published closed form $E_tA_t = C\,E_t\,d\,w$ with "
      r"$C=0.621$ at a $1.5d$ horizon and $C=0.102$ at $3.01d$ (2019 Appendix, Eq 6, at "
      r"$\nu=1/3$).")

    # --- constitutive inputs ---------------------------------------------------------------------
    def axis_row(label, axis):
        v = prm[axis]
        w(rf"{label} & {rs.axis_desc.get(axis, {}).get(v, tex(v))} \\")

    fy, es, b_spec = rs.steel
    w(r"\begin{factsheet}{Constitutive inputs}")
    w(rf"concrete & $f_c$ = {num(float(prm['fc']),1,' MPa')}, "
      rf"$f_t$ = {num(float(prm['ft']),2,' MPa')}, "
      rf"$G_f$ = {num(float(prm['gf']),4,' N/mm')} \\")
    for label, axis, where in rs.axis_rows:
        if where == "material":
            axis_row(label, axis)
    w(rf"crushing residual & {num(float(prm.get('concrete_residual',0.2)),2)} $\times f_c$ \\")
    b = float(prm['steel_b']) if prm.get('steel_b') is not None else b_spec
    w(rf"steel & $f_y$ = {num(fy,0,' MPa')}, $E_s$ = {num(es,0,' MPa')}, "
      rf"$b$ = {num(b,4)} \\")
    sr = float(prm.get('steel_rupture') or 0)
    w(rf"bar rupture & {(num(sr,3) + ' strain') if sr else 'none --- bars cannot break'} \\")
    for label, axis, where in rs.axis_rows:
        if where == "bond":
            axis_row(label, axis)
    w(r"\end{factsheet}")

    # --- materials -------------------------------------------------------------------------------
    mat_rows, fam_rows = inspect_model(spec, prm)
    w(r"\begin{factsheet}{OpenSees materials}")
    w(r"provenance & rebuilt from this run's own parameters; the build is deterministic, so these "
      r"are the records the run used \\")
    if filled:
        w(r"parameters filled from today's defaults & \val{" + tex(", ".join(filled)) +
          r"} --- this run predates them \\")
    w(r"\end{factsheet}")
    w(r"\begin{matsheet}")
    for m in mat_rows:
        w(f"{m['tag']} & \\val{{{tex(m['type'])}}} & {arglist(m['type'], m['args'])} & "
          f"{tex(m['where'])} \\\\")
    w(r"\end{matsheet}")

    # --- strut families: axial stiffness ---------------------------------------------------------
    w(r"\begin{factsheet}{Strut families --- axial stiffness}")
    w(r"note & $E$ is the material's INITIAL axial modulus; for Concrete02 that is $2f_c/"
      r"\varepsilon_{c0}$, which is the modulus the calibration was struck against, not an "
      r"argument of the record \\")
    w(r"\end{factsheet}")
    w(r"\begin{famsheet}")
    for f in fam_rows:
        w(f"{tex(f['name'])} & {num(f['count'])} & \\val{{{tex(f['etype'])}}} & "
          f"{num(f['L'],2)} & {num(f['A'],1)} & {num(f['E'],0)} & "
          f"{sci(f['EA'])} & {sci(f['EA_L'])} \\\\")
    w(r"\end{famsheet}")

    # --- solver ------------------------------------------------------------------------------------
    integ = dat.get("integrator")
    integ = ", ".join(str(x) for x in integ) if isinstance(integ, list) else str(integ)
    w(r"\begin{factsheet}{Analysis and solver}")
    w(rf"analysis & {tex(prm['analysis'])} \\")
    if prm["analysis"] == "cyclic":
        w(rf"protocol & {tex(prm.get('proto'))}, {num(int(prm.get('cycles',1)))} cycle(s) per "
          rf"level, largest $\pm${num(float(prm['drift'])*100,3,PCT)} drift \\")
        pro = dat.get("protocol") or {}
        lv = pro.get("levels") if isinstance(pro, dict) else None
        if lv:
            w(r"levels (\% drift) & \val{" + tex(", ".join(f"{v*100:g}" for v in lv)) + r"} \\")
        if isinstance(pro, dict) and pro.get("drive_path_mm"):
            w(rf"drive path & {num(pro['drive_path_mm'],1,' mm')} \\")
    else:
        w(rf"target & {num(float(prm['drift'])*100,3,PCT)} drift \\")
    w(rf"solution scheme & dynamic relaxation, quasi-static \\")
    w(rf"integrator & \val{{{tex(integ)}}} \\")
    w(rf"time step & {num((dat.get('dt') or 0)*1e6,3,' us')}, "
      rf"{num(dat.get('steps_per_period_used'))} steps per fundamental period \\")
    w(rf"first period $T_1$ & {num((dat.get('T1') or 0)*1e3,3,' ms')} \\")
    w(rf"drive rate & {num(float(prm['rate']),2,' mm/s')} \\")
    w(rf"damping & {num(float(prm['damping']),2)} critical, mass proportional \\")
    w(rf"steps & {num(dat.get('steps'))} \\")
    w(rf"wall clock & {num((dat.get('elapsed_s') or 0)/3600,2,' h')} \\")
    w(rf"converged & {tex(dat.get('converged'))} \\")
    w(r"\end{factsheet}")

    # --- what the run produced ---------------------------------------------------------------------
    pk = dat.get("peak_shear_smooth") or dat.get("peak_shear")
    cap = dat.get("drift_capacity")
    w(r"\begin{factsheet}{Recorded response}")
    w(rf"peak base shear & {num((pk or 0)/1e3,1,' kN')} "
      rf"({num(dat.get('peak_smoothing_window_s',0)*1e3,1,' ms')} moving average) \\")
    w(rf"raw sample maximum & {num((dat.get('peak_shear_raw') or 0)/1e3,1,' kN')} \\")
    w(rf"drift at peak & {num((dat.get('drift_at_peak_smooth') or 0)*100,3,PCT)} \\")
    w(rf"drift capacity ({num(dat.get('capacity_drop',0.2)*100,0,PCT)} drop) & "
      rf"{(num(cap*100,3,PCT) if cap else 'not reached')}"
      rf"{' --- basis: ' + tex(dat['capacity_basis']) if dat.get('capacity_basis') else ''} \\")
    w(rf"traced to & {num((dat.get('end_drift') or 0)*100,3,PCT)} drift \\")
    w(rf"inertia + damping in the reaction, ascending branch (p95) & "
      rf"{num((dat.get('residual_ascending_p95') or 0)/1e3,1,' kN')} \\")
    w(r"\end{factsheet}")

    # --- figure --------------------------------------------------------------------------------
    fig = f"figures/{slug(d)}.pdf"
    w(r"\begin{figure}[H]\centering")
    w(rf"\includegraphics[width=\textwidth]{{{fig}}}")
    cap_txt = rs.figure_caption_cyclic if prm["analysis"] == "cyclic" else rs.figure_caption_push
    w(rf"\caption{{{cap_txt}}}")
    w(r"\end{figure}")
    return "\n".join(L) + "\n"


def _sort_key(spec: StudySpec, prm: dict, d: Path) -> tuple:
    """Sibling runs adjacent: model cell first, then what differs within a cell, then date."""
    rs = spec.runsheet
    extra = rs.sort_key(prm) if rs.sort_key else ()
    return (*spec.cell(prm), float(prm.get("mesh", spec.registry.by_name["mesh"].default)),
            *extra, float(prm.get("drift") or 0), d.name)


def main(spec: StudySpec, here: Path, argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="+", metavar="SUBSTR",
                    help="one or more substrings of the run-directory name")
    ap.add_argument("--list", action="store_true", help="report what would be generated")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")

    src = here / "src"
    dirs = discover(spec, a.only)
    if a.list:
        for d in dirs:
            print(d.name)
        print(f"{len(dirs)} run(s)")
        return 0
    if not dirs:
        print("no runs matched", file=sys.stderr)
        return 1

    (src / "runs").mkdir(parents=True, exist_ok=True)
    (src / "figures").mkdir(parents=True, exist_ok=True)
    order = []
    for d in dirs:
        prm, _filled = load_params(spec, d)
        dat = json.loads((d / "data.json").read_text())
        make_figure(spec, d, prm, dat, src / "figures" / f"{slug(d)}.pdf")
        order.append((prm["analysis"], _sort_key(spec, prm, d), slug(d), d))
        print(f"  {d.name}")

    # the manifest report.tex inputs; cyclic runs first, then pushovers, each by date
    lines = ["% Generated by generate.py — do not edit.", ""]
    for kind, label in (("cyclic", "Cyclic analyses"), ("pushover", "Pushover analyses")):
        sel = [o for o in order if o[0] == kind]
        if not sel:
            continue
        lines.append(rf"\part{{{label}}}")
        ordered = sorted(sel, key=lambda o: o[1])
        # Two runs in the same section whose contents entries would read alike get the run's own
        # timestamp appended, so an entry always identifies exactly one run.
        n_axes = len(spec.axes)
        entries = {}
        for _k, key, sl, dd in ordered:
            prm2, _f = load_params(spec, dd)
            dat2 = json.loads((dd / "data.json").read_text())
            entries[sl] = (key[:n_axes], run_titles(spec, prm2, dat2)[1])
        seen = {}
        for sl, (cell, short) in entries.items():
            seen.setdefault((cell, short), []).append(sl)
        dupes = {sl for group in seen.values() if len(group) > 1 for sl in group}

        prev, first_section = None, True
        for _k, key, sl, dd in ordered:
            cell = key[:n_axes]
            new_section = cell != prev
            if new_section:
                if not first_section:
                    lines.append(r"\clearpage")
                lines.append(rf"\section{{{spec.runsheet.section_title(cell)}}}")
                prev, first_section = cell, False
            st = dd.name[:16].replace("_", " ") if sl in dupes else ""
            (src / "runs" / f"{sl}.tex").write_text(
                sheet(spec, dd, first_in_section=new_section, stamp=st))
            lines.append(rf"\input{{runs/{sl}}}")
        lines.append("")
    (src / "runs" / "_manifest.tex").write_text("\n".join(lines) + "\n")
    print(f"{len(dirs)} run sheet(s) -> {src}")
    return 0
