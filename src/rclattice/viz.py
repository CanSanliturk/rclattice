"""matplotlib visualizer for lattice / continuum models (D-frame).

Pure drawing: works from a `Model` plus a displacement/mode-shape dict ({node_id: [u...]}).
It does NOT run analysis or import openseespy — callers build/run, then pass results here.
Uses the Agg backend so it works headless (save to file).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402

from .model import Model  # noqa: E402


def _arrays(model: Model):
    ids = sorted(model.nodes)
    idx = {nid: i for i, nid in enumerate(ids)}
    pts = np.array([model.nodes[nid].coords for nid in ids], dtype=float)
    lines = [(idx[e.nodes[0]], idx[e.nodes[1]]) for e in model.elements if len(e.nodes) == 2]
    quads = [[idx[n] for n in e.nodes] for e in model.elements if len(e.nodes) == 4]
    return ids, idx, pts, lines, quads


def _disp_array(ids, idx, disp) -> np.ndarray:
    d = np.zeros((len(ids), 2))
    if disp:
        for nid in ids:
            v = disp.get(nid)
            if v is not None:
                d[idx[nid]] = (v[0], v[1])
    return d


def autoscale_factor(model: Model, disp, frac: float = 0.12) -> float:
    """Scale so the largest nodal displacement is ~`frac` of the model's diagonal size."""
    ids, idx, pts, *_ = _arrays(model)
    d = _disp_array(ids, idx, disp)
    dmax = float(np.linalg.norm(d, axis=1).max()) if len(d) else 0.0
    if dmax <= 0.0:
        return 1.0
    size = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))
    return frac * size / dmax


def _segments(pts, lines):
    return [[pts[a], pts[b]] for a, b in lines]


def _verts(pts, quads):
    return [pts[q] for q in quads]


def draw_model(ax, model: Model, disp=None, scale: float = 1.0, *, undeformed=True,
               color="C0", title=None):
    """Draw a model's deformed shape (lattice struts as lines, continuum quads as polygons)."""
    ids, idx, pts, lines, quads = _arrays(model)
    d = _disp_array(ids, idx, disp)
    dpts = pts + scale * d

    if undeformed:
        if lines:
            ax.add_collection(LineCollection(_segments(pts, lines), colors="0.8", linewidths=0.5, zorder=1))
        if quads:
            ax.add_collection(PolyCollection(_verts(pts, quads), facecolors="none", edgecolors="0.85",
                                             linewidths=0.4, zorder=1))
    if lines:
        ax.add_collection(LineCollection(_segments(dpts, lines), colors=color, linewidths=1.1, zorder=3))
    if quads:
        ax.add_collection(PolyCollection(_verts(dpts, quads), facecolors=color, alpha=0.25,
                                         edgecolors=color, linewidths=0.4, zorder=2))
    ax.set_aspect("equal")
    ax.autoscale()
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=9)


# Per-kind styling for the undeformed analysis-model drawing. Reinforcement is drawn in distinct
# colours (longitudinal vs stirrup/tie), with lineweights kept light so they don't swamp the
# concrete skeleton in the combined view.
_KIND_STYLE = {
    "concrete":     {"color": "0.6",  "lw": 0.4, "zorder": 2, "label": "concrete struts"},
    "longitudinal": {"color": "C3",   "lw": 1.1, "zorder": 5, "label": "longitudinal rebar"},
    "stirrup":      {"color": "C2",   "lw": 0.9, "zorder": 4, "label": "stirrup / tie"},
    "rebar":        {"color": "C1",   "lw": 1.0, "zorder": 4, "label": "rebar"},
}
_KIND_DEFAULT = {"color": "C0", "lw": 0.8, "zorder": 3, "label": "strut"}

_CONCRETE_KINDS = {"concrete", "quad"}            # concrete skeleton (struts + continuum quads)
_REBAR_KINDS = {"longitudinal", "stirrup", "rebar"}  # reinforcement


def draw_model_kinds(ax, model: Model, *, which="all", title=None, legend=True, lim=None):
    """Draw the UNDEFORMED analysis model, styling each element by its `kind`: concrete struts thin
    and grey, continuum quads as light polygons, and reinforcement in distinct colours (longitudinal
    vs stirrup). `which` selects what to show: "concrete" (skeleton only), "rebar" (reinforcement
    only), or "all" (combined). `lim` is an optional ((x0,x1),(y0,y1)) to force a shared extent so
    the three panels of `figure_model` line up. Used for the `--draw` analysis-model figures."""
    ids, idx, pts, _lines, quads = _arrays(model)
    show_concrete = which in ("concrete", "all")
    show_rebar = which in ("rebar", "all")

    if quads and show_concrete:
        ax.add_collection(PolyCollection(_verts(pts, quads), facecolors="0.92", edgecolors="0.78",
                                         linewidths=0.4, zorder=1))

    groups: dict[str, list] = {}
    for e in model.elements:
        if len(e.nodes) != 2:
            continue
        if e.kind in _CONCRETE_KINDS and not show_concrete:
            continue
        if e.kind in _REBAR_KINDS and not show_rebar:
            continue
        groups.setdefault(e.kind, []).append([pts[idx[e.nodes[0]]], pts[idx[e.nodes[1]]]])

    handles, labels = [], []
    for kind in sorted(groups, key=lambda k: _KIND_STYLE.get(k, _KIND_DEFAULT)["zorder"]):
        st = _KIND_STYLE.get(kind, _KIND_DEFAULT)
        ax.add_collection(LineCollection(groups[kind], colors=st["color"], linewidths=st["lw"],
                                         zorder=st["zorder"]))
        handles.append(plt.Line2D([], [], color=st["color"], lw=max(st["lw"], 1.2)))
        labels.append(st["label"])

    ax.set_aspect("equal")
    if lim is not None:
        ax.set_xlim(*lim[0])
        ax.set_ylim(*lim[1])
    else:
        ax.autoscale()
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=9)
    if legend and handles:
        ax.legend(handles, labels, fontsize=7, loc="upper right", framealpha=0.85)


def _model_lim(model: Model, *, margin=0.05):
    """Shared ((x0,x1),(y0,y1)) extent for a model, with a small relative margin."""
    pts = np.array([n.coords for n in model.nodes.values()], dtype=float)
    (x0, y0), (x1, y1) = pts.min(axis=0), pts.max(axis=0)
    mx = margin * max(x1 - x0, 1e-9)
    my = margin * max(y1 - y0, 1e-9)
    return ((x0 - mx, x1 + mx), (y0 - my, y1 + my))


def figure_model(panels, *, savepath=None, suptitle="Analysis model", dpi=170):
    """Draw one or more undeformed analysis models, one ROW per model and three columns:
    concrete skeleton only, reinforcement only, and combined. `panels` is a list of (title, Model).

    Each panel is sized to the model's true aspect ratio (so a tall, narrow column is drawn tall and
    narrow rather than squeezed into a wide axes), which spreads the dense lattice out enough to
    inspect the struts. High DPI + thin concrete lines keep the skeleton legible even when crowded."""
    cols = (("concrete", "concrete lattice"), ("rebar", "reinforcement"), ("all", "combined"))
    n = len(panels)
    lims = [_model_lim(model) for _, model in panels]
    # per-panel size from the tallest aspect (h/w) across models; clamp width so squat models stay sane
    aspect = max(((ly[1] - ly[0]) / max(lx[1] - lx[0], 1e-9)) for lx, ly in lims)
    panel_h = 7.5
    panel_w = float(np.clip(panel_h / max(aspect, 1e-9), 1.7, 6.5))
    fig, axes = plt.subplots(n, 3, figsize=(3 * panel_w + 0.6, n * panel_h + 0.7), squeeze=False)
    for r, ((title, model), lim) in enumerate(zip(panels, lims)):
        for c, (which, col_label) in enumerate(cols):
            ax = axes[r][c]
            draw_model_kinds(ax, model, which=which, lim=lim, legend=(which != "concrete"))
            if r == 0:
                ax.set_title(col_label, fontsize=10)
        axes[r][0].set_ylabel(title, fontsize=10)
    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=dpi)
    return fig


def align_sign(reference: dict, shape: dict) -> dict:
    """Flip `shape`'s sign to best match `reference` (both keyed by the same node ids)."""
    keys = reference.keys() & shape.keys()
    dot = sum(
        reference[k][0] * shape[k][0] + reference[k][1] * shape[k][1] for k in keys
    )
    if dot < 0:
        return {k: [-c for c in v] for k, v in shape.items()}
    return shape


def sign_fix(shape: dict) -> dict:
    """Flip a mode shape so its largest-magnitude translational component is positive.

    A deterministic, model-agnostic sign (unlike `align_sign`, which needs shared node ids) — so the
    reference and the lattice mode shapes, which live on DIFFERENT meshes, both lean the same way."""
    best = 0.0
    for v in shape.values():
        for c in (v[0], v[1]):
            if abs(c) > abs(best):
                best = c
    if best < 0.0:
        return {k: [-c for c in v] for k, v in shape.items()}
    return shape


def figure_static(lattice, lat_disp, continuum, cont_disp, *, savepath=None):
    """Side-by-side static deformed shapes (shared scale from the continuum reference)."""
    scale = autoscale_factor(continuum, cont_disp)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    draw_model(axes[0], lattice, lat_disp, scale, color="C0",
               title=f"Lattice — static  (deform x{scale:.0f})")
    draw_model(axes[1], continuum, cont_disp, scale, color="C3",
               title=f"Continuum — static  (deform x{scale:.0f})")
    fig.suptitle("Static deflected shape: lattice vs continuum", fontsize=11)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=130)
    return fig


def figure_modes(panels, *, savepath=None):
    """Grid of mode shapes: rows = modes, cols = [lattice, continuum]. Each panel is a dict:
    {"lattice": (model, shape), "continuum": (model, shape), "T_lat": float, "T_cont": float}."""
    n = len(panels)
    fig, axes = plt.subplots(n, 2, figsize=(12, 3.4 * n), squeeze=False)
    for i, p in enumerate(panels):
        lm, ls = p["lattice"]
        cm, cs = p["continuum"]
        ls = align_sign(cs, ls)
        scale = autoscale_factor(cm, cs)
        draw_model(axes[i][0], lm, ls, scale, color="C0",
                   title=f"Lattice — mode {i + 1}   T = {p['T_lat']:.4f} s")
        draw_model(axes[i][1], cm, cs, scale, color="C3",
                   title=f"Continuum — mode {i + 1}   T = {p['T_cont']:.4f} s")
    fig.suptitle("Mode shapes: lattice vs continuum", fontsize=11)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=130)
    return fig


def figure_modal_calibration(reference, lattice, table_rows, *, caption="", savepath=None,
                             title="Modal calibration: reference vs lattice", dpi=150):
    """Calibration mode-shape report (D35): an N-column mode-shape grid with the reference model on the
    TOP row and the lattice on the BOTTOM row, and a periods table underneath.

    `reference` / `lattice` are dicts {"model": Model, "shapes": [shape per mode], "label": str,
    "color": str} (a mode shape is {node_id: [u...]}). `table_rows` is a list of
    (mode, T_ref, T_lattice, delta_pct) — one row per mode (its length sets the number of columns).
    `caption` is the one-line calibration summary (areas / RMS / K0-match), drawn under the table.

    Each mode shape is sign-fixed (`sign_fix`) independently — the reference and lattice are on
    different meshes, so they cannot be cross-aligned by node id — and auto-scaled to its own model."""
    n = len(table_rows)
    rows = [reference, lattice]
    fig = plt.figure(figsize=(max(3.6 * n + 0.6, 6.0), 9.0))
    gs = fig.add_gridspec(3, n, height_ratios=[1.0, 1.0, 0.5], hspace=0.16, wspace=0.05,
                          left=0.06, right=0.98, top=0.93, bottom=0.04)
    for r, panel in enumerate(rows):
        model, shapes = panel["model"], panel["shapes"]
        for c in range(n):
            ax = fig.add_subplot(gs[r, c])
            shape = sign_fix(shapes[c])
            draw_model(ax, model, shape, autoscale_factor(model, shape), color=panel["color"])
            if r == 0:
                ax.set_title(f"mode {c + 1}", fontsize=10)
            if c == 0:
                ax.set_ylabel(panel["label"], fontsize=10)

    ax_t = fig.add_subplot(gs[2, :])
    ax_t.axis("off")
    col_labels = ["mode", "T_ref (s)", "T_lattice (s)", "Δ vs reference (%)"]
    cell_text = [[str(int(m)), f"{tr:.4f}", f"{tl:.4f}", f"{d:+.2f}"] for m, tr, tl, d in table_rows]
    tbl = ax_t.table(cellText=cell_text, colLabels=col_labels, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.5)

    fig.suptitle(title, fontsize=12)
    if caption:
        fig.text(0.5, 0.012, caption, ha="center", fontsize=8.5, style="italic")
    if savepath:
        fig.savefig(savepath, dpi=dpi)
    return fig


def figure_pushover(curves, *, savepath=None, xlabel="roof displacement", ylabel="base shear",
                    title="Pushover: base shear vs roof displacement"):
    """Plot one or more pushover curves on shared axes.

    `curves` is a list of dicts: {"disp": [...], "shear": [...], "label": str, optional
    "style": dict of plot kwargs}. Pure drawing — caller supplies the recorded arrays.
    """
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    for c in curves:
        ax.plot(c["disp"], c["shear"], label=c.get("label"), **c.get("style", {}))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)
    if any(c.get("label") for c in curves):
        ax.legend()
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=130)
    return fig


def figure_timehistory(series, *, savepath=None, drift_label="roof drift (in)",
                       shear_label="base shear (kip)", title="Seismic time-history",
                       include_hysteresis: bool = True):
    """Seismic comparison panels: roof drift vs time, base shear vs time, and (optionally) the
    overlaid hysteresis loop. `series` is a list of dicts with "t", "disp", "shear", "label",
    optional "style". With `include_hysteresis=False` only the two history panels are drawn (the
    hysteresis is then reported separately by `figure_hysteresis`). Pure drawing."""
    ncol = 3 if include_hysteresis else 2
    fig, axes = plt.subplots(1, ncol, figsize=(5.3 * ncol, 4.6), squeeze=False)
    ax_d, ax_v = axes[0][0], axes[0][1]
    ax_h = axes[0][2] if include_hysteresis else None
    for s in series:
        style = s.get("style", {})
        ax_d.plot(s["t"], s["disp"], label=s.get("label"), **style)
        ax_v.plot(s["t"], s["shear"], label=s.get("label"), **style)
        if ax_h is not None:
            ax_h.plot(s["disp"], s["shear"], label=s.get("label"), **style)
    ax_d.set_xlabel("time (s)"); ax_d.set_ylabel(drift_label); ax_d.set_title("roof drift history", fontsize=10)
    ax_v.set_xlabel("time (s)"); ax_v.set_ylabel(shear_label); ax_v.set_title("base shear history", fontsize=10)
    used = [ax_d, ax_v]
    if ax_h is not None:
        ax_h.set_xlabel(drift_label); ax_h.set_ylabel(shear_label); ax_h.set_title("hysteresis (overlaid)", fontsize=10)
        used.append(ax_h)
    for ax in used:
        ax.grid(True, alpha=0.3)
        if any(s.get("label") for s in series):
            ax.legend(fontsize=9)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=130)
    return fig


def figure_hysteresis(series, *, savepath=None, drift_label="roof drift (in)",
                      shear_label="base shear (kip)", title="Hysteresis (base shear vs roof drift)"):
    """Per-model hysteresis loops plus an overlaid panel. `series` is a list of dicts with "disp",
    "shear", "label", optional "style". Draws one panel per series (each autoscaled on its own
    axes, so a small loop is not crushed by a large one) followed by an overlaid panel — so each
    model's loop is legible on its own AND the two are shown together. Pure drawing."""
    n = len(series)
    fig, axes = plt.subplots(1, n + 1, figsize=(5.0 * (n + 1), 4.6), squeeze=False)
    ax = axes[0]
    for i, s in enumerate(series):
        ax[i].plot(s["disp"], s["shear"], **s.get("style", {}))
        ax[i].set_title(s.get("label", f"series {i + 1}"), fontsize=10)
    for s in series:
        ax[n].plot(s["disp"], s["shear"], label=s.get("label"), **s.get("style", {}))
    ax[n].set_title("overlaid", fontsize=10)
    for a in ax:
        a.set_xlabel(drift_label); a.set_ylabel(shear_label)
        a.grid(True, alpha=0.3)
    if any(s.get("label") for s in series):
        ax[n].legend(fontsize=9)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=130)
    return fig


def animate_modes(panels, *, savepath, frames: int = 36, fps: int = 18, scale_frac: float = 0.14):
    """Animate (oscillate) the mode shapes in a rows=modes x cols=[lattice, continuum] grid,
    saved as a GIF. `panels` has the same structure as figure_modes."""
    n = len(panels)
    fig, axes = plt.subplots(n, 2, figsize=(12, 3.4 * n), squeeze=False)
    cells = []
    for i, p in enumerate(panels):
        for j, key, color in ((0, "lattice", "C0"), (1, "continuum", "C3")):
            model, shape = p[key]
            if key == "lattice":
                shape = align_sign(p["continuum"][1], shape)
            ids, idx, pts, lines, quads = _arrays(model)
            d = _disp_array(ids, idx, shape)
            scale = autoscale_factor(model, shape, scale_frac)
            ax = axes[i][j]
            # fixed limits that fit the full oscillation
            ext = pts + scale * d
            lo = np.minimum(pts.min(0), ext.min(0))
            hi = np.maximum(pts.max(0), ext.max(0))
            pad = 0.08 * (hi - lo + 1e-9)
            ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0])
            ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            tkey = "T_lat" if key == "lattice" else "T_cont"
            ax.set_title(f"{key.capitalize()} — mode {i + 1}  T={p[tkey]:.4f}s", fontsize=9)
            if lines:
                ax.add_collection(LineCollection(_segments(pts, lines), colors="0.85", linewidths=0.4, zorder=1))
            if quads:
                ax.add_collection(PolyCollection(_verts(pts, quads), facecolors="none",
                                                 edgecolors="0.88", linewidths=0.3, zorder=1))
            lc = LineCollection([], colors=color, linewidths=1.1, zorder=3) if lines else None
            pc = PolyCollection([], facecolors=color, alpha=0.25, edgecolors=color,
                                linewidths=0.4, zorder=2) if quads else None
            if lc is not None:
                ax.add_collection(lc)
            if pc is not None:
                ax.add_collection(pc)
            cells.append((pts, d, scale, lines, quads, lc, pc))

    def update(frame):
        amp = np.sin(2.0 * np.pi * frame / frames)
        for pts, d, scale, lines, quads, lc, pc in cells:
            dpts = pts + amp * scale * d
            if lc is not None:
                lc.set_segments(_segments(dpts, lines))
            if pc is not None:
                pc.set_verts(_verts(dpts, quads))
        return []

    fig.suptitle("Animated mode shapes: lattice vs continuum", fontsize=11)
    fig.tight_layout()
    anim = FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False)
    anim.save(savepath, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return savepath
