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


def align_sign(reference: dict, shape: dict) -> dict:
    """Flip `shape`'s sign to best match `reference` (both keyed by the same node ids)."""
    keys = reference.keys() & shape.keys()
    dot = sum(
        reference[k][0] * shape[k][0] + reference[k][1] * shape[k][1] for k in keys
    )
    if dot < 0:
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
