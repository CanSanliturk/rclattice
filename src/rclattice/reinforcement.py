"""Reinforcement geometry (D13): map a rebar polyline onto the lattice node set.

Backend-agnostic. A `Rebar.path` is a polyline of (x, y) points; in the lattice the bar becomes
steel truss struts connecting the lattice nodes that lie on the path, in order (shared nodes =>
perfect bond, D5). The mesh is expected to align so nodes fall on the path; `tol` snaps nodes
that are within a small distance of it.
"""

from __future__ import annotations

import numpy as np


def _nodes_on_segment(coords: np.ndarray, p0, p1, tol: float) -> list[tuple[float, int]]:
    """(param t, node index) for nodes within `tol` of segment p0->p1, ordered along it."""
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    d = p1 - p0
    L2 = float(d @ d)
    if L2 == 0.0:
        return []
    out: list[tuple[float, int]] = []
    for i, c in enumerate(coords):
        w = c - p0
        t = float(w @ d) / L2
        if -1e-9 <= t <= 1.0 + 1e-9:
            proj = p0 + t * d
            if float(np.linalg.norm(c - proj)) <= tol:
                out.append((max(0.0, min(1.0, t)), i))
    out.sort()
    return out


def rebar_node_chain(coords: np.ndarray, path, tol: float = 1e-6) -> list[int]:
    """Ordered, de-duplicated lattice node indices (0-based) lying on a rebar polyline `path`.

    Walks each segment, collecting on-path nodes in order; consecutive duplicates (shared
    segment endpoints) are merged. Raises if a segment finds no nodes (mesh not aligned).
    """
    chain: list[int] = []
    for a, b in zip(path, path[1:]):
        seg = _nodes_on_segment(coords, a, b, tol)
        if not seg:
            raise ValueError(f"rebar segment {a}->{b} matched no lattice nodes (tol={tol}); "
                             "mesh is not aligned to the rebar path")
        for _t, i in seg:
            if not chain or chain[-1] != i:
                chain.append(i)
    seen: set[int] = set()
    ordered: list[int] = []
    for i in chain:
        if i not in seen:
            seen.add(i)
            ordered.append(i)
    return ordered
