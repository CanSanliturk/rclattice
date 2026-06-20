"""Lattice/continuum meshing: gmsh places nodes (D6/D10); a horizon rule builds struts (D9).

`mesh_rectangle_grid` returns BOTH node coordinates and quad connectivity from one structured
gmsh mesh, so the lattice (nodes + horizon struts) and the continuum (nodes + quads) share an
identical node set — making the verification comparison fair.

Backend-agnostic w.r.t. OpenSees — this module never imports openseespy.
"""

from __future__ import annotations

import gmsh
import numpy as np

_GMSH_QUAD4 = 3  # gmsh element type id for a 4-node quadrangle


def mesh_rectangle_grid(
    length: float,
    height: float,
    mesh_size: float,
    *,
    origin: tuple[float, float] = (0.0, 0.0),
    decimals: int = 9,
) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    """Structured (regular grid) mesh of a rectangle via gmsh transfinite recombination.

    Returns (coords, quads):
      - coords: (N, 2) array of node coordinates (row order is the node index).
      - quads:  list of 4-tuples of node indices (CCW), referring to rows of coords.
    """
    ox, oy = origin
    nx = max(1, round(length / mesh_size))
    ny = max(1, round(height / mesh_size))

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("rect")
        g = gmsh.model.geo

        p = [
            g.addPoint(ox, oy, 0.0),
            g.addPoint(ox + length, oy, 0.0),
            g.addPoint(ox + length, oy + height, 0.0),
            g.addPoint(ox, oy + height, 0.0),
        ]
        lines = [g.addLine(p[k], p[(k + 1) % 4]) for k in range(4)]
        loop = g.addCurveLoop(lines)
        surf = g.addPlaneSurface([loop])

        g.mesh.setTransfiniteCurve(lines[0], nx + 1)
        g.mesh.setTransfiniteCurve(lines[2], nx + 1)
        g.mesh.setTransfiniteCurve(lines[1], ny + 1)
        g.mesh.setTransfiniteCurve(lines[3], ny + 1)
        g.mesh.setTransfiniteSurface(surf)
        g.mesh.setRecombine(2, surf)
        g.synchronize()
        gmsh.model.mesh.generate(2)

        node_tags, coord_flat, _ = gmsh.model.mesh.getNodes()
        coords = np.round(np.array(coord_flat, dtype=float).reshape(-1, 3)[:, :2], decimals)
        tag2idx = {int(t): i for i, t in enumerate(node_tags)}

        etypes, _etags, enodes = gmsh.model.mesh.getElements(2, surf)
        quads: list[tuple[int, int, int, int]] = []
        for et, en in zip(etypes, enodes):
            if et == _GMSH_QUAD4:
                rows = np.array(en, dtype=int).reshape(-1, 4)
                for row in rows:
                    q = tuple(tag2idx[int(t)] for t in row)
                    quads.append(_ensure_ccw(coords, q))
    finally:
        gmsh.finalize()

    return coords, quads


def mesh_rectangle_nodes(length: float, height: float, mesh_size: float, **kw) -> np.ndarray:
    """Convenience: just the node coordinates of the structured rectangle grid."""
    coords, _quads = mesh_rectangle_grid(length, height, mesh_size, **kw)
    return coords


def mesh_compound_rectangles(
    rects: list[tuple[float, float, float, float]],
    mesh_size: float,
    *,
    decimals: int = 9,
) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    """Mesh several axis-aligned rectangles (ox, oy, w, h) and merge coincident nodes.

    Each rectangle is meshed structured (mesh_rectangle_grid); nodes shared on touching edges
    are merged by rounded coordinate, so members connect at their joints. Returns (coords,
    quads) with quad indices into the merged coords. A single rectangle is the trivial case.
    """
    parts: list[np.ndarray] = []
    quads: list[tuple[int, int, int, int]] = []
    offset = 0
    for ox, oy, w, h in rects:
        c, q = mesh_rectangle_grid(w, h, mesh_size, origin=(ox, oy), decimals=decimals)
        parts.append(c)
        quads.extend(tuple(int(i) + offset for i in quad) for quad in q)
        offset += len(c)

    coords = np.vstack(parts)
    unique, inverse = np.unique(np.round(coords, decimals), axis=0, return_inverse=True)
    inverse = inverse.ravel()
    merged_quads = [tuple(int(inverse[i]) for i in quad) for quad in quads]
    return unique, merged_quads


def _ensure_ccw(coords: np.ndarray, q: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    """Return the quad's node indices in counter-clockwise order (OpenSees expects CCW)."""
    pts = coords[list(q)]
    area2 = 0.0
    for i in range(4):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % 4]
        area2 += x1 * y2 - x2 * y1
    return q if area2 > 0 else (q[0], q[3], q[2], q[1])


def connect_horizon(
    coords: np.ndarray,
    mesh_size: float,
    horizon: float = 1.5,
    *,
    tol: float = 1e-9,
) -> list[tuple[int, int]]:
    """Connect every node pair within `horizon * mesh_size` of each other (D9).

    Returns a deduplicated list of (i, j) index pairs with i < j (one element per pair).
    """
    r = horizon * mesh_size + tol
    n = len(coords)
    pairs: list[tuple[int, int]] = []
    for i in range(n):
        diff = coords[i + 1 :] - coords[i]
        dist = np.sqrt((diff * diff).sum(axis=1))
        for off, d in enumerate(dist):
            if 0.0 < d <= r:
                pairs.append((i, i + 1 + off))
    return pairs
