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
    merge_rel: float = 1e-6,
) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    """Mesh several axis-aligned rectangles (ox, oy, w, h) and merge coincident nodes.

    Each rectangle is meshed structured (mesh_rectangle_grid); nodes shared on touching edges are
    merged, so members connect at their joints. Returns (coords, quads) with quad indices into the
    merged coords. A single rectangle is the trivial case.

    THE MERGE TOLERANCE IS RELATIVE TO THE MESH, NOT ABSOLUTE, and that distinction is not cosmetic.
    Merging on `round(coords, 9)` fails whenever gmsh returns a shared edge's node as, say,
    -749.999999996 from one rectangle and -750.0 from the other: 4e-9 apart is a relative difference
    of 1e-12 on a 750 mm coordinate, i.e. ordinary double-precision noise, but it rounds to two
    different 9-decimal keys and the nodes survive as a coincident PAIR.

    That is not a harmless duplicate. `connect_horizon` then joins the pair with a strut of length
    ~4e-9 mm, whose EA/L is some ten orders of magnitude above every real strut's — a near-rigid
    link that wrecks the conditioning of the stiffness matrix while (kinematically) doing what the
    merge should have done anyway. Merging on `coords / (mesh_size * merge_rel)` keys the comparison
    to the grid instead: at mesh 50 the tolerance is 5e-5 mm, which absorbs the noise while still
    separating genuinely distinct nodes by a factor of ~1e6.
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
    key = np.round(coords / (mesh_size * merge_rel))
    _key, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
    unique = np.round(coords[first], decimals)
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
    min_rel: float = 1e-4,
) -> list[tuple[int, int]]:
    """Connect every node pair within `horizon * mesh_size` of each other (D9).

    Returns a deduplicated list of (i, j) index pairs with i < j (one element per pair).

    Pairs closer than `min_rel * mesh_size` are REJECTED rather than connected. A strut of near-zero
    length has an EA/L orders of magnitude above every real strut, so it behaves as a rigid link and
    destroys the conditioning of the assembly; such a pair always means two nodes that should have
    been merged and were not. Callers that mesh a compound domain should merge properly (see
    `mesh_compound_rectangles`) — this is the backstop, not the cure.
    """
    r = horizon * mesh_size + tol
    floor = min_rel * mesh_size
    n = len(coords)
    pairs: list[tuple[int, int]] = []
    for i in range(n):
        diff = coords[i + 1 :] - coords[i]
        dist = np.sqrt((diff * diff).sum(axis=1))
        for off, d in enumerate(dist):
            if floor < d <= r:
                pairs.append((i, i + 1 + off))
    return pairs


def connect_between(
    from_coords: np.ndarray,
    to_coords: np.ndarray,
    radius: float,
    *,
    tol: float = 1e-9,
    min_abs: float = 0.0,
) -> list[tuple[int, int]]:
    """Connect every node of `from_coords` to every node of `to_coords` within `radius` (D72).

    Unlike `connect_horizon` this joins TWO node sets, so there is no i < j dedup — each returned
    (i, j) indexes `from_coords` and `to_coords` respectively. Used for bond elements, which link a
    bar's own steel nodes to the surrounding concrete nodes.

    Pairs closer than `min_abs` are skipped. For bond that matters: a steel node normally sits
    exactly ON a concrete node (the mesh is aligned to the bar path), and a zero-length truss is not
    an element. Pass `min_abs` a small fraction of the grid spacing to drop the coincident partner
    and keep the surrounding ring.
    """
    r = radius + tol
    out: list[tuple[int, int]] = []
    for i, c in enumerate(from_coords):
        diff = to_coords - c
        dist = np.sqrt((diff * diff).sum(axis=1))
        for j, d in enumerate(dist):
            if min_abs < d <= r:
                out.append((i, j))
    return out


def perturb_nodes(
    coords: np.ndarray,
    mesh_size: float,
    rmax_ratio: float,
    *,
    seed: int | None = None,
    fix_nodes: "set[int] | None" = None,
    tol: float = 1e-9,
) -> np.ndarray:
    """Aydin's grid perturbation: move each node by `R <= rmax_ratio*mesh_size` at a random angle.

    From Aydin, Binici & Tuncay (2021), Fig. 2(b). `R` is uniform on `[0, Rmax]` and the angle
    `theta` uniform on `[0, 2*pi)`, both independent per node — NOT a uniform sample of the disc,
    which would bias `R` outward. `Rmax/d = rmax_ratio` is the model's single calibration parameter;
    the paper's fitted values run 0.04-0.08 (Table 1).

    Why it exists: on a structured grid the vertical struts form a straight, perfectly aligned load
    path from platen to platen. Once the transverse struts split, that path carries the whole load
    with no kink to destabilize it — Aydin's "vertical locking", which makes a uniform-grid lattice
    read roughly twice its concrete's compressive strength (his Fig. 3(a)) or, with a compression
    cap on the struts, exactly the vertical struts' share of the face. Perturbing the grid replaces
    the straight columns with zigzag ones, so every kink demands transverse equilibrium from struts
    that have already cracked, and the column loses stability instead of locking.

    BOUNDARY HANDLING (a documented deviation — the paper does not say what it does). Nodes on an
    edge move only ALONG that edge, using the tangential component of the same random vector, and
    corner nodes do not move at all. Every boundary node therefore stays exactly on its face, so
    box selection of supports, driven nodes and reaction nodes remains exact, and the loaded faces
    are still flat. The alternative — perturbing boundary nodes off their face — would make the
    platen contact ragged and the specimen's own dimensions random.

    `fix_nodes` are additional row indices left exactly in place (e.g. a single node carrying a
    stability restraint that a box query must find). Returns a NEW array; `coords` is not modified.
    """
    c = np.array(coords, dtype=float, copy=True)
    if rmax_ratio <= 0.0:
        return c

    rng = np.random.default_rng(seed)
    rmax = rmax_ratio * mesh_size
    r = rng.uniform(0.0, rmax, size=len(c))
    theta = rng.uniform(0.0, 2.0 * np.pi, size=len(c))
    delta = np.column_stack((r * np.cos(theta), r * np.sin(theta)))

    x, y = c[:, 0], c[:, 1]
    xlo, xhi, ylo, yhi = x.min(), x.max(), y.min(), y.max()
    on_x = np.isclose(x, xlo, atol=tol) | np.isclose(x, xhi, atol=tol)   # left/right edge
    on_y = np.isclose(y, ylo, atol=tol) | np.isclose(y, yhi, atol=tol)   # bottom/top edge

    delta[on_x, 0] = 0.0          # left/right edge nodes slide vertically only
    delta[on_y, 1] = 0.0          # bottom/top edge nodes slide horizontally only
    delta[on_x & on_y] = 0.0      # corners are pinned
    if fix_nodes:
        delta[list(fix_nodes)] = 0.0

    out = c + delta
    # Keep every node inside the original bounding box: an edge node sliding past a corner would
    # lengthen the specimen. Rmax is a few per cent of the spacing, so this almost never bites.
    np.clip(out[:, 0], xlo, xhi, out=out[:, 0])
    np.clip(out[:, 1], ylo, yhi, out=out[:, 1])
    return out
