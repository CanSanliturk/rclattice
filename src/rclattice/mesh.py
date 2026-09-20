"""Lattice/continuum meshing: gmsh places nodes (D6/D10); a horizon rule builds struts (D9).

`mesh_rectangle_grid` returns BOTH node coordinates and quad connectivity from one structured
gmsh mesh, so the lattice (nodes + horizon struts) and the continuum (nodes + quads) share an
identical node set — making the verification comparison fair.

Backend-agnostic w.r.t. OpenSees — this module never imports openseespy.
"""

from __future__ import annotations

import math

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


# --- rebar-aligned (graded) grids, D104 ---------------------------------------------------------
#
# The uniform grid forces every bar axis onto a multiple of the mesh, which the Thomsen & Wallace
# RW2 layout (bars at 19, 70, 121, 172 mm; web bars at 323.5 + n*191) cannot satisfy at any usable
# spacing. The alternative is a STRUCTURED grid whose lines are the union of the bar axes and a
# regular fill: every interval between consecutive hard lines is split into `round(interval/mesh)`
# equal parts, so bars sit on nodes exactly and the spacing stays close to the target everywhere.
# gmsh remains the single node source (D6/D10): the boundary is split at the hard coordinates, each
# piece is transfinite with its own count, and one transfinite surface with the four corners named
# yields the tensor grid.


def graded_lines(extent: float, mesh_size: float, hard: "tuple[float, ...] | list[float]" = (),
                 *, origin: float = 0.0, tol: float = 1e-9) -> list[float]:
    """Grid-line coordinates across `extent`: the hard lines, each gap filled at ~`mesh_size`.

    Every interval between consecutive hard lines (the two ends included) is divided into
    `max(1, round(interval / mesh_size))` EQUAL parts, so the local spacing is the closest possible
    to the target while every hard line is honoured exactly. Duplicates and lines outside the
    extent are dropped.
    """
    pts = sorted({origin, origin + extent, *(float(h) for h in hard
                                              if -tol <= float(h) - origin <= extent + tol)})
    out = [pts[0]]
    for a, b in zip(pts, pts[1:]):
        n = max(1, int(round((b - a) / mesh_size)))
        out += [a + (b - a) * k / n for k in range(1, n + 1)]
    return out


def mesh_rectangle_lines(
    length: float,
    height: float,
    mesh_size: float,
    *,
    x_lines: "tuple[float, ...] | list[float]" = (),
    y_lines: "tuple[float, ...] | list[float]" = (),
    origin: tuple[float, float] = (0.0, 0.0),
    decimals: int = 9,
) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    """Structured gmsh mesh of a rectangle on the graded lines of `graded_lines` in x and y.

    Same return contract as `mesh_rectangle_grid`: `(coords, quads)`, quads CCW into `coords`.
    With no hard lines and an `extent` that `mesh_size` divides it reproduces the uniform grid.
    """
    ox, oy = origin
    xs = graded_lines(length, mesh_size, x_lines, origin=ox)
    ys = graded_lines(height, mesh_size, y_lines, origin=oy)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("graded")
        g = gmsh.model.geo
        # boundary points, counter-clockwise from the origin corner
        bottom = [g.addPoint(x, ys[0], 0.0) for x in xs]
        right = [bottom[-1]] + [g.addPoint(xs[-1], y, 0.0) for y in ys[1:]]
        top = [right[-1]] + [g.addPoint(x, ys[-1], 0.0) for x in xs[-2::-1]]
        left = [top[-1]] + [g.addPoint(xs[0], y, 0.0) for y in ys[-2:0:-1]] + [bottom[0]]
        loop_pts = bottom + right[1:] + top[1:] + left[1:-1]
        lines = []
        for a, b in zip(loop_pts, loop_pts[1:] + loop_pts[:1]):
            ln = g.addLine(a, b)
            g.mesh.setTransfiniteCurve(ln, 2)        # every boundary piece is ONE cell long
            lines.append(ln)
        loop = g.addCurveLoop(lines)
        surf = g.addPlaneSurface([loop])
        g.mesh.setTransfiniteSurface(surf, cornerTags=[bottom[0], bottom[-1], top[0], left[0]])
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
                for row in np.array(en, dtype=int).reshape(-1, 4):
                    quads.append(_ensure_ccw(coords, tuple(tag2idx[int(t)] for t in row)))
    finally:
        gmsh.finalize()

    # Transfinite interpolation on a rectangle with straight sides puts every interior node on the
    # tensor grid, but float noise can leave it a few 1e-12 off; snap to the exact lines so a bar
    # path at a hard coordinate matches its nodes at the reinforcement tolerance.
    xs_a, ys_a = np.asarray(xs), np.asarray(ys)
    coords[:, 0] = xs_a[np.abs(coords[:, 0][:, None] - xs_a[None, :]).argmin(axis=1)]
    coords[:, 1] = ys_a[np.abs(coords[:, 1][:, None] - ys_a[None, :]).argmin(axis=1)]
    return coords, quads


def grid_indices(coords: np.ndarray, *, tol: float = 1e-6) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """`(ix, iy, xs, ys)`: each node's column/row index on the structured lines it sits on."""
    xs = np.unique(np.round(coords[:, 0] / tol)) * tol
    ys = np.unique(np.round(coords[:, 1] / tol)) * tol
    ix = np.abs(coords[:, 0][:, None] - xs[None, :]).argmin(axis=1)
    iy = np.abs(coords[:, 1][:, None] - ys[None, :]).argmin(axis=1)
    return ix, iy, xs, ys


def connect_index_horizon(coords: np.ndarray, horizon: float = 1.5) -> list[tuple[int, int]]:
    """The horizon rule in INDEX space: connect (i, j) to (i+di, j+dj) whenever di^2+dj^2 <= h^2.

    On a uniform grid this is exactly `connect_horizon` (D9): at 1.5 the orthogonal and diagonal
    neighbours, at 3.01 the next rings too. On a GRADED grid the physical rule would wire a short
    interval's nodes to their second neighbours and skip a long interval's diagonals, making the
    topology a function of the grading; the index rule keeps it the regular-grid topology. Returns
    (i, j) pairs with i < j, one element per pair.
    """
    ix, iy, xs, ys = grid_indices(coords)
    lookup = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(ix, iy))}
    r = int(math.floor(horizon))
    offsets = [(di, dj) for di in range(-r, r + 1) for dj in range(0, r + 1)
               if (di, dj) != (0, 0) and (dj > 0 or di > 0) and di * di + dj * dj <= horizon * horizon + 1e-9]
    pairs: list[tuple[int, int]] = []
    for k, (a, b) in enumerate(zip(ix, iy)):
        for di, dj in offsets:
            m = lookup.get((int(a) + di, int(b) + dj))
            if m is not None:
                pairs.append((min(k, m), max(k, m)))
    return sorted(set(pairs))


def tributary_area_scale(coords: np.ndarray, pairs, mesh_size: float) -> np.ndarray:
    """Per-strut factor turning a uniform-grid area into the area a GRADED grid needs.

    The energy balance returns `A_t` for a grid of spacing `d`, and Aydin's closed form makes the
    dependence explicit: `E_t A_t = C E_t d w`, i.e. a strut's area is proportional to the grid
    spacing it represents. A strut's axial stiffness `EA/L` stands in for a strip of the continuum
    whose WIDTH is the perpendicular spacing and whose length is the strut's own, so on a graded
    grid an orthogonal strut takes `A_t * d_perp / d` with `d_perp` the mean spacing of the
    neighbouring lines in the perpendicular direction, and a diagonal takes
    `A_t * sqrt(dx dy) / d`, the isotropic scale of the cell it spans. Under a uniform rescaling
    of the whole grid both reduce to the scale factor, as they must; on a mildly graded grid they
    keep `EA/L` locally matched to the continuum, which the elastic gate then verifies.
    """
    ix, iy, xs, ys = grid_indices(coords)
    dx = np.diff(xs)
    dy = np.diff(ys)

    def spacing_at(lines_d, i, i2):
        """Mean line spacing over the index span [i, i2] (perpendicular width of a strut)."""
        lo, hi = min(i, i2), max(i, i2)
        # spacing "owned" by a node is the mean of the intervals either side of it
        left = lines_d[max(lo - 1, 0)] if lo > 0 else lines_d[0]
        right = lines_d[min(hi, len(lines_d) - 1)] if hi < len(lines_d) else lines_d[-1]
        return 0.5 * (left + right)

    out = np.empty(len(pairs))
    for k, (i, j) in enumerate(pairs):
        di, dj = int(ix[j]) - int(ix[i]), int(iy[j]) - int(iy[i])
        if dj == 0:                                   # horizontal: width is the y-spacing
            d_eff = spacing_at(dy, iy[i], iy[i])
        elif di == 0:                                 # vertical: width is the x-spacing
            d_eff = spacing_at(dx, ix[i], ix[i])
        else:                                         # inclined: the cell's isotropic scale
            cx = abs(xs[ix[j]] - xs[ix[i]]) / abs(di)
            cy = abs(ys[iy[j]] - ys[iy[i]]) / abs(dj)
            d_eff = math.sqrt(cx * cy)
        out[k] = d_eff / mesh_size
    return out
