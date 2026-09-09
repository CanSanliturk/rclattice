"""Builders (D12): translate a backend-agnostic `Problem` into an OpenSees FE `Model`.

Both builders mesh the SAME structured node grid (mesh_rectangle_grid), so the lattice and
the continuum share an identical node set and the verification comparison is fair. Each
builder also returns the edge->node-id sets so callers can query/compare results.

Builders are backend-agnostic: they produce a `Model`; they never import openseespy.
"""

from __future__ import annotations

from typing import Callable

import math

import numpy as np

from .materials import concrete_nd_elastic, concrete_uniaxial_elastic, steel_uniaxial
from .mesh import connect_between, connect_horizon, mesh_compound_rectangles
from .model import Load, Model, NDMaterial, Support, UniaxialMaterial
from .problem import BoxLoad, BoxSupport, EdgeLoad, EdgeSupport, Problem, Rebar
from .reinforcement import rebar_node_chain

EdgeNodes = dict[str, list[int]]


def build_lattice(
    problem: Problem,
    mesh_size: float,
    *,
    horizon: float = 1.5,
    strut_area: "float | Callable[[float], float]" = 1.0,
) -> tuple[Model, EdgeNodes]:
    """Lattice builder: gmsh nodes + horizon struts; uniaxial-elastic concrete struts.

    `strut_area` is a uniform float, or a callable mapping strut length -> area (used by the
    area-group calibration, D16: e.g. orthogonal vs diagonal struts).
    """
    coords, quads = _grid(problem, mesh_size)
    area_fn: Callable[[float], float] = strut_area if callable(strut_area) else (lambda _L: strut_area)

    model = Model(ndm=2, ndf=2)
    mat = concrete_uniaxial_elastic(problem.material, 1)
    model.uniaxial_materials.append(mat)

    for idx, (x, y) in enumerate(coords, start=1):
        model.add_node(idx, (float(x), float(y)))
    for eid, (i, j) in enumerate(connect_horizon(coords, mesh_size, horizon), start=1):
        length = float(np.linalg.norm(coords[j] - coords[i]))
        model.add_element(eid, "Truss", (i + 1, j + 1), (area_fn(length), mat.id))

    _assign_tributary_mass(model, coords, quads, problem)
    edges = _edges(coords)
    _apply_supports_loads(model, problem, coords, edges)
    return model, edges


def build_lattice_rc(
    problem: Problem,
    mesh_size: float,
    *,
    material_for: "Callable[[str, float], UniaxialMaterial]",
    zone_of: "Callable[[float, float], str]",
    rebars: "tuple[Rebar, ...]" = (),
    horizon: float = 1.5,
    strut_area: "float | Callable[[float], float]" = 1.0,
    rebar_tol: float = 1e-6,
    strut_element: str = "Truss",
    rebar_material: "Callable[[object, int], UniaxialMaterial]" = steel_uniaxial,
    grid: "tuple[np.ndarray, list] | None" = None,
    pairs: "list[tuple[int, int]] | None" = None,
    bond_material: "Callable[[str, float], UniaxialMaterial] | None" = None,
    bond_area: "float | Callable[[float], float] | None" = None,
    bond_horizon: "float | None" = None,
    bond_mass_share: float = 0.1,
    i_accept_the_known_bond_defect: bool = False,
) -> tuple[Model, EdgeNodes]:
    """RC lattice builder (D19/D20, Stage 2): per-zone, length-regularized concrete struts +
    steel rebar struts.

    `strut_element` is the OpenSees element for ALL struts (concrete + rebar): "Truss" (default,
    small-displacement) or "corotTruss" for geometric consistency at large drift (carries the
    P-Delta effect the beam-column reference shows; D22).

    `rebar_material(steel_grade, tag)` maps each rebar's steel grade to its OpenSees uniaxial
    material; defaults to `steel_uniaxial` (Steel02). Pass `steel_uniaxial_elastic` to keep the
    rebar linear-elastic for the linear-material verification (same topology, elastic constitutive).

    Each horizon strut is assigned a concrete material by `material_for(zone, length)`, where the
    zone comes from `zone_of(x, y)` at the strut midpoint and `length` is the strut length —
    enabling fracture-energy regularization that depends on strut length (D20). Materials are
    cached/deduplicated by (zone, rounded length), so the regular grid yields just a couple of
    materials per zone (orthogonal + diagonal). Each `Rebar` becomes steel truss struts on the
    lattice nodes lying on its path (perfect bond, D13); one Steel02 material per distinct steel
    grade. Same grid / mass / supports / loads plumbing as `build_lattice`.

    `grid` and `pairs` override the generated geometry, for a lattice whose node positions are not
    the plain structured grid (D60, Aydin's perturbed mesh). `grid` is the `(coords, quads)` pair
    `mesh_rectangle_grid` would return — quads are used only for tributary mass, so they may keep
    the unperturbed connectivity while `coords` carry the moved positions. `pairs` overrides the
    horizon connectivity; pass it whenever `coords` are perturbed, so the strut TOPOLOGY still comes
    from the regular grid. Re-running `connect_horizon` on moved nodes silently drops the struts
    whose length happens to cross `horizon*mesh_size` (5 of 420 at Rmax/d = 0.08), which would make
    the lattice's connectivity a random variable on top of its geometry.

    BOND (D72, optional — default OFF, i.e. perfect bond as before).

    **KNOWN DEFECT — DO NOT USE FOR RESULTS.** The horizon-star topology below makes the assembly
    STIFFER than perfect bond, which is physically impossible: an elastic Aldemir wall comes out at
    2.19x its transformed-section stiffness with a 12.4% secant spread in a linear model. Each bond
    link carries the concrete strut's EA, so the ring around every steel node is a stiff parallel
    load path between the surrounding concrete nodes, ADDED to a concrete lattice that is still
    fully present. It fails the defining limit test — as bond stiffness goes to infinity the model
    must reduce to perfect bond — and the star cannot satisfy that limit at any stiffness. The fix is
    a coincident-node interface spring (`zeroLength`: bond law along the bar axis, stiff transverse),
    not a smaller `bond_area`. Passing `bond_material` therefore raises unless
    `i_accept_the_known_bond_defect=True`.

    Passing `bond_material` switches the reinforcement from shared nodes to Aydin's (2019) scheme:

      * each `Rebar` gets its OWN steel nodes, duplicated at the coordinates of the lattice nodes on
        its path, so steel and concrete no longer share a DOF;
      * steel struts run between consecutive steel nodes as before;
      * BOND elements join every steel node to the concrete nodes within `bond_horizon * mesh_size`,
        EXCLUDING the concrete node it sits on (that link is zero-length). At the default horizon
        that is the surrounding ring of 8.

    `bond_material(zone, length)` mirrors `material_for` — see `materials.bond_elastic_brittle`.
    `bond_area` defaults to `strut_area`, and `bond_horizon` to `horizon`; the paper notes bond
    stiffens AND strengthens as the horizon grows, because the links get both more numerous and
    longer, so these are separated from the concrete strut settings deliberately.

    Bars that CROSS get separate steel nodes at the crossing and are coupled only through the
    concrete — which is what the physical bars do, sitting in different layers through the thickness.

    `bond_mass_share` is a numerical device, not a physical steel mass: each steel node takes that
    fraction of its coincident concrete node's tributary mass, and the concrete node keeps the rest,
    so the TOTAL is conserved exactly. It exists because a steel node has no tributary area of its
    own, and a zero-mass DOF makes the dynamic runners singular.
    """
    if bond_material is not None and not i_accept_the_known_bond_defect:
        raise ValueError(
            "bond_material is set, but the horizon-star bond scheme has a KNOWN DEFECT (D72): it "
            "makes the model STIFFER than perfect bond (2.19x the transformed section on the "
            "Aldemir wall) because each link carries the concrete strut's EA, forming a parallel "
            "load path beside the intact concrete lattice. It cannot satisfy the limit test "
            "'rigid bond == perfect bond' at any stiffness; the fix is a coincident-node zeroLength "
            "interface spring, not a smaller bond_area. Pass "
            "i_accept_the_known_bond_defect=True only to reproduce or debug the defect.")
    coords, quads = _grid(problem, mesh_size) if grid is None else grid
    if pairs is None:
        pairs = connect_horizon(coords, mesh_size, horizon)
    area_fn: Callable[[float], float] = strut_area if callable(strut_area) else (lambda _L: strut_area)

    model = Model(ndm=2, ndf=2)
    tag = 1

    for idx, (x, y) in enumerate(coords, start=1):
        model.add_node(idx, (float(x), float(y)))

    eid = 1
    mat_cache: dict[tuple[str, float], int] = {}  # (zone, rounded length) -> material tag
    for (i, j) in pairs:
        length = float(np.linalg.norm(coords[j] - coords[i]))
        mx, my = 0.5 * (coords[i] + coords[j])
        zone = zone_of(float(mx), float(my))
        key = (zone, round(length, 6))
        if key not in mat_cache:
            mat = material_for(zone, length)
            mat.id = tag
            model.uniaxial_materials.append(mat)
            mat_cache[key] = tag
            tag += 1
        model.add_element(eid, strut_element, (i + 1, j + 1), (area_fn(length), mat_cache[key]))
        eid += 1

    steel_tag: dict[int, int] = {}  # id(SteelGrade) -> material tag (one material per grade)
    bond_mat_cache: dict[tuple[str, float], int] = {}
    bond_area_fn: Callable[[float], float] = area_fn
    if bond_area is not None:
        bond_area_fn = bond_area if callable(bond_area) else (lambda _L: bond_area)
    bond_radius = (horizon if bond_horizon is None else bond_horizon) * mesh_size
    steel_of_concrete: dict[int, list[int]] = {}   # concrete idx -> steel node ids sitting on it
    next_node = len(coords) + 1

    for rb in rebars:
        gid = id(rb.steel)
        if gid not in steel_tag:
            # A factory may return ONE material or a SEQUENCE of them. The sequence form exists so
            # a wrapper can be built on top of a base law — `MinMax` giving `Steel02` a rupture
            # strain, say — where the wrapper must reference the base material's own tag. The LAST
            # entry is what elements use; earlier ones are its dependencies. A factory returning a
            # single material behaves exactly as before.
            made = rebar_material(rb.steel, tag)
            mats = list(made) if isinstance(made, (list, tuple)) else [made]
            model.uniaxial_materials.extend(mats)
            steel_tag[gid] = mats[-1].id
            tag = max(m.id for m in mats) + 1
        chain = rebar_node_chain(coords, rb.path, rebar_tol)

        if bond_material is None:
            # Perfect bond (D13): the bar's struts run on the concrete nodes themselves.
            bar_nodes = [c + 1 for c in chain]
        else:
            # Aydin (2019): the bar gets its own nodes, coincident but independent, and is tied
            # back to the concrete through bond elements (D72).
            bar_nodes = []
            for c in chain:
                model.add_node(next_node, (float(coords[c][0]), float(coords[c][1])))
                steel_of_concrete.setdefault(int(c), []).append(next_node)
                bar_nodes.append(next_node)
                next_node += 1

        for a, b in zip(bar_nodes, bar_nodes[1:]):
            model.add_element(eid, strut_element, (a, b), (rb.area, steel_tag[gid]), kind=rb.role)
            eid += 1

        if bond_material is None:
            continue
        # `min_abs` drops the coincident concrete node (a zero-length truss is not an element) while
        # keeping the surrounding ring; the floor mirrors `connect_horizon`'s own near-zero guard.
        links = connect_between(coords[chain], coords, bond_radius, min_abs=1e-4 * mesh_size)
        for k, cj in links:
            si = bar_nodes[k]
            length = float(np.linalg.norm(coords[cj] - coords[chain[k]]))
            mx, my = 0.5 * (coords[chain[k]] + coords[cj])
            zone = zone_of(float(mx), float(my))
            key = (zone, round(length, 6))
            if key not in bond_mat_cache:
                # A bond factory may return SEVERAL materials — a wrapper must reference its base
                # by tag, so a damaging law is a small stack (D95). Elements use the LAST, exactly
                # as they do for a wrapped rebar material.
                bmats = bond_material(zone, length)
                bmats = list(bmats) if isinstance(bmats, (list, tuple)) else [bmats]
                # Remap the factory's OWN tags onto the model's numbering. Build the whole map
                # first and rewrite args before touching any id: a wrapper references its base by
                # tag, so mutating ids while still reading them corrupts the references.
                id_map = {m.id: tag + i for i, m in enumerate(bmats)}
                for bmat in bmats:
                    bmat.args = tuple(id_map[v] if isinstance(v, int) and v in id_map else v
                                      for v in bmat.args)
                for bmat in bmats:
                    bmat.id = id_map[bmat.id]
                    model.uniaxial_materials.append(bmat)
                bond_mat_cache[key] = bmats[-1].id      # elements use the LAST
                tag += len(bmats)
            model.add_element(eid, strut_element, (si, cj + 1),
                              (bond_area_fn(length), bond_mat_cache[key]), kind="bond")
            eid += 1

    _assign_tributary_mass(model, coords, quads, problem)
    _share_mass_with_steel_nodes(model, steel_of_concrete, bond_mass_share)
    model.steel_nodes = {sid for ids in steel_of_concrete.values() for sid in ids}
    edges = _edges(coords)
    _apply_supports_loads(model, problem, coords, edges, steel_of_concrete=steel_of_concrete)
    return model, edges


def build_continuum(
    problem: Problem,
    mesh_size: float,
    *,
    plane: str = "PlaneStress",
) -> tuple[Model, EdgeNodes]:
    """2D continuum builder: structured quads with nD ElasticIsotropic material.

    `plane` is "PlaneStress" (default) or "PlaneStrain" (D14, configurable).
    """
    coords, quads = _grid(problem, mesh_size)
    model = Model(ndm=2, ndf=2)
    mat = concrete_nd_elastic(problem.material, 1)
    model.nd_materials.append(mat)

    for idx, (x, y) in enumerate(coords, start=1):
        model.add_node(idx, (float(x), float(y)))
    thickness = problem.domain.thickness
    for eid, q in enumerate(quads, start=1):
        nodes = tuple(i + 1 for i in q)
        model.add_element(eid, "quad", nodes, (thickness, plane, mat.id), kind="quad")

    _assign_tributary_mass(model, coords, quads, problem)
    edges = _edges(coords)
    _apply_supports_loads(model, problem, coords, edges)
    return model, edges


def build_continuum_rc(
    problem: Problem,
    mesh_size: float,
    *,
    nd_material_for: "Callable[[str], tuple[NDMaterial, NDMaterial]]",
    zone_of: "Callable[[float, float], str]",
    rebars: "tuple[Rebar, ...]" = (),
    plane: str = "PlaneStress",
    rebar_tol: float = 1e-6,
    rebar_material: "Callable[[object, int], UniaxialMaterial]" = steel_uniaxial,
) -> tuple[Model, EdgeNodes]:
    """RC continuum builder (D29): per-zone nonlinear nD-concrete quads + steel rebar struts.

    The continuum verification reference that matches the RC lattice (D12/D14): the SAME structured
    node grid, plane-stress `quad` elements with a nonlinear nD concrete (ASDConcrete3D + PlaneStress
    wrapper, length-regularized — see `materials.concrete_nd_nonlinear`), and reinforcement as steel
    truss struts on the shared quad nodes (perfect bond, D5/D13). Like-for-like with `build_lattice_rc`
    so the two pushovers are directly comparable.

    `nd_material_for(zone)` returns the (ASDConcrete3D base, PlaneStress wrapper) NDMaterial pair for a
    zone; each pair is emitted once (cached by zone — the structured grid has a single quad size, so
    no per-length split as in the lattice) with builder-assigned, namespace-separate nDMaterial tags,
    and the quad uses the wrapper tag. The zone comes from `zone_of(x, y)` at the quad centroid. Each
    `Rebar` becomes steel struts on its on-path nodes via `rebar_material` (one material per grade).
    Same grid / mass / supports / loads plumbing as the other builders.
    """
    coords, quads = _grid(problem, mesh_size)
    model = Model(ndm=2, ndf=2)
    thickness = problem.domain.thickness

    for idx, (x, y) in enumerate(coords, start=1):
        model.add_node(idx, (float(x), float(y)))

    nd_tag = 1
    zone_wrapper: dict[str, int] = {}  # zone -> PlaneStress wrapper tag (one nD material pair per zone)
    for eid, q in enumerate(quads, start=1):
        cx, cy = coords[list(q)].mean(axis=0)
        zone = zone_of(float(cx), float(cy))
        if zone not in zone_wrapper:
            base, wrapper = nd_material_for(zone)
            base.id, wrapper.id, wrapper.args = nd_tag, nd_tag + 1, (nd_tag,)
            model.nd_materials.extend((base, wrapper))
            zone_wrapper[zone] = nd_tag + 1
            nd_tag += 2
        nodes = tuple(i + 1 for i in q)
        model.add_element(eid, "quad", nodes, (thickness, plane, zone_wrapper[zone]), kind="quad")

    eid = len(quads) + 1
    steel_tag_counter = 1
    steel_tag: dict[int, int] = {}  # id(SteelGrade) -> uniaxial material tag (one per grade)
    for rb in rebars:
        gid = id(rb.steel)
        if gid not in steel_tag:
            mat = rebar_material(rb.steel, steel_tag_counter)
            if isinstance(mat, (list, tuple)):
                raise TypeError("build_continuum_rc does not support multi-material rebar "
                                "factories yet; pass a single-material factory")
            model.uniaxial_materials.append(mat)
            steel_tag[gid] = steel_tag_counter
            steel_tag_counter += 1
        chain = rebar_node_chain(coords, rb.path, rebar_tol)
        for a, b in zip(chain, chain[1:]):
            model.add_element(eid, "Truss", (a + 1, b + 1), (rb.area, steel_tag[gid]), kind=rb.role)
            eid += 1

    _assign_tributary_mass(model, coords, quads, problem)
    edges = _edges(coords)
    _apply_supports_loads(model, problem, coords, edges)
    return model, edges


# --- shared helpers ---------------------------------------------------------

def _grid(problem: Problem, mesh_size: float):
    return mesh_compound_rectangles(problem.domain.rectangles(), mesh_size)


def _quad_area(pts: np.ndarray) -> float:
    """Shoelace area of a 4-node polygon."""
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _share_mass_with_steel_nodes(model: Model, steel_of_concrete: dict, share: float) -> None:
    """Move `share` of each concrete node's lumped mass onto the steel nodes sitting on it (D72).

    Conserves the total exactly: if a concrete node carries n steel nodes, each takes `share` of its
    mass and the concrete node is left with `1 - n*share`. Raises if that would go negative rather
    than silently producing a node with no mass.
    """
    for ci, steel_ids in steel_of_concrete.items():
        mass = model.masses.get(ci + 1)
        if mass is None:
            continue
        taken = share * len(steel_ids)
        if taken >= 1.0:
            raise ValueError(f"bond_mass_share={share} x {len(steel_ids)} steel nodes on concrete "
                             f"node {ci + 1} would take all of its mass; lower bond_mass_share")
        for sid in steel_ids:
            model.masses[sid] = tuple(share * m for m in mass)
        model.masses[ci + 1] = tuple((1.0 - taken) * m for m in mass)


def _assign_tributary_mass(model: Model, coords: np.ndarray, quads, problem: Problem) -> None:
    """Lump mass to nodes by tributary volume (rho * thickness * area/4 per quad corner), D16.

    Computed from the shared quad grid, so the lattice and continuum get identical nodal mass
    (same total, same distribution) — mass is independent of strut area.
    """
    rho = problem.material.rho
    thickness = problem.domain.thickness
    m = np.zeros(len(coords))
    for q in quads:
        share = rho * thickness * _quad_area(coords[list(q)]) / 4.0
        for idx in q:
            m[idx] += share
    for idx, mass in enumerate(m):
        model.masses[idx + 1] = (float(mass), float(mass))  # 2D: both translational DOFs


def _edges(coords: np.ndarray, *, tol: float = 1e-9) -> EdgeNodes:
    x, y = coords[:, 0], coords[:, 1]
    sel = {
        "xmin": np.isclose(x, x.min(), atol=tol),
        "xmax": np.isclose(x, x.max(), atol=tol),
        "ymin": np.isclose(y, y.min(), atol=tol),
        "ymax": np.isclose(y, y.max(), atol=tol),
    }
    return {edge: [int(i) + 1 for i in np.where(mask)[0]] for edge, mask in sel.items()}


def select_nodes(model: Model, box: tuple[float, float, float, float], *, tol: float = 1e-6,
                 kind: str = "any") -> list[int]:
    """Node ids inside an axis-aligned box (xmin, xmax, ymin, ymax), queried off a built Model.

    A post-build selector (the spec selectors above run during the build, off the coord array).
    Used by the pushover example to pick the control node, lateral-load nodes, and base nodes
    from the assembled model. Returns ids sorted ascending.

    `kind` MATTERS ONLY UNDER THE BOND SCHEME, where a bar node is duplicated at a concrete node's
    coordinates and a box therefore returns BOTH (D80 item 4):

      * `"any"` (default) — every node in the box. Correct for a base REACTION set: the anchored
        steel nodes carry part of the base shear and dropping them loses the bars' contribution.
      * `"concrete"` — the lattice nodes only. Correct for a DRIVE set: the actuator loads the
        concrete, and the bars follow through their bond links. Imposing the drive on the duplicate
        as well pins zero slip at the loaded row, which is perfect bond exactly where the bond
        model is being asked a question.
      * `"steel"` — the duplicates only; a probe, mostly for checks like "is every base bar node
        anchored".

    On a model built without bond, `model.steel_nodes` is empty and all three agree.
    """
    if kind not in ("any", "concrete", "steel"):
        raise ValueError(f"kind must be 'any', 'concrete' or 'steel'; got {kind!r}")
    xmin, xmax, ymin, ymax = box
    out = [
        nid
        for nid, n in model.nodes.items()
        if xmin - tol <= n.coords[0] <= xmax + tol and ymin - tol <= n.coords[1] <= ymax + tol
    ]
    if kind == "concrete":
        out = [nid for nid in out if nid not in model.steel_nodes]
    elif kind == "steel":
        out = [nid for nid in out if nid in model.steel_nodes]
    return sorted(out)


def _select(spec, coords: np.ndarray, edges: EdgeNodes, *, tol: float = 1e-9) -> list[int]:
    """Resolve a support/load spec to 1-based node ids (edge-based or box-based selection)."""
    if isinstance(spec, (EdgeSupport, EdgeLoad)):
        return edges[spec.edge]
    if isinstance(spec, (BoxSupport, BoxLoad)):
        xmin, xmax, ymin, ymax = spec.box
        x, y = coords[:, 0], coords[:, 1]
        mask = (x >= xmin - tol) & (x <= xmax + tol) & (y >= ymin - tol) & (y <= ymax + tol)
        return [int(i) + 1 for i in np.where(mask)[0]]
    raise TypeError(f"unknown support/load spec: {type(spec).__name__}")


def _apply_supports_loads(model: Model, problem: Problem, coords: np.ndarray, edges: EdgeNodes,
                          *, steel_of_concrete: "dict[int, list[int]] | None" = None) -> None:
    """Resolve support and load specs onto node ids.

    SUPPORTS PROPAGATE TO COINCIDENT STEEL NODES (D80 item 4). Specs resolve against the CONCRETE
    coordinate array, so under the bond scheme (`bond_material`, D72) a duplicated bar node at a
    fixed row was never itself fixed: on the Aldemir wall 30 steel nodes sat on the fixed base held
    only by their own ring of bond links, which fail at `ft/E * mesh` = 0.0037 mm of slip. The bar
    was effectively unanchored, and every bonded result before 2026-09-05 carries it.

    Fixing the duplicate is also the physical reading: a bar crossing the support row is anchored
    into the foundation. The consequence for post-processing is that a base REACTION sum must now
    include those steel nodes, or the bars' share of the base shear goes missing — select base nodes
    with `kind="any"` (the default), never `kind="concrete"`.
    """
    for bc in problem.supports:
        for nid in _select(bc, coords, edges):
            model.supports.append(Support(nid, bc.fix))
            for sid in (steel_of_concrete or {}).get(nid - 1, ()):
                model.supports.append(Support(sid, bc.fix))
    for ld in problem.loads:
        ids = _select(ld, coords, edges)
        if not ids:
            raise ValueError(f"load spec selected no nodes: {ld}")
        per = tuple(v / len(ids) for v in ld.total)
        for nid in ids:
            model.loads.append(Load(nid, per))


def critical_time_step(model: Model, modulus_of: "Callable[[Element], float]") -> tuple[float, float, object]:
    """Stability limit for an EXPLICIT integrator: `(dt_crit, w_max, worst_element)` (D74).

    An explicit march is conditionally stable at `dt < 2 / w_max`, and `w_max` is set by the
    STIFFEST, LIGHTEST element — NOT by the fundamental period. That is the trap: `steps_per_period`
    sized off T1 (the implicit habit) is typically 15-20x too large for an explicit run, and the
    result diverges rather than converging slowly, so it fails loudly rather than silently.

    Estimated element-wise as `w_e = 2*sqrt(k_e / m_min)` with `k_e = E*A/L` and `m_min` the smaller
    of the two end masses — the standard two-mass bound, and conservative for a lattice because a
    node carries several struts.

    `modulus_of(element)` supplies E, since a `Model` stores only material TAGS: a caller typically
    branches on `element.kind` to separate rebar from concrete.

    NOTE the stability limit falls further if stiffness-proportional damping is used — see
    `run_pushover_dynamic`, which drops the `betaKinit` term under an explicit integrator for
    exactly this reason.
    """
    masses = {nid: m[0] for nid, m in model.masses.items()}
    w_max, worst = 0.0, None
    for e in model.elements:
        if len(e.nodes) != 2 or not e.args:
            continue
        a = np.asarray(model.nodes[e.nodes[0]].coords, dtype=float)
        b = np.asarray(model.nodes[e.nodes[1]].coords, dtype=float)
        length = float(np.linalg.norm(b - a))
        m_min = min(masses.get(e.nodes[0], 0.0), masses.get(e.nodes[1], 0.0))
        if length <= 0.0 or m_min <= 0.0:
            continue
        w = 2.0 * math.sqrt(modulus_of(e) * float(e.args[0]) / length / m_min)
        if w > w_max:
            w_max, worst = w, e
    if w_max <= 0.0:
        raise ValueError("no element gave a positive frequency — check masses and areas")
    return 2.0 / w_max, w_max, worst
