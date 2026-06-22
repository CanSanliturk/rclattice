"""Builders (D12): translate a backend-agnostic `Problem` into an OpenSees FE `Model`.

Both builders mesh the SAME structured node grid (mesh_rectangle_grid), so the lattice and
the continuum share an identical node set and the verification comparison is fair. Each
builder also returns the edge->node-id sets so callers can query/compare results.

Builders are backend-agnostic: they produce a `Model`; they never import openseespy.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .materials import concrete_nd_elastic, concrete_uniaxial_elastic, steel_uniaxial
from .mesh import connect_horizon, mesh_compound_rectangles
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
    """
    coords, quads = _grid(problem, mesh_size)
    area_fn: Callable[[float], float] = strut_area if callable(strut_area) else (lambda _L: strut_area)

    model = Model(ndm=2, ndf=2)
    tag = 1

    for idx, (x, y) in enumerate(coords, start=1):
        model.add_node(idx, (float(x), float(y)))

    eid = 1
    mat_cache: dict[tuple[str, float], int] = {}  # (zone, rounded length) -> material tag
    for (i, j) in connect_horizon(coords, mesh_size, horizon):
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
    for rb in rebars:
        gid = id(rb.steel)
        if gid not in steel_tag:
            mat = rebar_material(rb.steel, tag)
            model.uniaxial_materials.append(mat)
            steel_tag[gid] = tag
            tag += 1
        chain = rebar_node_chain(coords, rb.path, rebar_tol)
        for a, b in zip(chain, chain[1:]):
            model.add_element(eid, strut_element, (a + 1, b + 1), (rb.area, steel_tag[gid]), kind=rb.role)
            eid += 1

    _assign_tributary_mass(model, coords, quads, problem)
    edges = _edges(coords)
    _apply_supports_loads(model, problem, coords, edges)
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


def select_nodes(model: Model, box: tuple[float, float, float, float], *, tol: float = 1e-6) -> list[int]:
    """Node ids inside an axis-aligned box (xmin, xmax, ymin, ymax), queried off a built Model.

    A post-build selector (the spec selectors above run during the build, off the coord array).
    Used by the pushover example to pick the control node, lateral-load nodes, and base nodes
    from the assembled model. Returns ids sorted ascending.
    """
    xmin, xmax, ymin, ymax = box
    out = [
        nid
        for nid, n in model.nodes.items()
        if xmin - tol <= n.coords[0] <= xmax + tol and ymin - tol <= n.coords[1] <= ymax + tol
    ]
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


def _apply_supports_loads(model: Model, problem: Problem, coords: np.ndarray, edges: EdgeNodes) -> None:
    for bc in problem.supports:
        for nid in _select(bc, coords, edges):
            model.supports.append(Support(nid, bc.fix))
    for ld in problem.loads:
        ids = _select(ld, coords, edges)
        if not ids:
            raise ValueError(f"load spec selected no nodes: {ld}")
        per = tuple(v / len(ids) for v in ld.total)
        for nid in ids:
            model.loads.append(Load(nid, per))
