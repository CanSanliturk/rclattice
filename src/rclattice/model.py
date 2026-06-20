"""Backend-agnostic intermediate FE model produced by the builders (D12).

This is the assembly that the OpenSees backend translates. It is generic over element
type so it can hold truss struts (lattice), quads/bricks/shells (continuum) and
beam-columns alike. Builders populate it from a `Problem`; `opensees.py` emits it.

Units are the user's responsibility (use a consistent system, e.g. SI: N, m, Pa).
Object ids are reused directly as OpenSees tags.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Node:
    """A node. `coords` has length ndm (2 in 2D, 3 in 3D)."""

    id: int
    coords: tuple[float, ...]


@dataclass
class UniaxialMaterial:
    """A uniaxial material: ops.uniaxialMaterial(mtype, id, *args)."""

    id: int
    mtype: str
    args: tuple = ()


@dataclass
class NDMaterial:
    """An nD material: ops.nDMaterial(mtype, id, *args)."""

    id: int
    mtype: str
    args: tuple = ()


@dataclass
class Element:
    """A generic element: ops.element(etype, id, *nodes, *args)."""

    id: int
    etype: str
    nodes: tuple[int, ...]
    args: tuple = ()


@dataclass
class Support:
    """Boundary condition at a node. `fix` is per-DOF (1 = fixed, 0 = free), length ndf."""

    node: int
    fix: tuple[int, ...]


@dataclass
class Load:
    """A nodal load. `values` is per-DOF, length ndf."""

    node: int
    values: tuple[float, ...]


@dataclass
class Model:
    """An assembled FE model. Dimension-agnostic via ndm/ndf."""

    ndm: int
    ndf: int
    nodes: dict[int, Node] = field(default_factory=dict)
    elements: list[Element] = field(default_factory=list)
    uniaxial_materials: list[UniaxialMaterial] = field(default_factory=list)
    nd_materials: list[NDMaterial] = field(default_factory=list)
    supports: list[Support] = field(default_factory=list)
    loads: list[Load] = field(default_factory=list)
    masses: dict[int, tuple[float, ...]] = field(default_factory=dict)  # node id -> per-DOF mass

    def add_node(self, node_id: int, coords: tuple[float, ...]) -> Node:
        if len(coords) != self.ndm:
            raise ValueError(f"node {node_id}: expected {self.ndm} coords, got {len(coords)}")
        node = Node(node_id, coords)
        self.nodes[node_id] = node
        return node

    def add_element(self, eid: int, etype: str, nodes: tuple[int, ...], args: tuple = ()) -> Element:
        el = Element(eid, etype, tuple(nodes), tuple(args))
        self.elements.append(el)
        return el
