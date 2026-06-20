"""Backend-agnostic problem definition (D12).

A `Problem` describes a specimen ONCE — geometry, materials (physical grades), supports and
loads — independently of how it will be discretised. Builders (lattice / continuum /
beam-column) translate the same Problem into an OpenSees FE Model.

Reinforcement (free 3D curves, D13) is not part of this first elastic-verification slice and
will be added here later.
"""

from __future__ import annotations

from dataclasses import dataclass

Edge = str  # one of: "xmin", "xmax", "ymin", "ymax"


@dataclass
class ConcreteGrade:
    """A physical concrete grade.

    `E, nu, rho` drive the elastic / modal models. The remaining fields parameterize the
    nonlinear uniaxial law (Concrete02, D15/D19) and are stored as POSITIVE magnitudes — the
    material mapping applies the OpenSees sign convention (compression negative). `ft`/`Ets`
    default in the mapping (ft ~ 0.1*fc, Ets ~ 0.1*E) when left None.
    """

    name: str
    E: float
    nu: float
    rho: float = 2400.0  # mass density (for modal analysis)
    fc: float | None = None     # peak compressive strength (magnitude)
    ft: float | None = None     # tensile strength (magnitude); mapping defaults ~0.1*fc
    Gf: float | None = None      # fracture energy (reserved; D15)
    # Concrete02 extras (magnitudes; signs applied in materials.py):
    epsc0: float | None = None  # strain at peak compression
    fcu: float | None = None    # crushing (residual) strength
    epsU: float | None = None   # crushing strain
    Ets: float | None = None    # tension softening slope; mapping defaults ~0.1*E
    lam: float = 0.1            # Concrete02 unloading-slope ratio (lambda)


@dataclass
class SteelGrade:
    """A physical reinforcing-steel grade -> Steel02 (Giuffre-Menegotto-Pinto), D19.

    `fy, E0, b` are the bilinear backbone; `R0, cR1, cR2` shape the elastic-plastic transition
    (Steel02 defaults). Reuse Steel01's exact match by setting a large R0 if ever needed.
    """

    name: str
    fy: float
    E0: float
    b: float = 0.01
    R0: float = 18.0
    cR1: float = 0.925
    cR2: float = 0.15


@dataclass
class Rebar:
    """A reinforcing bar: a polyline `path` (list of (x, y)) with cross-section `area` and a
    steel grade. In the lattice it becomes steel truss struts on the lattice nodes lying on the
    path (shared nodes => perfect bond, D5/D13). Mesh must align so nodes fall on the path."""

    path: list[tuple[float, float]]
    area: float
    steel: SteelGrade


@dataclass
class RectangleDomain:
    """A 2D rectangular member. `thickness` is the out-of-plane dimension (plane models)."""

    length: float
    height: float
    thickness: float = 1.0
    origin: tuple[float, float] = (0.0, 0.0)

    def rectangles(self) -> list[tuple[float, float, float, float]]:
        ox, oy = self.origin
        return [(ox, oy, self.length, self.height)]


@dataclass
class CompoundRectangles:
    """A 2D domain made of axis-aligned rectangles (ox, oy, width, height) that may share
    edges (e.g. a portal frame: 2 columns + a beam). Meshed structured per rectangle with
    coincident nodes merged at the joints (D-frame)."""

    rects: list[tuple[float, float, float, float]]
    thickness: float = 1.0

    def rectangles(self) -> list[tuple[float, float, float, float]]:
        return list(self.rects)


def  portal_frame(
    *,
    height: float = 3.0,
    span: float = 4.0,
    col_depth: float = 0.3,
    beam_depth: float = 0.3,
    thickness: float = 0.3,
) -> CompoundRectangles:
    """A one-bay, one-storey portal frame as compound rectangles: two columns (centerlines at
    x=0 and x=span, base at y=0) and a beam on top spanning the column outer faces."""
    half = col_depth / 2.0
    rects = [
        (-half, 0.0, col_depth, height),                 # left column
        (span - half, 0.0, col_depth, height),           # right column
        (-half, height, span + col_depth, beam_depth),   # beam on top
    ]
    return CompoundRectangles(rects=rects, thickness=thickness)


@dataclass
class EdgeSupport:
    """Fix DOFs on all nodes of an edge. `fix` is per-DOF (1 = fixed), length ndf."""

    edge: Edge
    fix: tuple[int, ...]


@dataclass
class EdgeLoad:
    """A total load applied to an edge, distributed equally over its nodes.

    `total` is the per-DOF resultant (length ndf), e.g. (0.0, -10e3) for 10 kN down.
    """

    edge: Edge
    total: tuple[float, ...]


@dataclass
class BoxSupport:
    """Fix DOFs on all nodes inside an axis-aligned box (xmin, xmax, ymin, ymax)."""

    box: tuple[float, float, float, float]
    fix: tuple[int, ...]


@dataclass
class BoxLoad:
    """A total load distributed equally over all nodes inside an axis-aligned box."""

    box: tuple[float, float, float, float]
    total: tuple[float, ...]


@dataclass
class Problem:
    """A fully-defined specimen, agnostic of discretisation."""

    ndm: int
    ndf: int
    domain: RectangleDomain
    material: ConcreteGrade
    supports: list[EdgeSupport]
    loads: list[EdgeLoad]
