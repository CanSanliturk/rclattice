"""Shared specimen definition for the plain-concrete cube uniaxial-tension study.

A 0.1 m concrete cube modelled two ways and compared under uniaxial tension:
  - reference: a SINGLE displacement-based fiber beam-column (the 1D material coupon);
  - lattice:   the same specimen as gmsh nodes + horizon struts (Concrete02), with the strut
               area calibrated so the lattice's initial axial stiffness equals the beam-column's.

Because the lattice/continuum tooling is 2D (the gmsh mesher builds a plane), the cube is
represented as a 2D PLANE-STRESS SQUARE: a 100 x 100 mm face with a 100 mm out-of-plane
thickness (so the section perpendicular to the pull is L x THK = the cube face). Uniaxial tension
is fully captured this way; nothing here imports openseespy (D8).

Units: N, mm (so the 0.1 m edge is 100 mm, E = 30000 N/mm^2 = 30 GPa, stresses in MPa). This keeps
displacements O(0.01 mm) — well conditioned for the default solver tolerances. The specimen is
backend-agnostic (geometry, grade, supports); the model BUILDERS + calibration live in `build.py`,
the entry script in `pushover.py`. Run as `python examples/cube/pushover.py`.
"""

from __future__ import annotations

from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxSupport, ConcreteGrade, Problem, RectangleDomain

L, THK = 100.0, 100.0        # cube edge (0.1 m) and out-of-plane thickness (mm)
MESH, EPS = 10.0, 1e-6       # lattice node spacing (10 divisions across the cube) and geom tol
HORIZON = 1.5                # strut connectivity: pairs within HORIZON * MESH are connected (D9)
DU, TARGET = 1.0e-3, 0.12    # tension pushover: control-disp increment and target (mm)
PULL_REF = 1.0               # axial reference-load magnitude (shape only; DisplacementControl scales it)
A_FACE = L * THK             # cross-section perpendicular to the pull (mm^2) -> stress = N / A_FACE
OUT = Path(__file__).resolve().parent.parent / "output" / "cube"

# Plain concrete, generic C30 (N, mm). Compression params are present for Concrete02 but not
# exercised by the tension pull; ft = 3 MPa is the tensile strength that caps the coupon. `ft`/`Ets`
# left to the mapping defaults (ft as given; Ets ~ 0.1*E) -> a plain (unregularized) softening branch.
CONCRETE = ConcreteGrade("C30", E=30000.0, nu=0.2, rho=2.4e-9,
                         fc=30.0, epsc0=0.002, fcu=6.0, epsU=0.006, ft=3.0)


def cube_problem() -> Problem:
    """The backend-agnostic cube: a 100 x 100 mm plane-stress square (thickness THK) with a
    FULLY FIXED base (bottom edge fixed in X and Y) and no body loads — the tension pull is applied
    as a separate reference pattern by the runner (see `axial_loads`)."""
    domain = RectangleDomain(length=L, height=L, thickness=THK, origin=(0.0, 0.0))
    supports = [BoxSupport(box=(-EPS, L + EPS, -EPS, EPS), fix=(1, 1))]   # fully fixed base
    return Problem(ndm=2, ndf=2, domain=domain, material=CONCRETE, supports=supports, loads=[])


def axial_loads(model) -> list[Load]:
    """Uniaxial tension reference pattern: an equal upward (+Y) nodal force over the top edge,
    summing to `PULL_REF`. DisplacementControl drives its magnitude; the shape pulls the whole top
    face uniformly so the specimen stays in uniform axial tension."""
    ids = select_nodes(model, (-EPS, L + EPS, L - EPS, L + EPS))
    return [Load(nid, (0.0, PULL_REF / len(ids))) for nid in ids]
