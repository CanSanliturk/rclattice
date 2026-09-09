"""Aydin, Binici & Tuncay (2021) uniaxial compression cube — the specimen, verbatim (D60).

Reproduces the cube of *Lattice simulation of concrete compressive behaviour as indirect tension
failure*, Magazine of Concrete Research 73(8): 394-409. NOTE this is a DIFFERENT source from the
2017 METU MSc thesis that `examples/wall/` and `examples/compression_cube/` calibrate against; the
elastic energy balance (p. 396) is the same, but everything about compression here is new.

The point of the paper, and therefore of this example: the struts have NO compressive strength.
Compression is linear elastic forever. The cube's compressive strength is not a material input —
it emerges as the stability loss of a load path whose lateral restraint has been cracked away by
transverse tension. The single calibration parameter is geometric, the grid perturbation `Rmax/d`.

Paper values used here (his "first" study, p. 396):
    100 x 100 x 100 mm cube, grid d = 10 mm, dynamic explicit solution, uniform displacement
    imposed on the uppermost nodes; fc = 20 MPa and Gf = 50 N/m CHOSEN, with ft and Ec then
    computed from the ACI 318-19 / TS 500 expressions he cites:
        ft = 0.35*sqrt(fc)      (TS 500)
        Ec = 4700*sqrt(fc)      (ACI 318)
Two boundary conditions, his Fig. 2(a): LR (lateral restraint) and NR (no restraint).

Units: N, mm (stresses MPa). Backend-agnostic — nothing here imports openseespy (D8).
"""

from __future__ import annotations

from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.problem import BoxSupport, ConcreteGrade, Problem, RectangleDomain

# --- geometry (p. 396: "100 x 100 x 100 mm cube specimen", "a grid size of 10 mm") --------------
L, THK = 100.0, 100.0
MESH = 10.0
HORIZON = 1.5                # 8 neighbours: the orthogonal + diagonal cell of his Fig. 2(b)
EPS = 1e-6
A_FACE = L * THK
OUT = Path(__file__).resolve().parent.parent / "output" / "aydin_cube"

# --- material (p. 396) --------------------------------------------------------------------------
FC = 20.0                                   # chosen
GF = 0.050                                  # chosen, 50 N/m -> N/mm
FT = 0.35 * FC ** 0.5                       # TS 500          -> 1.5652 MPa
EC = 4700.0 * FC ** 0.5                     # ACI 318         -> 21019 MPa
EPSC0 = 2.0 * FC / EC                       # p. 397, "epsc0 = 2*fc/Ec" -> 1.903e-3
NU = 0.20                                   # not stated in the paper; used only by the balance

# `fc` and `epsc0` are carried on the grade for the RSM knee and for reporting only. NOTHING in
# this example gives a strut compressive strength — that is the whole idea (see `build.py`).
CONCRETE = ConcreteGrade("aydin-fc20", E=EC, nu=NU, rho=2.4e-9,
                         fc=FC, epsc0=EPSC0, ft=FT, Gf=GF)

# Trilinear tension-softening shape, his Fig. 1(a) caption and Table 1 (a1 = 1.5 for every
# specimen). `a2`/`a3` are NOT fixed here: they are solved per strut from Gf so the model stays
# mesh-objective — see `rclattice.materials.aydin_lattice_softening`.
A1, B1, B2 = 1.5, 0.6, 0.2
A3_OVER_A2 = 5.0                            # Jansen & Shah column of Table 1 (300/60)

# Reduced stiffness model, his Fig. 1(b): 40% of the modulus past a third of epsc0.
RSM_ALPHA, RSM_BETA = 1.0 / 3.0, 0.4

# Grid perturbation. p. 397 fits Rmax = 0.6 mm (NR) and 0.75 mm (LR) at d = 10 mm for this very
# specimen, i.e. Rmax/d = 0.060 and 0.075. Those are the paper's ANSWERS for fc = 20, Gf = 50;
# `--calibrate` re-derives them here rather than assuming they transfer.
RMAX_RATIO = {"LR": 0.075, "NR": 0.060}

# --- loading ------------------------------------------------------------------------------------
# His Fig. 3 runs to ~0.30 mm of shortening on a 100 mm cube; 0.4% strain covers peak and softening.
TARGET_STRAIN = 0.004
TARGET = -TARGET_STRAIN * L

# "Uniform displacement was imposed on the uppermost nodes at a sufficiently slow speed to simulate
# a quasistatic test" (p. 396). His Table 1 lists v = 0.1-20 mm/s for the validation specimens; no
# value is given for the cube. 1.0 mm/s is quasi-static here (`--rate` sweeps it; the run prints an
# inertia check so the claim is measured, not asserted).
QUASI_STATIC_RATE = 1.0

REALIZATIONS = 5             # p. 399: "five different simulation sets were run", results averaged


def cube_problem(coords, boundary: str = "LR") -> Problem:
    """The cube under one of the paper's two end conditions (his Fig. 2(a)).

    LR — *lateral restraint*: "pin supports at the bottom nodes and a roller support at the top
    nodes". The bottom face is fully fixed and the top face is held against lateral movement while
    being pushed down, so the platens confine the specimen. This is the rough/bonded-platen bound,
    and it is the case that develops the two diagonal cracks of his Fig. 3.

    NR — *no lateral restraint*: "only roller supports at the bottom nodes and leaving the top
    nodes free". Both faces may spread. One extra node — the bottom node nearest mid-width — also
    has ux fixed, which the paper does not need and we do: his explicit scheme integrates a mass
    matrix and tolerates a free-floating rigid-body mode, while the implicit Newmark solve here
    would see a singular stiffness. It removes exactly that mode and restrains nothing else.

    `coords` are the PERTURBED node positions, so the supports are selected against the geometry
    that is actually built. Boundary nodes only ever slide along their own face
    (`mesh.perturb_nodes`), so the faces stay flat and these boxes stay exact.
    """
    if boundary not in ("LR", "NR"):
        raise ValueError(f"boundary must be 'LR' or 'NR', got {boundary!r}")
    domain = RectangleDomain(length=L, height=L, thickness=THK, origin=(0.0, 0.0))

    if boundary == "LR":
        supports = [
            BoxSupport(box=(-EPS, L + EPS, -EPS, EPS), fix=(1, 1)),          # pinned base
            BoxSupport(box=(-EPS, L + EPS, L - EPS, L + EPS), fix=(1, 0)),   # top roller: ux held
        ]
    else:
        # Partition the bottom edge so no node is fixed twice (`build` emits one ops.fix per spec).
        xpin = _pin_x(coords)
        supports = [
            BoxSupport(box=(-EPS, xpin - EPS, -EPS, EPS), fix=(0, 1)),
            BoxSupport(box=(xpin - EPS, xpin + EPS, -EPS, EPS), fix=(1, 1)),
            BoxSupport(box=(xpin + EPS, L + EPS, -EPS, EPS), fix=(0, 1)),
        ]
    return Problem(ndm=2, ndf=2, domain=domain, material=CONCRETE, supports=supports, loads=[])


def _pin_x(coords) -> float:
    """x of the bottom-edge node nearest mid-width — the NR stability pin."""
    bottom = coords[abs(coords[:, 1]) < EPS]
    return float(bottom[abs(bottom[:, 0] - L / 2.0).argmin(), 0])


def top_nodes(model) -> list[int]:
    """The uppermost nodes — the loaded face, driven downward together."""
    return select_nodes(model, (-EPS, L + EPS, L - EPS, L + EPS))


def base_nodes(model) -> list[int]:
    """The supported bottom face — where the axial reaction is summed."""
    return select_nodes(model, (-EPS, L + EPS, -EPS, EPS))


def control_node(model) -> int:
    """Displacement-control node: the top-face node nearest mid-width."""
    top = top_nodes(model)
    return min(top, key=lambda n: abs(model.nodes[n].coords[0] - L / 2.0))


def boundary_note(boundary: str) -> str:
    """One-line description of the end condition, for figure titles and printed reports."""
    return ("LR: pinned base + laterally held top platen (rough platens, confined ends)"
            if boundary == "LR" else
            "NR: rollers under the base, top free to spread (smooth platens)")
