"""Perturbed lattice + energy-balance calibration for Aydin's compression cube (D60).

Three pieces, in the order the paper applies them:

  1. `perturbed_grid` — the paper's novelty. A structured grid whose nodes are displaced by
     `R <= Rmax` at a random angle (his Fig. 2(b)). The strut TOPOLOGY is taken from the
     unperturbed grid, so only the geometry is random, not the connectivity.
  2. `calibrate` — the elastic energy balance (his p. 396, and the same Eq. 2.1 the wall study
     uses): equate the strain energy a continuum stores under a uniform strain field to the
     lattice's, and solve for one uniform `E*A`. Closed form; no FE solve, no optimiser.
  3. `cube_lattice` — assembly, with `materials.concrete_lattice_aydin` on every strut.

The compressive strength is an OUTPUT of the analysis, not an input to any of the three.
"""

from __future__ import annotations

import numpy as np

from rclattice.builders import build_lattice_rc
from rclattice.calibration import EnergyBalanceResult, energy_balance_area
from rclattice.materials import aydin_lattice_softening, concrete_lattice_aydin
from rclattice.mesh import connect_horizon, mesh_rectangle_grid, perturb_nodes

from specimen import (
    A1, A3_OVER_A2, B1, B2, CONCRETE, EPS, GF, HORIZON, L, MESH, NU, RSM_ALPHA, RSM_BETA, THK,
    cube_problem,
)


def perturbed_grid(mesh_size: float = MESH, *, rmax_ratio: float, seed: int | None = None,
                   horizon: float = HORIZON, boundary: str = "LR"):
    """`(coords, quads, pairs)` for one realization of the perturbed cube.

    `pairs` come from the UNPERTURBED grid on purpose. Re-running `connect_horizon` on the moved
    nodes drops every strut whose length happens to cross `horizon*mesh_size` — 5 of 420 at
    `Rmax/d = 0.08` — which would make the lattice's connectivity a second random variable on top
    of its geometry, and would quietly change the number of load paths from realization to
    realization. Aydin perturbs a grid and keeps its cells (his Fig. 2(b) still shows both
    diagonals in every cell), so topology-from-the-regular-grid is the faithful reading.

    For NR the bottom node nearest mid-width is held exactly in place, because it carries the
    single ux restraint that keeps the implicit solve non-singular and a box query has to find it.
    """
    coords0, quads = mesh_rectangle_grid(L, L, mesh_size)
    pairs = connect_horizon(coords0, mesh_size, horizon)

    fix: set[int] = set()
    if boundary == "NR":
        bottom = np.flatnonzero(np.abs(coords0[:, 1]) < EPS)
        fix.add(int(bottom[np.abs(coords0[bottom, 0] - L / 2.0).argmin()]))

    coords = perturb_nodes(coords0, mesh_size, rmax_ratio, seed=seed, fix_nodes=fix)
    return coords, quads, pairs


def calibrate(coords, pairs) -> EnergyBalanceResult:
    """Aydin's elastic energy balance over this realization's own geometry.

    Recalibrated per realization rather than once on the regular grid: perturbation changes every
    strut length, so the lattice that stores the continuum's energy is a slightly different one
    each time. The spread is small (well under 1%) but computing it costs nothing and keeps the
    elastic modulus of every realization on target instead of only the mean.
    """
    return energy_balance_area(coords, pairs, E=CONCRETE.E, nu=NU, thickness=THK,
                               area_inplane=L * L)


def material_factory(*, rsm: bool = True, Gf: float = GF):
    """`material_for(zone, length)` hook emitting Aydin's tension-only strut law.

    Every strut gets its own `a2`/`a3` from its own length, so the fracture energy per unit crack
    area is `Gf` regardless of grid spacing or of how perturbation stretched that particular strut.
    """
    def material_for(_zone: str, length: float):
        return concrete_lattice_aydin(CONCRETE, 0, length, Gf=Gf, a1=A1, b1=B1, b2=B2,
                                      a3_over_a2=A3_OVER_A2, rsm=rsm,
                                      alpha=RSM_ALPHA, beta=RSM_BETA)
    return material_for


def cube_lattice(coords, quads, pairs, area: float, *, mesh_size: float = MESH,
                 horizon: float = HORIZON, boundary: str = "LR", rsm: bool = True,
                 Gf: float = GF, strut_element: str = "corotTruss"):
    """Assemble the analysis model for one realization.

    `strut_element` defaults to **corotTruss**, and that is a modelling decision, not a default
    carried over from elsewhere. Aydin's failure mode is the loss of stability of the vertical load
    path — a geometric event. A small-displacement `Truss` formulation cannot express it: the
    library's own D53 experiment ran this material's unbounded compression on a straight grid and
    the stress climbed monotonically to 6*fc with no failure at all, while the same run on
    corotational struts peaked and then collapsed. Aydin's explicit scheme integrates the nodes'
    actual positions, so he has this for free and never has to mention it.
    """
    return build_lattice_rc(
        cube_problem(coords, boundary), mesh_size,
        material_for=material_factory(rsm=rsm, Gf=Gf),
        zone_of=lambda _x, _y: "concrete",
        rebars=(),
        strut_area=area,
        horizon=horizon,
        strut_element=strut_element,
        grid=(coords, quads),
        pairs=pairs,
    )[0]


def softening_report(mesh_size: float = MESH, *, Gf: float = GF) -> str:
    """The derived trilinear tail for the two strut lengths of an unperturbed grid, vs Table 1."""
    rows = []
    for label, length in (("orthogonal", mesh_size), ("diagonal", mesh_size * 2.0 ** 0.5)):
        a2, a3 = aydin_lattice_softening(CONCRETE.ft, CONCRETE.E, length, Gf=Gf,
                                         a1=A1, b1=B1, b2=B2, a3_over_a2=A3_OVER_A2)
        rows.append(f"{label} L = {length:5.2f} mm -> a2 = {a2:5.1f}, a3 = {a3:6.1f}")
    return "; ".join(rows)
