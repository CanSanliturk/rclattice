"""Lattice-vs-continuum framework checks (D12): shared node grid + scalar calibration."""

import numpy as np

from rclattice import ConcreteGrade, EdgeLoad, EdgeSupport, Problem, RectangleDomain
from rclattice.builders import build_continuum, build_lattice
from rclattice.opensees import run_static


def _problem() -> Problem:
    return Problem(
        ndm=2,
        ndf=2,
        domain=RectangleDomain(length=1.0, height=0.2, thickness=0.1),
        material=ConcreteGrade(name="C30", E=30e9, nu=0.2),
        supports=[EdgeSupport(edge="xmin", fix=(1, 1))],
        loads=[EdgeLoad(edge="xmax", total=(0.0, -10e3))],
    )


def _tip(model, edges) -> float:
    res = run_static(model)
    assert res["ok"] == 0
    return float(np.mean([res["disps"][nid][1] for nid in edges["xmax"]]))


def test_lattice_and_continuum_share_node_grid():
    problem = _problem()
    cont, _ = build_continuum(problem, 0.1)
    lat, _ = build_lattice(problem, 0.1)
    assert len(cont.nodes) == len(lat.nodes)  # same gmsh grid -> fair comparison


def test_scalar_area_calibration_matches_continuum():
    problem = _problem()
    cont, ce = build_continuum(problem, 0.1)
    delta_c = _tip(cont, ce)

    lat0, le = build_lattice(problem, 0.1, strut_area=1.0)
    delta_l0 = _tip(lat0, le)
    area_star = delta_l0 / delta_c  # linear: deflection ~ 1/area

    lat, le = build_lattice(problem, 0.1, strut_area=area_star)
    assert np.isclose(_tip(lat, le), delta_c, rtol=1e-6)
