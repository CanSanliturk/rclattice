"""Portal-frame geometry, box selectors, mode shapes, and visualizer smoke tests (D-frame)."""

from rclattice import viz
from rclattice.builders import build_continuum, build_lattice
from rclattice.opensees import run_modal, run_static
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, portal_frame

H, SPAN = 1.0, 1.0
MESH = 0.1


def _frame() -> Problem:
    domain = portal_frame(height=H, span=SPAN, col_depth=0.2, beam_depth=0.2, thickness=0.2)
    return Problem(
        ndm=2,
        ndf=2,
        domain=domain,
        material=ConcreteGrade(name="C30", E=30e9, nu=0.2, rho=2400.0),
        supports=[BoxSupport(box=(-1.0, SPAN + 1.0, -1e-3, 1e-3), fix=(1, 1))],
        loads=[BoxLoad(box=(-1.0, SPAN + 1.0, H - 1e-3, H + 0.2 + 1e-3), total=(10e3, 0.0))],
    )


def test_frame_lattice_continuum_share_grid_and_run():
    problem = _frame()
    lat, _ = build_lattice(problem, MESH, strut_area=1e-2)
    cont, _ = build_continuum(problem, MESH)
    assert len(lat.nodes) == len(cont.nodes)        # compound mesh merged, shared grid
    assert any(len(e.nodes) == 2 for e in lat.elements)   # struts
    assert any(len(e.nodes) == 4 for e in cont.elements)  # quads
    assert run_static(lat)["ok"] == 0
    assert run_static(cont)["ok"] == 0


def test_frame_modes_have_shapes_and_positive_periods():
    problem = _frame()
    cont, _ = build_continuum(problem, MESH)
    modal = run_modal(cont, 3)
    assert len(modal["shapes"]) == 3
    assert all(p > 0 for p in modal["periods"])
    assert set(modal["shapes"][0]) == set(cont.nodes)  # one vector per node


def test_visualizer_renders(tmp_path):
    problem = _frame()
    lat, _ = build_lattice(problem, MESH, strut_area=1e-2)
    cont, _ = build_continuum(problem, MESH)
    lat_d = run_static(lat)["disps"]
    cont_d = run_static(cont)["disps"]
    out = tmp_path / "static.png"
    viz.figure_static(lat, lat_d, cont, cont_d, savepath=str(out))
    assert out.exists() and out.stat().st_size > 0
