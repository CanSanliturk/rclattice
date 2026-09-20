"""The rebar-aligned graded grid (D104): equivalence on a uniform panel, exact bar lines, scaling."""
import numpy as np

from rclattice.mesh import (connect_horizon, connect_index_horizon, graded_lines,
                            mesh_rectangle_grid, mesh_rectangle_lines, tributary_area_scale)


def _keys(c, pairs):
    return {tuple(sorted((tuple(np.round(c[i], 6)), tuple(np.round(c[j], 6))))) for i, j in pairs}


def test_graded_grid_reduces_to_the_uniform_grid():
    c1, q1 = mesh_rectangle_grid(200.0, 100.0, 25.0)
    c2, q2 = mesh_rectangle_lines(200.0, 100.0, 25.0)
    assert len(c1) == len(c2) == 45 and len(q1) == len(q2) == 32
    for h in (1.5, 3.01):
        assert _keys(c1, connect_horizon(c1, 25.0, h)) == _keys(c2, connect_index_horizon(c2, h))
    assert np.allclose(tributary_area_scale(c2, connect_index_horizon(c2, 1.5), 25.0), 1.0)


def test_hard_lines_land_on_nodes_and_spacing_stays_near_target():
    hard = [19.0, 70.0, 121.0, 172.0, 323.5]
    xs = graded_lines(400.0, 25.0, hard)
    assert all(any(abs(x - h) < 1e-9 for x in xs) for h in hard)
    d = np.diff(xs)
    assert d.min() >= 19.0 - 1e-9 and d.max() <= 25.5 + 1e-9
    c, q = mesh_rectangle_lines(400.0, 100.0, 25.0, x_lines=hard, y_lines=[38.0])
    for h in hard:
        assert np.isclose(c[:, 0], h).sum() == len(np.unique(np.round(c[:, 1], 6)))
    pairs = connect_index_horizon(c, 1.5)
    s = tributary_area_scale(c, pairs, 25.0)
    assert 0.7 < s.min() < 1.0 < s.max() + 0.05 and len(pairs) == len(s)
