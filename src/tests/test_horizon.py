"""Sanity checks for the horizon connectivity rule (D9)."""

import numpy as np

from rclattice.mesh import connect_horizon


def test_horizon_captures_orthogonal_and_diagonal_not_next_ring():
    # unit 3x3 grid, spacing s = 1.0
    coords = np.array([[x, y] for y in (0.0, 1.0, 2.0) for x in (0.0, 1.0, 2.0)])
    pairs = connect_horizon(coords, mesh_size=1.0, horizon=1.5)
    dists = sorted({round(float(np.linalg.norm(coords[i] - coords[j])), 6) for i, j in pairs})
    # horizon 1.5 -> includes 1.0 (orthogonal) and sqrt(2)~1.414 (diagonal), excludes 2.0
    assert dists == [1.0, round(float(np.sqrt(2)), 6)]


def test_horizon_one_element_per_pair():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    pairs = connect_horizon(coords, mesh_size=1.0, horizon=1.5)
    assert len(pairs) == len(set(pairs))
    assert all(i < j for i, j in pairs)
