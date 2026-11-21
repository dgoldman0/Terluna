"""The page builder's traced borders: chains between junctions, their smoothing and their simplification."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render as R                                     # noqa: E402


def test_chains_run_between_junctions_and_round_loops():
    # a T: three arms meeting at corner 0, and a separate square loop
    arms = [(0, 1), (1, 2), (0, 3), (3, 4), (0, 5)]
    loop = [(10, 11), (11, 12), (12, 13), (13, 10)]
    chains = R.chains_of(arms + loop)
    opened = sorted(tuple(c) for c, closed in chains if not closed)
    assert sorted(len(c) for c in opened) == [2, 3, 3] and all(c[0] == 0 or c[-1] == 0 for c in opened)
    closed = [c for c, closed in chains if closed]
    assert len(closed) == 1 and closed[0][0] == closed[0][-1] and len(closed[0]) == 5


def test_a_staircase_smooths_to_a_line_through_its_ends():
    steps = np.array([[i // 2 + (i % 2), i // 2] for i in range(41)], float)   # a 45 degree staircase of cell edges
    out = R.smooth_chain(steps, closed=False)
    assert np.allclose(out[0], steps[0]) and np.allclose(out[-1], steps[-1])   # the ends stay where borders meet
    inner = out[len(out) // 4: -len(out) // 4]
    off = np.abs(inner[:, 0] - inner[:, 1] - 0.5) / np.sqrt(2.0)               # distance from the staircase's diagonal
    assert off.max() < 0.1


def test_simplify_keeps_the_ends_and_drops_points_on_a_straight_run():
    x = np.linspace(0, 10, 21)
    line = np.stack([x, 1.0 - np.abs(x - 5.0) / 5.0], axis=1)                 # two straight runs and the corner
    out = R.simplify(line, 0.08)
    assert np.allclose(out[0], line[0]) and np.allclose(out[-1], line[-1])
    assert len(out) == 3 and np.allclose(out[1], line[10])
