"""The works' trenches: the ways they take on a small surface, and, once computed, the trenches against the routes
governed_drainage.py records (the same deepest cut, the branches' flows through them)."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import design_cuts as DC
from research.studies.atlas_round import governed_drainage as GD


def ridge():
    """Two basins, water at 0 m, either side of a ridge 300 m high whose lowest pass, 100 m, lies at column 25."""
    h = np.full((40, 50), -50.0)
    h[:, 20:30] = 300.0
    h[18:22, 20:30] = 100.0
    return h


def test_downhill_stops_at_the_water():
    h = np.array([[300.0, 200.0, 100.0, 50.0, -10.0, -20.0]])
    receiver = np.array([1, 2, 3, 4, 5, 5])
    assert DC.downhill(0, h.ravel(), receiver, 0.0) == [0, 1, 2, 3]


def test_the_way_between_passes_keeps_to_the_low_ground():
    h = ridge()
    cols = h.shape[1]
    lat = np.zeros(h.shape[0])
    way = DC.low_way(5 * cols + 20, 35 * cols + 29, h, 0.0, 1.9, lat)
    rows = np.array(way) // cols
    assert way[-1] == 5 * cols + 20 and way[0] == 35 * cols + 29
    assert np.mean(h.ravel()[way] <= 0.0) > 0.5                             # it leaves the ridge for the basins' floor


@pytest.fixture(scope='module')
def record():
    if not DC.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.design_cuts first')
    return json.loads(DC.OUT.read_text())


def test_the_trenches_are_the_governed_routes(record):
    s = record['summary']
    governed = json.loads(GD.OUT.read_text()) if GD.OUT.exists() else None
    recorded = [c['depth_m'] for k in (GD.HELD_BRANCH, GD.REST_BRANCH) for c in s['recorded_cuts'][k]]
    assert s['deepest_m'] == pytest.approx(max(recorded), abs=1.0)
    with np.load(DC.PRODUCT, allow_pickle=False) as z:
        flow = z['cut_flow_m3s']
    shares = s['shares_m3s']
    assert flow.max() == pytest.approx(shares[GD.HELD_BRANCH] + shares[GD.REST_BRANCH], rel=1e-4)
    allowed = np.array([shares[GD.HELD_BRANCH], shares[GD.REST_BRANCH], shares[GD.HELD_BRANCH] + shares[GD.REST_BRANCH]])
    assert np.all(np.abs(flow[:, None] - allowed[None, :]).min(axis=1) < 1.0)   # one branch's share, the other's, or both
    if governed is not None:
        assert governed['branches'][GD.HELD_BRANCH]['cuts'][0]['depth_m'] == pytest.approx(max(recorded), abs=1.0)
    assert record['product']['sha256'] == DC.sha256(DC.PRODUCT)


def test_each_stretch_runs_unbroken_over_the_cut_ground(record):
    """Every cut cell lies on a stretch's course, and each course steps from cell to neighbouring cell."""
    with np.load(DC.PRODUCT, allow_pickle=False) as z:
        cols = json.loads(str(z['metadata']))['grid']['cols']
        cut, course, stretch, level = z['cut_cell'], z['course_cell'], z['course_stretch'], z['course_level_m']
    assert set(cut.tolist()) <= set(course.tolist())
    assert level.size == stretch.max() + 1 == record['summary']['stretches']
    for s in range(level.size):
        cells = course[stretch == s].astype(np.int64)
        rows, columns = cells // cols, cells % cols
        steps = np.abs((np.diff(columns) + cols // 2) % cols - cols // 2)
        assert cells.size >= 1 and np.all(np.abs(np.diff(rows)) <= 1) and np.all(steps <= 1)
