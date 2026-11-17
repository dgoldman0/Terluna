"""Checks of the dryland rivers: a stretch's sum runs down the routing to the first stop and leaves the stop out, the
water surface takes each sea's and lake's level, a route costs the near-side sea its loss, or its whole flow where it
ends in a closed lake, and the record holds together once computed."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import dryland_rivers as DR


def test_a_stretch_sums_down_to_the_first_stop():
    # a chain 0 -> 1 -> 2 -> 3 -> 4 (an outlet), and a branch 5 -> 2; cell 3 is a stop
    receiver = np.array([1, 2, 3, 4, 4, 2])
    value = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    stop = np.array([False, False, False, True, True, False])
    total, end = DR.downstream_sum(receiver, value, stop)
    np.testing.assert_allclose(total, [7.0, 6.0, 4.0, 0.0, 0.0, 36.0])
    assert list(end) == [3, 3, 3, 3, 4, 3]


def test_a_long_chain_converges():
    n = 3000
    receiver = np.minimum(np.arange(n) + 1, n - 1)
    stop = np.zeros(n, bool)
    stop[-1] = True
    total, end = DR.downstream_sum(receiver, np.ones(n), stop)
    np.testing.assert_allclose(total, np.arange(n, 0, -1) - 1.0)
    assert (end == n - 1).all()


def test_the_water_surface_takes_seas_and_lakes():
    meta = dict(seas=['near_side', 'orientale'], levels_m={'near_side': -1615.0, 'orientale': 726.0})
    sea = np.array([0, -1, 1, -1, -1], np.int8)
    h = np.array([-2000.0, 50.0, 300.0, 10.0, 5.0])
    surface = DR.water_surface(meta, sea, np.array([3]), np.array([12.5]), h)
    assert surface[0] == -1615.0 and surface[2] == 726.0 and surface[3] == 12.5
    assert np.isnan(surface[1]) and np.isnan(surface[4])


def test_what_the_near_side_gives_up():
    names = ['near_side', 'south_pole_aitken', 'orientale']
    # cells: 0 near-side sea, 1 salt sea, 2 a river into the near side, 3 a river into the salt sea, 4 a closed lake
    receiver = np.array([0, 1, 0, 1, 4])
    r = dict(receiver=receiver, discharge=np.array([0.0, 0.0, 2000.0, 2000.0, 0.0]))
    sea = np.array([0, 1, -1, -1, -1], np.int8)
    term = np.array([0, 1, 0, 1, 4])
    rate = np.full(5, 1e-2)
    route = dict(free_cell=2, end_cell=0, cut_m=50.0, full_sea=True, source_sea=2, source_outflow_m3s=1e4,
                 source_level=700.0, pass_cell=2, shore_cell=2, source_terminal=0)
    dry = np.full(5, 100.0)
    length = np.full(5, 120.0)
    back = DR.evaluate(route, r, sea, names, term, dry, length, rate, 5, np.eye(3), 0.03)
    assert back['ends_in'] == 'near_side' and back['source_water_ends_in'] == 'near_side'
    assert back['near_side_gives_up_m3s'] == back['loss_m3s'] and back['loss_m3s']['wide'] > 0
    salt = DR.evaluate(dict(route, free_cell=3, end_cell=1), r, sea, names, term, dry, length, rate, 5, np.eye(3), 0.03)
    assert salt['ends_in'] == 'south_pole_aitken'
    assert salt['near_side_gives_up_m3s'] == salt['loss_m3s']                   # the weir tops the salt sea up
    closed = DR.evaluate(dict(route, free_cell=4, end_cell=4), r, sea, names, term, dry, length, rate, 5, np.eye(3), 0.03)
    assert closed['ends_in'] == 'a closed lake'
    assert closed['near_side_gives_up_m3s'] == {'narrow': DR.Q_ROUTE, 'wide': DR.Q_ROUTE}
    assert closed['near_side_south_part_falls_m']['narrow'] == pytest.approx(0.03 * DR.Q_ROUTE, abs=0.05)


@pytest.fixture(scope='module')
def record():
    if not DR.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.dryland_rivers first')
    return json.loads(DR.OUT.read_text())


def test_the_record_holds_together(record):
    routes = record['routes']
    ranks = [x['green_km2_per_m_of_cut'] for x in routes]
    assert ranks == sorted(ranks, reverse=True)
    for x in routes:
        assert 0 < x['cut_m'] <= DR.CUT_MAX_M
        assert x['dryland_km'] >= DR.DRY_MIN_KM and x['dryland_km'] <= x['new_river_km'] + 1
        assert x['green_km2']['narrow'] <= x['green_km2']['wide']
        assert x['loss_m3s']['wide'] <= x['flow_m3s'] + 1
    s = record['summary']
    assert s['arid_and_semi_arid_land_km2'] <= s['land_km2']
    assert sum(record['drylands_by_turned_latitude_km2'].values()) == pytest.approx(s['arid_and_semi_arid_land_km2'], rel=0.01)
