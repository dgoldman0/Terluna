"""Checks of the atlas round's satoyama candidate map: the side valleys, the methane rates, and the product once it has been
computed."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import satoyama_candidates as S


def test_side_valleys_end_where_they_reach_a_large_river_water_or_a_closed_hollow():
    #            0     1      2      3      4     5     6     7    8    9     10    11    12    13
    receiver = [0,    0,     1,     2,     2,    4,    5,    3,   7,   0,    9,    13,   11,   0]
    drained = [0.0, 500.0, 400.0, 300.0, 50.0, 30.0, 10.0, 5.0, 2.0, 60.0, 20.0, 40.0, 15.0, 5.0]
    water = np.zeros(14, bool)
    water[0] = True                                   # the sea; cells 1-3 a river draining more than 100 km2
    unit = S.units(np.array(receiver), np.array(drained), water)
    # a side valley into the river (4-6), one too small (7-8), one straight into the sea (9-10) and a closed hollow
    # (11-12) whose spill carries none of its land on (cell 13 drains less than the hollow)
    assert unit.tolist() == [-1, -1, -1, -1, 4, 4, 4, -1, -1, 9, 9, 11, 11, -1]


def test_the_paddy_rates_are_satoyama_mds():
    assert S.rate_t_km2() == pytest.approx(7.4, abs=0.05)                  # drained more than once, cropped once
    assert 2 * S.rate_t_km2(drained=1.0) == pytest.approx(26.9, abs=0.1)   # continuously flooded, cropped twice
    assert S.rate_t_km2(dry_before=S.DRY_BEFORE) == pytest.approx(4.36, abs=0.01)


@pytest.fixture(scope='module')
def product():
    if not S.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.satoyama_candidates first')
    return json.loads(S.OUT.read_text())


def test_the_candidates_hold_together(product):
    c = product['candidates']
    assert 0 < c['count'] <= product['valleys']['count']
    assert c['flooded_possible_km2'] <= min(c['valley_floor_km2'], S.IRRIGATED_SHARE_MAX * c['km2']) + 1
    assert sum(b['candidates'] for b in product['bands']) == c['count']
    m = product['methane']
    rate = m['paddy_t_km2_yr']['drained_more_than_once_cropped_once']
    assert m['at_the_possible_area_tg_yr'][1] == pytest.approx(c['flooded_possible_km2'] * rate / 1e6, abs=0.01)
