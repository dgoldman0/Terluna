"""The falls: steep steps joined into reaches on a small chain, and, once computed, the record's order and the river
over Moscoviense's pass."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import falls as F


def test_steep_steps_join_into_one_reach_with_its_power():
    # one river of 1,000 m3/s down a row: flat, then three steps of 100 m, then flat into the sea at the end
    cols, n = 6, 6
    receiver = np.array([1, 2, 3, 4, 5, 5])
    surface = np.array([400.0, 400.0, 300.0, 200.0, 100.0, 100.0])
    q = np.full(n, 1000.0)
    water = np.zeros(n, bool)
    water[5] = True
    cells, slope, found = F.reaches(receiver, q, surface, water, cols, ppd=1)
    assert len(found) == 1
    x = found[0]
    assert x['cells'] == [1, 2, 3, 4]
    assert x['drop_m'] == pytest.approx(300.0)
    assert x['power_gw'] == pytest.approx(1000.0 * 1.62 * 1000.0 * 300.0 / 1e9)


@pytest.fixture(scope='module')
def record():
    if not F.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.falls first')
    return json.loads(F.OUT.read_text())


def test_the_record_lists_the_largest_first_and_follows_moscovienses_river(record):
    powers = [x['power_gw'] for x in record['reaches']]
    assert powers == sorted(powers, reverse=True)
    great = record['summary']['moscoviense_river']
    assert great['into'] == 'South Pole–Aitken Sea'
    assert great['falls_m'] > 3000
    assert any(x['on_moscoviense_river'] for x in record['reaches'][:3])
    assert record['product']['sha256'] == F.sha256(F.PRODUCT)
