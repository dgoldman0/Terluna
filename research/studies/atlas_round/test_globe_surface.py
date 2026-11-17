"""Checks of the design surface the globe draws: the spill is put in on the ground beside the main part's water, and the
record holds together once computed (each sea at the level its source gives it, the shares adding up, the neck's river
reaching the southern seas, the grid on the research drive matching its hash)."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import globe_surface as GS


def test_the_spill_leaves_from_the_lowest_ground_beside_the_main_water():
    rows, cols = 64, 128
    h = np.full((rows, cols), 50.0, np.float32)
    main = np.zeros((rows, cols), bool)
    main[10:20, 40:60] = True                     # the main part's water
    h[main] = -10.0
    h[20, 50] = 5.0                               # the lowest ground beside it, south of its middle
    h[25, 50] = 4.0                               # lower ground further off is not beside the water
    centre = (90.0 - 20.5 / 16.0, 50.5 / 16.0)    # the throat's centre, at cell (20, 50)
    assert GS.spill_cell(h.ravel(), main.ravel(), cols, centre) == 20 * cols + 50


@pytest.fixture(scope='module')
def record():
    if not GS.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.globe_surface first')
    return json.loads(GS.OUT.read_text())


def test_the_record_holds_together(record):
    s = record['summary']
    seas = s['seas']
    neck = json.loads(GS.N.OUT.read_text())['budgets']['closed_budget']['natural_neck']['steady']
    assert seas['near-side sea']['level_m'] == pytest.approx(neck['main_level_m'], abs=0.02)
    assert seas['southern near-side seas']['level_m'] == pytest.approx(neck['south_level_m'], abs=0.02)
    for k in ('South Pole–Aitken Sea', 'Mare Humboldtianum'):
        assert abs(seas[k]['level_m']) < 0.02                                # held at the common level
    for k in ('Mare Moscoviense', 'Mare Orientale', 'Smythii–Marginis Sea'):
        assert seas[k]['level_m'] > 700.0                                     # full at their rims
    assert s['water_share'] == pytest.approx(s['sea_share'] + s['lake_share'], abs=2e-4)
    assert sum(v['share'] for v in seas.values()) == pytest.approx(s['sea_share'], abs=2e-4)
    assert s['neck_river']['reaches'] == 'southern near-side seas'
    assert s['neck_spill_m3s'] == round(neck['neck_flow_m3s'])
    product = record['product']
    if GS.PRODUCT.exists():
        assert GS.sha256(GS.PRODUCT) == product['sha256']
