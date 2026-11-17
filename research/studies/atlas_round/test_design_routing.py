"""Checks of the design surface's routing: its seas are the globe surface's bodies in order, the neck's spill is put in on
the ground beside the main part's water, and once computed the record holds together (each sea at the level its source
gives it, the product matching its hash) and the product keeps its lakes off the seas."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import near_side_neck as N


def test_the_seas_are_the_globe_surfaces_bodies_in_order():
    assert len(DR.SEAS) == len(GS.BODIES) - 1
    for name, body in GS.BODY_OF.items():
        assert DR.SEAS[body - 1] == name
    assert DR.SEAS[:2] == ('near_side_main', 'near_side_south')            # bodies 1 and 2, the parted near-side sea
    assert set(DR.UPPER_SEAS) <= set(DR.SEAS)


def test_the_spill_leaves_from_the_lowest_ground_beside_the_main_water():
    rows, cols = 64, 128
    h = np.full((rows, cols), 50.0, np.float32)
    main = np.zeros((rows, cols), bool)
    main[10:20, 40:60] = True                     # the main part's water
    h[main] = -10.0
    h[20, 50] = 5.0                               # the lowest ground beside it, south of its middle
    h[25, 50] = 4.0                               # lower ground further off is not beside the water
    centre = (90.0 - 20.5 / 16.0, 50.5 / 16.0)    # the throat's centre, at cell (20, 50)
    assert DR.spill_cell(h.ravel(), main.ravel(), cols, centre) == 20 * cols + 50


@pytest.fixture(scope='module')
def record():
    if not DR.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.design_routing first')
    return json.loads(DR.OUT.read_text())


def test_the_record_holds_together(record):
    s = record['summary']
    seas = s['seas']
    neck = json.loads(N.OUT.read_text())['budgets'][DR.BUDGET]['natural_neck']['steady']
    assert set(seas) == set(DR.SEAS)
    assert seas['near_side_main']['level_above_common_m'] == pytest.approx(neck['main_level_m'], abs=0.02)
    assert seas['near_side_south']['level_above_common_m'] == pytest.approx(neck['south_level_m'], abs=0.02)
    for k in ('south_pole_aitken', 'humboldtianum'):
        assert abs(seas[k]['level_above_common_m']) < 0.02                     # held at the common level
    for k in DR.UPPER_SEAS:
        assert seas[k]['level_above_common_m'] > 700.0                          # full at their rims
    assert s['neck_spill_m3s'] == round(neck['neck_flow_m3s'])
    assert 0.0 < s['closed_lake_share'] <= s['lake_share']
    assert set(record['producer']['files']) == set(DR.CODE)
    if DR.PRODUCT.exists():
        assert DR.sha256(DR.PRODUCT) == record['product']['sha256']


def test_the_product_keeps_its_lakes_off_the_seas(record):
    if not DR.PRODUCT.exists():
        pytest.skip('the routing product is not on the research drive')
    with np.load(DR.PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        sea, lake_cell, receiver = z['sea'], z['lake_cell'], z['receiver']
        lake_level = z['lake_level_m']
    assert meta['seas'] == list(DR.SEAS)
    assert set(np.unique(sea).tolist()) == set(range(-1, len(DR.SEAS)))
    assert (sea[lake_cell] < 0).all()
    assert np.isfinite(lake_level).all()
    assert 0 <= receiver.min() and receiver.max() < sea.size
    assert sea.size == meta['grid']['rows'] * meta['grid']['cols']
