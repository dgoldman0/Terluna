"""Checks of the design surface's routing: its seas are the globe surface's bodies in order, the neck's spill is put in on
the ground beside the main part's water, and once computed the record holds together (each sea at the level its source
gives it, the product matching its hash) and the product keeps its lakes off the seas."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import lake_rule_balance as LB
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
    crest = json.loads(N.OUT.read_text())['crest']['four_connected_m']
    assert set(seas) == set(DR.SEAS)
    for part in ('main', 'south'):                                             # each part where its loss meets its supply
        within = max(LB.PART_SETTLED_M3S, LB.PART_SETTLED_SHARE * seas[f'near_side_{part}']['receives_m3s'])
        assert abs(s['near_side']['left_out_of_balance_m3s'][part]) < within
    joined = s['near_side_sea_joined_across_its_neck']
    assert joined == (seas['near_side_south']['level_above_common_m'] > crest)
    assert seas['near_side_main']['level_above_common_m'] > crest or s['neck_spill_m3s'] == 0
    for k in ('south_pole_aitken', 'humboldtianum'):
        assert abs(seas[k]['level_above_common_m']) < 0.02                     # held at the common level
    for k in DR.UPPER_SEAS:
        assert seas[k]['level_above_common_m'] > 700.0                          # full at their rims
    weirs = s['weirs']
    assert sum(weirs['released_m3s'].values()) == pytest.approx(weirs['moscoviense_surplus_m3s'], abs=3)
    for k in ('south_pole_aitken', 'humboldtianum'):                           # the weirs meet each held sea's need
        assert abs(seas[k]['receives_m3s'] - seas[k]['loses_m3s']) <= 2 or k in weirs['held_seas_left_short_m3s']
    assert abs(sum(s['budget_m3s'][k] for k in ('runoff_of_dry_ground', 'put_in', 'lakes_rain_less_evaporation'))
               - s['budget_m3s']['reaches_the_seas'] - s['budget_m3s']['kept_by_lakes_not_full']) < 1.0
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
        outlet = z['outlet']
    assert meta['seas'] == list(DR.SEAS)
    assert set(np.unique(sea).tolist()) == set(range(-1, len(DR.SEAS)))
    assert (sea[lake_cell] < 0).all()
    assert np.isfinite(lake_level).all()
    assert 0 <= receiver.min() and receiver.max() < sea.size
    assert sea.size == meta['grid']['rows'] * meta['grid']['cols']
    at = receiver.astype(np.int64)                                             # every course leads through its lakes to a sea
    for _ in range(40):
        at = at[at]
    assert outlet[at].all()
    assert (meta['neck']['spill_cell'] is None) == meta['neck']['joined']
