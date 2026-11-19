"""Checks of the design surface the globe draws, once computed: it is drawn from the design routing as that routing
stands, and the record holds together (each sea at the level its source gives it, the shares adding up, the neck's river
reaching the southern seas, the grid on the research drive matching its hash). The spill's placement is checked with
the routing (test_design_routing.py)."""
import json

import pytest

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import near_side_neck as N


@pytest.fixture(scope='module')
def record():
    if not GS.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.globe_surface first')
    return json.loads(GS.OUT.read_text())


def test_the_record_holds_together(record):
    s = record['summary']
    seas = s['seas']
    routing = json.loads(DR.OUT.read_text())['summary']
    assert seas['near-side sea']['level_m'] == pytest.approx(routing['seas']['near_side_main']['level_above_common_m'], abs=0.02)
    assert seas['southern near-side seas']['level_m'] == pytest.approx(
        routing['seas']['near_side_south']['level_above_common_m'], abs=0.02)
    for k in ('South Pole–Aitken Sea', 'Mare Humboldtianum'):
        assert abs(seas[k]['level_m']) < 0.02                                # held at the common level
    for k in ('Mare Moscoviense', 'Mare Orientale', 'Smythii–Marginis Sea'):
        assert seas[k]['level_m'] > 700.0                                     # full at their rims
    assert s['water_share'] == pytest.approx(s['sea_share'] + s['lake_share'], abs=2e-4)
    assert sum(v['share'] for v in seas.values()) == pytest.approx(s['sea_share'], abs=2e-4)
    if routing['near_side_sea_joined_across_its_neck']:                      # one sea across the neck: no river over it
        assert s['neck_river'] is None
    else:
        assert s['neck_river']['reaches'] == 'southern near-side seas'
    assert s['neck_spill_m3s'] == routing['neck_spill_m3s']
    product = record['product']
    if GS.PRODUCT.exists():
        assert GS.sha256(GS.PRODUCT) == product['sha256']


def test_the_surface_is_drawn_from_the_design_routing_as_it_stands(record):
    assert record['inputs'][str(DR.OUT.relative_to(GS.ROOT))] == GS.sha256(DR.OUT)
    routing = json.loads(DR.OUT.read_text())
    assert record['ignored_inputs'][str(DR.PRODUCT)] == routing['product']['sha256']
    assert record['summary']['neck_spill_m3s'] == routing['summary']['neck_spill_m3s']
    assert record['summary']['sea_share'] == pytest.approx(routing['summary']['sea_share'], abs=2e-5)


def test_the_works_cuts_are_cut_into_the_ground_and_carry_their_flow(record):
    import numpy as np
    from research.studies.atlas_round import design_cuts as DC
    if not (DC.OUT.exists() and DC.PRODUCT.exists() and GS.PRODUCT.exists()):
        pytest.skip('run design_cuts and globe_surface first')
    assert record['inputs'][str(DC.OUT.relative_to(GS.ROOT))] == GS.sha256(DC.OUT)
    with np.load(DC.PRODUCT, allow_pickle=False) as z:
        cut, level, flow = z['cut_cell'], z['cut_level_m'], z['cut_flow_m3s']
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        h = z['height_m'].ravel()
        channel = dict(zip(z['channel_cell'].tolist(), z['channel_discharge_m3s'].tolist()))
        water = dict(zip(z['water_cell'].tolist(), z['water_level_m'].tolist()))
    assert cut.size > 0
    assert np.all(h[cut] <= level + 1e-3)                                   # the ground taken down to the water
    for c, lev, f in zip(cut.tolist(), level.tolist(), flow.tolist()):       # the flow runs through, or a lake it drained
        assert channel.get(c, 0.0) >= f - 1e-3 or water.get(c, np.inf) == pytest.approx(lev, abs=0.01)
    assert record['summary']['cuts']['deepest_m'] == pytest.approx(json.loads(DC.OUT.read_text())['summary']['deepest_m'],
                                                                   abs=0.2)
    pieces = record['summary']['cuts']['pieces']
    assert sum(x['cells'] for x in pieces) == cut.size and pieces[0]['deepest_m'] == record['summary']['cuts']['deepest_m']
