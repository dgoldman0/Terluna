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
    neck = json.loads(N.OUT.read_text())['budgets']['closed_budget']['natural_neck']['steady']
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


def test_the_surface_is_drawn_from_the_design_routing_as_it_stands(record):
    assert record['inputs'][str(DR.OUT.relative_to(GS.ROOT))] == GS.sha256(DR.OUT)
    routing = json.loads(DR.OUT.read_text())
    assert record['ignored_inputs'][str(DR.PRODUCT)] == routing['product']['sha256']
    assert record['summary']['neck_spill_m3s'] == routing['summary']['neck_spill_m3s']
    assert record['summary']['sea_share'] == pytest.approx(routing['summary']['sea_share'], abs=2e-5)
