"""Checks of the satoyama map on the final surface: its seas are the globe surface's bodies in order, and once computed
it reads the design routing as that routing stands, each sea stands where the globe surface has it, and the candidates
hold together as the first map's do."""
import json

import pytest

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import satoyama_candidates as SC
from research.studies.atlas_round import satoyama_final as S


def test_the_seas_are_the_globe_surfaces_bodies_in_order():
    assert len(S.SEAS) == len(GS.BODIES) - 1
    for name, body in GS.BODY_OF.items():
        assert S.SEAS[body - 1] == name
    assert S.SEAS[:2] == ('near_side_main', 'near_side_south')            # bodies 1 and 2, the parted near-side sea


@pytest.fixture(scope='module')
def product():
    if not S.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.satoyama_final first')
    return json.loads(S.OUT.read_text())


def test_the_map_reads_the_design_routing_as_it_stands(product):
    routing = product['routing']
    assert routing['record'] == str(DR.OUT.relative_to(S.ROOT))
    assert product['inputs'][routing['record']] == S.sha256(DR.OUT)
    assert product['ignored_inputs'][str(DR.PRODUCT)] == json.loads(DR.OUT.read_text())['product']['sha256']
    seas = json.loads(GS.OUT.read_text())['summary']['seas']
    for name, level in routing['levels_above_the_common_level_m'].items():
        assert level == pytest.approx(seas[GS.BODIES[S.SEAS.index(name) + 1]]['level_m'], abs=0.02)


def test_the_candidates_hold_together(product):
    c = product['candidates']
    assert 0 < c['count'] <= product['valleys']['count']
    assert c['flooded_possible_km2'] <= min(c['valley_floor_km2'], SC.IRRIGATED_SHARE_MAX * c['km2']) + 1
    assert sum(b['candidates'] for b in product['bands']) == c['count']
    m = product['methane']
    rate = m['paddy_t_km2_yr']['drained_more_than_once_cropped_once']
    assert m['at_the_possible_area_tg_yr'][1] == pytest.approx(c['flooded_possible_km2'] * rate / 1e6, abs=0.01)
    first = product['first_map']
    assert first['sha256'] == S.sha256(SC.OUT)
    assert first['candidates'] == json.loads(SC.OUT.read_text())['candidates']['count']
