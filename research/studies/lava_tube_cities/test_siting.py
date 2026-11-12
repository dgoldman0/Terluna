"""Checks of spatial identity, scale and reproducible source binding."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest

from shared.constants import MOON_RADIUS
from shared.provenance import constants_changed
from . import siting


def test_containing_cell_and_longitude_seam():
    atlas = {'lat_deg': np.array([67.5, 22.5, -22.5, -67.5]),
             'lon_deg': np.array([45., 135., 225., 315.]),
             'height_m': np.arange(16).reshape(4, 4),
             'water_label': np.zeros((4, 4), int),
             'metadata': {'sea_level_m': 0}}
    assert siting.atlas_point(atlas, 0, 0)['cell_index'] == [2, 0]
    assert siting.atlas_point(atlas, 0, 360) == siting.atlas_point(atlas, 0, 0)
    assert siting.atlas_point(atlas, 0, -1)['cell_index'] == [2, 3]
    assert siting.atlas_point(atlas, -90, 0)['cell_index'] == [3, 0]
    with pytest.raises(ValueError):
        siting.atlas_point(atlas, 91, 0)


@pytest.fixture(scope='module')
def result():
    spec = json.loads((siting.HERE / 'siting_inputs.json').read_text())
    path = siting.ROOT / spec['inputs']['atlas']['path']
    if not path.exists():
        pytest.skip('Ignored 28% atlas grid is not restored')
    return siting.run()


def test_population_and_floor_scope(result):
    ref = result['population_scales'][0]
    assert ref['residents_per_city'] == 50000
    assert ref['residents_four_cities'] == 200000
    assert ref['residents_cities_plus_towns'] == 250000
    assert ref['ordinary_floor_m2_per_city'] == 3500000
    assert result['illustrative_geometry']['ordinary_floor_m2'] == 3600000
    assert result['illustrative_geometry']['residents_at_inherited_floor_allowance'] == pytest.approx(3600000 / 70)


def test_water_screen_has_correct_sign_and_scale(result):
    rows = {r['key']: r for r in result['candidates']}
    for key in ('tranquillitatis', 'marius'):
        assert rows[key]['coarse_water_screen']['ground_above_common_sea_m'] > 0
        assert not rows[key]['coarse_water_screen']['common_level_water']
    for key in ('mairan', 'schroter_extension', 'ingenii', 'mortis'):
        q = rows[key]['coarse_water_screen']
        assert q['common_level_water']
        assert q['water_above_cell_ground_m'] == -q['ground_above_common_sea_m']
    assert rows['marius']['coarse_water_screen']['north_south_cell_km'] == pytest.approx(math.pi * MOON_RADIUS / 720000)


def test_anomaly_is_not_the_historical_valley_centre(result):
    row = next(r for r in result['candidates'] if r['key'] == 'schroter_extension')
    old = row['historical_point_comparison']
    assert [row['latitude_deg'], row['longitude_east_deg']] == [24, 306]
    assert not old['same_coordinate']
    assert old['ground_above_sea_m'] > 0
    assert row['coarse_water_screen']['ground_above_common_sea_m'] < 0


def test_input_hash_mismatch_is_not_silently_accepted():
    spec = json.loads((siting.HERE / 'siting_inputs.json').read_text())
    wrong = deepcopy(spec)
    wrong['inputs']['town_design']['sha256'] = '0' * 64
    with pytest.raises(ValueError, match='hash mismatch'):
        siting.load_bound(wrong, 'town_design')


def test_sources_resolve_and_proposed_dimensions_remain_separate(result):
    sources = json.loads((siting.HERE / 'sources.json').read_text())
    known = {s['key'] for s in sources['sources']}
    assert all(set(c['sources']) <= known for c in result['candidates'])
    assert all(not c['demonstrated_city_volume'] for c in result['candidates'])
    assert result['assessment']['selected_city_sites'] == []
    assert 'not an observed cavity' in result['illustrative_geometry']['status']


def test_committed_result_reproduces_and_binds_constants(result):
    saved = json.loads((siting.HERE / 'results/siting.json').read_text())
    assert not constants_changed(saved['producer']['constants'])
    assert saved == result
