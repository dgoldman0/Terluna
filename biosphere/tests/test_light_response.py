"""Checks on the canopy's light response: its binding, and the arithmetic and orderings it rests on."""
import hashlib
import json

import numpy as np
import pytest

from biosphere.canopy import light_response as lr
from biosphere.canopy import plant


@pytest.fixture(scope='module')
def product():
    path = lr.RESULTS / 'light_response.json'
    if not path.is_file():
        pytest.skip('results/light_response.json has not been generated')
    return json.loads(path.read_text())


def test_the_product_binds_its_producer_and_the_plant_model(product):
    assert product['schema'] == lr.SCHEMA
    for path, digest in product['producer']['files'].items():
        assert hashlib.sha256((lr.ROOT / path).read_bytes()).hexdigest()[:16] == digest, path
    plant_json = (lr.RESULTS / 'plant.json').read_bytes()
    assert hashlib.sha256(plant_json).hexdigest()[:16] == product['producer']['inputs']['plant_sha256']


def test_the_undimmed_stand_is_the_plant_models(product):
    stored = json.loads((lr.RESULTS / 'plant.json').read_text())['results'][lr.SIDE]
    for band, out in product['bands'].items():
        for name in lr.STRATEGIES:
            assert out[name]['base']['growth'] == pytest.approx(stored[band][name]['growth'], rel=1e-3)
            assert out[name]['base']['store'] == pytest.approx(stored[band][name]['store'], rel=1e-3)


def test_elasticities_follow_from_the_dimmed_runs(product):
    for out in product['bands'].values():
        for name in lr.STRATEGIES:
            base = out[name]['base']
            for d, dim in out[name]['dimmed'].items():
                for k in ('gpp', 'growth', 'store'):
                    expected = (base[k] - dim[k]) / base[k] / float(d)
                    assert out[name]['elasticity'][d][k] == pytest.approx(expected, rel=2e-2, abs=1e-3)


def test_growth_follows_the_light_more_than_photosynthesis_and_the_store_grows(product):
    """Upkeep does not shrink with the light, so growth loses a larger share than photosynthesis; the twilight that
    shortens the night is dimmed too, so the store grows."""
    for out in product['bands'].values():
        for name in lr.STRATEGIES:
            e = out[name]['elasticity']['0.05']
            assert 0.0 < e['gpp'] < e['growth']
            assert e['store'] < 0.0
    # poleward, where the Sun stays low, the light matters more
    gpp = [product['bands'][b]['earth_like']['elasticity']['0.05']['gpp'] for b in ('0', '30', '60')]
    assert np.all(np.diff(gpp) > 0)


def test_deeper_canopies_and_lower_suns_follow_the_light_more(product):
    for out in product['bands'].values():
        by_lai = [out['by_lai'][f'{lai:g}']['elasticity']['0.05'] for lai in plant.LAIS]
        assert np.all(np.diff(by_lai) > 0)
        parts = out['by_part_at_5_percent']
        lit = [p for p in parts.values() if p['share_of_photosynthesis'] > 0]
        assert sum(p['share_of_photosynthesis'] for p in parts.values()) == pytest.approx(1.0, abs=1e-3)
        assert sum(p['share_of_the_loss'] for p in parts.values()) == pytest.approx(1.0, abs=1e-3)
        assert np.all(np.diff([p['elasticity'] for p in lit]) < 0)
    s = product['by_sun_height']
    up = [e for h, e in zip(s['elevation_deg'], s['elasticity']['0.05']) if h >= 0]
    assert np.all(np.diff(up) < 0)
    assert all(0.0 < e <= 1.0 + 1e-6 for e in s['elasticity']['0.05'])
