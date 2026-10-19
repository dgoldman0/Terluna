"""Scenes read the arrivals product, keep one water budget, and carry the target's air."""
import numpy as np
import pytest

import scenes
from atmosphere.electricity import thunder

PRODUCT = scenes.load_arrivals()


@pytest.mark.parametrize('surface', ('concrete', 'puddles', 'grass-leaves'))
def test_every_scene_delivers_its_rain_once(surface):
    for target in scenes.TARGETS:
        for world in scenes.WORLDS:
            s = scenes.scene(PRODUCT, target, world, surface, rain_mm_h=7.0)
            assert abs(s.ledger()['balance_mm_h']) < 1e-12
            assert (s.drips is not None) == (surface == 'grass-leaves')
            if s.drips is not None:
                assert s.drips.rain_mm_h == pytest.approx(7.0 * scenes.DEFAULT_DRIPS['fraction'])
    heavy = scenes.scene(PRODUCT, 'lock', 'moon', surface, tail_fraction=0.02)
    assert heavy.tail.rain_mm_h == pytest.approx(0.02 * scenes.EARTH_REFERENCE_MM_H)
    assert abs(heavy.ledger()['balance_mm_h']) < 1e-12


def test_the_direct_rain_is_the_products_arriving_spectrum():
    side = PRODUCT['targets']['resort']['moon']
    row = next(r for r in side['spectra'] if (r['dm_mm'], r['shape_mu']) == (1.5, 3.0))
    p = scenes.direct_rain(PRODUCT, 'resort', 'moon', 2.0)
    assert p.arrivals_m2_s == pytest.approx(2.0 * row['arrivals_m2_s_per_mm_h'], rel=1e-12)
    survivors = [d for d, s in zip(side['classes']['surface_diameter_mm'], side['classes']['status']) if s == 'surface']
    assert np.allclose(np.sort(p.diameter_m * 1e3), np.sort(survivors))
    with pytest.raises(ValueError):
        scenes.scene(PRODUCT, 'resort', 'moon', 'concrete', rain_mm_h=0.0)


def test_the_moons_drops_arrive_slower_and_fewer():
    for target in scenes.TARGETS:
        moon, earth = (scenes.direct_rain(PRODUCT, target, w, 5.0).describe() for w in ('moon', 'earth'))
        assert moon['arrivals_m2_s'] < earth['arrivals_m2_s']
        assert moon['mean_speed_m_s_by_number'] < 0.5 * earth['mean_speed_m_s_by_number']
        assert moon['kinetic_power_w_m2'] < 0.25 * earth['kinetic_power_w_m2']


def test_the_air_carries_the_thunder_models_sound_speed():
    for target in scenes.TARGETS:
        for world in scenes.WORLDS:
            air = scenes.air(PRODUCT, target, world)
            assert 340 < air['sound_speed_m_s'] < 360
            dry = thunder.sound_speed(air['temperature_k'], 0.0, air['dry_molar_mass_kg_mol'])
            assert dry < air['sound_speed_m_s'] < dry * 1.01                     # vapour lightens the air
        moon, earth = (scenes.air(PRODUCT, target, w) for w in ('moon', 'earth'))
        assert moon['x_o2'] < earth['x_o2'] and moon['pressure_pa'] > earth['pressure_pa']


def test_the_moons_rates_are_the_gcms():
    rates = scenes.load_rates()
    for target in scenes.TARGETS:
        r = scenes.gcm_rates(rates, target)
        wet = rates['targets'][target]['six_and_a_half_hour_means']['wet']
        assert r == dict(median=wet['p50'], heavy=wet['p99'])
        assert rates['wet_threshold_mm_h'] < r['median'] < r['heavy'] < 20
    lock, hot = (scenes.gcm_rates(rates, t) for t in ('lock', 'hot_world'))
    assert hot['median'] > lock['median'] and hot['heavy'] > lock['heavy']
