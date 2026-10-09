"""The provisioning screen binds its inputs and its arithmetic stays checkable."""
import hashlib
import json
from pathlib import Path

import numpy as np

from provisioning import screen
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(screen.OUT.read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == screen.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_fleet_infrared_follows_the_view_factor():
    heat = PRODUCT['heat']
    area = 4 * np.pi * screen.MOON_RADIUS ** 2
    shares = (heat['moon_facing_share']['deep_stack'], heat['moon_facing_share']['one_layer'])
    for absorbed, share, view, w in zip(heat['tiles_absorbed_pw'], shares, heat['view_factor'], heat['fleet_infrared_w_m2']):
        assert np.isclose(absorbed * 1e15 * share * view / area, w, rtol=0.01)
    # The face turned to the Moon is the same all orbit; the films emit a little more there than sunward.
    assert 0.5 < shares[1] < 0.65


def test_two_layers_split_their_heat_between_the_even_split_and_one_layer():
    assert np.isclose(screen.two_layers(0.4, 0.4, 0.4)[0], 0.5)                # symmetric faces split evenly
    assert 0.5 < screen.two_layers(0.37, 0.44, 0.44)[0] < 0.44 / (0.44 + 0.37)
    mirror = PRODUCT['heat']['heat_mirror']
    assert 0.5 < mirror['two_layers_as_is'] < PRODUCT['heat']['moon_facing_share']['one_layer']
    for case in (v for k, v in mirror.items() if k.startswith('mirror_')):
        # A mirror on every tile traps the hidden layer between two weak emitters; on the nearest layer alone it does better.
        assert case['one_layer_share'] < case['nearest_layer_only']['share'] < case['every_tile']['share'] < mirror['two_layers_as_is']
        assert case['every_tile']['hottest_layer_K'] > case['nearest_layer_only']['hottest_layer_K']


def test_human_heat_is_the_array_studys_figure():
    array_heat = json.loads(screen.HEAT.read_text())
    assert np.isclose(PRODUCT['heat']['human_use']['w_m2_per_tw'], array_heat['placement']['moon_W_m2_per_TW_used_on_moon'], rtol=0.01)


def test_food_land_follows_the_harvest():
    for site in PRODUCT['food']['sites'].values():
        kcal = site['harvest_kg_per_cycle'] * 1000 * screen.FRUIT_KCAL_PER_G / screen.SYNODIC_MONTH_DAYS
        assert np.isclose(site['kcal_m2_day'], kcal, rtol=0.01)
        assert np.isclose(site['m2_per_person'], screen.PERSON_KCAL_PER_DAY / kcal, rtol=0.02)


def test_night_energy_is_the_metropolis_for_half_a_month():
    night = PRODUCT['night']
    assert np.isclose(night['metropolis_night_twh'], 200e9 * night['night_hours'] / 1e12, rtol=0.01)


def test_tidal_power_scales_with_gravity_and_frequency():
    ratio = PRODUCT['waters']['tidal_power_per_area_moon_over_earth']
    assert 1 / 400 < ratio < 1 / 250


def test_transit_scales_with_gravity():
    t = PRODUCT['transit']
    for e, m in zip(t['standing_acceleration_m_s2']['earth'], t['standing_acceleration_m_s2']['moon']):
        assert np.isclose(m / e, screen.MOON_SURFACE_GRAVITY / screen.STANDARD_GRAVITY, rtol=0.01)


def test_fusion_fuel_energy_per_deuteron():
    f = PRODUCT['fusion_fuel']
    joules_per_kg_d = f['energy_j'] / f['deuterium_kg']
    assert np.isclose(joules_per_kg_d, 7.2 * 1.602176634e-13 * screen.AVOGADRO / 2.014e-3, rtol=0.02)


def test_industry_ledger_follows_the_glow_share():
    heat = PRODUCT['heat']
    ind = heat['industry']
    area = 4 * np.pi * screen.MOON_RADIUS ** 2
    for w, a, s in zip(heat['fleet_infrared_w_m2'], heat['tiles_absorbed_pw'], ind['glow_share_of_absorbed']):
        assert np.isclose(w * area / (a * 1e15), s, rtol=0.01)
    # Halving the glow through computing means taking half the films' heat; today's data centres are a millionth of it.
    for c, a in zip(ind['computing_tw_to_halve_glow'], heat['tiles_absorbed_pw']):
        assert np.isclose(c, a * 1e3 / 2, rtol=0.01)
    assert ind['data_centres_tw'] / min(ind['computing_tw_to_halve_glow']) < 1e-5
    # An absorbing dimmer on the window rings returns the glow share of all the light those rings take in.
    for r, g, b in zip(ind['window_rings_intercept_per_moon_sunlight'], ind['glow_share_of_absorbed'], ind['absorbing_dimmer_gives_back']):
        assert np.isclose(r * g, b, rtol=0.02)
    assert ind['moon_w_m2_from_industry_tw_released_in_orbit']['1000'] < 0.1 < ind['moon_w_m2_from_industry_tw_used_on_moon']['100']
    # Computing in orbit brings its collectors' heat: 2.9-4.4 W in all for each watt computed.
    lo, hi = ind['moon_w_m2_from_computing_tw_near_fleet_with_collectors']['1000']
    assert np.isclose(lo, ind['moon_w_m2_from_industry_tw_released_in_orbit']['1000'] * (1 + ind['collectors_heat_per_watt'][0]), rtol=0.02)
    assert 0.14 < lo < hi < 0.23
