"""The population screen binds its inputs and its arithmetic stays checkable."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from provisioning import population as pop
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(pop.OUT.read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == pop.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_every_source_entry_is_complete():
    sources = json.loads(pop.SOURCES.read_text())['sources']
    for s in sources:
        assert set(s) == {'citation', 'doi_or_url', 'use', 'access'}, s['citation']
        assert all(isinstance(v, str) and v for v in s.values()), s['citation']


def test_growth_paths_reach_their_targets():
    for key, path in PRODUCT['horizons']['paths_from_2100'].items():
        target, year = key.split('_by_')
        reached = pop.POP_2100 * math.exp(path['percent_per_year'] / 100 * (int(year) - 2100))
        assert np.isclose(reached, float(target) * 1e9, rtol=1e-3), key
        assert np.isclose(path['doubling_years'] * path['percent_per_year'] / 100, math.log(2), rtol=1e-2), key
    history = PRODUCT['horizons']['history_percent_per_year']
    # Reaching 20-30 billion by the year 3000 takes growth of the pre-industrial world's order.
    by_3000 = [PRODUCT['horizons']['paths_from_2100'][f'{t:g}_by_3000']['percent_per_year'] for t in pop.TARGETS_BILLION]
    assert history['year0_to_1750'] < min(by_3000) < max(by_3000) < history['year1500_to_1750']


def test_heat_ceilings_follow_the_moons_kelvin_per_terawatt():
    heat = PRODUCT['heat']
    array_heat = json.loads(pop.HEAT.read_text())
    assert np.isclose(heat['moon_w_m2_per_tw'], array_heat['placement']['moon_W_m2_per_TW_used_on_moon'], rtol=0.01)
    lo, hi = heat['dimming_percent_per_tw']
    assert np.isclose(hi * heat['tw_per_percent_dimming'][0], 1.0, rtol=0.01)
    assert np.isclose(lo * heat['tw_per_percent_dimming'][1], 1.0, rtol=0.01)
    for band in heat['moon'].values():
        for row in band.values():
            h = row['kw_per_person']
            for tol, (b_lo, b_hi) in row['billion_at_tolerance'].items():
                t = float(tol.rstrip('K'))
                assert np.isclose(b_lo, t * heat['moon_tw_per_k'][0] / h, rtol=0.02)
                assert np.isclose(b_hi, t * heat['moon_tw_per_k'][1] / h, rtol=0.02)
            assert np.isclose(row['city_w_m2_at_40000_km2'], pop.METROPOLIS_DENSITY_KM2 * h / 1e3, rtol=0.01)
    # A TW heats the Moon about 13.5 times as much per square metre as Earth: their areas' ratio.
    assert np.isclose(heat['moon_over_earth_per_tw'], (pop.EARTH_RADIUS / pop.MOON_RADIUS) ** 2, rtol=0.01)
    # Moving industry to orbit removes more heat than beaming the rest of the power down.
    f = heat['placement_factor']
    assert f['industry_in_orbit_power_from_orbit'] < f['industry_in_orbit_thermal'] < f['all_on_moon_thermal']
    assert f['all_on_moon_thermal'] - f['industry_in_orbit_thermal'] > f['industry_in_orbit_thermal'] - f['industry_in_orbit_power_from_orbit']


def test_water_land_and_food_identities():
    water, food = PRODUCT['water'], PRODUCT['food']
    drainage = json.loads(pop.DRAINAGE.read_text())
    atlas = json.loads(pop.ATLAS.read_text())
    runoff = drainage['rivers']['discharge_to_sea_m3s'] * pop.HOURS_PER_YEAR * 3600 / 1e9
    assert np.isclose(water['moon_runoff_km3'], runoff, rtol=1e-3)
    for limit, billion in water['moon_billion_at'].items():
        assert np.isclose(billion * 1e9 * pop.FALKENMARK_M3[limit], runoff * 1e9, rtol=0.01)
    dry = 4 * math.pi * (pop.MOON_RADIUS / 1e3) ** 2 * (atlas['main_land_share'] + atlas['island_share'] - drainage['lakes']['area_share'])
    assert np.isclose(food['moon_dry_land_km2'], dry, rtol=1e-3)
    # The Moon's runoff is about a sixth of Earth's, close to its share of Earth's ice-free land.
    assert 0.12 < water['moon_share_of_earth'] < 0.2
    assert 0.12 < food['moon_dry_land_km2'] / pop.ICE_FREE_KM2 < 0.2


def test_array_spin_and_shield():
    arr = PRODUCT['array']
    for rpm, radii in arr['spin'].items():
        omega = float(rpm.rstrip('rpm')) * 2 * math.pi / 60
        assert np.isclose(radii['earth_g_radius_m'], pop.STANDARD_GRAVITY / omega ** 2, atol=1)
        assert np.isclose(radii['moon_g_radius_m'], pop.MOON_SURFACE_GRAVITY / omega ** 2, atol=0.1)
    assert np.isclose(arr['shield_t_per_person'], pop.TORUS_SHIELD_T / pop.TORUS_PEOPLE)
    assert arr['habitats']['1_billion']['shield_gt'][1] == round(arr['shield_t_per_person'], 1)


def test_commute_transfers():
    c = PRODUCT['commute']
    r0, r1 = pop.MOON_RADIUS, c['ring_radius_km'] * 1e3
    a = (r0 + r1) / 2
    assert np.isclose(c['transfers']['hohmann']['hours'], math.pi * math.sqrt(a ** 3 / pop.MOON_GM) / 3600, atol=0.05)
    # Faster transfers cost more; every one costs more than the vacuum minimum per kilogram.
    t = c['transfers']
    assert t['parabolic']['hours'] < t['apolune_twice_ring']['hours'] < t['hohmann']['hours']
    assert t['hohmann']['delta_v_vacuum_m_s'] < t['apolune_twice_ring']['delta_v_vacuum_m_s'] < t['parabolic']['delta_v_vacuum_m_s']
    assert np.isclose(c['minimum_kwh_per_kg'], (pop.MOON_GM / r0 - pop.MOON_GM / (2 * r1)) / 3.6e6, rtol=1e-3)
    mwh = t['hohmann']['round_trip_mwh_per_traveller'][0][0]
    assert mwh * 1e3 > c['minimum_kwh_per_traveller'][0]
    weekly = t['hohmann']['continuous_kw_per_commuter']['weekly'][0]
    assert np.isclose(weekly, pop.ROTATIONS['weekly'] * mwh * 1e3 / pop.HOURS_PER_YEAR, rtol=0.01)


def test_earth_service_geometry():
    e = PRODUCT['earth_service']
    # The 1978 reference system sits near tau = 2, so the same tau sizes its 1 km transmitter.
    assert 1.8 < e['tau_1978'] < 2.2
    assert np.isclose(e['transmitter_diameter']['2.45GHz_geo_km'], pop.SPS_1978['transmitter_km'], rtol=0.1)
    # The transmitter grows with distance at a fixed rectenna.
    ratio = e['transmitter_diameter']['2.45GHz_moon_km'] / e['transmitter_diameter']['2.45GHz_geo_km']
    assert np.isclose(ratio, pop.EARTH_MOON_DISTANCE / (e['geo_altitude_km'] * 1e3), rtol=0.01)
    assert np.isclose(e['rectenna_km2_per_tw'] * e['rectenna_mean_w_m2'] * 1e6, 1e12, rtol=0.01)
    assert np.isclose(e['light_time_s'], pop.EARTH_MOON_DISTANCE / pop.SPEED_OF_LIGHT, rtol=1e-3)


def test_scenarios_add_up():
    band = PRODUCT['affluence']['band_primary_kw']
    for key, s in PRODUCT['scenarios'].items():
        b = s['billion']
        assert np.isclose(b['earth'] + b['moon'] + b['array'] + b['elsewhere'], b['total'])
        assert pop.TARGETS_BILLION[0] <= b['total'] <= pop.TARGETS_BILLION[-1]
        for name, row in s['by_energy'].items():
            assert np.isclose(row['civilization_tw'], b['total'] * band[name], atol=0.6)
            assert np.isclose(row['moon_tw']['all_on_moon_thermal'], b['moon'] * band[name], rtol=0.01)
        assert np.isclose(s['moon_water_m3_per_person'] * b['moon'] * 1e9, PRODUCT['water']['moon_runoff_km3'] * 1e9, rtol=0.01)
    # The binding lists are ordered from the first limit reached.
    for body in ('moon', 'earth'):
        firsts = [row['billion'][0] for row in PRODUCT['binds_first'][body]]
        assert firsts == sorted(firsts)
