"""The computing screen binds its inputs, and its supply, demand and limits stay checkable."""
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np

from provisioning import computing
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(computing.OUT.read_text())
HEAT = json.loads(computing.HEAT.read_text())
MONTHS = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
          'November', 'December', 'Q4'}
CITATION = re.compile(r"([A-Z][A-Za-z0-9\-]+)(?: et al\.| & [A-Z][A-Za-z\-]+|, [A-Z][A-Za-z\-]+ & [A-Z][A-Za-z\-]+)? "
                      r"(?:19|20)\d\d[a-c]?\b")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == computing.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert set(PRODUCT['producer']['inputs']) == {str(p.relative_to(ROOT)) for p in computing.INPUTS}
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_supply_reads_the_array_study_and_the_night_half():
    sup = PRODUCT['supply']
    ledgers = json.loads(computing.LEDGERS.read_text())
    assert np.isclose(sup['night_half_intercepts_pw'] * 1e15, ledgers['energy']['fleet_night_half_intercepts_w'], rtol=1e-3)
    per = sup['per_tw']
    assert np.isclose(per['radiator_km2'], HEAT['computing']['radiator_panel_km2_per_TW'], rtol=1e-3)
    assert np.allclose(per['collector_km2'], sorted(HEAT['computing']['collector_km2_per_TW']), rtol=1e-3)
    assert np.allclose(per['radiator_mt'], HEAT['computing']['radiator_Mt_per_TW'], rtol=1e-3)
    # Today's computing per terawatt is the array study's accelerator figure.
    assert np.isclose(sup['scenarios']['1']['flop_s']['today'], HEAT['computing']['accelerator_flop_per_s_per_TW'], rtol=1e-3)
    # Heat radiated evenly from the ring radius reproduces the array study's placement, and falls as 1/d^2 farther out.
    d = sup['heat_by_distance']
    assert np.isclose(d['ring_radius']['moon_w_m2_per_tw_heat'], HEAT['placement']['moon_W_m2_per_TW_released_in_orbit'], rtol=5e-3)
    ratio = d['ring_radius']['moon_w_m2_per_tw_heat'] / d['at_100000_km']['moon_w_m2_per_tw_heat']
    assert np.isclose(ratio, (100000 / 20000) ** 2, rtol=2e-2)
    # Orbital computing releases its own heat and its collectors' heat: 2.9-4.4 W per W at 20-30% conversion.
    assert np.allclose(per['heat_released_in_orbit_tw'], [1 + 1.9333, 1 + 3.4], rtol=1e-3)


def test_scenarios_scale_linearly_with_power():
    s = PRODUCT['supply']['scenarios']
    for key in ('light_pw', 'mass_gt', 'orbit_moon_w_m2', 'collector_and_radiator_km2', 'if_on_moon_k'):
        assert np.allclose(np.array(s['1000'][key]), 1000 * np.array(s['1'][key]), rtol=1e-2), key
    assert np.allclose(np.array(s['1']['light_pw']) * 1e3, [1 / 0.3, 1 / 0.2], rtol=1e-3)
    per = PRODUCT['supply']['per_tw']
    assert np.allclose(per['mass_mt'], np.array(per['radiator_mt']) + per['collector_mt'] + np.array(per['hardware_mt']), rtol=1e-3)
    # The night half's light is never the bound below 10,000 TW; the fleet's mass is passed before it.
    assert s['10000']['share_of_night_half_light'][1] < 0.1 and s['10000']['share_of_fleet_mass'][0] > 1


def test_radiators_and_the_landauer_floor():
    r = PRODUCT['radiators']
    assert np.isclose(r['330K']['km2_per_tw'], HEAT['computing']['radiator_panel_km2_per_TW'], rtol=1e-3)
    # Area per terawatt grows as T^-4 as radiators cool; the Landauer-limited rate per square metre falls as T^3.
    assert np.isclose(r['150K']['km2_per_tw'] / r['330K']['km2_per_tw'], (330 / 150) ** 4, rtol=1e-3)
    assert np.isclose(r['500K']['per_m2_over_330K'], (500 / 330) ** 3, rtol=1e-3)
    e = PRODUCT['efficiency']
    assert np.isclose(e['landauer_bit_j']['330K'], computing.BOLTZMANN * 330 * math.log(2), rtol=1e-3)
    assert np.isclose(e['cmos_limit_flop_per_j'], 4.7e15 / 16, rtol=1e-3)
    assert np.isclose(e['h100_check'], HEAT['design']['accelerator_flop_per_J'], rtol=2e-2)
    lo, hi = e['landauer_floor_flop_per_j']['330K']
    assert np.isclose(hi / lo, 100, rtol=1e-3)
    assert e['today_flop_per_j'] < e['cmos_limit_flop_per_j'] < lo


def test_rates_and_campaigns_follow_the_efficiencies():
    e = PRODUCT['efficiency']
    night = PRODUCT['supply']['night_half_intercepts_pw'] * 1e15
    for key, r in PRODUCT['comparison'].items():
        if r['kind'] == 'rate':
            assert np.allclose(np.array(r['flop_s']) / e['today_flop_per_j'] / 1e12, r['power_today_tw'], rtol=2e-2), key
            assert np.allclose(np.array(r['power_today_tw']) / e['cmos_over_today'], r['power_cmos_limit_tw'], rtol=2e-2), key
            share = [r['power_today_tw'][0] * 1e12 / 0.3 / night, r['power_today_tw'][1] * 1e12 / 0.2 / night]
            assert np.allclose(share, r['share_of_night_half_light_today'], rtol=2e-2), key
        else:
            days = np.array(r['flop']) / e['today_flop_per_j'] / 1e12 / 86400
            assert np.allclose(days, r['days_at_1_tw_today'], rtol=2e-2), key


def test_weather_costs_a_thousandfold_per_decade_of_resolution():
    w = PRODUCT['weather']
    for coarse, fine in (('1000m', '100m'), ('100m', '10m'), ('10m', '1m')):
        ratio = np.array(w['grids'][fine]['real_time_w']) / np.array(w['grids'][coarse]['real_time_w'])
        assert np.all(ratio >= 1000 * 0.999), (coarse, fine)
    assert np.isclose(w['cosmo_cells'], 4 * math.pi * 6371e3 ** 2 * 0.984 / 930 ** 2 * 60, rtol=1e-3)
    # COSMO's 596 MWh per simulated year over its cells and 6-second steps.
    steps = 365.25 * 86400 / 6
    assert np.isclose(w['j_per_cell_step_p100'], 596 * 3.6e9 / (w['cosmo_cells'] * steps), rtol=1e-3)


def test_limits_and_minds_reproduce_the_literature():
    lim = PRODUCT['limits']
    assert np.isclose(lim['margolus_levitin_ops_per_s_per_kg'], 5.4258e50, rtol=1e-3)     # Lloyd 2000
    assert np.isclose(lim['bremermann_bits_per_s_per_kg'], 1.36e50, rtol=1e-2)
    assert 187 < lim['schneier']['counter_bits'] < 188                                     # Schneier 1996
    c = PRODUCT['comparison']
    assert np.allclose(c['ancestor_simulation']['product_of_terms'], [1.5e34, 1.5e37], rtol=1e-2)
    # A brain at the roadmap's spiking level is 1e18 FLOP/s, 0.71 MW at today's 1.4e12 FLOP/J.
    assert np.isclose(c['wbe_level_4']['power_today_tw'][0], 1e18 / 1.4e12 / 1e12, rtol=1e-2)


def test_latency_places_conversation_within_reach_of_the_array():
    lat = PRODUCT['latency']
    screen = json.loads(computing.SCREEN.read_text())
    assert np.isclose(lat['moon_array_round_trip_s'], screen['latency']['round_trip_s_to_ring_radius'], rtol=1e-2)
    assert lat['moon_array_round_trip_s'] < lat['conversation_gap_s'] < lat['earth_moon_round_trip_s']


def test_every_cited_source_is_listed_and_every_listed_source_cited():
    sources = json.loads(computing.SOURCES.read_text())['sources']
    for s in sources:
        assert set(s) == {'citation', 'doi_or_url', 'use', 'access'} and all(s.values()), s
    listed = ' '.join(s['citation'] for s in sources)
    code = Path(computing.__file__).read_text()
    for name in set(CITATION.findall(code)) - MONTHS:
        assert name in listed, name
    note = (computing.HERE / 'computing.md').read_text() + code
    for s in sources:
        first = re.split(r'[ ,&]', s['citation'])[0]
        assert first in note, s['citation']
