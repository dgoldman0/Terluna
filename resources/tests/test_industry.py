"""The industry screen binds its inputs and its arithmetic stays checkable."""
import hashlib
import json
from pathlib import Path

import numpy as np

from resources import industry
from shared.constants import AU, JULIAN_DAY, JULIAN_YEAR_DAYS, SUN_GM
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(industry.OUT.read_text())
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == industry.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert 'resources/results/first_screen.json' in PRODUCT['producer']['inputs']
    assert 'research/studies/solar_shield_array/results/array_heat.json' in PRODUCT['producer']['inputs']
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_every_source_is_registered_with_its_access():
    sources = json.loads(industry.SOURCES.read_text())['sources']
    for s in sources:
        assert set(s) == {'citation', 'doi_or_url', 'use', 'access'} and all(s.values()), s
    cited = ' '.join(s['citation'] for s in sources)
    for name in ('Krausmann', 'Pauliuk', 'Berglund-Brown', 'SP-413', 'Lodders', 'Meyer 2010', 'Meyers 1994', 'Metal Stocks',
                 'Recycling Rates', 'USGS', 'Chen', 'Roberts', 'DGX', 'Microsoft', 'Durante', 'Energy Institute',
                 'World Population Prospects', 'Mendonca', 'Provisioning branch'):
        assert name in cited, name


def test_annulus_tiles_carry_the_window_tiles_loading():
    c = PRODUCT['fleet_carriage']
    assert np.isclose(c['loading_g_m2_per_unit_pressure'] * c['annulus_pressure_coefficient'], c['annulus_tile_g_m2'], rtol=2e-3)
    assert np.isclose(c['window_tile_g_m2'], 26.0, rtol=2e-3)          # 1 um of titania on 10 um of silica
    assert np.isclose(c['annulus_tile_g_m2'] - c['annulus_film_g_m2'], c['allowance_g_m2'], rtol=5e-3)
    assert 0 < c['store_share_of_allowance'][0] < c['store_share_of_allowance'][1] < 1
    used = PRODUCT['places']['moon_surface']['heat_w_m2_per_tw_used']
    for glow, tw in zip(c['uniform_percent']['glow_w_m2_per_percent_absorbed'], c['uniform_percent']['equivalent_tw_used_on_moon']):
        assert np.isclose(glow / used, tw, rtol=0.01)


def test_tile_exchanges_carry_the_upkeep():
    t = PRODUCT['plant_traffic']
    for per_day, gt in zip(t['exchanges_per_day'], t['upkeep_gt_yr']):
        assert np.isclose(per_day * t['tile_t'] * 1e3 * JULIAN_YEAR_DAYS, gt * 1e12, rtol=0.01)
    one, eight = t['plants']['1'], t['plants']['8']
    assert eight['dv_median_m_s'] < one['dv_median_m_s'] < one['dv_max_m_s']
    for p in (one, eight):
        assert np.isclose(p['sail_days_per_leg'] * JULIAN_DAY * t['sail_acceleration_m_s2'], p['dv_median_m_s'], rtol=0.02)
        for transit, per_day in zip(p['tiles_in_transit'], t['exchanges_per_day']):
            assert np.isclose(transit, 2 * per_day * p['sail_days_per_leg'], rtol=0.02)


def test_heat_reaches_the_moon_less_the_farther_out_it_is_released():
    pl = PRODUCT['places']
    moon = pl['moon_surface']['heat_w_m2_per_tw_used']
    ring = pl['ring_platforms']['heat_w_m2_per_tw_even']
    hubs = pl['sun_earth_l1_l2']['heat_w_m2_per_tw_even']
    assert max(hubs) < ring < min(pl['ring_platforms']['heat_w_m2_per_tw_from_tile_faces']) < moon
    assert np.isclose(pl['ratios']['moon_over_ring'], moon / ring, rtol=0.01)
    # the array heat screen's own placement factor for heat released in orbit
    assert np.isclose(ring, 4.98e-5, rtol=0.01)


def test_arrivals_follow_vis_viva():
    for src in PRODUCT['places']['sources'].values():
        r1, r2 = AU, src['distance_au'] * AU
        a = (r1 + r2) / 2
        v_p = np.sqrt(SUN_GM * (2 / r1 - 1 / a))
        assert np.isclose(src['arrival_v_inf_km_s'] * 1e3 + np.sqrt(SUN_GM / r1), v_p, rtol=2e-3)
        assert np.isclose(src['braking_mj_kg'], 0.5 * (src['arrival_v_inf_km_s'] * 1e3) ** 2 / 1e6, rtol=0.01)
        assert np.isclose(src['transit_years'] * YEAR_S, np.pi * np.sqrt(a ** 3 / SUN_GM), rtol=0.01)


def test_cycle_times_are_inventory_over_flow():
    after = PRODUCT['after_build']
    air = after['air']
    for row in air['losses'].values():
        assert np.isclose(row['kg_s'] * row['cycle_years'] * YEAR_S, air['inventory_kg'], rtol=0.01)
    ox = after['oxygen_to_rock']
    o2 = ox['windows']['years_500_1000']['kg'][1] / ox['windows']['years_500_1000']['air_oxygens'][1]
    assert np.isclose(ox['kg_s'][1] * ox['cycle_years'][0] * YEAR_S, o2, rtol=0.02)
    co2 = after['carbon']
    for years, rate in zip(co2['cycle_years']['weathering'], co2['weathering_gt_co2_yr'][::-1]):
        assert np.isclose(years * rate * 1e12, co2['inventory_kg'], rtol=0.02)
    film = after['film']
    fleet_gt = PRODUCT['after_build']['pre_air_lift']['fleet_gt']
    for r, years in film['fleet_replaced_by_fresh_feed_years'].items():
        assert np.isclose(film['fresh_feed_kg_s'][r][1] * years[0] * YEAR_S, fleet_gt[0] * 1e12, rtol=0.02)


def test_near_horizons_barely_touch_the_air():
    after = PRODUCT['after_build']
    for flow in ('air', 'oxygen_to_rock', 'nitrogen'):
        w = after[flow]['windows']['years_1000_10000']
        share = next(v for k, v in w.items() if k != 'kg')
        assert max(share) < 0.01, flow
    assert min(after['oxygen_to_rock']['windows']['operation_1e9']['air_oxygens']) > 1


def test_habitats_follow_sp413_and_copper_follows_its_grade():
    sp = PRODUCT['people']['sp413']
    for s in PRODUCT['people']['scenarios'].values():
        per_person_t = (sp['shield_t'] + sp['structure_t']) / sp['people']
        assert np.isclose(s['array']['habitat_gt'][0] * 1e12, s['array_people'] * per_person_t * 1e3, rtol=0.01)
        cu = s['moon']['copper']
        for stock, ore in zip(cu['stock_gt'], cu['ore_gt']['ci_chondrite']):
            assert np.isclose(ore, stock / 127e-6, rtol=0.01)
        assert cu['ore_gt']['highland_soil'][0] > cu['ore_gt']['ci_chondrite'][0]


def test_the_build_dwarfs_what_follows():
    op = PRODUCT['operation']
    assert min(op['build_over_after']) > 1e3
    assert op['after_build_total_kg_s'][1] < op['build_kg_s'][0]
    streams = op['after_build_streams_kg_s']
    assert streams['air_escape'][1] < streams['film_fresh_feed'][1] < streams['array_habitats'][1]
    assert np.isclose(sum(v[0] for v in streams.values()), op['after_build_total_kg_s'][0], rtol=0.01)
