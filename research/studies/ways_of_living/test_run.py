"""Conservation, units, physical limits and bindings of the ways-of-living product."""
import json
import math

import numpy as np
import pytest

from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_changed
from research.studies.ways_of_living import run as m


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


def missing_inputs():
    missing = [p for p in m.IGNORED if not (m.ROOT / p).exists()]
    for key in m.BOUND:
        try:
            m.bound_bytes(key)
        except m.MissingInput as error:
            missing.append(str(error))
    return missing


@pytest.fixture(scope='module')
def fresh():
    missing = missing_inputs()
    if missing:
        pytest.skip(f'inputs absent: {missing}')
    return json.loads(m.dumps(m.results()))


def test_product_regenerates(product, fresh):
    assert product['schema'] == m.SCHEMA
    assert product == fresh


def test_producer_and_input_hashes(product):
    producer = product['producer']
    for path, h in {**producer['files'], **producer['inputs']}.items():
        assert m.digest(m.ROOT / path) == h, path
    for path, h in producer['ignored_inputs'].items():
        if (m.ROOT / path).exists():
            assert m.digest(m.ROOT / path) == h, path
    assert not constants_changed(producer['constants'])
    assert set(producer['ignored_inputs']) == set(m.IGNORED)


def test_bound_inputs_pinned(product):
    bound = product['producer']['bound_inputs']
    assert set(bound) == set(m.BOUND)
    for key, b in bound.items():
        pinned = m.BOUND[key]
        assert (b['path_at_integration'], b['branch'], b['commit'], b['sha256']) == \
            (pinned['path'], pinned['branch'], pinned['commit'], pinned['sha256'])
        try:
            raw = m.bound_bytes(key)
        except m.MissingInput:
            pytest.skip(f'{key} is not readable here')
        assert json.loads(raw)['schema'] == pinned['schema']


def test_rejects_wrong_schema(tmp_path, monkeypatch):
    bad = tmp_path / 'wrong.json'
    bad.write_text('{"schema": "wrong"}')
    monkeypatch.setattr(m, 'ROOT', tmp_path)
    with pytest.raises(ValueError):
        m.load_json('wrong.json', 'terluna.scenario.population/1')
    with pytest.raises(m.MissingInput):
        m.load_npz('absent.npz', 'terluna.geography.atlas-grid/1')


def test_people_conserved_in_every_mix(product):
    shared = json.loads((m.ROOT / 'shared/scenarios/population.json').read_text())
    for name, case in product['cases'].items():
        c = shared['cases'][name]
        lunar = c['lunar_surface'] + c['lunar_aerial']
        assert math.isclose(case['lunar_billion'], lunar)
        off = case['off_moon']
        assert math.isclose(off['earth']['people_billion'] + lunar + off['orbital']['people_billion'] +
                            off['elsewhere']['people_billion'], c['total'])
        for f, e in case['by_aerial_share'].items():
            assert math.isclose(e['lunar_surface_billion'] + e['lunar_aerial_billion'], lunar, rel_tol=1e-6)
            assert e['is_case_split'] == math.isclose(e['lunar_aerial_billion'], c['lunar_aerial'], rel_tol=1e-6)
            for mix, s in e['surface'].items():
                assert math.isclose(sum(s['people_billion'].values()), e['lunar_surface_billion'], rel_tol=1e-5)
                assert all(v >= 0 for v in s['people_billion'].values())
            for amix, a in e['aerial'].items():
                assert math.isclose(sum(a['people_billion'].values()), e['lunar_aerial_billion'], rel_tol=1e-5,
                                    abs_tol=1e-9)


def test_mix_rules():
    homes = {k: 1.0 for k in m.SURFACE_MIXED}
    fixed = {'summit_metropolis': 0.1, 'undersea_communities': 0.0}
    for mix, rule in m.MIXES.items():
        people = m.place_surface(mix, 5.0, homes, fixed)
        assert math.isclose(sum(people.values()), 5.0)
        favoured = [people[k] for k in rule['favoured']]
        others = [people[k] for k in m.SURFACE_MIXED if k not in rule['favoured']]
        if favoured and others:
            assert math.isclose(min(favoured) / max(others), m.MIX_WEIGHT)
    with pytest.raises(ValueError):
        m.place_surface('spread', 0.05, homes, fixed)
    for amix in m.AERIAL_MIXES:
        assert math.isclose(sum(m.place_aerial(amix, 3.0).values()), 3.0)


def test_twilight_reproduces_the_register(product):
    d = product['definitions']
    assert abs(d['practical_dusk_lux'] - 2.98) < 0.01
    assert abs(d['practical_dusk_sun_deg'] + 41.4) < 0.1
    assert abs(d['no_dark_night_poleward_of_deg'] - 48.6) < 0.1
    rows = {row['lat_deg']: row for row in product['light']['by_latitude']}
    assert abs(rows[0]['practical_dusk_ends_after_sunset_hours'] - 82.0) < 1.0
    assert abs(rows[30]['practical_dusk_ends_after_sunset_hours'] - 98.0) < 1.0
    assert abs(rows[0]['day_vision_ends_after_sunset_hours'] - 46.0) < 1.0
    assert rows[50]['practical_dusk_ends_after_sunset_hours'] is None
    assert rows[50]['beyond_practical_dusk_hours'] == 0.0
    month = SYNODIC_MONTH_DAYS * 24.0
    assert math.isclose(rows[0]['beyond_practical_dusk_hours'] + 2 * rows[0]['practical_dusk_ends_after_sunset_hours'],
                        month / 2, rel_tol=0.01)


def test_hours_to_depression_limits():
    month = SYNODIC_MONTH_DAYS * 24.0
    assert math.isclose(m.hours_to_depression(0.0, 30.0), 30.0 / 360.0 * month)
    assert m.hours_to_depression(60.0, 30.0) is not None
    assert m.hours_to_depression(60.0, 30.1) is None
    assert m.hours_to_depression(90.0, 1.0) is None
    hs = [m.hours_to_depression(20.0, x) for x in (5.0, 20.0, 40.0, 60.0)]
    assert all(a < b for a, b in zip(hs, hs[1:]))
    dip = math.degrees(math.acos(MOON_RADIUS / (MOON_RADIUS + 1e4)))
    assert math.isclose(m.sun_up_longer_aloft_hours(10.0), dip / 360.0 * month)


def test_units_per_billion(product):
    common = product['per_billion_common']
    pop = json.loads((m.ROOT / 'provisioning/results/population.json').read_text())
    placement = pop['heat']['placement_factor']['industry_in_orbit_thermal']
    per_tw = pop['heat']['dimming_percent_per_tw']
    for k, kw in common['energy_tw'].items():
        heat = common['heat_on_moon_tw'][k]
        assert math.isclose(heat, kw * placement, rel_tol=2e-3)
        assert math.isclose(common['added_dimming_percent'][k][0], heat * per_tw[0], rel_tol=5e-3)
    for sid, s in product['settings'].items():
        if 'settlement_m2_per_person' in s:
            for km2, m2 in zip(s['per_billion']['settlement_km2'], s['settlement_m2_per_person']):
                assert math.isclose(km2, 1e3 * m2, rel_tol=5e-3), sid
    assert math.isclose(common['lit_sky_km2'][0], 1e3 * 791.0, rel_tol=2e-3)
    assert all(a <= b for a, b in zip(common['domestic_water_makeup_km3_year'], common['domestic_water_gross_km3_year']))


def test_land_partition(product):
    moon = product['moon']
    regions = moon['regions']
    total = sum(v['dry_land_km2'] for v in regions.values() if v)
    assert math.isclose(total, moon['dry_land_km2'], rel_tol=2e-3)
    for v in regions.values():
        if v:
            assert math.isclose(sum(v['earth_up_shares'].values()), 1.0, abs_tol=5e-3)
            assert all(b >= 0 for b in v['bands_dry_land_km2'].values())
    homes = {}
    for sid, s in product['settings'].items():
        if 'home_dry_land_km2' in s:
            homes.setdefault(tuple(s['home_classes']), 0.0)
            homes[tuple(s['home_classes'])] += s['home_dry_land_km2']
    for classes, km2 in homes.items():
        assert math.isclose(km2, sum(regions[c]['dry_land_km2'] for c in classes), rel_tol=2e-3)


def test_food_fill_conserves():
    need = 3.0
    avail = {'a': 4.0, 'b': 10.0}
    yields = {'a': 0.5, 'b': 0.1}
    used, left = m.food_fill(need, avail, yields)
    assert math.isclose(used['a'] * 0.5 + used['b'] * 0.1 + left, need)
    assert math.isclose(used['a'], 4.0)                    # the wetter band fills first
    used, left = m.food_fill(0.0, avail, yields)
    assert left == 0.0 and all(v == 0.0 for v in used.values())
    used, left = m.food_fill(100.0, avail, yields)
    assert math.isclose(left, 100.0 - 4.0 * 0.5 - 10.0 * 0.1)


def test_food_shares_against_rain_fed_land(product):
    eco_rain_fed = product['moon']['rain_fed_land_ecology_km2']
    for case in product['cases'].values():
        e = case['by_aerial_share']['0.5']
        need = e['food']['kept_30']['cropland_equivalent_km2']
        share = e['food']['kept_30']['share_of_rain_fed_land_outside_kept']
        assert math.isclose(share[0], need[0] / (eco_rain_fed[1] * 0.7), rel_tol=5e-3)
        assert math.isclose(share[1], need[1] / (eco_rain_fed[0] * 0.7), rel_tol=5e-3)
        for s in e['surface'].values():
            for tag, use in s['rain_belt_land_use'].items():
                assert 0.0 <= use[0] <= use[1] <= 1.0 + 1e-9


def test_drag_and_holding(product):
    a = product['aerial_reference']
    assert a['hull']['cd_over_efficiency_on_volume_area'] < a['sphere_cd_over_efficiency_on_volume_area']
    assert math.isclose(a['hull']['drag_ratio_sphere_to_hull_same_volume'],
                        a['sphere_cd_over_efficiency_on_volume_area'] / a['hull']['cd_over_efficiency_on_volume_area'],
                        rel_tol=5e-3)
    follow = a['hull_follow_sun_kw_per_resident_by_lat']
    assert math.isclose(follow['60'] / follow['0'], math.cos(math.radians(60.0)) ** 3, rel_tol=1e-2)
    assert a['hull_hold_place_kw_per_resident'] < a['sphere_holding_kw_per_resident_5m_s']
    assert math.isclose(a['wind_effective_m_s'] ** 3, a['wind_mean_cube_m3_s3'], rel_tol=1e-2)
    duty = a['roaming_kw_per_resident']
    assert math.isclose(duty['0.1'] / duty['0.01'], 10.0, rel_tol=1e-2)


def test_mean_cube_speed_of_steady_wind():
    winds = {'height_km': [0.0, 10.0]}
    for key in ('speed_p10_m_s', 'speed_p50_m_s', 'speed_p90_m_s', 'speed_p99_m_s', 'speed_p99_9_m_s',
                'speed_max_m_s'):
        winds[key] = [4.0, 4.0]
    v3 = m.mean_cube_speed(winds, 5.0)
    # the quantile curve rises from zero to 4 m/s over the first tenth and then stays: 0.1 * 4^3 / 4 + 0.9 * 64
    assert math.isclose(v3, 0.1 * 64.0 / 4.0 + 0.9 * 64.0, rel_tol=1e-3)


def test_sky_spacing_and_scaling(product):
    kept = product['assumptions']['kept_shares']
    clear = {b: v['sky_clear_of_terrain_km2'] for b, v in product['moon']['bands'].items()}
    for case in product['cases'].values():
        rows = case['by_aerial_share']
        base = rows['0.25']['sky']
        for f, e in rows.items():
            sky = e['sky']
            if e['lunar_aerial_billion'] == 0:
                assert sky['envelopes'] == 0 and sky['by_band']['equatorial']['spacing_km_kept_30'] is None
                continue
            ratio = e['lunar_aerial_billion'] / rows['0.25']['lunar_aerial_billion']
            assert math.isclose(sky['mass_gt'], base['mass_gt'] * ratio, rel_tol=1e-2)
            assert math.isclose(sky['lifting_gas_mixture_gt'], base['lifting_gas_mixture_gt'] * ratio, rel_tol=1e-2)
            for b, row in sky['by_band'].items():
                for k in kept:
                    sp = row[f'spacing_km_kept_{int(k * 100)}']
                    assert math.isclose(sp, math.sqrt(clear[b] * (1 - k) / row['towns']), rel_tol=1e-2)


def test_hydrogen_upkeep(product):
    a = product['aerial_reference']['hydrogen_upkeep']
    per = a['per_billion_per_percent_a_year']
    assert math.isclose(per['air_hydrogen_ppm_held'],
                        per['released_mt_year'] * a['soil_uptake_years'] * a['air_hydrogen_ppm_per_mt_held'],
                        rel_tol=1e-2)
    e = product['cases']['total25']['by_aerial_share']['0.5']
    up = e['sky']['hydrogen_upkeep']
    aloft = e['lunar_aerial_billion']
    assert math.isclose(up['released_mt_year'][0], aloft * per['released_mt_year'] * 100 * a['loss_bracket_per_year'][0],
                        rel_tol=1e-2)


def test_grid_helpers():
    mask = np.zeros((4, 6), bool)
    mask[2, 0] = True
    n = m.neighbours(mask)
    assert n[2, 5] and n[1, 1] and n[3, 0] and not n[2, 0] and not n[0, 3]
    vals = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.array([1.0, 1.0, 1.0, 1.0])
    assert m.weighted_pct(vals, w, (0.5,)) == [2.0]
    area = m.pixel_areas(np.array([0.125]), 4, 1)
    assert math.isclose(float(area[0, 0]), (MOON_RADIUS / 1e3) ** 2 * math.radians(0.25) *
                        (math.sin(math.radians(0.25)) - 0.0))


def test_winds_source(product):
    winds = product['winds']
    if winds['source'] == m.ZONAL_WINDS:
        assert m.ZONAL_WINDS in product['producer']['inputs']
        assert math.isclose(sum(winds['occupancy_by_band'].values()), 1.0, abs_tol=2e-2)
        for row in winds['hold_time_of_day']:
            assert row['airspeed_m_s_best_hour'] <= row['airspeed_m_s_hour_mean'] + 1e-9
    else:
        assert 'placeholder' in winds['status']


def test_off_moon_figures(product):
    for case in product['cases'].values():
        earth = case['off_moon']['earth']
        share = earth['share_of_food_grown_off_the_land']
        assert 0.0 <= share < 1.0
        if earth['over_food_within_boundaries'] <= 1.0:
            assert share == 0.0
        orbital = case['off_moon']['orbital']
        assert orbital['shield_gt'][0] < orbital['shield_gt'][1]
