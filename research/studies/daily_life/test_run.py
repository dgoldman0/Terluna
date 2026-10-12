"""Physical limits, the light states' bookkeeping and the bindings of the daily-life product."""
import json
import math

import numpy as np
import pytest

from shared.provenance import constants_changed
from research.studies.daily_life import run as m


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


def test_stored_product_and_inputs(product):
    if not (m.ROOT / m.CLIMATOLOGY).exists():
        pytest.skip('The design run\'s stored climatology (ignored input) is absent; the stored product was made with it')
    if m.main_places()[0] is None:
        pytest.skip('Main\'s pinned places product is not readable from Git here')
    assert product['schema'] == m.SCHEMA
    assert product == json.loads(json.dumps(m.results()))
    for path, h in {**product['producer']['files'], **product['producer']['inputs']}.items():
        assert m.digest(m.ROOT / path) == h, path
    assert product['producer']['optional_inputs'][m.CLIMATOLOGY] == m.digest(m.ROOT / m.CLIMATOLOGY)
    assert not constants_changed(product['producer']['constants'])


def test_rejects_wrong_schema(tmp_path):
    p = tmp_path / 'wrong.json'
    p.write_text('{"schema": "wrong"}')
    with pytest.raises(ValueError):
        m.load_product(p, m.SCHEMA)


def test_main_places_bound_by_hash():
    data, source = m.main_places()
    if data is None:
        pytest.skip(f'Main\'s places product unavailable ({source})')
    assert data['schema'] == m.MAIN_PLACES['schema']
    assert source in ('in-tree', 'git object at the pinned commit')
    raw = json.dumps(data)
    assert 'natural_darkest_lux' in raw          # the corrected product carries the darkest hour; the old one did not


def test_states_partition_the_cycle(product):
    for name, c in product['light']['calendars'].items():
        assert math.isclose(c['day_h'] + c['dusk_h'] + c['night_h'] + c['dawn_h'], m.CYCLE_H, abs_tol=0.05), name
        assert math.isclose(c['below_dusk_h'] + c['earth_holds_above_dusk_h'], c['night_h'], abs_tol=0.05), name
        assert min(c['day_h'], c['dusk_h'], c['night_h'], c['dawn_h'], c['below_dusk_h']) >= 0, name
        assert c['earth_holds_above_dusk_h'] >= -0.05, name
        assert c['below_dark_h'] <= c['below_dusk_h'] + 0.05, name
        assert math.isclose(c['dusk_h'], c['dawn_h'], abs_tol=0.05), name            # the Sun's twilight is symmetric
        if c['dusk_runs_into_dawn']:
            assert c['night_h'] == 0 and c['darkest_lux'] >= product['levels']['dusk'] * 0.999, name
        assert c['darkest_lux'] <= c['midnight_lux'] + 1e-9, name


def test_register_light_and_time(product):
    cal, cycle = product['light']['calendars'], product['cycle']
    assert 81 < cal['far side 0']['dusk_h'] < 83                      # 82 h after an equatorial sunset
    assert 97 < cal['far side 30']['dusk_h'] < 99.5                   # 98 h at 30 degrees
    assert abs(cycle['never_falls_to_practical_dusk_poleward_of_deg'] - 48.6) < 0.3
    assert 45.5 < cal['far side 0']['day_vision_after_sunset_h'] < 47.5   # day vision for 46 h
    assert abs(product['levels']['dusk'] - 2.98) < 0.01
    assert abs(cal['far side 0']['below_dusk_h'] - 190.2) < 1.5       # the joint synthesis's far-side equator
    assert cal['nearside 0']['dusk_runs_into_dawn'] is False and cal['nearside 0']['below_dusk_h'] == 0
    assert cal['polar_plain']['dusk_runs_into_dawn'] and cal['far side 60']['dusk_runs_into_dawn']
    for lat in ('0', '15', '30', '45'):                                # the far side darkens; the high Earth holds the nearside
        assert cal[f'far side {lat}']['darkest_lux'] < cal[f'nearside {lat}']['darkest_lux']


def test_study_light_matches_main_corrected_earthlight(product):
    check = product['light']['cross_check_main_places']
    if not check['places']:
        pytest.skip('Main\'s places product was not readable when the product was made')
    assert len(check['places']) == 7
    for name, c in check['places'].items():
        for key in ('midnight_lux', 'darkest_lux'):
            ours, theirs = c[f'study_{key}'], c[f'main_{key}']
            assert abs(ours - theirs) <= max(0.003 * theirs, 2e-4), (name, key, ours, theirs)
        assert abs(c['study_earth_elevation_deg'] - c['main_earth_elevation_deg']) <= 0.1, name
        assert abs(c['study_hours_below_dark'] - c['main_hours_below_dark']) <= 1.0, name


def test_crossing_hours_on_an_exponential():
    t = np.linspace(0, 100, 2001)[:-1]
    lux = 1000 * np.exp(-0.1 * t)                    # falls to 1 lux at t = 10 ln(1000)
    below, falls, rises = m.crossing_hours(t, lux, 1.0)
    exact = 10 * math.log(1000)
    assert math.isclose(falls[0], exact, abs_tol=1e-6)
    assert math.isclose(below, (t[-1] + t[1] - t[0]) - exact, abs_tol=0.06)   # the wrap back to t=0 rises at the end
    assert len(rises) == 1


def test_twilight_curve_agrees_with_the_calendars(product):
    curve = {c['lat']: c['reaches']['practical dusk']['hours_after_sunset'] for c in product['light']['twilight']}
    cal = product['light']['calendars']
    for lat in (0.0, 15.0, 30.0, 45.0):
        assert abs(curve[lat] - cal[f'far side {lat:g}']['dusk_h']) < 0.2
    assert curve[60.0] is None and curve[75.0] is None
    levels = [r['lux'] for r in product['light']['twilight'][0]['light']]
    assert all(a > b for a, b in zip(levels, levels[1:]))                 # the equator's light only falls after sunset


def test_town_day_limits():
    h, lat = 10.0, 0.0
    v = m.sun_ground_speed(h, lat)
    assert math.isclose(v, 2 * math.pi * (m.MOON_RADIUS + 1e4) / (m.CYCLE_H * 3600))
    assert 4.25 < v < 4.35
    assert math.isclose(m.town_day_hours(h, lat, 0.0)[0], m.CYCLE_H)
    assert m.town_day_hours(h, lat, -v) == (None, False)                  # drifting west at the Sun's speed holds the hour
    day, reversed_ = m.town_day_hours(h, lat, -2 * v)
    assert math.isclose(day, m.CYCLE_H) and reversed_                      # the Sun runs backwards at the same pace
    assert math.isclose(m.town_day_hours(h, lat, -v / 2)[0], 2 * m.CYCLE_H)
    assert m.town_day_hours(h, lat, 1e6)[0] < 1.0
    days = [m.town_day_hours(h, lat, u)[0] for u in (-v * 0.9, -2.0, 0.0, 2.0, 8.0)]
    assert all(a > b for a, b in zip(days, days[1:]))                     # eastward drift shortens the day
    assert math.isclose(m.sun_ground_speed(0, 0) * 3.6, 15.4, abs_tol=0.05)
    assert m.sun_ground_speed(0, 60) < m.sun_ground_speed(0, 0)


def test_horizon_and_sunlit_hours():
    assert m.horizon_dip_deg(0) == 0
    assert math.isclose(m.sunlit_hours_at_height(0, 0), m.CYCLE_H / 2)
    dips = [m.horizon_dip_deg(h) for h in (0.1, 1, 10, 24)]
    assert all(a < b for a, b in zip(dips, dips[1:]))
    # Independent: at the equator the Sun sets once it is the dip below the level horizon.
    extra = 2 * m.horizon_dip_deg(10) / m.DEG_PER_H
    assert math.isclose(m.sunlit_hours_at_height(10, 0) - m.CYCLE_H / 2, extra, rel_tol=1e-9)


def test_seasonal_day_and_cycle(product):
    assert product['cycle']['seasonal_day']['0']['longest_day_h'] == product['cycle']['seasonal_day']['0']['shortest_day_h']
    for lat in ('30', '60', '75', '85'):
        s = product['cycle']['seasonal_day'][lat]
        assert math.isclose(s['longest_day_h'] + s['shortest_day_h'], m.CYCLE_H, abs_tol=0.11)
    swings = [product['cycle']['seasonal_day'][k]['longest_day_h'] for k in ('0', '30', '60', '75', '85')]
    assert all(a < b for a, b in zip(swings, swings[1:]))
    assert m.seasonal_day(89.0) is None


def test_oxygen_equivalent_heights(product):
    sea = m.COESA['o2_share'] * m.COESA['p0_pa']
    assert abs(m.earth_equivalent_m(sea)) < 1e-6
    exponent = m.STANDARD_GRAVITY * m.COESA['molar_mass_kg_mol'] / (m.GAS_CONSTANT * m.COESA['lapse_k_m'])
    p2000 = m.COESA['p0_pa'] * (1 - m.COESA['lapse_k_m'] * 2000 / m.COESA['t0_k']) ** exponent
    assert math.isclose(m.earth_equivalent_m(m.COESA['o2_share'] * p2000), 2000, abs_tol=1e-6)
    zones = product['volume']['open_air']
    heights = [z['earth_equivalent_m'] for z in zones]
    assert all(a < b for a, b in zip(heights, heights[1:]))
    # The summit brief's table: Mexico City at the summit, about La Paz at 10 km above it, El Alto at 15 km.
    assert 1700 < zones[0]['earth_equivalent_m'] < 2000
    assert 3300 < zones[2]['earth_equivalent_m'] < 3600
    assert 4000 < zones[4]['earth_equivalent_m'] < 4400
    assert abs(product['volume']['o2_share'] - 0.1746) < 0.0005


def test_kepler_transfers():
    mu, r_moon, r_earth = m.EARTH_GM, m.EARTH_MOON_DISTANCE, m.EARTH_RADIUS
    hohmann = m.kepler_hours(r_moon, r_earth, r_moon, mu)
    a = (r_moon + r_earth) / 2
    assert math.isclose(hohmann, math.pi * math.sqrt(a**3 / mu) / 3600, rel_tol=1e-9)
    faster = m.kepler_hours(r_moon, r_earth, 2 * r_moon, mu)
    parabola = m.kepler_hours(r_moon, r_earth, math.inf, mu)
    assert parabola < faster < hohmann
    # A very long ellipse approaches the parabola.
    assert math.isclose(m.kepler_hours(r_moon, r_earth, 1e6 * r_moon, mu), parabola, rel_tol=1e-3)


def test_falls_and_spin():
    small = m.fall(0.5, 17.9)
    assert math.isclose(small['impact_m_s'], math.sqrt(2 * m.MOON_SURFACE_GRAVITY * 0.5), rel_tol=0.01)
    assert math.isclose(small['like_earth_fall_m'], 0.5 * m.MOON_SURFACE_GRAVITY / m.STANDARD_GRAVITY, abs_tol=0.01)
    assert m.fall(1e4, 17.9)['impact_m_s'] <= 17.9
    for rpm in (1, 2, 4):
        omega = rpm * 2 * math.pi / 60
        r = m.centrifuge_radius_m(rpm, m.STANDARD_GRAVITY, m.MOON_SURFACE_GRAVITY)
        assert math.isclose(math.hypot(omega**2 * r, m.MOON_SURFACE_GRAVITY), m.STANDARD_GRAVITY)
        assert math.isclose(m.centrifuge_radius_m(rpm, m.STANDARD_GRAVITY, 0.0), m.STANDARD_GRAVITY / omega**2)


def test_earth_ties(product):
    e = product['earth']
    assert e['light_one_way_s']['perigee'] < e['light_one_way_s']['mean'] < e['light_one_way_s']['apogee']
    assert 24.7 < e['earth_turns_once_h'] < 25.0                       # the Earth's face turns once a tidal day
    assert all(e['earth_elevation_by_place'][k] < -8 for k in ('summit', 'ingenii_coast', 'highland_box'))
    assert e['earth_elevation_by_place']['nubium'] > 50
    t = e['transfer_days']
    assert 4.8 < t['apogee at the Moon (least energy)'] < 5.1


def test_travel_phase_offsets_are_antisymmetric(product):
    places = {k: v for s in product['settings'].values() for k, v in s['places'].items()}
    for row in product['travel']['pairs']:
        dlon = (places[row['b']]['lon'] - places[row['a']]['lon'] + 180) % 360 - 180
        assert math.isclose(row['phase_offset_h'], dlon / m.DEG_PER_H, abs_tol=0.06)
        assert abs(row['phase_offset_h']) <= m.CYCLE_H / 2 + 0.06
        assert row['winged_liner_h'] < row['long_haul_h']


@pytest.fixture(scope='module')
def zonal():
    return m.load_product(m.ROOT / m.ZONAL_WINDS, 'terluna.climate.gcm-zonal-winds/1')


def test_zonal_winds_sun_speed_is_this_studys(zonal):
    winds = m.Winds(zonal)
    for h in (0.0, 10.0, 40.0, 90.0):
        i = winds.row(h)
        for lat, v in zip(winds.lat, winds.arrays['v_sun'][i]):
            assert math.isclose(v, m.sun_ground_speed(h, lat), abs_tol=0.002), (h, lat)


def test_zonal_winds_read_as_the_note_reads_them(product, zonal):
    winds = m.Winds(zonal)
    # climate/gcm/zonal_winds.md: the median drifting day on the equator and at 60 degrees, and the mean eastward wind.
    for band, expected in (('equator', (21.8, 17.1, 9.3)), ('60', (17.3, 10.4, 4.3))):
        for h, day in zip((10.0, 20.0, 40.0), expected):
            assert abs(winds.band('day_p50', h, m.WIND_BANDS[band]) - day) <= 0.15, (band, h)
    for h, u in ((7.5, 0.85), (20.0, 3.5), (40.0, 10.7), (75.0, 21.9)):
        assert abs(winds.global_u(h) - u) <= 0.06, h
    profile = [r['u_m_s'] for r in product['volume']['wind_by_height']]
    assert all(a < b for a, b in zip(profile[1:], profile[2:]))          # superrotation grows with height above 2.5 km


def test_sky_town_days(product):
    towns = product['volume']['sky_towns']
    for t in towns:
        d = t['day_days']
        assert 0 < d['p10'] <= d['p50'] <= d['p90']
        assert 0 <= t['held_within_1_m_s'] <= 1 and 0 <= t['sun_turned_back'] <= 1
        v, u = t['sun_speed_m_s'], t['u_mean_m_s']
        assert math.isclose(t['day_at_mean_wind_days'], m.CYCLE_H / 24 * v / abs(v + u), abs_tol=0.06)
    for band in m.WIND_BANDS:
        days = [t['day_days']['p50'] for t in towns if t['latitudes'] == band and t['height_km'] >= 10]
        assert all(a > b for a, b in zip(days, days[1:])), band       # above 8 km the higher town lives the shorter day
        assert all(day < m.CYCLE_H / 24 for day in days)
    holds = {tuple(r['abs_lat_deg']): r for r in product['volume']['wind_holds_the_hour'] if r['height_km'] == [0.0, 5.0]}
    assert holds[(60.0, 75.0)]['held_within_1_m_s'] > holds[(40.0, 60.0)]['held_within_1_m_s'] > holds[(0.0, 20.0)]['held_within_1_m_s']


def test_routes_split_their_day_and_linger_from_evening_to_morning(product):
    routes = {r['height_km']: r for r in product['volume']['routes']}
    for r in routes.values():
        assert math.isclose(sum(r['days_by_quarter'].values()), r['mean_day_days'], abs_tol=0.25)
        assert r['years_5_10_poleward_of_70_deg'] <= r['years_5_10_poleward_of_50_deg'] + 1e-9
    assert routes[2.5]['lingers_at_hour_deg'] > 90                    # the evening near the ground
    assert abs(routes[10.0]['lingers_at_hour_deg']) < 30               # about noon at 10 km
    assert routes[15.0]['lingers_at_hour_deg'] < 0 and routes[20.0]['lingers_at_hour_deg'] < 0   # the morning higher up
    assert routes[2.5]['days_by_quarter']['evening_and_first_half_of_night'] > m.CYCLE_H / 24 / 4
    assert routes[30.0]['linger_strength'] < 1.1                        # an even pace at 30 km
