"""Limits, bindings and consistency of the land-cities product."""
import json
import math

import numpy as np
import pytest

from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_changed
from research.studies.land_cities import grid, groundwater as gw, highlands as hl, run as m, tubes as tb


@pytest.fixture(scope='module')
def product():
    return json.loads(m.OUT.read_text())


@pytest.fixture(scope='module')
def fresh():
    missing = m.missing_inputs()
    if missing:
        pytest.skip(f'inputs absent: {missing}')
    return json.loads(m.dumps(m.product()))


# ---- bindings ----

def test_product_regenerates(product, fresh):
    assert product['schema'] == m.SCHEMA == 'terluna.research.land-cities/1'
    assert product == fresh


def test_producer_and_input_hashes(product):
    producer = product['producer']
    assert set(producer['files']) == set(m.FILES)
    assert set(producer['inputs']) == set(m.INPUTS)
    for path, h in {**producer['files'], **producer['inputs']}.items():
        assert m.digest(m.ROOT / path) == h, path
    for path, h in producer['ignored_inputs'].items():
        if (m.ROOT / path).exists():
            assert m.digest(m.ROOT / path) == h, path
    assert set(producer['ignored_inputs']) <= set(m.IGNORED)
    assert not constants_changed(producer['constants'])
    assert 'MOON_SURFACE_GRAVITY' in producer['constants']


def test_bound_input_pinned(product):
    bound = product['producer']['bound_inputs']
    assert set(bound) == set(m.BOUND)
    for key, b in bound.items():
        pinned = m.BOUND[key]
        assert (b['path_at_integration'], b['branch'], b['commit'], b['sha256']) == \
            (pinned['path'], pinned['branch'], pinned['commit'], pinned['sha256'])
        try:
            data = m.load_bound(key)
        except m.MissingInput:
            pytest.skip(f'{key} is not readable here')
        assert data['schema'] == pinned['schema']
        assert 'natural_darkest_lux' in json.dumps(data)      # main's corrected product carries the darkest hour


def test_rejects_wrong_schema(tmp_path, monkeypatch):
    bad = tmp_path / 'wrong.json'
    bad.write_text('{"schema": "wrong"}')
    monkeypatch.setattr(m, 'ROOT', tmp_path)
    with pytest.raises(ValueError):
        m.load_json('wrong.json', 'terluna.geography.atlas/1')
    with pytest.raises(m.MissingInput):
        m.load_npz('absent.npz', 'terluna.geography.atlas-grid/1')


# ---- physics at lunar gravity ----

def test_pressure_and_conductivity_scale_with_gravity():
    assert math.isclose(float(tb.pressure_pa(100.0)), 1000.0 * MOON_SURFACE_GRAVITY * 100.0)
    ratio = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    assert math.isclose(float(tb.earth_head_m(600.0)), 600.0 * ratio)
    assert math.isclose(float(tb.k_moon(1e-4)), 1e-4 * ratio)
    assert math.isclose(tb.ITA_MEMBRANE_HEAD_EARTH_M / ratio, 15.0 * STANDARD_GRAVITY / MOON_SURFACE_GRAVITY)


def test_mound_limits():
    R, K, b, L = 1e-8, 1e-5, 500.0, 50e3
    assert tb.mound_m(0.0, L, R, K, b) == 0.0
    x = np.linspace(0.0, L, 11)
    h = tb.mound_m(x, L, R, K, b)
    assert np.all(np.diff(h) > 0)                          # rises from the shore to the divide
    assert float(tb.mound_m(2 * L, L, R, K, b)) == float(h[-1])    # flat beyond the divide
    thick = tb.mound_m(L, L, R, K, 1e6)                    # a thick aquifer: the linear limit R L^2 / (2 K b)
    assert math.isclose(float(thick), R * L ** 2 / (2 * K * 1e6), rel_tol=1e-3)
    thin = tb.mound_m(L, L, R, K, 0.0)                     # Woessner & Poeter's flat-base form with h at the shore 0
    assert math.isclose(float(thin), math.sqrt(R / K) * L, rel_tol=1e-9)
    limits = tb.mound_ratio_limits()
    ratio = STANDARD_GRAVITY / MOON_SURFACE_GRAVITY
    assert limits['mound_small_against_base'] == pytest.approx(ratio, rel=5e-3)
    assert limits['mound_large_against_base'] == pytest.approx(math.sqrt(ratio), rel=5e-3)


# ---- the two-dimensional water table ----

CELL_M = grid.RADIUS_KM * 1e3 * math.radians(0.25)


def strip(n, ground=1e9):
    """A one-row window on the equator: the sea in its first cell at level 0, land beyond with no flow at the end."""
    wet = np.zeros((1, n), bool)
    wet[0, 0] = True
    level = np.full((1, n), np.nan)
    level[0, 0] = 0.0
    g = np.full((1, n), float(ground))
    g[0, 0] = -10.0
    return gw.Window(wet, level, g, np.full((1, n), CELL_M ** 2), np.array([0.0]), (0, 0), 4)


def test_steady_strip_matches_the_dupuit_mound():
    R, K, b, n = 5e-9, 4e-5, 740.0, 40
    st = gw.steady(strip(n), np.full((1, n), R), K, b, cap=False)
    x = np.arange(n) * CELL_M
    L = (n - 0.5) * CELL_M                                 # the no-flow face beyond the last cell
    assert np.allclose(st['z'][0], tb.mound_m(x, L, R, K, b), atol=1e-6)


def test_capped_strip_seeps_beyond_its_free_boundary():
    R, K, b, n, G = 5e-9, 4e-5, 740.0, 60, 300.0
    st = gw.steady(strip(n, ground=G), np.full((1, n), R), K, b)
    z = st['z'][0]
    assert np.all(z[1:] <= G + 1e-6) and st['pinned'][0].any()
    # The free part is a mound with its crest at the ground: zero slope where the seepage face begins.
    xs = math.sqrt(2 * K * (0.5 * (b + G) ** 2 - 0.5 * b ** 2) / R)
    first = int(np.argmax(st['pinned'][0]))
    assert abs(first * CELL_M - xs) <= CELL_M
    x = np.arange(first) * CELL_M
    expect = np.sqrt(b * b + R / K * (2 * xs * x - x * x)) - b
    assert np.allclose(z[:first], expect, atol=2.0)


def nonlinear_strip(n, length, R, K, b, S, years, dt_years=1.0):
    """A fully implicit 1-D Boussinesq run (Picard on the transmissivity) on nodes 0..n, the shore held at node 0
    and no flow beyond node n; the water table above the shore after each year."""
    from scipy.linalg import solve_banded
    dx = length / n
    dt = dt_years * tb.SECONDS_PER_YEAR
    h = np.zeros(n + 1)
    out = []
    for _ in range(int(years / dt_years)):
        new = h.copy()
        for _ in range(4):
            T = K * (b + 0.5 * (new[1:] + new[:-1])) / dx ** 2     # faces between nodes
            ab = np.zeros((3, n + 1))
            ab[1, 0] = 1.0
            ab[1, 1:] = S / dt + T
            ab[1, 1:n] += T[1:]
            ab[0, 2:] = -T[1:]
            ab[2, :n] = -T
            rhs = S * h / dt + R
            rhs[0] = 0.0
            new = solve_banded((1, 1), ab, rhs)
        h = new
        out.append(h.copy())
    return np.array(out)


def test_rise_matches_a_nonlinear_strip_and_settles_on_the_steady_mound():
    R, K, b, S = 6.3e-9, 4.14e-5, 740.0, 0.07
    n, centre = 19, 6
    win = strip(n)
    win.centre = (0, centre)
    end = gw.steady(win, np.full((1, n), R), K, b, cap=False)
    series = gw.rise(win, np.full((1, n), R), K, S, 0.0, end)
    target = 0.9 * end['z'][0, centre]
    t2d = gw.first_reached(series, target)
    length = (n - 0.5) * CELL_M
    nodes = 370
    runs = nonlinear_strip(nodes, length, R, K, b, S, 1500.0)
    k = int(round(centre * CELL_M / (length / nodes)))
    t1d = float(np.argmax(runs[:, k] >= target) + 1)
    assert t2d == pytest.approx(t1d, rel=0.05)
    assert series['levels'][-1] == pytest.approx(end['z'][0, centre], abs=0.5)
    assert gw.first_reached(series, -1.0) == 0.0 and gw.first_reached(series, 1e6) is None


def test_tight_rock_rises_at_the_first_rate():
    R, K, b, S, n, G = 8e-9, 5e-8, 100.0, 0.07, 19, 600.0
    win = strip(n, ground=G)
    win.centre = (0, 12)
    end = gw.steady(win, np.full((1, n), R), K, b)
    series = gw.rise(win, np.full((1, n), R), K, S, 0.0, end)
    level = 400.0
    assert gw.first_reached(series, level) == pytest.approx(level / (R * tb.SECONDS_PER_YEAR / S), rel=0.02)


def test_land_water_tables(product):
    rows = [r_ for r_ in product['tubes']['on_the_atlas']['rows'] if 'water_table' in r_]
    assert len(rows) >= 10
    lo, hi = tb.mound_ratio_limits().values()
    for r_ in rows:
        w = r_['water_table']
        unc = [w[c]['uncapped_mound_above_sea_m'] for c in ('low', 'central', 'high')]
        assert unc == sorted(unc)                          # more rain into tighter, thinner rock: a higher mound
        for c in ('low', 'central', 'high'):
            assert w[c]['water_table_above_sea_m'] <= w[c]['uncapped_mound_above_sea_m'] + 0.1
            assert w[c]['water_table_above_sea_m'] <= r_['atlas_cell_above_sea_m'] + 0.1
            assert w[c]['depth_below_ground_m'] >= 0.0
        ratio = w['central']['rain_mound_moon_over_earth']
        assert ratio is None or hi - 0.01 <= ratio <= lo + 0.01
        rise = w['central'].get('rise')
        if rise and rise['years_to_90_percent']:
            a, b_ = rise['years_to_90_percent']
            assert a <= b_ and b_ / a == pytest.approx(tb.SPECIFIC_YIELD[2] / tb.SPECIFIC_YIELD[0], rel=1e-2)
        if r_.get('floor_below_water_table'):
            for c, verdict in r_['floor_below_water_table'].items():
                depth = w[c]['depth_below_ground_m']
                floor = r_['floor_below_ground_m']
                if min(abs(depth - f) for f in floor) <= 5e-3 * max(floor):
                    continue                                   # on a floor to the product's rounding
                assert verdict == ('yes' if depth < floor[0] else 'partly' if depth < floor[1] else 'no')


def test_vault_levels_and_the_estimate():
    widths = tb.level_widths(1000.0, 50.0, 5.0)
    assert widths.size == 9 and np.all(np.diff(widths) < 0) and widths.max() <= 1000.0 and widths.min() > 0.0
    grail = tb.holds_levels(1000.0, 10.0, 10, tb.HEIGHT_TO_WIDTH[0])
    assert grail['levels'] == [10, 10] and grail['usable_floor_km2'][1] < 1000.0 * 10.0 * 10 * 0.7 / 1e3
    est = tb.check_estimate()
    assert est['floor_km2'] == 500.0 and est['people_million'] == pytest.approx(7.14, rel=1e-3)
    vault = est['vault']['5_m_storeys']
    assert vault['levels'] == 9 and vault['full_width_floor_km2'] == 450.0 and vault['gross_floor_km2'] < 450.0
    circle = tb.envelope_m2_per_m(100.0, 50.0)            # half a circle of radius 50 over its floor
    assert math.isclose(circle, 100.0 + math.pi * 50.0, rel_tol=1e-9)


def test_seepage_linear_and_pumping():
    a = tb.seepage_m3_s(10.0, 1e6, 1e-8, 5.0)
    assert math.isclose(tb.seepage_m3_s(20.0, 1e6, 1e-8, 5.0), 2 * a)
    assert math.isclose(tb.seepage_m3_s(10.0, 1e6, 1e-7, 5.0), 10 * a)
    assert math.isclose(tb.pumping_mw(1.0, 100.0, 1.0), 1000.0 * MOON_SURFACE_GRAVITY * 100.0 / 1e6)


def test_first_phase_physics():
    hold = tb.roof_to_hold_air(101.325)
    lo, hi = hold['roof_m']
    assert lo < hi and math.isclose(lo, 101325.0 / (tb.BASALT_BULK_DENSITY[2] * MOON_SURFACE_GRAVITY), rel_tol=2e-3)
    col = tb.air_matching_earth_column()
    assert math.isclose(col['lunar_surface_pressure_kpa'], 101.325 * MOON_SURFACE_GRAVITY / STANDARD_GRAVITY,
                        rel_tol=2e-3)
    fast, slow = tb.warming_lag_years(100.0)
    assert fast < slow
    from scipy.special import erfc
    t = fast * tb.SECONDS_PER_YEAR
    assert erfc(100.0 / (2 * math.sqrt(tb.ROCK_DIFFUSIVITY[1] * t))) == pytest.approx(0.5, abs=0.01)


def test_air_and_floods_scale_with_gravity():
    moon = tb.stack_pressure_pa(100.0, 5.0, 1.2, 293.0)
    earth = tb.stack_pressure_pa(100.0, 5.0, 1.2, 293.0, gravity=STANDARD_GRAVITY)
    assert math.isclose(moon / earth, MOON_SURFACE_GRAVITY / STANDARD_GRAVITY)
    assert math.isclose(tb.breach_inflow_m3_s(100.0, 1.0, 1.0), math.sqrt(2 * MOON_SURFACE_GRAVITY * 100.0))


def test_beam_matches_the_haze_product(product):
    haze = m.load_json('illumination/aerosol/results/haze.json', m.INPUTS['illumination/aerosol/results/haze.json'])
    day = haze['regions']['seas']['cases']['clean']['day']['daylight']
    for z, v in day.items():
        assert tb.beam_transmission(float(z), day) == pytest.approx(v['clear']['direct'], rel=1e-6)
    rows = product['tubes']['lasting']['daylight']
    by_lat = {}
    for row in rows:
        by_lat.setdefault(row['target_lux'], []).append(row['mirror_m2_per_floor_m2'][0])
    for shares in by_lat.values():
        assert all(a <= b for a, b in zip(shares, shares[1:]))   # poleward the Sun stands lower: more mirror


# ---- the atlas and the inventory ----

def test_candidates_classified_once(product):
    rows = product['tubes']['on_the_atlas']['rows']
    keys = [r_['key'] for r_ in rows]
    assert len(keys) == len(set(keys))
    assert {r_['class'] for r_ in rows} <= {'under the sea', 'under a lake', 'coastal', 'dry'}
    counts = product['tubes']['on_the_atlas']['counts']
    assert sum(sum(v.values()) for v in counts.values()) == len(rows)
    for r_ in rows:
        if r_['class'] in ('under the sea', 'under a lake'):
            assert r_['water_depth_m'] > 0 and r_['ground_above_sea_m'] < 0 or r_['class'] == 'under a lake'
        else:
            assert r_['water_depth_m'] == 0.0 and r_['km_to_standing_water'] >= 0.0
            assert r_['ground_above_sea_m'] >= 0.0


def test_capacity_sums(product):
    cap = product['tubes']['capacity']
    for kind, t in cap['by_type'].items():
        members = [v for v in cap['by_candidate'].values() if v['type'] == kind]
        assert t['candidates'] == len(members)
        for size in ('small', 'large'):
            assert t[size][1] == pytest.approx(sum(v[size][1] for v in members), rel=1e-2)
    for label, tot in cap['totals'].items():
        for size in ('small', 'large'):
            assert tot['all'][size][1] == pytest.approx(sum(v[size][1] for v in tot['by_class'].values()), rel=1e-2)
        straddling = ' '.join(tot['by_class'].get('straddles the coast', {}).get('regions', []))
        assert 'marius_hills' in straddling
    assert cap['totals']['mapped']['all']['large'][1] < cap['totals']['with_gravity']['all']['large'][1]


def tiny_atlas():
    """A 1 px/deg atlas with the sea west of 10 E and a dry land elsewhere."""
    lat = 89.5 - np.arange(180.0)
    lon = 0.5 + np.arange(360.0)
    LA, LO = np.meshgrid(lat, lon, indexing='ij')
    area = grid.pixel_areas(lat, 1, 360)
    water = LO < 10.0
    zeros = np.zeros(LA.shape)
    return grid.Atlas(lat=lat, lon=lon, LA=LA, LO=LO, area=area, height_m=np.where(water, -500.0, 100.0),
                      sea_level_m=0.0, hsl_km=np.where(water, -0.5, 0.1), water=water,
                      water_label=water.astype(int), f_lake=zeros, lake_level_m=np.full(LA.shape, np.nan),
                      land=np.where(water, 0.0, area), dry=np.where(water, 0.0, area), side=zeros.astype(int), ppd=1)


def test_cell_takes_the_cell_south_and_east_of_an_edge():
    at = tiny_atlas()
    assert at.cell(10.0, 20.0) == (80, 20)                 # the cell spanning 9-10 N and 20-21 E
    assert at.cell(-89.9, 359.99) == (179, 359)
    assert at.cell(9.99, 20.01) == (80, 20)


def test_box_across_the_prime_meridian():
    at = tiny_atlas()
    box = tb.box_shares(at, ((-5.0, 5.0), (350.0, 20.0)))   # 350 E to 20 E: the sea from 0 to 10 E, land around it
    assert box['sea_share'] == pytest.approx(10.0 / 30.0, rel=1e-2)
    assert box['dry_share'] == pytest.approx(20.0 / 30.0, rel=1e-2)


def test_highland_fixes(product):
    hl_out = product['highlands']
    assert hl_out['scattered']['undersea_people'] == 48600.0
    boat = hl_out['access']['sky_boat']
    assert boat['range_with_climb_km'] == pytest.approx(
        boat['range_km'] - boat['climb_kwh'] / boat['cruise_kwh_per_km'], rel=1e-2)
    for trip in hl_out['access']['trips']:
        assert trip['sky_boat']['within_range_with_climb'] == (trip['km'] <= boat['range_with_climb_km'])
        if trip['sky_boat']['within_range_with_climb']:
            assert trip['sky_boat']['kwh_with_climb'] <= boat['mission_kwh'] * 1.01
    for row in hl_out['towers']['by_storeys']:
        c = row['canyon_at_coverage_0_2']
        assert 0.0 < c['sky_view_factor_at_least'] < 1.0 and 0.0 < c['open_share_of_street_side'] < 1.0


def test_inventory_bounds():
    for p in tb.PITS:
        assert -90 <= p['lat'] <= 90 and 0 <= p['lon'] < 360
        assert p['inner'][0] <= p['inner'][1] and p['floor_m'][0] <= p['floor_m'][1]
        assert p['floor_m'][1] >= p['depth']
    for x in tb.RILLES:
        for length, width, roof, rank in x['segments']:
            assert 0 < length < 30 and 0 < width < 2 and rank in 'ABC'
    for c in tb.CHAINS:
        assert 0 < c['intact_km'] <= c['chain_km']


# ---- the highlands ----

def test_highland_land_partitions(product):
    land = product['highlands']['land']
    total = land['all_highland']['dry_land_km2']
    assert sum(b['dry_land_km2'] for b in land['by_height']) == pytest.approx(total, rel=2e-3)
    assert sum(b['dry_land_km2'] for b in land['by_latitude']) == pytest.approx(total, rel=2e-3)
    halves = land['poleward_of_30']['dry_land_km2'] + land['equatorward_of_30']['dry_land_km2']
    assert halves == pytest.approx(total, rel=2e-3)
    for b in land['by_height']:
        assert sum(x['dry_land_km2'] for x in b['by_latitude']) == pytest.approx(b['dry_land_km2'], rel=2e-3)
    assert land['ways_of_living_home']['dry_land_km2'] < total < land['dry_land_km2']


def test_summit_has_the_most_flashes(product):
    land = product['highlands']['land']
    top = land['by_height'][-1]['flashes_per_km2_year']
    assert top == max(b['flashes_per_km2_year'] for b in land['by_height'])
    places = product['highlands']['lightning_at_places']
    assert places['summit port']['central'] > places['highland box']['central'] > places['polar plain']['central']


def test_towers_lighter_on_the_moon(product):
    t = product['highlands']['towers']
    for row in t['by_storeys']:
        if row['moon'] and row['earth']:
            assert row['moon']['frame_steel_kg_per_m2'] < row['earth']['frame_steel_kg_per_m2']
            assert row['moon']['sway_service_m'] <= hl.SWAY_LIMIT * row['height_m'] * 1.0001
        for c, d in row['people_per_km2_by_coverage'].items():
            assert d == pytest.approx(float(c) * row['storeys'] * 1e6 / hl.FLOOR_M2_PER_PERSON, rel=5e-3)
    assert t['wind']['wind_to_weight_against_earth'] > 1.0


def test_people_and_water_arithmetic(product):
    gp = product['highlands']['ground_people']
    for row in gp['rows']:
        expected = row['ground_billion'] * row['highland_share'] * 1e9 / row['people_per_km2']
        assert row['settled_km2'] == pytest.approx(expected, rel=5e-3)
    w = product['highlands']['water']
    per_km = w['gw_per_km3_year_per_km_of_lift']
    assert per_km == pytest.approx(1e9 * 1000.0 * MOON_SURFACE_GRAVITY * 1000.0 / w['pump_efficiency']
                                   / tb.SECONDS_PER_YEAR / 1e9, rel=5e-3)
    for k, v in w['per_billion_by_recovery'].items():
        assert v['from_sea_gw'][0] <= v['from_sea_gw'][1]


def test_grid_helpers():
    lat = np.array([0.125, -0.125])
    a = grid.pixel_areas(lat, 4, 1440)
    assert a.shape == (2, 1440) and np.allclose(a[0], a[1])
    assert grid.great_circle_km(0.0, 0.0, 0.0, 180.0) == pytest.approx(math.pi * grid.RADIUS_KM)
    field = np.arange(6.0).reshape(2, 3)
    out = grid.regrid(field, [10.0, -10.0], [0.0, 120.0, 240.0], np.array([10.0]), np.array([0.0, 120.0]))
    assert np.allclose(out, [[0.0, 1.0]])


def test_every_cited_source_is_listed():
    listed = {s_['key'] for s_ in json.loads((m.HERE / 'sources.json').read_text())['sources']}
    cited = set()

    def walk(x):
        if isinstance(x, dict):
            for k, v in x.items():
                if k in ('sources', 'source'):
                    cited.update([v] if isinstance(v, str) else v)
                else:
                    walk(v)
        elif isinstance(x, (list, tuple)):
            for v in x:
                walk(v)
    for name in ('PITS', 'GRAVITY', 'RADAR', 'TEMPERATURES', 'DOSE', 'QUAKES', 'NAMED_RILLE_NOTES'):
        walk(getattr(tb, name))
    cited.update(tb.RILLE_SOURCES + tb.CHAIN_SOURCES + tb.STABILITY_SOURCES + (hl.DESIGNED_LATERAL_SOURCE,))
    assert cited <= listed, sorted(cited - listed)


def test_missing_inputs_are_reported(monkeypatch):
    monkeypatch.setattr(m, 'IGNORED', {**m.IGNORED, 'geography/products/absent_grid.npz': None})
    assert 'geography/products/absent_grid.npz' in m.missing_inputs()
    monkeypatch.setattr(m, 'BOUND', {'places': dict(m.BOUND['places'], commit='0' * 40, sha256='0' * 64)})
    assert any('places.json' in x and 'not readable' in x for x in m.missing_inputs())
