"""Tests for the waters' warmth across the Moon: the product's binding and checks, the array form of the lake model's
bulk exchange against the lake model itself, the tables, windows and percentiles, the many-cell slab against the
forcing's one-site slab on the same synthetic climate, and the lake model's calm column beside the mixed layers."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from biosphere.ecology import lake as L
from research.studies.littoral_community import forcing as F
from research.studies.littoral_community import water_map as W
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import constants_changed

PRODUCT = W.OUT


@pytest.fixture(scope='module')
def product():
    if not PRODUCT.exists():
        pytest.skip('product not built')
    return json.loads(PRODUCT.read_text())


# ---- binding ---------------------------------------------------------------------------------------------------------
def test_product_schema_files_constants_and_inputs(product):
    assert product['schema'] == W.SCHEMA
    for path, digest in product['producer']['files'].items():
        assert F.sha256(W.ROOT / path) == digest, path
    assert constants_changed(product['producer']['constants']) == {}
    for path in (W.MOISTURE, W.FORCING):
        assert F.sha256(path) == product['inputs'][str(path.relative_to(W.ROOT))], path
    current = dict(zip(('air_gas_constant', 'air_cp', 'o2_fraction', 'salinity_g_kg'), F.air_constants()))
    assert product['design_air'] == pytest.approx(current, rel=1e-12)


def test_product_grid_cells_and_shares(product):
    g, c = product['grid'], product['cells']
    n = len(g['lat']) * len(g['lon'])
    assert n == 2048 and len(c['sea']) == n and len(g['lat_edges']) == len(g['lat']) + 1
    assert sum(c['area_km2']) == pytest.approx(4.0 * math.pi * 1737.4 ** 2, rel=1e-3)
    for name in [x['name'] for x in product['settings']['cases']] + ['sea_surface']:
        s = {k: np.array(v) for k, v in c[name].items()}
        assert np.all(s['p01_c'] <= s['mean_c']) and np.all(s['mean_c'] <= s['p99_c']), name
        assert np.all(s['p99_c'] <= s['max_c'] + 0.26), name
        assert np.all(s['hours_above_34'] >= s['hours_above_36']) and np.all(s['hours_above_36'] >= s['hours_above_38'])
        assert 10.0 < s['p01_c'].min() and s['max_c'].max() < 60.0, name
    for stat in product['shares'].values():
        for key, shares in stat.items():
            assert sum(shares) == pytest.approx(1.0, abs=1e-9), key
    assert len(product['classes']['names']) == len(product['classes']['edges_c']) + 1
    # the shares are those anyone counts from the cells: the seas' surface by the water cells' area
    sea = np.array(c['sea']) > 0
    area = np.array(c['area_km2'])
    for stat in ('mean_c', 'p99_c'):
        k = W.classify(np.array(c['sea_surface'][stat]))
        counted = [float(area[sea & (k == n)].sum() / area[sea].sum()) for n in range(len(W.CLASS_NAMES))]
        assert counted == pytest.approx(product['shares'][stat]['seas_surface'], abs=1e-6)
        still = W.classify(np.array(c['still_4m'][stat]))
        counted = [float(area[sea & (still == n)].sum() / area[sea].sum()) for n in range(len(W.CLASS_NAMES))]
        assert counted == pytest.approx(product['shares'][stat]['seas_still_4m'], abs=1e-6)
    bands = product['summary']['by_latitude']
    assert [b['band_deg'] for b in bands] == [list(b) for b in W.LATITUDE_BANDS]
    assert bands[0]['water_cells']['still_2m_p99_c'] > bands[-1]['water_cells']['still_2m_p99_c'] + 3.0


def test_product_physics_reads_as_expected(product):
    c = product['cells']
    sea = np.array(c['sea']) > 0
    area = np.array(c['area_km2'])
    lat = np.repeat(product['grid']['lat'], len(product['grid']['lon']))
    get = lambda name, key: np.array(c[name][key])
    # the seas the run carries stand near the design climate's 29-30 C, and barely swing
    assert 28.5 < product['summary']['seas']['surface_mean_c'] < 30.5
    assert np.median((get('sea_surface', 'p99_c') - get('sea_surface', 'mean_c'))[sea]) < 2.0
    # shallower still water peaks hotter, and trading water with the open water cools the peak, cell by cell
    assert np.all(get('still_1m', 'p99_c') >= get('still_2m', 'p99_c') - 0.26)
    assert np.all(get('still_2m', 'p99_c') >= get('still_4m', 'p99_c') - 0.26)
    assert np.mean((get('shore_2m_x1', 'p99_c') <= get('shore_2m_x0.3', 'p99_c') + 0.26)[sea]) > 0.99
    assert np.mean((get('shore_2m_x0.3', 'p99_c') <= get('still_2m', 'p99_c') + 0.26)[sea]) > 0.99
    # the peaks belong to the low latitudes: still 2 m shallows peak hotter within 30 degrees than poleward of 60
    low, high = np.abs(lat) < 30.0, np.abs(lat) > 60.0
    peak = get('still_2m', 'p99_c')
    assert np.average(peak[low], weights=area[low]) > np.average(peak[high], weights=area[high]) + 3.0


def test_product_agrees_with_the_five_sites(product):
    """One Sun for every cell and the per-site Suns of the forcing give the same shore water."""
    assert product['sun']['noon_spread_h'] < 1.0 and abs(product['sun']['period_h'] / F.SYNODIC_H - 1.0) < 1e-4
    for name, site in product['checks']['sites'].items():
        assert abs(site['noon_minus_site_fit_h']) < 1.5, name
        for case in ('shore_2m_x0.3', 'shore_2m_x1'):
            if case in site:
                assert site[case]['p99_c'] == pytest.approx(site[case]['forcing_p99_c'], abs=0.4), (name, case)
                assert site[case]['mean_c'] == pytest.approx(site[case]['forcing_mean_c'], abs=0.2), (name, case)
        if name != 'lake_0':                                     # the forcing's lake margin takes a sheltered wind
            assert site['still_2m_p99_c'] == pytest.approx(site['forcing_closed_2m_pan_p99_c'], abs=1.5), name


def test_a_calm_deep_waters_surface_lies_between_the_still_layers(product):
    """At the forcing's five sites the lake model's 20 m column, which resolves the day's warm layer, puts its surface
    in its hottest hundredth of hours above the open water's body and between the still 4 m and 2 m layers."""
    sites = product['checks']['sites']
    assert set(sites) == {'shore_0', 'shore_30', 'shore_45', 'shore_60', 'lake_0'}
    for name, site in sites.items():
        column = site['calm_column']
        assert column['depth_m'] == F.OPEN_LAKE_M
        assert site['still_4m_p99_c'] - 0.1 <= column['surface_p99_c'] <= site['still_2m_p99_c'] + 0.1, name
        assert column['surface_p99_c'] > site['body_p99_c'] + 1.0, name
        assert column['floor_mean_c'] < column['surface_mean_c'] < column['surface_p99_c'] <= column['surface_max_c'], name
        assert column['surface_hours_above_34'] >= column['surface_hours_above_36'] >= column['surface_hours_above_38']
    # the 10 m layer under the run's near-surface air stands within about a kelvin of the run's own sea surface
    apart = product['checks']['mixed_10m_minus_the_runs_sea_surface_means']
    assert -1.0 < apart['mean_c'] < 0.0 and -1.6 < apart['range_c'][0] <= apart['range_c'][1] < 0.5


# ---- the pieces --------------------------------------------------------------------------------------------------------
def test_bulk_exchange_for_arrays_is_the_lake_models():
    rng = np.random.default_rng(5)
    air = L.Air(288.38, 1008.9, MOON_SURFACE_GRAVITY)
    n = 400
    ts = rng.uniform(15.0, 45.0, n)
    air_k = rng.uniform(288.0, 308.0, n)
    q = rng.uniform(0.004, 0.03, n)
    wind = rng.uniform(0.0, 9.0, n)
    p = rng.uniform(9.0e4, 1.25e5, n)
    sens, latent = W.bulk_fluxes(ts, air_k, q, wind, p, air)
    stable = 0
    for i in range(n):
        s, l, _, _ = L.bulk_fluxes(float(ts[i]), float(air_k[i]), float(q[i]), float(wind[i]), float(p[i]), air)
        assert sens[i] == pytest.approx(s, rel=1e-10, abs=1e-9) and latent[i] == pytest.approx(l, rel=1e-10, abs=1e-9)
        stable += s < 0
    assert 20 < stable < n - 20                                    # both branches of the stability factor were met
    earth = L.Air(287.05, 1004.0, STANDARD_GRAVITY)
    s, l, _, _ = L.bulk_fluxes(30.0, 300.0, 0.015, 2.0, 101325.0, earth)
    sv, lv = W.bulk_fluxes(np.array([30.0]), np.array([300.0]), np.array([0.015]), np.array([2.0]), np.array([101325.0]), earth)
    assert sv[0] == pytest.approx(s, rel=1e-12) and lv[0] == pytest.approx(l, rel=1e-12)


def test_lookup_is_linear_between_the_centres_and_constant_beyond():
    rng = np.random.default_rng(2)
    tables = rng.uniform(0.0, 800.0, (5, W.CENTRES.size))
    elev = np.array([-90.0, -89.0, -12.3, 0.0, 88.7])
    got = W.lookup(tables, elev)
    for c in range(5):
        assert got[c] == pytest.approx(float(np.interp(elev[c], W.CENTRES, tables[c])), rel=1e-12)
    assert W.lookup(tables, np.full(5, 90.0)) == pytest.approx(tables[:, -1])


def test_between_weights_are_np_interps():
    n_w = 7
    centres = (np.arange(n_w) + 0.5) * F.STEPS_PER_WINDOW
    values = np.arange(n_w, dtype=float) ** 2
    for w in (0, 3, n_w - 1):
        for k in (0, 17, 35, 36, 71):
            a, b, wt = W.between(k, n_w, w)
            assert values[a] * (1 - wt) + values[b] * wt == pytest.approx(
                float(np.interp(w * F.STEPS_PER_WINDOW + k + 0.5, centres, values)), rel=1e-12, abs=1e-12)


def test_percentiles_from_counts():
    rng = np.random.default_rng(9)
    x = rng.normal(30.0, 3.0, 200000)
    counts = np.histogram(x, bins=W.BIN0 + np.arange(W.NBINS + 1) * W.BIN_W)[0][None, :]
    for q in (1.0, 50.0, 99.0):
        assert W.percentile(counts, q)[0] == pytest.approx(float(np.percentile(x, q)), abs=0.05)
    assert W.classify([27.99, 28.0, 39.99, 40.0, 55.0]).tolist() == [0, 1, 6, 7, 7]
    assert W.shares(np.array([0, 0, 7]), [1.0, 1.0, 2.0]) == [0.5, 0, 0, 0, 0, 0, 0, 0.5]


def synthetic(n_w=60, lat=(20.0, -40.0), lon=(0.0, 90.0), period=708.7, noon0=100.0):
    """A small made-up climate on a 2 x 2 grid, with windows' fields consistent with a Sun of known noon and period."""
    rng = np.random.default_rng(11)
    lat, lon = np.array(lat), np.array(lon)
    shape = (n_w, lat.size, lon.size)
    lat3 = lat[None, :, None] * np.ones(shape)
    noon = (noon0 - lon / 360.0 * period)[None, None, :] * np.ones(shape)
    steps = (np.arange(n_w * F.WINDOW_STEPS) + 0.5) * 0.5
    def window_mean(f):
        return f.reshape(n_w, F.WINDOW_STEPS, *f.shape[1:]).mean(1)
    cosz = np.cos(np.radians(lat))[None, :, None] * np.cos(2 * np.pi * (steps[:, None, None] - noon[0][None]) / period)
    czen = window_mean(np.maximum(cosz, 0.0))
    end = (np.arange(n_w) + 1.0) * F.WINDOW_H
    cos_end = np.maximum(np.cos(np.radians(lat3)) * np.cos(2 * np.pi * (end[:, None, None] - noon) / period), 0.0)
    alb = np.full(shape, 0.069)
    cloud = rng.uniform(0.5, 0.95, shape)
    data = dict(lat=lat, lon=lon, lsm=np.array([[0.0, 1.0], [0.0, 1.0]]), czen=czen.astype(np.float32),
                rsscs=(900.0 * cos_end ** 1.2 * (1 - alb)).astype(np.float32),
                rss=(900.0 * window_mean(np.maximum(cosz, 0.0) ** 1.2) * cloud * (1 - alb)).astype(np.float32),
                alb=alb.astype(np.float32), tas=(299.0 + rng.normal(0, 1.0, shape)).astype(np.float32),
                ts=(302.0 + rng.normal(0, 0.5, shape)).astype(np.float32), rls=np.full(shape, -70.0, np.float32),
                ps=np.full(shape, 1215.9, np.float32), uas=rng.uniform(0.2, 3.0, shape).astype(np.float32),
                vas=rng.uniform(-1.0, 1.0, shape).astype(np.float32),
                hus_low=rng.uniform(0.012, 0.02, shape).astype(np.float32))
    return data, dict(noon0_h=noon0, period_h=period)


def test_many_cells_at_once_match_the_forcings_one_site_slab():
    """The same made-up climate through the map's stepping and through forcing.site_forcing with forcing.slab."""
    data, sun = synthetic()
    win = W.windows(data, sun)
    sea = (data['lsm'] <= 0.5).ravel()
    cases = (('still_2m', 2.0, 0.0), ('unused', 1.0, 0.0), ('unused', 4.0, 0.0), ('mixed_10m', 10.0, 0.0),
             ('shore_2m_x0.3', 2.0, 0.3))
    out = W.run(win, sun, sea, cases=cases, spin_windows=2)
    stats = W.summarise(out)
    sky = dict(elevation=np.array([-90.0, 90.0]), shortwave=np.array([1.0, 1.0]), par0=np.array([0.0, 0.0]))
    for j in range(2):
        i = 0                                                             # the water cells
        c = j * data['lon'].size + i
        site = dict(j=j, i=i, lat=float(data['lat'][j]), lon=float(data['lon'][i]), wind_share=1.0)
        own = dict(noon_h=float(win['noon2'][c]), period_h=sun['period_h'])
        forcing, _, _ = F.site_forcing(data, site, own, sky)
        start = float(forcing.air_k[0] - 273.15)
        first = 2 * F.STEPS_PER_WINDOW
        still, _ = F.slab(forcing, 2.0, 0.0, np.zeros(len(forcing)), start)
        shore, _ = F.slab(forcing, 2.0, 0.3, F.open_sea_c(data, site, len(forcing), forcing.dt), start)
        assert stats['mean_c'][0, c] == pytest.approx(still[first:].mean(), abs=1e-6)
        assert stats['max_c'][0, c] == pytest.approx(still[first:].max(), abs=1e-6)
        assert stats['mean_c'][4, c] == pytest.approx(shore[first:].mean(), abs=1e-6)
        assert stats['min_c'][4, c] == pytest.approx(shore[first:].min(), abs=1e-6)
        assert stats['p99_c'][4, c] == pytest.approx(float(np.percentile(shore[first:], 99)), abs=W.BIN_W)
    # on a land cell the shore water trades with that cell's own 10 m water
    land = 1
    assert stats['mean_c'][4, land] != pytest.approx(stats['mean_c'][0, land], abs=1e-3)
    assert min(stats['mean_c'][0, land], stats['mean_c'][3, land]) - 0.5 < stats['mean_c'][4, land] < max(
        stats['mean_c'][0, land], stats['mean_c'][3, land]) + 0.5


def test_the_calm_column_resolves_a_warm_surface_over_a_cooler_floor():
    """The lake model's 20 m pan on the made-up climate: after the spin-up its surface swings above its floor, and in
    its hottest hours stands above a 10 m mixed layer under the same forcing."""
    data, sun = synthetic(n_w=60)
    win = W.windows(data, sun)
    site = dict(j=0, i=0, lat=float(data['lat'][0]), lon=float(data['lon'][0]),
                sun=dict(noon_h=float(win['noon2'][0]), period_h=sun['period_h']))
    column = W.calm_column(data, site)
    assert set(column) == {'depth_m', 'wind_mean_m_s', 'surface_mean_c', 'surface_p99_c', 'surface_max_c', 'floor_mean_c',
                           'surface_hours_above_34', 'surface_hours_above_36', 'surface_hours_above_38'}
    assert column['floor_mean_c'] < column['surface_mean_c'] < column['surface_p99_c'] <= column['surface_max_c']
    out = W.run(win, sun, (data['lsm'] <= 0.5).ravel(), cases=(('a', 2.0, 0.0), ('b', 4.0, 0.0), ('c', 6.0, 0.0), ('mixed_10m', 10.0, 0.0)),
                spin_windows=F.SPIN_WINDOWS)
    assert column['surface_p99_c'] > W.summarise(out)['p99_c'][3, 0]


def test_coast_takes_its_neighbouring_water_cells():
    sea2 = np.zeros((3, 6), bool)
    sea2[0, 0] = sea2[2, 5] = True
    values = np.arange(18, dtype=float).reshape(3, 6)
    got, beside = W.coast_values(values, sea2, [0, 1, 1, 2], [1, 3, 5, 0])
    assert beside.tolist() == [True, False, True, True]
    assert got[0] == 0.0                      # only the water cell to its west
    assert got[1] == 9.0                      # no water cell beside it: its own value
    assert got[2] == pytest.approx((0.0 + 17.0) / 2)     # the wrap finds (0, 0), and (2, 5) lies below it
    assert got[3] == 17.0                     # across the longitude's wrap
