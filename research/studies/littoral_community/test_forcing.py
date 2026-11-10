"""Tests for the littoral's physical setting: the product's binding and its own checks, the coast mask, the Sun fitted
to window means, the run's clear sky read from end-of-window values, and the exchanging slab against its limits."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from biosphere.ecology import lake as L
from biosphere.ecology import littoral as M
from research.studies.littoral_community import forcing as F
from shared.provenance import constants_changed

PRODUCT = F.OUT


@pytest.fixture(scope='module')
def product():
    if not PRODUCT.exists():
        pytest.skip('product not built')
    return json.loads(PRODUCT.read_text())


# ---- binding ---------------------------------------------------------------------------------------------------------
def test_product_schema_files_and_constants(product):
    assert product['schema'] == F.SCHEMA
    for path, digest in product['producer']['files'].items():
        assert F.sha256(F.ROOT / path) == digest, path
    assert constants_changed(product['producer']['constants']) == {}


def test_design_air_matches_the_waters_oxygen_study(product):
    current = dict(zip(('air_gas_constant', 'air_cp', 'o2_fraction', 'salinity_g_kg'), F.air_constants()))
    assert product['design_air'] == pytest.approx(current, rel=1e-12)


def test_series_on_the_drive_are_the_recorded_ones(product):
    for name, entry in product['series'].items():
        path = F.Path(entry['path'])
        if not path.exists():
            pytest.skip('series not on this machine')
        assert F.sha256(path) == entry['sha256'], name


def test_product_checks(product):
    assert set(product['sites']) == {'shore_0', 'shore_30', 'shore_45', 'shore_60', 'lake_0'}
    for name, site in product['sites'].items():
        assert site['sun']['rms'] < 0.02, name                       # the fitted Sun against the windows' cosines
        assert abs(site['sun']['period_h'] / F.SYNODIC_H - 1.0) < 1e-4, name
        c = site['light_checks']
        assert c['window_shortwave_kept'] < 1e-6, name               # each lit window keeps the run's sunlight
        # the clear-sky diagnostic is each window's last value: it follows the Sun at the window's end, falls to
        # twilight once the Sun has set by then, and a window mean would correlate far less
        assert c['clear_sky_follows_window_end'] > 0.99 > 0.95 > c['clear_sky_follows_window_mean'], name
        assert c['clear_sky_after_sunset_w'] < 20.0, name
        lo, hi = c['run_clear_over_solved_channels']
        assert 1.3 < lo <= hi < 1.9, name                            # the solved channels hold 0.5-0.75 of the shortwave
        assert 0.5 < c['cloud_factor_ten_years'] < 0.9 and c['lit_windows_above_cap'] < 0.05, name
        s = site['setting']
        assert 100.0 < s['dark_hours_per_lunar_day'] < 280.0, name
        assert 5.0 < s['par0_mol_m2_d'] < 35.0, name
        for key, case in site['cases'].items():
            assert 15.0 < case['water_min_c'] <= case['water_mean_c'] <= case['water_max_c'] < 50.0, (name, key)
            assert 0.1 < case['k_m_day_mean'] < 5.0, (name, key)
    # darkness shortens toward the pole, and the equatorial shore's dark hours stand near the clear sky's 255
    dark = [product['sites'][k]['setting']['dark_hours_per_lunar_day'] for k in ('shore_0', 'shore_30', 'shore_45', 'shore_60')]
    assert dark == sorted(dark, reverse=True) and 250.0 < dark[0] < 275.0


def test_a_site_series_drives_the_shelved_simulator(product):
    """The forcing's series feed the community model as they stand: two boxes at one shore through a lunar day, with
    carbon, phosphorus and oxygen closing."""
    path = F.Path(product['series']['shore_45']['path'])
    if not path.exists():
        pytest.skip('series not on this machine')
    with np.load(path) as z:
        dt = float(z['dt_s']) / 86400.0
        steps = int(round(29.53 / dt))
        cases = ('slab_2m_x0.3', 'slab_4m_x0.3')
        both = lambda key: np.repeat(z[key][:steps, None], 2, 1)
        forcing = dict(t=np.stack([z[c + '_t'][:steps] for c in cases], 1), t_open=both('open_c'), par0=both('par0'),
                       k=np.stack([z[c + '_k'][:steps] for c in cases], 1), wind=both('wind10'))
    setting = M.Setting(depth=np.array([2.0, 4.0]), exchange=0.3, tp_open=8.0,
                        o2_fraction=product['design_air']['o2_fraction'], salinity=product['design_air']['salinity_g_kg'])
    model = M.Littoral(M.example_community(idle_share=0.25), setting, 2)
    state = model.initial()
    start = model.totals({k: v.copy() for k, v in state.items()})
    seeds = [(0.0, 'phyto', 0.05), (0.0, 'peri', 0.5), (0.0, 'chara', 5.0), (0.0, 'snail', 0.3), (0.0, 'cray_a', 0.3),
             (0.0, 'fish_a', 0.1)]
    out = M.drive(model, state, forcing, dt, steps, seeds, record_every=24)
    a, brought = out['accounts'], out['seeds']
    carbon, phosphorus, _, oxygen = model.totals(state)
    assert np.allclose(carbon - start[0], brought['c'] + a['gpp_c'] - a['resp_c'] - a['burial_c'] - a['emerged_c']
                       - a['harvest_c'] - a['ch4_c'] + a['c_exchange'], atol=1e-9)
    assert np.allclose(phosphorus - start[1], brought['p'] - a['p_calcite'] - a['burial_p'] - a['emerged_p']
                       - a['harvest_p'] + a['p_exchange'], atol=1e-8)
    assert np.allclose(oxygen - start[3], a['o2_bio'] + a['o2_air'] - a['o2_bubbles'] + a['o2_exchange'], atol=1e-9)
    assert np.all(out['accounts']['gpp_c'] > 0) and np.all(out['o2_low'] >= 0.0)


def test_exchange_with_the_sea_tames_the_shore(product):
    """More exchange and more depth both narrow the water's swing toward the open water's."""
    for name in ('shore_0', 'shore_30', 'shore_45', 'shore_60'):
        cases = product['sites'][name]['cases']
        for depth in ('1', '2', '4'):
            peaks = [cases[f'slab_{depth}m_x{x}']['water_p99_c'] for x in ('0.1', '0.3', '1')]
            assert peaks == sorted(peaks, reverse=True), (name, depth)
        for x in ('0.1', '0.3', '1'):
            peaks = [cases[f'slab_{d}m_x{x}']['water_p99_c'] for d in ('1', '2', '4')]
            assert peaks == sorted(peaks, reverse=True), (name, x)


# ---- the pieces --------------------------------------------------------------------------------------------------------
def test_coast_is_sea_beside_land_and_wraps_in_longitude():
    lsm = np.zeros((4, 6))
    lsm[1, 0] = 1.0                                   # one land cell at the grid's western edge
    coast = F.coastal(lsm)
    assert not coast[1, 0]                            # land is no coast
    assert coast[1, 1] and coast[1, 5] and coast[0, 0] and coast[2, 5] and coast[0, 5]
    assert not coast[3, 0] and not coast[1, 3]
    assert coast.sum() == 8
    # a land cell in the top row has no neighbours across the pole
    top = np.zeros((4, 6))
    top[0, 2] = 1.0
    assert F.coastal(top)[3].sum() == 0 and F.coastal(top).sum() == 5


def test_median_cell_takes_the_lower_middle():
    ts = np.zeros((3, 2, 2), dtype=np.float32) + 273.15
    ts[:, 0, 0] += 10.0
    ts[:, 0, 1] += 30.0
    ts[:, 1, 0] += 20.0
    ts[:, 1, 1] += 40.0
    j, i, t, n = F.median_cell(dict(ts=ts), np.ones((2, 2), bool))
    assert (j, i, n) == (1, 0, 4) and t == pytest.approx(20.0, abs=1e-4)


def test_window_means_of_the_suns_cosine():
    # over whole lunar days the mean of max(cos, 0) is cos(lat) / pi
    n = 400
    lat = 30.0
    period = F.WINDOW_H * 10.0                         # ten windows to the day, so means repeat
    means = F.window_mean_cosine(lat, 0.0, period, n)
    assert means.mean() == pytest.approx(math.cos(math.radians(lat)) / math.pi, rel=1e-3)
    assert means.max() < math.cos(math.radians(lat)) and means.min() == 0.0


def test_the_sun_fit_recovers_a_known_noon_and_period():
    pytest.importorskip('scipy')
    lat, noon, period = 47.0, 301.7, 708.9
    czen = F.window_mean_cosine(lat, noon, period, 1200)
    rng = np.random.default_rng(3)
    fit = F.fit_sun(czen + rng.normal(0.0, 1e-3, czen.size), lat)
    assert fit['period_h'] == pytest.approx(period, abs=0.002)
    assert (fit['noon_h'] - noon + period / 2) % period - period / 2 == pytest.approx(0.0, abs=0.2)
    assert fit['rms'] < 2e-3


def test_elevation_follows_the_fitted_sun():
    sun = dict(noon_h=100.0, period_h=708.7)
    e, hour = F.elevation_at(20.0, np.array([100.0, 100.0 + 708.7 / 4, 100.0 + 708.7 / 2]), sun)
    assert e[0] == pytest.approx(70.0) and e[1] == pytest.approx(0.0, abs=1e-9) and e[2] == pytest.approx(-70.0)
    assert hour[0] == pytest.approx(0.0) and hour[1] == pytest.approx(90.0)


def test_run_clear_sky_reads_bin_medians_by_elevation():
    elevation = np.array([-40.0, -39.5, 10.2, 10.9, 11.5, 60.0])
    flux = np.array([0.0, 0.2, 100.0, 120.0, 500.0, 900.0])
    centres, medians = F.run_clear_sky(flux, elevation)
    assert centres.tolist() == [-39.0, 11.0, 61.0]
    assert medians.tolist() == [0.1, 120.0, 900.0]


def steady(steps, dt=3575.0, sw=0.0, lw=420.0, air_k=302.0, q=0.018, wind=2.0, p=121590.0):
    full = lambda x: np.full(steps, float(x))
    return L.Forcing(dt=dt, sw_down=full(sw), lw_down=full(lw), air_k=full(air_k), q_air=full(q), wind10=full(wind),
                     pressure=full(p), o2_fraction=0.1746)


def test_slab_follows_the_open_water_under_fast_exchange_and_its_own_balance_without():
    forcing = steady(2000, sw=300.0)
    open_c = np.full(2000, 29.0)
    fast, k = F.slab(forcing, 1.0, 200.0, open_c, 20.0)
    assert fast[-1] == pytest.approx(29.0, abs=0.05) and np.all(k > 0)
    alone, _ = F.slab(forcing, 1.0, 0.0, open_c, 29.0)
    some, _ = F.slab(forcing, 1.0, 0.3, open_c, 29.0)
    # the closed slab settles where its surface balance closes; exchange holds it between that and the open water
    assert abs(alone[-1] - alone[-200]) < 0.05
    assert min(alone[-1], 29.0) < some[-1] < max(alone[-1], 29.0)
    deep, _ = F.slab(forcing, 4.0, 0.0, open_c, 29.0)
    assert abs(deep[50] - 29.0) < abs(alone[50] - 29.0)            # a deeper slab warms more slowly


def test_slab_matches_the_lake_model_on_a_mixed_pan():
    """Under steady forcing a closed slab and the lake model's pan of the same depth settle at the same temperature."""
    forcing = steady(3000, sw=250.0, wind=4.0)
    alone, _ = F.slab(forcing, 1.0, 0.0, np.full(3000, 29.0), 29.0)
    out = F.run_water(forcing, 1.0, 29.0)
    pan = 0.5 * (out['surface_t'][-1] + out['bottom_t'][-1])
    assert alone[-1] == pytest.approx(pan, abs=0.6)
