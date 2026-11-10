"""Tests for the waters under a year's 6.5-hour records against its 3-day means: the product's binding and what it
shows, the record brought to windows, and the passes' winds on a made-up record."""
from __future__ import annotations

import json

import numpy as np
import pytest

from research.studies.littoral_community import forcing as F
from research.studies.littoral_community import inside_windows as I
from research.studies.littoral_community import water_map as W
from research.studies.littoral_community import wind_record as R
from shared.provenance import constants_changed

PRODUCT = I.OUT


@pytest.fixture(scope='module')
def product():
    if not PRODUCT.exists():
        pytest.skip('product not built')
    return json.loads(PRODUCT.read_text())


# ---- the product -----------------------------------------------------------------------------------------------------
def test_product_schema_files_constants_and_inputs(product):
    assert product['schema'] == I.SCHEMA
    for path, digest in product['producer']['files'].items():
        assert F.sha256(F.ROOT / path) == digest, path
    assert constants_changed(product['producer']['constants']) == {}
    assert F.sha256(W.FORCING) == product['inputs'][str(W.FORCING.relative_to(F.ROOT))]
    for path in (I.fine_file(), I.parent_file(), I.parent_file(I.PARENT_YEAR + 1)):
        if not path.exists():
            pytest.skip('the year files are not on this machine')
        assert F.sha256(path) == product['inputs'][str(path.relative_to(F.RUNS))], path


def test_the_record_is_the_design_years_own_weather(product):
    """The two runs end the year on one restart, and the record's fields, brought to windows as each field is, are
    the design year's 3-day records: the surface fields to rounding, the lowest level's air and wind to a few parts
    in a thousand of their largest values."""
    checks = product['checks']
    assert checks['restart_identical'] is True
    agree = checks['agreement']
    assert set(agree) == set(I.MEANS + I.LAST + I.FIRST)
    for k in ('tas', 'ts', 'rss', 'rls', 'alb', 'uas', 'vas', 'czen', 'rsscs'):
        assert agree[k]['largest'] < 1e-5, k
    for k in ('ps', 'ta_low', 'ua_low', 'va_low', 'hus_low'):
        assert agree[k]['rms'] < 1e-2 and agree[k]['largest'] < 0.2, k
    s = product['settings']
    assert s['record_h'] == 6.5 and s['records_per_window'] * s['record_h'] == s['window_h'] == F.WINDOW_H
    assert s['design_run'] == F.RUN and s['design_year'] == 10 and s['wind_raise'] == pytest.approx(F.wind_raise())


def test_a_year_files_first_humidity_is_empty_and_the_fill_stands_near_the_record(product):
    h = product['checks']['humidity_at_the_years_ends']
    assert 0.01 < h['as_written_over_the_next_record'] < 0.04
    for key in ('start_fill_against_6h5_inside', 'end_fill_against_6h5_inside'):
        assert h[key]['rms'] < 0.15 and abs(h[key]['mean']) < 0.05, key
        assert h[key]['rms'] < 2.0 * h['a_window_of_real_change']['rms'] + 0.02, key


def test_the_passes_stand_in_their_order(product):
    """The unraised means leave the water warmest, and each account of the wind's speed inside the windows cools it:
    the true window speeds bring the still water nearer the 6.5-hour forcing than the unraised means do."""
    assert set(product['passes']) == set(I.PASSES) - {I.REFERENCE}
    wind = product['wind_at_the_level']
    for kind in ('water', 'land'):
        assert wind['means_unraised'][kind] < wind['means_true_speed'][kind]
        assert wind['means_true_speed'][kind] == pytest.approx(wind['six_hourly'][kind], rel=1e-3)
        assert wind['means_cell_factor'][kind] == pytest.approx(wind['means_true_speed'][kind], rel=1e-3)
        assert wind['means_cell_excess'][kind] == pytest.approx(wind['means_true_speed'][kind], rel=1e-3)
        assert wind['means_raised'][kind] == pytest.approx(wind['means_unraised'][kind] * F.wind_raise(), rel=1e-3)
    for case in ('still_2m', 'still_4m'):
        for kind in ('water_cells', 'land_cells'):
            unraised = product['passes']['means_unraised'][case]['mean_c'][kind]['mean']
            true = product['passes']['means_true_speed'][case]['mean_c'][kind]['mean']
            assert unraised > true and abs(true) < abs(unraised), (case, kind)
    for name, rows in product['passes'].items():
        assert set(rows) == set(I.REPORT)
        for case in rows:
            for key in I.KEYS:
                for kind in ('water_cells', 'land_cells'):
                    e = rows[case][key][kind]
                    assert len(e['by_band']) == len(W.LATITUDE_BANDS) and e['rms'] >= abs(e['mean']) - 1e-3


def test_the_sites_carry_every_pass(product):
    forcing = json.loads(W.FORCING.read_text())
    assert set(product['sites']) == set(forcing['sites'])
    for site in product['sites'].values():
        assert set(site['passes']) == set(I.PASSES) and site['wind_factor'] > 1.0 and site['wind_excess_m_s'] > 0.0


# ---- the record brought to windows -------------------------------------------------------------------------------------
def made_up(n_w=24, lat=(20.0, -40.0), lon=(0.0, 90.0), period=708.7, noon0=100.0, seed=5):
    """A made-up record on a 2 x 2 grid, eleven records to a window, consistent with a Sun of known noon and period,
    with a wind that gusts inside its windows; and its 3-day record, each field aggregated as the field is."""
    rng = np.random.default_rng(seed)
    per, write = I.PER, I.FINE_STEPS
    n = n_w * per
    lat, lon = np.array(lat), np.array(lon)
    shape = (n, lat.size, lon.size)
    noon = (noon0 - lon / 360.0 * period)[None, None, :]
    steps = (np.arange(n * write) + 0.5) * 0.5
    cosz = np.maximum(np.cos(np.radians(lat))[None, :, None] * np.cos(2 * np.pi * (steps[:, None, None] - noon) / period), 0.0)
    mean = lambda f: f.reshape(n, write, *f.shape[1:]).mean(1)
    end = (np.arange(n) + 1.0) * I.FINE_H
    cos_end = np.maximum(np.cos(np.radians(lat))[None, :, None] * np.cos(2 * np.pi * (end[:, None, None] - noon) / period), 0.0)
    alb = np.full(shape, 0.069)
    cloud = np.repeat(rng.uniform(0.5, 0.95, (n_w,) + shape[1:]), per, 0)
    slow = lambda lo, hi: np.repeat(rng.uniform(lo, hi, (n_w,) + shape[1:]), per, 0)
    fine = dict(lat=lat, lon=lon, lsm=np.array([[0.0, 1.0], [0.0, 1.0]]), sigma_low=0.99415, year_starts=np.array([0]),
                czen=mean(cosz), rsscs=900.0 * cos_end ** 1.2 * (1 - alb), rss=900.0 * mean(cosz ** 1.2) * cloud * (1 - alb),
                alb=alb, tas=slow(298.0, 300.0), ts=slow(301.5, 302.5), rls=np.full(shape, -70.0),
                ps=np.full(shape, 1215.9), uas=slow(0.2, 3.0), vas=slow(-1.0, 1.0), hus_low=slow(0.012, 0.02),
                ta_low=slow(296.0, 298.0), ua_low=slow(0.2, 4.0) + rng.normal(0.0, 1.2, shape),
                va_low=slow(-1.5, 1.5) + rng.normal(0.0, 1.2, shape))
    coarse = dict(I.aggregate(fine), lat=lat, lon=lon, lsm=fine['lsm'], sigma_low=fine['sigma_low'],
                  year_starts=np.array([0]))
    return fine, coarse, dict(noon0_h=noon0, period_h=period)


def test_aggregate_takes_each_field_as_it_is_written():
    per = I.PER
    ramp = np.arange(3 * per, dtype=float)[:, None, None] * np.ones((1, 2, 2))
    fine = {k: ramp for k in I.MEANS + I.LAST + I.FIRST}
    out = I.aggregate(fine)
    assert out['tas'][:, 0, 0].tolist() == pytest.approx([5.0, 16.0, 27.0])            # a window's mean
    assert out['rsscs'][:, 0, 0].tolist() == [10.0, 21.0, 32.0]                         # its last record
    assert out['hus_low'][:, 0, 0].tolist() == [0.0, 11.0, 22.0]                        # its opening record
    # a record made of a 3-day record repeated agrees with it, whatever its first humidity
    coarse = {k: (v + 1.0) for k, v in out.items()}
    again = {k: np.repeat(v, per, 0) for k, v in coarse.items()}
    again['hus_low'][0] = 0.0
    assert max(v['largest'] for v in I.agreement(again, coarse).values()) == 0.0
    coarse['tas'] = coarse['tas'] * 1.001
    both = I.agreement(again, coarse)['tas']
    assert both['largest'] == pytest.approx(0.001 / 1.001, rel=1e-6) and 0.0 < both['rms'] < both['largest']
    assert both['scale'] == pytest.approx(float(np.abs(coarse['tas']).max()))


def test_the_humiditys_fill_is_read_against_the_record():
    fine, coarse, _ = made_up()
    q = np.array(coarse['hus_low'])
    fine['hus_low'] = np.array(fine['hus_low'])
    fine['hus_low'][0] *= 0.02                                                        # as the run writes a year's first
    coarse['hus_low'] = q.copy()
    coarse['hus_low'][0] *= 0.02
    h = I.humidity_at_the_years_ends(fine, coarse, q[1])
    assert h['as_written_over_the_next_record'] == pytest.approx(0.02 * (q[0] / q[1]).mean(), rel=1e-3)
    assert set(h) == {'as_written_over_the_next_record', 'start_fill_against_6h5_inside', 'end_fill_against_6h5_inside',
                      'a_window_of_real_change'}
    assert 0.0 < h['start_fill_against_6h5_inside']['rms'] < 0.5


# ---- the passes ------------------------------------------------------------------------------------------------------
def test_the_passes_take_the_winds_they_name():
    fine, coarse, sun = made_up()
    wins, wind = I.build_passes(fine, coarse, sun)
    assert set(wins) == set(I.PASSES)
    speed, vector = R.windows(fine['ua_low'], fine['va_low'], I.PER)
    inside, outside, factor = R.floored(speed, vector)
    flat = lambda a: a.reshape(a.shape[0], -1)
    assert wins['six_hourly']['n_w'] == I.PER * wins['means_raised']['n_w']
    assert wins['six_hourly']['window_h'] == 6.5 and wins['six_hourly']['steps_per_window'] == 13
    assert wins['means_raised']['window_h'] == F.WINDOW_H and wins['means_raised']['steps_per_window'] == F.STEPS_PER_WINDOW
    assert wins['means_unraised']['level'] == pytest.approx(flat(outside))
    assert wins['means_raised']['level'] == pytest.approx(flat(outside) * F.wind_raise())
    assert wins['means_true_speed']['level'] == pytest.approx(flat(inside))
    six = wins['six_hourly']['level'].reshape(-1, I.PER, flat(inside).shape[1]).mean(1)
    assert six == pytest.approx(flat(inside))                                      # the record's own floored speeds
    for name in ('means_cell_factor', 'means_cell_excess'):
        assert wins[name]['level'].mean(0) == pytest.approx(flat(inside).mean(0), rel=1e-9)
    assert np.all(wind['factor'] > 1.0) and np.all(wind['excess'] > 0.0)
    # every pass reads one Sun and one set of clear-sky tables, and the humidity's windows keep the same means
    for name in I.PASSES:
        assert wins[name]['tables'] is wins['six_hourly']['tables']
    assert wins['six_hourly']['q'][1:-1].mean() == pytest.approx(wins['means_raised']['q'][1:-1].mean(), rel=0.02)


def test_the_waters_under_the_record_and_under_its_means_stand_close_and_in_order():
    """On the made-up record the still water under the unraised means runs warmer than under the record, and the true
    window speeds bring it nearer."""
    fine, coarse, sun = made_up(n_w=40)
    wins, _ = I.build_passes(fine, coarse, sun)
    sea = (coarse['lsm'] <= 0.5).ravel()
    cases = (('still_2m', 2.0, 0.0), ('still_4m', 4.0, 0.0), ('unused', 1.0, 0.0), ('mixed_10m', 10.0, 0.0))
    mean = {}
    for name in ('six_hourly', 'means_unraised', 'means_true_speed'):
        spin = 4 * (I.PER if name == 'six_hourly' else 1)
        mean[name] = W.summarise(W.run(wins[name], sun, sea, cases=cases, spin_windows=spin))['mean_c'][0]
    assert np.all(np.isfinite(mean['six_hourly'])) and np.all(np.abs(mean['means_unraised'] - mean['six_hourly']) < 3.0)
    assert np.all(mean['means_unraised'] > mean['means_true_speed'])
    assert np.abs(mean['means_true_speed'] - mean['six_hourly']).mean() < np.abs(mean['means_unraised'] - mean['six_hourly']).mean()
