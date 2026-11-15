"""The cycle screen: a known sinusoid in red noise is found at its period, pure red noise is flagged near the nominal
rate within a series and across a family, the clock reads a known solar day, the aliasing and the multiple-testing
arithmetic are exact, and the two-timescale null recovers a known memory and holds a long-memory series that AR(1)
flags."""
import json
from pathlib import Path

import numpy as np
import pytest

from climate.gcm import cycles

SOLAR = cycles.solar_day(27.321661, 365.25636)
ANNUAL = cycles.forcing_year_days(365.25636)
PRODUCT = cycles.RESULTS / 'cycles_A28R_301k_design_l15_snow_s924_turned_refit_y10-49.json'


def times(years):
    return (np.arange(years * cycles.RECORDS) + 0.5) * cycles.DT


def test_the_forcing_year_and_solar_day_follow_the_configuration():
    assert ANNUAL == 365.25
    assert SOLAR == pytest.approx(29.530588, abs=1e-6)
    assert cycles.FILE_YEAR == 357.5


def test_a_sinusoid_in_red_noise_is_found_at_its_period():
    rng = np.random.default_rng(11)
    t = times(20)
    y = cycles.ar1(0.5, len(t), 1, rng)[0] + 0.8 * np.cos(2 * np.pi * t / 500.0 + 1.0)
    sc = cycles.screen(y[:, None], t, SOLAR, ANNUAL, 400, rng)
    p_series, _, rows = cycles.peaks(sc, [('x', 'test')], t)
    best = min(rows, key=lambda r: r['p_series'])
    assert best['period_days'] == pytest.approx(500.0, rel=0.03)
    assert best['p_series'] <= 1 / 401 + 1e-12
    assert p_series[0] <= 1 / 401 + 1e-12
    assert best['amplitude'] == pytest.approx(0.8, rel=0.25)


def test_pure_red_noise_is_flagged_near_the_nominal_rate():
    rng = np.random.default_rng(12)
    t = times(20)
    rhos = rng.uniform(0.1, 0.8, 300)
    Y = np.stack([cycles.ar1(r, len(t), 1, rng)[0] for r in rhos], axis=1)
    sc = cycles.screen(Y, t, SOLAR, ANNUAL, 200, rng, batch=200)
    p_series, _, _ = cycles.peaks(sc, [('x', str(i)) for i in range(len(rhos))], t)
    rate = float(np.mean(p_series <= cycles.ALPHA))
    assert 0.02 <= rate <= 0.09, rate
    # Across families of ten series, Holm keeps the chance of any detection near or below the level.
    families = p_series.reshape(30, 10)
    any_detection = sum(bool((cycles.holm(p) <= cycles.ALPHA).any()) for p in families)
    assert any_detection <= 4, any_detection


def test_the_clock_reads_a_known_solar_day():
    lat = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(32)[0]))[::-1]
    lon = np.arange(64) * 5.625
    t = times(3)
    sub = 40.0 - 360.0 * t / SOLAR
    czen = np.maximum(0, np.cos(np.radians(lat))[None, :, None]
                      * np.cos(np.radians(lon[None, None, :] - sub[:, None, None])))
    sun = czen[:, np.abs(lat) < 20, :].mean(axis=1) @ np.exp(1j * np.radians(lon))
    observed, rms = cycles.clock(t, sun)
    assert observed == pytest.approx(SOLAR, rel=1e-9)
    assert rms < 1e-6


def test_the_records_alias_the_solar_harmonics_where_the_screen_found_them():
    assert 1 / cycles.alias(8 / SOLAR) == pytest.approx(15.444, abs=0.01)
    assert 1 / cycles.alias(9 / SOLAR) == pytest.approx(32.37, abs=0.01)
    assert 1 / cycles.alias(10 / SOLAR) == pytest.approx(336.97, abs=0.05)
    assert 1 / cycles.alias(11 / SOLAR) == pytest.approx(27.151, abs=0.01)
    name, f, _ = cycles.nearest(1 / 15.46, cycles.catalogue(SOLAR, ANNUAL))
    assert name == 'solar 8 (aliased)'


def test_holm_and_benjamini_hochberg_adjust_as_published():
    p = np.array([0.01, 0.04, 0.03, 0.005])
    assert np.allclose(cycles.holm(p), [0.03, 0.06, 0.06, 0.02])
    assert np.allclose(cycles.bh(p), [0.02, 0.04, 0.04, 0.02])


def test_the_file_year_composite_sees_a_restart_spike():
    rng = np.random.default_rng(13)
    r = cycles.ar1(0.3, 40 * cycles.RECORDS, 2, rng)
    r[1, ::cycles.RECORDS] += 3.0           # one record a year, three standard deviations of the noise
    stat = cycles.file_year_stat(r)
    assert stat[1] > 3 * stat[0]


def test_the_two_timescale_fit_recovers_a_known_memory():
    rng = np.random.default_rng(14)
    x = cycles.two_ar1(0.6, 0.4, 0.99, 40 * cycles.RECORDS, 1, rng)[0]
    w, fast, slow, err = cycles.fit_memory(x)
    assert w == pytest.approx(0.6, abs=0.15)
    assert fast == pytest.approx(0.4, abs=0.15)
    assert slow == pytest.approx(0.99, abs=0.01)
    assert err < 0.05


def test_long_memory_passes_the_two_timescale_null_where_ar1_flags_it():
    rng = np.random.default_rng(15)
    t = times(40)
    y = cycles.two_ar1(0.5, 0.5, 0.995, len(t), 1, rng)[0]
    sc = cycles.screen(y[:, None], t, SOLAR, ANNUAL, 300, rng)
    p_ar1, observed_max, _ = cycles.peaks(sc, [('x', 'memory')], t)
    ms = cycles.memory_screen(sc, t, SOLAR, ANNUAL, 300, rng)
    assert p_ar1[0] < 0.01
    assert cycles.exceed(ms['smax'], observed_max)[0] > 0.05


@pytest.mark.skipif(not PRODUCT.exists(), reason='product not built')
def test_the_product_reads_the_runs_clock_and_records_its_inputs():
    product = json.loads(PRODUCT.read_text())
    assert product['schema'] == cycles.SCHEMA
    clock = product['clock']
    assert clock['solar_days_observed'] == pytest.approx(clock['solar_days_config'], rel=1e-5)
    assert clock['solar_phase_rms_degrees'] < 0.1
    assert clock['forcing_year_best_fit_days'] == pytest.approx(365.25, abs=0.05)
    first, last = product['years']
    assert sorted(int(y) for y in product['inputs']['year_files']) == list(range(first, last + 1))
    assert product['evidence'] and product['reading_rule']
    assert set(product['windows']) == {'main', 'compare', 'first half', 'second half'}
    assert 'memory_null' in product['windows']['main'] and 'memory_null' in product['windows']['second half']
    assert len(product['memory']['main']['fits']) == product['series']['count'] == 234
