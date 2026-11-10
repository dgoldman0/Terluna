"""Tests for the wind inside a 3-day mean: the product's binding, the ratios, the floor and the quadrature excess on
made-up winds, and what the two records show of the raise the forcing takes."""
from __future__ import annotations

import json

import numpy as np
import pytest

from research.studies.littoral_community import forcing as F
from research.studies.littoral_community import wind_record as R
from shared.provenance import constants_changed

PRODUCT = R.OUT


@pytest.fixture(scope='module')
def product():
    if not PRODUCT.exists():
        pytest.skip('product not built')
    return json.loads(PRODUCT.read_text())


def test_product_schema_files_constants_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, digest in product['producer']['files'].items():
        assert F.sha256(F.ROOT / path) == digest, path
    assert constants_changed(product['producer']['constants']) == {}
    assert F.sha256(F.LAKE_OXYGEN) == product['inputs'][str(F.LAKE_OXYGEN.relative_to(F.ROOT))]
    assert set(product['records']) == set(R.RECORDS) == {'corrected', 'earlier'}
    assert product['settings']['write_h'] == 6.5 and product['settings']['write_steps'] * R.WRITES_PER_WINDOW == F.WINDOW_STEPS
    for key in R.RECORDS:
        for path in R.record_files(key):
            if not path.exists():
                pytest.skip('the 6.5-hour records are not on this machine')
            assert F.sha256(path) == product['inputs'][str(path.relative_to(F.RUNS))], path


def test_the_corrected_record_is_a_design_year_at_the_level_the_exchange_is_paired_at(product):
    corrected, earlier = product['records']['corrected'], product['records']['earlier']
    assert corrected['run'] == F.RUN + '_eddies' and corrected['windows'] == 120 and earlier['windows'] == 240
    assert 250.0 < corrected['lowest_level_height_m'] < 400.0 < 700.0 < earlier['lowest_level_height_m'] < 1100.0
    assert 'at_10_m' in corrected and 'at_10_m' not in earlier


def test_the_forcings_raise_is_the_waters_oxygen_studys_and_sits_in_the_earlier_grids_range(product):
    factors = F.wind_factors()
    assert len(factors) == 6 and all(1.05 < x < 1.4 for x in factors)
    assert F.wind_raise() == pytest.approx(float(np.mean(factors)), rel=1e-12)
    assert product['the_forcings_raise'] == dict(factor=pytest.approx(F.wind_raise(), rel=1e-12), cells=factors)
    grid = product['records']['earlier']['lowest_level']['grid']
    cube, linear = grid['cube_ratio'], grid['mean_speed_over_vector']
    assert cube['cells']['p05'] < F.wind_raise() < cube['cells']['p95']
    assert abs(F.wind_raise() - cube['area_mean']) < 0.05
    # the cube's ratio stands at or above the mean speed's over nearly all the grid, and both above one
    assert 1.0 < linear['cells']['lowest'] and 1.0 < cube['cells']['lowest']
    assert linear['area_mean'] < cube['area_mean'] < 1.3


def test_calm_windows_hold_more_beyond_their_vector_than_one_factor_gives(product):
    """In every record the mean speed's ratio to the vector's falls from well above the forcing's factor in the calmest
    windows toward one in the windiest, so one factor raises calm windows too little and windy ones too much. The
    earlier record's excess, about 0.9 km up, is nearly one figure in every class; at the corrected run's 313 m less
    wind hides in a window, and it grows with the wind, though by less than the ratio falls."""
    shown = lambda key, wind: [r for r in product['records'][key][wind]['by_vector_speed'] if r['share_of_windows'] > 0.01]
    for key, wind, calmest in (('corrected', 'lowest_level', 2.0), ('corrected', 'at_10_m', 1.4), ('earlier', 'lowest_level', 2.0)):
        table = product['records'][key][wind]['by_vector_speed']
        assert len(shown(key, wind)) >= 4 and sum(r['share_of_windows'] for r in table) == pytest.approx(1.0, abs=1e-3)
        ratios = [r['ratio'] for r in shown(key, wind)]
        assert ratios == sorted(ratios, reverse=True) and ratios[0] > calmest and ratios[-1] < 1.1, (key, wind)
        assert ratios[0] > F.wind_raise() > ratios[-1], (key, wind)
    early = [r['excess_rms_m_s'] for r in shown('earlier', 'lowest_level')]
    assert max(early) - min(early) < 0.1 * np.mean(early)
    late = [r['excess_rms_m_s'] for r in shown('corrected', 'lowest_level')]
    assert late == sorted(late) and 1.3 < late[-1] / late[0] < 2.0
    assert max(late) < min(early)


def test_the_corrected_records_own_factor_stands_below_the_forcings_raise(product):
    """At the level the exchange is paired at, the design year's own windows need a smaller raise than the factor
    borrowed from the earlier run's level about 0.9 km up, over water and over land alike."""
    own = product['records']['corrected']['lowest_level']['floored']
    for kind in ('water', 'land'):
        assert 1.05 < own['factor'][kind] < F.wind_raise() - 0.05
        assert own['quadrature_excess_m_s'][kind] < product['records']['earlier']['lowest_level']['floored']['quadrature_excess_m_s'][kind]
    assert product['records']['earlier']['lowest_level']['floored']['factor']['water'] > own['factor']['water']
    grid = product['records']['corrected']['lowest_level']['grid']
    assert 1.1 < grid['mean_speed_over_vector']['water'] < grid['cube_ratio']['water'] < F.wind_raise()


def test_the_floored_factor_and_the_quadrature_excess_stand_recorded(product):
    for key in R.RECORDS:
        f = product['records'][key]['lowest_level']['floored']
        assert f['floor_m_s'] == F.LEVEL_FLOOR_MS and len(f['by_latitude']) == len(R.BANDS)
        for kind in ('water', 'land'):
            assert 1.0 < f['factor'][kind] < 1.6 and 0.0 < f['quadrature_excess_m_s'][kind] < 4.0
            assert 0.0 <= f['writes_under_the_floor'][kind] < 0.6
            cells = f['factor'][f'{kind}_cells']
            assert cells['p05'] <= cells['median'] <= cells['p95']


def test_ratios_and_excess_on_made_up_winds():
    # a steady wind has no excess; a wind that turns through the window keeps its speed and loses its vector
    steady = np.full((4 * R.WRITES_PER_WINDOW, 2), 3.0, np.float32)
    speed, vector = R.windows(steady, np.zeros_like(steady))
    linear, cube = R.ratios(speed, vector)
    assert linear == pytest.approx(1.0) and cube == pytest.approx(1.0) and R.excess(speed, vector) == pytest.approx(0.0, abs=1e-6)
    angle = np.tile(np.linspace(-0.6, 0.6, R.WRITES_PER_WINDOW), 4)[:, None] * np.ones((1, 2))
    u, v = (3.0 * np.cos(angle)).astype(np.float32), (3.0 * np.sin(angle)).astype(np.float32)
    speed, vector = R.windows(u, v)
    linear, cube = R.ratios(speed, vector)
    expected = 3.0 / float(np.hypot(u[:R.WRITES_PER_WINDOW, 0].mean(), v[:R.WRITES_PER_WINDOW, 0].mean()))
    assert linear == pytest.approx(expected, rel=1e-5) and cube == pytest.approx(expected, rel=1e-5)
    assert R.excess(speed, vector) == pytest.approx(float(np.sqrt(9.0 - (3.0 / expected) ** 2)), rel=1e-4)
    # an added gust of fixed size raises a calm window's ratio more than a windy one's
    rng = np.random.default_rng(4)
    gust = rng.normal(0.0, 1.5, (400 * R.WRITES_PER_WINDOW, 2, 2)).astype(np.float32)
    mean = np.array([0.5, 5.0], np.float32)[None, None, :]
    speed, vector = R.windows(gust[:, 0, :] + mean[0], gust[:, 1, :])
    linear, _ = R.ratios(speed, vector)
    assert linear[0] > 1.5 > 1.1 > linear[1] > 1.0


def test_the_floor_and_the_quadrature_excess_on_made_up_winds():
    """A steady wind over the floor needs no raise and holds no excess; under the floor both sides stand at the floor.
    A gusty wind's factor is the ratio of the floored means, and its excess, added in quadrature to each window's
    floored vector speed, gives the windows their mean speed."""
    per = R.WRITES_PER_WINDOW
    steady = np.stack([np.full(6 * per, 3.0), np.full(6 * per, 0.4)], 1).astype(np.float32)
    speed, vector = R.windows(steady, np.zeros_like(steady))
    inside, outside, factor = R.floored(speed, vector)
    assert factor == pytest.approx([1.0, 1.0]) and outside[0].tolist() == pytest.approx([3.0, F.LEVEL_FLOOR_MS])
    assert R.quadrature_excess(inside, outside) == pytest.approx([0.0, 0.0], abs=1e-9)
    rng = np.random.default_rng(7)
    gust = rng.normal(0.0, 1.5, (300 * per, 2, 3)).astype(np.float32)
    mean = np.array([0.3, 2.0, 6.0], np.float32)[None, :]
    speed, vector = R.windows(gust[:, 0, :] + mean, gust[:, 1, :])
    inside, outside, factor = R.floored(speed, vector)
    assert factor == pytest.approx(inside.mean(0) / outside.mean(0)) and factor[0] > factor[1] > factor[2] > 1.0
    e = R.quadrature_excess(inside, outside)
    assert np.sqrt(outside ** 2 + e ** 2).mean(0) == pytest.approx(inside.mean(0), rel=1e-9)
    assert np.all(e > 0.5) and np.all(e < 3.0)
    # the raise by one factor leaves the calmest windows short of their speed, the excess does not
    calm = outside[:, 1] < np.percentile(outside[:, 1], 20)
    by_factor = (outside[:, 1] * factor[1])[calm].mean()
    by_excess = np.sqrt(outside[:, 1] ** 2 + e[1] ** 2)[calm].mean()
    assert by_factor < by_excess and abs(by_excess - inside[calm, 1].mean()) < abs(by_factor - inside[calm, 1].mean())
