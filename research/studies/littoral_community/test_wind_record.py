"""Tests for the wind inside a 3-day mean: the record's binding, the ratios on made-up winds, and what the product
shows of the raise the forcing takes."""
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
    for path in R.record_files():
        if not path.exists():
            pytest.skip('the 6.4-hour record is not on this machine')
        assert F.sha256(path) == product['inputs'][str(path.relative_to(F.RUNS))], path


def test_the_forcings_raise_is_the_waters_oxygen_studys_and_sits_in_the_grids_range(product):
    factors = F.wind_factors()
    assert len(factors) == 6 and all(1.05 < x < 1.4 for x in factors)
    assert F.wind_raise() == pytest.approx(float(np.mean(factors)), rel=1e-12)
    assert product['the_forcings_raise'] == dict(factor=pytest.approx(F.wind_raise(), rel=1e-12), cells=factors)
    cube = product['grid']['cube_ratio']
    assert cube['cells']['p05'] < F.wind_raise() < cube['cells']['p95']
    assert abs(F.wind_raise() - cube['area_mean']) < 0.05
    # the cube's ratio stands at or above the mean speed's over nearly all the grid, and both above one
    linear = product['grid']['mean_speed_over_vector']
    assert 1.0 < linear['cells']['lowest'] and 1.0 < cube['cells']['lowest']
    assert linear['area_mean'] < cube['area_mean'] < 1.3
    for row in product['by_latitude']:
        for kind in ('water', 'land'):
            assert 1.1 < row[kind]['mean_speed_over_vector'] < 1.25 and 1.1 < row[kind]['cube_ratio'] < 1.25, row


def test_the_record_holds_a_like_excess_in_calm_and_windy_windows(product):
    """The windows' speed beyond their mean vector is nearly one figure at every vector speed, so the ratio falls
    from several in calm windows toward one in windy ones: one factor raises calm windows too little."""
    rows = [r for r in product['by_vector_speed'] if r['share_of_windows'] > 0.01]
    assert len(rows) >= 5 and sum(r['share_of_windows'] for r in product['by_vector_speed']) == pytest.approx(1.0, abs=1e-3)
    excess = [r['excess_rms_m_s'] for r in rows]
    assert max(excess) - min(excess) < 0.1 * np.mean(excess)
    ratios = [r['ratio'] for r in rows]
    assert ratios == sorted(ratios, reverse=True) and ratios[0] > 2.0 > 1.1 > ratios[-1]
    assert rows[0]['ratio'] > F.wind_raise() > rows[-1]['ratio']


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
