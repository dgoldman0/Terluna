import json

import numpy as np
import pytest

from climate.gcm import global_winds as gw


def test_weighted_percentiles_of_a_uniform_histogram():
    edges = np.linspace(0.0, 10.0, 101)
    counts = np.ones((2, 100))
    counts[1] = 0.0                                   # an empty row gives no percentile
    p = gw.weighted_percentiles(counts, edges, (10, 50, 99))
    assert p[:, 0] == pytest.approx([1.0, 5.0, 9.9])
    assert np.all(np.isnan(p[:, 1]))


def test_to_heights_interpolates_within_the_column_and_holds_below_it():
    z = np.array([[900.0, 4000.0, 10000.0]])
    speed = np.array([[2.0, 5.0, 11.0]])
    rho = np.array([[1.3, 1.1, 0.9]])
    h = np.array([0.0, 2450.0, 7000.0, 12000.0])
    s = gw.to_heights(z, speed, h)[:, 0]
    assert s[:3] == pytest.approx([2.0, 3.5, 8.0])
    assert np.isnan(s[3])
    r = gw.to_heights(z, rho, np.array([7000.0]), log=True)[0, 0]
    assert r == pytest.approx(np.sqrt(1.1 * 0.9))


def test_stored_design_run_product():
    path = gw.RESULTS / 'global_winds_A28_dim5.json'
    if not path.is_file():
        pytest.skip('the global wind product has not been generated')
    d = json.loads(path.read_text())
    assert d['schema'] == gw.SCHEMA
    z = np.array(d['height_km'])
    rho = np.array(d['density_mean_kg_m3'], dtype=float)
    assert np.all(np.diff(rho[z >= 5.0]) < 0)
    band = (z >= 35.0) & (z <= 45.0)
    assert 0.6 < rho[band].min() and rho[band].max() < 0.85
    for key in ('speed_p50_m_s', 'speed_p99_m_s', 'speed_max_m_s'):
        v = np.array(d[key], dtype=float)
        assert np.all(v[z >= 5.0] > 0) and np.all(v < 40.0)
    assert all(a <= b for a, b in zip(d['speed_p99_m_s'], d['speed_max_m_s']))
