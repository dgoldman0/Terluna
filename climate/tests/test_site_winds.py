import json

import numpy as np
import pytest

from climate.gcm import site_winds as sw


def test_bilinear_reproduces_a_linear_field_and_wraps_in_longitude():
    lat = np.array([8.31, 2.77, -2.77])            # north first, as in the model
    lon = np.arange(64) * 5.625
    field = 2.0 * lat[:, None] + 0.0 * lon[None, :] + 3.0
    w = sw.bilinear(lat, lon, 5.375, 201.375)
    assert sw.at_site(field, w) == pytest.approx(2.0 * 5.375 + 3.0)
    ring = np.cos(np.radians(lon))[None, :].repeat(3, axis=0)
    w = sw.bilinear(lat, lon, 0.0, 358.0)          # between the last column and the first
    assert sum(w[3]) == pytest.approx(1.0) and set(w[1]) == {63, 0}
    assert np.cos(np.radians(354.375)) <= sw.at_site(ring, w) <= 1.0


def test_layer_heights_of_an_isothermal_column():
    R, gm, radius, t = 288.0, 4.9e12, 1.7374e6, 280.0
    sigma = np.array([0.9, 0.6, 0.3])
    ps = np.array([1.2e5])
    z, rho = sw.layer_heights(sigma, ps, np.full((1, 3), t), np.zeros((1, 3)), R, gm, radius)
    phi = R * t * np.log(1.0 / sigma)
    assert z[0] == pytest.approx(gm / (gm / radius - phi) - radius)
    assert rho[0] == pytest.approx(sigma * ps[0] / (R * t))


def test_gumbel_return_level_rises_with_the_period():
    maxima = [5.1, 5.6, 4.9, 6.2, 5.4, 5.8, 5.0, 6.0, 5.3, 5.7]
    assert sw.gumbel_return(maxima, 10.0) < sw.gumbel_return(maxima, 50.0) < 9.0


def test_stored_summit_product():
    path = sw.RESULTS / 'site_winds_A28_dim5_summit.json'
    if not path.is_file():
        pytest.skip('the summit wind product has not been generated')
    d = json.loads(path.read_text())
    assert d['schema'] == sw.SCHEMA
    assert np.all(np.diff(d['height_m']) > 0) and np.all(np.diff(d['density_kg_m3']) < 0)
    assert 5000.0 < d['ground_height_m'] < 9000.0
    for key in ('speed_p50_m_s', 'speed_p99_m_s', 'speed_max_m_s'):
        assert all(v > 0 for v in d[key])
    assert all(a <= b for a, b in zip(d['mean_cube_u_m3_s3'], d['mean_cube_m3_s3']))
