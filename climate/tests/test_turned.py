"""The GCM's surface on the turned Moon (turned.py): the rotation as domain/geology converts places, the turned surface's
form and kept water on synthetic grids (no geography files needed), the runner's frame, and the committed record when
present."""
import hashlib
import json
import unittest

import numpy as np
import pytest

from climate.gcm import boundary, turned, exoplasim_run as er


def synthetic(nlat=36, nlon=72, seed=3):
    """A lakes surface in boundary.lakes_product's form on a regular grid, north first, from 0 east."""
    rng = np.random.default_rng(seed)
    dlat, dlon = 180.0 / nlat, 360.0 / nlon
    lat = 90.0 - dlat * (np.arange(nlat) + 0.5)
    lon = dlon * (np.arange(nlon) + 0.5)
    sea = rng.random((nlat, nlon)) < 0.3
    wet = np.where(sea, 1.0, np.where(rng.random((nlat, nlon)) < 0.2, rng.random((nlat, nlon)), 0.0))
    height = rng.normal(0.0, 2000.0, (nlat, nlon))
    depth = np.where(wet > 0, rng.random((nlat, nlon)) * 900.0, 0.0)
    surface = np.where(sea, -1654.2, np.where(wet > 0, height + depth, -1654.2))
    prod = dict(lat_deg=lat, lon_deg=lon, water_fraction=wet, mean_height_m=height, mean_water_depth_m=depth,
                water_surface_m=surface, meta=dict(level_m=-1654.2), sha256='a+b', drainage_sha256='b',
                routing_run='synthetic')
    return prod, sea


def about_z(degrees):
    """The frame turned about the spin axis: its first axis, toward Earth, stands at today's longitude `degrees`."""
    a = np.radians(degrees)
    return np.array([[np.cos(a), np.sin(a), 0.0], [-np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])


class RotationTests(unittest.TestCase):
    def test_the_frames_conversions_invert_each_other_and_turn_about_the_axis(self):
        r = about_z(30.0)
        lat, lon = turned.to_turned(10.0, 50.0, r)
        self.assertAlmostEqual(float(lat), 10.0, places=10)
        self.assertAlmostEqual(float(lon), 20.0, places=10)
        back = turned.to_today(lat, lon, r)
        self.assertAlmostEqual(float(back[0]), 10.0, places=10)
        self.assertAlmostEqual(float(back[1]), 50.0, places=10)

    def test_the_identity_keeps_the_surface(self):
        prod, sea = synthetic()
        t = turned.turned_product(prod, sea, np.eye(3))
        np.testing.assert_allclose(t['water_fraction'], prod['water_fraction'], rtol=0, atol=1e-12)
        np.testing.assert_allclose(t['mean_height_m'], prod['mean_height_m'], rtol=1e-12, atol=1e-9)
        wet = prod['water_fraction'] > 0
        np.testing.assert_allclose(t['mean_water_depth_m'][wet], prod['mean_water_depth_m'][wet], rtol=1e-12)
        np.testing.assert_allclose(t['water_surface_m'][wet], prod['water_surface_m'][wet], rtol=1e-12)
        self.assertTrue(np.all(t['mean_water_depth_m'][~wet] == 0.0))

    def test_a_turn_about_the_axis_moves_every_cell_by_its_longitude(self):
        prod, sea = synthetic()
        t = turned.turned_product(prod, sea, about_z(90.0))           # a whole number of 5-degree cells
        np.testing.assert_allclose(t['water_fraction'], np.roll(prod['water_fraction'], -18, axis=1), atol=1e-12)
        np.testing.assert_allclose(t['mean_height_m'], np.roll(prod['mean_height_m'], -18, axis=1), rtol=1e-12,
                                   atol=1e-9)

    def test_a_tilted_frame_keeps_the_water_and_the_seas(self):
        prod, sea = synthetic(nlat=180, nlon=360)
        axis = turned.unit(0.0, -77.0)                                  # about the equator, as the Moon's turn is
        angle = np.radians(8.05)
        cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
        r = np.cos(angle) * np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * np.outer(axis, axis)
        t = turned.turned_product(prod, sea, r)
        w = np.cos(np.radians(prod['lat_deg']))[:, None] * np.ones(prod['lon_deg'].size)
        share = lambda f: float((f * w).sum() / w.sum())
        self.assertAlmostEqual(share(t['water_fraction']), share(prod['water_fraction']), delta=2e-3)
        self.assertAlmostEqual(t['sea_share'], share(sea.astype(float)), delta=2e-3)
        self.assertTrue(np.all((t['water_fraction'] >= 0) & (t['water_fraction'] <= 1)))
        self.assertTrue(np.all(t['water_surface_m'][t['water_fraction'] == 0] == prod['meta']['level_m']))


def quiet(monkeypatch):
    """The configuration without its one side effect, patching the shared PlaSim source, and without the installed
    package's version."""
    import importlib.metadata
    monkeypatch.setattr(er, 'ensure_patched', lambda *a, **k: dict(stub='source'))
    monkeypatch.setattr(importlib.metadata, 'version', lambda name: 'stub')


def test_without_a_frame_the_configuration_is_todays(monkeypatch):
    quiet(monkeypatch)
    cfg = er.configuration('A28R', 4)
    assert 'frame' not in cfg
    planet = json.loads((er.PRODUCTS / 'moon_gcm_configuration.json').read_text())['planet']
    assert cfg['planet']['obliquity'] == planet['obliquity_deg']


def test_a_record_of_another_surface_is_refused(tmp_path):
    path = tmp_path / 'turned_surface_x.json'
    path.write_text(json.dumps(dict(schema=turned.SCHEMA, surface='moon_28pct_water_lakes_0000_turned_x_gaussian_T21.nc',
                                    files={}, frame=dict(state='elastic', tilt_to_ecliptic_deg=1.6))))
    with pytest.raises(ValueError, match='another surface'):
        er.frame_setting(path, er.EXPERIMENTS['A28R'])
    path.write_text(json.dumps(dict(schema='other/1')))
    with pytest.raises(ValueError, match='record'):
        er.frame_setting(path, er.EXPERIMENTS['A28R'])


def test_the_built_frame_runs_the_turned_surface_and_tilt(monkeypatch):
    records = sorted(boundary.RESULTS.glob('turned_surface_*.json'))
    if not records:
        pytest.skip('the turned surface has not been built')
    record = json.loads(records[0].read_text())
    if not (er.PRODUCTS / record['surface']).is_file():
        pytest.skip('the turned surface file is not in products')
    quiet(monkeypatch)
    cfg = er.configuration('A28R', 4, frame=str(records[0]))
    assert cfg['frame']['surface'] == record['surface']
    assert cfg['planet']['obliquity'] == record['frame']['tilt_to_ecliptic_deg']
    assert cfg['planet']['obliquity'] != er.configuration('A28R', 4)['planet']['obliquity']
    assert {k: v for k, v in cfg.items() if k not in ('frame', 'planet')} == \
        {k: v for k, v in er.configuration('A28R', 4).items() if k not in ('frame', 'planet')}


class RecordTests(unittest.TestCase):
    def test_the_record_keeps_the_water_and_follows_geologys_places(self):
        records = sorted(boundary.RESULTS.glob('turned_surface_*.json'))
        if not records:
            self.skipTest('the turned surface has not been built')
        for path in records:
            record = json.loads(path.read_text())
            self.assertEqual(record['schema'], turned.SCHEMA)
            today, now = record['shares']['today'], record['shares']['turned']
            self.assertAlmostEqual(now['water'], today['water'], delta=1e-4)
            self.assertAlmostEqual(now['seas'], today['seas'], delta=1e-4)
            self.assertAlmostEqual(record['t21']['turned']['water_share'], now['water'], places=6)
            self.assertEqual(record['t21']['turned']['cells'], 32 * 64)
            for place in record['places']:
                dlat = place['turned'][0] - place['product_turned'][0]
                dlon = (place['turned'][1] - place['product_turned'][1] + 180.0) % 360.0 - 180.0
                self.assertLess(max(abs(dlat), abs(dlon)), 0.01, place['place'])
            r = np.asarray(record['frame']['rotation_today_to_turned'])
            np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
            for name, sha in record['files'].items():
                local = boundary.HERE / 'products' / name
                if local.exists():
                    self.assertEqual(hashlib.sha256(local.read_bytes()).hexdigest(), sha)


if __name__ == '__main__':
    unittest.main()
