"""Checks of the atlas round's surface for the design climate's run: the 16 to 4 px/deg averaging keeps the dry ground's
height exactly; a frame made from a rotation alone places the pole and the point under Earth; the dry ground's height
averaged as it stands agrees with boundary.regrid wherever nothing lies below the water level; the turned surface's
files and record take the form the runner reads, and the runner's experiment A28S refuses another surface; and the
product once built: the settled frame under the final water, the water kept through the turn, the files as recorded and
the runner taking the record."""
import hashlib
import json

import numpy as np
import pytest

from climate.gcm import boundary as bd
from climate.gcm import exoplasim_run as er
from climate.gcm import turned
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import settled_turn as ST
from research.studies.atlas_round import turned_surface as TS


def synthetic(nlat=36, nlon=72, seed=3):
    """A lakes surface in boundary.lakes_product's form on a regular grid, north first, from 0 east, with dry ground
    above and below the water level."""
    rng = np.random.default_rng(seed)
    dlat, dlon = 180.0 / nlat, 360.0 / nlon
    lat = 90.0 - dlat * (np.arange(nlat) + 0.5)
    lon = dlon * (np.arange(nlon) + 0.5)
    sea = rng.random((nlat, nlon)) < 0.3
    wet = np.where(sea, 1.0, np.where(rng.random((nlat, nlon)) < 0.2, rng.random((nlat, nlon)), 0.0))
    height = rng.normal(0.0, 2000.0, (nlat, nlon))
    depth = np.where(wet > 0, rng.random((nlat, nlon)) * 900.0, 0.0)
    surface = np.where(sea, -1654.2, np.where(wet > 0, height + depth, -1654.2))
    w = np.cos(np.radians(lat))[:, None] * np.ones(nlon)
    share = lambda a: float((a * w).sum() / w.sum())
    prod = dict(lat_deg=lat, lon_deg=lon, water_fraction=wet, mean_height_m=height, mean_water_depth_m=depth,
                water_surface_m=surface, meta=dict(level_m=-1654.2), sha256='a+b', drainage_sha256='b',
                routing_run='synthetic', sea_share=share(sea.astype(float)), lake_share=share(np.where(sea, 0.0, wet)))
    return prod, sea


def about(axis_lat, axis_lon, degrees):
    axis = turned.unit(axis_lat, axis_lon)
    a = np.radians(degrees)
    cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    return np.cos(a) * np.eye(3) + np.sin(a) * cross + (1.0 - np.cos(a)) * np.outer(axis, axis)


def test_a_frame_from_a_rotation_alone_places_the_pole_and_the_point_under_earth():
    r = about(0.0, -60.8, 6.97)
    frame = TS.frame_from_rotation(r, 1.707, 'settled', 'ab' * 32)
    pole, earth = frame['frame']['north_pole'], frame['frame']['sub_earth_point']
    assert float(turned.to_turned(pole['lat_deg'], pole['lon_deg'], r)[0]) == pytest.approx(90.0, abs=1e-6)
    lat, lon = turned.to_turned(earth['lat_deg'], earth['lon_deg'], r)
    assert (float(lat), float(lon)) == pytest.approx((0.0, 0.0), abs=1e-9)
    assert frame['frame']['pole_move_deg'] == pytest.approx(6.97, abs=1e-3)
    assert (frame['state'], frame['tilt_deg'], frame['sha256']) == ('settled', 1.707, 'ab' * 32)
    with pytest.raises(ValueError, match='proper rotation'):
        TS.frame_from_rotation(-r, 1.707, 'settled', 'ab' * 32)


def test_the_dry_grounds_height_is_averaged_as_it_stands():
    prod, _ = synthetic()
    plain, signed = bd.regrid(prod, 'gaussian_T21'), TS.regrid_signed(prod, 'gaussian_T21')
    assert (plain['land_elevation_m'] >= 0).all() and (signed['land_elevation_m'] < 0).any()
    assert (signed['land_elevation_m'] <= plain['land_elevation_m'] + 1e-6).all()
    for name in ('land_fraction', 'ocean_depth_m', 'cell_area_fraction'):
        np.testing.assert_array_equal(signed[name], plain[name])
    np.testing.assert_allclose(signed['water_surface_m'], plain['water_surface_m'], atol=1e-6)
    low = dict(prod, meta=dict(level_m=prod['meta']['level_m'] - 10000.0))      # nothing below the water level
    a, b = bd.regrid(low, 'gaussian_T21'), TS.regrid_signed(low, 'gaussian_T21')
    for name in ('land_elevation_m', 'water_surface_m'):
        np.testing.assert_allclose(b[name], a[name], atol=1e-6)
    with pytest.raises(ValueError, match='below the water level'):
        TS.regrid_signed(prod, 'gaussian_T21', below=100.0)


def test_the_turned_surface_is_written_in_the_form_the_runner_reads(tmp_path):
    prod, sea = synthetic()
    turn = np.array([[np.cos(np.radians(30.0)), np.sin(np.radians(30.0)), 0.0],
                     [-np.sin(np.radians(30.0)), np.cos(np.radians(30.0)), 0.0], [0.0, 0.0, 1.0]])
    frame = TS.frame_from_rotation(turn, 1.7, 'settled', 'cd' * 32)
    describe = dict(producer='test', source_product='synthetic', file_evidence='a test', evidence='a test')
    path = tmp_path / 'record.json'
    record = TS.write_surface(prod, sea, frame, source='synthetic', out=tmp_path, inputs={}, path=path,
                              producer=dict(domain='test', files={}), describe=describe)
    assert json.loads(path.read_text()) == record and record['schema'] == turned.SCHEMA
    assert record['surface'] == turned.turned_name('b', 'cd' * 8, 'gaussian_T21')
    assert sorted(record['files']) == sorted(turned.turned_name('b', 'cd' * 8, grid) for grid in bd.GRIDS)
    for name, sha in record['files'].items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == sha
    assert record['frame']['state'] == 'settled' and record['frame']['tilt_to_ecliptic_deg'] == 1.7
    assert record['shares']['turned']['water'] == pytest.approx(record['shares']['today']['water'], abs=1e-12)
    assert record['t21']['turned']['water_share'] == pytest.approx(record['t21']['today']['water_share'], abs=2e-3)
    t21 = bd.read_t21(tmp_path / record['surface'])
    assert t21['water_share'] == pytest.approx(record['t21']['turned']['water_share'])


def test_the_rounds_experiment_takes_only_its_own_surface(tmp_path):
    lakes = er.EXPERIMENTS['A28S']['lakes']
    path = tmp_path / 'turned_surface_x.json'
    for surface, match in ((f'moon_28pct_water_lakes_{lakes}_turned_x_gaussian_T21.nc', 'missing'),
                           (f"moon_28pct_water_lakes_{er.EXPERIMENTS['A28R']['lakes']}_turned_x_gaussian_T21.nc",
                            'another surface')):
        path.write_text(json.dumps(dict(schema=turned.SCHEMA, surface=surface, files={surface: '0' * 64},
                                        frame=dict(state='settled', tilt_to_ecliptic_deg=1.7))))
        with pytest.raises(ValueError, match=match):
            er.frame_setting(path, er.EXPERIMENTS['A28S'])


def test_the_coarse_cells_keep_the_dry_grounds_height_exactly():
    rng = np.random.default_rng(5)
    h = rng.normal(0.0, 500.0, (8, 16))
    sea = np.full(h.shape, -1, np.int8)
    sea[:4, :4] = 0                                               # one coarse cell all sea, standing at 100 m
    h[:4, :4] = rng.uniform(-900.0, 50.0, (4, 4))
    sea[4:6, 8:12] = 1                                            # half of another, standing at -200 m
    h[4:6, 8:12] = rng.uniform(-900.0, -250.0, (2, 4))
    lake_cell = np.array([6 * 16 + 1, 7 * 16 + 2, 5 * 16 + 13])     # three lake cells on land
    lake_level = h.ravel()[lake_cell] + 30.0
    c = TS.coarse(h, sea, lake_cell, lake_level, [100.0, -200.0])
    assert c['water_fraction'][0, 0] == 1.0 and c['sea_fraction'][0, 0] == 1.0
    assert c['mean_water_depth_m'][0, 0] == pytest.approx(np.mean(100.0 - h[:4, :4]))
    assert c['water_surface_m'][0, 0] == pytest.approx(100.0)
    assert c['water_fraction'][1, 2] == 0.5 and c['water_surface_m'][1, 2] == pytest.approx(-200.0)
    assert c['water_fraction'][1, 0] == 2 / 16 and c['mean_water_depth_m'][1, 0] == pytest.approx(30.0)
    assert c['lakes_on_seas'] == 0 and c['water_at_or_below_ground'] == 0
    wet = np.zeros(h.size, bool)
    wet[(sea.ravel() >= 0)] = True
    wet[lake_cell] = True
    wet = wet.reshape(h.shape)
    for i in range(2):
        for j in range(4):
            block, dry = h[4 * i:4 * i + 4, 4 * j:4 * j + 4], ~wet[4 * i:4 * i + 4, 4 * j:4 * j + 4]
            w = c['water_fraction'][i, j]
            assert c['mean_height_m'][i, j] == pytest.approx(block.mean())
            if dry.any():                                         # the dry part's height as boundary.regrid derives it
                bed = c['water_surface_m'][i, j] - c['mean_water_depth_m'][i, j] if w > 0 else 0.0
                assert (c['mean_height_m'][i, j] - w * bed) / (1.0 - w) == pytest.approx(block[dry].mean())
    assert np.isnan(c['water_surface_m'][0, 1])                   # a dry cell has no water surface


@pytest.fixture(scope='module')
def record():
    if not ST.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.settled_turn first')
    path = TS.record_path(hashlib.sha256(ST.OUT.read_bytes()).hexdigest())
    if not path.exists():
        pytest.skip('run python -m research.studies.atlas_round.turned_surface first')
    return path, json.loads(path.read_text())


def test_the_record_turns_the_final_water_into_the_settled_frame(record):
    _, r = record
    st = json.loads(ST.OUT.read_text())
    settled = st['cases']['final_water']['settled']
    np.testing.assert_array_equal(r['frame']['rotation_today_to_turned'], st['frame_today_to_turned'])
    assert r['frame']['tilt_to_ecliptic_deg'] == settled['tilt_to_ecliptic_deg']
    assert r['frame']['pole_move_deg'] == pytest.approx(settled['pole_deg'], abs=1e-3)
    assert r['frame']['state'] == TS.STATE
    today, now = r['shares']['today'], r['shares']['turned']
    assert now['water'] == pytest.approx(today['water'], abs=1e-3)
    assert now['seas'] == pytest.approx(today['seas'], abs=1e-3)
    assert r['water']['lake_cells_on_seas'] == 0
    assert r['t21']['turned']['water_share'] == pytest.approx(now['water'], abs=1e-3)
    for name, sha in r['files'].items():
        local = TS.PRODUCTS / name
        if local.exists():
            assert hashlib.sha256(local.read_bytes()).hexdigest() == sha


def test_the_runner_takes_the_record_for_the_rounds_experiment(record):
    from climate.gcm import exoplasim_run as er
    path, r = record
    assert r['ignored_inputs'][str(GD.GOVERNED_ROUTING)][:16] == er.EXPERIMENTS['A28S']['lakes']
    if not (er.PRODUCTS / r['surface']).is_file():
        pytest.skip('the surface file is not in products')
    frame = er.frame_setting(path, er.EXPERIMENTS['A28S'])
    assert (frame['surface'], frame['tilt_deg'], frame['state']) == \
        (r['surface'], r['frame']['tilt_to_ecliptic_deg'], TS.STATE)
