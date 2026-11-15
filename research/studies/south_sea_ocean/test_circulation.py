import json

import numpy as np
import pytest

from climate.ocean import lunar, mitgcm
from climate.ocean import stratified as ST
from research.studies.south_sea_ocean import circulation as C


def test_a_cell_is_water_where_the_sea_covers_half_of_it():
    # pixels of 0.25 degrees on a 1 degree grid near the equator; each pixel's area as the grid counts it
    lat, lon = np.meshgrid(np.arange(0.125, 2.0, 0.25), np.arange(0.125, 3.0, 0.25), indexing='ij')
    r = lunar.RADIUS / 1000.0
    area = (r * np.radians(0.25)) ** 2 * np.cos(np.radians(lat))
    sea = (lon < 1.0) | ((lon < 1.5) & (lat < 1.0)) | ((lon > 2.5) & (lat > 1.5))       # full, half and a quarter
    depth = np.where(lon < 1.0, 100.0, 40.0)
    grid = C.floor_on_grid(lat[sea], lon[sea], depth[sea], area[sea], 1.0)
    j, i, inside = C.cells(np.array([0.5, 0.5, 1.5, 1.5]), np.array([0.5, 1.5, 1.5, 2.6]), grid.xg, grid.yg)
    assert inside.all()
    assert grid.depth[j, i] == pytest.approx([100.0, 40.0, 0.0, 0.0])
    assert grid.depth[0].sum() == 0 and grid.depth[:, 0].sum() == 0          # a land cell round the sea


def test_shallow_cells_count_as_land():
    lat, lon = np.meshgrid(np.arange(0.125, 1.0, 0.25), np.arange(0.125, 2.0, 0.25), indexing='ij')
    depth = np.where(lon < 1.0, 9.0, 11.0)
    grid = C.floor_on_grid(lat.ravel(), lon.ravel(), depth.ravel(), np.ones(lat.size) * 59.0, 1.0)
    assert sorted(grid.depth[grid.depth > 0].tolist()) == [11.0]


def test_the_record_is_carried_between_its_cells_and_round_the_seam():
    lat = np.array([60.0, 20.0, -20.0, -60.0])
    lon = np.arange(0.0, 360.0, 90.0)
    field = lat[:, None] + 0.0 * lon[None, :]
    index, weight = C.bilinear(lat, lon, np.array([0.0, -40.0, -80.0, 80.0]), np.array([10.0, 200.0, 359.0, 45.0]))
    assert weight.sum(axis=0) == pytest.approx(1.0)
    assert (field.ravel()[index] * weight).sum(axis=0) == pytest.approx([0.0, -40.0, -60.0, 60.0])   # held past the ends
    east = 0.0 * lat[:, None] + lon[None, :]
    index, weight = C.bilinear(lat, lon, np.array([20.0, 20.0]), np.array([45.0, 315.0]))
    assert (east.ravel()[index] * weight).sum(axis=0) == pytest.approx([45.0, 135.0])    # 315 lies between 270 and 0


def test_land_cells_take_their_nearest_water():
    lat = np.array([30.0, 0.0, -30.0])
    lon = np.arange(0.0, 360.0, 60.0)
    water = np.zeros((3, 6), bool)
    water[1, 0] = water[2, 3] = True
    fill = C.nearest_water(lat, lon, water)
    assert fill[1 * 6 + 0] == 6 and fill[2 * 6 + 3] == 15
    assert fill[0 * 6 + 0] == 6 and fill[0 * 6 + 3] == 15 and fill[1 * 6 + 5] == 6


def test_coriolis_follows_the_true_latitude():
    frame = ST.Frame(-50.0, 180.0)
    grid = mitgcm.Grid(xg=np.array([-1.0, 0.0, 1.0]), yg=np.array([-10.0, 0.0, 10.0]), depth=np.ones((2, 2)))
    f = C.coriolis(grid, frame)
    lat, _ = frame.to_true(*np.meshgrid(grid.yc, grid.xc, indexing='ij'))
    assert f['centre'] == pytest.approx(2.0 * lunar.OMEGA * np.sin(np.radians(lat)))
    assert f['corner'][1, 1] == pytest.approx(2.0 * lunar.OMEGA * np.sin(np.radians(-50.0)))
    assert (f['centre'][1] > f['centre'][0]).all() and (f['centre'] < 0).all()
    assert f['centre'] ** 2 + f['cosine'] ** 2 == pytest.approx((2.0 * lunar.OMEGA) ** 2)


def test_the_fluxes_keep_their_size_and_point_upward():
    lat = np.array([60.0, 20.0, -20.0, -60.0])
    lon = np.arange(0.0, 360.0, 90.0)
    shape = (2, 4, 4)
    record = dict(lat=lat, lon=lon, lsm=np.zeros((4, 4)), tauu=np.full(shape, 0.03), tauv=np.full(shape, 0.04),
                  rss=np.full(shape, 200.0), rls=np.full(shape, -50.0), hfss=np.full(shape, -10.0),
                  hfls=np.full(shape, -100.0), evap=np.full(shape, -3.0e-8), pr=np.full(shape, 1.0e-8),
                  ts=np.full(shape, 303.15))
    frame = ST.Frame(-40.0, 100.0)
    depth = np.ones((2, 3))
    depth[0, 0] = 0.0
    grid = mitgcm.Grid(xg=np.array([-20.0, -10.0, 0.0, 10.0]), yg=np.array([-10.0, 0.0, 10.0]), depth=depth)
    f = C.forcing_on_grid(grid, frame, record)
    wet = grid.wet
    assert np.hypot(f['ustress'], f['vstress'])[:, wet] == pytest.approx(0.05)
    assert f['hflux'][:, wet] == pytest.approx(-40.0)        # 40 W/m2 into the sea is a flux of -40 upward
    assert f['swflux'][:, wet] == pytest.approx(-200.0)
    assert f['sflux'][:, wet] == pytest.approx(2.0e-8)       # evaporation less rain, upward
    assert f['climsst'][:, wet] == pytest.approx(30.0)
    assert all((v[:, ~wet] == 0).all() for v in f.values())


def test_the_trial_is_recorded_with_its_setting():
    out = json.loads((C.HERE / 'results' / 'circulation_trial_1deg.json').read_text())
    assert out['schema'] == C.SCHEMA
    s = out['setting']
    assert abs(s['water_left_over_m3s']) < 1e-3
    assert s['rivers']['six_seas_m3s'] + s['rivers']['moscoviense_route_m3s'] == pytest.approx(s['forcing']['loss_m3s'])
    assert s['area_km2'] == pytest.approx(s['lola_pixels']['area_km2'], rel=0.03)
    assert s['volume_km3'] == pytest.approx(s['lola_pixels']['volume_km3'], rel=0.03)
    assert out['hours_a_model_year'] == pytest.approx(out['seconds_a_step'] * C.YEAR_S / s['deltat_s'] / 3600.0)
    m = out['monitor']
    assert m['time_secondsf'] == pytest.approx(out['model_days'] * 86400.0)
    assert 0.0 <= m['dynstat_salt_min'] <= m['dynstat_salt_max'] < 40.0
    assert 20.0 < m['dynstat_theta_min'] <= m['dynstat_theta_max'] < 40.0
