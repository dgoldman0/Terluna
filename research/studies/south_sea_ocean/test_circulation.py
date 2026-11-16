import json
from pathlib import Path

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


def little_sea(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, 'RUNS', tmp_path)
    depth = np.array([[0.0, 30.0], [12.0, 30.0]])                 # land, two deep columns and one of 12 m
    grid = mitgcm.Grid(xg=np.array([0.0, 1.0, 2.0]), yg=np.array([-1.0, 0.0, 1.0]), depth=depth)
    flat = np.zeros((2, 2))
    return ST.Sea(name='little', grid=grid, deltat=1800.0, delr=np.array([10.0, 20.0]),
                  coriolis=dict(centre=flat, corner=flat, cosine=flat))


def test_the_dye_and_the_warmth_are_weighed_on_the_models_own_cells(tmp_path, monkeypatch):
    sea = little_sea(tmp_path, monkeypatch)
    volume = C.cell_volumes(sea)                                  # no run yet: the floor's own shares
    assert volume[0, 0, 0] == 0 and volume[1, 1, 0] == pytest.approx(2.0 * sea.grid.area[1, 0])
    # the model rounds the 2 m left in the second layer of the 12 m column to rock; its own shares say so
    sea.folder.mkdir(parents=True)
    share = np.array([[[0.0, 1.0], [1.0, 1.0]], [[0.0, 1.0], [0.0, 1.0]]])
    share.astype('>f4').tofile(sea.folder / 'hFacC.data')
    (sea.folder / 'hFacC.meta').write_text(" nDims = [   3 ];\n dimList = [\n     2,    1,    2,\n     2,    1,    2,\n"
                                           "     2,    1,    2\n ];\n dataprec = [ 'float32' ];\n nrecords = [     1 ];\n")
    volume = C.cell_volumes(sea)
    assert volume[1, 1, 0] == 0.0
    dye = np.where(share > 0, 1.0, 0.0)                           # rock reads as zero in the model's output
    dye[0] *= 0.5
    old = C.old_water(sea, dye, volume)
    assert old['below 0 m'] == pytest.approx((0.5 * volume[0].sum() + volume[1].sum()) / volume.sum())
    assert C.by_depth(dye, volume) == pytest.approx([0.5, 1.0])
    theta = np.where(share > 0, 27.0, 0.0)
    theta[0] = np.where(share[0] > 0, [[0.0, 30.0], [29.0, 31.0]], 0.0)
    true_lat = np.array([[-20.0, -20.0], [-50.0, -50.0]])
    w = C.warmth(sea, theta, volume, true_lat)
    assert w['by_depth_c'][1] == pytest.approx(27.0)              # the rock cell's zero is left out
    assert w['bands']['14-30S']['top_layer_c']['mean'] == pytest.approx(30.0)
    assert w['bands']['45-60S']['top_layer_c'] == pytest.approx(dict(mean=30.0, coldest=29.0, warmest=31.0), rel=0.01)
    assert w['bands']['14-30S']['at_or_under_the_limit_from_m']['median'] == pytest.approx(20.0)
    assert w['bands']['45-60S']['at_or_under_the_limit_from_m']['share_of_area_with_any'] == pytest.approx(0.5, abs=0.01)


def test_a_patch_holds_fine_cells_inside_the_coarse_grids_span():
    xg, yg = np.arange(-33.0, 33.1, 1.0), np.arange(-37.0, 35.1, 1.0)
    x, y = C.patch_edges(xg, yg, (-19.0, -9.0, 24.0, 34.0), 0.1)
    assert x[0] == xg[0] and x[-1] == xg[-1] and y[0] == yg[0] and y[-1] == yg[-1]
    assert (np.diff(x) > 0).all() and (np.diff(y) > 0).all()
    inside = (x[:-1] >= -19.0 - 1e-9) & (x[1:] <= -9.0 + 1e-9)
    assert np.diff(x)[inside] == pytest.approx(0.1) and inside.sum() == 100
    assert np.diff(x).max() <= 1.0 + 1e-9
    grow = np.diff(x)[np.flatnonzero(inside)[-1] + 1:][:5]
    assert grow[1:] / grow[:-1] == pytest.approx(1.15, rel=1e-6)          # widening by the ratio away from the box


def test_a_start_is_carried_over_from_the_nearest_water_at_each_level():
    parent = mitgcm.Grid(xg=np.array([0.0, 1.0, 2.0]), yg=np.array([0.0, 1.0]), depth=np.array([[30.0, 10.0]]))
    share = np.array([[[1.0, 1.0]], [[1.0, 0.0]], [[0.0, 0.0]]])        # the east column holds one level, none the third
    fields = dict(theta=np.array([[[28.0, 30.0]], [[27.0, 0.0]], [[0.0, 0.0]]]))
    child = mitgcm.Grid(xg=np.linspace(0.0, 2.0, 5), yg=np.array([0.0, 0.5, 1.0]), depth=np.ones((2, 4)))
    out = C.carry_over(parent, share, fields, child)['theta']
    assert out[0, 0].tolist() == [28.0, 28.0, 30.0, 30.0]
    assert out[1, 0].tolist() == [27.0, 27.0, 27.0, 27.0]               # only the west column holds the second level
    assert (out[2] == out[1]).all()                                     # below the parent's water the level above serves


def test_the_way_over_water_goes_round_the_land():
    # a U of water round a tongue of land, at the equator: cells of 0.1 degrees
    depth = np.full((3, 5), 50.0)
    depth[0:2, 2] = 0.0
    grid = mitgcm.Grid(xg=np.linspace(0.0, 0.5, 6), yg=np.linspace(-0.15, 0.15, 4), depth=depth)
    d = C.water_distance_km(grid, (0, 0))
    cell = np.radians(0.1) * lunar.RADIUS / 1000.0
    assert d[0, 0] == 0.0 and np.isinf(d[0, 2])
    assert d[0, 1] == pytest.approx(cell, rel=1e-3)
    # to the far foot of the U: four steps corner to corner, over the tongue's head
    assert d[0, 4] == pytest.approx(4.0 * np.sqrt(2.0) * cell, rel=1e-3)
    # two cells that meet at a corner with land on both sides of it are not joined
    corner = mitgcm.Grid(xg=np.array([0.0, 0.1, 0.2]), yg=np.array([0.0, 0.1, 0.2]), depth=np.array([[50.0, 0.0], [0.0, 50.0]]))
    assert np.isinf(C.water_distance_km(corner, (0, 0))[1, 1])


def test_the_reach_is_read_by_distance_from_the_mouth(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, 'RUNS', tmp_path)
    grid = mitgcm.Grid(xg=np.linspace(0.0, 0.4, 5), yg=np.array([0.0, 0.1]), depth=np.full((1, 4), 30.0))
    sea = ST.Sea(name='arm', grid=grid, deltat=900.0, delr=np.array([10.0, 20.0]))
    volume = C.cell_volumes(sea)
    cell = np.radians(0.1) * lunar.RADIUS / 1000.0
    distance = C.water_distance_km(grid, (0, 0))
    salt = np.array([[[0.0, 10.0, 20.0, 35.0]], [[10.0, 35.0, 35.0, 35.0]]]) * C.R.practical_from_reference(
        C.R.sc.SEAWATER_G_KG) / C.R.sc.SEAWATER_G_KG
    out = C.reach(sea, salt, volume, distance, grid.wet, band_km=cell * 1.9)      # two cells to a band
    near, far = out['bands']
    assert near['area_km2'] == pytest.approx(grid.area[0, :2].sum() / 1e6) and near['depth_m'] == pytest.approx(30.0)
    assert near['top_layer_g_kg'] == pytest.approx(dict(mean=5.0, freshest=0.0, saltiest=10.0), abs=1e-6)
    assert near['floor_g_kg'] == pytest.approx(22.5) and far['floor_g_kg'] == pytest.approx(35.0)
    assert near['fresher_water_m'] == pytest.approx({'under 30': 20.0, 'under 18': 20.0, 'under 5': 5.0})
    assert far['fresher_water_m'] == pytest.approx({'under 30': 5.0, 'under 18': 0.0, 'under 5': 0.0})
    reaches = out['top_layer_reaches_km']
    assert reaches['under 30'] == pytest.approx(2.0 * cell, rel=1e-3) and reaches['under 18'] == pytest.approx(cell, rel=1e-3)
    assert reaches['under 5'] == 0.0 and reaches['under 0.5'] == 0.0


def test_the_tide_is_the_equilibrium_height_at_the_turned_places():
    from climate.ocean import forcing as F
    frame = ST.Frame(-47.8, -176.9)
    grid = mitgcm.Grid(xg=np.array([-1.0, 0.0, 1.0]), yg=np.array([-1.0, 0.0, 1.0]), depth=np.ones((2, 2)))
    h = C.tide_on_grid(grid, frame, 2.0)
    assert h.shape == (2, 2, 2)
    lat, lon = frame.to_true(np.array([grid.yc[1]]), np.array([grid.xc[0]]))
    want = F.ephemeris_series(np.array([86400.0]))[0] @ F.monomials(lat, lon)[0]
    assert h[1, 1, 0] == pytest.approx(want, abs=1e-9)
    assert 0.01 < np.abs(h).max() < 5.0                                   # metres of the Moon's monthly tide


def test_flows_on_faces_are_carried_to_the_centres():
    u = np.array([[[1.0, 3.0, 5.0]]])
    v = np.array([[[2.0], [4.0]]])
    uc, vc = C.at_centres(u, v)
    assert uc.tolist() == [[[2.0, 4.0, 5.0]]] and vc.tolist() == [[[3.0], [4.0]]]


def write_mds(stem, fields):
    shape = next(iter(fields.values())).shape
    dims = ",\n".join(f"{n:6d},{1:6d},{n:6d}" for n in shape[::-1])
    names = " ".join(f"'{k:<8s}'" for k in fields)
    Path(f"{stem}.meta").write_text(f" nDims = [ {len(shape)} ];\n dimList = [\n{dims}\n ];\n dataprec = [ 'float32' ];\n"
                                    f" nrecords = [ {len(fields)} ];\n fldList = {{\n {names}\n }};\n")
    np.concatenate([np.ravel(v) for v in fields.values()]).astype('>f4').tofile(f"{stem}.data")   # big-endian, as MITgcm writes


def test_the_tides_part_is_the_difference_of_two_runs(tmp_path, monkeypatch):
    monkeypatch.setattr(mitgcm, 'RUNS', tmp_path)
    grid = mitgcm.Grid(xg=np.array([0.0, 1.0, 2.0]), yg=np.array([-1.0, 0.0, 1.0]), depth=np.full((2, 2), 30.0))
    pattern = np.array([[-1.0, 1.0], [-1.0, 1.0]])                    # a tilt over two cells of equal area
    tide = np.array([0.0, 1.0, 0.0, -1.0, 0.0])[:, None, None] * pattern
    seas = {}
    for on in (False, True):
        sea = ST.Sea(name='on' if on else 'off', grid=grid, deltat=86400.0, delr=np.array([10.0, 20.0]),
                     tide_m=tide if on else 0.0 * tide, tide_period_s=86400.0)
        sea.folder.mkdir(parents=True)
        for day in range(1, 5):
            level = 0.5 * (tide[day - 1] + tide[day]) if on else 0.0 * pattern       # the model follows equilibrium
            flow = np.full((2, 2, 2), 0.1 if on else 0.0)
            write_mds(sea.folder / f'deep.{day:010d}', dict(UVEL=flow, VVEL=0.0 * flow))
            write_mds(sea.folder / f'surface.{day:010d}', dict(ETAN=level))
        seas[on] = sea
    out = C.tide_part(seas[False], seas[True], skip_days=0)
    assert out['days'] == [1.0, 4.0]
    assert out['level_range_m']['ratio'] == pytest.approx([1.0, 1.0, 1.0])
    assert out['level_range_m']['model'] == pytest.approx([1.0, 1.0, 1.0])
    assert out['level_range_m']['correlation'] == pytest.approx(1.0)
    assert out['depth_mean_tidal_flow_m_s']['largest_by_cell'] == pytest.approx([0.1] * 4, rel=1e-6)
    assert out['floor_tidal_flow_m_s']['largest_by_cell'][-1] == pytest.approx(0.1, rel=1e-6)
    power = seas[True].rho_const * seas[True].drag * 0.1 ** 3 * grid.area.sum() / 1e6
    assert out['bed_drag_power']['whole_sea_mw'] == pytest.approx(power, rel=1e-5)
    assert out['bed_drag_power']['under_200_m_of_water_mw'] == pytest.approx(power, rel=1e-5)


def test_the_background_mixing_spends_its_power_on_the_stratification(tmp_path, monkeypatch):
    sea = little_sea(tmp_path, monkeypatch)
    theta = np.array([[[0.0, 30.0], [30.0, 30.0]], [[0.0, 27.0], [27.0, 27.0]]])      # 3 K warmer at the top
    salt = np.full((2, 2, 2), 35.0)
    wet = sea.grid.wet
    drho = C.R.seawater_density(35.0, 27.0) - C.R.seawater_density(35.0, 30.0)
    want = sea.diff_r * lunar.GRAVITY * drho * sea.grid.area[wet].sum() / 1e6
    assert C.stirring_power_mw(sea, dict(THETA=theta, SALT=salt)) == pytest.approx(want)
    assert want > 0


def test_the_sea_is_the_lake_models_basin_and_keeps_its_demand(tmp_path, monkeypatch):
    sea = little_sea(tmp_path, monkeypatch)
    volume = C.cell_volumes(sea)
    basin = C.hypsometry(sea, volume)
    assert basin.z.tolist() == [0.0, 10.0, 30.0]
    assert basin.area[0] == pytest.approx(sea.grid.area[sea.grid.wet].sum()) and basin.area[-1] == 0.0
    s = C.central_scenarios()['oligotrophic']
    lake = C.demand_lake_form(basin, s)
    assert (lake['water'] >= 0).all() and (lake['floor'] >= 0).all() and lake['floor'][-1] > 0   # the floor lies below 10 m
    ocean = C.demand_ocean_form(basin, s)
    kept = lambda d, b=None: float(((d['water'] + d['floor']) * (b or basin).volume).sum())
    assert kept(ocean) == pytest.approx(kept(lake), rel=1e-9)               # the depth changes, the total does not
    assert ocean['export_g_m2_yr'] == 0.0 and ocean['from_m'] == 30.0       # a sea of 30 m holds no water below 100 m
    # on a deep basin the flux's water demand falls with depth below 100 m, and the total is kept
    from biosphere.ecology import lake as L
    deep = L.Basin(np.linspace(0.0, 1000.0, 101), np.linspace(1.0, 0.9, 101) * 1e12)
    ocean = C.demand_ocean_form(deep, s)
    lower = ocean['water'][np.flatnonzero(deep.z[:-1] >= 100.0)]
    assert (np.diff(lower) < 0).all() and lower[0] > 5.0 * lower[-1]
    assert ocean['export_g_m2_yr'] > 0 and ocean['from_m'] == 100.0
    assert kept(ocean, deep) == pytest.approx(kept(C.demand_lake_form(deep, s), deep), rel=1e-9)


def test_the_waters_age_is_the_area_under_its_dye():
    # a level whose dye falls by 1% a year holds water 100 years old on average; one cleared at once, none
    years = [dict(year=float(y), old_water_by_depth=[0.0, float(np.exp(-0.01 * (y - 0.5) * C.YEAR_S / C.R.YEAR_S))])
             for y in range(1, 41)]
    aged = C.water_age(years, 2)
    assert aged['renewed'].tolist() == [True, False] and aged['age'][0] == 0.0
    assert aged['rate'][1] == pytest.approx(0.01)
    assert aged['age'][1] == pytest.approx(100.0, rel=2e-3)


def test_oxygen_is_saturation_less_the_demand_over_the_waters_age(tmp_path, monkeypatch):
    sea = little_sea(tmp_path, monkeypatch)
    volume = C.cell_volumes(sea)
    band = dict(top_layer_c=dict(mean=28.0, coldest=27.5, warmest=28.5))
    years = [dict(year=float(y), old_water_by_depth=[0.0, float(np.exp(-0.01 * (y - 0.5) * C.YEAR_S / C.R.YEAR_S))],
                  salinity_by_depth_g_kg=[35.0, 35.1], warmth=dict(bands={'60-81S': band})) for y in range(1, 12)]
    out = C.oxygen_estimate(sea, years, volume, [5.0, 20.0])
    assert out['renewed_within_the_first_year'] == [True, False]
    assert out['age_years'][0] == 0.0 and out['age_years'][1] == pytest.approx(100.0, rel=2e-3)
    sat = out['saturation_mg_l']
    assert 6.0 < sat < 6.8 and out['at_c'] == 28.0                          # sea water at 28 C under the Moon's air
    assert out['left_to_the_air'] == [True, False]
    for forms in out['states'].values():
        for f in forms.values():
            d = np.array(f['demand_g_m3_yr']['water']) + np.array(f['demand_g_m3_yr']['floor'])
            assert f['oxygen_mg_l'][0] is None
            assert f['oxygen_mg_l'][1] == pytest.approx(max(0.0, sat - d[1] * out['age_years'][1]))


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
