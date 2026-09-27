"""Checks of the CM1 runner's pure pieces: patching, namelist editing, grids, surfaces and restarts."""
import numpy as np
import pytest

from climate.crm import cm1_run as c


def test_patch_applies_once_and_is_idempotent():
    text = 'a = 1\nb = 2\nb = 2\n'
    edits = [('', '! head\n', 1), ('a = 1\n', 'a = 1  ! Terluna: marked\n', 1), ('b = 2', 'b = 3', 2)]
    once = c.apply_patch(text, 'demo', 'Terluna: marked', edits)
    assert once == '! head\na = 1  ! Terluna: marked\nb = 3\nb = 3\n'
    assert c.apply_patch(once, 'demo', 'Terluna: marked', edits) == once


def test_patch_refuses_a_different_source():
    with pytest.raises(RuntimeError):
        c.apply_patch('b = 2\n', 'demo', 'Terluna: marked', [('b = 2', 'b = 3  ! Terluna: marked', 2)])
    with pytest.raises(RuntimeError):
        c.apply_patch('b = 2\n', 'demo', 'Terluna: never added', [('b = 2', 'b = 3', 1)])


def test_every_patch_adds_its_marker():
    for name, marker, edits in [*c.MAKEFILE.values(), *c.PATCHES]:
        assert any(marker in new for _, new, _ in edits), name


def test_morrison_fall_speeds_follow_gravity():
    # V ~ g^((B+1)/3): Stokes droplets (B = 2) fall g times slower, hail (B = 0.5) as sqrt(g)
    assert (1 / 6) ** ((2 + 1) / 3) == pytest.approx(1 / 6)
    assert (1 / 6) ** ((0.5 + 1) / 3) == pytest.approx((1 / 6) ** 0.5)
    assert 'AR = AR*(G/9.81)**((BR+1.)/3.)' in c.MORRISON_BLOCK


def test_set_namelist_changes_only_the_named_entry():
    text = ' &param1\n dx     =  3000.0,\n dx_inner = 5.0,\n /\n &param10\n dx     =  1.0,\n /\n'
    out = c.set_namelist(text, 'param1', 'dx', 6011.2)
    assert ' dx     =  6011.2,' in out and ' dx_inner = 5.0,' in out and out.endswith(' dx     =  1.0,\n /\n')
    assert '.true.' in c.set_namelist(' &p\n flag = .false.,\n /\n', 'p', 'flag', True)
    with pytest.raises(KeyError):
        c.set_namelist(text, 'param1', 'dy', 1.0)


def test_vertical_grid_is_fine_below_and_even_aloft():
    z = np.array(c.vertical_grid(dz0=100.0, z_fine=1000.0, stretch=1.08, dz_max=2000.0, ztop=150000.0))
    dz = np.diff(z)
    assert z[0] == 0.0 and z[-1] == pytest.approx(150000.0)
    assert np.all(dz > 0) and dz[:10] == pytest.approx(100.0)
    assert np.all(dz[1:] / dz[:-1] < 1.081) and dz.max() <= 2000.0 * 1.5


def test_ring_surface_segments_cover_the_ring(monkeypatch):
    lon = (np.arange(1440) + 0.5) * 0.25
    lat = 89.875 - 0.25 * np.arange(720)
    water = np.zeros((720, 1440))
    water[:, (lon > 90) & (lon < 180)] = 1.0                           # a quarter of the equator is sea
    from climate.gcm import boundary
    monkeypatch.setattr(boundary, 'lakes_product', lambda *a: dict(lat_deg=lat, lon_deg=lon, water_fraction=water,
                                                                   sha256='test'))
    segments, share, _ = c.ring_surface(nx=360, dx=1000.0, sea_k=299.0, land_k=295.0)
    assert share == pytest.approx(0.25, abs=0.01)
    assert segments[0][0] < 0 and segments[-1][1] > 360 * 1000.0
    assert [s[2] for s in segments] == [1, 2, 1]
    assert all(a[1] == b[0] for a, b in zip(segments, segments[1:]))


def test_latest_restart_needs_a_complete_set(tmp_path):
    for part in 'isuvwx':
        (tmp_path / f'cm1rst_t000001_{part}.dat').write_bytes(b'')
    for part in 'isuv':
        (tmp_path / f'cm1rst_t000002_{part}.dat').write_bytes(b'')
    assert c.latest_restart(tmp_path) == 1
    empty = tmp_path / 'empty'
    empty.mkdir()
    assert c.latest_restart(empty) == 0


def test_storms_wrap_around_the_ring():
    from climate.crm import ring_analysis as ra
    rain = np.zeros(20)
    rain[[0, 1, 18, 19]] = 5.0                                         # one storm across the seam
    rain[[8, 9, 10]] = 2.0
    found = sorted(ra.storms(rain, threshold=1.0))
    assert found == [(8, 3), (18, 4)]
    assert ra.storms(np.full(5, 3.0)) == [(0, 5)] and ra.storms(np.zeros(5)) == []


def test_hour_angle_grows_eastward_and_with_time():
    from climate.crm import ring_analysis as ra
    geo = dict(x=np.array([0.0, 250.0, 500.0]), sun=dict(h0=-90.0, length=1000.0))
    assert ra.hour_angle(geo, 0.0, 100.0) == pytest.approx([-90.0, 0.0, 90.0])
    assert ra.hour_angle(geo, 25.0, 100.0) == pytest.approx([0.0, 90.0, -180.0])


def test_read_snapshot_follows_the_grads_layout(tmp_path):
    from climate.crm import ring_analysis as ra
    ctl = ('dset ^cm1out_t00%y4_s.dat\nxdef 3 linear 0 1\nydef 1 linear 0 1\nzdef 2 levels\n1\n2\nvars 2\n'
           'rain 0 99 rain\nth 2 99 theta\nendvars\n')
    (tmp_path / 'cm1out_s.ctl').write_text(ctl)
    np.arange(9, dtype='<f4').tofile(tmp_path / 'cm1out_t000001_s.dat')
    d = ra.read_snapshot(tmp_path, 1)
    assert d['rain'].tolist() == [0, 1, 2] and d['th'].tolist() == [[3, 4, 5], [6, 7, 8]]


def test_wet_bulb_matches_the_psychrometric_chart():
    from climate.crm import ring_analysis as ra
    t = 303.15
    assert ra.wet_bulb(t, 0.5 * ra.saturation_pa(t), 101325.0) - 273.15 == pytest.approx(22.0, abs=0.2)
    assert ra.wet_bulb(298.15, ra.saturation_pa(298.15), 121600.0) == pytest.approx(298.15, abs=1e-6)
    assert ra.vapour_pa(0.622, 2000.0) == pytest.approx(1000.0)                  # r = 0.622 is half the air


def test_hour_angle_table_keeps_land_and_water_apart():
    from climate.crm import ring_analysis as ra
    bins = ra.HOUR_BINS
    comp = {k: np.zeros((2, bins)) for k in ('t2', 'tsk', 'rh2', 'tw2', 'prate', 'evap', 's10', 'hpbl', 'lcl', 'cape',
                                              'cloud', 'cloud_high', 'cloud_flight', 'cloud_top_m', 'pwat')}
    comp['t2'][0], comp['t2'][1] = 300.0, 290.0
    comp['wind_levels'] = np.ones((2, bins, len(ra.WIND_HEIGHTS_M)))
    table = ra.hour_angle_table(np.linspace(-175.0, 175.0, bins), comp)
    assert table['land']['air_2m_c'][0] == pytest.approx(26.9) and table['water']['air_2m_c'][0] == pytest.approx(16.9)
    assert len(table['land']['wind_by_height_m_s']) == len(ra.WIND_HEIGHTS_M) == len(table['wind_heights_km'])


def test_storm_tracks_follow_overlap_and_measure_drift():
    from climate.crm import ring_analysis as ra
    found = [[(10, 4)], [(12, 4), (40, 3)], [(14, 4)], []]            # one storm drifting east, one brief
    tracks = ra.track_storms(found, nx=50, reach=1)
    assert sorted(len(t) for t in tracks) == [1, 3]
    summary = ra.track_summary(tracks, output_s=3600.0, dx=1000.0, nx=50)
    assert summary['drift_m_s']['median'] == pytest.approx(4 * 1000.0 / (2 * 3600.0))   # 4 columns in 2 hours
    assert summary['lifetime_h']['max'] == pytest.approx(3.0) and summary['seen_once_share'] == pytest.approx(0.5)
    across = ra.track_storms([[(48, 3)], [(0, 3)]], nx=50, reach=1)   # across the seam of the ring
    assert len(across) == 1 and ra.track_summary(across, 3600.0, 1000.0, 50)['drift_m_s']['median'] > 0


def test_rain_spells_count_runs_at_each_place():
    from climate.crm import ring_analysis as ra
    prate = np.zeros((6, 2))
    prate[[1, 2], 0] = 5.0                                              # one 2-snapshot spell on land
    prate[[0, 3, 5], 1] = 5.0                                           # three single spells on water
    out = ra.rain_spells(prate, np.array([True, False]), output_s=3600.0, threshold=1.0)
    assert out['land']['duration_h']['median'] == pytest.approx(2.0)
    assert out['water']['spells_per_lunar_day'] is None                # six hours is not a lunar day
    month = ra.rain_spells(prate, np.array([True, False]), output_s=29.530589 * 86400.0 / 6, threshold=1.0)
    assert month['water']['spells_per_lunar_day'] == pytest.approx(3.0, rel=1e-3)


def test_a_case_ends_at_its_last_whole_restart_interval():
    two_lunar_days = dict(days=59.06, restart_s=43200.0)
    assert c.case_length_s(two_lunar_days, 1.0) == pytest.approx(59.0 * 86400.0)
    scaled = dict(days=20.0, restart_s=86400.0)                         # the Earth-gravity pair case
    s = 0.16556767345166287
    assert c.case_length_s(scaled, s) == pytest.approx(20.0 * 86400.0 * s)


def test_noon_lies_where_the_hour_angle_is_zero():
    from climate.crm import ring_analysis as ra
    geo = dict(x=np.array([0.0]), sun=dict(h0=-90.0, length=1000.0))
    assert ra.noon_x_km(geo, 0.0, 100.0) == pytest.approx(0.25)
    assert ra.noon_x_km(geo, 25.0, 100.0) == pytest.approx(0.0)                 # a quarter day later noon reached x = 0
    x = ra.noon_x_km(geo, 60.0, 100.0) * 1000.0
    assert ra.hour_angle(dict(x=np.array([x]), sun=geo['sun']), 60.0, 100.0)[0] == pytest.approx(0.0, abs=1e-9)


def test_field_product_pools_rain_and_keeps_sections_to_90_km(tmp_path):
    import json
    from climate.crm import ring_analysis as ra
    nz, nx, nt, w = 6, 8, 3, 2 * ra.STORM_HALF_WIDTH
    heights = np.array([1.0, 20.0, 60.0, 89.0, 95.0, 120.0]) * 1000.0
    rain = np.zeros((nt, nx)); rain[1, 5] = 7.0
    arrays = dict(heights_m=heights, hour_angle_deg=np.array([-5.0, 5.0]), x_km=np.arange(nx) * 6.0,
                  section_u=np.ones((2, nz)), section_w=np.zeros((2, nz)), section_cloud=np.zeros((2, nz)),
                  hov_prate_mm_h=rain, hov_time_day=np.arange(nt) * 0.125, hov_noon_x_km=np.zeros(nt),
                  land=np.array([1, 1, 0, 0, 1, 0, 1, 1], bool), storm_day=30.0, storm_hour_angle_deg=45.0,
                  storm_rain_max_mm_h=7.0, storm_condensate_g_kg=np.ones((nz, w)), storm_rain_water_g_kg=np.zeros((nz, w)),
                  storm_updraft_m_s=np.zeros((nz, w)), storm_rain_mm_h=np.zeros(w), storm_wind_10m_m_s=np.zeros(w),
                  storm_air_2m_c=np.zeros(w), storm_land=np.ones(w, bool))
    path = tmp_path / 'fields.npz'
    ra.write_fields(path, dict(case='t', span_days=[0, 1], evidence='e'), arrays)
    f = np.load(path)
    assert json.loads(str(f['metadata']))['schema'] == ra.FIELDS_SCHEMA
    assert f['section_eastward_m_s'].shape == (2, 4) and f['storm_condensate_g_kg'].shape == (5, w)
    assert f['rain_map_max_mm_h'].shape == (nt, 2) and float(f['rain_map_max_mm_h'][1, 1]) == 7.0
    assert f['rain_map_land_share'].tolist() == [0.5, 0.75]


def test_water_mean_profile_skips_land_columns():
    land = np.array([True, False, False])
    snap = dict(th=np.array([[300.0, 290.0, 292.0]]), qv=np.array([[0.02, 0.01, 0.012]]),
                psfc=np.array([1.0e5, 1.2e5, 1.22e5]), t2=np.array([305.0, 298.0, 300.0]), q2=np.array([0.02, 0.014, 0.016]))
    prof = c.water_mean_profile([snap, snap], land)
    assert prof['theta_k'][0] == pytest.approx(291.0) and prof['qv_kg_kg'][0] == pytest.approx(0.011)
    assert prof['surface_pa'] == pytest.approx(1.21e5) and prof['air_2m_k'] == pytest.approx(299.0)
