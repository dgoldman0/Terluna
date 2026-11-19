"""Checks of the sand study: the transport laws at lunar gravity, the turning of wind directions between the turned
Moon and today's, the drift sums on made-up winds, the carrying scheme's balance, and the written product."""
import hashlib
import json
import math

import numpy as np
import pytest

from research.studies.atlas_round import sand as S

G_E, RHO_E = S.G_EARTH, S.AIR_DENSITY_EARTH


def rotation_for_test():
    """The design surface's rotation where the drive holds it, else a turn of 6.7 degrees about a tilted axis."""
    if S.SURFACE.exists():
        z = np.load(S.SURFACE, allow_pickle=False)
        return np.asarray(json.loads(str(z['metadata']))['frame']['rotation_today_to_turned'], float)
    axis = np.array([0.55, -0.83, 0.05]); axis /= np.linalg.norm(axis)
    a = math.radians(6.7)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(a) * k + (1 - math.cos(a)) * (k @ k)


def bearing(lat1, lon1, lat2, lon2):
    """Initial bearing (degrees clockwise from north) of the great circle from place 1 to place 2."""
    p1, p2, dl = math.radians(lat1), math.radians(lat2), math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x)) % 360.0


def apart(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def test_threshold_scales_with_gravity_and_air():
    for rho_p in S.GRAINS.values():
        for gamma in S.SHAO_LU_GAMMA:
            d_m, u_m = S.easiest_grain(rho_p, gamma=gamma)
            d_e, u_e = S.easiest_grain(rho_p, RHO_E, G_E, gamma)
            # the lowest threshold goes as g^(1/4) rho_a^(-1/2), whatever the cohesion and the grain's density
            expected = (S.G_MOON / G_E) ** 0.25 * (S.AIR_DENSITY / RHO_E) ** -0.5
            assert u_m / u_e == pytest.approx(expected, rel=2e-3)
            assert d_m / d_e == pytest.approx(math.sqrt(G_E / S.G_MOON), rel=2e-3)
    d, u = S.easiest_grain(S.GRAINS['feldspar'])
    assert 240e-6 < d < 280e-6 and 0.13 < u < 0.15
    # the easiest grain is a minimum of the threshold curve
    assert S.fluid_threshold(d, S.GRAINS['feldspar']) < min(S.fluid_threshold(0.5 * d, S.GRAINS['feldspar']),
                                                              S.fluid_threshold(2.0 * d, S.GRAINS['feldspar']))
    # the record's figure: bare loose ground starts to move at 3.2-4.0 m/s at 10 m (quartz, 1.44 kg/m3, roughness 1-0.1 mm)
    uq = S.easiest_grain(2650.0, 1.44)[1]
    assert S.wind_at_10m(uq, 1e-3) == pytest.approx(3.2, abs=0.1) and S.wind_at_10m(uq, 1e-4) == pytest.approx(4.0, abs=0.1)


def test_flux_law():
    u_it = 0.12
    assert S.saltation_flux(0.119, u_it) == 0.0 and S.saltation_flux(u_it, u_it) == 0.0
    q = S.saltation_flux(np.array([0.13, 0.2, 0.3, 0.5]), u_it)
    assert np.all(np.diff(q) > 0)
    # at equal friction velocity and threshold the Moon carries (rho_M / rho_E)(g_E / g_M) times Earth's sand
    ratio = S.saltation_flux(0.3, u_it) / S.saltation_flux(0.3, u_it, RHO_E, G_E)
    assert ratio == pytest.approx(S.AIR_DENSITY / RHO_E * G_E / S.G_MOON, rel=1e-9)
    # Earth's 250 um sand at u* 0.4 m/s and a threshold of 0.2: about 0.015 kg/m/s
    assert S.saltation_flux(0.4, 0.2, RHO_E, G_E) == pytest.approx(0.015, rel=0.05)
    assert S.wind_at_10m(S.friction_velocity(4.0)) == pytest.approx(4.0)


def test_wind_factors():
    assert S.sand_wind_factor() == pytest.approx(S.TO_SAND[2], abs=0.01)                    # both profiles neutral
    assert S.sand_wind_factor(model_share=0.52) == pytest.approx(S.TO_SAND[0], abs=0.01)   # the run's unstable share
    assert S.TO_SAND[0] < S.TO_SAND[1] < S.TO_SAND[2]
    assert 1.1 < S.ring_to_sand() < S.RING_TO_SAND < S.ring_to_sand(level=313.0) < 1.3


def test_wind_along_a_turned_parallel_turns_with_its_place():
    r = rotation_for_test()
    assert np.allclose(r @ r.T, np.eye(3), atol=1e-9)
    pole_lat, pole_lon = S.to_today(90.0, 0.0, r)            # the turned Moon's north pole on today's Moon
    for lat_t, lon_t in ((0.0, 0.0), (0.0, 123.0), (45.0, -60.0), (-45.0, 200.0), (70.0, 10.0), (-80.0, 300.0),
                         (85.0, 95.0), (20.0, 181.0)):
        lat, lon = (float(x) for x in S.to_today(lat_t, lon_t, r))
        east = float(S.turn_azimuth(90.0, lat_t, lon_t, r))
        # a step east along the turned parallel, seen from today's Moon
        lat2, lon2 = (float(x) for x in S.to_today(lat_t, lon_t + 0.01, r))
        assert apart(east, bearing(lat, lon, lat2, lon2)) < 0.02
        # the turned pole lies a quarter turn to the left of a wind blowing east round it
        assert apart(east - 90.0, bearing(lat, lon, float(pole_lat), float(pole_lon))) < 1e-6
        # north, south and west turn by the same angle, and the turn is the angle between the two norths
        turn = (east - 90.0) % 360.0
        for az in (0.0, 180.0, 270.0, 37.0):
            assert apart(float(S.turn_azimuth(az, lat_t, lon_t, r)), az + turn) < 1e-6
    # the vector keeps its length, and no rotation leaves every direction alone
    v = S.vector_to_today(np.array([3.0]), np.array([-4.0]), np.array([33.0]), np.array([71.0]), r)
    assert np.linalg.norm(v) == pytest.approx(5.0)
    assert float(S.turn_azimuth(123.0, 12.0, 34.0, np.eye(3))) == pytest.approx(123.0)
    # the pole has moved by the stated 6.7 degrees, so directions near the old pole turn by tens of degrees
    assert 6.0 < 90.0 - float(pole_lat) < 7.5


def made_up(u, v, rain_mm_day=0.0, records=40, shape=(2, 3)):
    full = lambda x: np.full((records,) + shape, x, np.float32)
    return dict(uas=full(u), vas=full(v), pr=full(rain_mm_day / 86400.0 / 1000.0))


def test_drift_sums():
    calm = S.drift(made_up(0.5, 0.0), 'none')
    assert np.all(calm['dp_m2_yr'] == 0.0) and np.all(calm['time_above'] == 0.0)
    east = S.drift(made_up(3.0, 0.0), 'none')                # 4.8 m/s over sand, above the threshold all year
    assert np.all(east['time_above'] == 1.0)
    assert np.allclose(east['east_m2_yr'], east['dp_m2_yr']) and np.allclose(east['north_m2_yr'], 0.0, atol=1e-9)
    assert east['by_direction_m2_yr'][S.DIRECTIONS // 4].min() > 0 and east['by_direction_m2_yr'][0].max() == 0.0
    u = S.friction_velocity(3.0 * S.TO_SAND[1])
    expected = float(S.saltation_flux(u, east['ustar_it'])) * S.SECONDS_PER_YEAR / S.BULK_DENSITY
    assert east['dp_m2_yr'][0, 0] == pytest.approx(expected, rel=1e-6)
    north = S.drift(made_up(0.0, 3.0), 'none')
    assert north['by_direction_m2_yr'][0].min() > 0 and np.allclose(north['east_m2_yr'], 0.0, atol=1e-9)
    # gusts add to what a steady wind carries, and storms add drift without direction
    assert S.drift(made_up(1.5, 0.0), 'central')['dp_m2_yr'][0, 0] > S.drift(made_up(1.5, 0.0), 'none')['dp_m2_yr'][0, 0]
    wet = S.drift(made_up(0.5, 0.0, rain_mm_day=1.0), 'central')
    assert wet['storm_m2_yr'][0, 0] > 0 and np.allclose(np.hypot(wet['east_m2_yr'], wet['north_m2_yr']), 0.0, atol=1e-9)
    assert np.allclose(S.weibull_nodes(4.0).mean(), 1.0)


def test_cover_and_classes():
    assert S.plant_cover(np.array([0.028]))[0] < 0.04 < 0.1 < S.plant_cover(np.array([0.12]))[0] < 0.3
    assert S.plant_cover(np.array([2.0]))[0] > 0.97
    assert list(S.class_of(np.array([0.01, 0.1, 0.3, 0.6, 0.8, 1.5]))) == [0, 1, 2, 3, 4, 5]
    assert S.erosion_m_per_yr(np.array([0.0, 52.5, 63.75, -80.0])).tolist() == pytest.approx([79e-6, 35e-6, 30e-6, 25e-6])
    # the kept crust's share of the open ground rises through the divide between the hyper-arid and the arid land
    crust = 1.0 / (1.0 + (S.CRUST_EDGE / np.array([0.02, 0.05, 0.12])) ** S.CRUST_STEEP)
    assert crust[0] < 0.05 and crust[1] == pytest.approx(0.5) and crust[2] > 0.95


def test_carry_keeps_the_sand_and_moves_it_downwind():
    shape = (S.ROWS, S.COLS)
    free = np.zeros(shape, np.float32)
    row, col = S.ROWS // 2, 1000
    free[row, col] = 4.0                                         # 4 m of sand in one cell on the equator
    years = 1000.0
    q = np.full(shape, 3.5 * 8.0 * S.CELL_M / years, np.float32)   # enough to carry an 8 m thick body three cells
    az = np.full(shape, 90.0, np.float32)
    zero = np.zeros(shape, np.float32)
    water = np.zeros(shape, bool)
    bowl, cap = np.zeros(shape, np.float32), np.zeros(S.ROWS * S.COLS, np.float32)
    lost = S.carry(free, q, az, zero, water, bowl, cap, years, log=lambda *_: None)
    assert lost == 0.0 and free.sum() == pytest.approx(4.0, rel=1e-5)
    assert free[row, col + 3] == pytest.approx(4.0, rel=1e-5) and free[row, col] == 0.0
    # a body fans out where the winds vary, and keeps its sand
    fan = np.zeros(shape, np.float32); fan[row, col] = 4.0
    S.carry(fan, q, az, zero, water, bowl, cap, years, side=np.full(shape, 0.25, np.float32), log=lambda *_: None)
    assert fan.sum() == pytest.approx(4.0, rel=1e-3) and fan[row, col + 3] < 2.0 and np.count_nonzero(fan) > 3
    # a pile too deep for one cell is shared over the land round it, round on the ground, and keeps its sand
    ground = S.Ground(water)
    pile = np.zeros(shape, np.float32); pile[row, col] = 40.0
    settled = S.settle(pile, ground)
    assert settled.sum() == pytest.approx(40.0, rel=1e-3) and settled.max() < 4.0
    assert settled[row, col + 4] == pytest.approx(settled[row + 4, col], rel=0.05)
    # sand that meets water is lost
    free[:] = 0.0; free[row, col] = 4.0
    water[row, col + 2] = True
    lost = S.carry(free, q, az, zero, water, bowl, cap, years, log=lambda *_: None)
    assert free.sum() == pytest.approx(0.0, abs=1e-6)
    assert lost == pytest.approx(4.0 * S.cell_area_km2()[row] * 1e6 / 1e9, rel=1e-5)


def test_product_and_record():
    if not (S.PRODUCT.exists() and S.OUT.exists()):
        pytest.skip('the product is on the research drive')
    record = json.loads(S.OUT.read_text())
    assert record['schema'] == S.SCHEMA
    path = S.ROOT / S.CODE[0]
    if record['producer']['files'][S.CODE[0]] != hashlib.sha256(path.read_bytes()).hexdigest():
        pytest.skip('sand.py has changed since the record was written')
    z = np.load(S.PRODUCT, allow_pickle=False)
    meta = json.loads(str(z['metadata']))
    assert meta['schema'] == S.SCHEMA and meta['grid'] == dict(pixels_per_degree=16, rows=2880, cols=5760)
    assert record['product']['sha256'] == S.sha256(S.PRODUCT)
    kinds = dict(sand_share=np.float32, active_share=np.float32, form=np.uint8, crest_azimuth_deg=np.float32,
                 crest_spacing_m=np.float32, drift_azimuth_deg=np.float32)
    area = S.cell_area_km2()[:, None]
    water = None
    if S.SURFACE.exists():
        water = np.zeros(S.ROWS * S.COLS, bool)
        water[np.load(S.SURFACE, allow_pickle=False)['water_cell']] = True
        water = water.reshape(S.ROWS, S.COLS)
    for name in meta['variants']:
        suffix = '' if name == 'central' else f'__{name}'
        f = {k: z[k + suffix] for k in kinds}
        for k, kind in kinds.items():
            assert f[k].shape == (S.ROWS, S.COLS) and f[k].dtype == kind, (name, k)
        share, active, form = f['sand_share'], f['active_share'], f['form']
        assert share.min() >= 0.0 and share.max() <= 1.0 and active.min() >= 0.0 and active.max() <= 1.0
        assert form.max() <= 5 and np.all((form > 0) == (share > 0))
        if water is not None:
            assert share[water].max() == 0.0
        dunes = form >= 2
        assert np.all(f['crest_spacing_m'][~dunes] == 0.0) and (not dunes.any() or f['crest_spacing_m'][dunes].min() > 0.0)
        lined = np.isin(form, (2, 3, 5))
        assert np.all(np.isfinite(f['crest_azimuth_deg'][lined])) and np.all(np.isnan(f['crest_azimuth_deg'][~lined]))
        assert np.all(np.isnan(f['drift_azimuth_deg'][(share == 0) | (active == 0)]))
        moving = np.isfinite(f['drift_azimuth_deg'])
        assert not moving.any() or (f['drift_azimuth_deg'][moving].min() >= 0.0 and f['drift_azimuth_deg'][moving].max() < 360.0)
        assert float((share * area).sum()) == pytest.approx(record['variants'][name]['sand_km2'], rel=1e-3, abs=2.0)
        # away from water and from the edges of shelter, few neighbouring cells differ by more than a few hundredths
        assert record['variants'][name]['neighbour_steps']['share_of_sandy_pairs'] < 0.05
