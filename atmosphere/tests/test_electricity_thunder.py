"""Checks of thunder's propagation to the ground (atmosphere/electricity/thunder.py)."""
import numpy as np
import pytest

from atmosphere.electricity import thunder as th


def iso_h(rh_percent, t_k, p_pa=101325.0):
    """ISO 9613-1's molar concentration of water vapour (%) at a relative humidity."""
    c = -6.8346 * (th.ISO_T01 / t_k) ** 1.261 + 4.6151
    return rh_percent * 10.0 ** c / (p_pa / th.ISO_P)


def test_absorption_matches_iso_9613_s_table_at_20_c_and_70_percent():
    # ISO 9613-2 table 2 (20 C, 70 %, one atmosphere), dB/km at 500, 1000, 2000 and 4000 Hz
    h = iso_h(70.0, 293.15)
    got = [1000.0 * th.absorption_db_per_m(f, 293.15, 101325.0, h) for f in (500.0, 1000.0, 2000.0, 4000.0)]
    assert got == pytest.approx([2.8, 5.0, 9.0, 22.9], rel=0.07)
    assert 1000.0 * th.absorption_db_per_m(63.0, 293.15, 101325.0, h) < 0.15


def test_less_oxygen_and_nitrogen_absorb_less_by_relaxation():
    h = iso_h(70.0, 293.15)
    earth = th.absorption_db_per_m(125.0, 293.15, 101325.0, h)
    lunar = th.absorption_db_per_m(125.0, 293.15, 101325.0, h, x_o2=0.1746, x_n2=0.8153)
    assert 0.8 * earth < lunar < earth * 1.05


def test_band_shares_add_up_and_peak_at_the_few_frequency():
    shares = th.band_shares(63.0)
    assert shares.sum() == pytest.approx(1.0, abs=1e-6)
    assert int(np.argmax(shares)) == list(th.BANDS).index(63.0)
    assert th.few_peak_hz(5.0e5, 1.0e5, 340.0) == pytest.approx(96.0, abs=1.0)   # an Earth return stroke


def uniform(c=340.0, top=20000.0, dz=100.0):
    z = np.arange(0.0, top + dz, dz)
    t = np.full(z.size, (c / th.sound_speed(1.0)) ** 2)        # the temperature giving c in dry air
    return th.Atmosphere(z, t, np.full(z.size, 1.0e5), np.zeros(z.size), np.zeros(z.size), np.zeros(z.size))


def test_straight_rays_spread_as_the_inverse_square_and_travel_at_the_speed_of_sound():
    atm = uniform()
    atm.alpha[:] = 0.0
    h = 5000.0
    hits = th.trace(atm, h, n_rays=36001, max_bounces=1)
    edges = np.arange(0.0, 30001.0, 1000.0)
    lv = th.ground_levels([hits], [1.0], 63.0, atm, edges, time_bin=1.0, min_sin=0.0)
    r = np.hypot(lv['x_m'], h)
    doubling = np.where(th.BANDS <= th.COHERENT_HZ, 4.0, 2.0)
    expect = atm.rho_c_ground * doubling * th.band_shares(63.0) / (4.0 * np.pi * r[:, None] ** 2)
    sel = slice(2, 25)
    assert lv['exposure_pa2s'][sel] == pytest.approx(expect[sel], rel=0.03)
    assert lv['first_s'][sel] == pytest.approx(np.hypot(edges[2:25], h) / 340.0, abs=1.01)


def test_a_ray_rising_into_faster_air_comes_back_down_where_its_arc_says():
    z = np.arange(0.0, 20001.0, 100.0)
    c = 330.0 + 0.01 * z                                          # sound speed growing 10 m/s per km
    t = (c / th.sound_speed(1.0)) ** 2
    atm = th.Atmosphere(z, t, np.full(z.size, 1.0e5), np.zeros(z.size), np.zeros(z.size), np.zeros(z.size))
    hits = th.trace(atm, 0.0, n_rays=181, max_bounces=1)          # one ray every degree, from the ground
    up = (hits['bounce'] == 0) & (hits['x'] > 0.0)                # the rays launched upward
    theta = np.radians(np.arange(1.0, 90.0))
    a = np.cos(theta) / 330.0
    turns = a * c.max() >= 1.0                                    # rays that turn below the top
    expect = 2.0 * np.sin(theta[turns]) / (a[turns] * 0.01)
    assert np.sort(hits['x'][up]) == pytest.approx(np.sort(expect), rel=1e-6)


def test_air_cooling_with_height_leaves_a_silent_zone_beyond_the_grazing_ray():
    z = np.arange(0.0, 40001.0, 100.0)
    c = 340.0 - 0.004 * z                                         # cooling about 6.5 K per km
    t = (c / th.sound_speed(1.0)) ** 2
    atm = th.Atmosphere(z, t, np.full(z.size, 1.0e5), np.zeros(z.size), np.zeros(z.size), np.zeros(z.size))
    h = 5000.0
    hits = th.trace(atm, h, n_rays=18001, max_bounces=1)
    grazing = np.sqrt(1.0 - (c[50] / 340.0) ** 2) * 340.0 / 0.004    # the ray that just touches the ground
    assert hits['x'].max() == pytest.approx(grazing, rel=0.01)
    downwind = th.trace(th.Atmosphere(z, t, np.full(z.size, 1.0e5), np.zeros(z.size), 0.004 * z, np.zeros(z.size)),
                        h, azimuth_deg=90.0, n_rays=18001, max_bounces=1)
    assert downwind['x'].max() > 3.0 * grazing                    # a wind growing as fast as the cooling bends it back


def test_audible_range_reads_the_threshold_and_the_background():
    lv = dict(x_m=np.array([1000.0, 2000.0, 3000.0]),
              peak_db=np.array([[np.nan, 70.0, 40.0, 20.0, -20, -20, -20, -20, -20],
                                [np.nan, 50.0, 30.0, 20.0, -20, -20, -20, -20, -20],
                                [np.nan, 50.0, 30.0, 10.0, -20, -20, -20, -20, -20]]),
              peak_dba=np.array([35.0, 25.0, 10.0]))
    assert th.audible_range(lv) == 1000.0                          # 63 Hz at 40 dB beats its 37.5-dB threshold
    assert th.audible_range(lv, quiet_dba=30.0) == 1000.0
