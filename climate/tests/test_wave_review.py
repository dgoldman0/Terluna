"""The review checks' small calculations: sea drag, wave appearance and the first tide estimate."""
import math

import numpy as np

from climate.waves import review_checks as R
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY


def test_charnock_drag_solves_the_open_sea_law():
    for wind in (2.0, 4.0, 6.0, 10.0):
        cd = float(R.charnock_drag(wind, MOON_SURFACE_GRAVITY))
        ustar = math.sqrt(cd) * wind
        z0 = max(R.CHARNOCK * ustar**2 / MOON_SURFACE_GRAVITY, R.Z0_FLOOR)
        assert math.isclose(wind, ustar / R.KAPPA * math.log(10 / z0), rel_tol=1e-9)
        assert math.isclose(float(R.neutral_wind(ustar, MOON_SURFACE_GRAVITY)), wind, rel_tol=1e-9)


def test_lunar_gravity_roughens_the_sea_above_the_floor():
    lunar, earth = R.charnock_drag(np.array([4.0, 8.0]), MOON_SURFACE_GRAVITY), R.charnock_drag(np.array([4.0, 8.0]), STANDARD_GRAVITY)
    assert np.all(lunar > earth)
    # Below the floor's wind the two gravities share one roughness and one drag.
    calm_lunar, calm_earth = R.charnock_drag(0.5, MOON_SURFACE_GRAVITY), R.charnock_drag(0.5, STANDARD_GRAVITY)
    assert math.isclose(float(calm_lunar), float(calm_earth), rel_tol=1e-12)


def test_deep_water_appearance_scales_with_the_square_root_of_gravity():
    length, speed, group = R.deep_water(12.0, MOON_SURFACE_GRAVITY)
    assert math.isclose(length, MOON_SURFACE_GRAVITY * 144 / (2 * math.pi))
    assert math.isclose(group, speed / 2)
    period, earth_speed = R.earth_equivalent(length)
    assert math.isclose(speed / earth_speed, math.sqrt(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY), rel_tol=1e-12)
    assert math.isclose(period / 12.0, math.sqrt(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY), rel_tol=1e-12)


def test_first_tide_estimate_keeps_each_sea_level_on_average():
    lat, lon = np.meshgrid(np.arange(-20.0, 21.0, 5.0), np.arange(-30.0, 31.0, 5.0))
    weights = np.cos(np.radians(lat.ravel()))
    height, settings = R.first_tide_estimate(lat.ravel(), lon.ravel(), weights, days=60, step_days=0.5)
    assert np.allclose((height * weights).sum(axis=1), 0, atol=1e-9)
    assert 12 < settings["earth_coefficient_m"] < 14 and 0.05 < settings["sun_coefficient_m"] < 0.1
    # A sea of one point has no tide relative to itself.
    single, _ = R.first_tide_estimate(np.array([10.0]), np.array([20.0]), np.array([1.0]), days=30)
    assert np.allclose(single, 0)


def test_sea_labels_follow_the_atlas_node_layout():
    grid = np.zeros((720, 1440), dtype=int)
    grid[int((89.875 - 2.375) * 4), int((93.375 - 0.125) * 4)] = 874
    grid[int((89.875 + 10.125) * 4), int((350.125 - 0.125) * 4)] = 433
    assert R.sea_labels_at(np.array([2.375, -10.125]), np.array([93.375, -9.875]), grid).tolist() == [874, 433]
