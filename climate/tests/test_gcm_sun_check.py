"""Checks of the GCM Sun check's pure pieces: PlaSim's model day, its solar-position arithmetic and the annual means."""
import numpy as np

from climate.gcm import sun_check as s

MOON_ROTSPD = 0.03660099581793362                                        # the A runs' namelist: 1 / 27.3217 days


def test_model_day_is_the_solar_day_plasim_derives():
    assert s.plasim_day_steps(MOON_ROTSPD, 12, 30.0) == 1430             # 12/11 of the rotation, 29.8 days
    assert s.plasim_day_steps(1.0, 360, 30.0) == 48


def test_solang_sweeps_evenly_only_when_its_model_day_is_one_turn():
    moon = s.solang_mean_mu(0.0, 64, 1430, 30.0, MOON_ROTSPD)
    assert np.isclose(moon.mean(), 1 / np.pi, rtol=2e-3)                # each row's mean stays right
    assert np.isclose(moon.max() / moon.min(), 1.28, atol=0.01)          # a 33-degree band gets noon twice
    flat = np.isclose(moon, moon.min(), rtol=1e-3)
    assert flat[0] and not flat[28]                                      # 0 E in the flat part, 157.5 E in the band
    even = s.solang_mean_mu(0.0, 64, 48, 30.0, 1.0)
    assert np.ptp(even) / even.mean() < 5e-3


def test_the_patched_clocks_turn_the_sun_evenly():
    one_turn = s.solang_mean_mu(0.0, 64, 1430, 30.0, MOON_ROTSPD, clock='model_day')
    assert np.ptp(one_turn) / one_turn.mean() < 1e-3
    year = dict(sidereal_day_s=86400.0 / MOON_ROTSPD, sidereal_year_s=365.25636 * 86400.0)
    synodic_s = 1.0 / (1.0 / year['sidereal_day_s'] - 1.0 / year['sidereal_year_s'])
    assert np.isclose(synodic_s / 86400.0, 29.5306, atol=1e-3)            # the real lunar day
    steps = np.arange(int(20 * synodic_s / 1800.0))                      # 20 lunar days of half-hour steps
    synodic = s.solang_mean_mu(0.0, 64, 1430, 30.0, MOON_ROTSPD, clock='synodic', steps=steps, **year)
    assert np.ptp(synodic) / synodic.mean() < 1e-3


def test_annual_means_keep_the_global_sunlight_and_spread_it_by_tilt():
    lat = np.linspace(-89.5, 89.5, 180)
    weight = np.cos(np.radians(lat))
    for tilt in (0.0, 1.54, 23.7):
        q = s.annual_mean_insolation(lat, tilt, 0.0167, 1173.5)
        assert np.isclose(np.sum(q * weight) / weight.sum(), 1173.5 / 4 / np.sqrt(1 - 0.0167 ** 2), rtol=2e-3)
    equator = s.annual_mean_insolation(0.0, 23.7, 0.0, np.pi)[0]
    assert np.isclose(equator, 0.958, atol=0.002)                        # Earth's tilt takes 4% from the equator
    assert np.isclose(s.annual_mean_insolation(0.0, 0.0, 0.0, np.pi)[0], 1.0)
