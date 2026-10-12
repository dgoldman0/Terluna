"""Day lengths, hour attractors and sailing force balances for the sailing component."""
import json
import math

import pytest

from shared.constants import SYNODIC_MONTH_DAYS
from research.studies.aerophytes import sailing as s
from research.studies.aerophytes.photoperiod import sun_follow_speed


def test_relative_period_limits():
    assert s.relative_solar_period(0., 10000., 0.) == pytest.approx(SYNODIC_MONTH_DAYS)
    v = sun_follow_speed(0., 10000.)
    assert s.relative_solar_period(0., 10000., v) == pytest.approx(SYNODIC_MONTH_DAYS/2)
    assert s.relative_solar_period(0., 10000., -v) is None
    assert s.relative_solar_period(0., 10000., -3*v) == pytest.approx(SYNODIC_MONTH_DAYS/2)


def test_attractor_found_where_wind_matches_the_sun():
    v = sun_follow_speed(60., 1000.)
    hours = [-175+10*i for i in range(36)]
    wind = [-v - .5*math.sin(math.radians(h)) for h in hours]  # crossing at 0 and +-180
    found = s.hour_attractors(hours, wind, 60., 1000.)
    stable = [a for a in found if a['stable']]
    assert len(stable) == 1 and stable[0]['hour_angle_deg'] == pytest.approx(0., abs=1e-6)
    assert any(not a['stable'] for a in found)


def test_dwell_counts_time_near_the_slowest_hour():
    v = sun_follow_speed(0., 1000.)
    hours = [-175+10*i for i in range(36)]
    wind = [0.]*36
    row = s.hour_dwell(hours, wind, 0., 1000., 10.)
    assert row['dwell_days'] == pytest.approx(SYNODIC_MONTH_DAYS*20/360, rel=1e-6)


def test_sail_balance_closes_and_drag_only_stays_on_the_shear():
    r = s.tethered_sail(1000., 100., 2., 1.2, 1.3, 1., .1)
    assert r['residual_n'] < 1e-6
    assert r['across_m_s'] > 0
    drag = s.tethered_sail(1000., 100., 2., 1.2, 1.3, 0., 1.)
    assert drag['across_m_s'] == pytest.approx(0., abs=1e-9)
    assert 0 < drag['along_m_s'] < 2
    # Equal drag areas and densities: the body moves at half the shear.
    half = s.tethered_sail(100., 100., 2., 1.2, 1.2, 0., 1.)
    assert half['along_m_s'] == pytest.approx(1.)
    other = s.tethered_sail(1000., 100., 2., 1.2, 1.3, 1., .1, tack=-1)
    assert other['across_m_s'] == pytest.approx(-r['across_m_s'])
    assert s.tethered_sail(1000., 100., 0., 1.2, 1.3)['speed_m_s'] == 0.


def test_reach_and_tether_sizing():
    assert s.meridional_reach(1., 1.) == pytest.approx(86400/s.METRES_PER_DEGREE)
    d = s.tether_diameter(1000., 28e6)
    assert math.pi*d*d/4*28e6 == pytest.approx(1000.)


def test_evaluate_profile_and_rings():
    out = s.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out
    heights = [p['height_km'] for p in out['global_mean_profile'] if p['height_km'] is not None]
    assert heights == sorted(heights, reverse=True)
    assert out['westernmost_mean_wind_m_s'] < 0
    assert set(out['rings']) == set(s.RINGS)


def test_circuit_period_equals_the_steady_period_for_a_steady_wind():
    hours = list(range(-180, 180, 10))
    steady = s.circuit_period(hours, [2.]*len(hours), 0., 10000.)
    assert steady == pytest.approx(s.relative_solar_period(0., 10000., 2.), rel=1e-6)
    varying = s.circuit_period(hours, [2.+1.5*math.sin(math.radians(h)) for h in hours], 0., 10000.)
    assert varying > steady  # the harmonic mean of a varying pace is slower than its mean
    assert s.circuit_period(hours, [-5.]*len(hours), 0., 10000.) is None


def test_hanging_tether_reaches_its_depth_within_its_stress():
    hang = s.hanging_tether(1000., 5000., 28e6, 1500., 1.6)
    w = 1500*1.6*hang['area_m2']
    assert hang['depth_m'] >= .9*5000*(1-1e-9)
    assert hang['top_tension_n'] <= 28e6*hang['area_m2']*(1+1e-9)
    assert hang['top_tension_n'] == pytest.approx(math.hypot(1000., hang['ballast_kg']*1.6+w*5000))
    assert s.hanging_depth(1000., 1e12, w, 5000.) == pytest.approx(5000., rel=1e-6)
    assert s.hanging_tether(1000., 20000., 28e6, 1500., 1.6) is None  # longer than it can hang


def test_lowest_sigma_level_is_extrapolated_below_the_profile():
    h, p = [2.5, 5.], [115618., 110251.]
    z = s.sigma_heights([.9833], 121590., h, p, extrapolate_below=True)[0]
    slope = 2.5/(math.log(110251.)-math.log(115618.))
    assert z == pytest.approx(2.5+(math.log(.9833*121590.)-math.log(115618.))*slope)
    assert 0.5 < z < 1.
