"""Independent geometry, dimensional checks and evidence-boundary checks."""
import json
import math

import pytest

from shared.constants import JULIAN_DAY, MOON_RADIUS, SYNODIC_MONTH_DAYS
from research.studies.floater_viability import photoperiod as p


def test_sun_follow_traverses_one_latitude_circle_per_synodic_cycle():
    velocity=p.sun_follow_speed(60,20000)
    circumference=2*math.pi*(MOON_RADIUS+20000)*.5
    assert velocity*SYNODIC_MONTH_DAYS*JULIAN_DAY==pytest.approx(circumference)
    assert p.sun_follow_speed(-60,20000)==pytest.approx(velocity)
    assert p.sun_follow_speed(90)==0
    assert p.sun_follow_speed(0,20000)>p.sun_follow_speed(0,0)


def test_solar_geometry_matches_independent_limiting_cases():
    assert p.solar_height(60,0)==pytest.approx(30)
    assert p.solar_height(60,180)==pytest.approx(-30)
    assert p.solar_height(60,90)==pytest.approx(0,abs=1e-12)
    assert p.solar_height(90,73,1.5)==pytest.approx(1.5)
    assert p.solar_height(60,0,1.5)==pytest.approx(31.5)


def test_drift_symmetry_inverse_error_and_stationary_check():
    a=p.drift_window(30,20000,1)
    b=p.drift_window(30,20000,-2)
    assert a['centre_to_edge_days']==pytest.approx(2*b['centre_to_edge_days'])
    assert a['hour_angle_drift_deg_per_day']<0<b['hour_angle_drift_deg_per_day']
    assert p.drift_window(30,20000,0)['centre_to_edge_days'] is None
    # Stationary organism starts at noon; a90deg edge would take one quarter month.
    stopped=p.drift_window(30,20000,-p.sun_follow_speed(30,20000),90)
    assert stopped['centre_to_edge_days']==pytest.approx(SYNODIC_MONTH_DAYS/4)
    assert stopped['relative_solar_period_days']==pytest.approx(SYNODIC_MONTH_DAYS)
    slow=p.drift_window(30,20000,-p.sun_follow_speed(30,20000)/2)
    assert slow['relative_solar_period_days']==pytest.approx(2*SYNODIC_MONTH_DAYS)
    assert slow['equinox_geometric_night_days']==pytest.approx(SYNODIC_MONTH_DAYS)
    with pytest.raises(ValueError):p.drift_window(90,0,1)


def test_threshold_duty_distinguishes_daylight_and_extended_twilight():
    assert p.threshold_duty(60,0)['fraction_above']==pytest.approx(.5)
    assert p.threshold_duty(70,-30)['fraction_above']==1
    assert p.threshold_duty(70,30)['fraction_above']==0
    assert p.threshold_duty(90,0,1)['fraction_above']==1
    assert p.threshold_duty(90,0,-1)['fraction_above']==0
    assert p.threshold_duty(60,-20)['fraction_above']>.5


def test_ground_spectra_preserve_grid_values_and_exclude_bolometric_claim():
    light=p.ground_light(30)
    assert light['par_w_m2']==133.268
    assert light['par_umol_m2_s']==624.567
    assert light['optical_202_1000_w_m2']==pytest.approx(206.52927803)
    # Full stored band splits independently into PAR and all other bins.
    assert light['optical_202_1000_w_m2']-light['nonpar_202_1000_w_m2']==pytest.approx(light['spectral_par_w_m2'])
    # Spectral bins and pre-binning source summaries have slightly different reduction precision.
    assert light['spectral_par_w_m2']==pytest.approx(light['par_w_m2'],rel=2e-4)
    assert p.ground_light(20)['par_w_m2']<p.ground_light(25)['par_w_m2']<light['par_w_m2']
    with pytest.raises(ValueError):p.ground_light(0)


def test_twilight_threshold_is_not_compensation_threshold():
    h=p.twilight_threshold(1)
    assert -30<h<-25
    assert p.threshold_duty(60,h)['hours_below_per_stationary_synodic_cycle']<40
    assert p.threshold_duty(70,h)['fraction_above']==1
    with pytest.raises(ValueError):p.twilight_threshold(.01)


def test_twilight_energy_is_a_declared_spectral_shape_not_zero_night():
    row=p.twilight_energy(-10)
    assert row['par_umol_m2_s']==25.2881
    # Any400–700nm spectrum lies inside these independent photon-energy bounds.
    assert .170<row['par_w_m2']/row['par_umol_m2_s']<.300
    assert row['optical_202_1000_w_m2']>row['par_w_m2']>0
    assert p.twilight_energy(-35)['optical_202_1000_w_m2']>0
    assert p.required_irradiance(.07,.01)==pytest.approx(7)
    with pytest.raises(ValueError):p.twilight_energy(-40)


def test_evaluation_preserves_plant_night_labels_and_no_infinite_json():
    out=p.evaluate()
    assert json.loads(json.dumps(out,allow_nan=False))==out
    near=next(r for r in out['plant_reference'] if r['side']=='near' and r['latitude_deg']==60 and r['strategy']=='earth_like')
    far=next(r for r in out['plant_reference'] if r['side']=='far' and r['latitude_deg']==60 and r['strategy']=='earth_like')
    assert near['par_below_1_hours']==far['par_below_1_hours']==33.7141
    assert near['maintenance_deficit_hours']==308.718
    assert far['maintenance_deficit_hours']==307.487
    assert near['sun_below_horizon_respiration_g_c_m2_cycle']==49.2869
    assert near['maintenance_deficit_hours']>9*near['par_below_1_hours']
