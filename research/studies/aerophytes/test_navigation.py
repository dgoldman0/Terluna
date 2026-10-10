"""Vector, work, time-average and carbon-accounting checks; no flight claim."""
import hashlib
import json
import math
import pytest

from research.studies.aerophytes import navigation as n
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS


def test_velocity_is_galilean_invariant_and_crosswind_counts():
    a=n.relative_velocity(-4.,0.,-3.,2.)
    b=n.relative_velocity(6.,-7.,7.,-5.)
    assert a==b
    assert a['airspeed_m_s']==pytest.approx(math.sqrt(5))
    assert n.relative_velocity(-4.,0.,-4.,0.)['airspeed_m_s']==0
    assert n.relative_velocity(-4.,0.,4.,0.)['airspeed_m_s']==8


def test_drag_work_and_cubic_scaling():
    a=n.drag_power(2.,1.204)
    assert a['drag_n_m2']==pytest.approx(1.13176)
    assert a['mechanical_w_m2']==pytest.approx(2.26352)
    assert a['input_w_m2']==pytest.approx(9.05408)
    assert a['mechanical_w_m2']==a['drag_n_m2']*2
    assert n.drag_power(4.,1.204)['input_w_m2']==pytest.approx(8*a['input_w_m2'])
    assert n.drag_power(2.,1.204,drag_area_per_projected_area=2.)['input_w_m2']==pytest.approx(2*a['input_w_m2'])
    assert n.drag_power(0.,1.204)['input_w_m2']==0


def test_wind_tolerance_independently_closes_available_power():
    for cd in (.05,.47):
        for eta in (.25,.7):
            for q in (0.,.1,.5,1.,2.):
                speed=n.wind_match_tolerance(q,1.204,cd,eta)
                # Independently close useful work = efficiency * chemical input.
                assert .5*1.204*cd*speed**3==pytest.approx(eta*q)


def test_motion_cost_is_paid_once_from_already_net_margin():
    a=n.motion_budget(1.,1.,1.204,duty_fraction=.1)
    assert a['mean_input_w_m2']==pytest.approx(.113176)
    assert a['remaining_input_w_m2']==pytest.approx(.886824)
    assert a['annual_motion_carbon_equivalent_kg_m2']*40e6==pytest.approx(
        a['mean_input_w_m2']*JULIAN_DAY*JULIAN_YEAR_DAYS)
    assert not n.motion_budget(.1,2.,1.204)['margin_nonnegative']
    assert n.motion_budget(-1.,0.,1.204)['remaining_input_w_m2']==-1.


def test_two_layer_kinematic_match_and_unresolved_crosswind():
    a=n.two_layer_advection(-4.,0.,-2.,0.,-6.,0.)
    assert a['high_layer_time_fraction']==.5
    assert a['minimum_mean_vector_mismatch_m_s']==0
    b=n.two_layer_advection(-4.,0.,-2.,1.,-6.,1.)
    assert b['minimum_mean_vector_mismatch_m_s']==1
    # Beyond the accessible segment chooses its endpoint, not extrapolation.
    c=n.two_layer_advection(-8.,0.,-2.,0.,-6.,0.)
    assert c['high_layer_time_fraction']==1
    assert c['minimum_mean_vector_mismatch_m_s']==2
    same=n.two_layer_advection(-4.,0.,-2.,0.,-2.,0.)
    assert same['minimum_mean_vector_mismatch_m_s']==2


def test_pulsed_correction_conserves_displacement_but_costs_more():
    steady=n.intermittent_correction(.5,1.,1.204)
    pulsed=n.intermittent_correction(.5,.1,1.204)
    assert pulsed['active_airspeed_m_s']*.1==.5
    assert pulsed['mean_input_w_m2']==pytest.approx(100*steady['mean_input_w_m2'])
    assert n.intermittent_correction(0.,.1,1.204)['mean_input_w_m2']==0


def test_light_advantage_is_net_and_incident_light_is_not_created():
    a=n.sunlight_advantage(1.,.4,.2)
    assert a['net_advantage_w_m2']==pytest.approx(.4)
    assert a['energetically_advantageous']
    assert not n.sunlight_advantage(1.,1.)['energetically_advantageous']
    assert n.sunlight_advantage(-1.,0.)['net_advantage_w_m2']==-1


@pytest.mark.parametrize('call',[
    lambda:n.drag_power(-1.,1.), lambda:n.drag_power(1.,0.),
    lambda:n.drag_power(1.,1.,efficiency=1.1),
    lambda:n.drag_power(float('nan'),1.),
    lambda:n.wind_match_tolerance(-1.,1.),
    lambda:n.motion_budget(1.,1.,1.,duty_fraction=2.),
    lambda:n.intermittent_correction(1.,0.,1.),
    lambda:n.relative_velocity(0.,0.,float('inf'),0.)])
def test_nonphysical_inputs_rejected(call):
    with pytest.raises(ValueError):call()


def test_mean_solar_follow_uses_wind_vector_and_shared_geometry():
    from research.studies.aerophytes.photoperiod import sun_follow_speed
    from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS
    speed=sun_follow_speed(0.,10000.)
    assert speed*SYNODIC_MONTH_DAYS*JULIAN_DAY==pytest.approx(2*math.pi*(MOON_RADIUS+10000.))
    match=n.solar_follow_case(0.,10000.,1.204,-speed)
    assert match['airspeed_m_s']==0
    assert match['mean_input_w_m2']==0
    cross=n.solar_follow_case(0.,10000.,1.204,-speed,1.)
    assert cross['airspeed_m_s']==1.
    assert n.solar_follow_case(60.,10000.,1.204,0.)['airspeed_m_s']==pytest.approx(speed/2)


def test_evaluation_provenance_and_conditional_metadata():
    product=n.evaluate()
    assert product['schema']=='terluna.research.aerophyte-navigation/1'
    for path,digest in product['input_hashes'].items():
        assert digest==hashlib.sha256((n.ROOT/path).read_bytes()).hexdigest()[:16]
    assert product['air_10km']['density_kg_m3']==1.204
    assert next(c for c in product['imposed_vector_cases'] if c['name']=='exact_westward_match')['mean_input_w_m2']==0
    assert all(not c['vertical_control_cost_accounted'] for c in product['illustrative_layers'])
    json.dumps(product,allow_nan=False)
