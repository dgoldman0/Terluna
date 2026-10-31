import copy
import json
import math

import pytest

from shared.constants import BOLTZMANN, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from research.studies.world_economy import slow_settling as m


def test_result_reproduces_bound_product():
    assert m.results() == json.loads((m.HERE/'results/slow_settling.json').read_text())


def test_changed_source_rejected(tmp_path):
    s,_ = m.load_inputs()
    s['products'][next(iter(s['products']))]['sha256'] = '0'*64
    p = tmp_path/'inputs.json';p.write_text(json.dumps(s))
    with pytest.raises(ValueError,match='Input changed'):
        m.load_inputs(p)


@pytest.mark.parametrize('g',[0.,MOON_SURFACE_GRAVITY,STANDARD_GRAVITY])
def test_terminal_speed_closes_weight_drag_balance(g):
    rp,rf,d,mu = 3000.,1.44,1e-4,1.8e-5
    r = m.terminal_speed(rp,rf,d,mu,g);v=r['speed_m_s'];re=r['particle_reynolds']
    weight=(rp-rf)*math.pi*d**3/6*g
    drag=3*math.pi*mu*d*v*(1+.15*re**.687)
    assert drag == pytest.approx(weight,rel=1e-12,abs=1e-30)


def test_stokes_limit_and_finite_inertia_differ():
    e=m.terminal_speed(2500,1000,1e-6,.01,STANDARD_GRAVITY)
    l=m.terminal_speed(2500,1000,1e-6,.01,MOON_SURFACE_GRAVITY)
    assert e['speed_m_s']==pytest.approx(e['stokes_speed_m_s'],rel=1e-5)
    assert l['speed_m_s']/e['speed_m_s']==pytest.approx(MOON_SURFACE_GRAVITY/STANDARD_GRAVITY,rel=1e-5)
    e=m.terminal_speed(1500,1000,2e-4,.001,STANDARD_GRAVITY)
    assert not e['stokes_within_five_percent']


def test_neutral_buoyancy_and_creaming_direction():
    assert m.terminal_speed(1000,1000,1e-4,.001,STANDARD_GRAVITY)['speed_m_s']==0
    r=m.terminal_speed(900,1000,1e-4,.001,STANDARD_GRAVITY)
    assert r['direction']=='up' and r['speed_m_s']>0


@pytest.mark.parametrize('bad',[True,-1,float('nan'),float('inf')])
def test_invalid_transport_inputs_rejected(bad):
    with pytest.raises(ValueError):m.terminal_speed(2500,1000,1e-4,.001,bad)


def test_outside_drag_branch_is_rejected():
    with pytest.raises(ValueError,match='Re below 1000'):
        m.terminal_speed(3000,1,1,1e-5,STANDARD_GRAVITY)


def test_column_diffusivity_closes_zero_flux_and_contrast():
    s,_=m.load_inputs();c=s['particles'][0]
    r=m.column_screen(c,MOON_SURFACE_GRAVITY,301,1.2)
    k=r['required_total_vertical_diffusivity_m2_s']
    assert math.exp(r['speed_m_s']*c['height_m']/k)==pytest.approx(1.2)
    assert r['required_additional_vertical_diffusivity_m2_s']+r['brownian_diffusivity_m2_s']==pytest.approx(k)


def test_brownian_motion_survives_freefall_and_controls_small_colloid():
    s,_=m.load_inputs();c=s['particles'][-1]
    z=m.column_screen(c,0,301,1.2)
    assert z['brownian_diffusivity_m2_s']==pytest.approx(BOLTZMANN*301/(3*math.pi*c['viscosity']*c['diameter']))
    assert z['required_additional_vertical_diffusivity_m2_s']==0
    assert z['time_to_cross_height_s'] is None
    assert m.column_screen(c,STANDARD_GRAVITY,301,1.2)['brownian_alone_meets_contrast']


def test_colloid_size_scaling_preserves_brownian_gravity_ratio():
    s,_=m.load_inputs();c=s['particles'][-1]
    e=m.column_screen(c,STANDARD_GRAVITY,301,1.2)
    large=dict(c,diameter=c['diameter']*(STANDARD_GRAVITY/MOON_SURFACE_GRAVITY)**(1/3))
    l=m.column_screen(large,MOON_SURFACE_GRAVITY,301,1.2)
    assert l['brownian_only_log_bottom_top_contrast']==pytest.approx(e['brownian_only_log_bottom_top_contrast'],rel=1e-7)


def test_lunar_filler_pass_does_not_imply_all_particles_stay_suspended():
    r=m.results()['particle_columns']
    assert r[0]['worlds']['moon']['hold_within_drift_limit']
    assert not r[0]['worlds']['earth']['hold_within_drift_limit']
    assert not r[1]['worlds']['moon']['hold_within_drift_limit']


@pytest.mark.parametrize('floor',[1,1.2])
def test_transport_floor_removes_mixer_saving(floor):
    assert m.suspension_power(MOON_SURFACE_GRAVITY/STANDARD_GRAVITY,floor)['turbulent_power_ratio']==1


def test_mixing_extrapolation_counts_only_selected_power():
    g=MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    assert m.suspension_power(g,.25)['turbulent_power_ratio']==pytest.approx(g**1.35)
    assert m.suspension_power(g,.7)['turbulent_power_ratio']==pytest.approx(.7**3)
    assert m.suspension_power(0,.25)['turbulent_power_ratio']==.25**3


def test_passive_harvest_penalty_and_zero_gravity_failure():
    r=m.batch_cycle(86400,7200,900,MOON_SURFACE_GRAVITY/STANDARD_GRAVITY)
    assert r['additional_critical_path_s']==pytest.approx(900*(STANDARD_GRAVITY/MOON_SURFACE_GRAVITY-1))
    assert .95<r['throughput_relative_to_earth']<.96
    short=m.batch_cycle(3600,1800,900,MOON_SURFACE_GRAVITY/STANDARD_GRAVITY)
    assert .58<short['throughput_relative_to_earth']<.59
    assert m.batch_cycle(1,1,1,0)['throughput_relative_to_earth']==0


@pytest.mark.parametrize('g,rho_ratio',[(STANDARD_GRAVITY,1),(MOON_SURFACE_GRAVITY,1.2)])
def test_ergun_pressure_balance(g,rho_ratio):
    s,_=m.load_inputs();c=s['fluid_bed'];r=m.fluidized_bed(c,g,rho_ratio,0)
    u=r['minimum_fluidization_speed_m_s'];e=c['onset_voidage'];d=c['diameter_m']
    gradient=150*c['viscosity_Pa_s']*(1-e)**2*u/(e**3*d*d)+1.75*c['rho_gas']*rho_ratio*(1-e)*u*u/(e**3*d)
    assert gradient*c['onset_height_m']==pytest.approx(r['bed_pressure_drop_Pa'])


def test_reactant_flow_floor_can_destroy_lunar_bed_retention():
    s,_=m.load_inputs();c=s['fluid_bed']
    earth=m.fluidized_bed(c,STANDARD_GRAVITY,1,.1)
    moon=m.fluidized_bed(c,MOON_SURFACE_GRAVITY,1.2,.1)
    assert earth['feasible_flow_interval']
    assert not moon['feasible_flow_interval']
    assert not m.fluidized_bed(c,0,1,0)['feasible_flow_interval']


def test_non_bed_pressure_loss_is_not_multiplied_by_gravity():
    s,_=m.load_inputs();c=s['fluid_bed']
    earth=m.fluidized_bed(c,STANDARD_GRAVITY,1,.05)
    moon=m.fluidized_bed(c,MOON_SURFACE_GRAVITY,1.2,.05)
    ratio=moon['flow_pressure_power_W_m2']/earth['flow_pressure_power_W_m2']
    assert ratio>MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    assert earth['required_superficial_speed_m_s']==moon['required_superficial_speed_m_s']


def test_cost_threshold_counts_electricity_share_once_and_preserves_rotation_competitor():
    s,_=m.load_inputs();g=MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    rows=m.cost_screen(s['costs'],g)
    r=next(x for x in rows if x['mixing_electricity_cost_share']==.05 and x['other_speed_floor']==.25 and x['moon_yield']==.95)
    assert r['lunar_cost_saving_fraction_vs_earth']==pytest.approx(.05*(1-g**1.35))
    assert r['moon_cost_U_good']==r['array_rotating_before_extra_service_U_good']
    assert r['array_freefall_before_extra_service_U_good']<r['moon_cost_U_good']
    assert r['minimum_array_freefall_extra_cost_for_moon_parity_U_good']>0
    assert all(x['lunar_cost_saving_fraction_vs_earth']==0 for x in rows if x['other_speed_floor']==1 and x['moon_yield']==.95)


def test_hypothetical_quality_gain_is_separate_from_physical_mixer_gain():
    s,_=m.load_inputs();rows=m.cost_screen(s['costs'],MOON_SURFACE_GRAVITY/STANDARD_GRAVITY)
    r=next(x for x in rows if x['mixing_electricity_cost_share']==.05 and x['other_speed_floor']==1 and x['moon_yield']==.98)
    assert r['reference_electricity_saving_fraction_of_total']==0
    assert r['lunar_cost_saving_fraction_vs_earth']==pytest.approx(1-.95/.98)
    assert r['moon_cost_U_good']==r['array_rotating_before_extra_service_U_good']
    assert r['moon_yield']==r['rotating_array_yield']


def test_no_gwp_or_market_assignment_and_old_casting_evidence_preserved():
    s,p=m.load_inputs();r=m.results()
    assert r['casting_reference']==p['research/studies/world_economy/results/partial_gravity.json']['casting']
    assert 'GWP' in ' '.join(r['boundaries'])
    assert r['mixing_cost_share_needed_for_25pct_factory_saving']>.27


def test_unbounded_drift_flags_vessel_boundary():
    r=m.results()['particle_columns'][2]['worlds']['earth']
    assert r['unbounded_hold_drift_m']>.1
    assert r['reaches_boundary_during_hold']


def test_whole_cycle_counts_restored_capacity_and_energy_once():
    r=m.cycle_cost(.05,.3,.1,.5)
    assert r['mixing_electricity_saving_U']==pytest.approx(4.5)
    assert r['added_vessel_occupancy_cost_U']==pytest.approx(30)
    assert r['moon_attempt_cost_U']==pytest.approx(125.5)


def test_whole_cycle_recovery_can_reverse_positive_energy_case():
    rows=m.results()['cycle_adjusted_costs']
    long=next(r for r in rows if r['cycle']=='long_batch' and r['mixing_cost_share']==.05)
    short=next(r for r in rows if r['cycle']=='short_batch' and r['mixing_cost_share']==.05)
    assert 0<long['net_cost_saving_fraction']<.04
    assert short['net_cost_saving_fraction']<-.15


def test_immediate_forced_recovery_does_not_erase_its_cost_allowance():
    r=m.cycle_cost(.05,.3,.1,1)
    assert r['added_vessel_occupancy_cost_U']==0
    assert r['incremental_cost_allowance_to_restore_earth_cycle_U']==pytest.approx(4.5)


def test_whole_cycle_rejects_overlapping_shares_and_zero_throughput():
    with pytest.raises(ValueError):m.cycle_cost(.8,.3,.1,1)
    with pytest.raises(ValueError):m.cycle_cost(.05,.3,.1,0)
