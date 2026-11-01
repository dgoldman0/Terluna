"""Independent physical limits for the declared stagnant compact-crystal control."""
import math

import pytest

from atmosphere.snow import facet_growth as growth


def sphere_state(radius=20e-6, alpha=.05):
    # Explicit idealized state, independent of the external coefficient helper.
    return dict(capacitance_m=radius, boundary_capacitance_m=radius+60e-9,
        diffusivity_m2_s=2e-5, conductivity_W_m_K=.023,
        vapor_molecular_speed_m_s=550., alpha_a=alpha, alpha_c=alpha,
        thermal_jump_length_m=1e-7, vapor_gas_constant_J_kg_K=461.5,
        latent_sublimation_J_kg=2.83e6, ice_saturation_pressure_Pa=165.)


def test_sphere_recovers_scalar_kinetic_resistances():
    radius=20e-6
    state=sphere_state(radius)
    dv,ka=growth.effective_transport(state,radius,radius)
    ratio=radius/state['boundary_capacitance_m']
    # Resistances in series for a sphere, independently arranged.
    vapor_resistance=ratio/state['diffusivity_m2_s']+4/(state['alpha_a']*550*radius)
    heat_resistance=ratio/state['conductivity_W_m_K']+state['thermal_jump_length_m']/(.023*radius)
    assert dv==pytest.approx(1/vapor_resistance,rel=1e-14)
    assert ka==pytest.approx(1/heat_resistance,rel=1e-14)


def test_sphere_recovers_free_molecular_attachment_limit():
    radius=20e-6;alpha=1e-12;temperature=258.15;si=.05
    state=sphere_state(radius,alpha)
    state['latent_sublimation_J_kg']=2.83e6
    result=growth.mass_and_axis_rates(state,radius,radius,temperature,si)
    rhoeq=165/(461.5*temperature)
    # Molecular impingement v/4 onto a sphere's 4*pi*r^2 surface.
    expected=math.pi*radius**2*alpha*550*rhoeq*si
    assert result['mass_rate_kg_s']==pytest.approx(expected,rel=1e-8)


def test_zero_supersaturation_stops_mass_and_axes():
    result=growth.mass_and_axis_rates(sphere_state(),20e-6,20e-6,258.15,0.)
    assert result['mass_rate_kg_s']==0.
    assert result['equatorial_rate_m_s']==result['polar_rate_m_s']==0.


@pytest.mark.parametrize('a,c,aa,ac',[(40e-6,4e-6,.2,.01),(4e-6,40e-6,.01,.2)])
def test_mass_equals_density_times_geometrical_volume_rate(a,c,aa,ac):
    state=sphere_state(10e-6)
    state.update(alpha_a=aa,alpha_c=ac)
    density=919.
    rates=growth.mass_and_axis_rates(state,a,c,258.15,.05,density)
    da=rates['equatorial_rate_m_s'];dc=rates['polar_rate_m_s']
    volume_rate=(4*math.pi/3)*(2*a*c*da+a*a*dc)
    assert density*volume_rate==pytest.approx(rates['mass_rate_kg_s'],rel=1e-14)
    assert dc/da==pytest.approx(ac/aa,rel=1e-14)


def test_supersaturation_scales_fixed_state_mass_flux():
    state=sphere_state()
    low=growth.mass_and_axis_rates(state,20e-6,20e-6,258.15,.01)
    high=growth.mass_and_axis_rates(state,20e-6,20e-6,258.15,.04)
    assert high['mass_rate_kg_s']==pytest.approx(4*low['mass_rate_kg_s'],rel=1e-14)


@pytest.mark.parametrize('temperature',[0.,-1.,float('nan'),float('inf')])
def test_invalid_temperature_is_rejected(temperature):
    with pytest.raises(ValueError):
        growth.mass_and_axis_rates(sphere_state(),20e-6,20e-6,temperature,.01)


@pytest.mark.parametrize('field,value',[('alpha_a',0.),('alpha_c',-1.),
    ('alpha_a',1.01),('capacitance_m',float('nan')),('boundary_capacitance_m',0.),
    ('diffusivity_m2_s',-1.),('conductivity_W_m_K',0.),('thermal_jump_length_m',-1.)])
def test_invalid_transport_state_is_rejected(field,value):
    state=sphere_state();state[field]=value
    with pytest.raises(ValueError):
        growth.effective_transport(state,20e-6,20e-6)


def test_sublimation_is_outside_this_growth_control():
    with pytest.raises(ValueError):
        growth.mass_and_axis_rates(sphere_state(),20e-6,20e-6,258.15,-.01)
