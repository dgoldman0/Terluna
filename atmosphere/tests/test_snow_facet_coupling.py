"""Physical balances and limiting solutions for the implicit facet calculation."""
import math
import pytest
from atmosphere.snow import facet_coupling as model
from atmosphere.tests.test_snow_facet_growth import sphere_state


def state(radius=20e-6):
    result = sphere_state(radius)
    result.update(characteristic_supersaturation_a=.01, characteristic_supersaturation_c=.02)
    return result


@pytest.mark.parametrize('m', [1., 3., 10., 15.])
def test_attachment_limits_and_monotonicity(m):
    values = [model.efficiency(s, .01, m) for s in (0., 1e-9, 1e-5, .01, .1, 100.)]
    assert values == sorted(values)
    assert values[0] == 0 and 0 < values[1] < 1e-6
    assert values[-1] == pytest.approx(1., abs=4e-9)
    assert model.efficiency(.01, .01, m) == pytest.approx(math.tanh(1.))


@pytest.mark.parametrize('a,c', [(20e-6, 20e-6), (50e-6, 5e-6), (5e-6, 50e-6)])
@pytest.mark.parametrize('fv,ft', [(1., 1.), (2., 1.5), (100., 100.)])
def test_surface_and_heat_balances(a,c,fv,ft):
    source = state()
    solution = model.solve(source, a,c,258.15,.02, vapor_ventilation=fv,heat_ventilation=ft)
    assert 0 < solution['surface_supersaturation_a'] <= solution['thermally_reduced_supersaturation'] <= .02
    assert 0 < solution['surface_supersaturation_c'] <= solution['thermally_reduced_supersaturation']
    assert max(map(abs,solution['relative_surface_balance_residuals'])) < 2e-10
    assert abs(solution['relative_thermal_balance_residual']) < 2e-10
    # Kinetic flux and continuum diffusion must agree separately on each axis.
    for axis, area, dx in zip(('a','c'), (a*c,a*a),solution['axis_diffusivities_m2_s']):
        kinetic = solution['alpha_'+axis] * source['vapor_molecular_speed_m_s']/4 * solution['surface_supersaturation_'+axis]
        diffusive = source['capacitance_m']/area * dx * solution['thermally_reduced_supersaturation']
        assert kinetic == pytest.approx(diffusive,rel=2e-10)
    rates=model.axis_rates(solution,a,c,919.)
    volume_rate=4*math.pi/3*(2*a*c*rates['equatorial_rate_m_s']+a*a*rates['polar_rate_m_s'])
    assert 919*volume_rate == pytest.approx(solution['mass_rate_kg_s'])


def test_spherical_dislocation_limit_has_independent_quadratic_solution():
    r=20e-6;source=state(r);source['characteristic_supersaturation_c']=.01
    si=1e-6;t=258.15
    solution=model.solve(source,r,r,t,si,1.)
    dv=source['diffusivity_m2_s'];ka=source['conductivity_W_m_K'];cap=source['capacitance_m']
    ratio=cap/source['boundary_capacitance_m'];beta=4*dv*cap/(550*r*r)
    kstar=ka/(ratio+source['thermal_jump_length_m']*cap/(r*r))
    latent=source['latent_sublimation_J_kg'];rv=source['vapor_gas_constant_J_kg_K']
    thermal=165/(rv*t)*latent/t*(latent/(rv*t)-1)
    coefficient=(ratio/beta+thermal*dv/(kstar*beta))/.01
    expected=2*si/(1+math.sqrt(1+4*coefficient*si))
    assert solution['surface_supersaturation_a']==pytest.approx(expected,rel=2e-10)
    assert solution['alpha_a']==pytest.approx(expected/.01,rel=2e-10)


def test_spherical_symmetry_is_preserved():
    source=state();source['characteristic_supersaturation_c']=.01
    solution=model.solve(source,20e-6,20e-6,258.15,.02)
    assert solution['alpha_a']==solution['alpha_c']
    rates=model.axis_rates(solution,20e-6,20e-6,919.)
    assert rates['equatorial_rate_m_s']==rates['polar_rate_m_s']


def test_ventilation_cannot_remove_surface_attachment_resistance():
    source=state();source['characteristic_supersaturation_c']=.01
    si=1e-5;r=20e-6
    low=model.solve(source,r,r,258.15,si,10.)
    high=model.solve(source,r,r,258.15,si,10.,vapor_ventilation=1e8,heat_ventilation=1e8)
    rhoeq=165/(461.5*258.15)
    kinetic_bound=math.pi*r*r*550*rhoeq*si*model.efficiency(si,.01,10.)
    assert low['mass_rate_kg_s'] <= high['mass_rate_kg_s'] <= kinetic_bound*(1+1e-12)
    assert high['mass_rate_kg_s']==pytest.approx(kinetic_bound,rel=1e-10)


def test_zero_supersaturation_is_equilibrium():
    solution=model.solve(state(),20e-6,20e-6,258.15,0.)
    assert solution['mass_rate_kg_s']==solution['alpha_a']==solution['alpha_c']==0.
    assert model.axis_rates(solution,20e-6,20e-6,919.)['equatorial_rate_m_s']==0.


@pytest.mark.parametrize('value', [-1.,float('nan'),float('inf')])
def test_invalid_forcing_rejected(value):
    with pytest.raises(ValueError): model.solve(state(),20e-6,20e-6,258.15,value)


def test_sphere_large_crystal_diffusion_limit():
    source=state(1e-3)
    source['characteristic_supersaturation_a']=source['characteristic_supersaturation_c']=1e-10
    solution=model.solve(source,1e-3,1e-3,258.15,.02)
    assert solution['alpha_a']==pytest.approx(1.)
    # Independently compute the scalar series resistance including both jumps.
    r=1e-3;ratio=r/source['boundary_capacitance_m']
    d=1/(ratio/2e-5+4/(550*r));k=1/(ratio/.023+source['thermal_jump_length_m']/(.023*r))
    l=source['latent_sublimation_J_kg'];rho=165/(461.5*258.15)
    expected=4*math.pi*r*rho*.02/(1/d+rho*l/(k*258.15)*(l/(461.5*258.15)-1))
    assert solution['mass_rate_kg_s']==pytest.approx(expected,rel=2e-11)


def test_perfect_attachment_bracket_has_exact_endpoint_sign():
    for length in (1e-8,.1,1.,143.04,1e8):
        ss,aa=model.surface_state(.1,length,1e-30,10.)
        assert aa==1.
        assert ss==pytest.approx(.1/(1+length),rel=1e-13)
