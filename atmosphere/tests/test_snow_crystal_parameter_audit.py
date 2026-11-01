import math
import numpy as np
import pytest
from shared.constants import BOLTZMANN
from atmosphere.snow.crystal_parameter_audit import mapped_state, MAPPINGS, nucleation_barrier


def test_equal_area_mapping_matches_polygon_area_and_growth_independently():
    a,c,da,dc=40e-6,200e-6,3e-8,2e-7
    theta=np.arange(6)*math.pi/3
    def volume(a,c):
        x,y=a*np.cos(theta),a*np.sin(theta)
        return abs(np.dot(x,np.roll(y,-1))-np.dot(y,np.roll(x,-1)))*c
    f=MAPPINGS['equal_basal_area']
    assert 2*math.pi*(f*a)**2*c==pytest.approx(volume(a,c),rel=1e-14)
    dt=1e-3
    derivative=(volume(a+da*dt,c+dc*dt)-volume(a-da*dt,c-dc*dt))/(2*dt)
    mapped_derivative=2*math.pi*(2*(f*a)*(f*da)*c+(f*a)**2*dc)
    assert mapped_derivative==pytest.approx(derivative,rel=1e-8)


@pytest.mark.parametrize('convention',list(MAPPINGS))
def test_mapping_keeps_forcing_and_mass_target_and_converts_radial_derivative(convention):
    original=dict(radius_m=4.,half_length_m=20.,fitted_rates_a_c_mass=[.3,.8,1.2],
                  ventilation=1.7,air={'diffusivity':2.})
    mapped=mapped_state(original,convention)
    assert mapped['radius_m']/original['radius_m']==pytest.approx(mapped['fitted_rates_a_c_mass'][0]/.3)
    assert mapped['fitted_rates_a_c_mass'][1:]==[.8,1.2]
    assert mapped['half_length_m']==20.
    assert mapped['ventilation']==1.7
    mapped['air']['diffusivity']=100.
    assert original['air']=={'diffusivity':2.}


@pytest.mark.parametrize('T,kappa,area',[(250.,1e-12,8.3e-20),(268.,8.5e-10,8.3e-20)])
def test_barrier_matches_independent_critical_nucleus_free_energy(T,kappa,area):
    s=1e-5
    mu=BOLTZMANN*T*math.log1p(s)
    radius=kappa*area/mu
    free_energy=lambda r:2*math.pi*r*kappa-math.pi*r*r*mu/area
    barrier=free_energy(radius)/(BOLTZMANN*T)
    assert free_energy(radius)>free_energy(radius*.99)
    assert free_energy(radius)>free_energy(radius*1.01)
    assert nucleation_barrier(T,kappa,area,polynuclear=False)==pytest.approx(barrier*math.log1p(s),rel=1e-13)
    assert nucleation_barrier(T,kappa,area)==pytest.approx(barrier*math.log1p(s)/3,rel=1e-13)


@pytest.mark.parametrize('T,kappa,area',[(0.,1e-12,8e-20),(260.,-1e-12,8e-20),(260.,1e-12,float('nan'))])
def test_energy_comparison_rejects_nonphysical_inputs(T,kappa,area):
    with pytest.raises(ValueError):nucleation_barrier(T,kappa,area)
