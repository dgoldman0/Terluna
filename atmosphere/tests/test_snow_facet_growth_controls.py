"""Independent spherical growth solution for the trajectory-control integrator."""
import math
import pytest
from atmosphere.snow import facet_growth_controls as controls
from atmosphere.tests.test_snow_facet_growth import sphere_state


class ContinuumSphere:
    """Idealized transport state: unity attachment and negligible molecular jumps."""
    def evaluate(self,si,t,p,a,c,m):
        assert a==pytest.approx(c,rel=1e-11)
        source=sphere_state(a,alpha=1.)
        source.update(boundary_capacitance_m=a,vapor_molecular_speed_m_s=1e99,
            thermal_jump_length_m=1e-100,characteristic_supersaturation_a=1e-100,
            characteristic_supersaturation_c=1e-100)
        return source


@pytest.mark.parametrize('method,factor', [('semianalytic',1.),('implicit',1.),('implicit',2.)])
def test_integrated_radius_matches_diffusion_heat_solution(method,factor):
    r0=10e-6;t=258.15;si=.01;latent=2.83e6;rhoeq=165/(461.5*t)
    diffusion=1/(1/2e-5+rhoeq*latent/(.023*t)*(latent/(461.5*t)-1))
    expected=math.sqrt(r0*r0+2*rhoeq*si*diffusion*factor/controls.DENSITY*600)
    result=controls.integrate(ContinuumSphere(),t,101000.,si,r0,r0,10.,method,rtol=1e-9,ventilation=factor)
    assert result['complete']
    assert result['endpoint']['equatorial_dimension_m']==pytest.approx(2*expected,rel=2e-8)
    assert result['endpoint']['polar_dimension_m']==pytest.approx(2*expected,rel=2e-8)
    assert result['maximum_sampled_mass_volume_residual']<2e-8


def test_semianalytic_cannot_silently_ignore_ventilation():
    with pytest.raises(ValueError,match='no ventilation input'):
        controls.integrate(ContinuumSphere(),258.15,101000.,.01,10e-6,10e-6,10.,'semianalytic',ventilation=2.)


def test_growth_control_rejects_unknown_method():
    with pytest.raises(ValueError,match='Unknown attachment method'):
        controls.integrate(ContinuumSphere(),258.15,101000.,.01,10e-6,10e-6,10.,'unspecified')
