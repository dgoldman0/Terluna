import math
import pytest
from shared.constants import BOLTZMANN
from atmosphere.snow.step_humidity_audit import infer_under_depletion, mean_attachment_for_depletion


@pytest.mark.parametrize('s,f',[(.001,.0027),(.05,.0027),(.1,.5),(.01,1.)])
def test_inverse_preserves_measured_velocity_and_spiral_spacing(s,f):
    beta,kappa,T,area=2e-5,8.5e-10,258.,8.3e-20
    got=infer_under_depletion(beta,kappa,s,f)
    assert got['inferred_beta_m_s']*got['actual_surface_sigma']==pytest.approx(beta*s)
    before=19*area*kappa/(BOLTZMANN*T*math.log1p(s))
    after=19*area*got['inferred_energy_J_m']/(BOLTZMANN*T*math.log1p(got['actual_surface_sigma']))
    assert after==pytest.approx(before,rel=1e-14)


def test_planar_humidity_matches_independent_diffusion_and_attachment_flux():
    X,H,f=1.4e-7,.002,.0027
    alpha=mean_attachment_for_depletion(X,H,f)
    # Normalize by the kinetic flux at unit reservoir supersaturation.
    assert X/H*(1-f)==pytest.approx(alpha*f,rel=1e-14)
    assert mean_attachment_for_depletion(X,H,1.)==0


@pytest.mark.parametrize('fraction',[0.,1.1,float('nan')])
def test_inverse_rejects_undefined_or_enhanced_humidity(fraction):
    with pytest.raises(ValueError):infer_under_depletion(1e-5,1e-10,.01,fraction)


def test_required_path_obeys_planar_ficks_law_for_the_observed_growth():
    from atmosphere.snow.step_humidity_audit import diffusion_path_for_growth
    # rho_eq/rho_ice, D and molecular speed are independently supplied.
    density_ratio,D,vth=2e-6,2e-5,140.
    reservoir,local_fraction,normal=.01,.003,3e-11
    path=diffusion_path_for_growth(D/vth,density_ratio*vth,reservoir,normal,local_fraction)
    recovered=density_ratio*D*(reservoir-local_fraction*reservoir)/path
    assert recovered==pytest.approx(normal,rel=1e-14)
