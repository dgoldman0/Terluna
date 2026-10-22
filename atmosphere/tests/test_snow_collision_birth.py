"""Check finite-time variable capture against analytic and matrix references."""
import numpy as np
import pytest
from scipy.linalg import expm
from scipy.stats import poisson
from atmosphere.snow import collision_entropy as entropy
from atmosphere.snow import collision_birth as birth
from atmosphere.snow import collision_jumps as jumps


def test_variable_rate_birth_and_first_exit_match_dense_matrix_exponential():
    grid=entropy.MomentGrid([.5,.8,1.2,1.7])
    initial_mass=.8;donor=.2;dt=.7
    kernel=lambda m,d: .3+np.asarray(m)**2
    count=jumps.first_exit_count(initial_mass,donor,grid.bounds[-1])
    mass=initial_mass+donor*np.arange(count+1)
    rates=kernel(mass[:-1],donor)
    Q=np.diag(np.r_[-rates,0.])+np.diag(rates,-1)
    exact=expm(dt*Q)[:,0]
    op=birth.BirthCapture(grid,kernel,donor,1.)
    output,ledger,tail,cap=op.step(grid.from_particles([initial_mass],[1.]),dt)
    values=grid.moments(output)
    np.testing.assert_allclose(values,[exact[:-1]@mass[:-1]**k for k in range(3)],rtol=2e-13,atol=1e-14)
    np.testing.assert_allclose(ledger['upper_exit'],exact[-1]*np.array([1,mass[-1],mass[-1]**2]),rtol=2e-13)
    assert ledger['supplied_mass']==pytest.approx(exact@(mass-initial_mass),rel=2e-13)
    assert np.max(cap)==0


def test_constant_rate_uniformization_matches_poisson_moments_at_large_step():
    grid=entropy.MomentGrid(np.geomspace(.5,2.,20));initial=grid.from_particles([1.],[1.])
    op=birth.BirthCapture(grid,lambda m,d:np.ones_like(m),.001,50.)
    output,ledger,_,_=op.step(initial,1.)
    np.testing.assert_allclose(grid.moments(output),[1.,1.05,1.05**2+.00005],rtol=1e-13)
    assert grid.variance(output)==pytest.approx(.00005,rel=1e-10)
    assert ledger['supplied_mass']==pytest.approx(.05,rel=1e-13)


def test_linear_birth_rate_matches_negative_binomial_mean_and_variance():
    # lambda(k)=a*(k+r): negative-binomial birth increments from k=0.
    a=.2;r=3.;donor=.1;m0=1.;dt=1.
    grid=entropy.MomentGrid(np.geomspace(.5,100,40))
    kernel=lambda m,d:a*((np.asarray(m)-m0)/donor+r)
    op=birth.BirthCapture(grid,kernel,donor,1.)
    output,ledger,_,_=op.step(grid.from_particles([m0],[1.]),dt)
    expected_mean=m0+donor*r*np.expm1(a*dt)
    expected_variance=donor**2*r*np.expm1(a*dt)*np.exp(a*dt)
    assert grid.moments(output)[1]==pytest.approx(expected_mean,rel=2e-13)
    assert grid.variance(output)==pytest.approx(expected_variance,rel=2e-10)


@pytest.mark.parametrize('scale',[1e-100,1e-300])
def test_variable_rate_scaling_has_no_population_cutoff(scale):
    grid=entropy.MomentGrid(np.geomspace(.5,2.,12));initial=grid.from_particles([1.],[1.])
    kernel=lambda m,d: np.asarray(m)**2
    answers=[]
    for factor in (1.,scale):
        op=birth.BirthCapture(grid,kernel,.0005,30.)
        row=birth.integrate(op,initial*factor,[1.,2.],max_step_s=.5)['snapshots'][-1]
        answers.append(row)
    np.testing.assert_allclose(answers[1]['state']/scale,answers[0]['state'],rtol=5e-9,atol=2e-12)
    assert abs(answers[1]['mass_residual']/scale)<2e-12


def test_work_budget_rejects_before_admitting_an_incomplete_transition():
    grid=entropy.MomentGrid([1.,2.]);op=birth.BirthCapture(grid,lambda m,d:1.,.001,100.)
    with pytest.raises(birth.StepBudgetError):op.step(grid.from_particles([1.],[1.]),1.)


def test_continuous_sources_and_first_exit_against_uniform_overshoot_identity():
    grid=entropy.MomentGrid([1.,2.])
    state=np.zeros((2,3));u=.5-grid.offsets[1];state[1]=[1.,u,u*u+1/12]
    op=birth.BirthCapture(grid,lambda m,d:np.ones_like(m),.25,4.)
    output,ledger,_,_=op.step(state,2.)
    exit_probability=poisson.sf(np.arange(6),8.).sum()/6
    expected=exit_probability*np.array([1.,3.125,9.+.75+.25**2/3])
    np.testing.assert_allclose(ledger['upper_exit'],expected,rtol=3e-12)
    assert grid.moments(output)[0]+exit_probability==pytest.approx(1.,abs=3e-12)
    assert grid.moments(output)[1]+ledger['upper_exit'][1]-ledger['supplied_mass']==pytest.approx(2.25,abs=3e-12)
