"""Independent moment, distribution, exit and scaling controls for capture jumps."""
import math
import numpy as np
import pytest
from scipy.stats import poisson
from atmosphere.snow import collision_entropy as entropy
from atmosphere.snow import collision_jumps as jumps


@pytest.mark.parametrize('mean',[0.,.01,8.,100.,30000.])
def test_poisson_rule_preserves_probabilities_without_renormalization(mean):
    rule=jumps.poisson_rule(mean)
    k,w=rule.counts,rule.weights
    assert (w>=0).all()
    assert abs(rule.normalization_residual)<3e-14
    assert abs(w.sum()-1)<3e-14+rule.omitted_count_moments[0]
    assert abs(w@k-mean)<3e-12*max(mean,1)+rule.omitted_count_moments[1]
    assert abs(w@(k-mean)**2-mean)<1e-9*max(mean,1)
    if mean:
        indices=np.unique(np.linspace(0,len(k)-1,9).astype(int))
        np.testing.assert_allclose(w[indices],poisson.pmf(k[indices],mean),rtol=2e-10,atol=0)


def test_large_macrostep_has_physical_poisson_variance_and_number_cdf():
    grid=entropy.MomentGrid(np.geomspace(.1,10,31))
    initial=grid.from_particles([1.],[1.]);d=3e-6;rate=100.;dt=200.
    op=jumps.ConstantCapture(grid,d,rate)
    result=jumps.integrate(op,initial,[0,dt],max_step_s=dt)
    row=result['snapshots'][-1];mean=1+d*rate*dt;variance=d*d*rate*dt
    assert row['moments'][0]==pytest.approx(1.,abs=1e-13)
    assert row['moments'][1]==pytest.approx(mean,rel=1e-13)
    assert row['variance']==pytest.approx(variance,rel=1e-8)
    k=np.arange(int(rate*dt-10*np.sqrt(rate*dt)),int(rate*dt+10*np.sqrt(rate*dt)))
    query=np.sort(np.r_[1+d*k,np.nextafter(1+d*k,-np.inf)])
    exact=np.empty_like(query);exact[::2]=poisson.cdf(k-1,rate*dt);exact[1::2]=poisson.cdf(k,rate*dt)
    assert np.max(abs(grid.cdf(row['state'],query)-exact))<.005
    assert abs(row['mass_residual'])<1e-13
    assert result['diagnostics']['macrosteps']==1
    assert result['diagnostics']['maximum_poisson_mean']==20000


def test_continuous_source_moments_cross_faces_and_refine_quadrature():
    grid=entropy.MomentGrid([.8,1.,1.2,1.5,2.,3.])
    state=np.zeros((6,3));mean=.5;variance=.015;u=mean-grid.offsets[1]
    state[1]=[1,u,u*u+variance]
    initial=grid.moments(state);d=.004;mu=100.
    rows=[]
    for order in (8,16,32):
        op=jumps.ConstantCapture(grid,d,10.,quadrature_order=order)
        output,ledger,bound,storage=op.step(state,10.)
        m=grid.moments(output)
        expected=[initial[0],initial[1]+d*mu,initial[2]+2*d*mu*initial[1]+d*d*(mu*mu+mu)]
        np.testing.assert_allclose(m,expected,rtol=2e-11,atol=1e-14)
        assert np.count_nonzero(output[:,0])>1
        assert ledger['supplied_mass']==pytest.approx(d*mu,rel=2e-11)
        rows.append(output)
    np.testing.assert_allclose(rows[-1],rows[-2],rtol=2e-9,atol=1e-13)


def test_atomic_exit_stops_at_first_crossing_instead_of_final_poisson_mass():
    grid=entropy.MomentGrid([.8,1.,1.2])
    m0=1.;d=.1;mu=8.;upper=grid.bounds[-1]
    count=jumps.first_exit_count(m0,d,upper)
    op=jumps.ConstantCapture(grid,d,4.)
    row=jumps.integrate(op,grid.from_particles([m0],[1.]),[2.],max_step_s=2.)['snapshots'][0]
    probability=poisson.sf(count-1,mu);mass=m0+count*d
    np.testing.assert_allclose(row['ledgers']['upper_exit'],probability*np.array([1,mass,mass*mass]),rtol=1e-14)
    k=np.arange(count);expected_input=d*(poisson.pmf(k,mu)@k+count*probability)
    assert row['ledgers']['supplied_mass']==pytest.approx(expected_input,rel=1e-14)
    assert row['ledgers']['supplied_mass']<d*mu
    assert abs(row['number_residual'])<1e-13 and abs(row['mass_residual'])<1e-13


def test_continuous_first_exit_has_uniform_overshoot_in_this_exact_case():
    grid=entropy.MomentGrid([1.,2.])
    state=np.zeros((2,3));u=.5-grid.offsets[1];state[1]=[1,u,u*u+1/12]
    op=jumps.ConstantCapture(grid,.25,4.)
    output,ledger,_,_=op.step(state,2.)
    probability=poisson.sf(np.arange(6),8.).sum()/6
    expected=probability*np.array([1.,3.125,9+.75+.25**2/3])
    np.testing.assert_allclose(ledger['upper_exit'],expected,rtol=2e-12)
    assert grid.moments(output)[0]+probability==pytest.approx(1.,abs=2e-12)
    assert grid.moments(output)[1]+ledger['upper_exit'][1]-ledger['supplied_mass']==pytest.approx(2.25,abs=2e-12)


@pytest.mark.parametrize('population',[1e-100,1e-300])
def test_population_scaling_does_not_restore_an_event_timescale(population):
    grid=entropy.MomentGrid(np.geomspace(.5,2,12))
    initial=grid.from_particles([1.],[1.])
    rows=[]
    for factor in (1.,population):
        op=jumps.ConstantCapture(grid,.0005,100.)
        result=jumps.integrate(op,initial*factor,[1.,2.],max_step_s=1.)
        rows.append(result['snapshots'][-1])
        assert result['diagnostics']['macrosteps']==2
    np.testing.assert_allclose(rows[1]['state']/population,rows[0]['state'],rtol=2e-9,atol=2e-12)
    assert abs(rows[1]['mass_residual']/population)<2e-12


def test_zero_rate_is_exact_identity():
    grid=entropy.MomentGrid([1.,2.]);initial=grid.from_particles([1.],[1e-300])
    result=jumps.integrate(jumps.ConstantCapture(grid,.1,0),initial,[3.])
    np.testing.assert_array_equal(result['snapshots'][0]['state'],initial)


def test_tail_bounds_cover_probability_and_both_unbounded_count_moments():
    for mean in (.001,5.,100.,30000.):
        rule=jumps.poisson_rule(mean,tail_tolerance=1e-12)
        lo,hi=rule.counts[[0,-1]]
        q0=poisson.cdf(lo-1,mean)+poisson.sf(hi,mean)
        q1=mean*(poisson.cdf(lo-2,mean)+poisson.sf(hi-1,mean))
        q2=mean**2*(poisson.cdf(lo-3,mean)+poisson.sf(hi-2,mean))+q1
        assert np.all(rule.omitted_count_moments >= np.array([q0,q1,q2]))
        assert rule.omitted_count_moments[0]<=1e-12
    assert jumps.poisson_tail_bound(1000,.001,upper=True)>0


def test_atomic_cdf_uses_physical_mass_on_shifted_grid():
    # This grid exposes cancellation in (query-bound)/width at the exact atom.
    from atmosphere.snow import collision_entropy_benchmarks as reference
    import json
    product = json.loads((reference.ROOT/reference.control.REFERENCE).read_text())
    unit = product['capture']['initial_collector_mass_kg']
    grid = entropy.MomentGrid(reference.reference.nested_grid(2)*(1.02/unit))
    state = grid.from_particles([1.], [0.4])
    query = np.array([np.nextafter(1.,-np.inf),1.,np.nextafter(1.,np.inf)])
    np.testing.assert_array_equal(jumps.number_cdf(grid,state,query), [0.,.4,.4])
