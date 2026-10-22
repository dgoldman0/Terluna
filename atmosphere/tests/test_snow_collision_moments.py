"""Independent controls for the bounded local collision-moment closure."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import betaln

from atmosphere.snow import collision_moments as m


def test_full_beta_quadrature_retains_moments_at_narrow_and_broad_limits():
    for mean, variance in ((.4, .0004), (.4, 1e-14), (.4, .2), (.01, .009)):
        measure = m.BetaMeasure(mean, variance)
        for order in (2, 4):
            x, w = measure.quadrature(order)
            assert (w > 0).all() and (x >= 0).all() and (x <= 1).all()
            assert w.sum() == pytest.approx(1., abs=2e-14)
            assert w@x == pytest.approx(mean, abs=2e-14)
            assert w@((x-mean)**2) == pytest.approx(variance, rel=2e-9)


def test_atom_limits_and_adjacent_segments_do_not_duplicate_boundary_atoms():
    for measure in (m.BetaMeasure(.5, 0.), m.BetaMeasure(.4, .24)):
        nodes, weights = measure.quadrature(4)
        assert weights.sum() == 1
        assert weights@nodes == pytest.approx(measure.mean)
        parts = [measure.quadrature(4, a, b) for a,b in ((0.,.5),(.5,1.))]
        assert sum(w.sum() for _,w in parts) == 1
        assert measure.cdf(-1) == 0 and measure.cdf(2) == 1


def test_narrow_truncated_beta_moments_against_independent_smooth_integration():
    # This interval catastrophically cancels degree4/6 raw recurrence moments.
    a, b = 200., 300.
    mean, variance = a/(a+b), a*b/((a+b)**2*(a+b+1))
    sd = math.sqrt(variance)
    lower, upper = mean-.001*sd, mean+.001*sd
    measure = m.BetaMeasure(mean, variance)
    reference = []
    scale = .001*sd
    for degree in range(8):
        value = quad(lambda z: z**degree*math.exp((a-1)*math.log(mean+scale*z)
                     +(b-1)*math.log1p(-mean-scale*z)-betaln(a,b))*scale,
                     -1, 1, epsabs=1e-15)[0]
        reference.append(value)
    for order in (2,4):
        x, w = measure.quadrature(order, lower, upper)
        assert x.min() >= lower and x.max() <= upper
        for degree in range(2*order):
            observed = w@((x-mean)/scale)**degree
            assert observed == pytest.approx(reference[degree], abs=2e-13, rel=2e-8)


def test_exact_particle_injection_and_reconstruction_cdf():
    grid = m.MomentGrid([1., 2., 4., 8.])
    state = grid.from_particles([1., 1.1, 4.5], [2., 3., 5.])
    np.testing.assert_allclose(grid.moments(state), [10., 27.8, 106.88], rtol=1e-14)
    for weighting in ('number','mass'):
        cdf = grid.cdf(state, np.linspace(0,20,100), weighting)
        assert cdf[0] == 0 and cdf[-1] == pytest.approx(1)
        assert (np.diff(cdf) >= -1e-14).all()
    exact = grid.from_particles([2.], [1.])
    np.testing.assert_array_equal(grid.cdf(exact, [1.,2.,3.]), [0.,1.,1.])


def test_invalid_state_and_kernel_are_rejected_without_population_clipping():
    grid = m.MomentGrid([1.,2.,4.])
    initial = grid.from_particles([2.],[1.])
    invalid = []
    x = initial.copy(); x[1,0] = -.01; invalid.append(x)
    x = initial.copy(); x[1,1] = 2.; invalid.append(x)
    x = initial.copy(); x[1,2] = -1e-3; invalid.append(x)
    x = initial.copy(); x[0,2] = .001; invalid.append(x)
    x = initial.copy(); x[1,2] = np.nan; invalid.append(x)
    for x in invalid:
        with pytest.raises(m.RealizabilityError):
            grid.distributions(x)
    for kernel in (lambda a,b: -1., lambda a,b: np.nan, lambda a,b: a+2*b):
        with pytest.raises(ValueError):
            m.BathCapture(grid,kernel,.01,1.).evaluate(initial)


def test_bath_constant_rate_has_exact_poisson_moments_and_water_ledger():
    grid = m.MomentGrid(np.geomspace(.01,100,24))
    initial = grid.from_particles([1.], [2.])
    result = m.integrate(m.BathCapture(grid,lambda a,b: 1.,.001,3.), initial,
                         [0.,.1,1.], max_step_s=.04)
    for row in result['snapshots']:
        t = row['time_s']
        np.testing.assert_allclose(row['moments'][:2], [2.,2*(1+.003*t)], rtol=1e-12)
        assert row['variance'] == pytest.approx(3e-6*t, abs=2e-15)
        assert row['ledgers']['supplied_mass'] == pytest.approx(.006*t, abs=1e-14)
        assert abs(row['number_residual']) < 1e-12 and abs(row['mass_residual']) < 1e-12


def test_unrealizable_large_trial_is_rejected_and_successful_state_keeps_variance():
    grid = m.MomentGrid(np.geomspace(.1,100,24))
    initial = grid.from_particles([grid.pivots[8]], [1.])
    result = m.integrate(m.BathCapture(grid,lambda a,b: 1.,1e-4,100.), initial,
                         [1.], max_step_s=1.)
    assert result['diagnostics']['rejected_steps'] > 0
    assert result['snapshots'][0]['variance'] == pytest.approx(1e-6, rel=1e-8)
    assert np.min(result['snapshots'][0]['state'][:,0]) >= 0


def test_collision_reaction_moments_and_explicit_domain_exit():
    grid = m.MomentGrid([1.,2.,4.,8.])
    state = grid.from_particles([2.,8.],[3.,5.])
    operator = m.CollisionMoments(grid, lambda a,b: 1.)
    rhs, ledger = operator.evaluate(state)
    derivative = grid.moments(rhs)
    initial = grid.moments(state)
    np.testing.assert_allclose(derivative+ledger['upper_exit'],
                              [-.5*initial[0]**2, 0., initial[1]**2], atol=1e-12)
    assert ledger['collision_number_loss'] == pytest.approx(.5*initial[0]**2)
    assert ledger['upper_exit'][0] > 0
    assert (ledger['lower_exit'] == 0).all()


@pytest.mark.parametrize('kind', ['constant','additive'])
def test_constant_and_additive_kernel_analytic_moments(kind):
    grid = m.MomentGrid(np.geomspace(.25,2048,28))
    state = grid.from_particles([1.,2.],[.7,.3])
    initial = grid.moments(state)
    kernel = (lambda a,b: 1.) if kind == 'constant' else (lambda a,b: a+b)
    result = m.integrate(m.CollisionMoments(grid,kernel),state,[.01],max_step_s=.002)
    row = result['snapshots'][0]
    if kind == 'constant':
        expected = [initial[0]/(1+.5*initial[0]*.01), initial[1],initial[2]+initial[1]**2*.01]
    else:
        expected = [initial[0]*math.exp(-initial[1]*.01),initial[1],initial[2]*math.exp(2*initial[1]*.01)]
    np.testing.assert_allclose(row['moments'], expected, rtol=2e-6, atol=1e-12)
    assert abs(row['mass_residual']) < 1e-11


def test_endpoint_singularity_and_atomic_rounding_regressions():
    d = m.BetaMeasure(.9968295730174938,1.2462046563496227e-5)
    p = d._probability(.9999743119551164,1.)
    for order in (2,4):
        x,w = d.quadrature(order,.9999743119551164,1.)
        assert w.sum() == pytest.approx(p,rel=1e-11)
        assert x.min() >= .9999743119551164 and x.max() <= 1.
    grid = m.MomentGrid([1.,2.,4.,8.])
    state = np.zeros((4,3))
    u = np.array([-grid.offsets[1],1-grid.offsets[1]])
    state[1] = [1.,np.array([.98,.02])@u,np.array([.98,.02])@u**2]
    measure = grid.distributions(state)[1]
    assert measure.atoms is not None
    assert np.isfinite(measure.cdf(.5))


def test_skew_beta_small_segment_uses_complementary_probability():
    a = 1e-12
    mean = a/(a+1)
    variance = a/((a+1)**2*(a+2))
    d = m.BetaMeasure(mean,variance)
    lower,upper = .4*mean,.6*mean
    # Beta(a,1) has exact CDF x**a; expm1 avoids cancellation here.
    expected = lower**a*math.expm1(a*math.log(upper/lower))
    assert d._probability(lower,upper) == pytest.approx(expected,rel=1e-12)


@pytest.mark.parametrize('order',[2,4])
def test_segmented_continuous_parent_collision_generator_has_exact_polynomial_moments(order):
    grid = m.MomentGrid([1.,2.,4.,8.])
    state = np.zeros((4,3))
    for i,number,mean,variance in ((0,.4,.3,.03),(1,.6,.6,.02)):
        local = mean-grid.offsets[i]
        state[i] = number*np.array([1.,local,local**2+variance])
    number,mass,second = grid.moments(state)
    for kernel,expected in ((lambda a,b:1.,[-.5*number**2,0.,mass**2]),
                            (lambda a,b:a+b,[-number*mass,0.,2*mass*second])):
        rhs,ledger = m.CollisionMoments(grid,kernel,order).evaluate(state)
        observed = grid.moments(rhs)+ledger['lower_exit']+ledger['upper_exit']
        np.testing.assert_allclose(observed,expected,rtol=2e-8,atol=2e-10)


def test_progress_hook_and_zero_kernel_leave_population_exactly_unchanged():
    grid = m.MomentGrid([1.,2.,4.])
    state = grid.from_particles([1.,2.],[3.,4.])
    calls=[]
    result = m.integrate(m.CollisionMoments(grid,lambda a,b:0.),state,[1.],max_step_s=.25,
                         progress=lambda t,d:calls.append((t,d['steps'])),time_limit_s=10.)
    np.testing.assert_array_equal(result['snapshots'][0]['state'],state)
    assert calls == [(.25,1),(.5,2),(.75,3),(1.,4)]
    with pytest.raises(ValueError,match='wall time'):
        m.integrate(m.CollisionMoments(grid,lambda a,b:0.),state,[1.],time_limit_s=0.)
