"""Analytical controls for bounded quadratic maximum entropy."""
import math
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtr

from atmosphere.snow import collision_entropy as e


def test_untruncated_gaussian_limit_and_quadrature():
    d = e.EntropyMeasure(.4,.0004)
    assert d.A == pytest.approx(0,abs=1e-10)
    assert d.B == pytest.approx(-.5,abs=1e-10)
    x = np.linspace(.32,.48,33)
    np.testing.assert_allclose(d.cdf(x),ndtr((x-.4)/.02),atol=1e-11)
    for order in (2,4):
        nodes,weights = d.quadrature(order)
        assert weights.sum() == pytest.approx(1.,abs=1e-11)
        assert weights@nodes == pytest.approx(.4,abs=1e-11)
        assert weights@(nodes-.4)**2 == pytest.approx(.0004,rel=1e-10)


def test_half_gaussian_at_cell_boundary_is_representable_without_beta_shape():
    original_sd = .01
    mean = 1-original_sd*math.sqrt(2/math.pi)
    variance = original_sd**2*(1-2/math.pi)
    d = e.EntropyMeasure(mean,variance)
    x = np.linspace(.95,1.,101)
    expected = 2*ndtr((x-1)/original_sd)
    np.testing.assert_allclose(d.cdf(x),expected,atol=2e-10)
    for order in (2,4):
        nodes,weights = d.quadrature(order,.995,1.)
        assert weights.sum() == pytest.approx(1-2*ndtr(-.5),rel=2e-10)
        assert nodes.min() >= .995 and nodes.max() <= 1.


def test_positive_curvature_and_uniform_limits():
    for mean,variance in ((.4,.2),(.5,1/12),(.1,.04),(.01,.009)):
        d = e.EntropyMeasure(mean,variance)
        if variance > 1/12:
            assert d.B > 0
        for order in (2,4):
            nodes,w = d.quadrature(order)
            assert (w > 0).all() and (nodes >= 0).all() and (nodes <= 1).all()
            assert w.sum() == pytest.approx(1.,abs=2e-10)
            assert w@nodes == pytest.approx(mean,abs=2e-10)
            assert w@(nodes-mean)**2 == pytest.approx(variance,abs=2e-10)
        assert np.all(np.diff(d.cdf(np.linspace(0,1,51))) >= -1e-12)
    d = e.EntropyMeasure(.5,1/12)
    np.testing.assert_allclose(d.cdf([.1,.5,.9]),[.1,.5,.9],atol=1e-10)


def test_atomic_limits_and_exact_moment_event_bookkeeping():
    for mean,var in ((.3,0.),(.3,.21)):
        d = e.EntropyMeasure(mean,var)
        x,w = d.quadrature(4)
        assert w.sum() == 1 and w@x == pytest.approx(mean)
    grid = e.MomentGrid([1.,2.,4.,8.])
    state = grid.from_particles([2.,8.],[3.,5.])
    rhs,ledger = e.CollisionMoments(grid,lambda a,b:1.).evaluate(state)
    n,mass,_=grid.moments(state)
    np.testing.assert_allclose(grid.moments(rhs)+ledger['upper_exit'],[-n*n/2,0,mass*mass],atol=1e-12)


def test_constant_capture_preserves_physical_variance_and_water():
    grid = e.MomentGrid(np.geomspace(.01,100,24))
    state = grid.from_particles([1.],[1.])
    result = e.integrate(e.BathCapture(grid,lambda a,b:1.,.001,3.),state,[1.],max_step_s=.05,time_limit_s=10.)
    row = result['snapshots'][0]
    np.testing.assert_allclose(row['moments'][:2],[1.,1.003],rtol=1e-12)
    assert row['variance'] == pytest.approx(3e-6,rel=1e-9)
    assert abs(row['mass_residual']) < 1e-12


@pytest.mark.parametrize('linear,quadratic',[(-8.,-3.),(4.,-12.),(-20.,20.),(2.,0.),(0.,40.),(-10.,-.001),(-10.,.001)])
def test_inverse_density_and_cumulative_mass_against_independent_integrals(linear,quadratic):
    peak=max(0.,linear+quadratic)
    if quadratic < 0 and 0 < -linear/(2*quadratic) < 1:
        peak=max(peak,-linear**2/(4*quadratic))
    density=lambda x:math.exp(linear*x+quadratic*x*x-peak)
    moments=np.array([quad(lambda x:density(x)*x**r,0.,1.,epsabs=1e-13)[0] for r in range(3)])
    mean=moments[1]/moments[0]
    variance=moments[2]/moments[0]-mean**2
    d=e.EntropyMeasure(mean,variance)
    assert d.A == pytest.approx((linear+2*quadratic*mean)*math.sqrt(variance),abs=2e-9)
    assert d.B == pytest.approx(quadratic*variance,abs=2e-9)
    queries=np.array([.01,.1,.3,.7,.99])
    expected_number=np.array([quad(density,0.,x,epsabs=1e-13)[0]/moments[0] for x in queries])
    expected_mass=np.array([quad(lambda y:y*density(y),0.,x,epsabs=1e-13)[0]/moments[0] for x in queries])
    np.testing.assert_allclose(d.cdf(queries),expected_number,atol=3e-10)
    np.testing.assert_allclose(d.partial_first(queries),expected_mass,atol=3e-10)


def test_subnormal_storage_uncertainty_is_bounded_and_does_not_modify_state():
    # Exact last rejected state from the independently instrumented collapse.
    g=e.MomentGrid([.5,1.,1.9272963557631262,4.])
    # The archived cell's geometry fixes its centered coordinate convention.
    g.offsets[2]=.4672410861352382
    g.widths[2]=.2536332838613804
    state=np.zeros((4,3));state[2]=[4.74e-322,-2.17e-322,1e-322]
    original=state.copy()
    reconstructed=g.distributions(state)[2]
    np.testing.assert_array_equal(state,original)
    assert reconstructed.variance == 0.
    assert reconstructed.atoms is not None
    d=g.roundoff_diagnostics
    assert d['storage_reconstructed_cells'] == 1
    assert d['maximum_affected_number'] == state[2,0]
    assert d['maximum_absolute_reconstruction_second_moment_change'] <= np.nextafter(0.,1.)
    assert d['maximum_storage_number_rounding_bound'] == np.nextafter(0.,1.)
    # The first rejected subnormal trial had a real deficit larger than the
    # representation bound and must still trigger timestep rejection.
    state[2]=[4.0418849e-316,-1.8839404e-316,8.7810603e-317]
    with pytest.raises(e.RealizabilityError,match='beyond storage uncertainty'):
        g.distributions(state)


def test_reconstruction_does_not_erase_feasible_subnormal_variance():
    g=e.MomentGrid([1.,2.,4.]);state=np.zeros((3,3))
    state[1]=[1e-320,0.,1e-322]
    d=g.distributions(state)[1]
    assert d.variance > 0
    assert g.roundoff_diagnostics['storage_reconstructed_cells'] == 0


def test_stable_ssp_average_preserves_an_admissible_subnormal_stationary_state():
    # Exact final tuple from the second independently instrumented failure.
    grid=e.MomentGrid([.8,1.1402227451174085,1.3,1.5])
    state=np.zeros((4,3));state[2]=[2.03e-322,-9.4e-323,4.4e-323]
    # Match the archived geometric-grid centered coordinate.
    grid.offsets[2]=.4672410861352388
    grid.distributions(state)
    assert np.any(.5*state+.5*state != state)
    result=e.integrate(e.BathCapture(grid,lambda a,b:0.,1e-6,100.),state,[1.],max_step_s=.5)
    np.testing.assert_array_equal(result['snapshots'][0]['state'],state)
    assert result['diagnostics']['steps'] == 2
    assert result['diagnostics']['rejected_steps'] == 0


def test_new_ssp_sequence_matches_frozen_normal_range_moments_and_ledgers():
    from atmosphere.snow import collision_moments as frozen
    grid=e.MomentGrid(np.geomspace(.01,100,24))
    state=grid.from_particles([1.],[1.])
    old=frozen.integrate(e.BathCapture(grid,lambda a,b:1.,.001,3.),state,[1.],max_step_s=.05)
    new=e.integrate(e.BathCapture(grid,lambda a,b:1.,.001,3.),state,[1.],max_step_s=.05)
    np.testing.assert_allclose(new['snapshots'][0]['moments'],old['snapshots'][0]['moments'],rtol=1e-14)
    assert new['snapshots'][0]['variance'] == pytest.approx(old['snapshots'][0]['variance'],rel=1e-10)
    assert new['snapshots'][0]['ledgers']['supplied_mass'] == old['snapshots'][0]['ledgers']['supplied_mass']


def test_optimizer_plateau_is_admitted_by_unchanged_independent_accuracy_gate():
    # Archived third boundary-pilot failure: the optimizer exhausted its
    # iterations just above its tighter target, despite an accurate density.
    d=e.EntropyMeasure(5.4390656230423598e-6,1.6675232184383115e-10)
    assert d.fit_status == 'iteration_limit'
    assert d.fit_iterations == 80
    assert 2e-11 < d.fit_optimizer_residual < 2e-10
    assert d.fit_residual <= 2e-10
    assert d.fit_normalization_difference <= 2e-10
    nodes,weights=e._positive_rule(d.lower,d.upper,d.A,d.B,96)[:2]
    assert abs(weights@nodes) < 2e-10
    assert abs(weights@nodes**2-1.) < 2e-10


def test_stalled_inaccurate_reconstruction_is_rejected(monkeypatch):
    e._fit.cache_clear()
    monkeypatch.setattr(np.linalg,'solve',lambda a,b:np.zeros_like(b))
    with pytest.raises(e.QuadratureError,match='fail quadrature refinement after non_descent_direction'):
        e.EntropyMeasure(.12,.03)


def test_independent_normalization_gate_rejects_even_converged_optimizer(monkeypatch):
    e._fit.cache_clear()
    original=e._statistics
    def changed_normalization(*args,**kwargs):
        row=list(original(*args,**kwargs))
        if kwargs.get('order') == 48:
            row[-1] += 1e-8
        return tuple(row)
    monkeypatch.setattr(e,'_statistics',changed_normalization)
    with pytest.raises(e.QuadratureError,match='fail quadrature refinement after optimizer_target'):
        e.EntropyMeasure(.4,.0004)


def test_grid_reports_reconstruction_admission_separately_from_storage_changes():
    grid=e.MomentGrid([1.,2.,4.])
    state=grid.from_particles([1.9,2.1],[.5,.5])
    before=grid.fit_diagnostics['validating_calls']
    for _ in range(2):
        grid.distributions(state)
    d=grid.fit_diagnostics
    assert d['validating_calls'] == before+2
    assert d['reconstructed_finite_density_cells'] == before+2
    assert sum(d['optimizer_termination_counts'].values()) == before+2
    assert d['maximum_independently_refined_residual'] <= 2e-10
    assert d['maximum_normalization_log_difference'] <= 2e-10
    assert d['maximum_optimizer_iterations'] > 0
    assert grid.roundoff_diagnostics['storage_reconstructed_cells'] == 0


@pytest.mark.parametrize('mean,variance',[(.4,.239999),(.01,.009899),(.4,.0004)])
def test_positive_sample_quadrature_resolves_narrow_endpoint_modes(mean,variance):
    d=e.EntropyMeasure(mean,variance)
    for lower,upper in ((0.,1.),(0.,.5),(.5,1.)):
        lo,hi=(lower-mean)/d.sd,(upper-mean)/d.sd
        reference_nodes,reference_weights,log_part=e._positive_rule(lo,hi,d.A,d.B,96)
        probability=math.exp(log_part-d.logZ)
        if probability == 0:
            continue
        center=reference_weights@reference_nodes
        scale=math.sqrt(reference_weights@(reference_nodes-center)**2)
        reference_z=(reference_nodes-center)/scale
        for order in (2,4,6):
            nodes,weights=d.quadrature(order,lower,upper)
            assert (weights > 0).all()
            assert nodes.min() >= lower and nodes.max() <= upper
            assert weights.sum() == pytest.approx(probability,rel=2e-10)
            z=((nodes-mean)/d.sd-center)/scale
            # Compare directly integrated standardized moments; raw x powers
            # would conceal errors within a very narrow endpoint peak.
            expected=np.array([reference_weights@reference_z**k for k in range(2*order)])
            actual=np.array([(weights/probability)@z**k for k in range(2*order)])
            np.testing.assert_allclose(actual,expected,rtol=2e-8,atol=2e-8)


def test_positive_sample_quadrature_preserves_finite_atomic_measure():
    nodes=np.array([-2.,-.3,1.2]);weights=np.array([.05,.7,.25])
    q,w=e._gaussian_from_positive_samples(nodes,weights,3,-2.,1.2)
    np.testing.assert_allclose(q,nodes,atol=1e-14)
    np.testing.assert_allclose(w,weights,atol=1e-14)
    with pytest.raises(e.QuadratureError,match='positive sampled measure'):
        e._gaussian_from_positive_samples(nodes,[-.1,.7,.4],3,-2.,1.2)


def test_transformed_cell_warm_start_preserves_density_and_admission():
    e._fit.cache_clear()
    grid=e.MomentGrid([1.,2.,4.]);state=np.zeros((3,3))
    def set_cell(mean,variance):
        u=mean-grid.offsets[1]
        state[1]=[1.,u,u*u+variance]
    set_cell(.1,.04)
    old=grid.distributions(state)[1]
    set_cell(.10001,.039996)
    warm=grid.distributions(state)[1]
    assert warm.fit_warm_start and not warm.fit_cold_fallback
    e._fit.cache_clear()
    cold=e.EntropyMeasure(warm.mean,warm.variance)
    assert warm.fit_iterations < cold.fit_iterations
    np.testing.assert_allclose(warm.cdf(np.linspace(0,1,101)),cold.cdf(np.linspace(0,1,101)),atol=2e-10)
    for d in (warm,cold):
        assert d.fit_residual <= 2e-10
        assert d.fit_normalization_difference <= 2e-10
    assert grid.fit_diagnostics['warm_start_reconstructions'] >= 1


def test_exact_target_cache_ignores_initial_guess_and_bad_seed_falls_back():
    e._fit.cache_clear()
    bad=e.EntropyMeasure(.4,.0004,initial=(10000.,10000.))
    assert bad.fit_warm_start and bad.fit_cold_fallback
    assert bad.A == pytest.approx(0.,abs=1e-10)
    assert bad.B == pytest.approx(-.5,abs=1e-10)
    before=e._fit.cache_info()
    again=e.EntropyMeasure(.4,.0004,initial=(float('nan'),float('nan')))
    after=e._fit.cache_info()
    assert after.hits == before.hits+1 and after.misses == before.misses
    assert (again.A,again.B,again.logZ) == (bad.A,bad.B,bad.logZ)
