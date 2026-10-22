"""Independent checks of saved finite-time capture distributions and ledgers."""
import hashlib
import json
import numpy as np
import pytest
from scipy.special import gammaincc
from atmosphere.snow import collision_jump_benchmarks as b
from shared.provenance import constants_changed


def product(case):
    return json.loads((b.ROOT/f'atmosphere/snow/results/collision_jump_{case}.json').read_text())


@pytest.mark.parametrize('case',['boundary','bath'])
def test_saved_controls_bind_current_sources_and_complete_requested_variations(case):
    p=product(case)
    assert p['schema']==b.SCHEMA
    for group in ('files','inputs'):
        for name,digest in p['producer'][group].items():
            assert hashlib.sha256((b.ROOT/name).read_bytes()).hexdigest()==digest
    assert not constants_changed(p['producer']['constants'])
    assert p['requested_controls']==b.controls(case)
    assert [r['settings'] for r in p['runs']]==p['requested_controls']
    assert not p['execution_failures']
    assert p['acceptance']['passed']
    assert all(r['passed'] for r in p['acceptance']['initial_cdf_checks'])


@pytest.mark.parametrize('case',['boundary','bath'])
def test_saved_native_moments_accounting_and_cdf(case):
    p=product(case)
    unit=1. if case=='boundary' else p['settings']['initial_collector_mass_kg']
    for run in p['runs']:
        pivots=np.asarray(run['pivots_mass_ratio' if case=='boundary' else 'pivots_kg'],np.longdouble)
        width=np.diff(np.asarray(run['bounds_mass_ratio' if case=='boundary' else 'bounds_kg'],np.longdouble))
        for row in run['snapshots']:
            n,u,u2=np.asarray(row['state'],np.longdouble).T
            assert np.min(n)>=0
            native=[n.sum(),np.sum(n*pivots+u*width),np.sum(n*pivots**2+2*u*pivots*width+u2*width**2)]
            np.testing.assert_allclose(native,row['moments'],rtol=3e-13,atol=0)
            exits=np.asarray(row['ledgers']['lower_exit'])+row['ledgers']['upper_exit']
            assert abs(native[0]+exits[0]-1.)<1e-10
            assert abs(native[1]+exits[1]-row['ledgers']['supplied_mass']-unit)/unit<1e-10
            assert abs(row['number_residual'])<1e-10
            assert abs(row['mass_residual'])/unit<1e-10
            active=n>0
            mean=native[1]/native[0]
            local_mean=u[active]/n[active]
            variance=np.sum(n[active]*((pivots[active]+width[active]*local_mean-mean)**2+
                width[active]**2*np.maximum(0.,u2[active]/n[active]-local_mean**2)))/native[0]
            if row['time_s']:
                assert float(variance)==pytest.approx(row['variance'],rel=2e-7)
            q=np.asarray(row['query_mass_ratio' if case=='boundary' else 'number_cdf_query_mass_kg'])
            cdf=np.asarray(row['number_cdf'])
            assert np.isfinite(cdf).all() and np.all(np.diff(q)>0)
            assert cdf.min()>=-1e-10 and cdf.max()<=1+1e-10 and np.min(np.diff(cdf))>=-1e-10


def test_poisson_distributions_against_independent_gamma_identity():
    for run in product('boundary')['runs']:
        donor=run['donor_mass_ratio']
        for row in run['snapshots']:
            q=np.asarray(row['query_mass_ratio']);right=q[1::2]
            k=np.rint((right-1)/donor).astype(int)
            np.testing.assert_array_equal(q[::2],np.nextafter(right,-np.inf))
            exact=np.empty_like(q);exact[1::2]=gammaincc(k+1,row['poisson_mean'])
            exact[::2]=0.;positive=k>0
            exact[::2][positive]=gammaincc(k[positive],row['poisson_mean'])
            np.testing.assert_allclose(exact,row['exact_poisson_number_cdf'],rtol=0,atol=2e-15)
            error=float(np.max(abs(np.asarray(row['number_cdf'])-exact)))
            assert error==pytest.approx(row['number_cdf_supremum'],abs=2e-15)
            assert error<=.005
            assert row['expected_mean_mass_ratio']==pytest.approx(1+donor*row['poisson_mean'],rel=1e-14)
            assert row['expected_variance_mass_ratio']==pytest.approx(donor**2*row['poisson_mean'],rel=1e-14)
        # The shifted-face and alternate-grid controls sample their own crossings.
        assert run['snapshots'][2]['expected_mean_mass_ratio']==pytest.approx(run['reference_cell_face_mass_ratio'],rel=1e-14)


def test_physical_cdf_against_stored_independent_event_histories():
    source=json.loads((b.ROOT/b.old.control.REFERENCE).read_text())
    target=source['capture']['continuous_mass_references'][-1]
    for run in product('bath')['runs']:
        for row,ref in zip(run['snapshots'],target['snapshots']):
            mass=np.asarray(ref['masses_kg']);counts=np.asarray(ref['sample_counts'])
            cumulative=np.r_[0.,np.cumsum(counts)/counts.sum()]
            query=np.asarray(row['number_cdf_query_mass_kg'])
            empirical=cumulative[np.searchsorted(mass,query,side='right')]
            error=float(np.max(abs(np.asarray(row['number_cdf'])-empirical)))
            assert error==pytest.approx(row['empirical_cdf_supremum'],abs=2e-15)
            assert error<=target['number_cdf_uniform_95_percent_dkw_halfwidth']


@pytest.mark.parametrize('case',['boundary','bath'])
def test_truncation_bounds_resolve_centered_variance_and_are_not_total_error_claims(case):
    p=product(case)
    assert 'exclude reconstruction' in p['reading_rule']
    for run,reported in zip(p['runs'],p['acceptance']['truncation_budget_checks']):
        for row,check in zip(run['snapshots'][1:],reported['checks']):
            if case=='boundary':
                bounds=np.asarray(row['poisson_tail_moment_upper_bounds'])
                mean=row['expected_mean_mass_ratio'];variance=row['expected_variance_mass_ratio']
            else:
                bounds=np.asarray(row['uniformization_tail_moment_upper_bounds'])+row['finite_ladder_moment_upper_bounds']
                mean=row['reference_mean_mass_kg'];variance=row['reference_variance_mass_kg2']
            centered=(bounds[2]+2*mean*bounds[1]+mean**2*bounds[0])/variance
            assert centered==pytest.approx(check['centered_second_moment_bound_over_reference_variance'],rel=1e-14,abs=0)
            assert np.all(bounds>=0) and centered<=1e-4
            assert bounds[0]<=1e-12 and bounds[1]/mean<=1e-12


def test_one_and_two_macrosteps_preserve_distribution_within_declared_gate():
    check=product('boundary')['semigroup_control'];a,c=check['runs']
    error=float(np.max(abs(np.asarray(a['number_cdf'])-c['number_cdf'])))
    assert error==pytest.approx(check['number_cdf_difference'],abs=1e-15)
    assert error<=.005
    assert abs(a['variance']/c['variance']-1)<=1e-4
    assert a['diagnostics']['macrosteps']==1 and c['diagnostics']['macrosteps']==2


def test_completed_bath_preserves_prior_attempt_and_repeats_only_unfinished_controls():
    p=product('bath');history=p['execution_history']
    attempt=json.loads((b.ROOT/history['initial_attempt']).read_text())
    assert not attempt['acceptance']['passed']
    assert attempt['execution_failures']
    assert p['producer']['inputs'][history['initial_attempt']]==history['initial_attempt_sha256']
    for saved in attempt['runs']:
        assert next(r for r in p['runs'] if r['settings']==saved['settings'])==saved
    assert history['retained_completed_controls']==[r['settings'] for r in attempt['runs']]
    assert history['repeated_unfinished_controls']==[s for s in attempt['requested_controls'] if s not in history['retained_completed_controls']]
    assert history['continued_wall_limit_s']>history['initial_wall_limit_s']
    assert p['requested_controls']==attempt['requested_controls']
