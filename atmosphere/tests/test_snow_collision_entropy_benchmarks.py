"""Check the actual control wiring and saved independent distributions."""
from copy import deepcopy
import json
import numpy as np
import pytest
from atmosphere.snow import collision_entropy_benchmarks as b
from shared.baselines import recorded_file


def test_physical_bath_control_starts_from_exact_imposed_collector():
    product=deepcopy(json.loads((b.ROOT/b.control.REFERENCE).read_text()))
    product['capture']['times_s']=[0.]
    target=product['capture']['continuous_mass_references'][-1]
    target['snapshots']=target['snapshots'][:1]
    world,_=b.flakes.moon_column()
    run=b.bath_case(b.reference.SnowKernel(world),product)
    row=run['snapshots'][0]
    initial=product['capture']['initial_collector_mass_kg']
    np.testing.assert_allclose(row['moments'],[1.,initial,initial**2],rtol=1e-13)
    assert row['variance_mass_kg2']==0.
    assert row['number_residual']==0.
    assert row['mass_residual']==0.
    assert row['number_cdf'][0]==0. and row['number_cdf'][-1]==1.
    assert run['diagnostics']['ssp_average']=='state + 0.5*(second_euler_stage-state)'


@pytest.mark.parametrize('case', ['boundary', 'bath'])
def test_saved_entropy_control_identity_and_accounting(case):
    import hashlib
    from shared.provenance import constants_changed
    path=b.ROOT/f'atmosphere/snow/results/collision_entropy_{case}.json'
    product=json.loads(path.read_text())
    assert product['schema']==b.SCHEMA
    for kind in ('files','inputs'):
        for source,digest in product['producer'][kind].items():
            assert hashlib.sha256(recorded_file(source,digest).read_bytes()).hexdigest()==digest
    assert not constants_changed(product['producer']['constants'])
    complete=product['execution_failure'] is None and len(product['runs'])==len(product['requested_controls'])
    assert product['acceptance']['all_requested_controls_completed']==complete
    if not complete:
        assert product['acceptance']['passed'] is False
    assert product['acceptance']['accounting_and_cdf_checks']==b.accounting_checks(
        product['runs'],case,json.loads((b.ROOT/b.control.REFERENCE).read_text()))
    for run in product['runs']:
        for row in run['snapshots']:
            unit=1. if case=='boundary' else product['settings']['initial_collector_mass_kg']
            n,mass,_=row['moments']
            exits=np.asarray(row['ledgers']['lower_exit'])+row['ledgers']['upper_exit']
            assert abs(n+exits[0]-1.)<1e-10
            assert abs(mass+exits[1]-row['ledgers']['supplied_mass']-unit)/unit<1e-10
            assert abs(row['number_residual'])<1e-10
            assert abs(row['mass_residual'])/unit<1e-10
            cdf=np.asarray(row['number_cdf'])
            assert np.isfinite(cdf).all()
            assert cdf.min()>=-1e-10 and cdf.max()<=1+1e-10
            assert np.min(np.diff(cdf))>=-1e-10


def test_saved_entropy_poisson_cdf_against_gamma_identity():
    from scipy.special import gammaincc
    product=json.loads((b.ROOT/'atmosphere/snow/results/collision_entropy_boundary.json').read_text())
    for run in product['runs']:
        donor=run['donor_mass_ratio']
        for row in run['snapshots']:
            query=np.asarray(row['query_mass_ratio'])
            right=query[1::2]
            k=np.rint((right-1)/donor).astype(int)
            np.testing.assert_allclose(right,1+donor*k,rtol=0,atol=4e-16)
            np.testing.assert_array_equal(query[::2],np.nextafter(right,-np.inf))
            exact=np.empty_like(query)
            exact[1::2]=gammaincc(k+1,row['poisson_mean'])
            exact[::2]=0.
            positive=k>0
            exact[::2][positive]=gammaincc(k[positive],row['poisson_mean'])
            np.testing.assert_allclose(exact,row['exact_poisson_number_cdf'],rtol=0,atol=2e-15)
            assert np.max(np.abs(np.asarray(row['number_cdf'])-exact))==pytest.approx(
                row['number_cdf_supremum'],abs=2e-15)
            assert max(row['unqueried_left_probability'],row['unqueried_right_probability'])<1e-20
            expected_variance=row['poisson_mean']*donor**2
            assert row['expected_variance_mass_ratio']==pytest.approx(expected_variance,rel=1e-14)
            assert row['expected_mean_mass_ratio']==pytest.approx(1+donor*row['poisson_mean'],rel=1e-14)


def test_saved_entropy_bath_cdf_against_independent_event_samples():
    product=json.loads((b.ROOT/'atmosphere/snow/results/collision_entropy_bath.json').read_text())
    source=json.loads((b.ROOT/b.control.REFERENCE).read_text())
    targets=source['capture']['continuous_mass_references'][-1]['snapshots']
    for run in product['runs']:
        for row,target in zip(run['snapshots'],targets):
            sample=np.asarray(target['masses_kg'])
            counts=np.asarray(target['sample_counts'])
            cumulative=np.r_[0.,np.cumsum(counts)/sum(counts)]
            query=np.asarray(row['number_cdf_query_mass_kg'])
            empirical=cumulative[np.searchsorted(sample,query,side='right')]
            error=np.max(np.abs(np.asarray(row['number_cdf'])-empirical))
            assert error==pytest.approx(row['empirical_cdf_supremum'],abs=2e-15)
