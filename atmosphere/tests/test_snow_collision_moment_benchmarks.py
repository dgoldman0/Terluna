"""Retain the failed shape control and its independently checkable accounting."""
import hashlib
import json
import numpy as np
import pytest
from atmosphere.snow import collision_moment_benchmarks as b
from shared.baselines import recorded_file


def test_empirical_distance_checks_both_sides_of_atoms():
    masses=np.array([1.,2.]);counts=np.array([1,1])
    queries=np.sort(np.r_[masses,np.nextafter(masses,-np.inf)])
    assert b.empirical_distance(queries,[0.,.5,.5,1.],masses,counts)==0.
    assert b.empirical_distance(queries,[.5,.5,1.,1.],masses,counts)==.5


def test_saved_beta_boundary_failure_remains_distinct_from_moment_agreement():
    product=json.loads((b.ROOT/'atmosphere/snow/results/collision_moment_boundary.json').read_text())
    assert product['case']=='poisson_boundary'
    assert product['measure_name']=='bounded_beta'
    assert product['acceptance']['passed'] is False
    rows=product['runs'][0]['snapshots']
    center=rows[2]
    assert center['number_cdf_supremum']>.1
    assert center['mean_mass_ratio']==pytest.approx(center['expected_mean_mass_ratio'],rel=1e-10)
    assert center['variance_mass_ratio']==pytest.approx(center['expected_variance_mass_ratio'],rel=1e-7)
    assert rows[-1]['number_cdf_supremum']<.005
    for row in rows:
        q=np.asarray(row['query_mass_ratio']);cdf=np.asarray(row['number_cdf'])
        exact=np.asarray(row['exact_poisson_number_cdf'])
        assert (np.diff(q)>0).all()
        assert np.min(np.diff(cdf))>=-1e-11
        assert max(abs(row['number_residual']),abs(row['mass_residual']))<1e-11
        assert max(abs(cdf-exact))==pytest.approx(row['number_cdf_supremum'])
        assert row['unqueried_left_probability']<1e-20
        assert row['unqueried_right_probability']<1e-20


def test_saved_candidate_sources_and_closed_bath_ledgers():
    for case in ('bath','boundary'):
        product=json.loads((b.ROOT/f'atmosphere/snow/results/collision_moment_{case}.json').read_text())
        for kind in ('files','inputs'):
            for path,sha in product['producer'][kind].items():
                assert hashlib.sha256(recorded_file(path,sha).read_bytes()).hexdigest()==sha
        for run in product['runs']:
            for row in run['snapshots']:
                if case=='bath':
                    initial=product['settings']['initial_collector_mass_kg']
                    n,m1,_=row['moments']
                    exits=np.asarray(row['ledgers']['lower_exit'])+row['ledgers']['upper_exit']
                    assert n+exits[0]==pytest.approx(1.,abs=1e-11)
                    assert m1+exits[1]==pytest.approx(initial+row['ledgers']['supplied_mass'],rel=1e-11)
                    assert abs(row['mass_residual'])/initial<1e-11
