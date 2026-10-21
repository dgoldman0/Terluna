"""Check audit normalization, RNG independence and arithmetic counterexamples."""
import json

import numpy as np
import pytest

from atmosphere.snow import aggregation as ag, flakes as sf, numerical_audit as audit
from atmosphere.snow.habits import GEOMETRIES


def test_exhausted_pair_demonstrates_nonhomogeneity_without_losing_mass():
    rows=audit.arithmetic_examples()['exhausted_pair']
    assert rows[0]['normalized_weights']==[0.,1.5]
    assert rows[1]['normalized_weights']==[.5,1.]
    assert [r['positive_packets'] for r in rows]==[1,2,2]
    for r in rows:
        assert r['mass_before_kg_m2']==pytest.approx(r['mass_after_kg_m2'],rel=1e-15)


def test_fragment_count_cutoff_is_an_absolute_representation_threshold():
    rows=audit.arithmetic_examples()['fragment_threshold']
    assert [r['retained_count_per_m2'] for r in rows]==[0.,.5,.5]
    assert all(r['retained_count_per_m2']+r['discarded_count_per_m2']==.5 for r in rows)


def test_instrumentation_and_private_function_preserve_baseline_random_sequence():
    world=sf.earth_column()
    geometry=GEOMETRIES[audit.GEOMETRY]
    kw=dict(top_c=-.5,number_per_m3=1e4,snowfall_mm_h=3.,bins=2,inject=2,
            fill=1.,gather=.5,max_step_s=30.,seed=4,geometry=geometry)
    original=ag.collide
    expected=ag.column(world,**kw)
    actual,diagnostics=audit.scaled_column(world,area_m2=1.,**kw)
    assert ag.collide is original
    for key in ('masses','numbers'):
        np.testing.assert_array_equal(actual[key],expected[key])
    assert actual['mass_ledger_kg_m2']==expected['mass_ledger_kg_m2']
    assert actual['fragments']==expected['fragments']
    assert diagnostics['candidate_pairs']>0


def test_area_normalization_preserves_collision_free_columns_and_observers(monkeypatch):
    # Suppressing collisions isolates the wrapper's injection/supply/observer
    # normalization; the inherited collision arithmetic is tested separately.
    monkeypatch.setattr(ag,'collide',lambda *a,**kw: 0.)
    world=sf.earth_column()
    geometry=GEOMETRIES[audit.GEOMETRY]
    kw=dict(top_c=-.5,number_per_m3=1e4,snowfall_mm_h=1.,bins=2,inject=2,
            fill=1.,gather=.5,max_step_s=30.,seed=2,geometry=geometry)
    captures=[]
    results=[]
    for area in (1.,2.,256.):
        cloud=[]; arrivals=[]
        def observing(m,n,z,t,dt):
            cloud.append(float(np.dot(m,n)))
        def arriving(m,n,t,dt):
            arrivals.append(float(np.dot(m,n)))
        result,_=audit.scaled_column(world,area_m2=area,instrument=False,
                                     observer=observing,arrival_observer=arriving,**kw)
        captures.append((cloud,arrivals))
        results.append(result)
    assert results[0]['masses'].size>0
    for result,capture in zip(results[1:],captures[1:]):
        np.testing.assert_array_equal(result['masses'],results[0]['masses'])
        np.testing.assert_array_equal(result['numbers'],results[0]['numbers'])
        assert result['mass_ledger_kg_m2']==results[0]['mass_ledger_kg_m2']
        np.testing.assert_array_equal(capture[0],captures[0][0])
        np.testing.assert_array_equal(capture[1],captures[0][1])


def test_audit_product_identity_moments_and_nonempty_controls():
    if not audit.OUT.exists():
        pytest.skip('Numerical audit product not generated yet')
    from shared.provenance import constants_changed
    product=json.loads(audit.OUT.read_text())
    assert product['schema']==audit.SCHEMA
    assert not constants_changed(product['producer']['constants'])
    for key in ('files','inputs'):
        for path,digest in product['producer'][key].items():
            assert audit.sha(audit.ROOT/path)==digest,path
    assert len(product['runs'])==6
    for run in product['runs']:
        assert run['arriving_packet_count']>0
        assert abs(run['mass_ledger_kg_m2']['relative_residual'])<1e-10
        for location in ('cloud','freezing_level'):
            for key,total in run[location]['totals'].items():
                assert total==pytest.approx(sum(run[location][key]),rel=1e-13)
    assert len(product['comparisons'])==4


def test_heavy_replay_matches_recorded_reference_and_keeps_extreme_events():
    if not audit.OUT.exists():
        pytest.skip('Numerical audit product not generated yet')
    product=json.loads(audit.OUT.read_text())
    heavy=product.get('heavy_replay')
    if heavy is None:
        pytest.skip('Full heavy replay not generated yet')
    assert all(heavy['unchanged_replay_checks'].values())
    reference=heavy['reference']
    path=audit.ROOT/reference['path']
    if path.exists():
        assert audit.sha(path)==reference['sha256']
    diag=heavy['collision_diagnostics']
    assert diag['maximum_gt']>1
    assert diag['p_above_one_pairs']>0
    assert diag['first_threshold_events']
    for event in diag['first_threshold_events']:
        assert event['diameter_m']>=event['threshold_m']
        assert event['represented_mass_kg_m2']==pytest.approx(
            event['weight_per_m2']*event['particle_mass_kg'],rel=1e-14)
