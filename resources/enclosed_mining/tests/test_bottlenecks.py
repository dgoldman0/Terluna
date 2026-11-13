"""Diagnostic source assessment retains engineering failures and source bytes."""
from resources.enclosed_mining import bottlenecks as b


def test_empty_failed_physical_case_does_not_acquire_a_feed_path(monkeypatch):
    physical={'cases':[dict(name='failed',physical_scenario='capture_8',annual=[],convergence={'converged':False},first_capacity_shortfall=8)]}
    def accept_only_empty(p,*args):
        assert p['cases']==[]
        return {}
    monkeypatch.setattr(b.run,'assess_physical_cases',accept_only_empty)
    cases,unresolved=b.assess(physical,{},{});
    assert cases=={} and unresolved[0]['name']=='failed'
    assert not unresolved[0]['convergence']['converged']


def test_nonempty_capacity_failure_retains_requirement_and_intervention(monkeypatch):
    case=dict(name='test',physical_scenario='service_rollout',annual=[{'year':8}],convergence={'converged':True},first_capacity_shortfall=8,experiment={'variant':'service_rollout'},engineering_intervention={'applied_active_life_years':20.})
    def assessed(p,*args):
        assert p['cases'][0]['first_capacity_shortfall']==8
        return {'service_rollout:test':{'first_physical_capacity_shortfall':8}}
    monkeypatch.setattr(b.run,'assess_physical_cases',assessed)
    cases,unresolved=b.assess({'cases':[case]}, {}, {})
    assert unresolved==[]
    assert cases['service_rollout:test']['first_physical_capacity_shortfall']==8
    assert cases['service_rollout:test']['experiment']==case['experiment']
    assert 'engineering_intervention' not in case['annual'][0]


def test_published_baseline_resource_rows_reproduce_original_assessment():
    import gzip
    import hashlib
    import json
    from pathlib import Path
    root=Path(b.__file__).parent
    original=json.loads(gzip.decompress((root/'results/assessment.json.gz').read_bytes()))
    product=json.loads(gzip.decompress((root/'results/bottlenecks.json.gz').read_bytes()))
    for path,sha in product['dependencies'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==sha
    controls=[c for c in product['cases'].values() if c['experiment']['variant']=='baseline']
    assert len(controls)==4 and len(product['cases'])==28
    for case in controls:
        source=original['cases']['base:'+case['experiment']['baseline_name']]
        for field in ('annual','source_access','checkpoints','source_access_checkpoints','supply_driver_checkpoints'):
            assert case[field]==source[field]
    assert sum(c['first_physical_capacity_shortfall'] is not None for c in product['cases'].values())==7
