"""Published-product provenance, cash identities and scope controls."""
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import pytest

HERE=Path(__file__).resolve().parents[1]
ROOT=HERE.parents[1]

@lru_cache(maxsize=None)
def load(name):
    path=HERE/'results'/f'{name}.json.gz'
    manifest=json.loads(path.with_suffix('.manifest.json').read_text())
    raw=path.read_bytes();decoded=gzip.decompress(raw)
    assert hashlib.sha256(raw).hexdigest()==manifest['sha256']
    assert hashlib.sha256(decoded).hexdigest()==manifest['decompressed_sha256']
    product=json.loads(decoded)
    for field in ('sha256','model_sha256'):
        pathfield='path' if field=='sha256' else 'model_path'
        assert hashlib.sha256((ROOT/product['producer'][pathfield]).read_bytes()).hexdigest()==product['producer'][field]
    assert hashlib.sha256((ROOT/product['postprocessing']['path']).read_bytes()).hexdigest()==product['postprocessing']['sha256']
    for path,sha in product['dependencies'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    return product


def test_base_dimensions_and_completion():
    p=load('industry');assert len(p['cases'])==81
    assert len({c['name'] for c in p['cases']})==81
    for c in p['cases']:
        assert c['convergence']['converged']
        assert len(c['annual'])==213
        assert c['annual'][0]['year']==8 and c['annual'][-1]['year']==220
        cap=next(r for r in c['annual'] if r['year']==200)['population']
        assert cap==pytest.approx(6468289.238000949)
        assert all(r['population']==pytest.approx(cap) for r in c['annual'] if r['year']>=200)


def test_cash_reserves_and_payroll_count_once():
    for c in load('industry')['cases']:
        seed=c['seed'];assert seed['priced_initial_route_fuel']['total_USD2041']>0
        parts={k:v for k,v in seed['cost'].items() if k not in ('total_USD2041','identified_service_hardware_USD2041')}
        assert sum(parts.values())==pytest.approx(seed['cost']['total_USD2041'],rel=1e-12)
        assert seed['remaining_unpriced_requirements']
        for row in c['annual']:
            x=row['cash']
            assert abs(x['accounting_residual_USD2041'])<1.
            assert x['operating_resource_provision_USD2041']+x['capital_resource_provision_USD2041']==pytest.approx(x['resource_provision_USD2041'],rel=1e-12)
            assert x['capital_local_USD2041']+x['operating_local_payroll_USD2041']==pytest.approx(x['resident_payroll_USD2041'],rel=1e-12)
            assert x['operating_local_replacement_USD2041']==0
            assert x['identified_service_stock_value_USD2041']>=0


def test_adverse_inputs_and_unresolved_cases_remain_explicit():
    lunar=load('poor_lunar');source=load('poor_source');life=load('lifetimes')
    assert lunar['scenarios']['lunar']['metal_yield']==.05
    assert source['scenarios']['source']['water_grade']==.05
    assert source['scenarios']['source']['water_recovery']==.5
    assert life['scenarios']['lifetime_shape'] is None
    for product in (lunar,source,life):
        for case in product['cases']:
            if not case['convergence']['converged']:
                assert case['first_capacity_shortfall'] is not None
                assert not case['annual'] or all(row['actual_sales_capacity_fraction']==0 for row in case['annual'])
            for row in case['annual']:
                assert row['actual_sales_capacity_fraction']==(1. if row['feasible'] else 0.)


def test_all_values_finite():
    def walk(value):
        if isinstance(value,dict):
            for item in value.values():walk(item)
        elif isinstance(value,list):
            for item in value:walk(item)
        elif isinstance(value,float):assert math.isfinite(value)
    for name in ('industry','poor_lunar','poor_source','lifetimes'):walk(load(name))
