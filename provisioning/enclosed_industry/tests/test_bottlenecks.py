"""Intervention identity and input-isolation controls, not baseline restatement."""
from copy import deepcopy
import json
from pathlib import Path

import pytest
from provisioning.enclosed_industry import bottlenecks as b, model


@pytest.fixture(scope='module')
def engineering():
    return json.loads((Path(b.__file__).parent/'inputs/engineering.json').read_text())


def experiment(variant='baseline'):
    return dict(baseline_name='services:breakthrough:reference',variant=variant,
                engineering_overrides={'asteroid_tug_active_life_years':30.} if variant=='tug_life_30' else {})


def test_only_source_tug_active_life_changes_and_reference_is_preserved(engineering):
    original=deepcopy(engineering)
    changed,record=b.engineering_for_case(engineering,experiment('tug_life_30'))
    assert engineering==original
    expected=deepcopy(original);expected['templates']['asteroid_tugs']['life']=30.
    assert changed==expected
    assert record['original_active_life_years']==20.
    assert record['structural_life_years']==50.
    assert record['hardware_USD2026_kg']==1500.
    assert sum(model.survival(age,30.) for age in range(400))==pytest.approx(30.,rel=1e-12)


def test_variant_must_authorize_its_exact_engineering_override(engineering):
    bad=experiment();bad['engineering_overrides']={'asteroid_tug_active_life_years':30.}
    with pytest.raises(ValueError,match='disagree'):b.engineering_for_case(engineering,bad)
    bad=experiment('tug_life_30');bad['engineering_overrides']['unpriced_engine_efficiency']=2.
    with pytest.raises(ValueError,match='disagree'):b.engineering_for_case(engineering,bad)


def test_cache_changes_with_request_override_and_every_helper():
    case=dict(name='a',annual=[{'year':8,'goods':1.}],experiment=experiment())
    first=b.cache_identity(case,'engineering',{'model':'one','bootstrap':'two'})
    different=deepcopy(case);different['annual'][0]['goods']=2.
    assert first!=b.cache_identity(different,'engineering',{'model':'one','bootstrap':'two'})
    different=deepcopy(case);different['experiment']=experiment('tug_life_30')
    assert first!=b.cache_identity(different,'engineering',{'model':'one','bootstrap':'two'})
    assert first!=b.cache_identity(case,'engineering',{'model':'one','bootstrap':'three'})
    assert first!=b.cache_identity(case,'changed',{'model':'one','bootstrap':'two'})


def test_case_identity_and_repeated_names_are_rejected():
    ex=experiment();case=dict(name=ex['baseline_name']+':baseline',experiment=ex,market='breakthrough',market_growth_name='reference')
    book=dict(schema=b.ORDER_SCHEMA,cases=[case]);b.validate_book(book)
    with pytest.raises(ValueError,match='Duplicate'):b.validate_book(dict(book,cases=[case,case]))
    different=deepcopy(case);different['name']='silently-renamed'
    with pytest.raises(ValueError,match='disagree'):b.validate_book(dict(book,cases=[different]))


def test_failed_solution_remains_failure_not_empty_success(engineering,monkeypatch):
    case=dict(name='services:breakthrough:reference:baseline',portfolio='services',market='breakthrough',market_growth=.005,market_growth_name='reference',annual=[dict(year=8)],experiment=experiment())
    def fail(*args,**kwargs):raise ValueError('Deliberate unresolved stock requirement')
    monkeypatch.setattr(b.model,'solve',fail)
    result,cached=b.solve_case(engineering,case,'engineering',{})
    assert not cached and not result['convergence']['converged']
    assert result['first_capacity_shortfall']==8 and result['annual']==[]
    assert result['physical_scenario']=='baseline'


def test_published_baselines_equal_original_seed_and_every_annual_row():
    import gzip
    root=Path(b.__file__).parent
    diagnostic=json.loads(gzip.decompress((root/'results/bottlenecks.json.gz').read_bytes()))
    original=json.loads(gzip.decompress((root/'results/industry.json.gz').read_bytes()))
    original={c['name']:c for c in original['cases']}
    controls=[c for c in diagnostic['cases'] if c['experiment']['variant']=='baseline']
    assert len(controls)==4
    for case in controls:
        source=original[case['experiment']['baseline_name']]
        for field in ('seed','annual','convergence','first_capacity_shortfall'):
            assert case[field]==source[field]


def test_published_identity_and_failed_workforce_paths():
    import gzip
    import hashlib
    root=Path(b.__file__).parent;path=root/'results/bottlenecks.json.gz'
    manifest=json.loads(path.with_suffix('.manifest.json').read_text())
    raw=path.read_bytes();decoded=gzip.decompress(raw)
    assert hashlib.sha256(raw).hexdigest()==manifest['sha256']
    assert hashlib.sha256(decoded).hexdigest()==manifest['decompressed_sha256']
    product=json.loads(decoded)
    for name,sha in product['dependencies'].items():assert b.digest(b.ROOT/name)==sha
    assert len(product['cases'])==28
    failed=[c for c in product['cases'] if c['first_capacity_shortfall'] is not None]
    assert len(failed)==7
    for case in failed:
        first=next(r for r in case['annual'] if r['year']==case['first_capacity_shortfall'])
        assert first['shortfalls']['resident_workforce']>0
        assert first['actual_sales_capacity_fraction']==0
    for case in product['cases']:
        assert case['seed']['priced_initial_route_fuel']['total_USD2041']>0
        if case['experiment']['variant']=='tug_life_30':
            assert case['engineering_intervention']['applied_active_life_years']==30.
            assert case['engineering_intervention']['structural_life_years']==50.
            assert case['first_capacity_shortfall'] is None
