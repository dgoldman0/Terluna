"""Component stock, aging, construction and immutable producer controls."""
from copy import deepcopy
import gzip
import hashlib
import json
import math

import pytest
from provisioning.enclosed_industry import asset_export as a
from provisioning.enclosed_industry import model


def registry(life=4):
    return {'test:active':dict(asset='test',component='active',mass_share=1.,mean_life_years=life,
                               replacement_price_USD2041_kg=2.,functional_sector='industry')}


def cohort(birth,mass=100,replacement=0,origin='dated_order'):
    return dict(component='test:active',commission_year=birth,order_year=birth-2,
                original_mass_kg=mass,replacement_mass_kg=replacement,origin=origin)


def test_remaining_service_proxy_has_independent_deterministic_solution():
    survival,net=a.service_curves(4,None,8)
    assert survival[:6]==[1,1,1,1,0,0]
    assert net[:6]==[1,.75,.5,.25,0,0]
    rows=a.annual_stocks([cohort(0,origin='seed')],registry(),None,5)
    assert [r['totals']['gross_replacement_value_USD2041'] for r in rows]==[200,200,200,200,0,0]
    assert [r['totals']['net_remaining_service_value_USD2041'] for r in rows]==[200,150,100,50,0,0]
    assert a.check_stock_identities(rows)==0


def test_replacement_is_commissioned_once_and_two_year_wip_conserves_cost():
    rows=a.annual_stocks([cohort(0,origin='seed'),cohort(2,50,25)],registry(),None,6)
    assert rows[0]['totals']['construction_WIP_value_USD2041']==50
    assert rows[1]['totals']['construction_WIP_value_USD2041']==100
    assert rows[2]['totals']['construction_WIP_value_USD2041']==0
    assert rows[2]['totals']['gross_replacement_value_USD2041']==300
    assert rows[2]['totals']['net_remaining_service_value_USD2041']==200
    assert rows[2]['totals']['replacement_commissioned_value_USD2041']==50
    assert rows[4]['totals']['retired_replacement_value_USD2041']==200
    assert a.check_stock_identities(rows)==0


def test_seed_ages_gross_purchases_and_opening_stock_match_source_convention():
    result={'seed':{'stock_mass_kg':{'test':120}},'component_orders':[]}
    cohorts=a.export_cohorts(result,registry(),None,start=8)
    assert {r['commission_year'] for r in cohorts}=={4,5,6,7}
    assert sum(r['original_mass_kg'] for r in cohorts)==160 # one early cohort expires by opening
    annual=a.annual_stocks(cohorts,registry(),None,8)
    assert annual[8]['totals']['surviving_mass_kg']==120
    assert annual[2]['totals']['seed_WIP_value_USD2041']==40
    assert annual[7]['totals']['seed_WIP_value_USD2041']==0
    assert all(y>=0 for c in cohorts for y,_ in a.fabrication_schedule(c))
    assert a.fabrication_schedule(cohort(0,origin='seed'))==[(0,.5),(0,.5)]
    assert a.fabrication_schedule(cohort(1,origin='seed'))==[(0,.5),(0,.5)]
    assert a.check_stock_identities(annual)==0


def test_weibull_remaining_service_counts_retirement_without_double_depreciation():
    s,n=a.service_curves(10,3.,80)
    assert sum(s)==pytest.approx(10)
    assert n[0]==1 and all(0<=n[i]<=s[i]+1e-13 for i in range(len(s)))
    assert all(n[i+1]<=n[i] for i in range(len(n)-1))
    rows=a.annual_stocks([cohort(0,origin='seed')],registry(10),3.,30)
    assert a.check_stock_identities(rows)<1e-10
    total_depreciation=sum(r['totals']['consumption_of_fixed_capital_proxy_USD2041'] for r in rows)
    assert total_depreciation+rows[-1]['totals']['net_remaining_service_value_USD2041']==pytest.approx(200)


def test_uncommissioned_final_orders_remain_wip_and_bad_masses_fail():
    rows=a.annual_stocks([cohort(222,50)],registry(),None,220)
    assert rows[-1]['totals']['construction_WIP_value_USD2041']==50
    assert rows[-1]['totals']['gross_replacement_value_USD2041']==0
    assert a.check_stock_identities(rows)==0
    for c in [cohort(2,-1),cohort(2,2,3)]:
        with pytest.raises(ValueError):a.annual_stocks([c],registry(),None,4)


def test_full_component_registry_covers_housing_civil_services_and_prices_once():
    cfg,data,orders,physical=a.inputs();r=a.components(data);I=data['normalization']['USD2026_to_USD2041']
    assert len(r)==116
    assert {'residential_roads:structure','service_equipment:active','service_civil_bulk:structure',
            'visitor_equipment:active','visitor_civil_metal:structure','source_computer:active'}.issubset(r)
    assert r['residential_roads:structure']['replacement_price_USD2041_kg']==I*data['economy']['production']['local_bulk_value_USD2026_kg']
    assert r['asteroid_tugs:structure']['replacement_price_USD2041_kg']==I*data['economy']['production']['local_fabrication_value_USD2026_kg']
    assert r['asteroid_tugs:active']['replacement_price_USD2041_kg']==I*data['templates']['asteroid_tugs']['price']
    assert {v['functional_sector'] for v in r.values()}=={'industry','settlement','commerce','service','hospitality'}


def test_changed_frozen_record_is_rejected_before_valuation():
    cfg,data,orders,physical=a.inputs();name='propellant_logistics:breakthrough:wider'
    case=next(c for c in orders['cases'] if c['name']==name)
    expected=next(c for c in physical['cases'] if c['name']==name)
    fresh=model.solve(data,case['annual'],max_iterations=50)
    assert a.verify_physics(case,fresh,expected,data)==expected
    altered=deepcopy(expected);altered['annual'][0]['cash']['operating_total_USD2041']+=1
    with pytest.raises(ValueError,match='reproduce admitted'):
        a.verify_physics(case,fresh,altered,data)


def test_saved_product_identity_and_all_stock_sector_accounts():
    path=a.HERE/'results/asset_stock.json.gz'
    if not path.exists():pytest.skip('Asset export has not yet been generated')
    manifest=json.loads(path.with_suffix('.manifest.json').read_text());raw=path.read_bytes();decoded=gzip.decompress(raw)
    assert hashlib.sha256(raw).hexdigest()==manifest['sha256']
    assert hashlib.sha256(decoded).hexdigest()==manifest['decompressed_sha256']
    p=json.loads(decoded);assert p['producer']['sha256']==a.digest(a.HERE/'asset_export.py')
    for name,sha in p['dependencies'].items():assert a.digest(a.ROOT/name)==sha
    assert len(p['cases'])==18 and len({c['name'] for c in p['cases']})==18
    _,_,_,physical=a.inputs();source={c['name']:c for c in physical['cases']}
    for case in p['cases']:
        assert case['source_case_sha256']==hashlib.sha256(json.dumps(source[case['name']],sort_keys=True,separators=(',',':')).encode()).hexdigest()
        assert case['verification']['exact_compact_case_reproduced'] and case['verification']['available_stock_reproduced']
        assert len(case['annual'])==221 and [r['year'] for r in case['annual']]==list(range(221))
        assert case['annual_cost_registry']==[{'year':r['year'],**r['cash']} for r in source[case['name']]['annual']]
        assert case['seed_cost_registry']==source[case['name']]['seed']['cost']
        a.check_stock_identities(case['annual'])
        for row in case['annual']:
            for field in a.VALUES:
                assert math.fsum(s[field] for s in row['by_functional_sector'].values())==pytest.approx(row['totals'][field],rel=1e-12,abs=.1)
                assert row['totals'][field]>=0
            assert row['totals']['net_retirement_value_USD2041']+row['totals']['net_aging_value_USD2041']==pytest.approx(row['totals']['consumption_of_fixed_capital_proxy_USD2041'],rel=1e-12,abs=.1)
            if row['year'] in p['component_checkpoint_years']:
                for field in ('surviving_mass_kg','gross_replacement_value_USD2041','net_remaining_service_value_USD2041','construction_WIP_value_USD2041'):
                    assert math.fsum(s[field] for s in row['by_component'].values())==pytest.approx(row['totals'][field],rel=1e-12,abs=.1)
            else:assert 'by_component' not in row
        assert case['annual'][7]['totals']['seed_WIP_value_USD2041']==0
        assert case['cohort_account']['count']>case['cohort_account']['seed_count']>0
    assert p['unvalued_assets']
