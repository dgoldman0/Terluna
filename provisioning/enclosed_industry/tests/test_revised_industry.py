"""Revised order linkage and repair-aware asset reconciliation."""
from copy import deepcopy
import gzip,json
from pathlib import Path
import pytest
from provisioning.enclosed_industry import asset_export as a,revised_industry as r

HERE=Path(r.__file__).parent

@pytest.fixture(scope='module')
def inputs():return r.load()

@pytest.fixture(scope='module')
def physical():return json.loads(gzip.decompress((HERE/'results/revised_industry.json.gz').read_bytes()))

@pytest.fixture(scope='module')
def assets():return json.loads(gzip.decompress((HERE/'results/revised_asset_stock.json.gz').read_bytes()))


def registry():return {'machine:active':dict(asset='machine',component='active',mass_share=1.,mean_life_years=4,replacement_price_USD2041_kg=2.,functional_sector='industry')}

def cohort(birth,mass=100.,repaired=0.,origin='dated_order'):
    return dict(component='machine:active',commission_year=birth,order_year=birth-2,original_mass_kg=mass,replacement_mass_kg=0.,origin=origin,refurbished_mass_kg=repaired,new_unit_mass_kg=mass-repaired,repair_parts_mass_kg=repaired*.2,retained_core_mass_kg=repaired*.8)


def test_repaired_units_replace_retired_units_once_and_wip_prices_new_parts_only():
    annual,residual=r.revised_annual_stocks([cohort(0,origin='seed'),cohort(7,repaired=80)],registry(),None,8)
    assert residual==pytest.approx(0)
    assert annual[4]['totals']['gross_replacement_value_USD2041']==0
    assert annual[5]['totals']['construction_WIP_value_USD2041']==36
    assert annual[6]['totals']['construction_WIP_value_USD2041']==72
    final=annual[7]['totals']
    assert final['gross_replacement_value_USD2041']==200
    assert final['net_remaining_service_value_USD2041']==200
    assert final['construction_WIP_value_USD2041']==0
    assert final['new_input_commissioned_quote_USD2041']==72
    assert final['retained_core_service_restoration_USD2041']==128
    assert final['refurbished_commissioned_value_USD2041']==160
    assert final['refurbished_commissioned_mass_kg']==80
    assert final['new_units_commissioned_mass_kg']==20
    assert sum(x['totals']['repair_parts_purchase_USD2041'] for x in annual)==32
    assert annual[4]['totals']['retired_replacement_value_USD2041']==200


def test_new_only_stock_export_retains_original_asset_values():
    cohorts=[cohort(0,origin='seed'),cohort(6,50)]
    old=a.annual_stocks(cohorts,registry(),None,8);new,residual=r.revised_annual_stocks(cohorts,registry(),None,8)
    assert residual==0
    for x,y in zip(old,new):
        assert x['by_component']==y['by_component']
        for key in a.VALUES:assert x['totals'][key]==y['totals'][key]
        assert y['totals']['retained_core_service_restoration_USD2041']==0


def test_all_repaired_no_new_units_has_no_extra_complete_asset():
    rows,_=r.revised_annual_stocks([cohort(0,origin='seed'),cohort(7,repaired=100)],registry(),None,8)
    assert rows[6]['totals']['construction_WIP_value_USD2041']==40
    assert rows[7]['totals']['surviving_mass_kg']==100
    assert rows[7]['totals']['new_units_commissioned_mass_kg']==0
    assert rows[7]['totals']['new_input_commissioned_quote_USD2041']==40
    assert rows[7]['totals']['retained_core_service_restoration_USD2041']==160


def test_requested_comparison_keeps_primary_orders_and_separates_overrides(inputs):
    cfg,spec,data,orders,baseline,old,rollout=inputs
    planned=r.plan_cases(cfg,orders,rollout)
    assert len(planned)==93
    assert {g:sum(c['group']==g for c in planned) for g in ('primary','poor_source','rollout')}==dict(primary=81,poor_source=6,rollout=6)
    assert [c['case'] for c in planned[:81]]==orders['cases']
    assert all(not c['source_override'] for c in planned if c['group']!='poor_source')
    assert all(c['source_override']==cfg['poor_source'] for c in planned if c['group']=='poor_source')
    assert all(c['case']['experiment']['baseline_name']==c['baseline_case'] for c in planned if c['group']=='rollout')


def test_solver_failure_retains_unresolved_case_instead_of_invented_assets(inputs,monkeypatch):
    cfg,spec,data,orders,_,_,_=inputs
    planned=dict(case=orders['cases'][0],baseline_case=orders['cases'][0]['name'],variant='combined',group='primary',source_override={})
    def fail(*args,**kwargs):raise ValueError('synthetic unsatisfied physical solve')
    monkeypatch.setattr(r,'captured_case',fail)
    out=r.solve_planned(cfg,spec,data,planned,{})
    assert not out['case']['admitted_completed_path']
    assert out['case']['physical_case']['annual']==[]
    assert out['asset'] is None
    assert 'synthetic' in out['case']['physical_case']['convergence']['error']


def test_saved_products_share_physics_and_bind_all_sources(physical,assets):
    assert len(physical['cases'])==93 and len(assets['cases'])==18
    assert assets['source_physical_sha256']==r.digest(HERE/'results/revised_industry.json.gz')
    for product in (physical,assets):
        assert product['producer']['sha256']==r.digest(r.__file__)
        for name,expected in product['dependencies'].items():assert r.digest(r.ROOT/name)==expected
    byname={c['name']:c for c in physical['cases']}
    for asset in assets['cases']:
        source=byname[asset['name']]
        assert source['baseline_case']==asset['baseline_case'] and source['variant']=='combined'
        assert asset['source_case_sha256']==r.value_digest(source['physical_case'])
        assert asset['source_loop_case_sha256']==r.value_digest(source)
        assert asset['annual_cost_registry']==[dict(year=x['year'],**x['cash']) for x in source['physical_case']['annual']]


def test_all_primary_and_sensitivity_order_links_match_frozen_requests(physical,inputs):
    cfg,_,_,orders,baseline,old,rollout=inputs
    books={c['name']:c for c in orders['cases']+rollout['cases']};controls={c['name']:c for c in baseline['cases']}
    assert len(physical['baseline']['case_identities'])==81
    for name,h in physical['baseline']['case_identities'].items():assert h==r.value_digest(controls[name])
    for c in physical['cases']:
        order=books[c['order_case_name']]
        assert c['source_order_case_sha256']==r.value_digest(order)
        assert c['name']==c['baseline_case']+':'+c['variant']
        assert c['physical_case']['name']==c['order_case_name']
        if c['variant']=='combined_poor_source':
            for k,v in cfg['poor_source'].items():assert c['source_coefficients'][k]==v
        if c['physical_case']['first_capacity_shortfall'] is not None:assert not c['admitted_completed_path']


def test_asset_rollforward_has_no_saved_purchase_asset(assets):
    for c in assets['cases']:
        assert r.check_asset_identities(c['annual'])<=c['verification']['maximum_stock_identity_residual_USD2041']+1e-8
        assert [x['year'] for x in c['annual']]==list(range(221))
        assert c['cohort_account']['refurbished_original_mass_kg']>0
        for row in c['annual']:
            v=row['totals']
            assert v['net_remaining_service_value_USD2041']<=v['gross_replacement_value_USD2041']+.1
            assert v['commissioned_replacement_value_USD2041']==pytest.approx(v['new_input_commissioned_quote_USD2041']+v['retained_core_service_restoration_USD2041'],abs=.1)
            assert v['retained_core_service_restoration_USD2041']<=v['refurbished_commissioned_value_USD2041']+.1
            for key in a.VALUES+r.ASSET_EXTRA:
                assert sum(s[key] for s in row['by_functional_sector'].values())==pytest.approx(v[key],abs=.1)


def test_baseline_seed_bounds_preserve_feed_and_contain_existing_exact_seeds(physical,inputs):
    _,_,data,orders,baseline,previous,_=inputs;byname={c['name']:c for c in baseline['cases']};records={a['name']:a for a in physical['baseline']['seed_material_intervals']}
    assert len(records)==81
    for name,record in records.items():
        assert record['exact_source_feed_kg']==byname[name]['seed']['source_feed_kg']
        assert record['average_source_feed_kg']+record['late_source_feed_kg']==pytest.approx(record['exact_source_feed_kg'])
        for v in record['products'].values():assert 0<=v['total_lower_kg']<=v['total_upper_kg']
    for old in previous['cases']:
        if old['variant']!='baseline':continue
        bound=records[old['baseline_case']]
        for name,key in [('water','source_water_equivalent_kg'),('nitrogen','source_nitrogen_kg'),('virgin_metal','source_virgin_metal_kg')]:
            actual=sum(x[key] for x in old['annual'][:8]);v=bound['products'][name]
            assert v['total_lower_kg']-.01<=actual<=v['total_upper_kg']+.01


def test_one_primary_physical_and_asset_case_fresh_reproduction(physical,assets,inputs):
    cfg,spec,data,orders,baseline,old,rollout=inputs
    planned=next(p for p in r.plan_cases(cfg,orders,rollout) if p['baseline_case']=='services:breakthrough:reference' and p['variant']=='combined')
    previous={c['baseline_case']:c for c in old['cases'] if c['variant']=='combined'}
    fresh=r.solve_planned(cfg,spec,data,planned,previous)
    assert fresh['case']==next(c for c in physical['cases'] if c['name']==fresh['case']['name'])
    assert fresh['asset']==next(c for c in assets['cases'] if c['name']==fresh['asset']['name'])


def test_all_revised_material_accounts_cash_and_capacity_flags(physical):
    for case in physical['cases']:
        source=case['source_coefficients'];q=case['physical_case']
        if not q['annual']:
            assert not case['admitted_completed_path']
            continue
        assert [row['year'] for row in case['annual']]==list(range(221))
        for row in case['annual']:
            for field,grade,recovery in [('source_water_equivalent_kg','water_grade','water_recovery'),('source_nitrogen_kg','nitrogen_grade','nitrogen_recovery'),('source_virgin_metal_kg','metal_grade','metal_recovery')]:
                assert row[field]<=row['source_feed_kg']*source[grade]*source[recovery]*(1+1e-12)+.001
            assert row['chemical_exhaust_recovery_kg']==0
            if row['phase']!='operation':continue
            for values in row['materials']['refurbishment'].values():
                assert abs(values['mass_balance_residual_kg'])<.01
                assert values['qualified_inventory_kg']>=-.001
                assert values['repair_WIP_kg']>=-.001
            for values in row['materials']['structural'].values():assert abs(values['mass_balance_residual_kg'])<.01
        for row in q['annual']:
            cash=row['cash'];subtotal=sum(cash[k] for k in ('operating_external_replacement_USD2041','operating_external_logistics_materials_USD2041','operating_external_services_USD2041','operating_local_replacement_USD2041','operating_local_payroll_USD2041','operating_resource_provision_USD2041'))
            assert cash['operating_total_USD2041']==pytest.approx(subtotal,rel=2e-15,abs=.1)
            assert abs(cash['accounting_residual_USD2041'])<.1
            if row['maximum_capacity_use']>1+1e-6:
                assert not row['feasible'] and row['actual_sales_capacity_fraction']==0
                assert not case['admitted_completed_path']
