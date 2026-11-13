"""Conservation, timing, procurement and frozen-control tests for reclamation."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from provisioning.enclosed_industry import model,reclamation

HERE=Path(reclamation.__file__).parent

@pytest.fixture(scope='module')
def inputs():return reclamation.load()

@pytest.fixture(scope='module')
def saved():return json.loads(gzip.decompress((HERE/'results/reclamation.json.gz').read_bytes()))


def test_existing_net_water_and_goods_are_not_recycled_twice(inputs):
    spec,data,_,_=inputs
    account=reclamation.water_account(spec,spec['normalization_population'])
    assert account['operational_makeup_kg']==pytest.approx(data['demand']['yearly_water'])
    for v in account['streams'].values():
        assert v['recovered_kg']+v['makeup_kg']==pytest.approx(v['gross_kg'])
        assert v['loss_kg']==v['makeup_kg']
    modified=reclamation.configure(data,spec,'closure')
    g=spec['goods'];gross=g['gross_kg_person_year']*spec['normalization_population']
    assert modified['demand']['goods_virgin']+modified['settlement']['operating_import_kg_year']['everyday_carbon_feed']==pytest.approx(gross*(1-g['collection']*g['quality']))
    assert data['demand']['yearly_water']==pytest.approx(account['operational_makeup_kg'])


def synthetic_repair(spec,retirement_year):
    engine=SimpleNamespace(DESTS=model.DESTS,catalogue=lambda data:{'asteroid_tugs':dict(destination='source',owner='industry',price=100.)})
    q=deepcopy(spec);q['refurbishment']['eligible_assets']=['asteroid_tugs']
    rows=[dict(year=y) for y in range(8,16)];flows={y:model.blank_flow() for y in range(6,16)}
    events=[dict(asset='asteroid_tugs',component='active',mass_kg=100.,replacement_kg=100.,order_year=10,commission_year=12)]
    for y in (10,11):flows[y]['imports']['source']=50.;flows[y]['external_replacement']=5000.;flows[y]['industry_replacement_import_kg']=50.
    retired={retirement_year:{('asteroid_tugs','active'):100.}}
    audit=reclamation.refurbishment(engine,{},q,rows,flows,events,retired)
    return q,flows,events,audit


def test_retired_hardware_cannot_supply_same_year_orders(inputs):
    _,flows,events,audit=synthetic_repair(inputs[0],10)
    assert events[0]['refurbished_mass_kg']==0
    assert audit[10]['source']['pending_inspection_kg']==80
    assert audit[11]['source']['qualified_inventory_kg']==48
    assert audit[12]['source']['refurbished_commissioned_kg']==0
    assert sum(a['source']['imported_mass_avoided_kg'] for a in audit.values())==0


def test_accepted_repair_parts_and_rejects_close_mass_and_cash(inputs):
    q,flows,events,audit=synthetic_repair(inputs[0],9)
    assert events[0]['refurbished_mass_kg']==48
    assert events[0]['new_imported_mass_kg']==pytest.approx(61.6)
    assert audit[9]['source']['collected_kg']==80
    assert audit[10]['source']['quality_reject_kg']==32
    assert audit[10]['source']['discarded_old_parts_kg']==pytest.approx(9.6)
    assert audit[12]['source']['refurbished_commissioned_kg']==48
    assert audit[11]['source']['repair_WIP_kg']==pytest.approx(48)
    assert audit[12]['source']['repair_WIP_kg']==pytest.approx(0)
    assert sum(flows[y]['external_replacement'] for y in (10,11))==pytest.approx(6160)
    assert sum(a['source']['import_purchase_avoided_USD2026'] for a in audit.values())==pytest.approx(3840)
    # Collection, rejection and displaced repair parts are three different losses.
    assert sum(a['source']['uncollected_kg']+a['source']['quality_reject_kg']+a['source']['discarded_old_parts_kg'] for a in audit.values())==pytest.approx(61.6)
    for y,a in audit.items():
        assert a['source']['mass_balance_residual_kg']==pytest.approx(0,abs=1e-12)
        assert all(a[k]['withdrawn_for_repair_kg']==0 for k in ('surface','depot','leo'))
        assert all(flows[y]['imports'][k]==0 for k in ('surface','depot','leo'))
    assert sum(sum(f['loop_consume_kg'].values()) for f in flows.values())>0


def test_zero_recovery_removes_no_imports_or_purchase_cost(inputs):
    spec=deepcopy(inputs[0]);spec['refurbishment']['quality']=0
    _,flows,events,audit=synthetic_repair(spec,9)
    assert events[0]['refurbished_mass_kg']==0
    assert sum(flows[y]['external_replacement'] for y in (10,11))==10000
    assert sum(a['source']['quality_reject_kg'] for a in audit.values())==80


def test_no_legacy_allowance_endpoint_keeps_all_explicit_losses(inputs):
    spec,data,_,_=inputs
    e=reclamation.Experiment(data,spec,'combined_without_residual_allowance')
    a=reclamation.water_account(e.spec,250000,improved=True)
    b=reclamation.water_account(spec,250000,improved=True)
    assert a['streams']==b['streams']
    assert a['operational_makeup_kg']>0
    assert b['operational_makeup_kg']-a['operational_makeup_kg']==pytest.approx(spec['landscape_agriculture_allowance_kg_year'])
    assert e.data['demand']['yearly_water']==pytest.approx(a['operational_makeup_kg'])


def test_fresh_engine_does_not_mutate_baseline_module(inputs):
    spec,data,_,_=inputs;original=model.physical
    e=reclamation.Experiment(data,spec,'combined')
    assert e.engine.physical!=original
    assert model.physical is original
    assert 'loop_water_equipment' not in model.catalogue(data)
    assert 'loop_water_equipment' in e.engine.catalogue(e.data)


def test_saved_sources_and_six_complete_frozen_controls(saved,inputs):
    assert len(saved['cases'])==36
    for path,expected in saved['dependencies'].items():assert reclamation.sha(reclamation.ROOT/path)==expected
    assert saved['producer']['sha256']==reclamation.sha(reclamation.__file__)
    controls={c['name']:c for c in inputs[3]['cases']}
    selected=[c for c in saved['cases'] if c['variant']=='baseline']
    assert len(selected)==6
    for c in selected:
        assert c['baseline_exact_match']
        assert c['physical_case']==controls[c['baseline_case']]


def test_all_dated_material_balances_and_no_exhaust_credit(saved):
    for c in saved['cases']:
        assert [r['year'] for r in c['annual']]==list(range(221))
        for r in c['annual']:
            assert r['chemical_exhaust_recovery_kg']==0
            if r['phase']!='operation':continue
            w=r['water'];assert w['new_stock_and_construction_kg']>=0
            assert abs(w['rounding_residual_kg'])<=.001
            for stream in w['streams'].values():assert stream['gross_kg']==pytest.approx(stream['recovered_kg']+stream['loss_kg'])
            for v in r['materials']['structural'].values():
                assert v['retired_kg']==pytest.approx(v['recovered_kg']+v['collection_loss_kg']+v['quality_reject_kg'],abs=.001)
                assert v['opening_inventory_kg']+v['recovered_kg']==pytest.approx(v['used_kg']+v['closing_inventory_kg'],abs=.001)
                assert v['mass_balance_residual_kg']==pytest.approx(0,abs=.001)
            for v in r['materials']['refurbishment'].values():
                assert v['mass_balance_residual_kg']==pytest.approx(0,abs=.001)
                for field in ('pending_inspection_kg','qualified_inventory_kg','repair_WIP_kg'):
                    assert v[field]>=-.001
                assert v['collected_kg']+v['uncollected_kg']==pytest.approx(v['retired_kg'])
                assert v['accepted_kg']+v['quality_reject_kg']==pytest.approx(v['inspection_input_kg'])


def test_seed_products_and_feeds_count_once(saved,inputs):
    source=inputs[1]['industrial']['source']
    for c in saved['cases']:
        for r in c['annual']:
            for key,grade,recovery in [('source_water_equivalent_kg','water_grade','water_recovery'),('source_nitrogen_kg','nitrogen_grade','nitrogen_recovery'),('source_virgin_metal_kg','metal_grade','metal_recovery')]:
                assert r[key]<=r['source_feed_kg']*source[grade]*source[recovery]*(1+1e-12)+.001
        seed=c['physical_case']['seed'];a=c['seed_audit'];rows=c['annual'][:8]
        assert sum(r['source_feed_kg'] for r in rows)==pytest.approx(seed['source_feed_kg'])
        assert sum(r['lunar_raw_kg'] for r in rows)==pytest.approx(seed['lunar_raw_kg'])
        for new,old in [('source_water_equivalent_kg','source_water_kg_year'),('source_nitrogen_kg','source_nitrogen_kg_year'),('source_virgin_metal_kg','source_virgin_metal_kg_year')]:
            assert sum(r[new] for r in rows)==pytest.approx(8*a['physical_annual_average'][old]+a['late_preproduction_physical_total'][old])
        assert c['cumulative_source_feed_kg']==pytest.approx(sum(r['source_feed_kg'] for r in c['annual']))
        assert a['late_preproduction_physical_total']['source_feed_kg_year']>=0


def test_extra_work_is_priced_and_failed_capacity_is_not_erased(saved):
    for c in saved['cases']:
        for r in c['physical_case']['annual']:
            cash=r['cash'];op=sum(cash[k] for k in ('operating_external_replacement_USD2041','operating_external_logistics_materials_USD2041','operating_external_services_USD2041','operating_local_replacement_USD2041','operating_local_payroll_USD2041','operating_resource_provision_USD2041'))
            assert cash['operating_total_USD2041']==pytest.approx(op,abs=.01)
            assert abs(cash['accounting_residual_USD2041'])<.02
            if c['variant'] in ('combined','refurbishment','combined_without_residual_allowance') and r['year']>20:
                assert r['direct_refurbishment_jobs']>0
                assert r['direct_refurbishment_payroll_USD2041']>0
                assert sum(r['physical']['loop_extra_process_kWh_year'].values())>0
                assert sum(r['physical']['loop_repair_consumables_kg_year'].values())>0
            if not r['feasible']:assert r['actual_sales_capacity_fraction']==0
            if r['maximum_capacity_use']>1+1e-6:assert not r['feasible']


def test_one_complete_combined_case_reproduces_saved_product(inputs,saved):
    spec,data,book,_=inputs
    case=next(c for c in book['cases'] if c['name']=='services:breakthrough:reference')
    _,fresh=reclamation.solve_case(spec,data,case,'combined')
    old=next(c for c in saved['cases'] if c['baseline_case']==case['name'] and c['variant']=='combined')
    assert fresh==old
