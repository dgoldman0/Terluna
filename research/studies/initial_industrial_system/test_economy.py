"""Physical and economic controls for joint production and capacity reuse."""
from copy import deepcopy
import math
import pytest
from . import evolution as m
from . import economy as e
from . import run as base

@pytest.fixture(scope='module')
def data():
    return m.read_inputs()

@pytest.fixture(scope='module')
def result():
    return e.run()


def test_inherited_no_commerce_reproduces_transport_mass_balance(data):
    q,src,s,c=data;a=m.annual(*data)
    for k,old in [('source_water_kg_year','source_water_kg_year'),('source_feed_kg_year','source_ore_kg_year'),('lunar_metal_kg_year','metal_kg_year'),('lunar_bulk_kg_year','bulk_kg_year'),('LEO_kg_year','LEO_kg_year')]:
        assert a['physical'][k]==pytest.approx(c['steady_flows'][old],rel=1e-10)


def test_input_mismatch_is_an_error(data):
    q=deepcopy(data[0]);q['inputs']['settlement']['sha256']='0'*64
    with pytest.raises(ValueError,match='Input hash mismatch'):m.read_inputs(q)

@pytest.mark.parametrize('profile',['inherited','modular_imported','local_structures'])
def test_component_streams_and_local_share_conserve_mass(data,profile):
    a=m.annual(*data,profile=profile,growth_people=500,new_plant_scale=.002)
    for r in a['components']:
        expected=r['mass']*(r['structural_share']/r['structural_life']+(1-r['structural_share'])/r['life'])*r['stock_scale']
        assert r['local_renew_kg']+r['import_renew_kg']==pytest.approx(expected)
        assert 0<=r['local_share']<=r['structural_share']<=1
        if r['class'] in ['pressure','computer']:assert r['local_share']==0


def test_active_module_prices_change_cash_not_qualified_mass_or_capacity(data):
    q,src,s,c=data;v=deepcopy(q);v['production']['special_component_price_multiplier']=4
    a=m.annual(*data,profile='local_structures');b=m.annual(v,src,s,c,profile='local_structures')
    assert a['physical']==b['physical']
    assert b['cash_USD2026']['equipment_imports']>a['cash_USD2026']['equipment_imports']*3


def test_joint_feed_and_residue_close_without_selling_potential(data):
    q,src,s,c=data;a=m.annual(*data,profile='local_structures');f=a['physical'];co=m.coproducts(q,s,c,a)
    assert f['source_feed_kg_year']==pytest.approx(f['source_water_kg_year']+f['source_nitrogen_kg_year']+f['source_virgin_metal_kg_year']+f['source_retained_tails_kg_year'])
    assert co['raw_residue_revenue_USD2026']==0
    assert co['products']['aggregate_recoverable_metals']['unused_recovery_potential_kg_year']>20*co['current_extra_metal_refining_capacity_kg_year']


def test_refining_coproducts_uses_existing_joint_feed_and_capacity(data):
    a=m.annual(*data,profile='local_structures');b=m.annual(*data,profile='local_structures',extra_metal=100e6)
    assert b['physical']['source_metal_kg_year']-a['physical']['source_metal_kg_year']==pytest.approx(100e6*(1+data[0]['coproducts']['stock_containment_kg_per_kg']))
    assert b['physical']['source_feed_kg_year']<a['physical']['source_feed_kg_year']*1.02
    assert b['physical']['source_process_W']>a['physical']['source_process_W']
    assert b['physical']['imports_by_destination_kg_year']['source']>a['physical']['imports_by_destination_kg_year']['source']
    assert not m.annual(*data,profile='local_structures',extra_metal=1e9)['feasible']


def test_expanded_refinery_has_stock_power_delivery_and_new_capital(result):
    a=result['expanded_source_refinery'];commission=result['expanded_source_refinery_commissioning']
    assert a['engineering']['feasible'] and commission['feasible_overlap']
    assert commission['additional_external_purchases_USD2041']>70e9
    assert commission['local_civil_installation_allowance_USD2041']>0
    assert a['engineering']['physical']['source_metal_kg_year']>1e9
    assert result['expanded_source_refinery_incremental_external_cost_USD2041']>0
    assert a['external']['USD2041']['gross_external_sales']==0
    assert a['domestic']['USD2041']['inventory_change']>0


def test_nutrient_requirements_keep_element_and_compound_mass_distinct(result):
    n=result['nutrient_coproduct_screen']
    assert n['annual_element_requirement_kg']['carbon']==pytest.approx(5e6)
    assert n['upper_bound_displaceable_Earth_feed_kg_year']>sum(n['annual_element_requirement_kg'].values())
    assert n['Bennu_analogue_recovered_carbon_comparison_kg_year']>n['annual_element_requirement_kg']['carbon']
    assert result['coproducts']['raw_residue_revenue_USD2026']==0


def test_source_export_bypasses_surface_but_loads_both_space_routes(data):
    a=m.annual(*data,profile='local_structures');b=m.annual(*data,profile='local_structures',source_goods=50e6)
    assert b['physical']['Moon']['trips_year']==pytest.approx(a['physical']['Moon']['trips_year'])
    assert b['physical']['source']['trips_year']>a['physical']['source']['trips_year']
    assert b['physical']['Earth']['burned_kg_year']>a['physical']['Earth']['burned_kg_year']
    assert b['physical']['Earth']['Earth_fuel_kg_year']==0


def test_capsules_are_relaunched_and_not_free_return_freight(data):
    q,src,s,c=data;a=m.annual(*data,profile='local_structures',goods=100e6)
    assert a['cash_USD2026']['capsule_operations']>0
    assert a['cash_USD2026']['capsule_replacement']>0
    assert a['physical']['imports_by_destination_kg_year']['leo']>50e6
    assert a['saleable_export_kg_year']==pytest.approx(100e6*.99*.98)


def test_recovered_remote_frames_pay_for_return_routes(data):
    a=m.annual(*data,profile='local_structures')
    f=a['physical'];r=f['retired_components_kg_year']
    assert f['recycled_to_source_kg_year']==pytest.approx((r['depot']+r['leo'])*data[0]['production']['source_recovery'])
    assert f['recycled_from_LEO_kg_year']==pytest.approx(r['leo']*data[0]['production']['source_recovery'])
    assert f['source']['trips_year']>0


def test_operations_keep_maintenance_consumables_without_replacement(data):
    a=m.annual(*data,profile='local_structures',renewal_factor=0)
    assert a['cash_USD2026']['equipment_imports']==0
    assert a['physical']['imports_by_destination_kg_year']['source']>0
    assert a['cash_USD2026']['industrial_external_services']>0


def test_growth_cannot_use_unbuilt_capacity(data):
    q,src,s,c=data;f=m.frontier(*data,'growth_people',balanced_growth=True)
    a=f['at_limit'];assert a['new_plant_scale_year']>0
    assert a['capacity_use'][a['binding_capacity']]==pytest.approx(1,rel=1e-6)
    b=m.annual(*data,profile='local_structures',growth_people=f['maximum']*1.01,new_plant_scale=f['maximum']*1.01/c['demand']['population'])
    assert not b['feasible']


def test_growth_and_maximum_exports_compete(data,result):
    a=result['frontiers']['goods']['maximum'];b=result['balanced_growth_frontier']['local_structures']['new_capacity_people_year']
    assert not m.annual(*data,profile='local_structures',goods=a,growth_people=b,new_plant_scale=b/250000,service_portfolio=True)['feasible']
    assert result['high_value_capacity_with_growth']['engineering']['growth_people_year']<b


def test_frontier_is_independent_of_old_factory_budget_and_per_capita_cap(data):
    q,src,s,c=data;v=deepcopy(src);v['revenue_inputs']['manufacturing']['nameplate_finished_kg_per_resident_year']=1
    v['revenue_inputs']['markets']['reference']['factory_customer_budget_USD2026_year']=1
    a=m.frontier(*data,'goods')['maximum'];b=m.frontier(q,v,s,c,'goods')['maximum']
    assert a==pytest.approx(b) and a>200e6


def test_commissioning_pays_for_new_stock_and_qualification(data,result):
    a=result['high_value_capacity_commissioning'];assert a['additional_external_purchases_USD2041']>0
    assert a['local_civil_installation_allowance_USD2041']>0
    assert a['additional_source_feed_kg']>0
    assert a['feasible_overlap']


def test_growth_has_two_year_commissioning_delay_and_explicit_cap(result):
    path=result['growth']['None'];rows=path['years']
    assert rows[0]['occupied_capacity_people']==rows[1]['occupied_capacity_people']==250000
    assert rows[2]['occupied_capacity_people']==pytest.approx(250000+rows[0]['new_capacity_started_people'])
    assert all(r['feasible'] for r in rows)
    assert path['occupied_capacity_at_end_people']+path['pipeline_people']<=path['vacuum_architectural_reference_with_existing_cities_people']
    assert path['eventual_crater_architectural_reference_people']>500e6


def test_growth_zero_preserves_population_and_is_not_external_sales(result):
    assert result['growth']['0']['occupied_capacity_at_end_people']==250000
    a=result['high_value_capacity_case']['external']['USD2041'];b=result['high_value_capacity_with_growth']['external']['USD2041']
    assert a['gross_external_sales']==b['gross_external_sales']
    assert b['expansion_external_purchases']>0
    assert b['after_expansion_and_resource']<b['after_renewal_and_resource']

@pytest.mark.parametrize('name',['commodity','specialty','qualified_high_value'])
def test_supply_use_and_cash_accounts_reconcile(result,name):
    case=result['portfolios'][name];a=case['domestic']['USD2041'];rows=case['domestic']['sectors_USD2041'];t=case['external']['USD2041']
    assert sum(r['value_added'] for r in rows)==pytest.approx(a['GDP_GVA'])
    assert a['GDP_GVA']==pytest.approx(a['resident_consumption']+a['gross_capital_formation']+a['exports']-a['imports'])
    assert a['gross_domestic_output']==pytest.approx(a['GDP_GVA']+sum(r['domestic_intermediates']+r['imported_intermediates'] for r in rows))
    assert a['imports']==pytest.approx(t['total_imports'])
    assert a['net_domestic_product']==pytest.approx(a['GDP_GVA']-a['consumption_of_fixed_capital'])
    assert t['after_renewal_and_resource']==pytest.approx(t['gross_external_sales']-t['operating_imports']-t['external_renewal']-t['resource_provision'])
    assert case['external']['debt_service_USD2041'] is None


def test_imported_capital_does_not_become_domestic_output(result):
    case=result['portfolios']['specialty'];a=case['domestic']['USD2041'];rows=case['domestic']['sectors_USD2041']
    capital=rows[2]
    assert capital['gross_output']==pytest.approx(a['domestic_capital_formation'])
    assert a['gross_capital_formation']-capital['gross_output']==pytest.approx(a['capitalized_external_purchases'])


def test_internal_currency_revaluation_cannot_fix_external_cash(data):
    q,src,s,c=data;v=deepcopy(q);v['valuation']['inherited_U_service_normalization_USD2041_per_U']=60
    a=e.evaluate(*data,profile='local_structures');b=e.evaluate(v,src,s,c,profile='local_structures')
    assert a['external']==b['external']
    assert b['domestic']['USD2041']['GDP_GVA']>a['domestic']['USD2041']['GDP_GVA']


def test_compute_earth_comparator_includes_common_services(data):
    q=data[0];p=e.compute_price(q);ex=q['exports']
    assert p*(1-ex['compute_common_Earth_service_share'])==pytest.approx(47/5+1+.1*base.HOURS*.6/1000)
    a=m.annual(*data,profile='local_structures',compute=1e9);b=m.annual(*data,profile='local_structures')
    assert a['physical']['source_heat_W']-b['physical']['source_heat_W']>1e9
    assert a['cash_USD2026']['equipment_imports']>b['cash_USD2026']['equipment_imports']


def test_best_tested_cases_show_losses_and_conditional_surplus(result):
    assert result['best_tested_trade_cases']['3']['goods_kg_year']==0
    assert result['best_tested_trade_cases']['2000']['after_renewal_and_resource']>0
    assert result['high_value_capacity_with_growth']['external']['USD2041']['after_expansion_and_resource']>0
    assert result['high_value_break_even']['whole_system_USD2026_per_saleable_kg']>1000


def test_million_year_feed_is_not_declared_permitted_inventory(result):
    r=result['resource_horizon']
    assert r['one_million_year_constant_feed_kg']==pytest.approx(r['constant_local_structure_feed_kg_year']*1e6)
    assert result['currency']['adopted_exchange_rate'] is None


def test_coproduct_inventory_stays_inside_off_earth_economy(result):
    a=result['coproduct_inventory_comparison'];d=a['domestic']['USD2041']
    assert a['external']['USD2041']['gross_external_sales']==0
    assert d['inventory_change']>0
    assert d['gross_capital_formation']==pytest.approx(d['gross_fixed_capital_formation']+d['inventory_change'])
    assert d['GDP_GVA']==pytest.approx(d['GDP_from_final_expenditure'])
    assert all(x['credited_external_revenue_USD2041']==0 for x in result['trace_metal_sensitivity'])
