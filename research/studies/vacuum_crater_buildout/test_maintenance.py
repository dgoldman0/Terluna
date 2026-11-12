"""Annual-service boundaries and replacement-event accounting."""
from copy import deepcopy
import json

import pytest

from shared.provenance import constants_changed
from research.studies.vacuum_crater_buildout import maintenance as model


@pytest.fixture
def inputs():
    spec,sources=model.load_inputs()
    return spec,sources,sources['construction_result']['cases']['500']


def test_occupied_power_and_inventory_losses(inputs):
    spec,source,town=inputs
    r=model.account(town,source['construction_inputs'],spec)
    q=r['operating_quantities']
    assert q['electricity_MWh_year'] == pytest.approx(84738)
    assert q['air_makeup_kg_year'] == pytest.approx(18080.50891908222)
    assert q['processed_water_makeup_kg_year'] == pytest.approx(365250)
    assert q['additional_water_makeup_kg_year'] == pytest.approx(913125)
    assert q['air_inventory_cycle_years'] == 100
    assert r['maintained_installed_value_U']['firm_power'] == pytest.approx(53.75e6)
    assert r['variable_services_U_year']['utility_operation'] == pytest.approx(423690)


def test_capital_partition_and_growth_industry_exclusion(inputs):
    spec,source,town=inputs
    r=model.account(town,source['construction_inputs'],spec)
    a=r['maintained_installed_value_U']; initial=town['costs']['reference']['direct_U']
    assert a['ordinary_durable']+a['ordinary_equipment'] == pytest.approx(250e6)
    assert a['pressure_reinforcement']+a['pressure_barriers'] == pytest.approx(initial['pressure_system']/100)
    assert a['independent_support']+a['foundations'] == pytest.approx(initial['independent_support_and_foundations']/100)
    assert a['firm_power'] < initial['dedicated_firm_power']/100
    assert not set(a) & {'qualification','contingency','initial_air_and_water','shielding','removal_reserve','construction_plant'}


def test_annual_components_fund_renewal_once(inputs):
    spec,source,town=inputs
    r=model.account(town,source['construction_inputs'],spec)
    t=r['totals_U_year']
    assert t['routine_upkeep'] == pytest.approx(sum(r['maintained_installed_value_U'].values())*.02)
    assert t['operating_margin'] == pytest.approx(.2*(t['routine_upkeep']+t['variable_services']))
    assert t['total_annual_funding'] == pytest.approx(t['recurring_cash_before_renewal']+t['renewal_fund'])
    assert t['renewal_fund'] == pytest.approx(sum(r['full_replacement_value_U'][k]/life for k,life in r['asset_lives_years'].items()))


def test_replacement_freight_and_local_utility_are_included_once(inputs):
    spec,source,town=inputs
    r=model.account(town,source['construction_inputs'],spec)
    asset=r['maintained_installed_value_U']; replace=r['full_replacement_value_U']
    assert replace['ordinary_equipment']-asset['ordinary_equipment'] == pytest.approx(25e6)
    for k in ('independent_support','foundations','roads','radiators'):
        assert replace[k]-asset[k] == pytest.approx(r['local_replacement_utility_U_per_event'][k])
    for k in ('pressure_reinforcement','pressure_barriers','firm_power'):
        assert replace[k] == asset[k]


def test_one_time_construction_allowances_do_not_become_annual_upkeep(inputs):
    spec,source,town=inputs
    baseline=model.account(town,source['construction_inputs'],spec)
    changed=deepcopy(town)
    c=changed['costs']['reference']; c['per_community_U']*=2
    for k in ('incremental_qualification','regional_industrial_establishment_increment','survey_and_preservation'):
        c['direct_U'][k]*=10
    other=model.account(changed,source['construction_inputs'],spec)
    assert other['totals_U_year'] == baseline['totals_U_year']
    assert other['initial_build_share_U'] == 2*baseline['initial_build_share_U']


def test_no_replacement_at_or_after_removal_endpoint():
    assert model.replacement_count(50,50) == 0
    assert model.replacement_count(200,50) == 3
    assert model.replacement_count(200.1,50) == 4
    assert model.replacement_count(200,25) == 7
    with pytest.raises(ValueError):model.replacement_count(200,0)


def test_horizon_spending_uses_events_not_event_plus_fund(inputs):
    spec,source,town=inputs
    r=model.account(town,source['construction_inputs'],spec)
    h=r['horizons']['200']
    expected=r['totals_U_year']['recurring_cash_before_renewal']*200
    expected+=sum(r['full_replacement_value_U'][k]*count for k,count in h['component_replacements'].items())
    assert h['total_vacuum_expenditure_U'] == pytest.approx(expected)
    assert h['component_replacements']['pressure_reinforcement'] == 3
    assert h['component_replacements']['pressure_barriers'] == 7
    assert h['component_replacements']['foundations'] == 0


def test_longer_lives_reduce_renewal_without_erasing_routine_work(inputs):
    spec,source,town=inputs
    a=model.account(town,source['construction_inputs'],spec)
    changed=deepcopy(spec)
    changed['regimes']['reference']['life_years']={k:2*v for k,v in changed['regimes']['reference']['life_years'].items()}
    b=model.account(town,source['construction_inputs'],changed)
    assert b['totals_U_year']['renewal_fund'] == pytest.approx(a['totals_U_year']['renewal_fund']/2)
    assert b['totals_U_year']['routine_upkeep'] == a['totals_U_year']['routine_upkeep']


def test_operating_deliveries_are_additional_to_construction_capacity(inputs):
    spec,source,town=inputs
    r=model.account(town,source['construction_inputs'],spec)
    region=r['regional_100_town_requirements']; q=r['operating_quantities']
    assert region['steady_air_and_water_import_kg_year'] == pytest.approx(100*(q['air_makeup_kg_year']+q['processed_water_makeup_kg_year']+q['additional_water_makeup_kg_year']))
    assert region['construction_dedicated_volatile_import_kg_year'] == 500e6


def test_cost_price_and_life_cases_are_monotone(inputs):
    spec,source,town=inputs
    for price in ('lower','reference','upper'):
        annual=[model.account(town,source['construction_inputs'],spec,price,k)['totals_U_year']['total_annual_funding'] for k in ('longer_life','reference','higher_wear')]
        assert annual[0] < annual[1] < annual[2]


def test_committed_maintenance_reproduces_and_binds_constants():
    saved=json.loads((model.HERE/'results/maintenance.json').read_text())
    assert not constants_changed(saved['producer']['constants'])
    assert saved == model.run()
