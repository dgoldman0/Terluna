import math

import pytest

from research.studies.world_economy.mature_economics import component_cost, compare_locations, renewal_demand


SHARES = dict(feed=.35,gravity_step=.25,pressure_gas=.01,pressure_shell=.03,
              pressure_cycle=.03,other=.33)


def cost(**kw):
    settings = dict(base_attempt_cost=100,shares=SHARES,factors={},accepted_yield=1,annual_output=1000)
    settings.update(kw)
    return component_cost(**settings)


def test_pressure_penalty_applies_only_to_its_component():
    result=cost(factors={'pressure_gas':1.2})
    assert result['factory_cost_per_qualified_unit_U']==pytest.approx(100.2)
    assert result['qualified_component_costs_U']['pressure_gas']==1.2


def test_independent_step_saving_and_named_pressure_costs_close():
    result=cost(factors={'gravity_step':.8,'pressure_shell':math.sqrt(1.2)})
    assert result['factory_cost_per_qualified_unit_U']==pytest.approx(95+3*(math.sqrt(1.2)-1))
    assert sum(result['qualified_component_costs_U'].values())==pytest.approx(result['factory_cost_per_qualified_unit_U'])


def test_yield_changes_accepted_output_cost_not_delivered_good_units():
    low=cost(accepted_yield=.9);high=cost(accepted_yield=.95)
    assert high['factory_cost_per_qualified_unit_U']/low['factory_cost_per_qualified_unit_U']==pytest.approx(.9/.95)
    assert low['required_annual_attempts']-low['rejected_annual_attempts']==pytest.approx(1000)
    assert low['annual_qualified_output']==high['annual_qualified_output']==1000


def test_rotation_annual_service_and_delivery_charged_once():
    result=cost(annual_site_cost=2000,delivery_cost_per_good=3)
    assert result['factory_cost_per_qualified_unit_U']==102
    assert result['delivered_cost_per_qualified_unit_U']==105
    assert result['annual_levelized_production_cost_U']==102000


def test_break_even_component_reduction_reproduces_comparator():
    earth=cost(factors={'gravity_step':.9},delivery_cost_per_good=4)
    array=cost(annual_site_cost=2000,delivery_cost_per_good=1)
    moon=cost(factors={'pressure_gas':1.2},delivery_cost_per_good=5)
    comparison=compare_locations(dict(moon=moon,earth=earth,array=array),'gravity_step')
    required=comparison['total_reduction_fraction_of_reference_component_for_parity']
    amended=cost(factors={'pressure_gas':1.2,'gravity_step':1-required},delivery_cost_per_good=5)
    assert amended['delivered_cost_per_qualified_unit_U']==pytest.approx(min(earth['delivered_cost_per_qualified_unit_U'],array['delivered_cost_per_qualified_unit_U']))


def test_cheaper_orbital_recipe_blocks_false_earth_only_advantage():
    result=compare_locations(dict(moon=cost(),earth=cost(delivery_cost_per_good=10),
        array=cost(annual_site_cost=2000),array_freefall=cost(factors={'gravity_step':.5})),'gravity_step')
    assert result['cheapest_alternative']=='array_freefall'
    assert not result['moon_strictly_cheaper_than_all_alternatives']


def test_unattainable_threshold_is_retained():
    result=compare_locations(dict(moon=cost(delivery_cost_per_good=100),
        earth=cost(),array=cost()),'gravity_step')
    assert result['total_reduction_fraction_of_reference_component_for_parity']>1
    assert not result['component_reduction_alone_can_reach_parity']


def test_recurring_demand_closes_and_all_caps_apply():
    result=renewal_demand(140e12,40,.02,.5,63.7e9,50e9,3)
    assert result['total_annual_renewal_units']==3.5e12
    assert result['eligible_annual_orders_units']==70e9
    assert result['assumed_awarded_annual_orders_units']==35e9
    assert result['realized_annual_units']==pytest.approx(50e9/3)
    assert result['annual_customer_spending_U']+result['unused_customer_budget_U']==pytest.approx(50e9)


@pytest.mark.parametrize('options',[dict(shares={'a':.9}),dict(accepted_yield=0),
    dict(accepted_yield=1.1),dict(annual_output=0),dict(annual_site_cost=float('inf')),
    dict(factors={'unknown':1}),dict(factors={'gravity_step':-.1})])
def test_invalid_cost_inputs_rejected(options):
    with pytest.raises(ValueError):cost(**options)


def test_missing_competitor_rejected():
    with pytest.raises(ValueError):compare_locations(dict(moon=cost(),earth=cost()),'gravity_step')


def test_invalid_award_share_rejected():
    with pytest.raises(ValueError):renewal_demand(1,10,1,1.1,1,1,1)
