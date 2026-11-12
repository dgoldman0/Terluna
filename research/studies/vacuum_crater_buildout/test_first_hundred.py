"""Conserved quantities, monetary boundaries and limiting schedule cases."""
from copy import deepcopy
import json
import math

import pytest

from shared.provenance import constants_changed
from research.studies.vacuum_crater_buildout import first_hundred as model


@pytest.fixture
def inputs():
    return model.load_inputs()


def test_complete_town_geometry_independent_area_calculation(inputs):
    spec, products = inputs
    q = model.quantities(500,spec,products)
    assert q['population'] == 50000
    assert q['production_demand']['building_floor_m2'] == 3500000
    assert q['total_plan_area_m2'] == pytest.approx(3500000/(3*.35*.75))
    assert q['total_roof_area_m2'] == pytest.approx(q['total_plan_area_m2']*1.16)
    assert q['garden_area_m2']/q['population'] == pytest.approx(43+1/3)
    assert q['shielding_kg'] == pytest.approx(24746666666.66667)
    assert q['water_inventory']['lakes_kg'] == pytest.approx(2222222222.222222)
    assert 237 < q['dome_diameter_m'] < 239


def test_extraction_products_and_retained_tails_close(inputs):
    spec, products = inputs
    q = model.quantities(500,spec,products)
    assert q['raw_feed_kg'] == pytest.approx(q['local_finished_kg']+q['retained_tails_kg'])
    imported = sum(q['special_imports'].values())+q['production_demand']['volatile_import_kg']
    assert q['local_share_of_all_delivered_material'] == pytest.approx(q['local_finished_kg']/(q['local_finished_kg']+imported))
    assert q['production_demand']['handling_kg'] > q['raw_feed_kg']  # Separate grading pass, no new reserve credit.


def test_population_scaling_keeps_pressure_volume_cubic(inputs):
    spec, products = inputs
    small, large = [model.quantities(n,spec,products) for n in (500,1000)]
    assert large['total_plan_area_m2'] == pytest.approx(2*small['total_plan_area_m2'])
    assert large['pressure_reinforcement_kg'] == pytest.approx(small['pressure_reinforcement_kg']*2**1.5)
    assert large['air_including_reserve_kg'] == pytest.approx(small['air_including_reserve_kg']*2**1.5)
    assert large['shielding_kg'] == pytest.approx(2*small['shielding_kg'])


def test_better_strength_does_not_shrink_shield_or_fixed_barriers(inputs):
    spec, products = inputs
    a = model.quantities(500,spec,products)
    spec = deepcopy(spec)
    spec['quantity']['reference_allowable_Pa'] *= 2
    b = model.quantities(500,spec,products)
    assert b['pressure_reinforcement_kg'] == pytest.approx(a['pressure_reinforcement_kg']/2)
    assert b['barriers_kg'] == a['barriers_kg']
    assert b['shielding_kg'] == a['shielding_kg']
    assert b['independent_support_kg'] == a['independent_support_kg']


def test_lighting_peak_and_reserve_are_distinct(inputs):
    spec, products = inputs
    q = model.quantities(500,spec,products)
    s = q['services']
    assert s['garden_average_W'] == pytest.approx(216666666.6666667)
    assert s['garden_peak_W']/s['garden_average_W'] == pytest.approx(1.5)
    assert s['firm_capacity_W'] > s['peak_occupied_W'] > s['mean_occupied_W']
    assert s['reserve_Wh'] == pytest.approx(12e9)
    assert s['full_night_storage_Wh'] > 40*s['reserve_Wh']


def test_monetary_allowance_and_initial_funding_partition(inputs):
    spec, products = inputs
    q = model.quantities(500,spec,products)
    c = model.construction_cost(q,spec,products,'reference')
    assert c['direct_U']['ordinary_installed_capital'] == 25e9
    assert c['total_initial_funding_U'] == pytest.approx(sum(c['direct_U'].values())*1.4)
    assert c['per_community_U']*100 == pytest.approx(c['total_initial_funding_U'])
    changed = deepcopy(spec)
    changed['prices']['reference']['industrial_establishment_increment_U'] += 1e9
    d = model.construction_cost(q,changed,products,'reference')
    assert d['total_initial_funding_U']-c['total_initial_funding_U'] == pytest.approx(1.4e9)


def toy_case():
    return dict(lead_years={'survey':2,'plant':4,'qualification':3},effective_handling_systems=1,
                handling_tonnes_hour_system=1,handling_hours_year=.1,bulk_finished_kg_year=100,
                metal_finished_kg_year=100,pressure_system_kg_year=100,special_import_kg_year=100,
                volatile_import_kg_year=100,building_floor_m2_year=100,road_km_year=100,
                site_teams=100,site_cycle_years=1)


def settings():
    return dict(ramp_fractions=[],steps_per_year=12,maximum_years=50)


def test_constant_production_exact_capacity_and_parallel_leads():
    case=toy_case()
    demand=dict.fromkeys(model.production_rates(case),100)
    result=model.schedule(demand,100,case,settings())
    assert result['lead_years'] == 4  # Parallel prerequisites, not 2+4+3.
    assert result['total_years'] == pytest.approx(6)
    assert result['occupation_years']['1'] == pytest.approx(4+1/12+1)
    assert result['consumed'] == demand


def test_site_team_occupancy_limits_concurrency():
    case=toy_case()
    case['site_teams']=1
    demand=dict.fromkeys(model.production_rates(case),2)
    result=model.schedule(demand,2,case,settings(),ready_industry=True)
    assert result['total_years'] == pytest.approx(2+1/12)
    assert result['slowest_steady_stream'] == 'site_teams'


def test_ramp_is_integrated_and_ready_case_removes_only_lead_and_ramp():
    case=toy_case()
    demand=dict.fromkeys(model.production_rates(case),150)
    config=settings()
    config['ramp_fractions']=[.25,.5,.75]
    result=model.schedule(demand,100,case,config)
    ready=model.schedule(demand,100,case,config,ready_industry=True)
    assert result['total_years'] == pytest.approx(4+3+1)
    assert ready['total_years'] == pytest.approx(1.5+1)
    for key in demand:
        assert result['produced_before_last_start'][key] == pytest.approx(result['consumed'][key]+result['remaining_inventory_at_last_start'][key])


def test_slow_stream_controls_duration():
    case=toy_case()
    case['pressure_system_kg_year']=10
    demand=dict.fromkeys(model.production_rates(case),100)
    result=model.schedule(demand,100,case,settings(),ready_industry=True)
    assert result['total_years'] == pytest.approx(11)
    assert result['slowest_steady_stream'] == 'pressure_system_kg'
    for key in demand:
        assert result['produced_before_last_start'][key] == pytest.approx(demand[key])
        assert result['remaining_inventory_at_last_start'][key] == pytest.approx(0, abs=1e-8)


def test_invalid_capacity_fails_instead_of_inventing_completion():
    case=toy_case()
    demand=dict.fromkeys(model.production_rates(case),100)
    case['metal_finished_kg_year']=0
    with pytest.raises(ValueError,match='capacity'):
        model.schedule(demand,100,case,settings())


def test_smaller_timestep_preserves_required_resources(inputs):
    spec, products=inputs
    demand=model.quantities(500,spec,products)['production_demand']
    a=model.schedule(demand,100,spec['production']['reference'],spec['schedule'])
    refined=deepcopy(spec['schedule'])
    refined['steps_per_year']=120
    b=model.schedule(demand,100,spec['production']['reference'],refined)
    assert abs(a['total_years']-b['total_years']) <= .25
    for k in demand:
        assert b['consumed'][k] == pytest.approx(demand[k])
        assert b['remaining_inventory_at_last_start'][k] >= -1e-5


def test_committed_product_and_input_identities():
    committed=json.loads((model.HERE/'results/first_hundred.json').read_text())
    assert not constants_changed(committed['producer']['constants'])
    assert committed == model.run()


def test_establishment_arrives_before_production_and_is_counted_once(inputs):
    spec, products=inputs
    q=model.quantities(500,spec,products)
    assert q['production_demand']['special_import_kg']+q['establishment_import_kg'] == pytest.approx(sum(q['special_imports'].values()))
    case=toy_case()
    demand=dict.fromkeys(model.production_rates(case),100)
    result=model.schedule(demand,100,case,settings(),establishment_import_kg=250)
    assert result['lead_years'] == pytest.approx(10)
    assert result['total_years'] == pytest.approx(12)


def test_future_reservation_scales_with_usable_return_losses(inputs):
    spec, products=inputs
    q=model.quantities(500,spec,products)
    lower, middle, upper=q['future_raw_feed_reservation']
    assert lower['raw_feed_kg'] == pytest.approx(2*middle['raw_feed_kg'])
    assert upper['raw_feed_kg'] == pytest.approx(.2*middle['raw_feed_kg'])
    assert middle['processed_asset_inventory_cycle_years'] == pytest.approx(1000)
    assert middle['raw_feed_kg'] == pytest.approx(.2*(q['raw_feed_kg']-q['shielding_kg']))
