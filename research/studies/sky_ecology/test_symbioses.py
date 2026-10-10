"""Dimensional scaling, independent hand budgets, and allocation boundaries."""
import json
import math
from pathlib import Path

import pytest

from research.studies.sky_ecology import symbioses as m


TRAITS = dict(dry_fraction=.3, carbon_per_dry=.5, phosphorus_per_dry=.01,
              carbon_assimilation=.8, population_production_efficiency=.25,
              phosphorus_assimilation=.5)


def guild(name='patrol', share=1., density=.01):
    return dict(name=name, live_kg_per_footprint_m2=density, turnover_per_year=4.,
                host_budget_share=share,
                external_food_c_supply_kg_m2_year=0., external_diet_p_supply_kg_m2_year=0.)


def test_one_hectare_at_ten_grams_is_hundred_kg_not_tiny_individuals():
    p = m.population(.01, 100., 10000., .3, .2, 3.6)
    assert p['live_kg'] == 100.
    assert p['dry_kg'] == 30.
    assert p['individuals'] == 1.
    assert p['active_individuals'] == .2  # expectation, not a fraction of a real defender
    assert p['individuals_per_exterior_m2'] == pytest.approx(1/36000)


def test_population_area_individual_size_and_active_fraction_are_independent():
    a = m.population(.01, .001, 100., .3, .2, 3.6)
    b = m.population(.01, .002, 200., .3, .1, 7.2)
    assert b['individuals'] == a['individuals']
    assert b['live_kg'] == 2*a['live_kg']
    assert b['active_individuals'] == a['active_individuals']/2
    assert b['active_per_exterior_m2'] == a['active_per_exterior_m2']/8


def test_grid_response_hand_example_and_deadline_inverse():
    # A sqrt(2)*10 m square has a 10 m corner-to-centre path, doubled by obstacles.
    t = m.response_seconds(math.sqrt(2)*10, .5, 30., 2.)
    assert t == pytest.approx(70.)
    assert m.spacing_for_deadline(70., .5, 30., 2.) == pytest.approx(math.sqrt(2)*10)
    assert m.spacing_for_deadline(20., .5, 30., 2.) is None
    assert m.spacing_for_deadline(30., .5, 30., 2.) is None
    assert m.response_seconds(20., .2, 30., 1.)-30 == pytest.approx(
        2*(m.response_seconds(20., .4, 30., 1.)-30))


def test_independent_annual_carbon_and_phosphorus_hand_budget():
    # 10 g live -> 3 g dry, replaced four times: 12 g dry = 6 g C + .12 g P.
    # 80% C assimilation and 25% production efficiency require 30 g food C.
    r = m.annual_requirements(.01, turnover_per_year=4., **TRAITS)
    assert r['replacement_dry_kg_m2_year'] == pytest.approx(.012)
    assert r['replacement_c_kg_m2_year'] == pytest.approx(.006)
    assert r['replacement_p_kg_m2_year'] == pytest.approx(.00012)
    assert r['food_c_required_kg_m2_year'] == pytest.approx(.03)
    assert r['dietary_p_required_kg_m2_year'] == pytest.approx(.00024)
    assert r['respiratory_c_kg_m2_year'] == pytest.approx(.018)
    assert r['unassimilated_c_kg_m2_year'] == pytest.approx(.006)
    assert r['food_c_required_kg_m2_year'] == pytest.approx(sum(r[k] for k in (
        'replacement_c_kg_m2_year', 'respiratory_c_kg_m2_year', 'unassimilated_c_kg_m2_year')))


def test_density_and_turnover_scale_requirements_without_extra_respiration():
    a = m.annual_requirements(.01, turnover_per_year=2., **TRAITS)
    b = m.annual_requirements(.02, turnover_per_year=4., **TRAITS)
    for key in ('food_c_required_kg_m2_year', 'replacement_p_kg_m2_year', 'respiratory_c_kg_m2_year'):
        assert b[key] == pytest.approx(4*a[key])


def test_shared_host_budget_rejects_double_spending():
    with pytest.raises(ValueError, match='over-allocate'):
        m.allocate_food([guild('one', .7), guild('two', .7)], .1, TRAITS)
    out = m.allocate_food([guild('one', .1), guild('two', .2)], .1, TRAITS)
    assert out['host_c_allocated_kg_m2_year'] == pytest.approx(.03)
    assert out['host_c_unallocated_kg_m2_year'] == pytest.approx(.07)
    assert out['host_c_used_kg_m2_year'] == pytest.approx(.03)
    assert sum(r['food_c_shortfall_kg_m2_year'] for r in out['guilds']) == pytest.approx(.03)


def test_sugar_carbon_can_close_but_zero_phosphorus_cannot_replace_animals():
    r = m.allocate_food([guild()], .1, TRAITS)['guilds'][0]
    assert r['food_c_shortfall_kg_m2_year'] == 0
    assert r['assimilable_p_supply_kg_m2_year'] == 0
    assert r['incorporation_p_shortfall_kg_m2_year'] == pytest.approx(.00012)


def test_external_supply_is_separate_and_p_is_not_all_assimilated():
    g = guild(share=.1)
    g.update(external_food_c_supply_kg_m2_year=.02,
             external_diet_p_supply_kg_m2_year=.00012)
    r = m.allocate_food([g], .1, TRAITS)['guilds'][0]
    assert r['host_c_used_kg_m2_year'] == pytest.approx(.01)
    assert r['external_c_used_kg_m2_year'] == pytest.approx(.02)
    assert r['external_c_unused_supply_kg_m2_year'] == 0
    assert r['food_c_shortfall_kg_m2_year'] == pytest.approx(0., abs=1e-15)
    assert r['incorporation_p_shortfall_kg_m2_year'] == pytest.approx(.00006)


def test_unused_external_food_carries_its_p_and_both_conserve():
    g = guild()
    g.update(external_food_c_supply_kg_m2_year=.06, external_diet_p_supply_kg_m2_year=.00024)
    r = m.allocate_food([g], .1, TRAITS)['guilds'][0]
    assert r['external_c_used_kg_m2_year'] == pytest.approx(.03)
    assert r['host_c_used_kg_m2_year'] == 0
    assert r['external_diet_p_consumed_kg_m2_year'] == pytest.approx(.00012)
    assert r['assimilable_p_supply_kg_m2_year'] == pytest.approx(.00006)
    assert r['external_c_supply_kg_m2_year'] == pytest.approx(
        r['external_c_used_kg_m2_year']+r['external_c_unused_supply_kg_m2_year'])
    assert r['external_diet_p_supply_kg_m2_year'] == pytest.approx(
        r['external_diet_p_consumed_kg_m2_year']+r['external_diet_p_unconsumed_kg_m2_year'])
    assert r['external_diet_p_consumed_kg_m2_year'] == pytest.approx(
        r['assimilable_p_supply_kg_m2_year']+r['unassimilated_consumed_p_kg_m2_year'])


def test_unpaired_food_phosphorus_is_rejected():
    g = guild()
    g['external_diet_p_supply_kg_m2_year'] = .0001
    with pytest.raises(ValueError, match='paired'):
        m.allocate_food([g], .1, TRAITS)


def test_nutrient_import_boundary_conserves_p_and_can_be_negative():
    r = m.nutrient_import(10., .4, .5, .25, 3.)
    assert r['deposited_p_kg_m2_year'] == 4.
    assert r['retained_import_p_kg_m2_year'] == 2.
    assert r['deposition_loss_p_kg_m2_year'] == 2.
    assert r['initially_bioavailable_import_p_kg_m2_year'] == .5
    assert r['net_total_p_import_kg_m2_year'] == -1.
    assert r['deposited_p_kg_m2_year'] == r['retained_import_p_kg_m2_year']+r['deposition_loss_p_kg_m2_year']


@pytest.mark.parametrize('fn,args', [
    (m.population, (-.01, 1., 100., .3, .2, 3.6)),
    (m.population, (.01, 0., 100., .3, .2, 3.6)),
    (m.response_seconds, (20., 0., 30., 1.)),
    (m.response_seconds, (20., .2, 30., .5)),
    (m.response_seconds, (20., float('nan'), 30., 1.)),
    (m.nutrient_import, (1., 1.1, .5, .5, 0.)),
])
def test_invalid_inputs_fail(fn, args):
    with pytest.raises(ValueError):
        fn(*args)


def test_zero_turnover_cannot_claim_zero_maintenance_cost():
    with pytest.raises(ValueError, match='positive turnover'):
        m.annual_requirements(.01, turnover_per_year=0., **TRAITS)


def test_product_binds_source_register_and_uses_reef_surface_product(tmp_path):
    # Temporary register verifies binding only; it is not a scholarly source register.
    register = tmp_path/'sources.json'
    register.write_text('{"test_only": 1}')
    root = Path(m.st.__file__).parents[3]
    first = m.results(root=root, sources_path=register)
    register.write_text('{"test_only": 2}')
    second = m.results(root=root, sources_path=register)
    assert first['producer']['inputs']['symbioses_sources.json'] != second['producer']['inputs']['symbioses_sources.json']
    assert first['exterior_per_footprint'] == pytest.approx(2*math.pi/math.sqrt(3))
    assert first['producer']['model_sha256'] == m._digest(Path(m.__file__))
    assert first['producer']['dependencies']['stoichiometry.py'] == m._digest(Path(m.st.__file__))
    assert first['schema'] == 'terluna.research.sky-reef-symbioses/1'
    assert len(first['annual_requirements']) == 9


def test_zero_host_budget_rejected_before_product_fraction_division(tmp_path):
    scenario = json.loads((Path(m.__file__).parent/'symbioses_scenarios.json').read_text())
    scenario['host_food_c_kg_per_footprint_m2_year'] = 0.
    path = tmp_path/'scenario.json'
    path.write_text(json.dumps(scenario))
    with pytest.raises(ValueError, match='positive host budget'):
        m.results(scenarios_path=path)


def test_stored_symbiosis_product_regenerates_exactly():
    import json
    stored = json.loads((m.HERE/'results/symbioses.json').read_text())
    assert stored == m.results()


def test_reference_host_food_matches_parent_product():
    import json
    parent = json.loads((m.HERE/'results/sky_ecology.json').read_text())
    data = m.results()
    reference = parent['food_web']['sky']['kg_c_m2_year']['food_to_tenants']
    assert reference[0] == reference[1] == data['scenarios']['host_food_c_kg_per_footprint_m2_year']
