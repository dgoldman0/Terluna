"""Conservation, reserve floors and reproduction accounting checks."""
import math

import pytest

from research.studies.floater_viability.biology import (
    bud_attachment_budget, dark_storage, evaluate, fermentation_feed,
    gross_carbon_budget, net_surplus_budget, photosynthetic_hydrogen_allocation,
    reproduction_budget,
)


def test_gross_carbon_ledger_respects_growth_and_maintenance():
    case = gross_carbon_budget(3, 1, 0.003)
    assert case['carbon_residual_kg_m2_year'] == pytest.approx(0)
    assert case['npp_c_kg_m2_year'] == pytest.approx(case['new_dry_kg_m2_year']*.45)
    assert case['gpp_c_kg_m2_year'] == pytest.approx(
        case['maintenance_c_kg_m2_year'] + case['new_dry_kg_m2_year']*1.39*.4)
    assert case['growth_respiration_c_kg_m2_year'] > 0


def test_maintenance_failure_cannot_generate_positive_growth():
    case = gross_carbon_budget(1, 10, 0.031)
    assert not case['structural_maintenance_feasible']
    assert case['new_dry_kg_m2_year'] == 0
    assert case['npp_c_kg_m2_year'] < 0


def test_dark_storage_has_fixed_tissue_floor_and_separate_reserve_mass():
    case = dark_storage(0.01, 15, community_c_kg_m2_day=0.002)
    assert case['minimum_stored_c_kg_m2'] == pytest.approx(.18)
    assert case['minimum_reserve_dry_kg_m2'] == pytest.approx(.405)
    assert case['reserve_decay_loss_c_kg_m2'] == pytest.approx(0)
    doubled = dark_storage(.01, 30, community_c_kg_m2_day=.002)
    assert doubled['minimum_stored_c_kg_m2'] == pytest.approx(.36)


def test_dark_storage_with_decay_reaches_the_prescribed_final_reserve():
    case = dark_storage(.01, 15, reserve_decay_per_day=.02,
                        reserve_recovery_fraction=.8, final_usable_reserve_c_kg_m2=.03)
    final = (case['initial_usable_c_kg_m2']+.01/.02)*math.exp(-.02*15)-.01/.02
    assert final == pytest.approx(.03)
    assert case['minimum_stored_c_kg_m2']*.8 == pytest.approx(case['initial_usable_c_kg_m2'])
    assert case['reserve_decay_loss_c_kg_m2'] > 0


def test_net_npp_does_not_pay_baseline_respiration_again():
    assert net_surplus_budget(1)['remaining_c_kg_m2_year'] == 1
    cost = net_surplus_budget(1, hydrogen_cost_w_m2=1)
    assert cost['hydrogen_opportunity_c_kg_m2_year'] == pytest.approx(.78894)
    assert not net_surplus_budget(1, hydrogen_cost_w_m2=2)['allocation_feasible']
    assert net_surplus_budget(1, turnover_dry_kg_m2_year=1)['remaining_c_kg_m2_year'] == pytest.approx(.55)


def test_hydrogen_photon_allocation_cannot_also_fix_carbon():
    case = photosynthetic_hydrogen_allocation(1, 100, .02, 3)
    assert case['required_hydrogen_area_time_fraction'] == .5
    assert case['remaining_gpp_c_kg_m2_year'] == 1.5
    failed = photosynthetic_hydrogen_allocation(3, 100, .02, 3)
    assert not failed['area_time_feasible']
    assert failed['remaining_gpp_c_kg_m2_year'] == 0
    assert failed['unmet_hydrogen_w_m2'] == 1


def test_hydrogen_fermentation_conserves_molar_feed_ratio():
    # Four H2 moles from one mole of glucose in the idealised acetate route.
    case = fermentation_feed(4*.002016, 4)
    assert case['glucose_feed_kg'] == pytest.approx(.180156)
    assert fermentation_feed(1, 2)['glucose_feed_kg'] == pytest.approx(
        2*fermentation_feed(1, 4)['glucose_feed_kg'])


def test_bud_can_be_supported_before_independent_lift():
    case = bud_attachment_budget(100, 300, 250)
    assert case['attached_mass_feasible']
    assert not case['independent_mass_feasible']
    assert not bud_attachment_budget(40, 300, 250)['attached_mass_feasible']


def test_recruitment_and_construction_define_different_budget_stages():
    net = reproduction_budget(100, .5, 100, reproductive_lifetime_years=10,
                              recruitment_probability=.1)
    assert net['minimum_years_per_bud'] == pytest.approx(.9)
    assert net['lifetime_replacement_possible']
    gross = reproduction_budget(100, .5, 100, reproductive_lifetime_years=10,
                                recruitment_probability=.1, construction_glucose_kg_kg_dry=1.39)
    assert gross['minimum_years_per_bud'] == pytest.approx(1.112)
    assert not gross['lifetime_replacement_possible']
    assert reproduction_budget(100, 0, 100)['minimum_years_per_bud'] is None


def test_invalid_traits_and_complete_small_sensitivity():
    with pytest.raises(ValueError):
        gross_carbon_budget(1, 1, .001, construction_glucose_kg_kg_dry=.5)
    with pytest.raises(ValueError):
        dark_storage(.01, 15, reserve_recovery_fraction=0)
    with pytest.raises(ValueError):
        fermentation_feed(1, 0)
    assert len(evaluate()['maintenance']) == 36
