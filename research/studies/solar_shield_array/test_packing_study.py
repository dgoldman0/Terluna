"""Keep the local packing result, actual costs and unresolved fleet scope distinct."""
import json
from pathlib import Path

import pytest

from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):return json.loads((HERE/'results'/name).read_text())


def test_packing_products_bind_the_executed_models_and_inputs():
    for name in ['packing_screen.json','packing.json','packing_validation.json','packing_budget.json','packing_audit.json']:
        p=read(name)
        for path,h in p['producer']['source_hashes'].items():assert digest(ROOT/path)==h,path
        for path,h in p['producer']['inputs'].items():assert digest(HERE/'results'/path)==h,path
        assert not constants_changed(p['producer']['constants'])
        assert not p['accepted_fleet']


def test_selected_packing_preserves_inventory_and_replays_the_whole_local_passage():
    screen=read('packing_screen.json');candidate=read('packing.json');fine=read('packing_validation.json')
    assert len(screen['configs'])==6 and len(screen['cases'])==12 and screen['count']==361
    selected=candidate['selected_config']
    assert selected['depth_m']==12000. and not selected['swapped']
    assert selected['rows']==19 and selected['pitch_m']==8500.
    assert candidate['accepted_departure'] and fine['accepted_service_departure']
    assert len(fine['runs'])==2
    for entry,start,end in zip(fine['runs'],[0.,21600.],[21600.,43200.]):
        assert entry['start_s']==start and entry['end_s']==end and entry['completed']
        assert entry['local_geometry_pass'] and not entry['clearance_event']
        assert entry['mean_scheduled_translation_m_s']==0. and entry['burn_count']==0
        assert entry['comparison']['position_discrepancy_plus_reconstruction_bounds_m']<50.
    assert fine['runs'][0]['coverage']['passed']
    assert len(fine['runs'][0]['coverage']['samples'])==49
    assert not fine['full_return_executed'] and not fine['handover_executed']


def test_exact_distance_and_solar_source_search_preserve_the_constraints():
    audit=read('packing_audit.json')
    assert audit['completed'] and audit['passes'] and not audit['continuous_certificate']
    assert audit['coverage_adversary']['passed']
    assert audit['coverage_adversary']['total_probes']>audit['coverage_adversary']['grid_probes']
    assert audit['coverage_adversary']['worst']['margin_m']>0.
    for row in audit['runs']:
        assert row['minimum_sampled_surface_distance_m']>=row['minimum_conditional_all_pair_clearance_m']>100.
        assert row['intervals_below_100m']==0


def test_energy_and_optical_inventory_do_not_become_recurring_fleet_claims():
    p=read('packing_budget.json');old=read('energy_budget.json')['same_model_totals']['previous']
    assert p['completed'] and p['accepted_local_service_departure']
    total=p['totals'];assert 0.<total['service_plus_departure_m_s']<old['service_plus_departure_m_s']/2
    assert total['same_model_fractional_energy_reduction']==pytest.approx(1-total['propulsion_energy_J']/old['service_plus_departure_energy_J'])
    assert total['total_burn_count']==0 and total['provisional_individual_power_hardware_kg']>0
    for dataset in p['datasets'].values():
        for item in dataset['evaluations']:
            assert len(item['cost']['per_member_peak_propulsion_W'])==361
            assert item['cost']['translation_m_s']==0.
            light=item['optical']
            first=light['first_intercept_bolometric_equivalent_W']['energy_J']
            assert 0<light['redirected_band_optical_W']['energy_J']<first
            assert first==pytest.approx(light['front_first_intercept_W']['energy_J']+light['back_first_intercept_W']['energy_J'])
    assert not p['initial_preparation_reference']['maneuver_executed']
    assert not p['hardware_mass_feedback_propagated'] and not p['recurring_cost_measured']
    for field in ['measured_recurring_impulse_m_s_day','fleet_propulsion_W','collection_capacity_W','delivered_power_W']:
        assert p[field] is None
