"""Departure evidence preserves failed reductions and an incomplete fleet ledger."""
import json
from pathlib import Path

from shared.provenance import constants_changed
from .fleet_run import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(name):
    return json.loads((HERE/'results'/name).read_text())


def test_departure_products_bind_the_actual_producers_and_inputs():
    for name in ['departure_screen.json','departure_refine.json','departure.json',
                 'departure_validation.json','departure_budget.json']:
        product=read(name)
        for path,expected in product['producer']['source_hashes'].items():
            assert digest(ROOT/path)==expected,path
        assert not constants_changed(product['producer']['constants'])
        for path,expected in product['producer']['inputs'].items():
            assert digest(HERE/'results'/path)==expected,path
        assert product['accepted_fleet'] is False


def test_depth_only_failures_are_retained_and_refined_intervals_are_distinct():
    screen=read('departure_screen.json'); refinement=read('departure_refine.json')
    coarse=read('departure.json'); validation=read('departure_validation.json')
    assert len(screen['cases'])==96 and not screen['survivors']
    assert sum(c['rejected_constraint']=='acceleration' for c in screen['cases'])==24
    assert all('clearance' in c['rejected_constraint'] for c in screen['cases'] if c['rejected_constraint']!='acceleration')
    assert refinement['coupled_candidate']
    assert refinement['iterations'][-1]['optimization']['success']
    assert not coarse['accepted_departure']
    assert coarse['runs'][0]['completed_prefix']
    assert coarse['runs'][0]['minimum_sampled_clearance_m']>100
    assert coarse['runs'][0]['swept_separation']['possible_conflict_pairs']>0
    assert validation['accepted_departure']
    for geometry in [validation['coarse_two_second_audit'],validation['runs'][0]['geometry']]:
        assert geometry['step_s']<=2.
        assert geometry['swept_separation']['possible_conflict_pairs']==0
        assert geometry['swept_separation']['minimum_proven_axis_clearance_m']>100.
    assert validation['independent_replay']['maximum_position_difference_m']<50.
    assert not validation['continuous_certificate']
    fine=validation['runs'][0]
    assert not fine['clearance_event']
    assert fine['between_sample_crossings']['physical_intersections']==0
    assert fine['beam_and_rate']['minimum_sampled_beam_margin_deg']>0
    assert fine['beam_and_rate']['maximum_sampled_rate_deg_s']<.1
    assert fine['maximum_thrust_acceleration_m_s2']<.001
    inside=[r for r in fine['extra_moving_window_samples'] if r['window_inside_target']]
    assert len(inside)==4 and all(r['minimum_source_coverage']>=1-1e-10 for r in inside)
    assert not fine['extra_moving_window_samples'][-1]['window_inside_target']


def test_actuation_cost_includes_rotation_but_does_not_invent_recurrence():
    p=read('departure_budget.json')
    assert p['nominal_attitude_equivalent_delta_v_m_s']>p['translation_mean_delta_v_m_s']
    assert p['combined_mean_equivalent_delta_v_m_s']>4.
    assert p['service_plus_departure_equivalent_delta_v_m_s']>p['combined_mean_equivalent_delta_v_m_s']
    assert p['provisional_peak_power_hardware_kg']>0
    assert not p['mass_feedback_propagated'] and not p['structural_realizability_demonstrated']
    assert not p['full_return_executed'] and p['handovers_executed']==0
    assert p['derivative_convention_residual_rad_s']<1e-9
    for field in ['measured_recurring_delta_v_m_s_per_day','fleet_propulsion_W',
                  'collection_capacity_W','delivered_power_W']:
        assert p[field] is None
