"""Keep local simulated success separate from unmeasured fleet recurrence."""
import json
from pathlib import Path

import pytest

from shared.provenance import constants_changed
from .fleet_run import digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def test_patterns_products_bind_sources_and_preserve_evidence_boundaries():
    for name in ['patterns.json', 'patterns_budget.json']:
        p = json.loads((HERE/'results'/name).read_text())
        for path, expected in p['producer']['source_hashes'].items():
            assert digest(ROOT/path) == expected, path
        assert not constants_changed(p['producer']['constants'])
        for name, expected in p['producer'].get('inputs', {}).items():
            assert digest(HERE/'results'/name) == expected
        assert p['accepted_fleet'] is False
        assert p['collection_capacity_W'] is None
        assert p['delivered_power_W'] is None
    p = json.loads((HERE/'results/patterns.json').read_text())
    assert not p['recurring_cost_measured']
    assert not p['attitude']['continuous_certificate']
    assert p['attitude']['input_product_sha256'] == digest(HERE/'results/pilot_validation.json')
    assert len(p['affine_trials']) <= p['budget']['affine_cases']
    assert len(p['coupled']) <= p['budget']['coupled_cases']
    assert p['elapsed_s'] < p['budget']['wall_s']
    for r in p['coupled']+([p['replay']] if 'replay' in p else []):
        assert not r['closure_maneuver_executed'] and not r['recurring_cost_measured']
        assert not r['continuous_coverage_certificate']
        if r['sampled_local_gates_pass']:
            assert r['minimum_coverage'] >= 1-1e-9
            assert r['separation']['possible_conflict_pairs'] == 0
            assert r['guarded_nonstencil_margin_m'] > 0
            assert r['minimum_beam_margin_deg'] > 0
            assert r['max_sampled_acceleration_m_s2'] <= r['separation']['acceleration_bound_m_s2']
    if 'replay' in p:
        assert p['replay']['passes_50m'] == (p['replay']['max_position_difference_m'] < 50.)


def test_service_ledger_does_not_treat_terminal_mismatch_as_closed_maneuver():
    p = json.loads((HERE/'results/patterns_budget.json').read_text())
    assert p['service']['electric_delta_v_m_s'] == 0
    assert p['service']['propulsion_energy_J'] == 0
    assert p['preparation']['electric_energy_J'] > 0
    assert not p['preparation']['maneuver_executed']
    assert not p['endpoint']['closure_maneuver_executed']
    assert p['endpoint']['maximum_position_mismatch_m'] > 100.
    assert p['measured_recurring_delta_v_m_s_per_day'] is None
    assert p['recurring_maneuver_frequency_per_day'] is None
    assert p['recurrence_diagnostic']['velocity_only_two_impulse_proxy_m_s'] == pytest.approx(
        p['preparation']['mean_delta_v_m_s']+p['endpoint']['mean_velocity_mismatch_m_s'])
    assert p['handovers_executed'] == p['return_legs_executed'] == p['recovery_events_executed'] == 0
