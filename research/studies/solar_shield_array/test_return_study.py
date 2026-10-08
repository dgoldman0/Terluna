"""The return obstruction and executed-prefix ledger keep their evidence scope."""
import json
from pathlib import Path

import pytest

from shared.provenance import constants_changed
from .fleet_run import digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read(name):
    return json.loads((HERE/'results'/name).read_text())


def test_return_products_bind_sources_constants_and_inputs():
    for name in ['return_screen.json', 'return_prefix.json', 'return_budget.json']:
        p = read(name)
        for path, expected in p['producer']['source_hashes'].items():
            assert digest(ROOT/path) == expected, path
        assert not constants_changed(p['producer']['constants'])
        for name, expected in p['producer']['inputs'].items():
            assert digest(HERE/'results'/name) == expected
        assert p['accepted_fleet'] is False


def test_low_endpoint_impulse_does_not_erase_rejected_path_constraints():
    p = read('return_screen.json')
    assert len(p['cases']) == 24
    assert not p['accepted_return'] and not p['recurring_cost_measured']
    for c in p['cases']:
        if c['shooting']['success']:
            assert c['encounters']['sampled_pair_violations'] > 0
            assert c['beam_and_rate']['maximum_sampled_rate_deg_s'] > .1
    failed = p['refinements'][-1]
    assert failed['between_samples']['physical_intersections'] > 0
    assert not failed['correction']['success']
    assert failed['correction']['status'] == 2
    assert failed['correction']['cuts'] == 512
    assert p['proposal_mean_delta_v_m_s'] < 1.


def test_coupled_clearance_stop_repeats_and_ledger_counts_only_executed_prefix():
    p = read('return_prefix.json'); b = read('return_budget.json')
    assert len(p['runs']) == 2
    for r in p['runs']:
        assert r['clearance_event']
        assert r['minimum_clearance_m'] == pytest.approx(100., abs=1e-5)
        assert r['beam_and_rate']['minimum_sampled_beam_margin_deg'] > 0
        assert r['beam_and_rate']['maximum_sampled_rate_deg_s'] < .1
        assert r['max_acceleration_command_m_s2'] < .001
    assert p['independent_replay']['same_limiting_pair']
    assert p['independent_replay']['max_position_difference_m'] < .02
    assert not p['full_return_executed'] and not p['accepted_handover']
    assert b['prefix_between_sample_crossing_audit']['physical_intersections'] == 0
    assert b['frozen_state_diagnostic']['sun_facing_plane_minimum_clearance_m'] > 2000
    assert b['frozen_state_diagnostic']['failed_pair_centre_distance_m'] < 10000
    assert b['executed_prefix']['propulsion_energy_J'] > 0
    assert b['executed_prefix']['thrust_arcs_completed'] == 0
    assert b['executed_prefix']['mean_delta_v_m_s'] < b['endpoint_only_proposal']['mean_delta_v_m_s']
    assert not b['full_return_executed'] and not b['recovery_executed']
    assert b['handovers_executed'] == 0
    for field in ['measured_recurring_delta_v_m_s_per_day', 'fleet_propulsion_W',
                  'collection_capacity_W', 'delivered_power_W']:
        assert b[field] is None
