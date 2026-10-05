"""Keep local savings, failed returns and supply evidence distinct."""
import json
from pathlib import Path

import numpy as np

from shared.provenance import constants_changed
from .fleet_run import digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read(name):
    return json.loads((HERE / 'results' / (name + '.json')).read_text())


def checks():
    return json.loads((HERE / 'natural_checks.json').read_text())


def test_natural_products_bind_sources_parents_and_constants():
    manifest = checks()
    assert len(manifest['products']) == 16
    for entry in manifest['products']:
        path = ROOT / entry['path']
        assert digest(path) == entry['sha256']
        data = json.loads(path.read_text())
        assert data['completed'] == (path.name != 'natural_scout.json')
        for name, expected in data['producer']['source_hashes'].items():
            assert digest(ROOT / name) == expected, (path, name)
        for name, expected in data['producer']['inputs'].items():
            assert digest(HERE / 'results' / name) == expected, (path, name)
        assert not constants_changed(data['producer']['constants'])
        assert not data.get('accepted_return', False)
        assert not data.get('accepted_cycle', False)
    for entry in manifest['unchanged_pinned_requirements']:
        assert digest(ROOT / entry['path']) == entry['sha256']
        assert all(entry['unchanged_since'].values())


def test_slower_turn_passes_local_prefix_but_fails_return_with_uncertainty():
    fine = read('natural_free_slow_fine')
    audit = read('natural_audit')
    prefix = audit['first_twelve_hours']
    assert prefix['passed_local_prefix'] and prefix['zero_electric_translation']
    assert prefix['replayed_position_difference_m'] < 50
    assert prefix['clearance']['minimum_conditional_all_pair_clearance_m'] > 100
    assert .30 < prefix['energy_reduction_fraction'] < .31
    assert fine['runs'][0]['sampled_coverage']['passed']
    assert len(fine['runs']) == 2
    stop = fine['runs'][1]
    assert stop['clearance_event'] and not stop['completed']
    assert stop['witness']['pair'] == [325, 347]
    assert abs(stop['end_s'] / 3600 - 13.799825) < .00001
    assert fine['comparison'][1]['conditional_clearance_after_numerical_allowance_m'] < 100
    assert not any(row['actuator_event'] for row in fine['runs'])
    assert not audit['next_service_executed'] and audit['handovers_executed'] == 0
    assert audit['replacement_service_required_from_s'] == 21600


def test_false_positive_and_restricted_failures_remain_visible():
    face = read('natural_face_refined_coarse')
    assert not face['runs'][0]['sampled_coverage']['passed']
    assert face['runs'][1]['clearance_event']
    dynamic = read('natural_dynamic_rays')['cases']
    assert len(dynamic) == 4
    assert dynamic[0]['raw_sha256'] == dynamic[1]['raw_sha256']
    assert checks()['duplicate_variant']['distinct_variants'] == 3
    assert dynamic[0]['iterations'][-1]['fit']['success']
    assert dynamic[0]['iterations'][-1]['geometry']['violating_samples'] > 0
    refined = read('natural_refine')
    assert not refined['iterations'][0]['fit']['success']
    assert refined['final_next_service']['minimum_source_coverage'] < .3
    for row in read('natural_groups')['cases']:
        assert row['first_clearance_failure']['distance_m'] < 100
    force = read('natural_force_audit')['stages'][-1]
    assert force['maximum_reduced_vs_linear_proposal_m'] < 50
    assert force['maximum_coupled_vs_nonlinear_reduced_m'] > 5000


def test_prefix_cost_light_and_exhaust_are_complete_without_supply_credit():
    account = read('natural_account')
    for key in ['actual_generation_W', 'actual_bus_delivery_W',
                'actual_storage_state_J', 'full_cycle_energy_J']:
        assert account[key] is None
    for row in account['cases']:
        stages = row['stages']
        assert stages[0]['start_s'] == 0
        assert [s['start_s'] for s in stages][1:] == [s['end_s'] for s in stages][:-1]
        energy = sum(s['evaluations'][-1]['electrical_energy_J'] for s in stages)
        light = sum(s['evaluations'][-1]['first_intercept_J'] for s in stages)
        np.testing.assert_allclose(energy, row['executed_electrical_energy_J'], rtol=1e-14)
        np.testing.assert_allclose(light, row['first_intercept_J'], rtol=1e-14)
        np.testing.assert_allclose(row['ideal_exhaust_kg'], energy * 1.4 / 30000**2)
        np.testing.assert_allclose(row['executed_mean_power_W'] * row['duration_s'], energy)
        assert 0 < row['redirected_band_J'] < light
        assert row['overloaded_members'] == 0
        assert row['minimum_installed_to_peak'] > 1
        assert not row['complete_return'] and not row['actual_supply_validated']
    slow = next(x for x in account['cases'] if x['case'] == 'free_slow')
    assert slow['translation_mean_m_s'] == 0
    assert slow['first_twelve_hours_energy_J'] < slow['executed_electrical_energy_J']


def test_resource_bound_charges_interruption_once_and_long_coasts_need_inventory():
    manifest = checks()
    budget = manifest['budget']
    completed = sum(p['resources']['cpu_s'] for p in manifest['products'] if p['completed'])
    assert budget['completed_producer_cpu_s'] == completed
    assert budget['interrupted_attempt_charged_cpu_s'] == 400
    assert budget['total_charged_cpu_s'] == completed + 400 + budget['charged_ancillary_allowance_s']
    assert budget['total_charged_cpu_s'] <= budget['cpu_limit_s']
    assert budget['peak_recorded_rss_MiB'] < budget['address_space_limit_MiB']
    assert budget['new_raw_MiB'] < budget['new_raw_limit_MiB']
    assert not manifest['interrupted_attempt']['completed']
    assert digest(HERE / 'natural_scout.py') == manifest['interrupted_attempt']['producer_sha256']
    for case in read('natural_fast_scout')['cases']:
        cadences = [row['cadence'] for row in case['opportunities']]
        assert len(cadences) == 4
        for cadence in cadences:
            multiple = cadence['start_to_start_cadence_s'] / 21600
            np.testing.assert_allclose(cadence['time_capacity_group_lower_bound'], multiple)
            np.testing.assert_allclose(cadence['time_capacity_tile_lower_bound'], multiple * 361)
            assert not cadence['spatial_coverage_and_handover_demonstrated']
        assert cadences[-1]['time_capacity_mass_lower_bound_kg'] > 3.9 * cadences[0]['time_capacity_mass_lower_bound_kg']
