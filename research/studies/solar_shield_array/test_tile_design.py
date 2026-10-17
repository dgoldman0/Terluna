"""Evidence/provenance and independent integration checks for the tile study."""
import json

import numpy as np
import pytest

from shared import constants as K
from shared.provenance import constants_changed
from . import tile_design as td


def product():
    return json.loads(td.OUT.read_text())


def test_product_records_exact_producers_inputs_and_constants():
    p = product()
    assert p['schema'] == 'terluna.research.array-tile-design/1'
    assert set(p['producer']['files']) == set(td.FILES)
    for name, digest in p['producer']['files'].items():
        assert td.digest(td.ROOT/name) == digest, name
    for name, record in p['inputs'].items():
        assert td.digest(td.ROOT/name) == record['sha256'], name
        assert json.loads((td.ROOT/name).read_text())['schema'] == record['schema']
    assert not constants_changed(p['producer']['constants'])


def test_ideal_orbit_is_keplerian_and_roll_is_periodic():
    t, states, frames, omega, _ = td.ideal_attitude(20e6)
    speed = np.linalg.norm(states[:, 3:], axis=1)
    assert speed**2 == pytest.approx(K.MOON_GM/20e6)
    assert frames[0] == pytest.approx(frames[-1], abs=1e-12)
    assert np.max(np.linalg.norm(omega, axis=1))*np.max(np.diff(t)) < .03


def test_mass_and_momentum_iteration_is_consistent_but_not_certified():
    p = product()
    for name, v in p['hardware']['variants'].items():
        assert v['validated_design'] is False
        assert sum(v['components_kg'].values()) == pytest.approx(v['total_mass_kg'])
        assert v['total_mass_kg']/v['design']['side_m']**2 == pytest.approx(v['total_areal_mass_kg_m2'])
        assert v['rigid_sizing_iteration']['converged'], name
        assert v['momentum_capacity_utilization_in_screen'] <= 1.001/1.2
        assert v['motor_torque_utilization_in_screen'] <= 1.001/1.2
        assert np.linalg.eigvalsh(v['locked_inertia_kg_m2']).min() > 0
    assert p['hardware']['window_bare_film_excess_over_complete_allocation_g_m2'] > 0
    assert p['checks']['momentum_relative_peak_refinement'] < 3e-4


def test_recovery_does_not_claim_optical_repair_or_late_safety():
    p = product()['local_failures']
    early = p['cases']['backup_after_1800s']
    late = p['cases']['backup_after_3600s']
    assert early['returned_to_slot_by_horizon']
    assert early['first_safety_crossing_s'] is None
    assert not early['filter_function_restored']
    assert late['first_safety_crossing_s'] is not None
    assert late['first_envelope_contact_s'] is None
    assert p['timestep_check']['safety_time_difference_s'] < .1
    assert p['timestep_check']['contact_time_difference_s'] < .1


def test_population_delay_matches_low_hazard_little_law_and_keeps_correlation_explicit():
    p = product()['fleet_service']
    base = p['scenarios']['local_spares_1day']
    # Low-hazard expected unavailable share=lambda * response time. Allow
    # finite population draws, startup transient and the 3h queue timestep.
    expected = .01/K.JULIAN_YEAR_DAYS
    assert base['mean_missing_fraction'] == pytest.approx(expected, rel=.03)
    assert p['scenarios']['plant_delivery_27days']['mean_missing_fraction'] > 20*base['mean_missing_fraction']
    assert p['scenarios']['correlated_10000_slots']['ideal_independent_peak_missing_ray_fraction'] is None
    assert p['scenarios']['capacity_below_renewal']['final_renewal_waiting'] > 0
    assert all(all(r == 0 for r in s['balance_residual'].values()) for s in p['scenarios'].values())


def test_industry_comparison_charges_both_planned_and_failed_hardware():
    p = product()
    i = p['industry_comparison']
    s = p['fleet_service']['scenarios']['local_spares_1day']
    assert i['inventory_if_copied_Gt'] == pytest.approx(p['fleet_service']['population']*i['mean_reference_tile_mass_kg']/1e12)
    yearly_jobs = (s['planned_per_day']+s['expected_failure_per_day_at_full_population'])*K.JULIAN_YEAR_DAYS
    assert i['whole_tile_exchange_Gt_per_year'] == pytest.approx(yearly_jobs*i['mean_reference_tile_mass_kg']/1e12)
