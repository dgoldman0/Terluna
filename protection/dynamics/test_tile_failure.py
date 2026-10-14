from dataclasses import replace
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from shared.constants import MOON_GM
from protection.dynamics.tile_failure import (FailureConfig, free_hcw_state,
    hcw_acceleration, parallel_box_clearances, rotation_bound_window_s,
    rotation_clearance_lower_bound, simulate_failure, tile_frame)


def test_unforced_integration_matches_closed_form_hcw():
    n = math.sqrt(MOON_GM/(20e6)**3)
    initial = np.array([100., -50., 23., .002, -.003, .001])
    sol = solve_ivp(lambda t, y: np.r_[y[3:], hcw_acceleration(y, n)],
                    (0., 20_000.), initial, method="DOP853", rtol=1e-12,
                    atol=1e-10, dense_output=True, max_step=100.)
    times = np.linspace(0., 20_000., 101)
    assert np.allclose(sol.sol(times).T, free_hcw_state(initial, times, n), atol=1e-8, rtol=1e-10)


def test_constant_normal_force_matches_independent_oscillator_solution():
    force, duration = 3e-5, 1800.
    cfg = FailureConfig(duration_s=duration, residual_acceleration_m_s2=(0., 0., force), ring_step_m=None)
    result = simulate_failure(cfg)
    n = math.sqrt(MOON_GM/cfg.radius_m**3)
    expected_z = force/n**2*(1-math.cos(n*duration))
    expected_v = force/n*math.sin(n*duration)
    assert np.isclose(result["final_state_rtn_m_m_s"][2], expected_z, atol=1e-8)
    assert np.isclose(result["final_state_rtn_m_m_s"][5], expected_v, atol=1e-10)


def test_finite_square_edges_and_shingle_clearance_are_not_centre_distances():
    frame = tile_frame(1.5)
    clearance = parallel_box_clearances([0., 0., 0.], [[0., 8500., 0.]], frame, 10_000., 1.)[0]
    assert np.isclose(clearance, 8500*math.sin(math.radians(1.5))-1.)
    # A neighbour translated just beyond a square corner is farther than its
    # normal-plane separation; both finite in-plane dimensions enter.
    offset = frame @ np.array([10_003., 10_004., 13.])
    assert np.isclose(parallel_box_clearances([0., 0., 0.], [offset], frame, 10_000., 1.)[0], 13.)


def test_permanent_residual_stops_at_contact_and_has_converged_first_crossings():
    coarse = simulate_failure(FailureConfig(max_step_s=120.))
    fine = simulate_failure(FailureConfig(max_step_s=15.))
    assert coarse["first_safety_crossing_s"] > 0
    assert coarse["first_envelope_contact_s"] > coarse["first_safety_crossing_s"]
    assert coarse["physical_prefix_end_s"] == coarse["first_envelope_contact_s"]
    assert coarse["status"] == "envelope_contact"
    assert coarse["minimum_sampled_surface_clearance_m"] < 1e-6
    assert coarse["maximum_sampled_relative_speed_m_s"] > 0
    for key in ("first_safety_crossing_s", "first_envelope_contact_s"):
        assert abs(coarse[key]-fine[key]) < 1e-5
    assert abs(coarse["maximum_sampled_displacement_m"]-fine["maximum_sampled_displacement_m"]) < 1e-5


def test_thin_contact_between_output_steps_is_detected_by_face_roots():
    cfg = FailureConfig(shingle_deg=0., ring_step_m=None, envelope_thickness_m=1e-5,
                        duration_s=100., initial_state=(0., 0., 0., 1., 0., 0.),
                        residual_acceleration_m_s2=(0., 0., 0.), max_step_s=100.)
    result = simulate_failure(cfg, neighbour_offsets_m=[[20., 0., 0.]])
    assert result["status"] == "envelope_contact"
    assert abs(result["first_envelope_contact_s"]-20.) < .001


def test_early_backup_returns_continuously_and_late_backup_cannot_erase_collision():
    base = FailureConfig(response_delay_s=600., backup_acceleration_cap_m_s2=1e-4)
    early = simulate_failure(base)
    assert early["returned_to_slot_by_horizon"]
    assert early["first_safety_crossing_s"] is None
    assert early["backup_integrated_acceleration_m_s"] > 0
    assert early["maximum_sampled_backup_acceleration_m_s2"] <= base.backup_acceleration_cap_m_s2*(1+1e-12)
    assert early["maximum_sampled_displacement_m"] > 1.
    late = simulate_failure(replace(base, response_delay_s=10_000.))
    assert late["status"] == "envelope_contact"
    assert late["backup_integrated_acceleration_m_s"] == 0.
    assert not late["returned_to_slot_by_horizon"]
    assert not early["filter_function_restored"]


def test_underpowered_backup_does_not_get_unlimited_control_credit():
    cfg = FailureConfig(response_delay_s=0., backup_acceleration_cap_m_s2=1e-7)
    result = simulate_failure(cfg)
    assert result["status"] == "envelope_contact"
    assert result["maximum_sampled_backup_acceleration_m_s2"] <= 1.0000001e-7


def test_zero_residual_and_no_initial_error_stays_in_slot_without_control():
    result = simulate_failure(FailureConfig(residual_acceleration_m_s2=(0., 0., 0.), duration_s=600.))
    assert result["status"] == "clear_through_horizon"
    assert result["maximum_sampled_displacement_m"] == 0.
    assert result["backup_integrated_acceleration_m_s"] == 0.
    assert not result["returned_to_slot_by_horizon"]


def test_tumbling_clearance_is_a_conservative_bound_not_a_collision_claim():
    initial, gate, side, rate = 221.5, 150., 10_000., 1e-5
    window = rotation_bound_window_s(initial, gate, side, rate)
    assert window > 0.
    assert np.isclose(rotation_clearance_lower_bound(initial, side, rate, window), gate)
    assert rotation_bound_window_s(initial, gate, side, 0.) is None
    assert rotation_bound_window_s(initial, initial, side, rate) == 0.
    assert rotation_clearance_lower_bound(initial, side, rate, 2*window) < gate
    # Finite hardware corners are farther from the centre than a zero-depth sheet.
    thick = rotation_bound_window_s(initial, gate, side, rate, 1000.)
    assert thick < window
    assert np.isclose(rotation_clearance_lower_bound(initial, side, rate, thick, 1000.), gate)


def test_invalid_and_nonlocal_inputs_are_rejected():
    with pytest.raises(ValueError):
        simulate_failure(FailureConfig(duration_s=1e7))
    with pytest.raises(ValueError):
        simulate_failure(FailureConfig(), neighbour_offsets_m=[])
    with pytest.raises(ValueError):
        simulate_failure(FailureConfig(residual_acceleration_m_s2=(float("nan"), 0., 0.)))
