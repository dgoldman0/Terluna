"""Physical invariants and limiting cases of the conditional spinning disk."""
from math import pi, sqrt

import numpy as np
import pytest

from engineering.array_spin import ballasted_disk_screen, disk_grid_coverage


def case(**changes):
    kwargs = dict(radius_m=100., areal_mass_kg_m2=.01, spin_rad_s=.01,
                  solar_pressure_Pa=1e-6, allowable_sag_m=.1,
                  load_bearing_thickness_m=1e-5, admissible_stress_Pa=1e9,
                  times_s=[0., 1., 2.], omega_rad_s=np.zeros((3, 3)),
                  alpha_rad_s2=np.zeros((3, 3)),
                  radial_body=np.tile([0., 0., 1.], (3, 1)),
                  mean_motion_squared_s2=0.)
    kwargs.update(changes)
    return ballasted_disk_screen(**kwargs)


def test_finite_ballast_reacts_exact_edge_traction():
    p = case()
    inward_force_per_length = p['rim_line_density_kg_m']*.01**2*100
    assert inward_force_per_length == pytest.approx(p['edge_prestress_N_m'])
    assert p['ideal_film_plus_ballast_mass_kg'] == pytest.approx(
        p['film_mass_kg']+p['ideal_rim_ballast_mass_kg'])


def test_normal_spin_has_no_transverse_rigid_acceleration():
    p = case()
    assert p['peak_normal_acceleration_m_s2'] == 0
    assert p['peak_spin_precession_torque_N_m'] == 0
    assert p['peak_full_rigid_external_torque_N_m'] == 0


def test_spin_reduces_static_sag_at_fixed_prescribed_prestress():
    slow = case(spin_rad_s=1e-7)
    fast = case(spin_rad_s=.01)
    assert fast['edge_prestress_N_m'] == slow['edge_prestress_N_m']
    assert fast['exact_axisymmetric_static_sag_m'] < slow['exact_axisymmetric_static_sag_m']
    assert slow['exact_axisymmetric_static_sag_m'] == pytest.approx(
        slow['constant_tension_static_sag_bound_m'], rel=1e-7)


def test_axisymmetric_sag_matches_direct_radial_force_balance():
    p = case()
    radius = np.linspace(0, 100, 100001)
    a = (3+.3)*.01*.01**2/8
    tr = p['edge_prestress_N_m']+a*(100**2-radius**2)
    # Equilibrium integrates r*Tr*w'=-q*r^2/2 from the disk centre.
    sag = np.trapezoid(1e-6*radius/(2*tr), radius)
    assert sag == pytest.approx(p['exact_axisymmetric_static_sag_m'], rel=2e-8)


def test_roll_spin_coupling_equals_two_s_omega_not_one():
    omega = np.tile([.001, 0., 0.], (3, 1))
    p = case(omega_rad_s=omega)
    assert p['peak_normal_acceleration_m_s2'] == pytest.approx(100*2*.01*.001)
    assert p['peak_spin_precession_torque_N_m'] == pytest.approx(
        p['spin_angular_momentum_N_m_s']*.001)
    assert p['peak_full_rigid_external_torque_N_m'] == pytest.approx(
        p['peak_spin_precession_torque_N_m'])


def test_command_axial_rotation_changes_ballast_and_stress_concurrently():
    baseline = case()
    w = np.array([[0., 0., -.005], [0., 0., 0.], [0., 0., .005]])
    p = case(omega_rad_s=w)
    assert p['ideal_rim_ballast_mass_kg'] == pytest.approx(4*baseline['ideal_rim_ballast_mass_kg'])
    assert p['minimum_actual_axial_spin_rad_s'] == pytest.approx(.005)
    assert p['maximum_actual_axial_spin_rad_s'] == pytest.approx(.015)
    assert p['minimum_edge_prestress_N_m'] == pytest.approx(baseline['edge_prestress_N_m'])
    assert p['edge_prestress_N_m'] == pytest.approx(9*baseline['edge_prestress_N_m'])


def test_command_cannot_silently_pass_through_zero_spin():
    with pytest.raises(ValueError):
        case(omega_rad_s=np.tile([0., 0., -.01], (3, 1)))


def test_acceleration_and_tidal_gradient_can_cancel_before_norm():
    n2 = 1e-6
    u = np.tile([1/sqrt(2), 0, 1/sqrt(2)], (3, 1))
    alpha = np.tile([0, -1.5*n2, 0], (3, 1))
    p = case(radial_body=u, alpha_rad_s2=alpha, mean_motion_squared_s2=n2)
    assert p['peak_normal_acceleration_m_s2'] < 1e-18


def test_pressure_and_motion_use_concurrent_samples():
    alpha = np.array([[0., 0., 0.], [1e-5, 0., 0.], [0., 0., 0.]])
    p = case(alpha_rad_s2=alpha, solar_pressure_Pa=[2e-5, 0., 0.])
    assert p['peak_effective_pressure_Pa'] == pytest.approx(2e-5)
    assert p['pressure_peak_time_s'] == 0.


def test_stress_charges_only_load_bearing_thickness():
    p = case(load_bearing_thickness_m=1e-12)
    assert not p['static_stress_screen_passes']
    assert p['peak_membrane_stress_Pa'] == pytest.approx(
        p['peak_tension_resultant_N_m']/1e-12)
    assert not p['accepted_complete_design']


def test_filled_disk_and_rim_inertia():
    p = case()
    inertia = np.array(p['locked_axisymmetric_inertia_kg_m2'])
    assert inertia[2, 2] == pytest.approx(
        (.5*p['film_mass_kg']+p['ideal_rim_ballast_mass_kg'])*100**2)
    assert inertia[0, 0]+inertia[1, 1] == pytest.approx(inertia[2, 2])
    assert p['spin_kinetic_energy_J'] == pytest.approx(
        .5*p['spin_angular_momentum_N_m_s']*.01)


def test_gravity_gradient_is_external_torque_credit():
    u = np.tile([1/sqrt(2), 0, 1/sqrt(2)], (3, 1))
    p = case(radial_body=u, mean_motion_squared_s2=1e-6)
    assert p['peak_full_rigid_external_torque_N_m'] == 0
    inertia = p['locked_axisymmetric_inertia_kg_m2']
    expected = 1.5e-6*(inertia[2][2]-inertia[0][0])
    assert p['peak_full_rigid_nongravitational_torque_N_m'] == pytest.approx(expected)


def test_square_grid_circle_covering_has_pi_over_two_area_cost():
    p = disk_grid_coverage(1., sqrt(2), sqrt(2))
    assert p['original_grid_fully_covered']
    assert p['gross_film_area_per_covered_cell'] == pytest.approx(pi/2)


def test_smaller_disk_increases_count_for_equal_protected_area():
    p = disk_grid_coverage(1., 2., 2.)
    assert not p['original_grid_fully_covered']
    assert p['count_multiplier_at_fixed_protected_area'] == pytest.approx(2.)
    assert p['gross_film_area_per_covered_cell'] == pytest.approx(pi/2)


@pytest.mark.parametrize('key,value', [
    ('spin_rad_s', 0.), ('radius_m', -1.), ('areal_mass_kg_m2', float('nan')),
    ('load_bearing_thickness_m', 0.), ('poisson_ratio', .5),
    ('dynamic_factor', .9), ('mean_motion_squared_s2', -1.),
    ('times_s', [0., 0., 1.]), ('solar_pressure_Pa', [1., -1., 1.]),
    ('radial_body', np.zeros((3, 3))),
])
def test_invalid_inputs_are_rejected(key, value):
    with pytest.raises(ValueError):
        case(**{key: value})


def test_zero_residual_load_is_not_used_to_assert_taut_zero_ballast():
    with pytest.raises(ValueError):
        case(solar_pressure_Pa=0.)
