"""Analytic limits and conservation checks for the conditional tile screen."""
from dataclasses import replace
from math import sqrt

import pytest

from engineering.array_tile import (
    Material, TileDesign, hoop_store, membrane_coefficient, perimeter_truss,
    square_membrane, tile_budget,
)


def candidate(**changes):
    values = dict(side_m=10000., film_areal_mass_kg_m2=.01021,
                  target_total_areal_mass_kg_m2=.01307, pressure_Pa=6.3e-7,
                  angular_acceleration_rad_s2=2e-8,
                  required_momentum_N_m_s=1.5e10, required_torque_N_m=1e6,
                  orbit_period_s=250000.)
    values.update(changes)
    return TileDesign(**values)


def test_membrane_series_converges_to_square_poisson_solution():
    # Resolved square Poisson series, with the single (1,1) mode above the
    # converged centre displacement; this is numerical convergence, not data.
    c = membrane_coefficient(31)
    assert c == pytest.approx(.07367, abs=3e-6)
    assert abs(c-membrane_coefficient(61)) < 1e-6
    assert membrane_coefficient(1) > c


def test_membrane_scaling_and_zero_load():
    base = square_membrane(1000., 1e-6, 1., .01)
    twice = square_membrane(2000., 1e-6, 1., .01)
    assert twice['centre_sag_m'] == pytest.approx(4*base['centre_sag_m'])
    assert twice['fundamental_frequency_Hz'] == pytest.approx(base['fundamental_frequency_Hz']/2)
    higher = square_membrane(1000., 1e-6, 4., .01)
    assert higher['centre_sag_m'] == pytest.approx(base['centre_sag_m']/4)
    assert higher['fundamental_frequency_Hz'] == pytest.approx(2*base['fundamental_frequency_Hz'])
    zero = square_membrane(1000., 0., 1., .01)
    assert zero['centre_sag_m'] == 0
    assert zero['fundamental_frequency_Hz'] == pytest.approx(sqrt(100)/(sqrt(2)*1000))


def test_frame_reactions_and_stiffness_close_explicit_load_path():
    material = Material()
    result = perimeter_truss(10000., 5., material=material)
    assert result['side_axial_compression_N'] == 25000.
    assert result['max_bending_moment_N_m'] == 5*10000**2/8
    assert result['chord_compression_N'] == pytest.approx(5*10000**2/(8*200)+25000/4)
    assert sum(result['components_kg'].values()) == pytest.approx(result['mass_kg'])
    assert result['ideal_edge_bow_m'] <= 5.*(1+1e-12)
    assert result['edge_bow_bending_m'] == pytest.approx(2.5)
    # Stiffness sizing must exceed the misleading strength-only frame.
    assert result['chord_area_m2'] > result['chord_strength_and_buckling_area_m2']
    larger_allowance = perimeter_truss(10000., 5., allowable_edge_bow_m=50.)
    assert larger_allowance['mass_kg'] < result['mass_kg']


def test_frame_scales_linearly_with_membrane_tension():
    a = perimeter_truss(1000., 1., depth_m=20., bay_m=20., chord_radius_m=.2)
    b = perimeter_truss(1000., 3., depth_m=20., bay_m=20., chord_radius_m=.2)
    assert b['mass_kg'] == pytest.approx(3*a['mass_kg'])
    # A positive joint and support inventory exists, even for weak loads.
    assert all(x > 0 for x in a['components_kg'].values())


def test_hoop_conserves_momentum_energy_and_uses_three_axes():
    h, r, v = 1.5e10, 100., 1000.
    result = hoop_store(h, r, v)
    m = result['rim_mass_per_axis_kg']
    omega = result['angular_speed_rad_s']
    assert m*r*v == pytest.approx(h)
    assert result['axial_inertia_per_axis_kg_m2']*omega == pytest.approx(h)
    assert result['stored_energy_per_axis_J'] == pytest.approx(h*omega/2)
    assert result['total_rim_mass_kg'] == 3*m
    assert result['hoops'] == 6
    assert result['hoop_stress_Pa'] == 1600*v*v
    assert not hoop_store(h, r, v, admissible_stress_Pa=1e9)['stress_screen_passes']
    bigger = hoop_store(h, 2*r, v)
    assert bigger['total_rim_mass_kg'] == result['total_rim_mass_kg']/2
    assert bigger['total_maximum_stored_energy_J'] == result['total_maximum_stored_energy_J']/2


def test_zero_momentum_requires_no_rotor_energy_or_rim_mass():
    result = hoop_store(0., 100., 1000.)
    assert result['total_mass_kg'] == 0
    assert result['total_maximum_stored_energy_J'] == 0


def test_mass_ledger_and_motor_bus_energy_conserve():
    d = candidate()
    result = tile_budget(d)
    assert sum(result['components_kg'].values()) == pytest.approx(result['total_mass_kg'])
    assert result['mass_margin_kg'] == pytest.approx(d.side_m**2*d.target_total_areal_mass_kg_m2-result['total_mass_kg'])
    power = result['power']
    mechanical = result['rotors']['total_maximum_stored_energy_J']
    trim = result['optical_trim']['eclipse_device_energy_J']
    available = power['battery_nominal_J']*d.battery_depth_of_discharge
    assert available == pytest.approx(((mechanical+trim)/d.converter_efficiency+d.housekeeping_W*d.eclipse_s)/(1-d.wire_loss_fraction))
    assert power['peak_electrical_W']-power['peak_mechanical_W'] == pytest.approx(power['peak_waste_heat_W'])
    area = power['branch_wire_area_m2']
    resistance = d.wire_resistivity_ohm_m*d.side_m/area
    current = power['peak_electrical_W']/6/d.bus_voltage_V
    assert 6*current**2*resistance == pytest.approx(power['peak_electrical_W']*d.wire_loss_fraction)


def test_trim_switch_energy_enters_storage_once_and_respects_power_envelope():
    d = candidate()
    result = tile_budget(d)
    no_switch = tile_budget(replace(d, trim_switches_per_orbit=0))
    trim = result['optical_trim']
    energy = d.side_m**2*d.trim_area_fraction*d.trim_switch_J_m2*d.trim_switches_per_orbit
    assert trim['switch_energy_per_orbit_J'] == energy
    assert trim['average_switch_power_W'] == energy/d.orbit_period_s
    assert trim['peak_switch_power_W']*d.eclipse_s >= energy
    battery_delta = result['power']['battery_nominal_J']-no_switch['power']['battery_nominal_J']
    assert battery_delta == pytest.approx(energy/(d.converter_efficiency*(1-d.wire_loss_fraction)*d.battery_depth_of_discharge))
    # Recharging covers this same energy once; no extra orbit-average term.
    generation_delta = result['power']['generation_W']-no_switch['power']['generation_W']
    assert generation_delta == pytest.approx(energy/(d.recharge_s*d.converter_efficiency**2*(1-d.wire_loss_fraction)**2))
    assert trim['qualified_electrical_device'] is False


def test_nonzero_trim_hold_power_has_continuous_eclipse_and_thermal_cost():
    d = candidate()
    optimistic = tile_budget(d)
    held = tile_budget(replace(d, trim_hold_W_m2=.01))
    holding = d.side_m**2*d.trim_area_fraction*.01
    assert held['optical_trim']['continuous_hold_power_W'] == holding
    assert held['optical_trim']['eclipse_device_energy_J']-optimistic['optical_trim']['eclipse_device_energy_J'] == holding*d.eclipse_s
    assert held['power']['peak_waste_heat_W']-optimistic['power']['peak_waste_heat_W'] == pytest.approx(holding/(d.converter_efficiency*(1-d.wire_loss_fraction)))
    assert held['power']['generation_W'] > optimistic['power']['generation_W']
    assert held['total_mass_kg'] > optimistic['total_mass_kg']


def test_locked_inertia_preserves_products_of_inertia_and_positive_energy():
    result = tile_budget(candidate())
    tensor = result['locked_inertia_kg_m2']
    assert tensor[0][1] == tensor[1][0] < 0
    assert tensor[0][0] == tensor[1][1]
    assert tensor[0][0] > abs(tensor[0][1])
    assert tensor[2][2] > 0
    # Increasing centered mass through a perimeter allocation raises inertia.
    heavy = tile_budget(candidate(avionics_and_sensors_kg=50000.))
    assert heavy['locked_inertia_kg_m2'][0][0] > tensor[0][0]


def test_rotation_and_local_trim_loading_raise_structural_demand():
    d = candidate()
    still = tile_budget(replace(d, angular_acceleration_rad_s2=0))
    turning = tile_budget(d)
    no_trim = tile_budget(replace(d, trim_area_fraction=0))
    combined = tile_budget(replace(d, angular_speed_rad_s=.0002))
    assert turning['membrane']['tension_N_m'] > still['membrane']['tension_N_m']
    assert turning['membrane']['tension_N_m'] > no_trim['membrane']['tension_N_m']
    assert combined['membrane']['tension_N_m'] > turning['membrane']['tension_N_m']
    assert turning['total_mass_kg'] > turning['target_mass_kg']
    assert turning['validated_design'] is False


@pytest.mark.parametrize('kwargs', [dict(side_m=-1), dict(pressure_Pa=float('nan')),
                                   dict(converter_efficiency=0), dict(rotor_radius_m=1100),
                                   dict(trim_area_fraction=1.1), dict(wire_loss_fraction=0)])
def test_invalid_design_inputs_rejected(kwargs):
    with pytest.raises(ValueError):
        candidate(**kwargs)


def test_invalid_primitives_rejected():
    with pytest.raises(ValueError):
        Material(density_kg_m3=0)
    with pytest.raises(ValueError):
        square_membrane(1000, -1, 1, .01)
    with pytest.raises(ValueError):
        hoop_store(1, 100, 1000, axes=0)
    with pytest.raises(ValueError):
        membrane_coefficient(0)
