import numpy as np
import pytest
from engineering.large_coils import winding, upstream_width


def test_hoop_load_is_energy_derivative_at_fixed_current_and_section():
    a = 500e3; current = 3e8; delta = 10.
    w = winding(a, current, bundle_radius=20.)
    derivative = (winding(a+delta, current, bundle_radius=20.)['self_energy_J']-
                  winding(a-delta, current, bundle_radius=20.)['self_energy_J'])/(2*delta)
    assert w['hoop_tension_N'] == pytest.approx(derivative/(2*np.pi), rel=1e-8)


def test_smaller_current_changes_conductor_support_and_cold_budget():
    a = winding(500e3, 1e8, bundle_radius=100.)
    b = winding(500e3, 2e8, bundle_radius=100.)
    assert b['conductor_mass_kg']/a['conductor_mass_kg'] == pytest.approx(2.)
    assert b['ideal_tensile_mass_kg']/a['ideal_tensile_mass_kg'] == pytest.approx(4.)
    assert b['self_energy_J']/a['self_energy_J'] == pytest.approx(4.)
    assert a['cold_W'] > a['joints_cold_W'] > 0
    assert a['radiator']['heat_W'] == pytest.approx(a['cold_W']+a['refrigerator_W'])


def test_field_constraint_and_invalid_compact_coil_are_explicit():
    for current in [1e5, 1e8, 1e10]:
        w = winding(100e3, current)
        assert w['bundle_self_field_T'] <= 10.*(1+1e-12)
        assert w['bundle_radius_m'] >= w['conductor_pack_radius_m']
    assert not winding(100., 1e12)['admissible']


def test_upstream_zero_margin_and_distance_scaling():
    a = upstream_width(15e6, 7e6, lateral_speed=0., angle_deg=0.)
    assert a['required_obstacle_radius_m'] == 7e6
    b = upstream_width(30e6, 7e6)
    c = upstream_width(15e6, 7e6)
    assert b['required_obstacle_radius_m']-7e6 == pytest.approx(2*(c['required_obstacle_radius_m']-7e6))
