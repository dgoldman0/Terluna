"""Force closure, geometric integration and finite-stiffness support checks."""
from math import cos, pi, sin, sqrt
from itertools import product

import numpy as np
import pytest

from engineering.array_tile import Material
from engineering.array_support import (
    SupportInfeasible, cellular_support, scalloped_geometry, scalloped_support,
    normal_acceleration_coefficients, square_normal_acceleration_envelope,
    box_normal_acceleration_envelope,
)


def support(**updates):
    inputs = dict(side_m=10000., tension_N_m=5., edge_retreat_m=500.,
                  compression_material=Material(), cable_material=Material())
    inputs.update(updates)
    return scalloped_support(**inputs)


def test_normal_acceleration_matches_direct_cross_products_and_gravity_tensor():
    rng = np.random.default_rng(18491)
    omega = rng.normal(size=(31, 3))*.002
    alpha = rng.normal(size=(31, 3))*1e-6
    radial = rng.normal(size=(31, 3))
    radial /= np.linalg.norm(radial, axis=1)[:, None]
    position = rng.normal(size=(31, 3))*3000.
    position[:, 2] = 0.
    n2 = 8e-9
    coefficients = normal_acceleration_coefficients(omega, alpha, radial, n2)
    rigid = np.cross(alpha, position)+np.cross(omega, np.cross(omega, position))
    differential_gravity = n2*(3*radial*np.sum(radial*position, axis=1)[:, None]-position)
    expected = (rigid-differential_gravity)[:, 2]
    assert np.allclose(np.sum(coefficients*position[:, :2], axis=1), expected,
                       rtol=1e-12, atol=1e-15)


def test_gravity_gradient_converges_to_exact_spherical_gravity():
    radius = 2e7
    n2 = 1e-8
    gm = n2*radius**3
    radial = np.array([1., 2., 3.])/sqrt(14.)
    coefficient = normal_acceleration_coefficients([0., 0., 0.], [0., 0., 0.], radial, n2)
    errors = []
    for length in [20000., 10000., 5000.]:
        offset = np.array([length, -.4*length, 0.])
        plus, minus = radius*radial+offset, radius*radial-offset
        gravity_plus = -gm*plus/np.linalg.norm(plus)**3
        gravity_minus = -gm*minus/np.linalg.norm(minus)**3
        exact = -(gravity_plus[2]-gravity_minus[2])/2
        linear = np.dot(coefficient, offset[:2])
        errors.append(abs((exact-linear)/linear))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_pure_normal_spin_has_zero_normal_centrifugal_acceleration():
    coefficient = normal_acceleration_coefficients([0., 0., .02], [0., 0., 0.], [1., 0., 0.], 1e-8)
    assert np.array_equal(coefficient, [0., 0.])
    assert square_normal_acceleration_envelope(10000., coefficient) == 0.


def test_square_envelope_equals_explicit_corner_maximum_and_broadcasts():
    omega = np.array([[.001, .002, .003], [-.002, .001, .005]])
    alpha = np.array([1e-6, -2e-6, 0.])
    coefficients = normal_acceleration_coefficients(omega, alpha, [0., 1., 0.], 0.)
    corners = np.array([[-500., -500.], [-500., 500.], [500., -500.], [500., 500.]])
    explicit = np.max(np.abs(coefficients@corners.T))
    assert square_normal_acceleration_envelope(1000., coefficients) == pytest.approx(explicit)


def test_sampled_normal_envelope_has_expected_temporal_convergence():
    errors = []
    for sample_count in [10, 20, 40]:
        times = np.linspace(0., pi, sample_count)
        alpha = np.zeros((sample_count, 3))
        alpha[:, 1] = np.sin(times)
        coefficients = normal_acceleration_coefficients([0., 0., 0.], alpha, [1., 0., 0.], 0.)
        errors.append(1.-square_normal_acceleration_envelope(2., coefficients))
    assert errors[0] > errors[1] > errors[2] > 0.
    assert errors[0]/errors[1] > 4.
    assert errors[1]/errors[2] > 4.


@pytest.mark.parametrize('omega,alpha,radial,n2', [
    ([0., 0.], [0., 0., 0.], [1., 0., 0.], 0.),
    ([0., 0., float('nan')], [0., 0., 0.], [1., 0., 0.], 0.),
    ([0., 0., 0.], [0., 0., 0.], [2., 0., 0.], 0.),
    ([0., 0., 0.], [0., 0., 0.], [1., 0., 0.], -1.),
    ([0., 0., 0.], [0., 0., 0.], [1., 0., 0.], float('inf')),
    ([0., 0., 0.], [0., 0., 0.], [1., 0., 0.], [0.]),
    (np.zeros((2, 3)), np.zeros((4, 3)), [1., 0., 0.], 0.),
])
def test_invalid_normal_acceleration_inputs_rejected(omega, alpha, radial, n2):
    with pytest.raises(ValueError):
        normal_acceleration_coefficients(omega, alpha, radial, n2)


@pytest.mark.parametrize('side,coefficients', [(-1., [0., 0.]), (10., [0., 0., 0.]),
                                            (10., []), (10., [0., float('nan')])])
def test_invalid_square_acceleration_envelopes_rejected(side, coefficients):
    with pytest.raises(ValueError):
        square_normal_acceleration_envelope(side, coefficients)


def test_box_envelope_matches_cross_product_field_at_all_eight_corners():
    spans = np.array([5000., 5100., 70.])
    omega = np.array([[.001, -.002, .003], [-.003, .002, .001]])
    alpha = np.array([2e-7, -1e-7, 3e-7])
    radial = np.array([1., -2., 3.])/sqrt(14.)
    n2 = 1e-8
    corners = np.array(list(product([-1., 1.], repeat=3)))*spans
    values = []
    for velocity in omega:
        for position in corners:
            rigid = np.cross(alpha, position)+np.cross(velocity, np.cross(velocity, position))
            gravity = n2*(3*radial*np.dot(radial, position)-position)
            values.append(abs((rigid-gravity)[2]))
    assert box_normal_acceleration_envelope(spans, omega, alpha, radial, n2) == pytest.approx(max(values), rel=1e-12)


def test_flat_box_reduces_to_square_and_finite_depth_adds_normal_acceleration():
    omega, alpha, radial, n2 = [.002, 0., 0.], [0., 0., 0.], [1., 0., 0.], 0.
    coefficients = normal_acceleration_coefficients(omega, alpha, radial, n2)
    flat = box_normal_acceleration_envelope([5000., 5000., 0.], omega, alpha, radial, n2)
    assert flat == square_normal_acceleration_envelope(10000., coefficients) == 0.
    assert box_normal_acceleration_envelope([5000., 5000., 50.], omega, alpha, radial, n2) == pytest.approx(50.*.002**2)


def test_box_gravity_envelope_converges_to_spherical_corner_accelerations():
    orbital_radius = 2e7
    n2 = 1e-8
    gm = n2*orbital_radius**3
    radial = np.array([1., 2., 3.])/sqrt(14.)
    errors = []
    for scale in [1., .5, .25]:
        spans = np.array([20000., 10000., 3000.])*scale
        corners = np.array(list(product([-1., 1.], repeat=3)))*spans
        exact = []
        for position in corners:
            plus, minus = orbital_radius*radial+position, orbital_radius*radial-position
            gp, gm_ = -gm*plus/np.linalg.norm(plus)**3, -gm*minus/np.linalg.norm(minus)**3
            exact.append(abs((gp[2]-gm_[2])/2))
        linear = box_normal_acceleration_envelope(spans, [0., 0., 0.], [0., 0., 0.], radial, n2)
        errors.append(abs((max(exact)-linear)/linear))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


@pytest.mark.parametrize('spans', [[1., 2.], [1., 2., -1.], [1., 2., float('inf')]])
def test_invalid_box_spans_rejected(spans):
    with pytest.raises(ValueError):
        box_normal_acceleration_envelope(spans, [0., 0., 0.], [0., 0., 0.], [1., 0., 0.], 0.)


def test_retained_area_matches_independent_quadrature_and_shallow_limit():
    result = scalloped_geometry(10000., 500.)
    radius = result['cable_radius_m']
    x = np.linspace(-5000., 5000., 100001)
    inward = np.sqrt(radius**2-x**2)-(radius-500.)
    integrated = np.trapezoid(inward, x)
    assert result['lost_area_m2'] == pytest.approx(4*integrated, rel=1e-8)
    assert result['guaranteed_centred_square_side_m'] == 9000.
    shallow = scalloped_geometry(10000., .0001)
    assert shallow['lost_area_m2'] == pytest.approx(4*2*10000.*.0001/3, rel=1e-10)


@pytest.mark.parametrize('retreat', [.0001, 50., 500., 2000.])
def test_area_moments_converge_and_respect_square_symmetry(retreat):
    coarse = scalloped_geometry(10000., retreat, quadrature_order=64)
    fine = scalloped_geometry(10000., retreat, quadrature_order=128)
    assert coarse['quadrature_area_relative_error'] < 1e-12
    assert np.allclose(coarse['membrane_area_inertia_per_mass_m2'],
                       fine['membrane_area_inertia_per_mass_m2'], rtol=1e-12, atol=0.)
    tensor = np.array(coarse['membrane_area_inertia_per_mass_m2'])
    assert tensor[0, 0] == tensor[1, 1] > 0
    assert tensor[2, 2] == 2*tensor[0, 0]
    if retreat == .0001:
        assert tensor[0, 0] == pytest.approx(10000.**2/12, rel=1e-7)


def test_area_moments_match_independent_strip_integration():
    geometry = scalloped_geometry(10000., 1500.)
    x = np.linspace(-5000., 5000., 200001)
    radius = geometry['cable_radius_m']
    h = 1500.-x*x/(radius+np.sqrt(radius*radius-x*x))
    qx = np.trapezoid(x*x*h, x)
    qy = np.trapezoid((5000.**3-(5000.-h)**3)/3, x)
    expected = 10000.**4/12-2*(qx+qy)
    assert geometry['retained_area_second_moment_x_m4'] == pytest.approx(expected, rel=1e-9)


def test_cable_load_closes_at_each_corner_and_finite_hub():
    result = support()
    theta = result['geometry']['cable_half_angle_rad']
    force = result['cable_tension_N']
    corner_component = force*(cos(theta)+sin(theta))
    assert result['radial_boom']['axial_force_N']/sqrt(2) == pytest.approx(corner_component)
    assert result['hub_member']['axial_force_N'] == pytest.approx(corner_component)
    assert result['cable_tension_N']/result['geometry']['cable_radius_m'] == 5.
    assert result['components_kg']['finite_hub'] > 0
    assert sum(result['components_kg'].values()) == pytest.approx(result['mass_kg'])
    assert result['diagonal_anchor_retreat_m'] <= 5.*(1+1e-12)
    for truss in [result['radial_boom'], result['hub_member']]:
        assert all(mass > 0 for mass in truss['components_kg'].values())
        assert truss['global_buckling_utilization'] <= .5*(1+1e-12)


def test_cable_inward_retreat_reduces_prestress_support_demand_but_loses_area():
    shallow = support(edge_retreat_m=50.)
    deep = support(edge_retreat_m=1000.)
    assert deep['mass_kg'] < shallow['mass_kg']
    assert deep['geometry']['retained_area_m2'] < shallow['geometry']['retained_area_m2']
    assert deep['geometry']['guaranteed_centred_square_side_m'] < shallow['geometry']['guaranteed_centred_square_side_m']


def test_cantilever_sensitivity_has_finite_mass_and_checks_tip_displacement():
    unloaded = support()
    loaded = support(tip_transverse_force_N=200., tip_moment_N_m=1e6,
                     distributed_transverse_load_N_m=.01)
    boom = loaded['radial_boom']
    length = boom['length_m']
    assert boom['root_transverse_force_N'] == pytest.approx(200.+.01*length)
    assert boom['root_bending_moment_N_m'] == pytest.approx(1e6+200.*length+.01*length**2/2)
    assert boom['lateral_displacement_m'] <= 5.*(1+1e-12)
    assert loaded['hub_member']['lateral_displacement_m'] <= 5.*(1+1e-12)
    assert loaded['mass_kg'] > unloaded['mass_kg']
    assert loaded['hub_member']['mass_kg'] >= unloaded['hub_member']['mass_kg']
    assert loaded['hub_member']['root_bending_moment_N_m'] > unloaded['hub_member']['root_bending_moment_N_m']


def test_full_membrane_normal_resultant_is_charged_once_at_boom_tips():
    pressure = 1e-5
    result = support(tip_transverse_force_N=17., membrane_normal_pressure_Pa=pressure)
    total = pressure*result['geometry']['retained_area_m2']
    assert result['membrane_normal_resultant_N'] == total
    assert result['membrane_normal_force_per_boom_N'] == total/4
    assert result['radial_boom']['root_transverse_force_N'] == pytest.approx(17.+total/4)
    assert result['membrane_normal_force_transfer_residual_N'] == 0.
    assert 'unresolved' in result['normal_boundary_rule']


def test_structural_acceleration_sizes_own_mass_not_only_external_payload():
    unloaded = support()
    loaded = support(structural_acceleration_m_s2=5e-4)
    boom = loaded['radial_boom']
    assert boom['distributed_load_including_self_N_m'] == pytest.approx(
        5e-4*boom['mass_kg']/boom['length_m'], rel=1e-8)
    assert loaded['mass_kg'] >= unloaded['mass_kg']
    assert boom['self_load_iterations'] > 2
    with pytest.raises(SupportInfeasible):
        support(structural_acceleration_m_s2=100.)


def test_inertia_is_symmetric_positive_and_respects_square_symmetry():
    result = support()
    inertia = np.array(result['inertia_kg_m2'])
    assert np.allclose(inertia, inertia.T)
    assert inertia[0, 0] == pytest.approx(inertia[1, 1])
    assert np.linalg.eigvalsh(inertia).min() > 0
    assert np.max(np.abs(inertia-np.diag(np.diag(inertia)))) < 1.
    # Positive finite thickness implies Ixx+Iyy exceeds Izz by 2*integral(z² dm).
    assert inertia[0, 0]+inertia[1, 1] > inertia[2, 2]
    assert result['half_depth_m'] >= 52.
    assert result['outside_span_m'] > 10000.


def test_shallow_cable_inertia_approaches_square_perimeter_limit():
    result = support(edge_retreat_m=.0001)
    # For a square perimeter, Izz/M=L²/3. Truss/cable masses are independently
    # finite; their combined tensor must remain positive even as R grows.
    inertia = np.array(result['inertia_kg_m2'])
    assert np.isfinite(inertia).all()
    assert np.linalg.eigvalsh(inertia).min() > 0


def test_obscuration_is_a_bound_separate_from_geometric_coverage():
    result = support()
    retained = result['geometry']['retained_area_m2']
    assert 0 < result['projected_support_obscuration_upper_bound_m2'] < retained
    assert result['unobscured_membrane_area_lower_bound_m2'] + result['projected_support_obscuration_upper_bound_m2'] == retained
    assert 'overcounted' in result['obscuration_rule']


def test_finite_tube_radius_closes_false_shallow_frame_clearance_margin():
    result = support(boom_depth_m=50.)
    coverage_pitch = 8500.*9000./9890.
    centre_plane_spacing = coverage_pitch*sin(1.5*pi/180.)
    naive_centreline_gap = centre_plane_spacing-50.
    material_gap = centre_plane_spacing-2*result['half_depth_m']
    assert naive_centreline_gap > 150.
    assert material_gap < 150.
    assert result['maximum_member_outer_radius_m'] > 2.
    assert result['enclosing_truss_width_m'] > 54.
    assert result['outside_span_m'] > 10000.+50./sqrt(2)+4.


def test_independent_cells_preserve_area_bookkeeping_and_parallel_axis_inertia():
    cells = cellular_support(10000., 2, 5., .05, Material(), Material(),
                             boom_depth_m=50., hub_side_m=500., bay_m=25., chord_radius_m=1.)
    cell = cells['cell']
    assert cells['cells'] == 4
    assert cells['mass_kg'] == 4*cell['mass_kg']
    assert cells['retained_area_m2'] == pytest.approx(scalloped_geometry(10000., 500.)['retained_area_m2'])
    expected_ixx = 4*cell['inertia_kg_m2'][0][0]+4*cell['mass_kg']*2500.**2
    assert cells['instantaneous_inventory_inertia_kg_m2'][0][0] == pytest.approx(expected_ixx)
    assert 'independent' in cells['evidence']


@pytest.mark.parametrize('changes', [dict(edge_retreat_m=2500.), dict(edge_retreat_m=2200.),
                                   dict(boom_depth_m=2000.), dict(buckling_factor=1.),
                                   dict(tip_transverse_force_N=-1), dict(hub_side_m=10000.),
                                   dict(structural_acceleration_m_s2=float('nan'))])
def test_invalid_geometry_or_loads_are_rejected(changes):
    with pytest.raises(ValueError):
        support(**changes)
