"""Independent field identities, analytic benchmarks and reduced-model checks."""
from copy import deepcopy
from math import pi

import numpy as np
import pytest

from engineering.array_tile import Material, membrane_coefficient
from engineering.array_support import scalloped_support
from engineering.array_flexible import (
    fixed_square_membrane, model_from_support, scalloped_mesh, triangle_areas,
    _membrane_matrices,
)


def support():
    return scalloped_support(1000., 1., 50., Material(), Material(),
                             boom_depth_m=10., bay_m=5., chord_radius_m=.2,
                             hub_side_m=100.)


def model(divisions=6, modes=8, **kwargs):
    return model_from_support(support(), .01, mesh_divisions=divisions,
                              retained_modes=modes, **kwargs)


def test_triangle_patch_energy_reproduces_exact_constant_gradient():
    nodes, triangles, _, _ = scalloped_mesh(10., 0., 8)
    stiffness, mass, area = _membrane_matrices(nodes, triangles, 7., .1)
    displacement = .4+2*nodes[:, 0]-3*nodes[:, 1]
    energy = .5*displacement@(stiffness@displacement)
    assert energy == pytest.approx(.5*7.*100.*(2**2+3**2), rel=1e-13)
    assert mass.sum() == pytest.approx(.1*area.sum())


def test_fixed_square_frequency_and_sag_converge_to_independent_analytic_solutions():
    eigen_exact = 2*pi**2*2/.1/10**2
    sag_exact = membrane_coefficient(81)*.01*10**2/2
    eigen_errors, sag_errors = [], []
    for n in (8, 16, 24):
        result = fixed_square_membrane(side_m=10., tension_N_m=2., areal_mass_kg_m2=.1,
                                       mesh_divisions=n, pressure_Pa=.01)
        eigen_errors.append(abs(result['omega_squared_s2'][0]/eigen_exact-1))
        sag_errors.append(abs(result['centre_displacement_m']/sag_exact-1))
        # Square symmetry preserves the (1,2)/(2,1) frequency pair.
        assert result['omega_squared_s2'][1] == pytest.approx(result['omega_squared_s2'][2], rel=1e-12)
    assert eigen_errors[0] > eigen_errors[1] > eigen_errors[2]
    assert sag_errors[0] > sag_errors[1] > sag_errors[2]
    assert eigen_errors[-1] < .002
    assert sag_errors[-1] < .0005


def test_scalloped_mesh_has_positive_areas_and_converges_to_analytic_area():
    baseline = support()
    exact = baseline['geometry']['retained_area_m2']
    errors = []
    for n in (4, 8, 16):
        nodes, triangles, corners, edges = scalloped_mesh(1000., 50., n)
        areas = triangle_areas(nodes, triangles)
        assert (areas > 0).all()
        assert nodes[corners, :2] == pytest.approx(np.array([[-500., -500.], [500., -500.], [500., 500.], [-500., 500.]]))
        errors.append(abs(areas.sum()/exact-1))
        assert all(len(edge) == n+1 for edge in edges)
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_all_mass_is_finite_and_accounted_without_a_massless_hub_anchor():
    result = model(hub_equipment_mass_kg=2000.)
    expected = support()['mass_kg']+.01*support()['geometry']['retained_area_m2']+2000.
    assert result.nodal_masses_kg.sum() == pytest.approx(expected, rel=1e-13)
    assert (result.nodal_masses_kg > 0).all()
    assert result.nodal_masses_kg@result.nodes_m == pytest.approx([0., 0., 0.], abs=1e-8)
    assert result.metadata['hub_equipment_mass_kg'] == 2000.
    assert 'endpoint' in result.metadata['mass_rule']
    assert result.metadata['nodal_inertia_relative_frobenius_error'] > 0.
    inertia = np.array(result.metadata['nodal_reference_inertia_kg_m2'])
    assert inertia[2, 2] == pytest.approx(np.sum(result.nodal_masses_kg*np.sum(result.nodes_m[:, :2]**2, axis=1)))


def test_exactly_three_affine_normal_modes_have_zero_energy():
    result = model(divisions=4, modes=8)
    stiffness = result.stiffness_N_m.toarray()
    rigid = result.rigid_normal_modes
    relative = np.linalg.norm(stiffness@rigid)/(np.linalg.norm(stiffness)*np.linalg.norm(rigid))
    assert relative < 1e-14
    eigenvalues = np.linalg.eigvalsh(stiffness)
    tolerance = np.linalg.norm(stiffness)*1e-10
    assert np.max(np.abs(eigenvalues[:3])) < tolerance
    assert eigenvalues[3] > tolerance
    assert result.metadata['checks']['removed_rigid_normal_modes'] == 3


def test_modes_are_mass_orthogonal_to_rigid_motions_and_diagonalize_energy():
    result = model()
    shapes = result.scalar_modes
    assert shapes.T@(result.nodal_masses_kg[:, None]*shapes) == pytest.approx(np.eye(8), abs=1e-12)
    overlap = result.rigid_normal_modes.T@(result.nodal_masses_kg[:, None]*shapes)
    assert np.max(np.abs(overlap)) < 1e-9
    assert shapes.T@(result.stiffness_N_m@shapes) == pytest.approx(np.diag(result.omega_squared_s2), abs=1e-12)
    assert result.metadata['checks']['generalized_eigen_residual_relative'] < 1e-10


def test_internal_force_is_energy_gradient_and_has_no_net_force_or_torque():
    result = model()
    rng = np.random.default_rng(29018)
    w, direction = rng.normal(size=(2, len(result.nodes_m)))
    force = -(result.stiffness_N_m@w)
    epsilon = 1e-5
    energy = lambda displacement: .5*displacement@(result.stiffness_N_m@displacement)
    derivative = (energy(w+epsilon*direction)-energy(w-epsilon*direction))/(2*epsilon)
    assert derivative == pytest.approx(-force@direction, rel=1e-9)
    assert result.rigid_normal_modes.T@force == pytest.approx([0., 0., 0.], abs=1e-7)


def test_uniform_acceleration_of_all_mass_has_no_elastic_forcing():
    result = model()
    force = .02*result.nodal_masses_kg
    assert result.scalar_modes.T@force == pytest.approx(np.zeros(8), abs=1e-12)


def test_finite_boundary_moves_and_modal_static_compliance_converges():
    coarse = model(divisions=6, modes=4)
    medium = model(divisions=6, modes=12)
    full = model(divisions=6, modes=len(coarse.nodes_m)-3)
    force = np.zeros(len(coarse.nodes_m))
    force[coarse.corner_indices[0]] = 1.
    work = []
    for current in (coarse, medium, full):
        projection = current.scalar_modes.T@force
        coordinates = projection/current.omega_squared_s2
        displacement = current.scalar_modes@coordinates
        work.append(float(force@displacement))
    assert 0 < work[0] <= work[1] <= work[2]
    # Full static response after mass-balanced removal of rigid acceleration.
    rigid, mass = full.rigid_normal_modes, full.nodal_masses_kg
    rigid_acceleration = np.linalg.solve(rigid.T@(mass[:, None]*rigid), rigid.T@force)
    balanced_force = force-mass*(rigid@rigid_acceleration)
    assert full.stiffness_N_m@displacement == pytest.approx(balanced_force, abs=1e-10)
    # The deformed cable corner differs from the best-fit hub plane.
    fit = np.linalg.lstsq(rigid[full.hub_indices], displacement[full.hub_indices], rcond=None)[0]
    residual = displacement[full.corner_indices]-rigid[full.corner_indices]@fit
    assert np.max(np.abs(residual)) > 1e-5


def test_softer_support_and_heavier_hub_have_physical_modal_effects():
    base_support = support()
    softer = deepcopy(base_support)
    for name in ('radial_boom', 'hub_member'):
        softer[name]['beam_column_amplification_bound'] *= 4
    base = model_from_support(base_support, .01, mesh_divisions=6, retained_modes=4)
    soft = model_from_support(softer, .01, mesh_divisions=6, retained_modes=4)
    heavy = model_from_support(base_support, .01, mesh_divisions=6, retained_modes=4, hub_equipment_mass_kg=1e5)
    assert soft.omega_squared_s2[0] < base.omega_squared_s2[0]
    assert heavy.omega_squared_s2[0] < base.omega_squared_s2[0]


def test_low_modes_show_mesh_convergence_without_claiming_a_converged_structure():
    results = [model(divisions=n, modes=4) for n in (8, 12, 16)]
    first_difference = np.linalg.norm(results[1].omega_squared_s2-results[0].omega_squared_s2)
    second_difference = np.linalg.norm(results[2].omega_squared_s2-results[1].omega_squared_s2)
    assert second_difference < first_difference/2
    assert abs(results[2].omega_squared_s2[0]/results[1].omega_squared_s2[0]-1) < .005


def test_displaced_nodes_preserve_planar_coordinates_and_modal_shape():
    result = model()
    coordinates = np.arange(8)*.1
    displaced = result.displaced_nodes(coordinates)
    assert np.array_equal(displaced[:, :2], result.nodes_m[:, :2])
    assert displaced[:, 2] == pytest.approx(result.scalar_modes@coordinates)


@pytest.mark.parametrize('kwargs', [dict(mesh_divisions=1), dict(mesh_divisions=41),
                                   dict(retained_modes=0), dict(retained_modes=10000),
                                   dict(hub_equipment_mass_kg=-1)])
def test_invalid_model_requests_rejected(kwargs):
    with pytest.raises(ValueError):
        model_from_support(support(), .01, **kwargs)


def test_inverted_triangles_rejected():
    nodes = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    with pytest.raises(ValueError):
        triangle_areas(nodes, np.array([[0, 2, 1]]))
