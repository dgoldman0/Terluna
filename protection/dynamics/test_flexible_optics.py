"""Analytic momentum, ray topology and deformed-geometry regression tests."""
from math import cos, pi, sin, sqrt

import numpy as np
import pytest

from protection.dynamics.flexible_optics import (
    MeshTile, photon_forces, receiver_rays, solar_disk_quadrature,
)


def triangle(z=0., a=1., r=0., t=0., reverse=False):
    faces = [[0, 2, 1]] if reverse else [[0, 1, 2]]
    return MeshTile([[0., 0., z], [2., 0., z], [0., 1., z]], faces, a, r, t)


def square(z=0., a=1., r=0., t=0., half=1., offset_x=0.):
    return MeshTile(np.array([[-half, -half, z], [half, -half, z],
                              [half, half, z], [-half, half, z]])+[offset_x, 0., 0.],
                    [[0, 1, 2], [0, 2, 3]], a, r, t)


def force(tiles, s=(0., 0., 1.), **kwargs):
    return photon_forces(tiles, [s], [1.], pressure_Pa=2., **kwargs)


@pytest.mark.parametrize('a,r,t,expected', [(1., 0., 0., -2.),
                                           (0., 1., 0., -4.),
                                           (0., 0., 1., 0.),
                                           (.2, .3, .5, -1.6)])
def test_flat_triangle_analytic_normal_pressure(a, r, t, expected):
    p = force([triangle(a=a, r=r, t=t)])
    np.testing.assert_allclose(p['total_force'], [0, 0, expected], atol=1e-14)


def test_oblique_absorption_keeps_incoming_momentum_direction():
    s = np.array([sin(pi/3), 0., cos(pi/3)])
    p = force([triangle()], s=s)
    np.testing.assert_allclose(p['total_force'], -2*cos(pi/3)*s, atol=1e-14)


def test_oblique_specular_force_has_two_cosines():
    p = force([triangle(a=0., r=1., t=0.)], s=[sqrt(3)/2, 0., .5])
    np.testing.assert_allclose(p['total_force'], [0., 0., -1.], atol=1e-14)


def test_back_illumination_and_winding_invariance():
    p = force([triangle(a=.2, r=.8)], s=[0., 0., -1.])
    q = force([triangle(a=.2, r=.8, reverse=True)], s=[0., 0., -1.])
    np.testing.assert_allclose(p['total_force'], [0., 0., 3.6])
    np.testing.assert_allclose(p['nodal_forces'][0], q['nodal_forces'][0])


def test_nodal_force_torque_and_affine_virtual_work_consistency():
    mesh = triangle(a=.3, r=.7)
    p = force([mesh], s=[.6, 0., .8], torque_origin=[.2, .1, -.3])
    np.testing.assert_allclose(p['force_balance_residual_N'], 0., atol=1e-14)
    np.testing.assert_allclose(p['torque_balance_residual_N_m'], 0., atol=1e-14)
    velocity = np.array([.3, -.2, .4])+np.cross([.2, -.1, .3], mesh.vertices)
    nodal_power = np.sum(p['nodal_forces'][0]*velocity)
    face_power = p['face_forces'][0][0]@velocity.mean(axis=0)
    assert nodal_power == pytest.approx(face_power, abs=1e-14)


def test_rigid_rotation_and_translation_transform_wrench():
    mesh = triangle(a=.4, r=.6)
    angle = .7
    rotation = np.array([[cos(angle), 0, sin(angle)], [0, 1, 0],
                         [-sin(angle), 0, cos(angle)]])
    shift = np.array([10., -4., 7.])
    moved = MeshTile(mesh.vertices@rotation.T+shift, mesh.triangles, .4, .6, 0.)
    p = force([mesh])
    q = force([moved], s=rotation@[0., 0., 1.])
    np.testing.assert_allclose(q['total_force'], rotation@p['total_force'], atol=1e-13)
    np.testing.assert_allclose(q['total_torque'], rotation@p['total_torque']+
                               np.cross(shift, q['total_force']), atol=1e-12)


def test_opaque_upstream_tile_shadows_lower_tile_and_order_is_geometric():
    lower, upper = square(), square(z=1.)
    p = force([lower, upper])
    q = force([upper, lower])
    np.testing.assert_allclose(p['tile_forces'][0], 0., atol=1e-14)
    np.testing.assert_allclose(p['tile_forces'][1], [0, 0, -8.])
    np.testing.assert_allclose(p['total_force'], q['total_force'])
    np.testing.assert_allclose(p['face_direct_source_transmission'][0], 0.)


def test_multilayer_absorption_conserves_direct_momentum():
    lower = square(a=.5, t=.5)
    upper = square(z=1., a=.5, t=.5)
    p = force([lower, upper])
    # A=4, incident pressure=2; 75% of direct light is absorbed in two layers.
    np.testing.assert_allclose(p['total_force'], [0., 0., -6.])
    trace = receiver_rays([lower, upper], [[.1, .2, -1.]], [[0., 0., 1.]])
    np.testing.assert_allclose(trace['energy_transmittance'], [.25])
    assert trace['crossing_count'][0] == 2
    assert trace['first_hit_tile'][0] == 0


def test_shared_triangle_edge_is_one_crossing_not_two():
    mesh = square(a=.4, t=.6)
    p = receiver_rays([mesh], [[0., 0., -1.]], [[0., 0., 1.]], triangle_chunk_size=1)
    assert p['crossing_count'][0] == 1
    assert p['energy_transmittance'][0] == pytest.approx(.6)


def test_seam_with_different_coefficients_uses_lower_transmission():
    mesh = square(a=[.2, .7], t=[.8, .3])
    p = receiver_rays([mesh], [[0., 0., -1.]], [[0., 0., 1.]])
    assert p['crossing_count'][0] == 1
    assert p['energy_transmittance'][0] == pytest.approx(.3)


def test_folded_tile_distinct_crossings_are_both_counted():
    low, high = square(a=.5, t=.5), square(z=1., a=.5, t=.5)
    mesh = MeshTile(np.concatenate([low.vertices, high.vertices]),
                    np.concatenate([low.triangles, high.triangles+4]), .5, 0., .5)
    p = receiver_rays([mesh], [[0., 0., -1.]], [[0., 0., 1.]])
    assert p['crossing_count'][0] == 2
    assert p['energy_transmittance'][0] == pytest.approx(.25)
    shadow = force([mesh])
    no_shadow = force([mesh], self_shadow=False)
    np.testing.assert_allclose(shadow['total_force'], [0., 0., -6.])
    np.testing.assert_allclose(no_shadow['total_force'], [0., 0., -8.])


def test_parallel_near_parallel_miss_and_backward_hits_are_rejected():
    mesh = square()
    origins = [[0., 0., -1.], [0., 0., -1.], [0., 0., 1.], [3., 0., -1.]]
    directions = [[1., 0., 0.], [1., 0., 1e-15], [0., 0., 1.], [0., 0., 1.]]
    p = receiver_rays([mesh], origins, directions)
    assert not p['geometric_covered'].any()
    assert np.isinf(p['first_hit_distance_m']).all()
    np.testing.assert_allclose(p['energy_transmittance'], 1.)


def test_transparent_intersection_is_geometric_but_not_filtered_coverage():
    p = receiver_rays([square(a=0., t=1.)], [[0., 0., -1.]], [[0., 0., 1.]])
    assert p['geometric_covered'][0]
    assert not p['filtered_covered'][0]
    assert p['energy_transmittance'][0] == 1


def test_actual_warp_changes_force_and_receiver_hit():
    flat = triangle(a=0., r=1.)
    warped = MeshTile([[0., 0., 0.], [2., 0., 2.], [0., 1., 0.]],
                      [[0, 1, 2]], 0., 1., 0.)
    p, q = force([flat]), force([warped])
    assert abs(q['total_force'][0]) > 0
    assert not np.allclose(q['total_force'], p['total_force'])
    ray_o, ray_d = [[.5, .2, .25]], [[0., 0., 1.]]
    assert not receiver_rays([flat], ray_o, ray_d)['geometric_covered'][0]
    assert receiver_rays([warped], ray_o, ray_d)['geometric_covered'][0]


def test_updated_positions_are_not_cached_between_calls():
    mesh = square()
    before = receiver_rays([mesh], [[0., 0., -1.]], [[0., 0., 1.]])
    mesh.vertices[:, 0] += 10.
    after = receiver_rays([mesh], [[0., 0., -1.]], [[0., 0., 1.]])
    assert before['geometric_covered'][0] and not after['geometric_covered'][0]


def test_reflectivity_derivative_matches_local_force_difference():
    a, r, eps = .2, .3, 1e-6
    p = force([triangle(a=a, r=r, t=1-a-r)], s=[.6, 0., .8])
    q = force([triangle(a=a, r=r+eps, t=1-a-r-eps)], s=[.6, 0., .8])
    np.testing.assert_allclose((q['face_forces'][0]-p['face_forces'][0])/eps,
                               p['reflection_force_derivative'][0], rtol=1e-9, atol=1e-9)


def test_reflectivity_change_recomputes_downstream_illumination():
    lower = square()
    p = force([lower, square(z=1., a=.2, r=.2, t=.6)])
    q = force([lower, square(z=1., a=.2, r=.3, t=.5)])
    np.testing.assert_allclose(q['tile_forces'][0]-p['tile_forces'][0], [0, 0, .8])


def test_chunks_and_aabb_culling_preserve_intersections():
    meshes = [square(z=1., a=.5, t=.5), square(z=2., offset_x=10.),
              triangle(z=3., a=.25, t=.75)]
    rng = np.random.default_rng(4)
    origins = rng.uniform(-2., 2., size=(39, 3))
    origins[:, 2] = -1.
    p = receiver_rays(meshes, origins, [0., 0., 1.], chunk_size=1, triangle_chunk_size=1)
    q = receiver_rays(meshes, origins, [0., 0., 1.], chunk_size=128, triangle_chunk_size=256)
    for key in ('first_hit_distance_m', 'crossing_count', 'energy_transmittance'):
        np.testing.assert_allclose(p[key], q[key])


@pytest.mark.parametrize('rings,count', [(0, 1), (1, 7), (2, 19), (5, 91)])
def test_finite_source_counts_normalization_and_symmetry(rings, count):
    directions, weights = solar_disk_quadrature([0., 0., 1.], .2, rings=rings)
    assert directions.shape == (count, 3)
    np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1., atol=1e-14)
    assert weights.sum() == pytest.approx(1.)
    assert (weights > 0).all()
    np.testing.assert_allclose(weights@directions[:, :2], 0., atol=1e-15)


def test_finite_source_force_converges_to_uniform_cap_analytic_integral():
    alpha = .4
    errors = []
    for rings in (0, 1, 2, 5, 10):
        directions, weights = solar_disk_quadrature([0., 0., 1.], alpha, rings)
        # Force_z=-P*A*mean(cos(theta)^2) for absorption. Uniform cap has
        # cos(theta) uniform in [cos(alpha),1].
        p = photon_forces([triangle()], directions, weights, 2.)
        exact = -2*(1+cos(alpha)+cos(alpha)**2)/3
        errors.append(abs(p['total_force'][2]-exact))
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 2e-5


def test_finite_source_can_partly_illuminate_behind_edge():
    # A receiver at the lower triangle centroid sees an upper screen edge
    # through the extended source; a central source alone lies exactly on it.
    upper = square(z=10., half=1., offset_x=1.)
    directions, weights = solar_disk_quadrature([0., 0., 1.], .1, 5)
    traces = receiver_rays([upper], np.zeros((len(directions), 3)), directions)
    fraction = weights@traces['energy_transmittance']
    assert .3 < fraction < .7


def test_empty_scene_is_transparent_and_force_free():
    p = receiver_rays([], [[0., 0., 0.]], [0., 0., 1.])
    assert p['energy_transmittance'][0] == 1
    q = force([])
    np.testing.assert_allclose(q['total_force'], 0.)


@pytest.mark.parametrize('a,r,t', [(.2, .2, .2), (-.1, .1, 1.), (0., float('nan'), 1.)])
def test_invalid_optical_energy_partition_rejected(a, r, t):
    with pytest.raises(ValueError):
        triangle(a=a, r=r, t=t)


def test_invalid_rays_and_degenerate_geometry_rejected():
    with pytest.raises(ValueError):
        receiver_rays([square()], [[0., 0., 0.]], [0., 0., 2.])
    with pytest.raises(ValueError):
        MeshTile([[0., 0., 0.]], [[0., 0., 0.]], 1., 0., 0.)
    degenerate = MeshTile([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]], [[0, 1, 2]], 1., 0., 0.)
    with pytest.raises(ValueError):
        force([degenerate])


def test_external_attenuation_multiplies_mesh_layers_and_trim_derivative():
    lower, upper = square(a=.2, r=.3, t=.5), square(z=1., a=.2, r=.3, t=.5)
    baseline = force([lower, upper])
    result = force([lower, upper], external_transmission=lambda centres, s: np.full(len(centres), .25))
    np.testing.assert_allclose(result['total_force'], .25*baseline['total_force'])
    np.testing.assert_allclose(result['total_torque'], .25*baseline['total_torque'])
    for i in range(2):
        np.testing.assert_allclose(result['reflection_force_derivative'][i],
                                   .25*baseline['reflection_force_derivative'][i])
        np.testing.assert_allclose(result['face_external_source_transmission'][i], .25)
        np.testing.assert_allclose(result['face_direct_source_transmission'][i],
                                   .25*baseline['face_direct_source_transmission'][i])
        np.testing.assert_array_equal(result['face_unoccluded_source_weights'][i], 0.)
    assert result['external_transmission_applied']
    assert not baseline['external_transmission_applied']


def test_external_eclipse_uses_actual_centres_and_each_source_direction():
    mesh = square()
    directions = np.array([[0., 0., 1.], [.6, 0., .8]])
    calls = []
    def eclipse(centres, direction):
        calls.append((centres.copy(), direction.copy()))
        return ((centres[:, 0] > 0) & (direction[0] == 0)).astype(float)
    result = photon_forces([mesh], directions, [.25, .75], 2., external_transmission=eclipse)
    centres = mesh.vertices[mesh.triangles].mean(axis=1)
    assert len(calls) == 2
    for i, (points, direction) in enumerate(calls):
        np.testing.assert_allclose(points, centres)
        np.testing.assert_array_equal(direction, directions[i])
    np.testing.assert_allclose(result['face_external_source_transmission'][0], [.25, 0.])
    np.testing.assert_allclose(result['face_forces'][0], [[0., 0., -1.], [0., 0., 0.]])
    np.testing.assert_allclose(result['face_unoccluded_source_weights'][0], [.25, 0.])
    moved = MeshTile(mesh.vertices+[3., 0., 0.], mesh.triangles, 1., 0., 0.)
    after = photon_forces([moved], directions, [.25, .75], 2., external_transmission=eclipse)
    np.testing.assert_allclose(after['face_external_source_transmission'][0], [.25, .25])
    np.testing.assert_allclose(after['total_force'], 2*result['total_force'])


@pytest.mark.parametrize('callback', [0.5, lambda p, s: .5,
    lambda p, s: np.full(len(p), -1e-6), lambda p, s: np.full(len(p), 1.01),
    lambda p, s: np.full(len(p), np.nan), lambda p, s: np.ones((len(p), 1))])
def test_invalid_external_attenuation_rejected(callback):
    with pytest.raises(ValueError, match='[Ee]xternal'):
        force([triangle()], external_transmission=callback)
