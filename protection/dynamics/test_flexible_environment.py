"""Analytic third-body fields and forward-ray eclipse geometry, without files."""
import numpy as np
import pytest

from shared import constants as K
from protection.dynamics.flexible_environment import (
    differential_point_acceleration, lunar_environment, sphere_occlusion,
)


def test_differential_acceleration_matches_exact_axis_formula_without_cancellation():
    radius, gm = 1e12, 8e22
    for fraction in (1e-3, 1e-8, 1e-12):
        point = np.array([fraction*radius, 0., 0.])
        result = differential_point_acceleration(point, [dict(position_m=[radius, 0., 0.], gm_m3_s2=gm)])
        exact = gm/radius**2*fraction*(2-fraction)/(1-fraction)**2
        assert result == pytest.approx([exact, 0., 0.], rel=2e-14, abs=0.)


def test_general_field_matches_direct_point_gravity_and_superposes():
    points = np.array([[2., 3., -1.], [4., -.5, 6.]])
    sources = [dict(position_m=[100., 20., 50.], gm_m3_s2=400.),
               dict(position_m=[-40., -50., 70.], gm_m3_s2=200.)]
    expected = np.zeros_like(points)
    for source in sources:
        position, gm = np.asarray(source['position_m']), source['gm_m3_s2']
        offset = position-points
        expected += gm*(offset/np.linalg.norm(offset, axis=1)[:, None]**3-position/np.linalg.norm(position)**3)
    result = differential_point_acceleration(points, sources)
    assert result == pytest.approx(expected, rel=1e-12, abs=1e-15)
    separate = sum(differential_point_acceleration(points, [source]) for source in sources)
    assert result == pytest.approx(separate, rel=1e-14)


def test_origin_has_exactly_zero_acceleration_and_no_sources_has_zero_field():
    source = dict(position_m=[100., 20., 50.], gm_m3_s2=400.)
    assert np.array_equal(differential_point_acceleration(np.zeros((2, 3)), [source]), np.zeros((2, 3)))
    assert np.array_equal(differential_point_acceleration([[1., 2., 3.]], []), np.zeros((1, 3)))
    assert differential_point_acceleration(np.zeros((0, 3)), [source]).shape == (0, 3)


def test_differential_field_converges_to_independent_tidal_tensor():
    radial = np.array([1., 2., -3.])/np.sqrt(14.)
    radius, gm = 1e6, 1e12
    errors = []
    for length in (2000., 1000., 500.):
        point = length*np.array([1., -.4, .3])
        source = [dict(position_m=radius*radial, gm_m3_s2=gm)]
        odd_field = (differential_point_acceleration(point, source)-differential_point_acceleration(-point, source))/2
        linear = gm/radius**3*(3*radial*np.dot(radial, point)-point)
        errors.append(np.linalg.norm(odd_field-linear)/np.linalg.norm(linear))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_near_source_direct_branch_and_gravity_rotation_covariance():
    points = np.array([[8., 0., 0.], [3., -1., 2.]])
    source = dict(position_m=[10., 0., 0.], gm_m3_s2=100.)
    result = differential_point_acceleration(points, [source])
    assert result[0] == pytest.approx([24., 0., 0.])
    matrix, _ = np.linalg.qr(np.array([[1., 3., 2.], [2., 1., -3.], [3., -2., 1.]]))
    rotated = differential_point_acceleration(points@matrix.T,
                 [dict(position_m=matrix@np.asarray(source['position_m']), gm_m3_s2=100.)])
    assert rotated == pytest.approx(result@matrix.T, rel=1e-13, abs=1e-13)


def test_sphere_blocks_forward_rays_but_not_backward_rays():
    faces = np.array([[0., 0., 0.], [2., 0., 0.], [1., 0., 0.]])
    source = [[0., 0., 1.], [0., 0., -1.]]
    transmission = sphere_occlusion(faces, source, [dict(position_m=[0., 0., 5.], radius_m=1.)])
    assert np.array_equal(transmission, [[0., 1.], [1., 1.], [0., 1.]])


def test_inside_and_surface_limits_of_sphere_occlusion():
    source = [[0., 0., 1.], [0., 0., -1.]]
    sphere = [dict(position_m=[0., 0., 0.], radius_m=1.)]
    result = sphere_occlusion([[0., 0., 0.], [0., 0., 1.], [0., 0., -1.]], source, sphere)
    assert np.array_equal(result, [[0., 0.], [1., 0.], [0., 1.]])
    assert sphere_occlusion(np.empty((0, 3)), source, sphere).shape == (0, 2)
    assert np.array_equal(sphere_occlusion([[0., 0., 0.]], source, []), [[1., 1.]])


def test_external_grazing_and_two_occultors_are_resolved_separately():
    source = [[0., 0., 1.]]
    faces = [[1.+1e-8, 0., 0.], [3., 0., 0.]]
    one = dict(position_m=[0., 0., 5.], radius_m=1.)
    two = dict(position_m=[3., 0., 4.], radius_m=.5)
    assert np.array_equal(sphere_occlusion(faces, source, [one]), [[1.], [1.]])
    assert np.array_equal(sphere_occlusion(faces, source, [one, two]), [[1.], [0.]])


class SyntheticInterpolatedEphemeris:
    """Synthetic verification fixture; never substituted for a real ephemeris."""
    def __init__(self):
        self.calls = 0

    def at(self, t):
        self.calls += 1
        if not 0 <= t <= 100:
            raise ValueError('Outside synthetic test interval')
        return dict(positions=dict(sun=K.AU*np.array([np.cos(t/100), np.sin(t/100), 0.]),
                                   earth=np.array([0., 4e8, 0.])))


def test_ephemeris_wrapper_uses_supplied_time_state_and_shared_constants():
    ephemeris = SyntheticInterpolatedEphemeris()
    environment = lunar_environment(ephemeris, solar_disk_rings=1)
    first = environment(0.)
    assert len(first['directions']) == 7
    assert first['weights'].sum() == pytest.approx(1.)
    projected = first['pressure_Pa']*(first['weights']@(first['directions']@np.array([1., 0., 0.])))
    assert projected == pytest.approx(K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT)
    assert first['perturbing_bodies'][0]['gm_m3_s2'] == K.EARTH_GM
    assert first['perturbing_bodies'][1]['gm_m3_s2'] == K.SUN_GM
    assert first['occluders'][0]['radius_m'] == K.MOON_RADIUS
    assert first['occluders'][1]['radius_m'] == K.EARTH_RADIUS
    assert environment(0.) is first
    assert ephemeris.calls == 1
    changed = environment(20.)
    assert not np.allclose(changed['directions'], first['directions'])
    with pytest.raises(ValueError, match='Outside'):
        environment(101.)


def test_barycentric_state_wrapper_subtracts_moon_and_obeys_inverse_square_flux():
    class SyntheticStates:
        def state(self, name, t):
            moon = np.array([1e10+t, -2e10, 3e10])
            offset = dict(moon=np.zeros(3), sun=np.array([2*K.AU, 0., 0.]), earth=np.array([0., 4e8, 0.]))[name]
            return moon+offset, np.zeros(3)
    result = lunar_environment(SyntheticStates(), solar_disk_rings=0)(0.)
    assert result['pressure_Pa'] == pytest.approx(K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT/4)
    assert result['perturbing_bodies'][1]['position_m'] == pytest.approx([2*K.AU, 0., 0.])


@pytest.mark.parametrize('source', [dict(position_m=[0., 0., 0.], gm_m3_s2=1.),
                                  dict(position_m=[1., 0., 0.], gm_m3_s2=-1.),
                                  dict(position_m=[1., 0.], gm_m3_s2=1.),
                                  dict(position_m=[1., 0., 0.], gm_m3_s2=float('nan'))])
def test_invalid_perturbers_are_rejected(source):
    with pytest.raises(ValueError):
        differential_point_acceleration([[2., 0., 0.]], [source])


def test_point_gravity_singularity_and_nonunit_rays_are_rejected():
    with pytest.raises(ValueError, match='singular'):
        differential_point_acceleration([1., 0., 0.], [dict(position_m=[1., 0., 0.], gm_m3_s2=1.)])
    with pytest.raises(ValueError):
        sphere_occlusion([[0., 0., 0.]], [[0., 0., 2.]], [])
    with pytest.raises(ValueError):
        sphere_occlusion([[0., 0., 0.]], [[0., 0., 1.]], [dict(position_m=[0., 0., 1.], radius_m=0.)])
