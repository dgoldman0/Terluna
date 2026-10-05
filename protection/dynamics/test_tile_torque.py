from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from shapely.geometry import box
from shapely.ops import unary_union

from shared import constants as K
from .tile_torque import rectangle_moments, illuminated_moments, gravity_torque
from .parallel_shadow import source_illumination
from .active_formation import lattice


def test_union_area_and_first_moments_match_independent_polygon_geometry():
    lo = np.array([[[-2., -1.], [-1., -2.], [0., -.5], [100., 100.]]])
    hi = np.array([[[1., 1.], [2., .5], [3., 2.], [200., 200.]]])
    active = np.array([[True, True, True, False]])
    shape = unary_union([box(*np.r_[a, b]) for a, b in zip(lo[0, :3], hi[0, :3])])
    actual = rectangle_moments(lo, hi, active)[0]
    np.testing.assert_allclose(actual, [shape.area, shape.area*shape.centroid.x,
                                      shape.area*shape.centroid.y], atol=1e-12)


def test_illuminated_moments_preserve_the_force_shadow_union():
    offsets = lattice(5, 8500., 4000.)[0]
    q = offsets[:, [1, 2, 0]]+[1e6, 2e6, 15e6]
    source = np.array([0., 0., K.AU])
    for angle in [0., .4, 1.565, 1.575, 2.2]:
        frame = Rotation.from_rotvec([angle, .12, .4]).as_matrix()
        m = illuminated_moments(q, frame, source)
        np.testing.assert_allclose(m[:, 0]/9890.**2,
                                  source_illumination(q, frame, source), atol=1e-12)
        assert np.all(abs(m[:, 1:]) <= (9890/2)*m[:, :1]+1e-2)


def test_uniform_unshadowed_aperture_has_no_centre_ray_photon_torque():
    result = illuminated_moments(np.zeros((1, 3)), np.eye(3), np.array([0., 0., K.AU]))
    np.testing.assert_allclose(result, [[9890.**2, 0., 0.]])


def test_finite_gravity_torque_matches_the_gravity_gradient_limit():
    env = SimpleNamespace(gravity=lambda t, q: -K.MOON_GM*q/np.linalg.norm(q, axis=-1)[..., None]**3)
    q = np.array([[15e6, 0., 0.]])
    frame = Rotation.from_rotvec([0., .7, .2]).as_matrix()
    direction = q[0]@frame/np.linalg.norm(q[0])
    inertia = 10000.**2/12*np.array([1., 1., 2.])
    expected = 3*K.MOON_GM/np.linalg.norm(q[0])**3*np.cross(direction, inertia*direction)
    a = gravity_torque(env, 0., q, frame)[0]
    b = gravity_torque(env, 0., q, frame, order=5)[0]
    # The finite square also has higher spatial moments, including a small
    # normal-axis torque absent from the leading gravity-gradient formula.
    assert np.linalg.norm(a-expected) < 2e-6*np.linalg.norm(expected)
    np.testing.assert_allclose(a, b, rtol=1e-8, atol=1e-11)
