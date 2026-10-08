import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from shared import constants as K
from .collection import parallel_interception, require_uneclipsed_tiles, integrate_inventory
from .fleet import sun_points


def test_isolated_tile_cosine_power_and_specular_momentum():
    q = np.zeros((1, 3)); sources = np.array([[0., 0., K.AU]])
    for angle in [0., np.pi/3, 2*np.pi/3, np.pi]:
        frame = Rotation.from_rotvec([angle, 0., 0.]).as_matrix()
        out = parallel_interception(q, frame, sources, .2)
        expected = K.SOLAR_CONSTANT*9890.**2*abs(np.cos(angle))
        assert out['first_intercept_bolometric_equivalent_W'][0] == pytest.approx(expected)
        assert out['redirected_band_optical_W'][0] == pytest.approx(.2*expected)
        expected_force = -2*.2*expected*np.cos(angle)/K.SPEED_OF_LIGHT*frame[:, 2]
        np.testing.assert_allclose(out['reflected_force_N'][0], expected_force, atol=1e-12)
        np.testing.assert_allclose(out['front_first_intercept_W']+out['back_first_intercept_W'],
                                  out['first_intercept_bolometric_equivalent_W'])


def test_aligned_stack_counts_only_the_front_aperture():
    q = np.array([[0., 0., 0.], [0., 0., 4000.]])
    source = np.array([[0., 0., K.AU]])
    out = parallel_interception(q, np.eye(3), source, .138)
    expected = K.SOLAR_CONSTANT*9890.**2*(K.AU/(K.AU-4000.))**2
    assert out['first_intercept_bolometric_equivalent_W'][0] == 0.
    assert out['first_intercept_bolometric_equivalent_W'].sum() == pytest.approx(expected)
    assert out['unshadowed_aperture_W'].sum() > 1.99*expected


def test_finite_sun_feathered_power_is_small_and_illuminates_both_faces():
    q = np.zeros((1, 3)); frame = Rotation.from_rotvec([np.pi/2, 0., 0.]).as_matrix()
    sources = sun_points(np.array([0., 0., K.AU]), 512, .317)
    out = parallel_interception(q, frame, sources, .138)
    # Mean |y| on a uniform circular disk is 4R/(3*pi); small-angle limit.
    expected = K.SOLAR_CONSTANT*9890.**2*4*K.SUN_RADIUS/(3*np.pi*K.AU)
    assert out['first_intercept_bolometric_equivalent_W'][0] == pytest.approx(expected, rel=.01)
    assert out['front_first_intercept_W'][0] > 0
    assert out['back_first_intercept_W'][0] > 0


def test_body_eclipse_guard_rejects_occulted_prefixes():
    sun = np.array([[K.AU, 0., 0.]])
    earth = np.array([[0., 3e8, 0.]])
    require_uneclipsed_tiles(dict(sun=sun, earth=earth, moon=np.array([[-15e6, 0., 0.]])))
    with pytest.raises(ValueError, match='joint source/area'):
        require_uneclipsed_tiles(dict(sun=sun, earth=earth, moon=np.array([[15e6, 0., 0.]])))


def test_energy_inventory_uses_time_weights_and_preserves_member_totals():
    result = integrate_inventory([0., 1., 3.], [[2., 0.], [2., 2.], [2., 6.]])
    assert result['energy_J'] == 15.
    assert result['mean_W'] == 5.
    assert result['per_member_energy_J'] == [6., 9.]
    with pytest.raises(ValueError):
        integrate_inventory([0., 0.], [[1.], [1.]])
