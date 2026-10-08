"""Independent identities and limiting cases for the cycling experiment."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from shared import constants as K
from .cycling import (specular_sail, visible_sun, to_inertial, from_inertial,
                      dro_seed, cr3bp_derivative, jacobi_value, ray_geometry,
                      service_fraction, outgoing_clearance, service_coast_angles)
from .optical import optical_projection


def test_specular_force_obeys_reflected_momentum_and_projected_area():
    sun = np.array([K.AU, 0., 0.])
    angles = np.array([np.arctan(1/np.sqrt(2)), .37])
    a, area = specular_sail(sun, angles, .05, 1., 1.)
    n = a/np.linalg.norm(a)
    incoming = -sun/K.AU
    outgoing = incoming-2*np.dot(incoming, n)*n
    expected = K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*.05)*area*(incoming-outgoing)
    np.testing.assert_allclose(a, expected, rtol=1e-14)
    assert np.linalg.norm(outgoing) == pytest.approx(1)
    assert area == pytest.approx(np.sqrt(2/3))
    # The real tilted mirror lies inside the relaxed Sun-facing envelope.
    _, remainder = optical_projection(a, sun, .05)
    assert np.linalg.norm(remainder) < 1e-18
    assert a[0] < 0
    assert specular_sail(sun, angles, .05, .138, 0.)[0] == pytest.approx(np.zeros(3))
    with pytest.raises(ValueError):
        specular_sail(sun, np.array([1.6, 0]), .05, 1., 1.)


def test_moon_eclipse_and_double_occultation_are_not_double_counted():
    sun = np.array([K.AU, 0., 0.])
    behind = np.array([-K.EARTH_MOON_DISTANCE, 0., 0.])
    moon = np.array([50e6, 0., 0.])
    assert visible_sun(sun, behind, moon) == 0
    assert visible_sun(sun, -behind, moon) == 0
    assert visible_sun(sun, behind, -moon) == 1
    np.testing.assert_allclose(visible_sun(np.array([sun, sun]),
                                        np.array([behind, behind]), np.array([moon, -moon])), [0, 1])


def test_partial_three_disk_union_converges_and_exceeds_either_occultation():
    sun = np.array([K.AU, 0., 0.])
    earth = K.EARTH_RADIUS/np.sin(.003)*np.array([np.cos(.0025), np.sin(.0025), 0.])
    moon = K.MOON_RADIUS/np.sin(.002)*np.array([np.cos(.001), -np.sin(.001), 0.])
    v96 = visible_sun(sun, earth, moon, 96)
    v384 = visible_sun(sun, earth, moon, 384)
    assert v96 == pytest.approx(v384, abs=2e-5)
    assert 0 < v96 < visible_sun(sun, earth, -moon)
    assert v96 < visible_sun(sun, -earth, moon)


def test_moving_frame_transform_includes_rotation_velocity():
    omega = 2e-6
    frame = np.eye(3)
    derivative = np.array([[0, -omega, 0], [omega, 0, 0], [0, 0, 0]])
    local = np.array([50e6, 0, 0, 0, -470, 2.])
    inertial = to_inertial(local, frame, derivative)
    assert inertial[4] == pytest.approx(-370)
    np.testing.assert_allclose(from_inertial(inertial, frame, derivative), local, atol=1e-12)


def test_dro_seed_returns_and_jacobi_is_conserved_without_sail():
    local, duration = dro_seed(50e6)
    mu = K.MOON_GM/(K.EARTH_GM+K.MOON_GM)
    D = K.EARTH_MOON_DISTANCE
    T = np.sqrt(D**3/(K.EARTH_GM+K.MOON_GM))
    scale = np.array([D]*3+[D/T]*3)
    initial = local/scale
    initial[0] += 1-mu
    sol = solve_ivp(cr3bp_derivative, [0, duration/T], initial, args=(mu,),
                    rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(sol.y[:, -1], initial, atol=2e-9, rtol=0)
    states = sol.y.T.copy()
    states[:, 0] -= 1-mu
    values = jacobi_value(states*scale, D)
    assert np.ptp(values) < 1e-10
    # An independent centred derivative checks the sign of the feedback law.
    a = np.array([.00002, -.00004, .00001])
    step = .1
    plus, minus = local.copy(), local.copy()
    plus[3:] += a*step; minus[3:] -= a*step
    measured = (jacobi_value(plus, D)-jacobi_value(minus, D))/(2*step)
    assert measured == pytest.approx(-2*np.dot(local[3:], a)/(D/T)**2, rel=1e-6)


def static_sample():
    return dict(positions={"sun": np.array([K.AU, 0., 0.]),
                           "earth": np.array([-K.EARTH_MOON_DISTANCE, 0., 0.])},
                sun_v=np.zeros(3), earth_v=np.zeros(3), moon_v=np.zeros(3),
                sun_a=np.zeros(3), earth_a=np.zeros(3), moon_a=np.zeros(3))


def test_service_counts_projected_rays_and_rejects_nightside_and_bad_reflections():
    q = np.array([[50e6, 0, 0], [-50e6, 0, 0], [50e6, 0, 0]])
    rays = ray_geometry(q, static_sample())
    angles = np.array([[0., 0.], [0., 0.], [np.pi/2, 0.]])
    trace = dict(state=np.column_stack([q, np.zeros_like(q)]), angles=angles,
                 projected=np.cos(angles[:, 0]), visibility=np.ones(3), **rays)
    np.testing.assert_allclose(service_fraction(trace, 4*K.MOON_RADIUS), [1, 0, 0], atol=1e-15)
    assert rays["central"].tolist() == [True, False, True]
    assert outgoing_clearance(q, rays, angles, 4*K.MOON_RADIUS)[2] < 0


def test_feathering_keeps_the_entire_reflected_solar_cone_off_the_protected_sphere():
    rng = np.random.default_rng(20261004)
    q = rng.normal(size=(2000, 3))
    q *= 15e6/np.linalg.norm(q, axis=-1)[:, None]
    state = np.column_stack([q, np.zeros_like(q)])
    sample = static_sample()
    rays = ray_geometry(q, sample)
    angles = service_coast_angles(state, rays, sample, 4*K.MOON_RADIUS)
    active = np.cos(angles[:, 0]) > 1e-7
    clearance = outgoing_clearance(q, rays, angles, 4*K.MOON_RADIUS)
    assert np.count_nonzero(active) > 50
    assert np.min(clearance[active]) > 0
