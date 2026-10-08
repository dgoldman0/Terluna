import numpy as np

from shared import constants as K
from .relative_orbit import hcw_elements, hcw_state
from .ring_keeping import argument_of_latitude, plan

A = 15.0e6
N = np.sqrt(K.MOON_GM/A**3)


def test_planned_burns_remove_the_element_error_in_the_linear_model():
    start = np.array([2e-5, 1e-3, 3e-5, -1e-5, 4e-5, -2e-5])
    error = np.array([-1.5e-5, 0., -2e-5, 2.5e-5, 1e-5, 3e-5])
    elements, u = start.copy(), 0.3
    for phase, (dvr, dvt, dvn) in sorted(plan(error, A, N), key=lambda b: (b[0]-u) % (2*np.pi)):
        dt = ((phase-u) % (2*np.pi))/N
        state = hcw_state(elements, A, N, u, dt)
        u = u+N*dt
        state[3:] += [dvr, dvt, dvn]
        elements = hcw_elements(state, A, N, u)
    keep = [0, 2, 3, 4, 5]
    assert np.allclose(elements[keep], (start+error)[keep], rtol=0, atol=1e-12)


def test_argument_of_latitude_is_measured_toward_the_motion():
    r = A*np.array([np.cos(.4), np.sin(.4), 0.])
    v = np.sqrt(K.MOON_GM/A)*np.array([-np.sin(.4), np.cos(.4), 0.])
    assert np.isclose(argument_of_latitude(np.r_[r, v], np.array([1., 0., 0.])), .4)
    assert np.isclose(argument_of_latitude(np.r_[r, -v], np.array([1., 0., 0.])), -.4)


def test_geometric_elements_match_the_linear_ones_for_small_offsets_and_hold_for_large_ones():
    from .relative_orbit import rtn_frame
    from .ring_keeping import elements, geometric_elements
    chief = np.r_[A, 0., 0., 0., np.sqrt(K.MOON_GM/A), 0.]
    small = np.array([1e-6, 0., 2e-6, -1e-6, 1.5e-6, -.5e-6])
    rel = hcw_state(small, A, N, 0., 0.)
    frame = rtn_frame(chief)
    omega = np.cross(chief[:3], chief[3:])/np.dot(chief[:3], chief[:3])
    dr = rel[:3]@frame
    deputy = np.r_[chief[:3]+dr, chief[3:]+rel[3:]@frame+np.cross(omega, dr)]
    origin = np.array([1., 0., 0.])
    exact = geometric_elements(chief, np.array([chief, deputy]), origin)[0][1]
    linear = elements(chief, np.array([chief, deputy]), origin)[0][1]
    keep = [0, 2, 3, 4, 5]
    assert np.allclose(exact[keep], linear[keep], rtol=0, atol=2e-9)
    # A circular orbit tilted by 0.3 degrees about the origin direction has no eccentricity difference.
    tilt = np.radians(.3)
    v = np.sqrt(K.MOON_GM/A)*np.array([0., np.cos(tilt), np.sin(tilt)])
    tilted = np.r_[A, 0., 0., v]
    g = geometric_elements(chief, np.array([tilted]), origin)[0][0]
    assert abs(g[0]) < 1e-12 and np.hypot(g[2], g[3]) < 1e-12
    assert np.isclose(np.hypot(g[4], g[5]), np.sin(tilt), rtol=1e-9)
