import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from shared import constants as K
from .relative_orbit import (drift_per_revolution, energy_matching_burns, hcw_elements, hcw_state,
                             minimum_rn_separation, normal_thickness, period, rtn_frame, rtn_relative,
                             semi_major_axis)

A = 15.0e6
N = np.sqrt(K.MOON_GM/A**3)


def two_body(state, duration):
    def rhs(t, y):
        r = y[:3]
        return np.r_[y[3:], -K.MOON_GM*r/np.linalg.norm(r)**3]
    return solve_ivp(rhs, (0., duration), state, method='DOP853', rtol=1e-12, atol=1e-6).y[:, -1]


def circular(radius=A):
    return np.r_[radius, 0., 0., 0., np.sqrt(K.MOON_GM/radius), 0.]


def test_elements_round_trip_through_the_linear_state():
    rng = np.random.default_rng(1)
    elements = rng.normal(scale=1e-3, size=(20, 6))
    u0 = rng.uniform(-np.pi, np.pi, size=20)
    state = hcw_state(elements, A, N, u0, 0.)
    assert np.allclose(hcw_elements(state, A, N, u0), elements, rtol=0, atol=1e-15)


def test_linear_solution_satisfies_hill_equations():
    elements = np.array([2e-4, -1e-3, 3e-4, -2e-4, 1e-4, 5e-4])
    t = np.linspace(0., 2e5, 4001)
    s = hcw_state(elements, A, N, .7, t)
    h = t[1]-t[0]
    acc = (s[2:, :3]-2*s[1:-1, :3]+s[:-2, :3])/h**2
    x, z = s[1:-1, 0], s[1:-1, 2]
    xd, yd = s[1:-1, 3], s[1:-1, 4]
    scale = A*N**2*1e-3
    assert np.abs(acc[:, 0]-2*N*yd-3*N**2*x).max() < 1e-5*scale
    assert np.abs(acc[:, 1]+2*N*xd).max() < 1e-5*scale
    assert np.abs(acc[:, 2]+N**2*z).max() < 1e-5*scale
    assert np.allclose(np.gradient(s[:, 0], h)[5:-5], s[5:-5, 3], rtol=0, atol=1e-6*A*N*1e-3)


def test_minimum_separation_matches_brute_force_and_the_published_form():
    rng = np.random.default_rng(7)
    de = rng.normal(scale=1e-4, size=(200, 2))
    di = rng.normal(scale=1e-4, size=(200, 2))
    def distance(u, k):
        x = -A*(de[k, 0]*np.cos(u)+de[k, 1]*np.sin(u))
        z = A*(di[k, 0]*np.sin(u)-di[k, 1]*np.cos(u))
        return np.hypot(x, z)
    grid = np.linspace(0., 2*np.pi, 4001)
    sampled = []
    for k in range(len(de)):
        values = distance(grid, k)
        j = int(np.argmin(values))
        local = minimize_scalar(distance, bounds=(grid[max(j-1, 0)], grid[min(j+1, len(grid)-1)]),
                                args=(k,), method='bounded', options=dict(xatol=1e-12))
        sampled.append(min(values[j], local.fun))
    sampled = np.array(sampled)
    closed = minimum_rn_separation(de, di, A)
    assert np.all(closed <= sampled+1e-6)
    assert np.allclose(closed, sampled, rtol=0, atol=1e-3)
    dot = np.abs(np.sum(de*di, axis=1))
    published = np.sqrt(2)*A*dot/np.sqrt(np.sum(de**2, axis=1)+np.sum(di**2, axis=1)+
                                         np.linalg.norm(de+di, axis=1)*np.linalg.norm(de-di, axis=1))
    assert np.allclose(closed, published, rtol=1e-9, atol=1e-9)
    parallel = minimum_rn_separation([[2e-4, 0.]], [[3e-4, 0.]], A)
    assert np.isclose(parallel[0], 2e-4*A)
    assert minimum_rn_separation([[2e-4, 0.]], [[0., 3e-4]], A)[0] < 1e-6


def test_energy_matching_sets_the_semi_major_axis_along_the_velocity():
    rng = np.random.default_rng(3)
    states = np.tile(circular(), (50, 1))
    states[:, :3] += rng.normal(scale=20e3, size=(50, 3))
    states[:, 3:] += rng.normal(scale=2., size=(50, 3))
    burns = energy_matching_burns(states, A)
    after = states.copy()
    after[:, 3:] += burns
    assert np.abs(semi_major_axis(after)-A).max() < 1e-6
    cosine = np.sum(burns*states[:, 3:], axis=1)/np.linalg.norm(burns, axis=1)/np.linalg.norm(states[:, 3:], axis=1)
    assert np.all(np.abs(np.abs(cosine)-1) < 1e-12)


def test_drift_per_revolution_matches_two_body_propagation():
    chief = circular()
    deputy = chief.copy()
    deputy[3:] += energy_matching_burns(deputy[None], A+1000.)[0]
    T = period(A)
    c, d = two_body(chief, T), two_body(deputy, T)
    along = rtn_relative(c, d[None])[0, 1]
    expected = drift_per_revolution(semi_major_axis(deputy[None])[0], A)
    assert abs(expected+3*np.pi*1000.) < 0.01*3*np.pi*1000.
    assert abs(along-expected) < 0.01*abs(expected)


def test_linear_model_tracks_two_body_relative_motion_through_an_orbit():
    chief = circular()
    elements = np.array([0., 2e-4, 1e-4, -0.5e-4, 0.8e-4, 0.6e-4])
    relative = hcw_state(elements, A, N, 0., 0.)
    frame = rtn_frame(chief)
    omega = np.cross(chief[:3], chief[3:])/np.dot(chief[:3], chief[:3])
    dr = relative[:3]@frame
    deputy = np.r_[chief[:3]+dr, chief[3:]+relative[3:]@frame+np.cross(omega, dr)]
    for fraction in (.25, .5, 1.):
        t = fraction*period(A)
        c, d = two_body(chief, t), two_body(deputy, t)
        actual = rtn_relative(c, d[None])[0]
        predicted = hcw_state(elements, A, N, 0., t)
        scale = A*np.linalg.norm(elements)
        assert np.linalg.norm(actual[:3]-predicted[:3]) < 0.01*scale


def test_normal_thickness_measures_the_spread_along_the_orbit_normal():
    states = np.tile(circular(), (3, 1))
    states[:, 2] = [-500., 0., 700.]
    assert np.isclose(normal_thickness(states), 1200.)
