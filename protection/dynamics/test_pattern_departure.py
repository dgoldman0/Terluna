from types import SimpleNamespace

import numpy as np
import pytest

from .test_cycling import static_sample
from .pattern_departure import sun_frame, departure_command, expansion_maps, expansion_state
from .departure_control import control_modes, mode_impulses, minimize_impulse


def inertial_environment():
    return SimpleNamespace(at=lambda t: static_sample(), gravity=lambda t, q: np.zeros_like(q),
        central_fraction=0., sigma=.05)


def test_expansion_response_has_the_analytic_symmetric_pulse_moment():
    env = inertial_environment()
    initial = np.array([[15e6, -100., 30., 0., .02, .01], [15e6, 100., -30., 0., -.02, -.01]])
    end, burn = 7200., 3600.
    maps = expansion_maps(env, initial, 0., end, burn)
    impulse = np.array([[.2, -.3, .1], [-.2, .3, -.1]])
    final = expansion_state(maps, initial, impulse, [end])[0]
    expected = initial.copy()
    expected[:, :3] += initial[:, 3:]*end+impulse*(end-burn/2)
    expected[:, 3:] += impulse
    np.testing.assert_allclose(final, expected, atol=2e-6, rtol=1e-12)
    np.testing.assert_allclose(expansion_state(maps, initial, impulse, [0.])[0], initial, atol=1e-9)


def test_edge_turn_has_the_quintic_rate_limit_and_quiet_boundaries():
    env = inertial_environment(); end, delay, duration = 7200., 1800., 3600.
    command = departure_command(env, 0., end, delay, duration, 1, -1)
    frame = sun_frame(env, 0.)
    np.testing.assert_allclose(command(0.).as_matrix(), frame, atol=1e-14)
    assert abs(command(end).as_matrix()[:, 2]@frame[:, 2]) < 1e-14
    times = np.arange(0., end+1, 2.)
    rate = np.rad2deg(np.linalg.norm(command(times, 1), axis=1))
    assert rate.max() == pytest.approx(90*1.875/duration, rel=2e-5)
    assert np.linalg.norm(command(delay, 1)) < 1e-7
    assert np.linalg.norm(command(delay+duration, 1)) < 1e-7
    assert np.linalg.norm(command(delay, 2)) < 1e-8


def test_preparation_modes_have_zero_net_impulse_and_exact_vector_cap():
    modes = control_modes(inertial_environment(), 0.)
    impulse = mode_impulses(modes, np.array([1., .7, 1.3]))
    np.testing.assert_allclose(impulse.sum(axis=0), np.zeros(3), atol=1e-12)
    # These three orthogonal modes give the squared norm by summing squares,
    # including every tile's distinct layer and row/column coefficient.
    squared = np.sum(modes*modes, axis=1)@np.array([1., .7, 1.3])**2
    np.testing.assert_allclose(squared, np.sum(impulse*impulse, axis=1), atol=1e-14)


def test_local_convex_program_matches_an_independent_minimum_norm_problem():
    modes = np.eye(3)[None, :, :]
    target = np.array([.25, .5, .75])
    cuts = [dict(gradient=(10000*np.eye(3)[i]).tolist(), rhs=10000*target[i]) for i in range(3)]
    x, fit = minimize_impulse(modes, np.ones(3), cuts, 3600.)
    assert fit['success']
    np.testing.assert_allclose(x, target, atol=1e-8)
    assert fit['mean_delta_v_m_s'] == pytest.approx(np.linalg.norm(target), abs=1e-8)
    _, blocked = minimize_impulse(modes, target, cuts, 3600., cap=.0001)
    assert not blocked['success']
