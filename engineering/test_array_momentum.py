"""Conservation and frame-invariance tests for imposed-attitude momentum sizing."""
import numpy as np
import pytest

from engineering.array_momentum import angular_velocity_from_frames, integrate_rotor_momentum


def rz(angle):
    angle = np.asarray(angle)
    result = np.zeros(angle.shape + (3, 3))
    result[..., 0, 0] = result[..., 1, 1] = np.cos(angle)
    result[..., 0, 1] = -np.sin(angle)
    result[..., 1, 0] = np.sin(angle)
    result[..., 2, 2] = 1
    return result


def test_rotating_body_conserves_inertial_total_momentum():
    t = np.linspace(0, 20, 501)
    frames = rz(0.3 * t)
    omega = np.tile([0, 0, 0.3], (len(t), 1))
    inertia = np.diag([2., 3., 4.])
    row = integrate_rotor_momentum(t, frames, inertia, np.zeros((len(t), 3)),
                                   body_angular_velocity_rad_s=omega,
                                   initial_rotor_momentum_body_N_m_s=[8., 0., 0.])
    # With no external torque, a transverse rotor momentum cannot retain the
    # same body-axis components while its spacecraft turns.
    expected_rotor = np.c_[8*np.cos(0.3*t), -8*np.sin(0.3*t), np.zeros(len(t))]
    assert np.allclose(row['rotor_momentum_body_N_m_s'], expected_rotor, atol=1e-12)
    actual = np.einsum('tij,tj->ti', frames,
                       row['rotor_momentum_body_N_m_s'] + omega @ inertia)
    assert np.allclose(actual, np.tile([8, 0, 1.2], (len(t), 1)), atol=1e-12)
    assert row['peak_rotor_gyroscopic_torque_N_m'] == pytest.approx(2.4)
    assert row['max_derivative_balance_residual_N_m'] < 1.2e-4


def test_stationary_body_with_constant_external_torque():
    t = np.linspace(2, 12, 21)
    frames = np.broadcast_to(np.eye(3), (len(t), 3, 3))
    torque = np.tile([2., -3., 4.], (len(t), 1))
    row = integrate_rotor_momentum(t, frames, np.diag([2., 3., 4.]), torque)
    assert np.allclose(row['rotor_momentum_body_N_m_s'], (t-t[0])[:, None] * torque)
    assert np.allclose(row['rotor_reaction_torque_body_N_m'], -torque)
    assert row['max_derivative_balance_residual_N_m'] < 1e-12


def test_inertial_rotation_equivariance():
    t = np.linspace(0, 10, 101)
    frames = rz(0.1*t)
    rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    inertia = np.diag([2., 3., 4.])
    torque = np.c_[np.sin(t), np.cos(t), t/10]
    a = integrate_rotor_momentum(t, frames, inertia, torque)
    b = integrate_rotor_momentum(t, rotation @ frames, inertia, torque)
    assert np.allclose(a['rotor_momentum_body_N_m_s'], b['rotor_momentum_body_N_m_s'], atol=1e-12)
    assert np.allclose(b['total_momentum_inertial_N_m_s'],
                       a['total_momentum_inertial_N_m_s'] @ rotation.T, atol=1e-12)


def test_rotating_torque_quadrature_converges_at_second_order():
    errors = []
    for intervals in [40, 80, 160]:
        t = np.linspace(0, 1.7, intervals+1)
        frames = rz(t)
        row = integrate_rotor_momentum(t, frames, np.eye(3), np.tile([1., 0., 0.], (len(t), 1)),
                                       body_angular_velocity_rad_s=np.tile([0., 0., 1.], (len(t), 1)))
        exact_impulse = np.array([np.sin(t[-1]), 1-np.cos(t[-1]), 0.])
        errors.append(np.linalg.norm(row['external_angular_impulse_inertial_N_m_s'][-1] - exact_impulse))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_frame_angular_velocity_converges_and_recovers_axis():
    errors = []
    for intervals in [40, 80, 160]:
        t = np.linspace(0, 1.7, intervals+1)
        omega = angular_velocity_from_frames(t, rz(0.8*t))
        errors.append(np.max(np.abs(omega[:, 2] - 0.8)))
        assert np.allclose(omega[:, :2], 0)
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_changing_locked_inertia_requires_rotor_momentum():
    t = np.linspace(0, 1, 51)
    inertia = np.broadcast_to(np.eye(3), (len(t), 3, 3)).copy()
    inertia[:, 2, 2] = 1+t
    omega = np.tile([0., 0., 1.], (len(t), 1))
    row = integrate_rotor_momentum(t, rz(t), inertia, np.zeros((len(t), 3)),
                                   body_angular_velocity_rad_s=omega)
    assert np.allclose(row['rotor_momentum_body_N_m_s'][:, 2], -t)
    assert np.allclose(row['rotor_reaction_torque_body_N_m'][:, 2], 1)


def test_rejects_invalid_frames_and_inertia():
    t = np.arange(3.)
    frames = rz(t)
    with pytest.raises(ValueError, match='proper orthonormal'):
        integrate_rotor_momentum(t, frames*2, np.eye(3), np.zeros((3, 3)))
    with pytest.raises(ValueError, match='positive definite'):
        integrate_rotor_momentum(t, frames, np.diag([1., 0., 1.]), np.zeros((3, 3)))
