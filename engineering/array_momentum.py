"""Transport stored angular momentum through an imposed tile attitude history.

This is a rigid-body kinematic requirement calculation, not an attitude
controller or a flexible-body simulation. Frames map body vectors to inertial
vectors. The body inertia is the locked inertia of all hardware, including
rotors with their relative spins stopped. Relative rotor angular momentum is
then the residual between total inertial momentum and that locked-body motion.
External torque must include every applied environmental and unloading torque;
internal reaction-wheel torque cannot unload total angular momentum.

Unlike componentwise integration of body-axis torques, this calculation obeys
dH_inertial/dt = R @ torque_external_body. It does not reconstruct the old
attitude study's time histories from their maximum/median summary values.
"""
from __future__ import annotations

import numpy as np


def _sampled_inputs(times_s, frames_body_to_inertial):
    t = np.asarray(times_s, dtype=float)
    frames = np.asarray(frames_body_to_inertial, dtype=float)
    if t.ndim != 1 or len(t) < 3 or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError('times_s must contain at least three finite increasing samples')
    if frames.shape != (len(t), 3, 3) or not np.all(np.isfinite(frames)):
        raise ValueError('frames must have shape (samples, 3, 3) and be finite')
    gram = np.einsum('tji,tjk->tik', frames, frames)
    if not np.allclose(gram, np.eye(3), atol=1e-9, rtol=0) or not np.allclose(
            np.linalg.det(frames), 1, atol=1e-9, rtol=0):
        raise ValueError('frames must be proper orthonormal body-to-inertial rotations')
    return t, frames


def angular_velocity_from_frames(times_s, frames_body_to_inertial):
    """Body angular velocity from sampled proper rotations; second-order FD.

    The time grid must resolve rotations. This derivative loses accuracy for
    large inter-sample angles; callers with an analytic history should supply
    its analytic angular velocity to ``integrate_rotor_momentum`` instead.
    """
    t, frames = _sampled_inputs(times_s, frames_body_to_inertial)
    derivative = np.gradient(frames, t, axis=0, edge_order=2)
    spin = np.einsum('tji,tjk->tik', frames, derivative)
    spin = (spin - spin.transpose(0, 2, 1)) / 2
    return np.stack([spin[:, 2, 1], spin[:, 0, 2], spin[:, 1, 0]], axis=-1)


def integrate_rotor_momentum(
        times_s, frames_body_to_inertial, body_inertia_kg_m2,
        external_torque_body_N_m, *, body_angular_velocity_rad_s=None,
        initial_rotor_momentum_body_N_m_s=None):
    """Required net relative rotor momentum for the supplied body motion.

    ``body_inertia_kg_m2`` is a positive-definite 3x3 locked inertia, or a
    (samples, 3, 3) history. ``external_torque_body_N_m`` is (samples, 3).
    Initial relative rotor momentum defaults to zero, an explicit initial
    condition rather than an optimized initial bias. Numerical integration
    uses the composite trapezoid rule in inertial axes. Body derivatives use
    second-order finite differences. No wheel-axis allocation, torque limits,
    momentum limits, flexible response or initial spin-up is solved here.

    ``rotor_momentum_rate_body_N_m`` is dh/dt in moving body axes. The net
    reaction on the locked body is -(dh/dt + omega cross h); motor and bearing
    hardware must be allocated separately. A changing body inertia represents
    prescribed redistribution; its mechanism and energy are not designed.
    """
    t, frames = _sampled_inputs(times_s, frames_body_to_inertial)
    inertia = np.asarray(body_inertia_kg_m2, dtype=float)
    if inertia.shape == (3, 3):
        inertia = np.broadcast_to(inertia, (len(t), 3, 3))
    if inertia.shape != (len(t), 3, 3) or not np.all(np.isfinite(inertia)):
        raise ValueError('body inertia must be a finite 3x3 tensor or one tensor per sample')
    if not np.allclose(inertia, inertia.transpose(0, 2, 1), rtol=1e-12, atol=1e-12):
        raise ValueError('body inertia must be symmetric')
    if np.any(np.linalg.eigvalsh(inertia) <= 0):
        raise ValueError('body inertia must be positive definite')
    external = np.asarray(external_torque_body_N_m, dtype=float)
    if external.shape != (len(t), 3) or not np.all(np.isfinite(external)):
        raise ValueError('external torque must have shape (samples, 3) and be finite')
    omega = (angular_velocity_from_frames(t, frames) if body_angular_velocity_rad_s is None
             else np.asarray(body_angular_velocity_rad_s, dtype=float))
    if omega.shape != (len(t), 3) or not np.all(np.isfinite(omega)):
        raise ValueError('body angular velocity must have shape (samples, 3) and be finite')
    initial_rotor = (np.zeros(3) if initial_rotor_momentum_body_N_m_s is None
                     else np.asarray(initial_rotor_momentum_body_N_m_s, dtype=float))
    if initial_rotor.shape != (3,) or not np.all(np.isfinite(initial_rotor)):
        raise ValueError('initial rotor momentum must be a finite body-axis three-vector')

    locked = np.einsum('tij,tj->ti', inertia, omega)
    external_inertial = np.einsum('tij,tj->ti', frames, external)
    impulse = np.zeros_like(external_inertial)
    impulse[1:] = np.cumsum((external_inertial[:-1] + external_inertial[1:])
                            * np.diff(t)[:, None] / 2, axis=0)
    total_inertial = frames[0] @ (locked[0] + initial_rotor) + impulse
    total_body = np.einsum('tji,tj->ti', frames, total_inertial)
    rotor = total_body - locked
    rotor_rate = np.gradient(rotor, t, axis=0, edge_order=2)
    gyro = np.cross(omega, rotor)
    rotor_reaction = -rotor_rate - gyro
    # An independent derivative check converges with the time step; unlike the
    # algebraic momentum identity, it is not expected to vanish to roundoff.
    locked_required = np.gradient(locked, t, axis=0, edge_order=2) + np.cross(omega, locked)
    residual = locked_required - external - rotor_reaction
    return {
        'times_s': t,
        'body_angular_velocity_rad_s': omega,
        'locked_body_momentum_body_N_m_s': locked,
        'total_momentum_inertial_N_m_s': total_inertial,
        'external_angular_impulse_inertial_N_m_s': impulse,
        'rotor_momentum_body_N_m_s': rotor,
        'rotor_momentum_rate_body_N_m': rotor_rate,
        'rotor_gyroscopic_torque_body_N_m': gyro,
        'rotor_reaction_torque_body_N_m': rotor_reaction,
        'locked_body_required_torque_body_N_m': locked_required,
        'derivative_balance_residual_body_N_m': residual,
        'peak_rotor_momentum_N_m_s': float(np.max(np.linalg.norm(rotor, axis=1))),
        'peak_rotor_momentum_rate_N_m': float(np.max(np.linalg.norm(rotor_rate, axis=1))),
        'peak_rotor_gyroscopic_torque_N_m': float(np.max(np.linalg.norm(gyro, axis=1))),
        'peak_rotor_reaction_torque_N_m': float(np.max(np.linalg.norm(rotor_reaction, axis=1))),
        'max_derivative_balance_residual_N_m': float(np.max(np.linalg.norm(residual, axis=1))),
    }
