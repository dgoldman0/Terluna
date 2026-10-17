"""Independent uncondensed kinetic-energy checks of coupled accelerations.

The reference keeps three relative wheel-angle velocities, recovers the full
mass Hessian by evaluating kinetic energy, and differentiates generalized
momenta with respect to deformation. It does not reuse the implementation's
condensed inertia, cross-mass, centrifugal or Coriolis formulas.
"""
from itertools import product

import numpy as np
import pytest

from protection.dynamics.flexible_fleet import Body, make_state, mechanics


def energetic_body(rotors=True, damping=.03):
    # A small synthetic verification body, not an alternate physical sail.
    nodes = np.array([[x, y, 0.] for x, y in product((-1., 0., 1.), repeat=2)])
    masses = np.full(9, 2.)
    rigid = np.zeros((27, 6))
    for axis in range(3):
        rigid[axis::3, axis] = 1.
        rigid[:, axis+3] = np.cross(np.eye(3)[axis], nodes).ravel()
    rigid_basis, _ = np.linalg.qr(np.sqrt(np.repeat(masses, 3))[:, None]*rigid, mode='reduced')
    rng = np.random.default_rng(50329)
    shapes = rng.normal(size=(27, 3))
    shapes -= rigid_basis@(rigid_basis.T@shapes)
    shapes, _ = np.linalg.qr(shapes, mode='reduced')
    shapes /= np.sqrt(np.repeat(masses, 3))[:, None]
    options = dict(rotor_inertia_kg_m2=[.7, 1.1, .9],
                   torque_limit_N_m=[1., 1., 1.], wheel_capacity_N_m_s=[10., 10., 10.]) if rotors else {}
    return Body(nodes_m=nodes, triangles=np.array([[0, 3, 4], [0, 4, 1]]),
                nodal_masses_kg=masses, modes=shapes.reshape(9, 3, 3),
                omega_squared_s2=np.array([.3, .5, .9]), damping_ratio=damping,
                **options)


def uncondensed_kinetic_energy(body, q, generalized_speeds):
    k = body.k
    omega = generalized_speeds[:3]
    qdot = generalized_speeds[3:3+k]
    nodes = body.nodes_m+np.einsum('nik,k->ni', body.modes, q)
    shape_velocity = np.einsum('nik,k->ni', body.modes, qdot)
    velocity = np.cross(omega, nodes)+shape_velocity
    result = .5*np.sum(body.nodal_masses_kg[:, None]*velocity**2)
    if body.rotors:
        # Independently construct each thin rotor's transverse kinetic energy.
        relative_spin = generalized_speeds[3+k:]
        inertia = body.rotor_inertia_kg_m2
        for axis in range(3):
            transverse_square = np.sum(omega**2)-omega[axis]**2
            result += .25*inertia[axis]*transverse_square
            result += .5*inertia[axis]*(omega[axis]+relative_spin[axis])**2
    return float(result)


def kinetic_hessian(body, q):
    size = 3+body.k+(3 if body.rotors else 0)
    unit = np.eye(size)
    diagonal_energy = np.array([uncondensed_kinetic_energy(body, q, row) for row in unit])
    result = np.diag(2*diagonal_energy)
    for i in range(size):
        for j in range(i):
            value = uncondensed_kinetic_energy(body, q, unit[i]+unit[j])-diagonal_energy[i]-diagonal_energy[j]
            result[i, j] = result[j, i] = value
    return result


def kinetic_acceleration_reference(body, state, forces_body, motor_torque):
    omega, q, qdot = state['omega_rad_s'], state['q'], state['qdot']
    relative_spin = (state['wheel_momentum_N_m_s']/body.rotor_inertia_kg_m2-omega
                     if body.rotors else np.zeros(0))
    speed = np.r_[omega, qdot, relative_spin]
    matrix = kinetic_hessian(body, q)
    epsilon = 2e-5
    momentum_rate_at_fixed_speeds = (
        (kinetic_hessian(body, q+epsilon*qdot)-kinetic_hessian(body, q-epsilon*qdot))@speed/(2*epsilon))
    kinetic_gradient_q = np.empty(body.k)
    for j in range(body.k):
        step = np.eye(body.k)[j]*epsilon
        kinetic_gradient_q[j] = (uncondensed_kinetic_energy(body, q+step, speed)-
                                  uncondensed_kinetic_energy(body, q-step, speed))/(2*epsilon)
    momentum = matrix@speed
    bias = momentum_rate_at_fixed_speeds
    bias[:3] += np.cross(omega, momentum[:3])
    bias[3:3+body.k] -= kinetic_gradient_q
    nodes = body.nodes_m+np.einsum('nik,k->ni', body.modes, q)
    torque_external = np.sum(np.cross(nodes, forces_body), axis=0)
    force_modal = np.einsum('nik,ni->k', body.modes, forces_body)
    force_modal -= body.omega_squared_s2*q+body.damping*qdot
    # Motor action is internal to the body + wheels. In these uncondensed
    # coordinates it enters relative wheel angles only, with force -u.
    force = np.r_[torque_external, force_modal, -motor_torque] if body.rotors else np.r_[torque_external, force_modal]
    acceleration = np.linalg.solve(matrix, force-bias)
    return acceleration[:3], acceleration[3:3+body.k], acceleration[3+body.k:]


@pytest.mark.parametrize('rotors', [False, True])
def test_condensed_equations_match_independent_uncondensed_kinetic_energy(rotors):
    body = energetic_body(rotors=rotors)
    momentum = np.array([.3, -.25, .6]) if rotors else np.zeros(3)
    state = body.unpack(make_state(body, omega_rad_s=[.2, -.15, .3],
                                   q=[.04, -.03, .02], qdot=[-.02, .05, .03],
                                   wheel_momentum_N_m_s=momentum))
    forces = np.random.default_rng(50142).normal(size=(9, 3))*.07
    torque = np.array([.1, -.07, .03]) if rotors else np.zeros(3)
    expected_alpha, expected_qdd, relative_spin_acceleration = kinetic_acceleration_reference(body, state, forces, torque)
    alpha, qdd = mechanics(body, state, forces, torque)
    assert alpha == pytest.approx(expected_alpha, rel=2e-8, abs=2e-10)
    assert qdd == pytest.approx(expected_qdd, rel=2e-8, abs=2e-10)
    if rotors:
        assert body.rotor_inertia_kg_m2*(expected_alpha+relative_spin_acceleration) == pytest.approx(-torque, rel=2e-8, abs=2e-10)


def test_coupled_instantaneous_energy_balance_includes_motor_relative_power():
    body = energetic_body()
    state = body.unpack(make_state(body, omega_rad_s=[.2, -.15, .3],
                                   q=[.04, -.03, .02], qdot=[-.02, .05, .03],
                                   wheel_momentum_N_m_s=[.3, -.25, .6]))
    forces = np.random.default_rng(3425).normal(size=(9, 3))*.07
    torque = np.array([.1, -.07, .03])
    alpha, qdd = mechanics(body, state, forces, torque)
    omega, q, qdot = state['omega_rad_s'], state['q'], state['qdot']
    position = body.nodes_m+np.einsum('nik,k->ni', body.modes, q)
    shape_velocity = np.einsum('nik,k->ni', body.modes, qdot)
    velocity = np.cross(omega, position)+shape_velocity
    acceleration = (np.cross(alpha, position)+np.einsum('nik,k->ni', body.modes, qdd)+
                    np.cross(omega, np.cross(omega, position))+2*np.cross(omega, shape_velocity))
    energy_rate = np.sum(body.nodal_masses_kg[:, None]*velocity*acceleration)
    energy_rate += omega@body.extra_inertia@alpha
    energy_rate -= np.sum(state['wheel_momentum_N_m_s']/body.rotor_inertia_kg_m2*torque)
    energy_rate += np.sum(body.omega_squared_s2*q*qdot)
    expected = np.sum(forces*velocity)+torque@(omega-state['wheel_momentum_N_m_s']/body.rotor_inertia_kg_m2)
    expected -= np.sum(body.damping*qdot**2)
    assert energy_rate == pytest.approx(expected, rel=1e-12, abs=1e-13)


def test_combined_nodal_and_wheel_angular_momentum_receives_only_external_torque():
    body = energetic_body()
    state = body.unpack(make_state(body, omega_rad_s=[.2, -.15, .3],
                                   q=[.04, -.03, .02], qdot=[-.02, .05, .03],
                                   wheel_momentum_N_m_s=[.3, -.25, .6]))
    forces = np.random.default_rng(333).normal(size=(9, 3))*.07
    torque = np.array([.1, -.07, .03])
    alpha, qdd = mechanics(body, state, forces, torque)
    omega, q, qdot = state['omega_rad_s'], state['q'], state['qdot']
    position = body.nodes_m+np.einsum('nik,k->ni', body.modes, q)
    shape_velocity = np.einsum('nik,k->ni', body.modes, qdot)
    acceleration = (np.cross(alpha, position)+np.einsum('nik,k->ni', body.modes, qdd)+
                    np.cross(omega, np.cross(omega, position))+2*np.cross(omega, shape_velocity))
    momentum_rate = np.sum(np.cross(position, body.nodal_masses_kg[:, None]*acceleration), axis=0)
    momentum_rate += body.extra_inertia@alpha+np.cross(omega, body.extra_inertia@omega)
    momentum_rate += -torque+np.cross(omega, state['wheel_momentum_N_m_s'])
    assert momentum_rate == pytest.approx(np.sum(np.cross(position, forces), axis=0), rel=1e-12, abs=1e-13)
