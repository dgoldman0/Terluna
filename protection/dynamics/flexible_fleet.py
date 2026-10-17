"""Floating-frame reduced membrane dynamics with actual mesh photon loading.

The mechanical equations follow the kinetic energy of lumped moving points,
not a rigid attitude prescribed on a separately forced membrane. Positions are
x + R(r0 + Phi q). The full rotation/modal mass matrix and centrifugal and
Coriolis terms are retained. Elastic energy is the supplied linear modal
potential; finite rotation does not make that constitutive law nonlinear FEM.

Optional orthogonal momentum wheels attach to the mean reference frame. Their
local hub/attachment deformation is NOT resolved. Rotor point mass must already
be in the structural product; the additional tensor accounts only for rotor
self inertia. Photon trim is bounded grey specular reflection, not arbitrary
six-axis force authority. Neither actuator assumption is qualified hardware.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import lsq_linear

from .flexible_optics import MeshTile, photon_forces


def rotation(quaternion):
    """Body-to-world rotation, scalar-first quaternion, normalized for evaluation."""
    quaternion = np.asarray(quaternion, dtype=float)
    norm = np.linalg.norm(quaternion)
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError('Invalid orientation quaternion')
    w, x, y, z = quaternion/norm
    return np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
                     [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
                     [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])


def quaternion_rate(quaternion, omega):
    w, v = quaternion[0], quaternion[1:]
    return .5*np.r_[-v@omega, w*omega+np.cross(v, omega)]


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


@dataclass
class Body:
    nodes_m: object
    triangles: object
    nodal_masses_kg: object
    modes: object
    omega_squared_s2: object
    absorptivity: float = .1
    reflectivity: float = .1
    damping_ratio: float = 0.
    rotor_inertia_kg_m2: object = None
    torque_limit_N_m: object = None
    wheel_capacity_N_m_s: object = None
    trim_max_reflectivity: float | None = None
    trim_fraction: float = .25
    envelope_margin_m: float = 0.

    def __post_init__(self):
        for name in ('nodes_m', 'nodal_masses_kg', 'modes', 'omega_squared_s2'):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        self.triangles = np.asarray(self.triangles)
        if not np.issubdtype(self.triangles.dtype, np.integer):
            raise ValueError('Triangle indices must be integers')
        n, k = len(self.nodes_m), len(self.omega_squared_s2)
        if (self.nodes_m.shape != (n, 3) or self.nodal_masses_kg.shape != (n,)
                or self.modes.shape != (n, 3, k) or (self.nodal_masses_kg <= 0).any()
                or (self.omega_squared_s2 <= 0).any()
                or not all(np.isfinite(a).all() for a in (self.nodes_m, self.nodal_masses_kg,
                                                        self.modes, self.omega_squared_s2))):
            raise ValueError('Invalid reduced body arrays')
        self.mass_kg = float(self.nodal_masses_kg.sum())
        self.scale_m = float(np.max(np.linalg.norm(self.nodes_m, axis=1)))
        weighted = self.nodal_masses_kg[:, None]*self.nodes_m
        if np.linalg.norm(weighted.sum(axis=0))/self.mass_kg > 1e-9*max(1., self.scale_m):
            raise ValueError('Reference nodes must be mass centred')
        translation = np.einsum('n,nik->ik', self.nodal_masses_kg, self.modes)
        rotation_modes = np.sum(self.nodal_masses_kg[:, None, None]*
                                np.cross(self.nodes_m[:, None, :], self.modes.transpose(0, 2, 1)), axis=0)
        if (np.linalg.norm(translation) > 1e-8*np.sqrt(self.mass_kg)
                or np.linalg.norm(rotation_modes) > 1e-8*np.sqrt(self.mass_kg)*max(self.scale_m, 1.)):
            raise ValueError('Elastic modes must be mass orthogonal to rigid motion')
        self.modal_mass = np.einsum('n,nik,nil->kl', self.nodal_masses_kg, self.modes, self.modes)
        if not np.allclose(self.modal_mass, np.eye(k), atol=1e-8, rtol=1e-8):
            raise ValueError('Supplied modes must be mass normalized')
        if (not np.isfinite(self.damping_ratio) or self.damping_ratio < 0
                or not np.isfinite(self.envelope_margin_m) or self.envelope_margin_m < 0):
            raise ValueError('Damping and envelope margin must be finite and nonnegative')
        self.damping = 2*self.damping_ratio*np.sqrt(self.omega_squared_s2)
        MeshTile(self.nodes_m, self.triangles, self.absorptivity, self.reflectivity,
                 1-self.absorptivity-self.reflectivity)
        self.rotors = self.rotor_inertia_kg_m2 is not None
        if self.rotors:
            for name in ('rotor_inertia_kg_m2', 'torque_limit_N_m', 'wheel_capacity_N_m_s'):
                values = np.broadcast_to(np.asarray(getattr(self, name), dtype=float), (3,)).copy()
                if not np.isfinite(values).all() or (values <= 0).any():
                    raise ValueError('Rotor inertia, motor torque and relative momentum capacity must be positive')
                setattr(self, name, values)
            inertia = self.rotor_inertia_kg_m2
            self.extra_inertia = np.diag(.5*(inertia.sum()-inertia))
        else:
            self.extra_inertia = np.zeros((3, 3))
        if self.trim_max_reflectivity is not None:
            if not (self.reflectivity <= self.trim_max_reflectivity <= 1-self.absorptivity
                    and 0 < self.trim_fraction <= 1):
                raise ValueError('Invalid finite reflection trim bounds')
        self.trim_basis = self._trim_basis()
        self.k = k
        # x, v, quaternion, omega, q, qdot, absolute wheel momentum, work,
        # external angular impulse about the chosen nonrotating origin.
        self.state_size = 20+2*k

    @classmethod
    def from_product(cls, product, **options):
        if product.get('schema') != 'terluna.engineering.array-flexible/1':
            raise ValueError('Unsupported structural product schema')
        model = product['model']
        if canonical_hash(model) != product.get('model_sha256'):
            raise ValueError('Structural model hash mismatch')
        if not product.get('producer_hashes') or not product.get('input_hashes') or not product.get('evidence'):
            raise ValueError('Structural product needs producer/input provenance and evidence')
        names = ('nodes_m', 'triangles', 'nodal_masses_kg', 'modes', 'omega_squared_s2')
        for name in ('absorptivity', 'reflectivity', 'damping_ratio', 'envelope_margin_m'):
            if name in model:
                options.setdefault(name, model[name])
        return cls(**{name: model[name] for name in names}, **options)

    def _trim_basis(self):
        """Four radial sectors, outer faces first; partial last face homogenized."""
        faces = self.nodes_m[self.triangles]
        centres = faces.mean(axis=1)
        areas = np.linalg.norm(np.cross(faces[:, 1]-faces[:, 0], faces[:, 2]-faces[:, 0]), axis=1)/2
        basis = np.zeros((len(faces), 4))
        sector = (centres[:, 0] >= 0).astype(int)+2*(centres[:, 1] >= 0)
        for index in range(4):
            chosen = np.flatnonzero(sector == index)
            remaining = self.trim_fraction*areas[chosen].sum()
            for face in chosen[np.argsort(-np.sum(centres[chosen, :2]**2, axis=1))]:
                part = min(areas[face], remaining)
                basis[face, index] = part/areas[face]
                remaining -= part
        return basis

    def unpack(self, state):
        state = np.asarray(state, dtype=float)
        if state.shape != (self.state_size,):
            raise ValueError('Incorrect body state size')
        k = self.k
        return dict(position_m=state[:3], velocity_m_s=state[3:6],
                    quaternion_wxyz=state[6:10], omega_rad_s=state[10:13],
                    q=state[13:13+k], qdot=state[13+k:13+2*k],
                    wheel_momentum_N_m_s=state[13+2*k:16+2*k],
                    work_J=state[16+2*k], angular_impulse_N_m_s=state[17+2*k:20+2*k])


def make_state(body, *, position_m=(0., 0., 0.), velocity_m_s=(0., 0., 0.),
               quaternion_wxyz=(1., 0., 0., 0.), omega_rad_s=(0., 0., 0.),
               q=None, qdot=None, wheel_momentum_N_m_s=None):
    result = np.r_[position_m, velocity_m_s, quaternion_wxyz, omega_rad_s,
                   np.zeros(body.k) if q is None else q,
                   np.zeros(body.k) if qdot is None else qdot,
                   np.zeros(3) if wheel_momentum_N_m_s is None else wheel_momentum_N_m_s,
                   np.zeros(4)].astype(float)
    body.unpack(result)
    if not np.isfinite(result).all() or abs(np.linalg.norm(result[6:10])-1) > 1e-10:
        raise ValueError('Initial state must be finite with a unit quaternion')
    if not body.rotors and np.any(body.unpack(result)['wheel_momentum_N_m_s']):
        raise ValueError('Wheel momentum requires a rotor model')
    return result


def mechanics(body, state, forces_body, motor_torque=None):
    """Solve the kinetic-energy-consistent rotation/modal acceleration equations."""
    if not body.rotors and motor_torque is not None and np.any(motor_torque):
        raise ValueError('Motor torque requires a rotor model')
    q, qdot, omega = state['q'], state['qdot'], state['omega_rad_s']
    phi, mass = body.modes, body.nodal_masses_kg
    r = body.nodes_m+np.einsum('nik,k->ni', phi, q)
    b = np.einsum('nik,k->ni', phi, qdot)
    inertia = np.eye(3)*np.sum(mass*np.sum(r*r, axis=1))-np.einsum('n,ni,nj->ij', mass, r, r)
    cross_mass = np.einsum('n,nki->ik', mass, np.cross(r[:, None, :], phi.transpose(0, 2, 1)))
    matrix = np.block([[inertia+body.extra_inertia, cross_mass],
                       [cross_mass.T, body.modal_mass]])
    bias = np.cross(omega, np.cross(omega, r))+2*np.cross(omega, b)
    torque = np.sum(np.cross(r, forces_body-mass[:, None]*bias), axis=0)
    # The fixed state layout retains three inactive wheel coordinates when
    # rotors are absent. Implicit solvers perturb every coordinate to estimate
    # Jacobians; these padding coordinates must have no mechanical effect.
    wheel_momentum = state['wheel_momentum_N_m_s'] if body.rotors else np.zeros(3)
    torque -= np.cross(omega, body.extra_inertia@omega+wheel_momentum)
    if motor_torque is not None:
        torque += motor_torque
    force_q = np.einsum('nik,ni->k', phi, forces_body-mass[:, None]*bias)
    force_q -= body.omega_squared_s2*q+body.damping*qdot
    acceleration = np.linalg.solve(matrix, np.r_[torque, force_q])
    return acceleration[:3], acceleration[3:]


class Fleet:
    """Coupled propagation with pure controller and environment callbacks.

    A controller returns force_world_N (requested photon force) and/or
    torque_body_N_m (requested combined photon/motor torque). Gravity is
    outside this actuator wrench; controllers may compensate it explicitly.
    Reflection allocation has finite amplitude but instantaneous response:
    bandwidth, switching electricity and thermal transients are not modeled.
    Callbacks must not mutate hidden controller history on adaptive RHS calls.
    """
    def __init__(self, bodies, *, source_directions=((0., 0., 1.),), source_weights=(1.,),
                 pressure_Pa=0., gravitational_parameter_m3_s2=0., controller=None,
                 self_shadow=True, slope_limit=.2, source=None, perturbing_bodies=None,
                 occluders=None):
        self.bodies = list(bodies)
        if not self.bodies or not all(isinstance(b, Body) for b in self.bodies):
            raise ValueError('A fleet needs reduced bodies')
        self.offsets = np.r_[0, np.cumsum([b.state_size for b in self.bodies])]
        self.source_directions = np.asarray(source_directions, dtype=float)
        self.source_weights = np.asarray(source_weights, dtype=float)
        self.pressure_Pa = pressure_Pa
        self.gm = gravitational_parameter_m3_s2
        if not np.isfinite(self.gm) or self.gm < 0 or not np.isfinite(slope_limit) or slope_limit <= 0:
            raise ValueError('Invalid gravity parameter or slope limit')
        photon_forces([], self.source_directions, self.source_weights, self.pressure_Pa)
        self.controller, self.self_shadow, self.slope_limit = controller, self_shadow, slope_limit
        # Optional functions of elapsed seconds. A source returns directions,
        # weights and pressure_Pa in a dict. Perturbers/occluders return lists
        # using the explicit SI schemas in flexible_environment. No ephemeris
        # data are downloaded or silently synthesized here.
        self.source, self.perturbing_bodies, self.occluders = source, perturbing_bodies, occluders
        self.rhs_calls = 0

    def initial_state(self, states=None):
        if states is None:
            states = [{} for _ in self.bodies]
        if len(states) != len(self.bodies):
            raise ValueError('One initial state per body is required')
        return np.concatenate([make_state(body, **state) for body, state in zip(self.bodies, states)])

    def geometry(self, y):
        states, frames, local, world, meshes = [], [], [], [], []
        for index, body in enumerate(self.bodies):
            state = body.unpack(y[self.offsets[index]:self.offsets[index+1]])
            frame = rotation(state['quaternion_wxyz'])
            r = body.nodes_m+np.einsum('nik,k->ni', body.modes, state['q'])
            vertices = state['position_m']+r@frame.T
            states.append(state); frames.append(frame); local.append(r); world.append(vertices)
            meshes.append(MeshTile(vertices, body.triangles, body.absorptivity, body.reflectivity,
                                   1-body.absorptivity-body.reflectivity))
        return states, frames, local, world, meshes

    def evaluate(self, t, y):
        states, frames, local, world, meshes = self.geometry(y)
        environment = (dict(directions=self.source_directions, weights=self.source_weights,
                            pressure_Pa=self.pressure_Pa) if self.source is None else self.source(t))
        directions, weights, pressure = environment['directions'], environment['weights'], environment['pressure_Pa']
        external_transmission = None
        if self.occluders is not None:
            from .flexible_environment import sphere_occlusion
            spheres = self.occluders(t)
            external_transmission = lambda centres, direction: sphere_occlusion(centres, direction, spheres)[:, 0]
        def actual_photons():
            return photon_forces(meshes, directions, weights, pressure, self_shadow=self.self_shadow,
                                 external_transmission=external_transmission)
        photon = actual_photons()
        requests, controls = [], []
        for index, (body, state) in enumerate(zip(self.bodies, states)):
            request = {} if self.controller is None else self.controller(t, index, state)
            request = {} if request is None else request
            for name in ('force_world_N', 'torque_body_N_m'):
                if name in request and (np.shape(request[name]) != (3,) or not np.isfinite(request[name]).all()):
                    raise ValueError('Controller requests must be finite three-vectors')
            requests.append(request)
            command = np.zeros(4)
            allocation = None
            if request and body.trim_max_reflectivity is not None and pressure > 0:
                baseline = photon['nodal_forces'][index]
                lever = local[index]@frames[index].T
                force = baseline.sum(axis=0)
                torque_world = np.cross(lever, baseline).sum(axis=0)
                derivative = photon['reflection_force_derivative'][index]
                face_lever = lever[body.triangles].mean(axis=1)
                gain = body.trim_basis*(body.trim_max_reflectivity-body.reflectivity)
                force_gain = derivative.T@gain
                torque_gain = np.cross(face_lever, derivative).T@gain
                rows, targets = [], []
                if 'force_world_N' in request:
                    rows.append(force_gain)
                    targets.append(np.asarray(request['force_world_N'])-force)
                if 'torque_body_N_m' in request:
                    rows.append(torque_gain/max(body.scale_m, 1.))
                    targets.append((frames[index]@request['torque_body_N_m']-torque_world)/max(body.scale_m, 1.))
                if rows and np.linalg.norm(np.vstack(rows)) > 0:
                    allocation = lsq_linear(np.vstack(rows), np.concatenate(targets), bounds=(0., 1.),
                                            tol=1e-10, max_iter=50)
                    command = allocation.x
                reflection = body.reflectivity+gain@command
                meshes[index] = MeshTile(world[index], body.triangles, body.absorptivity,
                                         reflection, 1-body.absorptivity-reflection)
            controls.append(dict(trim_command=command,
                                 allocation_success=None if allocation is None else bool(allocation.success),
                                 allocation_status=None if allocation is None else int(allocation.status),
                                 allocation_scaled_residual_N=None if allocation is None else float(np.linalg.norm(allocation.fun))))
        # Actual simultaneous commands alter downstream illumination. Recompute
        # the whole scene; the local allocator's Jacobian is not assumed exact.
        if any(np.any(c['trim_command']) for c in controls):
            photon = actual_photons()
        rates = []
        perturbers = None if self.perturbing_bodies is None else self.perturbing_bodies(t)
        for index, (body, state, frame, r, points, request, control) in enumerate(
                zip(self.bodies, states, frames, local, world, requests, controls)):
            photons = photon['nodal_forces'][index]
            gravity = np.zeros_like(points)
            if self.gm:
                radius = np.linalg.norm(points, axis=1)
                if (radius <= 0).any():
                    raise ValueError('Point gravity is singular at the origin')
                gravity = -self.gm*body.nodal_masses_kg[:, None]*points/radius[:, None]**3
            external = photons.copy()
            if perturbers is not None:
                from .flexible_environment import differential_point_acceleration
                external += body.nodal_masses_kg[:, None]*differential_point_acceleration(points, perturbers)
            body_force = (external+gravity)@frame
            photon_torque = np.sum(np.cross(r, photons@frame), axis=0)
            motor = np.zeros(3)
            requested_motor = np.zeros(3)
            if body.rotors and 'torque_body_N_m' in request:
                requested_motor = np.asarray(request['torque_body_N_m'])-photon_torque
                motor = np.clip(requested_motor, -body.torque_limit_N_m, body.torque_limit_N_m)
            alpha, qdd = mechanics(body, state, body_force, motor)
            b = np.einsum('nik,k->ni', body.modes, state['qdot'])
            velocity = state['velocity_m_s']+(np.cross(state['omega_rad_s'], r)+b)@frame.T
            motor_power = 0.
            if body.rotors:
                motor_power = motor@(state['omega_rad_s']-state['wheel_momentum_N_m_s']/body.rotor_inertia_kg_m2)
            power = np.sum(external*velocity)+motor_power-np.sum(body.damping*state['qdot']**2)
            external_torque = np.cross(points, external).sum(axis=0)
            rates.append(np.r_[state['velocity_m_s'], (external+gravity).sum(axis=0)/body.mass_kg,
                               quaternion_rate(state['quaternion_wxyz'], state['omega_rad_s']),
                               alpha, state['qdot'], qdd, -motor, power, external_torque])
            control.update(motor_torque_N_m=motor, requested_motor_torque_N_m=requested_motor,
                           motor_saturated=bool(np.any(np.abs(requested_motor-motor) > 1e-10)),
                           motor_mechanical_power_W=float(motor_power),
                           force_residual_N=np.asarray(request.get('force_world_N', photons.sum(axis=0)))-photons.sum(axis=0),
                           torque_residual_N_m=np.asarray(request.get('torque_body_N_m', photon_torque+motor))-photon_torque-motor)
        return dict(states=states, frames=frames, local=local, world=world, meshes=meshes,
                    photon_result=photon, controls=controls, rates=np.concatenate(rates))

    def rhs(self, t, y):
        self.rhs_calls += 1
        return self.evaluate(t, y)['rates']

    def diagnostics(self, t, y):
        states, frames, local, world, meshes = self.geometry(y)
        energy, work = 0., 0.
        momentum, impulse = np.zeros(3), np.zeros(3)
        slopes, deformation, wheel_use, quaternion_error = [], [], [], []
        for body, state, frame, r, points in zip(self.bodies, states, frames, local, world):
            velocity = np.cross(state['omega_rad_s'], r)+np.einsum('nik,k->ni', body.modes, state['qdot'])
            energy += .5*body.mass_kg*np.sum(state['velocity_m_s']**2)+.5*np.sum(body.nodal_masses_kg[:, None]*velocity**2)
            energy += .5*np.sum(body.omega_squared_s2*state['q']**2)
            energy += .5*state['omega_rad_s']@body.extra_inertia@state['omega_rad_s']
            if body.rotors:
                energy += .5*np.sum(state['wheel_momentum_N_m_s']**2/body.rotor_inertia_kg_m2)
                relative = state['wheel_momentum_N_m_s']-body.rotor_inertia_kg_m2*state['omega_rad_s']
                wheel_use.append(float(np.max(np.abs(relative)/body.wheel_capacity_N_m_s)))
            if self.gm:
                energy -= self.gm*np.sum(body.nodal_masses_kg/np.linalg.norm(points, axis=1))
            spin = np.cross(r, body.nodal_masses_kg[:, None]*velocity).sum(axis=0)
            spin += body.extra_inertia@state['omega_rad_s']
            if body.rotors:
                spin += state['wheel_momentum_N_m_s']
            momentum += frame@spin+np.cross(state['position_m'], body.mass_kg*state['velocity_m_s'])
            work += state['work_J']; impulse += state['angular_impulse_N_m_s']
            faces = r[body.triangles]
            normals = np.cross(faces[:, 1]-faces[:, 0], faces[:, 2]-faces[:, 0])
            slopes.append(float(np.max(np.linalg.norm(normals[:, :2], axis=1)/np.maximum(np.abs(normals[:, 2]), 1e-300))))
            deformation.append(float(np.max(np.linalg.norm(r-body.nodes_m, axis=1))))
            quaternion_error.append(abs(float(np.linalg.norm(state['quaternion_wxyz']))-1))
        clearance = None
        for i in range(len(world)):
            for j in range(i):
                gap = np.maximum(np.maximum(world[i].min(axis=0)-world[j].max(axis=0),
                                             world[j].min(axis=0)-world[i].max(axis=0)), 0.)
                lower = max(0., float(np.linalg.norm(gap))-self.bodies[i].envelope_margin_m-self.bodies[j].envelope_margin_m)
                clearance = lower if clearance is None else min(clearance, lower)
        return dict(energy_J=float(energy), work_J=float(work), energy_minus_work_J=float(energy-work),
                    angular_momentum_N_m_s=momentum, angular_impulse_N_m_s=impulse,
                    angular_minus_impulse_N_m_s=momentum-impulse, max_slope=max(slopes),
                    max_deformation_m=max(deformation), wheel_utilization=max(wheel_use, default=0.),
                    max_quaternion_norm_error=max(quaternion_error), clearance_lower_bound_m=clearance,
                    clearance_rule='Instantaneous mesh AABB separation minus declared isotropic hardware margins; zero is inconclusive, positive is a lower bound. No between-sample guarantee.')

    def integrate(self, y0, duration_s, *, output_step_s=10., max_step_s=10., rtol=1e-8,
                  atol=1e-10, method='DOP853'):
        if not all(np.isfinite(x) and x > 0 for x in (duration_s, output_step_s, max_step_s)):
            raise ValueError('Integration duration and step bounds must be positive')
        for index, body in enumerate(self.bodies):
            initial = body.unpack(y0[self.offsets[index]:self.offsets[index+1]])
            if not body.rotors and np.any(initial['wheel_momentum_N_m_s']):
                raise ValueError('Initial wheel momentum requires a rotor model')
        def validity(t, y):
            diagnostic = self.diagnostics(t, y)
            return min(self.slope_limit-diagnostic['max_slope'], 1-diagnostic['wheel_utilization'])
        validity.terminal, validity.direction = True, -1
        if validity(0., y0) <= 0:
            raise ValueError('Initial state violates small-slope or momentum-capacity bound')
        times = np.linspace(0., duration_s, int(np.ceil(duration_s/output_step_s))+1)
        self.rhs_calls = 0
        result = solve_ivp(self.rhs, (0., duration_s), y0, method=method,
                           t_eval=times, max_step=max_step_s, rtol=rtol, atol=atol,
                           events=validity)
        result['rhs_calls'] = self.rhs_calls
        result['validity_stop'] = bool(result.status == 1)
        result['completed_duration'] = bool(result.success and result.status == 0 and result.t[-1] == duration_s)
        result['stop_reason'] = 'completed' if result['completed_duration'] else 'integrator_failure'
        if result['validity_stop']:
            last_time, last_state = result.t_events[0][-1], result.y_events[0][-1]
            terminal = self.diagnostics(last_time, last_state)
            result['stop_reason'] = ('small_slope_limit' if terminal['max_slope'] >= self.slope_limit*(1-1e-8)
                                     else 'wheel_momentum_capacity')
            # Preserve the crossing state even when it falls between outputs.
            if result.t[-1] < last_time:
                result.t = np.r_[result.t, last_time]
                result.y = np.column_stack((result.y, last_state))
        return result
