"""Independent limiting solutions and invariants for the coupled reduced solver."""
import copy

import numpy as np
import pytest

from protection.dynamics.flexible_fleet import Body, Fleet, canonical_hash, make_state, mechanics, rotation
from protection.dynamics.flexible_optics import receiver_rays


def toy_body(k=2, **options):
    nodes = np.array([[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.], [0., 0., 0.]])
    modes = np.zeros((5, 3, k))
    if k:
        modes[:, 2, 0] = np.array([-1., -1., -1., -1., 4.])/np.sqrt(20.)
    if k > 1:
        modes[:, 2, 1] = np.array([1., -1., 1., -1., 0.])/2
    return Body(nodes, [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
                np.ones(5), modes, np.arange(1, k+1, dtype=float)**2, **options)


def invariants(fleet, result):
    data = [fleet.diagnostics(t, y) for t, y in zip(result.t, result.y.T)]
    energy = np.array([d['energy_minus_work_J'] for d in data])
    angular = np.array([d['angular_minus_impulse_N_m_s'] for d in data])
    return np.max(np.abs(energy-energy[0])), np.max(np.linalg.norm(angular-angular[0], axis=1))


def test_modal_oscillator_matches_closed_solution():
    body = toy_body(1)
    fleet = Fleet([body])
    result = fleet.integrate(fleet.initial_state([dict(q=[.05])]), 9., max_step_s=.15, output_step_s=.2)
    assert result.completed_duration
    assert np.max(np.abs(result.y[13]-.05*np.cos(result.t))) < 2e-10
    assert invariants(fleet, result)[0] < 1e-11


def test_free_rotating_deforming_body_conserves_energy_and_angular_momentum():
    body = toy_body()
    fleet = Fleet([body])
    state = fleet.initial_state([dict(q=[.025, .01], qdot=[.01, -.02], omega_rad_s=[.2, .3, .1])])
    result = fleet.integrate(state, 4., max_step_s=.1, output_step_s=.2, rtol=2e-10, atol=1e-12)
    assert result.completed_duration
    energy_error, angular_error = invariants(fleet, result)
    assert energy_error < 2e-9
    assert angular_error < 2e-9
    assert np.max(np.abs(result.y[10:13, -1]-state[10:13])) > .01


def test_deformation_changes_photon_force_and_attitude_acceleration():
    body = toy_body(reflectivity=np.array([.1, .2, .3, .4]))
    fleet = Fleet([body], pressure_Pa=.01, source_directions=[[.6, 0., .8]])
    flat = fleet.initial_state()
    bent = fleet.initial_state([dict(q=[.08, .04])])
    a, b = fleet.evaluate(0., flat), fleet.evaluate(0., bent)
    assert np.linalg.norm(a['photon_result']['total_force']-b['photon_result']['total_force']) > 1e-4
    assert np.linalg.norm(a['rates'][10:13]-b['rates'][10:13]) > 1e-5
    assert np.linalg.norm(b['rates'][15:17]) > .01


def test_photon_work_and_angular_impulse_close_over_actual_trajectory():
    fleet = Fleet([toy_body()], pressure_Pa=.02, source_directions=[[.6, 0., .8]])
    state = fleet.initial_state([dict(q=[.02, .01], omega_rad_s=[.1, .05, .1], position_m=[2., 1., 3.])])
    result = fleet.integrate(state, 1., max_step_s=.04, output_step_s=.1, rtol=1e-10, atol=1e-12)
    assert result.completed_duration
    assert result.y[20, -1] != 0  # integrated work, k=2
    assert max(invariants(fleet, result)) < 2e-9


def test_central_gravity_node_forces_conserve_total_energy_and_angular_momentum():
    fleet = Fleet([toy_body()], gravitational_parameter_m3_s2=1000.)
    state = fleet.initial_state([dict(position_m=[100., 0., 0.], velocity_m_s=[0., np.sqrt(10), 0.],
                                     omega_rad_s=[.03, .01, .02], q=[.01, 0.])])
    result = fleet.integrate(state, 2., max_step_s=.1, output_step_s=.1, rtol=1e-10, atol=1e-12)
    assert result.completed_duration
    assert max(invariants(fleet, result)) < 1e-7
    assert result.y[1, -1] > 6.


def test_damping_loss_is_accounted_as_work():
    fleet = Fleet([toy_body(1, damping_ratio=.1)])
    state = fleet.initial_state([dict(q=[.05])])
    result = fleet.integrate(state, 5., max_step_s=.1, output_step_s=.2)
    first = fleet.diagnostics(0., state)
    last = fleet.diagnostics(result.t[-1], result.y[:, -1])
    assert last['energy_J'] < first['energy_J']*.5
    assert last['work_J'] < 0
    assert invariants(fleet, result)[0] < 1e-10


def test_wheels_exchange_momentum_and_account_motor_work():
    body = toy_body(rotor_inertia_kg_m2=[.2, .3, .4], torque_limit_N_m=.05, wheel_capacity_N_m_s=10.)
    fleet = Fleet([body], controller=lambda t, i, s: dict(torque_body_N_m=np.array([.1, -.02, .01])))
    state = fleet.initial_state([dict(q=[.01, .01], omega_rad_s=[.03, .02, .04],
                                     wheel_momentum_N_m_s=[.02, -.01, .03])])
    control = fleet.evaluate(0., state)['controls'][0]
    assert control['motor_saturated']
    assert np.allclose(control['motor_torque_N_m'], [.05, -.02, .01])
    result = fleet.integrate(state, 2., max_step_s=.1, output_step_s=.1, rtol=1e-10, atol=1e-12)
    assert result.completed_duration
    assert max(invariants(fleet, result)) < 2e-9
    assert result.y[17, -1] == pytest.approx(.02-.05*2, abs=1e-11)


def test_wheel_capacity_stops_integration_instead_of_clipping_state():
    body = toy_body(0, rotor_inertia_kg_m2=1., torque_limit_N_m=.1, wheel_capacity_N_m_s=.03)
    fleet = Fleet([body], controller=lambda t, i, s: dict(torque_body_N_m=[.1, 0., 0.]))
    result = fleet.integrate(fleet.initial_state(), 2., max_step_s=.025, output_step_s=.05)
    assert result.validity_stop and not result.completed_duration
    assert result.t_events[0][0] < .3
    at_event = fleet.diagnostics(result.t_events[0][0], result.y_events[0][0])
    assert at_event['wheel_utilization'] == pytest.approx(1., abs=1e-10)


def test_trim_is_bounded_photon_authority_and_recomputes_neighbor_transmission():
    bodies = [toy_body(0, trim_max_reflectivity=.7, trim_fraction=1.) for _ in range(2)]
    requested = np.array([10., 0., -100.])  # impossible lateral and normal command
    fleet = Fleet(bodies, pressure_Pa=.01,
                  controller=lambda t, i, s: dict(force_world_N=requested) if i == 0 else {})
    state = fleet.initial_state([dict(position_m=[0., 0., 2.]), {}])
    answer = fleet.evaluate(0., state)
    assert np.all(answer['controls'][0]['trim_command'] <= 1.)
    assert np.all(answer['controls'][0]['trim_command'] >= 0.)
    assert answer['controls'][0]['force_residual_N'][0] == pytest.approx(10.)
    assert answer['photon_result']['tile_forces'][0][0] == pytest.approx(0.)
    assert np.max(answer['photon_result']['face_direct_source_transmission'][1]) < .21
    assert np.linalg.norm(answer['controls'][0]['force_residual_N']) > 90.


def test_actual_motion_changes_receiver_coverage():
    fleet = Fleet([toy_body(0)])
    initial = fleet.initial_state([dict(position_m=[0., 0., 2.], velocity_m_s=[2., 0., 0.])])
    result = fleet.integrate(initial, 2., max_step_s=.5, output_step_s=1.)
    first = fleet.geometry(initial)[-1]
    last = fleet.geometry(result.y[:, -1])[-1]
    assert receiver_rays(first, [[0., 0., 0.]], [0., 0., 1.])['geometric_covered'][0]
    assert not receiver_rays(last, [[0., 0., 0.]], [0., 0., 1.])['geometric_covered'][0]


def test_clearance_bound_includes_declared_hardware_extent():
    fleet = Fleet([toy_body(0, envelope_margin_m=.2), toy_body(0, envelope_margin_m=.3)])
    state = fleet.initial_state([{}, dict(position_m=[4., 0., 0.])])
    assert fleet.diagnostics(0., state)['clearance_lower_bound_m'] == pytest.approx(1.5)
    state = fleet.initial_state([{}, dict(position_m=[1., 0., 0.])])
    assert fleet.diagnostics(0., state)['clearance_lower_bound_m'] == 0.


def test_structural_product_hash_and_evidence_are_enforced():
    body = toy_body()
    names = ('nodes_m', 'triangles', 'nodal_masses_kg', 'modes', 'omega_squared_s2')
    model = {name: getattr(body, name).tolist() for name in names}
    model.update(absorptivity=.2, reflectivity=.3, envelope_margin_m=.5)
    product = dict(schema='terluna.engineering.array-flexible/1', model=model,
                   model_sha256=canonical_hash(model), producer_hashes={'code': 'test'},
                   input_hashes={'support': 'test'}, evidence='Synthetic verification fixture')
    assert Body.from_product(product).mass_kg == 5.
    assert Body.from_product(product).absorptivity == .2
    assert Body.from_product(product).envelope_margin_m == .5
    assert Body.from_product(product, reflectivity=.4).reflectivity == .4
    changed = copy.deepcopy(product)
    changed['model']['nodes_m'][0][0] += 1
    with pytest.raises(ValueError, match='hash mismatch'):
        Body.from_product(changed)
    del product['evidence']
    with pytest.raises(ValueError, match='provenance and evidence'):
        Body.from_product(product)


def test_invalid_centres_modes_and_slope_are_rejected():
    body = toy_body()
    with pytest.raises(ValueError, match='mass centred'):
        Body(body.nodes_m+1., body.triangles, body.nodal_masses_kg, body.modes, body.omega_squared_s2)
    with pytest.raises(ValueError, match='rigid motion'):
        Body(body.nodes_m, body.triangles, body.nodal_masses_kg, body.modes+1., body.omega_squared_s2)
    fleet = Fleet([body])
    with pytest.raises(ValueError, match='small-slope'):
        fleet.integrate(fleet.initial_state([dict(q=[1., 1.])]), 1.)
    with pytest.raises(ValueError, match='rotor model'):
        make_state(body, wheel_momentum_N_m_s=[1., 0., 0.])
    with pytest.raises(ValueError, match='indices must be integers'):
        Body(body.nodes_m, body.triangles.astype(float)+.1, body.nodal_masses_kg,
             body.modes, body.omega_squared_s2)


def test_quaternion_world_rotation_convention():
    q = np.array([np.cos(np.pi/4), 0., 0., np.sin(np.pi/4)])
    assert np.allclose(rotation(q)@[1., 0., 0.], [0., 1., 0.], atol=1e-14)


def test_step_refinement_matches_same_endpoint():
    fleet = Fleet([toy_body()], pressure_Pa=.02, source_directions=[[.6, 0., .8]])
    state = fleet.initial_state([dict(q=[.02, .01], omega_rad_s=[.1, .05, .1])])
    results = [fleet.integrate(state, 2., max_step_s=step, output_step_s=2.,
                               method='RK45', rtol=.01, atol=1e-7).y[:, -1] for step in (.4, .2, .1)]
    assert np.linalg.norm(results[1]-results[2]) < np.linalg.norm(results[0]-results[2])/10


def test_time_dependent_sources_and_spherical_eclipse_use_actual_positions():
    body = toy_body(0)
    source = lambda t: dict(directions=[[0., 0., 1.]], weights=[1.], pressure_Pa=.01*(1+t))
    spheres = lambda t: [dict(position_m=[0., 0., 4.], radius_m=2.)]
    fleet = Fleet([body], source=source, occluders=spheres)
    covered = fleet.initial_state([dict(position_m=[0., 0., 1.])])
    clear = fleet.initial_state([dict(position_m=[5., 0., 1.])])
    assert np.allclose(fleet.evaluate(0., covered)['photon_result']['total_force'], 0.)
    first = fleet.evaluate(0., clear)['photon_result']['total_force']
    second = fleet.evaluate(1., clear)['photon_result']['total_force']
    assert np.linalg.norm(first) > .001
    assert np.allclose(second, 2*first)


def test_third_body_work_and_angular_impulse_close():
    perturb = lambda t: [dict(position_m=[100.+t, 200., 50.], gm_m3_s2=1e4)]
    fleet = Fleet([toy_body()], perturbing_bodies=perturb)
    state = fleet.initial_state([dict(position_m=[1., 2., 3.], velocity_m_s=[.1, -.2, .05],
                                     omega_rad_s=[.03, .02, .04], q=[.01, .02])])
    result = fleet.integrate(state, 2., max_step_s=.1, output_step_s=.2, rtol=1e-10, atol=1e-12)
    assert result.completed_duration
    assert max(invariants(fleet, result)) < 1e-8
    assert abs(fleet.diagnostics(result.t[-1], result.y[:, -1])['work_J']) > 1e-6


def test_inactive_wheel_coordinates_do_not_create_implicit_jacobian_forces():
    body = toy_body(1)
    fleet = Fleet([body])
    initial = fleet.initial_state([dict(q=[.05], omega_rad_s=[.03, .02, .04])])
    perturbed = initial.copy()
    perturbed[15:18] = [.01, -.02, .03]
    assert np.allclose(fleet.rhs(0., initial), fleet.rhs(0., perturbed), atol=0., rtol=0.)
    assert np.allclose(fleet.diagnostics(0., initial)['angular_momentum_N_m_s'],
                       fleet.diagnostics(0., perturbed)['angular_momentum_N_m_s'])
    with pytest.raises(ValueError, match='rotor model'):
        fleet.integrate(perturbed, 1.)
    explicit = fleet.integrate(initial, 1., max_step_s=.1, output_step_s=1., method='DOP853')
    implicit = fleet.integrate(initial, 1., max_step_s=.1, output_step_s=1., method='Radau')
    assert explicit.completed_duration and implicit.completed_duration
    assert np.allclose(explicit.y[:, -1], implicit.y[:, -1], atol=1e-9, rtol=1e-7)
