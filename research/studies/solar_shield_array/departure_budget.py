"""Charge translation and a realizable ideal rigid-body torque implementation."""
import json
import time

import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline

from shared.provenance import constants_used, constants_changed
from protection.dynamics.pattern_departure import departure_command, sun_frame
from protection.dynamics.attitude_load import rigid_square_load, disturbance_torque_bound
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO
from .return_screen import load_terminal


def main():
    before = time.monotonic()
    p = json.loads((HERE/'results/departure_validation.json').read_text())
    r = json.loads((HERE/'results/departure_refine.json').read_text())
    for product in [p, r]:
        for path, h in product['producer']['source_hashes'].items():
            if digest(ROOT/path) != h: raise ValueError('Changed source: '+path)
        if constants_changed(product['producer']['constants']): raise ValueError('Changed constants')
    last = p['runs'][-1]
    if digest(ROOT/last['raw_path']) != last['raw_sha256']: raise ValueError('Raw changed')
    raw = np.load(ROOT/last['raw_path']); env = environment(1.)
    if digest(ROOT/r['raw_path']) != r['raw_sha256']: raise ValueError('Proposal raw changed')
    proposal = np.load(ROOT/r['raw_path'])
    if not np.array_equal(raw['t'], proposal['t']): raise ValueError('Comparison dates differ')
    discrepancy = raw['state'][:, :, :3]-proposal['state'][:, :, :3]
    common_discrepancy = discrepancy.mean(axis=1)
    start, end = p['start_s'], last['end_s']; side = 10000.
    command = departure_command(env, start, p['intended_end_s'], r['delay_s'], r['slew_s'], r['axis'], r['sign'])
    times = np.linspace(start, end, int(np.ceil((end-start)/2))+1)
    load = rigid_square_load(command, times)
    nominal = float(np.trapezoid(load['force_per_mass'], times))
    disturbance = []
    for t, state in zip(raw['t'], raw['state']):
        disturbance.append(disturbance_torque_bound(env, t, state, command(t).as_matrix()))
    disturbance = np.asarray(disturbance)
    allowance = np.trapezoid(disturbance[:, 2], raw['t'], axis=0)
    mass_per_tile = SCENARIO['baseline_areal_mass_kg_m2']*side**2; mass = mass_per_tile*p['count']
    prop = SCENARIO['propulsion']; cant = np.cos(np.deg2rad(prop['cant_deg']))
    watts_per_newton = prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant)
    added_force = np.interp(times, raw['t'], disturbance[:, 2].mean(axis=-1))
    power = mass*(load['force_per_mass']+added_force)*watts_per_newton
    f = np.clip((times-start)/r['burn_s'], 0, 1)
    translation_shape = np.where(times <= start+r['burn_s'], 2/r['burn_s']*np.sin(np.pi*f)**2, 0.)
    power += mass*last['mean_delta_v_m_s']*translation_shape*watts_per_newton
    peak = float(power.max())
    # Count attitude support during the already executed six-hour service as
    # well. Its zero translational impulse never implied zero attitude cost.
    service_product, service_raw = load_terminal()
    service_t = service_raw['t']
    service_command = RotationSpline(service_t,
        Rotation.from_matrix([sun_frame(env, t) for t in service_t]))
    service_dense_t = np.arange(service_t[0], service_t[-1]+1., 2.)
    service_load = rigid_square_load(service_command, service_dense_t)
    service_nominal = float(np.trapezoid(service_load['force_per_mass'], service_dense_t))
    service_disturbance = np.array([disturbance_torque_bound(env, t, state,
        service_command(t).as_matrix())[2] for t, state in zip(service_t, service_raw['state'])])
    service_allowance = float(np.trapezoid(service_disturbance.mean(axis=1), service_t))
    service_power = mass*watts_per_newton*(rigid_square_load(service_command, service_t)['force_per_mass']+
        service_disturbance.mean(axis=1))
    peak = max(peak, float(service_power.max()))
    # Independently verify the library derivative convention against R^T dR.
    test_t = start+r['delay_s']+.37*r['slew_s']; h=.01
    frame = command(test_t).as_matrix()
    skew = frame.T@(command(test_t+h).as_matrix()-command(test_t-h).as_matrix())/(2*h)
    finite = np.array([skew[2, 1], skew[0, 2], skew[1, 0]])
    convention_error = float(length(finite-command(test_t, 1)))
    if convention_error > 1e-9: raise ValueError('Attitude derivative convention mismatch')
    sources = {**p['producer']['source_hashes'],
        'protection/dynamics/attitude_load.py': digest(ROOT/'protection/dynamics/attitude_load.py'),
        'research/studies/solar_shield_array/departure_budget.py': digest(__file__)}
    combined = last['mean_delta_v_m_s']+nominal+float(allowance.mean())
    out = dict(schema='terluna.research.departure-actuation-budget/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={n: digest(HERE/'results'/n) for n in ['departure_validation.json', 'departure_refine.json', 'patterns.json']}),
        count=p['count'], optical_mass_kg=mass, translation_mean_delta_v_m_s=last['mean_delta_v_m_s'],
        attitude_option='Uniform rigid square; three independently commanded opposite edge-midpoint electric force couples, each with 10 km lever arm and zero net force; same exhaust/efficiency/cant as translation.',
        nominal_attitude_equivalent_delta_v_m_s=nominal,
        sampled_disturbance_allowance_mean_m_s=float(allowance.mean()),
        sampled_disturbance_allowance_max_m_s=float(allowance.max()),
        combined_mean_equivalent_delta_v_m_s=combined,
        prior_service_attitude=dict(duration_s=float(service_t[-1]-service_t[0]),
            nominal_equivalent_delta_v_m_s=service_nominal,
            sampled_disturbance_allowance_mean_m_s=service_allowance,
            translation_delta_v_m_s=0., disturbance_sample_step_s=60.,
            interpretation='Ideal rigid-body support estimate on the existing coupled service trace; no hardware-mass replay.'),
        service_plus_departure_equivalent_delta_v_m_s=combined+service_nominal+service_allowance,
        service_plus_departure_propulsion_energy_J=mass*(combined+service_nominal+service_allowance)*watts_per_newton,
        maximum_nominal_torque_per_tile_N_m=float(length(load['torque_per_mass']).max()*mass_per_tile),
        maximum_nominal_total_couple_force_per_tile_N=float(load['force_per_mass'].max()*mass_per_tile),
        maximum_nominal_angular_momentum_per_tile_N_m_s=float(length(load['angular_momentum_per_mass']).max()*mass_per_tile),
        maximum_rigid_rotation_energy_per_tile_J=float(load['kinetic_energy_per_mass'].max()*mass_per_tile),
        derivative_convention_residual_rad_s=convention_error,
        nominal_quadrature_step_s=2., disturbance_quadrature_step_s=30.,
        combined_propulsion_energy_J=mass*combined*watts_per_newton,
        combined_propellant_kg=mass*combined/(prop['exhaust_velocity_m_s']*cant),
        peak_combined_power_with_allowance_W=peak,
        provisional_peak_power_hardware_kg=1.25*peak/300.,
        provisional_peak_power_hardware_fraction=1.25*peak/(300.*mass),
        mass_feedback_propagated=False, structural_realizability_demonstrated=False,
        maximum_linear_proposal_position_error_m=float(length(discrepancy).max()),
        maximum_linear_proposal_centre_error_m=float(length(common_discrepancy).max()),
        maximum_linear_proposal_relative_position_error_m=float(length(discrepancy-common_discrepancy[:, None]).max()),
        disturbance_scope='Pointwise norm envelopes for gravity-gradient and reflected-band torques, allowing arbitrary shadow distributions; their time integral is sampled, not continuously certified. No flexible, absorption, thermal or solar-wind torque.',
        control_scope='Explicit ideal torque authority and zero net force; actuator installation, power collection/storage, plume paths and film flexibility remain unclosed. The mass estimate changes sail acceleration and needs a new coupled run.',
        accepted_fleet=False, full_return_executed=False, handovers_executed=0,
        measured_recurring_delta_v_m_s_per_day=None, fleet_propulsion_W=None,
        collection_capacity_W=None, delivered_power_W=None,
        held_screen_benchmark_TW=238.35, elapsed_s=time.monotonic()-before)
    (HERE/'results/departure_budget.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
