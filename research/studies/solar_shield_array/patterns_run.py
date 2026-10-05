"""Bounded smooth-roll replay and naturally deforming local-pattern pilot."""
import json
import resource
import signal
import time
from itertools import product

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_formation import solar_frame, lattice, swept_square_conflicts
from protection.dynamics.active_retime import finite_tile_acceleration
from protection.dynamics.coverage_pilot import orientation_matrices
from protection.dynamics.natural_pattern import (NaturalPattern, transported_roll,
    finite_gravity, neighbour_exclusion, local_coverage)
from protection.dynamics.cycling import ray_geometry, visible_sun
from protection.dynamics.fleet import beam_margins
from protection.dynamics.optical import unit, length
from .fleet_run import environment, identities, digest
from .cycling_search import ROOT, HERE, SCENARIO

RUN = ROOT/'research/runs/solar_shield_array/patterns'


def source_ids():
    paths = ['protection/dynamics/active_formation.py', 'protection/dynamics/active_retime.py',
        'protection/dynamics/coverage_pilot.py', 'protection/dynamics/natural_pattern.py',
        'research/studies/solar_shield_array/patterns_run.py']
    return {**identities(), **{p: digest(ROOT/p) for p in paths}}


def attitude_replay(env):
    before = time.monotonic()
    product = json.loads((HERE/'results/pilot_validation.json').read_text())
    for p, h in product['producer']['source_hashes'].items():
        if digest(ROOT/p) != h:
            raise ValueError('Witness producer changed: '+p)
    if constants_changed(product['producer']['constants']):
        raise ValueError('Witness constants changed')
    raw = ROOT/product['raw_path']
    if digest(raw) != product['raw_sha256']:
        raise ValueError('Witness trace changed')
    trace = np.load(raw)
    times, states = trace['attitude_t'], trace['attitude_state'][:, 0]
    frames = []
    for t, state in zip(times, states):
        _, d, _ = finite_tile_acceleration(env, t, state, diagnostics=True)
        frames.append(orientation_matrices(state[None, :], d['normal'], d['sun'])[0])
    command = transported_roll(times, np.array(frames))

    def force(t, y):
        f = command(t).as_matrix()
        rays = ray_geometry(y[:3], env.at(t))
        n = f[:, 2]
        cosine = float(n@(-unit(rays['sun'])))
        if cosine < 0:
            n, cosine = -n, -cosine
        visibility = visible_sun(rays['sun'], rays['earth'], rays['moon'])
        scale = 2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)
        sail = scale*(K.AU/length(rays['sun']))**2*visibility*cosine**2*(9890/10000)**2*n
        gravity = finite_gravity(env, t, y[None, :3], f[:, 0], f[:, 1])[0]
        return np.r_[y[3:], gravity+sail]

    dense_t = np.arange(times[0], times[-1]+.1, .25)
    sol = solve_ivp(force, [times[0], times[-1]], states[0], method='DOP853',
        rtol=2e-12, atol=[1e-5]*3+[1e-9]*3, max_step=5., dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    dense_state = sol.sol(dense_t).T
    margins = []
    for t, y, frame in zip(dense_t[::4], dense_state[::4], command(dense_t[::4]).as_matrix()):
        rays = ray_geometry(y[None, :3], env.at(t))
        margins.append(beam_margins(y[None, :3], rays['sun'], frame[:, 2], env.at(t), 4*K.MOON_RADIUS, 10000.))
    rates = np.rad2deg(length(command(dense_t, 1)))
    acceleration = np.rad2deg(length(command(dense_t, 2)))
    path = RUN/'attitude.npz'
    np.savez_compressed(path, t=dense_t, state=dense_state, frames=command(dense_t).as_matrix(),
                        angular_rate_deg_s=rates, angular_acceleration_deg_s2=acceleration)
    return dict(member=product['attitude_witness']['member'], duration_s=float(np.ptp(times)),
        input_product_sha256=digest(HERE/'results/pilot_validation.json'),
        original_command_rate_deg_s=product['attitude_witness']['max_command_rate_deg_s'],
        transported_max_rate_deg_s=float(rates.max()), rate_cap_deg_s=.1,
        max_angular_acceleration_deg_s2=float(acceleration.max()), rate_sample_step_s=.25,
        minimum_beam_margin_deg=float(np.rad2deg(np.min(margins))), beam_sample_step_s=1.,
        maximum_position_change_m=float(length(sol.sol(times).T[:, :3]-states[:, :3]).max()),
        endpoint_velocity_change_m_s=float(length(sol.y[3:, -1]-states[-1, 3:])),
        electric_delta_v_m_s=0., continuous_certificate=False,
        scope='One 1800 s witness, two-sided rigid square, spline attitude with unmodelled torque/flexibility; no fleet-wide attitude acceptance.',
        raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path), nfev=sol.nfev,
        elapsed_s=time.monotonic()-before)


def reference_shooting(env, duration):
    """Gravity variational proposal along a sail-assisted finite-tile orbit."""
    frame = solar_frame(env.at(0)); u, b = frame[:, :2].T
    if np.dot(np.cross(u, b), env.frames(0)[:, 2]) > 0:
        b = -b
    radius, phase = 15e6, -np.pi/8
    q = radius*(u*np.cos(phase)+b*np.sin(phase))
    v = np.sqrt(K.MOON_GM/radius)*(-u*np.sin(phase)+b*np.cos(phase))
    _, neighbour, valid, _ = lattice(1)
    isolated = NaturalPattern(env, neighbour, valid)

    def rhs(t, y):
        state, phi = y[:6], y[6:].reshape(6, 6)
        eps = 10.
        g = env.gravity(t, state[:3]+np.r_[np.eye(3), -np.eye(3)]*eps)
        gradient = ((g[:3]-g[3:])/(2*eps)).T
        jacobian = np.block([[np.zeros((3, 3)), np.eye(3)], [gradient, np.zeros((3, 3))]])
        return np.r_[state[3:], isolated.forces(t, state[None, :])[0], (jacobian@phi).ravel()]

    sol = solve_ivp(rhs, [0, duration], np.r_[q, v, np.eye(6).ravel()], method='DOP853',
        rtol=2e-12, atol=1e-9, max_step=120., dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    phi = sol.y[6:, -1].reshape(6, 6)
    gradient = np.linalg.solve(phi[:3, 3:], solar_frame(env.at(duration))-phi[:3, :3]@frame)
    return sol, gradient


def affine_state(reference, t, offsets, gradient, frame0):
    y = reference.sol(t)
    phi = y[6:].reshape(6, 6)
    relative = np.c_[offsets@frame0.T, offsets@gradient.T]
    return y[:6]+relative@phi.T


def geometry_trials(env, reference, gradient, duration):
    rows = []
    times = np.linspace(0, duration, 13)
    frame0 = solar_frame(env.at(0))
    for count, pitch, depth, scale in product([19, 21], [8500., 8750., 9000.], [2000., 4000.], [.95, 1.]):
        offsets, neighbour, valid, _ = lattice(count, pitch, depth)
        states = np.array([affine_state(reference, t, offsets, gradient*scale, frame0) for t in times])
        coverage = [local_coverage(y, env.at(t), reference.sol(t)[:6]) for t, y in zip(times, states)]
        # Proposal-only snapshot separation. The accepted coupled run uses
        # swept interval exclusions at a much shorter cadence.
        sep = swept_square_conflicts(times, states, lambda t: solar_frame(env.at(t)),
            acceleration_bound=0., normal_rate_bound=0.)
        row = dict(rows=count, pitch_m=pitch, depth_m=depth, velocity_gradient_scale=scale,
            inventory=count**2, minimum_coverage=min(c['all_sampled_sun_coverage'] for c in coverage),
            coverage=coverage, proposal_chord_separation=sep)
        rows.append(row)
        print(json.dumps(dict(stage='affine_trial', index=len(rows),
            **{k: v for k, v in row.items() if k not in ['coverage', 'proposal_chord_separation']},
            clearance_m=sep['minimum_proven_axis_clearance_m'])), flush=True)
    # A sampled optical failure cannot outrank a passing proposal merely by
    # saving inventory; feasible proposals prefer fewer tiles, then clearance.
    def rank(r):
        margin = r['proposal_chord_separation']['minimum_proven_axis_clearance_m']
        passes = r['minimum_coverage'] >= 1-1e-9 and margin >= 100.
        return (not passes, r['inventory'] if passes else -r['minimum_coverage'], -margin, -r['pitch_m'])
    return rows, sorted(range(len(rows)), key=lambda k: rank(rows[k]))


def propagate_pattern(env, reference, gradient, config, duration, label, step=60., suns=8, rtol=2e-10):
    before = time.monotonic()
    offsets, neighbour, valid, _ = lattice(config['rows'], config['pitch_m'], config['depth_m'])
    frame0 = solar_frame(env.at(0))
    initial = affine_state(reference, 0., offsets, gradient*config['velocity_gradient_scale'], frame0)
    model = NaturalPattern(env, neighbour, valid, suns=suns)
    times = np.arange(0., duration+1, 60.)
    sol = solve_ivp(model.derivative, [0, duration], initial.ravel(), method='DOP853',
        rtol=rtol, atol=np.tile([1e-5]*3+[1e-9]*3, len(initial)), max_step=step,
        t_eval=times, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    states = sol.y.T.reshape(len(times), len(initial), 6)
    light, acceleration, margins, exclusion = [], [], [], []
    max_relative_speed = 0.
    for t, state in zip(times, states):
        a, lit, _ = model.forces(t, state, diagnostics=True)
        light.append(lit); acceleration.append(float(length(a).max()))
        frame = solar_frame(env.at(t)); rays = ray_geometry(state[:, :3], env.at(t))
        margins.append(float(beam_margins(state[:, :3], rays['sun'], -frame[:, 0], env.at(t), 4*K.MOON_RADIUS, 10000.).min()))
        exclusion.append(neighbour_exclusion(state, frame, env.at(t), neighbour, valid))
        max_relative_speed = max(max_relative_speed, float(2*length(state[:, 3:]-state[:, 3:].mean(axis=0)).max()))
    light = np.array(light)
    coverage_times = np.arange(0., duration+1, 900.)
    coverage = [local_coverage(sol.sol(t).reshape(-1, 6), env.at(t), reference.sol(t)[:6]) for t in coverage_times]
    worst = int(np.argmin([c['all_sampled_sun_coverage'] for c in coverage]))
    # Add half-cadence dates plus directions between old limb samples and
    # interior rays at the worst sampled time. Preserve concrete missed rays.
    extra_times = np.unique(np.r_[coverage_times[:-1]+450.,
        np.clip(coverage_times[worst]+np.arange(-300., 301., 60.), 0, duration)])
    adversarial = [local_coverage(sol.sol(t).reshape(-1, 6), env.at(t), reference.sol(t)[:6],
        limb_count=16, rotation=np.pi/16, interior=8) for t in extra_times]
    conflicts = swept_square_conflicts(times, states, lambda t: solar_frame(env.at(t)),
        acceleration_bound=.06, normal_rate_bound=3e-7)
    # Over each interval, the all-pair projected exclusion loses at most a
    # relative velocity bound, centre curvature and common-frame rotation.
    depth_span = float(np.max(np.ptp(states[:, :, :3], axis=1)))
    loss = max_relative_speed*60.+.06*60**2/4+3e-7*60*depth_span
    neighbour_margin = min(exclusion)-loss
    rate_times = np.arange(0., duration+1, 60.)
    frames = np.array([solar_frame(env.at(t)) for t in rate_times])
    angular_steps = np.arccos(np.clip((np.einsum('tij,tij->t', frames[:-1], frames[1:])-1)/2, -1, 1))/60
    predicted = np.array([affine_state(reference, t, offsets, gradient*config['velocity_gradient_scale'], frame0) for t in times])
    initial_relative_v = length(initial[:, 3:]-reference.sol(0)[None, 3:6])
    # Holding the initial Sun-oriented offsets gives a useful reference for
    # preparation/closure, not a claimed required recurring correction.
    frame_d = (solar_frame(env.at(1.))-solar_frame(env.at(-1.)))/2
    preparation = length(initial[:, 3:]-(reference.sol(0)[None, 3:6]+offsets@frame_d.T))
    final_frame_d = (solar_frame(env.at(duration+1))-solar_frame(env.at(duration-1)))/2
    end_target = reference.sol(duration)[:6]+np.c_[offsets@solar_frame(env.at(duration)).T, offsets@final_frame_d.T]
    endpoint_position = length(states[-1, :, :3]-end_target[:, :3])
    endpoint_velocity = length(states[-1, :, 3:]-end_target[:, 3:])
    minimum = min(c['all_sampled_sun_coverage'] for c in coverage+adversarial)
    window_outer_radius = max(c['window_centre_distance_m']+10000. for c in coverage+adversarial)
    ranges = length(states[:, :, :3])
    passed = (minimum >= 1-1e-9 and conflicts['possible_conflict_pairs'] == 0 and neighbour_margin > 0
        and min(margins) > 0 and max(acceleration) < .06 and max(angular_steps) < 3e-7
        and window_outer_radius < 4*K.MOON_RADIUS and ranges.min() > 4*K.MOON_RADIUS)
    path = RUN/(label+'.npz')
    np.savez_compressed(path, t=times, state=states, initial=initial, illumination=light)
    result = dict(config=config, solar_force_directions=suns, max_step_s=step, rtol=rtol,
        duration_s=duration, count=len(initial), optical_mass_kg=len(initial)*env.sigma*10000**2,
        window_outer_radius_m=float(window_outer_radius), range_m=[float(ranges.min()), float(ranges.max())],
        sampled_local_gates_pass=bool(passed), minimum_coverage=float(minimum),
        coverage_times_s=coverage_times.tolist(), coverage=coverage,
        adversarial_times_s=extra_times.tolist(), adversarial_coverage=adversarial,
        separation=conflicts, minimum_instant_nonstencil_margin_m=min(exclusion),
        interval_nonstencil_motion_allowance_m=loss, guarded_nonstencil_margin_m=neighbour_margin,
        max_sampled_acceleration_m_s2=max(acceleration), minimum_beam_margin_deg=float(np.rad2deg(min(margins))),
        max_frame_secant_rate_rad_s=float(max(angular_steps)),
        minimum_illumination=float(light.min()), mean_illumination=float(light.mean()),
        max_affine_position_defect_m=float(length(states[:, :, :3]-predicted[:, :, :3]).max()),
        max_initial_relative_speed_m_s=float(initial_relative_v.max()),
        preparation_delta_v_reference_mean_m_s=float(preparation.mean()),
        preparation_delta_v_reference_max_m_s=float(preparation.max()),
        endpoint_reference_position_defect_max_m=float(endpoint_position.max()),
        endpoint_reference_velocity_defect_mean_m_s=float(endpoint_velocity.mean()),
        endpoint_reference_velocity_defect_max_m_s=float(endpoint_velocity.max()),
        service_electric_delta_v_m_s=0., service_propulsion_energy_J=0.,
        closure_maneuver_executed=False, recurring_cost_measured=False, accepted_fleet=False,
        continuous_coverage_certificate=False, nfev=sol.nfev,
        raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path), elapsed_s=time.monotonic()-before)
    return result, states


def main():
    signal.alarm(1200)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    start = time.monotonic(); RUN.mkdir(parents=True, exist_ok=True)
    sources = source_ids(); env = environment(3.); duration = 6*3600.
    out = dict(schema='terluna.research.natural-pattern-pilot/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources)),
        epoch_tdb=SCENARIO['epoch_tdb'], budget=dict(wall_s=1200, affine_cases=24,
            coupled_cases=3, tighter_replays=1, maximum_tiles=441, memory_GiB=2, raw_GiB=1),
        collection_capacity_W=None, delivered_power_W=None, accepted_fleet=False,
        recurring_cost_measured=False)

    def save():
        out['elapsed_s'] = time.monotonic()-start
        out['max_rss_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        if source_ids() != sources:
            raise ValueError('Source changed during run')
        (HERE/'results/patterns.json').write_text(json.dumps(out, indent=2)+'\n')

    out['attitude'] = attitude_replay(env); save()
    print(json.dumps(dict(stage='attitude', **out['attitude'])), flush=True)
    reference, gradient = reference_shooting(env, duration)
    out['initial_velocity_gradient_s_inv'] = gradient.tolist()
    trials, order = geometry_trials(env, reference, gradient, duration)
    out['affine_trials'], out['candidate_order'] = trials, order
    out['coupled'] = []; save()
    for index in order[:3]:
        config = {k: trials[index][k] for k in ['rows', 'pitch_m', 'depth_m', 'velocity_gradient_scale']}
        result, states = propagate_pattern(env, reference, gradient, config, duration, f'candidate_{index}')
        result['affine_trial'] = index; out['coupled'].append(result); save()
        print(json.dumps(dict(stage='coupled', index=index, elapsed_s=result['elapsed_s'],
            minimum_coverage=result['minimum_coverage'], passed=result['sampled_local_gates_pass'],
            separation=result['separation'], nonstencil_margin=result['guarded_nonstencil_margin_m'])), flush=True)
        if result['sampled_local_gates_pass']:
            replay, independent = propagate_pattern(env, reference, gradient, config, duration,
                f'replay_{index}', step=30., suns=16, rtol=2e-12)
            replay['max_position_difference_m'] = float(length(independent[:, :, :3]-states[:, :, :3]).max())
            replay['max_velocity_difference_m_s'] = float(length(independent[:, :, 3:]-states[:, :, 3:]).max())
            replay['passes_50m'] = replay['max_position_difference_m'] < 50.
            out['replay'] = replay; save()
            print(json.dumps(dict(stage='replay', elapsed_s=replay['elapsed_s'],
                max_position_difference_m=replay['max_position_difference_m'],
                minimum_coverage=replay['minimum_coverage'], passed=replay['sampled_local_gates_pass'])), flush=True)
            break
    save()
    print(json.dumps(dict(stage='complete', elapsed_s=out['elapsed_s'], max_rss_MiB=out['max_rss_MiB'])), flush=True)


if __name__ == '__main__':
    main()
