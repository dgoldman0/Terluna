"""Coupled finite-Sun verification of a prepared departure and free coast."""
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_formation import swept_square_conflicts
from protection.dynamics.pattern_departure import departure_command
from protection.dynamics.parallel_shadow import ParallelPattern, minimum_clearance
from protection.dynamics.pattern_return import smooth_arc, beam_and_rate_screen, crossing_witnesses
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO
from .departure_screen import RUN


def main():
    signal.alarm(3000)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic()
    p = json.loads((HERE/'results/departure_refine.json').read_text())
    for path, h in p['producer']['source_hashes'].items():
        if digest(ROOT/path) != h: raise ValueError('Refinement source changed: '+path)
    if constants_changed(p['producer']['constants']): raise ValueError('Constants changed')
    if digest(ROOT/p['raw_path']) != p['raw_sha256']: raise ValueError('Raw proposal changed')
    if not p.get('coupled_candidate'): raise ValueError('No refined candidate')
    raw = np.load(ROOT/p['raw_path']); initial = raw['initial']; impulses = raw['impulses']
    start, end, burn = p['start_s'], p['end_s'], p['burn_s']; env = environment(1.)
    windows = np.array([[start, start+burn]])
    command = departure_command(env, start, end, p['delay_s'], p['slew_s'], p['axis'], p['sign'])
    sources = {**p['producer']['source_hashes'],
        'research/studies/solar_shield_array/departure_run.py': digest(__file__)}
    out = dict(schema='terluna.research.coupled-departure/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'departure_refine.json': digest(HERE/'results/departure_refine.json')}),
        start_s=start, intended_end_s=end, count=len(initial),
        force_scope='Every individual state; finite-area gravity, two-sided finite-Sun sail, dynamic exact parallel-square shadow unions; partial body eclipse rejects this prefix model.',
        accepted_fleet=False, accepted_return=False, accepted_handover=False,
        accepted_departure=False, continuous_certificate=False,
        recurring_cost_measured=False, measured_recurring_delta_v_m_s_per_day=None,
        fleet_propulsion_W=None, collection_capacity_W=None, delivered_power_W=None,
        handovers_executed=0, full_return_executed=False, recovery_executed=False,
        runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path) != h for path, h in sources.items()): raise ValueError('Source changed')
        (HERE/'results/departure.json').write_text(json.dumps(out, indent=2)+'\n')
    solutions = []
    for suns, step, rtol in [(8, 120., 2e-10), (16, 60., 2e-12)]:
        run_before = time.monotonic(); model = ParallelPattern(env, command, suns=suns)
        next_progress = [start+1800.]
        def rhs(t, flat):
            state = flat.reshape(-1, 6)
            if t >= next_progress[0]:
                print(json.dumps(dict(stage='propagating', sources=suns,
                    hours_after_service=(t-start)/3600, elapsed_s=time.monotonic()-run_before)), flush=True)
                next_progress[0] += 1800.
            thrust = smooth_arc(t, windows)[0]*impulses
            return np.c_[state[:, 3:], model.acceleration(t, state)+thrust].ravel()
        def clearance(t, flat):
            return minimum_clearance(flat.reshape(-1, 6), command(t).as_matrix())[0]-100.
        clearance.terminal = True; clearance.direction = -1
        sol = solve_ivp(rhs, [start, end], initial.ravel(), method='DOP853',
            max_step=step, rtol=rtol, atol=np.tile([1e-5]*3+[1e-9]*3, len(initial)),
            events=clearance, dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        times = np.linspace(start, sol.t[-1], int(np.ceil((sol.t[-1]-start)/30))+1)
        states = sol.sol(times).T.reshape(-1, len(initial), 6)
        audit_times = np.linspace(start, sol.t[-1], int(np.ceil((sol.t[-1]-start)/5))+1)
        audit = sol.sol(audit_times).T.reshape(-1, len(initial), 6)
        static = [minimum_clearance(y, command(t).as_matrix()) for t, y in zip(audit_times, audit)]
        worst = int(np.argmin([s[0] for s in static])); minimum, pair, projected = static[worst]
        # The earlier swept bound uses columns [normal, edge, edge]. The
        # departure command uses [edge, edge, normal]; reorder explicitly.
        swept = swept_square_conflicts(audit_times, audit,
            lambda t: command(t).as_matrix()[:, [2, 0, 1]],
            acceleration_bound=.15, normal_rate_bound=np.deg2rad(.1))
        crossings = crossing_witnesses(times, states, command,
            lambda t, ids: sol.sol(t).reshape(-1, 6)[ids], retain=8)
        gate = beam_and_rate_screen(env, audit_times, audit, command)
        endpoint = min(sol.t[-1], start+burn)
        fraction = (endpoint-start)/burn-np.sin(2*np.pi*(endpoint-start)/burn)/(2*np.pi)
        dv = length(impulses)*fraction
        light = []; acceleration = []
        for t, state in zip(times[::4], states[::4]):
            acc, lit, _ = model.acceleration(t, state, diagnostics=True)
            light.append([float(lit.min()), float(lit.mean())])
            acceleration.append(float(length(acc+smooth_arc(t, windows)[0]*impulses).max()))
        coverage = []
        for t in np.arange(start, min(start+p['delay_s'], sol.t[-1])+1, 300.):
            state = sol.sol(t).reshape(-1, 6)
            entry = local_coverage(state, env.at(t), state.mean(axis=0), limb_count=16, rotation=.137, interior=8)
            entry['time_s'] = float(t)
            entry['window_inside_target'] = entry['window_centre_distance_m']+10000. <= 4*K.MOON_RADIUS
            coverage.append(entry)
            if not entry['window_inside_target']: break
        row = dict(solar_sources=suns, max_step_s=step, rtol=rtol, nfev=sol.nfev,
            end_s=float(sol.t[-1]), duration_s=float(sol.t[-1]-start),
            clearance_event=bool(len(sol.t_events[0])), completed_prefix=bool(sol.t[-1] == end),
            transition_completed=bool(sol.t[-1] >= start+p['delay_s']+p['slew_s']),
            coast_after_transition_s=max(0., float(sol.t[-1]-start-p['delay_s']-p['slew_s'])),
            minimum_sampled_clearance_m=minimum, limiting_pair=pair.tolist(),
            worst_clearance_time_s=float(audit_times[worst]), pair_axis_projections_m=projected.tolist(),
            geometry_step_s=float(np.diff(audit_times).max()), swept_separation=swept,
            between_sample_crossings=crossings, beam_and_rate=gate,
            mean_delta_v_m_s=float(dv.mean()), maximum_delta_v_m_s=float(dv.max()),
            maximum_thrust_acceleration_m_s2=float(2*length(impulses).max()/burn),
            sampled_maximum_total_acceleration_m_s2=max(acceleration),
            illumination_diagnostics_step_s=120., minimum_sampled_illumination=min(x[0] for x in light),
            minimum_sampled_mean_illumination=min(x[1] for x in light),
            extra_moving_window_samples=coverage,
            service_scope='Moving local disk only while entirely inside four lunar radii; no replacement pattern or continuous handover is claimed.',
            thrust_arcs_started=1, thrust_arcs_completed=int(sol.t[-1] >= start+burn),
            elapsed_s=time.monotonic()-run_before)
        path = RUN/f'coupled_sun_{suns}.npz'
        np.savez_compressed(path, t=times, state=states, impulses=impulses,
            frames=command(times).as_matrix(), executed_delta_v=dv)
        row.update(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path))
        row['passed_prefix_constraints'] = bool(row['completed_prefix'] and minimum >= 100. and
            swept['possible_conflict_pairs'] == 0 and crossings['crossing_near_encounters'] == 0 and
            gate['minimum_sampled_beam_margin_deg'] > 0 and gate['maximum_sampled_rate_deg_s'] <= .1 and
            row['maximum_thrust_acceleration_m_s2'] <= .001)
        out['runs'].append(row); solutions.append(sol); save()
        print(json.dumps(dict(stage='coupled_complete', sources=suns,
            passed=row['passed_prefix_constraints'], minimum_clearance_m=minimum,
            swept=swept, mean_delta_v_m_s=float(dv.mean()), elapsed_s=row['elapsed_s'])), flush=True)
        if not row['passed_prefix_constraints']: break
    if len(solutions) == 2:
        times = np.arange(start, min(s.t[-1] for s in solutions)+1, 30.)
        a, b = [s.sol(times).T.reshape(-1, len(initial), 6) for s in solutions]
        out['independent_replay'] = dict(maximum_position_difference_m=float(length(a[:, :, :3]-b[:, :, :3]).max()),
            maximum_velocity_difference_m_s=float(length(a[:, :, 3:]-b[:, :, 3:]).max()),
            maximum_proposal_position_difference_m=float(length(b[:, :, :3]-raw['state'][:len(times), :, :3]).max()),
            comparison_step_s=30., placement_gate_m=50.)
        out['accepted_departure'] = bool(all(r['passed_prefix_constraints'] for r in out['runs']) and
            out['independent_replay']['maximum_position_difference_m'] <= 50.)
    fine = out['runs'][-1]; prop = SCENARIO['propulsion']; cant = np.cos(np.deg2rad(prop['cant_deg']))
    mass = len(initial)*SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    impulse = mass*fine['mean_delta_v_m_s']
    out['executed_optical_mass_budget'] = dict(mass_kg=mass, impulse_N_s=impulse,
        propulsion_energy_J=impulse*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant),
        propellant_kg=impulse/(prop['exhaust_velocity_m_s']*cant),
        mass_scope='Base optical inventory only; hardware, attitude actuation and propellant feedback remain unclosed.')
    out['stop_reason'] = ('Prepared departure and three-hour free coast pass the local numerical gates; return, repacking and handover remain open.' if out['accepted_departure'] else 'Candidate failed a coupled or replay constraint; inspect the retained witness.')
    save(); print(json.dumps(dict(stage='complete', accepted_departure=out['accepted_departure'],
        replay=out.get('independent_replay'), elapsed_s=out['elapsed_s'])), flush=True)


if __name__ == '__main__':
    main()
