"""Coupled geometry gates for at most two energy-first coast proposals."""
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.departure_control import control_modes, mode_impulses
from protection.dynamics.parallel_shadow import ParallelPattern, minimum_clearance
from protection.dynamics.pattern_return import smooth_arc, beam_and_rate_screen, crossing_witnesses
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .departure_validate import audit_geometry
from .energy_screen import RUN
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def main():
    signal.alarm(1800); resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic()
    p = json.loads((HERE/'results/energy_screen.json').read_text())
    for path, expected in p['producer']['source_hashes'].items():
        if digest(ROOT/path) != expected: raise ValueError('Changed source: '+path)
    if constants_changed(p['producer']['constants']): raise ValueError('Changed constants')
    if digest(ROOT/p['raw_path']) != p['raw_sha256']: raise ValueError('Changed raw')
    raw = np.load(ROOT/p['raw_path']); initial = raw['initial']; env = environment(1.)
    start, end, burn = p['start_s'], p['end_s'], p['burn_s']
    modes = control_modes(env, start); windows = np.array([[start, start+burn]])
    sources = dict(p['producer']['source_hashes'])
    for path in ['protection/dynamics/parallel_shadow.py',
                 'research/studies/solar_shield_array/departure_validate.py',
                 'research/studies/solar_shield_array/energy_run.py']:
        sources[path] = digest(ROOT/path)
    out = dict(schema='terluna.research.energy-coast-coupled/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'energy_screen.json': digest(HERE/'results/energy_screen.json')}),
        start_s=start, end_s=end, count=len(initial), runs=[],
        accepted_departure=False, accepted_fleet=False, full_return_executed=False,
        handover_executed=False, hardware_mass_feedback_propagated=False,
        continuous_certificate=False, recurring_cost_measured=False,
        collection_capacity_W=None, delivered_power_W=None)
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/energy.json').write_text(json.dumps(out, indent=2)+'\n')
    save()

    def propagate(case, suns, step, rtol):
        command = coast_command(env, start, end, **case['schedule'])
        impulses = mode_impulses(modes, np.array(case['parameters_m_s']))
        model = ParallelPattern(env, command, suns=suns)
        next_progress = [start+3600.]
        def rhs(t, flat):
            state = flat.reshape(-1, 6)
            if t >= next_progress[0]:
                print(json.dumps(dict(stage='coupled', case=case['index'], suns=suns,
                    hours=(t-start)/3600., elapsed_s=time.monotonic()-before)), flush=True)
                next_progress[0] += 3600.
            return np.c_[state[:, 3:], model.acceleration(t, state)+smooth_arc(t, windows)[0]*impulses].ravel()
        def event(t, flat):
            return minimum_clearance(flat.reshape(-1, 6), command(t).as_matrix())[0]-100.
        event.terminal = True; event.direction = -1
        sol = solve_ivp(rhs, [start, end], initial.ravel(), method='DOP853', max_step=step,
            rtol=rtol, atol=np.tile([1e-5]*3+[1e-9]*3, len(initial)), events=event, dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        times = np.linspace(start, sol.t[-1], int(np.ceil((sol.t[-1]-start)/30))+1)
        states = sol.sol(times).T.reshape(-1, len(initial), 6)
        audit_t = np.linspace(start, sol.t[-1], int(np.ceil((sol.t[-1]-start)/2))+1)
        audit = sol.sol(audit_t).T.reshape(-1, len(initial), 6)
        geometry = audit_geometry(audit_t, audit, command)
        beam = beam_and_rate_screen(env, audit_t, audit, command)
        roots = crossing_witnesses(times, states, command, lambda t, ids: sol.sol(t).reshape(-1, 6)[ids], retain=8)
        maximum_acceleration = max(float(length(rhs(t, state.ravel()).reshape(-1, 6)[:, 3:]).max())
                                   for t, state in zip(times[::4], states[::4]))
        coverage = []
        for t in np.arange(start, start+1201., 300.):
            if t > sol.t[-1] or (case['schedule']['angle_deg'] and t > start+case['schedule']['delay_s']): break
            state = sol.sol(t).reshape(-1, 6)
            entry = local_coverage(state, env.at(t), state.mean(axis=0), limb_count=16, interior=8)
            entry.update(time_s=float(t), window_inside_target=entry['window_centre_distance_m']+10000. <= 4*K.MOON_RADIUS)
            coverage.append(entry)
        path = RUN/f"candidate_{case['index']}_sun_{suns}.npz"
        np.savez_compressed(path, t=times, state=states, impulses=impulses, frames=command(times).as_matrix())
        complete = bool(abs(sol.t[-1]-end) < 1e-6)
        thrust = float(length(impulses).max()*2/burn)
        passed = (complete and not geometry['swept_separation']['possible_conflict_pairs']
            and not roots['crossing_near_encounters'] and beam['minimum_sampled_beam_margin_deg'] > 0
            and beam['maximum_sampled_rate_deg_s'] <= .1 and thrust <= .001
            and maximum_acceleration <= .15
            and all(c['minimum_source_coverage'] >= 1-1e-10 for c in coverage if c['window_inside_target']))
        result = dict(case=case['index'], schedule=case['schedule'], suns=suns, max_step_s=step, rtol=rtol,
            end_s=float(sol.t[-1]), completed=complete, clearance_event=bool(len(sol.t_events[0])),
            geometry=geometry, beam_and_rate=beam, between_sample_crossings=roots,
            maximum_sampled_acceleration_m_s2=maximum_acceleration,
            maximum_thrust_acceleration_m_s2=thrust, translation_mean_m_s=float(length(impulses).mean()),
            translation_account_scope='Full scheduled burn; all accepted traces must complete it.',
            extra_moving_window_samples=coverage, local_gates_pass=bool(passed),
            raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path))
        if not passed:
            result['failed_gates'] = [k for k, value in dict(completed=complete,
                interval_clearance=not geometry['swept_separation']['possible_conflict_pairs'],
                crossing_clearance=not roots['crossing_near_encounters'],
                beam=beam['minimum_sampled_beam_margin_deg'] > 0,
                rate=beam['maximum_sampled_rate_deg_s'] <= .1).items() if not value]
        return result, times, states, audit_t, audit

    for index in p['ranked_proposals'][:2]:
        case = p['cases'][index]
        try:
            row, times, states, _, audit = propagate(case, 8, 120., 2e-10)
        except ValueError as exc:
            out['runs'].append(dict(case=index, exception=str(exc), local_gates_pass=False)); save(); continue
        out['runs'].append(row); save(); del audit
        print(json.dumps(dict(stage='candidate_gates', case=index, passed=row['local_gates_pass'],
            geometry=row['geometry'])), flush=True)
        if not row['local_gates_pass']: continue
        coarse = CubicHermiteSpline(times, states[:, :, :3], states[:, :, 3:])
        fine, ft, fy, at, ay = propagate(case, 16, 60., 2e-12)
        out['runs'].append(fine)
        interpolated = CubicHermiteSpline(ft, fy[:, :, :3], fy[:, :, 3:])
        out['independent_replay'] = dict(maximum_position_difference_m=float(length(ay[:, :, :3]-coarse(at)).max()),
            maximum_velocity_difference_m_s=float(length(ay[:, :, 3:]-coarse(at, 1)).max()),
            reconstruction_position_error_m=float(length(ay[:, :, :3]-interpolated(at)).max()), comparison_step_s=2.)
        out['accepted_departure'] = bool(fine['local_gates_pass'] and out['independent_replay']['maximum_position_difference_m'] < 50.)
        out['selected_case'] = index; save(); break
    out['stop_reason'] = ('Local low-energy schedule passes the ideal-model gates; actual torque/collection account and full-cycle cost remain separate.'
        if out['accepted_departure'] else 'Coupled verification rejected the tested nominal-cost leaders; no cheaper departure accepted.')
    if any(digest(ROOT/p) != h for p, h in sources.items()): raise ValueError('Source changed')
    save(); print(json.dumps(dict(stage='complete', accepted=out['accepted_departure'], elapsed_s=out['elapsed_s'])), flush=True)


if __name__ == '__main__': main()
