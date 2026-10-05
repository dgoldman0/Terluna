"""Bounded joint coast-angle and preparation-impulse search from actual states."""
import json
import resource
import signal
import time

import numpy as np

from shared.provenance import constants_used
from protection.dynamics.energy_coast import coast_command, schedules
from protection.dynamics.pattern_departure import expansion_maps, expansion_state
from protection.dynamics.departure_control import (control_modes, mode_impulses,
    local_cuts, crossing_cuts, minimize_impulse)
from protection.dynamics.pattern_return import encounter_scan, crossing_witnesses, beam_and_rate_screen
from protection.dynamics.attitude_load import rigid_square_load
from .return_screen import load_terminal
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE

RUN = ROOT/'research/runs/solar_shield_array/energy'


def main():
    signal.alarm(1200); resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic(); RUN.mkdir(parents=True, exist_ok=True)
    parent, raw = load_terminal(); initial = raw['state'][-1].copy()
    start = float(raw['t'][-1]); end = start+21600.; burn = 7200.
    env = environment(1.); maps = expansion_maps(env, initial, start, end, burn)
    modes = control_modes(env, start)
    sources = {**parent['producer']['source_hashes']}
    for path in ['protection/dynamics/pattern_departure.py','protection/dynamics/departure_control.py',
                 'protection/dynamics/pattern_return.py','protection/dynamics/energy_coast.py',
                 'protection/dynamics/attitude_load.py',
                 'research/studies/solar_shield_array/return_screen.py',
                 'research/studies/solar_shield_array/energy_screen.py']:
        sources[path] = digest(ROOT/path)
    out = dict(schema='terluna.research.energy-coast-screen/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'patterns.json': digest(HERE/'results/patterns.json')}),
        early_checkpoint='374039df488132f67b1cb213499032100b5510e5',
        budget=dict(total_wall_s=2700, screen_wall_s=1200, maximum_schedules=48,
            corrections_per_schedule=4, address_space_GiB=2, numerical_threads=1),
        start_s=start, end_s=end, burn_s=burn, count=len(initial),
        hard_clearance_m=100., proposal_padding_m=400.,
        objective='Mean translation plus nominal rigid electric-couple impulse. Actual disturbance torque and optical inventory are charged on coupled traces; this proposal ranking is not a complete cycle-cost bound.',
        force_scope='Sun-facing point reference plus gravity variational response; coupled sail/shadow replay required for physical placement.',
        accepted_fleet=False, accepted_departure=False, full_return_executed=False,
        measured_recurring_impulse_m_s_day=None, collection_capacity_W=None,
        delivered_power_W=None, cases=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/energy_screen.json').write_text(json.dumps(out, indent=2)+'\n')
    times = np.arange(start, end+1., 120.); fine_times = np.arange(start, end+1., 30.)
    candidates = []
    for spec in schedules():
        case = dict(index=len(out['cases']), schedule=spec, iterations=[], proposal_pass=False)
        command = coast_command(env, start, end, **spec)
        parameters = np.array([1., 1.5, 1.5]); cuts = []; seen = set()
        for iteration in range(4):
            impulses = mode_impulses(modes, parameters)
            state = expansion_state(maps, initial, impulses, times)
            roots = crossing_witnesses(times, state, command,
                lambda t, ids: expansion_state(maps, initial, impulses, [t])[0, ids], retain=64)
            new = local_cuts(maps, initial, modes, command, times, parameters)+crossing_cuts(
                maps, modes, command, parameters, roots['witnesses'])
            for cut in new:
                key = (round(cut['time_s'], 1), cut['i'], cut['j'], cut['axis'], cut['sign'])
                if key not in seen: cuts.append(cut); seen.add(key)
            parameters, fit = minimize_impulse(modes, parameters, cuts, burn)
            row = dict(iteration=iteration, optimization=fit,
                preceding_crossings=roots['crossing_near_encounters'])
            case['iterations'].append(row)
            if not fit['success']:
                case['rejected_constraint'] = 'local_separation_half_spaces_or_acceleration'; break
            impulses = mode_impulses(modes, parameters)
            fine = expansion_state(maps, initial, impulses, fine_times)
            scan = encounter_scan(fine_times, fine, command, clearance=350., retain=4)
            crossings = crossing_witnesses(fine_times, fine, command,
                lambda t, ids: expansion_state(maps, initial, impulses, [t])[0, ids], retain=8)
            beam = beam_and_rate_screen(env, fine_times, fine, command)
            row.update(refined_encounters=scan, refined_crossings=crossings, beam_and_rate=beam)
            if beam['minimum_sampled_beam_margin_deg'] <= 0 or beam['maximum_sampled_rate_deg_s'] > .1:
                case['rejected_constraint'] = 'reflected_beam_or_attitude_rate'; break
            if not scan['sampled_pair_violations'] and not crossings['crossing_near_encounters']:
                dense = np.arange(start, end+1, 2.)
                nominal = float(np.trapezoid(rigid_square_load(command, dense)['force_per_mass'], dense))
                case.update(proposal_pass=True, parameters_m_s=parameters.tolist(),
                    translation_mean_m_s=fit['mean_delta_v_m_s'], nominal_attitude_m_s=nominal,
                    nominal_sum_m_s=fit['mean_delta_v_m_s']+nominal)
                candidates.append(case); break
        else:
            case['rejected_constraint'] = 'four_local_corrections_exhausted'
        out['cases'].append(case); save()
        print(json.dumps(dict(stage='schedule', index=case['index'], schedule=spec,
            passed=case['proposal_pass'], nominal_sum_m_s=case.get('nominal_sum_m_s'),
            rejected=case.get('rejected_constraint'), elapsed_s=out['elapsed_s'])), flush=True)
    if candidates:
        selected = min(candidates, key=lambda c: c['nominal_sum_m_s'])
        parameters = np.array(selected['parameters_m_s']); impulses = mode_impulses(modes, parameters)
        command = coast_command(env, start, end, **selected['schedule'])
        path = RUN/'screen.npz'
        np.savez_compressed(path, initial=initial, impulses=impulses, t=fine_times,
            state=expansion_state(maps, initial, impulses, fine_times), frames=command(fine_times).as_matrix())
        out.update(selected_case=selected['index'], raw_path=str(path.relative_to(ROOT)),
            raw_sha256=digest(path), ranked_proposals=[c['index'] for c in sorted(candidates,key=lambda c:c['nominal_sum_m_s'])],
            stop_reason='Lowest nominal-cost proposal requires coupled force, torque and collection checks.')
    else:
        out['stop_reason'] = 'No proposal survives this restricted local family.'
    if any(digest(ROOT/p) != h for p, h in sources.items()): raise ValueError('Source changed')
    save()
    print(json.dumps(dict(stage='complete', survivors=len(candidates), selected=out.get('selected_case'),
        elapsed_s=out['elapsed_s'], max_rss_MiB=out['max_rss_MiB'])), flush=True)


if __name__ == '__main__': main()
