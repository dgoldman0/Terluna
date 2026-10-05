"""Bounded actual-date return search from the saved coupled service terminal."""
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_formation import solar_frame, lattice
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.pattern_return import (natural_reference, shooting_maps,
    proposal_state, attitude_command, minimum_component_impulse,
    encounter_scan, crossing_witnesses, beam_and_rate_screen)
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE

RUN = ROOT/'research/runs/solar_shield_array/returns'


def load_terminal():
    product = json.loads((HERE/'results/patterns.json').read_text())
    for p, h in product['producer']['source_hashes'].items():
        if digest(ROOT/p) != h:
            raise ValueError('Pattern source changed: '+p)
    if constants_changed(product['producer']['constants']):
        raise ValueError('Pattern constants changed')
    raw = ROOT/product['replay']['raw_path']
    if digest(raw) != product['replay']['raw_sha256']:
        raise ValueError('Pattern terminal trace changed')
    return product, np.load(raw)


def next_passage(env, reference, original, start, end):
    frame0 = solar_frame(env.at(0))
    sign = np.sign(np.cross(original[:3], original[3:])@frame0[:, 2])
    def angle(t, state):
        frame = solar_frame(env.at(t))
        return np.arctan2(state[:3]@(sign*frame[:, 1]), state[:3]@frame[:, 0])
    target = angle(0, original)
    times = np.arange(start, end, 600.)
    phases = np.unwrap([angle(t, reference.sol(t)) for t in times])
    sought = target+2*np.pi
    k = np.flatnonzero((phases[:-1] < sought) & (phases[1:] >= sought))[0]
    def residual(t):
        a = angle(t, reference.sol(t))-target
        return np.arctan2(np.sin(a), np.cos(a))
    return float(brentq(residual, times[k], times[k+1]))


def main():
    signal.alarm(1800)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    start_wall = time.monotonic(); RUN.mkdir(parents=True, exist_ok=True)
    parent, raw = load_terminal()
    sources = {**parent['producer']['source_hashes'],
        'protection/dynamics/pattern_return.py': digest(ROOT/'protection/dynamics/pattern_return.py'),
        'research/studies/solar_shield_array/return_screen.py': digest(__file__)}
    env = environment(3.5)
    initial = raw['state'][-1].copy(); start = float(raw['t'][-1])
    reference = natural_reference(env, initial.mean(axis=0), start, 2.6*K.JULIAN_DAY)
    passage = next_passage(env, reference, raw['initial'].mean(axis=0), start, reference.t[-1])
    offsets, _, _, _ = lattice(19, 8500., 4000.); offsets -= offsets.mean(axis=0)
    gradient = np.array(parent['initial_velocity_gradient_s_inv'])*.95
    frame0 = solar_frame(env.at(0))
    out = dict(schema='terluna.research.local-return-screen/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'patterns.json': digest(HERE/'results/patterns.json')}),
        start_s=start, natural_next_passage_s=passage, initial_members=len(initial),
        initial_state_sha256=hashlib.sha256(initial.tobytes()).hexdigest(),
        budget=dict(total_continuation_wall_s=3600, this_screen_wall_s=1800,
            affine_cases=24, encounter_corrections=5, memory_GiB=2),
        reference_scope='Natural common-epoch point reference with unshadowed ideal centre-ray sail; member proposals use gravity variational maps. Full coupled force acceptance remains a separate gate.',
        accepted_return=False, accepted_fleet=False, recurring_cost_measured=False,
        collection_capacity_W=None, delivered_power_W=None, cases=[], refinements=[])

    def save():
        out['elapsed_s'] = time.monotonic()-start_wall
        out['max_rss_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        if any(digest(ROOT/p) != h for p, h in sources.items()):
            raise ValueError('Source changed during screen')
        (HERE/'results/return_screen.json').write_text(json.dumps(out, indent=2)+'\n')

    stored = []
    for date_shift in [-900., 0., 900.]:
        end = passage+date_shift
        maps, windows = shooting_maps(env, reference, start, end)
        times = np.linspace(start, end, int(np.ceil((end-start)/600))+1)
        for quarter in range(4):
            permutation = np.rot90(np.arange(361).reshape(19, 19), quarter).ravel()
            d = offsets[permutation]
            frame = solar_frame(env.at(end))
            target = reference.sol(end)+np.c_[d@frame.T, d@(frame@frame0.T@gradient).T]
            controls, fit = minimum_component_impulse(maps, reference, initial, target, windows)
            for roll in [0., np.pi/4]:
                case = dict(index=len(out['cases']), end_s=end, date_shift_s=date_shift,
                    target_quarter_turn=quarter, coast_roll_deg=float(np.rad2deg(roll)),
                    shooting=fit, windows_s=windows.tolist())
                if controls is not None:
                    command = attitude_command(env, reference, start, end, roll)
                    states = proposal_state(reference, maps, initial, controls, times)
                    case['encounters'] = encounter_scan(times, states, command, retain=12)
                    case['beam_and_rate'] = beam_and_rate_screen(env, times, states, command)
                    case['entry_coverage'] = local_coverage(states[-1], env.at(end), reference.sol(end))
                    case['range_m'] = [float(length(states[:, :, :3]).min()), float(length(states[:, :, :3]).max())]
                    stored.append((case['index'], maps, windows, times, target, controls, command))
                out['cases'].append(case); save()
                print(json.dumps(dict(stage='case', index=case['index'],
                    date_shift_s=date_shift, quarter=quarter, roll_deg=case['coast_roll_deg'],
                    shooting=fit, encounters=case.get('encounters', {}).get('sampled_pair_violations'),
                    beam_and_rate=case.get('beam_and_rate'))), flush=True)
    if not stored:
        out['stop_reason'] = 'Every restricted endpoint problem failed the thrust-cap LP.'; save(); return
    # Choose among beam/rate-admissible proposals, then prioritize encounters
    # and measured Euclidean impulse for the bounded local correction loop.
    def rank(item):
        c = out['cases'][item[0]]; gate = c['beam_and_rate']
        return (gate['minimum_sampled_beam_margin_deg'] <= 0 or gate['maximum_sampled_rate_deg_s'] > .1,
            c['encounters']['sampled_pair_violations'], c['shooting']['mean_delta_v_m_s'])
    index, maps, windows, times, target, controls, command = min(stored, key=rank)
    chosen = out['cases'][index]
    out['selected_case'] = index
    cuts = []; identities = set()
    for iteration in range(6):
        dense_times = np.linspace(start, chosen['end_s'], int(np.ceil((chosen['end_s']-start)/120))+1)
        states = proposal_state(reference, maps, initial, controls, dense_times)
        static = encounter_scan(dense_times, states, command)
        def state_at(t, ids):
            return proposal_state(reference, maps, initial[ids], controls[ids], [t])[0]
        crossings = crossing_witnesses(dense_times, states, command, state_at)
        row = dict(iteration=iteration, sampled=static, between_samples=crossings)
        out['refinements'].append(row)
        print(json.dumps(dict(stage='encounters', iteration=iteration,
            sampled=static['sampled_pair_violations'], crossings=crossings['physical_intersections'])), flush=True)
        if not static['sampled_pair_violations'] and not crossings['crossing_near_encounters']:
            out['stop_reason'] = 'No detected encounter in the refined linear proposal; coupled validation required.'
            break
        if iteration == 5:
            out['stop_reason'] = 'Encounter-correction budget reached; restricted proposal remains unresolved.'
            break
        for w in crossings['witnesses']+static['witnesses']:
            key = (round(w['time_s'], 1), w['i'], w['j'], tuple(np.round(w['axis'], 5)))
            if key in identities:
                continue
            identities.add(key)
            cuts.append((w['time_s'], w['i'], w['j'], np.array(w['axis']), w['required_projection_m']))
        corrected, fit = minimum_component_impulse(maps, reference, initial, target, windows, cuts)
        row['correction'] = fit; save()
        if corrected is None:
            out['stop_reason'] = 'Selected separating-axis cuts cannot be met in this three-arc capped family.'
            break
        controls = corrected
    final_states = proposal_state(reference, maps, initial, controls, dense_times)
    out['final_beam_and_rate'] = beam_and_rate_screen(env, dense_times, final_states, command)
    per_tile = length(controls.reshape(361, 3, 3)).sum(axis=1)
    out['proposal_mean_delta_v_m_s'] = float(per_tile.mean())
    out['proposal_maximum_delta_v_m_s'] = float(per_tile.max())
    out['proposal_mean_delta_v_m_s_per_day'] = float(per_tile.mean()/(chosen['end_s']/K.JULIAN_DAY))
    path = RUN/'screen.npz'
    np.savez_compressed(path, initial=initial, times=dense_times, state=final_states,
        controls=controls, target=target, windows=windows, frames=command(dense_times).as_matrix(),
        reference=reference.sol(dense_times).T, maps=maps.sol(dense_times).T,
        all_controls=np.array([s[5] for s in stored]), case_ids=np.array([s[0] for s in stored]))
    out['raw_path'] = str(path.relative_to(ROOT)); out['raw_sha256'] = digest(path)
    save()
    print(json.dumps(dict(stage='complete', selected=index, reason=out['stop_reason'],
        mean_delta_v_m_s=out['proposal_mean_delta_v_m_s'], elapsed_s=out['elapsed_s'],
        max_rss_MiB=out['max_rss_MiB'])), flush=True)


if __name__ == '__main__':
    main()
