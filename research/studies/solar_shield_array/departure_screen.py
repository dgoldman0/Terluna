"""Bounded layer-expansion and attitude-timing search at actual terminal states."""
import hashlib
import json
import resource
import signal
import time
from itertools import product

import numpy as np

from shared.provenance import constants_used
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.pattern_departure import (layer_labels, sun_frame,
    departure_command, expansion_maps, expansion_state, frozen_paths)
from protection.dynamics.pattern_return import encounter_scan, crossing_witnesses, beam_and_rate_screen
from protection.dynamics.optical import length
from .return_screen import load_terminal
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE

RUN = ROOT/'research/runs/solar_shield_array/departure'


def main():
    signal.alarm(600)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic(); RUN.mkdir(parents=True, exist_ok=True)
    parent, raw = load_terminal(); initial = raw['state'][-1].copy()
    start = float(raw['t'][-1]); end = start+21600.; env = environment(1.)
    labels = layer_labels(); u = solar_frame(env.at(start))[:, 0]
    paths = ['protection/dynamics/pattern_departure.py', 'protection/dynamics/pattern_return.py',
        'protection/dynamics/parallel_shadow.py', 'research/studies/solar_shield_array/return_screen.py',
        'research/studies/solar_shield_array/departure_screen.py']
    sources = {**parent['producer']['source_hashes'], **{p: digest(ROOT/p) for p in paths}}
    out = dict(schema='terluna.research.departure-screen/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'patterns.json': digest(HERE/'results/patterns.json')}),
        start_s=start, end_s=end, count=len(initial),
        initial_state_sha256=hashlib.sha256(initial.tobytes()).hexdigest(),
        accepted_fleet=False, accepted_return=False, accepted_departure=False,
        continuous_certificate=False, budget_wall_s=600, cases=[], survivors=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/p) != h for p, h in sources.items()):
            raise ValueError('Departure screen source changed')
        (HERE/'results/departure_screen.json').write_text(json.dumps(out, indent=2)+'\n')
    out['frozen_geometry'] = frozen_paths(initial, sun_frame(env, start), labels)
    # Rank axes by the smallest static spacing that passes, then its clearance.
    axes = []
    for axis in [0, 1]:
        for sign in [-1, 1]:
            good = [r for r in out['frozen_geometry'] if r['axis'] == axis and r['sign'] == sign and r['minimum_sampled_clearance_m'] >= 100.]
            axes.append((min((r['extra_layer_spacing_m'] for r in good), default=np.inf), axis, sign))
    axes = sorted(axes)[:2]; out['selected_axes'] = axes
    print(json.dumps(dict(stage='frozen', axes=axes, cases=out['frozen_geometry'])), flush=True)
    times = np.arange(start, end+1, 120.)
    maps_by_burn = {}; stored = []
    for burn in [3600., 7200.]:
        maps = expansion_maps(env, initial, start, end, burn); maps_by_burn[burn] = maps
        for delay, impulse, slew, (_, axis, sign) in product([7200., 10800., 14400.], [.5, 1., 1.5, 2.], [1800., 3600.], axes):
            controls = labels[:, None]*impulse*u
            cap = float(2*length(controls).max()/burn)
            row = dict(index=len(out['cases']), burn_s=burn, delay_s=delay, slew_s=slew,
                axis=axis, sign=sign, layer_impulse_m_s=impulse,
                mean_delta_v_m_s=float(length(controls).mean()), maximum_acceleration_m_s2=cap)
            if cap > .001:
                row['rejected_constraint'] = 'acceleration'; out['cases'].append(row); continue
            state = expansion_state(maps, initial, controls, times)
            command = departure_command(env, start, end, delay, slew, axis, sign)
            row['encounters'] = encounter_scan(times, state, command, retain=2)
            row['beam_and_rate'] = beam_and_rate_screen(env, times, state, command)
            gate = row['beam_and_rate']
            reasons = []
            if row['encounters']['sampled_pair_violations']: reasons.append('clearance')
            if gate['minimum_sampled_beam_margin_deg'] <= 0: reasons.append('beam')
            if gate['maximum_sampled_rate_deg_s'] > .1: reasons.append('attitude_rate')
            row['rejected_constraint'] = ','.join(reasons) if reasons else None
            out['cases'].append(row)
            if not reasons: stored.append((row, controls, command))
            if len(out['cases']) % 12 == 0:
                save(); print(json.dumps(dict(stage='screen', cases=len(out['cases']), survivors=len(stored))), flush=True)
    # Refine at most six strongest cases before requesting a coupled run.
    stored.sort(key=lambda r: (r[0]['mean_delta_v_m_s'], -r[0]['encounters']['minimum_sampled_separating_projection_m']))
    for row, controls, command in stored[:6]:
        maps = maps_by_burn[row['burn_s']]
        dense_times = np.arange(start, end+1, 30.)
        state = expansion_state(maps, initial, controls, dense_times)
        audit = crossing_witnesses(dense_times, state, command,
            lambda t, ids: expansion_state(maps, initial, controls, [t])[0, ids], retain=4)
        fine = encounter_scan(dense_times, state, command, retain=4)
        row['refined_encounters'] = fine; row['crossing_audit'] = audit
        if not fine['sampled_pair_violations'] and not audit['crossing_near_encounters']:
            out['survivors'].append(row['index'])
    path = RUN/'screen.npz'
    np.savez_compressed(path, initial=initial, t=times, labels=labels,
        maps_3600=maps_by_burn[3600.].sol(times).T, maps_7200=maps_by_burn[7200.].sol(times).T)
    out.update(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path),
        stop_reason='Coupled verification required for surviving proposals.' if out['survivors'] else 'Restricted layer-impulse family has no screened safe departure.')
    save(); print(json.dumps(dict(stage='complete', survivors=out['survivors'], elapsed_s=out['elapsed_s'])), flush=True)


if __name__ == '__main__':
    main()
