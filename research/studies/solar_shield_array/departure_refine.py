"""Use encounter bottlenecks to add transverse preparation to one departure."""
import json
import resource
import signal
import time

import numpy as np

from shared.provenance import constants_used, constants_changed
from protection.dynamics.pattern_departure import departure_command, expansion_maps, expansion_state
from protection.dynamics.departure_control import (control_modes, mode_impulses,
    local_cuts, crossing_cuts, minimize_impulse)
from protection.dynamics.pattern_return import encounter_scan, crossing_witnesses, beam_and_rate_screen
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE
from .departure_screen import RUN


def main():
    signal.alarm(600)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic()
    p = json.loads((HERE/'results/departure_screen.json').read_text())
    for path, h in p['producer']['source_hashes'].items():
        if digest(ROOT/path) != h: raise ValueError('Screen source changed: '+path)
    if constants_changed(p['producer']['constants']): raise ValueError('Screen constants changed')
    if digest(ROOT/p['raw_path']) != p['raw_sha256']: raise ValueError('Screen raw changed')
    initial = np.load(ROOT/p['raw_path'])['initial']; start, end = p['start_s'], p['end_s']
    env = environment(1.); burn, delay, slew, axis, sign = 7200., 7200., 3600., 1, -1
    maps = expansion_maps(env, initial, start, end, burn)
    modes = control_modes(env, start)
    command = departure_command(env, start, end, delay, slew, axis, sign)
    parameters = np.array([1., 1.5, 1.5])
    times = np.arange(start, end+1, 120.)
    sources = {**p['producer']['source_hashes'],
        'protection/dynamics/departure_control.py': digest(ROOT/'protection/dynamics/departure_control.py'),
        'research/studies/solar_shield_array/departure_refine.py': digest(__file__)}
    out = dict(schema='terluna.research.departure-refinement/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'departure_screen.json': digest(HERE/'results/departure_screen.json')}),
        start_s=start, end_s=end, count=len(initial), burn_s=burn, delay_s=delay,
        slew_s=slew, axis=axis, sign=sign, margin_m=400., initial_parameters_m_s=parameters.tolist(),
        accepted_fleet=False, accepted_return=False, accepted_departure=False,
        continuous_certificate=False, iterations=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/departure_refine.json').write_text(json.dumps(out, indent=2)+'\n')
    cuts = []; seen = set()
    for iteration in range(6):
        impulses = mode_impulses(modes, parameters)
        state = expansion_state(maps, initial, impulses, times)
        scan = encounter_scan(times, state, command, retain=8)
        roots = crossing_witnesses(times, state, command,
            lambda t, ids: expansion_state(maps, initial, impulses, [t])[0, ids], retain=128)
        new = local_cuts(maps, initial, modes, command, times, parameters)+crossing_cuts(maps, modes, command, parameters, roots['witnesses'])
        for c in new:
            key = (round(c['time_s'], 1), c['i'], c['j'], c['axis'], c['sign'])
            if key not in seen: cuts.append(c); seen.add(key)
        proposal, fit = minimize_impulse(modes, parameters, cuts, burn)
        row = dict(iteration=iteration, parameters_before_m_s=parameters.tolist(),
            encounters_before=scan, crossings_before=roots, optimization=fit)
        out['iterations'].append(row)
        if not fit['success']:
            out['stop_reason'] = 'Chosen separation half-spaces cannot be optimized in the bounded three-mode control family.'
            save(); break
        parameters = proposal
        impulses = mode_impulses(modes, parameters)
        fine_times = np.arange(start, end+1, 30.)
        fine = expansion_state(maps, initial, impulses, fine_times)
        row['refined_encounters'] = encounter_scan(fine_times, fine, command, clearance=350., retain=8)
        row['refined_crossings'] = crossing_witnesses(fine_times, fine, command,
            lambda t, ids: expansion_state(maps, initial, impulses, [t])[0, ids], retain=8)
        row['beam_and_rate'] = beam_and_rate_screen(env, fine_times, fine, command)
        save(); print(json.dumps(dict(stage='refinement', iteration=iteration, fit=fit,
            refined_violations=row['refined_encounters']['sampled_pair_violations'],
            min_clearance=row['refined_encounters']['minimum_sampled_separating_projection_m'],
            crossings=row['refined_crossings']['physical_intersections'])), flush=True)
        if (not row['refined_encounters']['sampled_pair_violations'] and
            not row['refined_crossings']['crossing_near_encounters'] and
            row['beam_and_rate']['minimum_sampled_beam_margin_deg'] > 0 and
            row['beam_and_rate']['maximum_sampled_rate_deg_s'] <= .1):
            out['stop_reason'] = 'Refined proposal passes sampled constraints; coupled verification required.'
            out['coupled_candidate'] = True; break
    else:
        out['stop_reason'] = 'Six local refinements exhausted; restricted family unresolved.'
    path = RUN/'refine.npz'
    impulses = mode_impulses(modes, parameters)
    final_times = np.arange(start, end+1, 30.)
    np.savez_compressed(path, initial=initial, t=final_times,
        state=expansion_state(maps, initial, impulses, final_times),
        impulses=impulses, parameters=parameters, frames=command(final_times).as_matrix())
    out.update(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path),
        parameters_m_s=parameters.tolist(), cuts=len(cuts))
    if any(digest(ROOT/path) != h for path, h in sources.items()): raise ValueError('Source changed')
    save(); print(json.dumps(dict(stage='complete', reason=out['stop_reason'], elapsed_s=out['elapsed_s'])), flush=True)


if __name__ == '__main__':
    main()
