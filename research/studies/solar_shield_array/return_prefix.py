"""Coupled validation of the first clearance obstruction in a return proposal."""
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation, RotationSpline

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.parallel_shadow import ParallelPattern, minimum_clearance
from protection.dynamics.pattern_return import smooth_arc, beam_and_rate_screen
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE

RUN = ROOT/'research/runs/solar_shield_array/returns'


def main():
    signal.alarm(1200)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    start_wall = time.monotonic()
    product = json.loads((HERE/'results/return_screen.json').read_text())
    for p, h in product['producer']['source_hashes'].items():
        if digest(ROOT/p) != h:
            raise ValueError('Return screen changed: '+p)
    if constants_changed(product['producer']['constants']):
        raise ValueError('Return screen constants changed')
    if digest(ROOT/product['raw_path']) != product['raw_sha256']:
        raise ValueError('Return raw changed')
    raw = np.load(ROOT/product['raw_path'])
    times, initial = raw['times'], raw['initial']
    controls, windows = raw['controls'].reshape(-1, 3, 3), raw['windows']
    command = RotationSpline(times, Rotation.from_matrix(raw['frames']))
    env = environment(3.5)
    sources = {**product['producer']['source_hashes'],
        'protection/dynamics/parallel_shadow.py': digest(ROOT/'protection/dynamics/parallel_shadow.py'),
        'research/studies/solar_shield_array/return_prefix.py': digest(__file__)}
    # First retained between-sample crossing is an upper limit for this focused
    # prefix. Events also inspect all tiles, so an earlier failure can stop it.
    witnesses = product['refinements'][-1]['between_samples']['witnesses']
    predicted = min(w['time_s'] for w in witnesses)
    end = min(predicted+600., times[-1])
    out = dict(schema='terluna.research.return-clearance-prefix/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={'return_screen.json': digest(HERE/'results/return_screen.json')}),
        selected_case=product['selected_case'], actual_initial_members=len(initial),
        prefix_end_limit_s=end, predicted_retained_crossing_s=predicted,
        force_scope='All individual states, finite-area gravity, dynamic exact parallel-square shadow unions, two-sided finite-Sun force; body eclipse explicitly excluded from this short prefix.',
        accepted_return=False, accepted_handover=False, accepted_fleet=False,
        recurring_cost_measured=False, full_return_executed=False, runs=[])
    dense_solutions = []
    for suns, step, rtol in [(8, 60., 2e-10), (16, 30., 2e-12)]:
        before = time.monotonic()
        model = ParallelPattern(env, command, suns=suns)
        def rhs(t, flat):
            state = flat.reshape(-1, 6)
            thrust = np.einsum('j,njk->nk', smooth_arc(t, windows), controls)
            return np.c_[state[:, 3:], model.acceleration(t, state)+thrust].ravel()
        def clearance(t, flat):
            return minimum_clearance(flat.reshape(-1, 6), command(t).as_matrix())[0]-100.
        clearance.terminal = True; clearance.direction = -1
        sol = solve_ivp(rhs, [times[0], end], initial.ravel(), method='DOP853',
            max_step=step, rtol=rtol, atol=np.tile([1e-5]*3+[1e-9]*3, len(initial)),
            events=clearance, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        ts = np.linspace(times[0], sol.t[-1], int(np.ceil((sol.t[-1]-times[0])/30))+1)
        state = sol.sol(ts).T.reshape(len(ts), len(initial), 6)
        minimum, pair, projected = minimum_clearance(state[-1], command(ts[-1]).as_matrix())
        electric = np.array([length(np.einsum('j,njk->nk', smooth_arc(t, windows), controls)) for t in ts])
        dv = np.trapezoid(electric, ts, axis=0)
        _, light, _ = model.acceleration(ts[-1], state[-1], diagnostics=True)
        row = dict(solar_sources=suns, max_step_s=step, rtol=rtol, nfev=sol.nfev,
            clearance_event=bool(len(sol.t_events[0])), event_time_s=float(ts[-1]),
            seconds_after_service=float(ts[-1]-times[0]), minimum_clearance_m=minimum,
            limiting_pair=None if pair is None else pair.tolist(),
            pair_axis_projections_m=None if projected is None else projected.tolist(),
            mean_executed_delta_v_m_s=float(dv.mean()), max_executed_delta_v_m_s=float(dv.max()),
            max_acceleration_command_m_s2=float(electric.max()),
            minimum_final_illumination=float(light.min()), mean_final_illumination=float(light.mean()),
            beam_and_rate=beam_and_rate_screen(env, ts, state, command),
            elapsed_s=time.monotonic()-before)
        path = RUN/f'prefix_sun_{suns}.npz'
        np.savez_compressed(path, t=ts, state=state, executed_delta_v=dv)
        row.update(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path))
        out['runs'].append(row); dense_solutions.append(sol)
        out['elapsed_s'] = time.monotonic()-start_wall
        out['max_rss_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        (HERE/'results/return_prefix.json').write_text(json.dumps(out, indent=2)+'\n')
        print(json.dumps(dict(stage='prefix', **row)), flush=True)
    common_times = np.linspace(times[0], min(s.t[-1] for s in dense_solutions), 201)
    a, b = [s.sol(common_times).T.reshape(-1, len(initial), 6) for s in dense_solutions]
    out['independent_replay'] = dict(max_position_difference_m=float(length(a[:, :, :3]-b[:, :, :3]).max()),
        max_velocity_difference_m_s=float(length(a[:, :, 3:]-b[:, :, 3:]).max()),
        clearance_event_time_difference_s=abs(out['runs'][0]['event_time_s']-out['runs'][1]['event_time_s']),
        same_limiting_pair=out['runs'][0]['limiting_pair'] == out['runs'][1]['limiting_pair'])
    if any(digest(ROOT/p) != h for p, h in sources.items()):
        raise ValueError('Source changed during prefix replay')
    out['elapsed_s'] = time.monotonic()-start_wall
    (HERE/'results/return_prefix.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(stage='complete', replay=out['independent_replay'], elapsed_s=out['elapsed_s'])), flush=True)


if __name__ == '__main__':
    main()
