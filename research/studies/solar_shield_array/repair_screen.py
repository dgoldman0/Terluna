"""Measure the actual encounter; compare bounded early-return freedoms."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command, independent_path
from protection.dynamics.return_escape import rolled_return, applied_impulse
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.pattern_return import crossing_witnesses, beam_and_rate_screen
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.optical import length
from .packing_screen import parents, load_raw
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE

RUN = ROOT/'research/runs/solar_shield_array/return_repair'


def source_parents(names, extra=()):
    p, sources = parents(names)
    for path in ['protection/dynamics/return_escape.py',
                 'research/studies/solar_shield_array/repair_screen.py', *extra]:
        sources[path] = digest(ROOT/path)
    return p, sources


def write_product(name, out, sources, before):
    out.update(elapsed_s=time.monotonic()-before,
        max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    if any(digest(ROOT/path) != h for path, h in sources.items()):
        raise ValueError('Source changed during computation')
    (HERE/'results'/name).write_text(json.dumps(out, indent=2)+'\n')


def min_geometry(t, state, command):
    best = (float('inf'), None)
    for ti, y, frame in zip(t, state, command(t).as_matrix()):
        value, pair, projection = closest_squares(y[:, :3], frame)
        if value < best[0]:
            best = (value, dict(time_s=float(ti), pair=pair.tolist(),
                projection_m=projection.tolist(), centre_distance_m=float(length(y[pair[0], :3]-y[pair[1], :3]))))
    return dict(minimum_sampled_surface_distance_m=best[0], witness=best[1])


def main():
    signal.alarm(600)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic(); RUN.mkdir(parents=True, exist_ok=True)
    names = ['closure_validation.json', 'closure_return.json', 'closure_sizing.json']
    p, sources = source_parents(names)
    entry = p['closure_validation.json']['runs'][-1]; raw = load_raw(entry)
    initial = raw['state'][0]; start = float(raw['t'][0]); end = entry['intended_end_s']
    ratio = raw['mass_ratio']; env = environment(3.)
    basecmd = return_command(env, start, end, 28800.)
    spline = CubicHermiteSpline(raw['t'], raw['state'][:, :, :3], raw['state'][:, :, 3:])
    out = dict(schema='terluna.research.return-repair-screen/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={n: digest(HERE/'results'/n) for n in names}),
        budget=dict(total_wall_s=3600, stage_wall_s=600, threads=1, address_space_GiB=2),
        accepted_return=False, accepted_cycle=False, accepted_fleet=False,
        scope='Independent unshadowed proposals from actual mass-loaded states. Physical geometry includes smooth attitude transitions. Coupled trials and load accounting decide acceptance.',
        actual_encounter=[], cases=[])
    for t in [start, 45000., 46800., 48600., float(raw['t'][-1])]:
        q, v, f = spline(t), spline(t, 1), basecmd(t).as_matrix()
        d = q[325]-q[347]; dv = v[325]-v[347]
        rate = (d@basecmd(t+1.).as_matrix()-d@basecmd(t-1.).as_matrix())/2+dv@f
        out['actual_encounter'].append(dict(time_s=t, pair=[325, 347],
            displacement_body_m=(d@f).tolist(), displacement_rate_body_m_s=rate.tolist(),
            global_geometry=min_geometry([t], np.c_[q,v][None], basecmd)))
    # The removed/scaled arcs keep later vectors solely to expose endpoint debt.
    specs = [dict(label='unchanged'), dict(label='omit_first', first_scale=0.),
        dict(label='half_first', first_scale=.5), dict(label='three_quarter_first', first_scale=.75),
        *[dict(label=f'delay_{d:g}h', delay_s=d*3600.) for d in [.5, 1., 2., 4.]],
        dict(label='roll_60', angle_deg=60., turn_s=5400.),
        dict(label='roll_75', angle_deg=75., turn_s=5400.),
        dict(label='roll_minus30', angle_deg=-30., turn_s=5400.)]
    for spec in specs:
        u = raw['controls'].copy(); windows = raw['windows'].copy()
        u[:, 0] *= spec.get('first_scale', 1.)
        windows[0] += spec.get('delay_s', 0.)
        command = rolled_return(basecmd, start, end, spec.get('angle_deg', 0.), spec.get('turn_s', 5400.))
        sol = independent_path(env, initial, start, end, command, u, windows, ratio)
        def state(t):
            return sol.sol(np.atleast_1d(t)).T.reshape(-1, len(initial), 6)
        early_t = np.arange(start, 64800.+1., 20.)
        full_t = np.unique(np.r_[np.arange(start, end, 120.), end])
        early = min_geometry(early_t, state(early_t), command)
        whole = min_geometry(full_t, state(full_t), command)
        roots = crossing_witnesses(early_t, state(early_t), command,
            lambda t, ids: state([t])[0, ids], retain=12)
        beam = beam_and_rate_screen(env, early_t, state(early_t), command)
        dense = np.arange(start, 64800.+1., 2.)
        nominal = float(np.trapezoid(rigid_square_load(command,dense)['force_per_mass'],dense))
        first_cost = float(applied_impulse(u, windows, start, 64800.).mean())
        if spec['label'] == 'unchanged':
            baseline_terminal = state([end])[0]
        terminal = state([end])[0]; difference = terminal-baseline_terminal
        delta_first = (u[:, 0]-raw['controls'][:, 0])@basecmd(raw['t'][-1]).as_matrix()
        row = dict(index=len(out['cases']), schedule=spec, early_geometry=early, full_proposal_geometry=whole,
            early_crossings=roots, beam_and_rate=beam,
            mean_first_six_hour_translation_m_s=first_cost, nominal_first_six_hour_attitude_m_s=nominal,
            nominal_first_six_hour_control_m_s=first_cost+nominal,
            scheduled_translation_m_s=float(length(u).sum(axis=1).mean()),
            endpoint_difference_from_unchanged_m=float(length(difference[:,:3]).max()),
            endpoint_velocity_difference_from_unchanged_m_s=float(length(difference[:,3:]).max()),
            limiting_pair_first_burn_change_body_m_s=(delta_first[325]-delta_first[347]).tolist(),
            early_proposal_pass=bool(early['minimum_sampled_surface_distance_m']>350. and
                not roots['crossing_near_encounters'] and beam['minimum_sampled_beam_margin_deg']>0 and
                beam['maximum_sampled_rate_deg_s']<=.1))
        path = RUN/(spec['label']+'_proposal.npz')
        np.savez_compressed(path,t=full_t,state=state(full_t),initial=initial,controls=u,
            windows=windows,mass_ratio=ratio,frames=command(full_t).as_matrix())
        row.update(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path))
        out['cases'].append(row)
        write_product('repair_screen.json',out,sources,before)
        print(json.dumps(dict(case=spec['label'],minimum_m=early['minimum_sampled_surface_distance_m'],
            crossing_count=roots['crossing_near_encounters'],cost_m_s=first_cost+nominal,
            endpoint_difference_m=row['endpoint_difference_from_unchanged_m'],elapsed_s=time.monotonic()-before)),flush=True)
    out['completed']=True
    write_product('repair_screen.json',out,sources,before)


if __name__=='__main__':main()
