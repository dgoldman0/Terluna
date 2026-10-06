"""Budgeted, resumable coupled search; every numerical stage writes a ledger.

Run --help for portable prepare, benchmark, sensitivity, solve and replay
commands. Only explicitly named modes are differentiated; absent columns
are never filled with an approximate force response.
"""
import argparse
import hashlib
import json
import os
import resource
import signal
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from protection.dynamics.ephemeris import Ephemeris, MANIFEST, BODY_GM
from protection.dynamics.coupled_control import contact_set, contact_values
from protection.dynamics.optical import length
from .coupled_common import limit, SPECIFIC
from .coupled_cycle import transplant
from .joint_common import spline, state_at
from .fleet_run import digest
from .cycling_search import SCENARIO
from .expanded_model import ExpandedExperiment, PACKAGE, ROOT, HERE, NAMES, DIM

RUN = ROOT/'research/runs/solar_shield_array/expanded_cycle'
SCALE = np.array([1e5]*3+[5.]*3)
SERVICE_J = 2.64802552665676e12  # replaced from source-bound account in prepare


def check_budget():
    from . import expanded_model
    return expanded_model.check_budget()


def native(x):
    if isinstance(x, dict): return {k: native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    return x


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(native(data), indent=2, allow_nan=False)+'\n'); tmp.replace(path)


def raw_read(path):
    with np.load(path, allow_pickle=False) as z: return {k: z[k] for k in z.files}


def save_case(name, row, raw, folder=None):
    check_budget()
    if folder is None: folder = RUN
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder/(name+'.npz'), **raw)
    row['raw_sha256'] = digest(folder/(name+'.npz'))
    write(folder/(name+'.json'), row)


def prepare():
    """Use pinned DE440 and the measured historical full-force benchmark once."""
    from .coupled_stable import StableExperiment
    PACKAGE.mkdir(exist_ok=True)
    eph = Ephemeris(epoch=SCENARIO['epoch_tdb'])
    samples = eph.sample(np.arange(-3600., 6*86400.+7200., 1800.)); eph.close()
    flat = {k: v for k, v in samples.items() if k not in ('positions', 'velocities')}
    flat.update({'p_'+n: samples['positions'][n] for n in BODY_GM})
    flat.update({'v_'+n: samples['velocities'][n] for n in BODY_GM})
    np.savez_compressed(PACKAGE/'ephemeris.npz', **flat)
    ex = StableExperiment(); base = raw_read(RUN/'benchmark.npz'); cmd = ex.command(np.zeros(8))
    dates = np.array([ex.arrival, ex.arrival+21600., ex.arrival+46800.])
    states = state_at(spline(base), dates)
    targets = [transplant(y, cmd(t).as_matrix(), cmd(t, 1), cmd(ex.arrival).as_matrix(),
                         cmd(ex.arrival, 1), states[0].mean(axis=0))
               for y, t in [(ex.arrays['initial_epoch_state'], 0.)]]
    residual = states[0]-targets[0]; f = cmd(ex.arrival).as_matrix()
    late = np.zeros((361, 3, 6))
    for group in range(2):
        r = residual[:, 3*group:3*group+3]@f
        for axis in range(3):
            p = -r[:, axis]; p -= np.average(p, weights=ex.mass)
            p /= max(abs(p))
            late[:, :, 3*group+axis] = p[:, None]*f[:, axis]
    np.savez_compressed(PACKAGE/'basis.npz', late=late,
                        reference_centre=states.mean(axis=1))
    service = json.loads((HERE/'results/coupled_validation.json').read_text())['runs'][0]
    config = dict(schema='terluna.expanded-cycle-inputs/1', seed=20261006, initial_commit='762ce4d810e07ddc0788f315212c5a6f6bdca5fe',
                  epoch_tdb=SCENARIO['epoch_tdb'], scenario_sha256=digest(HERE/'scenario.json'),
                  ephemeris_manifest=json.loads(MANIFEST.read_text()),
                  compact_ephemeris_sha256=digest(PACKAGE/'ephemeris.npz'),
                  basis_sha256=digest(PACKAGE/'basis.npz'), names=NAMES,
                  initial_service_energy_J=service['electrical_energy_J'],
                  limits=dict(prefix_energy_J=51.859e12, through_service_energy_J=100e12,
                              clearance_m=100., numerical_reserve_m=50., translation_m_s2=.001,
                              attitude_rate_deg_s=.1, recurrence_position_m=50.,
                              recurrence_velocity_m_s=.01, recurrence_attitude_rad=1e-4),
                  supply=dict(generation=None, storage=None, unmet_demand=None),
                  recurrence='Two subsequent actual-date cycles without state resets are required in addition to the finite sequence; not a proof of indefinite recurrence.',
                  hardware_mass_kg=float(ex.mass.sum()), hardware_resized=False,
                  reference='Same epoch and loaded seed; acquisition patterns follow baseline arrival component errors, normalized to 1 m/s per coefficient.')
    write(PACKAGE/'inputs.json', config)
    print('Compact inputs prepared', flush=True)


def aligned(ex, raw, x, times):
    shift = 1800*x[16]
    dates = times+shift*np.clip((times-30*3600)/(ex.arrival0-30*3600), 0., 1.)
    return state_at(spline(raw), dates), ex.command(x)(dates).as_matrix()


def energy_parts(ex, raw, x):
    times = raw['t']; dates = ex.dates(x)
    initial = json.loads((PACKAGE/'inputs.json').read_text())['initial_service_energy_J']
    out = []
    for end in [43200., dates[1], dates[-1]]:
        sel = times <= end+1e-7
        out.append(initial+float(np.trapezoid(raw['power_W'][sel].sum(axis=1), times[sel])))
    return np.array(out)


def feature(ex, raw, x, base, cuts):
    from protection.dynamics.natural_pattern import local_coverage
    y, f = aligned(ex, raw, x, base['t'])
    service = []
    for t in np.linspace(ex.dates(x)[0], ex.dates(x)[1], 13):
        check_budget()
        state = state_at(spline(raw), t)
        cover = local_coverage(state, ex.env.at(t), state.mean(axis=0),
                               limb_count=16, interior=8, rotation=.137)
        service.append(1.-cover['minimum_source_coverage'])
    return dict(terminal=(ex.terminal(raw, x)/SCALE).ravel(),
                service=np.array(service),
                gaps=contact_values(y, f, cuts),
                attitude_peak=(raw['attitude_acceleration']*ex.mass*SPECIFIC).max(axis=0),
                energy=energy_parts(ex, raw, x), state=y, frames=f)


def benchmark(ex, max_step):
    row, raw, _ = ex.evaluate(np.zeros(DIM), max_step=max_step)
    row['reference_120s_comparison'] = None
    if (RUN/'benchmark.npz').exists():
        reference = raw_read(RUN/'benchmark.npz')
        dates = reference['t']; delta = state_at(spline(raw), dates)-reference['state']
        row['reference_120s_comparison'] = dict(maximum_position_m=float(length(delta[:, :, :3]).max()),
            maximum_velocity_m_s=float(length(delta[:, :, 3:]).max()),
            first_event_time_difference_s=abs(row['clearance_events_s'][0]-json.loads((RUN/'benchmark.json').read_text())['clearance_events_s'][0]))
    row['cycle'] = ex.summary(raw, np.zeros(DIM)); row['energy_parts_J'] = energy_parts(ex, raw, np.zeros(DIM))
    save_case('base', row, raw)
    cuts = contact_set(raw['state'], raw['frames'], reserve=1000.)
    np.savez_compressed(RUN/'contacts.npz', cuts=cuts)
    print('base', row['cpu_s'], row['reference_120s_comparison'], 'cuts', len(cuts), flush=True)


def sensitivity(ex, columns, h, max_step, deadline, reserve):
    base = raw_read(RUN/'base.npz'); cuts = raw_read(RUN/'contacts.npz')['cuts']
    for j in columns:
        path = RUN/f'direction_{j:02d}.json'
        if path.exists(): continue
        cost = 160. if j < 10 else 65.
        if deadline-time.process_time() < reserve+cost:
            print('Reserved budget before column', j, flush=True); break
        x = np.zeros(DIM); x[j] = h
        # An exact baseline output is saved at each late-control onset.
        start = ex.onset(j); initial = state_at(spline(base), start)
        row, suffix, _ = ex.evaluate(x, initial, start, max_step=max_step)
        prefix = base['t'] < start
        raw = {k: np.r_[base[k][prefix], suffix[k]] if k in ['t', 'state', 'frames',
               'distance_m', 'power_W', 'attitude_acceleration', 'translation',
               'first_intercept_W', 'redirected_W', 'torque_N_m'] else v
               for k, v in suffix.items()}
        row.update(column=j, name=NAMES[j], h=h, restart_onset_s=start,
                   cycle=ex.summary(raw, x), energy_parts_J=energy_parts(ex, raw, x))
        if start > 21600.:
            row['inherited_first_failed_clearance_s'] = json.loads((RUN/'base.json').read_text())['clearance_events_s'][0]
            row['restoration_only'] = True
        save_case(f'direction_{j:02d}', row, raw)
        print('direction', j, NAMES[j], row['cpu_s'], flush=True)


def build_model(ex):
    available = [j for j in range(DIM) if (RUN/f'direction_{j:02d}.json').exists()]
    if (RUN/'model.npz').exists():
        cached = raw_read(RUN/'model.npz')
        if list(cached['columns']) == available and 'service_jac' in cached:
            return cached
    base = raw_read(RUN/'base.npz'); cuts = raw_read(RUN/'contacts.npz')['cuts']
    b = feature(ex, base, base['x'], base, cuts)
    columns = []; derivatives = {k: [] for k in ['terminal', 'gaps', 'service', 'attitude_peak', 'energy']}
    for j in range(DIM):
        p = RUN/f'direction_{j:02d}.json'
        if not p.exists(): continue
        row = json.loads(p.read_text()); raw = raw_read(p.with_suffix('.npz'))
        f = feature(ex, raw, raw['x'], base, cuts); columns.append(j)
        for k in derivatives: derivatives[k].append((f[k]-b[k])/row['h'])
    model = dict(columns=np.array(columns), cuts=cuts, cut_times_s=base['t'][cuts[:, 0].astype(int)], x=base['x'],
                 steps=np.array([json.loads((RUN/f'direction_{j:02d}.json').read_text())['h'] for j in columns]))
    for k in derivatives:
        model[k] = b[k]; model[k+'_jac'] = np.stack(derivatives[k], axis=-1)
    np.savez_compressed(RUN/'model.npz', **model)
    return model


def rank_report(model):
    out = {}
    for name, scale in [('terminal', 1.), ('gaps', 1000.), ('service', 1.)]:
        a = model[name+'_jac']/scale
        u, s, vt = np.linalg.svd(a, full_matrices=False)
        residual = np.minimum((model[name]-150.)/scale, 0.) if name == 'gaps' else model[name]/scale
        ranks = {str(r): int(np.count_nonzero(s > s[0]*r)) for r in [1e-2, 1e-3, 1e-4, 1e-6]}
        r = ranks['0.001']; projection = u[:, :r]@(u[:, :r].T@residual)
        out[name] = dict(singular_values=s, effective_ranks=ranks,
            residual_norm=float(np.linalg.norm(residual)),
            fraction_residual_norm_outside_rank_1e3=float(np.linalg.norm(residual-projection)/max(np.linalg.norm(residual), 1e-30)),
            column_norms=np.linalg.norm(a, axis=0))
    out['measured_columns'] = model['columns']; out['unmeasured_columns'] = sorted(set(range(DIM))-set(model['columns']))
    for label, rows, cols in [
        ('early_clearance', model['cut_times_s'] <= 15*3600., model['columns'] < 10),
        ('legacy_measured_terminal', None, model['columns'] < 6),
        ('expanded_early_terminal', None, model['columns'] < 10),
        ('arrival_terminal', None, model['columns'] >= 10)]:
        a = model['gaps_jac'][rows][:, cols]/1000. if rows is not None else model['terminal_jac'][:, cols]
        residual = np.minimum((model['gaps'][rows]-150.)/1000., 0.) if rows is not None else model['terminal']
        if not a.shape[1]: continue
        u, s, _ = np.linalg.svd(a, full_matrices=False)
        rank = int(np.count_nonzero(s > s[0]*1e-3)); projection = u[:, :rank]@(u[:, :rank].T@residual)
        out[label] = dict(singular_values=s, effective_rank_1e3=rank,
            fraction_residual_norm_outside_rank_1e3=float(np.linalg.norm(residual-projection)/max(np.linalg.norm(residual), 1e-30)))
    out['scope'] = 'Finite-difference local rank of the measured columns at the saved anchor; no global feasibility inference.'
    return out


def solve(ex, model, radius=1., energy_factor=1., capacity_factor=1., prefix_factor=1., warm_start=None):
    cols = model['columns']; n = len(cols); x0 = model['x']; initial = json.loads((PACKAGE/'inputs.json').read_text())['initial_service_energy_J']
    # The energy Jacobian includes tiny one-sided thrust norms. Remove their
    # linear term and add the actual vector norm cost at every trial instead.
    def control_cost(x):
        controls, windows = ex.controls(x)
        costs = []
        for end in [43200., ex.dates(x)[1], ex.dates(x)[-1]]:
            f = np.clip((end-windows[:, 0])/(windows[:, 1]-windows[:, 0]), 0., 1.)
            # Exact partial integral of smooth_arc's normalized sin² pulse.
            f = f-np.sin(2*np.pi*f)/(2*np.pi)
            costs.append(float(np.sum(length(controls)*f[None, :]*ex.mass[:, None]*SPECIFIC)))
        return np.array(costs)
    translation_jac = []
    for index, j in enumerate(cols):
        h = model['steps'][index] if 'steps' in model else json.loads((RUN/f'direction_{j:02d}.json').read_text())['h']
        x = x0.copy(); x[j] += h; translation_jac.append(control_cost(x)/h)
    att_jac = model['energy_jac']-np.array(translation_jac).T
    def quantities(d):
        check_budget()
        x = x0.copy(); x[cols] += d
        energy = model['energy']+att_jac@d+control_cost(x)
        u, windows = ex.controls(x)
        # Conservative non-overlapping burn peak plus whole-sequence attitude
        # peak. This never borrows aggregate capacity between members.
        peak_acc = (length(u)*2./(windows[:, 1]-windows[:, 0])[None, :]).max(axis=1)
        attitude = np.maximum(model['attitude_peak']+model['attitude_peak_jac']@d, 0.)
        power = attitude+peak_acc*ex.mass*SPECIFIC
        terminal = model['terminal']+model['terminal_jac']@d
        gaps = model['gaps']+model['gaps_jac']@d
        return energy, power, peak_acc, terminal, gaps
    def obj(z):
        e, _, _, r, _ = quantities(z[:n])
        service = model['service']+model['service_jac']@z[:n]
        return e[-1]/1e14+10*np.mean(r*r)+2*np.mean(service**2)+20*z[-2]+20*z[-1]
    def constraints(z):
        e, power, accel, r, g = quantities(z[:n])
        # Clearance slack is km; terminal slack is dimensionless (100 km/5 m/s).
        tolerance = np.tile(np.array([50./np.sqrt(3)/1e5]*3+[.01/np.sqrt(3)/5]*3), len(r)//6)
        return np.r_[(g-150.)/1000.+z[-2], tolerance+z[-1]-abs(r),
                     z[-1]-(model['service']+model['service_jac']@z[:n]),
                     1-power/(capacity_factor*ex.capacity), 1-accel/.001,
                     (51.859e12*prefix_factor-e[0])/1e14,
                     (100e12*energy_factor-e[1])/1e14]
    start = np.r_[np.zeros(n), max(0., (150-model['gaps'].min())/1000.),
                  max(abs(model['terminal']))]
    if warm_start is not None:
        start[:n] = np.asarray(warm_start['x'])[cols]-x0[cols]
        _, _, _, r, g = quantities(start[:n])
        start[-2:] = [max(0., (150-g.min())/1000.), max(abs(r))]
    bounds = []
    for j in cols:
        cap = 1. if j >= 16 else (5. if j >= 10 else 30.)
        bounds.append((max(-radius, -cap-x0[j]), min(radius, cap-x0[j])))
    fit = minimize(obj, start, method='SLSQP', bounds=bounds+[(0, None), (0, None)],
                   constraints=[dict(type='ineq', fun=constraints)],
                   options=dict(maxiter=80, ftol=1e-8))
    z = fit.x; e, power, acc, r, g = quantities(z[:n]); x = x0.copy(); x[cols] += z[:n]
    return dict(x=x, step=z[:n], columns=cols, success=bool(fit.success), message=fit.message,
                iterations=fit.nit, radius=radius, energy_factor=energy_factor,
                prefix_energy_factor=prefix_factor, capacity_factor=capacity_factor,
                minimum_scaled_constraint=float(constraints(z).min()),
                clearance_slack_m=z[-2]*1000, terminal_slack=z[-1],
                predicted_energy_J=e, predicted_terminal_rms=float(np.sqrt(np.mean(r*r))),
                predicted_minimum_gap_m=float(g.min()), predicted_max_power_ratio=float(np.max(power/ex.capacity)),
                predicted_minimum_service_coverage=float(1.-np.max(model['service']+model['service_jac']@z[:n])),
                predicted_max_acceleration_m_s2=float(acc.max()), objective=float(fit.fun),
                accepted_cycle=False, virtual_controls=0.)


def main():
    global RUN, check_budget
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage', choices=['prepare', 'benchmark', 'sensitivity', 'solve', 'replay', 'assess', 'validate', 'inspect'])
    ap.add_argument('--run-dir', type=Path, default=RUN)
    ap.add_argument('--cpu-budget', type=float, default=1800.)
    ap.add_argument('--reserve', type=float, default=400.)
    ap.add_argument('--columns', type=int, nargs='+', default=[6, 7, 0, 10, 11, 12, 13, 14, 15, 16, 17, 8, 9, 1, 2, 3, 4, 5])
    ap.add_argument('--h', type=float, default=.2)
    ap.add_argument('--max-step', type=float, default=480.)
    ap.add_argument('--radius', type=float, default=1.)
    ap.add_argument('--case', default='actual')
    ap.add_argument('--energy-factor', type=float, default=1.)
    ap.add_argument('--capacity-factor', type=float, default=1.)
    ap.add_argument('--prefix-factor', type=float, default=1.)
    ap.add_argument('--warm-start')
    ap.add_argument('--candidate', type=Path)
    ap.add_argument('--resume-model', type=Path)
    ap.add_argument('--fine', action='store_true')
    ap.add_argument('--stop', action='store_true')
    ap.add_argument('--end-hour', type=float)
    args = ap.parse_args(); RUN = args.run_dir; RUN.mkdir(parents=True, exist_ok=True)
    for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
        if os.environ.get(name) != '1': raise RuntimeError('Set '+name+'=1 before starting the numerical runner')
    limit(); cpu = time.process_time(); child = resource.getrusage(resource.RUSAGE_CHILDREN)
    ledger_path = RUN/'ledger.json'
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else dict(charged_cpu_s=0., stages=[])
    if not ledger['stages'] and (RUN/'benchmark.json').exists():
        prior = json.loads((RUN/'benchmark.json').read_text())
        ledger['charged_cpu_s'] = prior['total_cpu_s']+10.  # compilation and failed dependency probe allowance
    remaining = args.cpu_budget-ledger['charged_cpu_s']; deadline = cpu+remaining
    def explicit_guard():
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        consumed = time.process_time()+children.ru_utime+children.ru_stime-child.ru_utime-child.ru_stime
        if consumed >= deadline-5.:
            raise TimeoutError('Explicit cumulative CPU limit reached; completed checkpoints retained')
    check_budget = explicit_guard
    from . import expanded_model
    expanded_model.check_budget = explicit_guard
    def expired(signum, frame): raise TimeoutError('Numerical CPU budget reached; completed products retained')
    signal.signal(signal.SIGPROF, expired)
    if remaining <= 0: raise RuntimeError('Budget exhausted; set an explicitly larger local budget to continue')
    signal.setitimer(signal.ITIMER_PROF, max(.1, remaining-2.))
    record = dict(stage=args.stage, args=vars(args).copy(), completed=False)
    record['args']['run_dir'] = str(RUN)
    record['args']['candidate'] = str(args.candidate) if args.candidate else None
    record['args']['resume_model'] = str(args.resume_model) if args.resume_model else None
    try:
        check_budget()
        if args.stage == 'prepare': prepare()
        elif args.stage == 'inspect':
            ex = ExpandedExperiment(); print('Loaded', len(ex.mass), 'tiles;', DIM, 'controls; mass', ex.mass.sum())
            if (PACKAGE/'restart.json').exists():
                state = json.loads((PACKAGE/'restart.json').read_text())
                for field in ['source_hashes', 'input_hashes']:
                    for name, expected in state[field].items():
                        if digest(ROOT/name) != expected: raise ValueError('Restart identity changed: '+name)
                for name, expected in state['payload_hashes'].items():
                    if digest(PACKAGE/name) != expected: raise ValueError('Restart payload changed: '+name)
                model = raw_read(PACKAGE/'optimizer_model.npz')
                print('Checked optimizer model', model['terminal_jac'].shape, 'and contacts', model['cuts'].shape)
        else:
            ex = ExpandedExperiment()
            if args.candidate:
                write(RUN/(args.case+'_proposal.json'), json.loads(args.candidate.read_text()))
            if args.stage == 'benchmark': benchmark(ex, args.max_step)
            elif args.stage == 'sensitivity': sensitivity(ex, args.columns, args.h, args.max_step, deadline, args.reserve)
            elif args.stage == 'solve':
                model = raw_read(args.resume_model) if args.resume_model else build_model(ex)
                write(RUN/'rank.json', rank_report(model))
                warm = json.loads((RUN/(args.warm_start+'_proposal.json')).read_text()) if args.warm_start else None
                fit = solve(ex, model, args.radius, args.energy_factor, args.capacity_factor, args.prefix_factor, warm)
                write(RUN/(args.case+'_proposal.json'), fit)
                print(json.dumps(native(fit)), flush=True)
            elif args.stage in ('assess', 'validate'):
                from .expanded_validation import assess_trial, validate
                (assess_trial if args.stage == 'assess' else validate)(ex, RUN, args.case)
            elif args.stage == 'replay':
                proposal = json.loads((RUN/(args.case+'_proposal.json')).read_text())
                x = np.array(proposal['x'])
                end = args.end_hour*3600 if args.end_hour is not None else None
                row, raw, _ = ex.evaluate(x, fine=args.fine, stop=args.stop, end=end,
                                          max_step=args.max_step, output_step=60. if args.fine else 300.)
                if row['completed_horizon'] and end is None:
                    row['cycle'] = ex.summary(raw, x); row['energy_parts_J'] = energy_parts(ex, raw, x)
                save_case(args.case+('_fine' if args.fine else '_replay'), row, raw)
                print(json.dumps(native(row)), flush=True)
        record['completed'] = True
    except BaseException as exc:
        record['error'] = repr(exc); raise
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0.)
        expanded_model.check_budget = lambda: None
        check_budget = lambda: None
        endchild = resource.getrusage(resource.RUSAGE_CHILDREN)
        charge = time.process_time()-cpu+endchild.ru_utime+endchild.ru_stime-child.ru_utime-child.ru_stime
        record.update(cpu_s=charge, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        ledger['charged_cpu_s'] += charge; ledger['stages'].append(record)
        ledger['raw_MiB'] = sum(p.stat().st_size for p in RUN.rglob('*') if p.is_file())/1024**2
        write(ledger_path, ledger)


if __name__ == '__main__': main()
