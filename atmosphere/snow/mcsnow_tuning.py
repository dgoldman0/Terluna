"""Bounded diagnostic calibration of native McSnow growth; no adopted law.

Reuse the pinned model, preserve every evaluated vector, and keep validation
cases out of the fitting objective. Bounds are exploratory, not physical priors.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import differential_evolution, least_squares

from atmosphere.snow import mcsnow_shape_diagnosis as diagnosis

trace = diagnosis.trace
native = trace.native
model = trace.model
HERE = Path(__file__).resolve().parent
CACHE = trace.CACHE / 'tuning'
PRODUCT = HERE / 'results/mcsnow_tuning.json'
DRIVER = HERE / 'mcsnow_tune_driver.f90'
CORE_TEMPERATURES = (-14.4, -5.3)
TRAIN_AGES = (5., 10., 20.)
SEARCH_DT = 1 / 16
# Multipliers use logarithmic coordinates; exponent/enhancement use linear ones.
KNOBS = (
    ('transport', .25, 4., 'log', 'Common multiplier of capacitance in vapor and heat transfer'),
    ('axis_ratio', .25, 4., 'log', 'Multiplier of effective polar/equatorial logarithmic growth ratio'),
    ('directional', 0., 2., 'linear', 'Exponent of native directional ventilation ratio; zero removes it'),
    ('density', .25, 4., 'log', 'Multiplier of porous deposition density, bounded above by solid ice'),
    ('branch_radius', .25, 4., 'log', 'Multiplier of native 200 micrometre branching radius'),
    ('ventilation', 0., 2., 'linear', 'Multiplier of mean ventilation enhancement above unity'),
)


def physical(x):
    x = np.asarray(x, dtype=float)
    if x.shape != (len(KNOBS),) or not np.all(np.isfinite(x)) or np.any((x < 0) | (x > 1)):
        raise ValueError('Controls must be six finite coordinates in [0,1]')
    return np.array([lo*(hi/lo)**v if kind == 'log' else lo+(hi-lo)*v
                     for v, (_, lo, hi, kind, _) in zip(x, KNOBS)])


def coordinates(values):
    values = np.asarray(values, dtype=float)
    if values.shape != (len(KNOBS),) or not np.all(np.isfinite(values)):
        raise ValueError('Six finite physical controls required')
    x = np.array([math.log(v/lo)/math.log(hi/lo) if kind == 'log' else (v-lo)/(hi-lo)
                  for v, (_, lo, hi, kind, _) in zip(values, KNOBS)])
    physical(x)
    return x


BASE = coordinates(np.ones(6))


def deposition_source():
    body = trace.source_patch((model.BASELINE/'source/src/mo_depo.f90').read_text(), 'plate_base')
    replace = model.replace_once
    body = replace(body, '  REAL(wp), PUBLIC :: growth_trace(22)=0._wp',
        '  REAL(wp), PUBLIC :: growth_trace(22)=0._wp\n' +
        '  REAL(wp), PUBLIC :: tune(6)=1._wp\n  INTEGER, PUBLIC :: tune_cap_count=0')
    body = replace(body, '  C = capacitance( sp, Re )\n  growth_trace=0._wp',
                   '  C = capacitance( sp, Re ) * tune(1)\n  growth_trace=0._wp')
    start = body.index('SUBROUTINE vapourDiffusion( sp,')
    end = body.index('END SUBROUTINE vapourDiffusion', start)
    part = body[start:end]
    for var, arg in (('Dv', 'Sc'), ('Ka', 'Pr')):
        old = f'    {var} = {var}*ventilation_coef_i( sqrtRe*cubr{arg}, sp%phi )'
        part = replace(part, old, f'    {var} = {var}*(1._wp+tune(6)*(ventilation_coef_i( sqrtRe*cubr{arg}, sp%phi )-1._wp))')
    body = body[:start]+part+body[end:]
    body = replace(body, '(0.5_wp * sp%d) .GE. 200e-6_wp',
                   '(0.5_wp * sp%d) .GE. (200e-6_wp*tune(5))')
    body = replace(body, '  growth_trace(15)=rho_depo',
        '''  if (rho_depo<rhoi) then
    if (rho_depo*tune(4)>rhoi) tune_cap_count=tune_cap_count+1
    rho_depo=min(rhoi,rho_depo*tune(4))
  end if
  growth_trace(15)=rho_depo''')
    return replace(body, '  growth_trace(14)=IGF_vent',
        '  IGF_vent=IGF*tune(2)*(IGF_vent/IGF)**tune(3)\n  growth_trace(14)=IGF_vent')


def build():
    path = CACHE/'build'
    baseline = native.verify_build(model.BASELINE)
    if path.exists():
        return verify_build()
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source', path/'source')
    original = (model.BASELINE/'source/src/mo_depo.f90').read_text()
    changed = deposition_source()
    (path/'source/src/mo_depo.f90').write_text(changed)
    shutil.copyfile(DRIVER, path/DRIVER.name)
    for name in native.TABLE_NAMES:
        shutil.copyfile(model.BASELINE/name, path/name)
    (path/'tuning.patch').write_text(''.join(difflib.unified_diff(
        trace.source_patch(original, 'plate_base').splitlines(True), changed.splitlines(True),
        fromfile='improved/mo_depo.f90', tofile='diagnostic/mo_depo.f90')))
    compiler = baseline['compiler_path']
    selected = native.selected_sources(path/'source/src', path/DRIVER.name, compiler)
    commands, objects = [], []
    with (path/'build.log').open('w') as log:
        for i, src in enumerate(selected):
            obj = path/f'{i:02d}_{src.stem}.o'; objects.append(str(obj))
            cmd = [compiler, *native.FLAGS, '-I'+str(path/'source/src'), '-I'+str(path),
                   '-J'+str(path), '-c', str(src), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *native.FLAGS, *objects, '-o', str(path/'tuner')]
        commands.append(cmd)
        subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.mcsnow-tuning-build/1',
        baseline_sha256=native.sha(model.BASELINE/'build.json'), source_identity=baseline['source_identity'],
        compiled_sources={str(p): native.sha(p) for p in selected},
        compiler_path=compiler, compiler_sha256=native.sha(compiler),
        compiler_version=baseline['compiler_version'], flags=native.FLAGS, commands=commands,
        files={p.name: native.sha(p) for p in [path/'tuner', path/'tuning.patch', path/'build.log',
                                            *[path/n for n in native.TABLE_NAMES]]})
    model.save(path/'build.json', record)
    return record


def verify_build():
    path = CACHE/'build'
    record = json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if record['baseline_sha256'] != native.sha(model.BASELINE/'build.json'):
        raise ValueError('Baseline changed')
    for p, sha in record['compiled_sources'].items():
        if native.sha(p) != sha: raise ValueError('Compiled source changed: '+p)
    for name, sha in record['files'].items():
        if native.sha(path/name) != sha: raise ValueError('Build file changed: '+name)
    if (path/'source/src/mo_depo.f90').read_text() != deposition_source() or (path/DRIVER.name).read_bytes() != DRIVER.read_bytes():
        raise ValueError('Build differs from current declared source')
    return record


def cases():
    # Reuse the accepted experiment definitions and exact admissibility decisions.
    # The native air record identifies the solid density; this reads no additional simulation.
    rho = diagnosis.ice_density()
    return [e for e in diagnosis.cases(rho) if e['admissible'] and e['dt_s'] == min(diagnosis.TIMESTEPS)]


class Batch:
    """One native worker, persistent only within this context; deterministic cache."""
    def __init__(self, entries, dt=SEARCH_DT, *, cache=CACHE, decoder=physical, experiment=None):
        if decoder is not physical and experiment is None:
            raise ValueError('A custom decoder requires an explicit experiment identity')
        self.root = cache; self.decode = decoder
        self.entries = entries; self.dt = dt; self.proc = None; self.stderr = None
        self.rows = [(i, age) for i, e in enumerate(entries) for age in e.get('sample_ages_min', diagnosis.sample_ages(e['temperature_C']))
                     if age >= e['initial_age_min']]
        config = dict(entries=entries, dt=dt, build_sha256=native.sha(self.root/'build/build.json'))
        if experiment is not None: config['experiment'] = experiment
        encoded = json.dumps(config, sort_keys=True).encode()
        self.path = self.root/'evaluations'/hashlib.sha256(encoded).hexdigest()[:20]
        self.path.mkdir(parents=True, exist_ok=True)
        if (self.path/'configuration.json').exists():
            if json.loads((self.path/'configuration.json').read_text()) != config:
                raise ValueError('Batch configuration changed')
        else: model.save(self.path/'configuration.json', config)
        lines = [str(len(entries))]
        for e in entries:
            c = e['configuration']
            ages = [a for a in e.get('sample_ages_min', diagnosis.sample_ages(e['temperature_C'])) if a >= e['initial_age_min']]
            values = [c['temperature_K'], c['pressure_Pa'], 1., c['initial_mass_kg'],
                      c['initial_equatorial_m'], c['initial_polar_m'], dt, len(ages)]
            lines.extend([' '.join(map(str, values)), ' '.join(str((a-e['initial_age_min'])*60) for a in ages)])
        data = '\n'.join(lines)+'\n'
        p = self.path/'cases.txt'
        if p.exists() and p.read_text() != data: raise ValueError('Batch input changed')
        p.write_text(data)
        self.cache = {}
        journal = self.path/'evaluations.jsonl'
        if journal.exists():
            for line in journal.read_text().splitlines():
                r = json.loads(line); self.cache[tuple(r['coordinates'])] = r

    def start(self):
        self.stderr = (self.path/'stderr.log').open('a')
        self.proc = subprocess.Popen([str(self.root/'build/tuner'), str(self.path/'cases.txt')],
            cwd=self.root/'build', env=native.ENV, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=self.stderr, text=True, bufsize=1)
        while True:
            line = self.proc.stdout.readline()
            if not line: raise RuntimeError('Native worker failed to initialize')
            if line.strip() == 'READY': break

    def close(self):
        if self.proc:
            if self.proc.poll() is None:
                self.proc.stdin.close(); self.proc.wait(timeout=20)
            self.proc.stdout.close()
            self.proc = None
        if self.stderr: self.stderr.close(); self.stderr = None

    def __enter__(self): return self
    def __exit__(self, *_): self.close()

    def evaluate(self, x):
        values = self.decode(x); key = tuple(float(v) for v in x)
        if key in self.cache: return self.cache[key]
        started = time.monotonic(); output = []; error = None
        try:
            if self.proc is None: self.start()
            self.proc.stdin.write(' '.join(format(v, '.17g') for v in values)+'\n'); self.proc.stdin.flush()
            for _ in self.rows:
                line = self.proc.stdout.readline()
                if not line: raise RuntimeError('Native worker terminated during trajectory')
                row = list(map(float, line.split()))
                if len(row) != 8 or not np.all(np.isfinite(row)) or min(row[:3]) <= 0 or row[3] < 0:
                    raise RuntimeError('Invalid native output')
                if abs(row[4]) > 1e-10 or abs(row[5]) > 1e-10 or row[7] > 1+1e-10:
                    raise RuntimeError('Mass/volume/density invariant failed')
                output.append(row)
            if self.proc.stdout.readline().strip() != 'DONE': raise RuntimeError('Incomplete native batch')
        except (RuntimeError, ValueError, BrokenPipeError) as exc:
            error = str(exc)
            if self.proc and self.proc.poll() is None: self.proc.terminate()
            self.close()
        result = dict(coordinates=list(key), physical_values=values.tolist(), rows=output,
                      error=error, elapsed_s=time.monotonic()-started)
        with (self.path/'evaluations.jsonl').open('a') as stream:
            stream.write(json.dumps(result, allow_nan=False)+'\n'); stream.flush()
        self.cache[key] = result
        return result


def targets(batch):
    rows = []
    for i, age in batch.rows:
        e = batch.entries[i]
        obs = diagnosis.lab.observed(e['temperature_C'], age)
        rows.append([obs[k] for k in diagnosis.FIELDS[:3]] + [diagnosis.lab.fall_speed(e['temperature_C'], age)])
    return np.array(rows)


def reported(batch, result):
    if result['error']: raise ValueError(result['error'])
    values = np.array(result['rows'])[:, :4].copy()
    for j, (i, _) in enumerate(batch.rows): values[j, 1:3] /= batch.entries[i]['axis_scale']
    return values


def residual(batch, x, ages=TRAIN_AGES):
    selected = [j for j, (i, age) in enumerate(batch.rows) if age in ages and age > batch.entries[i]['initial_age_min']]
    counts = {i: sum(batch.rows[j][0] == i for j in selected) for i in range(len(batch.entries))}
    if not selected or any(n == 0 for n in counts.values()): raise ValueError('Every fit case needs targets')
    result = batch.evaluate(x)
    if result['error']: return np.full(3*len(selected), 100.)
    ratios = reported(batch, result)/targets(batch)
    return np.concatenate([np.log(ratios[j, :3])/math.sqrt(3*len(counts)*counts[batch.rows[j][0]]) for j in selected])


def optimize(name, entries, seed, generations, initial=BASE, fixed=None):
    """Global exploration followed by bounded local refinement; record both exits."""
    fixed = fixed or {}
    spec = dict(name=name, seed=seed, generations=generations, population_multiplier=6,
                initial=list(map(float, initial)), fixed={str(k): float(v) for k, v in fixed.items()},
                cases=[e['name'] for e in entries], training_ages_min=TRAIN_AGES,
                search_dt_s=SEARCH_DT, knobs=KNOBS, scipy_version=scipy.__version__,
                build_sha256=native.sha(CACHE/'build/build.json'))
    # Normalize tuples before exact comparisons with saved JSON.
    spec = json.loads(json.dumps(spec))
    path = CACHE/'fits'/name
    path.mkdir(parents=True, exist_ok=True)
    if (path/'specification.json').exists():
        if json.loads((path/'specification.json').read_text()) != spec:
            raise ValueError('Search specification changed: '+name)
    else: model.save(path/'specification.json', spec)
    if (path/'result.json').exists(): return json.loads((path/'result.json').read_text())
    free = [i for i in range(6) if i not in fixed]
    def expand(z):
        x = np.array(initial, dtype=float); x[free] = z
        for i, value in fixed.items(): x[i] = value
        return x
    with Batch(entries) as batch:
        started = time.monotonic(); last_report = started; generation = 0
        def fun(z): return residual(batch, expand(z))
        def objective(z): return float(np.sum(fun(z)**2))
        def callback(x, convergence):
            nonlocal last_report, generation
            generation += 1
            if time.monotonic()-last_report > 15:
                print(name, 'generation', generation, 'log-RMS', math.sqrt(objective(x)), flush=True)
                last_report = time.monotonic()
        de = differential_evolution(objective, [(0., 1.)]*len(free), seed=seed,
            maxiter=generations, popsize=6, tol=1e-5, atol=1e-9, workers=1,
            updating='immediate', polish=False, x0=np.asarray(initial)[free], callback=callback)
        local = least_squares(fun, de.x, bounds=(np.zeros(len(free)), np.ones(len(free))),
            diff_step=.005, x_scale='jac', max_nfev=160, ftol=1e-8, xtol=1e-7, gtol=1e-7)
        candidates = [np.asarray(initial)[free], de.x, local.x]
        best = min(candidates, key=objective); x = expand(best)
        result = dict(name=name, coordinates=x.tolist(), parameters=dict(zip([k[0] for k in KNOBS], physical(x))),
            training_loss=objective(best), training_log_rms=math.sqrt(objective(best)),
            baseline_loss=float(np.sum(residual(batch, BASE)**2)), specification=spec,
            differential_evolution=dict(success=bool(de.success), message=str(de.message),
                                        evaluations=int(de.nfev), generations=int(de.nit)),
            local_refinement=dict(success=bool(local.success), message=str(local.message), evaluations=int(local.nfev)),
            active_bounds=[KNOBS[i][0] for i in free if min(x[i], 1-x[i]) < .005],
            elapsed_s=time.monotonic()-started, evaluation_directory=str(batch.path.relative_to(model.ROOT)),
            population_coordinates=[expand(z).tolist() for z in de.population],
            population_losses=de.population_energies.tolist())
        model.save(path/'result.json', result)
        print(name, 'finished', result['training_log_rms'], result['parameters'], flush=True)
        return result


def search():
    build(); all_cases = cases()
    training = [e for e in all_cases if e['temperature_C'] in CORE_TEMPERATURES and e['initialization'] != 'volume']
    fits = [optimize('shared_'+str(seed), training, seed, 35) for seed in (20261023, 20261024)]
    best = min(fits, key=lambda r: r['training_loss'])
    for t in CORE_TEMPERATURES:
        subset = [e for e in training if e['temperature_C'] == t]
        fits.append(optimize('temperature_'+diagnosis.key(t), subset, 20261023, 25, initial=best['coordinates']))
        for init in ('sphere', 'axes'):
            fits.append(optimize('history_'+diagnosis.key(t)+'_'+init,
                                 [e for e in subset if e['initialization'] == init],
                                 20261023, 20, initial=best['coordinates']))
    for i, knob in enumerate(KNOBS):
        fits.append(optimize('lock_'+knob[0], training, 20261023, 25,
                             initial=best['coordinates'], fixed={i: BASE[i]}))
    model.save(CACHE/'searches.json', fits)
    return fits


def score_rows(rows):
    scores = {}
    for group in ('training', 'withheld_age', 'geometry_sensitivity', 'other_temperature'):
        rr = np.array([r['ratio'] for r in rows if r['group'] == group])
        if len(rr):
            scores[group] = dict(states=len(rr), log_rms_by_quantity=np.sqrt(np.mean(np.log(rr)**2,axis=0)).tolist(),
                median_absolute_percent_by_quantity=(100*np.median(abs(rr-1),axis=0)).tolist(),
                maximum_absolute_percent_by_quantity=(100*np.max(abs(rr-1),axis=0)).tolist())
    return scores


def assessment(batch, x):
    result = batch.evaluate(x)
    if result['error']: return dict(error=result['error'])
    values = reported(batch, result); obs = targets(batch); ratios = values/obs
    rows = []
    for j, (i, age) in enumerate(batch.rows):
        e = batch.entries[i]; temp=e['temperature_C']; init=e['initialization']
        group = ('other_temperature' if temp not in CORE_TEMPERATURES else
                 'geometry_sensitivity' if init == 'volume' else
                 'training' if age in TRAIN_AGES else 'withheld_age')
        # A fitted restart's initial state is prescribed, not a prediction.
        if age == e['initial_age_min']: group = 'prescribed_initial_state'
        rows.append(dict(temperature_C=temp, initialization=init, age_min=age, group=group,
                         values=values[j].tolist(), target=obs[j].tolist(), ratio=ratios[j].tolist()))
    return dict(rows=rows, scores=score_rows(rows), mass_balance_max=max(abs(r[4]) for r in result['rows']),
                volume_residual_max=max(abs(r[5]) for r in result['rows']),
                density_cap_max_step_count=max(r[6] for r in result['rows']),
                maximum_density_fraction=max(r[7] for r in result['rows']))


def suite(entries, x, dt):
    """Isolate held-out failures so that one case cannot hide the others."""
    rows=[]; failed=[]; ledger=volume=cap=density=0.
    for entry in entries:
        with Batch([entry],dt) as batch: result=assessment(batch,x)
        if 'error' in result:
            failed.append(dict(temperature_C=entry['temperature_C'],initialization=entry['initialization'],
                               error=result['error']))
            continue
        rows.extend(result['rows'])
        ledger=max(ledger,result['mass_balance_max']);volume=max(volume,result['volume_residual_max'])
        cap=max(cap,result['density_cap_max_step_count']);density=max(density,result['maximum_density_fraction'])
    return dict(rows=rows,scores=score_rows(rows),failed_cases=failed,
                mass_balance_max=ledger,volume_residual_max=volume,
                density_cap_max_step_count=cap,maximum_density_fraction=density)


def analyze():
    verify_build(); fits=json.loads((CACHE/'searches.json').read_text()); all_cases=cases()
    # Fixed-control refits are also feasible points in the common six-control space.
    # Select on the same training score only; no held-out prediction enters selection.
    best=min([r for r in fits if r['name'].startswith(('shared_', 'lock_'))], key=lambda r:r['training_loss'])
    training=[e for e in all_cases if e['temperature_C'] in CORE_TEMPERATURES and e['initialization']!='volume']
    checks=[]; assessments={}; reset=[]
    # Compare the new adapter with both preceding, independently executed controls.
    old=json.loads(diagnosis.PRODUCT.read_text())['trajectories']
    previous=json.loads((HERE/'results/mcsnow_shape_controls.json').read_text())['trajectories']
    for name,x in [('baseline',BASE),('shared',best['coordinates'])]:
        assessments[name]=suite(all_cases,x,1/128)
    maxdiff=0.
    for r in assessments['baseline']['rows']:
        target=next(e for e in old if e['temperature_C']==r['temperature_C'] and
                    e['initialization']==r['initialization'] and e['dt_s']==1/128)
        found=next(v for v in target['comparisons'] if v['age_min']==r['age_min'])
        ref=[found['age_matched']['values'][k] for k in diagnosis.FIELDS]
        maxdiff=max(maxdiff,float(np.max(abs(np.array(r['values'])/ref-1))))
    checks.append(dict(name='unchanged_adapter_parity',maximum_relative_difference=maxdiff,
                       passed=maxdiff<1e-10 and not assessments['baseline']['failed_cases']))
    with Batch(training,1/128) as batch:
        neutral=BASE.copy();neutral[2]=0.
        result=assessment(batch,neutral); maxdiff=0.
        for r in result['rows']:
            target=next(e for e in previous if e['temperature_C']==r['temperature_C'] and
                        e['initialization']==r['initialization'] and e['dt_s']==1/128 and e['shape_control']==1)
            found=next(v for v in target['comparisons'] if v['age_min']==r['age_min'])
            ref=[found['age_matched']['values'][k] for k in diagnosis.FIELDS]
            maxdiff=max(maxdiff,float(np.max(abs(np.array(r['values'])/ref-1))))
        checks.append(dict(name='unit_directional_control_parity',maximum_relative_difference=maxdiff,passed=maxdiff<1e-10))
    for fit in fits:
        selected=[e for e in all_cases if e['name'] in fit['specification']['cases']]
        with Batch(selected,1/128) as batch:
            fit['fine_assessment']=assessment(batch,fit['coordinates'])
            fit['fine_training_loss']=float(np.sum(residual(batch,fit['coordinates'])**2))
        with Batch(selected,SEARCH_DT) as batch:
            coarse=assessment(batch,fit['coordinates'])
        change=max(abs(a/b-1) for lo,hi in zip(coarse['rows'],fit['fine_assessment']['rows'])
                   for a,b in zip(lo['values'],hi['values']))
        checks.append(dict(name='timestep_'+fit['name'],maximum_relative_difference=change,passed=change<.01))
        if fit['name'].startswith(('temperature_', 'history_')):
            temp=selected[0]['temperature_C']
            crossed=[e for e in all_cases if e['temperature_C']==temp]
            cross=suite(crossed,fit['coordinates'],1/128)
            fit['cross_initialization_assessment']=cross
            if fit['name'].startswith('history_'):
                coarse_cross=suite(crossed,fit['coordinates'],1/64)
                failed=cross['failed_cases'] or coarse_cross['failed_cases']
                cross_change=None if failed else max(abs(a/b-1)
                    for lo,hi in zip(coarse_cross['rows'],cross['rows'])
                    for a,b in zip(lo['values'],hi['values']))
                checks.append(dict(name='cross_history_timestep_'+fit['name'],
                    maximum_relative_difference=cross_change,passed=not failed and cross_change<.01))
    with Batch(training) as batch:
        for i, knob in enumerate(KNOBS):
            x=np.array(best['coordinates']); x[i]=BASE[i]
            reset.append(dict(parameter=knob[0],reset_loss=float(np.sum(residual(batch,x)**2)),
                              refitted_loss=next(r['training_loss'] for r in fits if r['name']=='lock_'+knob[0])))
        # Finite-difference sensitivity in explicitly normalized search coordinates.
        jac=[]
        for i in range(6):
            lo=np.array(best['coordinates']);hi=lo.copy()
            lo[i]=max(0.,lo[i]-.005);hi[i]=min(1.,hi[i]+.005)
            jac.append((residual(batch,hi)-residual(batch,lo))/(hi[i]-lo[i]))
        _,singular,vh=np.linalg.svd(np.array(jac).T,full_matrices=False)
        sensitivity=dict(coordinates='Normalized search coordinates; finite difference width up to 0.01',
                         singular_values=singular.tolist(),weakest_combination=vh[-1].tolist(),
                         interpretation='Local sensitivity, not statistical confidence; branching can be nonsmooth.')
    manifests=[]; evaluations=0; failed=0; maxledger=0.; maxvolume=0.
    for p in sorted((CACHE/'evaluations').glob('*/evaluations.jsonl')):
        records=[json.loads(line) for line in p.read_text().splitlines()]
        evaluations+=len(records);failed+=sum(r['error'] is not None for r in records)
        for r in records:
            for row in r['rows']: maxledger=max(maxledger,abs(row[4]));maxvolume=max(maxvolume,abs(row[5]))
        manifests.append(dict(directory=str(p.parent.relative_to(model.ROOT)),evaluations=len(records),
                              files={f.name:native.sha(f) for f in sorted(p.parent.iterdir()) if f.is_file()}))
    return dict(schema='terluna.atmosphere.mcsnow-tuning/1',
        evidence='Diagnostic calibration to fitted Earth histories; best found parameters are not physical constants or an adopted lunar law.',
        quantities=diagnosis.FIELDS,
        raw_output_columns=[*diagnosis.FIELDS,'relative_mass_ledger_residual','relative_spheroid_volume_residual',
                            'solid_density_cap_step_count','mean_density_over_solid_density'],
        knobs=[dict(name=n,lower=lo,upper=hi,coordinate=k,meaning=m) for n,lo,hi,k,m in KNOBS],
        baseline_parameters=dict(zip([k[0] for k in KNOBS],np.ones(6))),selected_search=best['name'],
        selection_pool='Both shared searches and all neutral-control refits; same training cases and ages, no withheld data.',
        objective='Mean squared log mass/width/thickness ratios; equal weights per case, age and quantity. Speed is checked only.',
        training=dict(temperatures_C=CORE_TEMPERATURES,initializations=['sphere','axes'],ages_min=TRAIN_AGES),
        checks=checks,assessments=assessments,searches=fits,ablation=reset,sensitivity=sensitivity,
        evaluations=evaluations,failed_evaluations=failed,mass_balance_max=maxledger,volume_residual_max=maxvolume,
        producer={str(p.relative_to(model.ROOT)):native.sha(p) for p in
                  (Path(__file__),DRIVER,Path(diagnosis.__file__),Path(trace.__file__),Path(diagnosis.lab.__file__))},
        inputs={str(p.relative_to(model.ROOT)):native.sha(p) for p in
                (diagnosis.PRODUCT,HERE/'results/mcsnow_shape_controls.json')},
        laboratory_source_sha256=diagnosis.lab.SOURCE_SHA,
        build_record_sha256=native.sha(CACHE/'build/build.json'),manifest=manifests,
        search_artifacts={str(p.relative_to(model.ROOT)):native.sha(p) for p in
                          sorted((CACHE/'fits').glob('*/*.json'))},
        restrictions=['Bounds are exploratory search limits, not empirical uncertainty intervals.',
                      'Published fits summarize different crystals; planar thickness is partly inferred from mass and area.',
                      'Withheld ages share the fitted source; alternative envelope mappings are representation sensitivities.',
                      'Neither native drag nor intrinsic temperature tables are fitted; no humidity law is inferred.',
                      'Optimizer termination or parameter boundaries do not prove a global optimum or a unique mechanism.'])


def smoke():
    build()
    entries = [e for e in cases() if e['temperature_C'] in CORE_TEMPERATURES and e['initialization'] != 'volume']
    with Batch(entries) as batch:
        r = batch.evaluate(BASE)
        print(json.dumps(dict(elapsed_s=r['elapsed_s'], error=r['error'], rows=len(r['rows']),
                              objective=float(np.sum(residual(batch, BASE)**2))), indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('smoke','search','analyze'))
    args = p.parse_args()
    if args.action == 'smoke': smoke()
    elif args.action == 'search': search()
    else:
        result=analyze();model.save(PRODUCT,result)
        print(result['evaluations'],'parameter evaluations;',result['failed_evaluations'],'failed')


if __name__ == '__main__': main()
