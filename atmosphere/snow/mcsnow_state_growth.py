"""Test size-dependent growth allocation in a separate native McSnow build.

These are diagnostic closures inspired by observed morphology transitions, not
an implementation of a published facet-spreading or molecular attachment model.
The neutral limit preserves the six-control tuner and its pinned baseline.
"""
from __future__ import annotations

import argparse
import copy
import difflib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import differential_evolution, least_squares

from atmosphere.snow import mcsnow_tuning as tuning

diagnosis = tuning.diagnosis
native, model = tuning.native, tuning.model
HERE = Path(__file__).resolve().parent
CACHE = tuning.trace.CACHE/'state_growth'
PRODUCT = HERE/'results/mcsnow_state_growth.json'
MAX_RADIUS_UM = 150.
BASE = np.r_[tuning.BASE, 0.]
FAMILIES = ('axis', 'axis_density')


def stage_weight(radius_m, transition_radius_um):
    """Bounded diagnostic transition; the fourth power is a declared assumption."""
    if not all(map(math.isfinite, (radius_m, transition_radius_um))) or radius_m <= 0 or transition_radius_um < 0:
        raise ValueError('Positive radius and nonnegative transition scale required')
    return 1./(1.+(transition_radius_um*1e-6/radius_m)**4)


def decode(x, family='axis', seed_diameter_um=0.):
    x = np.asarray(x, dtype=float)
    if x.shape != (7,) or not np.all(np.isfinite(x)) or np.any((x < 0) | (x > 1)):
        raise ValueError('Seven finite coordinates in [0,1] required')
    if family not in FAMILIES or not math.isfinite(seed_diameter_um) or seed_diameter_um < 0:
        raise ValueError('Unknown family or invalid seed diameter')
    return np.r_[tuning.physical(x[:6]), MAX_RADIUS_UM*x[6], seed_diameter_um,
                 int(family == 'axis_density')]


def deposition_source():
    s = tuning.deposition_source()
    replace = model.replace_once
    s = replace(s, '  REAL(wp), PUBLIC :: tune(6)=1._wp',
        '  REAL(wp), PUBLIC :: tune(6)=1._wp\n  REAL(wp), PUBLIC :: stage(3)=0._wp')
    s = replace(s, '  LOGICAL  :: branch_crit', '  REAL(wp) :: stage_w\n  LOGICAL  :: branch_crit')
    s = replace(s, '  growth_trace(13)=IGF', '''  growth_trace(13)=IGF
  stage_w=1._wp
  if (stage(1)>0._wp) stage_w=1._wp/(1._wp+(stage(1)*1e-6_wp/(z3_4pi*V_tot)**z1_3)**4)''')
    s = replace(s, '  growth_trace(15)=rho_depo', '''  ! Diagnostic: compact deposition at small size, approaching the tuned density.
  if (stage(1)>0._wp.and.stage(3)>0.5_wp) rho_depo=rhoi*(rho_depo/rhoi)**stage_w
  if (sp%d<stage(2)*1e-6_wp) rho_depo=rhoi
  growth_trace(15)=rho_depo''')
    s = replace(s, 'IGF*tune(2)*(IGF_vent/IGF)**tune(3)',
                   'IGF*tune(2)**stage_w*(IGF_vent/IGF)**tune(3)')
    return replace(s, '  growth_trace(14)=IGF_vent',
        '  if (sp%d<stage(2)*1e-6_wp) IGF_vent=1._wp\n  growth_trace(14)=IGF_vent')


def driver_source():
    s = tuning.DRIVER.read_text()
    s = model.replace_once(s, 'sp_depo, tune, tune_cap_count', 'sp_depo, tune, stage, tune_cap_count')
    s = model.replace_once(s, 'read(input_unit,*,iostat=ios) tune', 'read(input_unit,*,iostat=ios) tune,stage')
    return s.replace('any(.not.ieee_is_finite(tune))', 'any(.not.ieee_is_finite([tune,stage]))')


def build():
    path = CACHE/'build'
    baseline = native.verify_build(model.BASELINE)
    if path.exists(): return verify_build()
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source', path/'source')
    (path/'source/src/mo_depo.f90').write_text(deposition_source())
    driver = path/'mcsnow_state_driver.f90'
    driver.write_text(driver_source())
    for name in native.TABLE_NAMES: shutil.copyfile(model.BASELINE/name, path/name)
    (path/'state.patch').write_text(''.join(difflib.unified_diff(
        tuning.deposition_source().splitlines(True), deposition_source().splitlines(True),
        fromfile='constant/mo_depo.f90', tofile='state/mo_depo.f90')))
    compiler = baseline['compiler_path']
    selected = native.selected_sources(path/'source/src', driver, compiler)
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
    record = dict(schema='terluna.atmosphere.mcsnow-state-build/1',
        baseline_sha256=native.sha(model.BASELINE/'build.json'), source_identity=baseline['source_identity'],
        compiled_sources={str(p): native.sha(p) for p in selected},
        compiler_path=compiler, compiler_sha256=native.sha(compiler),
        compiler_version=baseline['compiler_version'], flags=native.FLAGS, commands=commands,
        files={p.name: native.sha(p) for p in [path/'tuner', path/'state.patch', path/'build.log',
                                            *[path/n for n in native.TABLE_NAMES]]})
    model.save(path/'build.json', record)
    return record


def verify_build():
    path = CACHE/'build'
    r = json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if r['baseline_sha256'] != native.sha(model.BASELINE/'build.json'): raise ValueError('Baseline changed')
    for p, sha in r['compiled_sources'].items():
        if native.sha(p) != sha: raise ValueError('Compiled source changed: '+p)
    for name, sha in r['files'].items():
        if native.sha(path/name) != sha: raise ValueError('Build file changed: '+name)
    if (path/'source/src/mo_depo.f90').read_text() != deposition_source(): raise ValueError('Source changed')
    if (path/'mcsnow_state_driver.f90').read_text() != driver_source(): raise ValueError('Driver changed')
    return r


def batch(entries, dt=tuning.SEARCH_DT, family='axis', seed_diameter_um=0.):
    return tuning.Batch(entries, dt, cache=CACHE,
        decoder=lambda x: decode(x, family, seed_diameter_um),
        experiment=dict(family=family, seed_diameter_um=seed_diameter_um,
                        radius_upper_bound_um=MAX_RADIUS_UM, transition_power=4))


def training_cases():
    return [e for e in tuning.cases() if e['temperature_C'] in tuning.CORE_TEMPERATURES
            and e['initialization'] != 'volume']


def previous_best():
    p = json.loads(tuning.PRODUCT.read_text())
    return next(f for f in p['searches'] if f['name'] == p['selected_search'])


def fit(family, seed, generations=35, temperature=None, seed_diameter_um=0.):
    suffix = '' if temperature is None else diagnosis.key(temperature)+'_'
    if seed_diameter_um: suffix += f'seed{seed_diameter_um:g}_'
    name = f'{family}_{suffix}{seed}'
    path = CACHE/'fits'/name; path.mkdir(parents=True, exist_ok=True)
    entries = training_cases()
    if temperature is not None: entries = [e for e in entries if e['temperature_C'] == temperature]
    initial = np.r_[previous_best()['coordinates'], 30./MAX_RADIUS_UM]
    if temperature is not None:
        common = [json.loads((CACHE/'fits'/f'{family}_{sd}'/'result.json').read_text())
                  for sd in (20261024,20261025)]
        initial = np.array(min(common,key=lambda r:r['training_loss'])['coordinates'])
    spec = dict(family=family, seed=seed, generations=generations, population_multiplier=6,
                initial=initial.tolist(), dt=tuning.SEARCH_DT, training_ages=list(tuning.TRAIN_AGES),
                cases=[e['name'] for e in entries], radius_upper_bound_um=MAX_RADIUS_UM,
                scipy_version=scipy.__version__, build_sha256=native.sha(CACHE/'build/build.json'))
    if seed_diameter_um: spec['seed_diameter_um'] = seed_diameter_um
    if (path/'specification.json').exists():
        if json.loads((path/'specification.json').read_text()) != spec: raise ValueError('Search changed')
    else: model.save(path/'specification.json', spec)
    if (path/'result.json').exists(): return json.loads((path/'result.json').read_text())
    with batch(entries, family=family, seed_diameter_um=seed_diameter_um) as b:
        start = time.monotonic(); last = start; generation = 0
        def fun(x): return tuning.residual(b, x)
        def objective(x): return float(np.sum(fun(x)**2))
        def callback(x, convergence):
            nonlocal last, generation
            generation += 1
            if time.monotonic()-last > 15:
                print(name, generation, math.sqrt(objective(x)), flush=True); last = time.monotonic()
        de = differential_evolution(objective, [(0.,1.)]*7, seed=seed, maxiter=generations,
            popsize=6, tol=1e-5, atol=1e-9, workers=1, updating='immediate', polish=False,
            x0=initial, callback=callback)
        local = least_squares(fun, de.x, bounds=(np.zeros(7),np.ones(7)), diff_step=.005,
            x_scale='jac', max_nfev=200, ftol=1e-8, xtol=1e-7, gtol=1e-7)
        old = np.r_[previous_best()['coordinates'],0.]
        if temperature is not None:
            prior = json.loads(tuning.PRODUCT.read_text())
            previous = next(r for r in prior['searches'] if r['name']=='temperature_'+diagnosis.key(temperature))
            old = np.r_[previous['coordinates'],0.]
        x = min([old, initial, de.x, local.x], key=objective)
        r = dict(name=name, family=family, coordinates=x.tolist(), physical_controls=decode(x,family,seed_diameter_um).tolist(),
            training_loss=objective(x), training_log_rms=math.sqrt(objective(x)), specification=spec,
            differential_evolution=dict(success=bool(de.success), message=str(de.message),
                                        evaluations=int(de.nfev), generations=int(de.nit)),
            local_refinement=dict(success=bool(local.success), message=str(local.message), evaluations=int(local.nfev)),
            elapsed_s=time.monotonic()-start, evaluation_directory=str(b.path.relative_to(model.ROOT)),
            population_coordinates=de.population.tolist(), population_losses=de.population_energies.tolist())
        model.save(path/'result.json', r)
        print(name, 'complete', r['training_log_rms'], r['physical_controls'], flush=True)
        return r


def scan():
    build(); entries = training_cases(); rows = []
    # Isolate a restored compact, isotropic seed stage before refitting anything.
    for seed_um in (0., 10., 20., 40., 100.):
        with batch(entries, seed_diameter_um=seed_um) as b:
            for name, x in [('original', BASE), ('shared', np.r_[previous_best()['coordinates'],0.])]:
                r = tuning.assessment(b, x)
                rows.append(dict(name=name, seed_diameter_um=seed_um, assessment=r))
                print('scan', name, seed_um, r.get('scores',{}).get('training'), flush=True)
    model.save(CACHE/'seed_scan.json', rows)
    return rows


def search():
    build()
    return [fit(family, seed) for family in FAMILIES for seed in (20261024, 20261025)]


def temperatures():
    build()
    fits = [json.loads((CACHE/'fits'/f'{family}_{seed}'/'result.json').read_text())
            for family in FAMILIES for seed in (20261024,20261025)]
    selected = min(fits,key=lambda r:r['training_loss'])['family']
    return [fit(selected,20261024,25,temperature=t) for t in tuning.CORE_TEMPERATURES]


def seeds():
    build()
    fits = [json.loads((CACHE/'fits'/f'{family}_{seed}'/'result.json').read_text())
            for family in FAMILIES for seed in (20261024,20261025)]
    selected=min(fits,key=lambda r:r['training_loss'])['family']
    return [fit(selected,20261024,25,temperature=-5.3,seed_diameter_um=diameter)
            for diameter in (20.,40.)]


def suite(entries, x, family, dt, seed_diameter_um=0.):
    rows = []; failed = []; ledger = volume = density = 0.
    for e in entries:
        with batch([e], dt, family, seed_diameter_um) as b: r = tuning.assessment(b, x)
        if 'error' in r:
            failed.append(dict(case=e['name'], error=r['error'])); continue
        rows.extend(r['rows'])
        ledger = max(ledger, r['mass_balance_max'])
        volume = max(volume, r['volume_residual_max'])
        density = max(density, r['maximum_density_fraction'])
    return dict(rows=rows, scores=tuning.score_rows(rows), failed_cases=failed,
                mass_balance_max=ledger, volume_residual_max=volume, maximum_density_fraction=density)


def comparison_difference(a, b):
    def key(r): return r['temperature_C'], r['initialization'], r['age_min']
    aa = {key(r): r for r in a['rows']}; bb = {key(r): r for r in b['rows']}
    if aa.keys() != bb.keys() or a.get('failed_cases') or b.get('failed_cases'):
        raise ValueError('Different or incomplete comparison states')
    return max(abs(x/y-1) for k in aa for x,y in zip(aa[k]['values'],bb[k]['values']))


def probes(x, family, dt=1/128, duration=1/32, temperature=None, seed_diameter_um=0.):
    """Short independent starts at supplied fitted geometries; no history fitting."""
    result = []; rho = diagnosis.ice_density()
    for t in (tuning.CORE_TEMPERATURES if temperature is None else (temperature,)):
        for age in (3.5, 5.5, 8., 30.):
            for mapping in diagnosis.MAPPINGS:
                target = diagnosis.fit_rates(t, age, mapping)
                if target['density_kg_m3'] > rho: continue
                e = copy.deepcopy(next(e for e in tuning.cases()
                                       if e['temperature_C'] == t and e['initialization'] == mapping))
                e['initial_age_min'] = age
                e['sample_ages_min'] = [age, age+duration/60]
                e['configuration'].update(initial_mass_kg=target['mass_kg'],
                    initial_equatorial_m=target['model_a_m'], initial_polar_m=target['model_c_m'])
                with batch([e], dt, family, seed_diameter_um) as b: r = b.evaluate(x)
                if r['error']: raise ValueError('Fixed-state probe failed: '+r['error'])
                start, end = np.asarray(r['rows'])[:,:3]
                dm = (end[0]-start[0])/duration
                axes = np.log(end[1:]/start[1:])/duration
                v = math.pi/6*np.array([start[1]**2*start[2], end[1]**2*end[2]])
                rates = [dm, *axes]
                reference = [target['mass_rate_kg_s'], target['a_log_rate_s'], target['c_log_rate_s']]
                result.append(dict(temperature_C=t, age_min=age, mapping=mapping,
                    rates=rates, target_rates=reference, ratio=(np.array(rates)/reference).tolist(),
                    deposition_density_kg_m3=(end[0]-start[0])/(v[1]-v[0]),
                    target_deposition_density_kg_m3=target['deposition_density_kg_m3'],
                    stage_weight=stage_weight((3*v[0]/(4*math.pi))**(1/3), MAX_RADIUS_UM*x[6])))
    return result


def seed_sensitivity(x, family, seed_diameter_um=0., temperature=None):
    rows = []
    for t in (tuning.CORE_TEMPERATURES if temperature is None else (temperature,)):
        original = next(e for e in tuning.cases() if e['temperature_C']==t and e['initialization']=='sphere')
        for phi in (.25,.5,1.,2.,4.):
            e = copy.deepcopy(original)
            diameter = e['configuration']['initial_equatorial_m']
            e['configuration'].update(initial_equatorial_m=diameter*phi**(-1/3),
                                      initial_polar_m=diameter*phi**(2/3))
            with batch([e],1/128,family,seed_diameter_um) as b: r=tuning.assessment(b,x)
            rows.append(dict(temperature_C=t,initial_aspect_ratio=phi,
                initial_volume_equivalent_diameter_m=diameter,assessment=r))
    return rows


def equivalent_axis_ratio(initial_phi, final_phi, volume_ratio):
    """Invert the exact integrated native shape relation, without a time law.

    This is an equivalent constant ratio, not the arithmetic average of a
    variable ratio and not a measured early attachment coefficient.
    """
    if not all(map(math.isfinite,(initial_phi,final_phi,volume_ratio))):
        raise ValueError('Finite geometry required')
    if min(initial_phi,final_phi)<=0 or volume_ratio<=1:
        raise ValueError('Positive aspect ratios and growing volume required')
    exponent=math.log(final_phi/initial_phi)/math.log(volume_ratio)
    if not -.5 < exponent < 1:
        raise ValueError('Geometry requires a nonpositive or unbounded growth ratio')
    return (1+2*exponent)/(1-exponent)


def shape_budget():
    rows=[]
    for t in tuning.CORE_TEMPERATURES:
        config=model.configuration(model.native.scale.lunar.flakes.T0+t)
        volume0=math.pi/6*config['initial_equatorial_m']**2*config['initial_polar_m']
        phi0=config['initial_polar_m']/config['initial_equatorial_m']
        for mapping in diagnosis.MAPPINGS:
            target=diagnosis.fit_state(t,3.,mapping)
            phi=target['model_c_m']/target['model_a_m']
            volume_ratio=target['volume_m3']/volume0
            gamma=equivalent_axis_ratio(phi0,phi,volume_ratio)
            mature=diagnosis.fit_rates(t,3.5,mapping)['effective_growth_ratio']
            required_seed=phi/volume_ratio**((mature-1)/(mature+2))
            rows.append(dict(temperature_C=t,mapping=mapping,initial_phi=phi0,
                initial_envelope_volume_m3=volume0,final_phi=phi,volume_ratio=volume_ratio,
                equivalent_axis_ratio_before_3min=gamma,
                fitted_axis_ratio_at_3p5min=mature,
                fitted_axis_ratio_at_30min=diagnosis.fit_rates(t,30.,mapping)['effective_growth_ratio'],
                initial_phi_required_if_later_ratio_applied_throughout=required_seed,
                reconstruction_relative_error=phi0*volume_ratio**((gamma-1)/(gamma+2))/phi-1))
    return dict(evidence='Kinematic consequence of the assumed seed and fitted three-minute state; early history is unobserved.',
        equation='log(phi_final/phi_initial) = integral ((Gamma-1)/(Gamma+2)) dlog(V)',
        rows=rows,
        restriction='The inferred seed aspect ratio assumes the same initial envelope volume, freely changing its dimensions. It is a degeneracy demonstration, not a measured or adopted seed.')


def lateral_reference():
    """Read only explicitly identified complete rows from the pinned archive.

    Later rows contain empty series columns: whitespace splitting would shift
    their meanings. They are deliberately excluded, not padded or interpolated.
    """
    import io
    import zipfile
    source=json.loads((HERE/'mcsnow_state_sources.json').read_text())
    p=model.ROOT/source['data_archive']['cache_path']
    if native.sha(p)!=source['data_archive']['sha256']: raise ValueError('Lateral archive changed')
    with zipfile.ZipFile(p) as z:
        def read(suffix):
            name=next(n for n in z.namelist() if n.endswith('/'+suffix))
            data=z.read(name)
            import hashlib
            if hashlib.sha256(data).hexdigest()!=source['data_archive']['members'][name]:
                raise ValueError('Reference member changed')
            return data.decode()
        first=[list(map(float,l.split())) for l in read('fig3.output').splitlines()[:4]]
        if len(first)!=4 or any(len(row)!=40 for row in first):
            raise ValueError('Ambiguous or incomplete observation columns')
        # Metadata columns 17/18: measured normal a diameter, single pyramidal area.
        initial_d,initial_area=first[0][16:18]
        final_d,final_area=first[-1][16:18]
        output=np.loadtxt(io.StringIO(read('fig5b.output')))
        if output.shape!=(1001,4): raise ValueError('Author model output shape changed')
        row=output[np.flatnonzero(output[:,2]/output[:,1]>=.99)[0]]
    return dict(archive_sha256=native.sha(p),
        observed=dict(evidence='Author digitization of Gonda and Yamazaki (1984), not new measurements.',
            conditions='-15 C, 1.5% ice supersaturation; initial diameter about 23.8 micrometres.',
            complete_rows_used=4, columns_one_based=[17,18],
            first_dimension_um=initial_d,last_dimension_um=final_d,
            first_single_pyramidal_area_um2=initial_area,last_single_pyramidal_area_um2=final_area,
            dimension_ratio=final_d/initial_d,remaining_area_fraction=final_area/initial_area),
        author_model=dict(evidence='Published model output, not an independent experiment or our model reproduction.',
            rows=len(output),lateral_over_normal_threshold=.99,
            dimension_over_initial_at_threshold=float(row[1]/output[0,1])),
        scope='Confirms an observable facet-area transition at modest outer-size change. No facet-area predictions exist in the tested McSnow closures, so this is evidence for a missing observable, not their validation.')


def laboratory_plot_segments():
    result={}
    for t in tuning.CORE_TEMPERATURES:
        result[str(t)]={}
        for quantity,curve in zip(diagnosis.FIELDS[:3],diagnosis.lab.FITS[t][2:]):
            lower=3.;segments=[]
            for upper,_,_ in curve:
                end=min(float(upper),30.)
                if end>lower:
                    age=np.linspace(lower,np.nextafter(end,lower),100)
                    segments.append(dict(age_min=age.tolist(),
                        values=[diagnosis.lab.observed(t,float(a))[quantity] for a in age]))
                lower=end
            result[str(t)][quantity]=segments
    return result


def analyze():
    verify_build(); all_cases = tuning.cases(); checks = []; assessments = {}
    old = json.loads(tuning.PRODUCT.read_text())
    fits = [json.loads((CACHE/'fits'/f'{family}_{seed}'/'result.json').read_text())
            for family in FAMILIES for seed in (20261024,20261025)]
    settings = dict(original=(BASE,'axis'), constant=(np.r_[previous_best()['coordinates'],0.],'axis'))
    for family in FAMILIES:
        best = min((f for f in fits if f['family'] == family), key=lambda f:f['training_loss'])
        settings[family] = (np.asarray(best['coordinates']), family)
    for name, (x, family) in settings.items():
        fine = suite(all_cases, x, family, 1/128)
        assessments[name] = fine
        if name in ('original','constant'):
            ref = old['assessments']['baseline' if name == 'original' else 'shared']
            diff = comparison_difference(fine,ref)
            checks.append(dict(name='neutral_parity_'+name, maximum_relative_difference=diff, passed=diff<1e-10))
        else:
            coarse = suite(all_cases, x, family, 1/32)
            diff = comparison_difference(coarse,fine)
            checks.append(dict(name='all_case_timestep_'+name, maximum_relative_difference=diff, passed=diff<.01))
        print('assessment',name,fine['scores'],flush=True)
    for f in fits:
        with batch(training_cases(), 1/128, f['family']) as b:
            fine = tuning.assessment(b,f['coordinates'])
            f['fine_training_loss'] = float(np.sum(tuning.residual(b,f['coordinates'])**2))
        with batch(training_cases(), tuning.SEARCH_DT, f['family']) as b:
            coarse = tuning.assessment(b,f['coordinates'])
        diff = comparison_difference(coarse,fine)
        checks.append(dict(name='search_timestep_'+f['name'], maximum_relative_difference=diff, passed=diff<.01))
    temperature_fits = []
    for path in sorted((CACHE/'fits').glob('*_m*_20261024/result.json')):
        f = json.loads(path.read_text())
        entries = [e for e in all_cases if e['name'] in f['specification']['cases']]
        t = entries[0]['temperature_C']
        crossed = [e for e in all_cases if e['temperature_C']==t]
        seed_limit=f['specification'].get('seed_diameter_um',0.)
        fine = suite(crossed,f['coordinates'],f['family'],1/128,seed_limit)
        coarse = suite(crossed,f['coordinates'],f['family'],1/32,seed_limit)
        diff = comparison_difference(coarse,fine)
        checks.append(dict(name='temperature_timestep_'+f['name'],maximum_relative_difference=diff,passed=diff<.01))
        f['temperature_C'] = t
        f['seed_diameter_um'] = seed_limit
        f['assessment'] = fine
        f['fine_training_loss'] = float(np.mean([math.log(v)**2 for r in fine['rows']
            if r['group']=='training' for v in r['ratio'][:3]]))
        f['fixed_state_probes'] = probes(f['coordinates'],f['family'],temperature=t,seed_diameter_um=seed_limit)
        f['previous_constant_loss'] = next(r['training_loss'] for r in old['searches']
                                           if r['name']=='temperature_'+diagnosis.key(t))
        f['seed_shape_sensitivity']=seed_sensitivity(f['coordinates'],f['family'],seed_limit,t)
        if seed_limit:
            with batch(entries, family=f['family']) as b:
                output = b.evaluate(f['coordinates'])
                loss = None if output['error'] else float(np.sum(tuning.residual(b,f['coordinates'])**2))
            f['compact_stage_ablation'] = dict(training_loss=loss,error=output['error'],
                interpretation='Remove the compact stage with all fitted controls held fixed; compare also with the separately reoptimized zero-cutoff fit.')
        temperature_fits.append(f)
    warm = min((f for f in temperature_fits if f['temperature_C']==-5.3),key=lambda f:f['training_loss'])
    warm_cases = [e for e in training_cases() if e['temperature_C']==-5.3]
    with batch(warm_cases, family=warm['family'], seed_diameter_um=warm['seed_diameter_um']) as b:
        low=np.array(warm['coordinates']);high=low.copy();low[4]=0.;high[4]=1.
        low_result=tuning.assessment(b,low);high_result=tuning.assessment(b,high)
    branch_diff=comparison_difference(low_result,high_result)
    checks.append(dict(name='warm_branch_radius_inactive',maximum_relative_difference=branch_diff,
                       passed=branch_diff<1e-12))
    warm['inactive_control_check']=dict(control='branch_radius',
        physical_multipliers=[.25,4.],maximum_relative_difference=branch_diff,
        reason='Native columns hollow immediately; the branching-size criterion applies only when the underlying IGF is below one.')
    ablations = []
    for name in FAMILIES:
        x,family = settings[name]
        reset = x.copy(); reset[6] = 0.
        with batch(training_cases(),family=family) as b:
            output=b.evaluate(reset)
            loss=None if output['error'] else float(np.sum(tuning.residual(b,reset)**2))
        ablations.append(dict(family=name,reset_transition_loss=loss,error=output['error'],
            previous_reoptimized_constant_loss=previous_best()['training_loss']))
    fixed = {name:probes(x,family) for name,(x,family) in settings.items()}
    best_name = min(FAMILIES, key=lambda name:min(f['training_loss'] for f in fits if f['family']==name))
    x,family = settings[best_name]
    ref = probes(x,family,1/256,1/64)
    diff = float(max(abs(a/b-1) for l,h in zip(fixed[best_name],ref) for a,b in zip(l['rates'],h['rates'])))
    checks.append(dict(name='fixed_state_duration_and_timestep', maximum_relative_difference=diff, passed=diff<.01))
    # All preceding held-out comparisons are assessments, not selection criteria.
    scans = json.loads((CACHE/'seed_scan.json').read_text())
    seeds = {name:seed_sensitivity(*settings[name]) for name in ('original','constant',best_name)}
    manifest = []; evaluations = failures = 0; ledger = volume = 0.
    for p in sorted((CACHE/'evaluations').glob('*/evaluations.jsonl')):
        records = [json.loads(l) for l in p.read_text().splitlines()]
        evaluations += len(records); failures += sum(r['error'] is not None for r in records)
        for r in records:
            for row in r['rows']: ledger=max(ledger,abs(row[4]));volume=max(volume,abs(row[5]))
        manifest.append(dict(directory=str(p.parent.relative_to(model.ROOT)), evaluations=len(records),
            files={f.name:native.sha(f) for f in sorted(p.parent.iterdir()) if f.is_file()}))
    checks.append(dict(name='mass_ledger',maximum_relative_difference=ledger,passed=ledger<1e-10))
    checks.append(dict(name='spheroid_volume',maximum_relative_difference=volume,passed=volume<1e-10))
    result = dict(schema='terluna.atmosphere.mcsnow-state-growth/1',
        evidence='Diagnostic shared state-dependent closure; no adopted physical or lunar law.',
        quantities=diagnosis.FIELDS,
        units=dict(age_min='min',mass_kg='kg',equatorial_dimension_m='m',polar_dimension_m='m',terminal_speed_m_s='m/s'),
        reading_rule='Rows are sampled single-crystal states; laboratory plotting segments retain separate fitted branches. No interpolation across temperatures or initializations is defined. Ratios compare to fitted populations, not independently measured uncertainties.',
        laboratory_plot_segments=laboratory_plot_segments(),
        selected_temperature_searches={str(t):min([f for f in temperature_fits if f['temperature_C']==t],key=lambda f:f['training_loss'])['name'] for t in tuning.CORE_TEMPERATURES if any(f['temperature_C']==t for f in temperature_fits)},
        transition=dict(radius='enclosing-volume equivalent sphere radius',
            equation='w=1/(1+(r_transition/r_equivalent)^4)', radius_interval_um=[0,MAX_RADIUS_UM],
            axis_multiplier='axis_control**w', density_coupled='rho_deposition/rho_ice raised to w',
            limits='0 transition radius recovers the constant model; these equations are diagnostic assumptions.'),
        objective=old['objective'], training=old['training'], previous_product_path=str(tuning.PRODUCT.relative_to(model.ROOT)), previous_product_sha256=native.sha(tuning.PRODUCT),
        selected_family=best_name, selection='Minimum coarse training loss; no reserved measurements in the objective. Candidates must execute through the full requested horizon.',
        searches=fits, temperature_searches=temperature_fits, transition_ablation=ablations,
        assessments=assessments, fixed_state_probes=fixed, seed_scan=scans,
        seed_shape_sensitivity=seeds, shape_budget=shape_budget(), lateral_reference=lateral_reference(),
        source_record_sha256=native.sha(HERE/'mcsnow_state_sources.json'),
        checks=checks,evaluations=evaluations,failed_evaluations=failures,
        mass_balance_max=ledger,volume_residual_max=volume,manifest=manifest,
        producer={str(p.relative_to(model.ROOT)):native.sha(p) for p in [Path(__file__),Path(tuning.__file__),tuning.DRIVER]},
        fit_artifacts={str(p.relative_to(model.ROOT)):native.sha(p) for p in sorted((CACHE/'fits').glob('*/*.json'))},
        build_record_sha256=native.sha(CACHE/'build/build.json'),
        inputs={str(p.relative_to(model.ROOT)):native.sha(p) for p in [Path(diagnosis.__file__),Path(diagnosis.lab.__file__),Path(model.__file__)]},
        restrictions=['Fitted populations and unknown initial crystal shapes remain limitations.',
                      'Withheld ages share the fitting source; they are not independent experimental validation.',
                      'Reserved measurements do not enter the score, but a candidate must execute through the full requested horizon.',
                      'No reproduction of the 2021 lateral-growth model is claimed.',
                      'The original author executable/configuration remains unidentified.'])
    model.save(PRODUCT,result)
    print('checks',checks,flush=True)
    return result


def main():
    p = argparse.ArgumentParser(); p.add_argument('action', choices=['build','scan','search','temperatures','seeds','analyze'])
    a = p.parse_args(); globals()[a.action]()


if __name__ == '__main__': main()
