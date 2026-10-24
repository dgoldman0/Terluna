"""Independent early and later axis-growth limits in a separate McSnow build.

The preceding radius interpolation fixed its early multiplier to one. This
experiment releases that boundary without changing the source baseline or
adopting the fitted multipliers as physical attachment coefficients.
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

from atmosphere.snow import mcsnow_state_growth as state

tuning, diagnosis, native, model = state.tuning, state.diagnosis, state.native, state.model
HERE = Path(__file__).resolve().parent
CACHE = tuning.trace.CACHE/'growth_limits'
PRODUCT = HERE/'results/mcsnow_growth_limits.json'
EARLY_BOUNDS = (.25, 4.)


def early_coordinate(value,bounds=EARLY_BOUNDS):
    if len(bounds)!=2 or not all(map(math.isfinite,bounds)) or not 0<bounds[0]<bounds[1]:
        raise ValueError('Positive increasing bounds required')
    if not math.isfinite(value) or not bounds[0] <= value <= bounds[1]:
        raise ValueError('Early multiplier outside declared bounds')
    return math.log(value/bounds[0])/math.log(bounds[1]/bounds[0])


def decode(x, family='axis', seed_diameter_um=0., early_bounds=EARLY_BOUNDS, ventilation_max=2.):
    x=np.asarray(x,dtype=float)
    if x.shape!=(8,) or not np.all(np.isfinite(x)) or np.any((x<0)|(x>1)):
        raise ValueError('Eight finite coordinates in [0,1] required')
    early_coordinate(1.,early_bounds)
    if not math.isfinite(ventilation_max) or ventilation_max<=0:raise ValueError('Positive ventilation bound required')
    values=state.decode(x[:7],family,seed_diameter_um);values[5]=ventilation_max*x[5]
    return np.r_[values,early_bounds[0]*(early_bounds[1]/early_bounds[0])**x[7]]


def axis_multiplier(radius_m, transition_radius_um, early, late):
    if not all(map(math.isfinite,(early,late))) or min(early,late)<=0:
        raise ValueError('Positive finite axis multipliers required')
    w=state.stage_weight(radius_m,transition_radius_um)
    return early**(1-w)*late**w


def deposition_source():
    s=model.replace_once(state.deposition_source(),
        '  REAL(wp), PUBLIC :: stage(3)=0._wp',
        '  REAL(wp), PUBLIC :: stage(4)=[0._wp,0._wp,0._wp,1._wp]')
    return model.replace_once(s,'IGF*tune(2)**stage_w*(IGF_vent/IGF)**tune(3)',
        'IGF*stage(4)**(1._wp-stage_w)*tune(2)**stage_w*(IGF_vent/IGF)**tune(3)')


def build():
    state.verify_build();path=CACHE/'build'
    if path.exists(): return verify_build()
    baseline=native.verify_build(model.BASELINE);path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source',path/'source')
    (path/'source/src/mo_depo.f90').write_text(deposition_source())
    driver=path/'mcsnow_limits_driver.f90';driver.write_text(state.driver_source())
    for name in native.TABLE_NAMES: shutil.copyfile(model.BASELINE/name,path/name)
    patch=path/'limits.patch'
    patch.write_text(''.join(difflib.unified_diff(state.deposition_source().splitlines(True),
        deposition_source().splitlines(True),fromfile='state/mo_depo.f90',tofile='limits/mo_depo.f90')))
    compiler=baseline['compiler_path'];selected=native.selected_sources(path/'source/src',driver,compiler)
    commands=[];objects=[]
    with (path/'build.log').open('w') as log:
        for i,source in enumerate(selected):
            obj=path/f'{i:02d}_{source.stem}.o';objects.append(str(obj))
            cmd=[compiler,*native.FLAGS,'-I'+str(path/'source/src'),'-I'+str(path),
                 '-J'+str(path),'-c',str(source),'-o',str(obj)]
            commands.append(cmd)
            subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler,*native.FLAGS,*objects,'-o',str(path/'tuner')];commands.append(cmd)
        subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    record=dict(schema='terluna.atmosphere.mcsnow-limits-build/1',
        baseline_sha256=native.sha(model.BASELINE/'build.json'),source_identity=baseline['source_identity'],
        parent_build_sha256=native.sha(state.CACHE/'build/build.json'),
        compiled_sources={str(p):native.sha(p) for p in selected},
        compiler_path=compiler,compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],
        flags=native.FLAGS,commands=commands,
        files={p.name:native.sha(p) for p in [path/'tuner',patch,path/'build.log',*[path/n for n in native.TABLE_NAMES]]})
    model.save(path/'build.json',record);return record


def verify_build():
    state.verify_build();path=CACHE/'build';r=json.loads((path/'build.json').read_text())
    if r['baseline_sha256']!=native.sha(model.BASELINE/'build.json'): raise ValueError('Baseline changed')
    if r['parent_build_sha256']!=native.sha(state.CACHE/'build/build.json'): raise ValueError('Parent changed')
    for p,sha in r['compiled_sources'].items():
        if native.sha(p)!=sha: raise ValueError('Compiled source changed: '+p)
    for name,sha in r['files'].items():
        if native.sha(path/name)!=sha: raise ValueError('Build file changed: '+name)
    if (path/'source/src/mo_depo.f90').read_text()!=deposition_source(): raise ValueError('Source changed')
    if (path/'mcsnow_limits_driver.f90').read_text()!=state.driver_source(): raise ValueError('Driver changed')
    return r


def batch(entries,dt=tuning.SEARCH_DT,family='axis',seed_diameter_um=0.,early_bounds=EARLY_BOUNDS,ventilation_max=2.):
    experiment=dict(family=family,seed_diameter_um=seed_diameter_um,
        radius_upper_bound_um=state.MAX_RADIUS_UM,transition_power=4,
        early_axis_multiplier_bounds=list(early_bounds),axis_rule='early**(1-w)*late**w')
    if ventilation_max!=2.:experiment['ventilation_max']=ventilation_max
    return tuning.Batch(entries,dt,cache=CACHE,
        decoder=lambda x:decode(x,family,seed_diameter_um,early_bounds,ventilation_max),experiment=experiment)


def fit_options(fit):
    return dict(family=fit['family'],seed_diameter_um=fit['seed_diameter_um'],
        early_bounds=tuple(fit['specification']['early_axis_multiplier_bounds']),
        ventilation_max=fit['specification'].get('ventilation_max',2.))


def group_key(fit):
    options=fit_options(fit)
    return (fit['temperature_C'],*options.values(),tuple(sorted(fit['specification'].get('fixed',{}).items())))


def cases(temperature):
    return [e for e in state.training_cases() if e['temperature_C']==temperature]


def prior_fit(temperature,seed_diameter_um=0.):
    old=json.loads(state.PRODUCT.read_text())
    return min((f for f in old['temperature_searches'] if f['temperature_C']==temperature
                and f['seed_diameter_um']==seed_diameter_um),key=lambda f:f['training_loss'])


def fit(temperature,seed,family='axis',seed_diameter_um=0.,generations=35,early_bounds=EARLY_BOUNDS,ventilation_max=2.,fixed=None,label=None):
    fixed=fixed or {}
    if fixed and not label:raise ValueError('Fixed-control fits require a descriptive label')
    if any(i not in range(8) or not math.isfinite(v) or not 0<=v<=1 for i,v in fixed.items()):
        raise ValueError('Invalid fixed coordinates')
    build();tag=f'{family}_{diagnosis.key(temperature)}_seed{seed_diameter_um:g}_{seed}'
    extended=tuple(early_bounds)!=EARLY_BOUNDS or ventilation_max!=2.
    if extended:tag+=f'_early{early_bounds[0]:g}to{early_bounds[1]:g}_vent{ventilation_max:g}'
    if label:tag+='_'+label
    path=CACHE/'fits'/tag;path.mkdir(parents=True,exist_ok=True)
    previous=prior_fit(temperature,seed_diameter_um)
    fallback=np.r_[previous['coordinates'],early_coordinate(1.,early_bounds)]
    fallback[5]*=2/ventilation_max
    initial=fallback.copy()
    if temperature==-14.4:
        initial[1]=early_coordinate(1.9);initial[6]=35/state.MAX_RADIUS_UM;initial[7]=early_coordinate(.52,early_bounds)
    elif seed_diameter_um:
        initial[1]=early_coordinate(2.1);initial[6]=50/state.MAX_RADIUS_UM;initial[7]=early_coordinate(3.5,early_bounds)
    else:
        initial[1]=early_coordinate(2.5);initial[6]=30/state.MAX_RADIUS_UM;initial[7]=early_coordinate(.66,early_bounds)
    if extended or fixed:
        previous_trials=[json.loads(p.read_text()) for p in (CACHE/'fits').glob('*/result.json')]
        base=[f for f in previous_trials if f['temperature_C']==temperature and f['family']==family
              and f['seed_diameter_um']==seed_diameter_um and fit_options(f)['early_bounds']==EARLY_BOUNDS
              and fit_options(f)['ventilation_max']==2. and not f['specification'].get('fixed')]
        if not base:raise ValueError('Run the default-bound comparison first')
        best=min(base,key=lambda f:f['training_loss']);initial=np.array(best['coordinates'])
        initial[5]=best['physical_controls'][5]/ventilation_max
        initial[7]=early_coordinate(best['physical_controls'][-1],early_bounds)
    for i,value in fixed.items():initial[i]=value;fallback[i]=value
    free=[i for i in range(8) if i not in fixed]
    if not free:raise ValueError('At least one free coordinate required')
    def expand(z):
        x=initial.copy();x[free]=z;return x
    entries=cases(temperature)
    spec=dict(family=family,temperature_C=temperature,seed_diameter_um=seed_diameter_um,
        seed=seed,generations=generations,population_multiplier=6,initial=initial.tolist(),
        fallback=fallback.tolist(),early_axis_multiplier_bounds=list(early_bounds),
        cases=[e['name'] for e in entries],training_ages=list(tuning.TRAIN_AGES),dt=tuning.SEARCH_DT,
        build_sha256=native.sha(CACHE/'build/build.json'),previous_product_sha256=native.sha(state.PRODUCT),
        scipy_version=scipy.__version__)
    if ventilation_max!=2.:spec['ventilation_max']=ventilation_max
    if fixed:spec['fixed']={str(i):v for i,v in sorted(fixed.items())}
    if (path/'specification.json').exists():
        if json.loads((path/'specification.json').read_text())!=spec: raise ValueError('Search specification changed')
    else:model.save(path/'specification.json',spec)
    if (path/'result.json').exists():return json.loads((path/'result.json').read_text())
    with batch(entries,family=family,seed_diameter_um=seed_diameter_um,early_bounds=early_bounds,ventilation_max=ventilation_max) as b:
        started=time.monotonic();last=started;generation=0
        def fun(x):return tuning.residual(b,x)
        def objective(x):return float(np.sum(fun(x)**2))
        def callback(x,convergence):
            nonlocal last,generation
            generation+=1
            if time.monotonic()-last>15:
                print(tag,generation,math.sqrt(objective(expand(x))),flush=True);last=time.monotonic()
        de=differential_evolution(lambda z:objective(expand(z)),[(0.,1.)]*len(free),seed=seed,maxiter=generations,popsize=6,
            tol=1e-5,atol=1e-9,workers=1,updating='immediate',polish=False,x0=initial[free],callback=callback)
        local=least_squares(lambda z:fun(expand(z)),de.x,bounds=(np.zeros(len(free)),np.ones(len(free))),diff_step=.005,x_scale='jac',
            max_nfev=240,ftol=1e-8,xtol=1e-7,gtol=1e-7)
        x=min([fallback,initial,expand(de.x),expand(local.x)],key=objective)
        result=dict(name=tag,temperature_C=temperature,family=family,seed_diameter_um=seed_diameter_um,
            coordinates=x.tolist(),physical_controls=decode(x,family,seed_diameter_um,early_bounds,ventilation_max).tolist(),
            training_loss=objective(x),training_log_rms=math.sqrt(objective(x)),baseline_loss=objective(fallback),
            previous_training_loss=previous['training_loss'],specification=spec,
            differential_evolution=dict(success=bool(de.success),message=str(de.message),
                evaluations=int(de.nfev),generations=int(de.nit)),
            local_refinement=dict(success=bool(local.success),message=str(local.message),evaluations=int(local.nfev)),
            elapsed_s=time.monotonic()-started,evaluation_directory=str(b.path.relative_to(model.ROOT)),
            population_coordinates=[expand(z).tolist() for z in de.population],population_losses=de.population_energies.tolist())
        model.save(path/'result.json',result)
        print(tag,'complete',result['training_log_rms'],result['physical_controls'],flush=True)
        return result



def refine_from_candidate(target, source):
    """Test a better nested-family solution before interpreting a freed control."""
    fixed={int(k):v for k,v in target['specification'].get('fixed',{}).items()}
    if fit_options(target)!=fit_options(source) or target['temperature_C']!=source['temperature_C']:
        raise ValueError('Nested refit requires identical experiment conditions')
    initial=np.asarray(source['coordinates'])
    if any(initial[i]!=v for i,v in fixed.items()):
        raise ValueError('Source violates target fixed controls')
    free=[i for i in range(8) if i not in fixed]
    def expand(z):
        x=initial.copy();x[free]=z;return x
    name=target['name']+'_nested'
    path=CACHE/'fits'/name;path.mkdir(parents=True,exist_ok=True)
    spec=copy.deepcopy(target['specification'])
    spec.update(initial=initial.tolist(),fallback=initial.tolist(),generations=0,
        refinement='Local refinement from the better nested control family',
        refined_from=source['name'],refined_from_sha256=native.sha(CACHE/'fits'/source['name']/'result.json'))
    if (path/'specification.json').exists():
        if json.loads((path/'specification.json').read_text())!=spec:
            raise ValueError('Nested refinement specification changed')
    else:model.save(path/'specification.json',spec)
    if (path/'result.json').exists():return json.loads((path/'result.json').read_text())
    with batch(cases(target['temperature_C']),**fit_options(target)) as b:
        started=time.monotonic()
        def fun(z):return tuning.residual(b,expand(z))
        def objective(x):return float(np.sum(tuning.residual(b,x)**2))
        local=least_squares(fun,initial[free],bounds=(np.zeros(len(free)),np.ones(len(free))),
            diff_step=.005,x_scale='jac',max_nfev=240,ftol=1e-8,xtol=1e-7,gtol=1e-7)
        x=min([initial,expand(local.x)],key=objective)
        result=dict(name=name,temperature_C=target['temperature_C'],family=target['family'],
            seed_diameter_um=target['seed_diameter_um'],coordinates=x.tolist(),
            physical_controls=decode(x,**fit_options(target)).tolist(),
            training_loss=objective(x),training_log_rms=math.sqrt(objective(x)),baseline_loss=objective(initial),
            previous_training_loss=target['previous_training_loss'],specification=spec,
            differential_evolution=dict(skipped=True,evaluations=0,generations=0,
                message='This named followup refines an already evaluated nested-family candidate.'),
            local_refinement=dict(success=bool(local.success),message=str(local.message),evaluations=int(local.nfev)),
            elapsed_s=time.monotonic()-started,evaluation_directory=str(b.path.relative_to(model.ROOT)))
        model.save(path/'result.json',result)
        print(name,'complete',result['training_log_rms'],result['physical_controls'],flush=True)
        return result


def nested():
    """Preserve finite-search results and add explicit nested-family restarts."""
    fits=[json.loads(p.read_text()) for p in sorted((CACHE/'fits').glob('*/result.json'))]
    original=[f for f in fits if not f['specification'].get('refined_from')]
    for target in original:
        fixed=target['specification'].get('fixed',{})
        eligible=[f for f in original if f['temperature_C']==target['temperature_C']
            and fit_options(f)==fit_options(target)
            and len(f['specification'].get('fixed',{}))>len(fixed)
            and all(f['specification'].get('fixed',{}).get(k)==v for k,v in fixed.items())]
        if eligible:
            source=min(eligible,key=lambda f:f['training_loss'])
            if source['training_loss']<target['training_loss']:
                refine_from_candidate(target,source)


def suite(entries,x,dt=1/128,family='axis',seed_diameter_um=0.,early_bounds=EARLY_BOUNDS,ventilation_max=2.):
    rows=[];failed=[];ledger=volume=density=0.
    for e in entries:
        with batch([e],dt,family,seed_diameter_um,early_bounds,ventilation_max) as b:r=tuning.assessment(b,x)
        if 'error' in r:failed.append(dict(case=e['name'],error=r['error']));continue
        rows.extend(r['rows']);ledger=max(ledger,r['mass_balance_max'])
        volume=max(volume,r['volume_residual_max']);density=max(density,r['maximum_density_fraction'])
    return dict(rows=rows,scores=tuning.score_rows(rows),failed_cases=failed,
        mass_balance_max=ledger,volume_residual_max=volume,maximum_density_fraction=density)


def probes(fit,dt=1/128,duration=1/32):
    t=fit['temperature_C'];result=[];rho=diagnosis.ice_density()
    for age in (3.5,5.5,8.,30.):
        for mapping in diagnosis.MAPPINGS:
            target=diagnosis.fit_rates(t,age,mapping)
            if target['density_kg_m3']>rho:continue
            e=copy.deepcopy(next(e for e in tuning.cases() if e['temperature_C']==t and e['initialization']==mapping))
            e['initial_age_min']=age;e['sample_ages_min']=[age,age+duration/60]
            e['configuration'].update(initial_mass_kg=target['mass_kg'],
                initial_equatorial_m=target['model_a_m'],initial_polar_m=target['model_c_m'])
            with batch([e],dt,**fit_options(fit)) as b:r=b.evaluate(fit['coordinates'])
            if r['error']:raise ValueError('Fixed-state probe failed: '+r['error'])
            start,end=np.asarray(r['rows'])[:,:3]
            rates=np.r_[(end[0]-start[0])/duration,np.log(end[1:]/start[1:])/duration]
            reference=np.array([target['mass_rate_kg_s'],target['a_log_rate_s'],target['c_log_rate_s']])
            volumes=math.pi/6*np.array([start[1]**2*start[2],end[1]**2*end[2]])
            radius=(3*volumes[0]/(4*math.pi))**(1/3)
            controls=fit['physical_controls']
            result.append(dict(temperature_C=t,age_min=age,mapping=mapping,rates=rates.tolist(),
                target_rates=reference.tolist(),ratio=(rates/reference).tolist(),
                effective_axis_ratio=float(rates[2]/rates[1]),
                target_effective_axis_ratio=float(reference[2]/reference[1]),
                equivalent_radius_m=radius,
                axis_multiplier=axis_multiplier(radius,controls[6],controls[9],controls[1]),
                deposition_density_kg_m3=float((end[0]-start[0])/(volumes[1]-volumes[0])),
                target_deposition_density_kg_m3=target['deposition_density_kg_m3']))
    return result


def seed_sensitivity(fit):
    result=[];t=fit['temperature_C']
    original=next(e for e in cases(t) if e['initialization']=='sphere')
    for phi in (.25,.5,1.,2.,4.):
        e=copy.deepcopy(original);diameter=e['configuration']['initial_equatorial_m']
        e['configuration'].update(initial_equatorial_m=diameter*phi**(-1/3),initial_polar_m=diameter*phi**(2/3))
        with batch([e],1/128,**fit_options(fit)) as b:r=tuning.assessment(b,fit['coordinates'])
        result.append(dict(initial_aspect_ratio=phi,initial_volume_equivalent_diameter_m=diameter,assessment=r))
    return result


def analyze():
    verify_build();old=json.loads(state.PRODUCT.read_text());checks=[]
    fits=[json.loads(p.read_text()) for p in sorted((CACHE/'fits').glob('*/result.json'))]
    if not fits:raise ValueError('No completed searches')
    parity={}
    for name,previous in [('original',None),('constant',state.previous_best())]:
        x=np.r_[tuning.BASE if previous is None else previous['coordinates'],0.,early_coordinate(1.)]
        assessment=suite(tuning.cases(),x);diff=state.comparison_difference(assessment,old['assessments'][name])
        checks.append(dict(name='neutral_parity_'+name,maximum_relative_difference=diff,passed=diff<1e-10))
        parity[name]=dict(states=len(assessment['rows']),maximum_relative_difference=diff)
    for family in state.FAMILIES:
        previous=min((f for f in old['searches'] if f['family']==family),key=lambda f:f['training_loss'])
        assessment=suite(tuning.cases(),np.r_[previous['coordinates'],early_coordinate(1.)],family=family)
        diff=state.comparison_difference(assessment,old['assessments'][family])
        checks.append(dict(name='parent_state_parity_'+family,maximum_relative_difference=diff,passed=diff<1e-10))
        parity[family]=dict(states=len(assessment['rows']),maximum_relative_difference=diff)
    groups=sorted({group_key(f) for f in fits})
    selected=[]
    for group in groups:
        candidates=[f for f in fits if group_key(f)==group]
        best=min(candidates,key=lambda f:f['training_loss']);selected.append(best['name'])
        for f in candidates:
            with batch(cases(f['temperature_C']),1/128,**fit_options(f)) as b:
                f['fine_training_loss']=float(np.sum(tuning.residual(b,f['coordinates'])**2))
                fine=tuning.assessment(b,f['coordinates'])
            with batch(cases(f['temperature_C']),tuning.SEARCH_DT,**fit_options(f)) as b:
                coarse=tuning.assessment(b,f['coordinates'])
            diff=state.comparison_difference(coarse,fine)
            checks.append(dict(name='search_timestep_'+f['name'],maximum_relative_difference=diff,passed=diff<.01))
        entries=[e for e in tuning.cases() if e['temperature_C']==best['temperature_C']]
        best['assessment']=suite(entries,best['coordinates'],1/128,**fit_options(best))
        coarse=suite(entries,best['coordinates'],1/32,**fit_options(best))
        diff=state.comparison_difference(coarse,best['assessment'])
        checks.append(dict(name='all_start_timestep_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
        best['fixed_state_probes']=probes(best);ref=probes(best,1/256,1/64)
        diff=max(abs(a/b-1) for lo,hi in zip(best['fixed_state_probes'],ref) for a,b in zip(lo['rates'],hi['rates']))
        checks.append(dict(name='fixed_state_refinement_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
        reset=np.array(best['coordinates']);reset[7]=early_coordinate(1.,fit_options(best)['early_bounds'])
        with batch(cases(best['temperature_C']),**fit_options(best)) as b:
            r=b.evaluate(reset)
            best['early_limit_ablation']=dict(error=r['error'],training_loss=None if r['error'] else float(np.sum(tuning.residual(b,reset)**2)))
        if best['temperature_C']==-5.3 and not best['specification'].get('fixed'):
            with batch(cases(-5.3),1/128,**fit_options(best)) as b:
                reference=tuning.assessment(b,best['coordinates']);branch_checks=[]
                for value in (0.,1.):
                    x=np.array(best['coordinates']);x[4]=value
                    trial=tuning.assessment(b,x)
                    diff=state.comparison_difference(reference,trial)
                    branch_checks.append(dict(branch_multiplier=decode(x,**fit_options(best))[4],
                        maximum_relative_difference=diff))
            best['branch_radius_controls']=branch_checks
            diff=max(r['maximum_relative_difference'] for r in branch_checks)
            checks.append(dict(name='inactive_column_branch_'+best['name'],maximum_relative_difference=diff,passed=diff<1e-12))
            if best['seed_diameter_um']:
                options=fit_options(best);options['seed_diameter_um']=0.
                with batch(cases(-5.3),**options) as b:
                    r=b.evaluate(best['coordinates'])
                    best['compact_stage_ablation']=dict(error=r['error'],
                        training_loss=None if r['error'] else float(np.sum(tuning.residual(b,best['coordinates'])**2)))
        best['seed_shape_sensitivity']=seed_sensitivity(best)
    manifest=[];evaluations=failures=0;ledger=volume=0.
    for p in sorted((CACHE/'evaluations').glob('*/evaluations.jsonl')):
        rows=[json.loads(l) for l in p.read_text().splitlines()]
        evaluations+=len(rows);failures+=sum(r['error'] is not None for r in rows)
        for r in rows:
            for v in r['rows']:ledger=max(ledger,abs(v[4]));volume=max(volume,abs(v[5]))
        manifest.append(dict(directory=str(p.parent.relative_to(model.ROOT)),evaluations=len(rows),
            files={f.name:native.sha(f) for f in sorted(p.parent.iterdir()) if f.is_file()}))
    checks.extend([dict(name='mass_ledger',maximum_relative_difference=ledger,passed=ledger<1e-10),
                   dict(name='spheroid_volume',maximum_relative_difference=volume,passed=volume<1e-10)])
    result=dict(schema='terluna.atmosphere.mcsnow-growth-limits/1',
        evidence='Diagnostic early and later growth multipliers; no adopted physical or lunar growth law.',
        quantities=old['quantities'],units=old['units'],reading_rule=old['reading_rule'],
        equation='q(r)=q_early**(1-w)*q_late**w, w=1/(1+(R/r_equivalent)**4)',
        physical_control_order=[k[0] for k in tuning.KNOBS]+[
            'transition_radius_um','compact_stage_maximum_diameter_um','density_transition_flag','early_axis_ratio'],
        physical_control_meanings=dict(tuning_controls=[dict(name=k[0],lower=k[1],upper=k[2],
            coordinate=k[3],meaning=k[4]) for k in tuning.KNOBS],
            later_axis_ratio='The axis_ratio control is the later multiplier; early_axis_ratio is its independent early limit.',
            bounds='Per-search specifications supersede these default bounds.',
            transition='Enclosing-volume equivalent radius, fourth-power transition, upper search bound 150 micrometres.'),
        default_early_bounds=list(EARLY_BOUNDS),bounds_rule='Each search specification records its actual bounds; extended ranges are diagnostic controls.',previous_product_path=str(state.PRODUCT.relative_to(model.ROOT)),
        previous_product_sha256=native.sha(state.PRODUCT),
        selection='Minimum training loss within each declared temperature/family/compact-stage/bounds/fixed-control group; reserved ages excluded from the score but executions must reach them.',
        training=old['training'],objective=old['objective'],searches=fits,selected_searches=selected,
        parity=parity,checks=checks,evaluations=evaluations,failed_evaluations=failures,
        mass_balance_max=ledger,volume_residual_max=volume,manifest=manifest,
        producer={str(Path(__file__).relative_to(model.ROOT)):native.sha(Path(__file__))},
        inputs={str(p.relative_to(model.ROOT)):native.sha(p) for p in [Path(state.__file__),Path(tuning.__file__),
            tuning.DRIVER,Path(diagnosis.__file__),Path(diagnosis.lab.__file__),Path(model.__file__)]},
        fit_artifacts={str(p.relative_to(model.ROOT)):native.sha(p) for p in sorted((CACHE/'fits').glob('*/*.json'))},
        build_record_sha256=native.sha(CACHE/'build/build.json'),
        restrictions=old['restrictions']+['Releasing one early boundary tests a model-family limitation; it does not identify a molecular mechanism.'])
    model.save(PRODUCT,result)
    for f in fits:
        if f['name'] in selected:print('selected',f['name'],f['training_log_rms'],f['assessment']['scores'],flush=True)
    print('checks',checks,flush=True);return result


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','cold','warm','warm_compact','bounds','constraints','nested','analyze'])
    args=p.parse_args()
    if args.action=='cold':
        for seed in (20261024,20261025):fit(-14.4,seed)
    elif args.action=='warm':
        for seed in (20261024,20261025):fit(-5.3,seed)
    elif args.action=='warm_compact':
        for seed in (20261024,20261025):fit(-5.3,seed,seed_diameter_um=40.)
    elif args.action=='bounds':
        fit(-14.4,20261026,generations=25,early_bounds=(.0625,4.))
        fit(-14.4,20261026,generations=25,ventilation_max=4.)
    elif args.action=='constraints':
        for temperature,diameter in ((-14.4,0.),(-5.3,40.)):
            for label,indices in [('axis_only',(0,2,3,4,5)),('axis_density_branch',(0,2,5)),('native_transport_mean',(0,5))]:
                fit(temperature,20261027,seed_diameter_um=diameter,generations=25,
                    fixed={i:float(tuning.BASE[i]) for i in indices},label=label)
    else:globals()[args.action]()


if __name__=='__main__':main()
