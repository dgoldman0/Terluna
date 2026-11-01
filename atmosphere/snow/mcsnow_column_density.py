"""Test evolving column deposition density while preserving the cold benchmark.

The size transition is a diagnostic hypothesis, not a measured surface law.
Earlier products and native builds remain unchanged.
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

from atmosphere.snow import mcsnow_growth_limits as limits

state, tuning = limits.state, limits.tuning
native, model, diagnosis = limits.native, limits.model, limits.diagnosis
HERE = Path(__file__).resolve().parent
CACHE = tuning.trace.CACHE / 'column_density'
PRODUCT = HERE / 'results/mcsnow_column_density.json'
DENSITY_RADIUS_MAX_UM = 300.
DENSITY_BOUNDS = (.25,4.)
WARM_NAME = 'axis_m5p3_seed40_20261025'
COLD_NAME = 'axis_m14p4_seed0_20261025'
PREFERRED_NAME = 'joint_join_balanced'
JOIN_TRAIN_AGE = 6-tuning.SEARCH_DT/60


def previous_fit(name=WARM_NAME):
    p=json.loads(limits.PRODUCT.read_text())
    return next(f for f in p['searches'] if f['name']==name)


def decode(x, seed_diameter_um=40., density_bounds=DENSITY_BOUNDS):
    x=np.asarray(x,dtype=float)
    if x.shape not in ((10,),(11,)) or not np.all(np.isfinite(x)) or np.any((x<0)|(x>1)):
        raise ValueError('Ten or eleven finite coordinates in [0,1] required')
    limits.early_coordinate(1.,density_bounds)
    if len(x)==11:seed_diameter_um=80.*x[10]
    return np.r_[limits.decode(x[:8],seed_diameter_um=seed_diameter_um),
                 DENSITY_RADIUS_MAX_UM*x[8],density_bounds[0]*(density_bounds[1]/density_bounds[0])**x[9]]


def density_multiplier(radius_m, transition_radius_um, early_ratio, late_multiplier):
    if not all(map(math.isfinite,(early_ratio,late_multiplier))) or min(early_ratio,late_multiplier)<=0:
        raise ValueError('Positive finite density controls required')
    w=state.stage_weight(radius_m,transition_radius_um)
    return late_multiplier*early_ratio**(1-w)


def deposition_source():
    replace=model.replace_once
    s=replace(limits.deposition_source(),
        '  REAL(wp), PUBLIC :: stage(4)=[0._wp,0._wp,0._wp,1._wp]',
        '  REAL(wp), PUBLIC :: stage(6)=[0._wp,0._wp,0._wp,1._wp,0._wp,1._wp]')
    s=replace(s,'  REAL(wp) :: stage_w','  REAL(wp) :: stage_w, density_w, density_factor')
    return replace(s,'''  if (rho_depo<rhoi) then
    if (rho_depo*tune(4)>rhoi) tune_cap_count=tune_cap_count+1
    rho_depo=min(rhoi,rho_depo*tune(4))
  end if''','''  if (rho_depo<rhoi) then
    density_factor=tune(4)
    ! Test gradual change in the density of added ice for the column habit.
    ! IGF<=1 retains the preceding plate treatment exactly.
    if (IGF>1._wp.and.stage(5)>0._wp) then
      density_w=1._wp/(1._wp+(stage(5)*1e-6_wp/(z3_4pi*V_tot)**z1_3)**4)
      density_factor=density_factor*stage(6)**(1._wp-density_w)
    end if
    if (rho_depo*density_factor>rhoi) tune_cap_count=tune_cap_count+1
    rho_depo=min(rhoi,rho_depo*density_factor)
  end if''')


def driver_source():
    return state.driver_source()


def build():
    limits.verify_build();path=CACHE/'build'
    if path.exists():return verify_build()
    baseline=native.verify_build(model.BASELINE)
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source',path/'source')
    (path/'source/src/mo_depo.f90').write_text(deposition_source())
    driver=path/'mcsnow_column_density_driver.f90';driver.write_text(driver_source())
    for name in native.TABLE_NAMES:shutil.copyfile(model.BASELINE/name,path/name)
    patch=path/'column_density.patch'
    patch.write_text(''.join(difflib.unified_diff(limits.deposition_source().splitlines(True),
        deposition_source().splitlines(True),fromfile='limits/mo_depo.f90',tofile='column_density/mo_depo.f90')))
    compiler=baseline['compiler_path'];sources=native.selected_sources(path/'source/src',driver,compiler)
    commands=[];objects=[]
    with (path/'build.log').open('w') as log:
        for i,source in enumerate(sources):
            obj=path/f'{i:02d}_{source.stem}.o';objects.append(str(obj))
            cmd=[compiler,*native.FLAGS,'-I'+str(path/'source/src'),'-I'+str(path),'-J'+str(path),'-c',str(source),'-o',str(obj)]
            commands.append(cmd)
            subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler,*native.FLAGS,*objects,'-o',str(path/'tuner')];commands.append(cmd)
        subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    record=dict(schema='terluna.atmosphere.mcsnow-column-density-build/1',
        baseline_sha256=native.sha(model.BASELINE/'build.json'),source_identity=baseline['source_identity'],
        parent_build_sha256=native.sha(limits.CACHE/'build/build.json'),
        compiled_sources={str(p):native.sha(p) for p in sources},compiler_path=compiler,
        compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],
        flags=native.FLAGS,commands=commands,
        files={p.name:native.sha(p) for p in [path/'tuner',patch,path/'build.log',*[path/n for n in native.TABLE_NAMES]]})
    model.save(path/'build.json',record)
    return record


def verify_build():
    limits.verify_build();path=CACHE/'build';r=json.loads((path/'build.json').read_text())
    if r['parent_build_sha256']!=native.sha(limits.CACHE/'build/build.json'):raise ValueError('Parent build changed')
    if r['baseline_sha256']!=native.sha(model.BASELINE/'build.json'):raise ValueError('Baseline changed')
    for name,digest in r['compiled_sources'].items():
        if native.sha(name)!=digest:raise ValueError('Compiled source changed: '+name)
    for name,digest in r['files'].items():
        if native.sha(path/name)!=digest:raise ValueError('Build artifact changed: '+name)
    if (path/'source/src/mo_depo.f90').read_text()!=deposition_source():raise ValueError('Declared source changed')
    if (path/'mcsnow_column_density_driver.f90').read_text()!=driver_source():raise ValueError('Declared driver changed')
    return r


def batch(entries,dt=tuning.SEARCH_DT,seed_diameter_um=40.,free_compact=False,density_bounds=DENSITY_BOUNDS):
    experiment=dict(family='column_density',seed_diameter_um=seed_diameter_um,
            density_transition_radius_max_um=DENSITY_RADIUS_MAX_UM,density_early_ratio_bounds=list(density_bounds),
            density_rule='native_density * late_multiplier * early_ratio**(1-w_density), capped at solid ice',
            transition_power=4,column_gate='native IGF > 1; plate law unchanged')
    if free_compact:
        experiment['compact_limit_bounds_um']=[0.,80.]
        experiment['seed_diameter_um']=None
    def decoder(x):
        if len(x)!=(11 if free_compact else 10):raise ValueError('Coordinate count differs from experiment')
        return decode(x,seed_diameter_um,density_bounds)
    return tuning.Batch(entries,dt,cache=CACHE,decoder=decoder,experiment=experiment)


def evaluate(x,dt=1/128,temperature=-5.3,seed_diameter_um=40.,all_starts=True,density_bounds=DENSITY_BOUNDS,ages=None):
    entries=[e for e in tuning.cases() if e['temperature_C']==temperature and (all_starts or e['initialization']!='volume')]
    if ages is not None:
        entries=copy.deepcopy(entries)
        for e in entries:e['sample_ages_min']=list(ages)
    rows=[];failed=[];ledger=volume=density=0.
    for e in entries:
        with batch([e],dt,seed_diameter_um,free_compact=len(x)==11,density_bounds=density_bounds) as b:r=tuning.assessment(b,x)
        if 'error' in r:failed.append(dict(case=e['name'],error=r['error']));continue
        rows.extend(r['rows']);ledger=max(ledger,r['mass_balance_max'])
        volume=max(volume,r['volume_residual_max']);density=max(density,r['maximum_density_fraction'])
    return dict(rows=rows,scores=tuning.score_rows(rows),failed_cases=failed,
        mass_balance_max=ledger,volume_residual_max=volume,maximum_density_fraction=density)


def classify_assessment(assessment,training_ages):
    """Label fitting, supplied initial states and comparisons for this objective."""
    for row in assessment['rows']:
        if row['initialization']=='volume':row['group']='geometry_sensitivity'
        elif row['initialization']=='axes' and row['age_min']==3.:row['group']='initial_state'
        elif row['age_min'] in training_ages:row['group']='training'
        else:row['group']='withheld_age'
    assessment['scores']=tuning.score_rows(assessment['rows'])
    return assessment


def additional_ages(x,dt=1/128,density_bounds=DENSITY_BOUNDS,exclude_ages=()):
    entries=copy.deepcopy(limits.cases(-5.3))
    for entry in entries:entry['sample_ages_min']=[age for age in (4.,6.,8.,12.,25.) if age not in exclude_ages]
    with batch(entries,dt,free_compact=len(x)==11,density_bounds=density_bounds) as b:
        result=tuning.assessment(b,x)
    if 'error' in result:raise ValueError('Additional-age check failed: '+result['error'])
    return result



def fitting_cases(training_ages):
    entries=copy.deepcopy(limits.cases(-5.3))
    ages=sorted(set(diagnosis.sample_ages(-5.3))|set(training_ages))
    if ages!=list(diagnosis.sample_ages(-5.3)):
        for entry in entries:entry['sample_ages_min']=ages
    return entries


def interval_assessment(x,dt=1/128,density_bounds=DENSITY_BOUNDS,training_ages=tuning.TRAIN_AGES):
    """Quarter-minute states and both sides of the fitted length-history join."""
    ages=sorted(set(np.arange(3.,30.25,.25).tolist()+[6-1/32/60,6+1/32/60]))
    rows=[];ledger=volume=density=0.
    # The preserved native driver accepts at most eight output times per case.
    for j in range(0,len(ages),8):
        entries=copy.deepcopy(limits.cases(-5.3))
        for entry in entries:entry['sample_ages_min']=ages[j:j+8]
        with batch(entries,dt,free_compact=len(x)==11,density_bounds=density_bounds) as b:
            result=tuning.assessment(b,x)
        if 'error' in result:raise ValueError('Interval check failed: '+result['error'])
        rows.extend(result['rows']);ledger=max(ledger,result['mass_balance_max'])
        volume=max(volume,result['volume_residual_max']);density=max(density,result['maximum_density_fraction'])
    rows.sort(key=lambda r:(r['initialization'],r['age_min']))
    return classify_assessment(dict(rows=rows,mass_balance_max=ledger,volume_residual_max=volume,
        maximum_density_fraction=density,failed_cases=[],
        evidence='Sampled comparison with the same published fitted population histories; additional ages informed diagnosis and are not independent validation.'),training_ages)


def comparison_metrics(assessment,training_ages=tuning.TRAIN_AGES):
    rows=[r for r in assessment['rows'] if r['initialization'] in ('sphere','axes')]
    def worst(selected):return max(abs(v-1) for r in selected for v in r['ratio'][:3])
    return dict(maximum_30min_relative_error=worst([r for r in rows if r['age_min']==30.]),
        maximum_all_age_relative_error=worst(rows),
        maximum_training_relative_error=worst([r for r in rows if r['age_min'] in training_ages]))



def objective_residual(b,x,training_ages,age_weights=None):
    residual=tuning.residual(b,x,ages=training_ages)
    if not age_weights:return residual
    weights={float(age):float(weight) for age,weight in age_weights.items()}
    if any(age not in training_ages or not math.isfinite(weight) or weight<=0 for age,weight in weights.items()):
        raise ValueError('Weights require fitting ages and positive finite values')
    selected=[(i,age) for i,age in b.rows if age in training_ages and age>b.entries[i]['initial_age_min']]
    counts={i:sum(k==i for k,_ in selected) for i in range(len(b.entries))}
    totals={i:sum(weights.get(age,1.) for k,age in selected if k==i) for i in counts}
    factors=np.repeat([math.sqrt(weights.get(age,1.)*counts[i]/totals[i]) for i,age in selected],3)
    return residual*factors


def fit(seed=20261024,generations=35,local_only=False,initial=None,label=None,free_compact=False,density_bounds=DENSITY_BOUNDS,training_ages=tuning.TRAIN_AGES,fixed_late_axis=None,age_weights=None):
    build();prior=previous_fit()
    fallback=np.r_[prior['coordinates'],0.,limits.early_coordinate(1.,density_bounds)]
    if initial is None:
        initial=fallback.copy();initial[3]=limits.early_coordinate(.74)
        initial[8]=80/DENSITY_RADIUS_MAX_UM;initial[9]=limits.early_coordinate(.55,density_bounds)
    initial=np.asarray(initial)
    if free_compact:
        fallback=np.r_[fallback,.5]
        if initial.shape==(10,):initial=np.r_[initial,.5]
    elif initial.shape!=(10,):raise ValueError('Unexpected coordinate count')
    # The column's branching-radius control was independently shown to be inactive.
    fixed={4:float(initial[4])}
    if fixed_late_axis is not None:fixed[1]=limits.early_coordinate(fixed_late_axis)
    for i,value in fixed.items():initial[i]=value;fallback[i]=value
    free=[i for i in range(len(initial)) if i not in fixed]
    def expand(z):
        x=initial.copy();x[free]=z;return x
    name=label or (('compact_' if free_compact else '')+('local' if local_only else 'global')+f'_{seed}')
    path=CACHE/'fits'/name;path.mkdir(parents=True,exist_ok=True)
    spec=dict(name=name,family='column_density',seed=seed,generations=0 if local_only else generations,
        population_multiplier=6,initial=initial.tolist(),fallback=fallback.tolist(),fixed={str(i):v for i,v in sorted(fixed.items())},
        fitting_ages_min=list(training_ages),temperature_C=-5.3,seed_diameter_um=40.,
        dt_s=tuning.SEARCH_DT,density_radius_upper_um=DENSITY_RADIUS_MAX_UM,
        density_early_ratio_bounds=list(density_bounds),previous_product_sha256=native.sha(limits.PRODUCT),
        build_sha256=native.sha(CACHE/'build/build.json'),scipy_version=scipy.__version__)
    if free_compact:spec.update(compact_limit_fitted=True,compact_limit_bounds_um=[0.,80.])
    if age_weights:spec['age_weights']={str(age):weight for age,weight in sorted(age_weights.items())}
    if (path/'specification.json').exists():
        if json.loads((path/'specification.json').read_text())!=spec:raise ValueError('Search changed')
    else:model.save(path/'specification.json',spec)
    if (path/'result.json').exists():return json.loads((path/'result.json').read_text())
    with batch(fitting_cases(training_ages),free_compact=free_compact,density_bounds=density_bounds) as b:
        start=time.monotonic();last=start;generation=0
        def objective(x):return float(np.sum(objective_residual(b,x,training_ages,age_weights)**2))
        def callback(z,convergence):
            nonlocal last,generation
            generation+=1
            if time.monotonic()-last>15:
                print(name,'generation',generation,'rms',math.sqrt(objective(expand(z))),flush=True);last=time.monotonic()
        neutral=np.r_[prior['coordinates'],80/DENSITY_RADIUS_MAX_UM,limits.early_coordinate(1.,density_bounds)]
        if free_compact:neutral=np.r_[neutral,.5]
        for i,value in fixed.items():neutral[i]=value
        candidates=[neutral,fallback,initial]
        de_record=dict(skipped=True,evaluations=0,generations=0)
        if not local_only:
            de=differential_evolution(lambda z:objective(expand(z)),[(0.,1.)]*len(free),seed=seed,
                maxiter=generations,popsize=6,tol=1e-5,atol=1e-9,workers=1,updating='immediate',
                polish=False,x0=initial[free],callback=callback)
            candidates.append(expand(de.x))
            de_record=dict(success=bool(de.success),message=str(de.message),evaluations=int(de.nfev),generations=int(de.nit))
        start_local=min(candidates,key=objective)
        local=least_squares(lambda z:objective_residual(b,expand(z),training_ages,age_weights),start_local[free],bounds=(np.zeros(len(free)),np.ones(len(free))),
            diff_step=.005,x_scale='jac',max_nfev=400,ftol=1e-9,xtol=1e-8,gtol=1e-8)
        candidates.append(expand(local.x));x=min(candidates,key=objective)
        result=dict(name=name,coordinates=x.tolist(),physical_controls=decode(x,density_bounds=density_bounds).tolist(),specification=spec,
            training_loss=objective(x),training_log_rms=math.sqrt(objective(x)),
            differential_evolution=de_record,local_refinement=dict(success=bool(local.success),message=str(local.message),evaluations=int(local.nfev)),
            elapsed_s=time.monotonic()-start,evaluation_directory=str(b.path.relative_to(model.ROOT)))
    ages=sorted(set(diagnosis.sample_ages(-5.3))|set(training_ages))
    result['assessment']=classify_assessment(evaluate(x,density_bounds=density_bounds,ages=ages),training_ages)
    result['metrics']=comparison_metrics(result['assessment'],training_ages)
    model.save(path/'result.json',result)
    print(name,'completed',result['training_log_rms'],result['metrics'],result['physical_controls'],flush=True)
    return result


def parity():
    build();result=[]
    for name in (COLD_NAME,WARM_NAME):
        previous=previous_fit(name)
        x=np.r_[previous['coordinates'],0.,.5]
        assessment=evaluate(x,temperature=previous['temperature_C'],seed_diameter_um=previous['seed_diameter_um'])
        difference=state.comparison_difference(assessment,previous['assessment'])
        result.append(dict(name=name,states=len(assessment['rows']),maximum_relative_difference=difference,passed=difference<1e-12))
    model.save(CACHE/'parity.json',result);print('parity',result,flush=True)
    return result



def probes(x,dt=1/128,duration=1/32,density_bounds=DENSITY_BOUNDS):
    result=[]
    for age in (3.5,5.5,8.,30.):
        for mapping in diagnosis.MAPPINGS:
            target=diagnosis.fit_rates(-5.3,age,mapping)
            if target['density_kg_m3']>diagnosis.ice_density():continue
            e=copy.deepcopy(next(e for e in tuning.cases() if e['temperature_C']==-5.3 and e['initialization']==mapping))
            e['initial_age_min']=age;e['sample_ages_min']=[age,age+duration/60]
            e['configuration'].update(initial_mass_kg=target['mass_kg'],
                initial_equatorial_m=target['model_a_m'],initial_polar_m=target['model_c_m'])
            with batch([e],dt,free_compact=len(x)==11,density_bounds=density_bounds) as b:r=b.evaluate(x)
            if r['error']:raise ValueError('Supplied-state probe failed: '+r['error'])
            start,end=np.asarray(r['rows'])[:,:3]
            rates=np.r_[(end[0]-start[0])/duration,np.log(end[1:]/start[1:])/duration]
            reference=np.array([target['mass_rate_kg_s'],target['a_log_rate_s'],target['c_log_rate_s']])
            volumes=math.pi/6*np.array([start[1]**2*start[2],end[1]**2*end[2]])
            result.append(dict(age_min=age,mapping=mapping,rates=rates.tolist(),target_rates=reference.tolist(),
                ratio=(rates/reference).tolist(),effective_axis_ratio=float(rates[2]/rates[1]),
                target_effective_axis_ratio=float(reference[2]/reference[1]),
                deposition_density_kg_m3=float((end[0]-start[0])/(volumes[1]-volumes[0])),
                target_deposition_density_kg_m3=target['deposition_density_kg_m3']))
    return result


def early_history(x,dt=1/128,density_bounds=DENSITY_BOUNDS):
    e=copy.deepcopy(next(e for e in limits.cases(-5.3) if e['initialization']=='sphere'))
    ages=[0.,.05,.15,.3,.5,1.,2.,3.];e['sample_ages_min']=ages
    with batch([e],dt,free_compact=len(x)==11,density_bounds=density_bounds) as b:r=b.evaluate(x)
    if r['error']:raise ValueError('Early trajectory failed: '+r['error'])
    values=tuning.reported(b,r)
    return dict(evidence='Unobserved early model history from the assumed 4 micrometre sphere; not fitting targets.',
        rows=[dict(age_min=age,values=v.tolist(),aspect_ratio=float(v[2]/v[1])) for age,v in zip(ages,values)])


def seed_sensitivity(x,density_bounds=DENSITY_BOUNDS):
    original=next(e for e in limits.cases(-5.3) if e['initialization']=='sphere');result=[]
    for phi in (.25,.5,1.,2.,4.):
        e=copy.deepcopy(original);d=e['configuration']['initial_equatorial_m']
        e['configuration'].update(initial_equatorial_m=d*phi**(-1/3),initial_polar_m=d*phi**(2/3))
        with batch([e],1/128,free_compact=len(x)==11,density_bounds=density_bounds) as b:r=tuning.assessment(b,x)
        result.append(dict(initial_aspect_ratio=phi,initial_volume_equivalent_diameter_m=d,assessment=r))
    return result


def analyze():
    verify_build();checks=parity();previous_cold=previous_fit(COLD_NAME);previous_warm=previous_fit()
    fits=[json.loads(p.read_text()) for p in sorted((CACHE/'fits').glob('*/result.json'))]
    if not fits:raise ValueError('No completed fits')
    selected=[]
    groups=sorted({(len(f['coordinates']),tuple(f['specification']['density_early_ratio_bounds']),tuple(f['specification']['fitting_ages_min']),tuple(sorted(f['specification']['fixed'].items())),tuple(sorted(f['specification'].get('age_weights',{}).items()))) for f in fits})
    for count,density_bounds,training_ages,fixed,weights in groups:
        age_weights=dict(weights)
        group=[f for f in fits if len(f['coordinates'])==count and tuple(f['specification']['density_early_ratio_bounds'])==density_bounds
               and tuple(f['specification']['fitting_ages_min'])==training_ages and tuple(sorted(f['specification']['fixed'].items()))==fixed
               and tuple(sorted(f['specification'].get('age_weights',{}).items()))==weights]
        best=min(group,key=lambda f:f['training_loss']);selected.append(best['name'])
        ages=sorted(set(diagnosis.sample_ages(-5.3))|set(training_ages))
        for f in group:
            f['assessment']=classify_assessment(f['assessment'],training_ages)
            x=f['coordinates'];fine=evaluate(x,all_starts=False,density_bounds=density_bounds,ages=ages)
            coarse=evaluate(x,dt=tuning.SEARCH_DT,all_starts=False,density_bounds=density_bounds,ages=ages)
            diff=state.comparison_difference(coarse,fine)
            checks.append(dict(name='search_timestep_'+f['name'],maximum_relative_difference=diff,passed=diff<.01))
            with batch(fitting_cases(training_ages),1/128,free_compact=len(x)==11,density_bounds=density_bounds) as b:
                f['fine_training_loss']=float(np.sum(objective_residual(b,x,training_ages,age_weights)**2))
        x=best['coordinates'];fine=classify_assessment(evaluate(x,density_bounds=density_bounds,ages=ages),training_ages);best['assessment']=fine
        best['metrics']=comparison_metrics(fine,training_ages);coarse=evaluate(x,dt=1/32,density_bounds=density_bounds,ages=ages)
        diff=state.comparison_difference(coarse,fine)
        checks.append(dict(name='all_start_timestep_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
        cold_x=np.r_[previous_cold['coordinates'],x[8:10]]
        cold=evaluate(cold_x,temperature=-14.4,seed_diameter_um=0.,density_bounds=density_bounds)
        diff=state.comparison_difference(cold,previous_cold['assessment'])
        best['cold_control']=dict(coordinates=cold_x.tolist(),assessment=cold,maximum_relative_difference=diff,
            evidence='The frozen cold controls with the new density settings; the native IGF column condition leaves the cold plate law unchanged.')
        checks.append(dict(name='cold_unchanged_'+best['name'],maximum_relative_difference=diff,passed=diff==0.))
        best['fixed_state_probes']=probes(x,density_bounds=density_bounds);ref=probes(x,1/256,1/64,density_bounds=density_bounds)
        diff=max(abs(a/b-1) for lo,hi in zip(best['fixed_state_probes'],ref) for a,b in zip(lo['rates'],hi['rates']))
        checks.append(dict(name='supplied_state_refinement_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
        best['early_model_history']=early_history(x,density_bounds=density_bounds);coarse=early_history(x,1/32,density_bounds=density_bounds)
        diff=max(abs(a/b-1) for lo,hi in zip(coarse['rows'],best['early_model_history']['rows']) for a,b in zip(lo['values'],hi['values']))
        checks.append(dict(name='early_timestep_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
        endpoint=next(r for r in fine['rows'] if r['initialization']=='sphere' and r['age_min']==3.)
        diff=max(abs(a/b-1) for a,b in zip(best['early_model_history']['rows'][-1]['values'],endpoint['values']))
        checks.append(dict(name='early_endpoint_parity_'+best['name'],maximum_relative_difference=diff,passed=diff==0.))
        if 3. in training_ages:
            additional=additional_ages(x,density_bounds=density_bounds,exclude_ages=training_ages)
            additional_coarse=additional_ages(x,1/32,density_bounds=density_bounds,exclude_ages=training_ages)
            diff=state.comparison_difference(additional_coarse,additional)
            checks.append(dict(name='additional_age_timestep_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
            best['additional_age_assessment']=additional
            best['metrics']['maximum_additional_age_relative_error']=max(abs(v-1) for r in additional['rows'] for v in r['ratio'][:3])
            best['metrics']['maximum_combined_age_relative_error']=max(best['metrics']['maximum_all_age_relative_error'],best['metrics']['maximum_additional_age_relative_error'])
        if best['name'] in ('seed_state_density_floor','seed_state_early_interval','early_interval_density_profile','joint_early_constraints','joint_column_late_limit','joint_join_constraint','joint_join_weighted',PREFERRED_NAME):
            interval=interval_assessment(x,density_bounds=density_bounds,training_ages=training_ages)
            coarse=interval_assessment(x,1/32,density_bounds=density_bounds,training_ages=training_ages)
            diff=state.comparison_difference(coarse,interval)
            checks.append(dict(name='interval_timestep_'+best['name'],maximum_relative_difference=diff,passed=diff<.01))
            best['interval_assessment']=interval
            best['metrics']['maximum_interval_relative_error']=max(abs(v-1) for r in interval['rows'] for v in r['ratio'][:3])
        best['seed_shape_sensitivity']=seed_sensitivity(x,density_bounds=density_bounds)
        reset=np.array(x);reset[9]=limits.early_coordinate(1.,density_bounds)
        ablation=evaluate(reset,density_bounds=density_bounds,ages=ages)
        best['density_ablation']=dict(assessment=ablation,metrics=comparison_metrics(ablation,training_ages),
            evidence='New density ratio set to one with every other fitted control fixed; not a refitted restricted family.')
    manifest=[];evaluations=failures=0;ledger=volume=0.
    for p in sorted((CACHE/'evaluations').glob('*/evaluations.jsonl')):
        rows=[json.loads(line) for line in p.read_text().splitlines()]
        evaluations+=len(rows);failures+=sum(r['error'] is not None for r in rows)
        for r in rows:
            for v in r['rows']:ledger=max(ledger,abs(v[4]));volume=max(volume,abs(v[5]))
        manifest.append(dict(directory=str(p.parent.relative_to(model.ROOT)),evaluations=len(rows),
            files={f.name:native.sha(f) for f in sorted(p.parent.iterdir()) if f.is_file()}))
    checks.extend([dict(name='mass_ledger',maximum_relative_difference=ledger,passed=ledger<1e-10),
                   dict(name='spheroid_volume',maximum_relative_difference=volume,passed=volume<1e-10)])
    preferred=next(f for f in fits if f['name']==PREFERRED_NAME)
    worst=max(preferred['metrics']['maximum_all_age_relative_error'],preferred['metrics']['maximum_interval_relative_error'])
    result=dict(schema='terluna.atmosphere.mcsnow-column-density/1',
        evidence='Diagnostic fits of evolving deposition density for columns; no adopted temperature/humidity transfer or lunar aggregate distribution.',
        quantities=['mass_kg','equatorial_m','polar_m','fall_speed_m_s'],
        units=dict(age_min='minutes',values=['kg','m','m','m/s'],ratio='model / laboratory fitted history'),
        density_rule='For native IGF > 1: rho_depo=min(rho_ice,rho_native * late_density_multiplier * early_density_ratio**(1-w)); w=1/(1+(R_density/r_equivalent)**4). The existing compact stage overrides this with solid deposition.',
        physical_control_order=json.loads(limits.PRODUCT.read_text())['physical_control_order']+['density_transition_radius_um','early_density_ratio'],
        training=dict(temperature_C=-5.3,initializations=['sphere','axes'],ages_min=list(tuning.TRAIN_AGES),seed_constraint_ages_min=[3.,5.,10.,20.],interval_constraint_ages_min=[3.,5.,5.75,10.,20.],final_fitting_ages_min=[3.,4.,5.,JOIN_TRAIN_AGE,10.,20.],final_age_weights={str(JOIN_TRAIN_AGE):1.25},reserved_in_all_families=[15.,30.],additional_comparison_ages_min=[4.,6.,8.,12.,25.],quantities=['mass','width','polar']),
        selection='Minimum stated training loss within each fixed/free compact-stage, density-bound, fitting-age, fixed-control and age-weight family. The preferred diagnostic follows explicit early-age constraints and sensitivity comparisons; all inspected comparison ages have informed design.',
        target_assessment=dict(preferred_search=PREFERRED_NAME,maximum_checked_relative_error=worst,cold_maximum_change=preferred['cold_control']['maximum_relative_difference'],attained=worst<.1 and preferred['cold_control']['maximum_relative_difference']==0.),
        target=dict(maximum_30min_relative_error=.1,maximum_checked_relative_error=.1,quantities=['mass','width','polar'],starts=['sphere','axes'],cold='exactly preserve the preceding cold fit'),
        previous_product_path=str(limits.PRODUCT.relative_to(model.ROOT)),previous_product_sha256=native.sha(limits.PRODUCT),
        previous_warm=previous_warm,previous_cold=previous_cold,searches=fits,selected_searches=selected,preferred_search=PREFERRED_NAME,checks=checks,
        evaluations=evaluations,failed_evaluations=failures,mass_balance_max=ledger,volume_residual_max=volume,
        manifest=manifest,build_record_sha256=native.sha(CACHE/'build/build.json'),
        producer={str(Path(__file__).relative_to(model.ROOT)):native.sha(Path(__file__))},
        inputs={str(p.relative_to(model.ROOT)):native.sha(p) for p in [Path(limits.__file__),Path(state.__file__),Path(tuning.__file__),
            tuning.DRIVER,Path(diagnosis.__file__),Path(diagnosis.lab.__file__),Path(model.__file__)]},
        fit_artifacts={str(p.relative_to(model.ROOT)):native.sha(p) for p in sorted((CACHE/'fits').glob('*/*.json'))},
        restrictions=['Separate temperature-specific controls; no inferred interpolation.',
            'Fitted population histories share measurement and geometry assumptions; the reserved ages are not an untouched validation dataset.',
            'Early seed geometry and compact-stage threshold are unmeasured assumptions.',
            'Density and radius transitions are diagnostic controls, not measured attachment kinetics.',
            'Finite searches do not establish a global optimum; native failures remain in the journals.'])
    model.save(PRODUCT,result)
    print('Selected',[(f['name'],f['metrics']) for f in fits if f['name'] in selected],flush=True)
    print('Checks',len(checks),'passed',sum(c['passed'] for c in checks),flush=True)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','parity','local','global','compact','density_bounds','seed_state','late_control','density_floor','early_interval','interval_profile','joint_early','joint_late','joint_join','join_weights','analyze'])
    args=p.parse_args()
    if args.action=='analyze':analyze()
    elif args.action=='build':build()
    elif args.action=='parity':parity()
    elif args.action=='local':fit(local_only=True)
    elif args.action=='global':
        initial=json.loads((CACHE/'fits/local_20261024/result.json').read_text())['coordinates']
        for seed in (20261024,20261025):fit(seed=seed,initial=initial)
    elif args.action=='density_floor':
        previous=json.loads((CACHE/'fits/seed_state_native_late_axis/result.json').read_text())
        initial=np.array(previous['coordinates']);initial[9]=0.
        fit(local_only=True,initial=initial,free_compact=True,training_ages=(3.,5.,10.,20.),
            fixed_late_axis=1.,density_bounds=(.5,4.),label='seed_state_density_floor')
    elif args.action=='early_interval':
        previous=json.loads((CACHE/'fits/seed_state_density_floor/result.json').read_text())
        fit(local_only=True,initial=previous['coordinates'],free_compact=True,training_ages=(3.,5.,5.75,10.,20.),
            fixed_late_axis=1.,density_bounds=(.5,4.),label='seed_state_early_interval')
    elif args.action in ('interval_profile','joint_early','joint_late','joint_join','join_weights'):
        continuations={
            'interval_profile':('seed_state_early_interval','early_interval_density_profile',(3.,5.,5.75,10.,20.),(.4,4.),1.),
            'joint_early':('early_interval_density_profile','joint_early_constraints',(3.,4.,5.,5.75,10.,20.),DENSITY_BOUNDS,1.),
            'joint_late':('joint_early_constraints','joint_column_late_limit',(3.,4.,5.,5.75,10.,20.),DENSITY_BOUNDS,1.5),
            'joint_join':('joint_column_late_limit','joint_join_constraint',(3.,4.,5.,JOIN_TRAIN_AGE,10.,20.),DENSITY_BOUNDS,1.5),
            'join_weights':('joint_join_constraint',PREFERRED_NAME,(3.,4.,5.,JOIN_TRAIN_AGE,10.,20.),DENSITY_BOUNDS,1.5)}
        source,label,ages,bounds,late=continuations[args.action]
        previous=json.loads((CACHE/'fits'/source/'result.json').read_text())
        initial=np.array(previous['coordinates'])
        if list(bounds)!=previous['specification']['density_early_ratio_bounds']:
            initial[9]=limits.early_coordinate(previous['physical_controls'][-1],bounds)
        if args.action=='join_weights':
            for weight,label in [(1.5,'joint_join_weighted'),(1.25,PREFERRED_NAME)]:
                fit(local_only=True,initial=initial.copy(),free_compact=True,training_ages=ages,fixed_late_axis=late,
                    density_bounds=bounds,age_weights={JOIN_TRAIN_AGE:weight},label=label)
        else:fit(local_only=True,initial=initial,free_compact=True,training_ages=ages,fixed_late_axis=late,density_bounds=bounds,label=label)
    elif args.action=='late_control':
        previous=json.loads((CACHE/'fits/seed_state_local/result.json').read_text())
        fit(local_only=True,initial=previous['coordinates'],free_compact=True,training_ages=(3.,5.,10.,20.),
            fixed_late_axis=1.,label='seed_state_native_late_axis')
    elif args.action=='seed_state':
        for name,label in [('compact_local_20261024','seed_state_local'),('compact_local_declining','seed_state_declining')]:
            previous=json.loads((CACHE/'fits'/name/'result.json').read_text())
            fit(local_only=True,initial=previous['coordinates'],free_compact=True,training_ages=(3.,5.,10.,20.),label=label)
    elif args.action=='density_bounds':
        previous=json.loads((CACHE/'fits/compact_local_20261024/result.json').read_text())
        initial=np.array(previous['coordinates']);bounds=(.0625,4.)
        initial[9]=limits.early_coordinate(previous['physical_controls'][-1],bounds)
        fit(local_only=True,initial=initial,free_compact=True,density_bounds=bounds,label='compact_density_lower_bound')
    elif args.action=='compact':
        initial=json.loads((CACHE/'fits/local_20261024/result.json').read_text())['coordinates']
        fit(local_only=True,initial=initial,free_compact=True)
        declining=np.r_[initial,55/80]
        declining[1]=limits.early_coordinate(2.2)
        declining[3]=limits.early_coordinate(.75)
        declining[6]=60/state.MAX_RADIUS_UM
        declining[7]=limits.early_coordinate(3.5)
        declining[8]=70/DENSITY_RADIUS_MAX_UM
        declining[9]=limits.early_coordinate(.55)
        fit(local_only=True,initial=declining,free_compact=True,label='compact_local_declining')


if __name__=='__main__':main()
