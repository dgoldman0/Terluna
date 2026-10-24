"""Joint growth/settling diagnostics under a frozen pointwise acceptance contract.

Bulk-relative deposition density and equal normal compact-stage growth are
testable closures, not measured attachment laws. Rejected candidates are kept.
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
from scipy.optimize import least_squares

from atmosphere.snow import mcsnow_column_density as column
from atmosphere.snow import mcsnow_joint_acceptance as acceptance
from atmosphere.snow import mcsnow_joint_drag as drag

native,model,diagnosis,tuning=column.native,column.model,column.diagnosis,column.tuning
HERE=Path(__file__).resolve().parent
CACHE=drag.CACHE/'growth_sharp'
PRODUCT=HERE/'results/mcsnow_joint_growth.json'
NEUTRAL=dict(relative_density=False,normal_compact=False,density_power=4.,
             original_best=0.,span_width=0.,prism_width_fraction=0.)
STATE=dict(relative_density=True,normal_compact=True,density_power=16.,
           original_best=1.,span_width=0.,prism_width_fraction=3/math.pi)


def baseline_fit(temperature):
    p=acceptance.load_baseline()
    return next(f for f in p['searches'] if f['name']==p['preferred_search']) if temperature==-5.3 else p['previous_cold']


def baseline_coordinates(temperature):
    fit=baseline_fit(temperature);x=np.asarray(fit['coordinates'])
    if temperature==-14.4:x=np.r_[x,0.,.5,0.]
    return np.r_[x,0.,.5]


def decode(x,setup):
    x=np.asarray(x,dtype=float)
    if x.shape!=(13,) or not np.all(np.isfinite(x)) or np.any((x<0)|(x>1)):
        raise ValueError('Thirteen finite normalized controls required')
    if set(setup)!=set(NEUTRAL) or not 1<=setup['density_power']<=32:
        raise ValueError('Invalid joint experiment')
    if any(not 0<=setup[k]<=1 for k in ('original_best','span_width','prism_width_fraction')):
        raise ValueError('Invalid drag controls')
    return np.r_[column.decode(x[:11]),float(setup['relative_density']),float(setup['normal_compact']),
        setup['density_power'],150*x[11],.25*16**x[12],setup['original_best'],setup['span_width'],setup['prism_width_fraction']]


def deposition_source():
    replace=model.replace_once;s=column.deposition_source()
    s=replace(s,'  REAL(wp), PUBLIC :: tune(6)=1._wp',
        '  REAL(wp), PUBLIC :: tune(6)=1._wp\n  REAL(wp), PUBLIC :: joint(5)=[0._wp,0._wp,4._wp,0._wp,1._wp]')
    s=replace(s,'(stage(1)*1e-6_wp/(z3_4pi*V_tot)**z1_3)**4)', '(stage(1)*1e-6_wp/(z3_4pi*V_tot)**z1_3)**joint(3))')
    s=replace(s,'  C = capacitance( sp, Re ) * tune(1)', '''  C = capacitance( sp, Re ) * tune(1)
  ! Size-dependent transport is an explicit diagnostic control.
  if (joint(4)>0._wp) C=C*joint(5)**(1._wp- &
    1._wp/(1._wp+(joint(4)*1e-6_wp/(z3_4pi*sp%V_i)**z1_3)**4))''')
    s=replace(s,'  ! Diagnostic: compact deposition at small size, approaching the tuned density.', '''  ! Added-ice density responds to the present envelope density.
  ! Preserve solid deposition before the plate branching criterion.
  if (joint(1)>0.5_wp.and.(IGF>1._wp.or.branch_crit)) then
    density_w=1._wp
    if (stage(5)>0._wp) density_w=1._wp/(1._wp+ &
      (stage(5)*1e-6_wp/(z3_4pi*V_tot)**z1_3)**joint(3))
    density_factor=(sp%m_i+sp%m_r+sp%m_f)/V_tot*tune(4)*stage(6)**(1._wp-density_w)
    if (density_factor>rhoi) tune_cap_count=tune_cap_count+1
    rho_depo=min(rhoi,density_factor)
  end if
  ! Diagnostic: compact deposition at small size, approaching the tuned density.''')
    return replace(s,'  if (sp%d<stage(2)*1e-6_wp) IGF_vent=1._wp', '''  if (sp%d<stage(2)*1e-6_wp) then
    IGF_vent=1._wp
    ! Equal absolute normal face speeds round unequal axes as they grow.
    if (joint(2)>0.5_wp) IGF_vent=1._wp/sp%phi
  end if''')


def velocity_source():
    s=(model.BASELINE/'source/src/mo_velocity.f90').read_text();replace=model.replace_once
    s=replace(s,'CONTAINS','  REAL(wp), PUBLIC :: joint_drag(3)=0._wp\n\nCONTAINS')
    s=replace(s,'  REAL(wp) :: fx','  REAL(wp) :: fx, q_projected')
    s=replace(s,'      d = (4.0_wp * V_tot / pi / sp%phi)**(1.0_wp/3.0_wp)',
        '''      d = (4.0_wp * V_tot / pi / sp%phi)**(1.0_wp/3.0_wp)
      d = d*((sp%d/sp%phi)/d)**joint_drag(2)''')
    s=replace(s,'    ! prolate assumption q = A / A_ce ', '''    ! Prolate interpolation with an explicitly declared prism side projection.
    q_projected=sp%p_area/(pi_4*sp%d**2/sp%phi)
    if (joint_drag(3)>0._wp.and.sp%mm<1.01_wp) q_projected=z4_pi*joint_drag(3)
    ! prolate assumption q = A / A_ce ''')
    s=replace(s,'mtot, sp%p_area / (pi_4 * sp%d**2 / sp%phi),    &','mtot, q_projected,    &')
    return replace(s,'MAX(alpha,1.0_wp)*q**(1.0_wp/4.0_wp))                 ! porosity correction also for q > 1',
        '''MAX(alpha,1.0_wp)* &
      q**(1.0_wp/4.0_wp)*(MAX(q**(1.0_wp/4.0_wp),q)/q**(1.0_wp/4.0_wp))**joint_drag(1))''')


def driver_source():
    s=column.driver_source();replace=model.replace_once
    s=replace(s,'use mo_velocity, only: vterm_snow','use mo_velocity, only: vterm_snow, joint_drag')
    s=replace(s,'sp_depo, tune, stage, tune_cap_count','sp_depo, tune, stage, joint, tune_cap_count')
    s=replace(s,'read(input_unit,*,iostat=ios) tune,stage','read(input_unit,*,iostat=ios) tune,stage,joint,joint_drag')
    s=replace(s,'ieee_is_finite([tune,stage])','ieee_is_finite([tune,stage,joint,joint_drag])')
    return s.replace('times(8,ncase)','times(256,ncase)').replace('counts(i)>8','counts(i)>256')


def build():
    column.verify_build();path=CACHE/'build'
    if path.exists():return verify_build()
    baseline=native.verify_build(model.BASELINE);path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source',path/'source')
    sources_text={'mo_depo.f90':deposition_source(),'mo_velocity.f90':velocity_source()}
    for name,s in sources_text.items():(path/'source/src'/name).write_text(s)
    driver=path/'mcsnow_joint_driver.f90';driver.write_text(driver_source())
    for name in native.TABLE_NAMES:shutil.copyfile(model.BASELINE/name,path/name)
    patch=path/'joint.patch';patch.write_text(''.join(''.join(difflib.unified_diff(
        ((column.CACHE/'build/source/src'/name).read_text()).splitlines(True),s.splitlines(True),
        fromfile='column_density/'+name,tofile='joint/'+name)) for name,s in sources_text.items()))
    compiler=baseline['compiler_path'];sources=native.selected_sources(path/'source/src',driver,compiler)
    objects=[];commands=[]
    with (path/'build.log').open('w') as log:
        for i,src in enumerate(sources):
            obj=path/f'{i:02d}_{src.stem}.o';objects.append(str(obj))
            cmd=[compiler,*native.FLAGS,'-I'+str(path/'source/src'),'-I'+str(path),'-J'+str(path),'-c',str(src),'-o',str(obj)]
            commands.append(cmd);subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler,*native.FLAGS,*objects,'-o',str(path/'tuner')];commands.append(cmd)
        subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    r=dict(schema='terluna.atmosphere.mcsnow-joint-growth-build/1',
        parent_build_sha256=native.sha(column.CACHE/'build/build.json'),
        baseline_build_sha256=native.sha(model.BASELINE/'build.json'),source_identity=baseline['source_identity'],
        compiler_path=compiler,compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],
        flags=native.FLAGS,commands=commands,compiled_sources={str(p):native.sha(p) for p in sources},
        files={p.name:native.sha(p) for p in [path/'tuner',patch,path/'build.log',*[path/n for n in native.TABLE_NAMES]]})
    model.save(path/'build.json',r);return r


def verify_build():
    column.verify_build();path=CACHE/'build';r=json.loads((path/'build.json').read_text())
    if r['parent_build_sha256']!=native.sha(column.CACHE/'build/build.json') or r['baseline_build_sha256']!=native.sha(model.BASELINE/'build.json'):
        raise ValueError('Parent build changed')
    for src,digest in r['compiled_sources'].items():
        if native.sha(src)!=digest:raise ValueError('Compiled source changed: '+src)
    for name,digest in r['files'].items():
        if native.sha(path/name)!=digest:raise ValueError('Build artifact changed: '+name)
    for name,s in [('source/src/mo_depo.f90',deposition_source()),('source/src/mo_velocity.f90',velocity_source()),('mcsnow_joint_driver.f90',driver_source())]:
        if (path/name).read_text()!=s:raise ValueError('Declared source changed: '+name)
    return r


def batch(entries,x_setup,dt=1/128):
    return tuning.Batch(entries,dt,cache=CACHE,decoder=lambda x:decode(x,x_setup),
        experiment=dict(family='joint_growth',setup=x_setup,transport_radius_max_um=150.,transport_early_bounds=[.25,4.]))


def entries_for(temperature,ages=None,all_starts=True):
    entries=copy.deepcopy([e for e in tuning.cases() if e['temperature_C']==temperature and (all_starts or e['initialization']!='volume')])
    fit=baseline_fit(temperature)
    for e in entries:
        e['sample_ages_min']=list(ages) if ages is not None else sorted({r['age_min']
            for label in ('assessment','additional_age_assessment','interval_assessment') if label in fit
            for r in fit[label]['rows'] if r['initialization']==e['initialization']})
    return entries


def probe_entries(temperature):
    entries=[]
    for p in baseline_fit(temperature)['fixed_state_probes']:
        age=p['age_min'];mapping=p['mapping'];target=diagnosis.fit_rates(temperature,age,mapping)
        e=next(copy.deepcopy(e) for e in tuning.cases() if e['temperature_C']==temperature and e['initialization']==mapping)
        e['initial_age_min']=age;e['sample_ages_min']=[age,age+1/32/60]
        e['configuration'].update(initial_mass_kg=target['mass_kg'],initial_equatorial_m=target['model_a_m'],initial_polar_m=target['model_c_m'])
        entries.append(e)
    return entries


def probe_values(raw,temperature):
    if raw['error']:raise ValueError(raw['error'])
    result=[];rows=np.asarray(raw['rows'])
    for j,p in enumerate(baseline_fit(temperature)['fixed_state_probes']):
        start,end=rows[2*j:2*j+2,:3];duration=1/32
        rates=np.r_[(end[0]-start[0])/duration,np.log(end[1:]/start[1:])/duration]
        volumes=math.pi/6*np.array([start[1]**2*start[2],end[1]**2*end[2]])
        result.append(dict(age_min=p['age_min'],mapping=p['mapping'],rates=rates.tolist(),target_rates=p['target_rates'],
            ratio=(rates/p['target_rates']).tolist(),effective_axis_ratio=float(rates[2]/rates[1]),
            target_effective_axis_ratio=p['target_effective_axis_ratio'],deposition_density_kg_m3=float((end[0]-start[0])/(volumes[1]-volumes[0])),
            target_deposition_density_kg_m3=p['target_deposition_density_kg_m3']))
    return result


def evaluate_fit(x,temperature,setup,dt=1/128):
    build();prior=baseline_fit(temperature);out={}
    with batch(entries_for(temperature),setup,dt) as b:r=tuning.assessment(b,x)
    if 'error' in r:raise ValueError('Trajectory failure: '+r['error'])
    index={(p['initialization'],p['age_min']):p for p in r['rows']}
    for label in ('assessment','additional_age_assessment','interval_assessment'):
        if label not in prior:continue
        out[label]=dict(r,rows=[dict(index[(p['initialization'],p['age_min'])],group=p['group']) for p in prior[label]['rows']])
        out[label]['scores']=tuning.score_rows(out[label]['rows'])
    with batch(probe_entries(temperature),setup,dt) as b:out['fixed_state_probes']=probe_values(b.evaluate(x),temperature)
    seeds=[]
    for old in prior['seed_shape_sensitivity']:
        phi=old['initial_aspect_ratio'];e=entries_for(temperature,ages=[r['age_min'] for r in old['assessment']['rows']])[0]
        if e['initialization']!='sphere':raise ValueError('Expected spherical base seed')
        diameter=e['configuration']['initial_equatorial_m']
        e['configuration'].update(initial_equatorial_m=diameter*phi**(-1/3),initial_polar_m=diameter*phi**(2/3))
        with batch([e],setup,dt) as b:a=tuning.assessment(b,x)
        if 'error' in a:raise ValueError(f'Seed {phi} failed: '+a['error'])
        seeds.append(dict(initial_aspect_ratio=phi,initial_volume_equivalent_diameter_m=diameter,assessment=a))
    out.update(seed_shape_sensitivity=seeds,coordinates=np.asarray(x).tolist(),physical_controls=decode(x,setup).tolist(),setup=setup)
    return out


def local_probe_fit(temperature,setup,initial,label,max_nfev=250):
    """Fast screen of rate/state closure before expensive trajectories."""
    build();path=CACHE/'probe_fits'/label;path.mkdir(parents=True,exist_ok=True)
    spec=dict(temperature_C=temperature,setup=setup,initial=list(initial),max_nfev=max_nfev,
              objective='Pointwise excess above min(10%, baseline error), for each rate, effective axis ratio and deposition density.',
              build_sha256=native.sha(CACHE/'build/build.json'))
    if (path/'specification.json').exists() and json.loads((path/'specification.json').read_text())!=spec:
        raise ValueError('Probe fit identity changed')
    model.save(path/'specification.json',spec)
    if (path/'result.json').exists():return json.loads((path/'result.json').read_text())
    priors=baseline_fit(temperature)['fixed_state_probes']
    def matrix(probes):
        return np.array([[*p['rates'],p['effective_axis_ratio'],p['deposition_density_kg_m3']] for p in probes])
    target=np.array([[*p['target_rates'],p['target_effective_axis_ratio'],p['target_deposition_density_kg_m3']] for p in priors])
    allowed=np.minimum(.1,np.abs(matrix(priors)/target-1))
    initial=np.array(initial);free=[0,1,3,5,6,7,8,9,11,12]
    def expand(z):
        x=initial.copy();x[free]=z;return x
    with batch(probe_entries(temperature),setup) as b:
        def residual(z):
            raw=b.evaluate(expand(z))
            if raw['error']:return np.full(target.size,100.)
            relative=matrix(probe_values(raw,temperature))/target-1
            return (np.sign(relative)*np.maximum(abs(relative)-allowed,0)).ravel()
        started=time.monotonic()
        fit=least_squares(residual,np.clip(initial[free],1e-10,1-1e-10),bounds=(np.zeros(len(free)),np.ones(len(free))),
            diff_step=.0001,x_scale='jac',max_nfev=max_nfev,ftol=1e-11,xtol=1e-10,gtol=1e-10)
        x=expand(fit.x);probes=probe_values(b.evaluate(x),temperature)
        relative=matrix(probes)/target-1
        result=dict(label=label,specification=spec,coordinates=x.tolist(),physical_controls=decode(x,setup).tolist(),
            fixed_state_probes=probes,maximum_absolute_relative_error=float(np.max(abs(relative))),
            maximum_excess=float(np.max(np.maximum(abs(relative)-allowed,0))),regressions=int(np.sum(abs(relative)>abs(matrix(priors)/target-1)+acceptance.ROUNDING)),
            above_ten_percent=int(np.sum(abs(relative)>.1)),success=bool(fit.success),message=str(fit.message),
            evaluations=len(b.cache),elapsed_s=time.monotonic()-started,evaluation_directory=str(b.path))
    model.save(path/'result.json',result)
    print(label,{k:result[k] for k in ('maximum_absolute_relative_error','maximum_excess','regressions','evaluations','elapsed_s')},flush=True)
    return result


def initial_state_coordinates(temperature):
    x=baseline_coordinates(temperature)
    if temperature==-5.3:
        early_density=2/(2*.36+1.95);late_density=2/(2*.36+1.16)
        transition=diagnosis.fit_state(temperature,6,'axes')
        igf=1.513;early_axis=(1.95/.36)/igf;late_axis=(1.16/.36)/igf
    else:
        early_density=1.63/(2*.95+.52);late_density=2.37/(2*1.02+.52)
        transition=diagnosis.fit_state(temperature,4,'axes')
        igf=.279;early_axis=(.52/.95)/igf;late_axis=(.52/1.02)/igf
        x[10]=20/80
    log=column.limits.early_coordinate
    x[1]=log(late_axis);x[2]=0.;x[3]=log(late_density);x[7]=log(early_axis)
    radius=(3*transition['volume_m3']/(4*math.pi))**(1/3)*1e6
    x[6]=radius/150;x[8]=radius/300;x[9]=log(early_density/late_density)
    x[11]=50/150;x[12]=.5
    return x


def screen():
    """Reproduce the four full controls; preserve completed execution records."""
    build();results=[]
    setups=[('neutral',NEUTRAL),('prism_drag',dict(NEUTRAL,original_best=1.,prism_width_fraction=3/math.pi)),
            ('normal_compact',dict(NEUTRAL,normal_compact=True)),('state_probe_fit',STATE)]
    for label,setup in setups:
        path=CACHE/(label+'_screen.json');fits=[];started=time.monotonic()
        for temp in (-5.3,-14.4):
            x=baseline_coordinates(temp)
            if label=='state_probe_fit':
                x=json.loads((CACHE/'probe_fits'/('state_'+diagnosis.key(temp))/'result.json').read_text())['coordinates']
            fits.append(evaluate_fit(x,temp,setup))
        decision=acceptance.compare(acceptance.inventory(*fits))
        result=dict(label=label,setup=setup,warm=fits[0],cold=fits[1],acceptance=decision,
            elapsed_s=time.monotonic()-started,build_sha256=native.sha(CACHE/'build/build.json'),baseline_sha256=acceptance.BASELINE_SHA256)
        if path.exists():
            previous=json.loads(path.read_text())
            if any(previous[k]!=result[k] for k in result if k!='elapsed_s'):
                raise ValueError('Completed screen does not reproduce: '+label)
            result=previous
        else:model.save(path,result)
        results.append(result)
        print(label,{k:decision[k] for k in ('accepted','above_ten_percent','regression_count')},flush=True)
    return results


def candidate_records():
    records=[]
    for label in ('neutral','prism_drag','normal_compact','state_probe_fit'):
        path=CACHE/(label+'_screen.json');r=json.loads(path.read_text())
        records.append(dict(label=label,warm=r['warm'],cold=r['cold'],source_path=str(path),source_sha256=native.sha(path)))
    for path in sorted((CACHE/'trajectory_fits').glob('*/result.json')):
        r=json.loads(path.read_text());temp=r['specification']['temperature_C']
        warm=r['full_fit'] if temp==-5.3 else baseline_fit(-5.3)
        cold=r['full_fit'] if temp==-14.4 else baseline_fit(-14.4)
        records.append(dict(label=r['label'],warm=warm,cold=cold,source_path=str(path),source_sha256=native.sha(path),
            fitting_specification=r['specification'],optimizer={k:r[k] for k in ('objective','message','optimizer_success','calls','elapsed_s')}))
    return records


def refinement_checks():
    """Measure timestep effects separately from the physical-error contract."""
    checks=[]
    for r in candidate_records():
        if r['label']=='neutral':continue
        temperatures=(-5.3,-14.4)
        for temp in temperatures:
            fit=r['warm'] if temp==-5.3 else r['cold']
            if 'setup' not in fit:continue
            fine=evaluate_fit(fit['coordinates'],temp,fit['setup'],dt=1/256)
            other=baseline_fit(-14.4 if temp==-5.3 else -5.3)
            old=acceptance.inventory(fit,other) if temp==-5.3 else acceptance.inventory(other,fit)
            new=acceptance.inventory(fine,other) if temp==-5.3 else acceptance.inventory(other,fine)
            old_index={p['id']:p for p in old['points']}
            change=max(abs(p['value']/old_index[p['id']]['value']-1) for p in new['points'])
            record=dict(label=r['label'],temperature_C=temp,timesteps_s=[1/128,1/256],
                maximum_relative_change=change,tolerance=.01,passed=change<.01,
                evidence='Numerical refinement; not a physical acceptance decision.')
            checks.append(record);model.save(CACHE/'refinement_checks.json',checks)
            print('refinement',record,flush=True)
    return checks


def assemble():
    verify_build();baseline=acceptance.baseline_inventory();base_index={p['id']:p for p in baseline['points']}
    candidates=[];checks=[]
    for r in candidate_records():
        inv=acceptance.inventory(r['warm'],r['cold']);decision=acceptance.compare(inv,baseline)
        if r['label']=='neutral':
            difference=max(abs(p['value']/base_index[p['id']]['value']-1) for p in inv['points'])
            checks.append(dict(name='neutral_full_contract_parity',points=len(inv['points']),
                maximum_relative_difference=difference,tolerance=1e-12,passed=difference<1e-12))
        entry={k:v for k,v in r.items() if k not in ('warm','cold')}
        entry['decision']={k:v for k,v in decision.items() if k!='rows'}
        entry['temperatures']=[]
        for temp,fit in ((-5.3,r['warm']),(-14.4,r['cold'])):
            entry['temperatures'].append(dict(temperature_C=temp,coordinates=fit['coordinates'],
                coordinate_family='joint_growth_13' if len(fit['coordinates'])==13 else 'frozen_baseline_original',
                physical_controls=dict(zip(('transport_multiplier','late_log_axis_multiplier','directional_ventilation_exponent',
                    'late_density_factor','branch_radius_multiplier','mean_ventilation_multiplier','axis_transition_radius_um',
                    'compact_maximum_dimension_um','parent_density_flag','early_log_axis_multiplier','density_transition_radius_um',
                    'early_density_ratio','relative_density_flag','equal_normal_compact_flag','transition_power',
                    'transport_transition_radius_um','early_transport_ratio','original_best_weight','span_width_weight',
                    'prism_width_fraction'),fit['physical_controls'],strict=True)) if len(fit['coordinates'])==13 else None,
                setup=fit.get('setup',dict(NEUTRAL)),assessment=fit['assessment'],
                fixed_state_probes=fit['fixed_state_probes'],seed_shape_sensitivity=fit['seed_shape_sensitivity']))
        candidates.append(entry)
    fixed=json.loads(drag.PRODUCT.read_text())
    prism=next(r for r in candidate_records() if r['label']=='prism_drag')
    for temp,fit in ((-5.3,prism['warm']),(-14.4,prism['cold'])):
        predicted=next(r['values'][3] for r in fit['assessment']['rows'] if r['initialization']=='axes' and r['age_min']==3.)
        reference=next(r['speeds']['paper_prism_mean'] for r in fixed['rows'] if r['start']=='reference_axes' and r['age_min']==3. and r['temperature_C']==temp)
        difference=abs(predicted/reference-1)
        checks.append(dict(name='independent_python_prism_drag_'+diagnosis.key(temp),
            maximum_relative_difference=difference,tolerance=1e-11,passed=difference<1e-11))
    checks.extend(json.loads((CACHE/'refinement_checks.json').read_text()))
    evaluations=failed=0;raw_manifest={}
    for root in (drag.CACHE/'growth',CACHE):
        for path in sorted(root.glob('evaluations/*/evaluations.jsonl')):
            records=[json.loads(line) for line in path.read_text().splitlines()]
            evaluations+=len(records);failed+=sum(r['error'] is not None for r in records)
            raw_manifest[str(path)]=native.sha(path)
        for path in sorted(root.glob('probe_fits/*/result.json')):raw_manifest[str(path)]=native.sha(path)
        for path in (root/'build/build.json',root/'producer_initial.py'):
            if path.exists():raw_manifest[str(path)]=native.sha(path)
    for r in candidates:raw_manifest[r['source_path']]=r['source_sha256']
    sources=json.loads((drag.CACHE/'sources/retrieval.json').read_text())
    for src in sources:
        if src['status']=='retrieved' and native.sha(drag.CACHE/'sources'/src['name'])!=src['sha256']:
            raise ValueError('Retrieved source changed: '+src['name'])
    return dict(schema='terluna.atmosphere.mcsnow-joint-growth/1',
        evidence='Diagnostic model-family tests and constrained calibration against published fitted histories. No candidate meets every ten-percent and pointwise non-regression requirement. No replacement model or lunar transfer is adopted.',
        acceptance_product_sha256=native.sha(acceptance.PRODUCT),baseline_sha256=acceptance.BASELINE_SHA256,
        fixed_drag_product_sha256=native.sha(drag.PRODUCT),build_record_sha256=native.sha(CACHE/'build/build.json'),
        candidates=candidates,accepted_candidates=[c['label'] for c in candidates if c['decision']['accepted']],
        selected_model=None,checks=checks,evaluations=evaluations,failed_evaluations=failed,
        raw_manifest=raw_manifest,sources=sources,
        source_reading_scope='Westbrook small-crystal fall-speed preprint: sections 2 and 2.2, including viscous-limit and orientation qualifications. TRAIL part 2: sections 2 and 6 on geometric length scales and Reynolds-number dependence. Welss publisher retrieval failed (403); the original Best-number alternative is identified by the pinned native source comment, not a newly verified transcription of the original Boehm paper. The previously cached westbrook2008.pdf is the separate capacitance paper and was not substituted for the fall-speed paper.',
        closure_interpretation=dict(relative_density='Added-ice density equals current bulk density times a fitted, size-dependent ratio; solid ice is the upper bound. A diagnostic morphology closure, not measured porosity.',
            normal_compact='Equal normal face speed during the compact stage, versus the preceding equal relative-axis growth. This preserves the spherical trajectory but changes nonspherical seeds.',
            prism_drag='Original-documented q denominator plus the side projection of a horizontal regular hexagonal prism averaged over azimuth. This is a geometry sensitivity, not an established correction to McSnow.',
            transport='Optional early/late capacitance multiplier; same multiplier in vapor and heat transfer. Temperature-specific diagnostic bounds, not molecular attachment coefficients.'),
        producers={str(p.relative_to(model.ROOT)):native.sha(p) for p in (Path(__file__),HERE/'mcsnow_joint_search.py',HERE/'mcsnow_joint_acceptance.py',HERE/'mcsnow_joint_drag.py',HERE/'mcsnow_joint_drag_driver.f90')},
        restrictions=['All measured-fit comparisons, diagnostic supplied-state rates, and seed sensitivities retain distinct evidence labels.',
            'Searches are local and bounded; failure does not prove that every possible new growth law is infeasible.',
            'Optimizer convergence and smaller aggregate loss do not establish acceptance. Every point must pass.',
            'Only the two current calibrated temperatures are covered; existing lunar results remain conditional.'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('build','probe_fit','screen','refine','assemble'))
    args=parser.parse_args()
    if args.action=='build':build()
    elif args.action=='screen':screen()
    elif args.action=='refine':refinement_checks()
    elif args.action=='assemble':model.save(PRODUCT,assemble())
    else:
        for temp in (-5.3,-14.4):
            local_probe_fit(temp,STATE,initial_state_coordinates(temp),'state_'+diagnosis.key(temp))
