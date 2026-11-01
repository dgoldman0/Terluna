"""Isolate McSnow mean and directional ventilation with explicit source patches.

These are diagnostic interventions, not an adopted ice-growth parameterization.
Mode zero preserves the selected 2026 source expressions for parity checking.
"""
from __future__ import annotations
import argparse
import difflib
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time
from atmosphere.snow import mcsnow_crystal_growth as native

HERE=Path(__file__).resolve().parent
DRIVER=HERE/'mcsnow_ventilation_driver.f90'
TARGET='src/mo_depo.f90'
NATIVE_SHA='4c13bf06ba5010229a32097bb974b051cc9060aa57457b7b0bd73a7ac04491c0'
MEAN_MODES={0:'native mean ventilation',1:'remove only extra oblate shape term',2:'remove both extra shape terms; retain equivalent-size base',3:'unit mean ventilation; historical sensitivity'}
DIRECTIONAL_MODES={0:'native volume^(1/3) denominator',1:'unit directional factor ratio',2:'physical equivalent-volume sphere radius; source/paper sensitivity'}
TRACE_FIELDS=('mass_kg','equatorial_dimension_m','polar_dimension_m','aspect_ratio','capacitance_m',
    'reynolds_number','diffusivity_raw_m2_s','conductivity_raw_W_m_K','diffusivity_effective_m2_s',
    'conductivity_effective_W_m_K','ventilation_vapor','ventilation_heat','inherent_growth_factor',
    'effective_growth_factor','deposition_density_kg_m3','branch_criterion_met','terminal_speed_m_s',
    'X_vapor','X_heat','directional_radius_m','deposition_rate_kg_s','heat_resistance_ratio')


def producer():
    return {'native_harness':native.producer(),Path(__file__).name:native.sha(__file__),DRIVER.name:native.sha(DRIVER)}


def replace_once(body,old,new):
    count=body.count(old)
    if count!=1:raise ValueError(f'Expected one patch anchor, got {count}: {old[:80]}')
    return body.replace(old,new,1)


def patch_source(body):
    """Preserve native arithmetic in mode zero; add actual-source diagnostics."""
    if hashlib.sha256(body.encode()).hexdigest()!=NATIVE_SHA:raise ValueError('Unrecognized 2026 deposition source')
    body=replace_once(body,'  INTEGER, PARAMETER :: subcycles = 1', '''  INTEGER, PARAMETER :: subcycles = 1
  INTEGER, PUBLIC :: control_mean_mode=0, control_directional_mode=0
  REAL(wp), PUBLIC :: growth_trace(22)=0._wp''')
    start=body.index('SUBROUTINE vapourDiffusion( sp,')
    end=body.index('END SUBROUTINE vapourDiffusion',start)
    part=body[start:end]
    part=replace_once(part,'  C = capacitance( sp, Re )','''  C = capacitance( sp, Re )
  growth_trace=0._wp
  growth_trace(1)=sp%m_i
  growth_trace(4)=sp%phi
  growth_trace(5)=C
  growth_trace(6)=Re
  growth_trace(7)=Dv
  growth_trace(8)=Ka
  growth_trace(17)=vt
  growth_trace(18)=sqrtRe*cubrSc
  growth_trace(19)=sqrtRe*cubrPr''')
    part=part.replace('ventilation_coef_i( sqrtRe*cubrSc, sp%phi )','controlled_ventilation_i( sqrtRe*cubrSc, sp%phi )')
    part=part.replace('ventilation_coef_i( sqrtRe*cubrPr, sp%phi )','controlled_ventilation_i( sqrtRe*cubrPr, sp%phi )')
    part=replace_once(part,'  Hi = (L/(Rv*atmo%T) - 1.0_wp)*(L*Dv*es)/(Ka*Rv*atmo%T**2)','''  growth_trace(9)=Dv
  growth_trace(10)=Ka
  growth_trace(11)=Dv/growth_trace(7)
  growth_trace(12)=Ka/growth_trace(8)
  Hi = (L/(Rv*atmo%T) - 1.0_wp)*(L*Dv*es)/(Ka*Rv*atmo%T**2)
  growth_trace(22)=Hi''')
    part=replace_once(part,'  ! sum heat for possible melting at surface temperature T_3', '  growth_trace(21)=dMdt\n  ! sum heat for possible melting at surface temperature T_3')
    body=body[:start]+part+body[end:]
    start=body.index('SUBROUTINE habit_predict(')
    end=body.index('END SUBROUTINE habit_predict',start)
    part=body[start:end]
    part=replace_once(part,'  CALL spheroid_geometry(sp, r_a, r_c, V_tot)','''  CALL spheroid_geometry(sp, r_a, r_c, V_tot)
  growth_trace(2)=2._wp*r_a
  growth_trace(3)=2._wp*r_c
  growth_trace(13)=IGF''')
    part=replace_once(part,'  dVidt    = dMdt / rho_depo', '''  growth_trace(15)=rho_depo
  growth_trace(16)=MERGE(1._wp,0._wp,branch_crit)
  dVidt    = dMdt / rho_depo''')
    anchor='''  IGF_vent = IGF * ventilation_coef_IGF( ReSc, r_c, V_tot**z1_3, sp%phi ) /                       &
                   ventilation_coef_IGF( ReSc, r_a, V_tot**z1_3, sp%phi ) '''
    part=replace_once(part,anchor,anchor+'''
  growth_trace(20)=V_tot**z1_3
  IF (control_directional_mode==1) THEN
    IGF_vent=IGF
  ELSE IF (control_directional_mode==2) THEN
    growth_trace(20)=(z3_4pi*V_tot)**z1_3
    IGF_vent=IGF*ventilation_coef_IGF(ReSc,r_c,growth_trace(20),sp%phi)/ &
                 ventilation_coef_IGF(ReSc,r_a,growth_trace(20),sp%phi)
  END IF
  growth_trace(14)=IGF_vent''')
    body=body[:start]+part+body[end:]
    function='''! Explicit experimental mean-ventilation controls; native function remains unchanged.
FUNCTION controlled_ventilation_i(X,phi) RESULT(fv)
  REAL(wp), INTENT(in) :: X,phi
  REAL(wp) :: fv,X_equiv
  IF (control_mean_mode==0) THEN
    fv=ventilation_coef_i(X,phi)
    RETURN
  END IF
  IF (control_mean_mode==3) THEN
    fv=1._wp
    RETURN
  END IF
  IF (phi>=1._wp) THEN
    X_equiv=X*phi**(-1./3.)
  ELSE
    X_equiv=X*phi**(1./6.)
  END IF
  IF (X_equiv<1._wp) THEN
    fv=1._wp+0.14_wp*X_equiv**2
  ELSE
    fv=0.86_wp+0.28_wp*X_equiv
  END IF
  IF (control_mean_mode==1.AND.phi>1._wp) fv=fv+0.0280_wp*phi*X_equiv
END FUNCTION controlled_ventilation_i

'''
    body=replace_once(body,'PURE ELEMENTAL FUNCTION ventilation_coef_i( X, phi ) result ( fv )',function+'PURE ELEMENTAL FUNCTION ventilation_coef_i( X, phi ) result ( fv )')
    return body


def build(output,native_build_dir,compiler='gfortran'):
    output=Path(output).resolve();native_build_dir=Path(native_build_dir).resolve()
    if output.exists():raise FileExistsError(output)
    baseline=native.verify_build(native_build_dir)
    if baseline['source_identity']['version']!='2026':raise ValueError('Control patch is specific to 2026 source')
    own=producer();original=(native_build_dir/'source'/TARGET).read_text();patched=patch_source(original)
    output.mkdir(parents=True);source=output/'source';shutil.copytree(native_build_dir/'source',source)
    (source/TARGET).write_text(patched)
    diff=''.join(difflib.unified_diff(original.splitlines(keepends=True),patched.splitlines(keepends=True),fromfile=TARGET+'.native',tofile=TARGET+'.controls'))
    (output/'physics.patch').write_text(diff)
    for name in native.TABLE_NAMES:shutil.copyfile(native_build_dir/name,output/name)
    shutil.copyfile(DRIVER,output/DRIVER.name)
    compiler_path=str(Path(shutil.which(compiler)).resolve())
    order=native.selected_sources(source/'src',output/DRIVER.name,compiler_path)
    commands=[];objects=[];start=time.monotonic()
    with (output/'build.log').open('w') as log:
        for i,p in enumerate(order):
            obj=output/f'{i:02d}_{p.stem}.o';objects.append(str(obj))
            cmd=[compiler_path,*native.FLAGS,'-I'+str(source/'src'),'-I'+str(output),'-J'+str(output),'-c',str(p),'-o',str(obj)]
            commands.append(cmd);subprocess.run(cmd,cwd=output,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler_path,*native.FLAGS,*objects,'-o',str(output/'crystal')]
        commands.append(cmd);subprocess.run(cmd,cwd=output,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    if own!=producer() or native.verify_build(native_build_dir)!=baseline:raise RuntimeError('Producer or baseline changed during build')
    product=dict(schema='terluna.atmosphere.mcsnow-ventilation-control-build/1',producer=own,
        native_build_record_sha256=native.sha(native_build_dir/'build.json'),native_source_identity=baseline['source_identity'],
        source_files={str(p.relative_to(source)):native.sha(p) for p in sorted(source.rglob('*')) if p.is_file()},
        native_source_files=baseline['source_files'],compiled_source_files=[str(p.relative_to(output)) for p in order],
        table_source_revision=baseline['table_source_revision'],table_identity=baseline['table_identity'],
        physics_patches=[dict(path=TARGET,native_sha256=NATIVE_SHA,patched_sha256=native.sha(source/TARGET),patch_sha256=native.sha(output/'physics.patch'))],
        compiler_path=compiler_path,compiler_sha256=native.sha(compiler_path),
        compiler_version=subprocess.check_output([compiler_path,'--version'],text=True).splitlines()[0],flags=native.FLAGS,
        commands=commands,executable_sha256=native.sha(output/'crystal'),driver_filename=DRIVER.name,
        driver_sha256=native.sha(output/DRIVER.name),build_log_sha256=native.sha(output/'build.log'),elapsed_s=time.monotonic()-start,
        mean_modes=MEAN_MODES,directional_modes=DIRECTIONAL_MODES,
        evidence='Explicit source interventions for causal diagnosis; no intervention is calibrated or adopted.')
    (output/'build.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n');return product


def verify_build(build_dir):
    build_dir=Path(build_dir);record=json.loads((build_dir/'build.json').read_text())
    for filename,key in [('crystal','executable_sha256'),(record['driver_filename'],'driver_sha256'),('build.log','build_log_sha256')]:
        if native.sha(build_dir/filename)!=record[key]:raise ValueError('Control build changed: '+filename)
    for path,digest in record['source_files'].items():
        if native.sha(build_dir/'source'/path)!=digest:raise ValueError('Control source changed: '+path)
    for name,item in record['table_identity'].items():native.pinned.verify((build_dir/name).read_bytes(),item)
    if native.sha(build_dir/'physics.patch')!=record['physics_patches'][0]['patch_sha256']:raise ValueError('Control patch changed')
    return record


def validate_modes(mean_mode,directional_mode):
    if mean_mode not in MEAN_MODES or directional_mode not in DIRECTIONAL_MODES:raise ValueError('Unknown ventilation intervention')


def run_configuration(build_dir,output,config,*,metadata=None):
    """Run arbitrary pure-ice states, including one-step observed-state replays.

    Config uses native namelist names plus mean_mode/directional_mode and optional
    trace_every_step. Raw snapshot and diagnostic rows remain untampered by fits.
    """
    config=dict(config);config.setdefault('mean_mode',0);config.setdefault('directional_mode',0);config.setdefault('trace_every_step',False)
    validate_modes(config['mean_mode'],config['directional_mode'])
    for key in ('dt_s','duration_s','output_interval_s'):
        if not math.isfinite(config[key]) or config[key]<=0:raise ValueError('Invalid time: '+key)
    for key in ('duration_s','output_interval_s'):
        if not math.isclose(config[key]/config['dt_s'],round(config[key]/config['dt_s']),rel_tol=0,abs_tol=1e-8):raise ValueError('Time does not divide timestep: '+key)
    build_dir=Path(build_dir).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    own=producer();identity=verify_build(build_dir)
    if identity['producer']!=own:raise ValueError('Control producer differs')
    if native.sha(native.scale.lunar.CACHE/native.scale.SOURCE_PDF)!=native.lab.SOURCE_SHA:raise ValueError('Empirical source identity differs')
    output.mkdir(parents=True);(output/'crystal.nml').write_text(native.namelist(config))
    for name in native.TABLE_NAMES:shutil.copyfile(build_dir/name,output/name)
    inputs={p.name:native.sha(p) for p in output.iterdir() if p.is_file()};start_time=time.monotonic()
    with (output/'run.log').open('w') as log:
        completed=subprocess.run([str(build_dir/'crystal'),'crystal.nml'],cwd=output,env=native.ENV,stdout=log,stderr=subprocess.STDOUT)
    rows,errors=native.read_csv(output/'trajectory.csv',(*native.QUANTITIES,'elapsed_s','ice_mass_balance_residual_kg'))
    trace,trace_errors=native.read_csv(output/'growth_trace.csv',('evaluation_elapsed_s','state_elapsed_s',*TRACE_FIELDS));errors+=trace_errors
    air,air_errors=native.read_csv(output/'air.csv',('temperature_K','pressure_Pa','vapor_pressure_Pa'));errors+=air_errors
    steps,step_errors=native.read_csv(output/'step_diagnostic.csv',('maximum_step_mass_fraction','steps')) if completed.returncode==0 else ([],[])
    errors+=step_errors
    output_steps=round(config['output_interval_s']/config['dt_s'])
    nsteps=round(config['duration_s']/config['dt_s'])
    saved_steps=[0,*range(output_steps,nsteps+1,output_steps)]
    if saved_steps[-1]!=nsteps:saved_steps.append(nsteps)
    trace_steps=list(range(1,nsteps+1)) if config['trace_every_step'] else sorted(set([1,*saved_steps[1:]]))
    if completed.returncode==0 and len(rows)!=len(saved_steps):errors.append('Successful exit with incomplete raw trajectory')
    if completed.returncode==0 and len(trace)!=len(trace_steps):errors.append('Successful exit with incomplete growth trace')
    for row,step in zip(rows,saved_steps):
        if not math.isclose(row['elapsed_s'],step*config['dt_s'],rel_tol=0,abs_tol=1e-8):errors.append('Irregular raw trajectory time')
        if any(row[k]<=0 for k in native.QUANTITIES):errors.append('Nonpositive raw crystal state')
    for row,step in zip(trace,trace_steps):
        if not math.isclose(row['state_elapsed_s'],step*config['dt_s'],rel_tol=0,abs_tol=1e-8):errors.append('Irregular source trace time')
        if not math.isclose(row['evaluation_elapsed_s'],(step-1)*config['dt_s'],rel_tol=0,abs_tol=1e-8):errors.append('Incorrect source evaluation time')
        if any(row[k]<=0 for k in TRACE_FIELDS if k not in ('branch_criterion_met','reynolds_number','terminal_speed_m_s','X_vapor','X_heat')):
            errors.append('Nonpositive source growth diagnostic')
        if row['branch_criterion_met'] not in (0.,1.):errors.append('Invalid source branching flag')
    if len(air)!=1:errors.append('Expected one air state')
    for k in ('temperature_K','pressure_Pa','vapor_pressure_Pa'):
        if air and not math.isclose(air[0][k],config[k],rel_tol=1e-13,abs_tol=1e-10):errors.append('Air differs: '+k)
    if not rows or not trace:errors.append('No usable trajectory or actual growth trace')
    if rows:
        initial_keys={'initial_mass_kg':'mass_kg','initial_equatorial_m':'equatorial_dimension_m','initial_polar_m':'polar_dimension_m'}
        for k,q in initial_keys.items():
            if not math.isclose(rows[0][q],config[k],rel_tol=1e-12):errors.append('Initial mapping differs: '+q)
        if any(abs(r['ice_mass_balance_residual_kg'])/r['mass_kg']>1e-10 for r in rows):errors.append('Ice accounting failed')
    if own!=producer() or verify_build(build_dir)!=identity:raise RuntimeError('Control identities changed during run')
    if inputs!={n:native.sha(output/n) for n in inputs}:raise RuntimeError('Control input changed during run')
    product=dict(schema='terluna.atmosphere.mcsnow-ventilation-control-run/1',producer=own,
        build_record_sha256=native.sha(build_dir/'build.json'),configuration=config,metadata=metadata or {},
        complete=completed.returncode==0 and not errors,returncode=completed.returncode,output_validation_errors=errors,
        raw_snapshots=rows,growth_trace=trace,air=air[0] if air else {},step_diagnostic=steps[0] if steps else None,
        inputs=inputs,outputs={p.name:native.sha(p) for p in output.iterdir() if p.is_file() and p.name not in inputs},
        elapsed_s=time.monotonic()-start_time,source_pdf_sha256=native.lab.SOURCE_SHA,
        evidence='Declared diagnostic source intervention; prescribed Earth air with actual pre-step source terms.',
        limits=['Mean/directional modes are sensitivities, not accepted replacements or calibration.',
                'No water reservoir depletion, liquid collection, collisions, motion or atmospheric heat feedback.',
                'Trace evaluation time is the start of the step whose end is state_elapsed_s.'])
    (output/'run.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n');return product


def run(build_dir,output,temperature_C,dt,growth_mode,initialization='fit',mean_mode=0,directional_mode=0,trace_every_step=False):
    identity=verify_build(build_dir)
    body=(Path(build_dir)/'source/src/mo_atmo_types.f90').read_text()
    match=re.search(r'\brhoi\s*=\s*([0-9.]+)_wp',body,re.I)
    if not match:raise ValueError('Native ice density missing')
    density=float(match[1]);start=native.initial(temperature_C,initialization,density)
    config=native.configuration(temperature_C,dt,growth_mode,initialization,ice_density=density)
    config.update(mean_mode=mean_mode,directional_mode=directional_mode,trace_every_step=trace_every_step)
    product=run_configuration(build_dir,output,config,metadata=dict(temperature_C=temperature_C,initialization=initialization,initial=start))
    fit_analysis=native.analyze_output(output,config,start,temperature_C,growth_mode,product['returncode'])
    validation_errors=product['output_validation_errors']+fit_analysis['output_validation_errors']
    product.update(fit_analysis)
    product['output_validation_errors']=list(dict.fromkeys(validation_errors))
    product['complete']=product['complete'] and not product['output_validation_errors']
    product.update(temperature_C=temperature_C,initialization=initialization,initial=start,
                   mean_mode=MEAN_MODES[mean_mode],directional_mode=DIRECTIONAL_MODES[directional_mode])
    (Path(output)/'run.json').write_text(json.dumps(product,indent=2,allow_nan=False)+'\n');return product


def parity(native_output,control_output):
    """Compare unchanged-mode raw bytes, including matched native failures.

    A failed execution has no final step file. Such absence is reported rather
    than counted as a compared file, and is accepted only when both runs failed.
    """
    native_output=Path(native_output);control_output=Path(control_output)
    a_run=json.loads((native_output/'run.json').read_text())
    b_run=json.loads((control_output/'run.json').read_text())
    status_equal=(a_run['returncode']==b_run['returncode'] and a_run['complete']==b_run['complete'])
    checks={}
    for name in ('trajectory.csv','air.csv','branch_diagnostic.csv','step_diagnostic.csv'):
        a=native_output/name;b=control_output/name
        both_present=a.exists() and b.exists()
        matched_failure_absence=(name=='step_diagnostic.csv' and not a.exists() and not b.exists()
                                 and a_run['returncode']!=0 and b_run['returncode']!=0)
        checks[name]=dict(native_exists=a.exists(),control_exists=b.exists(),
            native_sha256=native.sha(a) if a.exists() else None,control_sha256=native.sha(b) if b.exists() else None,
            equal=a.read_bytes()==b.read_bytes() if both_present else None,
            matched_failure_absence=matched_failure_absence)
    return dict(all_equal=status_equal and all(c['equal'] is True or c['matched_failure_absence'] for c in checks.values()),
                status_equal=status_equal,compared_files=sum(c['equal'] is not None for c in checks.values()),files=checks)


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    b=s.add_parser('build');b.add_argument('--output',required=True,type=Path);b.add_argument('--native-build',required=True,type=Path)
    r=s.add_parser('run');r.add_argument('--build',required=True,type=Path);r.add_argument('--output',required=True,type=Path)
    r.add_argument('--temperature',required=True,type=float);r.add_argument('--dt',required=True,type=float)
    r.add_argument('--growth-mode',required=True,type=int,choices=(1,2));r.add_argument('--initialization',choices=('fit','sphere'),default='fit')
    r.add_argument('--mean-mode',type=int,choices=tuple(MEAN_MODES),default=0);r.add_argument('--directional-mode',type=int,choices=tuple(DIRECTIONAL_MODES),default=0)
    r.add_argument('--trace-every-step',action='store_true');a=p.parse_args()
    if a.command=='build':result=build(a.output,a.native_build);print(result['executable_sha256'])
    else:
        result=run(a.build,a.output,a.temperature,a.dt,a.growth_mode,a.initialization,a.mean_mode,a.directional_mode,a.trace_every_step)
        print(json.dumps({k:result[k] for k in ('complete','returncode','output_validation_errors')},indent=2))
if __name__=='__main__':main()
