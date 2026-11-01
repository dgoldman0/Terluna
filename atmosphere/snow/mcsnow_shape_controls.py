"""Causal shape controls: native uptake/density with an imposed axis ratio.

The fitted-ratio control deliberately uses the comparison fits. It diagnoses
feedback and cannot be counted as an independent validation or adopted law.
"""
from __future__ import annotations

import argparse
import difflib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess

from atmosphere.snow import mcsnow_shape_diagnosis as diagnosis

trace=diagnosis.trace
native=trace.native
model=trace.model
CACHE=diagnosis.CACHE/'shape_controls'
PRODUCT=diagnosis.HERE/'results/mcsnow_shape_controls.json'
MODES={0:'retained improved control',1:'unit directional ventilation; intrinsic growth ratio unchanged',
       2:'effective axis-growth ratio prescribed from the fitted history; native deposition density unchanged'}


def deposition_source():
    original=(model.BASELINE/'source/src/mo_depo.f90').read_text()
    body=trace.source_patch(original,'plate_base')
    body=model.replace_once(body,'  REAL(wp), PUBLIC :: growth_trace(22)=0._wp',
        '  REAL(wp), PUBLIC :: growth_trace(22)=0._wp\n  REAL(wp), PUBLIC :: diagnostic_growth_ratio=-1._wp')
    return model.replace_once(body,'  growth_trace(14)=IGF_vent',
        '  if (diagnostic_growth_ratio>0._wp) IGF_vent=diagnostic_growth_ratio\n  growth_trace(14)=IGF_vent')


def driver_source():
    body=trace.trace_driver()
    body=model.replace_once(body,'  use mo_depo, only: sp_depo, growth_trace',
        '  use mo_depo, only: sp_depo, growth_trace, diagnostic_growth_ratio')
    body=model.replace_once(body,'  integer :: growth_mode, nsteps, output_steps, k, unit',
        '  integer :: growth_mode, nsteps, output_steps, k, unit, shape_control\n'+
        '  real(wp) :: initial_age_s, ratio_before, ratio_after, ratio_join_s, ratio_start_s')
    body=model.replace_once(body,'initial_polar_m, growth_mode, prolate_interpolation, native_saturation',
        'initial_polar_m, &\n    growth_mode, prolate_interpolation, native_saturation, &\n'+
        '    shape_control, initial_age_s, ratio_before, ratio_after, ratio_join_s, ratio_start_s')
    body=model.replace_once(body,'    old_mass=particle%m_i',
        '''    old_mass=particle%m_i
    diagnostic_growth_ratio=-1._wp
    if (shape_control==1) diagnostic_growth_ratio=get_IGF(air%T)
    if (shape_control==2.and.initial_age_s+elapsed>=ratio_start_s) then
      diagnostic_growth_ratio=ratio_before
      if (initial_age_s+elapsed>=ratio_join_s) diagnostic_growth_ratio=ratio_after
    end if''')
    return body


def build():
    path=CACHE/'build'
    baseline=native.verify_build(model.BASELINE)
    if path.exists():
        return verify_build()
    changed=deposition_source(); driver=driver_source()
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE/'source',path/'source')
    (path/'source/src/mo_depo.f90').write_text(changed)
    (path/native.DRIVER.name).write_text(driver)
    for name in native.TABLE_NAMES:
        shutil.copyfile(model.BASELINE/name,path/name)
    patch=''.join(difflib.unified_diff(trace.source_patch((model.BASELINE/'source/src/mo_depo.f90').read_text(),'plate_base').splitlines(True),
        changed.splitlines(True),fromfile='improved/mo_depo.f90',tofile='diagnostic/mo_depo.f90'))
    patch+=''.join(difflib.unified_diff(trace.trace_driver().splitlines(True),driver.splitlines(True),
        fromfile='trace/driver.f90',tofile='diagnostic/driver.f90'))
    (path/'diagnostic.patch').write_text(patch)
    compiler=baseline['compiler_path']
    selected=native.selected_sources(path/'source/src',path/native.DRIVER.name,compiler)
    objects=[]; commands=[]
    with (path/'build.log').open('w') as log:
        for i,src in enumerate(selected):
            obj=path/f'{i:02d}_{src.stem}.o'; objects.append(str(obj))
            cmd=[compiler,*native.FLAGS,'-I'+str(path/'source/src'),'-I'+str(path),'-J'+str(path),
                 '-c',str(src),'-o',str(obj)]
            commands.append(cmd)
            subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
        cmd=[compiler,*native.FLAGS,*objects,'-o',str(path/'crystal')]; commands.append(cmd)
        subprocess.run(cmd,cwd=path,env=native.ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    record=dict(schema='terluna.atmosphere.mcsnow-shape-control-build/1',
        baseline_record_sha256=native.sha(model.BASELINE/'build.json'),source_identity=baseline['source_identity'],
        compiled_sources={str(p):native.sha(p) for p in selected},
        compiler_path=compiler,compiler_sha256=native.sha(compiler),compiler_version=baseline['compiler_version'],
        flags=native.FLAGS,commands=commands,
        file_sha256={p.name:native.sha(p) for p in
                     [path/'crystal',path/'build.log',path/'diagnostic.patch',*[path/n for n in native.TABLE_NAMES]]})
    model.save(path/'build.json',record)
    return record


def verify_build():
    path=CACHE/'build'; record=json.loads((path/'build.json').read_text())
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE/'build.json')!=record['baseline_record_sha256']:
        raise ValueError('Baseline record changed')
    for name,digest in record['compiled_sources'].items():
        if native.sha(name)!=digest:raise ValueError('Compiled source changed: '+name)
    for name,digest in record['file_sha256'].items():
        if native.sha(path/name)!=digest:raise ValueError('Build artifact changed: '+name)
    if (path/'source/src/mo_depo.f90').read_text()!=deposition_source() or (path/native.DRIVER.name).read_text()!=driver_source():
        raise ValueError('Diagnostic source changed')
    return record


def cases():
    for entry in diagnosis.cases(diagnosis.ice_density()):
        temperature=entry['temperature_C']
        if temperature not in (-5.3,-14.4):continue
        # Within each segment, the ratio of logarithmic axis rates is c exponent / a exponent.
        before=diagnosis.fit_rates(temperature,3.5,'axes')['effective_growth_ratio']
        after=diagnosis.fit_rates(temperature,10.,'axes')['effective_growth_ratio']
        join=6. if temperature==-5.3 else 4.
        for mode in MODES:
            config=dict(entry['configuration'],shape_control=mode,initial_age_s=entry['initial_age_min']*60,
                        ratio_before=before,ratio_after=after,ratio_join_s=join*60,ratio_start_s=180.)
            yield dict(**{k:v for k,v in entry.items() if k not in ('name','configuration')},
                       name=entry['name']+f'_control{mode}',shape_control=mode,
                       baseline_name=entry['name'],configuration=config)


def execute():
    build()
    for entry in cases():
        path=CACHE/'runs'/entry['name']; config=entry['configuration']
        if path.exists():
            verify_run(entry); continue
        path.mkdir(parents=True)
        (path/'crystal.nml').write_text(native.namelist(config))
        for name in native.TABLE_NAMES:shutil.copyfile(CACHE/'build'/name,path/name)
        with (path/'run.log').open('w') as log:
            proc=subprocess.run([str(CACHE/'build/crystal'),'crystal.nml'],cwd=path,env=native.ENV,
                                stdout=log,stderr=subprocess.STDOUT)
        record=dict(configuration=config,returncode=proc.returncode,
                    build_record_sha256=native.sha(CACHE/'build/build.json'),
                    **trace.read_run(path,config,proc.returncode),
                    file_sha256={p.name:native.sha(p) for p in path.iterdir() if p.is_file()})
        model.save(path/'run.json',record)
        print(entry['name'],'complete' if record['complete'] else 'FAILED',flush=True)


def verify_run(entry):
    path=CACHE/'runs'/entry['name']; record=json.loads((path/'run.json').read_text())
    if record['configuration']!=entry['configuration'] or record['build_record_sha256']!=native.sha(CACHE/'build/build.json'):
        raise ValueError('Control input or build changed')
    for name,digest in record['file_sha256'].items():
        if native.sha(path/name)!=digest:raise ValueError('Control output changed: '+name)
    if (path/'crystal.nml').read_text()!=native.namelist(entry['configuration']):
        raise ValueError('Control namelist changed')
    replay=trace.read_run(path,entry['configuration'],record['returncode'])
    if any(record[k]!=v for k,v in replay.items()):raise ValueError('Control raw replay differs')
    if not record['complete']:raise ValueError('Failed control retained at '+str(path))
    return record


def assemble():
    verify_build(); runs=[]; parity=[]; manifest=[]; source_checks=[]
    rho=diagnosis.ice_density()
    for entry in cases():
        record=verify_run(entry)
        runs.append(diagnosis.trajectory_record(entry,record))
        density_error=0.; ratio_error=0.
        for row in record['trace']:
            intrinsic=row['inherent_growth_factor']
            expected_density=rho/intrinsic if intrinsic>=1 else rho*(intrinsic if row['branch_criterion_met'] else 1.)
            density_error=max(density_error,abs(row['deposition_density_kg_m3']/expected_density-1))
            age=entry['configuration']['initial_age_s']+row['evaluation_elapsed_s']
            expected_ratio=None
            if entry['shape_control']==1:
                expected_ratio=intrinsic
            elif entry['shape_control']==2 and age>=entry['configuration']['ratio_start_s']:
                expected_ratio=entry['configuration']['ratio_before'] if age<entry['configuration']['ratio_join_s'] else entry['configuration']['ratio_after']
            if expected_ratio is not None:
                ratio_error=max(ratio_error,abs(row['effective_growth_factor']/expected_ratio-1))
        if density_error>5e-13 or ratio_error>5e-13:
            raise ValueError('Directional intervention changed the wrong source quantities')
        source_checks.append(dict(name=entry['name'],maximum_density_relative_residual=density_error,
                                  maximum_prescribed_ratio_relative_residual=ratio_error))
        path=CACHE/'runs'/entry['name']
        manifest.append(dict(name=entry['name'],path=str(path.relative_to(model.ROOT)),
                             record_sha256=native.sha(path/'run.json'),file_sha256=record['file_sha256']))
        if entry['shape_control']==0:
            old=trace.verify_run(trace.CACHE/'traces/shape_diagnosis'/entry['baseline_name'])
            parity.append(dict(name=entry['name'],exact_raw_equality=all(
                (path/n).read_bytes()==(trace.CACHE/'traces/shape_diagnosis'/entry['baseline_name']/n).read_bytes()
                for n in ('trajectory.csv','air.csv','growth_trace.csv','branch_diagnostic.csv','step_diagnostic.csv')),
                exact_state_equality=all(record[k]==old[k] for k in ('snapshots','air','trace'))))
    refinements=[]
    for fine in [r for r in runs if r['dt_s']==min(diagnosis.TIMESTEPS)]:
        coarse=next(r for r in runs if r['dt_s']==max(diagnosis.TIMESTEPS) and
                    (r['temperature_C'],r['initialization'],r['shape_control'])==
                    (fine['temperature_C'],fine['initialization'],fine['shape_control']))
        max_change=max(abs(lo['age_matched']['values'][k]/hi['age_matched']['values'][k]-1)
                       for lo,hi in zip(coarse['comparisons'],fine['comparisons']) for k in diagnosis.FIELDS)
        refinements.append(dict(temperature_C=fine['temperature_C'],initialization=fine['initialization'],
                                shape_control=fine['shape_control'],maximum_saved_relative_change=max_change))
    return dict(schema='terluna.atmosphere.mcsnow-shape-controls/1',
        evidence='Causal sensitivity using the comparison fits to prescribe effective axis growth. This is not independent validation or an adopted physical law.',
        producer={str(p.relative_to(model.ROOT)):native.sha(p) for p in (Path(__file__),Path(diagnosis.__file__),Path(trace.__file__),Path(diagnosis.lab.__file__))},
        source_product_sha256=native.sha(diagnosis.PRODUCT),modes={str(k):v for k,v in MODES.items()},trajectories=runs,
        unchanged_control_parity=parity,numerical_refinements=refinements,source_checks=source_checks,manifest=manifest,
        build_record_sha256=native.sha(CACHE/'build/build.json'),
        restrictions=['Only the effective directional growth ratio changes. Native IGF, deposition density, branching rule, drag and uptake calculations remain active.',
                      'The fitted-ratio control starts at three minutes, retaining native early growth for the spherical seed. It uses age and temperature-specific comparison data; it has no extrapolation rule for changing humidity, pressure or temperature.',
                      'Its smaller residuals would be conditional fit agreement, not independently validated accuracy.',
                      'Separately fitted mass and axis histories have no pointwise uncertainty and need not describe one tracked crystal.'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('run','assemble'));args=p.parse_args()
    if args.action=='run':execute()
    else:
        result=assemble();model.save(PRODUCT,result);print(len(result['trajectories']),'shape-control trajectories')


if __name__=='__main__':main()
