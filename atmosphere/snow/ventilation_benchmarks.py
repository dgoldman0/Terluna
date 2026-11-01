"""Factorial Earth growth controls for McSnow mean and directional ventilation.

Every control remains a diagnostic intervention. The matrix retains both primary
crystal histories and both axes, records native failures, and tests time steps.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import re
from atmosphere.snow import mcsnow_ventilation_controls as controls
from atmosphere.snow import mcsnow_growth_benchmarks as reference

native=controls.native
ROOT=native.ROOT
CACHE=ROOT/'research/runs/snow_scale/mcsnow_ventilation'
BASELINE=ROOT/'research/runs/snow_scale/mcsnow_growth'
RESULT=Path(__file__).with_name('results')/'ventilation_benchmarks.json'
TEMPERATURES=reference.TEMPERATURES
COMBINATIONS=((0,0),(1,0),(0,1),(1,1),(2,0),(2,1))
HISTORIES=(('fit',(1.,.2)),('sphere',(.05,.0125)))
QUANTITIES=native.QUANTITIES


def case_name(mean,directional,growth,temperature,initialization,dt):
    return f'm{mean}_d{directional}_igf{growth}_{initialization}_T{abs(temperature):g}_dt{dt:g}'.replace('.','p')


def matrix():
    for mean,directional in COMBINATIONS:
        for growth in (1,2):
            for temperature in TEMPERATURES:
                for initialization,timesteps in HISTORIES:
                    for dt in timesteps:yield mean,directional,growth,temperature,initialization,dt


def _close(a,b,label,tolerance=2e-12):
    if not math.isclose(a,b,rel_tol=tolerance,abs_tol=0.):raise ValueError(label)


def verify_run(path,build_dir):
    """Independently reread preserved source/input/output evidence and geometry."""
    path=Path(path);build_dir=Path(build_dir)
    record=json.loads((path/'run.json').read_text());identity=controls.verify_build(build_dir)
    if record['producer']!=controls.producer() or identity['producer']!=controls.producer():raise ValueError('Control producer differs')
    if record['build_record_sha256']!=native.sha(build_dir/'build.json'):raise ValueError('Control build identity differs')
    for group in ('inputs','outputs'):
        for name,digest in record[group].items():
            if native.sha(path/name)!=digest:raise ValueError('Changed '+group+' file: '+name)
    if record['source_pdf_sha256']!=native.lab.SOURCE_SHA or native.sha(native.scale.lunar.CACHE/native.scale.SOURCE_PDF)!=native.lab.SOURCE_SHA:
        raise ValueError('Empirical source differs')
    config=record['configuration']
    density_match=re.search(r'\brhoi\s*=\s*([0-9.]+)_wp',(build_dir/'source/src/mo_atmo_types.f90').read_text(),re.I)
    if not density_match:raise ValueError('Native density not found')
    density=float(density_match[1])
    regenerated=native.configuration(record['temperature_C'],config['dt_s'],config['growth_mode'],record['initialization'],ice_density=density)
    regenerated.update(mean_mode=config['mean_mode'],directional_mode=config['directional_mode'],trace_every_step=False)
    start=native.initial(record['temperature_C'],record['initialization'],density)
    if config!=regenerated or record['initial']!=start or (path/'crystal.nml').read_text()!=native.namelist(regenerated):raise ValueError('Physics inputs do not reproduce')
    raw,raw_errors=native.read_csv(path/'trajectory.csv',(*QUANTITIES,'elapsed_s','ice_mass_balance_residual_kg'))
    trace,trace_errors=native.read_csv(path/'growth_trace.csv',('evaluation_elapsed_s','state_elapsed_s',*controls.TRACE_FIELDS))
    if raw!=record['raw_snapshots'] or trace!=record['growth_trace']:raise ValueError('Raw control output differs from record')
    analyzed=native.analyze_output(path,config,start,record['temperature_C'],config['growth_mode'],record['returncode'])
    for key,value in analyzed.items():
        if key not in ('complete','output_validation_errors') and record[key]!=value:raise ValueError('Native output analysis differs: '+key)
    for error in raw_errors+trace_errors+analyzed['output_validation_errors']:
        if error not in record['output_validation_errors']:raise ValueError('Output error was discarded: '+error)
    raw_values=trace_values=0
    for row,snapshot in zip(raw,record['snapshots']):
        for key,value in row.items():
            if snapshot[key]!=value:raise ValueError('Raw snapshot differs: '+key)
            raw_values+=1
        volume=math.pi/6*row['equatorial_dimension_m']**2*row['polar_dimension_m']
        _close(volume,row['ice_volume_m3'],'Axis/volume relation differs')
        _close(max(row['equatorial_dimension_m'],row['polar_dimension_m']),row['diameter_m'],'Maximum dimension differs')
        _close(row['polar_dimension_m']/row['equatorial_dimension_m'],row['aspect_ratio'],'Aspect ratio differs')
        _close(row['mass_kg']/volume,row['density_kg_m3'],'Apparent density differs')
    if trace and raw:
        for key in ('mass_kg','equatorial_dimension_m','polar_dimension_m','aspect_ratio'):
            _close(trace[0][key],raw[0][key],'Initial source trace geometry differs: '+key)
    air=record['air']
    for row in trace:
        if not math.isclose(row['state_elapsed_s']-row['evaluation_elapsed_s'],config['dt_s'],abs_tol=1e-9,rel_tol=0.):raise ValueError('Actual-source trace time differs')
        _close(row['polar_dimension_m']/row['equatorial_dimension_m'],row['aspect_ratio'],'Source trace aspect ratio differs')
        _close(row['diffusivity_effective_m2_s']/row['diffusivity_raw_m2_s'],row['ventilation_vapor'],'Vapor ventilation diagnostic differs')
        _close(row['conductivity_effective_W_m_K']/row['conductivity_raw_W_m_K'],row['ventilation_heat'],'Heat ventilation diagnostic differs')
        reynolds=air['density_kg_m3']*row['terminal_speed_m_s']*max(row['equatorial_dimension_m'],row['polar_dimension_m'])/air['dynamic_viscosity_Pa_s']
        _close(reynolds,row['reynolds_number'],'Reynolds diagnostic differs')
        if row['branch_criterion_met'] not in (0.,1.):raise ValueError('Invalid actual-source branch flag')
        if config['directional_mode']==1:_close(row['inherent_growth_factor'],row['effective_growth_factor'],'Unit directional intervention differs')
        if not 0<row['deposition_density_kg_m3']<=density*(1+1e-12):raise ValueError('Unphysical deposition density')
        trace_values+=len(row)
    expected_end=native.lab.FITS[record['temperature_C']][0]
    complete=record['returncode']==0 and not record['output_validation_errors'] and bool(record['snapshots']) and record['snapshots'][-1]['age_min']==expected_end
    if record['complete']!=complete:raise ValueError('Incorrect endpoint classification')
    return dict(raw_snapshot_values_checked=raw_values,actual_source_trace_values_checked=trace_values,
        maximum_relative_ice_deposition_residual=record['maximum_relative_ice_deposition_residual'])


def summarize(record,path):
    rows=record['snapshots'];config=record['configuration']
    return dict(run=str(Path(path).relative_to(ROOT)),run_sha256=native.sha(Path(path)/'run.json'),
        mean_mode=config['mean_mode'],directional_mode=config['directional_mode'],growth_mode=config['growth_mode'],
        temperature_C=record['temperature_C'],initialization=record['initialization'],dt_s=config['dt_s'],
        complete=record['complete'],returncode=record['returncode'],output_validation_errors=record['output_validation_errors'],
        requested_end_age_min=native.lab.FITS[record['temperature_C']][0],
        endpoint=rows[-1] if record['complete'] else None,last_saved=rows[-1] if rows else None,
        age_three=next((row for row in rows if row['age_min']==3.),None),
        age_ten=next((row for row in rows if row['age_min']==10.),None),
        last_source_trace=record['growth_trace'][-1] if record['growth_trace'] else None,
        branch_diagnostic=record['branch_diagnostic'],step_diagnostic=record['step_diagnostic'],
        maximum_relative_ice_deposition_residual=record['maximum_relative_ice_deposition_residual'])


def history_comparison(fitted,spherical):
    """Report signed final ratios separately from full empirical axis/mass errors."""
    shared=reference.relative_difference(fitted,spherical)
    endpoints=None
    if shared['endpoint'] is not None:
        a=fitted['snapshots'][-1];b=spherical['snapshots'][-1]
        endpoints={q:a[q]/b[q] for q in QUANTITIES}
    return dict(common_age_min=shared['common_age_min'],both_complete=shared['both_complete'],
        maximum_relative_difference=shared['trajectory_maximum'],fitted_over_spherical_endpoint=endpoints)


def assemble(cache=CACHE,baseline=BASELINE):
    cache=Path(cache);baseline=Path(baseline);records={};cases=[];checks=[];parity=[]
    for args in matrix():
        mean,directional,growth,temperature,initialization,dt=args
        path=cache/'runs'/case_name(*args)
        checks.append(verify_run(path,cache/'build'))
        record=json.loads((path/'run.json').read_text());records[args]=record;cases.append(summarize(record,path))
        if (mean,directional)==(0,0):
            oldpath=baseline/'runs'/reference.case_name('2026',growth,temperature,initialization,dt)
            item=dict(case=case_name(*args),baseline_available=(oldpath/'run.json').exists())
            if item['baseline_available']:
                reference.verify_run(oldpath,baseline/'build_2026')
                comparison=controls.parity(oldpath,path)
                if not comparison['all_equal']:raise ValueError('Unchanged-mode parity failed: '+case_name(*args))
                item.update(baseline_run=str(oldpath.relative_to(ROOT)),comparison=comparison)
            parity.append(item)
    refinements=[];histories=[]
    for mean,directional in COMBINATIONS:
        for growth in (1,2):
            for temperature in TEMPERATURES:
                for initialization,steps in HISTORIES:
                    a,b=(records[mean,directional,growth,temperature,initialization,dt] for dt in steps)
                    refinements.append(dict(mean_mode=mean,directional_mode=directional,growth_mode=growth,
                        temperature_C=temperature,initialization=initialization,timesteps_s=steps,
                        **reference.relative_difference(a,b)))
                fitted=records[mean,directional,growth,temperature,'fit',.2]
                sphere=records[mean,directional,growth,temperature,'sphere',.0125]
                histories.append(dict(mean_mode=mean,directional_mode=directional,growth_mode=growth,
                    temperature_C=temperature,fit_dt_s=.2,sphere_dt_s=.0125,**history_comparison(fitted,sphere)))
    complete=sum(case['complete'] for case in cases)
    return dict(schema='terluna.atmosphere.ventilation-benchmarks/1',
        evidence='Causal ventilation interventions against Earth laboratory fits; numerical controls are not empirical validation or adopted lunar laws.',
        producer={**controls.producer(),Path(__file__).name:native.sha(__file__)},
        build=json.loads((cache/'build/build.json').read_text()),
        matrix=dict(combinations=[dict(mean_mode=m,directional_mode=d) for m,d in COMBINATIONS],
            mean_modes=controls.MEAN_MODES,directional_modes=controls.DIRECTIONAL_MODES,growth_modes=[1,2],
            temperatures_C=TEMPERATURES,fit_timesteps_s=[1.,.2],sphere_timesteps_s=[.05,.0125],
            executions=len(cases),complete=complete,failed=len(cases)-complete),
        verification=dict(executions=len(checks),raw_snapshot_values_checked=sum(x['raw_snapshot_values_checked'] for x in checks),
            actual_source_trace_values_checked=sum(x['actual_source_trace_values_checked'] for x in checks),
            source_input_output_hashes_verified=True,physics_inputs_regenerated=True,
            maximum_relative_ice_deposition_residual=max(x['maximum_relative_ice_deposition_residual'] for x in checks if x['maximum_relative_ice_deposition_residual'] is not None),
            parity_cases_available=sum(x['baseline_available'] for x in parity),
            parity_files_compared=sum(x['comparison']['compared_files'] for x in parity if x['baseline_available'])),
        cases=cases,refinements=refinements,initial_history_sensitivity=histories,unchanged_mode_parity=parity,
        limits=['Fixed prescribed Earth air; no finite-reservoir feedback, particle interactions, fracture or lunar adaptation.',
            'Mean and directional interventions are independent causal controls, not fitted replacements.',
            'The same laboratory evidence helped inform native growth packages; agreement is not independent validation.',
            'Maximum dimension, mass and both axes must be assessed jointly and over the full laboratory interval.',
            'Failed runs retain their return codes and partial trajectories but have no accepted endpoint.',
            'Timestep comparisons measure observed sensitivity; they do not establish asymptotic convergence order.',
            'The source trace describes pre-step evaluated quantities, while snapshots describe end-of-step states.',
            'Mode2 removes the prolate extra term too; its oblate trajectories should equal mean mode1 under identical histories.'])


def run_matrix(cache=CACHE):
    cache=Path(cache)
    if not (cache/'build/build.json').exists():raise FileNotFoundError('Completed diagnostic build required')
    controls.verify_build(cache/'build')
    for args in matrix():
        mean,directional,growth,temperature,initialization,dt=args
        path=cache/'runs'/case_name(*args)
        if path.exists():raise FileExistsError(path)
        record=controls.run(cache/'build',path,temperature,dt,growth,initialization,mean,directional)
        endpoint=record['snapshots'][-1] if record['snapshots'] else None
        print(case_name(*args),'complete' if record['complete'] else 'FAILED',
              'last_age',endpoint['age_min'] if endpoint else None,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true')
    p.add_argument('--cache',type=Path,default=CACHE);p.add_argument('--baseline',type=Path,default=BASELINE)
    p.add_argument('--output',type=Path,default=RESULT);a=p.parse_args()
    if a.run:run_matrix(a.cache)
    product=assemble(a.cache,a.baseline);a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'matrix':product['matrix'],'verification':product['verification']},indent=2))
if __name__=='__main__':main()
