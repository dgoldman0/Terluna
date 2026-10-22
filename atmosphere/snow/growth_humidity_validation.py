"""Independent 2018 humidity evidence for diagnostic McSnow growth interventions.

The published sizes were not used to choose or tune these controls. Initial
spheres bracket reported seed span only; unknown seed shape and nominal group
temperatures prevent treating this as an exact experiment reconstruction.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from atmosphere.snow import mcsnow_ventilation_controls as controls
from atmosphere.snow.flakes import ice_saturation_pressure, R_VAPOUR

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DATA=HERE/'growth_humidity_observations.json'
SOURCE=ROOT/'research/runs/snow_scale/growth_environment_sources/takahashi2018_ams.pdf'
SOURCE_SHA='06a87318db2cc1831add739d305e2923fd4f361896ba1b5d8ade49a943d14f49'
CACHE=ROOT/'research/runs/snow_scale/growth_humidity'
BUILD=ROOT/'research/runs/snow_scale/mcsnow_ventilation/build'
OUT=HERE/'results/growth_humidity_validation.json'
MODES={'native':(0,0),'no_extra_mean':(1,0),'unit_direction':(0,1),'both':(1,1)}
SEED_SPANS_UM=(5.,10.,20.)
BASE_DT=.05
FINE_DT=.0125


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def observations():
    data=json.loads(DATA.read_text())
    if data['source']['sha256']!=SOURCE_SHA or sha(SOURCE)!=SOURCE_SHA:
        raise ValueError('2018 empirical source identity differs')
    if len(data['cases'])!=13 or [r['id'] for r in data['cases']]!=list('ABCDEFGHIJKLM'):
        raise ValueError('Humidity observation identities differ')
    for row in data['cases']:
        if any(not math.isfinite(row[k]) for k in ('nominal_temperature_C','vapor_density_excess_g_m3','observed_diameter_mm')):
            raise ValueError('Nonfinite observed input')
        if min(row['vapor_density_excess_g_m3'],row['observed_diameter_mm'])<=0:raise ValueError('Invalid observed input')
    return data


def producer():
    return {'module_sha256':sha(__file__),'observation_sha256':sha(DATA),'source_pdf_sha256':sha(SOURCE),
            'thermodynamic_module_sha256':sha(HERE/'flakes.py'),'controls':controls.producer()}


def configuration(row,mode,seed_span_um,dt,ice_density,temperature_offset_C=0.):
    if mode not in MODES:raise ValueError('Unknown intervention')
    if seed_span_um not in SEED_SPANS_UM:raise ValueError('Unknown seed-span sensitivity')
    if not math.isfinite(dt) or dt<=0 or not math.isclose(60/dt,round(60/dt),rel_tol=0,abs_tol=1e-8):
        raise ValueError('Timestep must positively divide one minute')
    if not math.isfinite(ice_density) or ice_density<=0:raise ValueError('Invalid native bulk ice density')
    if not math.isfinite(temperature_offset_C):raise ValueError('Invalid temperature offset')
    t=row['nominal_temperature_C']+temperature_offset_C+273.15
    drho=row['vapor_density_excess_g_m3']*1e-3
    e=ice_saturation_pressure(t)+drho*R_VAPOUR*t
    diameter=seed_span_um*1e-6
    mean,directional=MODES[mode]
    return dict(temperature_K=t,pressure_Pa=101000.,vapor_pressure_Pa=e,dt_s=dt,duration_s=600.,output_interval_s=60.,
                initial_mass_kg=math.pi/6*diameter**3*ice_density,initial_equatorial_m=diameter,initial_polar_m=diameter,
                growth_mode=2,prolate_interpolation=True,mean_mode=mean,directional_mode=directional,trace_every_step=False)


def run_name(row,mode,seed_span_um,dt,temperature_offset_C=0.):
    return f"{row['id']}_{mode}_seed{seed_span_um:g}_dt{dt:g}_offset{temperature_offset_C:+g}".replace('.','p')


def verified_record(path,config,metadata):
    """Reread immutable run bytes and compare regenerated forcing and observations."""
    result=json.loads((path/'run.json').read_text())
    if result['configuration']!=config or result['metadata']!=metadata:raise ValueError(f'Run forcing differs: {path}')
    if result['producer']!=controls.producer():raise ValueError('Control producer changed')
    for name,digest in {**result['inputs'],**result['outputs']}.items():
        if sha(path/name)!=digest:raise ValueError(f'Raw run identity differs: {name}')
    raw,errors=controls.native.read_csv(path/'trajectory.csv',('elapsed_s',*controls.native.QUANTITIES))
    trace,trace_errors=controls.native.read_csv(path/'growth_trace.csv',('evaluation_elapsed_s',*controls.TRACE_FIELDS))
    if raw!=result['raw_snapshots'] or trace!=result['growth_trace']:raise ValueError('Saved raw analysis differs')
    if errors or trace_errors:raise ValueError('Raw humidity outputs malformed')
    return result


def compact(path,result,row,mode,seed,dt,offset):
    snapshots=result['raw_snapshots'];last=snapshots[-1]
    trace=result['growth_trace'];last_trace=trace[-1]
    return dict(id=row['id'],mode=mode,seed_span_um=seed,dt_s=dt,temperature_offset_C=offset,
        complete=result['complete'],returncode=result['returncode'],validation_errors=result['output_validation_errors'],
        raw_directory=str(path.relative_to(ROOT)),run_record_sha256=sha(path/'run.json'),
        nominal_temperature_C=row['nominal_temperature_C'],actual_temperature_C=result['configuration']['temperature_K']-273.15,
        vapor_density_excess_g_m3=row['vapor_density_excess_g_m3'],observed_diameter_mm=row['observed_diameter_mm'],observed_habit=row['habit'],
        endpoint_age_s=last['elapsed_s'],diameter_mm=last['diameter_m']*1e3,
        observed_diameter_ratio=last['diameter_m']*1e3/row['observed_diameter_mm'] if result['complete'] else None,
        mass_kg=last['mass_kg'],equatorial_dimension_mm=last['equatorial_dimension_m']*1e3,
        polar_dimension_mm=last['polar_dimension_m']*1e3,aspect_ratio=last['aspect_ratio'],
        predicted_primary_habit='oblate' if last['aspect_ratio']<1 else 'prolate' if last['aspect_ratio']>1 else 'sphere',
        habit_limit='Spheroid aspect does not identify fern, dendrite, stellar or sector branches.',
        terminal_speed_m_s=last['terminal_speed_m_s'],air=result['air'],
        endpoint_source_terms=last_trace,
        maximum_saved_ventilation_vapor=max(t['ventilation_vapor'] for t in trace),
        maximum_saved_mass_balance_relative=max(abs(s['ice_mass_balance_residual_kg'])/s['mass_kg'] for s in snapshots),
        step_diagnostic=result['step_diagnostic'],
        trajectory=[{k:s[k] for k in ('elapsed_s','mass_kg','diameter_m','equatorial_dimension_m','polar_dimension_m','aspect_ratio')} for s in snapshots])


def execute_one(build_dir,cache,row,mode,seed,dt,density,offset,do_run):
    config=configuration(row,mode,seed,dt,density,offset)
    metadata={'experiment':'Takahashi 2018 Figure 1','datum_id':row['id'],'source_pdf_sha256':SOURCE_SHA,
              'observations_sha256':sha(DATA),'nominal_temperature_C':row['nominal_temperature_C'],
              'temperature_offset_C':offset,'seed_span_um':seed,'seed_shape':'dense sphere sensitivity; actual seed shape unreported',
              'control_label':mode,'humidity_input_kind':'vapor density excess above ice saturation',
              'vapor_density_excess_g_m3':row['vapor_density_excess_g_m3'],'saturation_convention':'Murphy-Koop over ice',
              'conversion_Rv_J_kg_K':R_VAPOUR,'bulk_ice_density_kg_m3':density}
    path=cache/run_name(row,mode,seed,dt,offset)
    if not path.exists():
        if not do_run:raise FileNotFoundError(path)
        controls.run_configuration(build_dir,path,config,metadata=metadata)
    result=verified_record(path,config,metadata)
    if result['build_record_sha256']!=sha(build_dir/'build.json'):raise ValueError('Humidity run build differs')
    return compact(path,result,row,mode,seed,dt,offset)


def compare_refinement(coarse,fine):
    if not coarse['complete'] or not fine['complete']:return {'complete_pair':False}
    out={'complete_pair':True}
    for field in ('diameter_mm','mass_kg','equatorial_dimension_mm','polar_dimension_mm'):
        out[field+'_endpoint_relative_change']=fine[field]/coarse[field]-1
    for field in ('diameter_m','mass_kg','equatorial_dimension_m','polar_dimension_m'):
        out[field+'_max_trajectory_relative_change']=max(abs(a[field]/b[field]-1) for a,b in zip(coarse['trajectory'],fine['trajectory']))
    return out


def assemble(*,do_run=False,build_dir=BUILD,cache=CACHE,output=OUT):
    build_dir=Path(build_dir);cache=Path(cache);output=Path(output)
    data=observations();own=producer();build=controls.verify_build(build_dir)
    if do_run:cache.mkdir(parents=True,exist_ok=True)
    import re
    m=re.search(r'\brhoi\s*=\s*([0-9.]+)_wp',(build_dir/'source/src/mo_atmo_types.f90').read_text(),re.I)
    if not m:raise ValueError('Native ice density unavailable')
    density=float(m[1]);runs=[];index={}
    for row in data['cases']:
        for mode in MODES:
            for seed in SEED_SPANS_UM:
                r=execute_one(build_dir,cache,row,mode,seed,BASE_DT,density,0.,do_run);runs.append(r);index[(row['id'],mode,seed,BASE_DT,0.)]=r
        if do_run:print(f"completed nominal {row['id']}",flush=True)
    # All 10-micrometre cases are refined; endpoints with largest absolute log error
    # additionally refine the 5 and 20-micrometre starts. Selection does not tune physics.
    fine_keys={(row['id'],mode,10.) for row in data['cases'] for mode in MODES}
    for mode in MODES:
        admissible=[r for r in runs if r['mode']==mode and r['complete']]
        if admissible:
            worst=max(admissible,key=lambda r:abs(math.log(r['observed_diameter_ratio'])))['id']
            fine_keys.update((worst,mode,seed) for seed in SEED_SPANS_UM)
    byid={r['id']:r for r in data['cases']};refinements=[]
    for id_,mode,seed in sorted(fine_keys):
        r=execute_one(build_dir,cache,byid[id_],mode,seed,FINE_DT,density,0.,do_run);runs.append(r);index[(id_,mode,seed,FINE_DT,0.)]=r
        refinements.append(dict(id=id_,mode=mode,seed_span_um=seed,**compare_refinement(index[(id_,mode,seed,BASE_DT,0.)],r)))
    # Temperature uncertainty probes at three nominal groups, selected before results.
    # These offsets are sensitivity inputs, not claimed recovered actual temperatures.
    for id_ in ('B','G','K'):
        for mode in MODES:
            for offset in (-.3,.3):
                runs.append(execute_one(build_dir,cache,byid[id_],mode,10.,FINE_DT,density,offset,do_run))
    summaries=[]
    for row in data['cases']:
        for mode in MODES:
            members=[index[(row['id'],mode,seed,FINE_DT,0.)] if (row['id'],mode,seed,FINE_DT,0.) in index else index[(row['id'],mode,seed,BASE_DT,0.)] for seed in SEED_SPANS_UM]
            good=[r for r in members if r['complete']]
            sensitivity=[r for r in runs if r['id']==row['id'] and r['mode']==mode and r['seed_span_um']==10 and r['dt_s']==FINE_DT and r['temperature_offset_C']!=0]
            summaries.append(dict(id=row['id'],mode=mode,complete_seed_cases=len(good),
                predicted_seed_span_diameter_mm_range=[min(r['diameter_mm'] for r in good),max(r['diameter_mm'] for r in good)] if good else None,
                diameter_ratio_range=[min(r['observed_diameter_ratio'] for r in good),max(r['observed_diameter_ratio'] for r in good)] if good else None,
                temperature_probe_diameter_mm_range=[min(r['diameter_mm'] for r in sensitivity),max(r['diameter_mm'] for r in sensitivity)] if sensitivity else None))
    if own!=producer() or build!=controls.verify_build(build_dir):raise RuntimeError('Humidity study identity changed')
    result=dict(schema='terluna.atmosphere.growth-humidity-validation/1',producer=own,
        source_observations=data,build_record_sha256=sha(build_dir/'build.json'),native_source_identity=build['native_source_identity'],
        evidence='Independent 2018 experimental evidence relative to local 1991 comparisons, same apparatus/author; diagnostic controls evaluated without fitting; not blind validation or exact reproduction.',
        assumptions=['Reported seed span 5–20 micrometres is represented by dense spheres, not reported seed shapes.',
            'Figure 1 nominal group temperatures and site mean 1010hPa pressure are prescribed; individual runs differ.',
            'Reported density excess is converted with Murphy–Koop saturation and project water-vapor gas constant.',
            '2018 sizes are reserved checks; no growth parameter is fitted to them.',
            'Ranges across seed span and selected temperature probes are sensitivities, not statistical confidence intervals.',
            'Spheroid aspect does not resolve secondary habit; mass, thickness and speed lack matched observations.',
            'Water vapor and temperature are held fixed; no aggregation, riming, sedimentation or finite-water coupling.'],
        counts={'runs':len(runs),'completed':sum(r['complete'] for r in runs),'refinement_pairs':len(refinements),'temperature_probe_runs':24},
        runs=runs,refinements=refinements,summaries=summaries)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--build',type=Path,default=BUILD)
    p.add_argument('--cache',type=Path,default=CACHE);p.add_argument('--output',type=Path,default=OUT);a=p.parse_args()
    result=assemble(do_run=a.run,build_dir=a.build,cache=a.cache,output=a.output)
    print(json.dumps(result['counts'],indent=2))
if __name__=='__main__':main()
