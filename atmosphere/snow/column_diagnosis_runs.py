"""Controlled process removals and observed vapour competition in the full column.

Core matrix: four realizations per group, two volumes for each process removal,
three for the baseline, and an independent split-step control. At most two workers.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import copy
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from atmosphere.snow import column_diagnosis as model, event_column_benchmarks as old

ROOT=model.ROOT
CACHE=ROOT/'research/runs/snow_column_diagnosis'
OUT=ROOT/'atmosphere/snow/results/column_diagnosis.json'
BASELINE=old.OUT
BASELINE_SHA='54d6d2d4b3ffb449a77934f822e9addc682269bed10c603a73fee08da9a4fef7'
FIELDS=model.DIAGNOSTIC_COLUMNS


def identity(binding):
    result=model.identity(binding)
    for path in (Path(__file__).resolve(),Path(old.__file__).resolve()):
        result['files'][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def cases():
    groups=[]
    for family,changes in [('full',{}),('no_deposition',{'growth':False}),
                           ('no_aggregation',{'collisions':False}),('no_chipping',{'breakup':False})]:
        for scale in ([1,4,16] if family=='full' else [1,4]):
            groups.append((f'{family}_v{scale}',family,dict(area_m2=.0001*scale,**changes)))
    groups.append(('full_dt30_v4','full_dt30',dict(area_m2=.0004,max_step_s=30.)))
    return [dict(id=f'{name}_r{repeat}',group=name,family=family,
        parameters=dict(layers=30,max_step_s=60.,area_m2=.0001,breakup=True,seed=20261021+repeat,
                        settling_scheme='upwind')|changes)
        for name,family,changes in groups for repeat in range(4)]


def scientific(row):
    result=copy.deepcopy(row)
    result.pop('schema');result.pop('process_diagnostics',None)
    result['diagnostics'].pop('build');result['diagnostics'].pop('wall_seconds')
    return result


def original_case(spec):
    if hashlib.sha256(BASELINE.read_bytes()).hexdigest()!=BASELINE_SHA:raise RuntimeError('Baseline product changed')
    if spec['family'] not in ('full','full_dt30'):return None
    name=spec['id'].replace('full_dt30_v4','moon_dt30').replace('full_v','moon_v')
    return next(r['result'] for r in json.loads(BASELINE.read_text())['cases'] if r['id']==name)


def validate_trace(raw):
    a=np.array(raw['process_diagnostics']['blocks']);s=raw['settings'];area=s['area_m2']
    if raw['process_diagnostics']['columns']!=FIELDS or a.shape!=(5,s['layers'],len(FIELDS)):
        raise ValueError('Unexpected diagnostic shape or fields')
    if not np.isfinite(a).all() or (a<0).any():raise ValueError('Invalid diagnostic integrals')
    duration=s['total_time_s'];gather=s['gathering_time_s'];expected=np.r_[duration-gather,np.full(4,gather/4)]
    if not np.allclose(a[:,:,0],expected[:,None],rtol=2e-12,atol=1e-10):raise ValueError('Wrong sample time')
    cap=s['snowfall_mm_h']/3600/s['depth_m']*s['layer_volume_m3'] if s['growth'] else 0.
    if not np.allclose(a[:,:,5]+a[:,:,6]+a[:,:,7],cap*a[:,:,0],rtol=2e-12,atol=1e-20):raise ValueError('Vapour classification does not close')
    if not np.allclose(a[:,:,6],cap*a[:,:,1],rtol=2e-12,atol=1e-20):raise ValueError('Empty-layer loss does not match occupancy')
    if not np.isclose(a[:,:,5].sum()/area,raw['mass_ledger_kg_m2']['deposited'],rtol=2e-11,atol=1e-11):raise ValueError('Deposition observer disagrees with water ledger')
    if a[:,:,12].sum()/area!=raw['number_ledger_m2']['coagulation_loss']:raise ValueError('Sticking count mismatch')
    if a[:,:,13].sum()!=raw['diagnostics']['chipping_events']:raise ValueError('Chipping count mismatch')
    if a[:,:,14].sum()!=raw['diagnostics']['emitted_fragments']:raise ValueError('Fragment count mismatch')
    if not s['growth'] and np.any(a[:,:,3:8]):raise ValueError('Disabled growth has water demand')
    return True


def run_case(spec):
    world,binding=model.sf.moon_column();before=identity(binding);start=time.monotonic()
    print('starting '+spec['id'],flush=True)
    try:
        raw=model.column(world,wall_limit_s=1800.,**spec['parameters'])
        validate_trace(raw)
        reference=original_case(spec)
        parity=scientific(raw)==scientific(reference) if reference else None
        if parity is False:raise RuntimeError('Instrumented baseline differs from frozen case')
        if identity(binding)!=before:raise RuntimeError('Producer changed during run')
        row=dict(id=spec['id'],group=spec['group'],family=spec['family'],status='complete',
                 specification=spec,producer=before,baseline_scientific_parity=parity,result=raw)
    except Exception as error:
        row=dict(id=spec['id'],group=spec['group'],family=spec['family'],status='failed',
                 specification=spec,producer=before,error=str(error),wall_seconds=time.monotonic()-start)
    old.write(CACHE/(spec['id']+'.json'),row)
    print(json.dumps(dict(id=spec['id'],status=row['status'],elapsed_s=time.monotonic()-start,
        summary=row.get('result',{}).get('summary'),error=row.get('error')),default=old.plain),flush=True)
    return row


def metrics(raw):
    s=raw['settings'];a=np.array(raw['process_diagnostics']['blocks']);g=a[1:].sum(axis=0)
    time=g[:,0];area=s['area_m2'];hours=s['gathering_time_s']/3600
    def ratio(num,den):return np.divide(num,den,out=np.zeros_like(num),where=den>0)
    profile=dict(empty_fraction=g[:,1]/time,single_particle_fraction=g[:,2]/time,
        vapour_limited_fraction=g[:,3]/time,mean_particles=g[:,8]/time,
        particle_number_variance=np.maximum(0,g[:,9]/time-(g[:,8]/time)**2),
        mean_effective_growth_number=g[:,10]/time,
        largest_growth_share_when_occupied=ratio(g[:,11],time-g[:,1]),
        admitted_mm_h=g[:,5]/area/hours,empty_unused_mm_h=g[:,6]/area/hours,
        occupied_unused_mm_h=g[:,7]/area/hours,
        sticking_events_m2_h=g[:,12]/area/hours,chipping_events_m2_h=g[:,13]/area/hours)
    summary=dict(gather_deposition_mm_h=float(profile['admitted_mm_h'].sum()),
        gather_empty_unused_mm_h=float(profile['empty_unused_mm_h'].sum()),
        gather_occupied_unused_mm_h=float(profile['occupied_unused_mm_h'].sum()),
        whole_deposition_kg_m2=float(a[:,:,5].sum()/area),whole_empty_unused_kg_m2=float(a[:,:,6].sum()/area),
        whole_occupied_unused_kg_m2=float(a[:,:,7].sum()/area),
        bottom_empty_fraction=float(profile['empty_fraction'][0]),
        bottom_single_particle_fraction=float(profile['single_particle_fraction'][0]),
        bottom_vapour_limited_fraction=float(profile['vapour_limited_fraction'][0]),
        bottom_mean_particles=float(profile['mean_particles'][0]),
        bottom_effective_growth_number=float(profile['mean_effective_growth_number'][0]))
    return summary,profile


def summarize(rows):
    result=old.summarize_group(rows)
    if result['status']!='complete':return result
    allmetrics=[metrics(r['result']) for r in rows]
    for kind,index in [('process_summary',0),('process_profiles',1)]:
        means={};errors={}
        for name in allmetrics[0][index]:
            values=np.array([m[index][name] for m in allmetrics])
            means[name]=values.mean(axis=0);errors[name]=values.std(axis=0,ddof=1)/np.sqrt(len(values))
        result['mean_'+kind]=means;result['standard_error_'+kind]=errors
    result['mean_process_blocks']=np.mean([r['result']['process_diagnostics']['blocks'] for r in rows],axis=0)
    return result


def compare(groups,left,right,axis):
    a,b=groups[left],groups[right]
    if a['status']!='complete' or b['status']!='complete':return dict(left=left,right=right,axis=axis,status='incomplete',passed=False)
    cdf=float(np.max(np.abs(a['mean_mass_cdf']-b['mean_mass_cdf'])))
    changes={k:float(b['mean_summary'][k]/a['mean_summary'][k]-1) for k in ('snowfall_mm_h','mass_d50_m','mass_d90_m','mass_d99_m')}
    return dict(left=left,right=right,axis=axis,status='complete',maximum_mean_mass_cdf_difference=cdf,
        relative_summary_changes=changes,passed=cdf<=.05 and all(abs(v)<=.05 for v in changes.values()))


def assemble():
    rows=[];before=identity(model.sf.moon_column()[1])
    for spec in cases():
        path=CACHE/(spec['id']+'.json')
        if not path.exists():rows.append(dict(id=spec['id'],group=spec['group'],status='not_performed'));continue
        row=json.loads(path.read_text())
        if row['producer']!=before or row['specification']!=spec:raise RuntimeError('Stale case '+spec['id'])
        if row['status']=='complete':validate_trace(row['result'])
        rows.append(row)
    groups={name:summarize([r for r in rows if r['group']==name]) for name in dict.fromkeys(r['group'] for r in rows)}
    pairs=[('full_v1','full_v4','volume'),('full_v4','full_v16','volume'),
           ('full_v4','full_dt30_v4','time')]
    pairs += [(f'{family}_v1',f'{family}_v4','volume') for family in ('no_deposition','no_aggregation','no_chipping')]
    comparisons=[compare(groups,*pair) for pair in pairs]
    complete=all(r['status']=='complete' for r in rows)
    ledgers=complete and all(g['max_absolute_relative_water_residual']<1e-10 and g['max_absolute_relative_number_residual']<1e-10 for g in groups.values())
    result=dict(schema='terluna.atmosphere.snow-column-diagnosis-controls/1',producer=before,
        baseline=dict(path=str(BASELINE.relative_to(ROOT)),sha256=BASELINE_SHA),
        expected_cases=len(cases()),completed_cases=sum(r['status']=='complete' for r in rows),cases=rows,groups=groups,
        comparisons=comparisons,diagnostic_columns=FIELDS,
        checks=dict(complete=complete,ledgers_passed=ledgers,
                    all_instrumented_baselines_identical=complete and all(r.get('baseline_scientific_parity') is not False for r in rows)),
        criterion='Refinement screen: mean mass-CDF change <=0.05 and mean flux/d50/d90/d99 changes <=5%. These are practical thresholds, not sampling confidence bounds.',
        interventions=dict(full='All inherited processes.',no_deposition='Vapour growth disabled; collisions, chipping, injection and settling retained.',
            no_aggregation='All encounters disabled, including collision-generated chipping; vapour growth retained.',
            no_chipping='Sticking retained; nonsticking encounters cause no mass loss or fragment production.',
            full_dt30='All processes, half macrostep at the middle volume.'),
        evidence='Numerical process-removal experiments and passive accounting under the same prescribed lunar forcing. '
        'Removal changes the evolving particle state and does not isolate an additive contribution. Empty and under-demand '
        'supply losses account for instantaneous unused capacity, not a counterfactual prediction after restoring that water.',
        reading_rule='Whole realizations provide standard errors. Each trace has filling and four gathering blocks, then '
        'layers ordered from freezing level upward. Demand uses both Heun-stage quadrature weights; occupancy and limitation '
        'are substep observations. No tail deletion, physical cutoff, empirically calibrated fracture or storm optical forecast.')
    old.write(OUT,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=[r['id'] for r in cases()]);p.add_argument('--pilot',action='store_true')
    p.add_argument('--run-all',action='store_true');p.add_argument('--assemble',action='store_true')
    p.add_argument('--workers',type=int,choices=[1,2],default=2);args=p.parse_args()
    selected=cases() if args.run_all else []
    if args.pilot:selected=[r for r in cases() if r['id'] in ('full_v1_r0','no_deposition_v1_r0','no_aggregation_v1_r0','no_chipping_v1_r0')]
    if args.case:selected=[r for r in cases() if r['id']==args.case]
    pending=[];current=identity(model.sf.moon_column()[1])
    for spec in selected:
        path=CACHE/(spec['id']+'.json')
        if path.exists():
            row=json.loads(path.read_text())
            if row['status']=='complete' and row['producer']==current and row['specification']==spec:continue
        pending.append(spec)
    if pending:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(run_case,r) for r in pending]
            for future in as_completed(futures):future.result()
    if args.assemble or args.run_all:print(json.dumps(assemble()['checks']),flush=True)

if __name__=='__main__':main()
