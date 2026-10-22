"""Reproducible independent closed-box and complete Earth/lunar column matrix.

At most two numerical workers. Every expected case remains in the report,
including explicit failures. Numerical acceptance never means physical validation.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from atmosphere.snow import event_column as model, flakes as sf, collision_reference as ref

ROOT=model.ROOT
CACHE=ROOT/'research/runs/snow_event_column'
OUT=ROOT/'atmosphere/snow/results/event_column.json'
BOX=ROOT/'atmosphere/snow/results/event_closed_boxes.json'
REFERENCE=ROOT/'atmosphere/snow/results/collision_reference.json'
SCHEMA='terluna.atmosphere.snow-event-column-controls/1'
REFERENCE_SHA='ccfe6c2fc8ceaceefdda4814f3c122b9dd12f7954007e843d0392507239b74c9'


def plain(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=plain,allow_nan=False)+'\n')
    temporary.replace(path)


def identity(binding):
    result=model.identity(binding)
    path=Path(__file__).resolve()
    result['files'][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def cases():
    settings=[('moon_v1',dict(area_m2=.0001)),('moon_v4',dict(area_m2=.0004)),
              ('moon_v16',dict(area_m2=.0016)),
              ('moon_dt30',dict(area_m2=.0004,max_step_s=30.)),
              ('moon_z60',dict(area_m2=.0008,layers=60)),
              ('earth_v4',dict(area_m2=.0004,snowfall_mm_h=1.)),
              ('earth_v16',dict(area_m2=.0016,snowfall_mm_h=1.))]
    rows=[]
    for group,changes in settings:
        for repeat in range(4):
            rows.append(dict(id=f'{group}_r{repeat}',group=group,world='earth' if group.startswith('earth') else 'moon',
                             parameters=dict(layers=30,max_step_s=60.,breakup=True,seed=20261021+repeat,
                                             settling_scheme='upwind',**{k:v for k,v in changes.items() if k not in ('layers','max_step_s')})
                             | {k:v for k,v in changes.items() if k in ('layers','max_step_s')}))
    return rows


def run_case(spec):
    world,binding=sf.moon_column() if spec['world']=='moon' else (sf.earth_column(),{})
    before=identity(binding);start=time.monotonic()
    print('starting '+spec['id'],flush=True)
    try:
        raw=model.column(world,wall_limit_s=1800.,**spec['parameters'])
        if identity(binding)!=before:raise RuntimeError('Producer changed during calculation')
        row=dict(id=spec['id'],group=spec['group'],world=spec['world'],status='complete',specification=spec,
                 producer=before,result=raw)
    except Exception as error:
        row=dict(id=spec['id'],group=spec['group'],world=spec['world'],status='failed',specification=spec,
                 producer=before,error=str(error),wall_seconds=time.monotonic()-start)
    write(CACHE/(spec['id']+'.json'),row)
    print(json.dumps(dict(id=spec['id'],status=row['status'],wall_seconds=time.monotonic()-start,
          summary=row.get('result',{}).get('summary'),error=row.get('error'))),flush=True)
    return row


def boxes(repetitions=128):
    if hashlib.sha256(REFERENCE.read_bytes()).hexdigest()!=REFERENCE_SHA:raise RuntimeError('Independent reference hash mismatch')
    reference=json.loads(REFERENCE.read_text());world,binding=sf.moon_column();before=identity(binding)
    for path,digest in reference['producer']['files'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:raise RuntimeError('Reference producer changed: '+path)
    target=reference['ordinary_box'];query=model.GEOMETRY.mass(np.array(target['cdf_query_diameter_m']))
    times=np.array(target['times_s']);rows=[];start=time.monotonic()
    for original in target['exact_event_ensembles']:
        n=original['initial_particles'];species=ref.ordinary_initial();initial=np.repeat(species,n//len(species));volume=n/1e5
        sample=[]
        for repeat in range(repetitions):
            sample.append(model.box(world,initial,volume,times,query,seed=620000+n*1000+repeat))
        moments=np.array([r['moments'] for r in sample]);cdfs=np.array([r['mass_cdf'] for r in sample])
        comparisons=[]
        for j,old in enumerate(original['snapshots']):
            mean=moments[:,j].mean(axis=0);se=moments[:,j].std(axis=0,ddof=1)/math.sqrt(repetitions)
            oldmean=np.array([old['mean'][k] for k in ('number_m3','mass_kg_m3','second_mass_moment_kg2_m3')])
            oldse=np.array([old['standard_error'][k] for k in ('number_m3','mass_kg_m3','second_mass_moment_kg2_m3')])
            combined=np.sqrt(se*se+oldse*oldse)
            discrepancy=np.abs(mean-oldmean)
            cmean=cdfs[:,j].mean(axis=0);cse=cdfs[:,j].std(axis=0,ddof=1)/math.sqrt(repetitions)
            ccombined=np.sqrt(cse*cse+np.array(old['mass_cdf_standard_error'])**2)
            cdifference=np.abs(cmean-np.array(old['mass_cdf_mean']))
            moment_pass=bool(np.all(discrepancy[[0,2]]<=5*combined[[0,2]]+1e-12*np.maximum(np.abs(oldmean[[0,2]]),1e-30)))
            cdf_pass=bool(np.all(cdifference<=5*ccombined+1e-12))
            comparisons.append(dict(time_s=float(times[j]),moments_mean=mean,moments_standard_error=se,
                mass_cdf_mean=cmean,mass_cdf_standard_error=cse,
                largest_standardized_moment_difference=float(np.max(np.divide(discrepancy[[0,2]],combined[[0,2]],out=np.zeros(2),where=combined[[0,2]]>0))),
                maximum_mass_cdf_difference=float(cdifference.max()),moment_sampling_check=moment_pass,cdf_sampling_check=cdf_pass))
        rows.append(dict(initial_particles=n,volume_m3=volume,repetitions=repetitions,snapshots=comparisons,
            maximum_relative_mass_residual=max(float(np.max(np.abs(r['relative_mass_residual']))) for r in sample),
            total_proposals=sum(r['proposals'] for r in sample),
            maximum_particle_mass_kg_per_realization=[r['maximum_mass'].tolist() for r in sample]))
        print('closed boxes',n,'elapsed',time.monotonic()-start,flush=True)
    if identity(binding)!=before:raise RuntimeError('Producer changed during closed-box matrix')
    accepted=all(x['moment_sampling_check'] and x['cdf_sampling_check'] for r in rows for x in r['snapshots'])
    product=dict(schema='terluna.atmosphere.snow-event-closed-boxes/1',producer=before,
        input=dict(path=str(REFERENCE.relative_to(ROOT)),sha256=REFERENCE_SHA),runs=rows,
        query_diameter_m=target['cdf_query_diameter_m'],acceptance=dict(passed=accepted,
        criterion='Each exported number/second-moment mean and mass CDF query agrees within five combined realization standard errors plus 1e-12 relative/numerical tolerance; a compatibility screen, not a joint confidence guarantee.'),
        elapsed_wall_seconds=time.monotonic()-start,
        evidence='Independent positive-tree thinning of the finite-volume generator; matching direct-event ensembles. '
        'The twelve-hour reference remains unconverged in volume. No empirical validation or infinite-volume tail inference.')
    write(BOX,product);return product


def summarize_group(rows):
    completed=[r['result'] for r in rows if r['status']=='complete']
    if len(completed)!=len(rows):return dict(status='incomplete',complete=len(completed),expected=len(rows))
    arrays=np.array([r['arriving_cumulative_moments'] for r in completed]);totals=arrays[:,-1]
    masscdf=np.divide(arrays[:,:,1],totals[:,1,None])
    summary_names=('snowfall_mm_h','mass_d50_m','mass_d90_m','mass_d99_m','maximum_diameter_anywhere_m')
    summary=np.array([[r['summary'][k] for k in summary_names] for r in completed])
    return dict(status='complete',realizations=len(rows),settings=completed[0]['settings'],
        mean_summary=dict(zip(summary_names,summary.mean(axis=0))),
        standard_error_summary=dict(zip(summary_names,summary.std(axis=0,ddof=1)/math.sqrt(len(rows)))),
        mean_mass_cdf=masscdf.mean(axis=0),standard_error_mass_cdf=masscdf.std(axis=0,ddof=1)/math.sqrt(len(rows)),
        mean_cumulative_moments=arrays.mean(axis=0),standard_error_cumulative_moments=arrays.std(axis=0,ddof=1)/math.sqrt(len(rows)),
        mean_layer_moments=np.mean([r['time_mean_layer_moments'] for r in completed],axis=0),
        mean_time_block_moments=np.mean([r['arrival_time_blocks']['moments'] for r in completed],axis=0),
        individual_mass_d50_m=summary[:,1],individual_mass_d99_m=summary[:,3],
        max_absolute_relative_water_residual=max(abs(r['mass_ledger_kg_m2']['relative_residual']) for r in completed),
        max_absolute_relative_number_residual=max(abs(r['number_ledger_m2']['relative_residual']) for r in completed))


def assemble():
    rows=[]
    for spec in cases():
        path=CACHE/(spec['id']+'.json')
        if not path.exists():rows.append(dict(id=spec['id'],group=spec['group'],status='not_performed'));continue
        row=json.loads(path.read_text())
        world,binding=sf.moon_column() if spec['world']=='moon' else (sf.earth_column(),{})
        if row['producer']!=identity(binding) or row['specification']!=spec:raise RuntimeError('Stale case: '+spec['id'])
        rows.append(row)
    groups={name:summarize_group([r for r in rows if r['group']==name]) for name in dict.fromkeys(r['group'] for r in rows)}
    comparisons=[]
    for left,right,axis in [('moon_v1','moon_v4','volume'),('moon_v4','moon_v16','volume'),
                             ('moon_v4','moon_dt30','time'),('moon_v4','moon_z60','height at fixed layer volume'),
                             ('earth_v4','earth_v16','Earth volume')]:
        a,b=groups[left],groups[right]
        if a['status']!='complete' or b['status']!='complete':
            comparisons.append(dict(left=left,right=right,axis=axis,passed=False,status='incomplete'));continue
        delta=float(np.max(np.abs(a['mean_mass_cdf']-b['mean_mass_cdf'])))
        changes={k:float(b['mean_summary'][k]/a['mean_summary'][k]-1) for k in ('snowfall_mm_h','mass_d50_m','mass_d90_m','mass_d99_m')}
        passed=delta<=.05 and all(abs(v)<=.05 for v in changes.values())
        comparisons.append(dict(left=left,right=right,axis=axis,maximum_mean_mass_cdf_difference=delta,
                                relative_summary_changes=changes,passed=passed,status='complete'))
    complete=all(r['status']=='complete' for r in rows)
    ledgers=complete and all(g['max_absolute_relative_water_residual']<1e-10 and g['max_absolute_relative_number_residual']<1e-10 for g in groups.values())
    passed=complete and ledgers and all(r['passed'] for r in comparisons)
    product=dict(schema=SCHEMA,expected_cases=len(cases()),completed_cases=sum(r['status']=='complete' for r in rows),
        cases=rows,groups=groups,comparisons=comparisons,
        acceptance=dict(passed=passed,complete=complete,ledgers_passed=ledgers,
            criterion='All cases complete; relative water/number residuals below 1e-10; every listed refinement changes '
                      'mean mass CDF by at most 0.05 and mean flux/d50/d90/d99 by at most 5%. These practical numerical '
                      'screens do not bound sampling uncertainty or establish physical validity.'),
        evidence='Complete finite-volume prescribed-cloud calculations with explicit refinement and realization uncertainty. '
                 'No physical size cap, giant-aggregate fracture, evolving geometry, or climate-resolved moist trajectories.',
        reading_rule='Whole realizations supply standard errors. Correlated particles within a column are not independent '
                     'uncertainty samples. Mean quantiles are means of realization quantiles, not pooled-population quantiles. '
                     'Source flux and duration agree in expectation with the frozen column. Upwind settling matches its '
                     'height-space generator; continuous masses remove pivot allocation. Mass conservation does not imply '
                     'volume convergence. Failed refinements prevent use as snowfall or optical predictions.')
    write(OUT,product);return product


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--boxes',action='store_true')
    p.add_argument('--case',choices=[r['id'] for r in cases()]);p.add_argument('--run-all',action='store_true')
    p.add_argument('--assemble',action='store_true');p.add_argument('--workers',type=int,choices=(1,2),default=2)
    args=p.parse_args()
    if args.boxes: print('closed acceptance',boxes()['acceptance'],flush=True)
    if args.case:run_case(next(r for r in cases() if r['id']==args.case))
    if args.run_all:
        pending=[]
        for spec in cases():
            path=CACHE/(spec['id']+'.json')
            if path.exists():
                row=json.loads(path.read_text());binding=sf.moon_column()[1] if spec['world']=='moon' else {}
                if row['status']=='complete' and row['specification']==spec and row['producer']==identity(binding):continue
            pending.append(spec)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(run_case,spec) for spec in pending]
            for future in as_completed(futures):future.result()
        args.assemble=True
    if args.assemble:print('column acceptance',assemble()['acceptance'],flush=True)

if __name__=='__main__':main()
