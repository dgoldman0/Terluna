"""Follow-up: separate sticking aggregation from both deposition and chipping.

This is an additional eight-case study, motivated by the completed growth-off
comparison's residual upper-tail sensitivity. It reuses the unchanged instrumented
solver and core case runner; no baseline or primary matrix is replaced.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time
from atmosphere.snow import column_diagnosis_runs as runs

ROOT=runs.ROOT
OUT=ROOT/'atmosphere/snow/results/column_diagnosis_aggregation_only.json'


def cases():
    return [dict(id=f'aggregation_only_v{scale}_r{repeat}',group=f'aggregation_only_v{scale}',family='aggregation_only',
        parameters=dict(layers=30,max_step_s=60.,area_m2=.0001*scale,breakup=False,growth=False,
                        seed=20261021+repeat,settling_scheme='upwind'))
        for scale in (1,4) for repeat in range(4)]


def identity():
    result=runs.identity(runs.model.sf.moon_column()[1])
    source=Path(__file__).resolve();result['files'][str(source.relative_to(ROOT))]=hashlib.sha256(source.read_bytes()).hexdigest()
    return result


def assemble():
    rows=[];current=runs.identity(runs.model.sf.moon_column()[1])
    for spec in cases():
        path=runs.CACHE/(spec['id']+'.json')
        if not path.exists():rows.append(dict(id=spec['id'],group=spec['group'],status='not_performed'));continue
        row=json.loads(path.read_text())
        if row['producer']!=current or row['specification']!=spec:raise RuntimeError('Stale aggregation-only case')
        if row['status']=='complete':runs.validate_trace(row['result'])
        rows.append(row)
    groups={name:runs.summarize([r for r in rows if r['group']==name]) for name in ('aggregation_only_v1','aggregation_only_v4')}
    comparison=runs.compare(groups,'aggregation_only_v1','aggregation_only_v4','volume')
    complete=all(r['status']=='complete' for r in rows)
    result=dict(schema='terluna.atmosphere.snow-column-diagnosis-controls/1',study='aggregation_only_followup',producer=identity(),
        baseline=dict(path=str(runs.BASELINE.relative_to(ROOT)),sha256=runs.BASELINE_SHA),
        primary_diagnosis=dict(path=str(runs.OUT.relative_to(ROOT)),sha256=hashlib.sha256(runs.OUT.read_bytes()).hexdigest()),
        expected_cases=8,completed_cases=sum(r['status']=='complete' for r in rows),cases=rows,groups=groups,comparisons=[comparison],
        diagnostic_columns=runs.FIELDS,checks=dict(complete=complete,
            ledgers_passed=complete and all(g['max_absolute_relative_water_residual']<1e-10 and g['max_absolute_relative_number_residual']<1e-10 for g in groups.values())),
        criterion='Same volume screen as the primary matrix: sampled mean mass-CDF change <=0.05; mean flux/d50/d90/d99 changes <=5%.',
        intervention='Disable both vapour growth and collision chipping; retain sticking aggregation, source and settling. Same full duration, geometry and source as the primary matrix.',
        evidence='Follow-up numerical intervention. A complete stochastic column with source and transport remains; this is not an isolated stationary box or proof of an infinite-volume limit.',
        reading_rule='Four independent realizations at each volume. Mean quantiles, ensemble-mean CDF and their sampling limitations retain the primary study conventions.')
    runs.old.write(OUT,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wait-for-primary',action='store_true')
    p.add_argument('--workers',type=int,choices=(1,2),default=2);args=p.parse_args()
    if args.wait_for_primary:
        print('Waiting for the primary pool to finish and write its complete product.',flush=True)
        deadline=time.monotonic()+10800
        while not runs.OUT.exists():
            if time.monotonic()>deadline:raise RuntimeError('Primary matrix was not produced within the execution allowance')
            time.sleep(5)
    primary=json.loads(runs.OUT.read_text())
    if not primary['checks']['complete']:raise RuntimeError('Complete primary study required before this follow-up')
    before=identity();current=runs.identity(runs.model.sf.moon_column()[1]);pending=[]
    for spec in cases():
        path=runs.CACHE/(spec['id']+'.json')
        if path.exists():
            row=json.loads(path.read_text())
            if row['status']=='complete' and row['producer']==current and row['specification']==spec:continue
        pending.append(spec)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(runs.run_case,spec) for spec in pending]):future.result()
    if identity()!=before:raise RuntimeError('Follow-up producer changed during execution')
    print(json.dumps(assemble()['comparisons']),flush=True)

if __name__=='__main__':main()
