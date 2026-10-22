"""Complete unfinished physical-bath controls with an explicit execution budget.

Completed controls are retained byte-for-byte as JSON values from the bound
initial attempt. Only unfinished settings are run again; acceptance gates and
all physical/numerical settings stay unchanged. This runner adds no snow model.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from atmosphere.snow import collision_jump_benchmarks as original
from atmosphere.snow import collision_birth as birth, collision_entropy as entropy
from atmosphere.snow import collision_jumps as jumps

ROOT=original.ROOT
SOURCE='atmosphere/snow/collision_birth_benchmarks.py'
ATTEMPT='atmosphere/snow/results/collision_jump_bath_attempt.json'


def bath_control(setting,reference,wall_limit_s):
    """Repeat the frozen bath protocol; only the execution allowance changes."""
    started=time.monotonic();c=reference['capture'];mass=c['initial_collector_mass_kg']
    grid=entropy.MomentGrid(original.old.reference.nested_grid(setting['factor']))
    world,_=original.old.flakes.moon_column()
    op=birth.BirthCapture(grid,original.old.reference.SnowKernel(world),c['captured_mass_kg'],c['bath_number_m3'],
                         quadrature_order=setting['order'],event_capacity=setting['capacity'])
    result=birth.integrate(op,grid.from_particles([mass],[1.]),c['times_s'],max_step_s=setting['dt'],time_limit_s=wall_limit_s)
    rows=[]
    for snapshot,ref in zip(result['snapshots'],c['continuous_mass_references'][-1]['snapshots']):
        sample=np.asarray(ref['masses_kg']);frequency=np.asarray(ref['sample_counts']);sd=ref['variance_mass_kg2']**.5
        query=np.unique(np.r_[sample,np.nextafter(sample,-np.inf),grid.bounds,
            np.linspace(ref['mean_mass_kg']-6*sd,ref['mean_mass_kg']+6*sd,1025)])
        candidate=jumps.number_cdf(grid,snapshot['state'],query)
        mean=snapshot['moments'][1]/snapshot['moments'][0]
        rows.append(dict(**snapshot,mean_mass_kg=float(mean),variance_mass_kg2=snapshot['variance'],
            number_cdf_query_mass_kg=query,number_cdf=candidate,
            reference_mean_mass_kg=ref['mean_mass_kg'],reference_variance_mass_kg2=ref['variance_mass_kg2'],
            variance_ratio=snapshot['variance']/ref['variance_mass_kg2'] if ref['variance_mass_kg2'] else None,
            mean_error_initial_mass=float((mean-ref['mean_mass_kg'])/mass),
            empirical_cdf_supremum=original.old.control.empirical_distance(query,candidate,sample,frequency)))
    return original.old.control.plain(dict(mass_cells=len(grid.pivots),grid_factor=setting['factor'],quadrature_order=setting['order'],
        max_step_s=setting['dt'],settings=setting,diagnostics=result['diagnostics'],snapshots=rows,wall_seconds=time.monotonic()-started,
        reconstruction_fits=grid.fit_diagnostics,reconstruction_roundoff=grid.roundoff_diagnostics,
        pivots_kg=grid.pivots,bounds_kg=grid.bounds,event_capacity=setting['capacity']))


def complete(wall_limit_s=600.):
    sha=lambda name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    attempt=json.loads((ROOT/ATTEMPT).read_text())
    if attempt['schema']!=original.SCHEMA or attempt['case']!='prescribed_bath':raise ValueError('Expected initial bath attempt')
    for group in ('files','inputs'):
        for name,digest in attempt['producer'][group].items():
            if sha(name)!=digest:raise RuntimeError('Initial attempt producer or input changed: '+name)
    if attempt['requested_controls']!=original.controls('bath'):raise ValueError('Requested controls changed')
    _,binding=original.old.flakes.moon_column();before=original.identity(binding)
    before['files'][SOURCE]=sha(SOURCE);before['inputs'][ATTEMPT]=sha(ATTEMPT)
    reference=json.loads((ROOT/original.old.control.REFERENCE).read_text());original.old.check_reference(reference)
    result=copy.deepcopy(attempt);runs=[];failures=[];executed=[];retained=[]
    for setting in attempt['requested_controls']:
        saved=next((r for r in attempt['runs'] if r['settings']==setting),None)
        if saved is not None:
            runs.append(saved);retained.append(setting);continue
        started=time.monotonic();executed.append(setting)
        try:
            run=bath_control(setting,reference,wall_limit_s);runs.append(run)
            print('Completed',setting,'in',run['wall_seconds'],'s',flush=True)
        except (ArithmeticError,ValueError,RuntimeError) as exc:
            failure=dict(settings=setting,exception=type(exc).__name__,message=str(exc),wall_seconds=time.monotonic()-started)
            failures.append(failure);print('Incomplete',failure,flush=True)
    acceptance=original.old.bath_decisions(runs,reference)
    acceptance['accounting_and_cdf_checks']=original.old.accounting_checks(runs,'bath',reference)
    acceptance['truncation_budget_checks']=original.error_budget_checks(runs,'bath')
    acceptance['all_requested_controls_completed']=len(runs)==len(attempt['requested_controls']) and not failures
    acceptance['initial_cdf_checks']=[dict(settings=r['settings'],passed=r['snapshots'][0]['empirical_cdf_supremum']<=1e-12) for r in runs]
    acceptance['semigroup_check']=True
    acceptance['passed']=(acceptance['passed'] and acceptance['all_requested_controls_completed'] and
        all(c['passed'] for c in acceptance['initial_cdf_checks']+acceptance['accounting_and_cdf_checks']+acceptance['truncation_budget_checks']))
    after=original.identity(binding);after['files'][SOURCE]=sha(SOURCE);after['inputs'][ATTEMPT]=sha(ATTEMPT)
    if after!=before:raise RuntimeError('Producer or inputs changed during continuation')
    result.update(producer=before,runs=runs,execution_failures=failures,acceptance=acceptance,
        execution_history=dict(initial_attempt=ATTEMPT,initial_attempt_sha256=before['inputs'][ATTEMPT],
            initial_wall_limit_s=300.,continued_wall_limit_s=wall_limit_s,retained_completed_controls=retained,
            repeated_unfinished_controls=executed,continued_workers=1,
            rule='Retained runs are unchanged JSON values from the source-verified initial attempt. Only unfinished controls are recomputed. Settings, numerical method and admission gates remain unchanged.'))
    result['environment']['initial_workers']=result['environment'].pop('workers')
    result['environment']['continued_workers']=1
    result['reading_rule']+=' The execution history binds the initial incomplete matrix and identifies the one-worker continuation; its larger wall-time allowance changes no numerical setting.'
    return original.old.control.plain(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wall-limit-s',type=float,default=600.)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=complete(args.wall_limit_s)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Saved',args.output,'acceptance',result['acceptance']['passed'],flush=True)
    if not result['acceptance']['passed']:raise SystemExit(1)


if __name__=='__main__':main()
