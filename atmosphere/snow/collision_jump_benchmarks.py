"""Distribution and accounting controls for finite-time maintained-bath capture."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.stats import poisson
from atmosphere.snow import collision_jumps as jumps, collision_birth as birth
from atmosphere.snow import collision_entropy as entropy, collision_entropy_benchmarks as old
from shared.provenance import constants_used

ROOT=old.ROOT
SCHEMA='terluna.atmosphere.snow-finite-capture-benchmarks/1'
FILES=('atmosphere/snow/collision_jump_benchmarks.py','atmosphere/snow/collision_jumps.py',
       'atmosphere/snow/collision_birth.py',*old.FILES)


def identity(binding):
    baseline=old.reference.identity(binding)
    sha=lambda path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    return dict(files={p:sha(p) for p in FILES},inputs={old.control.REFERENCE:sha(old.control.REFERENCE),**baseline['inputs']},constants=constants_used(FILES))


def controls(case):
    # Nominal timestep, quadrature, cells, shifted faces and ladder capacity are
    # changed in separate controls. Each Poisson grid samples its own first face.
    if case=='boundary':
        return [dict(factor=f,order=q,dt=dt,phase=phase) for f,q,dt,phase in
                ((2,16,60.,1.),(2,16,30.,1.),(2,32,60.,1.),(1,16,60.,1.),(4,16,60.,1.),(2,16,60.,1.02))]
    return [dict(factor=f,order=q,dt=dt,capacity=cap) for f,q,dt,cap in
            ((2,16,.5,256),(2,16,.25,256),(2,32,.5,256),(1,16,.5,256),(4,16,.5,256),(2,16,.5,384))]


def run_control(task):
    case,setting=task;start=time.monotonic()
    try:
        p=json.loads((ROOT/old.control.REFERENCE).read_text());c=p['capture'];mass_unit=c['initial_collector_mass_kg']
        factor,order,dt=setting['factor'],setting['order'],setting['dt']
        grid=entropy.MomentGrid(old.reference.nested_grid(factor)*(setting.get('phase',1.)/mass_unit if case=='boundary' else 1.))
        if case=='boundary':
            donor=c['captured_mass_kg']/mass_unit;face=float(grid.bounds[grid.bounds>1][0]);expected=(face-1)/donor
            counts=np.r_[0.,expected+np.array([-2.,0.,2.])*np.sqrt(expected),1.5*expected];times=counts/100.
            op=jumps.ConstantCapture(grid,donor,100.,quadrature_order=order)
            result=jumps.integrate(op,grid.from_particles([1.],[1.]),times,max_step_s=dt,time_limit_s=300.)
            rows=[]
            for snapshot,count in zip(result['snapshots'],counts):
                k=np.arange(max(0,int(count-12*np.sqrt(count))),int(count+12*np.sqrt(count))+2)
                mass=1+donor*k;query=np.sort(np.r_[mass,np.nextafter(mass,-np.inf)])
                candidate=jumps.number_cdf(grid,snapshot['state'],query)
                exact=np.empty_like(query);exact[::2]=poisson.cdf(k-1,count);exact[1::2]=poisson.cdf(k,count)
                mean=snapshot['moments'][1]/snapshot['moments'][0]
                rows.append(dict(**snapshot,poisson_mean=float(count),mean_mass_ratio=float(mean),variance_mass_ratio=snapshot['variance'],
                    expected_mean_mass_ratio=float(1+donor*count),expected_variance_mass_ratio=float(count*donor**2),
                    query_mass_ratio=query,number_cdf=candidate,exact_poisson_number_cdf=exact,
                    number_cdf_supremum=float(np.max(abs(candidate-exact))),
                    unqueried_left_probability=float(poisson.cdf(k[0]-1,count)),unqueried_right_probability=float(poisson.sf(k[-1],count))))
            extra=dict(pivots_mass_ratio=grid.pivots,bounds_mass_ratio=grid.bounds,donor_mass_ratio=donor,
                       prescribed_capture_rate_s=100.,reference_cell_face_mass_ratio=face,face_sampling='own first upper face for this grid',grid_phase=setting['phase'])
        else:
            world,_=old.flakes.moon_column();kernel=old.reference.SnowKernel(world)
            op=birth.BirthCapture(grid,kernel,c['captured_mass_kg'],c['bath_number_m3'],quadrature_order=order,event_capacity=setting['capacity'])
            result=birth.integrate(op,grid.from_particles([mass_unit],[1.]),c['times_s'],max_step_s=dt,time_limit_s=300.)
            rows=[];target=c['continuous_mass_references'][-1]
            for snapshot,ref in zip(result['snapshots'],target['snapshots']):
                sample=np.asarray(ref['masses_kg']);frequency=np.asarray(ref['sample_counts'])
                sd=ref['variance_mass_kg2']**.5
                query=np.unique(np.r_[sample,np.nextafter(sample,-np.inf),grid.bounds,
                    np.linspace(ref['mean_mass_kg']-6*sd,ref['mean_mass_kg']+6*sd,1025)])
                candidate=jumps.number_cdf(grid,snapshot['state'],query)
                mean=snapshot['moments'][1]/snapshot['moments'][0]
                rows.append(dict(**snapshot,mean_mass_kg=float(mean),variance_mass_kg2=snapshot['variance'],
                    number_cdf_query_mass_kg=query,number_cdf=candidate,
                    reference_mean_mass_kg=ref['mean_mass_kg'],reference_variance_mass_kg2=ref['variance_mass_kg2'],
                    variance_ratio=snapshot['variance']/ref['variance_mass_kg2'] if ref['variance_mass_kg2'] else None,
                    mean_error_initial_mass=float((mean-ref['mean_mass_kg'])/mass_unit),
                    empirical_cdf_supremum=old.control.empirical_distance(query,candidate,sample,frequency)))
            extra=dict(pivots_kg=grid.pivots,bounds_kg=grid.bounds,event_capacity=setting['capacity'])
        run=dict(mass_cells=len(grid.pivots),grid_factor=factor,quadrature_order=order,max_step_s=dt,
            settings=setting,diagnostics=result['diagnostics'],snapshots=rows,wall_seconds=time.monotonic()-start,
            reconstruction_fits=grid.fit_diagnostics,reconstruction_roundoff=grid.roundoff_diagnostics,**extra)
        return dict(run=old.control.plain(run),failure=None)
    except (ArithmeticError,ValueError,RuntimeError) as exc:
        return dict(run=None,failure=dict(settings=setting,exception=type(exc).__name__,message=str(exc),wall_seconds=time.monotonic()-start))


def error_budget_checks(runs,case):
    checks=[]
    for run in runs:
        rows=[]
        for row in run['snapshots'][1:]:
            bound=(np.asarray(row['poisson_tail_moment_upper_bounds']) if case=='boundary' else
                   np.asarray(row['uniformization_tail_moment_upper_bounds'])+row['finite_ladder_moment_upper_bounds'])
            variance=row['expected_variance_mass_ratio'] if case=='boundary' else row['reference_variance_mass_kg2']
            mean=row['expected_mean_mass_ratio'] if case=='boundary' else row['reference_mean_mass_kg']
            # Bound the induced centered second moment, not just raw M2.
            centered=bound[2]+2*mean*bound[1]+mean**2*bound[0]
            rows.append(dict(time_s=row['time_s'],number_upper_bound=float(bound[0]),
                first_moment_relative_upper_bound=float(bound[1]/mean),centered_second_moment_bound_over_reference_variance=float(centered/variance),
                passed=bool(bound[0]<=1e-12 and bound[1]/mean<=1e-12 and centered/variance<=1e-4)))
        checks.append(dict(checks=rows,passed=all(r['passed'] for r in rows)))
    return checks


def semigroup_control(product):
    c=product['capture'];unit=c['initial_collector_mass_kg'];donor=c['captured_mass_kg']/unit
    grid=entropy.MomentGrid(old.reference.nested_grid(2)/unit);face=grid.bounds[grid.bounds>1][0]
    t=(face-1)/donor/100;initial=grid.from_particles([1.],[1.]);answers=[]
    query=np.linspace(face-8*donor*np.sqrt(100*t),face+8*donor*np.sqrt(100*t),4097)
    for parts in (1,2):
        op=jumps.ConstantCapture(grid,donor,100.)
        result=jumps.integrate(op,initial,[t],max_step_s=t/parts)
        row=result['snapshots'][0]
        answers.append(dict(subdivisions=parts,moments=row['moments'],variance=row['variance'],number_cdf=jumps.number_cdf(grid,row['state'],query),diagnostics=result['diagnostics']))
    return old.control.plain(dict(time_s=t,query_mass_ratio=query,runs=answers,
        number_cdf_difference=float(np.max(abs(answers[0]['number_cdf']-answers[1]['number_cdf']))),
        relative_variance_difference=float(abs(answers[0]['variance']/answers[1]['variance']-1)),
        reading_rule='One update versus two half updates at the crossing. Any difference includes intermediate reconstruction; the underlying constant-rate transition has the exact semigroup property.'))


def analyse(case,workers=1):
    _,binding=old.flakes.moon_column();before=identity(binding)
    reference=json.loads((ROOT/old.control.REFERENCE).read_text());old.check_reference(reference)
    settings=controls(case);runs=[];failures=[]
    jobs=[(case,s) for s in settings]
    if workers==1:
        results=map(run_control,jobs)
        for setting,result in zip(settings,results):
            print(case,setting,'completed' if result['run'] else result['failure'],flush=True)
            if result['run']:runs.append(result['run'])
            else:failures.append(result['failure'])
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for setting,result in zip(settings,pool.map(run_control,jobs)):
                print(case,setting,'completed' if result['run'] else result['failure'],flush=True)
                if result['run']:runs.append(result['run'])
                else:failures.append(result['failure'])
    semigroup=semigroup_control(reference) if case=='boundary' else None
    acceptance=old.boundary_decisions(runs) if case=='boundary' else old.bath_decisions(runs,reference)
    acceptance['accounting_and_cdf_checks']=old.accounting_checks(runs,case,reference)
    acceptance['truncation_budget_checks']=error_budget_checks(runs,case)
    acceptance['all_requested_controls_completed']=len(runs)==len(settings) and not failures
    acceptance['initial_cdf_checks']=[dict(settings=r['settings'],passed=r['snapshots'][0]['number_cdf_supremum' if case=='boundary' else 'empirical_cdf_supremum']<=1e-12) for r in runs]
    acceptance['semigroup_check']=semigroup is None or (semigroup['number_cdf_difference']<=.005 and semigroup['relative_variance_difference']<=1e-4)
    acceptance['passed']=(acceptance['passed'] and acceptance['all_requested_controls_completed'] and acceptance['semigroup_check']
        and all(c['passed'] for c in acceptance['initial_cdf_checks'])
        and all(c['passed'] for c in acceptance['accounting_and_cdf_checks']+acceptance['truncation_budget_checks']))
    if identity(binding)!=before:raise RuntimeError('Producer or inputs changed during finite-capture controls')
    result=dict(schema=SCHEMA,producer=before,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,workers=workers),
        case='poisson_boundary' if case=='boundary' else 'prescribed_bath',measure_name='bounded_quadratic_maximum_entropy',
        method='direct Poisson finite-time mixture' if case=='boundary' else 'state-dependent finite-time birth uniformization',
        requested_controls=settings,runs=runs,execution_failures=failures,acceptance=acceptance,semigroup_control=semigroup,
        evidence='Numerical maintained-bath capture controls. No closed aggregation, moisture depletion, full lunar column or empirical snow validation.',
        reading_rule='Distribution, moments, accounting and truncation are separate gates. Analytical tail and finite-ladder bounds exclude reconstruction, source quadrature and floating-point errors. Physical exits stop at first crossing. Poisson controls sample each grid at its own first face; the initial 191-cell setting matches the earlier crossing. The physical reference uses 8192 independent event histories and unbiased sample variance. Its 95% DKW band is uniform in mass at each time separately, not joint across controls or a numerical convergence bound.')
    if case=='bath':result['settings']={k:reference['capture'][k] for k in ('initial_collector_mass_kg','captured_mass_kg','initial_collector_diameter_m','captured_diameter_m','bath_number_m3','times_s')}
    return old.control.plain(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',choices=('boundary','bath'),required=True)
    parser.add_argument('--workers',type=int,choices=(1,2),default=1);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=analyse(args.case,args.workers)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Saved',args.output,'acceptance',result['acceptance']['passed'],flush=True)
    if not result['acceptance']['passed']:raise SystemExit(1)


if __name__=='__main__':main()
