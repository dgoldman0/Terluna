"""Check continuous within-cell collision moments against frozen references.

This producer changes the numerical collision representation only. It does not
change the imported snow kernel, or infer a physical lunar size distribution.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
from scipy.stats import poisson
import numpy as np
from atmosphere.snow import collision_reference as ref, flakes
from atmosphere.snow import collision_moments as cm
from shared.provenance import constants_used

ROOT=Path(__file__).resolve().parents[2]
REFERENCE='atmosphere/snow/results/collision_reference.json'
SCHEMA='terluna.atmosphere.snow-collision-moment-benchmarks/1'
FILES=('atmosphere/snow/collision_moment_benchmarks.py',
       'atmosphere/snow/collision_moments.py',*ref.FILES)

def identity(binding):
    sha=lambda p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    return dict(files={p:sha(p) for p in FILES},
                inputs={REFERENCE:sha(REFERENCE),**binding['inputs']},
                constants=constants_used(FILES))

def empirical_distance(query,cdf,sample_masses,sample_counts):
    """Supremum at either side of each empirical jump; absolute number CDF."""
    total=np.sum(sample_counts)
    empirical=np.cumsum(sample_counts)/total
    candidate=np.interp(sample_masses,query,cdf)
    return float(max(np.max(np.abs(candidate-empirical)),
                     np.max(np.abs(np.interp(np.nextafter(sample_masses,-np.inf),query,cdf)-np.r_[0.,empirical[:-1]]))))

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x

def bath_case(kernel,product,factor=1,order=2,max_step_s=1.,grid_class=cm.MomentGrid):
    capture=product['capture']; reference=capture['continuous_mass_references'][-1]
    grid=grid_class(ref.nested_grid(factor))
    initial=capture['initial_collector_mass_kg']; small=capture['captured_mass_kg']
    state=grid.from_particles([initial],[1.])
    operator=cm.BathCapture(grid,kernel,small,capture['bath_number_m3'],quadrature_order=order)
    started=time.monotonic()
    result=cm.integrate(operator,state,capture['times_s'],max_step_s=max_step_s,time_limit_s=600.)
    rows=[]
    for snapshot,target in zip(result['snapshots'],reference['snapshots']):
        target_mean=target['mean_mass_kg'];target_variance=target['variance_mass_kg2']
        sd=target_variance**.5
        sample=np.asarray(target['masses_kg']); counts=np.asarray(target['sample_counts'])
        query=np.unique(np.r_[sample,np.nextafter(sample,-np.inf),grid.bounds.ravel(),
            np.linspace(target_mean-6*sd,target_mean+6*sd,1025),initial*np.linspace(.8,1.8,201)])
        cdf=grid.cdf(snapshot['state'],query,weighting='number')*snapshot['moments'][0]
        n,m1,m2=snapshot['moments'];mean=m1/n;variance=snapshot['variance']
        row=dict(time_s=snapshot['time_s'],state=snapshot['state'],moments=snapshot['moments'],
                 ledgers=snapshot['ledgers'],number_residual=snapshot['number_residual'],mass_residual=snapshot['mass_residual'],mean_mass_kg=mean,variance_mass_kg2=variance,
                 number_cdf_query_mass_kg=query,number_cdf=cdf,
                 reference_mean_mass_kg=target_mean,reference_variance_mass_kg2=target_variance,
                 mean_error_initial_mass=(mean-target_mean)/initial,
                 variance_ratio=variance/target_variance if snapshot['time_s']>0 and target_variance else None,
                 empirical_cdf_supremum=empirical_distance(query,cdf,sample,counts))
        rows.append(row)
    return dict(mass_cells=len(grid.pivots),grid_factor=factor,quadrature_order=order,
                max_step_s=max_step_s,pivots_kg=grid.pivots,bounds_kg=grid.bounds,
                diagnostics=result['diagnostics'],wall_seconds=time.monotonic()-started,snapshots=rows)

def analyse_bath(*,short=False):
    world,binding=flakes.moon_column();before=identity(binding)
    product=json.loads((ROOT/REFERENCE).read_text())
    for kind in ('files','inputs'):
        for path,sha in product['producer'][kind].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=sha:
                raise RuntimeError('Frozen reference identity mismatch: '+path)
    kernel=ref.SnowKernel(world);runs=[]
    settings=[(1,2,1.),(2,2,1.),(4,2,1.),(2,4,1.),(2,2,.5)]
    for factor,order,dt in (settings[:2] if short else settings):
        run=bath_case(kernel,product,factor,order,dt);runs.append(run)
        print(f'bath cells={run["mass_cells"]}, q={order}, dt={dt}: '+
              f'{run["wall_seconds"]:.2f}s; final variance ratio='+
              str(run['snapshots'][-1]['variance_ratio']),flush=True)
    capture=product['capture'];reference=capture['continuous_mass_references'][-1]
    gates=dict(description='Predeclared comparison gates for high-count bath endpoints; not empirical validation.',
               mean_maximum_reference_standard_errors=5.,variance_ratio_tolerance=.05,
               number_cdf_maximum_distance=reference['number_cdf_uniform_95_percent_dkw_halfwidth'])
    decisions=[]
    for run in runs:
        checks=[]
        for row,target in zip(run['snapshots'][1:],reference['snapshots'][1:]):
            checks.append(dict(time_s=row['time_s'],
                mean=abs(row['mean_mass_kg']-target['mean_mass_kg'])<=5*target['mean_mass_standard_error_kg'],
                variance=abs(row['variance_ratio']-1)<=.05,
                cdf=row['empirical_cdf_supremum']<=gates['number_cdf_maximum_distance']))
        decisions.append(dict(cells=run['mass_cells'],order=run['quadrature_order'],dt=run['max_step_s'],checks=checks,
                              passed=all(all(r[k] for k in ('mean','variance','cdf')) for r in checks)))
    if identity(binding)!=before:raise RuntimeError('Sources changed during benchmark; result not admitted')
    return plain(dict(schema=SCHEMA,producer=before,case='prescribed_bath',measure_name='bounded_beta',control_scope='96/191-cell endpoints only; additional endpoint controls withheld after independent boundary-shape failure' if short else 'five endpoint controls',settings=dict(
        initial_collector_mass_kg=capture['initial_collector_mass_kg'],captured_mass_kg=capture['captured_mass_kg'],
        initial_collector_diameter_m=capture['initial_collector_diameter_m'],captured_diameter_m=capture['captured_diameter_m'],
        bath_number_m3=capture['bath_number_m3'],times_s=capture['times_s']),runs=runs,acceptance=dict(gates=gates,runs=decisions),
        evidence='Numerical collision closure comparison under the unchanged prescribed-bath snow kernel. Reference histories are independent exact events; candidate uses bounded within-cell beta reconstruction.',
        reading_rule='Number CDFs retain their original collector normalization. Reference DKW bands describe independent sampling only; matching them is an endpoint compatibility check, not proof of closure accuracy or tail convergence. The smooth closure cannot exactly represent low-count discrete capture distributions. No lunar column, moisture budget or aggregate fracture law is changed.'))

def box_case(kernel, product, factor=1, order=2, max_step_s=30.,grid_class=cm.MomentGrid):
    ordinary=product['ordinary_box']
    times=np.asarray(ordinary['times_s'][:-1])  # six hours; late reference is unresolved
    grid=grid_class(ref.nested_grid(factor))
    masses=ref.ordinary_initial()
    initial=grid.from_particles(masses,np.full(len(masses),1e5/len(masses)))
    operator=cm.CollisionMoments(grid,kernel,quadrature_order=order)
    started=time.monotonic()
    result=cm.integrate(operator,initial,times,max_step_s=max_step_s)
    target=ordinary['exact_event_ensembles'][-1]
    query_d=np.asarray(ordinary['cdf_query_diameter_m'])
    query_m=ref.GEOMETRY.mass(query_d)
    rows=[]
    for row,reference in zip(result['snapshots'],target['snapshots']):
        cdf=grid.cdf(row['state'],query_m,weighting='mass')*row['moments'][1]/grid.moments(initial)[1]
        values=dict(zip(('number_m3','mass_kg_m3','second_mass_moment_kg2_m3'),row['moments']))
        comparison={key:dict(candidate=float(value),reference_mean=reference['mean'][key],
            reference_standard_error=reference['standard_error'][key],
            relative_difference=float(value/reference['mean'][key]-1)) for key,value in values.items()}
        rows.append(dict(time_s=row['time_s'],state=row['state'],moments=row['moments'],
            ledgers=row['ledgers'],number_residual=row['number_residual'],mass_residual=row['mass_residual'],mass_cdf=cdf,comparison=comparison,
            maximum_absolute_mass_cdf_difference=float(np.max(np.abs(cdf-reference['mass_cdf_mean'])))))
    return dict(mass_cells=len(grid.pivots),grid_factor=factor,quadrature_order=order,
        max_step_s=max_step_s,pivots_kg=grid.pivots,bounds_kg=grid.bounds,
        diagnostics=result['diagnostics'],wall_seconds=time.monotonic()-started,snapshots=rows)

def analyse_box():
    world,binding=flakes.moon_column();before=identity(binding)
    product=json.loads((ROOT/REFERENCE).read_text())
    kernel=ref.SnowKernel(world);runs=[]
    for factor,order,dt in ((1,2,30.),(2,2,30.),(4,2,30.),(2,4,30.),(2,2,15.)):
        run=box_case(kernel,product,factor,order,dt);runs.append(run)
        print(f'box cells={run["mass_cells"]}, q={order}, dt={dt}: '+
              f'{run["wall_seconds"]:.2f}s',flush=True)
    if identity(binding)!=before:raise RuntimeError('Sources changed during benchmark; result not admitted')
    return plain(dict(schema=SCHEMA,producer=before,case='ordinary_box',measure_name='bounded_beta',
        settings=dict(temperature_c=-10.,initial_number_m3=1e5,aggregation=True,
            chipping=False,deposition=False,settling=False,
            initial_masses_kg=ref.ordinary_initial(),times_s=product['ordinary_box']['times_s'][:-1],
            cdf_query_diameter_m=product['ordinary_box']['cdf_query_diameter_m']),runs=runs,
        evidence='Deterministic within-cell collision closure compared with frozen independent finite-volume exact-event ensembles under identical discrete seed masses and collision law.',
        reading_rule='The largest-volume reference is not a proven infinite-volume solution. Standard errors use independent entire realizations. Earlier six-hour bulk convergence provides a diagnostic target; the twelve-hour reference remains unresolved and is excluded from acceptance. No deposition, breakup or column is tested here.'))

def boundary_case(product,*,grid_class=cm.MomentGrid,factor=2,order=2,max_step_s=.5):
    """Analytic Poisson distribution while its mean crosses a fixed cell face.

    Masses are normalized by the imposed lunar collector, not actual kilograms.
    The rate is prescribed100/s and independent of mass. Grid controls retain
    the same physical sampling times defined by the191-cell reference face.
    """
    capture=product['capture'];mass_unit=capture['initial_collector_mass_kg']
    donor=capture['captured_mass_kg']/mass_unit
    reference_grid=cm.MomentGrid(ref.nested_grid(2)/mass_unit)
    face=float(reference_grid.bounds[reference_grid.bounds>1][0])
    expected=(face-1)/donor
    counts=np.r_[0.,expected+np.array([-2.,0.,2.])*np.sqrt(expected),1.5*expected]
    rate=100.;times=counts/rate
    grid=grid_class(ref.nested_grid(factor)/mass_unit)
    initial=grid.from_particles([1.],[1.]);started=time.monotonic();last=[started]
    def progress(t,diagnostics):
        now=time.monotonic()
        if now-last[0]>30:
            print(f'boundary t={t:.4f}s, steps={diagnostics["steps"]}',flush=True);last[0]=now
    result=cm.integrate(cm.BathCapture(grid,lambda a,b:1.,donor,rate,quadrature_order=order),
        initial,times,max_step_s=max_step_s,time_limit_s=600.,progress=progress)
    rows=[]
    for row,count in zip(result['snapshots'],counts):
        k=np.arange(max(0,int(count-12*np.sqrt(count))),int(count+12*np.sqrt(count))+2)
        mass=1+donor*k;left=np.nextafter(mass,-np.inf)
        query=np.sort(np.r_[mass,left]);n=row['moments'][0]
        candidate=grid.cdf(row['state'],query)*n
        exact=np.empty_like(query);exact[::2]=poisson.cdf(k-1,count);exact[1::2]=poisson.cdf(k,count)
        supremum=float(np.max(np.abs(candidate-exact)))
        mean=row['moments'][1]/n;variance=row['variance'];expected_variance=count*donor*donor
        rows.append(dict(time_s=row['time_s'],poisson_mean=float(count),state=row['state'],moments=row['moments'],
            mean_mass_ratio=float(mean),variance_mass_ratio=float(variance),expected_mean_mass_ratio=float(1+donor*count),
            expected_variance_mass_ratio=float(expected_variance),
            query_mass_ratio=query,number_cdf=candidate,exact_poisson_number_cdf=exact,
            number_cdf_supremum=supremum,unqueried_left_probability=float(poisson.cdf(k[0]-1,count)),
            unqueried_right_probability=float(poisson.sf(k[-1],count)),
            ledgers=row['ledgers'],number_residual=row['number_residual'],mass_residual=row['mass_residual']))
    return dict(mass_cells=len(grid.pivots),grid_factor=factor,quadrature_order=order,max_step_s=max_step_s,
        pivots_mass_ratio=grid.pivots,bounds_mass_ratio=grid.bounds,donor_mass_ratio=donor,
        prescribed_capture_rate_s=rate,reference_cell_face_mass_ratio=face,
        diagnostics=result['diagnostics'],snapshots=rows,wall_seconds=time.monotonic()-started)

def analyse_boundary():
    _,binding=flakes.moon_column();before=identity(binding)
    product=json.loads((ROOT/REFERENCE).read_text())
    run=boundary_case(product)
    gates=dict(maximum_number_cdf_supremum=.005,maximum_relative_variance_error=1e-4,
               maximum_absolute_mean_error=1e-10)
    decision=[]
    for row in run['snapshots'][1:]:
        decision.append(dict(time_s=row['time_s'],cdf=row['number_cdf_supremum']<=.005,
            variance=abs(row['variance_mass_ratio']/row['expected_variance_mass_ratio']-1)<=1e-4,
            mean=abs(row['mean_mass_ratio']-row['expected_mean_mass_ratio'])<=1e-10))
    if identity(binding)!=before:raise RuntimeError('Sources changed during boundary benchmark')
    return plain(dict(schema=SCHEMA,producer=before,case='poisson_boundary',measure_name='bounded_beta',runs=[run],
        acceptance=dict(gates=gates,decisions=decision,passed=all(all(r[k] for k in ('cdf','variance','mean')) for r in decision)),
        evidence='Analytic constant-rate Poisson capture control. The numerical closure is tested before, at and after a fixed cell face; no Monte Carlo sampling uncertainty enters the exact CDF.',
        reading_rule='Mass is normalized by the imposed collector. This mathematical control changes the kernel to a constant rate and is not a lunar storm prediction. The CDF supremum is evaluated at both sides of every Poisson jump within12count standard deviations; omitted analytic probability is reported. Correct mean and variance do not clear a failed distribution-shape gate.'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=('bath','box','boundary'),default='bath')
    parser.add_argument('--short',action='store_true',help='Only96/191-cell bath endpoints after the boundary failure')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=analyse_bath(short=args.short) if args.case=='bath' else (analyse_box() if args.case=='box' else analyse_boundary())
    output=args.output or ROOT/f'atmosphere/snow/results/collision_moment_{args.case}.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Wrote '+str(output),flush=True)
if __name__=='__main__':main()
