"""Independent distribution controls for bounded entropy collision reconstruction.

The beta controls and collision accounting remain frozen. This runner repeats
their capture protocol with a new reconstruction and algebraically equivalent
SSPRK2 averaging that preserves stationary subnormal states. Both are bound.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.stats import poisson

from atmosphere.snow import collision_entropy as entropy
from atmosphere.snow import collision_moments as cm
from atmosphere.snow import collision_moment_benchmarks as control
from atmosphere.snow import collision_reference as reference, flakes
from shared.provenance import constants_used, constants_changed

ROOT = control.ROOT
FILES = ('atmosphere/snow/collision_entropy_benchmarks.py',
         'atmosphere/snow/collision_entropy.py', *control.FILES)
SCHEMA = 'terluna.atmosphere.snow-collision-entropy-benchmarks/1'


def identity(binding):
    result = control.identity(binding)
    result['files'] = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
    result['constants'] = constants_used(FILES)
    return result


# The two frozen capture protocols are repeated here so the repaired integrator
# is selected explicitly; frozen benchmark module globals are never patched.
def bath_case(kernel,product,factor=1,order=2,max_step_s=1.,grid_class=entropy.MomentGrid):
    capture=product['capture']; exact=capture['continuous_mass_references'][-1]
    grid=grid_class(reference.nested_grid(factor))
    initial=capture['initial_collector_mass_kg']; small=capture['captured_mass_kg']
    state=grid.from_particles([initial],[1.])
    operator=cm.BathCapture(grid,kernel,small,capture['bath_number_m3'],quadrature_order=order)
    started=time.monotonic()
    result=entropy.integrate(operator,state,capture['times_s'],max_step_s=max_step_s,time_limit_s=600.)
    rows=[]
    for snapshot,target in zip(result['snapshots'],exact['snapshots']):
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
                 empirical_cdf_supremum=control.empirical_distance(query,cdf,sample,counts))
        rows.append(row)
    return dict(mass_cells=len(grid.pivots),grid_factor=factor,quadrature_order=order,
                max_step_s=max_step_s,pivots_kg=grid.pivots,bounds_kg=grid.bounds,
                diagnostics=result['diagnostics'],wall_seconds=time.monotonic()-started,snapshots=rows)

def boundary_case(product,*,grid_class=entropy.MomentGrid,factor=2,order=2,max_step_s=.5):
    """Analytic Poisson distribution while its mean crosses a fixed cell face.

    Masses are normalized by the imposed lunar collector, not actual kilograms.
    The rate is prescribed 100/s and independent of mass. Grid controls retain
    the same physical sampling times defined by the 191-cell reference face.
    """
    capture=product['capture'];mass_unit=capture['initial_collector_mass_kg']
    donor=capture['captured_mass_kg']/mass_unit
    reference_grid=cm.MomentGrid(reference.nested_grid(2)/mass_unit)
    face=float(reference_grid.bounds[reference_grid.bounds>1][0])
    expected=(face-1)/donor
    counts=np.r_[0.,expected+np.array([-2.,0.,2.])*np.sqrt(expected),1.5*expected]
    rate=100.;times=counts/rate
    grid=grid_class(reference.nested_grid(factor)/mass_unit)
    initial=grid.from_particles([1.],[1.]);started=time.monotonic();last=[started]
    def progress(t,diagnostics):
        now=time.monotonic()
        if now-last[0]>30:
            print(f'boundary t={t:.4f}s, steps={diagnostics["steps"]}',flush=True);last[0]=now
    result=entropy.integrate(cm.BathCapture(grid,lambda a,b:1.,donor,rate,quadrature_order=order),
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


def check_reference(product):
    for kind in ('files', 'inputs'):
        for path, expected in product['producer'][kind].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
                raise RuntimeError('Frozen reference identity mismatch: '+path)
    if constants_changed(product['producer'].get('constants', {})):
        raise RuntimeError('Frozen reference constants changed')


def boundary_decisions(runs):
    gates = dict(maximum_number_cdf_supremum=.005,
                 maximum_relative_variance_error=1e-4,
                 maximum_absolute_mean_error=1e-10)
    decisions = []
    for run in runs:
        rows = [dict(time_s=r['time_s'],
                     cdf=r['number_cdf_supremum'] <= gates['maximum_number_cdf_supremum'],
                     variance=abs(r['variance_mass_ratio']/r['expected_variance_mass_ratio']-1)
                              <= gates['maximum_relative_variance_error'],
                     mean=abs(r['mean_mass_ratio']-r['expected_mean_mass_ratio'])
                          <= gates['maximum_absolute_mean_error'])
                for r in run['snapshots'][1:]]
        decisions.append(dict(cells=run['mass_cells'], order=run['quadrature_order'],
                              dt=run['max_step_s'], checks=rows,
                              passed=all(all(r[k] for k in ('cdf','variance','mean')) for r in rows)))
    return dict(gates=gates, runs=decisions, passed=bool(runs) and all(r['passed'] for r in decisions))


def bath_decisions(runs, product):
    target = product['capture']['continuous_mass_references'][-1]
    band = target['number_cdf_uniform_95_percent_dkw_halfwidth']
    gates = dict(mean_maximum_reference_standard_errors=5., variance_ratio_tolerance=.05,
                 number_cdf_maximum_distance=band)
    decisions = []
    for run in runs:
        rows = [dict(time_s=row['time_s'],
                     mean=abs(row['mean_mass_kg']-ref['mean_mass_kg']) <= 5*ref['mean_mass_standard_error_kg'],
                     variance=abs(row['variance_ratio']-1) <= .05,
                     cdf=row['empirical_cdf_supremum'] <= band)
                for row,ref in zip(run['snapshots'][1:],target['snapshots'][1:])]
        decisions.append(dict(cells=run['mass_cells'], order=run['quadrature_order'],
                              dt=run['max_step_s'], checks=rows,
                              passed=all(all(r[k] for k in ('cdf','variance','mean')) for r in rows)))
    return dict(gates=gates, runs=decisions, passed=bool(runs) and all(r['passed'] for r in decisions))


def accounting_checks(runs, case, product):
    """Conservation and CDF admission, separate from shape agreement."""
    mass_unit = 1. if case == 'boundary' else product['capture']['initial_collector_mass_kg']
    checks = []
    for run in runs:
        rows = []
        for row in run['snapshots']:
            query = np.asarray(row['query_mass_ratio'] if case == 'boundary' else row['number_cdf_query_mass_kg'])
            cdf = np.asarray(row['number_cdf'])
            n,mass,_ = row['moments']
            exits = np.asarray(row['ledgers']['lower_exit'])+row['ledgers']['upper_exit']
            rows.append(dict(time_s=row['time_s'],
                number=abs(n+exits[0]-1.) <= 1e-10 and abs(row['number_residual']) <= 1e-10,
                mass=abs(mass+exits[1]-row['ledgers']['supplied_mass']-mass_unit)/mass_unit <= 1e-10
                     and abs(row['mass_residual'])/mass_unit <= 1e-10,
                cdf_valid=bool(np.isfinite(cdf).all() and (np.diff(query)>0).all()
                               and np.min(np.diff(cdf)) >= -1e-10 and cdf.min() >= -1e-10 and cdf.max() <= 1+1e-10),
                omitted_poisson_probability=case != 'boundary' or
                    max(row['unqueried_left_probability'],row['unqueried_right_probability']) <= 1e-20))
        checks.append(dict(checks=rows,passed=all(all(r[k] for k in
            ('number','mass','cdf_valid','omitted_poisson_probability')) for r in rows)))
    return checks


def analyse(case, *, short=False):
    world,binding = flakes.moon_column()
    before = identity(binding)
    product = json.loads((ROOT/control.REFERENCE).read_text())
    check_reference(product)
    settings = ([(2,2,.5),(2,2,.25),(2,4,.5),(1,2,.5),(4,2,.5)] if case == 'boundary'
                else [(1,2,1.),(2,2,1.),(4,2,1.),(2,4,1.),(2,2,.5)])
    selected = settings[:1] if short else settings
    runs = []
    failure = None
    kernel = reference.SnowKernel(world)
    for factor,order,dt in selected:
        started = time.monotonic()
        created_grids = []
        def make_grid(pivots):
            grid = entropy.MomentGrid(pivots)
            created_grids.append(grid)
            return grid
        print(f'{case}: grid factor={factor}, quadrature={order}, maximum step={dt}s',flush=True)
        try:
            if case == 'boundary':
                run = boundary_case(product,grid_class=make_grid,
                                            factor=factor,order=order,max_step_s=dt)
            else:
                run = bath_case(kernel,product,factor,order,dt,grid_class=make_grid)
        except (ArithmeticError, ValueError, RuntimeError) as exc:
            failure = dict(grid_factor=factor,quadrature_order=order,max_step_s=dt,
                           exception=type(exc).__name__,message=str(exc),wall_seconds=time.monotonic()-started)
            if created_grids:
                failure['reconstruction_roundoff'] = dict(created_grids[0].roundoff_diagnostics)
                failure['reconstruction_fits'] = dict(created_grids[0].fit_diagnostics)
            print('Incomplete control: '+str(exc),flush=True)
            break
        run['reconstruction_roundoff'] = dict(created_grids[0].roundoff_diagnostics)
        run['reconstruction_fits'] = dict(created_grids[0].fit_diagnostics)
        runs.append(run)
        metric = 'number_cdf_supremum' if case == 'boundary' else 'empirical_cdf_supremum'
        print(f'Completed in {run["wall_seconds"]:.2f}s; CDF errors '+
              ', '.join(f'{r[metric]:.6g}' for r in run['snapshots']),flush=True)
    if identity(binding) != before:
        raise RuntimeError('Producer or inputs changed during entropy controls; no product admitted')
    acceptance = boundary_decisions(runs) if case == 'boundary' else bath_decisions(runs,product)
    acceptance['accounting_and_cdf_checks'] = accounting_checks(runs,case,product)
    acceptance['all_requested_controls_completed'] = failure is None and len(runs) == len(selected)
    acceptance['passed'] = (acceptance['passed'] and acceptance['all_requested_controls_completed']
                            and all(r['passed'] for r in acceptance['accounting_and_cdf_checks']))
    result = dict(schema=SCHEMA,producer=before,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),case='poisson_boundary' if case=='boundary' else 'prescribed_bath',
                  measure_name='bounded_quadratic_maximum_entropy',
                  control_scope='first control only; refinement not established' if short else 'grid, timestep and quadrature controls',
                  requested_controls=[dict(grid_factor=f,quadrature_order=q,max_step_s=t) for f,q,t in selected],
                  runs=runs,execution_failure=failure,acceptance=acceptance,
                  evidence=('Analytic constant-rate capture control using frozen event accounting and a new bounded within-cell reconstruction. The constant kernel is a mathematical test, not the physical snow law.' if case == 'boundary' else 'Numerical comparison with the frozen physical snow kernel and event accounting. The changes are bounded within-cell reconstruction and an algebraically equivalent, subnormal-stable SSPRK2 average. No physical size prediction or empirical validation.'),
                  reading_rule='Moment inversion, event accounting and distribution accuracy are separate checks. Reconstruction rounding bounds apply at validation time, including rejected trials; they do not bound subsequent nonlinear amplification. The SSPRK2 average preserves identical subnormal states algebraically. Smooth three-moment closures cannot reproduce arbitrary multimodal populations or discrete low-count capture. No full column, deposition, breakup, settling or moisture-history repair is exercised. Incomplete controls remain failures to establish adequacy.')
    if case == 'boundary':
        result['reading_rule'] += (' Mass is normalized by the imposed collector; rate 100/s is a mathematical control. '
            'Times remain fixed at the 191-cell reference boundary for grid comparisons. The exact Poisson CDF is tested '
            'on both sides of every atom within 12 count standard deviations, with omitted probability reported. No sampling band is used.')
    else:
        capture = product['capture']
        result['settings'] = {k:capture[k] for k in ('initial_collector_mass_kg','captured_mass_kg',
            'initial_collector_diameter_m','captured_diameter_m','bath_number_m3','times_s')}
        result['reading_rule'] += (' Number CDFs use original collector normalization. The reference variance is the '
            'unbiased sample variance of 8192 independent histories; candidate variance is its represented population variance. '
            'The DKW band is 95% uniform in mass at each saved time separately, not a joint 95% guarantee across times or controls. It describes sampling only and supplies endpoint compatibility, not closure convergence.')
    return control.plain(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=('boundary','bath'),default='boundary')
    parser.add_argument('--short',action='store_true')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = analyse(args.case,short=args.short)
    output = args.output or ROOT/f'atmosphere/snow/results/collision_entropy_{args.case}.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Wrote '+str(output),flush=True)
    if result['execution_failure'] is not None:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
