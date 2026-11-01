"""One declared compact-sphere capacitance sensitivity, not archive calibration.

Set C=D/2 only below the source's solid-sphere mass threshold. The frozen public
source evaluates C from a floating residual and picks D/pi in tested small cases.
This explicit geometric alternative is not adopted into any production model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from atmosphere.snow import welss_column as w
from atmosphere.snow import welss_comparison as suite
from atmosphere.snow import welss_archive_comparison as archive

HERE=Path(__file__).resolve().parent


def compact_geometry(mass,original):
    result=original(mass)
    result['C']=np.where(np.asarray(mass)<w.M_THRESHOLD,result['D']/2,result['C'])
    return result


def domain_checked_rhs(fun,n,counter):
    """Reject only unrepresentable internal trial states, never clip a mass."""
    def checked(z,y):
        with np.errstate(over='ignore',under='ignore',invalid='ignore'):
            mass=np.exp(y[:n])
        if not np.isfinite(y).all() or not np.isfinite(mass).all() or (mass<=0).any():
            counter['invalid_rk_stages']+=1
            return np.full_like(y,np.nan)
        return fun(z,y)
    return checked


def startup_limited_solver(original,fun,t_span,y0,*,counter=None,first_step_factor=1.,**kwargs):
    # A faster compact-seed growth rate can defeat SciPy's automatic first
    # log-mass step. Bound only the initial trial increment; later steps adapt.
    n=len(y0)//2
    counter=counter if counter is not None else dict(invalid_rk_stages=0)
    checked=domain_checked_rhs(fun,n,counter)
    log_rate=float(np.max(np.abs(checked(t_span[0],y0)[:n])))
    kwargs.setdefault('first_step',first_step_factor*min(1e-4,1e-3/max(log_rate,1e-300)))
    return original(checked,t_span,y0,**kwargs)


def integrate_with_intervention(*,first_step_factor=1.,**kwargs):
    original=w.geometry
    original_solver=w.solve_ivp
    counter=dict(invalid_rk_stages=0)
    try:
        w.geometry=lambda mass:compact_geometry(mass,original)
        w.solve_ivp=lambda fun,t_span,y0,**kw:startup_limited_solver(original_solver,fun,t_span,y0,counter=counter,first_step_factor=first_step_factor,**kw)
        result=w.integrate(**kwargs)
        result['diagnostics'].update(counter)
        return result
    finally:
        w.geometry=original
        w.solve_ivp=original_solver


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--physical-result',type=Path,default=HERE/'results/welss_column.json')
    parser.add_argument('--output',type=Path,default=HERE/'results/welss_compact_capacitance.json')
    args=parser.parse_args()
    paths=[Path(__file__),args.reference,args.physical_result,HERE/'welss_archive_comparison.py']
    identity={str(p):archive.sha(p) for p in paths};source_before=suite.identity()
    physical=json.loads(args.physical_result.read_text())
    if physical['producer']!=source_before:
        raise ValueError('baseline physical producer identity differs')
    reference=w.read_archive(args.reference)
    settings=dict(suite.cases()['base'],rtol=1e-9)
    result=integrate_with_intervention(**settings,progress=True,wall_limit_s=600.)
    if identity!={str(p):archive.sha(p) for p in paths} or source_before!=suite.identity():
        raise RuntimeError('intervention source or input changed during calculation')
    z=result['height'];rows={};profiles={};baseline=physical['cases']['base']
    for mode,profile in result['profiles'].items():
        profiles[mode]={key:profile[i].tolist() for i,key in enumerate(w.MOMENTS)}
        rows[mode]={key:archive.profile_metrics(profile[i],reference['arrays'][key],z)
                    for i,key in enumerate(w.MOMENTS[:4])}
    product=dict(schema='terluna.atmosphere.welss-compact-capacitance/1',
                 producer=dict(file='atmosphere/snow/welss_compact_capacitance.py',sha256=identity[str(Path(__file__))],
                               frozen_source=source_before,comparison_producer_sha256=identity[str(HERE/'welss_archive_comparison.py')]),
                 baseline_physical_product=dict(file='atmosphere/snow/results/welss_column.json',sha256=identity[str(args.physical_result)]),
                 reference={k:v for k,v in reference.items() if k!='arrays'},
                 evidence='One geometric capacitance sensitivity of a deterministic Earth source experiment, not a reconstruction of the archived executable or empirical validation.',
                 intervention=dict(condition='mass < M_THRESHOLD',threshold_kg=w.M_THRESHOLD,
                                   replacement='capacitance = maximum_dimension / 2',
                                   other_terms='All source geometry, fall speed, atmosphere, ventilation and source distribution unchanged.',
                                   selection='Compact-sphere electrostatic capacitance, selected by the existing solid-sphere transition; no fit to archived profiles.'),
                 numerical_settings=dict(result['settings'],first_step_rule='min(1e-4m,1e-3/max initial absolute d(log mass)/dz); subsequent adaptive controls unchanged',invalid_trial_rule='Unrepresentable internal log-mass trial states return nonfinite RHS so DOP853 rejects and shrinks; no accepted state or physical mass is clipped'),diagnostics=result['diagnostics'],height_m=z.tolist(),
                 units=baseline['units'],profiles=profiles,comparison=rows,
                 remaining='One intervention with baseline quadratures and tighter ODE tolerance; its full quadrature convergence is not separately established and no closure adoption follows.')
    args.output.write_text(json.dumps(product,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows['last5h'],indent=2),flush=True)


if __name__=='__main__':
    main()
