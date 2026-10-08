"""Seek missed finite-Sun service rays after the complete mass-feedback replay."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize

from shared.provenance import constants_used
from protection.dynamics.service_margin import service_margin
from .packing_screen import parents,load_raw
from .packing_run import full_command
from .patterns_run import reference_shooting
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE


def main():
    before=time.monotonic();names=['packing_validation.json','closure_screen.json','closure_hardware.json',
        'closure_sizing.json','closure_return.json','closure_validation.json','closure_budget.json']
    p,sources=parents(names);used=sum(p[n]['elapsed_s'] for n in names[1:])
    allowance=min(150,int(4800-used))
    if allowance<=0:raise ValueError('No audit budget remains')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    for path in ['protection/dynamics/service_margin.py','research/studies/solar_shield_array/closure_audit.py']:
        sources[path]=digest(ROOT/path)
    env=environment(1.);command=full_command(env,p['packing_validation.json']['schedule'])
    reference,_=reference_shooting(env,21600.)
    raw=load_raw(p['closure_validation.json']['runs'][0])
    spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:]);probes=[]
    def evaluate(x):
        hour,angle,rho=x;t=hour*3600.
        data=service_margin(np.c_[spline(t),spline(t,1)],command(t).as_matrix(),env.at(t),reference.sol(t)[:6],
                            rho*np.array([np.cos(angle),np.sin(angle)]))
        probes.append(dict(time_s=float(t),angle_rad=float(angle),solar_radius_fraction=float(rho),**data))
        return data['margin_m']
    for hour in np.linspace(0.,6.,13):
        evaluate([hour,0.,0.])
        for k in range(32):evaluate([hour,.137+np.pi/32+k*np.pi/16,1.])
        for k in range(8):evaluate([hour,.137+np.pi/16+k*np.pi/4,.5])
    grid=len(probes);seen=set();seeds=[]
    for row in sorted(probes,key=lambda r:r['margin_m']):
        sector=int((row['angle_rad']%(2*np.pi))/(np.pi/2))
        if sector not in seen:seeds.append(row);seen.add(sector)
        if len(seeds)==4:break
    fits=[]
    for row in seeds:
        x=np.array([row['time_s']/3600.,row['angle_rad'],row['solar_radius_fraction']]);simplex=np.tile(x,(4,1))
        simplex[1,0]+=.1 if x[0]<5.9 else -.1;simplex[2,1]+=np.pi/32;simplex[3,2]+=.05 if x[2]<.95 else -.05
        fit=minimize(evaluate,x,method='Nelder-Mead',bounds=[(0.,6.),(x[1]-np.pi/8,x[1]+np.pi/8),(0.,1.)],
            options={'maxfev':160,'xatol':1e-5,'fatol':1e-3,'initial_simplex':simplex})
        fits.append(dict(success=bool(fit.success),nfev=int(fit.nfev),minimum_margin_m=float(fit.fun)))
    worst=min(probes,key=lambda r:r['margin_m'])
    out=dict(schema='terluna.research.cycle-hardware-ray-audit/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,this_stage_wall_s=allowance),grid_probes=grid,total_probes=len(probes),
        fits=fits,worst=worst,passed=bool(worst['margin_m']>0 and min(r['covered_fraction'] for r in probes)>=1-1e-9),
        accepted_fleet=False,continuous_certificate=False,
        scope='Whole circular receiver for each tested source/date, with circle sagitta allowance; bounded searches over solar-disk position and service date. No global source/time or full-target certificate.',
        completed=True,elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    (HERE/'results/closure_audit.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
