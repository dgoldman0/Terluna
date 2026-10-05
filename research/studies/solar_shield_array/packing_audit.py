"""Independent all-pair surface-distance account for the replayed packing."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize

from shared.provenance import constants_used
from protection.dynamics.square_distance import distance_guard
from protection.dynamics.service_margin import service_margin
from .packing_screen import parents,load_raw
from .packing_run import full_command
from .patterns_run import reference_shooting
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE


def main():
    before=time.monotonic()
    names=['packing_screen.json','packing.json','packing_validation.json','packing_budget.json']
    p,sources=parents(names)
    used=sum(q['elapsed_s'] for q in p.values())
    restart=ROOT/'research/runs/solar_shield_array/packing_replay_restart.json'
    if restart.exists():used+=json.loads(restart.read_text()).get('charged_interrupted_wall_s',0)
    allowance=min(300,int(3600-used))
    if allowance<=0:raise ValueError('All-pair audit has no remaining numerical budget')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    validation=p['packing_validation.json'];env=environment(1.)
    if not validation['accepted_service_departure']:raise ValueError('No accepted replay')
    for path in ['protection/dynamics/square_distance.py','protection/dynamics/service_margin.py',
                 'research/studies/solar_shield_array/packing_audit.py']:
        sources[path]=digest(ROOT/path)
    out=dict(schema='terluna.research.packing-distance-audit/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_numerical_wall_s=used,this_stage_wall_s=allowance),
        scope='Exact closest surface distance for all parallel 10 km squares at two-second dates; interval guard assumes stated acceleration/rate bounds. Independent source/integrator replay is recorded separately. No full-model continuous certificate.',
        accepted_fleet=False,continuous_certificate=False,runs=[])
    command=full_command(env,validation['schedule'])
    for entry in validation['runs']:
        raw=load_raw(entry);interpolation=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        chunks=[]
        for start in np.arange(raw['t'][0],raw['t'][-1],3600.):
            times=np.arange(start,min(start+3600.,raw['t'][-1])+1,2.)
            state=np.concatenate([interpolation(times),interpolation(times,1)],axis=-1)
            chunks.append(distance_guard(times,state,command,point_error=entry['fine_interpolation_position_error_m']))
        minimum=min(chunks,key=lambda r:r['minimum_sampled_surface_distance_m'])
        result=dict(start_s=float(raw['t'][0]),end_s=float(raw['t'][-1]),chunks=chunks,
            minimum_sampled_surface_distance_m=minimum['minimum_sampled_surface_distance_m'],
            sampled_witness=minimum['sampled_witness'],
            minimum_conditional_all_pair_clearance_m=min(c['minimum_conditional_all_pair_clearance_m'] for c in chunks),
            intervals_below_100m=sum(c['intervals_below_100m'] for c in chunks))
        out['runs'].append(result)
        print(json.dumps(dict(stage='distance_audit',result={k:v for k,v in result.items() if k!='chunks'},elapsed_s=time.monotonic()-before)),flush=True)
    raw=load_raw(validation['runs'][0]);reference,_=reference_shooting(env,21600.)
    spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
    probes=[]
    def evaluate(x):
        t,angle,rho=x;t*=3600.
        state=np.c_[spline(t),spline(t,1)]
        result=service_margin(state,command(t).as_matrix(),env.at(t),reference.sol(t)[:6],
            rho*np.array([np.cos(angle),np.sin(angle)]))
        probes.append(dict(time_s=float(t),angle_rad=float(angle),solar_radius_fraction=float(rho),**result))
        return result['margin_m']
    for hour in np.linspace(0.,6.,13):
        evaluate([hour,0.,0.])
        for k in range(32):evaluate([hour,.137+np.pi/32+k*np.pi/16,1.])
        for k in range(8):evaluate([hour,.137+np.pi/16+k*np.pi/4,.5])
    grid_count=len(probes);seeds=[];sectors=set()
    for row in sorted(probes,key=lambda r:r['margin_m']):
        sector=int((row['angle_rad']%(2*np.pi))/(np.pi/2))
        if sector not in sectors:
            seeds.append(row);sectors.add(sector)
        if len(seeds)==4:break
    optimizations=[]
    for row in seeds:
        x=np.array([row['time_s']/3600.,row['angle_rad'],row['solar_radius_fraction']])
        simplex=np.tile(x,(4,1))
        simplex[1,0]+=.1 if x[0]<5.9 else -.1
        simplex[2,1]+=np.pi/32
        simplex[3,2]+=.05 if x[2]<.95 else -.05
        fit=minimize(evaluate,x,method='Nelder-Mead',bounds=[(0.,6.),(x[1]-np.pi/8,x[1]+np.pi/8),(0.,1.)],
            options={'maxfev':160,'xatol':1e-5,'fatol':1e-3,'initial_simplex':simplex})
        optimizations.append(dict(success=bool(fit.success),message=fit.message,nfev=int(fit.nfev),
            minimum_margin_m=float(fit.fun),parameters=fit.x.tolist()))
    worst=min(probes,key=lambda r:r['margin_m'])
    out['coverage_adversary']=dict(grid_probes=grid_count,total_probes=len(probes),
        optimizations=optimizations,worst=worst,minimum_fraction=min(r['covered_fraction'] for r in probes),
        passed=bool(worst['margin_m']>0 and all(r['covered_fraction']>=1-1e-9 for r in probes)),
        scope='Rotated 32-direction limb grid, interior/centre sources and four bounded local searches over date and solar-disk coordinates. Positive margin covers the exact circular receiver for each tested source; unsampled global source/time certification remains open.')
    out.update(completed=True,passes=not any(r['intervals_below_100m'] for r in out['runs']) and out['coverage_adversary']['passed'],
        elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    (HERE/'results/packing_audit.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
