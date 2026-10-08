"""Shared bounded-run mechanics for the joint departure/arrival investigation."""
import json
import resource
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial.transform import Rotation,RotationSpline
from shared.provenance import constants_used
from protection.dynamics.fast_parallel import FastParallelPattern
from protection.dynamics.pattern_return import smooth_arc,beam_and_rate_screen
from protection.dynamics.square_distance import closest_squares,distance_guard
from protection.dynamics.pattern_departure import sun_frame
from protection.dynamics.optical import length
from .packing_screen import parents,load_raw
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO

RUN=ROOT/'research/runs/solar_shield_array/joint_cycle'
NEW=['protection/dynamics/fast_parallel.py','protection/dynamics/fast_parallel.cpp',
     'research/studies/solar_shield_array/joint_common.py']


def setup(extra=()):
    RUN.mkdir(parents=True,exist_ok=True)
    names=['closure_validation.json','closure_sizing.json','closure_return.json',
           'packing_validation.json','closure_budget.json','repair_trial_local_refined_fine.json']
    p,sources=parents(names)
    sources.update({path:digest(ROOT/path) for path in [*NEW,*extra]})
    producer=dict(source_hashes=sources,constants=constants_used(sources),
        inputs={name:digest(HERE/'results'/name) for name in names})
    return p,producer


def write(name,out,before,cpu):
    out['resources']=dict(wall_s=time.monotonic()-before,cpu_s=time.process_time()-cpu,
        max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        numerical_threads=1,new_raw_MiB=sum(p.stat().st_size for p in RUN.glob('*'))/1024**2)
    (HERE/'results'/name).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


def command(env,end,arrival=None,tilt=90.,turn_start=7.,turn_hours=4.):
    times=np.unique(np.r_[np.arange(0.,end,120.),end,turn_start*3600.,(turn_start+turn_hours)*3600.])
    def step(x):
        f=np.clip(x,0,1);return 10*f**3-15*f**4+6*f**5
    angles=-np.deg2rad(tilt)*step((times-turn_start*3600.)/(turn_hours*3600.))
    if arrival is not None:
        angles*=1-step((times-arrival+8*3600.)/(8*3600.))
    frames=np.array([sun_frame(env,t) for t in times])@Rotation.from_rotvec(np.c_[angles,np.zeros((len(times),2))]).as_matrix()
    return RotationSpline(times,Rotation.from_matrix(frames))


def spline(raw):
    return CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])


def state_at(path,t):
    return np.concatenate([path(t),path(t,1)],axis=-1)


def propagate(env,initial,start,end,cmd,controls,windows,ratio,label,suns=8,max_step=120.,rtol=2e-10):
    before=time.monotonic();model=FastParallelPattern(env,cmd,ratio,suns)
    maxforce=[0.];progress=[start+4*3600.]
    def rhs(t,flat):
        y=flat.reshape(-1,6)
        a=model.acceleration(t,y)+np.einsum('k,nki->ni',smooth_arc(t,windows),controls)
        maxforce[0]=max(maxforce[0],float(length(a).max()))
        if t>=progress[0]:
            print(label,round(t/3600,3),round(time.monotonic()-before,2),flush=True);progress[0]+=4*3600
        return np.c_[y[:,3:],a].ravel()
    def event(t,flat):return closest_squares(flat.reshape(-1,6)[:,:3],cmd(t).as_matrix())[0]-100.
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=max_step,rtol=rtol,
        atol=np.tile([1e-5]*3+[1e-9]*3,len(initial)),events=event,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    finish=float(sol.t[-1]);dates=np.unique(np.r_[np.arange(start,finish,60.),finish])
    states=sol.sol(dates).T.reshape(-1,len(initial),6)
    path=RUN/(label+'.npz')
    np.savez_compressed(path,t=dates,state=states,controls=controls,windows=windows,mass_ratio=ratio,frames=cmd(dates).as_matrix())
    # Conditional adaptive distance proof: refine any interval whose speed/
    # corner-motion envelope cannot retain the 100 m floor plus 100 m reserve.
    cache={};minimum=[np.inf,None];all_low=[];unresolved=[]
    def at(t):
        if t not in cache:
            y=sol.sol(t).reshape(-1,6);f=cmd(t).as_matrix()
            d,pair,proj=closest_squares(y[:,:3],f)
            speed=float(length(y[:,3:]-y[:,3:].mean(axis=0)).max())
            cache[t]=(d,speed)
            if d<minimum[0]:minimum[:]=[d,dict(time_s=t,pair=pair.tolist(),projection_m=proj.tolist())]
        return cache[t]
    def guard(a,b):
        da,va=at(a);db,vb=at(b);dt=b-a
        low=min(da,db)-(max(va,vb)+10000/np.sqrt(2)*np.deg2rad(.1))*dt-.15*dt**2/2
        if low<200 and dt>0.5:
            mid=(a+b)/2;guard(a,mid);guard(mid,b)
        else:
            all_low.append(low)
            if low<100:unresolved.append([a,b,low])
    for a,b in zip(dates[:-1],dates[1:]):guard(float(a),float(b))
    beam=beam_and_rate_screen(env,dates,states,cmd)
    reconstruct=spline(dict(t=dates,state=states));mid=(dates[:-1]+dates[1:])/2
    error=float(length(reconstruct(mid)-sol.sol(mid).T.reshape(-1,len(initial),6)[:,:,:3]).max())
    row=dict(start_s=start,end_s=finish,intended_end_s=end,completed=abs(finish-end)<1e-6,
        clearance_event=bool(len(sol.t_events[0])),minimum_sampled_surface_distance_m=minimum[0],
        witness=minimum[1],conditional_interval_lower_bound_m=min(all_low),unresolved_intervals=unresolved,
        guard_evaluations=len(cache),guard_acceleration_bound_m_s2=.15,guard_rate_bound_deg_s=.1,
        reconstruction_error_m=error,beam_and_rate=beam,maximum_evaluated_acceleration_m_s2=maxforce[0],
        suns=suns,max_step_s=max_step,rtol=rtol,nfev=sol.nfev,elapsed_s=time.monotonic()-before,
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
    return row,sol
