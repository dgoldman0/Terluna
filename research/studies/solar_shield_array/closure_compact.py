"""Tighter replay using exact compact shadow unions and the same physical guards.

The earlier runner remains byte-pinned to its completed products. This explicit
replay variant records its different force backend in the producer identity.
"""
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline

from protection.dynamics.compact_shadow import CompactMassivePattern
from protection.dynamics.pattern_return import smooth_arc,beam_and_rate_screen,crossing_witnesses
from protection.dynamics.square_distance import closest_squares,distance_guard
from protection.dynamics.collection import require_uneclipsed_tiles
from protection.dynamics.cycling import ray_geometry
from protection.dynamics.optical import length
from .closure_screen import RUN
from .fleet_run import digest
from .cycling_search import ROOT


def propagate(env,initial,start,end,command,controls,windows,mass_ratio,label,
              suns=8,max_step=120.,rtol=2e-10):
    before=time.monotonic();model=CompactMassivePattern(env,command,mass_ratio,suns=suns)
    maximum=[0.];progress=[start+3600.]
    def rhs(t,flat):
        state=flat.reshape(-1,6)
        require_uneclipsed_tiles(ray_geometry(state[:,:3],env.at(t)))
        a=model.acceleration(t,state)+np.einsum('k,nki->ni',smooth_arc(t,windows),controls)
        maximum[0]=max(maximum[0],float(length(a).max()))
        if t>=progress[0]:
            print(f'{label}: {t/3600:.3f} h; {time.monotonic()-before:.2f} s',flush=True)
            progress[0]+=3600.
        return np.c_[state[:,3:],a].ravel()
    def event(t,flat):
        return closest_squares(flat.reshape(-1,6)[:,:3],command(t).as_matrix())[0]-100.
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=max_step,rtol=rtol,
        atol=np.tile([1e-5]*3+[1e-9]*3,len(initial)),events=event,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    finish=float(sol.t[-1]);times=np.linspace(start,finish,int(np.ceil((finish-start)/30.))+1)
    states=sol.sol(times).T.reshape(-1,len(initial),6)
    spline=CubicHermiteSpline(times,states[:,:,:3],states[:,:,3:])
    guard=[];beam=[];error=0.
    for t0 in np.arange(start,finish,3600.):
        t1=min(finish,t0+3600.);t=np.linspace(t0,t1,int(np.ceil((t1-t0)/2.))+1)
        y=sol.sol(t).T.reshape(-1,len(initial),6)
        error=max(error,float(length(y[:,:,:3]-spline(t)).max()))
        guard.append(distance_guard(t,y,command))
        beam.append(beam_and_rate_screen(env,t[::5],y[::5],command))
    crossing=crossing_witnesses(times,states,command,lambda t,ids:sol.sol(t).reshape(-1,6)[ids],retain=16)
    sampled=min(guard,key=lambda r:r['minimum_sampled_surface_distance_m'])
    peak=float((length(controls)*2/np.diff(windows,axis=1).ravel()).max())
    def primitive(t):
        x=np.clip((t-windows[:,0])/np.diff(windows,axis=1).ravel(),0.,1.)
        return x-np.sin(2*np.pi*x)/(2*np.pi)
    applied=length(controls)@(primitive(finish)-primitive(start))
    path=RUN/f'{label}.npz'
    np.savez_compressed(path,t=times,state=states,controls=controls,windows=windows,
        mass_ratio=mass_ratio,frames=command(times).as_matrix())
    complete=bool(abs(finish-end)<1e-6)
    result=dict(start_s=start,end_s=finish,intended_end_s=end,completed=complete,
        clearance_event=bool(len(sol.t_events[0])),suns=suns,max_step_s=max_step,rtol=rtol,
        minimum_sampled_surface_distance_m=sampled['minimum_sampled_surface_distance_m'],
        sampled_witness=sampled['sampled_witness'],
        minimum_conditional_surface_distance_m=min(g['minimum_conditional_all_pair_clearance_m'] for g in guard),
        flagged_intervals=sum(g['intervals_below_100m'] for g in guard),
        acceleration_guard_m_s2=.15,rate_guard_deg_s=.1,
        beam_margin_deg=min(b['minimum_sampled_beam_margin_deg'] for b in beam),
        rate_deg_s=max(b['maximum_sampled_rate_deg_s'] for b in beam),
        angular_acceleration_deg_s2=max(b['maximum_sampled_angular_acceleration_deg_s2'] for b in beam),
        maximum_evaluated_acceleration_m_s2=maximum[0],maximum_command_acceleration_m_s2=peak,
        scheduled_translation_m_s=float(length(controls).sum(axis=1).mean()),
        mean_applied_translation_m_s=float(applied.mean()),maximum_applied_translation_m_s=float(applied.max()),
        scheduled_burn_count=int(np.count_nonzero(length(controls)>1e-8)),
        started_burn_count=int(np.count_nonzero((length(controls)>1e-8)&(windows[None,:,0]<finish))),
        reconstruction_error_m=error,crossings=crossing,
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
        elapsed_s=time.monotonic()-before,nfev=sol.nfev)
    result['geometry_pass']=bool(complete and not result['flagged_intervals'] and not crossing['crossing_near_encounters']
        and result['beam_margin_deg']>0 and result['rate_deg_s']<=.1 and peak<=.001 and maximum[0]<=.15)
    return result,times,states
