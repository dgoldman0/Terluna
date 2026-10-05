"""All-361-tile nonlinear executions and independent tighter trajectory replay.

Power overloads do not get silently clipped or marked feasible: an imposed
command beyond a rating is an ideal-force diagnostic and fails acceptance.
"""
import argparse,time,signal,resource
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from protection.dynamics.eclipse_parallel import EclipseParallelPattern
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.cycle_control import independent_path
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.optical import length
from .joint_common import *
from .joint_bridge import EXTRA
from shared import constants as K


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fine',action='store_true');ap.add_argument('--case',default='corridor_108',choices=['corridor_108','early_same']);a=ap.parse_args()
    signal.alarm(850);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    p,producer=setup([*EXTRA,'research/studies/solar_shield_array/joint_execute.py'])
    product='joint_bridge.json' if a.case.startswith('corridor') else 'joint_screen.json'
    parent=json.loads((HERE/'results'/product).read_text());producer['inputs'][product]=digest(HERE/'results'/product)
    case=next(x for x in parent['cases'] if x['spec']['name']==a.case);raw=load_raw(case)
    env=environment(3.);ratio=raw['mass_ratio'];arrival=case['arrival_s'];end=arrival+21600
    cmd=command(env,end,arrival);model=EclipseParallelPattern(env,cmd,ratio,suns=16 if a.fine else 8,grid=32 if a.fine else 16)
    if a.case.startswith('corridor'):
        spline_u=CubicSpline(raw['t'],raw['acceleration'],axis=0)
        def thrust(t):return spline_u(t) if 21600<t<arrival else np.zeros((361,3))
    else:
        def thrust(t):return np.einsum('k,nki->ni',smooth_arc(t,raw['windows']),raw['controls']) if 21600<=t<arrival else np.zeros((361,3))
    initial=load_raw(p['closure_validation.json']['runs'][0])['state'][0]
    out=dict(schema='terluna.research.joint-execution/1',producer=producer,case=a.case,tighter_replay=a.fine,
        arrival_s=arrival,accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        actuator_and_supply_feasible=False,scope='Ideal commanded-acceleration diagnostic; acceptance also requires per-tile actuator and electrical supply, beyond dynamics/geometry.',runs=[])
    label=a.case+('_fine' if a.fine else '_coarse');current=initial;start=0.;solutions=[];lastprint=[-1.]
    for leg,finish in [('service',21600.),('transfer',arrival),('next_service',end)]:
        # Explicit boundaries keep a command interpolation outside a stage from
        # inventing an impulsive or skipped transition.
        def rhs(t,flat):
            y=flat.reshape(-1,6);acc=model.acceleration(t,y)+thrust(t)
            if t-lastprint[0]>=4*3600:
                print(label,leg,round(t/3600,3),round(time.monotonic()-before,2),flush=True);lastprint[0]=t
            return np.c_[y[:,3:],acc].ravel()
        def event(t,flat):return closest_squares(flat.reshape(-1,6)[:,:3],cmd(t).as_matrix())[0]-100.
        event.terminal=True;event.direction=-1
        sol=solve_ivp(rhs,[start,finish],current.ravel(),method='DOP853',max_step=60 if a.fine else 120,
            rtol=2e-12 if a.fine else 2e-10,atol=np.tile([1e-5]*3+[1e-9]*3,361),events=event,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        actual=float(sol.t[-1]);times=np.unique(np.r_[np.arange(start,actual,60),actual]);states=sol.sol(times).T.reshape(-1,361,6)
        path=RUN/(label+'_'+leg+'.npz')
        np.savez_compressed(path,t=times,state=states,mass_ratio=ratio,acceleration=np.array([thrust(t) for t in times]),frames=cmd(times).as_matrix())
        cache={};minimum=[np.inf,None];lower=[];unresolved=[]
        def sample(t):
            if t not in cache:
                y=sol.sol(t).reshape(-1,6);d,pair,proj=closest_squares(y[:,:3],cmd(t).as_matrix())
                speed=float(length(y[:,3:]-y[:,3:].mean(axis=0)).max());cache[t]=(d,speed)
                if d<minimum[0]:minimum[:]=[d,dict(time_s=t,pair=pair.tolist(),projection_m=proj.tolist())]
            return cache[t]
        def guard(lo,hi):
            d0,v0=sample(lo);d1,v1=sample(hi);dt=hi-lo
            bound=min(d0,d1)-(max(v0,v1)+10000/np.sqrt(2)*np.deg2rad(.1))*dt-.15*dt*dt/2
            if bound<300 and dt>.5:mid=(lo+hi)/2;guard(lo,mid);guard(mid,hi)
            else:
                lower.append(bound)
                if bound<100:unresolved.append([lo,hi,bound])
        for t0,t1 in zip(times[:-1],times[1:]):guard(float(t0),float(t1))
        fit=spline(dict(t=times,state=states));mid=(times[1:]+times[:-1])/2
        recon=float(length(fit(mid)-sol.sol(mid).T.reshape(-1,361,6)[:,:,:3]).max())
        row=dict(stage=leg,start_s=start,end_s=actual,intended_end_s=finish,completed=abs(actual-finish)<1e-6,
            clearance_event=bool(len(sol.t_events[0])),minimum_sampled_surface_distance_m=minimum[0],witness=minimum[1],
            conditional_interval_lower_bound_m=min(lower),unresolved_intervals=unresolved,guard_evaluations=len(cache),
            reconstruction_error_m=recon,beam_and_rate=beam_and_rate_screen(env,times,states,cmd),
            raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),nfev=sol.nfev)
        if leg=='transfer':
            target=raw['target'] if 'target' in raw else raw['state'][-1]
            row['terminal_position_error_m']=float(length(states[-1,:,:3]-target[:,:3]).max())
            row['terminal_velocity_error_m_s']=float(length(states[-1,:,3:]-target[:,3:]).max())
        if leg=='next_service':
            # Assigned receiver follows a fresh unshadowed reference beginning
            # at the planned terminal mean. Actual propagated tiles are tested.
            target=raw['target'] if 'target' in raw else raw['state'][-1]
            ref=independent_path(env,target.mean(axis=0)[None,:],arrival,end,cmd,mass_ratio=ratio.mean())
            rows=[]
            for t in np.linspace(arrival,end,49):
                y=sol.sol(t).reshape(-1,6);c=ref.sol(t).reshape(6)
                coverage=local_coverage(y,env.at(t),c,limb_count=16,interior=8,rotation=.137)
                coverage.update(time_s=float(t),window_inside_target=coverage['window_centre_distance_m']+10000<=4*K.MOON_RADIUS);rows.append(coverage)
            row['coverage']=dict(samples=rows,passed=all(x['minimum_source_coverage']>=1-1e-9 and x['window_inside_target'] for x in rows),
                minimum_coverage=min(x['minimum_source_coverage'] for x in rows),sources_per_date=25)
        out['runs'].append(row);write('joint_'+label+'.json',out,before,cpu)
        print('finished',label,leg,actual/3600,minimum[0],flush=True)
        current=states[-1];start=actual;solutions.append(sol)
        if not row['completed']:break
    if a.fine:
        coarse=json.loads((HERE/'results'/('joint_'+a.case+'_coarse.json')).read_text());comparison=[]
        producer['inputs']['joint_'+a.case+'_coarse.json']=digest(HERE/'results'/('joint_'+a.case+'_coarse.json'))
        for row,old,sol in zip(out['runs'],coarse['runs'],solutions):
            oldraw=load_raw(old);oldpath=spline(oldraw);stop=min(row['end_s'],old['end_s'])
            ts=np.unique(np.r_[np.arange(row['start_s'],stop,30.),stop]);y=sol.sol(ts).T.reshape(-1,361,6)
            err=float(length(y[:,:,:3]-oldpath(ts)).max())+row['reconstruction_error_m']+old['reconstruction_error_m']
            comparison.append(dict(stage=row['stage'],maximum_position_with_reconstruction_m=err,
                maximum_velocity_m_s=float(length(y[:,:,3:]-oldpath(ts,1)).max()),
                conditional_clearance_after_replay_allowance_m=row['conditional_interval_lower_bound_m']-2*err,
                time_difference_s=abs(row['end_s']-old['end_s'])))
        out['comparison']=comparison
    out['completed']=True;write('joint_'+label+'.json',out,before,cpu)
if __name__=='__main__':main()
