"""Coupled all-member replay with terminal clearance and actuator events."""
import argparse,signal,resource,time
import numpy as np
from scipy.integrate import solve_ivp
from protection.dynamics.service_search import free_command
from protection.dynamics.eclipse_parallel import EclipseParallelPattern,load
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.pattern_return import smooth_arc,beam_and_rate_screen
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .joint_common import spline,state_at
from .natural_common import *


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fine',action='store_true');ap.add_argument('--case',default='free_slow');a=ap.parse_args()
    signal.alarm(500);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    extra=['protection/dynamics/eclipse_parallel.py','protection/dynamics/eclipse_parallel.cpp','research/studies/solar_shield_array/natural_execute.py']
    arrays,producer=setup(extra);ratio=arrays['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300
    scout=json.loads((HERE/'results/natural_fast_scout.json').read_text());producer['inputs']['natural_fast_scout.json']=digest(HERE/'results/natural_fast_scout.json')
    source=next(x for x in scout['cases'] if x['spec']['name']=='edge_a');arrival=source['opportunities'][0]['arrival_s']
    tilt=90.;controls=np.zeros((361,0,3));windows=np.empty((0,2))
    receiver_base=None
    if a.case=='face_refined':
        parent=json.loads((HERE/'results/natural_refine.json').read_text());producer['inputs']['natural_refine.json']=digest(HERE/'results/natural_refine.json')
        raw=np.load(ROOT/parent['raw_path']);assert digest(ROOT/parent['raw_path'])==parent['raw_sha256']
        arrival=parent['arrival_s'];tilt=0.;controls=raw['controls'];windows=raw['windows']
        receiver_base=spline(dict(t=raw['t'],state=raw['base_state']))
    elif a.case!='free_slow':
        parent=json.loads((HERE/'results/natural_dynamic_rays.json').read_text());producer['inputs']['natural_dynamic_rays.json']=digest(HERE/'results/natural_dynamic_rays.json')
        row=next(x for x in parent['cases'] if x['name']==a.case);raw=np.load(ROOT/row['raw_path']);assert digest(ROOT/row['raw_path'])==row['raw_sha256']
        arrival=row['arrival_s'];tilt=row['tilt'];controls=raw['controls'];windows=raw['windows']
        receiver_base=spline(dict(t=raw['t'],state=raw['base_state']))
    original=json.loads((HERE/'results/joint_regression.json').read_text())['runs'][0]
    producer['inputs']['joint_regression.json']=digest(HERE/'results/joint_regression.json')
    service_raw=np.load(ROOT/original['raw_path']);assert digest(ROOT/original['raw_path'])==original['raw_sha256']
    original_receiver=spline(service_raw)
    env=environment(arrival/86400+1.);end=arrival+21600;cmd=free_command(env,end,tilt=tilt,arrival=arrival,turn_start=6.,turn_hours=6.)
    suns=16 if a.fine else 8;grid=32 if a.fine else 16
    model=EclipseParallelPattern(env,cmd,ratio,suns=suns,grid=grid)
    def thrust(t):return np.einsum('k,nki->ni',smooth_arc(t,windows),controls) if len(windows) else np.zeros((361,3))
    def power(t,y):
        value=load(env,t,y,cmd(t).as_matrix(),suns=suns,grid=grid,torque=True)
        nominal=rigid_square_load(cmd,np.array([t]))['torque_per_mass'][0]
        torque=nominal-value['gravity_torque_per_mass']-value['radiation_torque_N_m']/mass[:,None]
        return (length(thrust(t))+2*abs(torque).sum(axis=1)/10000)*mass*30000/(1.4*np.cos(np.pi/4))
    label=a.case+('_fine' if a.fine else '_coarse')
    out=dict(schema='terluna.research.natural-execution/1',producer=producer,completed=False,case=a.case,tighter_replay=a.fine,arrival_s=arrival,
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,actual_generation_W=None,actual_bus_supply_W=None,runs=[],
        settings=dict(suns=suns,area_grid=grid,rtol=2e-12 if a.fine else 2e-10,max_step_s=60 if a.fine else 120),
        scope='All 361 loaded members from the original epoch; no reset or imposed intermediate formation. Stop at physical 100 m floor or installed actuator rating. Electrical generation is uninstalled and earns no supply credit.')
    current=arrays['initial_epoch_state'];solutions=[]
    for stage,start,stop in [('service',0.,21600.),('departure_return',21600.,arrival),('next_service',arrival,end)]:
        def rhs(t,flat):
            y=flat.reshape(361,6);return np.c_[y[:,3:],model.acceleration(t,y)+thrust(t)].ravel()
        def clearance(t,flat):return closest_squares(flat.reshape(361,6)[:,:3],cmd(t).as_matrix())[0]-100.
        def rating(t,flat):return float(np.min(1-power(t,flat.reshape(361,6))/capacity))
        clearance.terminal=True;clearance.direction=-1;rating.terminal=True;rating.direction=-1
        sol=solve_ivp(rhs,[start,stop],current.ravel(),method='DOP853',max_step=60 if a.fine else 120,rtol=2e-12 if a.fine else 2e-10,
            atol=np.tile([1e-5]*3+[1e-9]*3,361),events=[clearance,rating],dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        finish=float(sol.t[-1]);ts=np.unique(np.r_[np.arange(start,finish,60.),finish]);ys=sol.sol(ts).T.reshape(-1,361,6)
        cache={};lower=[];minimum=[np.inf,None];unresolved=[]
        def sample(t):
            if t not in cache:
                y=sol.sol(t).reshape(361,6);d,pair,projection=closest_squares(y[:,:3],cmd(t).as_matrix())
                cache[t]=(d,float(length(y[:,3:]-y[:,3:].mean(axis=0)).max()))
                if d<minimum[0]:minimum[:]=[d,dict(time_s=t,pair=pair.tolist(),projection_m=projection.tolist())]
            return cache[t]
        def guard(lo,hi):
            d0,v0=sample(lo);d1,v1=sample(hi);dt=hi-lo
            bound=min(d0,d1)-(max(v0,v1)+10000/np.sqrt(2)*np.deg2rad(.1))*dt-.15*dt**2/2
            if bound<300 and dt>.5:mid=(lo+hi)/2;guard(lo,mid);guard(mid,hi)
            else:
                lower.append(bound)
                if bound<100:unresolved.append([lo,hi,bound])
        for lo,hi in zip(ts[:-1],ts[1:]):guard(float(lo),float(hi))
        fit=spline(dict(t=ts,state=ys));mid=(ts[1:]+ts[:-1])/2
        recon=float(length(fit(mid)-sol.sol(mid).T.reshape(-1,361,6)[:,:,:3]).max())
        row=dict(stage=stage,start_s=start,end_s=finish,intended_end_s=stop,completed=abs(finish-stop)<1e-6,
            clearance_event=bool(len(sol.t_events[0])),actuator_event=bool(len(sol.t_events[1])),minimum_sampled_surface_distance_m=minimum[0],witness=minimum[1],
            conditional_interval_lower_bound_m=min(lower),unresolved_intervals=unresolved,reconstruction_error_m=recon,
            beam_and_rate=beam_and_rate_screen(env,ts,ys,cmd),**save(label+'_'+stage,t=ts,state=ys,mass_ratio=ratio,controls=controls,windows=windows))
        if stage in ['service','next_service'] and row['completed']:
            # Initial reference is the validated passage's actual mean; a new
            # return must separately define and verify its actual receiver path.
            samples=[]
            for t in np.linspace(start,finish,25):
                y=sol.sol(t).reshape(361,6)
                reference=original_receiver if stage=='service' else receiver_base
                centre=state_at(reference,t).mean(axis=0) if reference is not None else y.mean(axis=0)
                samples.append(local_coverage(y,env.at(t),centre,limb_count=16,interior=8,rotation=.137))
            row['sampled_coverage']=dict(minimum=min(x['minimum_source_coverage'] for x in samples),passed=all(x['minimum_source_coverage']>=1-1e-9 for x in samples))
        out['runs'].append(row);solutions.append(sol);write('natural_'+label+'.json',out,before,cpu)
        print(label,stage,finish/3600,minimum[0],row['actuator_event'],time.monotonic()-before,flush=True)
        current=ys[-1]
        if not row['completed']:break
    if a.fine:
        name='natural_'+a.case+'_coarse.json';coarse=json.loads((HERE/'results'/name).read_text());producer['inputs'][name]=digest(HERE/'results'/name)
        out['comparison']=[]
        for row,old,sol in zip(out['runs'],coarse['runs'],solutions):
            raw=np.load(ROOT/old['raw_path']);assert digest(ROOT/old['raw_path'])==old['raw_sha256'];path=spline(raw)
            finish=min(row['end_s'],old['end_s']);times=np.unique(np.r_[np.arange(row['start_s'],finish,30.),finish]);y=sol.sol(times).T.reshape(-1,361,6)
            error=float(length(y[:,:,:3]-path(times)).max())+row['reconstruction_error_m']+old['reconstruction_error_m']
            out['comparison'].append(dict(stage=row['stage'],maximum_position_with_reconstruction_m=error,maximum_velocity_difference_m_s=float(length(y[:,:,3:]-path(times,1)).max()),event_time_difference_s=abs(row['end_s']-old['end_s']),conditional_clearance_after_numerical_allowance_m=row['conditional_interval_lower_bound_m']-2*error))
    out['completed']=True;write('natural_'+label+'.json',out,before,cpu)
if __name__=='__main__':main()
