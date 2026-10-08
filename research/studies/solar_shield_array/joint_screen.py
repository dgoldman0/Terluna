"""Bounded joint departure/arrival screen; retain every failed branch."""
import time,signal,resource
import numpy as np
from scipy.optimize import linear_sum_assignment
from protection.dynamics.cycle_control import independent_path
from protection.dynamics.packing_control import response_maps,matrices,packing_offsets
from protection.dynamics.joint_transfer import fit_controls,propose,violations
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from shared import constants as K
from .joint_common import *


def main():
    signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    p,producer=setup(['protection/dynamics/joint_transfer.py','research/studies/solar_shield_array/joint_screen.py'])
    env=environment(3.);ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    initial=load_raw(p['closure_validation.json']['runs'][0])['state'][-1];start=21600.
    original=load_raw(p['closure_validation.json']['runs'][0])['state'][0]
    mass=5e6*ratio;specific=30000/(2*.7*np.cos(np.pi/4));installed=(ratio-1.2)*5e6*300
    arrival0=p['closure_return.json']['runs'][0]['intended_end_s']
    out=dict(schema='terluna.research.joint-screen/1',producer=producer,cases=[],accepted_return=False,accepted_cycle=False,
        scope='Independent finite-source dynamics and gravity response proposals; hardware-rated endpoint fits with local contact cuts. Source shadows, nonlinear control response and next-service replay remain acceptance gates.')
    specs=[dict(name='early_same',depth=12000,assignment='identity',stagger=0,tilt=90,shift=0,slack=False),
           dict(name='deep_assignment',depth=24000,assignment='hungarian',stagger=0,tilt=90,shift=-1800,slack=False),
           dict(name='staggered_flexible',depth=24000,assignment='hungarian',stagger=1800,tilt=90,shift=-1800,slack=True),
           dict(name='early_steering_75',depth=24000,assignment='hungarian',stagger=1800,tilt=75,shift=-3600,slack=True),
           dict(name='early_steering_105',depth=24000,assignment='hungarian',stagger=1800,tilt=105,shift=-3600,slack=True),
           dict(name='deep_identity',depth=36000,assignment='identity',stagger=0,tilt=90,shift=-3600,slack=True)]
    for spec in specs:
        end=arrival0+spec['shift'];cmd=command(env,end+21600,end,spec['tilt'])
        # Shared windows allow independent tile vectors. Arrival slots switch
        # off arcs ending after a member's earlier hand-in time.
        windows=np.array([[6,10],[10,14],[18,24],[28,34],[end/3600-10,end/3600-6],[end/3600-4,end/3600]])*3600
        # Avoid overlap for earlier arrivals by shortening the middle arc.
        windows[3,1]=min(windows[3,1],windows[4,0])
        free=independent_path(env,initial,start,end,cmd,mass_ratio=ratio)
        base=lambda t:free.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
        maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),start,end,windows,np.eye(3))
        centre=base([end])[0].mean(axis=0)
        service=independent_path(env,centre[None,:],end,end+21600,cmd,mass_ratio=ratio.mean())
        sm=response_maps(env,lambda t:service.sol(t),end,end+21600,np.empty((0,2)),np.eye(3))
        phi=matrices(sm,[end+21600],0)[0];frame=solar_frame(env.at(end))
        gradient=.95*np.linalg.solve(phi[:3,3:6],solar_frame(env.at(end+21600))-phi[:3,:3]@frame)
        offsets=packing_offsets(spec['depth']);target=centre+np.c_[offsets@frame.T,offsets@gradient.T]
        assignment=np.arange(361)
        if spec['assignment']=='hungarian':
            # Assignment cost includes both position and velocity, on a 6 h
            # transfer scale. This is a proposal heuristic, not an optimum.
            predicted=base([end])[0]
            cost=length((predicted[:,None,:3]-target[None,:,:3])/21600)+length(predicted[:,None,3:]-target[None,:,3:])
            _,assignment=linear_sum_assignment(cost);target=target[assignment]
        labels=((np.arange(361)//19)%2*2+(np.arange(361)%19)%2)[assignment]
        endpoint_times=end-spec['stagger']*labels/3
        # Backcast each staggered terminal onto the same next-service orbit;
        # no target teleport and no coast impulse after this endpoint is free.
        for i,t in enumerate(endpoint_times):
            if t==end:continue
            dt=t-end;rel=target[i]-centre
            # Local Taylor backcast; its truncation is an explicit proposal
            # defect to be measured in the nonlinear re-entry/service replay.
            q=target[i,:3];a=env.gravity(end,q)
            target[i,:3]+=target[i,3:]*dt+.5*a*dt*dt
            target[i,3:]+=a*dt
        caps=np.zeros((361,len(windows)))
        for k,(a,b) in enumerate(windows):
            grid=np.linspace(a,b,65);att=float(rigid_square_load(cmd,grid)['force_per_mass'].max())
            caps[:,k]=np.maximum(0,np.minimum(.001,.95*installed/(mass*specific)-att-3e-5))
            caps[endpoint_times<b-1e-6,k]=0
        slack=np.array([1500]*3+[.03]*3) if spec['slack'] else None
        u,fit=fit_controls(base,maps,target,windows,caps,endpoint_times,slack=slack)
        row=dict(spec=spec,arrival_s=end,endpoint_slot_range_s=[float(endpoint_times.min()),float(endpoint_times.max())],
            windows_s=windows.tolist(),initial_fit=fit,steps=[],assignment_changed=int(np.count_nonzero(assignment!=np.arange(361))))
        if u is not None:
            dates=np.linspace(start,end,int(np.ceil((end-start)/180))+1);cuts={}
            for attempt in range(7):
                state=propose(base,maps,u,dates)
                additions,scan=violations(dates,state,cmd,history=initial)
                row['steps'].append(dict(attempt=attempt,scan=scan))
                if not additions:break
                if attempt==6:break
                for t,i,j,axis,needed in additions:
                    key=(round(t,1),i,j,tuple(np.round(axis,5)));cuts[key]=(t,i,j,axis,needed)
                v,f=fit_controls(base,maps,target,windows,caps,endpoint_times,cuts=list(cuts.values()),slack=slack,wall_s=8.)
                row['steps'][-1]['fit']=f
                if v is None:break
                u=v
            state=propose(base,maps,u,dates)
            row['final_scan']=violations(dates,state,cmd,history=initial)[1]
            row['mean_translation_m_s']=float(length(u).sum(axis=1).mean())
            row['beam_and_rate']=beam_and_rate_screen(env,dates[::3],state[::3],cmd)
            row['entry_coverage']=local_coverage(state[-1],env.at(end),centre,limb_count=16,interior=8)
            service_dates=np.linspace(end,end+21600,13)
            row['reference_next_service_inside_target']=all(length(service.sol(t)[:3]-(service.sol(t)[:3]@solar_frame(env.at(t))[:,0])*solar_frame(env.at(t))[:,0])+10000<=4*K.MOON_RADIUS for t in service_dates)
            path=RUN/(spec['name']+'_proposal.npz')
            np.savez_compressed(path,t=dates,state=state,controls=u,windows=windows,mass_ratio=ratio,
                initial=initial,target=target,assignment=assignment,endpoint_times=endpoint_times,
                caps=caps,free_state=base(dates),maps=maps.sol(dates).T)
            row.update(raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
        out['cases'].append(row);write('joint_screen.json',out,before,cpu)
        print(spec['name'],row.get('final_scan'),fit,round(time.monotonic()-before,2),flush=True)
    out['completed']=True;write('joint_screen.json',out,before,cpu)
if __name__=='__main__':main()
