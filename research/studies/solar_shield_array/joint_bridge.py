"""Separated return corridor: a charged inverse-control comparator.

Preserve the actual first service; capture its relative velocity gradually,
keep a spaced three-dimensional formation, then acquire the next service
position AND velocity. Polynomial geometry only proposes open-loop controls.
A separate all-tile IVP must test those controls without resetting states.
"""
import time,signal,resource
import numpy as np
from scipy.interpolate import BPoly,CubicSpline
from protection.dynamics.cycle_control import independent_path,independent_force
from protection.dynamics.packing_control import response_maps,matrices,packing_offsets
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.eclipse_parallel import load,EclipseParallelPattern
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.natural_pattern import local_coverage,finite_gravity
from protection.dynamics.optical import length
from .joint_common import *
from shared import constants as K

EXTRA=['protection/dynamics/eclipse_parallel.py','protection/dynamics/eclipse_parallel.cpp',
       'research/studies/solar_shield_array/joint_bridge.py']


def build(env,initial,ratio,arrival,spec):
    cmd=command(env,arrival+21600,arrival)
    ref=independent_path(env,initial.mean(axis=0)[None,:],21600.,arrival+21600.,cmd,mass_ratio=ratio.mean())
    centre=ref.sol(arrival).reshape(6)
    sm=response_maps(env,lambda t:ref.sol(t),arrival,arrival+21600.,np.empty((0,2)),np.eye(3))
    phi=matrices(sm,[arrival+21600.],0)[0];f=solar_frame(env.at(arrival))
    gradient=.95*np.linalg.solve(phi[:3,3:6],solar_frame(env.at(arrival+21600.))-phi[:3,:3]@f)
    offsets=packing_offsets(spec['depth'])
    terminal=centre+np.c_[offsets@f.T,offsets@gradient.T]
    model=EclipseParallelPattern(env,cmd,ratio)
    rel0=initial-initial.mean(axis=0);relt=terminal-centre
    def refacc(t):return independent_force(env,t,ref.sol(t).reshape(1,6),cmd,ratio.mean())[0]
    a0=model.acceleration(21600.,initial)-refacc(21600.)
    at=model.acceleration(arrival,terminal)-refacc(arrival)
    hold=spec['scale']*rel0[:,:3]
    dates=[21600.,spec['capture_hour']*3600.,arrival-spec['approach_hours']*3600.,arrival]
    values=[[rel0[:,:3],rel0[:,3:],a0],
            [hold,np.zeros_like(hold),np.zeros_like(hold)],
            [hold,np.zeros_like(hold),np.zeros_like(hold)],
            [relt[:,:3],relt[:,3:],at]]
    path=BPoly.from_derivatives(dates,values)
    def desired(t):
        ts=np.atleast_1d(t);c=ref.sol(ts).T
        return np.concatenate([c[:,None,:3]+path(ts),c[:,None,3:]+path(ts,1)],axis=-1)
    return cmd,ref,path,desired,terminal,refacc


def main():
    signal.alarm(200);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time();p,producer=setup(EXTRA)
    env=environment(3.);ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    initial=load_raw(p['closure_validation.json']['runs'][0])['state'][-1]
    arrival=p['closure_return.json']['runs'][0]['intended_end_s']-1800
    out=dict(schema='terluna.research.joint-corridor/1',producer=producer,cases=[],accepted_cycle=False,
        evidence='Inverse-control trajectory proposals; loads include charged departure capture and arrival preparation; no dynamics acceptance from kinematic paths.')
    specs=[dict(name='corridor_108',scale=1.08,capture_hour=14,approach_hours=8,depth=12000),
           dict(name='corridor_120',scale=1.2,capture_hour=12,approach_hours=10,depth=12000),
           dict(name='corridor_deep',scale=1.08,capture_hour=14,approach_hours=10,depth=18000)]
    for spec in specs[:2]:
        cmd,ref,poly,desired,target,refacc=build(env,initial,ratio,arrival,spec)
        times=np.unique(np.r_[np.arange(21600,arrival,120.),arrival])
        y=desired(times);control=[];powers=[];min_d=np.inf;witness=None;energy=[]
        nominal=rigid_square_load(cmd,times)['torque_per_mass'];specific=30000/(2*.7*np.cos(np.pi/4))
        for k,t in enumerate(times):
            f=cmd(t).as_matrix();l=load(env,t,y[k],f,suns=8,torque=True)
            g=finite_gravity(env,t,y[k,:,:3],f[:,0],f[:,1])
            a=refacc(t)+poly(t,2)-g-l['reflected_force_N']/(5e6*ratio[:,None]);control.append(a)
            disturbance=l['gravity_torque_per_mass']+l['radiation_torque_N_m']/(5e6*ratio[:,None])
            attitude=2*abs(nominal[k]-disturbance).sum(axis=1)/10000
            powers.append((length(a)+attitude)*5e6*ratio*specific)
            d,pair,projection=closest_squares(y[k,:,:3],f)
            if d<min_d:min_d=d;witness=dict(time_s=float(t),pair=pair.tolist(),projection_m=projection.tolist())
        control=np.array(control);powers=np.array(powers);capacity=(ratio-1.2)*5e6*300
        path=RUN/(spec['name']+'_proposal.npz')
        np.savez_compressed(path,t=times,state=y,acceleration=control,mass_ratio=ratio,target=target,
            power_W=powers,frames=cmd(times).as_matrix())
        inside=all(length(ref.sol(t)[:3]-(ref.sol(t)[:3]@solar_frame(env.at(t))[:,0])*solar_frame(env.at(t))[:,0])+10000<=4*K.MOON_RADIUS for t in np.linspace(arrival,arrival+21600,25))
        row=dict(spec=spec,arrival_s=arrival,minimum_sampled_distance_m=float(min_d),witness=witness,
            proposed_translation_m_s=float(np.trapezoid(length(control),times,axis=0).mean()),
            proposed_energy_J=float(np.trapezoid(powers.sum(axis=1),times)),
            peak_acceleration_m_s2=float(length(control).max()),
            overloaded_members=int(np.count_nonzero(powers.max(axis=0)>capacity)),
            minimum_installed_to_peak=float(np.min(capacity/powers.max(axis=0))),
            beam_and_rate=beam_and_rate_screen(env,times,y,cmd),
            target_coverage=local_coverage(y[-1],env.at(arrival),ref.sol(arrival),limb_count=16,interior=8),
            next_reference_inside_target=bool(inside),raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
        out['cases'].append(row);write('joint_bridge.json',out,before,cpu)
        print(spec['name'],min_d,row['proposed_translation_m_s'],row['overloaded_members'],time.monotonic()-before,flush=True)
    out['completed']=True;write('joint_bridge.json',out,before,cpu)
if __name__=='__main__':main()
