"""Coverage transitions, next-service margin and repeat-state diagnostics."""
import time,signal,resource
import numpy as np
from scipy.optimize import minimize_scalar
from protection.dynamics.service_margin import service_margin
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.cycle_control import independent_path
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.pattern_departure import sun_frame
from protection.dynamics.optical import length
from shared import constants as K
from .joint_common import *


def main():
    signal.alarm(110);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time();p,producer=setup(['research/studies/solar_shield_array/joint_audit.py'])
    fine=json.loads((HERE/'results/joint_corridor_108_fine.json').read_text());case=json.loads((HERE/'results/joint_bridge.json').read_text())['cases'][0]
    for name in ['joint_corridor_108_fine.json','joint_bridge.json']:producer['inputs'][name]=digest(HERE/'results'/name)
    env=environment(3);arrival=fine['arrival_s'];end=fine['runs'][-1]['end_s'];cmd=command(env,end,arrival)
    target=load_raw(case)['target'];ratio=load_raw(case)['mass_ratio']
    reference=independent_path(env,target.mean(axis=0)[None,:],arrival,end,cmd,mass_ratio=ratio.mean())
    next_raw=load_raw(fine['runs'][-1]);path=spline(next_raw);probes=[]
    for t in np.linspace(arrival,end,13):
        for angle,rho in [(0,0)]+[(.219+k*2*np.pi/32,1.) for k in range(32)]+[(.413+k*2*np.pi/8,.5) for k in range(8)]:
            row=service_margin(state_at(path,t),cmd(t).as_matrix(),env.at(t),reference.sol(t).reshape(6),rho*np.array([np.cos(angle),np.sin(angle)]))
            probes.append(dict(time_s=float(t),angle=float(angle),rho=rho,**row))
    worst=min(probes,key=lambda x:x['margin_m'])
    def objective(t):
        r=service_margin(state_at(path,t),cmd(t).as_matrix(),env.at(t),reference.sol(t).reshape(6),worst['rho']*np.array([np.cos(worst['angle']),np.sin(worst['angle'])]))
        probes.append(dict(time_s=float(t),angle=worst['angle'],rho=worst['rho'],**r));return r['margin_m']
    minimize_scalar(objective,bounds=(max(arrival,worst['time_s']-900),min(end,worst['time_s']+900)),method='bounded',options={'maxiter':30,'xatol':.1})
    worst=min(probes,key=lambda x:x['margin_m'])
    trans=load_raw(fine['runs'][1]);sp=spline(trans);coverage=[]
    dates=np.unique(np.r_[np.arange(21600,arrival,1800.),np.arange(6,12.1,.5)*3600,arrival])
    for t in dates:
        y=state_at(sp,t);centre=y.mean(axis=0);sun=env.at(t)['positions']['sun'];sun=sun/length(sun)
        upstream=bool(np.min(y[:,:3]@sun)>0)
        if upstream:
            cv=local_coverage(y,env.at(t),centre,limb_count=8,interior=4,rotation=.191)
        else:cv=dict(minimum_source_coverage=0.,window_centre_distance_m=float(length(centre[:3]-(centre[:3]@sun)*sun)))
        coverage.append(dict(time_s=float(t),upstream=upstream,minimum_coverage=cv['minimum_source_coverage'],
            window_inside_target=cv['window_centre_distance_m']+10000<=4*K.MOON_RADIUS,assigned_service=False))
    first=load_raw(fine['runs'][0])['state'][-1];last=next_raw['state'][-1]
    f0=sun_frame(env,21600);f1=sun_frame(env,end)
    dq=(first[:,:3]-first[:,:3].mean(axis=0))@f0-(last[:,:3]-last[:,:3].mean(axis=0))@f1
    dv=(first[:,3:]-first[:,3:].mean(axis=0))@f0-(last[:,3:]-last[:,3:].mean(axis=0))@f1
    geometry=[]
    # Refine each coarse minimum with a scalar-time search, retaining all pairs
    # at each evaluation instead of assuming the coarse witness remains closest.
    for row in fine['runs']:
        r=load_raw(row);s=spline(r);w=row['witness']['time_s']
        fun=lambda t:closest_squares(s(t),cmd(t).as_matrix())[0]
        lo=max(row['start_s'],w-120);hi=min(row['end_s'],w+120)
        fit=minimize_scalar(fun,bounds=(lo,hi),method='bounded',options={'xatol':.02,'maxiter':50})
        d,pair,_=closest_squares(s(fit.x),cmd(fit.x).as_matrix())
        geometry.append(dict(stage=row['stage'],time_s=float(fit.x),surface_distance_m=float(d),pair=pair.tolist(),evaluations=fit.nfev))
    out=dict(schema='terluna.research.joint-audit/1',producer=producer,completed=True,accepted_cycle=False,accepted_fleet=False,
        geometric_return_plus_service_replayed=bool(all(r['completed'] and r['beam_and_rate']['minimum_sampled_beam_margin_deg']>0 and r['beam_and_rate']['maximum_sampled_rate_deg_s']<=.1 for r in fine['runs']) and all(c['maximum_position_with_reconstruction_m']<50 and c['conditional_clearance_after_replay_allowance_m']>100 for c in fine['comparison']) and fine['runs'][-1]['coverage']['passed']),
        continuous_optical_certificate=False,handovers_executed=0,required_replacement_service_interval_s=[21600.,arrival],
        next_service_adversary=dict(probes=len(probes),worst=worst,passed=worst['margin_m']>0),
        off_service_coverage_samples=coverage,refined_geometry=geometry,
        terminal_attitude_error_rad=float(length((cmd(arrival).inv()*Rotation.from_matrix(sun_frame(env,arrival))).as_rotvec())),
        next_service_end_vs_first_service_end=dict(maximum_relative_position_difference_m=float(length(dq).max()),
            maximum_relative_velocity_difference_m_s=float(length(dv).max()),
            interpretation='Frame-normalized differences, not an executed reset. A further cycle and epoch-dependent acquisition controls remain untested.'))
    write('joint_audit.json',out,before,cpu);print(out['resources'],worst,flush=True)
if __name__=='__main__':main()
