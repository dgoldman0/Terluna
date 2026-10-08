"""Recompute loaded baseline and original encounter with exact fast shadows."""
import time,signal,resource
import numpy as np
from protection.dynamics.fast_parallel import load
from protection.dynamics.tile_torque import parallel_load
from protection.dynamics.optical import length
from protection.dynamics.cycle_control import return_command
from .packing_run import full_command,service_coverage
from .patterns_run import reference_shooting
from .joint_common import *


def main():
    signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    p,producer=setup(['research/studies/solar_shield_array/joint_regression.py'])
    env=environment(3.);ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    cmd=full_command(env,p['packing_validation.json']['schedule'])
    out=dict(schema='terluna.research.joint-regression/1',producer=producer,runs=[],backend_checks=[],
        accepted_cycle=False,accepted_fleet=False)
    initial=load_raw(p['closure_validation.json']['runs'][0])['state'][0]
    for k,entry in enumerate(p['closure_validation.json']['runs']):
        raw=load_raw(entry)
        if k==2:cmd=return_command(env,entry['start_s'],entry['intended_end_s'],8*3600.)
        for index in [0,-1]:
            t=float(raw['t'][index]);y=raw['state'][index];frame=cmd(t).as_matrix()
            a=load(env,t,y,frame,8,torque=True);b=parallel_load(env,t,y,frame,suns=8)
            out['backend_checks'].append(dict(leg=k,time_s=t,
                max_force_error_N=float(np.max(abs(a['reflected_force_N']-b['reflected_force_N']))),
                max_torque_error_N_m=float(np.max(abs(a['radiation_torque_N_m']-b['radiation_torque_per_mass']*5e6))),
                max_optical_relative_error=float(np.max(abs(a['optical']['first_intercept_bolometric_equivalent_W']-b['optical']['first_intercept_bolometric_equivalent_W']))/1e11)))
        row,sol=propagate(env,initial,entry['start_s'],entry['intended_end_s'],cmd,raw['controls'],raw['windows'],ratio,'regression_'+str(k))
        end=min(row['end_s'],entry['end_s']);t=np.linspace(entry['start_s'],end,301)
        y=sol.sol(t).T.reshape(-1,361,6);old=spline(raw)
        row['against_reviewed_tighter_run']=dict(max_position_m=float(length(y[:,:,:3]-old(t)).max()),
            max_velocity_m_s=float(length(y[:,:,3:]-old(t,1)).max()),event_time_difference_s=abs(row['end_s']-entry['end_s']))
        if k==0:
            reference,_=reference_shooting(env,21600.);r=load_raw(row)
            row['coverage']=service_coverage(env,reference,r['t'],r['state'])
        out['runs'].append(row);initial=sol.y[:,-1].reshape(-1,6)
        write('joint_regression.json',out,before,cpu)
        print('finished',k,row['end_s']/3600,row['elapsed_s'],flush=True)
    out['completed']=True;write('joint_regression.json',out,before,cpu)
if __name__=='__main__':main()
