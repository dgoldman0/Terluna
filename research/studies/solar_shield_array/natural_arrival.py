"""Bounded four-arc service-led continuation, no held intermediate geometry."""
import signal,resource,time
import numpy as np
from scipy.optimize import linear_sum_assignment
from protection.dynamics.service_search import free_command,cadence_account
from protection.dynamics.service_reduced import reduced_path
from protection.dynamics.service_constraints import service_rows
from protection.dynamics.service_arrival_controls import fit_service_arcs
from protection.dynamics.packing_control import response_maps,matrices,packing_offsets
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.joint_transfer import propose,violations
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .natural_flexible import estimate
from .natural_common import *


def main():
    signal.alarm(260);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    extra=['protection/dynamics/service_reduced.py','protection/dynamics/service_constraints.py','protection/dynamics/service_arrival_controls.py','research/studies/solar_shield_array/natural_flexible.py','research/studies/solar_shield_array/natural_arrival.py']
    arrays,producer=setup(extra);env=environment(7.);ratio=arrays['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300
    parent=json.loads((HERE/'results/natural_fast_scout.json').read_text());producer['inputs']['natural_fast_scout.json']=digest(HERE/'results/natural_fast_scout.json')
    out=dict(schema='terluna.research.natural-arrival/1',producer=producer,completed=False,cases=[],accepted_return=False,accepted_cycle=False,
        scope='Early and arrival smooth electric arcs with fixed total/prefix energy caps and individual ratings; freely propagated intervening coast and free terminal states. Grids seed ray ownership only. LP status is restricted to the selected assignments and separating branches.')
    for family,tilt,index,pitch in [('face',0.,0,8500.),('face',0.,0,9000.),('edge_a',90.,0,8500.),('edge_a',90.,2,8500.)]:
        scout=next(x for x in parent['cases'] if x['spec']['name']==family);arrival=scout['opportunities'][index]['arrival_s'];end=arrival+21600
        cmd=free_command(env,end,tilt=tilt,arrival=arrival,turn_start=6.,turn_hours=6.)
        sol=reduced_path(env,arrays['initial_epoch_state'],0.,end,cmd,mass_ratio=ratio)
        base=lambda ts:sol.sol(np.atleast_1d(ts)).T.reshape(-1,361,6)
        windows=np.array([[0.,21600.],[21600.,43200.],[arrival-28800.,arrival-14400.],[arrival-14400.,arrival]])
        maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),0.,end,windows,np.eye(3))
        attitude,caps=estimate(env,base,cmd,0.,end,windows,mass,capacity)
        total=1e14-attitude['estimated_attitude_J'];prefix=51.859e12-attitude['estimated_prefix_attitude_J']
        offsets=packing_offsets(12000.,pitch=pitch);natural=base([arrival])[0];frame=solar_frame(env.at(arrival))
        cost=length(natural[:,None,:3]-(natural[:,:3].mean(axis=0)+offsets@frame.T)[None,:,:]);_,j=linear_sum_assignment(cost);assigned=offsets[j]
        def seed(ts):
            ys=base(ts).copy()
            for k,t in enumerate(np.atleast_1d(ts)):ys[k,:,:3]=ys[k,:,:3].mean(axis=0)+assigned@solar_frame(env.at(t)).T
            return ys
        row=dict(name=f'{family}_{index+1}_{int(pitch)}',tilt=tilt,turn_start=6.,turn_hours=6.,arrival_s=arrival,seed_pitch_m=pitch,iterations=[],cadence=cadence_account(arrival,mass_kg=float(mass.sum())),**attitude)
        controls=np.zeros((361,4,3));cuts=[];dates=np.unique(np.r_[np.arange(0,end,900.),np.linspace(arrival,end,25),end])
        for iteration in range(5):
            rays=service_rows(env,base,maps,arrival,4,assignment=seed)+service_rows(env,base,maps,0.,4)
            new,fit=fit_service_arcs(base,maps,windows,mass,caps,rays,cuts,total,prefix)
            item=dict(iteration=iteration,fit=fit);row['iterations'].append(item)
            if new is None:row['rejected_at']='local_capped_LP';break
            controls=new;ys=propose(base,maps,controls,dates)
            cuts,geometry=violations(dates,ys,cmd,margin=400.,retain=1800,history=arrays['initial_epoch_state']);item['geometry']=geometry
            cover=[]
            for service_start in [0.,arrival]:
                for t in np.linspace(service_start,service_start+21600,7):
                    cover.append(local_coverage(propose(base,maps,controls,[t])[0],env.at(t),base([t])[0].mean(axis=0),limb_count=8,rotation=.137))
            item['minimum_sampled_service_coverage']=min(x['minimum_source_coverage'] for x in cover)
            item['estimated_sequence_J']=fit['translation_euclidean_J']+attitude['estimated_attitude_J']
            if not geometry['violating_samples'] and item['minimum_sampled_service_coverage']>=1-1e-9:row['proposal_survived']=True;break
        if np.any(controls):row.update(save('arrival_'+row['name'],t=dates,state=propose(base,maps,controls,dates),controls=controls,windows=windows,mass_ratio=ratio,base_state=base(dates),response=matrices(maps,dates,4)))
        out['cases'].append(row);write('natural_arrival.json',out,before,cpu)
        print(row['name'],row.get('rejected_at'),row['iterations'][-1],time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_arrival.json',out,before,cpu)
if __name__=='__main__':main()
