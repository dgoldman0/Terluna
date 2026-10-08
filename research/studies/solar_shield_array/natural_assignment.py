"""Coverage-feasible ray-bundle seeds; terminal states remain free variables."""
import signal,resource,time
import numpy as np
from scipy.optimize import linear_sum_assignment
from protection.dynamics.service_search import free_command,cadence_account
from protection.dynamics.service_reduced import reduced_path
from protection.dynamics.service_constraints import service_rows,fit_free_controls
from protection.dynamics.packing_control import response_maps,matrices,packing_offsets
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.joint_transfer import propose,violations
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .natural_flexible import estimate
from .natural_common import *


def main():
    signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['protection/dynamics/service_reduced.py','protection/dynamics/service_constraints.py','research/studies/solar_shield_array/natural_flexible.py','research/studies/solar_shield_array/natural_assignment.py'])
    env=environment(3.);ratio=arrays['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300
    parent=json.loads((HERE/'results/natural_fast_scout.json').read_text());producer['inputs']['natural_fast_scout.json']=digest(HERE/'results/natural_fast_scout.json')
    out=dict(schema='terluna.research.natural-assignment/1',producer=producer,completed=False,cases=[],accepted_return=False,accepted_cycle=False,
        scope='Use 8.5/9 km grids only to seed coherent ray ownership; no positions or velocities are imposed on the optimization. The 100 TJ and first-twelve-hour caps remain. A separate unbudgeted LP diagnoses a local branch and never becomes a candidate.')
    for family,tilt in [('face',0.),('edge_a',90.)]:
        scout=next(x for x in parent['cases'] if x['spec']['name']==family);arrival=scout['opportunities'][0]['arrival_s'];end=arrival+21600
        cmd=free_command(env,end,tilt=tilt,arrival=arrival,turn_start=6.,turn_hours=6.)
        sol=reduced_path(env,arrays['initial_epoch_state'],0.,end,cmd,mass_ratio=ratio)
        base=lambda ts:sol.sol(np.atleast_1d(ts)).T.reshape(-1,361,6)
        windows=np.array([[0.,6.],[6.,12.]])*3600
        maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),0.,end,windows,np.eye(3))
        attitude,caps=estimate(env,base,cmd,0.,end,windows,mass,capacity)
        total=1e14-attitude['estimated_attitude_J'];prefix=51.859e12-attitude['estimated_prefix_attitude_J']
        for pitch in [8500.,9000.]:
            offsets=packing_offsets(12000.,pitch=pitch);frame=solar_frame(env.at(arrival))
            natural=base([arrival])[0];centroid=natural.mean(axis=0)
            cost=length(natural[:,None,:3]-(centroid[:3]+offsets@frame.T)[None,:,:])
            i,j=linear_sum_assignment(cost);assigned=offsets[j]
            def seed(ts):
                ys=base(ts).copy()
                for k,t in enumerate(np.atleast_1d(ts)):
                    ys[k,:,:3]=ys[k,:,:3].mean(axis=0)+assigned@solar_frame(env.at(t)).T
                return ys
            controls=np.zeros((361,2,3));cuts=[];dates=np.unique(np.r_[np.arange(0,end,900.),np.linspace(arrival,end,25),end])
            row=dict(name=family+'_'+str(int(pitch)),tilt=tilt,turn_start=6.,turn_hours=6.,arrival_s=arrival,seed_pitch_m=pitch,iterations=[],**attitude)
            if min(total,prefix)<0:row['rejected_at']='estimated_attitude_budget'
            else:
                for iteration in range(4):
                    rays=service_rows(env,base,maps,arrival,2,assignment=seed)
                    rays+=service_rows(env,base,maps,0.,2)
                    new,fit=fit_free_controls(base,maps,windows,mass,caps,rays,cuts,min(total,prefix),prefix,wall_s=8.)
                    item=dict(iteration=iteration,fit=fit);row['iterations'].append(item)
                    if new is None:
                        row['rejected_at']='local_capped_LP'
                        if iteration==0:
                            unlimited,diag=fit_free_controls(base,maps,windows,mass,caps,rays,(),1e17,1e17,wall_s=8.)
                            row['unbudgeted_branch_diagnostic']=diag
                            if unlimited is not None:
                                row['unbudgeted_translation_euclidean_J']=float(np.sum(length(unlimited).sum(axis=1)*mass*30000/(1.4*np.cos(np.pi/4))))
                        break
                    controls=new;ys=propose(base,maps,controls,dates)
                    cuts,geometry=violations(dates,ys,cmd,margin=400.,retain=1600,history=arrays['initial_epoch_state']);item['geometry']=geometry
                    cover=[]
                    for service_start in [0.,arrival]:
                        for t in np.linspace(service_start,service_start+21600,7):
                            cover.append(local_coverage(propose(base,maps,controls,[t])[0],env.at(t),base([t])[0].mean(axis=0),limb_count=8,rotation=.137))
                    item['minimum_sampled_service_coverage']=min(x['minimum_source_coverage'] for x in cover)
                    item['weighted_translation_J']=float(np.sum(length(controls).sum(axis=1)*mass*30000/(1.4*np.cos(np.pi/4))))
                    item['estimated_sequence_J']=item['weighted_translation_J']+attitude['estimated_attitude_J']
                    if not geometry['violating_samples'] and item['minimum_sampled_service_coverage']>=1-1e-9:row['proposal_survived']=True;break
                if np.any(controls):row.update(save('assignment_'+row['name'],t=dates,state=propose(base,maps,controls,dates),controls=controls,windows=windows,mass_ratio=ratio,base_state=base(dates),response=matrices(maps,dates,2)))
            out['cases'].append(row);write('natural_assignment.json',out,before,cpu)
            print(row['name'],row.get('rejected_at'),row.get('unbudgeted_branch_diagnostic'),row['iterations'][-1] if row['iterations'] else {},time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_assignment.json',out,before,cpu)
if __name__=='__main__':main()
