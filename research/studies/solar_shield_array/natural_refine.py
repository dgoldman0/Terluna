"""Retain old encounter cuts and add missing receiver rays on one cheap branch."""
import signal,resource,time
import numpy as np
from protection.dynamics.service_search import free_command,cadence_account
from protection.dynamics.service_reduced import reduced_path
from protection.dynamics.service_arrival_controls import fit_service_arcs
from protection.dynamics.service_holes import uncovered_rows
from protection.dynamics.packing_control import response_maps,matrices
from protection.dynamics.joint_transfer import propose,violations
from .natural_flexible import estimate
from . import natural_dynamic_rays as dynamic
from .natural_common import *


def main():
    signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=dynamic.setup(['protection/dynamics/service_arrival_controls.py','protection/dynamics/service_holes.py','research/studies/solar_shield_array/natural_flexible.py','research/studies/solar_shield_array/natural_refine.py'])
    parent=json.loads((HERE/'results/natural_dynamic_rays.json').read_text());producer['inputs']['natural_dynamic_rays.json']=digest(HERE/'results/natural_dynamic_rays.json')
    row=parent['cases'][0];arrival=row['arrival_s'];end=arrival+21600;env=environment(3.);ratio=arrays['mass_ratio'];mass=5e6*ratio;capacity=(ratio-1.2)*5e6*300
    raw=np.load(ROOT/row['raw_path']);assert digest(ROOT/row['raw_path'])==row['raw_sha256'];controls=raw['controls'].copy();windows=raw['windows']
    cmd=free_command(env,end,tilt=0.,arrival=arrival,turn_start=6.,turn_hours=6.)
    sol=reduced_path(env,arrays['initial_epoch_state'],0.,end,cmd,mass_ratio=ratio)
    base=lambda ts:sol.sol(np.atleast_1d(ts)).T.reshape(-1,361,6)
    maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),0.,end,windows,np.eye(3))
    attitude,caps=estimate(env,base,cmd,0.,end,windows,mass,capacity)
    total=1e14-attitude['estimated_attitude_J'];prefix=51.859e12-attitude['estimated_prefix_attitude_J']
    dates=np.unique(np.r_[np.arange(0,end,600.),np.linspace(arrival,end,37),end]);cut_archive={}
    rays=dynamic.rows(env,base,maps,arrival,4)+dynamic.rows(env,base,maps,0.,4)
    out=dict(schema='terluna.research.natural-refine/1',producer=producer,completed=False,accepted_return=False,accepted_cycle=False,case='face_refined',arrival_s=arrival,tilt=0.,turn_start=6.,turn_hours=6.,iterations=[],**attitude,
        scope='One selected ray/encounter branch; retain every inserted separation cut, add largest receiver holes, keep both energy caps and individual ratings. No nonlinear dynamics acceptance from this linear proposal.')
    for iteration in range(14):
        current=lambda ts:propose(base,maps,controls,ts)
        candidates,geometry=violations(dates,current(dates),cmd,margin=400.,retain=100000,history=arrays['initial_epoch_state'])
        # Insert earliest encounters first, plus every previously imposed cut.
        candidates=sorted(candidates,key=lambda x:x[0])[:1200]
        for cut in candidates:
            t,i,j,axis,required=cut;key=(t,i,j,tuple(np.round(axis,10)));cut_archive[key]=cut
        newrays,coverage=uncovered_rows(env,base,current,maps,arrival,4)
        first_rays,first_coverage=uncovered_rows(env,base,current,maps,0.,4)
        rays.extend(newrays);rays.extend(first_rays)
        item=dict(iteration=iteration,geometry_before=geometry,next_service_before=coverage,first_service_before=first_coverage,retained_encounter_cuts=len(cut_archive),retained_ray_inequalities=len(rays))
        if not geometry['violating_samples'] and min(coverage['minimum_source_coverage'],first_coverage['minimum_source_coverage'])>=1-1e-9:
            out['proposal_survived']=True;out['iterations'].append(item);break
        new,fit=fit_service_arcs(base,maps,windows,mass,caps,rays,list(cut_archive.values()),total,prefix,wall_s=8.)
        item['fit']=fit;out['iterations'].append(item)
        if new is None:out['rejected_at']='retained_local_cuts_or_ray_assignments';break
        controls=new;out['last_estimated_sequence_J']=fit['translation_euclidean_J']+attitude['estimated_attitude_J']
        write('natural_refine.json',out,before,cpu)
        print(iteration,geometry['violating_samples'],coverage['minimum_source_coverage'],fit.get('translation_euclidean_J'),time.monotonic()-before,flush=True)
    out.update(save('refined_face',t=dates,state=propose(base,maps,controls,dates),controls=controls,windows=windows,mass_ratio=ratio,base_state=base(dates),response=matrices(maps,dates,4)))
    final=lambda ts:propose(base,maps,controls,ts)
    _,out['final_geometry']=violations(dates,final(dates),cmd,margin=400.,retain=1,history=arrays['initial_epoch_state'])
    _,out['final_next_service']=uncovered_rows(env,base,final,maps,arrival,4)
    out['completed']=True;write('natural_refine.json',out,before,cpu)
if __name__=='__main__':main()
