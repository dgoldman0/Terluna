"""Continue the local proposal with a compact set of encountered constraints."""
import json
import resource
import signal
import time

import numpy as np
from scipy.spatial import cKDTree

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command,independent_path,corrected_state
from protection.dynamics.packing_control import response_maps,matrices,independent_impulses
from protection.dynamics.return_escape import applied_impulse
from protection.dynamics.pattern_return import crossing_witnesses,beam_and_rate_screen
from protection.dynamics.optical import length
from .packing_screen import load_raw
from .repair_screen import RUN,source_parents,write_product,min_geometry
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE


def main():
    before=time.monotonic();names=['repair_screen.json','repair_local_screen.json']
    p,sources=source_parents(names,['research/studies/solar_shield_array/repair_refine.py'])
    allowance=int(600-sum(p[n]['elapsed_s'] for n in names));signal.alarm(allowance)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    old=load_raw(p['repair_screen.json']['cases'][0]);initial=old['initial'];start=old['t'][0];end=64800.
    env=environment(3.);cmd=return_command(env,start,old['t'][-1],28800.)
    free=independent_path(env,initial,start,end,cmd,mass_ratio=old['mass_ratio'])
    base=lambda t:free.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
    windows=np.array([[start,start+7200.],[start+7200.,start+14400.]])
    maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),start,end,windows,np.eye(3))
    u=np.zeros((361,2,3));u[:,0]=old['controls'][:,0]
    out=dict(schema='terluna.research.return-repair-local-refinement/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(stage_wall_s=allowance,threads=1,address_space_GiB=2),
        scope='Continuation of the twelfth local proposal, compacting redundant encounter dates. Initial separating branches and free terminal states. Full independent propagation and coupled validation still required.',
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,iterations=[],schedule=dict(label='local_refined'))
    times=np.arange(start,end+1,60.);cuts={};first_frame=cmd(start).as_matrix()
    for iteration in range(4):
        y=corrected_state(base,maps,u,times);new={}
        for t,state,frame,b in zip(times,y,cmd(times).as_matrix(),matrices(maps,times)[:,:3,6:]):
            pairs=cKDTree(state[:,:3]).query_pairs(15000.,output_type='ndarray')
            d=(state[pairs[:,0],:3]-state[pairs[:,1],:3])@frame
            for pair in pairs[np.max(abs(d)-[10000.,10000.,0.],axis=1)<600.]:
                i,j=map(int,pair);d0=(initial[i,:3]-initial[j,:3])@first_frame
                axis=int(np.argmax(abs(d0)-[10000.,10000.,0.]));sign=1. if d0[axis]>=0 else -1.
                g=sign*frame[:,axis]@b;actual=sign*frame[:,axis]@(state[i,:3]-state[j,:3])
                deficit=600.+(10000. if axis<2 else 0.)-actual;key=(i,j,axis,sign)
                if key not in new or new[key]['deficit']<deficit:
                    new[key]=dict(i=i,j=j,axis=axis,sign=sign,t=float(t),deficit=deficit,
                        gradient=g,rhs=deficit+g@(u[i]-u[j]).ravel())
        for key,c in new.items():cuts[(*key,c['t'])]=c
        previous=u.copy();u,fit=independent_impulses(u,list(cuts.values()),windows,trust=2.,wall_s=20.)
        row=dict(iteration=iteration,fit=fit);out['iterations'].append(row)
        if not fit['success']:u=previous;break
        row['geometry']=min_geometry(times,corrected_state(base,maps,u,times),cmd)
        write_product('repair_local_refined.json',out,sources,before)
        if row['geometry']['minimum_sampled_surface_distance_m']>450.:break
    # Retain the final feasible LP iterate even when a later solve fails.
    controls=np.concatenate([u,old['controls'][:,1:]],axis=1);w=np.vstack([windows,old['windows'][1:]])
    sol=independent_path(env,initial,start,old['t'][-1],cmd,controls,w,old['mass_ratio'])
    state=lambda t:sol.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
    fine=np.arange(start,end+1,20.)
    geometry=min_geometry(fine,state(fine),cmd)
    roots=crossing_witnesses(fine,state(fine),cmd,lambda t,ids:state([t])[0,ids],retain=12)
    beam=beam_and_rate_screen(env,fine,state(fine),cmd)
    path=RUN/'local_refined_proposal.npz';t=old['t']
    np.savez_compressed(path,t=t,state=state(t),initial=initial,controls=controls,windows=w,
        mass_ratio=old['mass_ratio'],frames=cmd(t).as_matrix())
    delta=state([t[-1]])[0]-old['state'][-1]
    out.update(early_geometry=geometry,early_crossings=roots,beam_and_rate=beam,
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
        mean_first_six_hour_translation_m_s=float(applied_impulse(controls,w,start,end).mean()),
        endpoint_difference_from_unchanged_m=float(length(delta[:,:3]).max()),
        endpoint_velocity_difference_from_unchanged_m_s=float(length(delta[:,3:]).max()),
        early_proposal_pass=bool(geometry['minimum_sampled_surface_distance_m']>350. and not roots['crossing_near_encounters']),completed=True)
    write_product('repair_local_refined.json',out,sources,before)
    print(json.dumps({k:v for k,v in out.items() if k not in ['producer','iterations','early_crossings']}),flush=True)


if __name__=='__main__':main()
