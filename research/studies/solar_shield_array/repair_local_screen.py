"""Twelfth proposal: preserve initial separating branches with local burns.

No terminal lattice is imposed. The proposed correction must be propagated,
and its terminal state debt remains an explicit obligation.
"""
import json
import resource
import signal
import time

import numpy as np
from scipy.spatial import cKDTree

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command, independent_path, corrected_state
from protection.dynamics.packing_control import response_maps, matrices, independent_impulses
from protection.dynamics.return_escape import applied_impulse
from protection.dynamics.pattern_return import crossing_witnesses, beam_and_rate_screen
from protection.dynamics.optical import length
from .packing_screen import load_raw
from .repair_screen import RUN, source_parents, write_product, min_geometry
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def main():
    before=time.monotonic();names=['repair_screen.json','closure_validation.json']
    p,sources=source_parents(names,['research/studies/solar_shield_array/repair_local_screen.py'])
    allowance=int(600-p['repair_screen.json']['elapsed_s'])
    if allowance<10:raise ValueError('No local-screen budget remains')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    old=load_raw(p['repair_screen.json']['cases'][0]);initial=old['initial'];start=old['t'][0]
    end=64800.;env=environment(3.);cmd=return_command(env,start,old['t'][-1],28800.)
    free=independent_path(env,initial,start,end,cmd,mass_ratio=old['mass_ratio'])
    base=lambda t:free.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
    windows=np.array([[start,start+7200.],[start+7200.,start+14400.]])
    maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),start,end,windows,np.eye(3))
    controls=np.zeros((361,2,3));controls[:,0]=old['controls'][:,0]
    out=dict(schema='terluna.research.return-repair-local-screen/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(stage_wall_s=allowance,threads=1,address_space_GiB=2),
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        scope='Two independent early arcs, initial separating-axis branches and no imposed terminal lattice. Gravity response is a proposal only. Later original burns remain as a comparison and require refitting.',
        iterations=[])
    times=np.arange(start,end+1,60.);cuts={};first_frame=cmd(start).as_matrix()
    for iteration in range(3):
        states=corrected_state(base,maps,controls,times)
        for t,state,frame,b in zip(times,states,cmd(times).as_matrix(),matrices(maps,times)[:,:3,6:]):
            pairs=cKDTree(state[:,:3]).query_pairs(15000.,output_type='ndarray')
            d=(state[pairs[:,0],:3]-state[pairs[:,1],:3])@frame
            gap=abs(d)-[10000.,10000.,0.]
            for pair in pairs[np.max(gap,axis=1)<500.]:
                i,j=map(int,pair)
                d0=(initial[i,:3]-initial[j,:3])@first_frame
                axis=int(np.argmax(abs(d0)-[10000.,10000.,0.]))
                sign=1. if d0[axis]>=0 else -1.
                gradient=sign*frame[:,axis]@b
                actual=sign*frame[:,axis]@(state[i,:3]-state[j,:3])
                required=500.+(10000. if axis<2 else 0.)
                cuts[(float(t),i,j,axis,sign)]=dict(i=i,j=j,gradient=gradient,
                    rhs=required-actual+gradient@(controls[i]-controls[j]).ravel())
        controls,fit=independent_impulses(controls,list(cuts.values()),windows,trust=2.,wall_s=15.)
        row=dict(iteration=iteration,fit=fit);out['iterations'].append(row)
        write_product('repair_local_screen.json',out,sources,before)
        if not fit['success']:break
        states=corrected_state(base,maps,controls,times)
        row['geometry']=min_geometry(times,states,cmd)
        if row['geometry']['minimum_sampled_surface_distance_m']>450.:break
    if fit['success']:
        # Nonlinear independent-member screen of the changed early burns.
        u=np.concatenate([controls,old['controls'][:,1:]],axis=1)
        w=np.vstack([windows,old['windows'][1:]])
        sol=independent_path(env,initial,start,old['t'][-1],cmd,u,w,old['mass_ratio'])
        state=lambda t:sol.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
        fine=np.arange(start,end+1,20.)
        geometry=min_geometry(fine,state(fine),cmd)
        roots=crossing_witnesses(fine,state(fine),cmd,lambda t,ids:state([t])[0,ids],retain=12)
        beam=beam_and_rate_screen(env,fine,state(fine),cmd)
        path=RUN/'local_proposal.npz';t=old['t']
        np.savez_compressed(path,t=t,state=state(t),initial=initial,controls=u,windows=w,
            mass_ratio=old['mass_ratio'],frames=cmd(t).as_matrix())
        delta=state([t[-1]])[0]-old['state'][-1]
        out.update(early_geometry=geometry,early_crossings=roots,beam_and_rate=beam,
            raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
            schedule=dict(label='local'),mean_first_six_hour_translation_m_s=float(applied_impulse(u,w,start,end).mean()),
            endpoint_difference_from_unchanged_m=float(length(delta[:,:3]).max()),
            endpoint_velocity_difference_from_unchanged_m_s=float(length(delta[:,3:]).max()),
            early_proposal_pass=bool(geometry['minimum_sampled_surface_distance_m']>350. and
                not roots['crossing_near_encounters'] and beam['minimum_sampled_beam_margin_deg']>0))
    out['completed']=True;write_product('repair_local_screen.json',out,sources,before)
    print(json.dumps({k:v for k,v in out.items() if k not in ['producer','iterations']}),flush=True)


if __name__=='__main__':main()
