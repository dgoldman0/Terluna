"""Actual stopping states and the first interval losing the numerical margin."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command
from protection.dynamics.return_escape import rolled_return,quintic_fraction
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.optical import length
from .packing_screen import load_raw
from .repair_screen import source_parents,write_product
from .fleet_run import environment,digest
from .cycling_search import HERE


def main():
    before=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    names=['repair_trial_local_refined_fine.json','repair_trial_roll_60.json',
        'repair_trial_roll_minus30.json','closure_validation.json']
    p,sources=source_parents(names,['research/studies/solar_shield_array/repair_audit.py'])
    old_stop=p['closure_validation.json']['runs'][-1]['end_s'];env=environment(3.)
    out=dict(schema='terluna.research.return-repair-audit/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        accepted_cycle=False,accepted_fleet=False,scope='Actual finite-square stopping states; conditional interval margins subtract both centres\' measured replay discrepancies. No collision is propagated and no arrival is accepted.',trials=[])
    for name in names[:-1]:
        j=p[name];raw=load_raw(j['run']);s=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        start,end=map(float,raw['t'][[0,-1]]);spec=j['schedule']
        cmd=rolled_return(return_command(env,start,j['arrival_s'],j['arrival_turn_s']),start,j['arrival_s'],
            spec.get('angle_deg',0.),spec.get('turn_s',5400.))
        row=dict(source=name,terminal_states=[])
        for ti in sorted(set([start,min(old_stop,end),end])):
            q=s(ti);v=s(ti,1);frame=cmd(ti).as_matrix();value,pair,d=closest_squares(q,frame)
            dv=(v[pair[0]]-v[pair[1]])@frame
            body_rate=((q[pair[0]]-q[pair[1]])@(cmd(ti+.1).as_matrix()-cmd(ti-.1).as_matrix()))/.2+dv
            old_d=(q[325]-q[347])@frame
            row['terminal_states'].append(dict(time_s=ti,global_minimum_m=value,pair=pair.tolist(),
                body_displacement_m=d.tolist(),body_displacement_rate_m_s=body_rate.tolist(),
                projected_centre_distance_m=float(length(d[:2])),
                projected_inscribed_disks_overlap_for_any_in_plane_roll=bool(length(d[:2])<10000.),
                pair_centres_m=q[pair].tolist(),pair_velocities_m_s=v[pair].tolist(),
                original_pair_325_347_surface_distance_m=float(length(np.maximum(abs(old_d)-[10000.,10000.,0.],0.)))))
        if spec.get('angle_deg'):
            row['roll_at_stop_deg']=float(spec['angle_deg']*quintic_fraction(end,start,spec['turn_s']))
        if 'comparison' in j:
            err=j['comparison']['position_difference_plus_reconstruction_m']
            first_bad=None;global_lower=float('inf');last=None
            for begin in np.arange(start,end,3600.):
                t=np.unique(np.r_[np.arange(begin,min(begin+3600.,end),2.),min(begin+3600.,end)])
                q=s(t);v=s(t,1);d=np.array([closest_squares(qi,f)[0] for qi,f in zip(q,cmd(t).as_matrix())])
                speed=length(v-v.mean(axis=1)[:,None,:]).max(axis=1);dt=np.diff(t)
                loss=(np.maximum(speed[:-1],speed[1:])+10000./np.sqrt(2)*np.deg2rad(.1))*dt+.15*dt**2/2+2*err
                bound=np.minimum(d[:-1],d[1:])-loss;global_lower=min(global_lower,float(bound.min()))
                bad=np.flatnonzero(bound<100.)
                if len(bad) and first_bad is None:first_bad=t[bad[0]:bad[0]+2].tolist()
                if begin<=old_stop<=t[-1]:
                    subset=np.flatnonzero(t[1:]<=old_stop)
                    last=float(bound[subset].min()) if len(subset) else None
            row['uncertainty_guard']=dict(per_centre_replay_allowance_m=err,global_conditional_lower_bound_m=global_lower,
                first_interval_below_100m_after_allowance_s=first_bad,
                acceleration_envelope_m_s2=.15,rate_envelope_deg_s=.1,
                scope='Conditional on bounded centre acceleration and rigid attitude rate; measured replay disagreement is an allowance rather than a rigorous physical uncertainty bound.')
        out['trials'].append(row)
    out['completed']=True;write_product('repair_audit.json',out,sources,before)
    print(json.dumps(out['trials']),flush=True)


if __name__=='__main__':main()
