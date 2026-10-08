"""Earlier service-preserving burns and slower sail turns under the same cost cap."""
import signal,resource,time
import numpy as np
from protection.dynamics.service_search import free_command,cadence_account
from protection.dynamics.service_reduced import reduced_path
from protection.dynamics.service_constraints import service_rows,fit_free_controls
from protection.dynamics.packing_control import response_maps,matrices
from protection.dynamics.joint_transfer import propose,violations
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.tile_torque import gravity_torque
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .natural_common import *


def estimate(env,base,cmd,start,end,windows,mass,capacity):
    times=np.unique(np.r_[np.arange(start,end,300.),end,43200.,windows.ravel()])
    nominal=rigid_square_load(cmd,times)['torque_per_mass']
    gravity=np.array([gravity_torque(env,t,base([t])[0,:,:3],cmd(t).as_matrix()) for t in times])
    specific=30000/(1.4*np.cos(np.pi/4));couple=2*abs(nominal[:,None,:]-gravity).sum(axis=-1)/10000
    powers=couple*mass*specific
    before=2.6480316453385626e12 if start else 0.
    total=float(np.trapezoid(powers.sum(axis=1),times))+before
    mask=times<=43200;prefix=float(np.trapezoid(powers[mask].sum(axis=1),times[mask]))+before
    caps=np.array([np.maximum(0,capacity/(mass*specific)-couple[(times>=lo)&(times<=hi)].max(axis=0)) for lo,hi in windows]).T
    return dict(estimated_attitude_J=total,estimated_prefix_attitude_J=prefix,estimated_attitude_overloads=int(np.count_nonzero(powers.max(axis=0)>capacity))),caps


def main():
    signal.alarm(330);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['protection/dynamics/service_reduced.py','protection/dynamics/service_constraints.py','research/studies/solar_shield_array/natural_flexible.py'])
    env=environment(10.5);ratio=arrays['mass_ratio'];mass=ratio*5e6;capacity=(ratio-1.2)*5e6*300
    parent=json.loads((HERE/'results/natural_fast_scout.json').read_text());producer['inputs']['natural_fast_scout.json']=digest(HERE/'results/natural_fast_scout.json')
    out=dict(schema='terluna.research.natural-flexible/1',producer=producer,completed=False,cases=[],accepted_return=False,accepted_cycle=False,
        scope='Same 100 TJ sequence and 51.859 TJ first-twelve-hour caps. Two early smooth arcs, including service-period control where specified, with ray constraints on both service passages and free terminal states. Slow feather turns free attitude budget. Selected assignments and separating branches remain local restrictions.')
    specs=[dict(name='face_epoch',family='face',index=0,start=0.,tilt=0.,turn_start=6.,turn_hours=6.),
        dict(name='slow_edge_epoch',family='edge_a',index=0,start=0.,tilt=90.,turn_start=6.,turn_hours=6.),
        dict(name='slow_edge_after_service',family='edge_a',index=0,start=21600.,tilt=90.,turn_start=6.,turn_hours=6.),
        dict(name='slow_edge_third',family='edge_a',index=2,start=0.,tilt=90.,turn_start=6.,turn_hours=6.)]
    for spec in specs:
        scout=next(x for x in parent['cases'] if x['spec']['name']==spec['family']);arrival=scout['opportunities'][spec['index']]['arrival_s'];end=arrival+21600;start=spec['start']
        cmd=free_command(env,end,tilt=spec['tilt'],arrival=arrival,turn_start=spec['turn_start'],turn_hours=spec['turn_hours'])
        initial=arrays['initial_epoch_state' if start==0 else 'first_service_end_state']
        sol=reduced_path(env,initial,start,end,cmd,mass_ratio=ratio)
        base=lambda ts:sol.sol(np.atleast_1d(ts)).T.reshape(-1,361,6)
        windows=np.array([[0,6],[6,12]])*3600 if start==0 else np.array([[6,9],[9,12]])*3600
        maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),start,end,windows,np.eye(3))
        attitude,caps=estimate(env,base,cmd,start,end,windows,mass,capacity)
        row=dict(spec=spec,arrival_s=arrival,cadence=cadence_account(arrival,mass_kg=float(mass.sum())),iterations=[],**attitude)
        total=1e14-attitude['estimated_attitude_J'];prefix=51.859e12-attitude['estimated_prefix_attitude_J']
        if min(total,prefix)<0 or attitude['estimated_attitude_overloads']:
            row['rejected_at']='estimated_attitude_energy_or_rating'
        else:
            controls=np.zeros((361,2,3));cuts=[];dates=np.unique(np.r_[np.arange(start,end,900.),np.linspace(arrival,end,25),end])
            for iteration in range(6):
                current=lambda ts:propose(base,maps,controls,ts)
                rays=service_rows(env,base,maps,arrival,2,assignment=current)
                if start==0:rays+=service_rows(env,base,maps,0.,2,assignment=current)
                # Both arcs end by hour 12, so all translation must satisfy both
                # available-energy ceilings; no early burn is hidden in a prefix.
                new,fit=fit_free_controls(base,maps,windows,mass,caps,rays,cuts,min(total,prefix),prefix,wall_s=8.)
                item=dict(iteration=iteration,fit=fit);row['iterations'].append(item)
                if new is None:row['rejected_at']='local_capped_LP';break
                controls=new;ys=propose(base,maps,controls,dates)
                cuts,geometry=violations(dates,ys,cmd,margin=400.,retain=1600,history=initial);item['geometry']=geometry
                cover=[]
                for service_start in ([0.,arrival] if start==0 else [arrival]):
                    for t in np.linspace(service_start,service_start+21600,7):
                        cover.append(local_coverage(propose(base,maps,controls,[t])[0],env.at(t),base([t])[0].mean(axis=0),limb_count=8,rotation=.137))
                item['minimum_sampled_service_coverage']=min(x['minimum_source_coverage'] for x in cover)
                item['weighted_translation_J']=float(np.sum(length(controls).sum(axis=1)*mass*30000/(1.4*np.cos(np.pi/4))))
                item['estimated_sequence_J']=item['weighted_translation_J']+attitude['estimated_attitude_J']
                if not geometry['violating_samples'] and item['minimum_sampled_service_coverage']>=1-1e-9:row['proposal_survived']=True;break
            if np.any(controls):row.update(save(spec['name'],t=dates,state=propose(base,maps,controls,dates),controls=controls,windows=windows,mass_ratio=ratio,base_state=base(dates),response=matrices(maps,dates,2)))
        out['cases'].append(row);write('natural_flexible.json',out,before,cpu)
        print(spec['name'],attitude,row.get('rejected_at'),row['iterations'][-1] if row['iterations'] else {},time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_flexible.json',out,before,cpu)
if __name__=='__main__':main()
