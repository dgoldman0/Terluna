"""Energy/rating-capped ray acquisition and clearance, with free terminal states."""
import signal,resource,time
import numpy as np
from scipy.interpolate import CubicHermiteSpline
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


def base_path(env,raw,cmd,arrival,tilt):
    path=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
    switch=arrival-8*3600
    def coast(t):return np.concatenate([path(t),path(t,1)],axis=-1)
    tail=reduced_path(env,coast(switch),switch,arrival+21600,cmd,mass_ratio=raw['mass_ratio'])
    def base(t):
        t=np.atleast_1d(t);out=coast(t);after=t>=switch
        if np.any(after):out[after]=tail.sol(t[after]).T.reshape(-1,361,6)
        return out
    return base


def attitude_estimate(env,base,cmd,end,windows,mass,capacity):
    times=np.unique(np.r_[np.arange(21600,end,600.),end,43200.,windows.ravel()]);times=times[times<=end]
    nominal=rigid_square_load(cmd,times)['torque_per_mass']
    gravity=np.array([gravity_torque(env,t,base([t])[0,:,:3],cmd(t).as_matrix()) for t in times])
    specific=30000/(1.4*np.cos(np.pi/4));couple=2*abs(nominal[:,None,:]-gravity).sum(axis=-1)/10000
    powers=couple*mass*specific
    total=float(np.trapezoid(powers.sum(axis=1),times))+2.6480316453385626e12
    mask=times<=43200;prefix=float(np.trapezoid(powers[mask].sum(axis=1),times[mask]))+2.6480316453385626e12
    caps=np.array([np.maximum(0,capacity/(mass*specific)-couple[(times>=lo)&(times<=hi)].max(axis=0)) for lo,hi in windows]).T
    return dict(estimated_attitude_J=total,estimated_prefix_attitude_J=prefix,estimated_attitude_overloads=int(np.count_nonzero(powers.max(axis=0)>capacity))),caps


def main():
    signal.alarm(350);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();cpu=time.process_time()
    extra=['protection/dynamics/service_reduced.py','protection/dynamics/service_constraints.py','research/studies/solar_shield_array/natural_fit.py']
    arrays,producer=setup(extra);env=environment(10.5);ratio=arrays['mass_ratio'];mass=ratio*5e6;capacity=(ratio-1.2)*5e6*300
    parent=json.loads((HERE/'results/natural_fast_scout.json').read_text());producer['inputs']['natural_fast_scout.json']=digest(HERE/'results/natural_fast_scout.json')
    out=dict(schema='terluna.research.natural-fit/1',producer=producer,completed=False,cases=[],accepted_return=False,accepted_cycle=False,
        scope='Early independent smooth burns, no terminal position/velocity lattice. Sparse assigned-ray inequalities and selected collision branches constrain a linear response model. Estimated attitude omits shadow/eclipse torque; actual replay must enforce ratings and energy. Infeasibility only excludes this selected local assignment/branch.')
    # Widely separated opportunities in each passive attitude family.
    for family,index in [('face',0),('face',2),('edge_a',0),('edge_a',2),('edge_b',0),('edge_b',2)]:
        scout=next(x for x in parent['cases'] if x['spec']['name']==family)
        if len(scout['opportunities'])<=index:continue
        spec=scout['spec'];arrival=scout['opportunities'][index]['arrival_s'];end=arrival+21600
        row=dict(name=f'{family}_{index+1}',spec=spec,arrival_s=arrival,cadence=cadence_account(arrival,mass_kg=float(mass.sum())),iterations=[])
        raw=np.load(ROOT/scout['raw_path']);assert digest(ROOT/scout['raw_path'])==scout['raw_sha256']
        cmd=free_command(env,end,tilt=spec['tilt'],axis=spec['axis'],arrival=arrival)
        base=base_path(env,raw,cmd,arrival,spec['tilt']);windows=np.array([[6,10],[10,14]])*3600.
        maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),21600.,end,windows,np.eye(3))
        attitude,caps=attitude_estimate(env,base,cmd,end,windows,mass,capacity);row.update(attitude)
        total=1e14-attitude['estimated_attitude_J'];prefix=51.859e12-attitude['estimated_prefix_attitude_J']
        if min(total,prefix)<0 or attitude['estimated_attitude_overloads']:
            row['rejected_at']='estimated_attitude_energy_or_rating';out['cases'].append(row);write('natural_fit.json',out,before,cpu);continue
        controls=np.zeros((361,2,3));cuts=[];dates=np.unique(np.r_[np.arange(21600,end,900.),np.linspace(arrival,end,25),end])
        for iteration in range(5):
            previous=controls.copy();current=lambda ts:propose(base,maps,controls,ts)
            rays=service_rows(env,base,maps,arrival,2,assignment=current)
            new,fit=fit_free_controls(base,maps,windows,mass,caps,rays,cuts,total,prefix,wall_s=8.)
            item=dict(iteration=iteration,fit=fit);row['iterations'].append(item)
            if new is None:row['rejected_at']='local_capped_LP';break
            controls=new;ys=propose(base,maps,controls,dates)
            cuts,geometry=violations(dates,ys,cmd,margin=400.,retain=1600,history=arrays['first_service_end_state'])
            item['geometry']=geometry
            cover=[]
            for t in np.linspace(arrival,end,7):
                cover.append(local_coverage(propose(base,maps,controls,[t])[0],env.at(t),base([t])[0].mean(axis=0),limb_count=8,rotation=.137))
            item['minimum_sampled_service_coverage']=min(x['minimum_source_coverage'] for x in cover)
            item['weighted_translation_J']=float(np.sum(length(controls).sum(axis=1)*mass*30000/(1.4*np.cos(np.pi/4))))
            item['estimated_sequence_J']=item['weighted_translation_J']+attitude['estimated_attitude_J']
            if not geometry['violating_samples'] and item['minimum_sampled_service_coverage']>=1-1e-9:
                row['proposal_survived']=True;break
        if np.any(controls):
            row.update(save('fit_'+row['name'],t=dates,state=propose(base,maps,controls,dates),controls=controls,windows=windows,mass_ratio=ratio,
                base_state=base(dates),response=matrices(maps,dates,2)))
        out['cases'].append(row);write('natural_fit.json',out,before,cpu)
        print(row['name'],arrival/3600,attitude,row.get('rejected_at'),row['iterations'][-1] if row['iterations'] else {},time.monotonic()-before,flush=True)
    out['completed']=True;write('natural_fit.json',out,before,cpu)
if __name__=='__main__':main()
