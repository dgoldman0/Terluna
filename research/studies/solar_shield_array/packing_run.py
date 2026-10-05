"""Coupled service and departure from a changed packing, with individual burns."""
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial import cKDTree

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_formation import swept_square_conflicts
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.parallel_shadow import ParallelPattern, minimum_clearance
from protection.dynamics.packing_control import response_maps, matrices, measured_base, controls_state
from protection.dynamics.pattern_return import smooth_arc, beam_and_rate_screen, crossing_witnesses
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.collection import require_uneclipsed_tiles
from protection.dynamics.cycling import ray_geometry
from protection.dynamics.optical import length
from .patterns_run import reference_shooting
from .packing_screen import RUN, parents, load_raw, refine
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def full_command(env,spec):
    return coast_command(env,0.,43200.,**{**spec,'delay_s':21600.+spec['delay_s']})


def geometry_audit(times,states,command):
    smallest=np.inf;witness=None
    for t,state,frame in zip(times,states,command(times).as_matrix()):
        value,pair,projection=minimum_clearance(state,frame)
        kind='parallel_square_axis'
        if pair is None:
            distances,neighbours=cKDTree(state[:,:3]).query(state[:,:3],k=2)
            i=int(np.argmin(distances[:,1]));pair=np.array([i,neighbours[i,1]])
            value=float(distances[i,1]-np.sqrt(2)*10000.)
            projection=(state[pair[0],:3]-state[pair[1],:3])@frame;kind='circumsphere'
        if value<smallest:
            smallest=value;witness=dict(time_s=float(t),pair=pair.tolist(),
                projection_m=projection.tolist(),bound_kind=kind)
    swept=swept_square_conflicts(times,states,lambda t:command(t).as_matrix()[:,[2,0,1]],
        acceleration_bound=.15,normal_rate_bound=np.deg2rad(.1))
    if not np.isfinite(swept['minimum_proven_axis_clearance_m']):
        swept['minimum_proven_axis_clearance_m']=None
        swept['no_nearby_swept_pairs']=True
    return dict(step_s=float(np.diff(times).max()),minimum_sampled_clearance_bound_m=float(smallest),
        witness=witness,swept_separation=swept)


def propagate(env,initial,command,controls,basis,windows,start,end,label,
              suns=8,step=120.,rtol=2e-10):
    before=time.monotonic();model=ParallelPattern(env,command,suns=suns)
    impulse=np.einsum('nkj,ij->nki',controls,basis);progress=[start+3600.]
    maximum=[0.]
    def rhs(t,flat):
        state=flat.reshape(-1,6)
        require_uneclipsed_tiles(ray_geometry(state[:,:3],env.at(t)))
        force=model.acceleration(t,state)+np.einsum('k,nki->ni',smooth_arc(t,windows),impulse)
        maximum[0]=max(maximum[0],float(length(force).max()))
        if t>=progress[0]:
            print(json.dumps(dict(stage=label,hours=(t-start)/3600.,elapsed_s=time.monotonic()-before)),flush=True)
            progress[0]+=3600.
        return np.c_[state[:,3:],force].ravel()
    def event(t,flat):return minimum_clearance(flat.reshape(-1,6),command(t).as_matrix())[0]-100.
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=step,rtol=rtol,
        atol=np.tile([1e-5]*3+[1e-9]*3,len(initial)),events=event,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    times=np.linspace(start,sol.t[-1],int(np.ceil((sol.t[-1]-start)/30))+1)
    states=sol.sol(times).T.reshape(-1,len(initial),6)
    audit_t=np.linspace(start,sol.t[-1],int(np.ceil((sol.t[-1]-start)/2))+1)
    audit=sol.sol(audit_t).T.reshape(-1,len(initial),6)
    geometry=geometry_audit(audit_t,audit,command)
    beam=beam_and_rate_screen(env,audit_t,audit,command)
    roots=crossing_witnesses(times,states,command,lambda t,ids:sol.sol(t).reshape(-1,6)[ids],retain=16)
    interpolated=CubicHermiteSpline(times,states[:,:,:3],states[:,:,3:])
    interpolation_error=float(length(audit[:,:,:3]-interpolated(audit_t)).max())
    thrust=float(np.max(length(controls)*2/np.diff(windows,axis=1).ravel()))
    path=RUN/(label+'.npz')
    np.savez_compressed(path,t=times,state=states,controls=controls,basis=basis,windows=windows,
        frames=command(times).as_matrix())
    complete=bool(abs(sol.t[-1]-end)<1e-6)
    passed=bool(complete and not geometry['swept_separation']['possible_conflict_pairs']
        and not roots['crossing_near_encounters'] and beam['minimum_sampled_beam_margin_deg']>0
        and beam['maximum_sampled_rate_deg_s']<=.1 and thrust<=.001 and maximum[0]<=.15)
    result=dict(start_s=start,intended_end_s=end,end_s=float(sol.t[-1]),completed=complete,
        clearance_event=bool(len(sol.t_events[0])),suns=suns,max_step_s=step,rtol=rtol,
        geometry=geometry,beam_and_rate=beam,between_sample_crossings=roots,
        maximum_evaluated_acceleration_m_s2=maximum[0],maximum_thrust_acceleration_m_s2=thrust,
        mean_scheduled_translation_m_s=float(length(controls).sum(axis=1).mean()),
        burn_count=int(np.count_nonzero(length(controls)>1e-8)),local_geometry_pass=passed,
        fine_interpolation_position_error_m=interpolation_error,
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
        elapsed_s=time.monotonic()-before,nfev=sol.nfev)
    return result,times,states


def service_coverage(env,reference,times,states):
    spline=CubicHermiteSpline(times,states[:,:,:3],states[:,:,3:])
    dates=np.unique(np.r_[np.arange(0.,21600.+1,900.),np.arange(450.,21600.,900.)])
    rows=[]
    for t in dates[dates<=times[-1]]:
        state=np.c_[spline(t),spline(t,1)]
        data=local_coverage(state,env.at(t),reference.sol(t)[:6],limb_count=16,rotation=.137,interior=8)
        data.update(time_s=float(t),window_inside_target=data['window_centre_distance_m']+10000.<=4*K.MOON_RADIUS)
        rows.append(data)
    return dict(samples=rows,minimum_coverage=min(r['minimum_source_coverage'] for r in rows),
        passed=bool(times[-1]>=21600.-1e-6 and all(r['minimum_source_coverage']>=1-1e-9 and r['window_inside_target'] for r in rows)),
        sources_per_date=25,spatial_scope='Complete polygonal receiver disk at each sampled solar source; continuous source/time certificate open.')


def main():
    p,sources=parents(['packing_screen.json']);screen=p['packing_screen.json']
    allowance=int(3600-screen['elapsed_s']-600)
    signal.alarm(max(1,allowance));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();env=environment(1.);reference,_=reference_shooting(env,21600.)
    for path in ['research/studies/solar_shield_array/packing_run.py','protection/dynamics/collection.py']:
        sources[path]=digest(ROOT/path)
    out=dict(schema='terluna.research.packing-coupled/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={'packing_screen.json':digest(HERE/'results/packing_screen.json')}),
        budget=dict(total_wall_s=3600,this_stage_allowance_s=allowance,reserved_accounting_s=600,
            maximum_coupled_departures=4,maximum_new_services=5,address_space_GiB=2),
        accepted_departure=False,accepted_fleet=False,independent_replay_pending=True,
        full_return_executed=False,handover_executed=False,recurring_cost_measured=False,
        collection_capacity_W=None,delivered_power_W=None,services=[],departures=[],corrections=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/packing.json').write_text(json.dumps(out,indent=2)+'\n')
    for index in screen['ranked_proposals']:
        if len(out['departures'])>=4 or len(out['services'])>=5:break
        case=screen['cases'][index];raw=load_raw(case);config=screen['configs'][case['config']]
        command=full_command(env,case['schedule']);controls=raw['controls'];basis=raw['basis'];windows=raw['windows']
        service,st,sy=propagate(env,raw['initial'],command,np.zeros_like(controls),basis,windows,0.,21600.,f"service_{config['index']}")
        service['config']=config;service['coverage']=service_coverage(env,reference,st,sy)
        service['service_pass']=service['local_geometry_pass'] and service['coverage']['passed']
        out['services'].append(service);save()
        print(json.dumps(dict(stage='service_gates',config=config['index'],passed=service['service_pass'],
            coverage=service['coverage']['minimum_coverage'],geometry=service['geometry'])),flush=True)
        if not service['service_pass']:continue
        predicted=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        previous=lambda t:np.concatenate([predicted(np.atleast_1d(t)),predicted(np.atleast_1d(t),1)],axis=-1)
        maps=response_maps(env,lambda t:previous([t])[0].mean(axis=0),21600.,43200.,windows,basis)
        delta=sy[-1]-previous([21600.])[0]
        base=lambda t,pred=previous,d=delta,m=maps:pred(t)+np.einsum('tij,nj->tni',matrices(m,t)[:,:,:6],d)
        for retry in range(4):
            if len(out['departures'])>=4:break
            anchor=controls.copy()
            controls,proposal,pt,fit=refine(env,base,maps,command,anchor,windows,trust=.75)
            proposal_path=RUN/f"corrected_proposal_{len(out['departures'])}.npz"
            np.savez_compressed(proposal_path,t=pt,state=proposal,controls=controls,basis=basis,windows=windows)
            fit.update(config=config['index'],schedule_case=case['schedule_case'],retry=retry,
                raw_path=str(proposal_path.relative_to(ROOT)),raw_sha256=digest(proposal_path))
            out['corrections'].append(fit);save()
            if not fit['proposal_pass']:break
            row,times,states=propagate(env,sy[-1],command,controls,basis,windows,21600.,43200.,f"departure_{len(out['departures'])}")
            row.update(config=config['index'],schedule_case=case['schedule_case'],schedule=case['schedule'])
            predict=lambda t,b=base,m=maps,u=controls.copy(),a=anchor:controls_state(b,m,t,u,a)
            new_base,error=measured_base(predict,maps,times,states)
            row['maximum_proposal_position_defect_m']=error
            out['departures'].append(row);save()
            print(json.dumps(dict(stage='departure_gates',passed=row['local_geometry_pass'],
                defect_m=error,geometry=row['geometry'])),flush=True)
            if row['local_geometry_pass']:
                out.update(accepted_departure=True,selected_screen_case=index,selected_config=config,
                    schedule=case['schedule'],stop_reason='Local service/departure gates pass; independent replay and cost accounting pending.')
                break
            base=new_base
        if out['accepted_departure']:break
    out['completed']=True
    if not out['accepted_departure']:out['stop_reason']='No coupled low-energy departure accepted in the tested local branches.'
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    save()


if __name__=='__main__':main()
