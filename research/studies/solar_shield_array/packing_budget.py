"""Actual torque, individual maneuver totals and optical inventory of the pilot."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline, interp1d
from scipy.optimize import brentq

from shared.provenance import constants_used
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.packing_control import packing_offsets
from protection.dynamics.tile_torque import parallel_load, gravity_torque
from protection.dynamics.collection import integrate_inventory, uneclipsed_collection
from protection.dynamics.attitude_load import rigid_square_load, disturbance_torque_bound
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.coverage_pilot import conditional_cost
from protection.dynamics.optical import length
from .packing_screen import parents, load_raw, RUN
from .packing_run import full_command
from .patterns_run import reference_shooting
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO


def account(times,gravity,radiation,command,controls,windows):
    dense=np.unique(np.r_[np.arange(times[0],times[-1],2.),times[-1]])
    nominal=rigid_square_load(command,dense)
    disturbance=interp1d(times,gravity+radiation,axis=0,assume_sorted=True)(dense)
    couples=2*np.abs(nominal['torque_per_mass'][:,None,:]-disturbance).sum(axis=-1)/10000.
    shape=np.array([smooth_arc(t,windows) for t in dense])
    norms=length(controls)
    # Burns do not overlap: each vector's norm scales with its nonnegative arc.
    force=couples+shape@norms.T
    def primitive(t):
        x=np.clip((t-windows[:,0])/np.diff(windows,axis=1).ravel(),0.,1.)
        return x-np.sin(2*np.pi*x)/(2*np.pi)
    translation=float((norms@(primitive(times[-1])-primitive(times[0]))).mean())
    attitude=float(np.trapezoid(couples.mean(axis=1),dense))
    prop=SCENARIO['propulsion'];mass=SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    conversion=mass*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*np.cos(np.deg2rad(prop['cant_deg'])))
    peak_members=force.max(axis=0)*conversion
    return dict(translation_m_s=translation,attitude_equivalent_m_s=attitude,
        total_equivalent_m_s=translation+attitude,
        nominal_attitude_m_s=float(np.trapezoid(nominal['force_per_mass'],dense)),
        propulsion_energy_J=(translation+attitude)*len(controls)*conversion,
        sampled_peak_fleet_propulsion_W=float(force.sum(axis=1).max()*conversion),
        sum_individual_peak_propulsion_W=float(peak_members.sum()),
        maximum_individual_peak_propulsion_W=float(peak_members.max()),
        per_member_peak_propulsion_W=peak_members.tolist(),
        scheduled_burn_count=int(np.count_nonzero(norms>1e-8)),
        gravity_only_equivalent_m_s=float(np.trapezoid((2*abs(gravity).sum(axis=-1)/10000.).mean(axis=1),times)),
        radiation_only_equivalent_m_s=float(np.trapezoid((2*abs(radiation).sum(axis=-1)/10000.).mean(axis=1),times)))


def main():
    names=['packing_screen.json','packing.json','packing_validation.json','energy_budget.json','pilot_budget.json']
    p,sources=parents(names);validation=p['packing_validation.json']
    if not validation['accepted_service_departure']:raise ValueError('Requires independently replayed local candidate')
    restart=ROOT/'research/runs/solar_shield_array/packing_replay_restart.json'
    interrupted=validation['budget']['interrupted_partial_replay_wall_s']
    if restart.exists():interrupted=max(interrupted,json.loads(restart.read_text()).get('charged_interrupted_wall_s',0))
    used=sum(p[n]['elapsed_s'] for n in names[:3])+interrupted
    allowance=min(600,int(3600-used))
    if allowance<=0:raise ValueError('No accounting budget remains')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();env=environment(1.);command=full_command(env,validation['schedule'])
    for path in ['research/studies/solar_shield_array/packing_budget.py','protection/dynamics/tile_torque.py']:
        sources[path]=digest(ROOT/path)
    out=dict(schema='terluna.research.packing-budget/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names},raw_inputs={}),
        budget=dict(prior_numerical_wall_s=used,this_stage_wall_s=allowance,numerical_threads=1,address_space_GiB=2),
        accepted_local_service_departure=True,accepted_fleet=False,hardware_mass_feedback_propagated=False,
        full_return_executed=False,handover_executed=False,recurring_cost_measured=False,
        measured_recurring_impulse_m_s_day=None,fleet_propulsion_W=None,collection_capacity_W=None,delivered_power_W=None,
        scope='Tighter replay of the same 361 members through six-hour service and six-hour departure. Uniform rigid 50 g/m2 optical mass, ideal electric edge couples, finite-area gravity and reflected-band shadow moments. Acquisition, full return, handover, storage, flexibility, recovery and recurrence remain unpaid.',
        datasets={},rejected_optical_prefixes=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/packing_budget.json').write_text(json.dumps(out,indent=2)+'\n')
    for label,entry in zip(['service','departure'],validation['runs']):
        raw=load_raw(entry);out['producer']['raw_inputs'][entry['raw_path']]=entry['raw_sha256']
        spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        start,end=map(float,raw['t'][[0,-1]]);evaluations=[]
        agreements=dict(maximum_collection_relative_error=0.,maximum_gravity_quadrature_error_N_m_kg=0.,maximum_envelope_excess_N_m_kg=0.)
        for suns,step,rotation in [(16,60.,.317),(32,120.,.317+np.pi/32)]:
            times=np.unique(np.r_[np.arange(start,end,step),end]);gravity=[];radiation=[];light={}
            for k,t in enumerate(times):
                state=np.c_[spline(t),spline(t,1)];frame=command(t).as_matrix()
                loads=parallel_load(env,t,state,frame,suns=suns,rotation=rotation)
                gravity.append(loads['gravity_torque_per_mass']);radiation.append(loads['radiation_torque_per_mass'])
                for name,value in loads['optical'].items():light.setdefault(name,[]).append(value)
                if suns==16 and k%60==0:
                    check=uneclipsed_collection(env,t,state,frame,suns=16)
                    value=loads['optical']['first_intercept_bolometric_equivalent_W']
                    agreements['maximum_collection_relative_error']=max(agreements['maximum_collection_relative_error'],
                        float(np.max(abs(value-check['first_intercept_bolometric_equivalent_W'])))/max(float(value.max()),1.))
                    g5=gravity_torque(env,t,state[:,:3],frame,order=5)
                    agreements['maximum_gravity_quadrature_error_N_m_kg']=max(agreements['maximum_gravity_quadrature_error_N_m_kg'],
                        float(length(g5-loads['gravity_torque_per_mass']).max()))
                    gb,rb,_=disturbance_torque_bound(env,t,state,frame)
                    agreements['maximum_envelope_excess_N_m_kg']=max(agreements['maximum_envelope_excess_N_m_kg'],
                        float((length(loads['gravity_torque_per_mass'])-gb).max()),float((length(loads['radiation_torque_per_mass'])-rb).max()))
            gravity=np.array(gravity);radiation=np.array(radiation);light={k:np.array(v) for k,v in light.items()}
            cost=account(times,gravity,radiation,command,raw['controls'],raw['windows'])
            summary={}
            for name,value in light.items():
                data=integrate_inventory(times,value);data.pop('per_member_energy_J');summary[name]=data
            path=RUN/f'{label}_loads_{suns}.npz'
            np.savez_compressed(path,t=times,member_id=np.arange(len(raw['controls'])),gravity_torque_N_m_kg=gravity,
                radiation_torque_N_m_kg=radiation,**light)
            result=dict(sources=suns,time_step_s=step,cost=cost,optical=summary,
                raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
            if suns==16:
                coarse=account(times[::2],gravity[::2],radiation[::2],command,raw['controls'],raw['windows'])
                result['time_coarsening_relative_attitude_change']=abs(coarse['attitude_equivalent_m_s']-cost['attitude_equivalent_m_s'])/cost['attitude_equivalent_m_s']
                first=light['first_intercept_bolometric_equivalent_W']
                coarse_light=integrate_inventory(times[::2],first[::2])['energy_J']
                result['time_coarsening_relative_optical_energy_change']=abs(coarse_light-summary['first_intercept_bolometric_equivalent_W']['energy_J'])/coarse_light
            else:
                result['source_refinement_relative_attitude_change']=abs(coarse['attitude_equivalent_m_s']-cost['attitude_equivalent_m_s'])/cost['attitude_equivalent_m_s']
                result['source_refinement_relative_optical_energy_change']=abs(coarse_light-summary['first_intercept_bolometric_equivalent_W']['energy_J'])/coarse_light
            evaluations.append(result)
            out['datasets'][label]=dict(start_s=start,end_s=end,evaluations=evaluations,agreements=agreements)
            save();print(json.dumps(dict(stage='account',dataset=label,suns=suns,
                equivalent_m_s=cost['total_equivalent_m_s'],elapsed_s=time.monotonic()-before)),flush=True)
        if agreements['maximum_collection_relative_error']>1e-10 or agreements['maximum_envelope_excess_N_m_kg']>1e-10:
            raise ValueError('Load consistency failure')
    # Earlier failed variants remain alternative experiments, never extra inventory.
    for k,entry in enumerate(p['packing.json']['departures']):
        if entry['local_geometry_pass']:continue
        raw=load_raw(entry);spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        times=np.unique(np.r_[np.arange(raw['t'][0],raw['t'][-1],120.),raw['t'][-1]]);light={}
        trial_command=full_command(env,entry['schedule'])
        for t in times:
            values=uneclipsed_collection(env,t,np.c_[spline(t),spline(t,1)],trial_command(t).as_matrix(),suns=16)
            for name,value in values.items():
                if name!='reflected_force_N':light.setdefault(name,[]).append(value)
        light={name:np.array(value) for name,value in light.items()};summary={}
        for name,value in light.items():
            item=integrate_inventory(times,value);item.pop('per_member_energy_J');summary[name]=item
        path=RUN/f'rejected_{k}_collection.npz';np.savez_compressed(path,t=times,**light)
        out['rejected_optical_prefixes'].append(dict(trial=k,end_s=float(raw['t'][-1]),optical=summary,
            raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path)))
    chosen=[out['datasets'][k]['evaluations'][0] for k in ['service','departure']]
    impulse=sum(x['cost']['total_equivalent_m_s'] for x in chosen)
    energy=sum(x['cost']['propulsion_energy_J'] for x in chosen)
    peaks=np.maximum.reduce([np.array(x['cost']['per_member_peak_propulsion_W']) for x in chosen])
    prop=SCENARIO['propulsion'];mass=361*SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    old=p['energy_budget.json']['same_model_totals']['previous']
    out['totals']=dict(service_plus_departure_m_s=impulse,propulsion_energy_J=energy,
        same_model_fractional_energy_reduction=1-energy/old['service_plus_departure_energy_J'],
        baseline_service_plus_departure_m_s=old['service_plus_departure_m_s'],
        sum_member_peak_power_W=float(peaks.sum()),maximum_member_peak_power_W=float(peaks.max()),
        provisional_individual_power_hardware_kg=float(peaks.sum()*prop['peak_margin_factor']/prop['specific_power_W_kg']),
        provisional_power_hardware_fraction=float(peaks.sum()*prop['peak_margin_factor']/prop['specific_power_W_kg']/mass),
        total_burn_count=sum(x['cost']['scheduled_burn_count'] for x in chosen))
    reference,_=reference_shooting(env,21600.)
    config=validation['selected_config'];offsets=packing_offsets(config['depth_m'],config['swapped'])
    initial=load_raw(validation['runs'][0])['state'][0]
    frame_rate=(solar_frame(env.at(1.))-solar_frame(env.at(-1.)))/2
    preparation=length(initial[:,3:]-(reference.sol(0.)[3:6]+offsets@frame_rate.T))
    out['initial_preparation_reference']=dict(mean_velocity_difference_m_s=float(preparation.mean()),
        maximum_velocity_difference_m_s=float(preparation.max()),maneuver_executed=False,
        scope='Velocity difference from an initially Sun-tracking lattice translating with the reference. Positions are specified initial conditions; this is no acquisition trajectory or recurring impulse measurement.')
    out['conditional_inventory_trade']=[]
    for row in p['pilot_budget.json']['rows']:
        n=row['tiles'];limit=2*brentq(lambda dv:conditional_cost(n,dv,prop,hardware_multiplier=1.2,duty=.1)['propulsion_TW']-23.835,0.,2.)
        out['conditional_inventory_trade'].append(dict(tiles=n,
            assumed_two_day_full_cycle_allowance_m_s=limit,measured_partial_impulse_m_s=impulse,
            remaining_for_return_handover_recovery_m_s=limit-impulse,
            prefix_only_repeated_every_two_days=conditional_cost(n,impulse/2,prop,hardware_multiplier=1.2,duty=.1)))
    out['conditional_scope']='Hypothetical two-day recurrence, 20% additional fixed mass and 10% thrust duty; scalar power-hardware feedback only. The 17.29-million inventory is a capacity relaxation. No repeating cycle, acquisition, handover, return or fleet coverage has been demonstrated.'
    out['completed']=True
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    save();print(json.dumps(dict(stage='budget_complete',totals=out['totals'],elapsed_s=out['elapsed_s'])),flush=True)


if __name__=='__main__':main()
