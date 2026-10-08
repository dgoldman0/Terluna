"""Executed impulse, mass allocation and optical accounts; no fictitious cycle."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline,interp1d
from scipy.spatial.transform import Rotation

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.tile_torque import parallel_load
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.cycle_control import return_command
from protection.dynamics.collection import integrate_inventory
from protection.dynamics.coverage_pilot import conditional_cost
from protection.dynamics.optical import length
from protection.dynamics.pattern_departure import sun_frame
from protection.dynamics.square_distance import closest_squares
from protection.dynamics.cycling import ray_geometry
from protection.dynamics.fleet import beam_margins
from .packing_screen import parents,load_raw
from .packing_run import full_command
from .closure_screen import RUN
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO


def cost_account(times,gravity,radiation,command,controls,windows,ratio):
    dense=np.unique(np.r_[np.arange(times[0],times[-1],2.),times[-1]])
    nominal=rigid_square_load(command,dense)['torque_per_mass']
    disturbance=interp1d(times,gravity+radiation/ratio[None,:,None],axis=0)(dense)
    couples=2*abs(nominal[:,None,:]-disturbance).sum(axis=-1)/10000.
    norms=length(controls);shapes=np.array([smooth_arc(t,windows) for t in dense])
    force=couples+shapes@norms.T
    def primitive(t):
        x=np.clip((t-windows[:,0])/np.diff(windows,axis=1).ravel(),0.,1.)
        return x-np.sin(2*np.pi*x)/(2*np.pi)
    translation=norms@(primitive(times[-1])-primitive(times[0]))
    attitude=np.trapezoid(couples,dense,axis=0);dv=translation+attitude
    prop=SCENARIO['propulsion'];mass=SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    specific=prop['exhaust_velocity_m_s']/(2*prop['efficiency']*np.cos(np.deg2rad(prop['cant_deg'])))
    energy=float(mass*specific*np.dot(ratio,dv));peak=force.max(axis=0)*mass*ratio*specific
    propellant=2*prop['efficiency']*energy/prop['exhaust_velocity_m_s']**2
    return dict(mean_applied_translation_m_s=float(translation.mean()),mean_attitude_equivalent_m_s=float(attitude.mean()),
        mass_weighted_equivalent_m_s=float(np.dot(ratio,dv)/ratio.sum()),
        optical_mass_equivalent_m_s=float(np.dot(ratio,dv)/len(ratio)),
        propulsion_energy_J=energy,per_member_energy_J=(mass*specific*ratio*dv).tolist(),
        ideal_exhaust_propellant_kg=propellant,fraction_of_carried_mass_expended=propellant/(mass*ratio.sum()),
        per_member_peak_W=peak.tolist(),sum_individual_peak_W=float(peak.sum()),
        started_burn_count=int(np.count_nonzero((norms>1e-8)&(windows[None,:,0]<times[-1]))))


def main():
    before=time.monotonic();names=['packing_budget.json','packing_validation.json','closure_screen.json','closure_hardware.json','closure_sizing.json','closure_return.json','closure_validation.json']
    p,sources=parents(names);used=sum(p[n]['elapsed_s'] for n in names[2:])
    allowance=min(600,int(4800-used))
    if allowance<=0:raise ValueError('No accounting budget remains')
    signal.alarm(allowance)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    sources['research/studies/solar_shield_array/closure_budget.py']=digest(__file__)
    out=dict(schema='terluna.research.cycle-closure-budget/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,this_stage_wall_s=allowance,threads=1,address_space_GiB=2),
        accepted_return=False,accepted_fleet=False,recurring_cost_measured=False,
        fleet_propulsion_W=None,collection_capacity_W=None,delivered_power_W=None,
        scope='Continuous hardware-loaded service, departure and rejected return prefix, all on the tighter replay. Count each member and executed interval; proposed full returns receive no operational cost or collection credit.',datasets={})
    env=environment(3.)
    command=full_command(env,p['packing_validation.json']['schedule'])
    cases=[('hardware_service',r,command) for r in p['closure_validation.json']['runs'][:1]]
    cases += [('hardware_departure',r,command) for r in p['closure_validation.json']['runs'][1:2]]
    r=p['closure_validation.json']['runs'][-1]
    retcommand=return_command(env,r['start_s'],r['intended_end_s'],p['closure_return.json']['arrival_turn_hours']*3600.)
    cases.append(('rejected_return',r,retcommand))
    # Preserve the distinct, preliminary departure-only hardware experiment.
    # Repeated numerical replays of the same case are not additional inventory.
    cases += [('prefix_power_service',r,command) for r in p['closure_hardware.json']['runs'][:1]]
    cases += [('prefix_power_departure',r,command) for r in p['closure_hardware.json']['runs'][1:]]
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
        (HERE/'results/closure_budget.json').write_text(json.dumps(out,indent=2)+'\n')
    for label,entry,cmd in cases:
        raw=load_raw(entry);spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        evaluations=[]
        settings=[(16,120.,.317)] if label.startswith('prefix_power_') else [(16,120.,.317),(32,240.,.317+np.pi/32)]
        for suns,spacing,rotation in settings:
            times=np.unique(np.r_[np.arange(raw['t'][0],raw['t'][-1],spacing),raw['t'][-1]])
            gravity=[];radiation=[];light={}
            for t in times:
                y=np.c_[spline(t),spline(t,1)];load=parallel_load(env,t,y,cmd(t).as_matrix(),suns=suns,rotation=rotation)
                gravity.append(load['gravity_torque_per_mass']);radiation.append(load['radiation_torque_per_mass'])
                for k,value in load['optical'].items():light.setdefault(k,[]).append(value)
            gravity=np.array(gravity);radiation=np.array(radiation);light={k:np.array(v) for k,v in light.items()}
            cost=cost_account(times,gravity,radiation,cmd,raw['controls'],raw['windows'],raw['mass_ratio'])
            optical={}
            for k,value in light.items():
                summary=integrate_inventory(times,value);summary.pop('per_member_energy_J');optical[k]=summary
            path=RUN/f'{label}_loads_{suns}.npz'
            np.savez_compressed(path,t=times,member_id=np.arange(361),mass_ratio=raw['mass_ratio'],
                gravity_torque_N_m_kg=gravity,radiation_torque_optical_N_m_kg=radiation,**light)
            evaluation=dict(sources=suns,time_step_s=spacing,cost=cost,optical=optical,
                raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
            if suns==16 and not label.startswith('prefix_power_'):
                indices=np.unique(np.r_[np.arange(0,len(times),2),len(times)-1])
                coarser=cost_account(times[indices],gravity[indices],radiation[indices],cmd,
                                    raw['controls'],raw['windows'],raw['mass_ratio'])
                coarse_light=integrate_inventory(times[indices],light['first_intercept_bolometric_equivalent_W'][indices])['energy_J']
                evaluation['time_coarsening_relative_energy_change']=abs(coarser['propulsion_energy_J']-cost['propulsion_energy_J'])/cost['propulsion_energy_J']
                evaluation['time_coarsening_relative_optical_energy_change']=abs(coarse_light-optical['first_intercept_bolometric_equivalent_W']['energy_J'])/coarse_light
            elif suns==32:
                evaluation['source_refinement_relative_energy_change']=abs(coarser['propulsion_energy_J']-cost['propulsion_energy_J'])/cost['propulsion_energy_J']
                evaluation['source_refinement_relative_optical_energy_change']=abs(coarse_light-optical['first_intercept_bolometric_equivalent_W']['energy_J'])/coarse_light
            evaluations.append(evaluation)
            out['datasets'][label]=dict(start_s=float(times[0]),end_s=float(times[-1]),evaluations=evaluations)
            save();print(json.dumps(dict(stage='account',dataset=label,suns=suns,energy_J=cost['propulsion_energy_J'],elapsed_s=time.monotonic()-before)),flush=True)
        if len(evaluations)>1:
            a,b=[e['cost']['propulsion_energy_J'] for e in evaluations]
            out['datasets'][label]['combined_time_source_relative_energy_change']=abs(a-b)/b
        else:out['datasets'][label]['alternative_mass_experiment']=True
    hardware=[out['datasets'][k]['evaluations'][0]['cost'] for k in ['hardware_service','hardware_departure'] if k in out['datasets']]
    all_costs=[e['cost'] for k in ['hardware_service','hardware_departure','rejected_return']
               for e in out['datasets'][k]['evaluations']]
    peaks=np.maximum.reduce([x['per_member_peak_W'] for x in all_costs])
    mass=SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    capacity=(ratio-1.2)*mass*SCENARIO['propulsion']['specific_power_W_kg']
    out['hardware_account']=dict(executed_energy_J=sum(x['propulsion_energy_J'] for x in hardware),
        complete_twelve_hour_passage=p['closure_validation.json']['accepted_hardware_passage'],
        previous_optical_only_energy_J=p['packing_budget.json']['totals']['propulsion_energy_J'],
        optical_mass_equivalent_m_s=sum(x['optical_mass_equivalent_m_s'] for x in hardware),
        total_mass_weighted_equivalent_m_s=sum(x['mass_weighted_equivalent_m_s'] for x in hardware),
        allocated_installed_power_W=float(capacity.sum()),sum_required_member_peak_W=float(peaks.sum()),
        minimum_individual_capacity_to_peak=float(np.min(capacity/peaks)),
        maximum_25_percent_margin_allocation_ratio=float(np.max(SCENARIO['propulsion']['peak_margin_factor']*peaks/capacity)),
        all_members_fit_installed_capacity=bool(np.all(peaks<=capacity)),
        full_margin_target_met=bool(np.all(SCENARIO['propulsion']['peak_margin_factor']*peaks<=capacity)),
        peak_sampling='Maximum over both 16-source/120-second and rotated 32-source/240-second load evaluations; nominal attitude and smooth burn power sampled every two seconds.',
        scope='Sampled individual peaks through the executed return prefix. Complete return, structural mass layout, deployment and repeated-cycle hardware still unclosed.')
    executed=[out['datasets'][k]['evaluations'][0]['cost']
              for k in ['hardware_service','hardware_departure','rejected_return']]
    out['executed_prefix_conditional_cost']=dict(assumed_cycle_days=2.,
        energy_through_stop_J=sum(x['propulsion_energy_J'] for x in executed),
        scope='Repeat only the measured expenditure through the failed-return stop once per two days, with the propagated mass allocation. No return completion, handover, recovery or global coverage is supplied by this arithmetic. Inventory remains hypothetical.',
        rows=[dict(tiles=row['tiles'],
            mean_propulsion_TW=sum(x['propulsion_energy_J'] for x in executed)*row['tiles']/361/(2*K.JULIAN_DAY)/1e12,
            installed_propulsion_capacity_TW=float(capacity.sum())*row['tiles']/361/1e12,
            propellant_kg_s=sum(x['ideal_exhaust_propellant_kg'] for x in executed)*row['tiles']/361/(2*K.JULIAN_DAY))
            for row in p['packing_budget.json']['conditional_inventory_trade']])
    final=load_raw(p['closure_validation.json']['runs'][-1]);t=float(final['t'][-1]);state=final['state'][-1]
    rays=ray_geometry(state[:,:3],env.at(t));orientations=[]
    for angle in [0.,15.,30.,45.,60.,90.]:
        frame=sun_frame(env,t)@Rotation.from_rotvec([-np.deg2rad(angle),0,0]).as_matrix()
        value,pair,projection=closest_squares(state[:,:3],frame)
        margin=beam_margins(state[:,:3],rays['sun'],frame[:,2],env.at(t),4*K.MOON_RADIUS,10000.)
        orientations.append(dict(angle_deg=angle,minimum_surface_distance_m=value,pair=pair.tolist(),
            minimum_beam_margin_deg=float(np.rad2deg(margin.min()))))
    out['frozen_event_attitudes']=dict(time_s=t,cases=orientations,
        scope='Same terminal centres and velocities with hypothetical common square orientations. No instantaneous slew, alternative force history, safe transition or operating saving is implied.')
    rolls=[];feathered=sun_frame(env,t)@Rotation.from_rotvec([-np.pi/2,0,0]).as_matrix()
    for angle in np.arange(0.,90.1,15.):
        frame=feathered@Rotation.from_rotvec([0,0,np.deg2rad(angle)]).as_matrix()
        value,pair,_=closest_squares(state[:,:3],frame)
        rolls.append(dict(roll_deg=float(angle),minimum_surface_distance_m=value,pair=pair.tolist()))
    out['frozen_event_attitudes']['same_normal_rolls']=rolls
    out['proposal_cost_scenarios']=[]
    for case in p['closure_screen.json']['cases']:
        if not case['shooting']['success']:continue
        proposed=case['nominal_prefix_plus_return_m_s'];rate=proposed/(case['end_s']/K.JULIAN_DAY)
        out['proposal_cost_scenarios'].append(dict(case=case['index'],nominal_proposed_equivalent_m_s=proposed,
            hypothetical_cycle_days=case['end_s']/K.JULIAN_DAY,includes_actual_return=False,
            excludes='Actual return disturbance loads, encounter repairs, handover, recovery and repeated-cycle changes.',
            rows=[dict(tiles=n,**conditional_cost(n,rate,SCENARIO['propulsion'],hardware_multiplier=1.2,duty=.1))
                  for n in [r['tiles'] for r in p['packing_budget.json']['conditional_inventory_trade']]]))
    out['completed']=True;save()


if __name__=='__main__':main()
