"""Torque and collection on accepted historical and rejected new prefixes."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline, interp1d
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation, RotationSpline

from shared.provenance import constants_used, constants_changed
from protection.dynamics.tile_torque import parallel_load, gravity_torque
from protection.dynamics.collection import uneclipsed_collection, integrate_inventory
from protection.dynamics.attitude_load import rigid_square_load, disturbance_torque_bound
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.pattern_departure import departure_command, sun_frame
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.coverage_pilot import conditional_cost
from protection.dynamics.optical import length
from .energy_screen import RUN
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO


def actuator_integral(times, gravity, radiation, command, translation, burn_s=7200.):
    """Two-second command, linearly interpolated disturbance, explicit couples."""
    dense = np.unique(np.r_[np.arange(times[0], times[-1], 2.), times[-1]])
    nominal = rigid_square_load(command, dense)
    disturbance = interp1d(times, gravity+radiation, axis=0, assume_sorted=True)(dense)
    required = nominal['torque_per_mass'][:, None, :]-disturbance
    force = 2*np.abs(required).sum(axis=-1)/10000.
    attitude = float(np.trapezoid(force.mean(axis=1), dense))
    nominal_impulse = float(np.trapezoid(nominal['force_per_mass'], dense))
    shape = np.array([smooth_arc(t, [[times[0], times[0]+burn_s]])[0] for t in dense])
    mean_force = force.mean(axis=1)+translation*shape
    fraction = min((times[-1]-times[0])/burn_s, 1.)
    executed_translation = translation*(fraction-np.sin(2*np.pi*fraction)/(2*np.pi))
    mass = force.shape[1]*SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    prop = SCENARIO['propulsion']
    conversion = prop['exhaust_velocity_m_s']/(2*prop['efficiency']*np.cos(np.deg2rad(prop['cant_deg'])))
    return dict(translation_m_s=executed_translation, attitude_equivalent_m_s=attitude,
        nominal_attitude_equivalent_m_s=nominal_impulse,
        total_equivalent_m_s=executed_translation+attitude,
        propulsion_energy_J=mass*(executed_translation+attitude)*conversion,
        sampled_peak_propulsion_W=float(mean_force.max()*mass*conversion),
        mean_gravity_only_couple_equivalent_m_s=float(np.trapezoid(
            (2*np.abs(gravity).sum(axis=-1)/10000.).mean(axis=1), times)),
        mean_radiation_only_couple_equivalent_m_s=float(np.trapezoid(
            (2*np.abs(radiation).sum(axis=-1)/10000.).mean(axis=1), times)))


def main():
    signal.alarm(900); resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic(); env = environment(1.)
    names = ['patterns.json', 'departure_validation.json', 'departure_refine.json',
             'departure_budget.json', 'energy_screen.json', 'energy.json',
             'energy_correction.json', 'energy_replay.json', 'pilot_budget.json']
    parents = {n: json.loads((HERE/'results'/n).read_text()) for n in names}
    accepted = parents['energy_replay.json']['accepted_departure']
    sources = {}
    for p in parents.values():
        for path, expected in p['producer']['source_hashes'].items():
            if digest(ROOT/path) != expected: raise ValueError('Changed source: '+path)
        if constants_changed(p['producer']['constants']): raise ValueError('Changed constants')
        sources.update(p['producer']['source_hashes'])
    for path in ['protection/dynamics/tile_torque.py','protection/dynamics/collection.py',
                 'research/studies/solar_shield_array/energy_budget.py']:
        sources[path] = digest(ROOT/path)
    selected = parents['energy_replay.json']
    entries = dict(service=parents['patterns.json']['replay'],
        previous_departure=parents['departure_validation.json']['runs'][-1],
        new_departure=parents['energy_replay.json']['runs'][-1])
    out = dict(schema='terluna.research.energy-coast-budget/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}, raw_inputs={}),
        budget=dict(wall_s=900, address_space_GiB=2, numerical_threads=1),
        torque_scope='Uniform rigid 10 km square, finite-area gravity and first spatial moments of reflected-band shadow unions. Euler torque is supplied by ideal opposite edge-midpoint electric couples. No absorption, thermal, flexible, plume or hardware realization model.',
        integration_scope='Actual command at 2 s; linear interpolation of calculated disturbance torques sampled at 60 s. Independent solar quadrature uses 32 rotated sources at 120 s. No new mass or force repropagation.',
        accepted_new_departure=accepted, accepted_fleet=False, full_return_executed=False, handover_executed=False,
        new_candidate_scope=('Complete ideal-model local departure' if accepted else
            'Rejected partial prefix ending at the clearance stop; no favorable full-interval energy comparison or recurring-power extrapolation is permitted.'),
        measured_recurring_impulse_m_s_day=None, fleet_propulsion_W=None,
        collection_capacity_W=None, delivered_power_W=None,
        hardware_mass_feedback_propagated=False, structural_realizability_demonstrated=False,
        datasets={})
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/energy_budget.json').write_text(json.dumps(out, indent=2)+'\n')
    for label, entry in entries.items():
        if digest(ROOT/entry['raw_path']) != entry['raw_sha256']: raise ValueError('Changed raw')
        out['producer']['raw_inputs'][entry['raw_path']] = entry['raw_sha256']
        raw = np.load(ROOT/entry['raw_path']); start, end = map(float, raw['t'][[0,-1]])
        interpolation = CubicHermiteSpline(raw['t'], raw['state'][:, :, :3], raw['state'][:, :, 3:])
        if label == 'service':
            command = RotationSpline(raw['t'], Rotation.from_matrix([sun_frame(env,t) for t in raw['t']]))
            translation = 0.
        elif label == 'previous_departure':
            r = parents['departure_refine.json']
            command = departure_command(env,start,end,r['delay_s'],r['slew_s'],r['axis'],r['sign'])
            translation = parents['departure_budget.json']['translation_mean_delta_v_m_s']
        else:
            # Rejected traces end before the commanded turn finishes. Preserve
            # its original duration and evaluate only the completed prefix.
            command = coast_command(env,start,parents['energy_replay.json']['end_s'],**selected['schedule'])
            translation = entry['translation_mean_m_s']
        evaluated = []
        agreements = dict(maximum_collection_relative_error=0., maximum_gravity_quadrature_torque_error_N_m_kg=0.,
            maximum_bound_excess_N_m_kg=0.)
        for suns, step, rotation in [(16,60.,.317),(32,120.,.317+np.pi/32)]:
            times = np.unique(np.r_[np.arange(start,end,step),end]); gravity=[]; radiation=[]; optical={}
            for k,t in enumerate(times):
                state = np.c_[interpolation(t),interpolation(t,1)]; frame=command(t).as_matrix()
                load=parallel_load(env,t,state,frame,suns=suns,rotation=rotation)
                gravity.append(load['gravity_torque_per_mass']); radiation.append(load['radiation_torque_per_mass'])
                for name,value in load['optical'].items():optical.setdefault(name,[]).append(value)
                if suns==16 and k%60==0:
                    check=uneclipsed_collection(env,t,state,frame,suns=16)
                    scale=max(float(check['first_intercept_bolometric_equivalent_W'].sum()),1.)
                    error=abs(float((check['first_intercept_bolometric_equivalent_W']-load['optical']['first_intercept_bolometric_equivalent_W']).sum()))/scale
                    agreements['maximum_collection_relative_error']=max(agreements['maximum_collection_relative_error'],error)
                    g5=gravity_torque(env,t,state[:,:3],frame,order=5)
                    agreements['maximum_gravity_quadrature_torque_error_N_m_kg']=max(
                        agreements['maximum_gravity_quadrature_torque_error_N_m_kg'],float(length(g5-load['gravity_torque_per_mass']).max()))
                    gb,rb,_=disturbance_torque_bound(env,t,state,frame)
                    agreements['maximum_bound_excess_N_m_kg']=max(agreements['maximum_bound_excess_N_m_kg'],
                        float((length(load['gravity_torque_per_mass'])-gb).max()),float((length(load['radiation_torque_per_mass'])-rb).max()))
                if k%180==0:
                    print(json.dumps(dict(stage='torque_collection',dataset=label,suns=suns,
                        hours=(t-start)/3600.,elapsed_s=time.monotonic()-before)),flush=True)
            gravity=np.array(gravity);radiation=np.array(radiation)
            optical={k:np.array(v) for k,v in optical.items()}
            cost=actuator_integral(times,gravity,radiation,command,translation)
            path=RUN/f'{label}_loads_sun_{suns}.npz'
            np.savez_compressed(path,t=times,gravity_torque_N_m_kg=gravity,radiation_torque_N_m_kg=radiation,**optical)
            light={}
            for name,value in optical.items():
                result=integrate_inventory(times,value);result.pop('per_member_energy_J');light[name]=result
            item=dict(sources=suns,step_s=step,cost=cost,optical=light,raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
            if suns==16:
                coarse_indices=np.unique(np.r_[np.arange(0,len(times),2),len(times)-1])
                coarse=actuator_integral(times[coarse_indices],gravity[coarse_indices],radiation[coarse_indices],command,translation)
                item['time_coarsening_relative_attitude_change']=abs(coarse['attitude_equivalent_m_s']-cost['attitude_equivalent_m_s'])/max(cost['attitude_equivalent_m_s'],1e-12)
                cost16_coarse=coarse
            else:
                item['source_refinement_relative_attitude_change']=abs(cost['attitude_equivalent_m_s']-cost16_coarse['attitude_equivalent_m_s'])/max(cost['attitude_equivalent_m_s'],1e-12)
            evaluated.append(item);out['datasets'][label]=dict(evaluations=evaluated,agreements=agreements,
                start_s=start,end_s=end,completed_intended_interval=bool(label!='new_departure' or accepted));save()
        if agreements['maximum_collection_relative_error']>1e-10 or agreements['maximum_bound_excess_N_m_kg']>1e-10:
            raise ValueError('Torque/collection consistency or envelope failure')
    chosen={k:v['evaluations'][0] for k,v in out['datasets'].items()}
    service=chosen['service']['cost']; old=chosen['previous_departure']['cost']; new=chosen['new_departure']['cost']
    totals={}
    prop=SCENARIO['propulsion']
    for label,cost in [('previous',old),('new' if accepted else 'rejected_partial_prefix',new)]:
        peak=max(service['sampled_peak_propulsion_W'],cost['sampled_peak_propulsion_W'])
        totals[label]=dict(service_plus_departure_m_s=service['total_equivalent_m_s']+cost['total_equivalent_m_s'],
            service_plus_departure_energy_J=service['propulsion_energy_J']+cost['propulsion_energy_J'],
            peak_propulsion_W=peak, provisional_power_hardware_kg=prop['peak_margin_factor']*peak/prop['specific_power_W_kg'])
    out['same_model_totals']=totals
    out['same_model_fractional_energy_reduction']=(1-totals['new']['service_plus_departure_energy_J']/totals['previous']['service_plus_departure_energy_J'] if accepted else None)
    out['previous_envelope_total_m_s']=parents['departure_budget.json']['service_plus_departure_equivalent_delta_v_m_s']
    out['conditional_two_day_comparisons']=[]
    for row in (parents['pilot_budget.json']['rows'] if accepted else []):
        n=row['tiles']
        limit=brentq(lambda dv:conditional_cost(n,dv,prop,hardware_multiplier=1.2,duty=.1)['propulsion_TW']-23.835,0.,2.)
        dv=totals['new']['service_plus_departure_m_s']
        out['conditional_two_day_comparisons'].append(dict(tiles=n,
            complete_two_day_impulse_allowance_at_10_percent_m_s=2*limit,
            partial_prefix_m_s=dv, remaining_impulse_for_every_other_stage_m_s=2*limit-dv,
            partial_prefix_repeated_once_per_two_days=conditional_cost(n,dv/2,prop,hardware_multiplier=1.2,duty=.1)))
    out['conditional_scope']='Assumed two-day repetition, 20% extra fixed mass and 10% thrust duty. Only the measured partial prefix is charged; return, acquisition, handover and recovery remain unpaid. The 17.29-million inventory is a capacity relaxation, not a sufficient fleet. Scalar mass feedback has not been propagated.'
    out['completed']=True
    if any(digest(ROOT/p)!=h for p,h in sources.items()):raise ValueError('Source changed')
    save();print(json.dumps(dict(stage='complete',totals=totals,reduction=out['same_model_fractional_energy_reduction'],
        elapsed_s=out['elapsed_s'],max_rss_MiB=out['max_rss_MiB'])),flush=True)


if __name__=='__main__':main()
