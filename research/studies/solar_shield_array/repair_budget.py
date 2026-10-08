"""Actual loads and collection for every executed correction, including stops."""
import argparse
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command
from protection.dynamics.return_escape import rolled_return
from protection.dynamics.collection import integrate_inventory
from protection.dynamics.compact_collection import compact_collection
from protection.dynamics.tile_torque import parallel_load
from .closure_budget import cost_account
from .packing_screen import load_raw
from .repair_screen import RUN,source_parents,write_product
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume-manifest',type=Path)
    args=parser.parse_args();before=time.monotonic()
    signal.alarm(660 if args.resume_manifest else 1200)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    trial_names=[]
    for path in sorted((HERE/'results').glob('repair_*.json')):
        j=json.loads(path.read_text())
        if j.get('schema')=='terluna.research.return-repair-coupled/1':trial_names.append(path.name)
    # Replays are evidence for the same physical experiment, not extra tiles.
    selected=[n for n in trial_names if n.endswith('_fine.json') or n.replace('.json','_fine.json') not in trial_names]
    names=selected+['closure_budget.json','closure_sizing.json','pilot_budget.json']
    p,sources=source_parents(names,['protection/dynamics/compact_collection.py',
        'research/studies/solar_shield_array/repair_budget.py'])
    cached={};cache_record=None
    if args.resume_manifest:
        manifest=json.loads(args.resume_manifest.read_text())
        for key in ['original_runner','partial_product']:
            if digest(ROOT/manifest[key]['path'])!=manifest[key]['sha256']:
                raise ValueError('Changed accounting continuation input')
        previous=json.loads((ROOT/manifest['partial_product']['path']).read_text())
        old_sources=previous['producer']['source_hashes']
        own='research/studies/solar_shield_array/repair_budget.py'
        if manifest['original_runner']['sha256']!=old_sources[own]:
            raise ValueError('Original accounting source identity missing')
        if {k:v for k,v in old_sources.items() if k!=own}!={k:v for k,v in sources.items() if k!=own}:
            raise ValueError('Accounting physics or parent sources changed')
        if previous['producer']['inputs']!={n:digest(HERE/'results'/n) for n in names}:
            raise ValueError('Accounting parent results changed')
        for trial in previous['trials']:
            for e in trial['evaluations']+trial['grazing_refinement']:
                if manifest['files'].get(e['raw_path'])!=e['raw_sha256'] or digest(ROOT/e['raw_path'])!=e['raw_sha256']:
                    raise ValueError('Completed accounting data changed')
            cached[trial['source']]=trial
        cache_record=dict(path=str(args.resume_manifest.relative_to(ROOT)),sha256=digest(args.resume_manifest),
            original_runner=manifest['original_runner'],partial_product=manifest['partial_product'],
            prior_attempt_wall_s=manifest['charged_first_attempt_wall_s'],completed_trials=list(cached))
    env=environment(3.);out=dict(schema='terluna.research.return-repair-budget/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        full_cycle_energy_J=None,full_cycle_average_W=None,delivered_electricity_W=None,
        scope='Each alternative retains the unchanged twelve-hour prefix once, followed by its own executed return through the stopping point. Tighter replay replaces the coarse trace in cost totals. No unexecuted remainder, recovery, handover or recurrence is credited.',trials=[])
    out['accounting_continuation']=cache_record
    ratio=np.array(p['closure_sizing.json']['per_member_total_to_optical_mass_ratio'])
    capacity=(ratio-1.2)*.05*10000.**2*SCENARIO['propulsion']['specific_power_W_kg']
    prefix=[p['closure_budget.json']['datasets'][k]['evaluations'][-1] for k in ['hardware_service','hardware_departure']]
    for name in selected:
        if name in cached:
            out['trials'].append(cached[name])
            write_product('repair_budget.json',out,sources,before)
            print(json.dumps(dict(stage='reuse_completed_account',trial=name)),flush=True)
            continue
        parent=p[name];entry=parent['run'];raw=load_raw(entry)
        spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        start,end=map(float,raw['t'][[0,-1]]);spec=parent['schedule']
        base=return_command(env,start,parent['arrival_s'],parent['arrival_turn_s'])
        command=rolled_return(base,start,parent['arrival_s'],spec.get('angle_deg',0.),spec.get('turn_s',5400.))
        trial=dict(source=name,start_s=start,end_s=end,evaluations=[],grazing_refinement=[])
        for suns,spacing,rotation in [(16,120.,.317),(32,240.,.317+np.pi/32)]:
            t=np.unique(np.r_[np.arange(start,end,spacing),end]);gravity=[];radiation=[];light={}
            for ti in t:
                y=np.c_[spline(ti),spline(ti,1)]
                data=parallel_load(env,ti,y,command(ti).as_matrix(),suns=suns,rotation=rotation)
                gravity.append(data['gravity_torque_per_mass']);radiation.append(data['radiation_torque_per_mass'])
                for key,value in data['optical'].items():light.setdefault(key,[]).append(value)
            gravity=np.array(gravity);radiation=np.array(radiation);light={k:np.array(v) for k,v in light.items()}
            cost=cost_account(t,gravity,radiation,command,raw['controls'],raw['windows'],ratio)
            optical={k:integrate_inventory(t,v) for k,v in light.items()}
            path=RUN/(name[:-5]+f'_loads_{suns}.npz')
            np.savez_compressed(path,t=t,member_id=np.arange(361),mass_ratio=ratio,
                gravity_torque_N_m_kg=gravity,radiation_torque_optical_N_m_kg=radiation,**light)
            ev=dict(sources=suns,time_step_s=spacing,cost=cost,optical=optical,
                raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
            if suns==16:
                idx=np.unique(np.r_[np.arange(0,len(t),2),len(t)-1])
                coarse=cost_account(t[idx],gravity[idx],radiation[idx],command,raw['controls'],raw['windows'],ratio)
                coarse_opt=integrate_inventory(t[idx],light['first_intercept_bolometric_equivalent_W'][idx])
                ev['time_coarsening_relative_energy_change']=abs(coarse['propulsion_energy_J']-cost['propulsion_energy_J'])/cost['propulsion_energy_J']
            else:
                ev['source_refinement_relative_energy_change']=abs(coarse['propulsion_energy_J']-cost['propulsion_energy_J'])/cost['propulsion_energy_J']
                ev['matched_time_source_relative_optical_change']=abs(coarse_opt['energy_J']-optical['first_intercept_bolometric_equivalent_W']['energy_J'])/optical['first_intercept_bolometric_equivalent_W']['energy_J']
            trial['evaluations'].append(ev)
        # More source directions for grazing interception, retaining per-member data.
        for suns in [64,128]:
            light={}
            for ti in t:
                d=compact_collection(env,ti,np.c_[spline(ti),spline(ti,1)],command(ti).as_matrix(),suns=suns,rotation=.317+np.pi/suns)
                for key,value in d.items():
                    if key!='reflected_force_N':light.setdefault(key,[]).append(value)
            light={k:np.array(v) for k,v in light.items()}
            optical={k:integrate_inventory(t,v) for k,v in light.items()}
            path=RUN/(name[:-5]+f'_optical_{suns}.npz')
            np.savez_compressed(path,t=t,member_id=np.arange(361),**light)
            trial['grazing_refinement'].append(dict(sources=suns,time_step_s=240.,optical=optical,
                raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path)))
        low,high=trial['grazing_refinement']
        trial['grazing_64_128_relative_optical_change']=abs(low['optical']['first_intercept_bolometric_equivalent_W']['energy_J']-
            high['optical']['first_intercept_bolometric_equivalent_W']['energy_J'])/high['optical']['first_intercept_bolometric_equivalent_W']['energy_J']
        costs=[e['cost'] for e in trial['evaluations']]
        peaks=np.maximum.reduce([np.array(e['cost']['per_member_peak_W']) for e in prefix]+[np.array(c['per_member_peak_W']) for c in costs])
        trial['hardware']=dict(all_member_peaks_fit=bool(np.all(peaks<=capacity)),
            all_members_keep_25_percent_margin=bool(np.all(1.25*peaks<=capacity)),
            members_over_capacity=int(np.count_nonzero(peaks>capacity)),
            minimum_installed_to_peak=float(np.min(capacity/peaks)),sum_member_peaks_W=float(peaks.sum()),
            required_extra_power_hardware_kg=float(np.maximum(1.25*peaks-capacity,0.).sum()/SCENARIO['propulsion']['specific_power_W_kg']),
            mass_feedback='Existing allocation propagated. Any indicated extra hardware has not been added; no hardware acceptance if its margin fails.')
        actual=costs[0];energy=sum(e['cost']['propulsion_energy_J'] for e in prefix)+actual['propulsion_energy_J']
        fuel=sum(e['cost']['ideal_exhaust_propellant_kg'] for e in prefix)+actual['ideal_exhaust_propellant_kg']
        optical_prefix=[e['optical'] for e in prefix];optical_return=high['optical']
        trial['whole_executed_prefix']=dict(duration_s=end,propulsion_energy_J=energy,mean_propulsion_W=energy/end,
            ideal_exhaust_propellant_kg=fuel,
            first_intercept_energy_J=sum(e['first_intercept_bolometric_equivalent_W']['energy_J'] for e in optical_prefix)+optical_return['first_intercept_bolometric_equivalent_W']['energy_J'],
            redirected_band_energy_J=sum(e['redirected_band_optical_W']['energy_J'] for e in optical_prefix)+optical_return['redirected_band_optical_W']['energy_J'])
        for key in ['first_intercept','redirected_band']:
            trial['whole_executed_prefix'][key+'_mean_W']=trial['whole_executed_prefix'][key+'_energy_J']/end
        trial['conditional_two_day_repetition']=[dict(tiles=r['tiles'],mean_propulsion_W=energy*r['tiles']/361/(2*K.JULIAN_DAY),
            exhaust_kg_s=fuel*r['tiles']/361/(2*K.JULIAN_DAY)) for r in p['pilot_budget.json']['rows']]
        out['trials'].append(trial);write_product('repair_budget.json',out,sources,before)
        print(json.dumps(dict(trial=name,return_energy_J=actual['propulsion_energy_J'],whole_energy_J=energy,
            hardware=trial['hardware'],grazing_change=trial['grazing_64_128_relative_optical_change'],elapsed_s=time.monotonic()-before)),flush=True)
    out['completed']=True;write_product('repair_budget.json',out,sources,before)


if __name__=='__main__':main()
