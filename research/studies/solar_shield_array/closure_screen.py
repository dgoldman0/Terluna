"""Bounded return/arrival proposals from the deeper pattern's actual states."""
import json
import resource
import signal
import time

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.cycle_control import (return_command, independent_path,
    cycle_windows, corrected_state, endpoint_controls)
from protection.dynamics.packing_control import packing_offsets, response_maps, matrices
from protection.dynamics.pattern_return import encounter_scan, crossing_witnesses, beam_and_rate_screen
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.attitude_load import rigid_square_load
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.optical import length
from .packing_screen import parents, load_raw
from .return_screen import next_passage
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE

RUN=ROOT/'research/runs/solar_shield_array/closure'


def sources_for(names):
    p,sources=parents(names)
    for path in ['protection/dynamics/cycle_control.py',
                 'research/studies/solar_shield_array/closure_screen.py']:
        sources[path]=digest(ROOT/path)
    return p,sources


def main():
    signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();RUN.mkdir(parents=True,exist_ok=True)
    names=['packing_validation.json','packing_budget.json','packing_audit.json']
    p,sources=sources_for(names);parent=p['packing_validation.json']
    raw=load_raw(parent['runs'][1]);initial=raw['state'][-1];start=float(raw['t'][-1])
    first=load_raw(parent['runs'][0])['state'][0];env=environment(3.)
    coast=return_command(env,start,2.6*K.JULIAN_DAY,14400.)
    reference=independent_path(env,initial.mean(axis=0)[None,:],start,2.6*K.JULIAN_DAY,coast)
    passage=next_passage(env,reference,first.mean(axis=0),start,reference.t[-1])
    out=dict(schema='terluna.research.cycle-closure-screen/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names},
            raw_inputs={parent['runs'][1]['raw_path']:parent['runs'][1]['raw_sha256']}),
        budget=dict(total_wall_s=4800,this_stage_wall_s=600,threads=1,address_space_GiB=2,
                    maximum_cases=12,maximum_refined_candidates=3,maximum_corrections_each=4),
        start_s=start,next_feathered_reference_passage_s=passage,
        prior_executed_impulse_m_s=p['packing_budget.json']['totals']['service_plus_departure_m_s'],
        accepted_return=False,accepted_fleet=False,hardware_mass_feedback_propagated=False,
        scope='Independent-member finite-source proposals and common gravity response maps; mutual shadows omitted until coupled propagation. Exact endpoint packing is one tested service assignment. No recurrence or fleet acceptance.',
        cases=[],refinements=[])
    stored={}
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
        (HERE/'results/closure_screen.json').write_text(json.dumps(out,indent=2)+'\n')
    for shift in [-1800.,0.,1800.]:
        end=passage+shift;windows=cycle_windows(start,end)
        for hours in [4.,8.]:
            command=return_command(env,start,end,hours*3600.)
            free=independent_path(env,initial,start,end,command)
            base=lambda t,f=free:f.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
            maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),start,end,windows,np.eye(3))
            centre=base([end])[0].mean(axis=0)
            # The pure Sun frame continues after the arrival command.
            service_command=coast_command(env,end,end+21600.,0.,0.,3600.)
            service=independent_path(env,centre[None,:],end,end+21600.,service_command)
            service_maps=response_maps(env,lambda t:service.sol(t),end,end+21600.,np.empty((0,2)),np.eye(3))
            phi=matrices(service_maps,[end+21600.],0)[0]
            frame=solar_frame(env.at(end))
            gradient=.95*np.linalg.solve(phi[:3,3:6],solar_frame(env.at(end+21600.))-phi[:3,:3]@frame)
            service_dates=np.linspace(end,end+21600.,13)
            inside=all(float(length(service.sol(t)[:3]-(service.sol(t)[:3]@solar_frame(env.at(t))[:,0])*solar_frame(env.at(t))[:,0]))+10000.<=4*K.MOON_RADIUS for t in service_dates)
            # Repeat the same 36 km depth packing, allowing identity reassignment.
            offsets=packing_offsets(12000.)
            times=np.linspace(start,end,int(np.ceil((end-start)/600.))+1)
            dt=np.linspace(start,end,int(np.ceil((end-start)/20.))+1)
            nominal=float(np.trapezoid(rigid_square_load(command,dt)['force_per_mass'],dt))
            for quarter in [0,2]:
                permutation=np.rot90(np.arange(361).reshape(19,19),quarter).ravel()
                d=offsets[permutation]
                target=centre+np.c_[d@frame.T,d@gradient.T]
                u,fit=endpoint_controls(base,maps,target,windows)
                row=dict(index=len(out['cases']),end_s=end,date_shift_s=shift,arrival_turn_hours=hours,
                    quarter_turn_assignment=quarter,windows_s=windows.tolist(),shooting=fit,
                    next_reference_six_hour_window_inside_target=inside,nominal_arrival_attitude_m_s=nominal)
                if u is not None:
                    state=corrected_state(base,maps,u,times)
                    row['encounters']=encounter_scan(times,state,command,clearance=350.,retain=8)
                    row['beam_and_rate']=beam_and_rate_screen(env,times,state,command)
                    row['entry_coverage']=local_coverage(state[-1],env.at(end),centre,limb_count=16,interior=8,rotation=.137)
                    row['nominal_prefix_plus_return_m_s']=out['prior_executed_impulse_m_s']+nominal+fit['mean_delta_v_m_s']
                    row['proposal_gate']=bool(inside and not row['encounters']['sampled_pair_violations']
                        and row['beam_and_rate']['minimum_sampled_beam_margin_deg']>0
                        and row['beam_and_rate']['maximum_sampled_rate_deg_s']<=.1)
                    stored[row['index']]=(free,maps,target,u,command,windows)
                out['cases'].append(row);save()
                print(json.dumps(dict(stage='return_case',index=row['index'],hours=hours,shift_s=shift,
                    assignment=quarter,shooting=fit,conflicts=row.get('encounters',{}).get('sampled_pair_violations'),
                    inside=inside,nominal_total_m_s=row.get('nominal_prefix_plus_return_m_s'),elapsed_s=time.monotonic()-before)),flush=True)
    def rank(index):
        c=out['cases'][index];b=c['beam_and_rate']
        return (not c['next_reference_six_hour_window_inside_target'],b['minimum_sampled_beam_margin_deg']<=0,
                c['nominal_prefix_plus_return_m_s'])
    ranked=sorted(stored,key=rank)
    for index in ranked[:3]:
        free,maps,target,u,command,windows=stored[index];end=out['cases'][index]['end_s']
        base=lambda t,f=free:f.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
        times=np.linspace(start,end,int(np.ceil((end-start)/120.))+1)
        cuts=[];seen=set();accepted=False;record=dict(case=index,steps=[])
        for attempt in range(5):
            state=corrected_state(base,maps,u,times)
            static=encounter_scan(times,state,command,clearance=350.,retain=128)
            roots=crossing_witnesses(times,state,command,lambda t,ids:corrected_state(base,maps,u,[t])[0,ids],retain=128)
            step=dict(attempt=attempt,encounters=static,crossings=roots);record['steps'].append(step)
            if not static['sampled_pair_violations'] and not roots['crossing_near_encounters']:
                accepted=True;break
            if attempt==4:break
            for w in roots['witnesses']+static['witnesses']:
                key=(round(w['time_s'],2),w['i'],w['j'],tuple(np.round(w['axis'],5)))
                if key in seen:continue
                seen.add(key);cuts.append((w['time_s'],w['i'],w['j'],np.array(w['axis']),max(w['required_projection_m'],350.)))
            changed,fit=endpoint_controls(base,maps,target,windows,cuts)
            step['correction']=fit
            if changed is None:break
            u=changed
        state=corrected_state(base,maps,u,times)
        path=RUN/f'return_proposal_{index}.npz'
        np.savez_compressed(path,t=times,state=state,initial=initial,controls=u,windows=windows,target=target,
                            free_state=base(times),maps=maps.sol(times).T,frames=command(times).as_matrix())
        record.update(proposal_geometry_pass=accepted,raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path),
            final_mean_translation_m_s=float(length(u).sum(axis=1).mean()),
            beam_and_rate=beam_and_rate_screen(env,times,state,command))
        out['refinements'].append(record);save()
        print(json.dumps(dict(stage='return_refined',case=index,passed=accepted,
            translation_m_s=record['final_mean_translation_m_s'],elapsed_s=time.monotonic()-before)),flush=True)
        if accepted and out['cases'][index]['next_reference_six_hour_window_inside_target'] and record['beam_and_rate']['minimum_sampled_beam_margin_deg']>0:
            out['selected_case']=index;break
    if 'selected_case' not in out:
        out['diagnostic_case']=out['refinements'][0]['case'] if out['refinements'] else None
    out['completed']=True;save()


if __name__=='__main__':main()
