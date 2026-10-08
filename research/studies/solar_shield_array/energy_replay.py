"""Coupled retry and independent replay of one measured-defect correction."""
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.departure_control import control_modes, mode_impulses
from protection.dynamics.parallel_shadow import ParallelPattern, minimum_clearance
from protection.dynamics.pattern_return import smooth_arc, beam_and_rate_screen, crossing_witnesses
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .departure_validate import audit_geometry
from .energy_screen import RUN
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def main():
    used=sum(json.loads((HERE/'results'/name).read_text())['elapsed_s']
             for name in ['energy_screen.json','energy.json','energy_correction.json'])
    allowance=min(1200,int(2700-used-900))
    if allowance<=0:raise ValueError('No remaining numerical budget with torque reserve')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();p=json.loads((HERE/'results/energy_correction.json').read_text())
    sources={**p['producer']['source_hashes'],
        'research/studies/solar_shield_array/energy_replay.py':digest(__file__)}
    for path,h in sources.items():
        if digest(ROOT/path)!=h:raise ValueError('Changed source: '+path)
    if constants_changed(p['producer']['constants']):raise ValueError('Changed constants')
    if digest(ROOT/p['raw_path'])!=p['raw_sha256']:raise ValueError('Changed raw')
    raw=np.load(ROOT/p['raw_path']);initial=raw['initial'];impulses=raw['impulses']
    start,end,burn=p['start_s'],p['end_s'],p['burn_s'];env=environment(1.)
    origin='measured_defect_correction'
    if not p['proposal_pass']:
        screen=json.loads((HERE/'results/energy_screen.json').read_text())
        failed=json.loads((HERE/'results/energy.json').read_text())
        tried={r['case'] for r in failed['runs']}
        index=next(i for i in screen['ranked_proposals'] if i not in tried)
        case=screen['cases'][index]
        p={**p,'case':index,'schedule':case['schedule'],'parameters_m_s':case['parameters_m_s']}
        impulses=mode_impulses(control_modes(env,start),np.array(case['parameters_m_s']))
        origin='next_screened_family_after_infeasible_measured_defect_correction'
    command=coast_command(env,start,end,**p['schedule']);windows=[[start,start+burn]]
    out=dict(schema='terluna.research.energy-coast-corrected-replay/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in ['energy_correction.json','energy_screen.json','energy.json']}),
        start_s=start,end_s=end,count=len(initial),schedule=p['schedule'],case=p['case'],
        candidate_origin=origin,parameters_m_s=p['parameters_m_s'],
        budget=dict(this_stage_wall_s=allowance,prior_numerical_wall_s=used,total_wall_s=2700,
            torque_stage_reserved_s=900,address_space_GiB=2),runs=[],
        accepted_departure=False,accepted_fleet=False,full_return_executed=False,
        handover_executed=False,hardware_mass_feedback_propagated=False,
        continuous_certificate=False,recurring_cost_measured=False,
        collection_capacity_W=None,delivered_power_W=None)
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/energy_replay.json').write_text(json.dumps(out,indent=2)+'\n')
    save();coarse=None
    for suns,step,rtol in [(8,120.,2e-10),(16,60.,2e-12)]:
        model=ParallelPattern(env,command,suns=suns);progress=[start+3600.]
        def rhs(t,flat):
            state=flat.reshape(-1,6)
            if t>=progress[0]:
                print(json.dumps(dict(stage='corrected_coupled',suns=suns,hours=(t-start)/3600,
                    elapsed_s=time.monotonic()-before)),flush=True);progress[0]+=3600.
            return np.c_[state[:,3:],model.acceleration(t,state)+smooth_arc(t,windows)[0]*impulses].ravel()
        def event(t,flat):return minimum_clearance(flat.reshape(-1,6),command(t).as_matrix())[0]-100.
        event.terminal=True;event.direction=-1
        sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=step,rtol=rtol,
            atol=np.tile([1e-5]*3+[1e-9]*3,len(initial)),events=event,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        ts=np.linspace(start,sol.t[-1],int(np.ceil((sol.t[-1]-start)/30))+1)
        ys=sol.sol(ts).T.reshape(-1,len(initial),6)
        at=np.linspace(start,sol.t[-1],int(np.ceil((sol.t[-1]-start)/2))+1)
        ay=sol.sol(at).T.reshape(-1,len(initial),6)
        geometry=audit_geometry(at,ay,command);beam=beam_and_rate_screen(env,at,ay,command)
        crossings=crossing_witnesses(ts,ys,command,lambda t,ids:sol.sol(t).reshape(-1,6)[ids],retain=8)
        maximum=max(float(length(rhs(t,y.ravel()).reshape(-1,6)[:,3:]).max()) for t,y in zip(ts[::4],ys[::4]))
        coverage=[]
        for t in np.arange(start,min(sol.t[-1],start+1200.)+1,300.):
            if t>start+p['schedule']['delay_s'] and p['schedule']['angle_deg']:break
            state=sol.sol(t).reshape(-1,6)
            entry=local_coverage(state,env.at(t),state.mean(axis=0),limb_count=16,interior=8)
            entry.update(time_s=float(t),window_inside_target=entry['window_centre_distance_m']+10000.<=4*K.MOON_RADIUS)
            coverage.append(entry)
        path=RUN/f'corrected_sun_{suns}.npz'
        np.savez_compressed(path,t=ts,state=ys,impulses=impulses,frames=command(ts).as_matrix())
        complete=bool(abs(sol.t[-1]-end)<1e-6)
        passed=bool(complete and not geometry['swept_separation']['possible_conflict_pairs']
            and not crossings['crossing_near_encounters'] and beam['minimum_sampled_beam_margin_deg']>0
            and beam['maximum_sampled_rate_deg_s']<=.1 and maximum<=.15
            and length(impulses).max()*2/burn<=.001
            and all(c['minimum_source_coverage']>=1-1e-10 for c in coverage if c['window_inside_target']))
        row=dict(suns=suns,max_step_s=step,rtol=rtol,end_s=float(sol.t[-1]),completed=complete,
            clearance_event=bool(len(sol.t_events[0])),geometry=geometry,beam_and_rate=beam,
            between_sample_crossings=crossings,maximum_sampled_acceleration_m_s2=maximum,
            maximum_thrust_acceleration_m_s2=float(length(impulses).max()*2/burn),
            translation_mean_m_s=float(length(impulses).mean()),extra_moving_window_samples=coverage,
            local_gates_pass=passed,raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
        interpolation=CubicHermiteSpline(ts,ys[:,:,:3],ys[:,:,3:])
        if coarse is None:coarse=interpolation
        else:
            out['independent_replay']=dict(maximum_position_difference_m=float(length(ay[:,:,:3]-coarse(at)).max()),
                maximum_velocity_difference_m_s=float(length(ay[:,:,3:]-coarse(at,1)).max()),
                reconstruction_position_error_m=float(length(ay[:,:,:3]-interpolation(at)).max()),comparison_step_s=2.)
            out['accepted_departure']=bool(passed and out['independent_replay']['maximum_position_difference_m']<50.)
        out['runs'].append(row);save();del ay
        print(json.dumps(dict(stage='corrected_gates',suns=suns,passed=passed,geometry=geometry)),flush=True)
        if not passed:break
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    save();print(json.dumps(dict(stage='complete',accepted=out['accepted_departure'],elapsed_s=out['elapsed_s'])),flush=True)


if __name__=='__main__':main()
