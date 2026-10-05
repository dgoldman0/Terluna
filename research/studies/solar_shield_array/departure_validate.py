"""Refine the conservative interval guard and independently replay departure."""
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_formation import swept_square_conflicts
from protection.dynamics.pattern_departure import departure_command
from protection.dynamics.parallel_shadow import ParallelPattern, minimum_clearance
from protection.dynamics.pattern_return import smooth_arc, beam_and_rate_screen, crossing_witnesses
from protection.dynamics.natural_pattern import local_coverage
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE
from .departure_screen import RUN


def audit_geometry(times, states, command):
    static = [minimum_clearance(y, command(t).as_matrix()) for t, y in zip(times, states)]
    k = int(np.argmin([x[0] for x in static])); minimum, pair, projections = static[k]
    swept = swept_square_conflicts(times, states, lambda t: command(t).as_matrix()[:, [2, 0, 1]],
        acceleration_bound=.15, normal_rate_bound=np.deg2rad(.1))
    return dict(step_s=float(np.diff(times).max()), minimum_sampled_clearance_m=minimum,
        limiting_pair=pair.tolist(), worst_time_s=float(times[k]),
        pair_axis_projections_m=projections.tolist(), swept_separation=swept)


def main():
    signal.alarm(1800)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic()
    p = json.loads((HERE/'results/departure.json').read_text())
    r = json.loads((HERE/'results/departure_refine.json').read_text())
    for parent in [p, r]:
        for path, h in parent['producer']['source_hashes'].items():
            if digest(ROOT/path) != h: raise ValueError('Source changed: '+path)
        if constants_changed(parent['producer']['constants']): raise ValueError('Constants changed')
    coarse = p['runs'][0]
    if digest(ROOT/coarse['raw_path']) != coarse['raw_sha256']: raise ValueError('Raw changed')
    raw = np.load(ROOT/coarse['raw_path']); initial=raw['state'][0]; impulses=raw['impulses']
    if not coarse['completed_prefix']: raise ValueError('Only a completed prefix can refine its interval guard')
    start, end = p['start_s'], p['intended_end_s']; env = environment(1.)
    command = departure_command(env, start, end, r['delay_s'], r['slew_s'], r['axis'], r['sign'])
    times = np.arange(start, end+1, 2.)
    interpolation = CubicHermiteSpline(raw['t'], raw['state'][:, :, :3], raw['state'][:, :, 3:])
    state = np.concatenate([interpolation(times), interpolation(times, 1)], axis=-1)
    sources = {**p['producer']['source_hashes'],
        'research/studies/solar_shield_array/departure_validate.py': digest(__file__)}
    out = dict(schema='terluna.research.departure-validation/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={n: digest(HERE/'results'/n) for n in ['departure.json', 'departure_refine.json']}),
        start_s=start, intended_end_s=end, count=len(initial), accepted_departure=False,
        accepted_fleet=False, accepted_return=False, accepted_handover=False, continuous_certificate=False,
        coarse_two_second_audit=audit_geometry(times, state, command),
        coarse_interpolation_scope='Cubic Hermite reconstruction of the stored 30 s positions and velocities; independent dense-output replay below checks the reconstruction and dynamics together.',
        runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/departure_validation.json').write_text(json.dumps(out, indent=2)+'\n')
    save(); del state
    print(json.dumps(dict(stage='refined_coarse_guard', result=out['coarse_two_second_audit'])), flush=True)
    model = ParallelPattern(env, command, suns=16); windows=np.array([[start,start+r['burn_s']]])
    next_progress=[start+1800.]
    def rhs(t, flat):
        state=flat.reshape(-1,6)
        if t >= next_progress[0]:
            print(json.dumps(dict(stage='replay', hours_after_service=(t-start)/3600,
                elapsed_s=time.monotonic()-before)), flush=True); next_progress[0]+=1800.
        return np.c_[state[:,3:], model.acceleration(t,state)+smooth_arc(t,windows)[0]*impulses].ravel()
    def event(t, flat):
        return minimum_clearance(flat.reshape(-1,6),command(t).as_matrix())[0]-100.
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,[start,end],initial.ravel(),method='DOP853',max_step=60.,rtol=2e-12,
        atol=np.tile([1e-5]*3+[1e-9]*3,len(initial)),events=event,dense_output=True)
    if not sol.success: raise RuntimeError(sol.message)
    ts=np.linspace(start,sol.t[-1],int(np.ceil((sol.t[-1]-start)/30))+1)
    ys=sol.sol(ts).T.reshape(-1,len(initial),6)
    audit_t=np.linspace(start,sol.t[-1],int(np.ceil((sol.t[-1]-start)/2))+1)
    audit=sol.sol(audit_t).T.reshape(-1,len(initial),6)
    geometry=audit_geometry(audit_t,audit,command)
    beam=beam_and_rate_screen(env,audit_t,audit,command)
    roots=crossing_witnesses(ts,ys,command,lambda t,ids:sol.sol(t).reshape(-1,6)[ids],retain=8)
    position_error=float(length(audit[:,:,:3]-interpolation(audit_t)).max())
    velocity_error=float(length(audit[:,:,3:]-interpolation(audit_t,1)).max())
    fine_interp=CubicHermiteSpline(ts,ys[:,:,:3],ys[:,:,3:])
    interpolation_error=float(length(audit[:,:,:3]-fine_interp(audit_t)).max())
    coverage=[]
    for t in np.arange(start,min(start+r['delay_s'],sol.t[-1])+1,300.):
        y=sol.sol(t).reshape(-1,6)
        entry=local_coverage(y,env.at(t),y.mean(axis=0),limb_count=16,rotation=.137,interior=8)
        entry.update(time_s=float(t),window_inside_target=entry['window_centre_distance_m']+10000.<=4*K.MOON_RADIUS)
        coverage.append(entry)
        if not entry['window_inside_target']:break
    path=RUN/'coupled_sun_16.npz'
    np.savez_compressed(path,t=ts,state=ys,impulses=impulses,frames=command(ts).as_matrix())
    fine=dict(solar_sources=16,max_step_s=60.,rtol=2e-12,nfev=sol.nfev,end_s=float(sol.t[-1]),
        duration_s=float(sol.t[-1]-start),completed_prefix=bool(sol.t[-1]==end),
        clearance_event=bool(len(sol.t_events[0])),geometry=geometry,beam_and_rate=beam,
        between_sample_crossings=roots,extra_moving_window_samples=coverage,
        mean_delta_v_m_s=float(length(impulses).mean()),
        maximum_delta_v_m_s=float(length(impulses).max()),
        maximum_thrust_acceleration_m_s2=float(2*length(impulses).max()/r['burn_s']),
        raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path))
    out['runs'].append(fine)
    out['independent_replay']=dict(maximum_position_difference_m=position_error,
        maximum_velocity_difference_m_s=velocity_error,comparison_step_s=2.,
        fine_30s_interpolation_position_error_m=interpolation_error,placement_gate_m=50.)
    out['accepted_departure']=bool(fine['completed_prefix'] and
        geometry['swept_separation']['possible_conflict_pairs']==0 and
        out['coarse_two_second_audit']['swept_separation']['possible_conflict_pairs']==0 and
        geometry['minimum_sampled_clearance_m']>=100. and roots['crossing_near_encounters']==0 and
        beam['minimum_sampled_beam_margin_deg']>0 and beam['maximum_sampled_rate_deg_s']<=.1 and
        fine['maximum_thrust_acceleration_m_s2']<=.001 and position_error<=50.)
    out['guard_scope']='Swept square exclusions conditional on 0.15 m/s2 centre acceleration and 0.1 degree/s command-rate bounds, with finite geometry and an independent replay; not a continuous full-model or engineering certificate.'
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    save();print(json.dumps(dict(stage='complete',accepted_departure=out['accepted_departure'],
        replay=out['independent_replay'],geometry=geometry,elapsed_s=out['elapsed_s'])),flush=True)


if __name__=='__main__':main()
