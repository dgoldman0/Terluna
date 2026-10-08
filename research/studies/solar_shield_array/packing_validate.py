"""Replay the complete changed-packing service/departure with tighter dynamics."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.optical import length
from .packing_screen import parents, load_raw
from .packing_run import propagate, full_command, service_coverage
from .patterns_run import reference_shooting
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def comparison(coarse_t,coarse,fine_t,fine):
    a=CubicHermiteSpline(coarse_t,coarse[:,:,:3],coarse[:,:,3:])
    b=CubicHermiteSpline(fine_t,fine[:,:,:3],fine[:,:,3:])
    position=velocity=0.
    start,end=max(coarse_t[0],fine_t[0]),min(coarse_t[-1],fine_t[-1])
    for left in np.arange(start,end,3600.):
        times=np.r_[np.arange(left,min(left+3600.,end),2.),min(left+3600.,end)]
        position=max(position,float(length(a(times)-b(times)).max()))
        velocity=max(velocity,float(length(a(times,1)-b(times,1)).max()))
    return dict(maximum_saved_interpolant_position_difference_m=position,
        maximum_saved_interpolant_velocity_difference_m_s=velocity,sample_step_s=2.)


def main():
    p,sources=parents(['packing_screen.json','packing.json']);candidate=p['packing.json']
    if not candidate['accepted_departure']:raise ValueError('No accepted local candidate to replay')
    used=sum(q['elapsed_s'] for q in p.values())
    restart=ROOT/'research/runs/solar_shield_array/packing_replay_restart.json'
    interrupted_wall=(json.loads(restart.read_text())['estimated_partial_wall_s']+1 if restart.exists() else 0)
    allowance=int(3600-used-interrupted_wall-600)
    if allowance<=0:raise ValueError('Replay has no remaining resource budget')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();env=environment(1.);reference,_=reference_shooting(env,21600.)
    sources['research/studies/solar_shield_array/packing_validate.py']=digest(__file__)
    service=next(r for r in candidate['services'] if r['config']==candidate['selected_config'])
    departure=candidate['departures'][-1];sr=load_raw(service);dr=load_raw(departure)
    controls,basis,windows=dr['controls'],dr['basis'],dr['windows']
    command=full_command(env,candidate['schedule'])
    out=dict(schema='terluna.research.packing-validation/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in p}),
        selected_config=candidate['selected_config'],schedule=candidate['schedule'],
        budget=dict(prior_numerical_wall_s=used,interrupted_partial_replay_wall_s=interrupted_wall,
            this_stage_wall_s=allowance,accounting_reserved_s=600),
        accepted_service_departure=False,accepted_fleet=False,full_return_executed=False,
        handover_executed=False,hardware_mass_feedback_propagated=False,
        recurring_cost_measured=False,continuous_certificate=False,
        collection_capacity_W=None,delivered_power_W=None,runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/packing_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    row,st,sy=propagate(env,sr['state'][0],command,np.zeros_like(controls),basis,windows,
        0.,21600.,'replay_service',suns=16,step=60.,rtol=2e-12)
    row['coverage']=service_coverage(env,reference,st,sy)
    row['comparison']=comparison(sr['t'],sr['state'],st,sy)
    row['comparison']['position_discrepancy_plus_reconstruction_bounds_m']=(
        row['comparison']['maximum_saved_interpolant_position_difference_m']+
        row['fine_interpolation_position_error_m']+service['fine_interpolation_position_error_m'])
    out['runs'].append(row);save()
    if row['local_geometry_pass'] and row['coverage']['passed']:
        fine,dt,dy=propagate(env,sy[-1],command,controls,basis,windows,
            21600.,43200.,'replay_departure',suns=16,step=60.,rtol=2e-12)
        fine['comparison']=comparison(dr['t'],dr['state'],dt,dy)
        fine['comparison']['position_discrepancy_plus_reconstruction_bounds_m']=(
            fine['comparison']['maximum_saved_interpolant_position_difference_m']+
            fine['fine_interpolation_position_error_m']+departure['fine_interpolation_position_error_m'])
        out['runs'].append(fine)
        out['accepted_service_departure']=bool(fine['local_geometry_pass'] and
            all(r['comparison']['position_discrepancy_plus_reconstruction_bounds_m']<=50. for r in out['runs']))
    out['completed']=True
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    save()
    print(json.dumps(dict(stage='validation_complete',accepted=out['accepted_service_departure'],
        elapsed_s=out['elapsed_s'],max_rss_MiB=out['max_rss_MiB'])),flush=True)


if __name__=='__main__':main()
