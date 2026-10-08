"""Coupled verification of the strongest return proposal's limiting event."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command,independent_path,endpoint_controls,corrected_state
from protection.dynamics.packing_control import response_maps
from protection.dynamics.pattern_return import encounter_scan
from protection.dynamics.optical import length
from .packing_screen import parents,load_raw
from .closure_propagate import propagate
from .closure_screen import RUN
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE


def main():
    before=time.monotonic();names=['closure_screen.json','closure_hardware.json','closure_sizing.json','packing_validation.json']
    p,sources=parents(names);used=sum(p[n]['elapsed_s'] for n in names[:3])
    allowance=min(1200,int(4800-used-600));signal.alarm(allowance)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    sources['research/studies/solar_shield_array/closure_return.py']=digest(__file__)
    screen=p['closure_screen.json'];index=screen.get('selected_case',screen.get('diagnostic_case'))
    if index is None:raise ValueError('No endpoint-feasible proposal to diagnose')
    selected=next(r for r in screen['refinements'] if r['case']==index)
    raw=load_raw(selected);case=screen['cases'][index]
    hardware=p['closure_sizing.json']
    if not hardware['accepted_passage']:raise ValueError('Hardware passage did not reach the return departure state')
    actual=load_raw(hardware['runs'][-1]);initial=actual['state'][-1];ratio=actual['mass_ratio']
    env=environment(3.);start=screen['start_s'];end=case['end_s']
    command=return_command(env,start,end,case['arrival_turn_hours']*3600.)
    free=independent_path(env,initial,start,end,command,mass_ratio=ratio)
    base=lambda t:free.sol(np.atleast_1d(t)).T.reshape(-1,361,6)
    maps=response_maps(env,lambda t:base([t])[0].mean(axis=0),start,end,raw['windows'],np.eye(3))
    target=raw['target']-raw['target'].mean(axis=0)+base([end])[0].mean(axis=0)
    controls,fit=endpoint_controls(base,maps,target,raw['windows'])
    if controls is None:raise ValueError('Hardware-state endpoint proposal infeasible')
    proposal_t=raw['t'];proposal=corrected_state(base,maps,controls,proposal_t)
    path=RUN/'hardware_return_proposal.npz'
    np.savez_compressed(path,t=proposal_t,state=proposal,initial=initial,controls=controls,
        windows=raw['windows'],mass_ratio=ratio,target=target)
    out=dict(schema='terluna.research.cycle-closure-return/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,this_stage_wall_s=allowance,threads=1,address_space_GiB=2),
        case=index,proposal_geometry_pass=selected['proposal_geometry_pass'],arrival_turn_hours=case['arrival_turn_hours'],
        accepted_return=False,accepted_fleet=False,handover_executed=False,second_service_executed=False,
        scope='Diagnostic return from the propagated hardware-loaded twelve-hour terminal states. Refit the selected endpoint proposal to the new actual initial states and freely evolving mean; stop at 100 m physical surface clearance. The original target offsets and relative velocity pattern remain a restricted assignment.',
        hardware_refit=fit,hardware_proposal_encounters=encounter_scan(proposal_t,proposal,command,retain=16),
        hardware_proposal_raw_path=str(path.relative_to(ROOT)),hardware_proposal_raw_sha256=digest(path),
        runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
        (HERE/'results/closure_return.json').write_text(json.dumps(out,indent=2)+'\n')
    for label,suns,step,rtol in [('return_coupled',8,120.,2e-10)]:
        row,times,states=propagate(env,initial,start,end,command,controls,raw['windows'],ratio,
                                   label,suns=suns,max_step=step,rtol=rtol)
        predictor=CubicHermiteSpline(proposal_t,proposal[:,:,:3],proposal[:,:,3:])
        row['maximum_proposal_position_defect_m']=float(length(states[:,:,:3]-predictor(times)).max())
        out['runs'].append(row);save()
        print(json.dumps(dict(stage=label,end_s=row['end_s'],surface_m=row['minimum_sampled_surface_distance_m'],
            witness=row['sampled_witness'],elapsed_s=time.monotonic()-before)),flush=True)
    out.update(completed=True,accepted_return=all(r['geometry_pass'] for r in out['runs']))
    save()


if __name__=='__main__':main()
