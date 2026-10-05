"""Coupled return corrections and tighter replay from corresponding prefixes."""
import argparse
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command
from protection.dynamics.return_escape import rolled_return
from protection.dynamics.optical import length
from .closure_compact import propagate
from .packing_screen import load_raw
from .repair_screen import RUN, source_parents, write_product
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',required=True)
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--end-hours',type=float,default=18.)
    a=parser.parse_args();before=time.monotonic()
    names=['repair_screen.json','closure_validation.json','closure_sizing.json']
    if a.case.startswith('local'):names.append('repair_local_refined.json' if a.case=='local_refined' else 'repair_local_screen.json')
    if a.replay:names.append('repair_trial_'+a.case+'.json')
    p,sources=source_parents(names,['research/studies/solar_shield_array/repair_run.py'])
    used=10.+p['repair_screen.json']['elapsed_s']
    for path in (HERE/'results').glob('repair_*.json'):
        if path.name not in ['repair_screen.json','repair_trade.json','repair_budget.json']:
            used+=json.loads(path.read_text()).get('elapsed_s',0.)
    allowance=min(1600,int(3600-used-600))
    if allowance<=0:raise ValueError('No coupled-propagation budget remains')
    signal.alarm(allowance);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    if a.case.startswith('local'):
        case=p['repair_local_refined.json' if a.case=='local_refined' else 'repair_local_screen.json']
    else:
        case=next(r for r in p['repair_screen.json']['cases'] if r['schedule']['label']==a.case)
    raw=load_raw(case);start=float(raw['t'][0]);arrival=float(raw['t'][-1]);end=min(a.end_hours*3600.,arrival)
    parent=p['closure_validation.json' if a.replay else 'closure_sizing.json']['runs'][1]
    prefix=load_raw(parent);initial=prefix['state'][-1]
    if abs(prefix['t'][-1]-start)>1e-9:raise ValueError('Discontinuous prefix date')
    env=environment(3.);base=return_command(env,start,arrival,28800.);spec=case['schedule']
    command=rolled_return(base,start,arrival,spec.get('angle_deg',0.),spec.get('turn_s',5400.))
    label='repair_trial_'+a.case+('_fine' if a.replay else '')
    out=dict(schema='terluna.research.return-repair-coupled/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,stage_wall_s=allowance,threads=1,address_space_GiB=2),
        schedule=spec,arrival_s=arrival,arrival_turn_s=28800.,
        accepted_local_escape=False,accepted_return=False,accepted_cycle=False,accepted_fleet=False,
        scope='All 361 mass-loaded tiles from the corresponding saved hour-12 states; finite-source dynamic mutual shadows and unchanged optical mode. Coupled propagation stops at 100 m. A completed six-hour escape leaves arrival and recovery unpaid.',
        preceding_passage_source=parent['raw_path'],preceding_passage_unchanged=True)
    row,times,states=propagate(env,initial,start,end,command,raw['controls'],raw['windows'],raw['mass_ratio'],
        label,suns=16 if a.replay else 8,max_step=60. if a.replay else 120.,rtol=2e-12 if a.replay else 2e-10)
    path=RUN/(label+'.npz');(ROOT/row['raw_path']).rename(path)
    row['raw_path']=str(path.relative_to(ROOT))
    proposal=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
    row['maximum_independent_proposal_position_defect_m']=float(length(states[:,:,:3]-proposal(times)).max())
    out['run']=row
    if a.replay:
        old=p['repair_trial_'+a.case+'.json']['run'];coarse_raw=load_raw(old)
        c=CubicHermiteSpline(coarse_raw['t'],coarse_raw['state'][:,:,:3],coarse_raw['state'][:,:,3:])
        f=CubicHermiteSpline(times,states[:,:,:3],states[:,:,3:]);stop=min(times[-1],coarse_raw['t'][-1])
        maximum=0.;velocity=0.
        for t0 in np.arange(start,stop,3600.):
            t=np.linspace(t0,min(t0+3600.,stop),1801)
            maximum=max(maximum,float(length(c(t)-f(t)).max()))
            velocity=max(velocity,float(length(c(t,1)-f(t,1)).max()))
        error=maximum+row['reconstruction_error_m']+old['reconstruction_error_m']
        bound=min(row['minimum_conditional_surface_distance_m'],old['minimum_conditional_surface_distance_m'])-2*error
        out['comparison']=dict(maximum_position_difference_m=maximum,
            position_difference_plus_reconstruction_m=error,maximum_velocity_difference_m_s=velocity,
            same_witness_pair=row['sampled_witness']['pair']==old['sampled_witness']['pair'],
            event_time_difference_s=abs(row['end_s']-old['end_s']),
            conditional_clearance_after_two_centre_discrepancies_m=bound,
            scope='Source/integrator comparison is a measured uncertainty allowance, not a rigorous physical error bound.')
        out['accepted_local_escape']=bool(row['geometry_pass'] and old['geometry_pass'] and error<50. and bound>=100.)
    out['completed']=True;write_product(label+'.json',out,sources,before)
    print(json.dumps(dict(label=label,end_hours=row['end_s']/3600.,minimum_m=row['minimum_sampled_surface_distance_m'],
        witness=row['sampled_witness'],comparison=out.get('comparison'),elapsed_s=out['elapsed_s'])),flush=True)


if __name__=='__main__':main()
