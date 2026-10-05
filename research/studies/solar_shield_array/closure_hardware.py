"""Repropagate the complete existing passage with explicit added hardware mass."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.mass_feedback import allocated_mass_ratio
from protection.dynamics.optical import length
from .packing_screen import parents,load_raw
from .packing_run import full_command,service_coverage
from .patterns_run import reference_shooting
from .closure_propagate import propagate
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO


def main():
    before=time.monotonic();names=['packing_validation.json','packing_budget.json','closure_screen.json']
    p,sources=parents(names);used=p['closure_screen.json']['elapsed_s']
    allowance=min(1800,int(4800-used-600));signal.alarm(allowance)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    for path in ['protection/dynamics/mass_feedback.py','research/studies/solar_shield_array/closure_propagate.py',
                 'research/studies/solar_shield_array/closure_hardware.py']:
        sources[path]=digest(ROOT/path)
    old=p['packing_validation.json'];env=environment(1.);command=full_command(env,old['schedule'])
    data=p['packing_budget.json']['datasets'];mass=SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    peak=np.maximum.reduce([data[k]['evaluations'][0]['cost']['per_member_peak_propulsion_W'] for k in ['service','departure']])
    ratio=allocated_mass_ratio(peak/mass,SCENARIO['propulsion'])
    initial=load_raw(old['runs'][0])['state'][0];reference,_=reference_shooting(env,21600.)
    out=dict(schema='terluna.research.cycle-hardware-feedback/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,this_stage_wall_s=allowance,threads=1,address_space_GiB=2),
        fixed_extra_mass_ratio=.2,mass_distribution='All fixed and allocated power mass distributed uniformly on each rigid square; optical area and optical base mass unchanged.',
        per_member_total_to_optical_mass_ratio=ratio.tolist(),
        optical_mass_kg=float(mass*len(ratio)),total_allocated_mass_kg=float(mass*ratio.sum()),
        power_hardware_kg=float(mass*(ratio-1.2).sum()),
        scope='Same initial positions/velocities and command as the accepted optical-mass passage, with mass feedback in sail acceleration. Power mass is an allocation awaiting actual-load verification; no return or full hardware design.',
        accepted_passage=False,accepted_return=False,accepted_fleet=False,runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
        (HERE/'results/closure_hardware.json').write_text(json.dumps(out,indent=2)+'\n')
    controls=np.zeros((361,2,3));windows=np.array([[21600.,28800.],[32400.,39600.]])
    for label,start,end,old_entry in zip(['hardware_service','hardware_departure'],[0.,21600.],[21600.,43200.],old['runs']):
        row,times,states=propagate(env,initial,start,end,command,controls,windows,ratio,label)
        raw=load_raw(old_entry);spline=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        row['maximum_change_from_optical_mass_trace_m']=float(length(states[:,:,:3]-spline(times)).max())
        if start==0:
            row['coverage']=service_coverage(env,reference,times,states)
            row['geometry_pass']=row['geometry_pass'] and row['coverage']['passed']
        out['runs'].append(row);save()
        print(json.dumps(dict(stage=label,gates=row['geometry_pass'],end_s=row['end_s'],
            minimum_surface_m=row['minimum_sampled_surface_distance_m'],mass_effect_m=row['maximum_change_from_optical_mass_trace_m'],elapsed_s=time.monotonic()-before)),flush=True)
        if not row['geometry_pass']:break
        initial=states[-1]
    out.update(completed=True,accepted_passage=len(out['runs'])==2 and all(r['geometry_pass'] for r in out['runs']))
    save()


if __name__=='__main__':main()
