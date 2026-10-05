"""Size proposed return power before propagating its added mass from the epoch."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.mass_feedback import allocated_mass_ratio
from protection.dynamics.cycle_control import return_command
from protection.dynamics.attitude_load import rigid_square_load,disturbance_torque_bound
from protection.dynamics.pattern_return import smooth_arc
from protection.dynamics.optical import length
from .packing_screen import parents,load_raw
from .packing_run import full_command,service_coverage
from .patterns_run import reference_shooting
from .closure_propagate import propagate
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE,SCENARIO


def main():
    before=time.monotonic();names=['packing_validation.json','packing_budget.json','closure_screen.json','closure_hardware.json']
    p,sources=parents(names);used=sum(p[n]['elapsed_s'] for n in names[2:])
    allowance=min(1300,int(4800-used-600));signal.alarm(allowance)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    sources['research/studies/solar_shield_array/closure_sizing.py']=digest(__file__)
    screen=p['closure_screen.json'];index=screen.get('selected_case',screen.get('diagnostic_case'))
    proposal=load_raw(next(r for r in screen['refinements'] if r['case']==index));case=screen['cases'][index]
    env=environment(3.);cmd=return_command(env,screen['start_s'],case['end_s'],case['arrival_turn_hours']*3600.)
    prop=SCENARIO['propulsion'];optical_mass=SCENARIO['baseline_areal_mass_kg_m2']*10000.**2
    conversion=prop['exhaust_velocity_m_s']/(2*prop['efficiency']*np.cos(np.deg2rad(prop['cant_deg'])))
    dt=proposal['t'];nominal=rigid_square_load(cmd,dt)['force_per_mass'];envelope=[]
    for t,state in zip(dt,proposal['state']):
        envelope.append(disturbance_torque_bound(env,t,state,cmd(t).as_matrix())[2])
    thrust=np.array([smooth_arc(t,proposal['windows']) for t in dt])@length(proposal['controls']).T
    return_peak=(nominal[:,None]+np.array(envelope)+thrust).max(axis=0)*conversion
    data=p['packing_budget.json']['datasets']
    previous_peak=np.maximum.reduce([data[k]['evaluations'][0]['cost']['per_member_peak_propulsion_W'] for k in ['service','departure']])/optical_mass
    # Preserve the stated 25% installation margin and allow 10% refit headroom.
    design_specific=1.1*np.maximum(return_peak,previous_peak)
    ratio=allocated_mass_ratio(design_specific,prop)
    old=p['packing_validation.json'];command=full_command(env,old['schedule'])
    initial=load_raw(old['runs'][0])['state'][0];reference,_=reference_shooting(env,21600.)
    out=dict(schema='terluna.research.cycle-return-hardware-sizing/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,this_stage_wall_s=allowance,threads=1,address_space_GiB=2),
        selected_return_case=index,fixed_extra_mass_ratio=.2,refit_power_headroom=1.1,
        sizing_method='Maximum proposed return translation plus nominal attitude and conservative gravity/radiation couple envelope, or prior actual prefix peak; 10% refit headroom and scenario 25% installation margin. Scalar proportional-load mass feedback; sampled actual loads checked later.',
        mass_distribution='All fixed and power hardware mass uniformly distributed over each rigid square. Optical mass and area unchanged. Constant carried mass; propellant depletion and detailed hardware layout remain outside this passage.',
        per_member_total_to_optical_mass_ratio=ratio.tolist(),per_member_design_specific_peak_W_kg=design_specific.tolist(),
        optical_mass_kg=float(optical_mass*361),total_allocated_mass_kg=float(optical_mass*ratio.sum()),
        power_hardware_kg=float(optical_mass*(ratio-1.2).sum()),
        allocated_installed_power_W=float(optical_mass*(ratio-1.2).sum()*prop['specific_power_W_kg']),
        accepted_passage=False,accepted_return=False,accepted_fleet=False,runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
        (HERE/'results/closure_sizing.json').write_text(json.dumps(out,indent=2)+'\n')
    save()
    controls=np.zeros((361,2,3));windows=np.array([[21600.,28800.],[32400.,39600.]])
    for label,start,end in [('sized_service',0.,21600.),('sized_departure',21600.,43200.)]:
        row,times,states=propagate(env,initial,start,end,command,controls,windows,ratio,label)
        if start==0:
            row['coverage']=service_coverage(env,reference,times,states)
            row['geometry_pass']=row['geometry_pass'] and row['coverage']['passed']
        out['runs'].append(row);save()
        print(json.dumps(dict(stage=label,gates=row['geometry_pass'],end_s=row['end_s'],
            minimum_surface_m=row['minimum_sampled_surface_distance_m'],elapsed_s=time.monotonic()-before)),flush=True)
        if not row['geometry_pass']:break
        initial=states[-1]
    out.update(completed=True,accepted_passage=len(out['runs'])==2 and all(r['geometry_pass'] for r in out['runs']));save()


if __name__=='__main__':main()
