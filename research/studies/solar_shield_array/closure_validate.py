"""One complete tighter hardware-loaded replay, including the return witness."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used
from protection.dynamics.cycle_control import return_command
from protection.dynamics.optical import length
from protection.dynamics.mass_feedback import MassiveParallelPattern
from protection.dynamics.compact_shadow import CompactMassivePattern
from .packing_screen import parents,load_raw
from .packing_run import full_command,service_coverage
from .patterns_run import reference_shooting
from .closure_compact import propagate
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE


def main():
    before=time.monotonic();names=['packing_validation.json','closure_screen.json','closure_hardware.json','closure_sizing.json','closure_return.json']
    p,sources=parents(names);used=sum(p[n]['elapsed_s'] for n in names[1:])
    allowance=min(3100,int(4800-used-600));signal.alarm(allowance)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    sources['research/studies/solar_shield_array/closure_validate.py']=digest(__file__)
    for path in ['protection/dynamics/compact_shadow.py','research/studies/solar_shield_array/closure_compact.py']:
        sources[path]=digest(ROOT/path)
    hardware=p['closure_sizing.json'];ret=p['closure_return.json'];env=environment(3.)
    first=load_raw(p['packing_validation.json']['runs'][0]);initial=first['state'][0]
    command=full_command(env,p['packing_validation.json']['schedule']);reference,_=reference_shooting(env,21600.)
    old=hardware['runs']+ret['runs'];ratio=np.array(hardware['per_member_total_to_optical_mass_ratio'])
    out=dict(schema='terluna.research.cycle-closure-validation/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={n:digest(HERE/'results'/n) for n in names}),
        budget=dict(prior_wall_s=used,this_stage_wall_s=allowance,threads=1,address_space_GiB=2),
        scope='Same hardware masses, initial epoch state and controls; complete sequential replay with twice the solar force sources, half the integration step and tighter tolerance. Exact shadow-union compaction is checked against the original backend. No coarse terminal state is substituted between legs.',
        accepted_hardware_passage=False,accepted_return=False,accepted_fleet=False,continuous_certificate=False,runs=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
        (HERE/'results/closure_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    comparisons=[]
    for k,entry in enumerate(old):
        raw=load_raw(entry);t=float(raw['t'][-1]);state=raw['state'][-1]
        cmd=command if k<2 else return_command(env,entry['start_s'],entry['intended_end_s'],ret['arrival_turn_hours']*3600.)
        tick=time.monotonic();a=MassiveParallelPattern(env,cmd,ratio,suns=8).acceleration(t,state,diagnostics=True);old_s=time.monotonic()-tick
        tick=time.monotonic();b=CompactMassivePattern(env,cmd,ratio,suns=8).acceleration(t,state,diagnostics=True);new_s=time.monotonic()-tick
        test=dict(leg=k,time_s=t,maximum_acceleration_difference_m_s2=float(np.max(abs(a[0]-b[0]))),
            maximum_visibility_difference=float(np.max(abs(a[1]-b[1]))),maximum_sail_difference_m_s2=float(np.max(abs(a[2]-b[2]))),
            original_force_wall_s=old_s,compact_force_wall_s=new_s)
        if test['maximum_acceleration_difference_m_s2']>1e-14 or test['maximum_visibility_difference']>1e-10:
            raise ValueError('Compact shadow model changes force')
        comparisons.append(test)
    out['force_backend_comparison']=comparisons;save()
    print(json.dumps(dict(stage='force_compaction',comparisons=comparisons)),flush=True)
    for k,entry in enumerate(old):
        raw=load_raw(entry)
        if k==2:command=return_command(env,entry['start_s'],entry['intended_end_s'],ret['arrival_turn_hours']*3600.)
        row,times,states=propagate(env,initial,entry['start_s'],entry['intended_end_s'],command,
            raw['controls'],raw['windows'],ratio,f'closure_fine_{k}',suns=16,max_step=60.,rtol=2e-12)
        coarse=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        fine=CubicHermiteSpline(times,states[:,:,:3],states[:,:,3:])
        end=min(times[-1],raw['t'][-1]);maximum=0.;velocity=0.
        for begin in np.arange(times[0],end,3600.):
            common=np.linspace(begin,min(begin+3600.,end),int(np.ceil((min(begin+3600.,end)-begin)/2))+1)
            maximum=max(maximum,float(length(coarse(common)-fine(common)).max()))
            velocity=max(velocity,float(length(coarse(common,1)-fine(common,1)).max()))
        row['comparison']=dict(maximum_position_difference_m=maximum,
            position_difference_plus_reconstruction_m=maximum+entry['reconstruction_error_m']+row['reconstruction_error_m'],
            maximum_velocity_difference_m_s=velocity,event_time_difference_s=abs(times[-1]-raw['t'][-1]),
            same_witness_pair=entry['sampled_witness']['pair']==row['sampled_witness']['pair'])
        if k==0:
            row['coverage']=service_coverage(env,reference,times,states)
            row['geometry_pass']=row['geometry_pass'] and row['coverage']['passed']
        out['runs'].append(row);save()
        print(json.dumps(dict(stage='fine_leg',leg=k,end_s=row['end_s'],comparison=row['comparison'],elapsed_s=time.monotonic()-before)),flush=True)
        if k<2 and not row['geometry_pass']:break
        initial=states[-1]
    out.update(completed=True,accepted_hardware_passage=len(out['runs'])>=2 and all(r['geometry_pass'] and
        r['comparison']['position_difference_plus_reconstruction_m']<50. for r in out['runs'][:2]),
        accepted_return=len(out['runs'])==3 and all(r['geometry_pass'] and r['comparison']['position_difference_plus_reconstruction_m']<50. for r in out['runs']))
    save()


if __name__=='__main__':main()
