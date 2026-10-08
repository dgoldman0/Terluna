"""Retain a per-tile optical inventory for every rejected coupled energy trial."""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from shared.provenance import constants_used, constants_changed
from protection.dynamics.energy_coast import coast_command
from protection.dynamics.collection import uneclipsed_collection, integrate_inventory
from .energy_screen import RUN
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE


def main():
    signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    before=time.monotonic();env=environment(1.)
    names=['energy_screen.json','energy.json','energy_replay.json']
    p={n:json.loads((HERE/'results'/n).read_text()) for n in names}
    sources={}
    for parent in p.values():
        for path,h in parent['producer']['source_hashes'].items():
            if digest(ROOT/path)!=h:raise ValueError('Source changed: '+path)
        if constants_changed(parent['producer']['constants']):raise ValueError('Constants changed')
        sources.update(parent['producer']['source_hashes'])
    for path in ['protection/dynamics/collection.py','research/studies/solar_shield_array/energy_collection.py']:
        sources[path]=digest(ROOT/path)
    rows=[(r,r['schedule'],r['case']) for r in p['energy.json']['runs']]
    rows += [(r,p['energy_replay.json']['schedule'],p['energy_replay.json']['case']) for r in p['energy_replay.json']['runs']]
    out=dict(schema='terluna.research.energy-trial-collection/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={n:digest(HERE/'results'/n) for n in names},raw_inputs={}),
        evidence='Collection on every saved coupled trial through its clearance stop; each alternative is a separate local experiment, not a simultaneous fleet. No favorable comparison across unequal durations is inferred.',
        accepted_fleet=False,collection_capacity_W=None,delivered_power_W=None,
        sample_step_s=60.,sources=16,trials=[])
    for entry,spec,index in rows:
        if digest(ROOT/entry['raw_path'])!=entry['raw_sha256']:raise ValueError('Raw changed')
        out['producer']['raw_inputs'][entry['raw_path']]=entry['raw_sha256']
        raw=np.load(ROOT/entry['raw_path']);start,end=map(float,raw['t'][[0,-1]])
        command=coast_command(env,start,p['energy_screen.json']['end_s'],**spec)
        state_at=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
        times=np.unique(np.r_[np.arange(start,end,60.),end]);ledger={}
        for t in times:
            state=np.c_[state_at(t),state_at(t,1)]
            data=uneclipsed_collection(env,t,state,command(t).as_matrix(),suns=16)
            for key,value in data.items():
                if key!='reflected_force_N':ledger.setdefault(key,[]).append(value)
        ledger={k:np.array(v) for k,v in ledger.items()}
        coarse=np.unique(np.r_[np.arange(0,len(times),2),len(times)-1]);summary={}
        for key,value in ledger.items():
            result=integrate_inventory(times,value);result.pop('per_member_energy_J')
            older=integrate_inventory(times[coarse],value[coarse])
            result['time_coarsening_relative_energy_change']=abs(older['energy_J']-result['energy_J'])/max(result['energy_J'],1.)
            summary[key]=result
        path=RUN/f'trial_{index}_collection.npz'
        np.savez_compressed(path,t=times,member_id=np.arange(raw['state'].shape[1]),**ledger)
        out['trials'].append(dict(case=index,schedule=spec,start_s=start,end_s=end,duration_s=end-start,
            accepted_departure=False,original_force_sources=entry['suns'],optical=summary,
            raw_path=str(path.relative_to(ROOT)),raw_sha256=digest(path)))
        print(json.dumps(dict(stage='trial_collection',case=index,elapsed_s=time.monotonic()-before)),flush=True)
    out.update(completed=True,elapsed_s=time.monotonic()-before,
        max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        refinement_scope='60 versus 120 s integration on the same rejected traces; solar-force dynamics were not independently replayed for rejected candidates.')
    if any(digest(ROOT/path)!=h for path,h in sources.items()):raise ValueError('Source changed')
    (HERE/'results/energy_collection.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
