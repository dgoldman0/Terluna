"""Separate controlled reduced dynamics from the LP and coupled replay defects."""
import time
import numpy as np
from protection.dynamics.service_reduced import reduced_path
from protection.dynamics.service_search import free_command
from protection.dynamics.optical import length
from .joint_common import spline
from .natural_common import *


def main():
    before=time.monotonic();cpu=time.process_time()
    arrays,producer=setup(['protection/dynamics/service_reduced.py','research/studies/solar_shield_array/natural_force_audit.py'])
    names=['natural_face_refined_coarse.json','natural_refine.json'];parents=[]
    for name in names:
        p=HERE/'results'/name;producer['inputs'][name]=digest(p);parents.append(json.loads(p.read_text()))
    coupled,proposal=parents;data=np.load(ROOT/proposal['raw_path']);assert digest(ROOT/proposal['raw_path'])==proposal['raw_sha256']
    env=environment(3.);finish=coupled['runs'][-1]['end_s'];cmd=free_command(env,proposal['arrival_s']+21600,tilt=0.,arrival=proposal['arrival_s'],turn_start=6.,turn_hours=6.)
    sol=reduced_path(env,arrays['initial_epoch_state'],0.,finish,cmd,controls=data['controls'],windows=data['windows'],mass_ratio=arrays['mass_ratio'])
    proposed=spline(data);rows=[]
    for entry in coupled['runs']:
        raw=np.load(ROOT/entry['raw_path']);assert digest(ROOT/entry['raw_path'])==entry['raw_sha256'];ts=raw['t'];state=sol.sol(ts).T.reshape(-1,361,6)
        rows.append(dict(stage=entry['stage'],maximum_reduced_vs_linear_proposal_m=float(length(state[:,:,:3]-proposed(ts)).max()),
            maximum_coupled_vs_nonlinear_reduced_m=float(length(raw['state'][:,:,:3]-state[:,:,:3]).max())))
    ts=np.unique(np.r_[np.arange(0.,finish,300.),finish]);states=sol.sol(ts).T.reshape(-1,361,6)
    out=dict(schema='terluna.research.natural-force-audit/1',producer=producer,completed=True,accepted_return=False,accepted_cycle=False,stages=rows,
        scope='Identical executed controls and initial loaded states. Coupled-minus-reduced includes mutual shadows, finite-area gravity and finite-Sun force differences. The passive finite-area/source comparison is a separate scale check; do not label this an exact isolated-shadow error.',
        **save('controlled_reduced_face',t=ts,state=states))
    write('natural_force_audit.json',out,before,cpu);print(rows,out['resources'],flush=True)
if __name__=='__main__':main()
