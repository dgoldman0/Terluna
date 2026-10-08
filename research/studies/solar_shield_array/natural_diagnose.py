"""Identify relative-motion compression at the independently replayed encounter."""
import time
import numpy as np
from protection.dynamics.service_search import free_command
from .joint_common import spline,state_at
from .natural_common import *


def main():
    before=time.monotonic();cpu=time.process_time();arrays,producer=setup(['research/studies/solar_shield_array/natural_diagnose.py'])
    name='natural_free_slow_fine.json';p=HERE/'results'/name;parent=json.loads(p.read_text());producer['inputs'][name]=digest(p)
    env=environment(3.);cmd=free_command(env,parent['arrival_s']+21600,arrival=parent['arrival_s'],turn_start=6.,turn_hours=6.)
    paths=[]
    for r in parent['runs']:
        path=ROOT/r['raw_path'];assert digest(path)==r['raw_sha256'];paths.append(spline(np.load(path)))
    X=arrays['initial_epoch_state'][:,:3];X-=X.mean(axis=0);rows=[]
    for t in [0.,21600.,43200.,parent['runs'][-1]['end_s']]:
        y=state_at(paths[0] if t<=21600 else paths[1],t);q=y[:,:3]-y[:,:3].mean(axis=0)
        A=np.linalg.lstsq(X,q,rcond=None)[0];residual=q-X@A
        rows.append(dict(time_s=t,affine_singular_values=np.linalg.svd(A,compute_uv=False).tolist(),
            affine_rms_residual_m=float(np.sqrt(np.mean(np.sum(residual**2,axis=1)))),
            body_frame_position_standard_deviations_m=np.std(q@cmd(t).as_matrix(),axis=0).tolist()))
    t=parent['runs'][-1]['end_s'];y=state_at(paths[1],t);i,j=parent['runs'][-1]['witness']['pair'];f=cmd(t).as_matrix()
    displacement=(y[i,:3]-y[j,:3])@f;velocity=(y[i,3:]-y[j,3:])@f-np.cross(cmd(t,1),displacement)
    out=dict(schema='terluna.research.natural-diagnosis/1',producer=producer,completed=True,accepted_return=False,accepted_cycle=False,
        formation=rows,encounter=dict(time_s=t,pair=[i,j],body_displacement_m=displacement.tolist(),body_relative_velocity_m_s=velocity.tolist()),
        scope='Measured affine compression and encounter kinematics, not a proof of impossibility or a tested avoidance maneuver.')
    write('natural_diagnose.json',out,before,cpu);print(out,flush=True)
if __name__=='__main__':main()
