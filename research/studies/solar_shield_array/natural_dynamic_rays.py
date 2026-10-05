"""Initialize ray ownership with a deforming, actually executed service passage.

Earlier static or greedy ray bundles can themselves overconstrain a free-state
LP. This branch transfers only ownership from the validated six-hour pattern;
it still imposes ray inequalities, not terminal lattice positions/velocities.
"""
import numpy as np
import json
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import linear_sum_assignment
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.service_constraints import service_rows as original_rows
from protection.dynamics.optical import length
from . import natural_arrival as runner
from .natural_common import HERE,ROOT,digest

original_setup=runner.setup
original_write=runner.write
path=None


def setup(extra):
    global path
    arrays,producer=original_setup([*extra,'research/studies/solar_shield_array/natural_dynamic_rays.py'])
    parent=json.loads((HERE/'results/joint_regression.json').read_text());entry=parent['runs'][0]
    assert digest(ROOT/entry['raw_path'])==entry['raw_sha256']
    raw=np.load(ROOT/entry['raw_path']);path=CubicHermiteSpline(raw['t'],raw['state'][:,:,:3],raw['state'][:,:,3:])
    producer['inputs']['joint_regression.json']=digest(HERE/'results/joint_regression.json')
    return arrays,producer


def rows(env,base,maps,arrival,arcs,**kwargs):
    if arrival==0:return original_rows(env,base,maps,arrival,arcs,**kwargs)
    initial=path(0.);local=(initial-initial.mean(axis=0))@solar_frame(env.at(0.))
    natural=base([arrival])[0];target=natural[:,:3].mean(axis=0)+local@solar_frame(env.at(arrival)).T
    _,assignment=linear_sum_assignment(length(natural[:,None,:3]-target[None,:,:]))
    def seed(ts):
        out=base(ts).copy()
        for k,t in enumerate(np.atleast_1d(ts)):
            tau=np.clip(t-arrival,0,21600);q=path(tau)
            rel=(q-q.mean(axis=0))@solar_frame(env.at(tau))
            out[k,:,:3]=out[k,:,:3].mean(axis=0)+rel[assignment]@solar_frame(env.at(t)).T
        return out
    kwargs['assignment']=seed
    return original_rows(env,base,maps,arrival,arcs,**kwargs)


def write(name,data,before,cpu):
    data['scope']='Four-arc rated and energy-capped fits; receiver-ray ownership is initialized from the actually deforming validated first passage with Hungarian member assignment, not from a frozen grid. No endpoint positions or velocities are constrained. This remains one selected assignment branch per case.'
    data['ray_seed']='Transported actual first six-hour service deformation'
    original_write('natural_dynamic_rays.json',data,before,cpu)


if __name__=='__main__':
    runner.setup=setup;runner.service_rows=rows;runner.write=write
    original_save=runner.save
    runner.save=lambda name,**arrays:original_save('dynamic_'+name,**arrays)
    runner.main()
