"""Bounded natural-return producers, rooted in exact committed loaded states."""
import json,resource,time
from pathlib import Path
import numpy as np
from shared.provenance import constants_used,constants_changed
from .fleet_run import environment,digest
from .cycling_search import ROOT,HERE

RUN=ROOT/'research/runs/solar_shield_array/natural_return'
BASE=['protection/dynamics/service_search.py','research/studies/solar_shield_array/natural_common.py']


def setup(extra=()):
    RUN.mkdir(parents=True,exist_ok=True)
    seed=json.loads((HERE/'results/joint_seed.json').read_text())
    assert not constants_changed(seed['producer']['constants'])
    for name,h in seed['producer']['source_hashes'].items():
        assert digest(ROOT/name)==h,(name,h)
    sources=dict(seed['producer']['source_hashes'])
    sources.update({p:digest(ROOT/p) for p in [*BASE,*extra]})
    producer=dict(source_hashes=sources,constants=constants_used(sources),inputs={'joint_seed.json':digest(HERE/'results/joint_seed.json')})
    return {k:np.asarray(v) for k,v in seed['arrays'].items()},producer


def write(name,data,before,cpu):
    data['resources']=dict(wall_s=time.monotonic()-before,cpu_s=time.process_time()-cpu,
        max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,numerical_threads=1,
        new_raw_MiB=sum(p.stat().st_size for p in RUN.glob('*') if p.is_file())/1024**2)
    (HERE/'results'/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def save(name,**arrays):
    p=RUN/(name+'.npz');np.savez_compressed(p,**arrays)
    return dict(raw_path=str(p.relative_to(ROOT)),raw_sha256=digest(p))
