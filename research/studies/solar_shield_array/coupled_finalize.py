"""Typed result export and source-checked recovery of completed validation runs.

The original validator completed both IVPs and accounts, then JSON rejected a
NumPy boolean in its interval audit. Recover those hashed trajectories; audit
their Hermite reconstructions with their recorded reconstruction allowance.
--fresh performs the actual integrations with the corrected writer instead.
"""
import argparse
import copy
import json
from types import SimpleNamespace
import numpy as np
from . import coupled_validate as validation
from .coupled_common import HERE,ROOT,read_raw,spline,state_at,digest,inputs,write as write_result


def native(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:native(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [native(v) for v in value]
    return value


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fresh',action='store_true')
    ap.add_argument('--extra-finalize-cpu',type=float,default=0.);args=ap.parse_args()
    recovery=None
    if not args.fresh:
        parent_path=HERE/'results/coupled_validation_attempt.json'
        parent=json.loads(parent_path.read_text())
        for name,sha in parent['producer']['source_hashes'].items():
            if digest(ROOT/name)!=sha:raise ValueError('Validation source changed: '+name)
        for name,sha in parent['producer']['inputs'].items():
            if digest(HERE/'results'/name)!=sha:raise ValueError('Validation input changed: '+name)
        rows=[parent['coarse'],*parent['runs']]
        def cached(self,x,initial,start,end,fine=False,**kwargs):
            row=next(r for r in rows if r['start_s']==start and (r['suns']==16)==fine)
            raw=read_raw(row)
            np.testing.assert_allclose(x,raw['x'],atol=0,rtol=0)
            np.testing.assert_allclose(initial,raw['state'][0],atol=0,rtol=0)
            path=spline(raw)
            def sol(t):
                state=state_at(path,t)
                return state.ravel() if np.ndim(t)==0 else state.reshape(len(t),-1).T
            return copy.deepcopy(row),raw,SimpleNamespace(sol=sol)
        validation.StableExperiment.evaluate=cached
        recovery=dict(parent='coupled_validation_attempt.json',parent_sha256=digest(parent_path),
            completed_execution_cpu_s=parent['resources']['cpu_s'],additional_failed_export_allowance_s=60.,
            additional_finalize_attempt_cpu_s=args.extra_finalize_cpu,
            additional_failed_export_allowance_reason='Conservative charge for completed audit/coverage work after the last saved resource entry and before NumPy-boolean serialization failed.',
            audit_trajectory='60 s Hermite reconstruction of the completed fine IVPs; recorded dense-output reconstruction error is included in every point allowance.')
    def export(name,out,cpu,wall):
        if recovery:
            out['recovery']=recovery
            out['producer']['inputs'][recovery['parent']]=recovery['parent_sha256']
        out['producer']['source_hashes']['research/studies/solar_shield_array/coupled_finalize.py']=digest(HERE/'coupled_finalize.py')
        out['export_mode']='fresh integrations' if args.fresh else 'source-checked recovery of completed integrations'
        clean=native(out);write_result(name,clean,cpu,wall);out['resources']=clean['resources']
    validation.write=export
    validation.main()


if __name__=='__main__':main()
