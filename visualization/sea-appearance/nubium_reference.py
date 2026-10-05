"""Render the supplemental Nubium midnight frame with the existing sea tracer.

Keeps its XYZ radiance and Monte Carlo errors. A saved build cache permits
sample-count refinements without repeating the spectral sky/coast preparation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--width',type=int,default=768)
    p.add_argument('--spp-side',type=int,default=4)
    args=p.parse_args()
    root=(args.root or Path(__file__).resolve().parents[2]).resolve()
    sys.path.insert(0,str(root))
    tracer=root/'visualization/reference-renderer/seas'
    sys.path.insert(0,str(tracer))
    import scene
    import render
    source=root/'research/studies/sea_appearance/results/nubium-midnight.json'
    product=json.loads(source.read_text()); s=product['scene']
    common=scene.Common()
    common.regimes['scenes'].append(dict(coast=s['coast'],moment=s['moment'],hour=s['hour']))
    common.inputs[str(source.relative_to(root))]=hashlib.sha256(source.read_bytes()).hexdigest()
    common.inputs[str(Path(__file__).relative_to(root))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    scene.VIEWS[(s['coast'],s['moment'])]=(s['camera']['view_azimuth_deg'],s['camera']['pitch_deg'])
    args.out.mkdir(parents=True,exist_ok=True)
    cache=args.out/'scene-build.npz'
    identity=hashlib.sha256(json.dumps(dict(product=product,width=args.width,
        sources={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),tracer/'scene.py',tracer/'kernels.py']}),sort_keys=True).encode()).hexdigest()
    original=scene.build
    for name in ['light','sea_state','coast']:
        fn=getattr(scene,name)
        def timed(*a,_fn=fn,_name=name,**kw):
            print('Preparing',_name,flush=True); t=time.monotonic()
            r=_fn(*a,**kw); print(_name,'ready in',round(time.monotonic()-t,1),'s',flush=True)
            return r
        setattr(scene,name,timed)
    def cached(*a,**kw):
        if cache.exists():
            with np.load(cache,allow_pickle=False) as z:
                if str(z['identity'])!=identity: raise ValueError('Scene build cache identity differs')
                arrays={k[6:]:z[k] for k in z.files if k.startswith('array/')}
                light={k[6:]:z[k] for k in z.files if k.startswith('light/')}
                return z['camera'],arrays,json.loads(str(z['record'])),json.loads(str(z['geometry'])),light
        camera,arrays,record,geometry,light=original(*a,**kw)
        np.savez(cache,identity=identity,camera=camera,record=json.dumps(record),geometry=json.dumps(geometry),
                 **{'array/'+k:v for k,v in arrays.items()},**{'light/'+k:v for k,v in light.items()})
        return camera,arrays,record,geometry,light
    scene.build=cached
    render.render_frame(common,s['coast'],s['moment'],args.width,args.width*3//2,args.spp_side,20380307,
                        s['camera']['horizontal_fov_deg'],args.out)


if __name__=='__main__':
    main()
