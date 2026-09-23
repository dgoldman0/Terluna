"""Bake the clear-sky atlases from illumination/sky into the immersion's sky asset.

Repacks previously generated optical outputs; never recalculates the atmosphere.
"""
import argparse,base64,gzip,hashlib,json
from pathlib import Path
import numpy as np

def pack(source,dest):
    out={'version':1,'units':'photopic-weighted linear sRGB; sky is cd/m2-equivalent, irradiance lux-equivalent','provenance':[],'worlds':{}}
    for name in ['moon','earth','moon_no_ozone']:
        p=source/f'{name}_atlas.npz';d=np.load(p)
        rgb=d['rgb'][:,::2,::2,:].astype(np.float64)
        # Preserve luminance while bringing signed spectral RGB inside the positive octant.
        Y=np.maximum(rgb@np.array([.2126,.7152,.0722]),0)
        mn=np.min(rgb,axis=-1);factor=np.where(mn<0,np.minimum(1,Y/(Y-mn+1e-30)),1)
        rgb=np.maximum(0,Y[...,None]+factor[...,None]*(rgb-Y[...,None]))
        scale=np.maximum(np.max(rgb,axis=(1,2,3)),1e-25)
        packed=(rgb/scale[:,None,None,None]).astype('<f2')
        rebuilt=packed.astype(float)*scale[:,None,None,None]
        maxerr=float(np.max(np.abs(rebuilt-rgb)/scale[:,None,None,None]))
        world={k:np.array(d[k]).tolist() for k in ['suns','direct','direct_horizontal','diffuse','cloud_direct','cloud_diffuse']}
        world.update(width=rgb.shape[2],height=rgb.shape[1],scale=scale.tolist(),radius=1737400 if name!='earth' else 6371000,
            data=base64.b64encode(gzip.compress(packed.tobytes(),compresslevel=9,mtime=0)).decode(),quantization_error_relative_to_frame_max=maxerr)
        out['worlds'][name]=world
        out['provenance'].append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'shape':list(d['rgb'].shape),'derived_shape':list(rgb.shape),'new_atmospheric_calculation':False})
        print(name,rgb.shape,len(world['data']),maxerr)
    dest.write_text(json.dumps(out,separators=(',',':')))
REPO=Path(__file__).resolve().parents[2]
if __name__=='__main__':
    # Source: the illumination domain's clear-sky atlases; output: this experience's baked sky.
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive_data',type=Path,nargs='?',default=REPO/'illumination/sky/data')
    p.add_argument('--out',type=Path,default=REPO/'immersion/assets/sky/atmosphere.json')
    a=p.parse_args();pack(a.archive_data,a.out)
