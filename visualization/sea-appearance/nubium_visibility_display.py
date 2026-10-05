"""Project the external CIE99 Earth-optics diagnostic into the sea view.

Only Earth's direct light receives the ocular kernel. A 4--5 degree taper omits
small computed halo tails, whose boundary magnitude is reported. This is a
display approximation, not a full retinal or neural colour-appearance model.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat
from scipy.ndimage import map_coordinates


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--build',type=Path,required=True)
    args=p.parse_args(); out=args.directory
    sys.path.insert(0,str(args.root))
    sys.path.insert(0,str(args.root/'visualization/reference-renderer/seas'))
    sys.path.append(str(args.root/'visualization/sea-appearance'))
    from nubium_visibility_view import tone,save_rgb
    from nubium_human_view import digest
    from nubium_guides import Camera
    import kernels
    record=json.loads((out/'nubium-visibility-guide.json').read_text())
    with np.load(out/'nubium-visibility-xyz.npz') as z:
        xyz=z['xyz'].astype(float)
    with np.load(args.build) as z:
        sky=z['array/sky']; elevations=z['array/sky_el']
        az_step=z['array/params'][kernels.P_SKY_AZ_STEP]
    data=loadmat(out/'nubium-earth-optics.mat')
    optical=data['optical']; foreground=data['foreground']
    cam=Camera(record['camera']); focus=Camera(record['focus_camera'])
    # The small-angle CIE kernel's sampling distortion is under 1% at 5 degrees.
    # Project rays, rather than stretching a circular halo in image coordinates.
    corner=focus.direction(np.array([0,0,1024,1024]),np.array([0,1024,0,1024]))
    px,py=cam.project(corner)
    x0=max(0,math.floor(px.min())-2); x1=min(xyz.shape[1],math.ceil(px.max())+2)
    y0=max(0,math.floor(py.min())-2); y1=min(xyz.shape[0],math.ceil(py.max())+2)
    for row in range(y0,y1,32):
        end=min(row+32,y1); ss=2
        yy,xx=np.mgrid[row*ss:end*ss,x0*ss:x1*ss]
        rays=cam.direction((xx.ravel()+.5)/ss,(yy.ravel()+.5)/ss)
        fx,fy=focus.project(rays)
        theta=np.degrees(np.arccos(np.clip(rays@focus.f,-1,1)))
        weight=.5*(1+np.cos(np.pi*np.clip(theta-4,0,1)))
        direct=np.stack([map_coordinates(optical[...,c],[fy-.5,fx-.5],order=1,
                         mode='constant',cval=0) for c in range(3)],-1)*weight[:,None]
        background=np.stack([kernels.sky_xyz(sky,elevations,az_step,d) for d in rays])
        samples=(direct+background).reshape((end-row)*ss,(x1-x0)*ss,3)
        xyz[row:end,x0:x1]=samples.reshape(end-row,ss,x1-x0,ss,3).mean(axis=(1,3))
    white=record['display_curve']['white_world_cd_m2']
    wide_rgb,_=tone(xyz,white); focus_rgb,_=tone(foreground+optical,white)
    save_rgb(out/'nubium-midnight-visibility-view.png',wide_rgb)
    save_rgb(out/'nubium-earth-visibility-view.png',focus_rgb)
    # A crop retains the same physical field and tone curve. It is not a larger
    # Earth in the wide composition and does not include the distant horizon.
    yy,xx=np.mgrid[:1024,:1024]
    rays=focus.direction(xx.ravel()+.5,yy.ravel()+.5)
    theta=np.degrees(np.arccos(np.clip(rays@focus.f,-1,1))).reshape(1024,1024)
    ring=(theta>=4)&(theta<=5)
    core=theta<.65
    result=dict(schema='terluna.visualization.nubium-visibility-display/1',
                producer_sha256=digest(__file__),reference=record,
                eye_optics=dict(model='Official HDR-VDP-3.0.7 CIE99 OTF',age=24,
                                spatial_colour='Same achromatic optical kernel on XYZ; no neural colour adaptation',
                                scope='Direct Earth only; background and sea remain unblurred',
                                fft_negative_energy_fraction_zeroed=float(data['negative_fraction'][0,0]),
                                approximation='Uniform angular sampling over 10-degree field; halo tapered from 4 to 5 degrees in the wide view',
                                maximum_taper_region_halo_cd_m2=float(optical[...,1][ring].max()),
                                core_luminance_p10_p50_p90=np.percentile((foreground+optical)[...,1][core],[10,50,90]).tolist()),
                outputs={},
                evidence='Visibility-informed display reference; global monotone tone curve and gamut compression, with conditional ocular optics. Not an exact naked-eye prediction.',
                inputs={p.name:digest(p) for p in [out/'nubium-earth-optics.mat',out/'nubium-visibility-guide.json',args.build]})
    for filename in ['nubium-midnight-visibility-view.png','nubium-earth-visibility-view.png']:
        result['outputs'][filename]=digest(out/filename)
    (out/'nubium-visibility-display.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['eye_optics'],indent=2))


if __name__=='__main__':main()
