"""Prepare matched-angle tests of the corrected Earth and its display image.

This checks whether the display retains detectable interior structure. It does
not test equality of appearance, colour, or recognizability of particular land.
Run the study's visibility_backend.m on the emitted backend.json afterwards.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image
from scipy.io import savemat


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--backend-template',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(args.root))
    from research.studies.sea_appearance.visibility import featureless_inner,digest
    record=json.loads((args.directory/'nubium-visibility-guide.json').read_text())
    world=args.directory/'nubium-earth-xyz.npy'
    display_file=args.directory/'nubium-earth-visibility-view.png'
    physical=np.load(world)[...,1].astype(float)
    png=np.asarray(Image.open(display_file).convert('RGB'))/255
    rgb=np.where(png<=.04045,png/12.92,((png+.055)/1.055)**2.4)
    display=.1+99.9*(rgb@np.array([.2126729,.7151522,.0721750]))
    n=physical.shape[0]
    focal=n/2/np.tan(np.radians(record['focus_camera']['horizontal_fov_deg']/2))
    ppd=focal*np.pi/180
    q=(np.arange(n)+.5-n/2)/focal; x,y=np.meshgrid(q,q)
    ct=1/np.sqrt(1+x*x+y*y); omega=ct**3/focal**2
    r=np.hypot(x,y)/np.tan(np.radians(record['geometry']['diameter_deg']/2))
    cases=[]
    for name,data,age in [('reprojected_world',physical,24),
                          ('display_age24',display,24),('display_age70',display,70)]:
        counter,roi=featureless_inner(data,r,omega,ct)
        path=args.out/(name+'.mat')
        bg=float(np.median(data[(r>3)&(r<4)]))
        savemat(path,dict(reference=data,comparison=counter,roi=roi.astype(np.uint8),
                         radius_fraction=r,ppd=ppd,background=bg),do_compression=True)
        cases.append(dict(name=name,file=str(path),age=age,identical=False))
    cfg=json.loads(args.backend_template.read_text())
    cfg.update(out=str(args.out),cases=cases)
    (args.out/'backend.json').write_text(json.dumps(cfg,indent=2)+'\n')
    product=dict(schema='terluna.visualization.display-visibility-inputs/1',
                 producer_sha256=digest(__file__),pixels_per_degree=ppd,
                 display=dict(peak_cd_m2=100,black_cd_m2=.1,encoding='sRGB/D65'),
                 inputs={p.name:digest(p) for p in [world,display_file]},
                 cases=[{**c,'file':Path(c['file']).name} for c in cases],
                 reading_rule='Detection of interior-texture removal at matched angular size. World stimuli have no pre-applied eye kernel. Display stimuli include the illustrative ocular halo; HDR-VDP then predicts looking at that display. Passing is not perceptual equivalence or colour validation.')
    (args.out/'stimuli.json').write_text(json.dumps(product,indent=2)+'\n')


if __name__=='__main__':main()
