"""Unlabelled optical-depth geometry guides, separate from radiance displays.

The white fraction is 1-exp(-cloud optical depth), not cloud brightness.
The upper-layer guide retains only cells centred above 20 km. No turbulence,
photographic texture, tone adjustment or invented 3-D morphology is added.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from numba import njit, prange
from PIL import Image


@njit(parallel=True, cache=True)
def depths(origin, rays, edges, theta0, dtheta, beta, half_width):
    result = np.zeros(len(rays))
    for i in prange(len(rays)):
        if rays[i,2] > 0:
            result[i] = cloud_depth(origin,rays[i],edges,theta0,dtheta,beta,half_width)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path)
    p.add_argument('--scene', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--width', type=int, default=768)
    a = p.parse_args()
    root = (a.root or Path(__file__).resolve().parents[2]).resolve()
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(root/'visualization/sea-appearance'))
    from nubium_guides import Camera
    from illumination.cloud_light.microphysics import extinction
    global cloud_depth
    from illumination.cloud_light.volume import cloud_depth
    from research.studies.sea_appearance.nobili import checked, digest
    from shared.constants import MOON_RADIUS
    product = json.loads(a.scene.read_text())
    source = checked(root/product['cloud']['file'], product['cloud']['sha256'])
    with np.load(source) as z:
        meta = json.loads(str(z['metadata']))
        section = {k:z[k] for k in z.files if k != 'metadata'}
    parts, _ = extinction(section, section['density_kg_m3'], meta['droplets_cm3'])
    beta = np.stack([parts['qc']+parts['qr'], parts['qi']+parts['qs']+parts['qg']], -1)
    edges = MOON_RADIUS+section['height_edges_m']
    origin = np.array([0.,0.,MOON_RADIUS+product['camera']['eye_height_m']])
    dt = meta['dx_m']/MOON_RADIUS
    t0 = (section['column_offset'][0]-.5)*dt
    width = a.width
    height = round(width*product['camera']['frame'][1]/product['camera']['frame'][0])
    cam = Camera({**product['camera'],'frame':[width,height]})
    yy, xx = np.mgrid[:height,:width]
    rays = cam.direction(xx.ravel()+.5,yy.ravel()+.5)
    a.out.mkdir(parents=True,exist_ok=True)
    outputs = {}
    for label, minimum in [('all-condensates',0.),('upper-clouds',20000.)]:
        b = beta.copy()
        centre = (section['height_edges_m'][:-1]+section['height_edges_m'][1:])/2
        b[centre < minimum] = 0
        tau = depths(origin,rays,edges,t0,dt,b,product['cloud']['cross_ring_width_m']/2).reshape(height,width)
        image = np.rint(255*(-np.expm1(-tau))).astype('uint8')
        target = a.out/f'{label}.png'
        Image.fromarray(image).save(target)
        outputs[target.name] = digest(target)
        np.save(a.out/f'{label}-tau.npy',tau.astype('f4'))
    record = dict(schema='terluna.visualization.cloud-geometry/1',scene_sha256=digest(a.scene),
        camera={**product['camera'],'frame':[width,height]},
        reading_rule='Greyscale optical opacity: black is transparent, white is opaque. Upper-clouds retains cells centred above 20 km. Neither guide represents radiance or photographic appearance.',
        evidence='Exact piecewise-constant rays through the saved cloud extrusion; native 6-km cells, unresolved cross-ring shape.',
        producer_sha256=digest(__file__),outputs=outputs)
    (a.out/'geometry.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
