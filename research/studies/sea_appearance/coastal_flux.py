"""Refine the ground-flux estimate independently of a saved angular cloud map."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from research.studies.sea_appearance.nobili import ROOT, digest, checked
from research.studies.sea_appearance.coastal_clouds import direction
from research.studies.sea_appearance.lighting import Sky, Earthlight
from illumination.cloud_light.finite_source import estimate
from illumination.sky.solved_sky import load
from illumination.sky.solved_transport import SUN_RADIUS_RAD


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scene',type=Path,required=True)
    p.add_argument('--cloud',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--photons',type=int,default=65536)
    a=p.parse_args()
    if a.photons < 16 or a.photons%16:
        raise ValueError('Use a positive multiple of sixteen photons')
    product=json.loads(a.scene.read_text())
    with np.load(a.cloud) as z:
        meta=json.loads(str(z['metadata']))
        beta=z['cloud_beta'];height=z['height_edges_m'];theta0=float(z['theta0']);dtheta=float(z['dtheta']);half=float(z['half_width'])
    if meta['input_scene_sha256']!=digest(a.scene):raise ValueError('Cloud and scene differ')
    air_path=checked(ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz',meta['molecular_fields_sha256'])
    air,_=load(air_path);radius=air['edges'][0];origin=np.array([0.,0.,radius+product['camera']['eye_height_m']])
    light,earth=product['light'],product['earth'];components={}
    for index,name in enumerate(meta['source_map_scope']):
        if name=='sun':
            source=direction(light['sun_elevation_deg'],light['sun_azimuth_deg']);xyz=air['xyz'];r=SUN_RADIUS_RAD
        else:
            sky=Sky();w=Earthlight(sky).weights(np.array([earth['phase_angle_deg']]),np.array([light['earth_distance_m']]),np.array([light['earth_sunlight_factor']]))[0]
            source=direction(earth['elevation_deg'],earth['azimuth_deg']);xyz=air['xyz']*w[:,None];r=math.radians(earth['diameter_deg']/2)
        groups,cuts=estimate(origin,np.array([[0.,0.,1.]]),source,air['edges'],air['scattering'],air['absorption'],xyz,
            radius+height,theta0,dtheta,beta,half,.1,.85,.8,photons=a.photons,groups=16,seed=203306261+index,
            flux=True,source_radius_rad=r)
        if cuts.sum() or not np.isfinite(groups).all() or groups.min()<0:raise ValueError('Invalid flux refinement')
        direct=meta['components'][name]['direct_horizontal_lux'];diffuse=groups[0,:,1]
        components[name]=dict(diffuse_horizontal_lux=float(diffuse.mean()),diffuse_standard_error_lux=float(diffuse.std(ddof=1)/4),
            direct_horizontal_lux=direct,total_horizontal_lux=float(diffuse.mean()+direct),photon_groups=diffuse.tolist())
    record=dict(schema='terluna.research.coastal-ground-light/1',input_scene_sha256=digest(a.scene),cloud_map_sha256=digest(a.cloud),
        photons_per_source=a.photons,groups=16,seed=203306261,truncated_paths=0,components=components,
        retained_formal_sources=meta['retained_formal_sources'],
        total_horizontal_lux=sum(c['total_horizontal_lux'] for c in components.values())+sum(meta['retained_formal_sources'].values()),
        evidence='A larger independent Monte Carlo ground-flux estimate for the same saved cloud field and source; angular sky radiance is unchanged.',
        producer={str(f.relative_to(ROOT)):digest(f) for f in [Path(__file__),ROOT/'illumination/cloud_light/finite_source.py',ROOT/'illumination/cloud_light/volume.py']})
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':main()
