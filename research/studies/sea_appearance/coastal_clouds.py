"""Conditional Sun-plus-Earth cloud transport for the selected coastal scenes.

This consumes a saved equatorial CM1 section. Its cross-ring width is an
explicit scenario, and the Lambertian ground has no resolved terrain shadows.
"""
from __future__ import annotations
import argparse
import json
import math
import time
from pathlib import Path
import numpy as np

from research.studies.sea_appearance.nobili import ROOT, checked, digest
from research.studies.sea_appearance.lighting import Sky, Earthlight
from illumination.cloud_light.microphysics import extinction
from illumination.cloud_light.finite_source import estimate, direct_horizontal
from illumination.cloud_light.volume import cloud_depth
from illumination.sky.solved_sky import load
from illumination.sky.solved_transport import SUN_RADIUS_RAD


def direction(elevation, azimuth):
    e, a = np.radians([elevation, azimuth])
    return np.array([np.cos(e)*np.sin(a), np.cos(e)*np.cos(a), np.sin(e)])


def build(scene_path, out, photons, sun_photons, fine=False):
    product = json.loads(scene_path.read_text())
    cloud = product['cloud']
    section_path = checked(ROOT/cloud['file'], cloud['sha256'])
    with np.load(section_path) as f:
        section = {k: f[k] for k in f.files if k != 'metadata'}
        meta = json.loads(str(f['metadata']))
    if abs(meta['observer_latitude_deg']) > .01:
        raise ValueError('This adapter requires the equatorial east/north/up basis')
    air_path = ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz'
    air, settings = load(air_path)
    light, earth = product['light'], product['earth']
    sky = Sky()
    weights = Earthlight(sky).weights(np.array([earth['phase_angle_deg']]),
        np.array([light['earth_distance_m']]), np.array([light['earth_sunlight_factor']]))[0]
    parts, _ = extinction(section, section['density_kg_m3'], meta['droplets_cm3'])
    beta = np.stack([parts['qc']+parts['qr'], parts['qi']+parts['qs']+parts['qg']], -1)
    radius = air['edges'][0]
    edges = radius+section['height_edges_m']
    origin = np.array([0., 0., radius+product['camera']['eye_height_m']])
    dtheta = meta['dx_m']/radius
    theta0 = (section['column_offset'][0]-.5)*dtheta
    halfwidth = cloud['cross_ring_width_m']/2
    elev = np.array([.25,1,2,4,6,8,10,12,14,16,18,20,23,26,30,35,40,45,50,55,60,70,80,89])
    if fine:
        elev = np.unique(np.r_[elev, np.arange(10., 56., 1.5)])
    azimuth = np.unique(np.r_[np.arange(0,360,15), np.arange(240,306,2.5)])
    e, a = np.radians(np.meshgrid(elev,azimuth,indexing='ij'))
    directions = np.stack([np.cos(e)*np.sin(a),np.cos(e)*np.cos(a),np.sin(e)],-1).reshape(-1,3)
    sources = [('earth', direction(earth['elevation_deg'],earth['azimuth_deg']),
                air['xyz']*weights[:,None], math.radians(earth['diameter_deg']/2), photons),
               ('sun', direction(light['sun_elevation_deg'],light['sun_azimuth_deg']),
                air['xyz'], SUN_RADIUS_RAD, sun_photons)]
    summed = np.zeros((len(directions),16,3))
    clear_sum = summed.copy()
    components = {}
    cuts_total = 0
    start = time.monotonic()
    for index, (name, source_dir, xyz, source_radius, count) in enumerate(sources):
        if count == 0:
            formal=light['earth_lux' if name=='earth' else 'sun_lux']
            components[name]=dict(mode='Unmodified formal clear-sky contribution; cloud modulation is unresolved',photons_per_direction=0,total_horizontal_lux=formal)
            print(name, 'retaining formal clear-sky contribution',formal,'lux; no cloud-modulation claim',flush=True)
            continue
        print(f'{name}: {len(directions)} directions, {count} photons/direction', flush=True)
        args = (origin,directions,source_dir,air['edges'],air['scattering'],air['absorption'],xyz,edges,theta0,dtheta)
        seed = 20330626+index*100000
        kwargs = dict(photons=count,groups=16,seed=seed,source_radius_rad=source_radius)
        groups, cuts = estimate(*args,beta,halfwidth,.1,.85,.8,**kwargs)
        clear, cc = estimate(*args,np.zeros_like(beta),halfwidth,.1,.85,.8,**kwargs)
        direct = direct_horizontal(origin,source_dir,air['edges'],air['scattering']+air['absorption'],xyz,
            edges,theta0,dtheta,beta,halfwidth,source_radius_rad=source_radius)
        flux, fc = estimate(origin,np.array([[0.,0.,1.]]),*args[2:],beta,halfwidth,.1,.85,.8,
            photons=max(2048,count*8),groups=16,seed=seed+1,flux=True,source_radius_rad=source_radius)
        summed += groups
        clear_sum += clear
        cuts_total += int(cuts.sum()+cc.sum()+fc.sum())
        components[name] = dict(mode='Monte Carlo through molecular air and clouds',photons_per_direction=count,diffuse_horizontal_lux=float(flux[0,:,1].mean()),
            diffuse_standard_error_lux=float(flux[0,:,1].std(ddof=1)/4),direct_horizontal_lux=float(direct[1]),
            total_horizontal_lux=float(flux[0,:,1].mean()+direct[1]),source_radius_deg=math.degrees(source_radius))
        print(name, components[name], 'elapsed',round(time.monotonic()-start,1),flush=True)
    if cuts_total or not np.isfinite(summed).all() or summed.min() < 0:
        raise ValueError('Incomplete or invalid cloud transport; no final product written')
    tau = np.array([cloud_depth(origin,d,edges,theta0,dtheta,beta,halfwidth) for d in directions])
    disk_e = np.linspace(earth['elevation_deg']-earth['diameter_deg']/2,earth['elevation_deg']+earth['diameter_deg']/2,21)
    disk_tau = np.array([cloud_depth(origin,direction(el,earth['azimuth_deg']),edges,theta0,dtheta,beta,halfwidth) for el in disk_e])
    record = dict(schema='terluna.research.coastal-cloud-light/1',input_scene_sha256=digest(scene_path),
        cloud_section_sha256=cloud['sha256'],molecular_fields_sha256=digest(air_path),components=components,
        total_horizontal_lux=sum(c['total_horizontal_lux'] for c in components.values()),
        clear_total_horizontal_lux=light['sun_lux']+light['earth_lux'],
        earth_total_horizontal_lux=components['earth']['total_horizontal_lux'],
        cross_ring_width_m=2*halfwidth,groups=16,truncated_paths=cuts_total,
        source_map_scope=[name for name,c in components.items() if c['photons_per_direction']>0],
        retained_formal_sources={name:c['total_horizontal_lux'] for name,c in components.items() if c['photons_per_direction']==0},
        disk_cloud_transmission_centre=float(np.exp(-disk_tau[10])),
        disk_cloud_transmission_minmax=[float(np.exp(-disk_tau.max())),float(np.exp(-disk_tau.min()))],
        median_pixel_relative_standard_error=float(np.median(summed[:,:,1].std(axis=1,ddof=1)/4/np.maximum(summed[:,:,1].mean(axis=1),1e-15))),
        evidence='Multi-source Monte Carlo radiative transfer through one saved, extruded 2-D experimental cloud field. No local terrain-weather prediction.',
        limits=['Finite Earth is uniform for scattered light; historical texture is used only in the direct visual disk.',
            'Cross-ring extent and particle phase functions are declared scenarios. Fine 3-D cloud structure is unresolved.',
            'Ground is a Lambertian sphere; local cloud/terrain shadows and shoreline circulation are not solved.',
            'Sources with zero requested photons retain their formal clear-sky term in the display; their cloud modulation is not solved. Zero-hit low-photon twilight runs are not evidence of zero light.'],
        seconds=time.monotonic()-start,
        producer={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'illumination/cloud_light/finite_source.py',ROOT/'illumination/cloud_light/volume.py',ROOT/'illumination/cloud_light/microphysics.py']})
    out.mkdir(parents=True,exist_ok=True)
    shape=(len(elev),len(azimuth),16,3)
    np.savez_compressed(out/'cloud-light.npz',elevation_deg=elev,azimuth_deg=azimuth,
        cloudy_groups=summed.reshape(shape),clear_groups=clear_sum.reshape(shape),cloud_tau=tau.reshape(shape[:2]),
        cloud_beta=beta,height_edges_m=section['height_edges_m'],theta0=theta0,dtheta=dtheta,half_width=halfwidth,metadata=json.dumps(record))
    (out/'cloud-light.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scene',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--photons',type=int,default=512,help='Earth-source photons per direction; zero retains the formal background')
    p.add_argument('--sun-photons',type=int,default=128,help='Solar photons per direction; zero retains the formal background')
    p.add_argument('--fine',action='store_true')
    a=p.parse_args()
    build(a.scene.absolute(),a.out,a.photons,a.sun_photons,a.fine)


if __name__ == '__main__':
    main()
