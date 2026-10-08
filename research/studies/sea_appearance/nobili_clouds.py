"""Conditional Earth-lit cloud context for Nobili's selected flat-ring snapshot.

The saved 2-D section is extruded through a stated 20-km width. Molecular and
cloud scattering are solved together. Ground is a flat Lambertian sphere here;
crater shadows and terrain-weather interaction are not computed.
"""
from __future__ import annotations
import argparse,json,math,time
from pathlib import Path
import numpy as np
from research.studies.sea_appearance.nobili import ROOT,HERE,checked,digest
from research.studies.sea_appearance.lighting import Sky,Earthlight
from illumination.cloud_light.microphysics import extinction
from illumination.cloud_light.finite_source import estimate,direct_horizontal
from illumination.cloud_light.volume import cloud_depth
from illumination.sky.solved_sky import load


def build(out,photons):
    product=json.loads((HERE/'results/nobili-night.json').read_text());c=product['cloud'];path=checked(ROOT/c['file'],c['sha256'])
    with np.load(path) as f:section={k:f[k] for k in f.files if k!='metadata'};meta=json.loads(str(f['metadata']))
    air,settings=load(ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz')
    source=product['light'];earth=product['earth'];sky=Sky();model=Earthlight(sky)
    weights=model.weights(np.array([earth['phase_angle_deg']]),np.array([source['earth_distance_m']]),np.array([source['earth_sunlight_factor']]))[0]
    beta_parts,_=extinction(section,section['density_kg_m3'],meta['droplets_cm3'])
    beta=np.stack([beta_parts['qc']+beta_parts['qr'],beta_parts['qi']+beta_parts['qs']+beta_parts['qg']],-1)
    radius=air['edges'][0];edges=radius+section['height_edges_m'];origin=np.array([0.,0.,radius+2.])
    dtheta=meta['dx_m']/radius;theta0=(section['column_offset'][0]-.5)*dtheta;halfwidth=c['cross_ring_width_m']/2
    # Equatorial ring basis x=east, y=north, z=up. Azimuth measured from north.
    el,az=np.radians([earth['elevation_deg'],earth['azimuth_deg']]);source_dir=np.array([np.cos(el)*np.sin(az),np.cos(el)*np.cos(az),np.sin(el)])
    elev=np.array([.25,1,2,4,6,8,10,12,14,16,18,20,23,26,30,35,45,60,80])
    azimuth=np.unique(np.r_[np.arange(0,360,15),np.arange(240,301,2.5)])
    e,a=np.radians(np.meshgrid(elev,azimuth,indexing='ij'));directions=np.stack([np.cos(e)*np.sin(a),np.cos(e)*np.cos(a),np.sin(e)],-1).reshape(-1,3)
    radius_rad=math.radians(earth['diameter_deg']/2)
    print('Tracing',len(directions),'cloud sightlines',photons,'photons each',flush=True);start=time.monotonic()
    args=(origin,directions,source_dir,air['edges'],air['scattering'],air['absorption'],air['xyz']*weights[:,None],edges,theta0,dtheta)
    groups,cuts=estimate(*args,beta,halfwidth,.1,.85,.8,photons=photons,groups=16,seed=20330621,source_radius_rad=radius_rad)
    print('Cloud radiance ready in',round(time.monotonic()-start,1),'s',flush=True)
    clear,cc=estimate(*args,np.zeros_like(beta),halfwidth,.1,.85,.8,photons=photons,groups=16,seed=20330621,source_radius_rad=radius_rad)
    tau=np.array([cloud_depth(origin,d,edges,theta0,dtheta,beta,halfwidth) for d in directions]).reshape(len(elev),len(azimuth))
    direct=direct_horizontal(origin,source_dir,air['edges'],air['scattering']+air['absorption'],air['xyz']*weights[:,None],edges,theta0,dtheta,beta,halfwidth,source_radius_rad=radius_rad)
    flux,fc=estimate(origin,np.array([[0.,0.,1.]]),*args[2:],beta,halfwidth,.1,.85,.8,photons=8192,groups=16,seed=20330622,flux=True,source_radius_rad=radius_rad)
    # Disk-centre and limb optical depths, including all five condensate species.
    ee=np.linspace(earth['elevation_deg']-earth['diameter_deg']/2,earth['elevation_deg']+earth['diameter_deg']/2,21)
    disk_tau=[cloud_depth(origin,np.array([np.cos(np.radians(e))*np.sin(az),np.cos(np.radians(e))*np.cos(az),np.sin(np.radians(e))]),edges,theta0,dtheta,beta,halfwidth) for e in ee]
    record=dict(schema='terluna.research.nobili-cloud-light/1',input_scene_sha256=digest(HERE/'results/nobili-night.json'),
        cloud_section_sha256=c['sha256'],photons_per_pixel=photons,groups=16,source_radius_deg=earth['diameter_deg']/2,
        source='Uniform finite Earth disk, normalized to observed phase-curve irradiance. Surface spatial texture is not used for cloud lighting.',
        cross_ring_width_m=2*halfwidth,earth_diffuse_horizontal_lux=float(flux[0,:,1].mean()),
        earth_diffuse_standard_error_lux=float(flux[0,:,1].std(ddof=1)/4),earth_direct_horizontal_lux=float(direct[1]),
        earth_total_horizontal_lux=float(flux[0,:,1].mean()+direct[1]),truncated_paths=int(cuts.sum()+cc.sum()+fc.sum()),
        disk_cloud_transmission_centre=float(np.exp(-disk_tau[10])),disk_cloud_transmission_minmax=[float(np.exp(-max(disk_tau))),float(np.exp(-min(disk_tau)))],
        median_pixel_relative_standard_error=float(np.median(groups[:,:,1].std(axis=1,ddof=1)/4/np.maximum(groups[:,:,1].mean(axis=1),1e-15))),
        evidence='Monte Carlo radiative transfer through one extruded experimental cloud field; no crater terrain, cross-ring dynamics, or cloud forecast.',
        limits=['Local cloud shadows on the crater rim are not solved.','Solar diffuse light is only 0.0012 lux in the clear reference and is not included in this Earth-only cloud calculation.',
            'The uniform source disk approximates the gibbous Earth in scattering; the direct visual Earth retains its phase and geographic texture.',
            'The 20-km cross-ring extent is an assumed geometry.'],
        seconds=time.monotonic()-start,producer={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'illumination/cloud_light/volume.py',ROOT/'illumination/cloud_light/microphysics.py']})
    out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out/'cloud-light.npz',elevation_deg=elev,azimuth_deg=azimuth,cloudy_groups=groups.reshape(len(elev),len(azimuth),16,3),clear_groups=clear.reshape(len(elev),len(azimuth),16,3),cloud_tau=tau,
        cloud_beta=beta,height_edges_m=section['height_edges_m'],theta0=theta0,dtheta=dtheta,half_width=halfwidth,metadata=json.dumps(record))
    (out/'cloud-light.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT/'research/runs/sea_appearance/nobili');p.add_argument('--photons',type=int,default=1024);a=p.parse_args();build(a.out,a.photons)
if __name__=='__main__':main()
