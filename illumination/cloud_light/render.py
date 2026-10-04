"""Absolute cloud radiance and surface-light estimates for an admitted CM1 section."""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
from illumination.cloud_light.microphysics import extinction
from illumination.cloud_light.volume import estimate, direct_horizontal
from illumination.sky.solved_sky import load,digest,photometry,horizontal_integral
from illumination.sky.solved_transport import evaluate

ROOT=Path(__file__).resolve().parents[2]
SCHEMA='terluna.illumination.cloud-evening-scene/1'


def render(section_path,output,camera_azimuth_deg,photons=256,shape=(24,48),width_km=200.,ice_g=.8,
           seed=819,observer_height_m=1.6,surface_photons=16384):
    start=time.monotonic()
    with np.load(section_path,allow_pickle=False) as f:
        section={k:f[k].copy() for k in f.files if k!='metadata'}
        meta=json.loads(str(f['metadata']))
    if meta['schema']!='terluna.climate.cloud-section/1':
        raise ValueError('Unsupported cloud section')
    source=ROOT/'illumination/sky/results/solved_sky.json'
    product=json.loads(source.read_text())
    path=ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz'
    if digest(path)!=product['additional_checks']['numerical_checks.json']['inputs'][path.name]:
        raise ValueError('Molecular archive differs')
    air,settings=load(path)
    beta_species,radii=extinction(section,section['density_kg_m3'],meta['droplets_cm3'])
    beta=np.stack((beta_species['qc']+beta_species['qr'],beta_species['qi']+beta_species['qs']+beta_species['qg']),axis=-1)
    radius=air['edges'][0]
    edges=radius+section['height_edges_m']
    dtheta=meta['dx_m']/radius
    theta0=(section['column_offset'][0]-.5)*dtheta
    origin=np.array([0.,0.,radius+observer_height_m])
    elevations=np.linspace(2.,65.,shape[0])
    azimuths=camera_azimuth_deg+np.linspace(-50.,50.,shape[1])
    e,a=np.radians(np.meshgrid(elevations,azimuths,indexing='ij'))
    directions=np.stack((np.cos(e)*np.cos(a),np.cos(e)*np.sin(a),np.sin(e)),axis=-1).reshape(-1,3)
    args=(origin,directions,section['sun_local'],air['edges'],air['scattering'],air['absorption'],air['xyz'],edges,theta0,dtheta)
    options=dict(photons=photons,groups=16,seed=seed)
    cloudy,cut=estimate(*args,beta,width_km*500,.1,.85,ice_g,**options)
    print('cloud paths',round(time.monotonic()-start,2),'s',flush=True)
    clear,clear_cut=estimate(*args,np.zeros_like(beta),width_km*500,.1,.85,ice_g,**options)
    sun_alt=np.degrees(np.arcsin(section['sun_local'][2]))
    sun_az=np.degrees(np.arctan2(section['sun_local'][1],section['sun_local'][0]))
    relative_az=np.abs((azimuths-sun_az+180)%360-180)
    background=evaluate(air,settings,[sun_alt],elevations,relative_az)[0]@air['xyz']
    # Positive unbiased estimates remain stable for opaque cloud paths.
    # Paired clear paths are retained as an independent molecular-sky check.
    means=cloudy.mean(axis=1).reshape((*shape,3))
    errors=(cloudy.std(axis=1,ddof=1)/4).reshape((*shape,3))
    groups=cloudy.reshape((*shape,16,3))
    if observer_height_m!=1.6:
        background=clear.mean(axis=1).reshape((*shape,3))
    # Independent cosine-sampled hemispheres include cloud attenuation and cloud reflection.
    surface_args=(origin,np.array([[0.,0.,1.]]),section['sun_local'],air['edges'],air['scattering'],air['absorption'],air['xyz'],edges,theta0,dtheta)
    light,light_cut=estimate(*surface_args,beta,width_km*500,.1,.85,ice_g,photons=surface_photons,groups=16,seed=seed+991,flux=True)
    clear_light,clear_light_cut=estimate(*surface_args,np.zeros_like(beta),width_km*500,.1,.85,ice_g,
                                       photons=surface_photons,groups=16,seed=seed+991,flux=True)
    grid_e=90*np.linspace(0,1,41)**2;grid_a=np.linspace(0,180,65)
    spectrum=evaluate(air,settings,[sun_alt],grid_e,grid_a)
    base_flux=float(horizontal_integral(spectrum@air['xyz'][:,1],grid_e,grid_a)[0])
    _,direct,_=photometry(air,np.array([sun_alt]))
    direct_cloud=direct_horizontal(origin,section['sun_local'],air['edges'],air['scattering']+air['absorption'],
                                   air['xyz'],edges,theta0,dtheta,beta,width_km*500)
    surface=light[0,:,1].mean()+direct_cloud[1]
    surface_se=light[0,:,1].std(ddof=1)/4
    clear_surface=base_flux+float(direct[0,1]) if observer_height_m==1.6 else clear_light[0,:,1].mean()
    contrast=(means[:,:,1]-background[:,:,1])/background[:,:,1]
    # A black viewing recess leaves an explicit angular opening toward the cloud.
    aperture_options=dict(photons=surface_photons,groups=16,seed=seed+1991,flux=True,
                          flux_mu_min=np.sin(np.radians(5)),flux_mu_max=np.sin(np.radians(25)),
                          flux_az_min=np.radians(camera_azimuth_deg-10),flux_az_max=np.radians(camera_azimuth_deg+10))
    aperture,aperture_cut=estimate(*surface_args,beta,width_km*500,.1,.85,ice_g,**aperture_options)
    aperture_clear,aperture_clear_cut=estimate(*surface_args,np.zeros_like(beta),width_km*500,.1,.85,ice_g,**aperture_options)
    ae=np.linspace(5,25,21);aa=camera_azimuth_deg-sun_az+np.linspace(-10,10,41)
    ar=evaluate(air,settings,[sun_alt],ae,aa)@air['xyz'][:,1]
    aperture_base=float(horizontal_integral(ar,ae,aa)[0]/2)
    recess=aperture[0,:,1]
    truncations=int(cut.sum()+clear_cut.sum()+light_cut.sum()+clear_light_cut.sum()+aperture_cut.sum()+aperture_clear_cut.sum())
    records=dict(schema=SCHEMA,case=meta['case'],snapshot=meta['snapshot'],time_day=meta['time_day'],
                 observer_latitude_deg=meta['observer_latitude_deg'],observer_longitude_deg=meta['observer_longitude_deg'],
                 observer_sun_elevation_deg=float(sun_alt),observer_height_m=observer_height_m,
                 camera_azimuth_in_ring_basis_deg=camera_azimuth_deg,cross_ring_width_km=width_km,
                 phase_parameters=dict(liquid_g=.85,frozen_g=ice_g),photons_per_pixel=photons,groups=16,
                 surface_photons=surface_photons,shape=list(shape),seed=seed,
                 surface_lux=float(surface),surface_standard_error_lux=float(surface_se),
                 clear_surface_lux=float(clear_surface),
                 cloudy_direct_horizontal_lux=float(direct_cloud[1]),
                 clear_surface_monte_carlo_lux=float(clear_light[0,:,1].mean()+direct[0,1]),
                 clear_surface_monte_carlo_standard_error_lux=float(clear_light[0,:,1].std(ddof=1)/4),
                 clear_recess_formal_lux=aperture_base,
                 clear_recess_monte_carlo_lux=float(aperture_clear[0,:,1].mean()),
                 recess_lux=float(recess.mean()),recess_standard_error_lux=float(recess.std(ddof=1)/4),
                 recess_geometry='Black recess with a sky opening 20 degrees wide, elevations 5–25 degrees, centred on the camera; excluded sky/terrain surfaces have zero radiance',
                 sampled_truncated_paths=truncations,seconds=time.monotonic()-start,
                 negative_Y_pixels=int(np.sum(means[:,:,1]<0)),
                 median_pixel_relative_standard_error=float(np.median(errors[:,:,1]/np.maximum(means[:,:,1],1))),
                 raw_clear_monte_carlo_mean_cd_m2=float(clear[:,:,1].mean()),
                 formal_clear_mean_cd_m2=float(background[:,:,1].mean()),
                 source=dict(cloud_section=str(section_path.relative_to(ROOT)),cloud_sha256=digest(section_path),
                             molecular_product_sha256=digest(source),molecular_fields_sha256=digest(path)),
                 producer={f:digest(ROOT/f) for f in ('illumination/cloud_light/render.py','illumination/cloud_light/volume.py',
                           'illumination/cloud_light/microphysics.py','illumination/sky/solved_transport.py')},
                 evidence='Spectral Monte Carlo with molecular and particle multiple scattering, finite-Sun illumination, '
                          'cloud self-shadowing and foreground atmospheric scattering. Equivalent-sphere particle areas '
                          'follow the recorded Morrison size distributions; nonabsorbing HG phase functions and cross-ring '
                          'extrusion are explicit optical and geometric scenarios.',
                 reading_rule='XYZ uses the photopic scale, Y in cd m^-2. Independent photon blocks estimate '
                              'absolute cloudy radiance and sampling uncertainty; paired clear histories check the '
                              'formal clear solution at sea level. The scene uses the common solved molecular column and Lambertian '
                              'ground albedo 0.1. Earthlight, aerosols, terrain beyond spherical ground and ice optical '
                              'features require additional inputs. RGB display mapping belongs to visualization.')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,metadata=json.dumps(records),xyz=means,standard_error_xyz=errors,group_xyz=groups,
                        clear_xyz=background,raw_cloud_xyz=cloudy.reshape((*shape,16,3)),
                        raw_clear_xyz=clear.reshape((*shape,16,3)),elevations_deg=elevations,azimuths_deg=azimuths,
                        contrast=contrast,cloud_extinction_m1=beta,height_edges_m=section['height_edges_m'],
                        column_offset=section['column_offset'],surface_group_xyz=light[0],
                        recess_group_xyz=aperture[0],clear_surface_group_xyz=clear_light[0],
                        clear_recess_group_xyz=aperture_clear[0])
    if truncations:
        raise RuntimeError(f'{truncations} photon paths reached a cap; inspect the saved diagnostic')
    records['archive']=dict(path=str(output.relative_to(ROOT)),sha256=digest(output))
    return records
