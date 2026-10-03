"""Refine the low-horizon view of an admitted cloud scene and its viewing opening."""
from __future__ import annotations
import json,time
from pathlib import Path
import numpy as np
from .volume import estimate
from illumination.sky.solved_sky import load,digest,horizontal_integral
from illumination.sky.solved_transport import evaluate

ROOT=Path(__file__).resolve().parents[2]


def low_horizon(parent,output,photons=2048,shape=(16,48)):
    if parent['observer_height_m']!=1.6:
        raise ValueError('The detailed clear reference requires the sea-level observer')
    liquid_g=parent['phase_parameters']['liquid_g'];ice_g=parent['phase_parameters']['frozen_g']
    start=time.monotonic();source=ROOT/parent['archive']['path']
    if digest(source)!=parent['archive']['sha256']:
        raise ValueError('Parent cloud scene differs')
    for file,value in parent['producer'].items():
        if digest(ROOT/file)!=value:
            raise ValueError(f'Parent producer differs: {file}')
    path=ROOT/'research/runs/optical_comfort/spherical/moon_1.2atm_standard.npz'
    if digest(path)!=parent['source']['molecular_fields_sha256']:
        raise ValueError('Molecular fields differ')
    air,settings=load(path)
    section_path=ROOT/parent['source']['cloud_section']
    if digest(section_path)!=parent['source']['cloud_sha256']:
        raise ValueError('Cloud section differs')
    with np.load(section_path) as f:
        header=json.loads(str(f['metadata']));sun=f['sun_local']
    with np.load(source) as f:
        beta=f['cloud_extinction_m1'];height=f['height_edges_m'];offset=f['column_offset']
    radius=air['edges'][0];edges=radius+height
    step=header['dx_m']/radius;theta0=(offset[0]-.5)*step
    origin=np.array([0.,0.,radius+parent['observer_height_m']])
    camera=parent['camera_azimuth_in_ring_basis_deg'];width=parent['cross_ring_width_km']*500
    e=np.linspace(.5,15,shape[0]);a=camera+np.linspace(-30,30,shape[1])
    ee,aa=np.radians(np.meshgrid(e,a,indexing='ij'))
    directions=np.stack((np.cos(ee)*np.cos(aa),np.cos(ee)*np.sin(aa),np.sin(ee)),axis=-1).reshape(-1,3)
    args=(origin,directions,sun,air['edges'],air['scattering'],air['absorption'],air['xyz'],edges,theta0,step)
    cloudy,cuts=estimate(*args,beta,width,.1,liquid_g,ice_g,photons=photons,groups=16,seed=9819)
    clear,clear_cuts=estimate(*args,np.zeros_like(beta),width,.1,liquid_g,ice_g,photons=128,groups=16,seed=9819)
    sun_alt=parent['observer_sun_elevation_deg'];sun_az=np.degrees(np.arctan2(sun[1],sun[0]))
    relative=np.abs((a-sun_az+180)%360-180)
    background=evaluate(air,settings,[sun_alt],e,relative)[0]@air['xyz']
    options=dict(photons=65536,groups=16,seed=10819,flux=True,flux_mu_min=np.sin(np.radians(.5)),
                 flux_mu_max=np.sin(np.radians(5)),flux_az_min=np.radians(camera-10),flux_az_max=np.radians(camera+10))
    args=(origin,np.array([[0.,0.,1.]]),sun,air['edges'],air['scattering'],air['absorption'],air['xyz'],edges,theta0,step)
    recess,recess_cuts=estimate(*args,beta,width,.1,liquid_g,ice_g,**options)
    clear_recess,clear_recess_cuts=estimate(*args,np.zeros_like(beta),width,.1,liquid_g,ice_g,**options)
    ae=np.linspace(.5,5,21);az=camera-sun_az+np.linspace(-10,10,41)
    formal_recess=float(horizontal_integral(evaluate(air,settings,[sun_alt],ae,az)@air['xyz'][:,1],ae,az)[0]/2)
    groups=cloudy.reshape((*shape,16,3));means=groups.mean(axis=2);error=groups.std(axis=2,ddof=1)/4
    count=int(cuts.sum()+clear_cuts.sum()+recess_cuts.sum()+clear_recess_cuts.sum())
    record=dict(parent)
    record.update(parent_scene=parent['archive'],shape=list(shape),photons_per_pixel=photons,clear_photons_per_pixel=128,
                  view_elevation_limits_deg=[.5,15.],view_azimuth_offset_limits_deg=[-30.,30.],
                  seed=9819,seconds=time.monotonic()-start,recess_elevation_limits_deg=[.5,5.],
                  recess_lux=float(recess[0,:,1].mean()),recess_standard_error_lux=float(recess[0,:,1].std(ddof=1)/4),
                  recess_geometry='Black recess with a sky opening 20 degrees wide, elevations 0.5–5 degrees, centred on the camera; excluded surrounding directions have zero radiance',
                  clear_recess_formal_lux=formal_recess,clear_recess_monte_carlo_lux=float(clear_recess[0,:,1].mean()),
                  raw_clear_monte_carlo_mean_cd_m2=float(clear[:,:,1].mean()),formal_clear_mean_cd_m2=float(background[:,:,1].mean()),
                  sampled_truncated_paths=count,negative_Y_pixels=int(np.sum(means[:,:,1]<0)),
                  median_pixel_relative_standard_error=float(np.median(error[:,:,1][means[:,:,1]>0]/means[:,:,1][means[:,:,1]>0])),
                  producer={**parent['producer'],str(Path(__file__).relative_to(ROOT)):digest(__file__)},
                  reading_rule=parent['reading_rule']+' This product refines the low-horizon angular field. Open-ground flux is inherited from the admitted parent scene; its independent aperture probe covers elevations 0.5–5 degrees.')
    record.pop('archive',None)
    output=Path(output)
    with np.load(source) as f:
        surface=f['surface_group_xyz'];clear_surface=f['clear_surface_group_xyz']
    np.savez_compressed(output,metadata=json.dumps(record),xyz=means,standard_error_xyz=error,group_xyz=groups,clear_xyz=background,
        raw_cloud_xyz=groups,raw_clear_xyz=clear.reshape((*shape,16,3)),elevations_deg=e,azimuths_deg=a,
        contrast=means[:,:,1]/background[:,:,1]-1,cloud_extinction_m1=beta,height_edges_m=height,column_offset=offset,
        surface_group_xyz=surface,clear_surface_group_xyz=clear_surface,recess_group_xyz=recess[0],clear_recess_group_xyz=clear_recess[0])
    if count:
        raise RuntimeError(f'{count} paths reached a cap in the detailed scene')
    record['archive']=dict(path=str(output.relative_to(ROOT)),sha256=digest(output))
    return record
