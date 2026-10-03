"""Measure regional scene products and their fixed-observer histories.

python -m research.studies.cloud_twilight.appearance
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from climate.crm.cloud_columns import sha256
from illumination.cloud_light.volume import cloud_depth
from illumination.cloud_light.earthlight import context as earthlight_context
from shared.constants import MOON_RADIUS

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
SCHEMA='terluna.research.cloud-evening-appearance/1'
MORRISON_SHA256='8427670b64e6bad03b68d9aa5532f85e685807671b344d6d4480cf6d0db0ab8c'


def patch(groups,clear,mask,weights,clear_groups=None):
    """A deterministic angular patch, with uncertainty across independent blocks."""
    w=weights*mask
    if w.sum()==0:
        return None
    blocks=(groups*w[:,:,None,None]).sum(axis=(0,1))/w.sum()
    background=(clear*w[:,:,None]).sum(axis=(0,1))/w.sum()
    mean=blocks.mean(axis=0);se=blocks.std(axis=0,ddof=1)/np.sqrt(blocks.shape[0])
    contrast_se=se[1]/background[1]
    if clear_groups is not None:
        clear_blocks=(clear_groups*w[:,:,None,None]).sum(axis=(0,1))/w.sum()
        residual=(blocks[:,1]-mean[1]/background[1]*clear_blocks[:,1])/background[1]
        contrast_se=residual.std(ddof=1)/np.sqrt(len(residual))
    return dict(pixels=int(mask.sum()),xyz_cd_m2=mean.tolist(),standard_error_xyz=se.tolist(),
                clear_xyz_cd_m2=background.tolist(),contrast=float(mean[1]/background[1]-1),
                contrast_standard_error=float(contrast_se),
                chromaticity_xy=(mean[:2]/mean.sum()).tolist())


def analyse(scene):
    archive=ROOT/scene['archive']['path']
    if sha256(archive)!=scene['archive']['sha256']:
        raise ValueError('Cloud scene archive differs')
    for file,expected in scene['producer'].items():
        if sha256(ROOT/file)!=expected:
            raise ValueError(f'Cloud scene producer changed: {file}')
    section=ROOT/scene['source']['cloud_section']
    if sha256(section)!=scene['source']['cloud_sha256']:
        raise ValueError('Cloud section differs')
    with np.load(section) as f:
        header=json.loads(str(f['metadata']))
    if header['thermodynamic_source_hashes']['morrison.F']!=MORRISON_SHA256:
        raise ValueError('Microphysical closure requires the inspected Morrison source')
    with np.load(archive) as f:
        groups=f['group_xyz'];clear=f['clear_xyz'];e=f['elevations_deg'];a=f['azimuths_deg']
        clear_groups=f['raw_clear_xyz'] if scene['observer_height_m']!=1.6 else None
        elevations,azimuths=np.radians(np.meshgrid(e,a,indexing='ij'))
        directions=np.stack((np.cos(elevations)*np.cos(azimuths),np.cos(elevations)*np.sin(azimuths),
                             np.sin(elevations)),axis=-1)
        edges=MOON_RADIUS+f['height_edges_m'];beta=f['cloud_extinction_m1']
        dtheta=header['dx_m']/MOON_RADIUS
        theta0=(f['column_offset'][0]-.5)*dtheta
        origin=np.array([0.,0.,MOON_RADIUS+scene['observer_height_m']])
        tau=np.array([cloud_depth(origin,d,edges,theta0,dtheta,beta,scene['cross_ring_width_km']*500)
                      for d in directions.reshape(-1,3)]).reshape(elevations.shape)
    means=groups.mean(axis=2);errors=groups.std(axis=2,ddof=1)/4
    relative_error=float(np.median(errors[:,:,1][means[:,:,1]>0]/means[:,:,1][means[:,:,1]>0]))
    weights=np.cos(elevations)
    camera=scene['camera_azimuth_in_ring_basis_deg']
    lower,upper=scene.get('recess_elevation_limits_deg',(5,25))
    aperture=(e[:,None]>=lower)&(e[:,None]<=upper)&(abs(a[None,:]-camera)<=10)
    adjacent=(e[:,None]>=lower)&(e[:,None]<=upper)&(abs(a[None,:]-camera)>=20)&(abs(a[None,:]-camera)<=40)
    # Predefined around the selected cloud top, retained through the fixed-site history.
    top=scene['initial_pick']['view_elevation_deg']
    feature=(e[:,None]>=max(e[0],top-15))&(e[:,None]<=min(e[-1],top+5))&(abs(a[None,:]-camera)<=10)
    def y_blocks(mask):
        w=weights*mask
        return (groups[:,:,:,1]*w[:,:,None]).sum(axis=(0,1))/w.sum()
    central_blocks=y_blocks(aperture);adjacent_blocks=y_blocks(adjacent)
    ratio=central_blocks.mean()/adjacent_blocks.mean()
    ratio_se=(central_blocks-ratio*adjacent_blocks).std(ddof=1)/(4*adjacent_blocks.mean())
    r=dict(recess_elevation_limits_deg=[lower,upper],case=scene['case'],snapshot=scene['snapshot'],time_day=scene['time_day'],
           latitude_deg=scene['observer_latitude_deg'],longitude_deg=scene['observer_longitude_deg'],
           sun_deg=scene['observer_sun_elevation_deg'],
           hours_after_sunset=scene['initial_pick']['hours_after_sunset']+scene['delta_hours'],
           cloud_top_at_selection_km=scene['initial_pick']['cloud_top_km'],
           width_km=scene['cross_ring_width_km'],ice_g=scene['phase_parameters']['frozen_g'],
           observer_height_m=scene['observer_height_m'],surface_lux=scene['surface_lux'],
           surface_standard_error_lux=scene['surface_standard_error_lux'],
           clear_surface_lux=scene['clear_surface_lux'],recess_lux=scene['recess_lux'],
           recess_standard_error_lux=scene['recess_standard_error_lux'],
           median_pixel_relative_standard_error=relative_error,
           field=patch(groups,clear,np.ones(tau.shape,bool),weights,clear_groups),
           cloud=patch(groups,clear,tau>=.1,weights,clear_groups),
           thin_cloud=patch(groups,clear,(tau>=.01)&(tau<1),weights,clear_groups),
           feature=patch(groups,clear,feature,weights,clear_groups),aperture=patch(groups,clear,aperture,weights,clear_groups),
           adjacent=patch(groups,clear,adjacent,weights,clear_groups),
           aperture_vs_adjacent_contrast=float(ratio-1),aperture_vs_adjacent_standard_error=float(ratio_se),
           cloud_ray_fraction=float(np.mean(tau>=.1)),
           clear_surface_mc_lux=scene['clear_surface_monte_carlo_lux'],
           clear_surface_mc_standard_error_lux=scene['clear_surface_monte_carlo_standard_error_lux'],
           clear_mc_relative_error=float(scene['raw_clear_monte_carlo_mean_cd_m2']/scene['formal_clear_mean_cd_m2']-1),
           source=scene['archive'])
    sky=json.loads((ROOT/'illumination/sky/results/solved_sky.json').read_text())
    r['earthlight_context']=earthlight_context(r['latitude_deg'],r['longitude_deg'],r['hours_after_sunset'],
                                                sky['worlds']['moon_1.2atm']['samples'])
    return r


def main():
    products={};inputs={}
    sky=json.loads((ROOT/'illumination/sky/results/solved_sky.json').read_text())
    for file,expected in sky['producer']['files'].items():
        if sha256(ROOT/file)!=expected:
            raise ValueError(f'Molecular dependency changed: {file}')
    for mode in ('selected','history','sensitivity','deep'):
        path=HERE/'results'/f'evening_scenes_{mode}.json'
        if not path.exists():
            raise FileNotFoundError(f'Complete {mode} scenes first')
        inputs[str(path.relative_to(ROOT))]=sha256(path)
        product=json.loads(path.read_text())
        products[mode]=[analyse(s) for s in product['scenes']]
    for file in ('illumination/cloud_light/earthlight.py','illumination/ephemeris.py','shared/constants.json','illumination/sky/results/solved_sky.json'):
        inputs[file]=sha256(ROOT/file)
    result=dict(schema=SCHEMA,producer=dict(file=str(Path(__file__).relative_to(ROOT)),sha256=sha256(__file__),inputs=inputs),
        evidence='Radiance, colour and contrast from the computed cloud scenes, on sampled second-cycle CM1 fields.',
        reading_rule='Contrast is cloudy radiance divided by the same clear sightline minus one. Angular patches '
                     'are fixed geometrically or by line-of-sight cloud optical depth, independently of sampled brightness. '
                     'Errors are one standard error across 16 independent photon blocks; model and geometry uncertainty '
                     'are separate. The feature patch stays at its initial viewing angles through the fixed-site history. '
                     'The adjacent-view comparison uses the same elevations and offsets 20–40 degrees on either side, '
                     'including any cloud there. The black-recess probe admits 20 degrees of azimuth and the recorded elevation range '
                     '(5–25 degrees in the regional views, 0.5–5 degrees in the deep horizon scene); surface reflection '
                     'inside a real shelter adds light. Elevated clear backgrounds use their own photon estimate.',
        assumptions=dict(admitted_morrison_source_sha256=MORRISON_SHA256),**products)
    target=HERE/'results/evening_appearance.json'
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(target)


if __name__=='__main__':
    main()
