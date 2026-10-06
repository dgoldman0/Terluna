"""A conditional Nobili lake night: native terrain, dated geometry and a saved cloud experiment.

This is a scene construction, not a weather forecast. The flat CM1 ring does not
resolve the crater's circulation. Its solar longitude is matched to a real
orbital configuration; its cloud field and local wind remain scenario inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.interpolate import RegularGridInterpolator
from geography.coast_terrain import destination
from geography.topography import Grid, geoid
from research.studies.sea_appearance import lighting
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY, EARTH_RADIUS

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
INPUTS=HERE/'nobili_inputs.json'
SITE=(75.6278,0.)
JD=2463769.879982721
SNAPSHOT,COLUMN=367,381


def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def checked(path,expected):
    path=Path(path)
    if digest(path)!=expected:raise ValueError(f'Input hash mismatch: {path}')
    return path


def terrain(manifest,out,spacing=100.):
    """Join admitted equatorial LOLA rows before interpolation across their seam."""
    pieces=[]
    ppd=manifest['terrain']['pixels_per_degree']
    c0,c1=int(73*ppd),int(79*ppd)
    lons=(np.arange(c0,c1)+.5)/ppd
    latitude=[]
    for key in ['north','south']:
        s=manifest['terrain'][key]
        path=checked(ROOT/s['local_image'],s['sha256'])
        raw=np.memmap(path,dtype='<i2',mode='r',shape=(s['row_count'],s['samples']))
        la=(s['line_projection_offset']-s['first_row']-np.arange(s['row_count']))/ppd
        use=np.flatnonzero(np.abs(la)<1.95)
        pieces.append(np.asarray(raw[use[0]:use[-1]+1,c0:c1],float)*s['scale_m'])
        latitude.extend(la[use])
    latitude=np.array(latitude);radial=np.vstack(pieces)
    f=RegularGridInterpolator((latitude[::-1],lons),radial[::-1],bounds_error=True)
    axis=np.arange(-40000,40000+spacing/2,spacing)
    x,y=np.meshgrid(axis,axis);la,lo=destination(*SITE,x,y)
    h=f(np.stack([la,lo],-1))
    glat=np.arange(la.max()+.05,la.min()-.05,-.05)
    glon=np.arange(lo.min()-.05,lo.max()+.05,.05)
    datum=geoid(Grid(glat,glon,20))
    df=RegularGridInterpolator((glat[::-1],glon),datum[::-1])
    level=json.loads((ROOT/'geography/results/atlas.json').read_text())['sea_level_m']
    h-=df(np.stack([la,lo],-1))+level
    metadata=dict(schema='terluna.research.nobili-terrain/1',centre_lon_lat=SITE,spacing_m=spacing,
        sea_level_above_geoid_m=level,evidence='NASA LOLA 256 pixels/degree, sampled onto a local metre grid above the degree-200 GRAIL geoid and the selected atlas water level.',
        reading_rule='Heights above nominal atlas water; no monthly lake-level evolution, erosion or water-load subsidence. Native sampling is about 118 m, including altimetric interpolation.',
        inputs={str(ROOT/s['local_image']):s['sha256'] for s in [manifest['terrain']['north'],manifest['terrain']['south']]})
    np.savez_compressed(out/'terrain.npz',x_m=axis,y_m=axis,height_m=h.astype('f4'),metadata=json.dumps(metadata))
    distance=np.arange(50,40000,50.);az=np.arange(0.,360.,.25)
    angles=[];ranges=[];heights=[];first=[]
    for a in az:
        xx=distance*np.sin(np.radians(a));yy=distance*np.cos(np.radians(a))
        z=map_coordinates(h,[(yy-axis[0])/spacing,(xx-axis[0])/spacing],order=1,mode='nearest')
        q=distance/MOON_RADIUS
        elevation=np.degrees(np.arctan2((MOON_RADIUS+np.maximum(z,0))*np.cos(q)-(MOON_RADIUS+2),(MOON_RADIUS+np.maximum(z,0))*np.sin(q)))
        k=int(np.argmax(elevation));land=np.flatnonzero(z>0)
        angles.append(float(elevation[k]));ranges.append(float(distance[k]));heights.append(float(z[k]));first.append(float(distance[land[0]]) if len(land) else None)
    profile=dict(azimuth_deg=az.tolist(),elevation_deg=angles,range_m=ranges,height_m=heights,first_land_m=first)
    return dict(file=str((out/'terrain.npz').relative_to(ROOT)),sha256=digest(out/'terrain.npz'),
        water_depth_at_observer_m=-float(h[len(axis)//2,len(axis)//2]),profile=profile,metadata=metadata)


def horizons(manifest):
    row=manifest['horizons'];path=checked(ROOT/row['local_file'],row['sha256'])
    line=path.read_text().split('$$SOE')[1].split('$$EOE')[0].strip()
    v=[x.strip() for x in line.split(',')]
    return dict(time_tt=v[0],jd_tt=JD,azimuth_deg=float(v[5]),elevation_deg=float(v[6]),
        lit_fraction=float(v[7])/100,diameter_deg=float(v[8])/3600,
        subobserver_lon_east_deg=float(v[9]),subobserver_lat_deg=float(v[10]),
        subsolar_lon_east_deg=float(v[11]),subsolar_lat_deg=float(v[12]),
        phase_angle_deg=float(v[17]),north_pole_ra_deg=float(v[18]),north_pole_dec_deg=float(v[19]))


def build(out,crm_case):
    from climate.crm import cloud_scene
    from climate.crm.ring_analysis import read_snapshot,case_geometry
    from illumination.water_surface import short_waves
    from illumination.water_surface.full_spectrum import full_curvature
    manifest=json.loads(INPUTS.read_text());out.mkdir(parents=True,exist_ok=True)
    topo=terrain(manifest,out);print('Native terrain prepared',flush=True)
    h=horizons(manifest)
    g=lighting.geometry(np.array([JD]),*SITE)
    sky=lighting.Sky();earth=lighting.Earthlight(sky)
    distance=EARTH_RADIUS/np.sin(np.radians(h['diameter_deg']/2))
    weights=earth.weights(np.array([h['phase_angle_deg']]),np.array([distance]),g['earth_sunlight_factor'])[0]
    sun=float(g['sun_elevation'][0]);e=h['elevation_deg']
    from research.studies.sea_appearance.regimes import normal_beam
    beam=normal_beam(sky,e)*weights
    ground_earth=float((sky.ground_channels([e])[0]*weights)@sky.channel_y)
    light=dict(sun_lux=float(sky.ground(sun)),earth_lux=ground_earth,placeholder_lux=lighting.FLOOR_LUX,
        total_clear_sky_lux=float(sky.ground(sun))+ground_earth+lighting.FLOOR_LUX,
        earth_direct_normal_lux=float(beam@sky.channel_y),earth_above_air_normal_lux=float(weights@sky.channel_y),
        sun_elevation_deg=sun,sun_azimuth_deg=float(g['sun_azimuth'][0]),
        earth_sunlight_factor=float(g['earth_sunlight_factor'][0]),earth_distance_m=float(distance),
        scope='Clear molecular atmosphere. The cloud experiment is a separate conditional calculation; this is not cloudy ground illumination.')
    source=manifest['cloud'];case=Path(crm_case)
    checked(case/'case.json',source['case_sha256']);checked(case/source['file'],source['sha256'])
    old=cloud_scene.RUNS
    try:
        cloud_scene.RUNS=case.parent
        cloud=cloud_scene.export(case.name,SNAPSHOT,COLUMN,source['sha256'],out/'cloud-section.npz',half_columns=80)
    finally:cloud_scene.RUNS=old
    d=read_snapshot(case,SNAPSHOT);u10=float(d['s10'][COLUMN]);ust=float(d['ust'][COLUMN]);u=float(d['u10'][COLUMN]);v=float(d['v10'][COLUMN])
    wind=dict(speed_10m_m_s=u10,friction_velocity_m_s=ust,east_m_s=u,north_m_s=v,toward_cartesian_deg=float(np.degrees(np.arctan2(v,u))%360))
    # The existing unified spectrum supplies conditional wave-age sensitivities.
    # Neither is a terrain-resolving wave forecast for this lake.
    waves=[]
    for omega in [.84,1.5]:
        kp=MOON_SURFACE_GRAVITY*omega**2/u10**2;km,_=short_waves.capillary_scales(MOON_SURFACE_GRAVITY)
        k=np.geomspace(kp/30,15*km,12000);long,short=full_curvature(k,u10,ust,omega,MOON_SURFACE_GRAVITY)
        variance=float(np.trapezoid((long+short)/k**2,np.log(k)))
        waves.append(dict(inverse_wave_age=omega,significant_height_m=4*np.sqrt(variance),peak_wavelength_m=2*np.pi/kp))
    # Match the experiment's Sun longitude, not an unrelated clock time.
    sg=cloud['metadata'];model_sun_lon=(SITE[0]-sg['observer_hour_angle_deg'])%360
    actual=__import__('geography.lunar_ephemeris',fromlist=['earth_and_sun']).earth_and_sun(np.array([JD]))
    actual_lon=float(np.degrees(np.arctan2(actual['sun'][0,1],actual['sun'][0,0]))%360)
    record=dict(schema='terluna.research.nobili-night/1',jd_tt=JD,site_lon_lat_deg=SITE,
        camera=dict(frame=[2560,1920],horizontal_fov_deg=50.,view_azimuth_deg=270.,pitch_deg=3.,eye_height_m=2.,
            lon_lat_deg=list(SITE),note='A water-level offshore viewpoint; no boat is visible. No solid standing platform is claimed.'),
        earth=h,light=light,terrain=topo,cloud=dict(file=str((out/'cloud-section.npz').relative_to(ROOT)),sha256=digest(out/'cloud-section.npz'),
            case=case.name,snapshot=SNAPSHOT,column=COLUMN,time_day=sg['time_day'],model_sun_longitude_deg=model_sun_lon,
            ephemeris_sun_longitude_deg=actual_lon,cross_ring_width_m=20000.,
            evidence='Saved CM1 2-D flat equatorial ring; 20-km cross-ring extrusion is an explicit scene assumption. No terrain-coupled weather or dated forecast.'),
        wind=wind,wave_sensitivity=waves,selected_inverse_wave_age=.84,
        evidence='Native terrain and a dated orbital configuration, with a solar-phase-matched experimental cloud/wind section. A physically constrained illustration scenario, not an integrated forecast.',
        limits=['The nominal atlas water level floods Nobili; local monthly tides and shoreline evolution are not solved.',
            'Cloud and wind fields do not resolve this crater or cross-ring weather. Cloud placement remains conditional.',
            'Lake waves use an equilibrium unified spectrum at the saved wind and stress, with wave-age sensitivity; no lake wave-history or breaking simulation.',
            'Historical NASA Earth maps illustrate surface/cloud structure; they do not predict 2033 weather or full spatial spectra.',
            'Land reflectance and sub-118-m detail are placeholders; mature vegetation is not established by this optical scene.'],
        producer=dict(path=str(Path(__file__).relative_to(ROOT)),sha256=digest(__file__)),
        inputs={str(p.relative_to(ROOT)):digest(p) for p in [INPUTS,ROOT/'geography/topography.py',ROOT/'geography/coast_terrain.py',ROOT/'geography/results/atlas.json',ROOT/'shared/constants.json']})
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--crm-case',type=Path,required=True)
    p.add_argument('--build-dir',type=Path,default=ROOT/'research/runs/sea_appearance/nobili')
    p.add_argument('--out',type=Path,default=HERE/'results/nobili-night.json');a=p.parse_args()
    r=build(a.build_dir,a.crm_case);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['earth','light','wind','wave_sensitivity']},indent=2))

if __name__=='__main__':main()
