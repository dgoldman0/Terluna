"""Unlabelled WP-6 composition and wave guides from the sea study's products.

These are display references, not the spectral path tracer. The mean guide uses
stored XYZ colours and the domain's uniform-disk glitter model. The wave guide
projects a reduced domain realization and uses approximate reflected-sky shading.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, map_coordinates


def xyz(point):
    x, y = point['chromaticity_xy']
    return np.array([x/y, 1, (1-x-y)/y]) * point['luminance_cd_m2']


def interpolate(points, values, width, height):
    columns, xs = [], []
    for side in ('left edge', 'centre', 'right edge'):
        selected = sorted((i for i,p in enumerate(points) if p['label'].endswith(side)),
                          key=lambda i:points[i]['y'])
        xs.append(points[selected[0]]['x'])
        ys = [points[i]['y'] for i in selected]
        columns.append(np.stack([np.interp(np.arange(height), ys, values[selected,c]) for c in range(3)], -1))
    result = np.empty((height,width,3))
    for y in range(height):
        for c in range(3):
            result[y,:,c] = np.interp(np.arange(width), xs, [a[y,c] for a in columns])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--water-structure', action='store_true')
    args = parser.parse_args()
    root = (args.root or Path(__file__).resolve().parents[2]).resolve()
    sys.path.insert(0,str(root))
    sys.path.insert(0,str(root/'visualization/sea-appearance'))
    from image_guides import digest, rgb, linear, encoded
    from research.studies.sea_appearance.scenes import Frame, unit, XYZ_TO_SRGB, SHOULDER
    from illumination.water_surface import realization, reflection
    from illumination.water_surface.hotfile import read_hotfile
    from illumination.water_surface.model import wavenumber
    from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY

    source = root/'research/studies/sea_appearance/results/scenes.json'
    scene = next(s for s in json.loads(source.read_text())['scenes']
                 if s['coast']=='W Procellarum' and s['moment']=='darkest')
    out = args.out or root/'visualization/sea-appearance/results/illustrations'
    out.mkdir(parents=True,exist_ok=True)
    width,height = scene['camera']['frame']
    horizon = scene['camera']['horizon_row_at_centre']
    frame = Frame(scene['camera']['view_azimuth_deg'],scene['camera']['pitch_deg'])
    eye,white = scene['camera']['eye_height_m'],scene['light']['white_cd_m2']
    earth = scene['earth']
    earth_direction = unit(earth['elevation_deg'],earth['azimuth_deg'])
    radius = math.radians(earth['diameter_deg']/2)
    covariance = np.array(scene['waves']['slope_covariance_east_north'])
    index = float(reflection.refractive_index(550.0))
    x,y = earth['colour_chromaticity']
    earth_xyz = np.array([x/y,1,(1-x-y)/y])*earth['normal_illuminance_lux']

    def display(values):
        relative = np.maximum(values[...,1]/white,0)
        rolled = np.where(relative<=SHOULDER, relative,
                          SHOULDER+(1-SHOULDER)*(1-np.exp(-np.maximum(relative-SHOULDER,0)/(1-SHOULDER))))
        scale = np.divide(rolled,relative,out=np.ones_like(relative),where=relative>0)
        return np.clip(values/white*scale[...,None] @ XYZ_TO_SRGB.T,0,1)

    # Subtract the domain glitter at the stored samples before interpolating the
    # background. Re-add it over the full frame, retaining the narrow physical path.
    sea_points = scene['sea']
    sea_xyz = np.array([xyz(p) for p in sea_points])
    views = -frame.direction(np.array([p['x'] for p in sea_points]),np.array([p['y'] for p in sea_points]))
    kernels = reflection.disk_glint(earth_direction,radius,views,covariance,index)
    background_samples = np.maximum(sea_xyz-kernels[:,None]*earth_xyz,0)
    background_xyz = interpolate(sea_points,background_samples,width,height)
    sky_xyz = interpolate(scene['sky'],np.array([xyz(p) for p in scene['sky']]),width,height)
    sky = display(sky_xyz)
    average = background_xyz.copy()
    for start in range(math.ceil(horizon),height,16):
        yy,xx = np.mgrid[start:min(start+16,height),0:width]
        views = -frame.direction(xx.ravel()+.5,yy.ravel()+.5)
        kernel = reflection.disk_glint(earth_direction,radius,views,covariance,index)
        average[start:start+len(yy)] += kernel.reshape(yy.shape)[...,None]*earth_xyz
    sea = display(average)
    pixels = sky.copy()
    pixels[math.ceil(horizon):] = sea[math.ceil(horizon):]
    region = np.zeros((height,width),np.uint8)
    region[math.ceil(horizon):] = 2

    terrain = root/'research/runs/sea_appearance/terrain/w_procellarum.npz'
    with np.load(terrain,allow_pickle=False) as data:
        z = data['height_m'].astype(float)-scene['tide_m']
        axis = data['x_m']
    spacing = float(axis[1]-axis[0])
    distances = np.arange(25,72000,25.0)
    skyline = np.full(width,horizon)
    land = linear(rgb(scene['land'][0]['apparent_srgb']))
    for col in range(width):
        direction = frame.direction(np.array([col+.5]),np.array([horizon]))[0]
        azimuth = math.atan2(direction[0],direction[1])
        east,north = distances*math.sin(azimuth),distances*math.cos(azimuth)
        heights = map_coordinates(z,[(north-axis[0])/spacing,(east-axis[0])/spacing],order=1,mode='nearest')
        elev = np.degrees(np.arctan2(heights-eye-distances**2/(2*MOON_RADIUS),distances))
        _,rows = frame.project(unit(elev,np.full_like(elev,math.degrees(azimuth))))
        top = math.ceil(horizon)
        for j in np.flatnonzero(heights>0):
            row = max(0,int(math.floor(rows[j])))
            if row<top:
                pixels[row:top,col] = land
                region[row:top,col] = 1
                top = row
        skyline[col] = top

    # The phase is orthographic sphere geometry, using the stored projected
    # illumination direction. Bloom is only a small illustrative optical finish.
    yy,xx = np.mgrid[0:height,0:width]
    dx=(xx+.5-earth['x'])/(earth['width_px']/2)
    dy=(earth['y']-yy-.5)/(earth['height_px']/2)
    radial=dx*dx+dy*dy
    a,b=math.radians(earth['phase_angle_deg']),math.radians(earth['lit_side_toward_deg'])
    illumination=dx*math.sin(a)*math.cos(b)+dy*math.sin(a)*math.sin(b)+np.sqrt(np.maximum(0,1-radial))*math.cos(a)
    lit=(radial<=1)&(illumination>0)
    bloom=gaussian_filter(lit.astype(float),2.0)*.18
    pixels=np.minimum(pixels+bloom[...,None]*linear(rgb('#fff8dc')),1)
    pixels[lit]=linear(rgb('#fffef3'))
    region[lit]=3

    def save(name,array):
        Image.fromarray(np.uint8(np.clip(np.rint(encoded(array)*255),0,255))).save(out/name)

    save('wp-6-calculated-guide.png',pixels)
    Image.fromarray(np.where(region==1,255,0).astype('uint8')).save(out/'wp-6-land-mask.png')
    Image.fromarray(np.where(lit,255,0).astype('uint8')).save(out/'wp-6-earth-phase-mask.png')
    errors=[{'x':p['x'],'difference_px':round(float(skyline[p['x']]-(horizon-p['px_above_horizon'])),2)}
            for p in scene['land_profile'] if p['px_above_horizon'] is not None]
    record={
        'schema':'terluna.visualization.sea-image-guide/1','scene_id':'WP-6',
        'producer':{'script':'visualization/sea-appearance/wp6_guides.py','sha256':digest(__file__)},
        'inputs':{str(p.relative_to(root)):digest(p) for p in [source,terrain]},
        'source_files':{str(p.relative_to(root)):digest(p) for p in [root/'visualization/sea-appearance/image_guides.py',root/'research/studies/sea_appearance/scenes.py',root/'illumination/water_surface/reflection.py',root/'shared/constants.json']},
        'evidence':'Projected LOLA coast and stored scene XYZ colours. Domain uniform-disk reflection reconstructs the mean glitter path. Disk phase uses the stored illumination angle; no textured Earth map is introduced.',
        'reading_rule':'An image-generation guide, not the full radiance renderer. Background colours interpolate sparse samples; land appearance is constant, shadows are not solved, and the small display bloom is illustrative.',
        'native_size':[width,height],'horizon_y':horizon,
        'earth':earth,'earth_mask_lit_fraction':float(lit.sum()/np.sum(radial<=1)),
        'maximum_skyline_difference_px':max(abs(e['difference_px']) for e in errors),
        'skyline_comparison':errors,
        'region_fractions':{str(i):float(np.mean(region==i)) for i in range(4)},
    }
    (out/'wp-6-calculated-guide.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['earth_mask_lit_fraction','maximum_skyline_difference_px','region_fractions']},indent=2),flush=True)

    if not args.water_structure:
        return
    waves=scene['waves']
    wave_source=root/waves['file']
    raw=read_hotfile(wave_source,nodes=[scene['camera']['lon_lat_deg']])
    f,d,e=raw['frequency'],raw['direction'],raw['variance'][0]
    density=realization.swan_density(f,d,e,waves['depth_m'],MOON_SURFACE_GRAVITY)
    seed=20380211
    tiles=realization.realize([density],waves['depth_m'],MOON_SURFACE_GRAVITY,seed,cascades=((512.0,512),(31.7,512)))
    ocean=pixels.copy()
    distance=np.geomspace(2,3000,6000)
    for col in range(width):
        direction=frame.direction(np.array([col+.5]),np.array([horizon]))[0]
        azimuth=math.atan2(direction[0],direction[1])
        east,north=distance*math.sin(azimuth),distance*math.cos(azimuth)
        z,sx,sy=np.zeros_like(distance),np.zeros_like(distance),np.zeros_like(distance)
        for tile in tiles:
            coords=[north/tile.spacing_m,east/tile.spacing_m]
            if tile.height is not None:
                z+=map_coordinates(tile.height,coords,order=1,mode='grid-wrap')
            sx+=map_coordinates(tile.slope[0],coords,order=1,mode='grid-wrap')
            sy+=map_coordinates(tile.slope[1],coords,order=1,mode='grid-wrap')
        elev=np.degrees(np.arctan2(z-eye-distance**2/(2*MOON_RADIUS),distance))
        _,rows=frame.project(unit(elev,np.full_like(elev,math.degrees(azimuth))))
        v=np.stack([-east,-north,eye-z],-1)
        v/=np.linalg.norm(v,axis=1,keepdims=True)
        normal=np.stack([-sx,-sy,np.ones_like(sx)],-1)
        normal/=np.linalg.norm(normal,axis=1,keepdims=True)
        cos_v=np.maximum(np.sum(v*normal,axis=1),0)
        reflected=2*cos_v[:,None]*normal-v
        dot=reflected@earth_direction
        hit=dot>=math.cos(radius)
        rho=reflection.fresnel(cos_v,index)
        contrast=np.clip(rho/np.maximum(reflection.fresnel(v[:,2],index),1e-6),.4,1.8)
        top=height
        for j in range(len(distance)):
            row=max(math.ceil(horizon),int(math.floor(rows[j])))
            if row<top:
                values=background_xyz[row:top,col]*contrast[j]
                if hit[j]:
                    values=values+rho[j]*(earth_xyz/earth['normal_illuminance_lux'])*earth['disk_luminance_cd_m2']
                ocean[row:top,col]=display(values)
                top=row
        ocean[region[:,col]!=2,col]=pixels[region[:,col]!=2,col]
    save('wp-6-wave-structure-guide.png',ocean)
    kp=float(wavenumber(np.array([f[np.argmax(e.sum(axis=1))]]),waves['depth_m'],MOON_SURFACE_GRAVITY)[0])
    record={
        'schema':'terluna.visualization.sea-wave-guide/1','scene_id':'WP-6',
        'producer_sha256':digest(__file__),
        'inputs':{str(wave_source.relative_to(root)):digest(wave_source),str(source.relative_to(root)):digest(source)},
        'seed':seed,'cascades':[[t.size_m,t.cells] for t in tiles],
        'source_hs_m':waves['significant_height_m'],
        'guide_hs_m':float(4*np.sqrt(sum(np.var(t.height) for t in tiles if t.height is not None))),
        'source_peak_wavelength_m':waves['peak_wavelength_m'],'guide_peak_wavelength_m':float(2*np.pi/kp),
        'source_mss':waves['mean_square_slope']['total'],'guide_mss':float(sum(np.trace(realization.sample_covariance(t)) for t in tiles)),
        'evidence':'Domain realization of the stored local W Procellarum SWAN spectrum, at lunar gravity, with no short wind-wave addition. Projected wave heights and normals with a uniform-disk mirror-hit reference.',
        'reading_rule':'Use for broad wave scale and glassy reflection structure. It is a single instantaneous realization; its bright facets need not match the statistical glitter guide.',
        'omissions':['spectral reflected-sky transport','terrain reflections','light-source visibility behind terrain','facet shadowing','phase-dependent Earth reflectance','height variations of the short-wave tile','wave visibility above the nominal flat horizon','time integration and camera motion of a long exposure'],
    }
    (out/'wp-6-wave-structure-guide.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['source_hs_m','guide_hs_m','source_peak_wavelength_m','guide_peak_wavelength_m','source_mss','guide_mss']},indent=2),flush=True)


if __name__=='__main__':
    main()
