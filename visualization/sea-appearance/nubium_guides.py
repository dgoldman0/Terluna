"""Portrait and textured-Earth display guides for the supplemental Nubium scene.

The scene product owns the physics. This renderer projects its colours, terrain
and a reduced domain wave realization. NASA maps and a saved Horizons observer
record supply an illustrative Earth texture, separately exposed to retain detail.
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


class Camera:
    """Rectilinear display projection, in the study's east/north/up frame."""
    def __init__(self, spec):
        from research.studies.sea_appearance.scenes import unit
        self.width, self.height = spec['frame']
        self.f = unit(spec['pitch_deg'], spec['view_azimuth_deg'])
        a = math.radians(spec['view_azimuth_deg'])
        self.r = np.array([math.cos(a), -math.sin(a), 0.0])
        self.u = np.cross(self.r, self.f)
        self.focal = self.width/(2*math.tan(math.radians(spec['horizontal_fov_deg'])/2))

    def direction(self, x, y):
        d = self.f + ((x-self.width/2)/self.focal)[:,None]*self.r + ((self.height/2-y)/self.focal)[:,None]*self.u
        return d/np.linalg.norm(d,axis=1,keepdims=True)

    def project(self, d):
        d = np.atleast_2d(d)
        return self.width/2+self.focal*(d@self.r)/(d@self.f), self.height/2-self.focal*(d@self.u)/(d@self.f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path)
    parser.add_argument('--scene', type=Path)
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--waves', action='store_true')
    args = parser.parse_args()
    root = (args.root or Path(__file__).resolve().parents[2]).resolve()
    sys.path.insert(0,str(root))
    sys.path.insert(0,str(root/'visualization/sea-appearance'))
    from image_guides import digest, rgb, linear, encoded
    from wp6_guides import interpolate, xyz
    from research.studies.sea_appearance import scenes, lighting
    from illumination.stars import sky as stars
    from illumination.water_surface import reflection, realization
    from illumination.water_surface.hotfile import read_hotfile
    from illumination.water_surface.model import wavenumber
    from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY

    source = args.scene or root/'research/studies/sea_appearance/results/nubium-midnight.json'
    product = json.loads(source.read_text())
    s = product['scene']; earth = s['earth']; cam = Camera(s['camera'])
    manifest_path = root/'visualization/sea-appearance/nubium-earth-inputs.json'
    manifest = json.loads(manifest_path.read_text())
    for entry in manifest['files']:
        if digest(args.assets/entry['filename']) != entry['sha256']:
            raise ValueError('Source bytes differ: '+entry['filename'])
    h = manifest['horizons']
    width,height = s['camera']['frame']; horizon = s['camera']['horizon_row_at_centre']
    out = args.out; out.mkdir(parents=True,exist_ok=True)
    # A chosen photographic exposure, not a model of mesopic colour perception.
    exposure = 0.2
    white = s['light']['white_cd_m2']/exposure
    def display(v):
        rel = np.maximum(v[...,1]/white,0)
        rolled = np.where(rel<=scenes.SHOULDER,rel,scenes.SHOULDER+(1-scenes.SHOULDER)*(1-np.exp(-np.maximum(rel-scenes.SHOULDER,0)/(1-scenes.SHOULDER))))
        return np.clip(v/white*np.divide(rolled,rel,out=np.ones_like(rel),where=rel>0)[...,None]@scenes.XYZ_TO_SRGB.T,0,1)
    def save(name,v):
        Image.fromarray(np.uint8(np.rint(np.clip(encoded(v),0,1)*255))).save(out/name)

    ed = scenes.unit(earth['elevation_deg'],earth['azimuth_deg'])
    radius = math.radians(earth['diameter_deg']/2)
    covariance = np.array(s['waves']['slope_covariance_east_north'])
    index = float(reflection.refractive_index(550))
    ex,ey = earth['colour_chromaticity']
    earth_xyz = np.array([ex/ey,1,(1-ex-ey)/ey])*earth['normal_illuminance_lux']
    sky_xyz = interpolate(s['sky'],np.array([xyz(p) for p in s['sky']]),width,height)
    sp = s['sea']
    sv = -cam.direction(np.array([p['x'] for p in sp]),np.array([p['y'] for p in sp]))
    kernels = reflection.disk_glint(ed,radius,sv,covariance,index)
    background = interpolate(sp,np.maximum(np.array([xyz(p) for p in sp])-kernels[:,None]*earth_xyz,0),width,height)
    sea_xyz = background.copy()
    for start in range(math.ceil(horizon),height,16):
        yy,xx = np.mgrid[start:min(start+16,height),:width]
        k = reflection.disk_glint(ed,radius,-cam.direction(xx.ravel()+.5,yy.ravel()+.5),covariance,index)
        sea_xyz[start:start+len(yy)] += k.reshape(yy.shape)[...,None]*earth_xyz
    pixels = display(sky_xyz)
    pixels[math.ceil(horizon):] = display(sea_xyz[math.ceil(horizon):])
    region = np.zeros((height,width),np.uint8); region[math.ceil(horizon):] = 2
    terrain = root/'research/runs/sea_appearance/terrain/s_nubium.npz'
    with np.load(terrain,allow_pickle=False) as data:
        z = data['height_m'].astype(float)-s['tide_m']; axis = data['x_m']
    distance = np.arange(25,105000,25.0); spacing=float(axis[1]-axis[0])
    eye=s['camera']['eye_height_m']; cx,cy=scenes.CAMERAS['S Nubium']['offset_m']
    skyline=np.full(width,horizon)
    land_rgb=linear(rgb(s['land'][0]['apparent_srgb']))*exposure
    for col in range(width):
        ray=cam.direction(np.array([col+.5]),np.array([horizon]))[0]
        az=math.atan2(ray[0],ray[1])
        east,north=cx+distance*math.sin(az),cy+distance*math.cos(az)
        inside=(east>=axis[0])&(east<=axis[-1])&(north>=axis[0])&(north<=axis[-1])
        heights=map_coordinates(z,[(north-axis[0])/spacing,(east-axis[0])/spacing],order=1,mode='nearest')
        elev=np.degrees(np.arctan2(heights-eye-distance**2/(2*MOON_RADIUS),distance))
        _,rows=cam.project(scenes.unit(elev,np.full_like(elev,math.degrees(az))))
        valid=inside&(heights>0)
        if np.any(valid):
            top=max(0,min(math.ceil(horizon),int(math.floor(np.min(rows[valid])))))
            pixels[top:math.ceil(horizon),col]=land_rgb
            region[top:math.ceil(horizon),col]=1; skyline[col]=top

    # Horizons gives the dated Earth-facing longitude/latitude. The domain's
    # existing star transform places its ICRF north pole in the lunar local frame.
    pole_body=stars.body_directions(np.radians(h['north_pole_ra_deg']),np.radians(h['north_pole_dec_deg']),product['jd_tt'])[0]
    east,north,up=lighting.local_frame(*s['camera']['lon_lat_deg'])
    pole=np.array([pole_body@east,pole_body@north,pole_body@up])
    tangent_r=cam.r-(cam.r@ed)*ed; tangent_r/=np.linalg.norm(tangent_r)
    tangent_u=np.cross(tangent_r,ed)
    pole_t=pole-(pole@ed)*ed; pole_t/=np.linalg.norm(pole_t)
    nr,nu=pole_t@tangent_r,pole_t@tangent_u
    lon,lat=np.radians([h['subobserver_lon_east_deg'],h['subobserver_lat_deg']])
    centre=np.array([math.cos(lat)*math.cos(lon),math.cos(lat)*math.sin(lon),math.sin(lat)])
    geographic_north=np.array([-math.sin(lat)*math.cos(lon),-math.sin(lat)*math.sin(lon),math.cos(lat)])
    geographic_east=np.array([-math.sin(lon),math.cos(lon),0])
    rb=nu*geographic_east+nr*geographic_north
    ub=-nr*geographic_east+nu*geographic_north
    slon,slat=np.radians([h['subsolar_lon_east_deg'],h['subsolar_lat_deg']])
    sun=np.array([math.cos(slat)*math.cos(slon),math.cos(slat)*math.sin(slon),math.sin(slat)])
    # Preserve the NASA maps' geography; clouds are a historical illustrative
    # prior. Opacity, Lambert shading and the RGB atmospheric transfer are display
    # approximations, not spatially resolved spectral Earth radiance.
    surface=linear(np.asarray(Image.open(args.assets/'nubium-earth-bluemarble-2048.png').convert('RGB'))/255)
    cloud=np.asarray(Image.open(args.assets/'nubium-earth-clouds-2048.jpg').convert('L'))/255
    inputs=scenes.Inputs(); g=lighting.geometry(np.array([product['jd_tt']]),*s['camera']['lon_lat_deg'])
    weights=inputs.earth_light.weights(g['earth_phase_angle'],g['earth_distance_m'],g['earth_sunlight_factor'])[0]
    above=weights@inputs.xyz; below=(inputs.beam(earth['elevation_deg'])*weights)@inputs.xyz
    transfer=(below/below[1]@scenes.XYZ_TO_SRGB.T)/(above/above[1]@scenes.XYZ_TO_SRGB.T)
    transfer/=max(transfer)
    def globe(dx,dy):
        rr=dx*dx+dy*dy; inside=rr<=1
        zz=np.sqrt(np.maximum(0,1-rr))
        normals=dx[...,None]*rb+dy[...,None]*ub+zz[...,None]*centre
        glon=np.arctan2(normals[...,1],normals[...,0]); glat=np.arcsin(np.clip(normals[...,2],-1,1))
        tx=(glon+math.pi)/(2*math.pi)*surface.shape[1]-.5
        ty=(math.pi/2-glat)/math.pi*surface.shape[0]-.5
        ground=np.stack([map_coordinates(surface[...,c],[ty,tx],order=1,mode='grid-wrap') for c in range(3)],-1)
        clouds=map_coordinates(cloud,[ty,tx],order=1,mode='nearest')[...,None]*.96
        lit=np.maximum(normals@sun,0)
        texture=(ground*(1-clouds)+clouds*.9)*lit[...,None]*transfer
        # One fixed shorter exposure for the disk, normalized to diffuse cloud
        # white. Surrounding sky keeps the scene exposure: a deliberate HDR blend.
        texture=np.clip(texture*1.08,0,.98)
        return texture,inside,lit>0

    # Evaluate only the disk's small bounding box at scene resolution.
    x0=max(0,int(earth['x']-earth['width_px']/2)-2); x1=min(width,int(earth['x']+earth['width_px']/2)+3)
    y0=max(0,int(earth['y']-earth['height_px']/2)-2); y1=min(height,int(earth['y']+earth['height_px']/2)+3)
    yy,xx=np.mgrid[y0:y1,x0:x1]
    # Local orthographic sphere, stretched by the rectilinear projection at its
    # off-axis position. Finite-distance perspective and oblateness are omitted.
    disk,inside,lit=globe((xx+.5-earth['x'])/(earth['width_px']/2),(earth['y']-yy-.5)/(earth['height_px']/2))
    pixels[y0:y1,x0:x1][inside]=disk[inside]
    full_lit=np.zeros((height,width),bool); full_lit[y0:y1,x0:x1]=inside&lit
    bloom=gaussian_filter(full_lit.astype(float),1.25)*.045
    pixels=np.clip(pixels+bloom[...,None]*np.array([1,.8,.4]),0,1)
    region[y0:y1,x0:x1][inside]=3
    save('nubium-midnight-calculated-guide.png',pixels)
    Image.fromarray(region).save(out/'nubium-midnight-regions.png')
    Image.fromarray(np.uint8(full_lit*255)).save(out/'nubium-midnight-earth-mask.png')
    # Separate enlarged component, with the same map orientation and filter.
    yy,xx=np.mgrid[:1024,:1024]
    globe_image,mask,_=globe((xx+.5-512)/370,(512-yy-.5)/370)
    close=np.broadcast_to(display(sky_xyz[184,601]),(1024,1024,3)).copy()
    close[mask]=globe_image[mask]
    save('nubium-midnight-earth-detail.png',close)
    # Display projection checks, independent of the generated image.
    px=np.array([0,width/2,width-1.0]); py=np.array([0,height/2,height-1.0])
    qx,qy=cam.project(cam.direction(px,py))
    record=dict(schema='terluna.visualization.nubium-guides/1',
        producer={'script':'visualization/sea-appearance/nubium_guides.py','sha256':digest(__file__)},
        inputs={str(p):digest(p) for p in [source,terrain,manifest_path,root/'visualization/sea-appearance/wp6_guides.py',root/'visualization/sea-appearance/image_guides.py',root/'illumination/stars/sky.py',root/'illumination/water_surface/realization.py',root/'shared/constants.json']},
        evidence='Scene XYZ readouts and statistical Earth glitter; tide-adjusted LOLA silhouette; NASA map projection oriented by JPL Horizons and the existing stellar transform.',
        reading_rule='An image-generation reference, not a radiance validation. Terrain colour is constant and shadowing unsolved. Sparse sky/sea colours are interpolated. Earth has a separate display exposure, an RGB atmospheric tint and approximate Lambert shading; cloud maps are historical. No prediction of 2038 cloud cover.',
        frame=[width,height],scene_exposure_multiplier=exposure,
        exposure_stops_from_scene=float(math.log2(exposure)),earth_disk_short_exposure=True,
        earth_rgb_transfer=transfer.tolist(),earth_north_screen_angle_deg=float(np.degrees(np.arctan2(nu,nr))),
        earth_north_pole_subobserver_lat_check_deg=float(np.degrees(np.arcsin(-pole@ed))),
        horizons_geometry_differences={'elevation_deg':earth['elevation_deg']-h['elevation_deg'],
             'azimuth_deg':earth['azimuth_deg']-h['azimuth_deg'],
             'phase_fraction':earth['lit_fraction']-h['lit_fraction']},
        projection_roundtrip_max_px=float(max(np.max(abs(qx-px)),np.max(abs(qy-py)))),
        horizon_y=horizon,earth=earth,region_fractions={str(i):float(np.mean(region==i)) for i in range(4)})
    (out/'nubium-midnight-calculated-guide.json').write_text(json.dumps(record,indent=2)+'\n')
    if not args.waves:
        return
    waves=s['waves']; wave_source=root/waves['file']
    # Use the same station as scene.sea_spectrum, rather than the offshore camera.
    station=next(t['lon_lat'] for t in inputs.calendar['sites'] if t['short']=='S Nubium')
    raw=read_hotfile(wave_source,nodes=[station]); f,d,e=raw['frequency'],raw['direction'],raw['variance'][0]
    depth=waves['depth_m']; gravity=MOON_SURFACE_GRAVITY
    fp=float(f[np.argmax(e.sum(axis=1))]); kp=float(wavenumber(np.array([fp]),depth,gravity)[0])
    k_from=float(wavenumber(np.array([f[-1]]),depth,gravity)[0])
    densities=[realization.swan_density(f,d,e,depth,gravity),
               realization.short_density(waves['friction_velocity_m_s'],(90-waves['wind_toward_compass_deg'])%360,gravity,k_from,2*math.pi*fp/kp)]
    seed=20380307
    tiles=realization.realize(densities,depth,gravity,seed,cascades=((512.,512),(31.7,512),(1.4,256)))
    pyramids=[realization.lean_pyramid(t.slope) for t in tiles]
    ocean=pixels.copy(); dist=np.geomspace(2,3000,5000)
    footprint=np.maximum(dist/cam.focal,dist**2/(eye*cam.focal))
    for col in range(width):
        ray=cam.direction(np.array([col+.5]),np.array([horizon]))[0]; az=math.atan2(ray[0],ray[1])
        east,north=dist*math.sin(az),dist*math.cos(az)
        z,sx,sy=np.zeros_like(dist),np.zeros_like(dist),np.zeros_like(dist)
        for tile,pyramid in zip(tiles,pyramids):
            if tile.height is not None:
                z+=map_coordinates(tile.height,[north/tile.spacing_m,east/tile.spacing_m],order=1,mode='grid-wrap')
            levels=np.clip(np.floor(np.log2(np.maximum(1,footprint/tile.spacing_m))).astype(int),0,len(pyramid)-1)
            for li in np.unique(levels):
                sel=levels==li; cell=tile.spacing_m*2**int(li); coord=[north[sel]/cell,east[sel]/cell]
                sx[sel]+=map_coordinates(pyramid[li][0],coord,order=1,mode='grid-wrap')
                sy[sel]+=map_coordinates(pyramid[li][1],coord,order=1,mode='grid-wrap')
        elev=np.degrees(np.arctan2(z-eye-dist**2/(2*MOON_RADIUS),dist)); _,rows=cam.project(scenes.unit(elev,np.full_like(elev,math.degrees(az))))
        v=np.stack([-east,-north,eye-z],-1); v/=np.linalg.norm(v,axis=1,keepdims=True)
        normal=np.stack([-sx,-sy,np.ones_like(sx)],-1);normal/=np.linalg.norm(normal,axis=1,keepdims=True)
        cos_v=np.maximum(np.sum(v*normal,axis=1),0);rho=reflection.fresnel(cos_v,index)
        reflected=2*cos_v[:,None]*normal-v;hit=(reflected@ed)>=math.cos(radius)
        contrast=np.clip(rho/np.maximum(reflection.fresnel(v[:,2],index),1e-6),.35,1.8)
        top=height
        for j in range(len(dist)):
            row=max(math.ceil(horizon),int(math.floor(rows[j])))
            if row<top:
                value=background[row:top,col]*contrast[j]
                if hit[j]: value=value+rho[j]*earth_xyz/earth['normal_illuminance_lux']*earth['disk_luminance_cd_m2']
                ocean[row:top,col]=display(value); top=row
        ocean[region[:,col]!=2,col]=pixels[region[:,col]!=2,col]
    save('nubium-midnight-wave-guide.png',ocean)
    wave_record=dict(schema='terluna.visualization.nubium-wave-guide/1',
        producer_sha256=digest(__file__),input_sha256=digest(wave_source),seed=seed,
        source_hs_m=waves['significant_height_m'],guide_hs_m=float(4*np.sqrt(sum(np.var(t.height) for t in tiles if t.height is not None))),
        source_mss=waves['mean_square_slope']['total'],guide_mss=float(sum(np.trace(realization.sample_covariance(t)) for t in tiles)),
        peak_wavelength_m=2*math.pi/kp,cascades=[[t.size_m,t.cells] for t in tiles],
        evidence='Domain SWAN realization plus lunar-gravity short-wave spectrum, with filtered slope pyramids. Wave snapshot is 11 hours after the illumination time; hourly wind sets the short waves.',
        reading_rule='Reduced display guide with approximate Fresnel contrast and an instantaneous uniform-Earth mirror-hit check. No spectral facet transport, terrain reflections, facet shadowing, finite-exposure averaging or validated breaking-wave treatment. Far-field subpixel normal variance and heights of short tiles are omitted.')
    (out/'nubium-midnight-wave-guide.json').write_text(json.dumps(wave_record,indent=2)+'\n')
    print(json.dumps(wave_record,indent=2))


if __name__=='__main__':
    main()
