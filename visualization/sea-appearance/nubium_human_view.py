"""Apply Radiance's human-response display operator to the Nubium XYZ frame.

The high-dynamic-range input retains physical units. A dated NASA-textured
Earth replaces only the tracer's uniform direct disk, with its integral matched
to the domain's direct earthlight. No independent Earth or landscape exposure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def condition(xyz, camera, out, bins, *, dynamic_range=1000, fixation=None):
    """XYZ in cd/m² through official pvalue and pcond, then exact sRGB encoding."""
    from image_guides import encoded
    height,width=xyz.shape[:2]
    view=f"VIEW= -vtv -vh {camera['horizontal_fov_deg']} -vv {camera['vertical_fov_deg']}\n"
    binary=np.ascontiguousarray(xyz,dtype=np.float32).tobytes()
    raw=subprocess.run([str(bins/'pvalue'),'-r','-h','-df','-pXYZ','-y',str(height),'+x',str(width)],
                       input=binary,capture_output=True,check=True).stdout
    raw=raw.replace(b'\n\n',b'\n'+view.encode()+b'\n',1)
    hdr=out.with_suffix('.xyze'); hdr.write_bytes(raw)
    # An independent round trip through Radiance catches unit/order/format errors.
    recovered=subprocess.run([str(bins/'pvalue'),'-h','-H','-df',str(hdr)],capture_output=True,check=True).stdout
    recovered=np.frombuffer(recovered,dtype=np.float32).reshape(xyz.shape)
    roundtrip=float(np.max(np.abs(recovered[...,1]-xyz[...,1])/np.maximum(xyz[...,1],1e-8)))
    # XYZE uses a shared exponent for three 8-bit mantissas. A faint channel
    # can lose more than 1% while remaining within the format's quantization.
    exponent_error=float(np.max(np.abs(recovered[...,1]-xyz[...,1])/
                                np.maximum(np.max(np.abs(xyz),axis=2),1e-8)))
    if exponent_error>1/256+1e-7:
        raise ValueError(f'XYZ roundtrip exceeds shared-exponent precision: {exponent_error}')
    output_hdr=out.with_suffix('.hdr'); mapping=out.with_suffix('.curve.txt')
    command=[str(bins/'pcond'),'-h+','-u','100','-d',str(dynamic_range),
             '-p','.64','.33','.30','.60','.15','.06','.3127','.3290',
             '-x',str(mapping),str(hdr)]
    fixation_input=None
    if fixation is not None:
        # Radiance fixation coordinates have their origin at the lower left.
        command[1:1]=['-i','1']
        fixation_input=f'{round(fixation[0])} {height-1-round(fixation[1])}\n'.encode()
    env=os.environ.copy(); env['PATH']=str(bins)+os.pathsep+env.get('PATH','')
    completed=subprocess.run(command,input=fixation_input,capture_output=True,check=True,env=env)
    output_hdr.write_bytes(completed.stdout)
    converted=subprocess.run([str(bins/'pvalue'),'-h','-H','-df',str(output_hdr)],capture_output=True,check=True).stdout
    rgb=np.frombuffer(converted,dtype=np.float32).reshape(xyz.shape)
    image=np.uint8(np.rint(np.clip(encoded(rgb),0,1)*255))
    Image.fromarray(image).save(out.with_suffix('.png'))
    return dict(command=command,stderr=completed.stderr.decode(),display_peak_cd_m2=100,
                fixation_input=None if fixation_input is None else fixation_input.decode(),
                adaptation='Uniform image sampling' if fixation is None else 'Earth-centred fixation only; sensitivity case, not a measured gaze history',
                operator_chose_linear_mapping=b'EXPOSURE=' in completed.stdout.split(b'\n\n')[0],
                display_dynamic_range=dynamic_range,display_primaries='sRGB/D65',
                transfer_function='exact piecewise sRGB',input_xyz_luminance_roundtrip_max_relative=roundtrip,
                input_xyz_roundtrip_error_per_largest_component=exponent_error,
                display_linear_rgb_range=[float(rgb.min()),float(rgb.max())],
                radiance_binaries={k:digest(bins/k) for k in ['pvalue','pcond','pfilt']},
                output_sha256=digest(out.with_suffix('.png')))


def texture_earth(base, spec, product, h, assets, sky_map, sky_el, az_step, earth_source, geometry, *, ss=8, sampler=None):
    """Project a historical Earth in physical XYZ units into a supplied camera."""
    from nubium_guides import Camera
    from image_guides import linear
    from research.studies.sea_appearance import scenes, lighting
    from illumination.stars import sky as stars
    from shared.constants import EARTH_RADIUS
    import kernels
    height,width=base.shape[:2]
    cam=Camera(spec)
    ed=earth_source[:3];radius=earth_source[3];distance=EARTH_RADIUS/math.sin(radius)
    cx,cy=cam.project(ed);cx,cy=float(cx[0]),float(cy[0])
    side=cam.r-(cam.r@ed)*ed;side/=np.linalg.norm(side)
    up=np.cross(ed,side)
    angle=np.linspace(0,2*math.pi,721)
    limb=math.cos(radius)*ed+math.sin(radius)*(np.cos(angle)[:,None]*side+np.sin(angle)[:,None]*up)
    lx,ly=cam.project(limb)
    half_width=float(max(abs(lx-cx)));half_height=float(max(abs(ly-cy)))
    x0=max(0,int(cx-half_width)-3);x1=min(width,int(cx+half_width)+4)
    y0=max(0,int(cy-half_height)-3);y1=min(height,int(cy+half_height)+4)
    # Subpixel angular sampling resolves the limb and the very thin terminator.
    yy,xx=np.meshgrid(y0+(np.arange((y1-y0)*ss)+.5)/ss,x0+(np.arange((x1-x0)*ss)+.5)/ss,indexing='ij')
    shape=xx.shape;directions=cam.direction(xx.ravel(),yy.ravel());ct=directions@ed
    discriminant=EARTH_RADIUS**2-distance**2*(1-ct*ct)
    inside=discriminant>=0
    t=distance*ct-np.sqrt(np.maximum(0,discriminant))
    normal=(t[:,None]*directions-distance*ed)/EARTH_RADIUS
    pole_body=stars.body_directions(np.radians(h['north_pole_ra_deg']),np.radians(h['north_pole_dec_deg']),product['jd_tt'])[0]
    e,n,u=lighting.local_frame(*spec['lon_lat_deg']);pole=np.array([pole_body@e,pole_body@n,pole_body@u])
    tr=cam.r-(cam.r@ed)*ed;tr/=np.linalg.norm(tr);tu=np.cross(tr,ed)
    north_projected=pole-(pole@ed)*ed;north_projected/=np.linalg.norm(north_projected)
    nr,nu=north_projected@tr,north_projected@tu
    lon,lat=np.radians([h['subobserver_lon_east_deg'],h['subobserver_lat_deg']])
    c=np.array([math.cos(lat)*math.cos(lon),math.cos(lat)*math.sin(lon),math.sin(lat)])
    gn=np.array([-math.sin(lat)*math.cos(lon),-math.sin(lat)*math.sin(lon),math.cos(lat)])
    ge=np.array([-math.sin(lon),math.cos(lon),0])
    rb=nu*ge+nr*gn;ub=-nr*ge+nu*gn
    body=(normal@tr)[:,None]*rb+(normal@tu)[:,None]*ub+(normal@(-ed))[:,None]*c
    glon=np.arctan2(body[:,1],body[:,0]);glat=np.arcsin(np.clip(body[:,2],-1,1))
    if sampler is None:
        surface=linear(np.asarray(Image.open(assets/'nubium-earth-bluemarble-2048.png').convert('RGB'))/255)
        clouds=np.asarray(Image.open(assets/'nubium-earth-clouds-2048.jpg').convert('L'))/255
        tx=(glon+math.pi)/(2*math.pi)*surface.shape[1]-.5;ty=(math.pi/2-glat)/math.pi*surface.shape[0]-.5
        texture=np.stack([map_coordinates(surface[...,ch],[ty,tx],order=1,mode='grid-wrap') for ch in range(3)],-1)
        alpha=map_coordinates(clouds,[ty,tx],order=1,mode='nearest')[:,None]*.96
        texture=texture*(1-alpha)+.9*alpha
    else:
        texture=sampler(body,inside)
    slon,slat=np.radians([h['subsolar_lon_east_deg'],h['subsolar_lat_deg']])
    sun=np.array([math.cos(slat)*math.cos(slon),math.cos(slat)*math.sin(slon),math.sin(slat)])
    incident=np.maximum(body@sun,0); texture*=incident[:,None]*inside[:,None]
    omega=(directions@cam.f)**3/cam.focal**2/ss**2
    # Calibrate each broad RGB channel to the modeled beam, then return to XYZ.
    # This constrains the integrated colour and illuminance but does not validate
    # spatially resolved spectra, cloud BRDF, or historical cloud forecasts.
    target_xyz=earth_source[7:10]*earth_source[10]*geometry['earth_lit_fraction']
    target_rgb=target_xyz@scenes.XYZ_TO_SRGB.T
    integral=np.sum(texture*(omega*ct)[:,None],axis=0)
    texture*=target_rgb/integral
    direct_xyz=texture@np.linalg.inv(scenes.XYZ_TO_SRGB).T
    foreground=np.stack([kernels.sky_xyz(sky_map,sky_el,az_step,d) for d in directions])
    patch=(direct_xyz+foreground).reshape(shape+(3,)).reshape(y1-y0,ss,x1-x0,ss,3).mean(axis=(1,3))
    base[y0:y1,x0:x1]=patch
    expected=np.sum(direct_xyz*(omega*ct)[:,None],axis=0)
    return dict(centre=[cx,cy],target_xyz=target_xyz,integrated_xyz=expected,
                diameter_deg=math.degrees(2*radius),disk_dimensions_px=[2*half_width,2*half_height])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path)
    parser.add_argument('--render',type=Path,required=True)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--assets',type=Path,required=True)
    parser.add_argument('--radiance-bin',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    root=(args.root or Path(__file__).resolve().parents[2]).resolve()
    sys.path.insert(0,str(root));sys.path.insert(0,str(root/'visualization/sea-appearance'))
    sys.path.insert(0,str(root/'visualization/reference-renderer/seas'))
    from nubium_guides import Camera
    import kernels
    source=root/'research/studies/sea_appearance/results/nubium-midnight.json'
    product=json.loads(source.read_text());s=product['scene'];earth=s['earth']
    manifest_path=root/'visualization/sea-appearance/nubium-earth-inputs.json'
    manifest=json.loads(manifest_path.read_text());h=manifest['horizons']
    for entry in manifest['files']:
        if digest(args.assets/entry['filename'])!=entry['sha256']:
            raise ValueError('Source hash mismatch: '+entry['filename'])
    with np.load(args.render,allow_pickle=False) as data:
        base=data['xyz'].astype(float)+data['stars'].astype(float)
        error=data['error'];kind=data['kind'];render_record=json.loads(str(data['record']))
    height,width=base.shape[:2]
    spec={**s['camera'],'frame':[width,height],
          'horizon_row_at_centre':s['camera']['horizon_row_at_centre']*height/s['camera']['frame'][1]}
    cam=Camera(spec)
    with np.load(args.build,allow_pickle=False) as data:
        sky_map=data['array/sky'];sky_el=data['array/sky_el'];params=data['array/params'];earth_source=data['array/earth']
        geometry=json.loads(str(data['geometry']))
    projection=texture_earth(base,spec,product,h,args.assets,sky_map,sky_el,params[kernels.P_SKY_AZ_STEP],earth_source,geometry)
    cx,cy=projection['centre']
    target_xyz=projection['target_xyz'];expected=projection['integrated_xyz']
    args.out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.out/'nubium-midnight-textured-xyz.npz',xyz=base.astype(np.float32),kind=kind,error=error)
    main_record=condition(base,spec,args.out/'nubium-midnight-human-guide',args.radiance_bin)
    # The display uncertainty is checked explicitly, without selecting an
    # arbitrary camera exposure to force a desired nocturnal colour palette.
    sensitivity=condition(base,spec,args.out/'nubium-midnight-human-contrast100',args.radiance_bin,dynamic_range=100)
    earth_fixation=condition(base,spec,args.out/'nubium-midnight-earth-fixation',args.radiance_bin,fixation=(cx,cy))
    # An explicitly narrower angular window explains apparent size without
    # enlarging Earth inside the wide sea frame. The horizon is outside it.
    focus_spec={**spec,'frame':[1024,1024],'horizontal_fov_deg':10.,'vertical_fov_deg':10.,
                'view_azimuth_deg':geometry['earth_az'],'pitch_deg':geometry['earth_el']}
    focus_camera=Camera(focus_spec)
    focus_spec['horizon_row_at_centre']=512+focus_camera.focal*math.tan(
        math.radians(focus_spec['pitch_deg']+focus_spec['horizon_dip_deg']))
    yy,xx=np.mgrid[:1024,:1024]
    directions=focus_camera.direction(xx.ravel()+.5,yy.ravel()+.5)
    focus=np.stack([kernels.sky_xyz(sky_map,sky_el,params[kernels.P_SKY_AZ_STEP],d)
                    for d in directions]).reshape(1024,1024,3)
    focus_projection=texture_earth(focus,focus_spec,product,h,args.assets,sky_map,sky_el,
                                  params[kernels.P_SKY_AZ_STEP],earth_source,geometry,ss=4)
    focus_display=condition(focus,focus_spec,args.out/'nubium-midnight-earth-focused',args.radiance_bin,fixation=(512,512))
    focus_record=dict(camera=focus_spec,display=focus_display,
                      earth_diameter_deg=focus_projection['diameter_deg'],
                      earth_disk_dimensions_px=focus_projection['disk_dimensions_px'],
                      matched_viewing_distance_per_image_width=1/(2*math.tan(math.radians(5))),
                      reading_rule='10-degree square angular crop centred on Earth, with Earth-centred display adaptation. Not a telescope or a larger physical Earth. The sea horizon is 53 degrees below Earth and must not appear in this view. Historical clouds are illustrative; resolving all rendered detail with an unaided eye is not established.')
    (args.out/'nubium-midnight-earth-focused.json').write_text(json.dumps(focus_record,indent=2)+'\n')
    sea=kind[...,0]>.99
    record=dict(schema='terluna.visualization.nubium-human-view/1',
        producer={'script':'visualization/sea-appearance/nubium_human_view.py','sha256':digest(__file__)},
        inputs={str(f):digest(f) for f in [args.render,args.build,source,manifest_path]},
        evidence='Existing sea path tracer XYZ with explicit lunar-gravity wave geometry, filtered facets, terrain shadows, water colour and aerial transport. Direct Earth is a geometrically projected historical NASA texture calibrated to the domain beam. Official Radiance pcond applies contrast sensitivity, acuity, veiling glare and mesopic colour response to the common physical-radiance input.',
        reading_rule='One documented human-vision display approximation, not empirical validation of an Open Moon or a guarantee of a particular observer/device. The main wide view loses bright Earth detail in display clipping; this is not evidence that a person cannot resolve Earth. The Earth-fixation sensitivity retains more disk detail while darkening the sea. No single frame reproduces looking around and adapting. Fixed 2-m eye camera represents an offshore observer, not a solved standing platform. Clear molecular atmosphere and productive-water design guess; no lunar clouds or mature terrestrial ecology supplied. Spatial Earth spectra are a normalized RGB approximation; reflected Earth still uses the original uniform-disk model. User display luminance, surround, ocular characteristics and adaptation history are unknown.',
        sources={'method':'https://graphics.cs.yale.edu/sites/default/files/1997tvcg_hdr_tonemapping.pdf',
                 'manual':'https://radsite.lbl.gov/radiance/man_html/pcond.1.html',
                 'implementation':'https://github.com/LBNL-ETA/Radiance',
                 'implementation_commit':'bcffc2b52d99b908adfdc2543e4cddf83268a125'},
        display=main_record,contrast100_sensitivity=sensitivity,earth_fixation_sensitivity=earth_fixation,
        earth_focused_view=focus_record,
        camera=spec,matched_viewing_distance_per_image_width=1/(2*math.tan(math.radians(spec['horizontal_fov_deg'])/2)),
        earth_normal_illuminance_xyz=target_xyz.tolist(),integrated_texture_xyz=expected.tolist(),
        earth_integral_max_relative_error=float(np.max(abs(expected-target_xyz)/target_xyz)),
        sea_luminance_percentiles=np.percentile(base[...,1][sea],[5,50,95]).tolist(),
        sea_standard_error_percentiles=np.percentile(error[sea],[50,90,99]).tolist(),
        render=render_record)
    (args.out/'nubium-midnight-human-guide.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['earth_integral_max_relative_error','sea_luminance_percentiles','sea_standard_error_percentiles','display']},indent=2))


if __name__=='__main__':
    main()
