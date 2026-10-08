"""Display a coastal radiance reference and its conditional cloud treatment.

One global tone curve preserves the calculated Earth detail. Historical NASA
maps supply geographic structure, not future weather. Optional cloud maps
supply angular radiance and extinction; their application to terrain and water
uses an explicit approximate illumination ratio, not a coupled terrain solver.
"""
from __future__ import annotations
import argparse,json,math,sys
from pathlib import Path
import numpy as np
from scipy.io import savemat,loadmat
from scipy.ndimage import map_coordinates,gaussian_filter
from scipy.interpolate import RegularGridInterpolator


def hidden_earth_projection(cam, source):
    """Keep geometry metadata without painting a planet over the sea."""
    ed = source[:3]
    radius = source[3]
    x, y = cam.project(ed)
    side = cam.r-(cam.r@ed)*ed
    side /= np.linalg.norm(side)
    up = np.cross(ed, side)
    angle = np.linspace(0, 2*np.pi, 721)
    limb = np.cos(radius)*ed+np.sin(radius)*(np.cos(angle)[:,None]*side+np.sin(angle)[:,None]*up)
    lx, ly = cam.project(limb)
    return dict(centre=[float(x[0]),float(y[0])], target_xyz=[0.,0.,0.], integrated_xyz=[0.,0.,0.],
                diameter_deg=float(np.degrees(2*radius)),
                disk_dimensions_px=[float(2*max(abs(lx-x[0]))),float(2*max(abs(ly-y[0])))], visible=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--scene',type=Path);p.add_argument('--prefix',default='nobili');p.add_argument('--sky-display',type=float,default=.025);p.add_argument('--linear-display',action='store_true');p.add_argument('--out',type=Path,required=True);p.add_argument('--cloud',type=Path);p.add_argument('--finish',action='store_true');p.add_argument('--display-reference',type=Path);p.add_argument('--sky-reference',type=Path);a=p.parse_args()
    root=(a.root or Path(__file__).resolve().parents[2]).resolve();sys.path.insert(0,str(root));sys.path.insert(0,str(root/'visualization/sea-appearance'));sys.path.insert(0,str(root/'visualization/reference-renderer/seas'))
    from nubium_guides import Camera
    from nubium_human_view import texture_earth,digest
    from nubium_visibility_view import tone,save_rgb
    import kernels
    a.out.mkdir(parents=True,exist_ok=True)
    source=(a.scene or root/'research/studies/sea_appearance/results/nobili-night.json').absolute();product=json.loads(source.read_text())
    with np.load(a.reference/'scene-build.npz') as z:
        sky=z['array/sky'];sky_el=z['array/sky_el'];az_step=float(z['array/params'][kernels.P_SKY_AZ_STEP]);earth=z['array/earth'];geometry=json.loads(str(z['geometry']))
    frames=list(a.reference.glob('*__night.npz'))
    if len(frames)!=1:raise ValueError('One completed scene frame is required')
    with np.load(frames[0]) as z:xyz=z['xyz'].astype(float)+z['stars'].astype(float);kind=z['kind'];render_record=json.loads(str(z['record']))
    if render_record['scenario_source_sha256']!=digest(source):raise ValueError('Radiance frame belongs to a different scene product')
    height,width=xyz.shape[:2];spec={**product['camera'],'frame':[width,height]};cam=Camera(spec)
    focus_spec={**spec,'frame':[1024,1024],'horizontal_fov_deg':10.,'view_azimuth_deg':geometry['earth_az'],'pitch_deg':geometry['earth_el']};focus_cam=Camera(focus_spec)
    assets=root/'visualization/sea-appearance/results/earth-inputs';manifest=json.loads((root/'visualization/sea-appearance/nubium-earth-inputs.json').read_text())
    for entry in manifest['files'][:2]:
        if digest(assets/entry['filename'])!=entry['sha256']:raise ValueError('NASA Earth map hash differs')
    from shared.constants import MOON_RADIUS
    h=product['earth']
    horizon=-np.degrees(np.arccos(MOON_RADIUS/(MOON_RADIUS+spec['eye_height_m'])))
    below_horizon=bool(h['elevation_deg']+h['diameter_deg']/2<horizon)
    projection=hidden_earth_projection(cam,earth) if below_horizon else texture_earth(xyz,spec,product,h,assets,sky,sky_el,az_step,earth,geometry,ss=6,normalize_xyz=True)
    yy,xx=np.mgrid[:1024,:1024];frays=focus_cam.direction(xx.ravel()+.5,yy.ravel()+.5)
    foreground=np.stack([kernels.sky_xyz(sky,sky_el,az_step,d) for d in frays]).reshape(1024,1024,3)
    focus=foreground.copy()
    if not below_horizon:texture_earth(focus,focus_spec,product,h,assets,sky,sky_el,az_step,earth,geometry,ss=2,normalize_xyz=True)
    direct=focus-foreground
    cloud_record=None
    if a.cloud:
        with np.load(a.cloud) as z:
            ce=z['elevation_deg'];ca=z['azimuth_deg'];cloud=z['cloudy_groups'].mean(axis=2);clear=z['clear_groups'].mean(axis=2);tau=z['cloud_tau'];cloud_record=json.loads(str(z['metadata']))
        # A smooth low-resolution correction from actual cloudy-vs-clear rays;
        # Monte Carlo errors and this interpolation remain separate evidence.
        from cloud_display import delta_interpolator
        smoothing=render_record.get('cloud_conditioning',{}).get('angular_smoothing_deg',0.)
        correction=delta_interpolator(ce,ca,cloud,clear,smoothing)
        aa=np.r_[ca[-1]-360,ca,ca[0]+360]
        transmission=RegularGridInterpolator((ce,aa),np.concatenate([tau[:,-1:],tau,tau[:,:1]],axis=1),bounds_error=False,fill_value=None)
        ratio=cloud_record.get('total_horizontal_lux',cloud_record['earth_total_horizontal_lux'])/cloud_record.get('clear_total_horizontal_lux',product['light']['earth_lux'])
        already_cloudy='cloud_conditioning' in render_record
        def angles(rays):
            return np.stack([np.clip(np.degrees(np.arcsin(rays[:,2])),ce[0],ce[-1]),np.degrees(np.arctan2(rays[:,0],rays[:,1]))%360],-1)
        # The terrain and sea's cloud-shadow modulation is an approximation;
        # their clear-sky paths and geometry remain in the saved control.
        if not already_cloudy:xyz*=1+(np.clip(ratio,0,2)-1)*np.minimum(kind.sum(axis=-1),1)[...,None]
        for row in range(0,height,32):
            end=min(row+32,height);yy,xx=np.mgrid[row:end,:width];rays=cam.direction(xx.ravel()+.5,yy.ravel()+.5)
            sky_share=1-np.clip(kind[row:end].sum(axis=-1),0,1)
            if not already_cloudy:xyz[row:end]+=correction(angles(rays)).reshape(end-row,width,3)*sky_share[...,None]
        # Replace the Earth patch with the cloud-attenuated direct component;
        # use an enlarged direct field to retain phase and geographical texture.
        fc=angles(frays)
        if not already_cloudy:foreground=np.maximum(foreground+correction(fc).reshape(1024,1024,3),0)
        direct*=np.exp(-np.maximum(transmission(fc),0)).reshape(1024,1024,1)
        if already_cloudy:direct/=render_record['cloud_conditioning']['disk_centre_transmission']
        if not below_horizon:
            x0=max(0,int(projection['centre'][0]-projection['disk_dimensions_px'][0]/2)-4);x1=min(width,int(projection['centre'][0]+projection['disk_dimensions_px'][0]/2)+5)
            y0=max(0,int(projection['centre'][1]-projection['disk_dimensions_px'][1]/2)-4);y1=min(height,int(projection['centre'][1]+projection['disk_dimensions_px'][1]/2)+5)
            yy,xx=np.mgrid[y0:y1,x0:x1];rays=cam.direction(xx.ravel()+.5,yy.ravel()+.5);fx,fy=focus_cam.project(rays)
            fg=np.stack([kernels.sky_xyz(sky,sky_el,az_step,d) for d in rays])
            if not already_cloudy:fg+=correction(angles(rays))
            projected=np.stack([map_coordinates(direct[...,c],[fy-.5,fx-.5],order=1,mode='constant',cval=0) for c in range(3)],-1)
            xyz[y0:y1,x0:x1]=np.maximum((fg+projected).reshape(y1-y0,x1-x0,3),0)
    xyz=np.maximum(xyz,0)
    direct_absent=not np.any(direct)
    if a.finish:
        data=loadmat(a.out/'earth-optics.mat');optical=data['optical']
        # Add the ocular change over all surfaces, including the rim beneath
        # Earth, instead of replacing those surfaces with a sky background.
        delta=optical-direct
        for row in range(0,height,32):
            end=min(row+32,height);yy,xx=np.mgrid[row:end,:width];rays=cam.direction(xx.ravel()+.5,yy.ravel()+.5)
            theta=np.degrees(np.arccos(np.clip(rays@focus_cam.f,-1,1)));use=theta<5
            if not np.any(use):continue
            fx,fy=focus_cam.project(rays[use]);weight=.5*(1+np.cos(np.pi*np.clip(theta[use]-4,0,1)))
            value=np.stack([map_coordinates(delta[...,c],[fy-.5,fx-.5],order=1,mode='constant',cval=0) for c in range(3)],-1)*weight[:,None]
            block=xyz[row:end].reshape(-1,3);block[use]+=value
        focus=foreground+optical;xyz=np.maximum(xyz,0)
    else:
        focus=foreground+direct
        savemat(a.out/'earth-optics-input.mat',dict(direct=direct,foreground=foreground,ppd=focus_cam.focal*np.pi/180),do_compression=True)
        cfg=json.loads((root/'visualization/sea-appearance/results/nubium-visibility/optics-config.json').read_text())
        cfg.update(input=str(a.out.resolve()/'earth-optics-input.mat'),output=str(a.out.resolve()/'earth-optics.mat'),negative_energy_tolerance=1e-7)
        (a.out/'optics-config.json').write_text(json.dumps(cfg,indent=2)+'\n')
        if direct_absent:
            # Complete cloud occultation leaves no resolved direct signal;
            # applying FFT optics to a zero field would produce 0/0 diagnostics.
            savemat(a.out/'earth-optics.mat',dict(optical=np.zeros_like(direct),foreground=foreground,ppd=focus_cam.focal*np.pi/180,negative_fraction=0.))
    # A declared global display mapping, not a guarantee of naked-eye appearance.
    white=float(max(xyz[...,1].max(),0. if below_horizon else focus[...,1].max())*1.1)
    sky_anchor=float(np.median(xyz[:max(1,height//8),:,1]));knee=.005;sky_display=a.sky_display
    if a.sky_reference:
        anchor=json.loads(a.sky_reference.read_text())['display_curve']
        sky_anchor=anchor['sky_world_cd_m2'];knee=anchor['knee_world_cd_m2'];sky_display=anchor['sky_display_fraction']
    if a.display_reference:
        fixed=json.loads(a.display_reference.read_text())['display_curve']
        white=fixed['white_world_cd_m2'];sky_anchor=fixed['sky_world_cd_m2'];knee=fixed['knee_world_cd_m2'];sky_display=fixed['sky_display_fraction']
    linear_mode=a.linear_display or (bool(a.display_reference) and fixed.get('mode')=='linear-exposure')
    if linear_mode:
        from research.studies.sea_appearance.scenes import XYZ_TO_SRGB
        scale=sky_display/sky_anchor
        def linear_display(field):
            y=np.maximum(field[...,1]*scale,0)
            mapped=np.where(y<=.8,y,.8+.2*(-np.expm1(-(y-.8)/.2)))
            ratio=np.divide(mapped,y,out=np.ones_like(y),where=y>0)
            rgb=(field*scale*ratio[...,None])@XYZ_TO_SRGB.T
            delta=rgb-mapped[...,None];amount=np.ones_like(mapped)
            for c in range(3):
                d=delta[...,c]
                amount=np.minimum(amount,np.where(d>0,(1-mapped)/np.maximum(d,1e-30),mapped/np.maximum(-d,1e-30)))
            return np.clip(mapped[...,None]+np.clip(amount,0,1)[...,None]*delta,0,1)
        rgb=linear_display(xyz);frgb=linear_display(focus)
        curve=dict(mode='linear-exposure',white_world_cd_m2=1/scale,knee_world_cd_m2=knee,
            sky_world_cd_m2=sky_anchor,sky_display_fraction=sky_display,
            equation='XYZ scaled by sky_display/sky_world; luminance shoulder starts at 0.8; encoded as sRGB',
            scope='Single photographic exposure, not a calibrated human colour-appearance model.')
    else:
        rgb,curve=tone(xyz,white,knee=knee,sky_world=sky_anchor,sky_display=sky_display)
        frgb,_=tone(focus,white,knee=knee,sky_world=sky_anchor,sky_display=sky_display)
    suffix='view' if a.finish else 'guide'
    save_rgb(a.out/f'{a.prefix}-night-{suffix}.png',rgb);save_rgb(a.out/f'{a.prefix}-earth-{suffix}.png',frgb)
    np.savez_compressed(a.out/'display-xyz.npz',xyz=xyz.astype('f4'),kind=kind)
    projection={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in projection.items()}
    record=dict(schema='terluna.visualization.nobili-view/1',scene_sha256=digest(source),
        camera=spec,earth_projection=projection,display_curve=curve,display_reference_sha256=digest(a.display_reference) if a.display_reference else None,display_reference=str(a.display_reference) if a.display_reference else None,sky_reference=str(a.sky_reference) if a.sky_reference else None,cloud=cloud_record,
        eye_optics=('No resolved direct Earth; zero-field optics bypass' if direct_absent else 'CIE99 age 24, direct Earth only; halo taper 4–5 degrees') if a.finish else 'Not yet applied',
        direct_earth_numerically_absent=direct_absent,earth_below_horizon=below_horizon,
        fft_negative_energy_fraction_zeroed=float(data['negative_fraction'][0,0]) if a.finish else None,
        fft_negative_energy_tolerance=1e-7,
        earth_texture='Historical NASA global Blue Marble surface and cloud composite with JPL geography and Lambert illumination, normalized to the calculated Earth beam. Spatial spectra and future weather are unvalidated.',
        caveats=['Global tone mapping and ocular optics are display choices, not exact human colour appearance.',
            'Lunar cloud radiance is a 2-D extrusion scenario. Water reflects the conditional angular sky and attenuated Earth; terrain diffuse spectral shape, local cloud shadows and aerial perspective remain approximations.'] if a.cloud else ['Clear molecular atmosphere; terrain reflectance and equilibrium lake wave age are scenario inputs.'],
        render_record=render_record,producer_sha256=digest(__file__),outputs={name:digest(a.out/name) for name in [f'{a.prefix}-night-{suffix}.png',f'{a.prefix}-earth-{suffix}.png']})
    (a.out/'display.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(projection=projection,curve=curve,cloud=bool(a.cloud)),indent=2))
if __name__=='__main__':main()
