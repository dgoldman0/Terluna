"""Prepare a display reference retaining visible Earth structure.

The view consumes physical scene XYZ; its tone curve and spectral spatial
approximation are display choices. It is not an exact naked-eye simulation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from PIL import Image
from scipy.io import savemat
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

XYZ_RGB=np.array([[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],
                  [.0556434,-.2040259,1.0572252]])


def tone(xyz, white, knee=.05, sky_world=.5291728, sky_display=.035):
    """One monotone log-power luminance curve for the entire scene.

    The sky anchor is a declared display choice; not a perceptual measurement.
    Gamut is compressed toward neutral at fixed display Y, avoiding RGB clips.
    There is no Earth-specific exposure or pasted-in brighter/darker disk.
    """
    y=np.maximum(xyz[...,1],0)
    power=math.log(sky_display)/math.log(math.log1p(sky_world/knee)/math.log1p(white/knee))
    display=np.minimum(np.log1p(y/knee)/math.log1p(white/knee),1)**power
    rgb=xyz@XYZ_RGB.T/np.maximum(y[...,None],1e-15)*display[...,None]
    delta=rgb-display[...,None];amount=np.ones_like(display)
    for c in range(3):
        d=delta[...,c]
        amount=np.minimum(amount,np.where(d>0,(1-display)/np.maximum(d,1e-30),
                                         display/np.maximum(-d,1e-30)))
    rgb=display[...,None]+np.clip(amount,0,1)[...,None]*delta
    return np.clip(rgb,0,1),dict(white_world_cd_m2=white,knee_world_cd_m2=knee,power=power,
                               sky_world_cd_m2=sky_world,sky_display_fraction=sky_display,
                               equation='[log(1+Y/knee)/log(1+white/knee)]^power',
                               scope='Same curve at every pixel; monotone display approximation with retained detail, not calibrated colour appearance.')


def save_rgb(path,rgb):
    srgb=np.where(rgb<=.0031308,12.92*rgb,1.055*np.maximum(rgb,0)**(1/2.4)-.055)
    Image.fromarray(np.uint8(np.rint(np.clip(srgb,0,1)*255))).save(path)


class EpicSampler:
    """Reproject a historical observed Earth onto the scene's JPL orientation.

    Three narrow bands supply relative spatial RGB structure, individually
    normalized by the caller to the solved XYZ beam. This does not turn them
    into measured full-spectrum spatial radiance. A Lambert incidence ratio
    approximates the small phase change. The unobserved sliver uses nearest
    observed surface-direction samples, explicitly recorded as an illustration.
    """
    def __init__(self,path):
        import h5py
        with h5py.File(path) as f:
            lat=f['Geolocation/Earth/Latitude'][:];lon=f['Geolocation/Earth/Longitude'][:]
            self.valid=np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)
            yy,xx=np.mgrid[:2048,:2048]
            def unit(la,lo):
                la,lo=np.radians(la),np.radians(lo)
                return np.stack([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)],-1)
            sample=self.valid[::12,::12]
            bodies=unit(lat[::12,::12][sample],lon[::12,::12][sample])
            px=(xx[::12,::12][sample]-1023.5)/1024
            py=(yy[::12,::12][sample]-1023.5)/1024
            b=np.c_[bodies,np.ones(len(bodies))];zero=np.zeros_like(b)
            matrix=np.r_[np.c_[b,zero,-px[:,None]*bodies],np.c_[zero,b,-py[:,None]*bodies]]
            fit=np.linalg.lstsq(matrix,np.r_[px,py],rcond=None)[0]
            self.fx,self.fy,self.den=fit[:4],fit[4:8],fit[8:]
            predicted_x=(b@self.fx)/(1+bodies@self.den)
            predicted_y=(b@self.fy)/(1+bodies@self.den)
            error=1024*np.hypot(predicted_x-px,predicted_y-py)
            source_mu=np.cos(np.radians(f['Geolocation/Earth/SunAngleZenith'][::12,::12][sample]))
            self.sun=np.linalg.lstsq(bodies,source_mu,rcond=None)[0]
            la=float(f.attrs['centroid_mean_latitude'][0]);lo=float(f.attrs['centroid_mean_longitude'][0])
            self.view=unit(la,lo)
            self.bands=[]
            for band in [680,551,443]:
                raw=f[f'Band{band}nm'][:]
                raw=np.where(self.valid&np.isfinite(raw),np.maximum(raw,0),0)
                self.bands.append(raw/np.mean(raw[self.valid]))
            # Sparse native samples suffice only for the explicitly unobserved rim.
            sparse=self.valid[::4,::4]
            self.tree=cKDTree(unit(lat[::4,::4][sparse],lon[::4,::4][sparse]))
            self.nearest_xy=np.c_[xx[::4,::4][sparse],yy[::4,::4][sparse]]
        if np.percentile(error,99)>2:
            raise ValueError('EPIC geographic projection fit is too inaccurate')
        self.record=dict(projection_fit_rms_source_pixels=float(np.sqrt(np.mean(error**2))),
                         projection_fit_p99_source_pixels=float(np.percentile(error,99)),
                         source='epic_1b_20151117002712_00.h5',bands_rgb_nm=[680,551,443],
                         source_centre_lon_lat_deg=[lo,la],
                         reading_rule='Historical observed narrow-band spatial structure, reprojected and normalized to the scene beam. Illustrative 2038 weather; Lambert phase transfer and nearest observed unmeasured rim; spatial colours not spectrally validated.')
        self.calls=[]

    def __call__(self,body,inside):
        b=np.c_[body,np.ones(len(body))]
        den=1+body@self.den
        x=(b@self.fx)/den*1024+1023.5;y=(b@self.fy)/den*1024+1023.5
        valid=map_coordinates(self.valid.astype(float),[y,x],order=1,mode='constant',cval=0)>.999
        good=valid&(body@self.view>.02)
        missing=inside&~good
        if np.any(missing):
            _,idx=self.tree.query(body[missing])
            x[missing],y[missing]=self.nearest_xy[idx].T
        texture=np.stack([map_coordinates(a,[y,x],order=1,mode='constant',cval=0) for a in self.bands],-1)
        texture/=np.maximum(body@self.sun,.1)[:,None]
        self.calls.append(dict(inside_samples=int(inside.sum()),unobserved_rim_samples=int(missing.sum()),
                               unobserved_rim_fraction=float(missing.sum()/inside.sum())))
        return texture


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--epic',type=Path,required=True)
    p.add_argument('--render',type=Path,required=True);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--assets',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(args.root));sys.path.insert(0,str(args.root/'visualization/reference-renderer/seas'))
    sys.path.append(str(args.root/'visualization/sea-appearance'))
    from nubium_human_view import texture_earth,digest
    from nubium_guides import Camera
    import kernels
    product=json.loads((args.root/'research/studies/sea_appearance/results/nubium-midnight.json').read_text())
    h=json.loads((args.root/'visualization/sea-appearance/nubium-earth-inputs.json').read_text())['horizons']
    with np.load(args.build) as z:
        sky=z['array/sky'];sky_el=z['array/sky_el'];az_step=z['array/params'][kernels.P_SKY_AZ_STEP]
        earth=z['array/earth'];geometry=json.loads(str(z['geometry']))
    with np.load(args.render) as z:
        xyz=z['xyz'].astype(float)+z['stars'].astype(float);kind=z['kind'];error=z['error']
    height,width=xyz.shape[:2]
    s=product['scene'];spec={**s['camera'],'frame':[width,height],
                           'horizon_row_at_centre':s['camera']['horizon_row_at_centre']*height/s['camera']['frame'][1]}
    manifest_path=args.root/'research/studies/sea_appearance/visibility_inputs.json'
    if not manifest_path.exists(): manifest_path=Path(__file__).with_name('visibility_inputs.json')
    source_manifest=json.loads(manifest_path.read_text())
    expected=next(row['sha256'] for row in source_manifest['files'] if row['filename']==args.epic.name)
    if digest(args.epic)!=expected: raise ValueError('EPIC source hash mismatch')
    sampler=EpicSampler(args.epic)
    projection=texture_earth(xyz,spec,product,h,args.assets,sky,sky_el,az_step,earth,geometry,sampler=sampler)
    focus_spec={**spec,'frame':[1024,1024],'horizontal_fov_deg':10.,'vertical_fov_deg':10.,
                'view_azimuth_deg':geometry['earth_az'],'pitch_deg':geometry['earth_el']}
    cam=Camera(focus_spec);yy,xx=np.mgrid[:1024,:1024]
    rays=cam.direction(xx.ravel()+.5,yy.ravel()+.5)
    focus=np.stack([kernels.sky_xyz(sky,sky_el,az_step,d) for d in rays]).reshape(1024,1024,3)
    foreground=focus.copy()
    fp=texture_earth(focus,focus_spec,product,h,args.assets,sky,sky_el,az_step,earth,geometry,ss=4,sampler=sampler)
    white=float(max(xyz[...,1].max(),focus[...,1].max())*1.2)
    wide_rgb,curve=tone(xyz,white);focus_rgb,_=tone(focus,white)
    save_rgb(args.out/'nubium-midnight-visibility-guide.png',wide_rgb)
    save_rgb(args.out/'nubium-earth-visibility-guide.png',focus_rgb)
    np.savez_compressed(args.out/'nubium-visibility-xyz.npz',xyz=xyz.astype(np.float32),kind=kind,error=error)
    np.save(args.out/'nubium-earth-xyz.npy',focus.astype(np.float32))
    savemat(args.out/'nubium-earth-optics-input.mat',dict(direct=focus-foreground, foreground=foreground, ppd=cam.focal*np.pi/180),do_compression=True)
    projection={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in projection.items()}
    record=dict(schema='terluna.visualization.nubium-visibility-view/1',
                producer_sha256=digest(__file__),geometry=projection,camera=spec,focus_camera=focus_spec,
                display_curve=curve,earth_texture={**sampler.record,'sampling':sampler.calls},
                evidence='Physical sea-tracer XYZ and observed historical Earth structure normalized to the existing direct beam. Global display curve retains detail; no claim of exact naked-eye appearance.',
                limitations=['These display guides omit extra ocular glare and neural colour adaptation; visibility is assessed separately using HDR-VDP-3 and CIE99.',
                             'Earth spatial RGB uses three narrow bands, not resolved full-spectrum colour; historical cloud pattern is illustrative.',
                             'Sea geometry, illumination and reflection are inherited from the existing render, including its sampling uncertainty and uniform-disk Earth reflection.'],
                inputs={str(f):digest(f) for f in [args.epic,args.render,args.build,
                          args.root/'research/studies/sea_appearance/results/nubium-midnight.json',
                          args.root/'visualization/sea-appearance/nubium-earth-inputs.json',
                          args.root/'visualization/sea-appearance/nubium_guides.py',
                          args.root/'visualization/reference-renderer/seas/kernels.py',
                          Path(__file__).with_name('nubium_human_view.py')]})
    (args.out/'nubium-visibility-guide.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(curve=curve,earth_texture=record['earth_texture'],geometry=projection),indent=2))


if __name__=='__main__':main()
