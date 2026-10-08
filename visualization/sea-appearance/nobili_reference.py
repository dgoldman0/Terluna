"""Render the admitted Nobili scenario using the existing spectral sea tracer.

The scene's dated geometry, native terrain and experimental wind are read from
its research product. The equilibrium lake-wave realization is conditional;
optional cloud and fetch products constrain its angular sky and wave age.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sys,time
from pathlib import Path
import numpy as np


def wave_state(product,seed):
    import scene
    from illumination.water_surface import realization as r
    from illumination.water_surface.full_spectrum import unified_density
    from illumination.water_surface.short_waves import phase_speed
    from shared.constants import MOON_SURFACE_GRAVITY as gravity
    w=product['wind'];u=w['speed_10m_m_s'];ust=w['friction_velocity_m_s'];toward=w['toward_cartesian_deg'];omega=product['selected_inverse_wave_age']
    depth=product['terrain']['water_depth_at_observer_m']
    peak_wavelength=2*np.pi*u*u/(gravity*omega*omega)
    outer_size=max(256.,float(2**math.ceil(math.log2(8*peak_wavelength))))
    tiles=r.realize([unified_density(u,ust,omega,toward,gravity)],depth,gravity,seed,
        cascades=((outer_size,1024),(17.3,512),(1.43,512)))
    peak_speed=float(phase_speed(np.array([gravity*omega**2/u**2]),gravity)[0])
    residual=r.residual_covariance(ust,toward,gravity,tiles[-1].k_high,peak_speed)
    pyramids=[r.lean_pyramid(t.slope) for t in tiles];levels=max(map(len,pyramids))
    offsets=np.zeros((len(tiles),levels),np.int64);sizes=np.ones_like(offsets);chunks=[];position=0
    for i,pyramid in enumerate(pyramids):
        for j,level in enumerate(pyramid):
            offsets[i,j]=position;sizes[i,j]=level.shape[1];chunks.append(level.ravel());position+=level.size
    heights=tiles[0].height;n=heights.shape[0];b=scene.SEA_BLOCK
    blocks=heights.reshape(n//b,b,n//b,b).max(axis=(1,3));first=heights[::b,::b]
    blocks=np.maximum.reduce([blocks,np.roll(first,-1,0),np.roll(first,-1,1),np.roll(np.roll(first,-1,0),-1,1),
        np.roll(heights[::b,:].reshape(n//b,n//b,b).max(axis=2),-1,0),
        np.roll(heights[:,::b].reshape(n//b,b,n//b).max(axis=1),-1,1)]).astype('f4')
    record=dict(spectrum='Conditional full-range ECKV equilibrium spectrum at saved CM1 wind/stress; corrected equation 41.',
        inverse_wave_age=omega,depth_m=depth,wind_10m_m_s=u,friction_velocity_m_s=ust,wind_toward_deg=toward,
        significant_height_m=float(4*np.std(heights)),expected_significant_height_m=float(4*np.sqrt(sum(t.expected_variance for t in tiles))),
        mean_square_slope=float(sum(np.trace(r.sample_covariance(t)) for t in tiles)+np.trace(residual)),seed=seed,
        tiles=[dict(size_m=t.size_m,cells=t.cells,height_variance_expected=t.expected_variance,
            slope_covariance_expected=t.expected_covariance.tolist(),slope_covariance_drawn=r.sample_covariance(t).tolist()) for t in tiles])
    return dict(heights=heights,blocks=blocks,hmin=float(heights.min()),hmax=float(heights.max()),dx=tiles[0].spacing_m,
        lean=np.concatenate(chunks).astype('f4'),offsets=offsets,sizes=sizes,spacing=np.array([t.spacing_m for t in tiles]),
        levels=np.array(list(map(len,pyramids)),np.int64),residual=np.array([residual[0,0],residual[1,1],residual[0,1]]),record=record)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--scene',type=Path);p.add_argument('--waves',type=Path);p.add_argument('--ground',type=Path);p.add_argument('--cloud-smoothing-deg',type=float,default=0.);p.add_argument('--cloud',type=Path);p.add_argument('--width',type=int,default=1280);p.add_argument('--spp-side',type=int,default=2);a=p.parse_args()
    root=(a.root or Path(__file__).resolve().parents[2]).resolve();sys.path.insert(0,str(root));sys.path.insert(0,str(root/'visualization/reference-renderer/seas'))
    import scene,render
    source=(a.scene or root/'research/studies/sea_appearance/results/nobili-night.json').absolute();product=json.loads(source.read_text());s=product['camera'];h=product['earth']
    if a.waves:
        waves=json.loads(a.waves.read_text())
        if waves['scene_sha256']!=scene.sha256(source):raise ValueError('Wave product belongs to a different scene')
        product['selected_inverse_wave_age']=waves['selected_inverse_wave_age']
        product['wave_selection']=waves
    coast_name=product.get('name','Nobili')
    a.out.mkdir(parents=True,exist_ok=True);common=scene.Common();common.jd=lambda hour:product['jd_tt'];common.tide=lambda *args:0.
    common.regimes['scenes'].append(dict(coast=coast_name,moment='night',hour=0.))
    common.calendar['sites'].append(dict(short=coast_name,lon_lat=product['site_lon_lat_deg']))
    scene.CAMERAS[coast_name]=dict(offset_m=(0.,0.),note=s['note']);scene.VIEWS[(coast_name,'night')]=(s['view_azimuth_deg'],s['pitch_deg'])
    scene.TERRAIN=(root/product['terrain']['file']).parent;scene.TERRAIN_NAMES[coast_name]='terrain';scene.LAND_SOIL[coast_name]='62231'
    # Shading is sampled every 400 m; the native terrain geometry remains 100 m.
    scene.TEXTURE_STRIDE=4
    scene.SKY_EL=np.unique(np.r_[np.arange(-.12,.1,.01),np.arange(.1,2,.1),np.arange(2,20,.4),np.arange(20,90.01,1.)])
    scene.SKY_REL_AZ=np.unique(np.r_[np.arange(0,5,.25),np.arange(5,20,1.),np.arange(20,180.01,3.)])
    original_geometry=scene.lighting.geometry
    def geometry(jd,lon,lat):
        g=original_geometry(jd,lon,lat)
        for key,value in dict(earth_elevation=h['elevation_deg'],earth_azimuth=h['azimuth_deg'],earth_phase_angle=h['phase_angle_deg'],
            earth_lit_fraction=h['lit_fraction'],earth_distance_m=product['light']['earth_distance_m']).items():g[key]=np.array([value])
        return g
    scene.lighting.geometry=geometry
    scene.sea_state=lambda common,coast,hour,seed:wave_state(product,seed)
    for path in [source,Path(__file__),root/'research/studies/sea_appearance/nobili.py',root/'illumination/water_surface/short_waves.py',root/'illumination/water_surface/realization.py',root/'illumination/water_surface/full_spectrum.py']:
        common.inputs[str(path.relative_to(root))]=scene.sha256(path)
    if a.waves:common.inputs[str(a.waves)]=scene.sha256(a.waves)
    if a.ground:common.inputs[str(a.ground)]=scene.sha256(a.ground)
    if a.cloud:
        from scipy.ndimage import gaussian_filter
        from scipy.interpolate import RegularGridInterpolator
        with np.load(a.cloud) as z:
            ce=z['elevation_deg'];ca=z['azimuth_deg'];cloud=z['cloudy_groups'].mean(axis=2);clear=z['clear_groups'].mean(axis=2);cloud_meta=json.loads(str(z['metadata']))
        if cloud_meta['input_scene_sha256']!=scene.sha256(source):raise ValueError('Cloud sky belongs to a different scene product')
        from cloud_display import delta_interpolator
        delta_fn=delta_interpolator(ce,ca,cloud,clear,a.cloud_smoothing_deg)
        common.inputs['visualization/sea-appearance/cloud_display.py']=scene.sha256(root/'visualization/sea-appearance/cloud_display.py')
        ground_record=None
        if a.ground:
            ground_record=json.loads(a.ground.read_text())
            if ground_record['input_scene_sha256']!=scene.sha256(source):raise ValueError('Ground light belongs to a different scene')
        old_light=scene.light
        def cloudy_light(common,geometry,weights):
            light=old_light(common,geometry,weights)
            ee,aa_=np.meshgrid(scene.SKY_EL,light['azimuths'],indexing='ij');points=np.stack([np.clip(ee,ce[0],ce[-1]),aa_],-1)
            light['sky']=np.maximum(light['sky']+delta_fn(points),0)
            # The terrain uses the angular cloudy/clear diffuse ratio with its
            # existing band spectra; local cloud shadows at the rim are open.
            ee,aa_=np.meshgrid(light['coarse_el'],light['coarse_az'],indexing='ij')
            baseline_y=light['sky_bands'][...,1].sum(axis=-1)
            change_y=delta_fn(np.stack([np.clip(ee,ce[0],ce[-1]),aa_],-1))[...,1]
            ratios=np.clip((baseline_y+change_y)/np.maximum(baseline_y,1e-12),0,4)
            light['sky_bands']*=ratios[:,:,None,None]
            trans=cloud_meta['disk_cloud_transmission_centre']
            light['earth'][7:10]*=trans;light['earth_bands']*=trans
            cloudy_lux=cloud_meta.get('total_horizontal_lux',cloud_meta['earth_total_horizontal_lux'])
            if ground_record:cloudy_lux=ground_record['total_horizontal_lux']
            clear_lux=cloud_meta.get('clear_total_horizontal_lux',product['light']['earth_lux'])
            light['water']*=cloudy_lux/clear_lux
            return light
        scene.light=cloudy_light
        common.inputs[str(a.cloud.relative_to(root)) if a.cloud.is_absolute() else str(a.cloud)]=scene.sha256(a.cloud)
    original_build=scene.build;cache=a.out/'scene-build.npz'
    ident=hashlib.sha256(json.dumps(dict(product=product,width=a.width,producer=scene.sha256(__file__),scene=scene.sha256(Path(scene.__file__)),cloud=scene.sha256(a.cloud) if a.cloud else None,smoothing_deg=a.cloud_smoothing_deg,ground=scene.sha256(a.ground) if a.ground else None),sort_keys=True).encode()).hexdigest()
    def cached(*args,**kwargs):
        if cache.exists():
            with np.load(cache) as z:
                if str(z['identity'])!=ident:raise ValueError('Scene cache identity differs; choose a fresh build directory')
                return z['camera'],{k[6:]:z[k] for k in z.files if k.startswith('array/')},json.loads(str(z['record'])),json.loads(str(z['geometry'])),{k[6:]:z[k] for k in z.files if k.startswith('light/')}
        for name in ['light','sea_state','coast']:
            fn=getattr(scene,name)
            def timed(*aa,_fn=fn,_name=name,**kw):
                print('Preparing',_name,flush=True);t=time.monotonic();r=_fn(*aa,**kw);
                if _name=='light':np.savez(a.out/'light-cache.npz',**r)
                print(_name,'ready',round(time.monotonic()-t,1),'s',flush=True);return r
            setattr(scene,name,timed)
        camera,arrays,record,g,light=original_build(*args,**kwargs)
        record['scenario_source_sha256']=scene.sha256(source)
        if a.waves:record['wave_selection']=waves
        record['sky']='Cloud-conditioned reference' if a.cloud else 'Clear molecular reference; lunar cloud optics are separate.'
        if a.cloud:record['cloud_conditioning']=dict(source_sha256=scene.sha256(a.cloud),disk_centre_transmission=cloud_meta['disk_cloud_transmission_centre'],angular_smoothing_deg=a.cloud_smoothing_deg,ground_refinement=ground_record,reading_rule='The tracer reflects the conditional cloudy angular sky and attenuated finite Earth. Terrain diffuse spectral shape and aerial perspective retain clear-reference approximations; cloud attenuation uses the observer value, not a terrain-coupled shadow solve.')
        np.savez(cache,identity=ident,camera=camera,record=json.dumps(record),geometry=json.dumps(g),
            **{'array/'+k:v for k,v in arrays.items()},**{'light/'+k:v for k,v in light.items()})
        return camera,arrays,record,g,light
    scene.build=cached
    # New output directories prevent accidental reuse of row blocks at a new sample count.
    work=a.out/'partial'/f"{coast_name.replace(' ', '_').lower()}__night"
    if work.exists() and any(work.glob('rows_*.npz')):
        raise ValueError('Retained partial rows: inspect before restarting this sample configuration')
    render.render_frame(common,coast_name,'night',a.width,round(a.width*s['frame'][1]/s['frame'][0]),a.spp_side,20330621,s['horizontal_fov_deg'],a.out)

if __name__=='__main__':main()
