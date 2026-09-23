"""Offline B1 scene frames, spherical hemispheres and explicit source identities."""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,json,math,os,time
from pathlib import Path
import numpy as np
from PIL import Image
import geometry as g
import scene_transport as b
import primary
ROOT=b.ROOT

SCENARIOS={
 'clear':dict(label='Clear air',height_shift_m=0.,extinction_factor=0.),
 'fair':dict(label='Broken cloud · 9.3–12.5 km',height_shift_m=0.,extinction_factor=1.),
 'high':dict(label='Raised same field · 29.3–32.5 km',height_shift_m=20000.,extinction_factor=1.),
 'thin':dict(label='Thin same field · 35% extinction',height_shift_m=0.,extinction_factor=.35),
 'thick':dict(label='Thick same field · 2× extinction',height_shift_m=0.,extinction_factor=2.),
}

def spectra(groups=12):
    s=json.loads(b.locate('data/a1/spectral-inputs.json').read_text())
    c=b.t.load_inputs();w=np.full(48,10.);w[[0,-1]]=5.
    W=683*np.array(s['cie1931_xyz'])*w[:,None]*c['solar'][:,None]
    blocks=np.array_split(np.arange(48),groups)
    # Integrate the supplied solar spectrum and CIE weights exactly in each
    # group; evaluate the unit-irradiance transfer at its central wavelength.
    lam=np.array([np.mean(c['lam'][ix]) for ix in blocks])
    weights=np.array([W[ix].sum(axis=0) for ix in blocks])
    beta=1.24062e-6*(lam/1000.)**-4
    sigma=np.interp(lam,c['lam'],c['sigma'])
    return c,lam,beta,sigma,weights,np.array(s['xyz_to_linear_srgb'])

def prepare(scenario,groups=12):
    c,lam,beta,sigma,W,M=spectra(groups);s=SCENARIOS[scenario]
    c['lo']=c['lo'].copy();c['hi']=c['hi'].copy()
    c['lo'][1]+=s['height_shift_m'];c['hi'][1]+=s['height_shift_m']
    c['ext']=c['ext']*s['extinction_factor'];c['cmax']=float(c['ext'].max())*(1+1e-9)
    verts,glo,ghi,step=g.make_coast(c['R'])
    return c,lam,beta,sigma,W,M,verts,glo,ghi,step

def raw_identity():
    inherited=['tools/a3/coupled_transport.py','tools/a1/cloud_reference.py','tools/a2/molecular_multiple_scattering.py','data/a1/spectral-inputs.json']
    own=['tools/b1/geometry.py','tools/b1/scene_transport.py','tools/b1/render_scene.py','tools/b1/primary.py']
    return {x:hashlib.sha256(b.locate(x).read_bytes()).hexdigest() for x in inherited+own}

def display(rgb,ev=0.,scale=10000.):
    v=1-np.exp(-np.maximum(0.,rgb)*2**ev/scale)
    srgb=np.where(v<=.0031308,12.92*v,1.055*v**(1/2.4)-.055)
    return np.uint8(np.clip(srgb*255+.5,0,255))

def filtered(rgb,gbuffer,rad=2):
    """Optional spatial presentation filter; raw radiance is always retained.
    Terrain class and relative surface depth prevent shoreline edge mixing.
    Clouds have no resolved G-buffer; their edges can soften, explicitly.
    """
    out=np.zeros_like(rgb);ws=np.zeros(rgb.shape[:2]);h,w=ws.shape
    material=gbuffer[:,:,1];depth=np.minimum(gbuffer[:,:,0],1e7)
    for y in range(-rad,rad+1):
        for x in range(-rad,rad+1):
            val=np.roll(rgb,(y,x),(0,1));mt=np.roll(material,(y,x),(0,1));dp=np.roll(depth,(y,x),(0,1))
            a=np.exp(-(x*x+y*y)/max(1.,rad*rad))*(mt==material)*np.exp(-abs(dp-depth)/np.maximum(30.,depth*.12))
            if y>0:a[:y]=0
            if y<0:a[y:]=0
            if x>0:a[:,:x]=0
            if x<0:a[:,x:]=0
            out+=a[:,:,None]*val;ws+=a
    return out/np.maximum(ws[:,:,None],1e-30)

def render(scenario,sun_elevation,width,height,spp,workers,outdir,projection='perspective',groups=12,seed=562399,black=False):
    c,lam,beta,sigma,W,M,verts,glo,ghi,step=prepare(scenario,groups)
    origin=np.array([-2200.,5000.,-6500.]);hit,n,mat=g.surface(origin,np.array([0.,-1.,0.]),c['R'],verts,glo,ghi,step)
    origin[1]-=hit;origin+=n*2.
    az=math.radians(25.)
    a=math.radians(sun_elevation);sun=np.array([math.sin(az)*math.cos(a),math.sin(a),math.cos(az)*math.cos(a)])
    radius=math.radians(.26667);clouds=scenario!='clear'
    fov,yaw,pitch=104.,10.,16.
    ident=dict(schema='open-moon-b1-frame/1',scenario=scenario,scenario_inputs=SCENARIOS[scenario],sun_elevation_deg=sun_elevation,sun_azimuth_deg=25.,sun_radius_deg=.26667,origin_m=origin.tolist(),camera=dict(projection=projection,horizontal_fov_deg=fov,yaw_deg=yaw,pitch_deg=pitch),width=width,height=height,samples_per_pixel_per_band=spp,wavelength_groups=groups,representative_wavelength_nm=lam.tolist(),profile_sha256=c['profile_sha256'],base_cloud_sha256=c['field_sha256'],cloud_bounds_min_m=c['lo'].tolist(),cloud_bounds_max_m=c['hi'].tolist(),cloud_extinction_sha256=hashlib.sha256(c['ext'].tobytes()).hexdigest(),geometry_sha256=hashlib.sha256(verts.tobytes()).hexdigest(),seed=seed,surface_boundary='black control' if black else 'authored spectral Lambertian terrain and Fresnel-reflecting deep water',sources=raw_identity())
    identity=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest();name=f'{scenario}_sun{sun_elevation:+g}_{projection}'+('_black' if black else '')
    dest=outdir/(name+'.npz');meta=outdir/(name+'.json')
    if dest.exists() and meta.exists() and json.loads(meta.read_text()).get('identity_sha256')==identity:
        print('Reuse',name,flush=True);return meta
    start=time.monotonic();outdir.mkdir(parents=True,exist_ok=True)
    def args(ja,jb):return (origin,width,height,ja,jb,fov,yaw,pitch,0 if projection=='perspective' else 1,sun,radius,lam,beta,sigma,W,c['R'],c['top'],c['rows'],c['oz'],c['edges'],c['air_bounds'],c['oz_bounds'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['cmax'],verts,glo,ghi,step,clouds,True,0. if black else 1.,not black,spp,seed)
    # Compile on an empty tile before worker dispatch.
    b.pixel_batch(*args(0,0))
    spec=np.zeros((height,width,len(lam)));specse=np.zeros_like(spec);se=np.zeros((height,width,3));sur=np.zeros_like(se);fail=0
    chunks=[(j,min(height,j+3)) for j in range(0,height,3)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for (ja,jb),r in zip(chunks,pool.map(lambda z:b.pixel_batch(*args(*z)),chunks)):
            spec[ja:jb],specse[ja:jb],se[ja:jb],sur[ja:jb],f=r;fail+=f
    disk,disk_surface=primary.primary_frame(origin,width,height,fov,yaw,pitch,0 if projection=='perspective' else 1,sun,radius,beta,sigma,c['R'],c['top'],c['rows'],c['oz'],c['shells'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],verts,glo,ghi,step,not black,8)
    scattered_rgb=(spec@W)@M.T
    spec+=disk;sur+=disk_surface@W
    XYZ=spec@W;rgb=XYZ@M.T
    if projection=='perspective':gb=g.camera_maps(origin,width,height,c['R'],verts,glo,ghi,step,fov,yaw,pitch)
    else:gb=np.zeros((height,width,5));gb[:,:,0]=1e30
    sm=filtered(scattered_rgb,gb,3)+(disk@W)@M.T
    np.savez_compressed(dest,unit_transfer=spec,unit_transfer_se=specse,wavelength_nm=lam,solar_CIE_XYZ_weights=W,XYZ=XYZ,linear_srgb=rgb,XYZ_standard_error=se,surface_path_XYZ=sur,geometry=gb,filtered_linear_srgb=sm,primary_XYZ=disk@W)
    Image.fromarray(display(sm)).resize((width*4,height*4),Image.Resampling.BILINEAR).save(outdir/(name+'.png'))
    Image.fromarray(display(rgb)).save(outdir/(name+'_raw.png'))
    ident.update(identity_sha256=identity,unresolved_paths=int(fail),elapsed_s=time.monotonic()-start,total_paths=int(width*height*spp*len(lam)),radiance_units='XYZ photopic; Y cd m^-2; unit_transfer is sr^-1 at representative wavelengths',display=dict(scale=10000.,ev=0.,filter='optional radius-3 material/depth-gated spatial average in linear RGB; deterministic primary solar disk excluded',white_balance='none',gamut='negative linear sRGB channels clipped for display only'),sampling=dict(median_Y_standard_error=float(np.median(se[:,:,1])),median_Y=float(np.median(XYZ[:,:,1]))),outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [dest,outdir/(name+'.png'),outdir/(name+'_raw.png')]})
    meta.write_text(json.dumps(ident,indent=2)+'\n');print(json.dumps(dict(frame=name,seconds=ident['elapsed_s'],unresolved=fail,median_Y=ident['sampling']['median_Y'])),flush=True)
    return meta

def main():
    a=argparse.ArgumentParser();a.add_argument('--scenarios',nargs='+',choices=list(SCENARIOS),default=['clear','fair']);a.add_argument('--suns',nargs='+',type=float,default=[60,30,12,6,2,0,-2,-6]);a.add_argument('--width',type=int,default=192);a.add_argument('--height',type=int,default=112);a.add_argument('--spp',type=int,default=24);a.add_argument('--workers',type=int,default=4);a.add_argument('--groups',type=int,choices=[12,24,48],default=12);a.add_argument('--projection',choices=['perspective','hemisphere'],default='perspective');a.add_argument('--black',action='store_true');a.add_argument('--seed',type=int,default=562399);a.add_argument('--output',type=Path,default=ROOT/'data/b1/generated/frames');x=a.parse_args()
    if min(x.width,x.height,x.spp,x.workers)<1: a.error('Positive dimensions, samples and workers required')
    for s in x.scenarios:
        for sun in x.suns:render(s,sun,x.width,x.height,x.spp,x.workers,x.output,x.projection,x.groups,x.seed,x.black)
if __name__=='__main__':main()
