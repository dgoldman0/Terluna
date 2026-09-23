"""Ground incident irradiance and direct-beam visibility from the B1 scene.

Direct illumination uses deterministic gas/cloud integrals and finite-disk
quadrature. Diffuse irradiance uses cosine-weighted B1 transport histories,
including subsequent water, terrain, gas and cloud interactions. It excludes
the unscattered camera-to-Sun term, which is recorded separately as direct.
"""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,json,math,time
from pathlib import Path
import numpy as np
from numba import njit
import render_scene as r
import scene_transport as b
import geometry as g
import deterministic as dt

@njit(cache=True,nogil=True)
def tile(j0,j1,N,extent,sun,lam,beta,sigma,W,R,top,rows,oz,shells,edges,airb,ozb,lo,hi,dims,ext,ice,cmax,V,glo,ghi,step,cloud_on,spp,seed):
    direct=np.zeros((j1-j0,N,3));diffuse=np.zeros_like(direct);stderr=np.zeros_like(direct);visibility=np.zeros((j1-j0,N));positions=np.zeros_like(direct);normals=np.zeros_like(direct);failures=0
    nodes=np.array([-.9602898564975363,-.7966664774136267,-.525532409916329,-.1834346424956498,.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363]);weights=np.array([.1012285362903763,.2223810344533745,.3137066458778873,.362683783378362,.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763])
    radius=math.radians(.26667)
    for j in range(j0,j1):
        for i in range(N):
            p=np.array([-extent+(i+.5)*2*extent/N,5000.,-extent+(j+.5)*2*extent/N]);d=np.array([0.,-1.,0.]);h,n,mat=g.surface(p,d,R,V,glo,ghi,step)
            p=p+d*h+n*.04;positions[j-j0,i]=p;normals[j-j0,i]=n
            # Eight equal-solid-angle disk samples. This deterministic quadrature
            # has no Monte Carlo variance; its discretization error remains open.
            for k in range(8):
                sd=b.t.rotate(sun,1-(k+.5)/8*(1-math.cos(radius)),k*2.399963229728653)
                cs=max(0.,np.dot(n,sd))
                if cs<=0 or b.terrain_blocked(p,sd,R,V,glo,ghi,step,True):continue
                tau=dt.cloud_depth(p,sd,lo,hi,dims,ext,ice) if cloud_on else 0.
                ct=math.exp(-tau);visibility[j-j0,i]+=ct/8.
                air,ozone,blocked=dt.gas_columns(p,sd,R,top,rows,oz,shells,nodes,weights)
                if blocked:continue
                for kband in range(len(lam)):
                    tr=math.exp(-beta[kband]*air-sigma[kband]*b.t.DU*ozone)*ct
                    direct[j-j0,i]+=W[kband]*cs*tr*(2/(1+math.cos(radius)))/8.
            sx=np.zeros(3);sx2=np.zeros(3)
            for sample in range(spp):
                ss=(seed+104729*(j*N+i)+196613*sample)%2147483647;np.random.seed(ss);direction=b.cosine_direction(n);X=np.zeros(3)
                for k in range(len(lam)):
                    np.random.seed(ss+17)
                    value,sv,f=b.path_sample(p,direction,sun,radius,lam[k],R,top,rows,oz,beta[k],sigma[k],edges,airb,ozb,lo,hi,dims,ext,ice,cmax,V,glo,ghi,step,True,cloud_on,True,1.,True,.95,0.,4096,False)
                    X+=math.pi*value*W[k];failures+=int(f!=0)
                sx+=X;sx2+=X*X
            diffuse[j-j0,i]=sx/spp;stderr[j-j0,i]=np.sqrt(np.maximum(0.,(sx2-sx*sx/spp)/(max(1,spp-1)*spp)))
    return direct,diffuse,stderr,visibility,positions,normals,failures

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--scenarios',nargs='+',default=['clear','fair']);ap.add_argument('--suns',nargs='+',type=float,default=[60,30,12,6,0,-2]);ap.add_argument('--size',type=int,default=32);ap.add_argument('--spp',type=int,default=16);ap.add_argument('--workers',type=int,default=4);a=ap.parse_args();dest=b.ROOT/'data/b1/generated/maps';dest.mkdir(parents=True,exist_ok=True)
    for case in a.scenarios:
        c,lam,beta,sigma,W,M,V,glo,ghi,step=r.prepare(case)
        for s in a.suns:
            t0=time.monotonic();su=math.radians(s);az=math.radians(25.);sun=np.array([math.sin(az)*math.cos(su),math.sin(su),math.cos(az)*math.cos(su)])
            def args(j0,j1):return(j0,j1,a.size,18000.,sun,lam,beta,sigma,W,c['R'],c['top'],c['rows'],c['oz'],c['shells'],c['edges'],c['air_bounds'],c['oz_bounds'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['cmax'],V,glo,ghi,step,case!='clear',a.spp,726451)
            tile(*args(0,0));arrays=[np.zeros((a.size,a.size,3)) for _ in range(6)];vis=np.zeros((a.size,a.size));fail=0
            with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
                chunks=[(j,min(a.size,j+2)) for j in range(0,a.size,2)]
                for (j0,j1),values in zip(chunks,pool.map(lambda z:tile(*args(*z)),chunks)):
                    di,df,se,vi,po,no,ff=values;arrays[0][j0:j1]=di;arrays[1][j0:j1]=df;arrays[2][j0:j1]=se;arrays[3][j0:j1]=po;arrays[4][j0:j1]=no;vis[j0:j1]=vi;fail+=ff
            f=dest/f'{case}_sun{s:+g}.npz'
            np.savez_compressed(f,direct_XYZ=arrays[0],diffuse_XYZ=arrays[1],total_XYZ=arrays[0]+arrays[1],diffuse_XYZ_standard_error=arrays[2],position_m=arrays[3],normal=arrays[4],direct_beam_visibility=vis)
            meta=dict(schema='open-moon-b1-ground-irradiance/1',scenario=case,sun_elevation_deg=s,extent_m=18000.,size=a.size,samples_per_band=a.spp,spectral_groups=12,solar_disk_quadrature=8,units='incident XYZ; Y lux',visibility='cloud + terrain/horizon visibility; gas attenuation excluded from this diagnostic',unresolved_paths=fail,elapsed_s=time.monotonic()-t0,profile_sha256=c['profile_sha256'],base_cloud_sha256=c['field_sha256'],cloud_extinction_sha256=hashlib.sha256(c['ext'].tobytes()).hexdigest(),geometry_sha256=hashlib.sha256(V.tobytes()).hexdigest(),cloud_bounds_min_m=c['lo'].tolist(),cloud_bounds_max_m=c['hi'].tolist(),representative_wavelength_nm=lam.tolist(),surface_boundary='authored spectral Lambertian terrain and Fresnel-reflecting deep water',source_hashes={**r.raw_identity(),'tools/b1/diagnostics.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},output_sha256=hashlib.sha256(f.read_bytes()).hexdigest())
            f.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps({k:meta[k] for k in ['scenario','sun_elevation_deg','elapsed_s','unresolved_paths']}),flush=True)
if __name__=='__main__':main()
