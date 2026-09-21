"""Evaluate fixed-viewpoint spectral panoramas from the solved source field."""
from __future__ import annotations
import json,math,time,argparse
from pathlib import Path
import numpy as np
from numba import njit,prange
from solver import Atmosphere,ray_path,weights_path,interpolate,solar_columns,bracket,ROOT

# Analytic approximation to the CIE 1931 2-degree colour matching functions.
# Wyman, Sloan & Shirley (2013), JCGT 2(2), cited in SOURCES.md.
def cmf(w):
    def g(mu,left,right):
        t=(w-mu)*np.where(w<mu,left,right);return np.exp(-.5*t*t)
    return np.array([1.056*g(599.8,.0264,.0323)+.362*g(442.,.0624,.0374)-.065*g(501.1,.0490,.0382),
      .821*g(568.8,.0213,.0247)+.286*g(530.9,.0613,.0322),
      1.217*g(437.,.0845,.0278)+.681*g(459.,.0385,.0725)]).T
M_XYZ_RGB=np.array([[3.2406,-1.5372,-.4986],[-.9689,1.8758,.0415],[.0557,-.2040,1.0570]])
def colour_weights(lam):
    dw=np.zeros_like(lam);dw[1:-1]=(lam[2:]-lam[:-2])*.5;dw[0]=(lam[1]-lam[0])*.5;dw[-1]=(lam[-1]-lam[-2])*.5
    return 683.*dw[:,None]*cmf(lam)@M_XYZ_RGB.T

def to_rgb(spec,lam):return spec@colour_weights(lam)

@njit(cache=True,inline='always')
def lut(tab,i,j,u,v,l,c):
    return (1-u)*((1-v)*tab[i,j,l,c]+v*tab[i,j+1,l,c])+u*((1-v)*tab[i+1,j,l,c]+v*tab[i+1,j+1,l,c])

@njit(parallel=True,cache=True)
def evaluate(suns,elevs,azs,rg,ag,srg,sag,mom,beam,beta,absorb,R,top,H,ozscale,albedo,nshell=100):
    nl=len(beta);ns=len(suns);ne=len(elevs);na=len(azs)
    out=np.zeros((ns,ne,na,nl))
    shells=R+H*10*np.linspace(0.,1.,nshell+1)**2
    fac=3./(16*math.pi)
    for loopi in prange(ne):
        ie=loopi//2 if loopi%2==0 else ne-1-loopi//2
        mu=math.sin(elevs[ie]);st=math.cos(elevs[ie])
        p,g,ce,se=ray_path(R+1.7,mu,R,top,shells,H,ozscale)
        ww,T=weights_path(p,beta,absorb)
        for iss in range(ns):
            ss=math.sin(suns[iss]);cs=math.cos(suns[iss])
            for ia in range(na):
                cp=math.cos(azs[ia]);nu=ss*mu+cs*st*cp
                if g:
                    mg=max(-1.,min(1.,ss*ce+cs*se*cp));ground_angle=math.asin(mg)
                    for l in range(nl):
                        flux=interpolate(beam,srg,sag,R,ground_angle,l,1)+interpolate(mom,rg,ag,R,ground_angle,l,5)
                        out[iss,ie,ia,l]=T[l]*albedo/math.pi*flux
                for k in range(len(p)):
                    r=p[k,0];qr=p[k,1];ms=max(-1.,min(1.,ss*p[k,2]+cs*p[k,3]*cp));aa=math.asin(ms)
                    qt=(nu-ms*qr)/math.sqrt(max(1e-16,1-ms*ms));qp2=max(0.,1-qr*qr-qt*qt)
                    i,u=bracket(rg,r);j,v=bracket(ag,aa);bi,bu=bracket(srg,r);bj,bv=bracket(sag,aa)
                    for l in range(nl):
                        F=lut(beam,bi,bj,bu,bv,l,0)
                        J=lut(mom,i,j,u,v,l,0);Krr=lut(mom,i,j,u,v,l,1);Ktt=lut(mom,i,j,u,v,l,2);Kpp=lut(mom,i,j,u,v,l,3);Krt=lut(mom,i,j,u,v,l,4)
                        src=F*(1+nu*nu)+J+qr*qr*Krr+qt*qt*Ktt+qp2*Kpp+2*qr*qt*Krt
                        out[iss,ie,ia,l]+=ww[k,l]*fac*max(0.,src)
    return out

def evaluate_file(path,suns,elevs,azs,nshell=100,single=False):
    d=np.load(path);meta=json.loads(str(d['meta']));atm=Atmosphere(**meta['atmosphere'])
    mom=np.zeros_like(d['moments']) if single else d['moments']-d['last_order']
    return evaluate(np.deg2rad(suns),np.deg2rad(elevs),np.deg2rad(azs),d['r'],d['a'],d['sr'],d['sa'],mom,d['beam'],d['beta'],d['absorb'],atm.radius_m,atm.top,atm.scale_height_m,atm.ozone_scale,atm.ground_albedo,nshell)


def spherical_harmonics(x,y,z):
    # Real orthonormal SH, vertical axis z; source scene uses x toward Sun.
    return np.stack([np.ones_like(x)*.282095,.488603*y,.488603*z,.488603*x,
        1.092548*x*y,1.092548*y*z,.315392*(3*z*z-1),1.092548*x*z,.546274*(x*x-y*y)],axis=-1)


def build_atlas(path,ne=81,na=65):
    d=np.load(path);meta=json.loads(str(d['meta']));atm=Atmosphere(**meta['atmosphere'])
    lam=d['lam']; beam=d['beam']; sr=d['sr']; sa=d['sa']
    suns=np.unique(np.r_[np.arange(-90.,-35.999,3.),np.arange(-35.,12.001,.5),np.arange(-1.,1.001,.125),np.arange(13.,30.001,1.),np.arange(32.,90.001,2.)])
    q=np.linspace(-1,1,ne);elevs=np.sign(q)*q*q*90
    azs=np.linspace(0,180,na)
    print('panorama',atm.name,len(suns),ne,na,flush=True);t=time.time()
    spectra=evaluate_file(path,suns,elevs,azs,nshell=64)
    rgb=to_rgb(spectra,lam)
    # Preserve raw signed linear sRGB in archive. Display gamut mapping is separate.
    direct=[];horizontal=[];diffuse=[];sh=[];direct_grid=[];cloud_direct=[];cloud_diffuse=[]
    # Complete an azimuth ring using reflection symmetry, without double-counting ends.
    fullaz=np.linspace(0,2*np.pi,2*(na-1),endpoint=False)
    sky_start=ne//2
    e=np.deg2rad(elevs[sky_start:]);mu=np.sin(e)
    integ_weights=np.zeros(len(e));integ_weights[1:-1]=(mu[2:]-mu[:-2])/2;integ_weights[0]=(mu[1]-mu[0])/2;integ_weights[-1]=(mu[-1]-mu[-2])/2
    X=np.cos(e)[:,None]*np.cos(fullaz)[None,:];Y=np.cos(e)[:,None]*np.sin(fullaz)[None,:];Z=np.sin(e)[:,None]+np.zeros_like(X)
    basis=spherical_harmonics(X,Y,Z);dw=integ_weights[:,None]*2*np.pi/len(fullaz)
    for k,sun in enumerate(suns):
        tab=np.zeros((len(lam),2))
        for l in range(len(lam)):
            for c in range(2):tab[l,c]=interpolate(beam,sr,sa,atm.radius_m+1.7,np.deg2rad(sun),l,c)
        dd=to_rgb(tab[:,0],lam);dh=to_rgb(tab[:,1],lam)
        # Prescribed optical weather layer at 2500 m: lighting inputs, not weather prediction.
        cc=np.zeros((len(lam),2))
        for l in range(len(lam)):
            cc[l,0]=interpolate(beam,sr,sa,atm.radius_m+2500,np.deg2rad(sun),l,0)
            cc[l,1]=interpolate(d['moments'],d['r'],d['a'],atm.radius_m+2500,np.deg2rad(sun),l,5)
        cloud_direct.append(to_rgb(cc[:,0],lam).tolist());cloud_diffuse.append(to_rgb(cc[:,1],lam).tolist())
        ring=np.concatenate([rgb[k,sky_start:],rgb[k,sky_start:, -2:0:-1]],axis=1)
        flux=np.sum(ring*(Z*dw)[...,None],axis=(0,1))
        coeff=np.einsum('eac,eak,ea->kc',ring,basis,dw+np.zeros_like(X))
        # Ground radiance in the lower hemisphere is approximated as homogeneous.
        gl=atm.ground_albedo*(flux+dh)/np.pi
        baselow=spherical_harmonics(X,Y,-Z)
        coeff+=np.einsum('k,c->kc',np.sum(baselow*dw[...,None],axis=(0,1)),gl)
        # Lambertian irradiance convolution.
        coeff*=np.array([np.pi]+[2*np.pi/3]*3+[np.pi/4]*5)[:,None]
        direct.append(dd.tolist());horizontal.append(dh.tolist());diffuse.append(flux.tolist());sh.append(coeff.tolist())
    dest=ROOT/'data';dest.mkdir(exist_ok=True)
    np.savez_compressed(dest/f'{atm.name}_atlas.npz',rgb=rgb.astype(np.float32),suns=suns,elevs=elevs,azs=azs,lam=d['lam'],direct=direct,direct_horizontal=horizontal,diffuse=diffuse,sh=sh,
      cloud_direct=cloud_direct,cloud_diffuse=cloud_diffuse,metadata=json.dumps(meta))
    print('done',atm.name,round(time.time()-t,1),'seconds',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('paths',nargs='+');args=p.parse_args()
    for path in args.paths:build_atlas(path)
