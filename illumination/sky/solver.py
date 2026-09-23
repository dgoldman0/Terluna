"""Spherical scalar Rayleigh transport by successive scattering orders.

A research prototype, not a calibrated lunar appearance model.
Straight rays, prescribed exponential gas and triangular ozone profiles,
Lambertian ground, finite solar disk. All transport lengths are in metres.

The exact scalar Rayleigh angular source is evaluated from its zeroth and
second angular moments: 3/(16*pi) * (J + d_i K_ij d_j). Axial symmetry about
the Sun eliminates two off-diagonal moments, without isotropising the source.
Each ray is integrated through spherical density shells. Angular quadrature
and bilinear interpolation remain numerical approximations, tested separately.
"""
from __future__ import annotations
import argparse, json, math, time
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
from numba import njit, prange

ROOT=Path(__file__).resolve().parent
PI=np.pi
# Binned ASTM G173 extraterrestrial irradiance and 233 K ozone cross sections,
# transcribed from Eric Bruneton's BSD-3-Clause 2017/2018 reference demo.
# Source details and complete upstream notice are in SOURCES.md/LICENSES.txt.
SW=np.arange(360.,831.,10.)
SF=np.array([1.11776,1.14259,1.01249,1.14716,1.72765,1.73054,1.6887,1.61253,1.91198,2.03474,2.02042,2.02212,1.93377,1.95809,1.91686,1.8298,1.8685,1.8931,1.85149,1.8504,1.8341,1.8345,1.8147,1.78158,1.7533,1.6965,1.68194,1.64654,1.6048,1.52143,1.55622,1.5113,1.474,1.4482,1.41018,1.36775,1.34188,1.31429,1.28303,1.26758,1.2367,1.2082,1.18737,1.14683,1.12362,1.1058,1.07124,1.04992])
O3=np.array([1.18e-27,2.182e-28,2.818e-28,6.636e-28,1.527e-27,2.763e-27,5.52e-27,8.451e-27,1.582e-26,2.316e-26,3.669e-26,4.924e-26,7.752e-26,9.016e-26,1.48e-25,1.602e-25,2.139e-25,2.755e-25,3.091e-25,3.5e-25,4.266e-25,4.672e-25,4.398e-25,4.701e-25,5.019e-25,4.305e-25,3.74e-25,3.215e-25,2.662e-25,2.238e-25,1.852e-25,1.473e-25,1.209e-25,9.423e-26,7.455e-26,6.566e-26,5.105e-26,4.15e-26,4.228e-26,3.237e-26,2.451e-26,2.801e-26,2.534e-26,1.624e-26,1.465e-26,2.078e-26,1.383e-26,7.105e-27])

@dataclass(frozen=True)
class Atmosphere:
    name:str
    radius_m:float
    scale_height_m:float
    density_scale:float
    ozone_du:float=300.
    ground_albedo:float=.1
    top_scale_heights:float=10.
    sun_radius_deg:float=32./120.
    visible_filter:float=1.
    @property
    def top(self): return self.radius_m+self.scale_height_m*self.top_scale_heights
    @property
    def ozone_scale(self): return self.scale_height_m/8000.

EARTH=Atmosphere('earth',6371000.,8000.,1.)
# Optical profile proxy: terrestrial effective scale height stretched by g_E/g_M.
# At the same profile temperature, 1.2x surface number density gives 7.26x column.
MOON=Atmosphere('moon',1737400.,8000.*9.80665/1.62,1.2)
MOON_ZERO=Atmosphere('moon_no_ozone',MOON.radius_m,MOON.scale_height_m,1.2,0.)

@njit(cache=True)
def bracket(grid,x):
    j=np.searchsorted(grid,x)-1
    j=max(0,min(j,len(grid)-2))
    return j,max(0.,min(1.,(x-grid[j])/(grid[j+1]-grid[j])))

@njit(cache=True)
def densities(h,H,ozscale):
    gas=math.exp(-max(h,0.)/H)
    q=h/ozscale
    ozone=max(0.,min((q-10000.)/15000.,(40000.-q)/15000.))
    return gas,ozone

@njit(cache=True)
def ray_path(r,mu,R,top,shells,H,ozscale):
    """Mid-shell formal integration; split every radial shell crossing."""
    disc=r*r*(mu*mu-1.)+R*R
    ground=mu<0. and disc>0.
    if ground: end=-r*mu-math.sqrt(disc)
    else: end=-r*mu+math.sqrt(max(0.,r*r*(mu*mu-1.)+top*top))
    end=max(0.,end)
    ts=np.empty(2*len(shells)+2);n=1;ts[0]=0.
    for rad in shells:
        dd=r*r*(mu*mu-1.)+rad*rad
        if dd>0.:
            rt=math.sqrt(dd)
            t1=-r*mu-rt;t2=-r*mu+rt
            if t1>1e-7 and t1<end-1e-7:ts[n]=t1;n+=1
            if t2>1e-7 and t2<end-1e-7:ts[n]=t2;n+=1
    ts[n]=end;n+=1;ts=np.sort(ts[:n]);count=n-1
    out=np.zeros((count,7))
    for k in range(count):
        t=(ts[k]+ts[k+1])*.5
        rad=math.sqrt(r*r+2*r*mu*t+t*t)
        gas,oz=densities(rad-R,H,ozscale)
        out[k,0]=rad
        out[k,1]=(r*mu+t)/rad # view direction dot local outward normal
        out[k,2]=(r+mu*t)/rad # rotation cosine of normals
        out[k,3]=t*math.sqrt(max(0.,1-mu*mu))/rad # rotation sine
        out[k,4]=gas
        out[k,5]=oz
        out[k,6]=ts[k+1]-ts[k]
    cend=(r+mu*end)/(R if ground else top)
    send=end*math.sqrt(max(0.,1-mu*mu))/(R if ground else top)
    return out,ground,cend,send

@njit(cache=True)
def weights_path(path,beta,absorb):
    n=len(path);nl=len(beta);w=np.zeros((n,nl));T=np.ones(nl)
    for k in range(n):
        for l in range(nl):
            bs=path[k,4]*beta[l];chi=bs+path[k,5]*absorb[l]
            dt=chi*path[k,6]
            v=-math.expm1(-dt)
            w[k,l]=T[l]*v*(bs/chi if chi>0 else 0.)
            T[l]*=math.exp(-dt)
    return w,T

@njit(cache=True)
def solar_columns(r,mu,R,top,H,ozscale,n=240):
    """Independent high-order Gauss-like midpoint slant-column quadrature.
    Uniform distance with squared spacing near observer; actual spherical rays.
    A ground-intersecting ray transmits zero.
    """
    if mu<0 and r*r*(1-mu*mu)<R*R-1e-4:return -1.,-1.
    d=-r*mu+math.sqrt(max(0.,top*top-r*r*(1-mu*mu)))
    cg=0.;co=0.
    for j in range(n):
        t0=d*(j/n)**2;t1=d*((j+1)/n)**2
        t=(t0+t1)/2.
        rad=math.sqrt(r*r+2*r*mu*t+t*t)
        ga,oz=densities(rad-R,H,ozscale)
        cg+=ga*(t1-t0);co+=oz*(t1-t0)
    return cg,co

@njit(parallel=True,cache=True)
def solar_table(rgrid,agrid,R,top,H,ozscale,beta,absorb,solar,sunrad,nsun=12,npath=240):
    """Finite disk sampled in elevation; uniform disk strip area weights."""
    nr=len(rgrid);na=len(agrid);nl=len(beta)
    out=np.zeros((nr,na,nl,3)) # direct normal beam, radial/tangent directions unused approximation, horizontal irradiance
    for p in prange(nr*na):
        ir=p//na;ia=p%na;r=rgrid[ir];a=agrid[ia]
        # Integrate only the visible circular segment, with its exact area.
        # This avoids stepwise disk contact from a few fixed solar samples.
        horizon=-math.acos(max(-1.,min(1.,R/r)))
        lower=max(-1.,min(1.,(horizon-a)/sunrad))
        fraction=.5-(math.asin(lower)+lower*math.sqrt(max(0.,1-lower*lower)))/math.pi
        if fraction<1e-15:continue
        norm=0.
        for k in range(nsun):
            off=lower+(1-lower)*(k+.5)/nsun
            wt=math.sqrt(max(0.,1.-off*off));norm+=wt
            aa=a+off*sunrad
            mu=math.sin(aa)
            cg,co=solar_columns(r,mu,R,top,H,ozscale,npath)
            if cg<0:continue
            for l in range(nl):
                F=solar[l]*math.exp(-beta[l]*cg-absorb[l]*co)
                out[ir,ia,l,0]+=wt*F
                out[ir,ia,l,1]+=wt*F*max(0.,mu)
        for l in range(nl):
            out[ir,ia,l,0]*=fraction/norm;out[ir,ia,l,1]*=fraction/norm
    return out

@njit(cache=True)
def interpolate(tab,rg,ag,r,alpha,l,c):
    i,u=bracket(rg,r);j,v=bracket(ag,alpha)
    return ((1-u)*((1-v)*tab[i,j,l,c]+v*tab[i,j+1,l,c])+
            u*((1-v)*tab[i+1,j,l,c]+v*tab[i+1,j+1,l,c]))


def angular_quadrature(nmu=16,nphi=16):
    # Split at horizon to improve near-horizontal integration.
    g,w=np.polynomial.legendre.leggauss(nmu//4)
    mus=[];ws=[]
    for lo,hi in [(-1.,-.12),(-.12,0.),(0.,.12),(.12,1.)]:
        mus.extend((g+1)*(hi-lo)/2+lo);ws.extend(w*(hi-lo)/2)
    phi=(np.arange(nphi)+.5)*2*PI/nphi
    return np.array(mus),np.array(ws),np.cos(phi),np.sin(phi),2*PI/nphi


def build_paths(atm,rgrid,mus,beta,absorb,nshell=56):
    H=atm.scale_height_m;R=atm.radius_m
    hh=np.unique(np.r_[H*atm.top_scale_heights*np.linspace(0,1,nshell+1)**2,
                       np.array([10000.,25000.,40000.])*atm.ozone_scale])
    shells=R+hh
    items=[];K=0
    for r in rgrid:
        row=[]
        for mu in mus:
            p,g,c,s=ray_path(r,mu,R,atm.top,shells,H,atm.ozone_scale)
            w,T=weights_path(p,beta,absorb);row.append((p,w,T,g,c,s));K=max(K,len(p))
        items.append(row)
    nr=len(rgrid);nm=len(mus);nl=len(beta)
    geom=np.zeros((nr,nm,K,4));weights=np.zeros((nr,nm,K,nl));ends=np.zeros((nr,nm,nl))
    nseg=np.zeros((nr,nm),np.int32);gb=np.zeros((nr,nm,3))
    for i in range(nr):
        for j in range(nm):
            p,w,T,g,c,s=items[i][j];n=len(p)
            geom[i,j,:n]=p[:,:4];weights[i,j,:n]=w;ends[i,j]=T;nseg[i,j]=n;gb[i,j]=[g,c,s]
    return geom,weights,ends,nseg,gb

@njit(parallel=True,cache=True)
def transport(prev,beam,rg,ag,srg,sag,mus,muw,cp,sp,dphi,geom,weights,ends,nseg,gb,albedo,include_beam):
    nr=len(rg);na=len(ag);nl=prev.shape[2];nm=len(mus);np_=len(cp)
    out=np.zeros_like(prev)
    fac=3./(16*math.pi)
    for idx in prange(nr*na):
        ir=idx//na;ia=idx%na
        ca=math.cos(ag[ia]);sa=math.sin(ag[ia])
        acc=np.zeros((nl,prev.shape[3]))
        for im in range(nm):
            mu=mus[im];st=math.sqrt(max(0.,1-mu*mu))
            for ip in range(np_):
                vx=st*cp[ip];vy=st*sp[ip]
                nu=sa*mu+ca*vx
                I=np.zeros(nl)
                if gb[ir,im,0]>.5:
                    sg=sa*gb[ir,im,1]+ca*gb[ir,im,2]*cp[ip]
                    aa=math.asin(max(-1.,min(1.,sg)))
                    for l in range(nl):
                        flux=interpolate(prev,rg,ag,rg[0],aa,l,5)
                        if include_beam:flux+=interpolate(beam,srg,sag,rg[0],aa,l,1)
                        I[l]=ends[ir,im,l]*albedo/math.pi*flux
                for k in range(nseg[ir,im]):
                    r=geom[ir,im,k,0];qr=geom[ir,im,k,1]
                    ss=sa*geom[ir,im,k,2]+ca*geom[ir,im,k,3]*cp[ip]
                    ss=max(-1.,min(1.,ss));aa=math.asin(ss)
                    qt=(nu-ss*qr)/math.sqrt(max(1e-16,1-ss*ss))
                    qp2=max(0.,1-qr*qr-qt*qt)
                    ii,u=bracket(rg,r);jj,v=bracket(ag,aa)
                    if include_beam:bi,bu=bracket(srg,r);bj,bv=bracket(sag,aa)
                    for l in range(nl):
                        M=np.empty(5)
                        for c in range(5):
                            M[c]=(1-u)*((1-v)*prev[ii,jj,l,c]+v*prev[ii,jj+1,l,c])+u*((1-v)*prev[ii+1,jj,l,c]+v*prev[ii+1,jj+1,l,c])
                        src=M[0]+qr*qr*M[1]+qt*qt*M[2]+qp2*M[3]+2*qr*qt*M[4]
                        if include_beam:
                            F=(1-bu)*((1-bv)*beam[bi,bj,l,0]+bv*beam[bi,bj+1,l,0])+bu*((1-bv)*beam[bi+1,bj,l,0]+bv*beam[bi+1,bj+1,l,0])
                            src+=F*(1+nu*nu)
                        I[l]+=weights[ir,im,k,l]*fac*max(0.,src)
                dw=muw[im]*dphi
                for l in range(nl):
                    x=I[l]*dw
                    acc[l,0]+=x;acc[l,1]+=x*mu*mu;acc[l,2]+=x*vx*vx
                    acc[l,3]+=x*vy*vy;acc[l,4]+=x*mu*vx;acc[l,5]+=x*max(mu,0.)
                    if prev.shape[3]>6:acc[l,6]+=x*max(-mu,0.)
        out[ir,ia]=acc
    return out


def grids(atm,quality):
    nr,da,nmu,nphi,nshell= {'draft':(25,1.,12,12,40),
      'standard':(35,.75,20,16,56),'fine':(47,.5,28,24,80)}[quality]
    r=atm.radius_m+atm.scale_height_m*atm.top_scale_heights*np.linspace(0,1,nr)**2
    a=np.unique(np.r_[np.arange(-90,-35,5.),np.arange(-35,30.0001,da),np.arange(35,91,5.)])*PI/180
    # Solar attenuation is interpolated on a finer grid than diffuse moments.
    sr=atm.radius_m+atm.scale_height_m*atm.top_scale_heights*np.linspace(0,1,101)**2
    sa=np.unique(np.r_[np.arange(-90,-40,2.),np.arange(-40,20.0001,.125),np.arange(20.5,90.1,.5)])*PI/180
    return r,a,sr,sa,nmu,nphi,nshell


def solve(atm,quality='standard',wave_step=20.,orders=12,wavelengths=None):
    t=time.time();r,a,sr,sa,nmu,nphi,nshell=grids(atm,quality)
    lam=np.arange(380.,800.0001,wave_step) if wavelengths is None else np.asarray(wavelengths,dtype=float)
    if np.any(np.diff(lam)<=0) or lam[0]<360 or lam[-1]>830: raise ValueError('Wavelengths must increase within 360–830 nm')
    sol=np.interp(lam,SW,SF)*atm.visible_filter
    beta=1.24062e-6*(lam/1000.)**-4*atm.density_scale
    absorb=np.interp(lam,SW,O3)*(atm.ozone_du*2.687e20/(15000*atm.ozone_scale))
    print(atm.name,'solar',len(r),len(a),len(lam),flush=True)
    beam=solar_table(sr,sa,atm.radius_m,atm.top,atm.scale_height_m,atm.ozone_scale,beta,absorb,sol,atm.sun_radius_deg*PI/180)
    mus,muw,cp,sp,dphi=angular_quadrature(nmu,nphi)
    pp=build_paths(atm,r,mus,beta,absorb,nshell)
    print(atm.name,'paths',pp[0].shape,'elapsed',time.time()-t,flush=True)
    prev=np.zeros((len(r),len(a),len(lam),6));history=[];first=None
    # Successive orders. Each transport call is linear in the previous order.
    total=np.zeros_like(prev)
    for order in range(1,orders+1):
        cur=transport(prev,beam,r,a,sr,sa,mus,muw,cp,sp,dphi,*pp,atm.ground_albedo,order==1)
        total+=cur
        if order==1:first=cur.copy()
        frac=float(np.max(cur[...,0]) / max(np.max(total[...,0]),1e-30))
        # Track ground flux relative remainder for each sampled solar angle.
        roi=(a>=np.deg2rad(-24))&(a<=np.deg2rad(90))
        rel=np.max(cur[0,roi,:,5]/np.maximum(total[0,roi,:,5],1e-12))
        history.append(dict(order=order,maximum_increment=frac,ground_max_fraction=float(rel),seconds=time.time()-t))
        print(atm.name,'order',order,'max delta',round(frac,8),'twilight relative',round(rel,5),'seconds',round(time.time()-t),flush=True)
        prev=cur
        if order>=6 and rel<.0005:break
    dest=ROOT/'work'/f'{atm.name}_{quality}_{int(wave_step)}.npz'
    temp=dest.with_name(dest.stem+'.tmp.npz')
    np.savez_compressed(temp,r=r,a=a,sr=sr,sa=sa,beam=beam,moments=total,first=first,
      last_order=prev,lam=lam,beta=beta,absorb=absorb,solar=sol,
      meta=json.dumps(dict(atmosphere=asdict(atm),quality=quality,history=history,wave_step=wave_step,elapsed_seconds=time.time()-t)))
    temp.replace(dest)
    return dest

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--world',choices=['earth','moon','moon_no_ozone','all'],default='all');p.add_argument('--quality',choices=['draft','standard','fine'],default='standard');p.add_argument('--wave-step',type=float,default=20.);p.add_argument('--orders',type=int,default=16);p.add_argument('--wavelengths',nargs='+',type=float)
    args=p.parse_args();(ROOT/'work').mkdir(exist_ok=True)
    for atm in [EARTH,MOON,MOON_ZERO]:
        if args.world in ('all',atm.name):print(solve(atm,args.quality,args.wave_step,args.orders,args.wavelengths),flush=True)
