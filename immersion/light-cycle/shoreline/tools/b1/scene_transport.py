"""B1 additive coupled gas/cloud/surface path tracer.

A3's validated collision and transmittance primitives are imported unchanged.
Nearest triangle/sea intersections truncate free flight. The finite solar disk
is explicitly sampled at diffuse vertices; delta paths alone see its emitted
radiance on escape. Thus direct sunlight is counted exactly once. Deep water
reflects Fresnel's fraction; transmitted energy enters an absorbing reservoir.
"""
from __future__ import annotations
import math,sys
from pathlib import Path
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/a3'))
import coupled_transport as t
import geometry as g

@njit(cache=True)
def fresnel(cos_i,eta=1.333):
    c=min(1.,max(0.,cos_i));st2=(1-c*c)/(eta*eta)
    ct=math.sqrt(max(0.,1-st2))
    rs=(c-eta*ct)/(c+eta*ct);rp=(eta*c-ct)/(eta*c+ct)
    return .5*(rs*rs+rp*rp)

@njit(cache=True)
def cosine_direction(n):
    return t.rotate(n,math.sqrt(np.random.random()),2*math.pi*np.random.random())

@njit(cache=True)
def sample_sun(sun,radius):
    if radius<=0:return sun.copy()
    return t.rotate(sun,1-np.random.random()*(1-math.cos(radius)),2*math.pi*np.random.random())

@njit(cache=True)
def reflectance(material,lam,scale):
    if material==2:r=.37+.06*(lam-550)/300 # authored pale mineral shore
    elif material==3:r=.045+.11*math.exp(-((lam-550)/55)**2)+.22/(1+math.exp(-(lam-710)/15))
    elif material==4:r=.20+.025*(lam-550)/300
    else:r=.14
    return max(0.,min(.95,r*scale))

@njit(cache=True)
def terrain_blocked(p,d,R,verts,glo,ghi,step,terrain):
    limit,planet=t.domain_end(p,d,R,600000.)
    if planet:return True
    h,n,mat=g.surface(p,d,R,verts,glo,ghi,step,terrain,limit)
    return mat!=0 and h<limit-1.e-4

@njit(cache=True)
def path_sample(origin,direction,sun,sun_radius,lam,R,top,rows,oz,beta,sigma,edges,airb,ozb,lo,hi,dims,ext,ice,cmax,verts,glo,ghi,step,gas_on,cloud_on,terrain,surface_scale,water,roulette=.95,uniform_boundary=0.,max_bounces=4096,initial_delta=True,suppress_primary=False):
    """Unit normal solar irradiance; returns L, surface-path L, failure count.

    uniform_boundary is an analytic-test isotropic sky replacing the Sun.
    The scientific image path always sets it to zero. Water=False makes sea
    use the same neutral Lambertian boundary as the far-field land.
    """
    p=origin.copy();d=direction.copy();weight=1.;L=0.;surface_L=0.;touched=False;delta=initial_delta;has_diffuse=False
    if sun_radius>0:
        sun_L=1/(math.pi*math.sin(sun_radius)**2)
        sun_integral=2/(1+math.cos(sun_radius))
    else:sun_L=0.;sun_integral=1.
    for bounce in range(max_bounces):
        start=p.copy()
        q,kind,ic,ground,trials=t.next_collision(p,d,R,top,rows,oz,beta,sigma,edges,airb,ozb,lo,hi,dims,ext,ice,cmax,1.,gas_on,cloud_on)
        if kind<0:return L,surface_L,1
        distance=math.sqrt(np.dot(q-p,q-p))
        limit=distance+1e-3
        dist,n,mat=g.surface(p,d,R,verts,glo,ghi,step,terrain,limit)
        if mat!=0 and dist<=distance+1e-4:
            p=p+d*dist
            # Physical normals; no normal-map changes to the scattering measure.
            if np.dot(n,d)>1e-8:return L,surface_L,2
            touched=True
            p=p+n*.02
            if mat==1 and water:
                weight*=fresnel(-np.dot(d,n))
                if weight==0:return L,surface_L,0
                d=d-2*np.dot(d,n)*n;d/=math.sqrt(np.dot(d,d));delta=True
            else:
                rho=reflectance(mat,lam,surface_scale)
                if rho==0:return L,surface_L,0
                if uniform_boundary==0:
                    sd=sample_sun(sun,sun_radius);cosine=max(0.,np.dot(n,sd))
                    if cosine>0 and not terrain_blocked(p,sd,R,verts,glo,ghi,step,terrain):
                        tr,failed=t.transmission(p,sd,R,top,rows,oz,beta,sigma,edges,airb,ozb,lo,hi,dims,ext,ice,cmax,gas_on,cloud_on)
                        if failed:return L,surface_L,3
                        value=weight*rho/math.pi*cosine*tr*sun_integral
                        L+=value;surface_L+=value
                weight*=rho;d=cosine_direction(n);delta=False;has_diffuse=True
        else:
            p=q
            if kind==0:
                value=weight*uniform_boundary
                if delta and (not suppress_primary or has_diffuse) and sun_radius>0 and uniform_boundary==0 and np.dot(d,sun)>=math.cos(sun_radius):value+=weight*sun_L
                L+=value
                if touched:surface_L+=value
                return L,surface_L,0
            if kind==3:return L,surface_L,0
            if uniform_boundary==0:
                sd=sample_sun(sun,sun_radius)
                if not terrain_blocked(p,sd,R,verts,glo,ghi,step,terrain):
                    tr,failed=t.transmission(p,sd,R,top,rows,oz,beta,sigma,edges,airb,ozb,lo,hi,dims,ext,ice,cmax,gas_on,cloud_on)
                    if failed:return L,surface_L,3
                    ph=t.rayleigh_phase(np.dot(d,sd)) if kind==1 else t.cloud.phase(np.dot(d,sd),ic)
                    value=weight*ph*tr*sun_integral;L+=value
                    if touched:surface_L+=value
            d=t.rayleigh_scatter(d) if kind==1 else t.cloud.scatter(d,ic);delta=False;has_diffuse=True
        if bounce>=11:
            # The added surface throughput already accounts for absorbed energy.
            survive=min(roulette,max(.05,weight))
            if np.random.random()>=survive:return L,surface_L,0
            weight/=survive
    return L,surface_L,4

@njit(cache=True,nogil=True)
def pixel_batch(origin,w,h,jstart,jend,fov,yaw,pitch,projection,sun,sun_radius,lam,beta,sigma,xyz_weights,R,top,rows,oz,edges,airb,ozb,lo,hi,dims,ext,ice,cmax,verts,glo,ghi,step,cloud_on,terrain,surface_scale,water,spp,seed):
    B=len(lam);out=np.zeros((jend-jstart,w,B));errors=np.zeros_like(out);xyz_se=np.zeros((jend-jstart,w,3));surface_xyz=np.zeros_like(xyz_se);fails=0
    for j in range(jstart,jend):
        for i in range(w):
            sums=np.zeros(B);sums2=np.zeros(B);xyz2=np.zeros(3);surf=np.zeros(3)
            for sample in range(spp):
                sseed=(seed+104729*(j*w+i)+196613*sample)%2147483647
                np.random.seed(sseed)
                u=(i+np.random.random())/w;v=(j+np.random.random())/h
                if projection==0:d=g.ray_direction(u,v,w/h,fov,yaw,pitch)
                else:
                    az=2*math.pi*u;el=(1-v)*math.pi/2
                    d=np.array([math.sin(az)*math.cos(el),math.sin(el),math.cos(az)*math.cos(el)])
                X=np.zeros(3)
                for b in range(B):
                    # Common random numbers across wavelengths preserve their
                    # covariance. Per-sample XYZ records give correct color SE.
                    np.random.seed(sseed+17)
                    value,sv,fail=path_sample(origin,d,sun,sun_radius,lam[b],R,top,rows,oz,beta[b],sigma[b],edges,airb,ozb,lo,hi,dims,ext,ice,cmax,verts,glo,ghi,step,True,cloud_on,terrain,surface_scale,water,.95,0.,4096,True,True)
                    sums[b]+=value;sums2[b]+=value*value;fails+=int(fail!=0)
                    X+=value*xyz_weights[b];surf+=sv*xyz_weights[b]
                xyz2+=X*X
            out[j-jstart,i]=sums/spp
            errors[j-jstart,i]=np.sqrt(np.maximum(0.,(sums2-sums*sums/spp)/(max(1,spp-1)*spp)))
            Xmean=(sums/spp)@xyz_weights
            xyz_se[j-jstart,i]=np.sqrt(np.maximum(0.,(xyz2-Xmean*Xmean*spp)/(max(1,spp-1)*spp)))
            surface_xyz[j-jstart,i]=surf/spp
    return out,errors,xyz_se,surface_xyz,fails
