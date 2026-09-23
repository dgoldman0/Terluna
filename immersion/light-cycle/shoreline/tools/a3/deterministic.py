"""Independent deterministic checks and matched-input production-closure replay.

Cloud line integrals are exact for the frozen trilinear grid: split at cell
planes, then use two-point Gauss on the resulting cubic polynomial. Molecular
columns use A1's independent shell-split Gauss quadrature, not MC majorants.
"""
from __future__ import annotations
import math
import numpy as np
from numba import njit
import coupled_transport as t

@njit(cache=True)
def cloud_intervals(p,d,lo,hi,dims,length=1e99):
    a,b=t.cloud.box(p,d,lo,hi);b=min(b,length)
    ts=np.empty(int(np.sum(dims))+2)
    if b<=a:return ts[:0]
    ts[0]=a;ts[1]=b;n=2
    for k in range(3):
        if abs(d[k])>1e-14:
            for i in range(1,dims[k]-1):
                x=lo[k]+(hi[k]-lo[k])*i/(dims[k]-1)
                s=(x-p[k])/d[k]
                if a+1e-8<s<b-1e-8:ts[n]=s;n+=1
    return np.sort(ts[:n])

@njit(cache=True)
def cloud_depth(p,d,lo,hi,dims,ext,ice,length=1e99):
    ts=cloud_intervals(p,d,lo,hi,dims,length);tau=0.;g=1/math.sqrt(3.)
    for i in range(1,len(ts)):
        mid=(ts[i]+ts[i-1])/2;half=(ts[i]-ts[i-1])/2
        sa,_=t.cloud.voxel(p+d*(mid-half*g),lo,hi,dims,ext,ice)
        sb,_=t.cloud.voxel(p+d*(mid+half*g),lo,hi,dims,ext,ice)
        tau+=half*(sa+sb)
    return tau

@njit(cache=True)
def gas_columns(p,d,cR,ctop,rows,oz,shells,nodes,weights,length=1e99):
    r,h=t.radial(p,cR);q=np.array([p[0],p[1]+cR,p[2]]);mu=np.dot(q,d)/r
    return t.optics.columns(cR,ctop,h,mu,length,rows,oz,shells,nodes,weights)

@njit(cache=True)
def view_intervals(p,d,R,top,lo,hi,dims,extra_edges):
    ts,g=t.segments(p,d,R,top,extra_edges,lo,hi,True)
    cs=cloud_intervals(p,d,lo,hi,dims,ts[-1]);a=np.empty(len(ts)+len(cs));a[:len(ts)]=ts;a[len(ts):]=cs
    return np.sort(a),g

@njit(cache=True,nogil=True)
def single_ray(p,d,sun,R,top,rows,oz,shells,edges,lo,hi,dims,ext,ice,calbedo,beta,sigma,solar,view_nodes,view_weights,col_nodes,col_weights,cloud_on=True,gas_on=True):
    ts,ground=view_intervals(p,d,R,top,lo,hi,dims,edges);L=np.zeros(len(beta));rphase=t.rayleigh_phase(np.dot(d,sun))
    if ground:return L
    for i in range(1,len(ts)):
        mid=(ts[i]+ts[i-1])/2;half=(ts[i]-ts[i-1])/2
        if half<1e-7:continue
        for j in range(len(view_nodes)):
            s=mid+half*view_nodes[j];q=p+d*s;r,h=t.radial(q,R)
            air=t.molecular.sample_air(h,rows) if gas_on else 0.
            vga=vgo=sga=sgo=0.;blocked=False
            if gas_on:
                vga,vgo,_=gas_columns(p,d,R,top,rows,oz,shells,col_nodes,col_weights,s)
                sga,sgo,blocked=gas_columns(q,sun,R,top,rows,oz,shells,col_nodes,col_weights)
            else:
                _,blocked=t.domain_end(q,sun,R,top)
            if blocked:continue
            c=ic=ctau=0.
            if cloud_on:
                c,ic=t.cloud.voxel(q,lo,hi,dims,ext,ice)
                ctau=cloud_depth(p,d,lo,hi,dims,ext,ice,s)+cloud_depth(q,sun,lo,hi,dims,ext,ice)
            ph=t.cloud.phase(np.dot(d,sun),ic)
            for k in range(len(beta)):
                tr=math.exp(-beta[k]*(vga+sga)-sigma[k]*t.DU*(vgo+sgo)-ctau)
                L[k]+=half*view_weights[j]*solar[k]*(beta[k]*air*rphase+c*calbedo*ph)*tr
    return L

def independent_single(c,p,d,sun,view_order=6,gas_order=12,cloud_on=True,gas_on=True):
    vn,vw=np.polynomial.legendre.leggauss(view_order);cn,cw=np.polynomial.legendre.leggauss(gas_order)
    edges=np.unique(np.r_[np.linspace(0,1,65)**2*c['top'],c['shells']])
    return single_ray(np.array(p,float),np.array(d,float),np.array(sun,float),c['R'],c['top'],c['rows'],c['oz'],c['shells'],edges,c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['calbedo'],c['beta'] if gas_on else c['beta']*0,c['sigma'] if gas_on else c['sigma']*0,c['solar'],vn,vw,cn,cw,cloud_on,gas_on)

def a2_clear(c,p,d,sun):
    p,d,sun=map(lambda a:np.array(a,float),(p,d,sun));q=p+np.array([0.,c['R'],0.]);n=q/np.linalg.norm(q)
    ve=math.asin(np.clip(np.dot(n,d),-1,1));se=math.asin(np.clip(np.dot(n,sun),-1,1))
    ca=(np.dot(d,sun)-math.sin(ve)*math.sin(se))/max(1e-15,math.cos(ve)*math.cos(se))
    az=math.acos(np.clip(ca,-1,1))
    r,h=t.radial(p,c['R'])
    if abs(h-2)>1e-4:raise ValueError('A2 eval_ray reference expects radial altitude 2 m')
    shells=t.molecular.shell_altitudes(c['top'],128,c['oz'],c['profile']['sounding']['end_m'])
    return t.molecular.eval_ray(ve,az,se,c['R'],c['top'],c['rows'],c['oz'],c['beta'],c['sigma'],c['beam'],c['rg'],c['ag'],c['srg'],c['sag'],c['mom'],shells)

@njit(cache=True,nogil=True)
def closure_replay(p,d,sun,R,top,rows,oz,shells,lo,hi,dims,ext,ice,calbedo,beta,sigma,solar,clear,diffuse,H,vn,vw,cn,cw):
    """Revision-04 illumination/composite expression, with shared A2 inputs.

    Kept spectral here to separate its formula from legacy RGB/gas mismatch.
    Ground-reflection term is zero to match A2's black boundary. Molecular eye
    attenuation and solar beam use the shared profile. Explicit empirical terms
    (.24, .34, .05, .09, .16) retain their production values. This is an equation
    replay, not measurement of the unchanged production GPU/cache output.
    """
    ts=cloud_intervals(p,d,lo,hi,dims);out=np.zeros(len(beta))
    if len(ts)==0:return clear.copy()
    total_tau=cloud_depth(p,d,lo,hi,dims,ext,ice);out=clear*math.exp(-total_tau)
    for i in range(1,len(ts)):
        mid=(ts[i]+ts[i-1])/2;half=(ts[i]-ts[i-1])/2
        if half<1e-8:continue
        for j in range(len(vn)):
            s=mid+half*vn[j];q=p+d*s;co,ic=t.cloud.voxel(q,lo,hi,dims,ext,ice)
            if co==0:continue
            tau=cloud_depth(p,d,lo,hi,dims,ext,ice,s);sun_tau=cloud_depth(q,sun,lo,hi,dims,ext,ice)
            vga,vgo,_=gas_columns(p,d,R,top,rows,oz,shells,cn,cw,s)
            sga,sgo,blocked=gas_columns(q,sun,R,top,rows,oz,shells,cn,cw)
            _,h=t.radial(q,R);ph=t.cloud.phase(np.dot(d,sun),ic)
            for k in range(len(beta)):
                eye=math.exp(-beta[k]*vga-sigma[k]*t.DU*vgo)
                direct=0. if blocked else solar[k]*math.exp(-beta[k]*sga-sigma[k]*t.DU*sgo)
                sky=diffuse[k]/math.pi*math.exp(-h/(2*H))*(.24+.34*math.exp(-sun_tau*.05))
                source=sky+direct*(ph*math.exp(-sun_tau)+.09*math.exp(-sun_tau*.16))
                out[k]+=half*vw[j]*co*math.exp(-tau)*(source*eye+clear[k]*(1-eye))
    return out

def production_replay(c,p,d,sun,clear,order=6):
    p,d,sun=map(lambda a:np.array(a,float),(p,d,sun));q=p+np.array([0.,c['R'],0.]);n=q/np.linalg.norm(q);a=math.asin(np.clip(np.dot(n,sun),-1,1))
    diffuse=np.array([t.molecular.interpolate(c['mom'],c['rg'],c['ag'],c['R'],a,k,5) for k in range(len(c['lam']))])
    H=float(c['profile']['planet']['g0']) if 'g0' in c['profile']['planet'] else None
    # Use the sounding's physically derived local scale height as the empirical
    # diffuse falloff length, rather than reverting to its old gas-density proxy.
    row=c['profile']['rows'][0];H=row['p']/(row['rho']*row['g'])
    vn,vw=np.polynomial.legendre.leggauss(order);cn,cw=np.polynomial.legendre.leggauss(12)
    L=closure_replay(p,d,sun,c['R'],c['top'],c['rows'],c['oz'],c['shells'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['calbedo'],c['beta'],c['sigma'],c['solar'],clear,diffuse,H,vn,vw,cn,cw)
    return L,diffuse,H
