"""B1 triangle-heightfield coast over a spherical sea, metres in A3 coordinates.

The coast is authored test geometry, not lunar geography. Every transport and
shadow ray uses the same triangles. A 2-D DDA visits candidate cells; independent
brute-force tests validate it. Water is the smooth analytic sea-level sphere.
"""
from __future__ import annotations
import math
import numpy as np
from numba import njit

@njit(cache=True)
def sphere_hit(p,d,R):
    q=np.array([p[0],p[1]+R,p[2]])
    pd=np.dot(q,d);c=np.dot(q,q)-R*R;disc=pd*pd-c
    if disc<0:return 1.e30
    rt=math.sqrt(disc);a=-pd-rt;b=-pd+rt
    if a>1.e-4:return a
    if b>1.e-4 and c<0:return b
    return 1.e30

@njit(cache=True)
def triangle(p,d,a,b,c,limit):
    e=b-a;f=c-a;v=np.cross(d,f);det=np.dot(e,v)
    if abs(det)<1.e-12:return limit,np.zeros(3)
    inv=1./det;q=p-a;u=np.dot(q,v)*inv
    if u< -1.e-9 or u>1.+1.e-9:return limit,np.zeros(3)
    z=np.cross(q,e);w=np.dot(d,z)*inv
    if w< -1.e-9 or u+w>1.+1.e-9:return limit,np.zeros(3)
    t=np.dot(f,z)*inv
    if t<1.e-4 or t>=limit:return limit,np.zeros(3)
    n=np.cross(e,f)
    if n[1]<0:n=-n
    n/=math.sqrt(np.dot(n,n))
    return t,n

@njit(cache=True)
def aabb(p,d,lo,hi,limit):
    a=0.;b=limit
    for k in range(3):
        if abs(d[k])<1e-14:
            if p[k]<lo[k] or p[k]>hi[k]:return 1.,0.
        else:
            x=(lo[k]-p[k])/d[k];y=(hi[k]-p[k])/d[k]
            a=max(a,min(x,y));b=min(b,max(x,y))
    return a,b

@njit(cache=True)
def material_at(p,n,R):
    h=math.sqrt(p[0]**2+(p[1]+R)**2+p[2]**2)-R
    if h<28.:return 2 # pale shore
    if h>920. or n[1]<.69:return 4 # rock
    return 3 # vegetation-coloured diffuse terrain, no trees

@njit(cache=True)
def surface(p,d,R,verts,lo,hi,step,enabled=True,limit=1e30):
    hit=min(limit,sphere_hit(p,d,R));n=np.zeros(3);mat=0
    if hit<limit:
        q=p+d*hit;n=np.array([q[0],q[1]+R,q[2]])/R
        # The authored water basin extends through the local patch; far terrain
        # is a specified neutral diffuse spherical boundary.
        mat=1 if abs(q[0])<24000. and abs(q[2])<24000. else 5
    if not enabled:return hit,n,mat
    enter,leave=aabb(p,d,lo,hi,hit)
    if leave<=enter:return hit,n,mat
    N=verts.shape[0]-1;q=p+d*(enter+1.e-5)
    ix=min(N-1,max(0,int(math.floor((q[0]-lo[0])/step))))
    iz=min(N-1,max(0,int(math.floor((q[2]-lo[2])/step))))
    sx=1 if d[0]>0 else -1;sz=1 if d[2]>0 else -1
    tx=(lo[0]+(ix+(1 if sx>0 else 0))*step-p[0])/d[0] if abs(d[0])>1e-14 else 1e30
    tz=(lo[2]+(iz+(1 if sz>0 else 0))*step-p[2])/d[2] if abs(d[2])>1e-14 else 1e30
    dx=step/abs(d[0]) if abs(d[0])>1e-14 else 1e30
    dz=step/abs(d[2]) if abs(d[2])>1e-14 else 1e30
    for _ in range(2*N+4):
        a=verts[iz,ix];b=verts[iz,ix+1];c=verts[iz+1,ix];e=verts[iz+1,ix+1]
        nextt=min(tx,tz)
        ya=p[1]+d[1]*enter;yb=p[1]+d[1]*min(nextt,leave)
        ymin=min(a[1],b[1],c[1],e[1]);ymax=max(a[1],b[1],c[1],e[1])
        if max(ya,yb)>=ymin-1e-5 and min(ya,yb)<=ymax+1e-5:
            t,nt=triangle(p,d,a,b,c,hit)
            if t<hit:hit=t;n=nt;mat=material_at(p+d*t,n,R)
            t,nt=triangle(p,d,b,e,c,hit)
            if t<hit:hit=t;n=nt;mat=material_at(p+d*t,n,R)
        if hit<=nextt+1e-5 or nextt>leave:break
        enter=nextt
        if tx<tz:ix+=sx;tx+=dx
        else:iz+=sz;tz+=dz
        if ix<0 or ix>=N or iz<0 or iz>=N:break
    return hit,n,mat

@njit(cache=True)
def sea_y(x,z,R):return -(x*x+z*z)/(math.sqrt(max(0.,R*R-x*x-z*z))+R)

def make_coast(R=1737400.,N=97,extent=24000.):
    x=np.linspace(-extent,extent,N);X,Z=np.meshgrid(x,x)
    # Ridges framing a branching inlet. Multiscale relief has a fixed analytic
    # seed, allowing exact mesh hashes without an external texture dependency.
    hill=lambda x,z,sx,sz,amp:amp*np.exp(-((X-x)/sx)**2-((Z-z)/sz)**2)
    H=-170.+hill(-7300,1500,5700,12000,1750)+hill(8700,6500,6000,8600,1400)+hill(2600,17300,7000,5700,1450)
    H+=hill(-2700,-6500,3300,4000,340)+hill(2000,5400,1000,1800,240)
    texture=.52*np.sin(X/2300+Z/3300)+.27*np.sin(X/720-Z/1350)+.13*np.cos(X/310+Z/480)+.08*np.sin(X/170-Z/240)
    H-=650*np.exp(-((X-900*np.sin(Z/6300))/3000)**2)*np.exp(-(Z/23000)**4)
    H+=np.maximum(H,0)*.27*texture
    edge=np.clip((extent-np.maximum(abs(X),abs(Z)))/4000.,0.,1.)
    edge=edge*edge*(3-2*edge)
    H=H*edge-170*(1-edge)
    Y=-(X*X+Z*Z)/(np.sqrt(R*R-X*X-Z*Z)+R)+H
    verts=np.stack([X,Y,Z],axis=-1).astype(np.float64)
    lo=np.array([-extent,float(Y.min())-1.,-extent]);hi=np.array([extent,float(Y.max())+1.,extent])
    return verts,lo,hi,float(x[1]-x[0])

@njit(cache=True)
def ray_direction(u,v,aspect,fov_deg,yaw_deg,pitch_deg):
    yaw=math.radians(yaw_deg);pitch=math.radians(pitch_deg)
    f=np.array([math.sin(yaw)*math.cos(pitch),math.sin(pitch),math.cos(yaw)*math.cos(pitch)])
    r=np.array([math.cos(yaw),0.,-math.sin(yaw)]);up=np.cross(f,r)
    scale=math.tan(math.radians(fov_deg)*.5)
    d=f+(2*u-1)*scale*r+(1-2*v)*scale/aspect*up
    return d/math.sqrt(np.dot(d,d))

@njit(cache=True)
def camera_maps(origin,w,h,R,verts,lo,hi,step,fov=100.,yaw=0.,pitch=16.):
    out=np.zeros((h,w,5))
    for j in range(h):
        for i in range(w):
            d=ray_direction((i+.5)/w,(j+.5)/h,w/h,fov,yaw,pitch)
            t,n,mat=surface(origin,d,R,verts,lo,hi,step)
            out[j,i,0]=t;out[j,i,1]=mat;out[j,i,2:]=n
    return out
