"""Deterministic camera/Sun and camera/water/Sun paths.

These terms are excluded from Monte Carlo camera estimates. Independent angular
subpixel quadrature avoids rare, enormous Bernoulli Sun hits and keeps the finite
solar image outside the optional spatial filter. Other specular light paths
(e.g. gas -> water -> Sun) remain in the full transport estimator.
"""
import math
import numpy as np
from numba import njit
import geometry as g
import scene_transport as b
import deterministic as dt

@njit(cache=True)
def direction(u,v,w,h,fov,yaw,pitch,projection):
    if projection==0:return g.ray_direction(u,v,w/h,fov,yaw,pitch)
    az=u*2*math.pi;el=(1-v)*math.pi/2
    return np.array([math.sin(az)*math.cos(el),math.sin(el),math.cos(az)*math.cos(el)])

@njit(cache=True)
def leg(origin,d,R,top,V,glo,ghi,step,water):
    limit,ground=b.t.domain_end(origin,d,R,top);dist,n,mat=g.surface(origin,d,R,V,glo,ghi,step,True,limit+1e-3)
    if mat:
        if mat!=1 or not water:return origin,d,0.,0.,False
        p=origin+d*dist+n*.02;dd=d-2*np.dot(d,n)*n
        return p,dd,b.fresnel(-np.dot(d,n)),dist,True
    return origin,d,1.,0.,True

@njit(cache=True,nogil=True)
def primary_frame(origin,w,h,fov,yaw,pitch,projection,sun,radius,beta,sigma,R,top,rows,oz,shells,lo,hi,dims,ext,ice,V,glo,ghi,step,water,grid=8):
    out=np.zeros((h,w,len(beta)));surf=np.zeros_like(out)
    nodes=np.array([-.9602898564975363,-.7966664774136267,-.525532409916329,-.1834346424956498,.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363]);weights=np.array([.1012285362903763,.2223810344533745,.3137066458778873,.362683783378362,.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763])
    margin=3*max(math.radians(fov)/w,2*math.pi/w if projection else 0.)
    source=1/(math.pi*math.sin(radius)**2)
    for j in range(h):
        for i in range(w):
            dc=direction((i+.5)/w,(j+.5)/h,w,h,fov,yaw,pitch,projection)
            p,dd,F,dist,ok=leg(origin,dc,R,top,V,glo,ghi,step,water)
            if not ok or np.dot(dd,sun)<math.cos(radius+margin):continue
            for sy in range(grid):
                for sx in range(grid):
                    d=direction((i+(sx+.5)/grid)/w,(j+(sy+.5)/grid)/h,w,h,fov,yaw,pitch,projection)
                    p,dd,F,dist,ok=leg(origin,d,R,top,V,glo,ghi,step,water)
                    if not ok or np.dot(dd,sun)<math.cos(radius) or b.terrain_blocked(p,dd,R,V,glo,ghi,step,True):continue
                    air,ozone,block=dt.gas_columns(p,dd,R,top,rows,oz,shells,nodes,weights)
                    if block:continue
                    tau=dt.cloud_depth(p,dd,lo,hi,dims,ext,ice)
                    if dist>0:
                        a,o,_=dt.gas_columns(origin,d,R,top,rows,oz,shells,nodes,weights,max(0.,dist-.03));air+=a;ozone+=o
                        tau+=dt.cloud_depth(origin,d,lo,hi,dims,ext,ice,dist)
                    for k in range(len(beta)):
                        val=F*source*math.exp(-beta[k]*air-sigma[k]*b.t.DU*ozone-tau)/(grid*grid)
                        out[j,i,k]+=val
                        if dist>0:surf[j,i,k]+=val
    return out,surf
