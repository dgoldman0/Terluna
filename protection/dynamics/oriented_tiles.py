"""Finite independently oriented squares: geometric shadow, force and contact.

The receiving square is clipped against each blocker's perspective shadow
cone AND its upstream plane. This handles a blocker crossing the receiving
plane; a centre-depth ordering alone would be wrong. Polygon unions count
first interceptions once. The force uses the inherited centre-ray irradiance
approximation and exact visible-area first moments. Body eclipses are rejected.
This reference backend prioritizes correctness for bounded grouped-steering
checks; it is not used as an unchecked substitute for parallel geometry.
"""
import numpy as np
from shapely.geometry import Polygon
from shapely import union_all
from shared import constants as K
from .parallel_shadow import source_neighbours
from .natural_pattern import retarded_sun,finite_gravity
from .fleet import sun_points
from .collection import require_uneclipsed_tiles
from .cycling import ray_geometry
from .optical import unit,length

SIGNS=np.array([[-1.,-1.],[1,-1],[1,1],[-1,1]])


def vertices(q,frames,side=10000.):
    return q[:,None,:]+side/2*np.einsum('va,nba->nvb',SIGNS,frames[:,:,:2])


def clip(poly,ab,offset):
    """Clip a convex 2-D polygon to ab dot x + offset >= 0."""
    if len(poly)==0:return poly
    val=poly@ab+offset;out=[]
    for i in range(len(poly)):
        j=(i+1)%len(poly);inside=val[i]>=0;next_inside=val[j]>=0
        if inside:out.append(poly[i])
        if inside!=next_inside:
            out.append(poly[i]+(poly[j]-poly[i])*val[i]/(val[i]-val[j]))
    return np.array(out).reshape(-1,2)


def shadow_polygon(qi,fi,qj,fj,source,clear_side=9890.):
    corners=qj+SIGNS@(fj[:,:2].T*clear_side/2)
    poly=SIGNS*clear_side/2
    for k in range(4):
        a=corners[k];edge=corners[(k+1)%4]-a
        normal=np.cross(a-source,edge)
        normal/=max(length(normal),1e-100)
        if normal@(qj-a)<0:normal=-normal
        poly=clip(poly,normal@fi[:,:2],normal@(qi-a))
        if len(poly)<3:return Polygon()
    n=fj[:,2]*(-1 if fj[:,2]@(source-qj)>=0 else 1)
    poly=clip(poly,n@fi[:,:2],n@(qi-qj))
    return Polygon(poly) if len(poly)>=3 else Polygon()


def moments(q,frames,source,clear_side=9890.,every=False):
    frames=np.broadcast_to(frames,(len(q),3,3))
    neighbour,valid,_=source_neighbours(q,source,every=every)
    out=np.zeros((len(q),3));out[:,0]=clear_side**2
    for i in range(len(q)):
        shadows=[]
        for j in neighbour[i,valid[i]]:
            poly=shadow_polygon(q[i],frames[i],q[j],frames[j],source,clear_side)
            if not poly.is_empty and poly.area>0:shadows.append(poly)
        if shadows:
            blocked=union_all(shadows);area=blocked.area
            out[i,0]=max(0.,clear_side**2-area)
            if area:out[i,1:]=-area*np.array([blocked.centroid.x,blocked.centroid.y])
    return out


def photon_load(env,t,state,frames,suns=8,rotation=.317):
    q=state[:,:3];frames=np.broadcast_to(frames,(len(q),3,3));sample=env.at(t)
    require_uneclipsed_tiles(ray_geometry(q,sample),10000.)
    force=np.zeros_like(q);torque=np.zeros_like(q);first=np.zeros(len(q));bare=first.copy()
    front=first.copy();back=first.copy();normals=frames[:,:,2]
    for source in sun_points(retarded_sun(sample),suns,rotation):
        area,mx,my=moments(q,frames,source).T
        ray=source-q;cosine=np.sum(unit(ray)*normals,axis=1)
        density=K.SOLAR_CONSTANT*(K.AU/length(ray))**2/suns*abs(cosine)
        intercepted=density*area;first+=intercepted;bare+=density*9890**2
        front+=np.where(cosine>=0,intercepted,0);back+=np.where(cosine<0,intercepted,0)
        pressure=-2*env.central_fraction/K.SPEED_OF_LIGHT*density*cosine
        force+=(pressure*area)[:,None]*normals
        torque+=pressure[:,None]*np.c_[my,-mx,np.zeros(len(q))]
    return dict(reflected_force_N=force,radiation_torque_N_m=torque,
        optical=dict(first_intercept_bolometric_equivalent_W=first,
            redirected_band_optical_W=first*env.central_fraction,
            unshadowed_aperture_W=bare,mutual_shadow_equivalent_loss_W=bare-first,
            front_first_intercept_W=front,back_first_intercept_W=back))


def square_distance(q1,f1,q2,f2,side=10000.):
    """Euclidean distance: vertex/face, segment/segment and edge/face crossings."""
    h=side/2;v1=q1+SIGNS@(f1[:,:2].T*h);v2=q2+SIGNS@(f2[:,:2].T*h)
    best=np.inf
    for va,q,f in [(v1,q2,f2),(v2,q1,f1)]:
        local=(va-q)@f;d=local.copy();d[:,:2]-=np.clip(d[:,:2],-h,h)
        best=min(best,float(length(d).min()))
        for k in range(4):
            a,b=local[k],local[(k+1)%4]
            if a[2]*b[2]<0:
                hit=a+(b-a)*a[2]/(a[2]-b[2])
                if np.max(abs(hit[:2]))<=h:return 0.
    for i in range(4):
        a=v1[i];u=v1[(i+1)%4]-a
        for j in range(4):
            b=v2[j];v=v2[(j+1)%4]-b;w=a-b
            aa=u@u;bb=u@v;cc=v@v;dd=u@w;ee=v@w;det=aa*cc-bb*bb
            candidates=[(0,np.clip(ee/cc,0,1)),(1,np.clip((ee+bb)/cc,0,1)),
                        (np.clip(-dd/aa,0,1),0),(np.clip((bb-dd)/aa,0,1),1)]
            if det>1e-12*aa*cc:
                s=(bb*ee-cc*dd)/det;t=(aa*ee-bb*dd)/det
                if 0<=s<=1 and 0<=t<=1:candidates.append((s,t))
            best=min(best,min(float(length(w+s*u-t*v)) for s,t in candidates))
    return best
