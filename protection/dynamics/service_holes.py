"""Adaptive missing-ray constraints from full receiver polygon differences."""
import numpy as np
import shapely
from shapely.geometry import Point
from .service_constraints import projection
from .packing_control import matrices
from .natural_pattern import retarded_sun
from .cycling import sail_basis
from .fleet import projected_polygons
from .optical import unit
from shared import constants as K


def uncovered_rows(env,base,current,maps,arrival,arcs,margin=200.,per_source=3):
    constraints=[];worst=1.;records=[]
    signs=np.array([[-1,-1],[1,-1],[1,1],[-1,1]])*9890/2
    for t in np.linspace(arrival,arrival+21600,9):
        y=current([t])[0];original=base([t])[0];sun=retarded_sun(env.at(t));u=unit(sun);_,b,c=sail_basis(sun)
        target=projection(original[:,:3].mean(axis=0)[None,:],sun,u,b,c)[0][0]
        disk=Point(*target).buffer(10000,quad_segs=128)
        vertices=y[:,None,:3]+b*signs[None,:,:1]+c*signs[None,:,1:]
        theta=np.arange(12)*2*np.pi/12+.137
        sources=np.r_[sun[None,:],sun+K.SUN_RADIUS*(np.cos(theta)[:,None]*b+np.sin(theta)[:,None]*c)]
        B=matrices(maps,[t],arcs)[0,:3,6:]
        for s,source in enumerate(sources):
            union=shapely.union_all(projected_polygons(vertices,source,u,b,c));missing=disk.difference(union)
            fraction=1-missing.area/disk.area;worst=min(worst,fraction)
            if missing.area<1e-4:continue
            pieces=list(missing.geoms) if hasattr(missing,'geoms') else [missing]
            xy,jac,mag=projection(original[:,:3],source,u,b,c);actual=projection(y[:,:3],source,u,b,c)[0]
            for piece in sorted(pieces,key=lambda x:x.area,reverse=True)[:per_source]:
                point=piece.representative_point();target=np.array([point.x,point.y])
                i=int(np.argmin(np.max(abs(actual-target),axis=1)))
                for axis in range(2):
                    for sign in [-1,1]:constraints.append((i,sign*(jac[i,axis]@B),float(9890/2*mag[i]-margin+sign*(target[axis]-xy[i,axis]))))
            records.append(dict(time_s=float(t),source=s,covered_fraction=float(fraction),missed_area_m2=float(missing.area)))
    return constraints,dict(minimum_source_coverage=float(worst),uncovered_source_dates=len(records),worst=sorted(records,key=lambda x:x['covered_fraction'])[:3])
