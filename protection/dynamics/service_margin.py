"""Coverage margin for a specified finite-Sun ray and physical service frame."""
import numpy as np
import shapely
from shapely.geometry import Point

from shared import constants as K
from .natural_pattern import retarded_sun
from .cycling import sail_basis
from .fleet import projected_polygons
from .optical import unit


def service_margin(state,frame,sample,centre,solar_xy,radius=10000.):
    """Signed geometric margin around the whole receiver disk for one source.

    Positive margin is reduced by the circle polygon's sagitta. Negative
    values use missed-area equivalent radius as a search objective.
    Motion during light flight matches the existing local-coverage model.
    """
    sun=retarded_sun(sample);u=unit(sun);_,b,c=sail_basis(sun)
    source=sun+K.SUN_RADIUS*(solar_xy[0]*b+solar_xy[1]*c)
    q=state[:,:3]-(state[:,3:]+sample['moon_v'])*np.maximum(state[:,:3]@u,0.)[:,None]/K.SPEED_OF_LIGHT
    qc=centre[:3]-(centre[3:]+sample['moon_v'])*max(centre[:3]@u,0.)/K.SPEED_OF_LIGHT
    intercept=sun+(qc-sun)*(sun@u)/(sun@u-qc@u)
    disk=Point(intercept@b,intercept@c).buffer(radius,quad_segs=128)
    signs=np.array([[-1,-1],[1,-1],[1,1],[-1,1]])*9890/2
    vertices=q[:,None,:]+signs[None,:,:1]*frame[:,0]+signs[None,:,1:]*frame[:,1]
    union=shapely.union_all(projected_polygons(vertices,source,u,b,c))
    missing=disk.difference(union).area
    sagitta=radius*(1-np.cos(np.pi/512))
    margin=(-np.sqrt(missing/np.pi) if missing>1e-6 else disk.distance(union.boundary)-sagitta)
    return dict(margin_m=float(margin),missed_area_m2=float(missing),
        covered_fraction=float(1-missing/disk.area),circle_sagitta_allowance_m=float(sagitta))
