"""Finite-source disk quadrature with a phase-dependent illuminated hemisphere.

The Lambert incidence pattern controls spatial weighting only. The empirical
Earth phase curve already sets the total spectral irradiance; normalized weights
avoid applying a second phase law. Surface/cloud structure is unresolved.
"""
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=12)
def gauss(order):
    return np.polynomial.legendre.leggauss(order)


def samples(elevation_deg, radius_deg, phase_deg=None, bright_limb_rad=0., order=12):
    node, weight=gauss(order)
    y=node[:,None]
    half=np.sqrt(1-y*y)
    alpha=0. if phase_deg is None else np.radians(phase_deg)
    lower=-half if phase_deg is None else -np.cos(alpha)*half
    width=half-lower
    x=lower+(node[None,:]+1)*width/2
    z=np.sqrt(np.maximum(0.,1-x*x-y*y))
    w=weight[:,None]*weight[None,:]*width/2
    if phase_deg is not None:
        w=w*np.maximum(0.,x*np.sin(alpha)+z*np.cos(alpha))
    sr=np.sin(np.radians(radius_deg))
    direction_cos=np.sqrt(1-sr*sr*(x*x+y*y))
    w=w/direction_cos
    vertical=x*np.cos(bright_limb_rad)+y*np.sin(bright_limb_rad)
    e=np.radians(elevation_deg)
    rays=np.degrees(np.arcsin(np.clip(np.sin(e)*direction_cos+np.cos(e)*sr*vertical,-1,1)))
    norm=w.sum()
    return rays.ravel(), (w/norm).ravel() if norm>0 else np.zeros(w.size)


def integrate(angles, values, elevation_deg, radius_deg, phase_deg=None, bright_limb_rad=0., order=12):
    rays,weights=samples(elevation_deg,radius_deg,phase_deg,bright_limb_rad,order)
    return float(weights@np.interp(rays,angles,values))
