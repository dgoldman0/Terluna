"""Point-source driver for the established spherical Rayleigh operator.

The original finite-Sun producer remains a reproducible historical baseline.
This driver reuses its paths, spatial grids, quadrature and scattering operator;
it replaces the source with a pencil beam and checks orders over the full night.
"""
import json
import math
import time
from pathlib import Path

import numpy as np
from numba import njit, prange
from illumination.sky.solved_transport import layer_lengths, grids, quadrature, build_paths
from illumination.sky.solver import transport


@njit(parallel=True,cache=True)
def point_table(rgrid,agrid,edges,extinction):
    nr,na,nl=len(rgrid),len(agrid),extinction.shape[1]
    out=np.zeros((nr,na,nl,3))
    for p in prange(nr*na):
        ir,ia=p//na,p%na
        mu=math.sin(agrid[ia])
        lengths,blocked=layer_lengths(rgrid[ir],mu,edges)
        if blocked:continue
        for l in range(nl):
            tau=0.
            for j in range(len(lengths)):tau+=lengths[j]*extinction[j,l]
            out[ir,ia,l,0]=math.exp(-tau)
            out[ir,ia,l,1]=math.exp(-tau)*max(0.,mu)
    return out


def solve(optical,radius,quality='standard',albedo=.1,max_orders=240,tolerance=1e-6,path=None):
    start=time.monotonic()
    edges=radius+optical['height']
    sca=optical['scattering']/np.diff(edges)[:,None]
    absorb=optical['absorption']/np.diff(edges)[:,None]
    r,a,sr,sa,angular,nphi,shells=grids(edges,quality)
    beam=point_table(sr,sa,edges,sca+absorb)
    mus,muw,cp,sp,dphi=quadrature(r,radius,angular,nphi)
    paths=build_paths(r,mus,edges,sca,absorb,shells)
    print('point transport grids',len(r),len(a),len(optical['energy']),flush=True)
    prev=np.zeros((len(r),len(a),len(optical['energy']),7))
    total=np.zeros_like(prev);last=np.zeros_like(prev)
    active=np.arange(len(optical['energy']));orders=np.zeros(len(active),dtype=int)
    history=[];active_paths,active_beam=paths,beam
    for order in range(1,max_orders+1):
        current=transport(np.ascontiguousarray(prev[:,:,active,:]),active_beam,r,a,sr,sa,
                          mus,muw,cp,sp,dphi,*active_paths,albedo,order==1)
        total[:,:,active,:]+=current
        last[:,:,active,:]=current;orders[active]=order
        # Every retained source angle, including the antisolar point.
        denominator=np.maximum(total[0][:,active,5],1e-12)
        per_channel=np.max(current[0,:,:,5]/denominator,axis=0)
        per_maximum=np.max(current[:,:,:,0],axis=(0,1))/np.maximum(np.max(total[:,:,active,0],axis=(0,1)),1e-30)
        history.append(dict(order=order,ground_increment_fraction=float(per_channel.max()),
                            maximum_increment_fraction=float(per_maximum.max()),active_channels=len(active),
                            elapsed_seconds=time.monotonic()-start))
        if order==1 or order%5==0:
            print('order',order,'ground increment',float(per_channel.max()),'seconds',round(time.monotonic()-start,1),flush=True)
        prev[:,:,active,:]=current
        if order>=8:
            keep=(per_channel>=tolerance)|(per_maximum>=tolerance)
            if not np.any(keep):
                active=active[:0];break
            if np.any(~keep):
                active=active[keep]
                active_paths=(paths[0],np.ascontiguousarray(paths[1][:,:,:,active]),
                              np.ascontiguousarray(paths[2][:,:,active]),paths[3],paths[4])
                active_beam=np.ascontiguousarray(beam[:,:,active,:])
    meta=dict(quality=quality,albedo=albedo,history=history,increment_tolerance=tolerance,
              order_cap=max_orders,order_converged=len(active)==0,channel_orders=orders.tolist(),
              solar_disk_radius_deg=0.,convergence_min_elevation_deg=-90.,
              angular_order_per_interval=angular,azimuth_samples=nphi,seconds=time.monotonic()-start)
    result=dict(r=r,a=a,sr=sr,sa=sa,beam=beam,moments=total,last_order=last,
                edges=edges,scattering=sca,absorption=absorb,shells=shells,
                **{k:optical[k] for k in ('energy','xyz','wavelength','band')})
    if path:
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
        temporary=path.with_suffix('.tmp.npz')
        np.savez_compressed(temporary,**result,meta=json.dumps(meta));temporary.replace(path)
    return result,meta
