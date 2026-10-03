"""Spherical successive-orders transport through tabulated molecular layers.

Reuses the scalar Rayleigh moment operator and shell geometry from solver.py.
Solar attenuation and formal-ray weights integrate each homogeneous layer
exactly. Diffuse sources retain spatial, angular and scattering-order errors;
the exporter records their numerical checks alongside the resulting sky.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
from numba import njit, prange

from illumination.sky.solver import bracket, interpolate, ray_path, transport
from shared.constants import AU, SUN_RADIUS

SUN_RADIUS_RAD = math.asin(SUN_RADIUS / AU)


@njit(cache=True)
def layer_lengths(r, mu, edges):
    """Length in each nested spherical shell along an outward half-line."""
    b2 = r*r*(1-mu*mu)
    centre = -r*mu
    previous = 0.
    out = np.zeros(len(edges)-1)
    for j in range(len(edges)):
        d = edges[j]*edges[j]-b2
        length = 0.
        if d > 0:
            root = math.sqrt(d)
            length = max(0., centre+root) - max(0., centre-root)
        if j == 0:
            if mu < 0 and length > 1e-7:
                return out, True
            length = 0.
        else:
            out[j-1] = max(0., length-previous)
        previous = length
    return out, False


@njit(parallel=True, cache=True)
def solar_table(rgrid, agrid, edges, extinction, sunrad=SUN_RADIUS_RAD, nsun=12):
    """Unit incident beam; finite uniform disk integrated over its visible area."""
    nr, na, nl = len(rgrid), len(agrid), extinction.shape[1]
    out = np.zeros((nr, na, nl, 3))
    for p in prange(nr*na):
        ir, ia = p//na, p%na
        r, a = rgrid[ir], agrid[ia]
        horizon = -math.acos(min(1., edges[0]/r))
        lower = max(-1., min(1., (horizon-a)/sunrad))
        fraction = .5-(math.asin(lower)+lower*math.sqrt(max(0.,1-lower*lower)))/math.pi
        if fraction < 1e-15:
            continue
        norm = 0.
        for k in range(nsun):
            off = lower+(1-lower)*(k+.5)/nsun
            w = math.sqrt(max(0., 1-off*off))
            norm += w
            mu = math.sin(a+off*sunrad)
            lengths, blocked = layer_lengths(r, mu, edges)
            if blocked:
                continue
            for l in range(nl):
                tau = 0.
                for j in range(len(lengths)):
                    tau += lengths[j]*extinction[j,l]
                f = math.exp(-tau)*w
                out[ir,ia,l,0] += f
                out[ir,ia,l,1] += f*max(0., mu)
        out[ir,ia] *= fraction/norm
    return out


@njit(cache=True)
def weights_path(path, edges, scattering, absorption):
    nl = scattering.shape[1]
    w = np.zeros((len(path), nl))
    trans = np.ones(nl)
    for k in range(len(path)):
        layer = min(len(edges)-2, max(0, np.searchsorted(edges, path[k,0])-1))
        for l in range(nl):
            bs = scattering[layer,l]
            chi = bs+absorption[layer,l]
            dt = chi*path[k,6]
            w[k,l] = trans[l]*(-math.expm1(-dt))*(bs/chi if chi > 0 else 0.)
            trans[l] *= math.exp(-dt)
    return w, trans


def quadrature(rgrid, radius, order, nphi):
    """Split each altitude's angular integral at its geometric ground horizon."""
    g, w = np.polynomial.legendre.leggauss(order)
    mu, mw = [], []
    for r in rgrid:
        horizon = -np.sqrt(max(0., 1-(radius/r)**2))
        # At ground level horizon coincides with zero. Keep fixed array sizes
        # by retaining the repeated bound; its angular weight is zero.
        bounds = np.sort(np.array([-1., horizon, -.12, 0., .12, 1.]))
        mu.append(np.concatenate([(g+1)*(hi-lo)/2+lo for lo,hi in zip(bounds[:-1],bounds[1:])]))
        mw.append(np.concatenate([w*(hi-lo)/2 for lo,hi in zip(bounds[:-1],bounds[1:])]))
    phi = (np.arange(nphi)+.5)*2*np.pi/nphi
    return np.array(mu), np.array(mw), np.cos(phi), np.sin(phi), 2*np.pi/nphi


def build_paths(rgrid, mus, edges, scattering, absorption, shells):
    items, maximum = [], 0
    for ir, r in enumerate(rgrid):
        row = []
        for mu in mus[ir]:
            path, ground, c, s = ray_path(r, mu, edges[0], edges[-1], shells, 1., 1.)
            w, t = weights_path(path, edges, scattering, absorption)
            row.append((path, w, t, ground, c, s))
            maximum = max(maximum, len(path))
        items.append(row)
    nr, nm, nl = len(rgrid), mus.shape[1], scattering.shape[1]
    geom = np.zeros((nr,nm,maximum,4))
    weights = np.zeros((nr,nm,maximum,nl))
    ends = np.zeros((nr,nm,nl))
    nseg = np.zeros((nr,nm), np.int32)
    gb = np.zeros((nr,nm,3))
    for i,row in enumerate(items):
        for j,(p,w,t,g,c,s) in enumerate(row):
            n = len(p)
            geom[i,j,:n] = p[:,:4]
            weights[i,j,:n] = w
            ends[i,j] = t
            nseg[i,j] = n
            gb[i,j] = [g,c,s]
    return geom, weights, ends, nseg, gb


def grids(edges, quality):
    nr, da, angular, nphi, shell_count = {
        'draft': (29, 2., 3, 12, 40),
        'standard': (45, 1., 5, 16, 72),
        'refined': (57, .75, 6, 20, 92),
        'fine': (65, .5, 7, 24, 112),
    }[quality]
    radius, height = edges[0], edges[-1]-edges[0]
    r = radius + height*np.linspace(0,1,nr)**2
    a = np.unique(np.r_[np.arange(-90,-40,5), np.arange(-40,30+da/2,da), np.arange(35,91,5)])*np.pi/180
    sr = radius + height*np.linspace(0,1,161)**2
    sa = np.unique(np.r_[np.arange(-90,-40,1.), np.arange(-40,20.0001,.125), np.arange(20.5,90.1,.5)])*np.pi/180
    shells = np.unique(np.r_[edges, radius+height*np.linspace(0,1,shell_count+1)**2])
    return r,a,sr,sa,angular,nphi,shells


def solve(optical, radius, quality='standard', albedo=.1, max_orders=100, tolerance=2e-5, path=None):
    start = time.monotonic()
    edges = radius+optical['height']
    sca = optical['scattering']/np.diff(edges)[:,None]
    absorb = optical['absorption']/np.diff(edges)[:,None]
    r,a,sr,sa,angular,nphi,shells = grids(edges,quality)
    beam = solar_table(sr,sa,edges,sca+absorb)
    mus,muw,cp,sp,dphi = quadrature(r,radius,angular,nphi)
    paths = build_paths(r,mus,edges,sca,absorb,shells)
    print('transport grids',len(r),len(a),len(optical['energy']),paths[0].shape,flush=True)
    prev = np.zeros((len(r),len(a),len(optical['energy']),7))
    total = np.zeros_like(prev)
    history = []
    roi = a >= np.deg2rad(-24)
    active = np.arange(len(optical['energy']))
    channel_orders = np.zeros(len(active),dtype=int)
    last = np.zeros_like(prev)
    active_paths,active_beam = paths,beam
    for order in range(1,max_orders+1):
        current = transport(np.ascontiguousarray(prev[:,:,active,:]),active_beam,r,a,sr,sa,
                            mus,muw,cp,sp,dphi,*active_paths,albedo,order==1)
        total[:,:,active,:] += current
        last[:,:,active,:] = current
        channel_orders[active] = order
        denominator = np.maximum(total[0,roi][:,active,5],1e-12)
        per_channel = np.max(current[0,roi,:,5]/denominator,axis=0)
        per_maximum = np.max(current[:,:,:,0],axis=(0,1))/np.maximum(
            np.max(total[:,:,active,0],axis=(0,1)),1e-30)
        rel,largest = float(per_channel.max()),float(per_maximum.max())
        history.append(dict(order=order,ground_increment_fraction=rel,maximum_increment_fraction=largest,
                            active_channels=len(active),elapsed_seconds=time.monotonic()-start))
        if order == 1 or order % 5 == 0:
            print('order',order,'ground increment',rel,'seconds',round(time.monotonic()-start,1),flush=True)
        prev[:,:,active,:] = current
        if order >= 8:
            keep = (per_channel >= tolerance)|(per_maximum >= tolerance)
            if not np.any(keep):
                active = active[:0]
                break
            if np.any(~keep):
                active = active[keep]
                active_paths = (paths[0],np.ascontiguousarray(paths[1][:,:,:,active]),
                                np.ascontiguousarray(paths[2][:,:,active]),paths[3],paths[4])
                active_beam = np.ascontiguousarray(beam[:,:,active,:])
    meta = dict(quality=quality,albedo=albedo,history=history,increment_tolerance=tolerance,
                order_cap=max_orders,order_converged=len(active)==0,channel_orders=channel_orders.tolist(),
                solar_disk_radius_deg=math.degrees(SUN_RADIUS_RAD),sun_strips=12,
                angular_order_per_interval=angular,azimuth_samples=nphi,
                seconds=time.monotonic()-start)
    result = dict(r=r,a=a,sr=sr,sa=sa,beam=beam,moments=total,last_order=last,
                  edges=edges,scattering=sca,absorption=absorb,shells=shells,**{
                      k: optical[k] for k in ('energy','xyz','wavelength','band')})
    if path:
        path = Path(path)
        path.parent.mkdir(parents=True,exist_ok=True)
        temporary = path.with_suffix('.tmp.npz')
        np.savez_compressed(temporary,**result,meta=json.dumps(meta))
        temporary.replace(path)
    return result,meta


@njit(cache=True,inline='always')
def sample_at(tab,i,j,u,v,l,c):
    return ((1-u)*((1-v)*tab[i,j,l,c]+v*tab[i,j+1,l,c])+
            u*((1-v)*tab[i+1,j,l,c]+v*tab[i+1,j+1,l,c]))


@njit(parallel=True,cache=True)
def radiance(data, beam, rg, ag, sr, sa, edges, sca, absorb, shells, albedo, suns, elevations, azimuths):
    """Formal integration for a ground observer, per unit source in each channel."""
    out = np.zeros((len(suns),len(elevations),len(azimuths),sca.shape[1]))
    for iv in prange(len(elevations)):
        mu = math.sin(elevations[iv])
        path,ground,cg,sg = ray_path(edges[0],mu,edges[0],edges[-1],shells,1.,1.)
        w,t = weights_path(path,edges,sca,absorb)
        for i,a in enumerate(suns):
            ca,ss = math.cos(a),math.sin(a)
            for j,phi in enumerate(azimuths):
                cp = math.cos(phi)
                nu = ss*mu+ca*math.cos(elevations[iv])*cp
                if ground:
                    aa = math.asin(max(-1.,min(1.,ss*cg+ca*sg*cp)))
                    for l in range(sca.shape[1]):
                        flux = interpolate(data,rg,ag,rg[0],aa,l,5)+interpolate(beam,sr,sa,rg[0],aa,l,1)
                        out[i,iv,j,l] = t[l]*albedo/math.pi*flux
                for k in range(len(path)):
                    r,qr = path[k,0],path[k,1]
                    ms = max(-1.,min(1.,ss*path[k,2]+ca*path[k,3]*cp))
                    aa = math.asin(ms)
                    qt = (nu-ms*qr)/math.sqrt(max(1e-16,1-ms*ms))
                    qp2 = max(0.,1-qr*qr-qt*qt)
                    ri,u = bracket(rg,r)
                    ai,v = bracket(ag,aa)
                    bi,bu = bracket(sr,r)
                    bj,bv = bracket(sa,aa)
                    for l in range(sca.shape[1]):
                        source = (sample_at(data,ri,ai,u,v,l,0)+
                                  qr*qr*sample_at(data,ri,ai,u,v,l,1)+
                                  qt*qt*sample_at(data,ri,ai,u,v,l,2)+
                                  qp2*sample_at(data,ri,ai,u,v,l,3)+
                                  2*qr*qt*sample_at(data,ri,ai,u,v,l,4)+
                                  (1+nu*nu)*sample_at(beam,bi,bj,bu,bv,l,0))
                        out[i,iv,j,l] += w[k,l]*3/(16*math.pi)*max(0.,source)
    return out


def evaluate(solution, meta, suns, elevations, azimuths):
    args = [solution[k] for k in ('beam','r','a','sr','sa','edges','scattering','absorption','shells')]
    return radiance(solution['moments']-solution['last_order'],*args,meta['albedo'],
                    np.radians(suns),np.radians(elevations),np.radians(azimuths))


def energy_audit(solution, meta):
    """Whole-sphere ledger normalized to the incident disk power in each channel.

    Layer absorption uses Gauss integration of r² times interpolated actinic
    flux within each homogeneous layer. Top and ground fluxes use the retained
    moments, with solar-elevation integration refined around the terminator.
    """
    r,a,sr,sa = [solution[k] for k in ('r','a','sr','sa')]
    angles = np.unique(np.r_[np.linspace(-np.pi/2,np.pi/2,1441),a])
    mu = np.sin(angles)
    moments = solution['moments']
    direct = solar_table(np.array([r[0],r[-1]]),angles,solution['edges'],
                         solution['scattering']+solution['absorption'])
    up = np.array([np.interp(angles,a,moments[-1,:,l,6]) for l in range(moments.shape[2])]).T
    down = np.array([np.interp(angles,a,moments[0,:,l,5]) for l in range(moments.shape[2])]).T
    escape = 2*np.trapezoid(direct[1,:,:,0]*np.maximum(-mu,0)[:,None],mu,axis=0)
    reflected = 2*np.trapezoid(up,mu,axis=0)
    ground = 2*(r[0]/r[-1])**2*(1-meta['albedo'])*np.trapezoid(down+direct[0,:,:,1],mu,axis=0)
    # Integrate in altitude first, splitting at both layer and moment-grid boundaries.
    gauss,weights = np.polynomial.legendre.leggauss(3)
    splits = np.unique(np.r_[solution['edges'],r])
    radii = ((splits[1:,None]+splits[:-1,None])/2+
             np.diff(splits)[:,None]*gauss/2).ravel()
    rw = (np.diff(splits)[:,None]*weights/2).ravel()
    absorbed = np.zeros(moments.shape[2])
    for rad,w in zip(radii,rw):
        ir,u = bracket(r,rad)
        layer = min(len(solution['edges'])-2,np.searchsorted(solution['edges'],rad)-1)
        j = (1-u)*moments[ir,:,:,0]+u*moments[ir+1,:,:,0]
        bi,bu = bracket(sr,rad)
        f = (1-bu)*solution['beam'][bi,:,:,0]+bu*solution['beam'][bi+1,:,:,0]
        diffuse_angular = np.trapezoid(j,np.sin(a),axis=0)
        beam_angular = np.trapezoid(f,np.sin(sa),axis=0)
        absorbed += 2*w*(rad/r[-1])**2*solution['absorption'][layer]*(diffuse_angular+beam_angular)
    residual = 1-escape-reflected-ground-absorbed
    return dict(direct_escape=escape.tolist(),diffuse_escape=reflected.tolist(),ground_absorption=ground.tolist(),
                atmospheric_absorption=absorbed.tolist(),residual=residual.tolist(),
                maximum_absolute_channel_residual=float(abs(residual).max()),
                spectral_energy_weighted_residual=float(residual@solution['energy']/solution['energy'].sum()),
                photopic_weighted_residual=float(residual@solution['xyz'][:,1]/solution['xyz'][:,1].sum()))
