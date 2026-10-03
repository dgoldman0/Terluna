"""Backward Monte Carlo checks of the solved layered molecular atmosphere.

Null collisions sample the continuous spherical paths. Next-event sunlight
uses exact shell columns, independent of the deterministic source grids. A
point Sun is appropriate for the chosen spots away from solar-disk contact.
"""
import math

import numpy as np
from numba import njit, prange

from illumination.sky.reference_mc import rotate_direction
from illumination.sky.solved_transport import layer_lengths


@njit(cache=True)
def solar_transmission(r,mu,edges,extinction):
    lengths,blocked = layer_lengths(r,mu,edges)
    if blocked:
        return 0.
    return math.exp(-np.dot(lengths,extinction))


@njit(parallel=True,cache=True)
def monte_carlo(edges,scatter,absorb,albedo,sun,elevation,azimuth,groups,per_group,seed,flux=False,cdf=None):
    means = np.zeros(groups)
    truncated = np.zeros(groups,dtype=np.int64)
    all_extinction = scatter+absorb
    radius,top = edges[0],edges[-1]
    sx,sz = math.cos(sun),math.sin(sun)
    for group in prange(groups):
        np.random.seed(seed+103*group)
        for photon in range(per_group):
            if scatter.ndim == 2:
                channel = min(len(cdf)-1,np.searchsorted(cdf,np.random.random()))
                extinction,scattering = all_extinction[channel],scatter[channel]
            else:
                extinction,scattering = all_extinction,scatter
            px,py,pz = 0.,0.,radius+.01
            if flux:
                dz = math.sqrt(np.random.random())
                phi = 2*math.pi*np.random.random()
                dx,dy = math.sqrt(1-dz*dz)*math.cos(phi),math.sqrt(1-dz*dz)*math.sin(phi)
            else:
                dx,dy,dz = math.cos(elevation)*math.cos(azimuth),math.cos(elevation)*math.sin(azimuth),math.sin(elevation)
            weight,value,events = 1.,0.,0
            while weight>1e-12 and events<1000:
                events += 1
                rr = px*px+py*py+pz*pz
                r = math.sqrt(rr)
                rd = px*dx+py*dy+pz*dz
                disc = rd*rd-rr+radius*radius
                ground = rd<0 and disc>0
                end = max(0.,-rd-math.sqrt(disc)) if ground else max(0.,-rd+math.sqrt(max(0.,rd*rd-rr+top*top)))
                rmin = radius if ground else (math.sqrt(max(radius*radius,rr-rd*rd)) if rd<0 else r)
                first = max(0,min(len(scattering)-1,np.searchsorted(edges,rmin)-1))
                majorant = np.max(extinction[first:])
                t,hit = 0.,False
                while majorant>0 and t<end:
                    t += -math.log(max(1e-16,np.random.random()))/majorant
                    if t>=end:
                        break
                    qx,qy,qz = px+t*dx,py+t*dy,pz+t*dz
                    qr = math.sqrt(qx*qx+qy*qy+qz*qz)
                    layer = max(0,min(len(scattering)-1,np.searchsorted(edges,qr)-1))
                    if np.random.random()<extinction[layer]/majorant:
                        hit = True
                        px,py,pz,r = qx,qy,qz,qr
                        omega = scattering[layer]/extinction[layer]
                        mu_sun = (px*sx+pz*sz)/r
                        nu = dx*sx+dz*sz
                        value += weight*omega*3/(16*math.pi)*(1+nu*nu)*solar_transmission(r,mu_sun,edges,extinction)
                        weight *= omega
                        while True:
                            cosine = 2*np.random.random()-1
                            if np.random.random()<.5*(1+cosine*cosine):
                                break
                        dx,dy,dz = rotate_direction(dx,dy,dz,cosine,2*math.pi*np.random.random())
                        break
                if not hit:
                    if not ground:
                        break
                    px,py,pz = px+dx*end,py+dy*end,pz+dz*end
                    r = math.sqrt(px*px+py*py+pz*pz)
                    nx,ny,nz = px/r,py/r,pz/r
                    mu_sun = nx*sx+nz*sz
                    if mu_sun>0:
                        value += weight*albedo/math.pi*mu_sun*solar_transmission(radius+.01,mu_sun,edges,extinction)
                    weight *= albedo
                    dx,dy,dz = rotate_direction(nx,ny,nz,math.sqrt(np.random.random()),2*math.pi*np.random.random())
                    px,py,pz = nx*(radius+.01),ny*(radius+.01),nz*(radius+.01)
                if weight<.02:
                    if np.random.random()>.2:
                        break
                    weight /= .2
            if events>=1000:
                truncated[group] += 1
            means[group] += value*(math.pi if flux else 1.)/per_group
    return means,truncated


def check(solution,meta,channel,sun_deg,elevation_deg=90.,azimuth_deg=90.,photons=80000,seed=81,flux=False):
    groups = 40
    if photons<groups or photons%groups:
        raise ValueError('Photon count must be a positive multiple of 40')
    means,truncated = monte_carlo(solution['edges'],solution['scattering'][:,channel],
                                 solution['absorption'][:,channel],meta['albedo'],math.radians(sun_deg),
                                 math.radians(elevation_deg),math.radians(azimuth_deg),groups,photons//groups,seed,flux)
    return dict(channel=channel,wavelength_nm=float(solution['wavelength'][channel]),sun_deg=sun_deg,
                elevation_deg=elevation_deg,azimuth_deg=azimuth_deg,quantity='diffuse_horizontal' if flux else 'radiance',
                photons=photons,seed=seed,mean=float(means.mean()),standard_error=float(means.std(ddof=1)/math.sqrt(groups)),
                truncated_paths=int(truncated.sum()))


def spectral_flux_check(fine,radius,sun_deg,albedo=.1,photons=240000,seed=117):
    """Sample wavelength from the fine-grid photopic source for every path.

    This measures the broadband diffuse light with the unreduced line spectra;
    it bypasses both the angular/source grids and the reduced spectral bins.
    """
    if photons<40 or photons%40:
        raise ValueError('Photon count must be a positive multiple of 40')
    source = fine['xyz'][:,1]
    cdf = np.cumsum(source)/source.sum()
    dz = np.diff(fine['height'])
    sca = np.ascontiguousarray((fine['scattering']/dz[:,None]).T)
    absorb = np.ascontiguousarray((fine['absorption']/dz[:,None]).T)
    means,truncated = monte_carlo(radius+fine['height'],sca,absorb,albedo,math.radians(sun_deg),
                                 0.,0.,40,photons//40,seed,True,cdf)
    means *= source.sum()
    return dict(sun_deg=sun_deg,photons=photons,seed=seed,mean_diffuse_lux=float(means.mean()),
                standard_error_lux=float(means.std(ddof=1)/math.sqrt(40)),truncated_paths=int(truncated.sum()),
                fine_spectral_samples=len(source),quantity='fine_grid_photopic_diffuse_horizontal')
