"""Finite angular-source drivers for the existing cloud transport kernel.

The historical solar-only volume.py remains byte-stable. This extension reuses
its complete photon-path integrator and changes only source-disk sampling in
the radiance/illuminance drivers. The seeded driver loop follows that baseline.
"""
from __future__ import annotations
import math
import numpy as np
from numba import njit, prange
from illumination.cloud_light.volume import trace, cloud_depth
from illumination.sky.reference_mc import rotate_direction
from illumination.sky.solved_transport import layer_lengths, SUN_RADIUS_RAD

@njit(parallel=True,cache=True)
def estimate(origin,directions,sun,edges,scatter,absorb,xyz,cloud_edges,theta0,dtheta,beta,half_width,
             albedo,liquid_g,ice_g,photons=1024,groups=16,seed=77,flux=False,max_events=100000,
             flux_mu_min=0.,flux_mu_max=1.,flux_az_min=0.,flux_az_max=2*math.pi,
             source_radius_rad=SUN_RADIUS_RAD):
    """Radiance XYZ (Y cd/m2), or cosine-sampled horizontal illuminance XYZ."""
    if not math.isfinite(source_radius_rad) or source_radius_rad < 0 or source_radius_rad >= math.pi/2:
        raise ValueError('Source angular radius must be finite and in [0, pi/2)')
    if photons < groups or photons % groups:
        raise ValueError("Photons must be a positive multiple of the block count")
    if flux and (flux_mu_min < 0 or flux_mu_max > 1 or flux_mu_min >= flux_mu_max or flux_az_min >= flux_az_max):
        raise ValueError("Invalid horizontal aperture")
    counts = photons//groups
    n = len(directions)
    result = np.zeros((n,groups,3))
    truncated = np.zeros((n,groups),np.int64)
    source = xyz.sum(axis=1)
    probability = source/source.sum()
    cdf = np.cumsum(probability)
    for index in prange(n*groups):
        i,g = index//groups,index%groups
        for photon in range(counts):
            # Independent path seeds also pair cloud and clear control histories.
            np.random.seed(((seed+i*photons+g*counts+photon)*2654435761) % 4294967296)
            l = min(len(cdf)-1,np.searchsorted(cdf,np.random.random()))
            d = directions[i].copy()
            if flux:
                mu = math.sqrt(flux_mu_min**2+(flux_mu_max**2-flux_mu_min**2)*np.random.random())
                phi = flux_az_min+(flux_az_max-flux_az_min)*np.random.random()
                d = np.array([math.sqrt(1-mu*mu)*math.cos(phi),math.sqrt(1-mu*mu)*math.sin(phi),mu])
            value,cut = trace(origin,d,sun,edges,scatter[:,l],absorb[:,l],cloud_edges,theta0,dtheta,beta,half_width,
                              albedo,liquid_g,ice_g,source_radius_rad,max_events)
            for c in range(3):
                result[i,g,c] += value*xyz[l,c]/probability[l]*((flux_az_max-flux_az_min)*(flux_mu_max**2-flux_mu_min**2)/2 if flux else 1)/counts
            truncated[i,g] += cut
    return result,truncated


@njit(cache=True)
def direct_horizontal(origin,sun,edges,extinction,xyz,cloud_edges,theta0,dtheta,beta,half_width,
                      radial_samples=8,azimuth_samples=32,source_radius_rad=SUN_RADIUS_RAD):
    """Finite solar disk incident on a horizontal plane, including cloud shadows."""
    result=np.zeros(3)
    if not math.isfinite(source_radius_rad) or source_radius_rad < 0 or source_radius_rad >= math.pi/2:
        raise ValueError('Source angular radius must be finite and in [0, pi/2)')
    if sun[2]+math.sin(source_radius_rad)<=0:
        return result
    for radial in range(radial_samples):
        cosine=math.sqrt(1-math.sin(source_radius_rad)**2*(radial+.5)/radial_samples)
        for azimuth in range(azimuth_samples):
            s=np.array(rotate_direction(sun[0],sun[1],sun[2],cosine,
                                       2*math.pi*(azimuth+.5)/azimuth_samples))
            mu=s[2]
            if mu<=0:
                continue
            radius=math.sqrt(np.dot(origin,origin))
            lengths,blocked=layer_lengths(radius,np.dot(origin,s)/radius,edges)
            if blocked:
                continue
            tau_cloud=cloud_depth(origin,s,cloud_edges,theta0,dtheta,beta,half_width)
            transmission=np.exp(-(lengths@extinction+tau_cloud))
            result += mu*(transmission@xyz)/(radial_samples*azimuth_samples)
    return result
