"""Spectral Monte Carlo through molecular air and a resolved cloud section.

Cloud cells occupy spherical shells and great-circle angular columns, extruded
through a stated cross-section width. Both molecular and cloud multiple
scattering, cloud shadows and foreground air are sampled. A finite uniform Sun
is sampled at each next-event estimate; the ground is Lambertian.
"""
from __future__ import annotations

import math
import numpy as np
from numba import njit, prange

from illumination.sky.reference_mc import rotate_direction
from illumination.sky.solved_transport import layer_lengths, SUN_RADIUS_RAD

EPS = .002


@njit(cache=True)
def phase_hg(cosine,g):
    return (1-g*g)/(4*math.pi*(1+g*g-2*g*cosine)**1.5)


@njit(cache=True)
def sample_hg(u,g):
    if abs(g)<1e-8:
        return 2*u-1
    return max(-1.,min(1.,(1+g*g-((1-g*g)/(1-g+2*g*u))**2)/(2*g)))


@njit(cache=True)
def sphere_distance(p,d,radius):
    rd = p[0]*d[0]+p[1]*d[1]+p[2]*d[2]
    discriminant = rd*rd-(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])+radius*radius
    if discriminant<=0:
        return 1e30
    root = math.sqrt(discriminant)
    a,b = -rd-root,-rd+root
    return a if a>EPS*.25 else b if b>EPS*.25 else 1e30


@njit(cache=True)
def cell(p,edges,cloud_edges,theta0,dtheta,beta,half_width):
    radius = math.sqrt(np.dot(p,p))
    m = max(0,min(len(edges)-2,np.searchsorted(edges,radius)-1))
    k = np.searchsorted(cloud_edges,radius)-1
    j = int(math.floor((math.atan2(p[0],p[2])-theta0)/dtheta))
    liquid,ice = 0.,0.
    if p[2]>0 and abs(p[1])<half_width and 0<=k<beta.shape[0] and 0<=j<beta.shape[1]:
        liquid,ice = beta[k,j,0],beta[k,j,1]
    return m,k,j,liquid,ice


@njit(cache=True)
def boundary(p,d,edges,cloud_edges,theta0,dtheta,nx,half_width):
    # Probe the outgoing side at an interface; every segment has constant optics.
    q = p+EPS*d
    r = math.sqrt(np.dot(q,q))
    m = max(0,min(len(edges)-2,np.searchsorted(edges,r)-1))
    k = max(0,min(len(cloud_edges)-2,np.searchsorted(cloud_edges,r)-1))
    distance = min(sphere_distance(p,d,edges[m]),sphere_distance(p,d,edges[m+1]),
                   sphere_distance(p,d,cloud_edges[k]),sphere_distance(p,d,cloud_edges[k+1]))
    angle = math.atan2(q[0],q[2])
    j = max(0,min(nx-1,int(math.floor((angle-theta0)/dtheta))))
    for edge in (j,j+1):
        a = theta0+edge*dtheta
        denominator = d[0]*math.cos(a)-d[2]*math.sin(a)
        if abs(denominator)>1e-15:
            t = -(p[0]*math.cos(a)-p[2]*math.sin(a))/denominator
            if t>EPS*.25:
                distance = min(distance,t)
    if abs(d[1])>1e-15:
        for y in (-half_width,half_width):
            t = (y-p[1])/d[1]
            if t>EPS*.25:
                distance = min(distance,t)
    return distance


@njit(cache=True)
def cloud_depth(p,d,edges,theta0,dtheta,beta,half_width):
    """Exact piecewise-constant cloud optical depth along a straight ray."""
    q = p.copy()
    tau = 0.
    limits = np.array([edges[0],edges[-1]])
    # Tangency and view blockage are handled by the molecular shell calculation.
    for segment in range(10000):
        r = math.sqrt(np.dot(q,q))
        if r<edges[0]-.01 or (r>edges[-1]+EPS and np.dot(q,d)>=0):
            break
        _,_,_,liquid,ice = cell(q+EPS*d,limits,edges,theta0,dtheta,beta,half_width)
        length = boundary(q,d,limits,edges,theta0,dtheta,beta.shape[1],half_width)
        if length>1e20:
            break
        tau += (liquid+ice)*length
        if tau>40:
            return tau
        q += (length+EPS)*d
    else:
        raise RuntimeError("Cloud optical-depth ray exceeded the segment limit")
    return tau


@njit(cache=True)
def sunlight(p,sun,edges,extinction,cloud_edges,theta0,dtheta,beta,half_width):
    r = math.sqrt(np.dot(p,p))
    lengths,blocked = layer_lengths(r,np.dot(p,sun)/r,edges)
    if blocked:
        return 0.
    tau = np.dot(lengths,extinction)
    if tau<40:
        tau += cloud_depth(p,sun,cloud_edges,theta0,dtheta,beta,half_width)
    return math.exp(-tau)


@njit(cache=True)
def disk_direction(sun,radius):
    cosine = math.sqrt(max(0.,1-math.sin(radius)**2*np.random.random()))
    return np.array(rotate_direction(sun[0],sun[1],sun[2],cosine,2*math.pi*np.random.random()))


@njit(cache=True)
def trace(origin,direction,sun,edges,scatter,absorb,cloud_edges,theta0,dtheta,beta,half_width,
          albedo,liquid_g,ice_g,sun_radius,max_events):
    p,d = origin.copy(),direction.copy()
    weight,value = 1.,0.
    ext = scatter+absorb
    for event in range(max_events):
        optical = -math.log(max(1e-16,np.random.random()))
        hit_ground = False
        escaped = False
        liquid,ice,m = 0.,0.,0
        for segment in range(10000):
            probe = p+EPS*d
            radius = math.sqrt(np.dot(probe,probe))
            if radius<=edges[0]:
                hit_ground = True
                break
            if radius>=edges[-1] and np.dot(p,d)>=0:
                escaped = True
                break
            m,_,_,liquid,ice = cell(probe,edges,cloud_edges,theta0,dtheta,beta,half_width)
            length = boundary(p,d,edges,cloud_edges,theta0,dtheta,beta.shape[1],half_width)
            coefficient = ext[m]+liquid+ice
            if coefficient>0 and optical<coefficient*length:
                p += optical/coefficient*d
                break
            optical -= coefficient*length
            p += (length+EPS)*d
        else:
            return value,2
        if escaped:
            return value,0
        s = disk_direction(sun,sun_radius)
        if hit_ground:
            r = math.sqrt(np.dot(p,p))
            normal = p/r
            p = normal*(edges[0]+.01)
            mu = np.dot(normal,s)
            if mu>0:
                value += weight*albedo/math.pi*mu*sunlight(p,s,edges,ext,cloud_edges,theta0,dtheta,beta,half_width)
            weight *= albedo
            d = np.array(rotate_direction(normal[0],normal[1],normal[2],math.sqrt(np.random.random()),2*math.pi*np.random.random()))
        else:
            total = scatter[m]+liquid+ice
            coefficient = ext[m]+liquid+ice
            cosine = np.dot(d,s)
            phase = scatter[m]*3/(16*math.pi)*(1+cosine*cosine)+liquid*phase_hg(cosine,liquid_g)+ice*phase_hg(cosine,ice_g)
            value += weight*phase/coefficient*sunlight(p,s,edges,ext,cloud_edges,theta0,dtheta,beta,half_width)
            weight *= total/coefficient
            if total<=0:
                return value,0
            component = np.random.random()*total
            if component<scatter[m]:
                while True:
                    cosine = 2*np.random.random()-1
                    if np.random.random()<.5*(1+cosine*cosine):
                        break
            else:
                cosine = sample_hg(np.random.random(),liquid_g if component<scatter[m]+liquid else ice_g)
            d = np.array(rotate_direction(d[0],d[1],d[2],cosine,2*math.pi*np.random.random()))
        if weight<.05:
            if np.random.random()>.2:
                return value,0
            weight /= .2
    return value,1


@njit(parallel=True,cache=True)
def estimate(origin,directions,sun,edges,scatter,absorb,xyz,cloud_edges,theta0,dtheta,beta,half_width,
             albedo,liquid_g,ice_g,photons=1024,groups=16,seed=77,flux=False,max_events=100000,
             flux_mu_min=0.,flux_mu_max=1.,flux_az_min=0.,flux_az_max=2*math.pi):
    """Radiance XYZ (Y cd/m2), or cosine-sampled horizontal illuminance XYZ."""
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
                              albedo,liquid_g,ice_g,SUN_RADIUS_RAD,max_events)
            for c in range(3):
                result[i,g,c] += value*xyz[l,c]/probability[l]*((flux_az_max-flux_az_min)*(flux_mu_max**2-flux_mu_min**2)/2 if flux else 1)/counts
            truncated[i,g] += cut
    return result,truncated


@njit(cache=True)
def direct_horizontal(origin,sun,edges,extinction,xyz,cloud_edges,theta0,dtheta,beta,half_width,
                      radial_samples=8,azimuth_samples=32):
    """Finite solar disk incident on a horizontal plane, including cloud shadows."""
    result=np.zeros(3)
    if sun[2]+math.sin(SUN_RADIUS_RAD)<=0:
        return result
    for radial in range(radial_samples):
        cosine=math.sqrt(1-math.sin(SUN_RADIUS_RAD)**2*(radial+.5)/radial_samples)
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
