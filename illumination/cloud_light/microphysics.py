"""Visible extinction from the recorded Morrison mass and number distributions.

The gamma/exponential distributions and size bounds follow CM1 r22.0's
morrison.F. Q_ext=2 and nonabsorbing particles are visible geometric-optics
approximations. Ice particles use the scheme's equivalent spheres; phase
functions remain explicit sensitivity choices.
"""
from __future__ import annotations

import numpy as np

SPECIES = ('qc', 'qr', 'qi', 'qs', 'qg')
# Model coefficients from the admitted CM1 Morrison source, in SI units.
PARTICLE_DENSITY = dict(qc=997., qr=997., qi=500., qs=100., qg=400.)
SLOPE_LENGTH_LIMITS = dict(qr=(20e-6,2800e-6), qi=(1e-6,350e-6),
                           qs=(10e-6,2000e-6), qg=(20e-6,2000e-6))
NUMBERS = dict(qr='ncr', qi='nci', qs='ncs', qg='ncg')


def effective_radius(mass, number_per_kg, density, species, droplets_cm3=100.):
    """Area/volume moment radius, after the producing scheme's slope limits."""
    q, rho = np.broadcast_arrays(np.asarray(mass,float), np.asarray(density,float))
    if species not in SPECIES or np.any(q < 0) or np.any(rho <= 0):
        raise ValueError('Invalid hydrometeor column')
    particle = PARTICLE_DENSITY[species]
    if species == 'qc':
        number = droplets_cm3*1e6/rho
        shape = np.clip(1/(.0005714*droplets_cm3+.2714)**2-1,2,10)
        slope = (np.pi*particle/6*number*(shape+1)*(shape+2)*(shape+3)/np.maximum(q,1e-30))**(1/3)
        slope = np.clip(slope,(shape+1)/60e-6,(shape+1)/1e-6)
        radius = (shape+3)/(2*slope)
    else:
        number = np.broadcast_to(np.asarray(number_per_kg,float),q.shape)
        if np.any(number < 0):
            raise ValueError('Negative hydrometeor number')
        short,long = SLOPE_LENGTH_LIMITS[species]
        slope = np.clip((np.pi*particle*number/np.maximum(q,1e-30))**(1/3),1/long,1/short)
        radius = 1.5/slope
    return radius


def extinction(fields,density,droplets_cm3=100.):
    result,radii = {},{}
    for species in SPECIES:
        q = np.asarray(fields[species],float)
        radius = effective_radius(q, fields.get(NUMBERS.get(species,'')),density,species,droplets_cm3)
        result[species] = np.where(q>=1e-14,3*density*q/(2*PARTICLE_DENSITY[species]*radius),0.)
        radii[species] = radius
    return result,radii
