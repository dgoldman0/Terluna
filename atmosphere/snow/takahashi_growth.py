"""Takahashi et al. (1991) Table 2: empirical single-crystal growth fits.

Original printed p.24, full a/c lengths in mm, mass in g, age in min.
These are fits without pointwise errors, not raw measurements. Planar c was
inferred from mass and basal area. Independent fits need not map to a spheroid
of admissible density. Reject that mapping explicitly; never clip it silently.
"""
import math
from atmosphere.snow import scale_crystal_growth as original

SOURCE_SHA = original.SOURCE_SHA
# temperature: (end age, habit, mass segments, equatorial segments, polar segments)
# A segment is (upper age, coefficient, exponent); the next fit applies at a join.
FITS = {
    -3.7: (20, 'thick plate', [(20,1.3e-8,1.40)], [(20,.031,.46)], [(20,.033,.42)]),
    -5.3: (30, 'hollow column to sheaths/needle', [(30,8.4e-9,2.)], [(30,.040,.36)], [(6,.0097,1.95),(30,.036,1.16)]),
    -8.6: (15, 'thick plate', [(15,2.5e-8,1.52)], [(15,.037,.54)], [(15,.029,.59)]),
    -10.6: (12, 'thick plate', [(12,3.3e-8,1.51)], [(12,.061,.48)], [(12,.023,.53)]),
    -12.2: (20, 'plate to sector', [(7,6.3e-8,1.28),(20,1.1e-8,2.16)], [(7,.12,.37),(20,.041,.92)], [(20,.0085,.44)]),
    -14.4: (30, 'broad branch to dendrite', [(4,6.4e-8,1.63),(30,2.3e-8,2.37)], [(4,.14,.95),(30,.13,1.02)], [(30,.0044,.52)]),
    -16.5: (20, 'sector', [(5,8.1e-8,1.36),(20,3.5e-8,1.91)], [(5,.13,.54),(20,.10,.74)], [(20,.0058,.48)]),
    -18.2: (15, 'plate to plate with simple extensions', [(15,4.2e-8,1.49)], [(15,.067,.64)], [(15,.017,.21)]),
    -20.1: (12, 'thick plate', [(12,3.8e-8,1.44)], [(12,.054,.54)], [(12,.026,.44)]),
    -22.: (12, 'thick plate', [(12,3.4e-8,1.42)], [(12,.040,.49)], [(12,.036,.48)]),
}


def observed(temperature_C, age_min):
    if temperature_C not in FITS:
        raise ValueError('No fit at this temperature; interpolation is not measured')
    end,habit,ms,aa,cc=FITS[temperature_C]
    if not math.isfinite(age_min) or not 3 <= age_min <= end:
        raise ValueError('Outside the selected measured-age interval')
    def value(segments):
        for upper,coefficient,power in segments:
            if age_min < upper:
                return coefficient*age_min**power
        return segments[-1][1]*age_min**segments[-1][2]
    mass,a,c=(value(v)*1e-3 for v in (ms,aa,cc))
    return dict(age_min=age_min,mass_kg=mass,equatorial_dimension_m=a,
                polar_dimension_m=c,diameter_m=max(a,c),habit=habit)


def initial(temperature_C, bulk_ice_density):
    obs=observed(temperature_C,3.)
    a=obs['equatorial_dimension_m']/2;c=obs['polar_dimension_m']/2
    rho=obs['mass_kg']/(4*math.pi/3*a*a*c)
    if rho > bulk_ice_density:
        raise ValueError(f'Fitted mass/axes require {rho:.6f} kg/m3 in a spheroid, above bulk ice {bulk_ice_density}')
    return dict(**obs,equatorial_radius_m=a,polar_radius_m=c,density_kg_m3=rho)


# Table 4, p.25: w in cm/s; each relationship represents ~4 measurements.
# Segments store upper age, coefficient, exponent, additive constant.
FALL_SPEED = {
    -3.7: [(11,1.5,1,0),(20,.89,1,6.4)],
    -5.3: [(7,1.7,1,0),(30,.73,1,6.7)],
    -8.6: [(8,2.5,1,0),(15,1.3,1,9.4)],
    -10.6: [(12,2.4,1,0)],
    -12.2: [(7,1.8,1,0),(20,1.5,1,2.2)],
    -14.4: [(6,2.8,.72,0),(30,.82,1,5.1)],
    -16.5: [(9,1.9,1,0),(20,1.2,1,5.9)],
    -18.2: [(6,2.2,1,0),(15,1.3,1,5.3)],
    -20.1: [(7,2.1,1,0),(12,1.3,1,5.9)],
    -22.: [(6,2.2,1,0),(12,1.3,1,5.0)],
}


def fall_speed(temperature_C,age_min):
    observed(temperature_C,age_min)  # common validity check
    for upper,coefficient,power,offset in FALL_SPEED[temperature_C]:
        if age_min < upper:break
    return (coefficient*age_min**power+offset)*.01
