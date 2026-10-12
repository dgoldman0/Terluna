"""Aerophytes as a resource and a hazard to people: membrane, food, fibre, shade, falls and tethers.

Pure functions plus evaluate(). Sky-ship hull sizes come from the sky-ship
product; harvest shares and yields are design guesses stated with their range.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from shared.constants import AU, JULIAN_YEAR_DAYS, MOON_GM, MOON_RADIUS, SUN_RADIUS

ROOT = Path(__file__).resolve().parents[3]
INPUT_FILES = ('research/studies/sky_ships/results/sky_ships.json',
               'research/studies/aerophytes/resources_sources.json')
SUN_ANGULAR_DIAMETER_RAD = 2*math.asin(SUN_RADIUS/AU)
FOOD_KCAL_PER_KG_DRY = 4000.  # aerial_ecology convention
PERSON_KCAL_PER_DAY = 2500.   # aerial_ecology convention
GOLDBEATER_FABRIC_KG_M2 = (0.13, 0.15)  # finished fabric, Chollet 1922 (envelope_sources.json)
SKY_SHIP_CELLS = 16  # the sky-ship study's gas-cell count (one of sixteen deflated)


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


def shadow(body_diameter_m, height_m, direct_share=.49,
           sun_diameter_rad=SUN_ANGULAR_DIAMETER_RAD):
    """Ground shadow of a body overhead: umbra width, penumbra and the light it removes.

    The umbra narrows by height x the Sun's angular diameter; a body narrower
    than that casts only a soft dimming. In the umbra the direct beam (49% of
    the light at the ground with the Sun overhead) is gone and the sky remains.
    """
    _positive(diameter=body_diameter_m, height=height_m, sun=sun_diameter_rad)
    if not 0 <= direct_share <= 1:
        raise ValueError('direct share in [0, 1]')
    umbra = body_diameter_m - height_m*sun_diameter_rad
    penumbra_outer = body_diameter_m + height_m*sun_diameter_rad
    peak = direct_share*min(1., (body_diameter_m/(height_m*sun_diameter_rad))**2)
    return dict(umbra_diameter_m=max(0., umbra), penumbra_diameter_m=penumbra_outer,
                umbra_length_m=body_diameter_m/sun_diameter_rad,
                full_shadow=umbra > 0, peak_light_removed=peak)


def gas_cell_area(length_m, diameter_m, cells=SKY_SHIP_CELLS):
    """Film area of a rigid ship's gas cells: cylinders between bulkheads."""
    _positive(length=length_m, diameter=diameter_m)
    if cells < 1:
        raise ValueError('at least one cell')
    r = diameter_m/2
    one = 2*math.pi*r*r + 2*math.pi*r*length_m/cells
    return cells*one


def membrane_harvest(envelope_area_m2, renewal_years, harvest_share):
    """Envelope film area an aerophyte sheds or yields each year at a harvest share of its renewal."""
    _positive(area=envelope_area_m2, renewal=renewal_years)
    if not 0 <= harvest_share <= 1:
        raise ValueError('harvest share in [0, 1]')
    return envelope_area_m2/renewal_years*harvest_share


def food_people(projected_area_m2, active_dry_kg_m2, turnover_per_year, harvest_share):
    """People fed a year by harvesting a share of an aerophyte's renewed active tissue."""
    _nonnegative(area=projected_area_m2, active=active_dry_kg_m2, turnover=turnover_per_year)
    if not 0 <= harvest_share <= 1:
        raise ValueError('harvest share in [0, 1]')
    dry = projected_area_m2*active_dry_kg_m2*turnover_per_year*harvest_share
    return dict(dry_kg_year=dry, people=dry*FOOD_KCAL_PER_KG_DRY/(PERSON_KCAL_PER_DAY*JULIAN_YEAR_DAYS))


def draped_fall(areal_mass_kg_m2, air_density_kg_m3, drag_coefficient=1.3, height_m=0.):
    """Sink speed of a deflated envelope falling flat, like a canopy, and its fall time from height.

    One density for the whole fall; fall_minutes integrates through a density profile.
    """
    _positive(areal=areal_mass_kg_m2, density=air_density_kg_m3, cd=drag_coefficient)
    g = MOON_GM/(MOON_RADIUS+height_m)**2
    speed = math.sqrt(2*areal_mass_kg_m2*g/(air_density_kg_m3*drag_coefficient))
    return dict(sink_m_s=speed, minutes_from_height=height_m/speed/60 if height_m else 0.)


def fall_minutes(areal_mass_kg_m2, height_m, heights_km, densities_kg_m3, drag_coefficient=1.3, steps=400):
    """Minutes for a flat-falling envelope to reach the ground through a density profile.

    The sink speed sqrt(2 m g/(rho Cd)) follows the local density (log-linear in
    height between the profile's points) and gravity.
    """
    _positive(areal=areal_mass_kg_m2, height=height_m, cd=drag_coefficient)
    zs = [z*1000 for z in heights_km]
    lr = [math.log(r) for r in densities_kg_m3]
    total = 0.
    for i in range(steps):
        z = (i+.5)*height_m/steps
        rho = math.exp(float(np.interp(z, zs, lr)))
        g = MOON_GM/(MOON_RADIUS+z)**2
        total += height_m/steps/math.sqrt(2*areal_mass_kg_m2*g/(rho*drag_coefficient))
    return total/60


def evaluate(root=ROOT):
    """Shadows, gas-cell supply against sky ships, food and fibre per body, falls and tethers."""
    root = Path(root)
    ships = json.loads((root/INPUT_FILES[0]).read_text())
    if ships.get('schema') != 'terluna.research.sky-ships/1':
        raise ValueError('unexpected sky-ships product')
    rho = {a['height_km']: a['density_kg_m3'] for a in ships['air']}
    shadows = [dict(diameter_m=d, height_km=h, **shadow(d, h*1000)) for d in (40., 200., 1000., 2000.) for h in (10., 20.)]
    cells = [dict(length_m=s['length_m'], diameter_m=s['diameter_m'],
                  gas_cell_area_m2=gas_cell_area(s['length_m'], s['diameter_m']))
             for s in ships['operations']['ships']]
    bodies = []
    for name, envelope, projected in (('module_40m', 4*math.pi*20**2, math.pi*20**2),
                                      ('round_200m', 1.1*4*math.pi*100**2, math.pi*100**2),
                                      ('colony_2km', 3.63*math.pi*1000**2, math.pi*1000**2)):
        for renewal in (10., 30.):
            for share in (.1, .5):
                yearly = membrane_harvest(envelope, renewal, share)
                bodies.append(dict(body=name, envelope_area_m2=envelope, renewal_years=renewal,
                                   harvest_share=share, film_m2_year=yearly,
                                   ships_of_largest_class_per_year=yearly/cells[-1]['gas_cell_area_m2'],
                                   ships_of_smallest_class_per_year=yearly/cells[0]['gas_cell_area_m2']))
    food = [dict(body=n, harvest_share=s, **food_people(a, 1.1, .5, s))
            for n, a in (('module_40m', math.pi*20**2), ('round_200m', math.pi*100**2), ('colony_2km', math.pi*1000**2))
            for s in (.05, .2)]
    levels = sorted(rho)
    falls = [dict(areal_mass_kg_m2=m, height_km=h, sink_at_ground_m_s=draped_fall(m, rho[0])['sink_m_s'],
                  minutes_from_height=fall_minutes(m, h*1000, levels, [rho[z] for z in levels]))
             for m in (3., 10., 25.) for h in (10., 20.)]
    return dict(
        evidence='Geometry of shadows and gas cells, sky-ship hull sizes from the sky-ship product, food at the aerial-ecology convention. Harvest shares and renewal intervals are design guesses.',
        reading_rule='Shadows for a body overhead with the Sun overhead; light removed is the share of ground light lost at the shadow centre. Film areas in m² a year; ship equivalents count the gas-cell film of one ship.',
        sun_angular_diameter_deg=math.degrees(SUN_ANGULAR_DIAMETER_RAD),
        shadows=shadows, sky_ship_cells=cells, membrane=bodies, food=food, falls=falls,
        goldbeater_fabric_kg_m2=list(GOLDBEATER_FABRIC_KG_M2),
        traffic_bands_km=dict(wing_band=[.03, .3], sky_boats=[.3, 1.5], regional_ships=[22., 27.],
                              liners=[35., 45.], winged_liners=[55., 55.], platforms=[70., 70.],
                              summit_tower_top=[35.6, 35.6]),
        unresolved=['Processing a living film into gas-cell fabric, and its hydrogen permeability after harvest.',
                    'Sustainable harvest from a wild population; this needs the population and recruitment model.',
                    'Collision rules between aerophytes, their tethers and sky traffic.'])
