"""How people live through the Moon's 709-hour cycle, setting by setting: the physical situation.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.daily_life.run    # results/daily_life.json, about half a minute

The study computes what people would live in and leaves the ways of living to the write-up (README.md).
- Light comes from the calendar's solved clear sky (illumination/calendar, read with its own Python evaluator)
  in the joint synthesis's mean geometry: the Sun in the Moon's equatorial plane, the Earth over the sub-Earth
  point at its mean distance, no libration. The practical-dusk level is the solver's Earth control at the end of
  civil twilight (illumination/sky).
- Main's corrected per-place Earthlight (joint synthesis, commit 78020e0) is read from Git at its pinned commit
  and checked against this study's own light at the same places; it binds in-tree once main is merged.
- Air by height and its oxygen come from the design run's global means (climate/results/gcm), the flyers from the
  sky-fleet study, the aerial district from the inhabited-volume screen, the array's commute and every latency
  from the provisioning products, evening cloud from the cloud-twilight study.
- The design run's zonal-wind product (climate/results/gcm/zonal_winds_A28_dim5_moon.json) gives the wind a drifting
  town rides by height and latitude, the day it lives, how often the wind holds its hour, and the routes passive towns
  follow. Its stored climatology (ignored input, climate/gcm/products) gives the land air and cloud through the dusk.
Everything else is closed-form arithmetic on shared constants and cited values.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from illumination.calendar.geometry import geometry_from_vectors
from illumination.calendar.validate import Response
from shared.constants import (AU, EARTH_GM, EARTH_MOON_DISTANCE, EARTH_RADIUS, EARTH_SIDEREAL_DAY_S, GAS_CONSTANT,
                              JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_EQUATOR_TO_ECLIPTIC_DEG, MOON_ORBIT_ECCENTRICITY,
                              MOON_RADIUS, MOON_SURFACE_GRAVITY, SIDEREAL_MONTH_DAYS, SPEED_OF_LIGHT,
                              STANDARD_GRAVITY, SUN_GM, SUN_RADIUS, SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.daily-life/1'
OUT = HERE / 'results' / 'daily_life.json'
CODE = ('illumination/calendar/geometry.py', 'illumination/calendar/validate.py', 'illumination/calendar/disk.py')
ZONAL_WINDS = 'climate/results/gcm/zonal_winds_A28_dim5_moon.json'
INPUTS = {
    'illumination/calendar/results/transfer.json': 'terluna.illumination.calendar-transfer/1',
    'illumination/sky/results/solved_sky.json': 'terluna.illumination.solved-spherical-sky/1',
    'climate/results/gcm/global_winds_A28_dim5_moon.json': 'terluna.climate.gcm-global-winds/1',
    'research/studies/sky_ships/results/sky_ships.json': 'terluna.research.sky-ships/1',
    'research/studies/sky_fleet/results/sky_fleet.json': 'terluna.research.sky-fleet/1',
    'research/studies/summit_tower/results/summit_tower.json': 'terluna.research.summit-tower/1',
    'research/studies/inhabited_volume/results/inhabited_volume.json': 'terluna.research.inhabited-volume/1',
    'research/studies/cloud_twilight/results/regional_evenings.json': 'terluna.research.regional-cloud-evenings/1',
    'research/studies/conservation/results/heritage_register.json': 'terluna.research.conservation-register/1',
    'provisioning/results/population.json': 'terluna.provisioning.population/1',
    'provisioning/results/computing.json': 'terluna.provisioning.computing/1',
    'shared/scenarios/population.json': 'terluna.scenario.population/1',
    'geography/results/atlas.json': 'terluna.geography.atlas/1',
    'biosphere/canopy/results/canopy.json': 'terluna.biosphere.canopy-photosynthesis/1',
    ZONAL_WINDS: 'terluna.climate.gcm-zonal-winds/1',
}
# The design run's stored climatology lives on the research drive (ignored by Git); the study runs without it.
CLIMATOLOGY = 'climate/gcm/products/climatology_A28_dim5_moon.npz'
CLIMATOLOGY_SCHEMA = 'terluna.climate.gcm-climatology/1'
# Latitudes read from the zonal-wind product: model rows in bands of |latitude|, both hemispheres, area-weighted.
WIND_BANDS = {'equator': (0.0, 10.0), '30': (25.0, 35.0), '60': (55.0, 65.0)}
TOWN_HEIGHTS_KM = (2.5, 5.0, 10.0, 20.0, 30.0, 40.0)
# Main's corrected per-place Earthlight, bound by commit until the branch merges main (then read in-tree).
MAIN_PLACES = dict(path='research/studies/joint_synthesis/results/places.json',
                   commit='78020e0f3273cd98b6b7ac68554b162cf3a96239',
                   blob='3aa6d0939faf502433c03105b831a9994ab3315d',
                   sha256='9707e4d1c7332ba352b90dbf77c66fde080ec5e0fe24e57b54067c61ecdcfc8d',
                   schema='terluna.research.joint-synthesis-places/1')

CYCLE_H = SYNODIC_MONTH_DAYS * JULIAN_DAY / 3600.0       # 708.73 h from noon to noon
DEG_PER_H = 360.0 / CYCLE_H
SUN_RADIUS_DEG = math.degrees(math.asin(SUN_RADIUS / AU))
SAMPLES = 14400                                           # 0.049 h apart through one cycle
EARTH_SAMPLES = 1440                                      # the Earth's light changes slowly; interpolated between
# Light levels. Practical dusk is read from the solver's Earth control at the end of civil twilight (the Sun 6
# degrees down, USNO); the rest are conventions with their sources.
CIVIL_TWILIGHT_END_DEG = -6.0                             # USNO's civil twilight
GREY = 0.18                                               # the register's 18% grey surface (CIE 191:2010 modes)
PHOTOPIC_CD_M2, SCOTOPIC_CD_M2 = 5.0, 0.005               # CIE 191:2010 mesopic range limits
DARK_LUX = 0.1                                            # the joint synthesis's dark level (night_sky.py)
CIRCADIAN_LUX = dict(daytime_at_least=250.0, evening_at_most=10.0, sleep_at_most=1.0)   # Brown et al. 2022, melanopic EDI at the eye
# U.S. Standard Atmosphere 1976 (COESA): table 2's P0 and T0, table 4's first-layer gradient, the sea-level mean
# molecular weight from table 3's composition, and table 3's oxygen. Heights come out in the standard's geopotential
# metres (within 0.1% of geometric below 6 km); its own R* (8.31432) differs from the shared constant by 3 in 10^6.
COESA = dict(p0_pa=101325.0, t0_k=288.15, lapse_k_m=0.0065, molar_mass_kg_mol=0.0289644, o2_share=0.209476)
# Measured walk-run transitions: Kram et al. 1997 (1 g, harness simulation) and De Witt et al. 2014 (lunar gravity in
# parabolic flight); Froude numbers rise steeply below 0.4 g, so no Earth-fitted Froude rule is carried to the Moon.
GAIT = dict(earth_transition_m_s=1.98, earth_froude=0.45, lunar_transition_m_s=1.42, lunar_froude=1.39)
FRICTION = 0.6                                            # assumed tyre or shoe friction on a dry path
RESOURCE_AU = {'main belt': 2.7, 'Centaurs': 15.0, 'Kuiper belt': 40.0}   # representative distances (assumption)
RING_0_KM = 3.1                                           # the port's lowest ring above the summit (metropolis brief)
# Indoor daylight through the long night (a sensitivity): Brown et al.'s daytime level up to twice it on the lit
# floor, a utilization of the luminaires' light of 1 to 0.5, 16 waking hours a day, and the DOE's 2022 indoor
# troffer at 126 lm/W (Pattison et al. 2022, table 2.3).
INDOOR = dict(lux=(250.0, 500.0), utilization=(1.0, 0.5), waking_h_per_day=16.0, efficacy_lm_w=126.0)

SETTING_NOTES = {  # keys follow research/studies/ways_of_living; name, kind, where the light is read, what is fixed
    'summit_metropolis': ('The summit metropolis', 'place', ['summit'],
                          'Far side, 5.4 N 158.6 W, ground 11.6 km above sea level; about 100 million people at about '
                          '40,000 per km2 within 32 km of the port (register, 2026-09-27).'),
    'coastal_and_lake_towns': ('Coastal and lake towns', 'places', ['korolev_shore', 'procellarum_russell'],
                               'A lake town on Korolev\'s shore 84 km from the summit, and a nearside sea coast under a low Earth.'),
    'wetland_edge_towns': ('Wetland-edge towns', 'places', ['korolev_shore', 'far side 0', 'nearside 0'],
                           'Lake shores and wetlands (lunar-cycle ecology F3); half the land lies within 50 km of standing water.'),
    'canopy_towns': ('Canopy towns at forest edges', 'latitude rows', ['far side 0', 'far side 15', 'nearside 0', 'nearside 15'],
                     'Forest margins in the rain belt; megaforest trees 100-500 m tall (megaforest_wind); the floor\'s light '
                     'under the canopy model\'s stand is in volume.under_canopy.'),
    'farming_country': ('Rain-belt farming country', 'latitude rows', ['far side 0', 'far side 15', 'nearside 0', 'nearside 15'],
                        'Within 15 degrees of the equator, where rain-fed lakes cover 24% of the ground.'),
    'highland_towns': ('Highland and ridge towns', 'place', ['highland_box'],
                       'The climate work\'s highland box at 44.7 S on the far side, dry and sunny, comfortable from 4 km up.'),
    'elevated_districts': ('Districts held up by towers or terrain', 'heights', ['summit'],
                           'The port tower\'s rings from 3.1 km and disks from 12.3 km above the summit; the horizon and air by '
                           'height are in volume.horizon and volume.open_air.'),
    'nearside_towns': ('Nearside towns under the fixed Earth', 'places', ['nubium', 'nearside 0', 'nearside 15', 'nearside 30', 'nearside 45'],
                       'The Earth fixed in the sky, near full at local midnight.'),
    'limb_towns': ('Limb towns where the Earth rises and sets', 'places', ['smythii_headland', 'limb 0'],
                   'The libration lifts the Earth above the horizon and lowers it again once a month.'),
    'far_side_places': ('The far side\'s twilit-night places', 'places', ['ingenii_coast', 'far side 30', 'far side 45'],
                        'The far side away from its equator, where twilight wraps round the Moon.'),
    'high_latitude_towns': ('High-latitude towns', 'latitude rows', ['polar_plain', 'far side 60', 'far side 75', 'nearside 60', 'nearside 75'],
                            'Poleward of about 48.5 degrees the clear night never falls to practical dusk\'s level.'),
    'undersea_communities': ('The undersea communities', 'place', ['tranquility'],
                             'Tranquility Base under 467 m of water, an enclosed community at 1 atm (conservation D6); its '
                             'harbour and the sea above take the nearside\'s light.'),
    'roaming_sky_towns': ('Roaming sky towns', 'heights', [],
                          'Buoyant towns drifting with the wind near 10 km; the author lets much floating habitation roam '
                          '(2026-10-09). Their days are in volume.sky_towns.'),
    'moored_or_streamlined_sky_towns': ('Moored or streamlined sky towns', 'heights', [],
                                        'Towns that hold a place or follow the Sun to keep an hour; the airspeed and the '
                                        'sphere\'s power to hold the hour are in volume.hold_the_hour.'),
    'high_platforms': ('High platforms near 70 km', 'heights', [],
                       'Long-endurance platforms for far-side astronomy and nearside earthshine photometry (July canon).'),
    'array_habitats': ('The array\'s habitats', 'orbit', [],
                       'Spinning habitats near the ring radius, 20,000 km out; 990 t of shield a person (Stanford torus).'),
    'earth': ('Earth', 'reference', [], 'A 24-hour day; most of the 20-30 billion live there.'),
    'elsewhere': ('Elsewhere in the Solar System', 'reference', [], 'The unmodelled remainder; the resource operation\'s '
                  'crews are in work.resource_operation.'),
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def load_product(path, schema):
    data = json.loads(Path(path).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema for {path}')
    return data


def climatology():
    """The design run's stored climatology, or None where the research drive's file is absent."""
    path = ROOT / CLIMATOLOGY
    if not path.exists():
        return None
    data = np.load(path)
    if json.loads(str(data['metadata']))['schema'] != CLIMATOLOGY_SCHEMA:
        raise ValueError(f'Unexpected schema for {CLIMATOLOGY}')
    return data


def main_places():
    """Main's places product at its pinned commit: the in-tree file once merged, else the Git object."""
    local = ROOT / MAIN_PLACES['path']
    raw = local.read_bytes() if local.exists() else b''
    source = 'in-tree'
    if hashlib.sha256(raw).hexdigest() != MAIN_PLACES['sha256']:
        try:
            raw = subprocess.run(['git', '-C', str(ROOT), 'show', f"{MAIN_PLACES['commit']}:{MAIN_PLACES['path']}"],
                                 check=True, capture_output=True).stdout
            source = 'git object at the pinned commit'
        except (OSError, subprocess.CalledProcessError):
            return None, 'unavailable'
    if hashlib.sha256(raw).hexdigest() != MAIN_PLACES['sha256']:
        return None, 'hash mismatch'
    data = json.loads(raw)
    if data.get('schema') != MAIN_PLACES['schema']:
        raise ValueError('Unexpected schema for main places')
    return data, source


def slerp(a, b, distance_km):
    """The point distance_km from a toward b along their great circle (lat, east lon in degrees)."""
    va, vb = (unit(*p) for p in (a, b))
    angle = math.acos(max(-1.0, min(1.0, float(va @ vb))))
    f = distance_km * 1e3 / MOON_RADIUS / angle
    v = (math.sin((1 - f) * angle) * va + math.sin(f * angle) * vb) / math.sin(angle)
    return round(math.degrees(math.asin(v[2])), 3), round(math.degrees(math.atan2(v[1], v[0])) % 360, 3)


def unit(lat, lon):
    la, lo = math.radians(lat), math.radians(lon)
    return np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])


def great_circle_km(a, b):
    return math.acos(max(-1.0, min(1.0, float(unit(*a) @ unit(*b))))) * MOON_RADIUS / 1e3


class Sky:
    """Clear-sky light on level ground in the mean geometry, through the calendar's evaluator."""

    def __init__(self, transfer):
        self.response = Response(transfer)
        self.elevation = np.linspace(-90.0, 90.0, 9001)
        base = dict(sun_radius_deg=SUN_RADIUS_DEG, sun_distance_factor=1.0, earth_phase_deg=180.0,
                    source_separation_deg=180.0, earth_radius_deg=0.0)
        parts = [self.response.light({**base, 'sun_elevation_deg': float(e)}, include_earth=False) for e in self.elevation]
        self.solar_lux = np.array([p['solar'] for p in parts])
        self.direct_lux = np.array([p['solar_direct'] for p in parts])
        if np.any(np.diff(self.solar_lux) <= 0):
            raise ValueError('The Sun\'s light must rise with its elevation')

    def sun(self, elevation_deg):
        return np.interp(elevation_deg, self.elevation, self.solar_lux)

    def depression_for(self, lux):
        """The Sun's elevation (negative) at which its sky gives this light on level ground."""
        return float(np.interp(math.log(lux), np.log(self.solar_lux), self.elevation))

    def earth(self, lat, lon, hour_angle_deg):
        """The Earth's light at a place for the Sun at the given hour angles from local noon, and the Earth's elevation."""
        n = hour_angle_deg.size
        sun_lon = np.radians(lon - hour_angle_deg)            # the subsolar point moves west through the cycle
        g = dict(earth=np.tile([1.0, 0.0, 0.0], (n, 1)), earth_distance_m=np.full(n, EARTH_MOON_DISTANCE),
                 sun=np.stack([np.cos(sun_lon), np.sin(sun_lon), np.zeros(n)], axis=-1), sun_distance_m=np.full(n, AU))
        geo = geometry_from_vectors(g, lon, lat)
        light = np.array([self.response.light({k: float(v[i]) for k, v in geo.items()})['earth'] for i in range(n)])
        return light, float(geo['earth_elevation_deg'][0])


def crossing_hours(t, lux, level):
    """Hours below a level and the times the light crosses it, by log-linear interpolation between samples."""
    y, line = np.log(np.maximum(lux, 1e-30)), math.log(level)
    y0, y1 = y, np.roll(y, -1)
    dt = t[1] - t[0]
    below = np.where((y0 < line) & (y1 < line), dt, 0.0)
    fall = (y0 >= line) & (y1 < line)
    rise = (y0 < line) & (y1 >= line)
    frac = np.where(fall | rise, (y0 - line) / np.where(y0 != y1, y0 - y1, 1.0), 0.0)
    below = below + np.where(fall, (1 - frac) * dt, 0.0) + np.where(rise, frac * dt, 0.0)
    return float(below.sum()), t[fall] + frac[fall] * dt, t[rise] + frac[rise] * dt


def calendar(sky, lat, lon, levels):
    """One place's light through the cycle: the four states and the hours at each reference level."""
    hour_angle = np.arange(SAMPLES) * 360.0 / SAMPLES       # the Sun's hour angle from local noon
    t = hour_angle / DEG_PER_H
    elevation = np.degrees(np.arcsin(np.cos(math.radians(lat)) * np.cos(np.radians(hour_angle))))
    sun = sky.sun(elevation)
    coarse = np.arange(EARTH_SAMPLES) * 360.0 / EARTH_SAMPLES
    earth_coarse, earth_elevation = sky.earth(lat, lon, coarse)
    earth = np.interp(hour_angle, np.append(coarse, 360.0), np.append(earth_coarse, earth_coarse[0]))
    total = sun + earth
    sunset, sunrise = 90.0 / DEG_PER_H, 270.0 / DEG_PER_H
    night = (t >= sunset) & (t < sunrise)
    out = dict(lat=lat, lon=lon, earth_elevation_deg=round(earth_elevation, 1),
               side='nearside' if earth_elevation > 8 else ('far side' if earth_elevation < -8 else 'limb'))
    # Dusk and dawn are the Sun's twilight: from sunset until its sky falls to the level, and the mirror before
    # sunrise. The Earth's light is then counted by the night hours it holds above the level.
    _, falls, rises = crossing_hours(t, sun, levels['dusk'])
    first = [x for x in falls if sunset <= x < sunrise]
    last = [x for x in rises if sunset <= x < sunrise]
    sun_down_h = sunrise - sunset
    if first:
        dusk_h, dawn_h = first[0] - sunset, sunrise - last[-1]
        night_h = last[-1] - first[0]
    else:
        dusk_h = dawn_h = sun_down_h / 2
        night_h = 0.0
    below_dusk = crossing_hours(t, total, levels['dusk'])[0]
    out.update(day_h=round(CYCLE_H - sun_down_h, 2), dusk_h=round(dusk_h, 2), night_h=round(night_h, 2),
               dawn_h=round(dawn_h, 2), dusk_runs_into_dawn=not first,
               below_dusk_h=round(below_dusk, 2), earth_holds_above_dusk_h=round(night_h - below_dusk, 2))
    out['below_dark_h'] = round(crossing_hours(t, total, levels['dark'])[0], 2)
    out['below_scotopic_h'] = round(crossing_hours(t, total, levels['scotopic'])[0], 2)
    day_vision = [x for x in crossing_hours(t, total, levels['photopic'])[1] if sunset <= x < sunrise]
    out['day_vision_after_sunset_h'] = round(day_vision[0] - sunset, 2) if day_vision else None
    mid = SAMPLES // 2
    out.update(midnight_lux=round(float(total[mid]), 4), darkest_lux=round(float(total[night].min()), 4),
               earth_midnight_lux=round(float(earth[mid]), 4),
               earth_outshines_sky_h=round(float(((earth > sun) & night).sum() * (t[1] - t[0])), 1),
               noon_lux=round(float(total[0]), 0))
    circ = {}
    for key, level in CIRCADIAN_LUX.items():
        hours_below = crossing_hours(t, total, level)[0]
        night_below = crossing_hours(t, np.where(night, total, 1e9), level)[0]
        circ[key] = dict(level_lux=level, night_hours_below=round(night_below, 1),
                         night_hours_at_or_above=round(sun_down_h - night_below, 1), cycle_hours_below=round(hours_below, 1))
    out['circadian_reference'] = circ
    out['sleeps'] = {k: round(out[f'{k}_h'] / 24.0, 1) for k in ('day', 'dusk', 'night', 'dawn')}
    return out


def twilight_curve(sky, lat, levels):
    """Hours after sunset at which the Sun's sky falls to each level, and its light at set hours."""
    out = dict(lat=lat)
    reach = {}
    for name, lux in levels.items():
        e = sky.depression_for(lux)
        c = math.sin(math.radians(e)) / math.cos(math.radians(lat))
        reach[name] = dict(lux=round(lux, 3), sun_deg=round(e, 2),
                           hours_after_sunset=round((math.degrees(math.acos(c)) - 90.0) / DEG_PER_H, 1) if c >= -1 else None)
    out['reaches'] = reach
    hours = np.array([0, 6, 12, 24, 36, 48, 60, 72, 84, 96, 120, 150, 90.0 / DEG_PER_H])
    e = np.degrees(np.arcsin(-math.cos(math.radians(lat)) * np.sin(np.radians(hours * DEG_PER_H))))
    out['light'] = [dict(hours_after_sunset=round(float(h), 1), sun_deg=round(float(x), 2), lux=float(f'{sky.sun(x):.4g}'))
                    for h, x in zip(hours, e)]
    return out


def seasonal_day(lat):
    """Day length through the Moon's year: the Sun's declination swings by the equator's tilt to the ecliptic."""
    if abs(lat) >= 90 - MOON_EQUATOR_TO_ECLIPTIC_DEG:
        return None
    x = math.tan(math.radians(abs(lat))) * math.tan(math.radians(MOON_EQUATOR_TO_ECLIPTIC_DEG))
    longest = math.degrees(math.acos(-x)) / 180.0 * CYCLE_H
    return dict(longest_day_h=round(longest, 1), shortest_day_h=round(CYCLE_H - longest, 1))


class Air:
    """The design run's global-mean air by height above its sea level, with the oxygen it holds."""

    def __init__(self, winds):
        self.z = np.array(winds['height_km'])
        self.lnp = np.log(np.array(winds['pressure_pa']))
        self.t = np.array(winds['temperature_k'])
        gases = winds['configuration']['gases_bar']
        self.o2_share = gases['pO2'] / sum(gases.values())
        self.surface_pa = winds['configuration']['pressure_pa']

    def at(self, z_km):
        p = float(np.exp(np.interp(z_km, self.z, self.lnp)))
        t = float(np.interp(z_km, self.z, self.t))
        po2 = p * self.o2_share
        return dict(height_km=round(z_km, 3), pressure_atm=round(p / COESA['p0_pa'], 3),
                    temperature_c=round(t - 273.15, 1), o2_kpa=round(po2 / 1e3, 2),
                    earth_equivalent_m=round(earth_equivalent_m(po2), -1))

    def height_of(self, p_pa):
        return float(np.interp(-math.log(p_pa), -self.lnp, self.z))


def earth_equivalent_m(o2_pa):
    """The height in the 1976 standard atmosphere whose dry air holds the same oxygen partial pressure."""
    exponent = STANDARD_GRAVITY * COESA['molar_mass_kg_mol'] / (GAS_CONSTANT * COESA['lapse_k_m'])
    p = o2_pa / COESA['o2_share']
    return COESA['t0_k'] / COESA['lapse_k_m'] * (1.0 - (p / COESA['p0_pa']) ** (1.0 / exponent))


def horizon_dip_deg(h_km):
    return math.degrees(math.acos(MOON_RADIUS / (MOON_RADIUS + h_km * 1e3)))


def sunlit_hours_at_height(h_km, lat):
    """Hours of direct sun at a height above level ground, geometric horizon without refraction."""
    c = -math.sin(math.radians(horizon_dip_deg(h_km))) / math.cos(math.radians(lat))
    return math.degrees(math.acos(max(-1.0, c))) / 180.0 * CYCLE_H


def sun_ground_speed(h_km, lat):
    """The westward speed of the subsolar point's meridian at a height and latitude, m/s."""
    return 2 * math.pi * (MOON_RADIUS + h_km * 1e3) * math.cos(math.radians(lat)) / (CYCLE_H * 3600.0)


def canopy_floor(canopy):
    """Share of the photosynthetic light reaching the floor under the canopy model's stand, by Sun height.

    The stand's energy balance: incident = absorbed by leaves + reflected out + absorbed by the soil, and the soil
    absorbs (1 - its albedo) of what reaches it. Photosynthetic photons stand in for illuminance.
    """
    soil = canopy['settings']['stand']['soil']
    return {r['sun_deg']: (r['incident_par_umol_m2_s'] * (1 - r['canopy_albedo_par']) - r['absorbed_umol_m2_s'])
            / (1 - soil) / r['incident_par_umol_m2_s'] for r in canopy['by_sun_height']['moon']}


def town_day_hours(h_km, lat, u, v=None):
    """A drifting town's solar day: T x v_sun / |v_sun + u|, u eastward; None where it holds the hour."""
    v = sun_ground_speed(h_km, lat) if v is None else v
    s = v + u
    if abs(s) < 1e-9:
        return None, False
    return CYCLE_H * v / abs(s), s < 0


class Winds:
    """The zonal-wind product's arrays at its heights, read in bands of |latitude| over both hemispheres."""

    NAMES = dict(u=('ground_frame', 'u_mean_m_s'), coverage=('ground_frame', 'coverage'), v_sun=(None, 'sun_speed_m_s'),
                 day_p10=('solar_day', 'p10_days'), day_p50=('solar_day', 'p50_days'), day_p90=('solar_day', 'p90_days'),
                 hold_1=('holding', 'hold_1'), reversed=('holding', 'reversed'))

    def __init__(self, zonal):
        self.z = np.array(zonal['axes']['height_km'], dtype=float)
        self.lat = np.array(zonal['axes']['lat_deg'], dtype=float)
        self.arrays = {k: np.array(zonal[a][b] if a else zonal[b], dtype=float) for k, (a, b) in self.NAMES.items()}

    def row(self, h_km):
        i = np.flatnonzero(np.isclose(self.z, h_km))
        if i.size != 1:
            raise ValueError(f'No zonal-wind level at {h_km} km')
        return int(i[0])

    def band(self, name, h_km, band):
        """Area-weighted mean over the model rows whose |latitude| lies in the band (null rows left out)."""
        x = self.arrays[name][self.row(h_km)]
        keep = (np.abs(self.lat) >= band[0]) & (np.abs(self.lat) <= band[1]) & np.isfinite(x)
        w = np.cos(np.radians(self.lat[keep]))
        return float((x[keep] * w).sum() / w.sum())

    def global_u(self, h_km):
        """The area-weighted eastward wind over the ground lying below the height."""
        i = self.row(h_km)
        u, w = self.arrays['u'][i], np.cos(np.radians(self.lat)) * np.nan_to_num(self.arrays['coverage'][i])
        keep = np.isfinite(u) & (w > 0)
        return float((u[keep] * w[keep]).sum() / w[keep].sum())


QUARTERS = {'after_midnight': (-180.0, -90.0), 'morning': (-90.0, 0.0), 'afternoon': (0.0, 90.0),
            'evening_and_first_half_of_night': (90.0, 180.0)}


def routes(zonal):
    """Passive routes at fixed heights through the 3-day means: their day, where in it they linger, where they go."""
    t = zonal['trajectories']
    long_run, first_year = t['sequence']['long_run'], t['sequence']['first_year']
    hour_edges, lat_edges = np.array(t['occupancy_hour_edges_deg']), np.array(t['occupancy_lat_edges_deg'])
    hour_c, lat_c = (hour_edges[:-1] + hour_edges[1:]) / 2, (lat_edges[:-1] + lat_edges[1:]) / 2
    w = np.cos(np.radians(np.array(t['release_lat_deg'], dtype=float)))
    days = long_run['gathering']['days']
    out = []
    for i, h in enumerate(t['heights_km']):
        mean_day = float(np.nanmedian(np.array(long_run['mean_day_days'][i], dtype=float)))
        occupancy = np.array(long_run['occupancy'][i], dtype=float)          # years 5-10, [lat band][hour sector]
        by_hour, by_lat = occupancy.sum(0) / occupancy.sum(), occupancy.sum(1) / occupancy.sum()
        linger = float(hour_c[int(np.argmax(by_hour))])
        western = np.array(first_year['western_sunrise_share'][i], dtype=float)
        out.append(dict(
            height_km=h, mean_day_days=round(mean_day, 1),
            first_pass_days=round(float(np.nanmedian(np.array(long_run['first_pass_days'][i], dtype=float))), 1),
            days_by_quarter={k: round(float(by_hour[(hour_c > lo) & (hour_c < hi)].sum()) * mean_day, 1)
                             for k, (lo, hi) in QUARTERS.items()},
            lingers_at_hour_deg=linger, lingers_at_clock=f'{int((12 + linger / 15) % 24):02d}:00',
            linger_strength=round(float(by_hour.max() * by_hour.size), 2),
            sun_up_share=round(float(np.nanmedian(np.array(long_run['sun_up_share'][i], dtype=float))), 3),
            light_ratio=round(float(np.nanmedian(np.array(long_run['light_ratio'][i], dtype=float))), 3),
            years_5_10_within_20_deg_of_equator=round(float(by_lat[np.abs(lat_c) < 20].sum()), 3),
            years_5_10_between_20_and_70_deg=round(float(by_lat[(np.abs(lat_c) > 20) & (np.abs(lat_c) < 70)].sum()), 3),
            years_5_10_poleward_of_50_deg=round(float(by_lat[np.abs(lat_c) > 50].sum()), 3),
            years_5_10_poleward_of_70_deg=round(float(by_lat[np.abs(lat_c) > 70].sum()), 3),
            nearest_neighbour_km={f'{d} d': long_run['gathering']['nearest_km_p50'][k][i]
                                  for k, d in enumerate(days) if d in (30, 357, 3575)},
            first_year_western_sunrise_share=round(float((western * w).sum() / w.sum()), 3)))
    return out


def kepler_hours(r_start, r_peri, r_apo, mu):
    """Time from r_start inbound to periapsis on a two-body conic (r_apo=inf for a parabola), hours."""
    if math.isinf(r_apo):
        cos_nu = 2 * r_peri / r_start - 1
        d = math.tan(math.acos(cos_nu) / 2)
        return math.sqrt(2 * r_peri**3 / mu) * (d + d**3 / 3) / 3600.0
    a = (r_peri + r_apo) / 2
    e = (r_apo - r_peri) / (r_apo + r_peri)
    big_e = math.acos(max(-1.0, min(1.0, (1 - r_start / a) / e)))
    return (big_e - e * math.sin(big_e)) / math.sqrt(mu / a**3) / 3600.0


def fall(h_m, terminal):
    """Impact speed from a fall of h_m with quadratic drag, and the Earth height that gives it without drag."""
    v = terminal * math.sqrt(1 - math.exp(-2 * MOON_SURFACE_GRAVITY * h_m / terminal**2))
    return dict(height_m=h_m, impact_m_s=round(v, 2), like_earth_fall_m=round(v**2 / (2 * STANDARD_GRAVITY), 2),
                seconds_without_drag=round(math.sqrt(2 * h_m / MOON_SURFACE_GRAVITY), 2))


def centrifuge_radius_m(rpm, g_target, floor_gravity):
    """Radius at which spin and the floor's own gravity give a resultant of g_target (0 for orbit)."""
    omega = rpm * 2 * math.pi / 60
    return math.sqrt(g_target**2 - floor_gravity**2) / omega**2


def results():
    data = {p: load_product(ROOT / p, s) for p, s in INPUTS.items()}
    (transfer, solved, global_winds, ships, fleet, tower, volume, evenings, register, pop, comp, scenario, atlas, canopy,
     zonal) = data.values()
    clim = climatology()
    main, main_source = main_places()
    sky, air = Sky(transfer), Air(global_winds)
    earth_control = next(r for r in solved['worlds']['earth_control']['samples'] if r['sun_deg'] == CIVIL_TWILIGHT_END_DEG)
    levels = dict(dusk=earth_control['total_horizontal_lux'], dark=DARK_LUX,
                  photopic=PHOTOPIC_CD_M2 * math.pi / GREY, scotopic=SCOTOPIC_CD_M2 * math.pi / GREY)
    site = tower['site']
    korolev = next(l for l in site['lakes'] if l['area_km2'] > 100000)
    tranquility = next(s for s in register['sites'] if s['id'] == 'apollo-11')
    places = {
        'summit': (5.4125, 201.3665),
        'korolev_shore': slerp((5.4125, 201.3665), (korolev['centre'][0], korolev['centre'][1] % 360), korolev['shore_km']),
        'procellarum_russell': (32.9, 284.6),
        'smythii_headland': (-1.5, 90.0),
        'nubium': (-25.0, 345.0),
        'ingenii_coast': (-34.0, 163.0),
        'highland_box': (-44.7, 246.1),
        'polar_plain': (75.0, 0.0),
        'tranquility': (tranquility['lat'], tranquility['lon'] % 360),
    }
    rows = {f'{side} {lat:g}': (lat, lon) for side, lon in (('far side', 180.0), ('nearside', 0.0))
            for lat in (0.0, 15.0, 30.0, 45.0, 60.0, 75.0)}
    rows['limb 0'] = (0.0, 90.0)
    calendars = {name: calendar(sky, lat, lon, levels) for name, (lat, lon) in {**places, **rows}.items()}

    # Main's corrected Earthlight at its seven places against this study's own light there.
    names = {'summit port': 'summit', 'Oceanus Procellarum by Russell': 'procellarum_russell',
             'eastern Smythii headland': 'smythii_headland', 'southern Mare Nubium': 'nubium',
             'South Pole-Aitken coast by Mare Ingenii': 'ingenii_coast', 'highland box': 'highland_box',
             'polar plain': 'polar_plain'}
    check = dict(path=MAIN_PLACES['path'], commit=MAIN_PLACES['commit'], blob=MAIN_PLACES['blob'],
                 sha256=MAIN_PLACES['sha256'], read_from=main_source, places={})
    if main:
        for theirs, ours in names.items():
            m, c = main['places'][theirs]['night'], calendars[ours]
            check['places'][ours] = dict(main_midnight_lux=m['natural_midnight_lux'], study_midnight_lux=c['midnight_lux'],
                                         main_darkest_lux=m['natural_darkest_lux'], study_darkest_lux=c['darkest_lux'],
                                         main_earth_elevation_deg=m['earth_elevation_deg'],
                                         study_earth_elevation_deg=c['earth_elevation_deg'],
                                         main_hours_below_dark=m['hours_below_0_1_lux_natural'],
                                         study_hours_below_dark=c['below_dark_h'])

    latitudes = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0)
    curve_levels = {'bright, 1,000 lux': 1000.0, 'daytime circadian reference, 250 lux': 250.0,
                    'day vision, 87 lux': levels['photopic'], 'evening circadian reference, 10 lux': 10.0,
                    'practical dusk': levels['dusk'], 'sleep reference, 1 lux': 1.0, 'dark, 0.1 lux': DARK_LUX,
                    'night vision, 0.087 lux': levels['scotopic']}
    dusk_deg = sky.depression_for(levels['dusk'])
    terminator = []
    for la in latitudes:
        dusk_hours = twilight_curve(sky, la, {'d': levels['dusk']})['reaches']['d']['hours_after_sunset']
        width = None if dusk_hours is None else round(
            math.radians(dusk_hours * DEG_PER_H) * MOON_RADIUS * math.cos(math.radians(la)) / 1e3)
        terminator.append(dict(lat=la, km_h=round(sun_ground_speed(0, la) * 3.6, 1), practical_dusk_band_km=width))

    # Weather through the dusk: evening cloud on the cloud-resolving rings, rain on the equatorial ring,
    # and the design run's land air and cloud cover by hour angle where its climatology is present.
    soaring = fleet['people_on_wings']['soaring']
    ha = np.array(soaring['hour_angle_deg'])
    dusk_bins = (ha > 90) & (ha < 90 + abs(dusk_deg))
    weather = dict(
        evening_cloudy_columns={c['case']: dict(by_latitude=[dict(abs_lat=b['abs_latitude_deg'], share=round(b['cloud_fraction'], 4)) for b in c['latitude_bands']],
                                                by_hours_after_sunset=[dict(hours=b['hours_after_sunset'], share=round(b['cloud_fraction'], 4)) for b in c['time_bins']])
                                for c in evenings['case_summaries']},
        ring_rain_mm_h_through_dusk=[dict(hours_after_sunset=round((h - 90) / DEG_PER_H, 1), rain_mm_h=r)
                                     for h, r in zip(ha[dusk_bins], np.array(soaring['rain_mm_h'])[dusk_bins])],
        ring_rain_spell_hours=soaring['rain_spell_hours'], ring_land_wet_share=soaring['land_wet_share'],
        night_mixed_layer_km=soaring['night_mixed_layer_km'])
    if clim is not None:
        lat_axis, hours_axis = np.array(clim['lat']), np.array(clim['hour_angle_deg'])
        dusk_cols = (hours_axis > 90) & (hours_axis < 90 + abs(dusk_deg))
        day_cols = np.abs(hours_axis) < 90
        night_cols = np.abs(hours_axis) >= 90 + abs(dusk_deg)
        weather['gcm_by_latitude'] = []
        for la in latitudes:
            rows_ = np.abs(np.abs(lat_axis) - la) <= 3.0
            tas, clt = clim['tas_land_sun'][rows_], clim['clt_sun'][rows_]
            weather['gcm_by_latitude'].append(dict(
                lat=la, land_air_c=dict(day=round(float(np.nanmean(tas[:, day_cols])) - 273.15, 1),
                                        dusk=round(float(np.nanmean(tas[:, dusk_cols])) - 273.15, 1),
                                        night=round(float(np.nanmean(tas[:, night_cols])) - 273.15, 1)),
                cloud_cover=dict(day=round(float(clt[:, day_cols].mean()), 2), dusk=round(float(clt[:, dusk_cols].mean()), 2),
                                 night=round(float(clt[:, night_cols].mean()), 2))))

    ferries = next(f for f in fleet['sky_ferries'] if f['length_m'] == 250.0)
    line = ferries['line_per_hour_per_direction']
    gathering = []
    for la in (0.0, 30.0):
        d = twilight_curve(sky, la, {'d': levels['dusk']})['reaches']['d']['hours_after_sunset']
        for n in (1e4, 1e5, 1e6):
            gathering.append(dict(lat=la, practical_dusk_h=d, attendance=int(n), stay_h=4.0,
                                  arrivals_per_hour=round(n / d), present_at_once=round(n * 4.0 / d),
                                  share_of_one_ferry_line=[round(n / d / line['3 min'], 3), round(n / d / line['2 min'], 3)]))

    # The inhabited volume: air by height, the summit's open-air zones, the horizon from above, drifting days.
    ground = site['ground_above_sea_m'] / 1e3
    zones = [('summit ground, the metropolis', 0.0), ('top of the parks for everyone', 3.0),
             ('top of the parks for residents and acclimatised visitors', 10.0), ('winter parks begin', 13.5),
             ('top of the open-air high gardens; sealed above', 15.0), ('the crown', 24.0)]
    open_air = [dict(label=label, above_summit_km=k, **air.at(ground + k)) for label, k in zones]
    heights = [dict(label=label, **air.at(z)) for label, z in
               (('model sea level', 0.0), ('a lake town on Korolev\'s shore', korolev['level_above_sea_m'] / 1e3),
                ('a roaming sky town', 10.0), ('regional ships', 25.0), ('the flight band', 40.0),
                ('winged liners', 55.0), ('high platforms', 70.0))]
    horizon = [dict(above_ground_km=h, dip_deg=round(horizon_dip_deg(h), 2),
                    extra_sun_h={f'{la:g}': round(sunlit_hours_at_height(h, la) - CYCLE_H / 2, 1) for la in (0.0, 30.0, 60.0)})
               for h in (0.1, 0.5, 1.0, 3.1, 10.0, 12.3, 24.0, 70.0)]
    # A drifting town's day, read from the zonal-wind product by height and latitude, beside the day at its mean wind.
    winds = Winds(zonal)
    towns = []
    for h in TOWN_HEIGHTS_KM:
        for name, band in WIND_BANDS.items():
            v, u = winds.band('v_sun', h, band), winds.band('u', h, band)
            day, reversed_ = town_day_hours(h, None, u, v)
            towns.append(dict(height_km=h, latitudes=name, abs_lat_deg=list(band), sun_speed_m_s=round(v, 2),
                              u_mean_m_s=round(u, 2), day_at_mean_wind_days=round(day / 24, 1),
                              sun_turned_back_at_mean_wind=reversed_,
                              day_days={q: round(winds.band(f'day_{q}', h, band), 1) for q in ('p10', 'p50', 'p90')},
                              held_within_1_m_s=round(winds.band('hold_1', h, band), 3),
                              sun_turned_back=round(winds.band('reversed', h, band), 3)))
    drift_formula = []
    for la in (0.0, 60.0):
        v = sun_ground_speed(10.0, la)
        for u in [-15.0, -10.0, -6.0, -2.0, 0.0, 2.0, 4.0, 8.0, 15.0, -v, -2 * v, -0.5 * v]:
            day, reversed_ = town_day_hours(10.0, la, u)
            drift_formula.append(dict(height_km=10.0, lat=la, u_m_s=round(u, 3), day_h=None if day is None else round(day, 1),
                                      holds_the_hour=day is None, sun_rises_in_west=reversed_))
    wind_holds = [dict(height_km=b['height_km'], abs_lat_deg=b['abs_lat_deg'], held_within_1_m_s=b['hold_1'],
                       best_hour_deg=b['hold_1_peak_hour_deg'], held_at_best_hour=b['hold_1_peak_share'],
                       sun_turned_back=b['reversed'])
                  for b in zonal['bands'] if b['height_km'][1] <= 15.0]
    paths = routes(zonal)
    held = next(r for r in volume['held_city_drag'] if r['relative_air_speed_m_s'] == 5.0)
    reference = volume['cases']['total25']['aerial_reference']
    kw_per_resident_5 = held['input_power_w'] / reference['residents'] / 1e3
    hold_hour = []
    for name, band in WIND_BANDS.items():
        rel = abs(-winds.band('v_sun', 10.0, band) - winds.band('u', 10.0, band))
        hold_hour.append(dict(latitudes=name, airspeed_m_s=round(rel, 2),
                              sphere_kw_per_resident=round(kw_per_resident_5 * (rel / 5.0)**3, 1)))
    solar_noon = {e: dict(total_lux=round(float(sky.sun(e))), diffuse_lux=round(float(sky.sun(e) - np.interp(e, sky.elevation, sky.direct_lux))))
                  for e in (90.0, 45.0, 20.0)}
    omega = math.radians(DEG_PER_H)
    shadow = dict(envelope_diameter_m=2 * reference['radius_m'], height_km=reference['height_km'],
                  ground_light_under_it=solar_noon,
                  held_hours_in_shadow_at_noon_equator=round(2 * reference['radius_m'] / (reference['height_km'] * 1e3 * omega), 1),
                  roaming_minutes_in_shadow={f'{u:g} m/s': round(2 * reference['radius_m'] / u / 60, 1) for u in (1.0, 1.75, 5.0)},
                  sky_blocked_share=round(reference['projected_area_km2'] * 1e6 / (math.pi * (reference['height_km'] * 1e3)**2), 4))
    floor = canopy_floor(canopy)
    diffuse_share = floor[min(floor)]                     # the lowest Sun's sky is nearly all diffuse, as in twilight
    under_canopy = dict(lai=canopy['settings']['stand']['lai'],
                        floor_share_by_sun_deg={f'{k:g}': round(v, 4) for k, v in sorted(floor.items(), reverse=True)},
                        noon_floor_lux=round(float(sky.sun(90.0)) * floor[max(floor)]),
                        twilight_floor_share=round(diffuse_share, 4))
    for name, level in (('practical_dusk', levels['dusk']), ('day_vision', levels['photopic'])):
        e = sky.depression_for(level / diffuse_share)
        under_canopy[f'{name}_after_sunset_h'] = {f'{la:g}': round((math.degrees(math.acos(math.sin(math.radians(e)) / math.cos(math.radians(la)))) - 90) / DEG_PER_H, 1)
                                                  for la in (0.0, 15.0)}
    boat = next(b for b in fleet['electric']['moon'] if b['seats'] == 4)
    wing = fleet['people_on_wings']['lunar_designs']['personal_wing']['moon_summit']
    glider = fleet['people_on_wings']['lunar_designs']['micro_glider']['moon_summit']
    efficiency = fleet['electric']['technology']['efficiency']
    vertical = dict(
        lift_kwh_per_100kg_km=round(100 * MOON_SURFACE_GRAVITY * 1e3 / 3.6e6, 4),
        lift_kwh_per_100kg_km_earth=round(100 * STANDARD_GRAVITY * 1e3 / 3.6e6, 4),
        stairs_earth_metres_per_100_m=round(100 * MOON_SURFACE_GRAVITY / STANDARD_GRAVITY, 1),
        sky_boat_climb_10km_kwh=round(boat['mass_kg'] * MOON_SURFACE_GRAVITY * 1e4 / efficiency / 3.6e6, 2),
        sky_boat_trip_10km_kwh=boat['city_trip_kwh'],
        pedalled_wing_climb_m_s={k: v['climb_m_s'] for k, v in wing.items() if isinstance(v, dict)},
        glider_range_from_10km_km=round(10 * glider['best_glide'], 0),
        glider_hours_down_from_10km=round(1e4 / glider['least_sink_m_s'] / 3600, 1),
        canopy_descent_minutes=fleet['people_on_wings']['canopy']['descent_minutes'],
        glide_from_port_ring_to_korolev=dict(
            ring_above_summit_km=RING_0_KM,
            height_km=round(RING_0_KM + korolev['summit_above_lake_m'] / 1e3, 2),
            reach_km=round((RING_0_KM + korolev['summit_above_lake_m'] / 1e3) * glider['best_glide'], 0),
            shore_km=korolev['shore_km'], hours=round(korolev['shore_km'] * 1e3 / glider['best_glide_speed_m_s'] / 3600, 1)))

    # Bodies and generations: gravity, falls, gait, traction, spin.
    terminal = fleet['people_on_wings']['canopy']['free_fall_m_s'][0]
    gait = dict(GAIT, minutes_per_km_at_lunar_transition=round(1e3 / GAIT['lunar_transition_m_s'] / 60, 1),
                minutes_per_km_at_earth_transition=round(1e3 / GAIT['earth_transition_m_s'] / 60, 1))
    bodies = dict(
        gravity_m_s2=round(MOON_SURFACE_GRAVITY, 4), gravity_g=round(MOON_SURFACE_GRAVITY / STANDARD_GRAVITY, 4),
        falls=[fall(h, terminal) for h in (0.5, 1.0, 3.0, 10.0, 30.0, 100.0)],
        terminal_m_s=terminal, terminal_like_earth_fall_m=round(terminal**2 / (2 * STANDARD_GRAVITY), 1),
        jump_height_factor=round(STANDARD_GRAVITY / MOON_SURFACE_GRAVITY, 2), gait=gait,
        traction=dict(friction=FRICTION, braking_m_s2=round(FRICTION * MOON_SURFACE_GRAVITY, 2),
                      stop_from_5_m_s_m=dict(moon=round(25 / (2 * FRICTION * MOON_SURFACE_GRAVITY), 1),
                                              earth=round(25 / (2 * FRICTION * STANDARD_GRAVITY), 1))),
        centrifuge=[dict(rpm=rpm, lunar_surface_1g_radius_m=round(centrifuge_radius_m(rpm, STANDARD_GRAVITY, MOON_SURFACE_GRAVITY), 1),
                         floor_bank_deg=round(math.degrees(math.atan(math.sqrt(STANDARD_GRAVITY**2 - MOON_SURFACE_GRAVITY**2) / MOON_SURFACE_GRAVITY)), 1),
                         orbit_1g_radius_m=round(centrifuge_radius_m(rpm, STANDARD_GRAVITY, 0.0), 1),
                         orbit_lunar_g_radius_m=round(centrifuge_radius_m(rpm, MOON_SURFACE_GRAVITY, 0.0), 1))
                    for rpm in (1, 2, 3, 4)],
        array_spin=pop['array']['spin'])

    # Work, the array, the resource operation, Earth.
    commute = pop['commute']
    light_s = lambda m: m / SPEED_OF_LIGHT
    a, e = EARTH_MOON_DISTANCE, MOON_ORBIT_ECCENTRICITY
    earth_turn_h = 1 / (1 / EARTH_SIDEREAL_DAY_S - 1 / (SIDEREAL_MONTH_DAYS * JULIAN_DAY)) / 3600
    transfers = {name: round(kepler_hours(EARTH_MOON_DISTANCE, EARTH_RADIUS, apo, EARTH_GM) / 24, 2)
                 for name, apo in (('apogee at the Moon (least energy)', EARTH_MOON_DISTANCE),
                                   ('apogee at twice the Moon\'s distance', 2 * EARTH_MOON_DISTANCE),
                                   ('parabolic', math.inf))}
    resource = {name: dict(au=r, light_one_way_min=[round(light_s((r - 1) * AU) / 60, 1), round(light_s((r + 1) * AU) / 60, 1)],
                           least_energy_transfer_years=round(math.pi * math.sqrt(((1 + r) / 2 * AU)**3 / SUN_GM) / (JULIAN_YEAR_DAYS * JULIAN_DAY), 2))
                for name, r in RESOURCE_AU.items()}
    floor_m2 = pop['affluence']['living']['decent_floor_m2']
    night_days = CYCLE_H / 2 / 24.0
    lighting = [round(lux * floor_m2 / (INDOOR['efficacy_lm_w'] * u) * INDOOR['waking_h_per_day'] * night_days / 1e3, 1)
                for lux, u in zip(INDOOR['lux'], INDOOR['utilization'])]
    europe_night = pop['night']['final_kwh_per_person']['european']
    work = dict(
        indoor_daylight_through_the_night=dict(assumptions=INDOOR, floor_m2=floor_m2, kwh_per_person=lighting,
                                               share_of_european_night=[round(x / europe_night, 3) for x in lighting]),
        primary_kw_band=pop['affluence']['band_primary_kw'], night_final_kwh_per_person=pop['night']['final_kwh_per_person'],
        travel_hours_per_day=pop['affluence']['living']['travel_hours_per_day'],
        array=dict(transfer_hours={k: v['hours'] for k, v in commute['transfers'].items()},
                   commuter_kw=commute['transfers']['hohmann']['continuous_kw_per_commuter'],
                   shield_t_per_person=pop['array']['shield_t_per_person'], spin=pop['array']['spin'],
                   round_trip_light_s=comp['latency']['moon_array_round_trip_s']),
        resource_operation=resource)
    earth = dict(
        light_one_way_s=dict(mean=comp['latency']['earth_moon_one_way_s'],
                             perigee=round(light_s(a * (1 - e)), 3), apogee=round(light_s(a * (1 + e)), 3)),
        conversation_gap_s=comp['latency']['conversation_gap_s'],
        round_trip_over_gap=round(comp['latency']['earth_moon_round_trip_s'] / comp['latency']['conversation_gap_s'], 1),
        across_the_moon_light_ms=round(light_s(math.pi * MOON_RADIUS) * 1e3, 1),
        transfer_days=transfers,
        earth_diameter_deg=dict(mean=round(2 * math.degrees(math.asin(EARTH_RADIUS / a)), 2),
                                perigee=round(2 * math.degrees(math.asin(EARTH_RADIUS / (a * (1 - e)))), 2),
                                apogee=round(2 * math.degrees(math.asin(EARTH_RADIUS / (a * (1 + e)))), 2)),
        earth_turns_once_h=round(earth_turn_h, 2),
        earth_elevation_by_place={k: calendars[k]['earth_elevation_deg'] for k in places})

    # Shared air and water.
    lunar = scenario['cases']['total25']
    people = (lunar['lunar_surface'] + lunar['lunar_aerial']) * 1e9
    air_kg = air.surface_pa * 4 * math.pi * MOON_RADIUS**2 / MOON_SURFACE_GRAVITY
    shared = dict(air_kg=float(f'{air_kg:.3g}'), people_case='total25 lunar', people=people,
                  air_t_per_person=round(air_kg / people / 1e3, -3),
                  water_kg=atlas['water_mass_kg'], water_m3_per_person=round(atlas['water_mass_kg'] / 1e3 / people, -3),
                  runoff_m3_per_person_year=round(pop['water']['moon_runoff_km3'] * 1e9 / people, 0),
                  heat_w_m2_per_tw=pop['heat']['moon_w_m2_per_tw'], tw_per_kelvin=pop['heat']['moon_tw_per_k'])

    # Travel between the places, and how far apart their phases of the cycle lie.
    liner = fleet['long_haul']['fast_winged']
    long_haul = next(o for o in fleet['long_haul']['options'] if o['length_m'] == 750.0)
    regional = next(s for s in fleet['regional']['ships'] if s['length_m'] == 350.0)
    pairs = []
    keys = list(places)
    for i, p in enumerate(keys):
        for q in keys[i + 1:]:
            km = great_circle_km(places[p], places[q])
            dlon = (places[q][1] - places[p][1] + 180) % 360 - 180
            row = dict(a=p, b=q, km=round(km), phase_offset_h=round(dlon / DEG_PER_H, 1),
                       winged_liner_h=round(km * 1e3 / liner['cruise_m_s'] / 3600, 1),
                       winged_liner_kwh=round(km * liner['wh_per_passenger_km'] / 1e3, 1),
                       long_haul_h=round(km * 1e3 / long_haul['cruise_m_s'] / 3600, 1),
                       long_haul_kwh=round(km * long_haul['wh_per_passenger_km'] / 1e3, 1))
            if km <= regional['route_km']:
                row.update(regional_h=round(km * 1e3 / regional['cruise_m_s'] / 3600, 1),
                           regional_kwh=round(km * regional['wh_per_passenger_km'] / 1e3, 1))
            if km <= boat['range_km']:
                row.update(sky_boat_h=round(km * 1e3 / boat['cruise_m_s'] / 3600, 2),
                           sky_boat_kwh=round(km * boat['city_trip_wh_per_km'] / 1e3, 2))
            pairs.append(row)
    local = [dict(mode=m, m_s=round(v, 2), minutes_for_10km=round(1e4 / v / 60, 0)) for m, v in
             (('walking at the measured lunar walk-run transition', GAIT['lunar_transition_m_s']),
              ('pedalled wing', wing['best_glide_speed_m_s']), ('micro glider', glider['best_glide_speed_m_s']),
              ('sky ferry with stops', ferries['commercial_speed_km_h'] / 3.6), ('sky boat', boat['cruise_m_s']))]
    twilight_reach = [dict(lat=la, longest_km_to_a_dusk_or_dawn=round(math.radians(90) * MOON_RADIUS * math.cos(math.radians(la)) / 1e3),
                           winged_liner_h=round(math.radians(90) * MOON_RADIUS * math.cos(math.radians(la)) / liner['cruise_m_s'] / 3600, 1))
                      for la in (0.0, 30.0, 60.0)]

    settings = {}
    for key, (name, kind, light, note) in SETTING_NOTES.items():
        settings[key] = dict(name=name, kind=kind, light=[l for l in light], note=note,
                             places={l: dict(zip(('lat', 'lon'), places[l])) for l in light if l in places})
    settings['summit_metropolis']['ground_above_sea_km'] = ground
    settings['undersea_communities'].update(depth_m=tranquility['depth_m'], outside_pressure_atm=tranquility['pressure_atm'],
                                          nearest_land_km=tranquility['nearest_land_km'])
    settings['roaming_sky_towns']['air'] = air.at(10.0)
    settings['high_platforms']['air'] = air.at(70.0)
    settings['array_habitats']['light'] = ['chosen: continuous sunlight shaped by mirrors and shades']
    settings['earth']['light'] = ['a 24-hour day']

    files = {str(Path(__file__).relative_to(ROOT)): digest(__file__)}
    files.update({c: digest(ROOT / c) for c in CODE})
    optional = {CLIMATOLOGY: digest(ROOT / CLIMATOLOGY) if clim is not None else None}
    aerial_day = next(r for r in paths if r['height_km'] == reference['height_km'])
    return dict(
        schema=SCHEMA,
        producer=dict(study='daily_life', files=files, constants=constants_used([__file__]),
                      inputs={p: digest(ROOT / p) for p in INPUTS}, optional_inputs=optional,
                      main_inputs={MAIN_PLACES['path']: dict(commit=MAIN_PLACES['commit'], sha256=MAIN_PLACES['sha256'][:16],
                                                             read_from=main_source,
                                                             bind='Read in-tree once the branch merges main; until then from Git at the commit.')}),
        evidence=('Clear-sky light on level ground from the calendar\'s solved column in the mean geometry (no libration, '
                  'refraction, terrain, cloud or haze), checked against main\'s corrected per-place Earthlight; global-mean '
                  'air by height and the zonal-wind product\'s 3-day-mean winds, drifting days and passive routes from the '
                  'design run (large-scale flow at about 170 km, no steering); closed-form arithmetic on the sky-fleet, '
                  'inhabited-volume, summit-tower, provisioning and conservation products, shared constants and cited '
                  'values. No settlement, schedule, health outcome or social arrangement is modelled.'),
        reading_rule=('Hours are Earth hours (3,600 s) in the 708.73-hour synodic cycle; light in lux on an unobstructed level '
                      'surface; states partition the cycle by the Sun: day (its centre up), dusk (sunset until its sky falls to '
                      'the practical-dusk level), night (its sky below the level), dawn (the mirror of dusk); where its sky never '
                      'falls to the level, dusk runs into dawn at local midnight. The Earth\'s light is counted within the night: '
                      'below_dusk_h are the hours the total light stays below the level, earth_holds_above_dusk_h the night hours '
                      'the Earth keeps above it. Places are (lat, east lon); far-side and nearside rows lie '
                      'on the 180 and 0 meridians. phase_offset_h is how many hours the second place\'s cycle runs ahead of the '
                      'first\'s. Circadian levels are Brown et al.\'s melanopic EDI at the eye, compared here with horizontal '
                      'illuminance as a reference only. Sky-town days are Earth days, read from the zonal-wind product\'s '
                      'arrays as area-weighted means over the model rows in each band of |latitude| (both hemispheres); '
                      'day_at_mean_wind_days applies T v_sun / |v_sun + u| to the band\'s mean wind. Routes are passive, at fixed '
                      'heights, carried by the 3-day means; days_by_quarter splits the ten-year mean day by the time spent in '
                      'each quarter of hour angle over years 5-10, lingers_at_hour_deg is the centre of the 30-degree sector '
                      'holding most of that time (0 noon, +90 sunset; on a 24-hour clock of the lunar day, 12:00 + hour/15) and '
                      'linger_strength its share against an even pace (1 is even).'),
        levels=dict(levels, dusk_basis=f"Earth control's light at the end of civil twilight (Sun {CIVIL_TWILIGHT_END_DEG:g} deg), solved_sky.json",
                    dusk_sun_deg=round(dusk_deg, 2), dark_sun_deg=round(sky.depression_for(DARK_LUX), 2),
                    circadian_reference=CIRCADIAN_LUX),
        cycle=dict(hours=round(CYCLE_H, 3), sun_deg_per_hour=round(DEG_PER_H, 5), days_of_24h=round(CYCLE_H / 24, 2),
                   terminator=terminator,
                   never_falls_to_practical_dusk_poleward_of_deg=round(90 + dusk_deg, 2),
                   never_falls_to_dark_poleward_of_deg=round(90 + sky.depression_for(DARK_LUX), 2),
                   seasonal_day={f'{la:g}': seasonal_day(la) for la in (0.0, 30.0, 60.0, 75.0, 85.0)}),
        light=dict(calendars=calendars, twilight=[twilight_curve(sky, la, curve_levels) for la in latitudes],
                   cross_check_main_places=check),
        settings=settings,
        twilight=dict(weather=weather, gathering=gathering, ferry_line_per_hour=line, reach=twilight_reach),
        volume=dict(o2_share=round(air.o2_share, 4), open_air=open_air, heights=heights, horizon=horizon,
                    sky_towns=towns, drift_formula=drift_formula, hold_the_hour=hold_hour, wind_holds_the_hour=wind_holds,
                    routes=paths,
                    wind_by_height=[dict(height_km=h, u_m_s=round(winds.global_u(h), 2))
                                    for h in (2.5, 5.0, 7.5, 10.0, 20.0, 30.0, 40.0, 55.0, 75.0)],
                    aerial_district=dict(residents=round(reference['residents']), radius_m=reference['radius_m'],
                                         height_km=reference['height_km'], floor_m2=round(reference['residents'] * volume['assumptions']['floor_m2_per_person']),
                                         held_5_m_s_kw_per_resident=round(kw_per_resident_5, 1),
                                         drifting_day=dict(routes=aerial_day,
                                                           by_latitude=[r for r in towns if r['height_km'] == reference['height_km']]),
                                         wind_10km=next(r for r in ships['air'] if r['height_km'] == 10.0)['gcm_3day'],
                                         ring_wind_10km=next(r for r in ships['air'] if r['height_km'] == 10.0)['ring_3hourly']),
                    shadow=shadow, vertical=vertical, under_canopy=under_canopy),
        bodies=bodies, work=work, earth=earth, shared=shared,
        travel=dict(pairs=pairs, local=local),
        unresolved=['Light at height above the ground and under cloud; refraction in the tall air; terrain horizons.',
                    'Libration moves the Earth by up to about 8 degrees in each place\'s sky; the dated calendar gives it.',
                    'How the circadian system, sleep and children\'s development respond to these light cycles and a sixth of '
                    'Earth\'s gravity: no human data exist for long stays below 1 g.',
                    'Weather downtime, rescue and traffic for elevated and aerial districts; storms and gusts inside the 3-day '
                    'mean winds; towns that steer, sail on the shear or change height.'])


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results(), indent=1) + '\n')
    print(OUT.relative_to(ROOT))


if __name__ == '__main__':
    main()
