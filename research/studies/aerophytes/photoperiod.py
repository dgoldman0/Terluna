"""Mean Sun-following kinematics and distinct optical/biological light measures.

Reads committed surface products. No atmospheric trajectory, altitude radiative
transfer, wind forecast, propulsion or continuous-light physiology is solved.
"""
from __future__ import annotations

from bisect import bisect_left
import json
import math
from pathlib import Path

from shared.constants import (AVOGADRO, JULIAN_DAY, MOON_EQUATOR_TO_ECLIPTIC_DEG,
                              MOON_RADIUS, PLANCK, SPEED_OF_LIGHT, SYNODIC_MONTH_DAYS)

ROOT = Path(__file__).resolve().parents[3]
SYNODIC_SECONDS = SYNODIC_MONTH_DAYS * JULIAN_DAY
INPUT_FILES = (
    'illumination/surface_light/results/surface_light.json',
    'biosphere/canopy/results/plant.json',
    'biosphere/canopy/results/canopy.json',
    'illumination/calendar/results/transfer.json',
    'research/studies/aerophytes/photoperiod_sources.json',
)


def _latitude(latitude_deg):
    if not -90 <= latitude_deg <= 90:
        raise ValueError('latitude must be in [-90,90] degrees')
    return math.radians(latitude_deg)


def sun_follow_speed(latitude_deg, altitude_m=0.):
    """Westward ground-relative m/s to hold mean solar hour angle at fixed latitude.

    Uses a circular latitude track at radius R+h and the mean synodic period.
    Solar declination and nonuniform ephemeris motion still change illumination.
    At a geographic pole longitude/hour angle is degenerate; speed is zero.
    """
    latitude = _latitude(latitude_deg)
    if altitude_m < 0:
        raise ValueError('nonnegative altitude required')
    if abs(latitude_deg) == 90:
        return 0.
    return 2*math.pi*(MOON_RADIUS+altitude_m)*math.cos(latitude)/SYNODIC_SECONDS


def solar_height(latitude_deg, hour_angle_deg, declination_deg=0.):
    """Geometric solar-centre elevation above the local horizontal, degrees.

    H=0 at local solar noon; no refraction, terrain, finite disk or horizon dip.
    Observer altitude does not add to this local-horizontal elevation.
    """
    latitude = _latitude(latitude_deg)
    if not -90 <= declination_deg <= 90:
        raise ValueError('declination must be in [-90,90] degrees')
    d, h = math.radians(declination_deg), math.radians(hour_angle_deg)
    sine = math.sin(latitude)*math.sin(d)+math.cos(latitude)*math.cos(d)*math.cos(h)
    return math.degrees(math.asin(max(-1.,min(1.,sine))))


def drift_window(latitude_deg, altitude_m, speed_error_m_s,
                 half_width_hour_angle_deg=30.):
    """Time from the centre to an hour-angle window edge at constant speed error.

    Error = actual westward ground speed minus Sun-follow target, m/s. Positive
    error moves toward earlier solar hours. This is residence in an angular
    window, not a biological light threshold or a weather residence time.
    None means no phase drift in this ideal mean-motion calculation.
    """
    latitude = _latitude(latitude_deg)
    if abs(latitude_deg) == 90 or altitude_m < 0 or not 0 < half_width_hour_angle_deg < 180:
        raise ValueError('nonpolar latitude, nonnegative altitude and window (0,180) required')
    radius = (MOON_RADIUS+altitude_m)*math.cos(latitude)
    rate = -math.degrees(speed_error_m_s/radius)*JULIAN_DAY
    exit_days = None if speed_error_m_s == 0 else half_width_hour_angle_deg/abs(rate)
    period = None if speed_error_m_s == 0 else 360/abs(rate)
    return dict(speed_error_m_s=speed_error_m_s,hour_angle_drift_deg_per_day=rate,
                half_width_hour_angle_deg=half_width_hour_angle_deg,
                centre_to_edge_days=exit_days,
                full_window_crossing_days=None if exit_days is None else 2*exit_days,
                relative_solar_period_days=period,
                equinox_geometric_day_days=None if period is None else period/2,
                equinox_geometric_night_days=None if period is None else period/2)


def threshold_duty(latitude_deg, minimum_solar_height_deg, declination_deg=0.):
    """Fraction of a stationary mean solar cycle above a geometric elevation.

    A negative threshold can describe a separately calibrated twilight level.
    It is never assumed to be photosynthetic compensation or H2 light duty.
    """
    latitude = _latitude(latitude_deg)
    if not -90 <= minimum_solar_height_deg <= 90 or not -90 <= declination_deg <= 90:
        raise ValueError('solar height and declination must lie in [-90,90]')
    d = math.radians(declination_deg)
    a, b = math.sin(latitude)*math.sin(d), math.cos(latitude)*math.cos(d)
    threshold = math.sin(math.radians(minimum_solar_height_deg))
    if abs(b) < 1e-14:
        fraction = float(a >= threshold)
    else:
        x = (threshold-a)/b
        fraction = 1. if x <= -1 else 0. if x >= 1 else math.acos(x)/math.pi
    return dict(fraction_above=fraction,
                hours_below_per_stationary_synodic_cycle=(1-fraction)*SYNODIC_SECONDS/3600)


def _interpolate(x, xs, ys):
    if not xs[0] <= x <= xs[-1]:
        raise ValueError('outside stored interpolation interval; no extrapolation')
    k = bisect_left(xs,x)
    if k == 0 or (k < len(xs) and x == xs[k]):
        return float(ys[k])
    f = (x-xs[k-1])/(xs[k]-xs[k-1])
    return (1-f)*ys[k-1]+f*ys[k]


def _product(root, index, schema):
    data = json.loads((Path(root)/INPUT_FILES[index]).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected product schema {data.get("schema")}')
    return data


def _ground_rows(product):
    case = product['cases']['moon_1.2atm/design']
    reflected = product['columns'][case['column']]['reflectance_from_below']
    centres = product['centres_nm']
    if len(centres) != len(reflected) or any(abs(b-a-1)>1e-8 for a,b in zip(centres,centres[1:])):
        raise ValueError('expected stored one-nanometre spectral bins')
    rows = []
    for key, spectrum in case['spectra'].items():
        summary = case['summaries'][key]['0.1']
        direct, diffuse = spectrum['direct_w_m2_nm'], spectrum['diffuse_black_w_m2_nm']
        if len(direct) != len(centres) or len(diffuse) != len(centres):
            raise ValueError('inconsistent stored spectrum lengths')
        # The surface-light product defines this exact ground-return formula.
        total = [(d+f)/(1-.1*r) for d,f,r in zip(direct,diffuse,reflected)]
        rows.append(dict(sun_deg=float(key),par_w_m2=summary['par_w_m2'],
            par_umol_m2_s=summary['photons_umol_m2_s']['par'],
            illuminance_lux=summary['illuminance_lux'],
            optical_202_1000_w_m2=sum(total),
            spectral_par_w_m2=sum(e for w,e in zip(centres,total) if 400<=w<700),
            nonpar_202_1000_w_m2=sum(e for w,e in zip(centres,total) if not 400<=w<700)))
    return sorted(rows,key=lambda row:row['sun_deg'])


def ground_light(sun_deg, root=ROOT, *, rows=None):
    """Stored clear-ground daylight at neutral albedo .1; interpolate in sin(height).

    Valid only for solar heights1..90 degrees. Optical power covers202..1000nm,
    not total solar energy. This older daylight solution underestimates low Sun
    relative to the spherical atlas used by plant.py. No altitude gain is applied.
    """
    if rows is None:
        rows = _ground_rows(_product(root,0,'terluna.illumination.surface-light/1'))
    if not rows[0]['sun_deg'] <= sun_deg <= rows[-1]['sun_deg']:
        raise ValueError('daylight solar height outside stored range')
    x = math.sin(math.radians(sun_deg))
    xs = [math.sin(math.radians(row['sun_deg'])) for row in rows]
    return dict(sun_deg=sun_deg,**{key:_interpolate(x,xs,[row[key] for row in rows])
                for key in rows[0] if key != 'sun_deg'})


def twilight_threshold(par_umol_m2_s, root=ROOT, *, table=None):
    """Invert the stored plant's clear-ground twilight PAR table, no extrapolation."""
    if table is None:
        table = _product(root,1,'terluna.biosphere.plant-carbon-cycle/1')['twilight']['moon']
    return _interpolate(par_umol_m2_s,table['par_umol_m2_s'],table['elevation_deg'])


def _twilight_shape(product):
    """Energy per PAR micromole for plant.py's assumed1deg diffuse-black spectrum."""
    wavelengths = product['centres_nm']
    spectrum = product['cases']['moon_1.2atm/design']['spectra']['1.0']['diffuse_black_w_m2_nm']
    par = [(w,e) for w,e in zip(wavelengths,spectrum) if 400<=w<700]
    photons = sum(e*w*1e-9/(PLANCK*SPEED_OF_LIGHT*AVOGADRO)*1e6 for w,e in par)
    return dict(par_w_per_umol=sum(e for _,e in par)/photons,
                optical_202_1000_w_per_par_umol=sum(spectrum)/photons)


def twilight_energy(sun_deg, root=ROOT, *, table=None, shape=None):
    """Inherited plant spectral approximation, NOT a newly solved twilight spectrum.

    Scales its1deg diffuse-black-ground spectrum by the stored twilight PAR.
    Ground-only, -35..0deg,202..1000nm; cannot be rescaled from calendar lux or
    treated as bolometric/altitude light. H2 efficiency is supplied externally.
    """
    if table is None:
        table = _product(root,1,'terluna.biosphere.plant-carbon-cycle/1')['twilight']['moon']
    if shape is None:
        shape = _twilight_shape(_product(root,0,'terluna.illumination.surface-light/1'))
    par = _interpolate(sun_deg,table['elevation_deg'],table['par_umol_m2_s'])
    return dict(sun_deg=sun_deg,par_umol_m2_s=par,
                par_w_m2=par*shape['par_w_per_umol'],
                optical_202_1000_w_m2=par*shape['optical_202_1000_w_per_par_umol'])


def required_irradiance(chemical_power_w_m2, efficiency):
    """Required incident power in the SAME spectral/area convention as efficiency."""
    if chemical_power_w_m2<0 or not 0<efficiency<=1:
        raise ValueError('nonnegative chemical power and efficiency (0,1] required')
    return chemical_power_w_m2/efficiency


def evaluate(root=ROOT):
    light = _product(root,0,'terluna.illumination.surface-light/1')
    plant = _product(root,1,'terluna.biosphere.plant-carbon-cycle/1')
    canopy = _product(root,2,'terluna.biosphere.canopy-photosynthesis/1')
    calendar = _product(root,3,'terluna.illumination.calendar-transfer/1')
    rows = _ground_rows(light)
    twilight = plant['twilight']['moon']
    shape = _twilight_shape(light)
    latitudes = (0.,30.,60.,70.,80.,85.)
    altitudes = (0.,10000.,20000.,40000.)
    plant_rows = []
    for side,bands in plant['results'].items():
        for latitude,strategies in bands.items():
            for strategy in ('earth_like','idle_quarter','idle_tenth'):
                p = strategies[strategy]
                plant_rows.append(dict(side=side,latitude_deg=float(latitude),strategy=strategy,
                    par_below_1_hours=p['dark_hours'],maintenance_deficit_hours=p['deficit_hours'],
                    reserve_g_c_m2=p['store'],gpp_g_c_m2_cycle=p['gpp'],
                    leaf_respiration_g_c_m2_cycle=p['leaf_respiration'],
                    stem_root_respiration_g_c_m2_cycle=p['stem_root_respiration'],
                    sun_below_horizon_respiration_g_c_m2_cycle=p['respiration_in_the_dark'],
                    twilight_gpp_g_c_m2_cycle=p['twilight_gpp']))
    locks=[]
    for latitude in latitudes:
        for hour in (0.,60.):
            for declination in (-MOON_EQUATOR_TO_ECLIPTIC_DEG,0.,MOON_EQUATOR_TO_ECLIPTIC_DEG):
                height=solar_height(latitude,hour,declination)
                locks.append(dict(latitude_deg=latitude,hour_angle_deg=hour,declination_deg=declination,
                    solar_height_deg=height,
                    ground_light=ground_light(height,rows=rows) if height>=rows[0]['sun_deg'] else None))
    return dict(
        evidence='Mean-motion kinematics and inherited clear-ground optical/plant calculations. No controllable Sun-following weather route or aerial carbon/H2 budget has been demonstrated.',
        reading_rule='Westward ground speed is distinct from air-relative propulsion. PAR darkness, geometric night, whole-plant maintenance deficit, visual illumination and H2 production thresholds are separate quantities. Ground PAR and202–1000nm irradiance are not altitude or bolometric irradiance.',
        assumptions=dict(synodic_days=SYNODIC_MONTH_DAYS,moon_radius_m=MOON_RADIUS,
            declination_bound_deg=MOON_EQUATOR_TO_ECLIPTIC_DEG,
            note='Fixed latitude, mean solar motion and steady speed-error scenarios. Real ephemeris, weather and solar declination vary. No terrain, refraction, eclipse, cloud, horizon-dip or altitude transport correction.'),
        sun_follow=[dict(latitude_deg=latitude,altitude_m=altitude,
                        westward_ground_speed_m_s=sun_follow_speed(latitude,altitude))
                    for latitude in latitudes for altitude in altitudes],
        mismatch=[dict(latitude_deg=latitude,altitude_m=20000.,
                       **drift_window(latitude,20000.,error))
                  for latitude in latitudes for error in (-5.,-2.,-1.,0.,1.,2.,5.)],
        passive_wind_windows=[dict(latitude_deg=latitude,altitude_m=20000.,
            wind_west_m_s=wind,target_west_m_s=sun_follow_speed(latitude,20000.),
            **drift_window(latitude,20000.,wind-sun_follow_speed(latitude,20000.)))
            for latitude in latitudes for wind in (-5.,-2.,-1.,0.,1.,2.,5.)],
        sun_locked_ground_reference=locks,ground_daylight_table=rows,
        ground_twilight=twilight,
        twilight_energy_assumption=dict(
            reading_rule='The inherited plant model assumes the1deg diffuse-black spectrum at every twilight angle. These spectral energy proxies are not the newer calendar solution, altitude predictions or measured H2 action spectra. No flux is forced to zero below a trait threshold.',
            shape=shape,
            rows=[twilight_energy(h,table=twilight,shape=shape) for h in twilight['elevation_deg']]),
        calendar_ground_visual_reference=dict(
            reading_rule='Solved ground point-source diffuse solar illuminance, lux; finite solar disk and Earthlight not integrated here. Calendar stores no ground spectral irradiance for energy conversion.',
            rows=[dict(sun_deg=h,diffuse_lux=_interpolate(h,calendar['diffuse_elevation_deg'],calendar['solar']['diffuse_lux']))
                  for h in (-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.)]),
        threshold_example=dict(chemical_power_w_m2=.07,efficiency=.01,
            required_incident_w_m2=required_irradiance(.07,.01),
            note='The7W/m² threshold applies only if1% efficiency is valid on the same input spectrum/area at twilight. Full-solar laboratory efficiency is not automatically a twilight or202–1000nm efficiency. This is not a default H2 duty.'),
        stationary_twilight_duty=[dict(latitude_deg=latitude,par_threshold_umol_m2_s=threshold,
            minimum_solar_height_deg=twilight_threshold(threshold,table=plant['twilight']['moon']),
            **threshold_duty(latitude,twilight_threshold(threshold,table=plant['twilight']['moon'])))
            for latitude in latitudes for threshold in (1.,10.,50.)],
        plant_reference=plant_rows,
        plant_reference_limits='LAI5 ground stand, fixed zero-declination geometry, terrestrial C3 model and near/far land temperature traces. near/far PAR darkness is identical because Earthlight is absent from plant.py light inputs. Respiration is not an aerophyte calibration.',
        reference_canopy_by_sun_height=canopy['by_sun_height']['moon'],
        unresolved=['Air-relative velocity vector, propulsion and controllable wind-shear path; see navigation component.',
            'Altitude-resolved spectral and angular light, cloud/shadow encounters, eclipse response and seasons.',
            'C3 compensation at actual tissue temperature, water/CO2 supply and constant-light acclimation.',
            'Separate wavelength/action-spectrum and minimum-flux response for any photobiological H2 pathway.',
            'Actual moving-trajectory carbon, water, heat, nutrient and gas-reserve integrals. No GCM run made.'])
