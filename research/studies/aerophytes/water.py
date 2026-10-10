"""Water from the air: rain, vapour, cloud droplets and dew, against an aerophyte's need.

Pure functions plus evaluate(). Water in kg per projected m² of organism (1 kg/m²
is 1 mm). Rain comes from the design GCM climatology; the clear-air humidity
profile from the equatorial CM1 box; rain timing and storms from the CM1 ring.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from shared.constants import GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS
from research.studies.aerophytes.environment import (
    LATENT_J_KG, WATER_MOLAR_KG, saturation_pa, water_for_assimilation)

ROOT = Path(__file__).resolve().parents[3]
INPUT_FILES = (
    'climate/gcm/products/climatology_A28_dim5_moon.npz',
    'climate/results/crm/mixed_phase_box_0e.json',
    'climate/results/crm/ring_ring_equator.json',
    'climate/results/crm/ring_ring_equator_fields.npz',
    'research/studies/sky_ships/results/sky_ships.json',
    'research/studies/aerophytes/water_sources.json',
)
AIR_CP_J_KG_K = 1005.  # dry-air engineering value for the psychrometric balance
LEWIS_NUMBER = .85     # water vapour in air, thermal over mass diffusivity (Chilton-Colburn analogy)
AIR_VISCOSITY_PA_S = 1.8e-5  # near 15 C
WATER_DENSITY_KG_M3 = 1000.


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


def _unit(name, value):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{name} must lie in [0, 1]')


def vapour_density(temperature_c, relative_humidity=1.):
    """kg water vapour per m³ of air."""
    _unit('relative humidity', relative_humidity)
    return relative_humidity*saturation_pa(temperature_c)*WATER_MOLAR_KG/(GAS_CONSTANT*(temperature_c+273.15))


# --- the need ---------------------------------------------------------------------

def daily_need(net_leaf_c_kg_m2_year, air_c, leaf_above_air_k, relative_humidity, pressure_pa,
               ci_ca=.7, co2_ppm=400.):
    """Transpiration that buys a year's net leaf carbon, as kg/m² a day and a year."""
    out = water_for_assimilation(net_leaf_c_kg_m2_year*1000, air_c, air_c+leaf_above_air_k,
                                 relative_humidity, pressure_pa, co2_ppm=co2_ppm, ci_ca=ci_ca)
    return dict(water_kg_m2_year=out['water_kg_m2_year'],
                water_kg_m2_day=out['water_kg_m2_year']/JULIAN_YEAR_DAYS,
                water_kg_per_kg_c=out['water_kg_m2_year']/net_leaf_c_kg_m2_year if net_leaf_c_kg_m2_year else None,
                vpd_pa=out['vpd_pa'], ci_ca=ci_ca, relative_humidity=relative_humidity,
                leaf_above_air_k=leaf_above_air_k)


# --- rain -------------------------------------------------------------------------

def zonal_rain(pr_mm_day, latitudes_deg):
    """Zonal-mean rain per model row and its cos-latitude-weighted Moon mean."""
    pr = np.asarray(pr_mm_day, dtype=float)
    lat = np.asarray(latitudes_deg, dtype=float)
    if pr.ndim != 2 or pr.shape[0] != lat.size:
        raise ValueError('rain must be [lat][lon] with matching latitudes')
    zonal = pr.mean(axis=1)
    weights = np.cos(np.radians(lat))
    return dict(latitude_deg=[float(x) for x in lat], rain_mm_day=[float(x) for x in zonal],
                moon_mean_mm_day=float((zonal*weights).sum()/weights.sum()))


def rain_share_meeting(pr_mm_day, latitudes_deg, need_mm_day, catch_share=1., height_factor=1.):
    """Area share of the Moon where mean rain caught aloft meets a daily need."""
    _nonnegative(need=need_mm_day, height=height_factor)
    _unit('catch share', catch_share)
    pr = np.asarray(pr_mm_day, dtype=float)*catch_share*height_factor
    weights = np.cos(np.radians(np.asarray(latitudes_deg, dtype=float)))[:, None]*np.ones_like(pr)
    return float((weights*(pr >= need_mm_day)).sum()/weights.sum())


def rain_closing_latitude(zonal, need_mm_day, catch_share=1., height_factor=1.):
    """Poleward edge, each hemisphere, where the zonal rain falls to the need.

    Takes the most poleward row that meets the need and interpolates linearly in
    latitude to the next row poleward; a hemisphere met to its last row returns it.
    """
    pairs = sorted(zip(zonal['latitude_deg'], zonal['rain_mm_day']))
    out = {}
    for key, rows in (('north_deg', [x for x in pairs if x[0] > 0]),
                      ('south_deg', sorted((x for x in pairs if x[0] < 0), reverse=True))):
        supply = [(l, r*catch_share*height_factor) for l, r in rows]
        meeting = [i for i, (_, r) in enumerate(supply) if r >= need_mm_day]
        if not meeting:
            out[key] = None
            continue
        i = meeting[-1]
        if i == len(supply)-1:
            out[key] = supply[i][0]
            continue
        (l0, r0), (l1, r1) = supply[i], supply[i+1]
        out[key] = l0 + (r0-need_mm_day)/(r0-r1)*(l1-l0)
    return out


def rain_height_factor(fields):
    """Downward rain flux at height over the surface flux, in the ring's heaviest storm.

    One instantaneous storm: rain water mixing ratio times air density times a
    fall speed scaled as (rho0/rho)^0.54, summed over the storm's columns.
    """
    z = np.asarray(fields['storm_height_km'], dtype=float)
    qr = np.asarray(fields['storm_rain_water_g_kg'], dtype=float)
    rho = np.asarray(fields['air_density_kg_m3'], dtype=float)
    if qr.shape[0] != z.size or rho.size != z.size:
        raise ValueError('inconsistent storm section')
    flux = (rho[:, None]*qr*((rho[0]/rho)**.54)[:, None]).sum(axis=1)
    if flux[0] <= 0:
        raise ValueError('no surface rain flux in the storm section')
    return dict(height_km=[float(x) for x in z], flux_over_surface=[float(x) for x in flux/flux[0]])


def storage_between_rains(need_kg_m2_day, dry_days):
    """Free water a body must hold to bridge a dry interval at its daily need."""
    _nonnegative(need=need_kg_m2_day, dry=dry_days)
    return need_kg_m2_day*dry_days


def daylight_in_gap(sector_start_deg, sector_end_deg):
    """Degrees of daylight hour angle (-90 to 90) in the dry gap after a rain sector.

    The gap runs from the sector's end, past midnight, to its start a cycle later.
    """
    gap = (sector_end_deg % 360, sector_start_deg % 360 + 360 if sector_start_deg % 360 <= sector_end_deg % 360
           else sector_start_deg % 360)
    lit = 0.
    for a, b in ((-90., 90.), (270., 450.)):
        lit += max(0., min(b, gap[1]) - max(a, gap[0]))
    return lit


def storage_daylight(need_kg_m2_day, relative_period_days, sector_start_deg, sector_end_deg):
    """Water a body transpires between rains when it transpires by day only.

    The need is a cycle mean, so the daylight rate is twice it; the stored water is
    that rate times the daylight days in the dry gap.
    """
    _nonnegative(need=need_kg_m2_day, period=relative_period_days)
    lit_days = daylight_in_gap(sector_start_deg, sector_end_deg)/360*relative_period_days
    return dict(daylight_dry_days=lit_days, storage_kg_m2=2*need_kg_m2_day*lit_days)


def dry_interval(relative_period_days, rain_sector_fraction):
    """Days between rain sectors for a body crossing them once per relative solar period."""
    _nonnegative(period=relative_period_days)
    _unit('rain sector', rain_sector_fraction)
    return relative_period_days*(1-rain_sector_fraction)


# --- the trim rule ----------------------------------------------------------------

def trim_pressure_per_kg(gas_column_m, air_density_kg_m3, ambient_pressure_pa):
    """Superpressure that holds the trim through a 1 kg/m² water swing with all gas kept.

    A fixed outer volume V (gas column V/A per projected m²) that takes in or
    loses dense air in a ballonet, or rides to a new density level, needs
    dp = dM p/(rho V/A). Round giants have tall columns; flat colonies short ones.
    """
    _positive(column=gas_column_m, density=air_density_kg_m3, pressure=ambient_pressure_pa)
    return ambient_pressure_pa/(air_density_kg_m3*gas_column_m)


def allowed_swing(spare_pressure_pa, gas_column_m, air_density_kg_m3, ambient_pressure_pa):
    """Largest water swing, kg/m², a spare superpressure absorbs without venting gas."""
    _nonnegative(spare=spare_pressure_pa)
    return spare_pressure_pa/trim_pressure_per_kg(gas_column_m, air_density_kg_m3, ambient_pressure_pa)


# --- vapour -----------------------------------------------------------------------

def hygroscopic_uptake(air_c, relative_humidity, water_activity, pressure_pa,
                       convective_w_m2_k=5., radiative_w_m2_k=5.4, air_density_kg_m3=1.2):
    """Steady vapour uptake by a sorbing surface held at water activity a_w.

    Uptake m = h_m (rho_v,air - a_w rho_sat(T_s)) with h_m = h_c/(rho c_p Le^(2/3))
    (Chilton-Colburn, Le = 0.85); the latent heat released warms the surface until
    (h_c + h_r)(T_s - T_a) = L m. Uptake stops where a_w reaches the ambient RH.
    """
    _unit('relative humidity', relative_humidity)
    _unit('water activity', water_activity)
    _positive(pressure=pressure_pa, hc=convective_w_m2_k, density=air_density_kg_m3)
    _nonnegative(hr=radiative_w_m2_k)
    hm = convective_w_m2_k/(air_density_kg_m3*AIR_CP_J_KG_K*LEWIS_NUMBER**(2/3))
    rv = vapour_density(air_c, relative_humidity)

    def residual(ts):
        m = hm*(rv - water_activity*vapour_density(ts, 1.))
        return (convective_w_m2_k+radiative_w_m2_k)*(ts-air_c) - LATENT_J_KG*m

    if water_activity >= relative_humidity:
        return dict(surface_c=air_c, uptake_kg_m2_day=0., latent_w_m2=0., warming_k=0.)
    lo, hi = air_c, air_c+40.
    for _ in range(80):
        mid = (lo+hi)/2
        if residual(mid) < 0:
            lo = mid
        else:
            hi = mid
    ts = (lo+hi)/2
    m = hm*(rv - water_activity*vapour_density(ts, 1.))
    return dict(surface_c=ts, warming_k=ts-air_c, uptake_kg_m2_day=m*JULIAN_DAY,
                latent_w_m2=LATENT_J_KG*m)


def vapour_work(relative_humidity, tissue_activity, temperature_c):
    """Least work, J per kg, to move water from air at RH into tissue at activity a."""
    _unit('relative humidity', relative_humidity)
    _unit('tissue activity', tissue_activity)
    if relative_humidity == 0:
        raise ValueError('no vapour to take')
    t = temperature_c+273.15
    return max(0., GAS_CONSTANT*t/WATER_MOLAR_KG*math.log(tissue_activity/relative_humidity))


def solar_regeneration_heat(uptake_kg_m2_day):
    """Solar heat, W/m², that drives absorbed water back out of a sorbent for use."""
    _nonnegative(uptake=uptake_kg_m2_day)
    return uptake_kg_m2_day*LATENT_J_KG/JULIAN_DAY


def dew_threshold(air_c, surface_depression_k):
    """Least ambient RH at which a surface that much colder than the air gathers dew."""
    _nonnegative(depression=surface_depression_k)
    return saturation_pa(air_c-surface_depression_k)/saturation_pa(air_c)


def radiative_depression(net_cooling_w_m2, convective_w_m2_k, air_c, emissivity=.96):
    """Steady sub-air temperature of an upward face losing net longwave to the sky."""
    _nonnegative(cooling=net_cooling_w_m2)
    _positive(hc=convective_w_m2_k)
    hr = 4*emissivity*STEFAN_BOLTZMANN*(air_c+273.15)**3
    return net_cooling_w_m2/(convective_w_m2_k+hr)


# --- cloud droplets ---------------------------------------------------------------

def stokes_number(droplet_diameter_m, collector_diameter_m, relative_speed_m_s,
                  viscosity_pa_s=AIR_VISCOSITY_PA_S):
    """St = 2 rho_w r² v/(9 mu R), droplet radius r and fibre radius R (Park et al. 2013)."""
    _positive(droplet=droplet_diameter_m, collector=collector_diameter_m,
              speed=relative_speed_m_s, viscosity=viscosity_pa_s)
    r, big = droplet_diameter_m/2, collector_diameter_m/2
    return 2*WATER_DENSITY_KG_M3*r*r*relative_speed_m_s/(9*viscosity_pa_s*big)


def impaction_efficiency(stokes):
    """Droplet deposition on a long cylinder, St/(St + pi/2): Langmuir and Blodgett's fit."""
    _nonnegative(stokes=stokes)
    return stokes/(stokes+math.pi/2)


def droplet_capture(liquid_g_m3, relative_speed_m_s, efficiency, collector_area_ratio=1., duty=1.):
    """Cloud water gathered per projected m², kg a day, by collectors in a moving stream."""
    _nonnegative(lwc=liquid_g_m3, speed=relative_speed_m_s, area=collector_area_ratio)
    _unit('efficiency', efficiency)
    _unit('duty', duty)
    return liquid_g_m3/1000*relative_speed_m_s*efficiency*collector_area_ratio*duty*JULIAN_DAY


def evaluate(root=ROOT, net_leaf_c_kg_m2_year=1.83, relative_periods_days=None):
    """Needs by trait set, rain by latitude, vapour and droplet capacities, trim limits.

    relative_periods_days: noon-to-noon periods an aerophyte sees (label -> days),
    from the sailing component; a fixed place's synodic month is always included.
    """
    root = Path(root)
    clim = np.load(root/INPUT_FILES[0], allow_pickle=False)
    meta = json.loads(str(clim['metadata']))
    if meta.get('schema') != 'terluna.climate.gcm-climatology/1':
        raise ValueError('unexpected climatology')
    box = json.loads((root/INPUT_FILES[1]).read_text())
    ring = json.loads((root/INPUT_FILES[2]).read_text())
    if ring.get('schema') != 'terluna.climate.crm-ring/1':
        raise ValueError('unexpected ring product')
    fields = np.load(root/INPUT_FILES[3], allow_pickle=False)
    ships = json.loads((root/INPUT_FILES[4]).read_text())
    air = {a['height_km']: a for a in ships['air']}
    profile = box['mean_profile']
    z = profile['z_km']

    def at(h, key):
        return float(np.interp(h, z, profile[key]))

    rh = {h: at(h, 'relative_humidity_clear') for h in (5, 10, 15, 20, 25)}
    rh_out = {str(h): v for h, v in rh.items()}
    zz = np.asarray(profile['z_km'])*1000
    above = zz >= 10000
    vapour_column = float(np.trapezoid((np.asarray(profile['air_density_kg_m3'])*np.asarray(profile['vapour_g_kg'])/1000)[above], zz[above]))
    t10, p10 = air[10]['temperature_c'], air[10]['pressure_atm']*101325.
    traits = [('reference', .6, 5., .7), ('box_humidity_cool_leaves', rh[10], 0., .7),
              ('humid_layer_cool_leaves', .9, 0., .7), ('concentrating_cool_leaves', rh[10], 0., .3),
              ('concentrating_humid_layer', .9, 0., .3), ('concentrating_saturated_layer', .98, 0., .3)]
    needs = [dict(name=n, **daily_need(net_leaf_c_kg_m2_year, t10, d, r, p10, ci_ca=c))
             for n, r, d, c in traits]
    zonal = zonal_rain(clim['pr_mm_day'], clim['lat'])
    storm = dict(storm_height_km=fields['storm_height_km'],
                 storm_rain_water_g_kg=fields['storm_rain_water_g_kg'],
                 air_density_kg_m3=np.interp(fields['storm_height_km'], z, profile['air_density_kg_m3']))
    height = rain_height_factor(storm)
    h10 = float(np.interp(10., height['height_km'], height['flux_over_surface']))
    h20 = float(np.interp(20., height['height_km'], height['flux_over_surface']))
    closures = []
    for need in needs:
        for catch in (.5, 1.):
            for factor in (1., h10):
                closures.append(dict(name=need['name'], need_mm_day=need['water_kg_m2_day'],
                    catch_share=catch, height_factor=factor,
                    moon_share=rain_share_meeting(clim['pr_mm_day'], clim['lat'], need['water_kg_m2_day'], catch, factor),
                    **rain_closing_latitude(zonal, need['water_kg_m2_day'], catch, factor)))
    hours = ring['by_hour_angle']['hour_angle_deg']
    land_rain = ring['by_hour_angle']['land']['rain_mm_h']
    sea_rain = ring['by_hour_angle']['water']['rain_mm_h']
    sector = [h for h, r in zip(hours, land_rain) if r >= .05]
    sector_fraction = len(sector)/len(hours)
    half_bin = (hours[1]-hours[0])/2  # the sector runs between its bins' outer edges
    periods = dict(relative_periods_days or {})
    periods['fixed_place'] = SYNODIC_MONTH_DAYS
    start, end = (min(sector)-half_bin, max(sector)+half_bin) if sector else (0., 0.)
    intervals = []
    for label, period in periods.items():
        lit = storage_daylight(1., period, start, end)['daylight_dry_days']
        intervals.append(dict(path=label, relative_period_days=period, rain_sector_fraction=sector_fraction,
                              rain_every_days=period, dry_days=dry_interval(period, sector_fraction),
                              daylight_dry_days=lit,
                              storage_kg_m2={n['name']: storage_daylight(n['water_kg_m2_day'], period, start, end)['storage_kg_m2']
                                             for n in needs}))
    trim = [dict(shape=s, gas_column_m=c, pa_per_kg_m2=trim_pressure_per_kg(c, air[10]['density_kg_m3'], p10),
                 **{f'swing_kg_m2_at_{q:g}_pa': allowed_swing(q, c, air[10]['density_kg_m3'], p10) for q in (100., 300., 1000.)})
            for s, c in (('raft_of_40m_modules', 4/3*20*math.pi/(2*math.sqrt(3))), ('quilt_60m', 60.), ('sphere_40m', 26.7),
                         ('sphere_200m', 133.3), ('sphere_500m', 333.3), ('sphere_2km', 1333.3))]
    vapour = [dict(height_km=h, relative_humidity=r, water_activity=a, convective_w_m2_k=hc,
                   **hygroscopic_uptake(air[h]['temperature_c'], r, a, air[h]['pressure_atm']*101325.,
                                        hc, air_density_kg_m3=air[h]['density_kg_m3']))
              for h, r in ((10, rh[10]), (10, .9), (20, rh[20])) for a in (.3, .5) for hc in (2., 5., 10.)]
    work = [dict(relative_humidity=r, tissue_activity=.99, j_per_kg=vapour_work(r, .99, t10),
                 w_m2_per_kg_m2_day=vapour_work(r, .99, t10)/JULIAN_DAY) for r in (rh[10], .8, .9, .98)]
    droplets = []
    for d in (5e-6, 10e-6, 20e-6):
        for collector in (1e-4, 5e-4, 2e-3):
            for v in (.3, 1., 3.):
                st = stokes_number(d, collector, v)
                droplets.append(dict(droplet_um=d*1e6, collector_mm=collector*1e3, relative_speed_m_s=v,
                                     stokes=st, efficiency=impaction_efficiency(st)))
    capture = [dict(liquid_g_m3=c, relative_speed_m_s=v, efficiency=e, duty=f,
                    kg_m2_day=droplet_capture(c, v, e, 1., f))
               for c in (.05, .1, .2) for v in (.3, 1., 2.) for e in (.1, .3) for f in (.05, .2, 1.)]
    cloud = ring['circulation']
    cloud_rows = [dict(height_km=hz, peak_cloud_fraction=max(row),
                       hour_angle_of_peak_deg=cloud['hour_angle_deg'][row.index(max(row))],
                       mean_cloud_fraction=sum(row)/len(row))
                  for hz, row in zip(cloud['heights_km'], cloud['cloud_fraction'])]
    dew = [dict(air_c=t10, depression_k=d, threshold_rh=dew_threshold(t10, d)) for d in (.5, 1., 2., 5.)]
    radiative = [dict(net_cooling_w_m2=q, convective_w_m2_k=hc, depression_k=radiative_depression(q, hc, t10),
                      threshold_rh=dew_threshold(t10, radiative_depression(q, hc, t10)),
                      energy_ceiling_kg_m2_day=q*JULIAN_DAY/LATENT_J_KG)
                 for q in (5., 10., 20.) for hc in (2., 5.)]
    return dict(
        evidence='Rain from the design GCM climatology (10-year means of 3-day outputs, T21) and the CM1 equatorial ring and box; transpiration from the stomatal diffusion model; vapour, dew and droplet capacities from heat and mass balances with stated transfer coefficients. No trajectory or along-path humidity is simulated.',
        reading_rule='kg water per projected m² (1 kg/m² = 1 mm). Needs buy the reference net leaf carbon and are cycle means; transpiration runs by day at twice the mean. Rain closes where zonal-mean rain x catch share x height factor meets the need, interpolated in latitude. Storage bridges the daylight hours of the dry gap after the afternoon rain sector. The trim rule gives the swing a spare superpressure absorbs with all gas kept.',
        net_leaf_c_kg_m2_year=net_leaf_c_kg_m2_year, box_clear_air_rh=rh_out,
        box_vapour_above_10_km_kg_m2=vapour_column,
        air_10km=dict(temperature_c=t10, pressure_pa=p10, density_kg_m3=air[10]['density_kg_m3']),
        needs=needs, zonal_rain=zonal, rain_height=dict(**height, at_10_km=h10, at_20_km=h20,
            note='One instantaneous storm, the heaviest of the ring: rain still aloft exceeds the surface flux, so 1 is the conservative factor and the storm value an upper sensitivity.'),
        rain_closure=closures,
        ring_rain_by_hour=dict(hour_angle_deg=hours, land_mm_h=land_rain, sea_mm_h=sea_rain,
                               land_sector_hours_deg=[min(sector)-half_bin, max(sector)+half_bin] if sector else None,
                               sector_fraction=sector_fraction,
                               wet_share_land=ring['rain_rate_mm_h']['land']['wet_share'],
                               storm_rain_max_mm_h=ring['storms']['rain_max_mm_h']),
        dry_intervals=intervals, trim=trim, vapour=vapour, vapour_work=work,
        droplet_efficiency=droplets, droplet_capture=capture, cloud_by_height=cloud_rows,
        dew_threshold=dew, radiative_dew=radiative,
        unresolved=['Humidity aloft away from the equatorial box, and along an aerophyte path.',
                    'The longwave sky seen from 10-20 km: the cooling used for dew is a stated range.',
                    'Cloud liquid water outside storm cores, and how often a route meets it.',
                    'A latitude and solar-hour wind and rain product (being built on another branch) for along-path rain timing.'])
