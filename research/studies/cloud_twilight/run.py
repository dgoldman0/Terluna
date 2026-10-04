"""Cloud-twilight mechanism screen using only committed Terluna inputs.

Run: OPENBLAS_NUM_THREADS=1 python -m research.studies.cloud_twilight.run
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from illumination.cloud_light.model import (
    SUN_RADIUS_DEG, atlas_proxy, cloud_view, photopic_weights, solved_moon,
    sunset_depression_deg,
)
from research.studies.optical_comfort.angular import PackedAtlas, Y_WEIGHTS
from shared.constants import EARTH_RADIUS, MOON_RADIUS, SYNODIC_MONTH_DAYS

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CRM = ROOT / 'climate/results/crm/ring_ring_equator.json'
FIELDS = CRM.with_name('ring_ring_equator_fields.npz')
ATLAS = ROOT / 'immersion/assets/sky/atmosphere.json'
SCHEMA = 'terluna.research.cloud-twilight/1'
HOURS_PER_DEGREE = SYNODIC_MONTH_DAYS * 24 / 360
HEIGHTS_KM = (20., 40., 60., 80.)
DEPRESSION_DEG = np.arange(0., 20.0001, .25)
RANGES_KM = (-200., -100., -50., 0., 25., 50., 75., 100., 150., 200., 300., 400.)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_clouds():
    summary = json.loads(CRM.read_text())
    with np.load(FIELDS, allow_pickle=False) as fields:
        metadata = json.loads(str(fields['metadata']))
        if (summary['schema'] != 'terluna.climate.crm-ring/1' or
                metadata['schema'] != 'terluna.climate.crm-ring-fields/1' or
                summary['case'] != 'ring_equator' or metadata['case'] != summary['case'] or
                summary['gcm']['run'] != 'A28_dim5_moon' or
                summary['span_days'] != metadata['span_days'] or summary['tilt_deg'] != 0):
            raise ValueError('Study requires matching corrected equatorial-ring products')
        height = fields['section_height_km'].astype(float)
        hour = fields['section_hour_angle_deg'].astype(float)
        cloud = fields['section_cloud_fraction'].astype(float)
        if (cloud.shape != (len(hour), len(height)) or not np.isfinite(cloud).all() or
                np.any((cloud < 0) | (cloud > 1)) or np.any(np.diff(height) <= 0) or
                np.any(np.diff(hour) <= 0)):
            raise ValueError('Invalid cloud occurrence grid')
        section = dict(height_km=height.tolist(), hour_angle_deg=hour.tolist(), fraction=cloud.tolist(),
                       rule='At each height/local-time bin, fraction of sampled columns with qc+qi >= 1e-5 kg/kg. '
                            'Not sky coverage, whole-cloud probability, optical depth or event duration.')
        context = dict(peak_rain_storm_day=float(fields['storm_day']),
                       peak_rain_storm_hour_angle_deg=float(fields['storm_hour_angle_deg']),
                       rule='The committed morphology slice is a daytime maximum-rain snapshot. '
                            'Its condensate combines qc+qi+qs+qg; it is not a twilight optical volume.')
    return summary, section, context


def path_row(column, weights, height_km, range_km, depression_deg, order=32):
    h = height_km * 1000
    view = cloud_view(column.radius_m, h, range_km * 1000)
    local_sun = view['arc_angle_deg'] - depression_deg
    solar = column.sun_transmission(h, local_sun, order)
    sight = column.transmission(1.6, view['mu'], view['ray_length_m'])
    _, blocked = column.path(1.6, view['mu'], view['ray_length_m'])
    # Incoming propagation is -sun; outgoing propagation is -view, so their
    # dot product is sun dot view. Zero degrees is forward scattering.
    dep = math.radians(depression_deg)
    cos_scatter = math.cos(dep) * view['ray_x'] - math.sin(dep) * view['ray_z']
    total = float(weights.sum())
    return dict(height_km=height_km, westward_surface_range_km=range_km,
                observer_sun_depression_deg=depression_deg,
                hours_after_surface_solar_centre_set=depression_deg * HOURS_PER_DEGREE,
                cloud_local_sun_elevation_deg=local_sun,
                cloud_hour_angle_deg=90 - local_sun,
                apparent_elevation_deg=view['apparent_elevation_deg'],
                scattering_angle_deg=math.degrees(math.acos(float(np.clip(cos_scatter, -1, 1)))),
                line_of_sight_blocked=blocked,
                cloud_direct_normal_lux=float(solar @ weights),
                cloud_solar_transmission=float(solar @ weights / total),
                observer_path_source_weighted_transmission=float(sight @ weights / total),
                two_path_source_fraction=float((solar * sight) @ weights / total),
                solar_transmission_450_550_650nm=[float(solar[list(column.wavelength_nm).index(w)])
                                                for w in (450, 550, 650)],
                two_path_transmission_450_550_650nm=[float((solar * sight)[list(column.wavelength_nm).index(w)])
                                                   for w in (450, 550, 650)])


def last_reference_level(column, weights, height_km, level_lux):
    """Post-sunset crossing of an arbitrary incident-light reference level."""
    end = sunset_depression_deg(column.radius_m, height_km * 1000) + SUN_RADIUS_DEG + .01
    beam = lambda d: float(column.sun_transmission(height_km * 1000, -d) @ weights)
    if beam(0) < level_lux:
        return None
    depression = brentq(lambda d: beam(d) - level_lux, 0, end, xtol=1e-7)
    return dict(depression_deg=depression, hours_after_surface_solar_centre_set=depression * HOURS_PER_DEGREE)


def provenance():
    sources = [Path(__file__), HERE / 'sources.json', ROOT / 'illumination/cloud_light/model.py',
               ROOT / 'illumination/surface_light/model.py', ROOT / 'illumination/sky/atmospheres.py',
               ROOT / 'illumination/sky/colour_matching.py', ROOT / 'research/studies/optical_comfort/angular.py',
               ROOT / 'shared/constants.py', ROOT / 'shared/constants.json',
               ROOT / 'atmosphere/middle_atmosphere/equilibrium.py',
               ROOT / 'atmosphere/radiative_convective/thermodynamics.py',
               ROOT / 'atmosphere/radiative_convective/optics.py',
               ROOT / 'atmosphere/radiative_convective/spectroscopy.py',
               ROOT / 'atmosphere/radiative_convective/climate.py']
    inputs = [CRM, FIELDS, ATLAS, ROOT / 'protection/spectra/shield_transmission.json',
              ROOT / 'atmosphere/middle_atmosphere/results/middle_atmosphere.json']
    for name in ('moon_1.2atm_titania_stack', 'moon_1.2atm_titania_stack_ts298'):
        inputs += [ROOT / f'atmosphere/middle_atmosphere/results/profiles/{name}.csv',
                   ROOT / f'atmosphere/middle_atmosphere/results/cases/{name}.json']
    return dict(lane='research', runner=str(Path(__file__).relative_to(ROOT)),
                files={str(p.relative_to(ROOT)): digest(p) for p in sources},
                inputs={str(p.relative_to(ROOT)): dict(sha256=digest(p)) for p in inputs})


def build_product():
    crm, occurrence, morphology = read_clouds()
    current, facts = solved_moon()
    proxy = atlas_proxy()
    unfiltered, design = photopic_weights(False), photopic_weights(True)
    configurations = dict(current_design=(current, design), current_unfiltered=(current, unfiltered),
                          atlas_unfiltered=(proxy, unfiltered))
    models = dict(current_design=dict(**facts, sunlight='titania_stack times 0.95',
                                      top_of_air_normal_lux=float(design.sum())),
                  current_unfiltered=dict(same_atmosphere_as='current_design', sunlight='unfiltered',
                                          top_of_air_normal_lux=float(unfiltered.sum())),
                  atlas_unfiltered=dict(density='Existing zero-ozone exponential optical proxy, 1-km shells',
                                        sunlight='unfiltered', top_of_air_normal_lux=float(unfiltered.sum())))
    overhead, views = [], []
    for name, (column, weights) in configurations.items():
        for height in HEIGHTS_KM:
            for depression in DEPRESSION_DEG:
                overhead.append(dict(model=name, **path_row(column, weights, height, 0., float(depression))))
            for depression in (0., 5., 10., 15.):
                for distance in RANGES_KM:
                    views.append(dict(model=name, **path_row(column, weights, height, distance, depression)))
    height_stats = crm['storms']['cloud_top_km']
    geometry, reference_levels = [], []
    for label, height in [('20km', 20.), ('40km', 40.), ('60km', 60.), ('80km', 80.), *height_stats.items()]:
        depression = sunset_depression_deg(MOON_RADIUS, height * 1000)
        geometry.append(dict(label=label, height_km=height, overhead_solar_centre_set_depression_deg=depression,
                             overhead_solar_centre_set_hours=depression * HOURS_PER_DEGREE,
                             overhead_last_limb_set_hours=(depression + SUN_RADIUS_DEG) * HOURS_PER_DEGREE,
                             same_height_earth_centre_set_minutes=
                             sunset_depression_deg(EARTH_RADIUS, height * 1000) * 24 / 360 * 60))
        reference_levels.append(dict(label=label, height_km=height, model='current_design',
                                     crossings={str(l): last_reference_level(current, design, height, l)
                                                for l in (10000, 1000, 100)}))

    # Preserve the atlas's own source/profile. Never ratio its background with
    # the new solved-column beam and label that result a cloud/sky contrast.
    atlas = PackedAtlas()
    world = atlas.data['worlds']['moon_no_ozone']
    ground = []
    for depression in np.arange(0., 25.0001, .5):
        index = world['suns'].index(-float(depression))
        diffuse = float(np.asarray(world['diffuse'][index]) @ Y_WEIGHTS)
        direct = float(np.asarray(world['direct_horizontal'][index]) @ Y_WEIGHTS)
        ground.append(dict(depression_deg=float(depression), hours=depression * HOURS_PER_DEGREE,
                           horizontal_diffuse_lux=diffuse, horizontal_direct_lux=direct,
                           horizontal_total_lux=diffuse + direct))
    # Two-path product has to be formed spectrally: multiplying two separately
    # photopic-weighted transmissions would lose their shared reddening.
    spectral_example = next(r for r in views if r['model'] == 'current_design' and r['height_km'] == 40 and
                            r['observer_sun_depression_deg'] == 5 and r['westward_surface_range_km'] == 100)
    fine_proxy = atlas_proxy(500.)
    disk_errors, shell_errors = [], []
    for height in HEIGHTS_KM:
        for depression in (0., 5., 10., 15., sunset_depression_deg(MOON_RADIUS, height * 1000)):
            s32 = current.sun_transmission(height * 1000, -depression, 32)
            s64 = current.sun_transmission(height * 1000, -depression, 64)
            disk_errors.append(abs(float((s32 - s64) @ design)))
            p1 = float(proxy.sun_transmission(height * 1000, -depression) @ unfiltered)
            p2 = float(fine_proxy.sun_transmission(height * 1000, -depression) @ unfiltered)
            if p2 >= 1:
                shell_errors.append(abs(p1 / p2 - 1))
    return dict(schema=SCHEMA, producer=provenance(),
                evidence='Repo-only mechanism screen using corrected 2-D CM1 cloud-height/occurrence products, '
                         'current solved-column molecular extinction, and separately labelled old sky-proxy context. '
                         'No cloud scene radiance, visibility probability or validated appearance is computed.',
                reading_rule='Equatorial equinox, sea-level observer at 1.6 m, smooth spherical ground. Positive '
                             'range is west toward sunset. Zero time is ground-level geometric solar-centre sunset. '
                             'Cloud direct normal lux is incident illumination on a plane normal to the beam, '
                             'not a cloud-face luminance. two_path_source_fraction integrates solar transmission '
                             'times observer-path transmission at each wavelength before photopic weighting; '
                             'it excludes the cloud scattering law and sky/path radiance. Reference lux levels '
                             'are arbitrary diagnostics, not visibility or comfort thresholds. Geometry windows '
                             'are not cloud lifetimes. Earth geometry uses the same height, not an Earth cloud forecast.',
                units=dict(length='km unless explicitly _m', angle='degrees', time='hours unless explicitly minutes',
                           illuminance='lux', transmission='dimensionless', occurrence='fraction per height/time cell'),
                assumptions=dict(synodic_month_days=SYNODIC_MONTH_DAYS, hours_per_degree=HOURS_PER_DEGREE,
                                 observer_height_m=1.6, solar_radius_deg=SUN_RADIUS_DEG, solar_disk_order=32,
                                 spectral_grid_nm=current.wavelength_nm.tolist(),
                                 solar_spectrum='Existing binned ASTM/Bruneton reference, not the WHI spectrum of surface_light.',
                                 omitted='Refraction, aerosols, other gas absorption, cloud self-shadowing/extinction, '
                                         'particle phase functions, multiple scattering by clouds, terrain, '
                                         'Earthlight, and a coupled current-column twilight sky.'),
                models=models, cloud_input=dict(case=crm['case'], gcm_run=crm['gcm']['run'],
                                                span_days=crm['span_days'], snapshots=crm['snapshots'],
                                                evidence=crm['evidence'], storms=crm['storms'],
                                                storm_tracks=crm['storm_tracks'], morphology=morphology),
                occurrence=occurrence, geometry=geometry, incident_reference_levels=reference_levels,
                overhead=overhead, views=views,
                ground_context=dict(model='Existing unfiltered moon_no_ozone atlas; albedo 0.1',
                                    rule='Different atmosphere/source from current_design. Context only; no '
                                         'quantitative current cloud-to-ground contrast can be derived. The '
                                         'atlas_unfiltered beam provides a matched-proxy comparison. Cloud light '
                                         'stored in the atlas is at 2.5 km and is not used for high clouds.',
                                    source_archive_sha256=next(p['sha256'] for p in atlas.data['provenance']
                                                              if p['file'] == 'moon_no_ozone_atlas.npz'), rows=ground),
                checks=dict(solar_disk_max_absolute_lux_32_vs_64=max(disk_errors),
                            proxy_shell_max_relative_lux_1000_vs_500m_above_1lux=max(shell_errors),
                            spectral_example=spectral_example,
                            naive_product_for_example=spectral_example['cloud_solar_transmission'] *
                            spectral_example['observer_path_source_weighted_transmission']))


def main():
    data = build_product()
    out = HERE / 'results/cloud_twilight.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    print(out)
    print(json.dumps(data['checks'], indent=2))


if __name__ == '__main__':
    main()
