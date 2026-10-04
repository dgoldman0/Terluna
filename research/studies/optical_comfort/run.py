"""Optical comfort screening from the committed surface-light product.

Run from the repository root with python -m research.studies.optical_comfort.run.
The sky's angular shapes, surface geometry and observer are explicit scenarios.
No outdoor discomfort probability or safety threshold is inferred from lux.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import AU, SUN_RADIUS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
INPUT = ROOT / 'illumination/surface_light/results/surface_light.json'
SCHEMA = 'terluna.research.optical-comfort/1'
INPUT_SCHEMA = 'terluna.illumination.surface-light/1'
CASES = {'moon': 'moon_1.2atm/design', 'earth': 'earth_control/unfiltered',
         'moon_1atm': 'moon_1.0atm/design'}
# Angular sensitivity shapes L(mu) proportional to a + b*mu, mu = sin(elevation).
# These are deliberately simple shapes, not predictions or uncertainty bounds.
SKY_SHAPES = {'uniform': (1.0, 0.0), 'zenith_bright': (1.0, 2.0),
              'horizon_bright': (3.0, -2.0)}
PATCH_REFLECTANCES = (0.1, 0.2, 0.5, 0.8)
AGES = (30.0, 50.0, 70.0)
GAZES_DEG = (-90.0, -60.0, -30.0, 0.0, 30.0)
WATER_REFRACTIVE_INDEX = 1.333  # prescribed clear-water, wavelength-independent scenario


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_input(path=INPUT):
    data = json.loads(path.read_text())
    if data.get('schema') != INPUT_SCHEMA:
        raise ValueError(f'{path}: expected {INPUT_SCHEMA}, got {data.get("schema")}')
    return data


def sky_vertical_factor(a, b):
    """Vertical sky illuminance / horizontal diffuse illuminance, axisymmetric sky.

    Horizontal integral = 2*pi*C*(a/2+b/3); vertical = C*(a*pi/2+2*b/3).
    """
    if min(a, a + b) < 0 or a / 2 + b / 3 <= 0:
        raise ValueError('Sky must be nonnegative with positive horizontal flux')
    return (a * math.pi / 2 + 2 * b / 3) / (2 * math.pi * (a / 2 + b / 3))


def eye_illuminance(direct_h, diffuse_h, albedo, sun_deg, gaze_deg=0.0,
                    azimuth_from_sun_deg=90.0):
    """Unoccluded eye plane normal to gaze, isotropic sky and flat Lambertian ground.

    A level gaze gives vertical illuminance. At other gaze elevations the plane
    tilts. The Sun is a collimated beam in this plane-projection calculation.
    This is incident photopic illuminance, before pupil/eyelid optics.
    """
    if not all(math.isfinite(v) for v in (direct_h, diffuse_h, albedo, sun_deg, gaze_deg, azimuth_from_sun_deg)):
        raise ValueError('Illuminances and geometry must be finite')
    if direct_h < 0 or diffuse_h < 0:
        raise ValueError('Illuminances must be nonnegative')
    if not (0 < sun_deg <= 90 and -90 <= gaze_deg <= 90 and 0 <= albedo <= 1):
        raise ValueError('Invalid solar elevation, gaze or albedo')
    h, e, az = map(math.radians, (sun_deg, gaze_deg, azimuth_from_sun_deg))
    dot = math.sin(e) * math.sin(h) + math.cos(e) * math.cos(h) * math.cos(az)
    direct = direct_h / math.sin(h) * max(0.0, dot)
    sky = diffuse_h * (1 + math.sin(e)) / 2
    ground = albedo * (direct_h + diffuse_h) * (1 - math.sin(e)) / 2
    return {'direct_lux': direct, 'sky_lux': sky, 'ground_lux': ground,
            'total_lux': direct + sky + ground}


def surface_rows(data):
    rows = []
    for world, case in CASES.items():
        for height, summaries in data['cases'][case]['summaries'].items():
            h = float(height)
            for albedo, s in summaries.items():
                a = float(albedo)
                total, direct = s['illuminance_lux'], s['direct_illuminance_lux']
                diffuse = total - direct
                if min(total, direct, diffuse) < 0:
                    raise ValueError('Nonphysical stored illuminance')
                ambient = (diffuse + a * total) / 2
                toward = eye_illuminance(direct, diffuse, a, h, azimuth_from_sun_deg=0)
                rows.append(dict(
                    world=world, sun_deg=h, landscape_albedo=a,
                    horizontal_total_lux=total, horizontal_direct_lux=direct,
                    horizontal_diffuse_lux=diffuse,
                    ground_luminance_cd_m2=a * total / math.pi,
                    uniform_sky_luminance_cd_m2=diffuse / math.pi,
                    vertical_ambient_lux=ambient,
                    vertical_sunward_lux=toward['total_lux'],
                    vertical_sun_direct_lux=toward['direct_lux'],
                    vertical_sky_lux=diffuse / 2, vertical_ground_lux=a * total / 2,
                    vertical_zenith_bright_lux=sky_vertical_factor(1, 2) * diffuse + a * total / 2,
                    vertical_horizon_bright_lux=sky_vertical_factor(3, -2) * diffuse + a * total / 2,
                    disk_only_shadow_fraction=diffuse / total,
                    sun_shadow_michelson=direct / (total + diffuse),
                    angular_model='uniform_sky_flat_ground_with_named_shape_sensitivities',
                    transfer_quality='two_stream_low_sun_caution' if h <= 20 else 'two_stream_screen',
                ))
    return rows


def gaze_scenes(rows):
    """Eye-plane orientation sensitivity; eye illuminance is not a glare score.

    Downward views weight the dimmer lunar ground more strongly. This can
    reverse the Moon/Earth ordering seen in a level gaze under the same sky.
    """
    out = []
    for r in rows:
        if r['world'] not in ('moon', 'earth') or r['sun_deg'] not in (90, 30) or r['landscape_albedo'] not in (0.1, 0.8):
            continue
        for gaze in GAZES_DEG:
            for azimuth in (0.0, 180.0):
                light = eye_illuminance(r['horizontal_direct_lux'], r['horizontal_diffuse_lux'],
                                        r['landscape_albedo'], r['sun_deg'], gaze, azimuth)
                out.append(dict(world=r['world'], sun_deg=r['sun_deg'],
                                landscape_albedo=r['landscape_albedo'], gaze_deg=gaze,
                                azimuth_from_sun_deg=azimuth, angular_model='uniform_sky_flat_ground',
                                ambient_lux=light['sky_lux'] + light['ground_lux'], **light))
    return out


def local_surfaces(rows):
    """Small matte patches do not reset the landscape albedo or its sky feedback."""
    out = []
    for row in rows:
        if row['sun_deg'] not in (90.0, 30.0) or row['landscape_albedo'] not in (0.1, 0.8):
            continue
        for rho in PATCH_REFLECTANCES:
            out.append(dict(world=row['world'], sun_deg=row['sun_deg'],
                            landscape_albedo=row['landscape_albedo'], patch_reflectance=rho,
                            horizontal_cd_m2=rho * row['horizontal_total_lux'] / math.pi,
                            vertical_sunward_cd_m2=rho * row['vertical_sunward_lux'] / math.pi,
                            vertical_away_cd_m2=rho * row['vertical_ambient_lux'] / math.pi,
                            horizontal_disk_shadow_cd_m2=rho * row['horizontal_diffuse_lux'] / math.pi))
    return out


def sky_fraction_below(elevation_deg):
    """Share of uniform sky's vertical illuminance below an elevation cutoff."""
    if not 0 <= elevation_deg <= 90:
        raise ValueError('Cutoff must lie in [0, 90]')
    e = math.radians(elevation_deg)
    return 2 * (e + math.sin(e) * math.cos(e)) / math.pi


def masks(rows):
    """Component removal scenarios, with unchanged distant sunlit ground.

    All masks remove the solar disk. They are angular masks, not built awnings:
    reflected roof light and changes to nearby ground illumination are omitted.
    """
    settings = [('disk_only', 1.0, 1.0),
                ('disk_and_sky_above_45deg', sky_fraction_below(45), 1.0),
                ('disk_and_sky_above_15deg', sky_fraction_below(15), 1.0),
                ('disk_and_all_sky', 0.0, 1.0),
                ('disk_and_90pct_each_hemisphere', 0.1, 0.1)]
    out = []
    for r in rows:
        if r['world'] == 'moon_1atm' or r['sun_deg'] not in (90, 30) or r['landscape_albedo'] not in (0.1, 0.8):
            continue
        for name, sky_share, ground_share in settings:
            value = sky_share * r['vertical_sky_lux'] + ground_share * r['vertical_ground_lux']
            out.append(dict(world=r['world'], sun_deg=r['sun_deg'],
                            landscape_albedo=r['landscape_albedo'], mask=name,
                            retained_sky_share=sky_share, retained_ground_share=ground_share,
                            eye_lux=value, fraction_of_open_sunward=value / r['vertical_sunward_lux']))
    return out


def veiling_point(eye_lux, separation_deg, age_years):
    """CIE age-adjusted Stiles-Holladay (Vos 2003), strictly 1 < theta < 30 deg."""
    theta = np.asarray(separation_deg, dtype=float)
    if np.any((theta <= 1) | (theta >= 30)) or eye_lux < 0 or age_years < 0:
        raise ValueError('CIE point screen requires nonnegative inputs and 1 < theta < 30 deg')
    return 10 * (1 + (age_years / 70) ** 4) * eye_lux / theta ** 2


def contrast_retained(background_cd_m2, veil_cd_m2):
    """Weber contrast after an additive uniform veil over target and background."""
    if background_cd_m2 <= 0 or veil_cd_m2 < 0:
        raise ValueError('Positive background and nonnegative veil required')
    return background_cd_m2 / (background_cd_m2 + veil_cd_m2)


def annular_veil_per_luminance(age_years, inner_deg=5.0, outer_deg=25.0, order=64):
    """Integrate a uniform excess-luminance annulus around fixation.

    dE = delta_L*cos(theta)*2*pi*sin(theta)*dtheta. Integration uses radians;
    the CIE kernel uses degrees. Baseline uniform background is subtracted.
    """
    if not 1 < inner_deg < outer_deg < 30:
        raise ValueError('Entire annulus must lie inside the CIE angular domain')
    x, w = np.polynomial.legendre.leggauss(order)
    low, high = map(math.radians, (inner_deg, outer_deg))
    theta = (x + 1) * (high - low) / 2 + low
    kernel = veiling_point(1.0, np.degrees(theta), age_years)
    return float(np.sum(w * kernel * 2 * math.pi * np.sin(theta) * np.cos(theta)) * (high - low) / 2)


def glare_scenarios(rows):
    solar, surface = [], []
    for r in rows:
        if r['world'] not in ('moon', 'earth') or r['sun_deg'] not in (30, 10) or r['landscape_albedo'] != 0.1:
            continue
        dni = r['horizontal_direct_lux'] / math.sin(math.radians(r['sun_deg']))
        for age in AGES:
            for theta in (5.0, 10.0, 20.0):
                eye = dni * math.cos(math.radians(theta))
                veil = float(veiling_point(eye, theta, age))
                solar.append(dict(world=r['world'], sun_deg=r['sun_deg'], age_years=age,
                                  separation_deg=theta, source_eye_lux=eye, veil_cd_m2=veil,
                                  contrast_retained_at_background={str(int(b)): contrast_retained(b, veil)
                                                                  for b in (100, 1000, 10000)}))
    # An ideal dark matte task in a bright annular surround, both lit equally.
    # Ratios match on both worlds; illumination amplitude and neutral tint cancel.
    for age in AGES:
        k = annular_veil_per_luminance(age)
        for ratio in (2.0, 5.0, 8.0):
            surface.append(dict(age_years=age, surround_to_task_luminance=ratio,
                                inner_deg=5.0, outer_deg=25.0,
                                excess_veil_over_background=k * (ratio - 1),
                                contrast_retained=contrast_retained(1, k * (ratio - 1))))
    return dict(direct_sun=solar, bright_surface_annulus=surface)


def fresnel_water(cos_incidence, n=WATER_REFRACTIVE_INDEX):
    """Unpolarized air-to-water dielectric Fresnel reflectance; n_air = 1."""
    if not 0 <= cos_incidence <= 1 or n <= 1:
        raise ValueError('Expected nonnegative incidence cosine and water index > 1')
    ct = math.sqrt(1 - (1 - cos_incidence ** 2) / n ** 2)
    rs = (cos_incidence - n * ct) / (cos_incidence + n * ct)
    rp = (n * cos_incidence - ct) / (n * cos_incidence + ct)
    return (rs * rs + rp * rp) / 2


def water_glint(rows):
    alpha = math.asin(SUN_RADIUS / AU)
    omega = 2 * math.pi * (1 - math.cos(alpha))
    projected_omega = math.pi * math.sin(alpha) ** 2
    out = []
    for r in rows:
        if r['world'] not in ('moon', 'earth') or r['landscape_albedo'] != 0.1:
            continue
        h = math.radians(r['sun_deg'])
        dni = r['horizontal_direct_lux'] / math.sin(h)
        reflectance = fresnel_water(math.sin(h))
        solar_luminance = dni / projected_omega
        out.append(dict(world=r['world'], sun_deg=r['sun_deg'],
                        reflectance_unpolarized=reflectance,
                        mean_solar_disk_luminance_cd_m2=solar_luminance,
                        resolved_glint_luminance_cd_m2=reflectance * solar_luminance))
    return dict(refractive_index=WATER_REFRACTIVE_INDEX, solar_angular_radius_deg=math.degrees(alpha),
                solar_solid_angle_sr=omega, solar_projected_solid_angle_sr=projected_omega,
                evidence='Resolved reflection of a uniform solar disk in a horizontal water patch, '
                         'observed in the specular direction at a downward gaze equal to the solar elevation. '
                         'Pointwise incident Fresnel factor at the disk centre; no rough-water extent, '
                         'wave statistics, sky reflection, polarization transport, or observer-path extinction.', rows=out)


def low_sun_proxy(data):
    """Transfer stored unfiltered atlas component ratios to the design case, A=0.1.

    This is a sensitivity, not a new shielded spherical solve. No high-albedo
    extension or spectral/angular radiance is invented from the flux ratios.
    """
    out = []
    checks = data['checks']['atlas_comparison']
    for world, unfiltered in [('moon', 'moon_1.2atm/unfiltered'), ('earth', CASES['earth'])]:
        for check in checks[unfiltered]['rows']:
            h = check['sun_deg']
            if h not in (20, 10, 7, 5, 3, 2, 1):
                continue
            if check['atlas_sun_deg'] != h:
                raise ValueError('Proxy requires identical solar elevation')
            s = data['cases'][CASES[world]]['summaries'][str(h)]['0.1']
            total, direct = s['illuminance_lux'], s['direct_illuminance_lux']
            b_ratio = check['atlas_direct_lux'] / check['direct_lux']
            d_ratio = check['atlas_diffuse_lux'] / check['diffuse_lux']
            b = direct * b_ratio
            d = (total - direct) * d_ratio
            out.append(dict(world=world, sun_deg=h, landscape_albedo=0.1,
                            source_atlas_hash=checks[unfiltered]['atlas_sha256'],
                            direct_component_ratio=b_ratio, diffuse_component_ratio=d_ratio,
                            two_stream_total_lux=total, proxy_total_lux=b + d,
                            proxy_diffuse_lux=d, proxy_direct_lux=b,
                            proxy_vertical_ambient_lux=(d + 0.1 * (b + d)) / 2))
    return out


def interval_hours(low_deg, high_deg, latitude_deg=0.0, period_hours=SYNODIC_MONTH_DAYS * 24):
    """One rising OR setting passage, declination zero, flat geometric horizon."""
    if not 0 <= low_deg < high_deg <= 90 or not 0 <= abs(latitude_deg) < 90:
        raise ValueError('Invalid solar interval or latitude')
    c = math.cos(math.radians(latitude_deg))
    hour_angle = lambda h: math.acos(min(1.0, math.sin(math.radians(h)) / c))
    return (hour_angle(low_deg) - hour_angle(high_deg)) * period_hours / (2 * math.pi)


def build(data, input_path=INPUT):
    rows = surface_rows(data)
    sources = [Path(__file__), HERE / 'sources.json']
    return dict(
        schema=SCHEMA,
        producer=dict(lane='research', runner='research/studies/optical_comfort/run.py',
                      files={str(p.relative_to(ROOT)): digest(p) for p in sources},
                      constants=constants_used(sources),
                      inputs={str(input_path.relative_to(ROOT)) if input_path.is_relative_to(ROOT) else str(input_path):
                              dict(schema=data['schema'], sha256=digest(input_path))}),
        evidence='Conditional clear-sky screening from committed illuminance summaries. '
                 'Angular skies, local matte surfaces, dark masks, glare tasks and smooth water are imposed scenarios. '
                 'No validated outdoor discomfort probability, resolved landscape radiance, clinical prediction '
                 'or eye-safety assessment. Low-Sun atlas transfer is labelled separately.',
        reading_rule='Rows are sampled solar elevations and spectrally neutral landscape albedos, without interpolation. '
                     'vertical_ambient excludes the solar disk; vertical_sunward adds its projected flux. '
                     'Local patch reflectance never changes the landscape feedback. Luminance uses photopic weighting '
                     'already applied in the input. gaze_scenes use the eye plane normal to the stated gaze; '
                     'zero is horizontal and negative angles look down. Angular sensitivities are scenarios, '
                     'not statistical bounds. Source and input hashes identify the calculation independently of Git HEAD.',
        units=dict(illuminance='lux at the named surface or eye plane', luminance='cd m^-2',
                   angles='degrees except quadrature internals', time='hours', reflectance='dimensionless'),
        assumptions=dict(sky_shapes=SKY_SHAPES, sky_vertical_factors={k: sky_vertical_factor(*v) for k, v in SKY_SHAPES.items()},
                         ground='unobstructed flat Lambertian landscape, albedo applied to the entire upstream column',
                         patch_reflectances=PATCH_REFLECTANCES, observer_ages_years=AGES, gazes_deg=GAZES_DEG,
                         plane_projection='Collimated Sun; the finite disk is used only in the water-radiance calculation',
                         source_model='CIE age-adjusted Stiles-Holladay, 1 < theta < 30 deg; ages are sensitivity cases; '
                                      'photopic achromatic screen without wavelength-dependent ocular scatter'),
        surface_scenes=rows, gaze_scenes=gaze_scenes(rows),
        local_matte_surfaces=local_surfaces(rows), angular_masks=masks(rows),
        disability_glare=glare_scenarios(rows), smooth_water=water_glint(rows), low_sun_atlas_proxy=low_sun_proxy(data),
        solar_residence=[dict(latitude_deg=lat, solar_interval_deg=[0, upper],
                              moon_one_pass_hours=interval_hours(0, upper, lat),
                              earth_one_pass_hours=interval_hours(0, upper, lat, 24))
                         for lat in (0, 30, 60) for upper in (10, 30)],
        checks=dict(annular_quadrature_max_abs_64_vs_128=max(
                        abs(annular_veil_per_luminance(age, order=64) - annular_veil_per_luminance(age, order=128))
                        for age in AGES)),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, default=HERE / 'results')
    args = parser.parse_args()
    data = build(read_input(args.input), args.input.resolve())
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'optical_comfort.json').write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    with (args.output / 'surface_scenes.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data['surface_scenes'][0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(data['surface_scenes'])
    print(f'Wrote {len(data["surface_scenes"])} scenes to {args.output}')


if __name__ == '__main__':
    main()
