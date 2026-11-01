"""Snow-only extinction and conditional light in the aerial winter.

Reads versioned snow populations and the existing solved-sky product; never
imports the illumination solver. The stored angular moments give light at the
snow's height. Directional thin-snow radiance uses a declared HG sensitivity and
an isotropic equivalent of the diffuse field. A separate conservative slab
Monte Carlo includes every scattering order and tests the thick-snow limit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from atmosphere.rain.drops import WATER_DENSITY
from shared.constants import AU, SUN_RADIUS
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.snow-optics/1'
POPULATION = HERE / 'results/population.json'
GROWTH_DEPTH = HERE / 'results/growth_depth.json'
NUMERICAL_AUDIT = HERE / 'results/numerical_audit.json'
SKY = ROOT / 'illumination/sky/results/solved_sky.json'
OUT = HERE / 'results/optics.json'
SOURCE_NOTES = HERE / 'optics_sources.json'


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def summarize_distribution(diameter_m, number_m3, projected_area_m2, *, q_ext=2.0):
    """Bin-integrated concentrations and particle areas, never density per dD.

    Area in each bin must be number weighted, so sum(n*A) retains the original
    particles' projected area. Diameter is used only to group optical shares.
    Q_ext=2 is the large-particle geometric-optics approximation, not an exact
    electromagnetic solution for the extrapolated lunar aggregates.
    """
    d, n, area = [np.asarray(x, float) for x in (diameter_m, number_m3, projected_area_m2)]
    if d.ndim != 1 or n.shape != d.shape or area.shape != d.shape:
        raise ValueError('diameter, number and area require matching one-dimensional arrays')
    if not all(np.isfinite(x).all() for x in (d, n, area)) or np.any(d <= 0) or np.any(n < 0) or np.any(area < 0):
        raise ValueError('finite positive diameters and nonnegative concentrations/areas required')
    if not math.isfinite(q_ext) or q_ext <= 0:
        raise ValueError('positive Q_ext required')
    contribution = q_ext * n * area
    beta = float(contribution.sum())
    below_one = float(contribution[d < 1.].sum())
    shares = []
    for lo, hi in ((0., .001), (.001, .01), (.01, .1), (.1, 1.), (1., 10.), (10., float('inf'))):
        selected = (d >= lo) & (d < hi)
        shares.append(dict(diameter_min_m=lo, diameter_max_m=hi if math.isfinite(hi) else None,
                           number_m3=float(n[selected].sum()),
                           extinction_share=float(contribution[selected].sum() / beta) if beta else 0.))
    return dict(number_m3=float(n.sum()), projected_area_per_volume_m_inverse=float(np.sum(n * area)),
                extinction_m_inverse=beta, e_folding_m=1 / beta if beta else None,
                beam_5_percent_distance_m=-math.log(.05) / beta if beta else None,
                q_ext_one_sensitivity=dict(extinction_m_inverse=beta/q_ext,
                    beam_5_percent_distance_m=-math.log(.05)*q_ext/beta if beta else None),
                thin_snow_max_path_m=.1 / beta if beta else None,
                transmission=[dict(path_m=x, fraction=math.exp(-beta*x)) for x in (10., 100., 1000.)],
                size_contributions=shares,
                below_one_metre_component=dict(extinction_m_inverse=below_one,
                    beam_5_percent_distance_m=-math.log(.05)/below_one if below_one else None,
                    fraction_of_total_extinction=below_one/beta if beta else 0.,
                    evidence='Only the existing sub-metre component; excluded tail mass is not redistributed. This is not a breakup solution.'))


def phase_hg(cosine, g):
    """Unit-solid-angle phase density; positive g favours forward scattering."""
    if not -1 < g < 1:
        raise ValueError('HG asymmetry must lie inside (-1, 1)')
    c = np.asarray(cosine, float)
    if np.any(np.abs(c) > 1 + 1e-12):
        raise ValueError('scattering cosine outside [-1, 1]')
    return (1-g*g) / (4*np.pi*(1+g*g-2*g*np.clip(c, -1, 1))**1.5)


def thin_radiance(direct_xyz, diffuse_scalar_xyz, scattering_angle_deg, *, tau=.05, g=.8):
    """Radiance added along a thin, homogeneous line of sight, CIE XYZ.

    Direct XYZ is normal illuminance; diffuse scalar XYZ is integral L dOmega,
    not horizontal irradiance. The diffuse angular distribution is replaced by
    its isotropic equivalent. Incident light is constant along this short path.
    Snow attenuation of that incident field and its repeated scattering are
    outside this thin calculation. tau > .1 is rejected, rather than presented
    as an optically thick storm image.
    """
    if not 0 <= tau <= .1:
        raise ValueError('thin-snow radiance requires 0 <= tau <= 0.1')
    if not 0 <= scattering_angle_deg <= 180:
        raise ValueError('scattering angle requires 0–180 degrees')
    direct, diffuse = np.asarray(direct_xyz, float), np.asarray(diffuse_scalar_xyz, float)
    if direct.shape != (3,) or diffuse.shape != (3,) or not np.isfinite(direct).all() or not np.isfinite(diffuse).all():
        raise ValueError('finite XYZ triplets required')
    if np.any(direct < 0) or np.any(diffuse < 0):
        raise ValueError('negative incident CIE XYZ')
    source = direct * phase_hg(math.cos(math.radians(scattering_angle_deg)), g) + diffuse/(4*np.pi)
    return source * -math.expm1(-tau)


def slab_flux(tau, *, g=.8, mu0=1., photons=50000, seed=0, albedo=1., max_events=100000):
    """All-order forward Monte Carlo in a horizontally infinite uniform slab.

    Unit-energy photons enter the top with downward cosine mu0. Exiting either
    surface ends a history. There is no reflecting ground, molecular air,
    surrounding cloud or finite storm edge. Fractions and binomial standard
    errors refer to integrated flux, not a directional radiance or image.
    """
    if not math.isfinite(tau) or tau < 0 or not 0 < mu0 <= 1 or not 0 <= albedo <= 1 or not -1 < g < 1:
        raise ValueError('invalid slab parameters')
    if photons < 1 or int(photons) != photons:
        raise ValueError('positive integer photon count required')
    rng = np.random.default_rng(seed)
    z, mu, events = np.zeros(photons), np.full(photons, mu0), np.zeros(photons, dtype=np.int32)
    reflected = transmitted = absorbed = unscattered = event_sum = 0
    for iteration in range(max_events):
        new_z = z - np.log1p(-rng.random(z.size)) * mu
        top, bottom = new_z < 0, new_z >= tau
        reflected += int(top.sum())
        transmitted += int(bottom.sum())
        unscattered += int(np.sum(bottom & (events == 0)))
        event_sum += int(events[top | bottom].sum())
        inside = ~(top | bottom)
        z, mu, events = new_z[inside], mu[inside], events[inside]
        if z.size == 0:
            break
        if albedo < 1:
            dead = rng.random(z.size) >= albedo
            absorbed += int(dead.sum())
            event_sum += int(events[dead].sum())
            z, mu, events = z[~dead], mu[~dead], events[~dead]
        if z.size == 0:
            break
        events += 1
        u = rng.random(z.size)
        ct = 2*u-1 if abs(g) < 1e-12 else (1+g*g-((1-g*g)/(1-g+2*g*u))**2)/(2*g)
        ct = np.clip(ct, -1, 1)
        mu = np.clip(mu*ct + np.sqrt(np.maximum(0., 1-mu*mu))*np.sqrt(1-ct*ct)*np.cos(2*np.pi*rng.random(z.size)), -1, 1)
    else:
        raise RuntimeError('unfinished photon histories; increase max_events explicitly')
    result = dict(tau=float(tau), g=float(g), incident_downward_cosine=float(mu0), single_scattering_albedo=float(albedo),
                  photons=int(photons), seed=int(seed), mean_scatter_events=event_sum/photons)
    for name, count in [('reflected', reflected), ('transmitted', transmitted), ('absorbed', absorbed), ('unscattered', unscattered)]:
        p = count/photons
        result[name] = p
        result[name+'_standard_error'] = math.sqrt(p*(1-p)/photons)
        # Wilson interval stays nonzero when a rare event has no samples.
        z95 = 1.959963984540054
        centre = (p+z95*z95/(2*photons))/(1+z95*z95/photons)
        half = z95*math.sqrt(p*(1-p)/photons+z95*z95/(4*photons**2))/(1+z95*z95/photons)
        result[name+'_interval95'] = [max(0., centre-half), min(1., centre+half)]
    result['energy_residual'] = 1-(reflected+transmitted+absorbed)/photons
    result['unscattered_exact'] = math.exp(-tau/mu0)
    return result


def shell_path_lengths(radius, direction_cosine, edges):
    """Straight half-ray lengths through concentric optical shells; None if blocked."""
    edges = np.asarray(edges, float)
    if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0) or not edges[0] <= radius <= edges[-1]:
        raise ValueError('ordered shells and an interior starting point required')
    if not -1 <= direction_cosine <= 1:
        raise ValueError('invalid direction cosine')
    impact2 = radius*radius*(1-direction_cosine**2)
    if direction_cosine < 0 and impact2 < edges[0]**2:
        return None
    root = np.sqrt(np.maximum(0., edges*edges-impact2))
    centre = -radius*direction_cosine
    lengths_inside = np.maximum(0., centre+root)-np.maximum(0., centre-root)
    return np.maximum(0., np.diff(lengths_inside))


def direct_transmission(height_m, sun_deg, edges, extinction, *, disk_strips=32):
    """Finite uniform solar disk through the product's homogeneous optical shells."""
    radius, alpha, disk = edges[0]+height_m, math.radians(sun_deg), math.asin(SUN_RADIUS/AU)
    horizon = -math.acos(edges[0]/radius)
    lower = np.clip((horizon-alpha)/disk, -1, 1)
    fraction = .5-(math.asin(lower)+lower*math.sqrt(max(0., 1-lower*lower)))/math.pi
    out = np.zeros(extinction.shape[1])
    if fraction <= 0:
        return out
    offsets = lower + (1-lower)*(np.arange(disk_strips)+.5)/disk_strips
    weights = np.sqrt(1-offsets*offsets)
    for offset, weight in zip(offsets, weights):
        lengths = shell_path_lengths(radius, math.sin(alpha+offset*disk), edges)
        if lengths is not None:
            out += weight*np.exp(-lengths@extinction)
    return out * fraction/weights.sum()


def _bilinear(rgrid, angles, values, radius, angle):
    if not rgrid[0] <= radius <= rgrid[-1] or not angles[0] <= angle <= angles[-1]:
        raise ValueError('sampling outside stored optical field')
    i = min(np.searchsorted(rgrid, radius, side='right')-1, len(rgrid)-2)
    j = min(np.searchsorted(angles, angle, side='right')-1, len(angles)-2)
    u, v = (radius-rgrid[i])/(rgrid[i+1]-rgrid[i]), (angle-angles[j])/(angles[j+1]-angles[j])
    return (1-u)*((1-v)*values[i,j]+v*values[i,j+1])+u*((1-v)*values[i+1,j]+v*values[i+1,j+1])


def read_sky_boundary(sky_path=SKY, *, heights_m=(29200., 35000., 40000., 45500.)):
    """Read saved moments at height, authenticated through the versioned atlas.

    The source solution has no independent schema; its hash is pinned by the
    atlas/1 header and its channel and moment meanings by the recorded producer.
    Only moments and small arrays are loaded; the much larger beam array is not.
    """
    source = json.loads(Path(sky_path).read_text())
    if source['schema'] != 'terluna.illumination.solved-spherical-sky/1':
        raise ValueError('unsupported solved sky schema')
    world = source['worlds']['moon_1.2atm']
    atlas = ROOT / world['atlas']['path']
    if digest(atlas) != world['atlas']['sha256']:
        raise ValueError('sky atlas hash mismatch')
    with np.load(atlas, allow_pickle=False) as saved:
        head = json.loads(str(saved['meta']))
    if head['schema'] != 'terluna.illumination.solved-spherical-sky-atlas/1':
        raise ValueError('unsupported sky atlas schema')
    field = atlas.with_name('moon_1.2atm_standard.npz')
    field_hash = digest(field)
    if field_hash != head['source_solution_sha256']:
        raise ValueError('sky moments do not match the admitted atlas')
    with np.load(field, allow_pickle=False) as saved:
        arrays = {k: saved[k].copy() for k in ('r', 'a', 'edges', 'scattering', 'absorption', 'xyz', 'moments', 'wavelength')}
        settings = json.loads(str(saved['meta']))
    if not settings['order_converged'] or arrays['moments'].shape[-1] != 7:
        raise ValueError('a converged seven-moment field is required')
    moment = arrays['moments']
    if not np.allclose(moment[...,0], moment[...,1:4].sum(-1), rtol=1e-8, atol=1e-12):
        raise ValueError('diffuse scalar irradiance differs from moment trace')
    # The evening product pins both completed spatial grids to this same sky.
    # Its coarse field provides a numerical sensitivity at the actual sampled
    # winter heights, beyond the surface checks in the original sky product.
    context_path = ROOT/'illumination/cloud_light/results/evening_context.json'
    context = json.loads(context_path.read_text())
    if context['schema'] != 'terluna.illumination.evening-context/1' or context['source']['sha256'] != digest(sky_path):
        raise ValueError('coarse-field manifest does not bind this sky product')
    coarse_path = field.with_name('moon_1.2atm_draft.npz')
    coarse_key = str(coarse_path.relative_to(ROOT))
    coarse_hash = digest(coarse_path)
    if coarse_hash != context['producer']['inputs'][coarse_key]:
        raise ValueError('coarse moment field hash mismatch')
    with np.load(coarse_path, allow_pickle=False) as saved:
        coarse = {k: saved[k].copy() for k in ('r', 'a', 'xyz', 'moments')}
        if not json.loads(str(saved['meta']))['order_converged']:
            raise ValueError('coarse moment field did not converge')
    rows = []
    coarse_controls = []
    for height in heights_m:
        for sun in (60., 10., 0., -5., -10., -20., -40.):
            direct = direct_transmission(height, sun, arrays['edges'], arrays['scattering']+arrays['absorption'])@arrays['xyz']
            moments = _bilinear(arrays['r'], arrays['a'], moment, arrays['edges'][0]+height, math.radians(sun))
            xyz = moments.T@arrays['xyz']
            c_moments = _bilinear(coarse['r'], coarse['a'], coarse['moments'], arrays['edges'][0]+height, math.radians(sun))
            c_xyz = c_moments.T@coarse['xyz']
            coarse_controls.append(dict(altitude_m=height, source_elevation_deg=sun,
                diffuse_scalar_relative_difference=float(c_xyz[0,1]/xyz[0,1]-1),
                diffuse_down_relative_difference=float(c_xyz[5,1]/xyz[5,1]-1)))
            rows.append(dict(source='shielded_sun', altitude_m=height, source_elevation_deg=sun,
                             direct_normal_xyz=direct.tolist(), diffuse_scalar_xyz=xyz[0].tolist(),
                             diffuse_down_horizontal_xyz=xyz[5].tolist(), diffuse_up_horizontal_xyz=xyz[6].tolist()))
        # At 25 degrees above the horizon the source-size difference is small.
        # This retains solar spectral weighting; it does not claim Earth's blue spectrum.
        earth_sun = 25.
        earth_above_air_lux = 7.8  # admitted rounded full-Earth value, earthlight/README.md
        scale = earth_above_air_lux/float(arrays['xyz'][:,1].sum())
        direct = direct_transmission(height, earth_sun, arrays['edges'], arrays['scattering']+arrays['absorption'])@arrays['xyz']*scale
        moments = _bilinear(arrays['r'], arrays['a'], moment, arrays['edges'][0]+height, math.radians(earth_sun))
        xyz = moments.T@arrays['xyz']*scale
        rows.append(dict(source='full_earth_solar_colour_proxy', altitude_m=height, source_elevation_deg=earth_sun,
                         above_air_normal_lux=earth_above_air_lux, direct_normal_xyz=direct.tolist(),
                         diffuse_scalar_xyz=xyz[0].tolist(), diffuse_down_horizontal_xyz=xyz[5].tolist(),
                         diffuse_up_horizontal_xyz=xyz[6].tolist()))
    # Match several published surface direct-beam values: independent shell readout.
    checks = []
    for sun in (0., 10., 60., 90.):
        expected = next(r['direct_normal_lux'] for r in world['samples'] if r['sun_deg'] == sun)
        actual = float((direct_transmission(0., sun, arrays['edges'], arrays['scattering']+arrays['absorption'])@arrays['xyz'])[1])
        checks.append(dict(sun_deg=sun, stored_lux=expected, independent_lux=actual, relative_difference=actual/expected-1))
    if max(abs(r['relative_difference']) for r in checks) > .005:
        raise ValueError('independent optical-shell readout failed the stored surface control')
    disk_controls = []
    for height, sun in ((29200., -10.), (35000., -10.), (45500., -12.), (35000., 0.)):
        normal = direct_transmission(height, sun, arrays['edges'], arrays['scattering']+arrays['absorption'])@arrays['xyz']
        refined = direct_transmission(height, sun, arrays['edges'], arrays['scattering']+arrays['absorption'], disk_strips=128)@arrays['xyz']
        disk_controls.append(dict(altitude_m=height, source_elevation_deg=sun, disk_strips=[32, 128],
                                  direct_normal_relative_difference=float(normal[1]/refined[1]-1)))
    return dict(rows=rows, input_files={str(Path(sky_path).relative_to(ROOT)): digest(sky_path),
                    str(atlas.relative_to(ROOT)): digest(atlas), str(field.relative_to(ROOT)): field_hash,
                    str(context_path.relative_to(ROOT)): digest(context_path), coarse_key: coarse_hash,
                    'illumination/earthlight/README.md': digest(ROOT/'illumination/earthlight/README.md')},
                column=world['optical_column']['column'], solar_shield='titania-stack transmission × 0.95',
                optical_key=world['optical_column']['key'], shield_sha256=world['optical_column']['shield_sha256'],
                surface_direct_controls=checks, coarse_moment_controls=coarse_controls, finite_disk_controls=disk_controls,
                reading_rule='Bilinear interpolation in physical radius and geometric solar elevation of the saved spectral moments. '
                             'Moment 0 is integral L dOmega; moments 5 and 6 are downward/upward horizontal flux. '
                             'Direct light uses independent analytic shell intersections and a finite uniform Sun. '
                             'XYZ weights already contain the shielded TSIS-1 source and 360–830 nm photopic response.',
                limits=['Molecular clear sky at a 294.9 K surface; snow population column uses the 301 K resort interpolation.',
                        'Snow and supercooled cloud do not feed back on this incident clear-sky field.',
                        'Earth proxy uses the admitted 7.8 lux full-Earth normalization with shielded-solar colour; '
                        'Earth spectral shape, phase, rotating continents and its larger finite disk are not modeled; '
                        'the proxy is a separate source contribution, not a combined sky at a dated orbital configuration.'])


def droplet_extinction(liquid_water_kg_m3, effective_radius_m):
    if liquid_water_kg_m3 < 0 or effective_radius_m <= 0:
        raise ValueError('nonnegative liquid water and positive radius required')
    return 3*liquid_water_kg_m3/(2*WATER_DENSITY*effective_radius_m)


def build(population_path=POPULATION, *, growth_depth_path=GROWTH_DEPTH, numerical_audit_path=NUMERICAL_AUDIT, photons=50000):
    population = json.loads(Path(population_path).read_text())
    if population['schema'] != 'terluna.atmosphere.snow-population/1':
        raise ValueError('unsupported snow population schema')
    growth_depth = json.loads(Path(growth_depth_path).read_text())
    if growth_depth['schema'] != 'terluna.atmosphere.snow-growth-depth/1':
        raise ValueError('unsupported snow growth-depth schema')
    numerical_audit = json.loads(Path(numerical_audit_path).read_text())
    if numerical_audit['schema'] != 'terluna.atmosphere.snow-numerical-audit/1':
        raise ValueError('unsupported snow numerical-audit schema')
    audit_identity = dict(path=str(Path(numerical_audit_path).relative_to(ROOT)), sha256=digest(numerical_audit_path))
    for source in (population, growth_depth):
        assessment = source.get('numerical_assessment', {})
        if (assessment.get('schema') != 'terluna.atmosphere.snow-assessment/1'
                or assessment.get('audit') != audit_identity):
            raise ValueError('Snow population requires a current assessment bound to the supplied numerical audit')
    runs = [(run, 'population') for run in population['runs']] + [(run, 'growth_depth') for run in growth_depth['runs']]
    if len({run['id'] for run, _ in runs}) != len(runs):
        raise ValueError('population and growth-depth run identifiers must be unique')
    rows, columns = [], []
    for run, source_product in runs:
        layer_tau, below_one_layer_tau = 0., 0.
        for level_name, distribution in [('freezing_level', run['freezing_level'])] + [
                (f'layer_{i}', layer) for i, layer in enumerate(run.get('layers', []))]:
            optics = summarize_distribution(distribution['diameter_m'], distribution['number_m3'], distribution['projected_area_m2'])
            ice_mass = float(np.sum(distribution['mass_kg_m3']))
            optics['ice_water_kg_m3'] = ice_mass
            sizes, masses = np.asarray(distribution['diameter_m']), np.asarray(distribution['mass_kg_m3'])
            for contribution in optics['size_contributions']:
                lo, hi = contribution['diameter_min_m'], contribution['diameter_max_m']
                selected = (sizes >= lo) & (sizes < hi if hi is not None else True)
                contribution['ice_water_kg_m3'] = float(masses[selected].sum())
                contribution['ice_mass_share'] = float(masses[selected].sum()/ice_mass) if ice_mass else 0.
            optics['mass_specific_extinction_m2_kg'] = optics['extinction_m_inverse']/ice_mass if ice_mass else None
            rows.append(dict(run_id=run['id'], source_product=source_product, cloud_top_c=run['top_c'],
                             world=run['world'], geometry=run['geometry'],
                             snowfall_set_mm_h=run['snowfall_set_mm_h'], snowfall_arriving_mm_h=run['summary']['snowfall_mm_h'],
                             crystals_per_litre=run['crystals_per_litre'],
                             seed=run['seed'], level=level_name, altitude_m=distribution.get('altitude_m'), **optics))
            if level_name != 'freezing_level':
                layer_tau += optics['extinction_m_inverse']*distribution['depth_m']
                below_one_layer_tau += optics['below_one_metre_component']['extinction_m_inverse']*distribution['depth_m']
        columns.append(dict(run_id=run['id'], source_product=source_product, cloud_top_c=run['top_c'], snow_vertical_optical_depth=layer_tau,
                            vertical_unscattered_transmission=math.exp(-layer_tau),
                            below_one_metre_component_optical_depth=below_one_layer_tau,
                            evidence='Integral of time-mean snow projected area only; excludes cloud liquid and molecular air.'))
    light = read_sky_boundary()
    thin = []
    for row in light['rows']:
        for g in (.7, .8, .9):
            directions = []
            for theta in (20., 90., 160.):
                xyz = thin_radiance(row['direct_normal_xyz'], row['diffuse_scalar_xyz'], theta, tau=.05, g=g)
                # Explicit comparison scene: a black target and a background
                # whose unobstructed luminance is the diffuse all-sky mean.
                background = row['diffuse_scalar_xyz'][1]/(4*np.pi)
                transmitted_background = background*math.exp(-.05)
                contrast = -transmitted_background/(transmitted_background+xyz[1]) if transmitted_background+xyz[1] else None
                directions.append(dict(scattering_angle_deg=theta, luminance_cd_m2=float(xyz[1]),
                                       chromaticity_xy=(xyz[:2]/xyz.sum()).tolist() if xyz.sum() else None,
                                       imposed_background_luminance_cd_m2=background,
                                       black_target_contrast=contrast))
            thin.append(dict(source=row['source'], altitude_m=row['altitude_m'], source_elevation_deg=row['source_elevation_deg'],
                             g=g, sightline_optical_depth=.05, unscattered_background_fraction=math.exp(-.05), directions=directions))
    slabs = [slab_flux(tau, g=g, mu0=mu, photons=photons, seed=2180+i*100+j*10+k)
             for i, tau in enumerate((.1, 1., 10., 100.)) for j, g in enumerate((.7, .8, .9)) for k, mu in enumerate((1., .2))]
    droplets = []
    for water in (.01, .1, .5):
        beta = droplet_extinction(water*1e-3, 10e-6)
        droplets.append(dict(liquid_water_g_m3=water, effective_radius_um=10., extinction_m_inverse=beta,
                             beam_5_percent_distance_m=-math.log(.05)/beta))
    model_files = ['atmosphere/snow/optics.py', 'atmosphere/rain/drops.py']
    product = dict(schema=SCHEMA, producer=dict(lane='atmosphere', runner='atmosphere/snow/optics.py',
                  files={p:digest(ROOT/p) for p in model_files}, constants=constants_used(model_files)),
                  inputs={str(Path(population_path).relative_to(ROOT)):digest(population_path),
                          str(Path(growth_depth_path).relative_to(ROOT)):digest(growth_depth_path),
                          str(Path(numerical_audit_path).relative_to(ROOT)):digest(numerical_audit_path),
                          str(SOURCE_NOTES.relative_to(ROOT)):digest(SOURCE_NOTES), **light['input_files']},
                  evidence='Optical diagnostics on exploratory, unvalidated inherited snow populations, with representation-dependent '
                           'particle arithmetic and uncapped extreme tails. Snow-only geometric-optics extinction retains explicit number/area moments; '
                           'height-resolved illumination read from an existing molecular-sky solution; thin-path HG radiance '
                           'and independent all-order slab flux sensitivities. No lunar observation or resolved storm image.',
                  reading_rule='Use a run and level together. Number_m3 is bin integrated; area is its number-weighted particle mean. '
                               'Q_ext=2 is wavelength independent. The 5% beam distance is a snow-only optical range, '
                               'not a prediction of human visibility under twilight or Earthlight. '
                               'Extinction is evaluated from the time-mean population as a homogeneous medium: '
                               'exp(-mean optical depth) is not the time mean of transmission through intermittent snow. '
                               'Thin directional outputs apply at tau=.05 in exposed snow with approximately clear incident illumination '
                               'and isotropized diffuse light; a short sightline deep inside an opaque cloud does not meet that condition. '
                               'slab fluxes use a distinct infinite-slab geometry and cannot be substituted for those radiances '
                               'or used alone to infer image detail or a visible horizon.',
                  assumptions=dict(q_ext=2., single_scattering_albedo=1., hg_asymmetry_sensitivity=[.7,.8,.9],
                                   size_parameter_regime='dimensions and monomers much larger than visible wavelengths',
                                   particle_area='population settling area used for every ray direction; no spherical envelope substituted or optical orientation distribution solved'),
                  population_numerical_audit=dict(path=str(Path(numerical_audit_path).relative_to(ROOT)),
                      schema=numerical_audit['schema'], evidence=numerical_audit['evidence'], reading_rule=numerical_audit['reading_rule']),
                  distributions=rows, column_snow_optical_depths=columns, lighting=light, thin_path_radiance=thin, multiple_scattering_slabs=slabs,
                  cloud_droplet_sensitivity=droplets,
                  limits=['Aggregate morphology, orientation, transparency and breakup remain population/optical assumptions. '
                          'The inherited uncapped aggregation law can generate nonphysical extreme tails without an upper size bound. '
                          'Their area remains in the total and is exposed separately in 1–10 m and >=10 m size contributions; '
                          'optical totals with substantial such contributions are runaway-model diagnostics, not credible storm predictions. '
                          'Sub-metre component values still inherit the population arithmetic and sampling limitations.',
                          'HG has no ice-halo peaks or specular facets; neither halo nor sparkle appearance is predicted.',
                          'The large-particle extinction paradox includes narrow forward diffraction. Instrument aperture and '
                          'finite visual targets may recover some scattered light; no human observer model is included. '
                          'The Q_ext=1 sensitivity counts geometric interception only; it is not a confidence bound or '
                          'a substitute for a phase matrix plus a specified receiver acceptance angle.',
                          'Small supercooled cloud droplets can dominate visibility; water content is not supplied by the snow population.',
                          'Conservative grey snow neglects weak visible ice absorption and impurities; no snow-induced spectral colour shift is solved.',
                          'No joint three-dimensional snow/cloud/molecular transport or calibrated crystal phase matrix.'])
    OUT.write_text(json.dumps(product, indent=2, allow_nan=False)+'\n')
    return product


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--population', type=Path, default=POPULATION)
    parser.add_argument('--growth-depth', type=Path, default=GROWTH_DEPTH)
    parser.add_argument('--photons', type=int, default=50000)
    args = parser.parse_args()
    result = build(args.population, growth_depth_path=args.growth_depth, photons=args.photons)
    print(OUT.relative_to(ROOT), len(result['distributions']), 'population levels;', len(result['lighting']['rows']), 'lighting cases')
