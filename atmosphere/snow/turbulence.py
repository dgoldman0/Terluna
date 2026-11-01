"""Snow-conditioned motions in archived CM1 fields and independent flake-scale sensitivity.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.snow.turbulence

Reads existing climate output as data; starts no model runs.  Coarse resolved
motions are not a measurement of centimetre-scale turbulence.  The epsilon
cases are uncalibrated sensitivities anchored to terrestrial cloud observations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from atmosphere.snow import flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / 'climate/crm/runs'
OUT = HERE / 'results/turbulence.json'
SCHEMA = 'terluna.atmosphere.snow-turbulence/1'
FIELDS = ('th', 'prs', 'qv', 'qc', 'qi', 'qs', 'uinterp', 'vinterp', 'winterp')
CODE = ('atmosphere/snow/turbulence.py', 'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py',
        'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py',
        'atmosphere/snow/habits.py', 'atmosphere/snow/habits_sources.json')
CASES = ('box_0e_small', 'box_0e_elec_corrected', 'box_0e_elec_corrected_storms_second')
WINTER_TEMPERATURE_C = (-20.0, 0.0)
HEIGHT_BAND_M = (20000.0, 60000.0)
SNOW_THRESHOLDS = (1e-6, 1e-5, 1e-4)  # kg snow / kg dry air; sensitivity, not universal snow thresholds
EPSILONS = (1e-5, 1e-4, 1e-3, 1e-2)
DIAMETERS_M = (0.001, 0.01, 0.05, 0.15, 0.30, 0.60)
CM1_CP = 1005.7  # J/kg/K: CM1 r22 constants.F; model closure, not design air's cp
CM1_P00 = 1e5  # Pa: CM1 r22 reference pressure


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def layout(path):
    """GrADS little-endian scalar layout, retaining record offsets and grid shape."""
    lines = Path(path).read_text().splitlines()
    dims = {a: int(next(l for l in lines if l.startswith(a + 'def')).split()[1]) for a in 'xyz'}
    start = next(i for i, line in enumerate(lines) if line.lower().startswith('vars'))
    count = int(lines[start].split()[1])
    fields, offset = {}, 0
    for line in lines[start + 1:start + 1 + count]:
        name, n = line.split()[:2]
        levels = max(1, int(n))
        fields[name] = (offset, levels)
        offset += levels * dims['x'] * dims['y']
    return dims, fields, offset


def read_selected(path, shape, fields, length, levels, digest):
    """Read/hash only named fields at declared heights, never the multi-GB run."""
    if Path(path).stat().st_size != length * 4:
        raise ValueError(f'GrADS size mismatch: {path}')
    raw = np.memmap(path, dtype='<f4', mode='r', shape=(length,))
    out = {}
    for name in FIELDS:
        offset, nz = fields[name]  # missing fields are errors, not substitutions
        field = np.asarray(raw[offset:offset + nz * shape['y'] * shape['x']].reshape(
            nz, shape['y'], shape['x'])[levels], dtype='<f4').copy()
        if not np.all(np.isfinite(field)) or np.any(field < -9e5):
            raise ValueError(f'Invalid {name}: {path}')
        digest.update(name.encode() + b'\0' + field.tobytes(order='C'))
        out[name] = field.astype(float)
    return out


def temperature(fields):
    if np.any(fields['prs'] <= 0):
        raise ValueError('Nonpositive pressure')
    return fields['th'] * (fields['prs'] / CM1_P00) ** (CM1_DRY_AIR_GAS_CONSTANT / CM1_CP)


def resolved_quantities(fields, heights_m, dx_m):
    """Grid-scale kinematics. A uniform translating wind has zero deformation."""
    u, v, w = (fields[k] for k in ('uinterp', 'vinterp', 'winterp'))
    shear = np.hypot(np.gradient(u, heights_m, axis=0), np.gradient(v, heights_m, axis=0))
    deformation2 = np.zeros_like(w)
    for f in (u, v, w):
        for axis in (1, 2):
            if f.shape[axis] > 1:
                derivative = (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2 * dx_m)
                deformation2 += derivative**2
    return {'vertical_wind_m_s': w, 'horizontal_wind_m_s': np.hypot(u, v),
            'vertical_horizontal_wind_shear_s-1': shear,
            'horizontal_velocity_gradient_norm_s-1': np.sqrt(deformation2)}


def winter_masks(t_kelvin, snow_kg_kg, relative_humidity_ice=None, cloud_water_kg_kg=None):
    winter = ((t_kelvin >= sf.T0 + WINTER_TEMPERATURE_C[0]) & (t_kelvin <= sf.T0))
    masks = {'all_winter_air': winter}
    masks.update({f'snow_ge_{threshold:g}_kg_kg': winter & (snow_kg_kg >= threshold)
                  for threshold in SNOW_THRESHOLDS})
    if relative_humidity_ice is not None and cloud_water_kg_kg is not None:
        snow = winter & (snow_kg_kg >= 1e-5)
        saturated = relative_humidity_ice >= 1.0
        liquid = cloud_water_kg_kg >= 1e-5
        masks.update(snow_and_ice_supersaturation=snow & saturated,
                     snow_and_liquid_cloud=snow & liquid,
                     snow_and_liquid_and_ice_supersaturation=snow & liquid & saturated)
    return masks


def longest_layer_depth(mask, layer_depth_m):
    """Maximum contiguous flagged depth at each column (no gap interpolation)."""
    mask = np.asarray(mask, bool)
    if mask.shape[0] != len(layer_depth_m):
        raise ValueError('Layer depths must match the height dimension')
    running = np.zeros(mask.shape[1:])
    longest = running.copy()
    for wet, dz in zip(mask, layer_depth_m):
        running = np.where(wet, running + dz, 0)
        longest = np.maximum(longest, running)
    return longest


def statistics(values):
    a = np.asarray(values, float)
    if a.size == 0:
        return {'count': 0, 'mean': None, 'p01': None, 'p50': None, 'p90': None, 'p99': None, 'min': None, 'max': None}
    q = np.percentile(a, (1, 50, 90, 99))
    return dict(count=int(a.size), mean=float(a.mean()), p01=float(q[0]), p50=float(q[1]), p90=float(q[2]),
                p99=float(q[3]), min=float(a.min()), max=float(a.max()))


def namelist_options(path):
    text = Path(path).read_text()
    names = ('cm1setup', 'sgsmodel', 'ipbl', 'horizturb', 'output_tke', 'output_km', 'output_kh')
    return {k: int(re.search(r'^\s*' + k + r'\s*=\s*(-?\d+)', text, re.M).group(1)) for k in names}


def analyse_case(name, max_snapshots=512, from_day=29.5):
    """Uniformly sample saved times, conditioning on temperature and local snow.

    Equal sample weights describe sampled cell-times (not particle trajectories,
    global weather probabilities, or snowfall at a community). All selected winter
    levels have equal 2-km spacing in the archived cases.
    """
    case = RUNS / name
    record = json.loads((case / 'case.json').read_text())
    shape, fields, length = layout(case / 'cm1out_s.ctl')
    if shape['x'] < 3 or shape['y'] < 3:
        raise ValueError('This diagnostic requires 3-D periodic boxes')
    opts = namelist_options(case / 'namelist.input')
    rawmeta = np.fromfile(case / 'cm1out_metadata.dat', dtype='<f4')
    _, mlayout, mlen = layout(case / 'cm1out_metadata.ctl')
    meta = rawmeta.reshape(-1, mlen)
    nindex, _ = mlayout['nwrite']
    tindex, _ = mlayout['mtime']
    times = {int(row[nindex]): float(row[tindex]) for row in meta}
    available = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    available = [n for n in available if times.get(n, -1) >= from_day * 86400]
    if not available:
        raise ValueError(f'No saved snapshots in {name} after {from_day} days')
    selection = np.unique(np.linspace(0, len(available) - 1, min(max_snapshots, len(available))).round().astype(int))
    use = [available[i] for i in selection]
    interfaces = np.loadtxt(case / 'input_grid_z')
    zall = 0.5 * (interfaces[1:] + interfaces[:-1])
    levels = np.where((zall >= HEIGHT_BAND_M[0]) & (zall <= HEIGHT_BAND_M[1]))[0]
    # Include neighbouring levels before differentiation; discard them from statistics.
    read_levels = np.arange(max(0, levels[0] - 1), min(len(zall), levels[-1] + 2))
    inner = np.isin(read_levels, levels)
    z = zall[read_levels]
    quantities = ('vertical_wind_m_s', 'horizontal_wind_m_s', 'vertical_horizontal_wind_shear_s-1',
                  'horizontal_velocity_gradient_norm_s-1', 'temperature_c', 'height_m', 'snow_g_kg',
                  'cloud_water_g_kg', 'relative_humidity_ice')
    groups = {key: {q: [] for q in quantities} for key in winter_masks(
        np.array([260.]), np.array([0.]), np.array([1.]), np.array([0.]))}
    depth_samples = {key: [] for key in ('snow', 'ice_supersaturated', 'near_ice_saturation',
                                       'liquid_cloud', 'liquid_and_ice_supersaturated')}
    source_rows = []
    for n in use:
        digest = hashlib.sha256()
        digest.update(np.asarray(read_levels, dtype='<i4').tobytes())
        d = read_selected(case / f'cm1out_t{n:06d}_s.dat', shape, fields, length, read_levels, digest)
        t = temperature(d)
        kinematics = resolved_quantities(d, z, record['grid']['dx_m'])
        vapour_pressure = d['prs'] * d['qv'] / (CM1_DRY_AIR_GAS_CONSTANT / CM1_VAPOUR_GAS_CONSTANT + d['qv'])
        ice_pressure = np.exp(9.550426 - 5723.265 / t + 3.53068 * np.log(t) - 0.00728332 * t)
        data = dict(kinematics, temperature_c=t - sf.T0, height_m=np.broadcast_to(z[:, None, None], t.shape),
                    snow_g_kg=d['qs'] * 1000, cloud_water_g_kg=d['qc'] * 1000,
                    relative_humidity_ice=vapour_pressure / ice_pressure)
        masks = winter_masks(t[inner], d['qs'][inner], (vapour_pressure / ice_pressure)[inner], d['qc'][inner])
        for key, mask in masks.items():
            for quantity, value in data.items():
                groups[key][quantity].append(value[inner][mask].astype(np.float32))
        snow = masks['snow_ge_1e-05_kg_kg']
        snowy_columns = snow.any(axis=0)
        winter = masks['all_winter_air']
        rhi, liquid = (vapour_pressure / ice_pressure)[inner], d['qc'][inner] >= 1e-5
        environments = dict(snow=snow, ice_supersaturated=winter & (rhi >= 1),
                            near_ice_saturation=winter & (rhi >= .95), liquid_cloud=winter & liquid,
                            liquid_and_ice_supersaturated=winter & liquid & (rhi >= 1))
        for key, environment in environments.items():
            depth_samples[key].append(longest_layer_depth(environment, np.diff(interfaces)[levels])[snowy_columns])
        source_rows.append(dict(output=n, model_day=times[n] / 86400, selected_fields_sha256=digest.hexdigest()))
    summaries = {}
    winter_count = sum(len(x) for x in groups['all_winter_air']['vertical_wind_m_s'])
    for group, data in groups.items():
        joined = {key: np.concatenate(value) for key, value in data.items()}
        w = joined['vertical_wind_m_s']
        n = len(w)
        summaries[group] = dict(sample_count=n, share_of_winter_cells=n / winter_count if winter_count else None,
            quantities={key: statistics(value) for key, value in joined.items()},
            fractions={f'abs_w_gt_{threshold:g}_m_s': float(np.mean(np.abs(w) > threshold)) if n else None
                       for threshold in (0.6, 2.0, 5.0)},
            upward_gt_0_6_m_s=float(np.mean(w > 0.6)) if n else None,
            ice_subsaturated_share=float(np.mean(joined['relative_humidity_ice'] < 1)) if n else None,
            ice_rh_below_0_95_share=float(np.mean(joined['relative_humidity_ice'] < .95)) if n else None,
            liquid_cloud_ge_0_01_g_kg_share=float(np.mean(joined['cloud_water_g_kg'] >= .01)) if n else None)
    return dict(case=name, source=f'climate/crm/runs/{name}', build=record['build']['label'],
        forcing_design=record['build'].get('air', {}).get('run'), grid=record['grid'], turbulence_options=opts,
        source_files={p: sha256(case / p) for p in ('case.json', 'namelist.input', 'cm1out_s.ctl',
                     'cm1out_metadata.ctl', 'cm1out_metadata.dat', 'input_grid_z')},
        input_fields=list(FIELDS), read_level_indices=read_levels.tolist(), sampled_height_m=zall[levels].tolist(),
        selected_byte_rule='For each output: little-endian int32 level indices, then each FIELDS name plus NUL '
                          'and its little-endian float32 array in C order (z,y,x) at those levels.',
        sampling=dict(method='Nearest equally spaced indices across available outputs after from_day; no storm selection '
                             'within a case. The storms_second case itself is a preselected active window.',
                      from_day=from_day, available_outputs=len(available), snapshots=len(use), records=source_rows),
        conditioning=dict(temperature_c=list(WINTER_TEMPERATURE_C), height_m=list(HEIGHT_BAND_M),
                          snow_thresholds_kg_kg=list(SNOW_THRESHOLDS), weighting='Equal sampled cell-times'),
        summaries=summaries,
        layer_depths_in_snowy_columns=dict(reading_rule='For each sampled column holding cold snow >=0.01 g/kg, '
            'maximum contiguous depth among flagged -20..0 C cells at 20-60 km. Environment masks need not '
            'themselves contain snow. Cell boundaries define depth; 2-km cells cannot locate thinner layers, '
            'and persistence along a flake trajectory is not measured. Zeros retained for snowy columns lacking '
            'that environment.', maximum_contiguous_depth_m={k: statistics(np.concatenate(v)) for k, v in depth_samples.items()}))


def increment_scale(epsilon_m2_s3, separation_m, viscosity_m2_s):
    """Order-of-magnitude velocity increment; prefactor 1, not an RMS prediction.

    Match a smooth-gradient estimate below eta to inertial scaling above it.
    The latter also requires separation well below the turbulence integral scale.
    """
    if epsilon_m2_s3 < 0 or separation_m <= 0 or viscosity_m2_s <= 0:
        raise ValueError('Need nonnegative dissipation, positive separation and viscosity')
    if epsilon_m2_s3 == 0:
        return 0.0
    eta = (viscosity_m2_s**3 / epsilon_m2_s3)**0.25
    return (separation_m * math.sqrt(epsilon_m2_s3 / viscosity_m2_s) if separation_m < eta
            else (epsilon_m2_s3 * separation_m)**(1 / 3))


def terminal_speed(diameter, geometry, gravity, state):
    """Existing Heymsfield-Westbrook law evaluated with the named geometry."""
    mass = float(geometry.mass(diameter))
    best = 8 * mass * gravity * state['rho'] / (math.pi * math.sqrt(geometry.area_ratio) * state['eta']**2)
    re = sf.DELTA0**2 / 4 * (math.sqrt(1 + 4 * math.sqrt(best) / (sf.DELTA0**2 * math.sqrt(sf.C0))) - 1)**2
    return state['eta'] * re / (state['rho'] * diameter)


def flake_scales():
    """Snow-column air and hypothetical local dissipation, not fitted to CM1."""
    moon, binding = sf.moon_column()
    height = moon.height_of(sf.T0 - 10)
    state = moon.state(height)
    nu = state['eta'] / state['rho']
    geometry = GEOMETRIES['bf95_dmax_dendrites_area']
    rows = []
    for epsilon in EPSILONS:
        flakes = []
        for diameter in DIAMETERS_M:
            du = increment_scale(epsilon, diameter, nu)
            mass = float(geometry.mass(diameter))
            speed = terminal_speed(diameter, geometry, moon.gravity, state)
            flakes.append(dict(diameter_m=diameter, mass_kg=mass, fall_speed_m_s=speed,
                increment_scale_m_s=du, eddy_time_s=diameter / du,
                increment_to_fall_speed_ratio=du / speed,
                dynamic_pressure_scale_pa=0.5 * state['rho'] * du**2,
                settling_dynamic_pressure_pa=0.5 * state['rho'] * speed**2,
                aligned_slip_pressure_increment_pa=state['rho'] * speed * du + 0.5 * state['rho'] * du**2,
                weight_per_projected_area_pa=mass * moon.gravity / float(geometry.area(diameter))))
        rows.append(dict(epsilon_m2_s3=epsilon, kolmogorov_length_m=(nu**3 / epsilon)**.25,
                         kolmogorov_time_s=math.sqrt(nu / epsilon), flakes=flakes))
    return dict(evidence='Uncalibrated scenarios. Epsilon is imposed, not recovered from CM1; no probabilities '
                        'or lunar turbulence classes are assigned. No bond strength or breakup criterion is fitted.',
        reference_air=dict(temperature_c=-10, height_m=height, density_kg_m3=state['rho'],
                           viscosity_m2_s=nu, gravity_m_s2=moon.gravity),
        column_input=binding,
        reference_geometry=geometry.to_dict(),
        formulas=dict(inertial_increment='(epsilon D)^(1/3), order-one prefactor fixed to 1',
                      viscous_increment='D sqrt(epsilon / nu)',
                      dynamic_pressure='rho (delta_u)^2 / 2; loading scale only, not tensile stress on a branch',
                      aligned_slip_pressure_increment='rho V delta_u + rho (delta_u)^2 / 2, from changing '
                          'settling slip V to V+delta_u. An adverse aligned fluctuation envelope, not actual '
                          'spatial stress or an intermittent maximum. Assumes unchanged drag coefficient.',
                      eddy_time='D / delta_u'),
        reading_rule='D must be within an approximately isotropic locally equilibrated cascade. At D below '
                     'the Kolmogorov length the smooth-gradient expression is used. Mean translation does not '
                     'enter these increments. No survival claim follows by comparing pressure with weight: '
                     'contacts, drag fluctuations, flutter and collisions need separate mechanical laws.', rows=rows)


def analyse(max_snapshots=512):
    sources = json.loads((HERE / 'turbulence_sources.json').read_text())
    return dict(schema=SCHEMA, producer=dict(domain='atmosphere', model='atmosphere/snow/turbulence.py',
                files={p: sha256(ROOT / p) for p in CODE}, constants=constants_used(CODE)),
        source_register_sha256=sha256(HERE / 'turbulence_sources.json'),
        evidence='Diagnostic of archived lunar CM1 simulations at kilometre resolution, conditioned on cold '
                 'snow-bearing cells, plus independent sub-metre turbulence scenarios. Neither empirical lunar '
                 'validation nor a prediction of giant-flake survival or resort weather.',
        reading_rule='Keep each run separate. The two full-cycle boxes use different microphysics and the third '
                     'is a selected storm window, not an occurrence sample. Summaries weight sampled cell-times '
                     'equally. The snow mask describes local model snow mass, not flakes crossing a surface. '
                     'Grid-scale velocities and gradients are not dissipation estimates. Epsilon scenarios '
                     'cannot be assigned to sampled cloud cells without an additional closure.',
        sources=sources, cases=[analyse_case(name, max_snapshots) for name in CASES], scenarios=flake_scales())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-snapshots', type=int, default=512)
    args = parser.parse_args()
    if args.max_snapshots < 1:
        parser.error('--max-snapshots must be positive')
    result = analyse(args.max_snapshots)
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    for case in result['cases']:
        print(case['case'], {k: dict(count=v['sample_count'], abs_w_over_06=v['fractions']['abs_w_gt_0.6_m_s'],
                                    w=v['quantities']['vertical_wind_m_s']['p99'])
                            for k, v in case['summaries'].items()})


if __name__ == '__main__':
    main()
