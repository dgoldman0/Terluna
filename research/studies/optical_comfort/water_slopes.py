"""Connect the local coastal wave history to resolved water-surface slopes.

OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.water_slopes
"""
from __future__ import annotations

import json
from pathlib import Path
import re

import numpy as np

from illumination.water_surface.model import fraction_above, iter_swan, slope_moments
from shared.constants import MOON_SURFACE_GRAVITY, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used
from .run import HERE, ROOT, digest

SCHEMA = 'terluna.research.optical-water-slopes/1'
SOURCE = ROOT / 'research/runs/waves/shore_history/coast_s4_dt900_0_1419'
MANIFEST = HERE / 'local_inputs.json'


def refine(f, e):
    x = np.empty(2 * len(f) - 1); x[::2] = f; x[1::2] = (f[1:] + f[:-1]) / 2
    y = np.empty((len(x), e.shape[1])); y[::2] = e; y[1::2] = (e[1:] + e[:-1]) / 2
    return x, y


def quantiles(values):
    return dict(zip(('minimum', 'p10', 'p50', 'p90', 'maximum'),
                    map(float, np.quantile(values, [0, .1, .5, .9, 1]))))


def build():
    admitted = json.loads(MANIFEST.read_text())['coastal_wave_history']
    if digest(SOURCE / 'product.json') != admitted['product_sha256']:
        raise ValueError('Coastal product differs from the admitted lunar-seas result')
    source = json.loads((SOURCE / 'product.json').read_text())
    if source['schema'] != 'terluna.climate.coastal-history-case/1':
        raise ValueError('Unsupported coastal history schema')
    input_hashes = source['run']['identity']['inputs']
    hashes = {}
    for name in ('sites.spc', 'sites.tbl', 'INPUT'):
        actual = digest(SOURCE / name)
        expected = input_hashes[name] if name == 'INPUT' else source['run']['output_sha256'][name]
        if actual != expected:
            raise ValueError(f'Coastal source hash mismatch: {name}')
        hashes[name] = actual
    commands = (SOURCE / 'INPUT').read_text()
    gravity = float(re.search(r'SET GRAV\s+([0-9.eE+\-]+)', commands)[1])
    if not np.isclose(gravity, MOON_SURFACE_GRAVITY, rtol=1e-10):
        raise ValueError('Lunar gravity convention differs')
    if re.search(r'(?im)^\s*(?:INPGRID\s+CUR|READINP\s+CUR)', commands):
        raise ValueError('Absolute-to-intrinsic frequency conversion needs the current field')
    names = list(source['sites'])
    locations = np.array(list(source['sites'].values()))
    bulk = np.loadtxt(SOURCE / 'sites.tbl').reshape(-1, len(names), len(source['columns']))
    ti, di, hi = (source['columns'].index(k) for k in ('time_s', 'depth_m', 'hs_m'))
    lookup = {float(v[0, ti]): v for v in bulk}
    bounds = np.array([1, 2]) * SYNODIC_MONTH_DAYS * 24
    rows, errors, coarse_errors, height_differences = [], [], [], []
    origin = None
    frequency = None
    for record in iter_swan(SOURCE / 'sites.spc'):
        if origin is None:
            origin = record['time']
        hour = (record['time'] - origin).total_seconds() / 3600 + source['start_hour']
        if not bounds[0] <= hour <= bounds[1]:
            continue
        if not np.allclose(record['locations'], locations, rtol=0, atol=1e-6):
            raise ValueError('Coastal station ordering differs')
        if frequency is not None and not np.array_equal(frequency, record['frequency']):
            raise ValueError('Frequency grid changes within the history')
        frequency = record['frequency']
        if not record['available'].all():
            raise ValueError('A requested station lacks a directional spectrum')
        station = lookup[hour * 3600]
        for i, name in enumerate(names):
            depth = float(station[i, di])
            f, e = frequency, record['variance'][i]
            levels = []
            for level in range(3):
                levels.append(slope_moments(f, record['direction'], e, depth, gravity))
                if level < 2:
                    f, e = refine(f, e)
            m = levels[-1]
            cov = m['covariance']
            total = float(np.trace(cov))
            hs = float(4 * np.sqrt(m['variance_m2']))
            if hs >= .1 and total > 0:
                coarse_errors.append(abs(np.trace(levels[0]['covariance']) / total - 1))
                errors.append(abs(np.trace(levels[1]['covariance']) / total - 1))
                height_differences.append(abs(hs / station[i, hi] - 1))
            rows.append(dict(site=name, hour=hour, depth_m=depth, resolved_hs_m=hs,
                             rms_gradient=float(np.sqrt(total)), covariance_east_north=cov.tolist(),
                             principal_variances=np.linalg.eigvalsh(cov).tolist(),
                             slope_fraction_above_025hz=fraction_above(f, m['slope_density'], .25),
                             elevation_variance_fraction_above_025hz=fraction_above(f, m['elevation_density'], .25),
                             slope_fraction_in_last_two_native_intervals=fraction_above(f, m['slope_density'], frequency[-3]),
                             dispersion_relative_residual=m['dispersion_relative_residual']))
    if not rows:
        raise ValueError('Second-cycle spectra absent')
    summaries = {}
    for name in names:
        selected = [r for r in rows if r['site'] == name]
        summaries[name] = dict(location_deg=source['sites'][name], depth_m=selected[0]['depth_m'],
                              snapshots=len(selected),
                              **{key: quantiles([r[key] for r in selected]) for key in
                                 ('resolved_hs_m', 'rms_gradient', 'slope_fraction_above_025hz',
                                  'elevation_variance_fraction_above_025hz', 'slope_fraction_in_last_two_native_intervals')})
    history = ROOT / 'research/runs/optical_comfort/water_slope_history.json'
    history.write_text(json.dumps(dict(schema=SCHEMA, reading_rule='One row per hourly spectrum and station; '
                                       'frequency integration uses four subdivisions per native interval.', rows=rows),
                                      separators=(',', ':'), allow_nan=False) + '\n')
    files = [Path(__file__), MANIFEST, HERE / 'local_sources.json', ROOT / 'illumination/water_surface/model.py']
    return dict(schema=SCHEMA, producer=dict(lane='research', runner=str(Path(__file__).relative_to(ROOT)),
                files={str(p.relative_to(ROOT)): digest(p) for p in files}, constants=constants_used(files),
                input=dict(file=str((SOURCE / 'product.json').relative_to(ROOT)), schema=source['schema'],
                           sha256=admitted['product_sha256'], verified_files=hashes,
                           branch='research/lunar-seas', source_commit=admitted['source_commit']),
                history=dict(file=str(history.relative_to(ROOT)), sha256=digest(history))),
                evidence='Resolved gravity-wave slope covariance from the second-cycle SWAN coastal history. '
                         'Wave variance is weighted by wavenumber squared at the simulated lunar gravity and local depth. '
                         'The spectral upper cutoff leaves shorter-wave roughness open; this result supplies part '
                         'of a water-reflection model.',
                reading_rule='Cartesian east/north slope covariance. The source spectrum is m2/Hz/degree; periodic '
                             'direction integration uses the printed bin width. Linear frequency interpolation of '
                             'variance is integrated with four subdivisions per native interval. Dispersion is '
                             'omega²=g*k*tanh(k*depth). The run prescribes zero current, so AFREQ is intrinsic here. '
                             'Quantiles weight saved hourly samples equally. Spin-up precedes the first requested '
                             'time. Spectral energy beyond the printed frequency range is omitted. '
                             'RMS gradient is sqrt(trace(covariance)); a probability distribution of facets and '
                             'a glitter image require additional assumptions. Earth wind-to-slope fits are unused.',
                units=dict(frequency='Hz', depth='m', wave_height='m', slope_covariance='dimensionless',
                           rms_gradient='dimensionless', time='hours since wave-run start'),
                source_evidence=source['evidence'], source_reading_rule=source['reading_rule'],
                requested_span_hours=bounds.tolist(), sampled_span_hours=[rows[0]['hour'], rows[-1]['hour']],
                gravity_m_s2=gravity, frequency_range_hz=[float(frequency[0]), float(frequency[-1])],
                native_frequency_count=len(frequency), summaries=summaries,
                checks=dict(maximum_slope_integral_relative_change_1_to_4=max(coarse_errors),
                            maximum_slope_integral_relative_change_2_to_4=max(errors),
                            maximum_resolved_hs_relative_difference_from_bulk=max(height_differences),
                            maximum_dispersion_relative_residual=max(r['dispersion_relative_residual'] for r in rows),
                            comparison_scope='Quadrature and Hs diagnostics use samples with resolved Hs >= 0.1 m. '
                                             'Bulk SWAN Hs also includes its diagnostic high-frequency tail.'))


def main():
    product = build()
    path = HERE / 'results/water_slopes.json'
    path.write_text(json.dumps(product, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(checks=product['checks'], summaries=product['summaries']), indent=2))


if __name__ == '__main__':
    main()
