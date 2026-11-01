"""Audit published geometry and extract checks from measured hollow-crystal series.

The 2026 archive supplies shape functions, surrogate 1991 observations and two
substrate time series. This runner does not implement the paper's growth theory.
Original external bytes stay in the ignored cache; compact derived diagnostics
retain conditions, uncertainty and source identities.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import numpy as np
from atmosphere.snow import hollow_shape as shape, takahashi_growth as lab

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = HERE / 'hollow_growth_sources.json'
RESULT = HERE / 'results/hollow_growth_assessment.json'


def read_sources():
    record = json.loads(SOURCES.read_text())
    for item in record['retrieved']:
        path = ROOT / item['cache_path']
        if shape.sha(path) != item['sha256'] or path.stat().st_size != item['bytes']:
            raise ValueError('Hollow-growth retrieval identity differs: ' + str(path))
    return record


def line_fit(time_s, values, absolute_error):
    """Descriptive OLS slope with a deterministic measurement-error envelope.

    Reported errors are bounds, not independent Gaussian standard deviations.
    The slope is a linear combination of measurements, so sum |weight|*error
    encloses its error for all perturbations within the supplied intervals.
    This does not bound curvature, forcing errors or unknown systematic bias.
    """
    time_s, values, absolute_error = map(lambda x: np.asarray(x, dtype=float),
                                        (time_s, values, absolute_error))
    if any(x.ndim != 1 for x in (time_s, values, absolute_error)) or not (
            time_s.shape == values.shape == absolute_error.shape) or len(time_s) < 3:
        raise ValueError('At least three matching one-dimensional measurements required')
    if not all(np.all(np.isfinite(x)) for x in (time_s, values, absolute_error)):
        raise ValueError('Finite measurements required')
    if np.any(np.diff(time_s) <= 0) or np.any(absolute_error < 0):
        raise ValueError('Increasing times and nonnegative errors required')
    centered = time_s - time_s.mean()
    weights = centered / np.dot(centered, centered)
    slope = float(np.dot(weights, values))
    error = float(np.dot(np.abs(weights), absolute_error))
    residual = values - (values.mean() + slope * centered)
    return dict(slope=slope, measurement_error_bound=error, interval=[slope-error, slope+error],
                residual_rms=float(np.sqrt(np.mean(residual**2))), measurements=len(time_s),
                time_range_s=[float(time_s[0]), float(time_s[-1])])


def parse_chamber(raw):
    lines = raw.decode().splitlines()
    if len(lines) < 8 or not lines[0].startswith('Experiment ') or 'time' not in lines[6]:
        raise ValueError('Unexpected chamber metadata/header')
    data = np.loadtxt(io.StringIO('\n'.join(lines[7:])))
    if data.ndim != 2 or data.shape[1] != 10 or not np.all(np.isfinite(data)):
        raise ValueError('Expected ten finite chamber columns')
    if np.any(np.diff(data[:, 0]) <= 0) or np.any(data[:, [1, 2, 7]] <= 0):
        raise ValueError('Nonincreasing time or nonpositive dimensions')
    if np.any(data[:, [3, 4, 5, 6, 8, 9]] < 0):
        raise ValueError('Negative chamber uncertainty')
    if np.any(data[:, [5, 6, 9]] < data[:, [3, 4, 8]]):
        raise ValueError('Total error smaller than camera difference')
    if np.any(data[:, 7] > data[:, 1]):
        raise ValueError('Reported rim wider than a semidimension; inspect convention')
    return lines[:7], data


def window(data, begin_min, end_min):
    selected = data[(data[:, 0] >= begin_min * 60) & (data[:, 0] < end_min * 60)]
    fits = {key: line_fit(selected[:, 0], selected[:, col], selected[:, err])
            for key, col, err in (('a_um_s', 1, 5), ('c_um_s', 2, 6), ('rim_um_s', 7, 9))}
    fit_a, fit_c = fits['a_um_s'], fits['c_um_s']
    ratio = fit_c['slope'] / fit_a['slope'] if fit_a['slope'] > 0 else None
    ratio_interval = None
    if fit_a['interval'][0] > 0 and fit_c['interval'][0] > 0:
        ratio_interval = [fit_c['interval'][0] / fit_a['interval'][1],
                          fit_c['interval'][1] / fit_a['interval'][0]]
    rim_interval = [float(np.median(selected[:, 7] - selected[:, 9])),
                    float(np.median(selected[:, 7] + selected[:, 9]))]
    return dict(requested_interval_min=[begin_min, end_min], measurements=len(selected),
        fits=fits, fitted_directional_growth_ratio=ratio, directional_ratio_interval=ratio_interval,
        rim_width_median_um=float(np.median(selected[:, 7])),
        rim_median_measurement_interval_um=rim_interval,
        rim_fraction_median=float(np.median(selected[:, 7] / selected[:, 1])),
        a_median_um=float(np.median(selected[:, 1])), c_median_um=float(np.median(selected[:, 2])),
        uncertainty='Deterministic envelopes from the source total measurement errors; not statistical confidence intervals. Null ratio interval means a positive denominator is not resolved.')


def chamber_controls(members):
    cases = []
    for name, pressure, step in (('dimensions-20231128.dat', 97000., False),
                                  ('dimensions-20240814.dat', 97200., True)):
        header, data = parse_chamber(members[name])
        if header[1] != 'Temperature: -50C' or header[3] != f'Pressure: {pressure/100:g} hPa':
            raise ValueError('Unexpected chamber forcing')
        expected = 'Supersaturation: 48%, switch to 20% at 230 min' if step else 'Supersaturation: 48%'
        if header[2] != expected:
            raise ValueError('Unexpected chamber humidity history')
        windows = [window(data, 150., 230.), window(data, 270., 336.)] if step else [window(data, 0., data[-1,0]/60+1.)]
        cases.append(dict(member=name, member_sha256=hashlib.sha256(members[name]).hexdigest(),
            measurements=len(data), temperature_K=223.15, pressure_Pa=pressure,
            ice_supersaturation_before=.48, ice_supersaturation_after=.20 if step else None,
            switch_elapsed_s=230*60 if step else None,
            elapsed_time_range_s=[float(data[0,0]), float(data[-1,0])], windows=windows,
            convention='a and c treated as semidimensions following the companion 2025 methods, equations 2-3; the archive header does not separately restate the factor of two. Rim is the reported physical width.',
            source_header=header,
            limits=['Substrate crystal, not free fall; no measured mass in these files.',
                    'Elapsed time starts with the recorded experiment, not necessarily nucleation.',
                    'Nominal forcing from the archive; no invented temperature/humidity error bars.',
                    'OLS windows summarize trends and do not reconstruct instantaneous facet rates.',
                    'The step-case windows exclude the first forty minutes after the humidity switch; no dynamic law was fitted.']))
    before, after = cases[1]['windows']
    pre_bounds = before['rim_median_measurement_interval_um']
    post_bounds = after['rim_median_measurement_interval_um']
    return dict(cases=cases, step_response=dict(
        median_rim_width_ratio=after['rim_width_median_um']/before['rim_width_median_um'],
        median_rim_width_ratio_measurement_interval=[post_bounds[0]/pre_bounds[1], post_bounds[1]/pre_bounds[0]],
        median_rim_fraction_ratio=after['rim_fraction_median']/before['rim_fraction_median'],
        c_slope_after_to_before=after['fits']['c_um_s']['slope']/before['fits']['c_um_s']['slope'],
        a_slope_after_to_before=after['fits']['a_um_s']['slope']/before['fits']['a_um_s']['slope']),
        evidence='Calculated summaries of two published substrate measurements; no growth-model validation or mass estimate.')


def free_fall_controls(members):
    cases = []
    for axis in ('a', 'c'):
        name = axis + '-dimension-TF1991-raw.csv'
        text = members[name].decode()
        if 'max dimension (mm)' not in text.splitlines()[0]:
            raise ValueError('Unrecognized free-fall dimension convention')
        data = np.genfromtxt(io.StringIO(text), delimiter=',', skip_header=3)
        if data.shape[1] != 5:
            raise ValueError('Unexpected surrogate observation columns')
        for col, temperature in enumerate((-5.3, -3.7, -8.6, -10.6), 1):
            compared = []; excluded = []
            for age, dimension in data[:, [0, col]]:
                if not np.isfinite(dimension):
                    continue
                if not 3 <= age <= lab.FITS[temperature][0]:
                    excluded.append(dict(age_min=float(age), reason='Outside existing selected fit interval'))
                    continue
                observed = lab.observed(temperature, float(age))
                fitted = 1000 * observed['equatorial_dimension_m' if axis == 'a' else 'polar_dimension_m']
                compared.append(float(dimension / fitted))
            cases.append(dict(axis=axis, temperature_C=temperature, member=name,
                compared_points=len(compared), excluded_points=excluded,
                digitization_to_existing_fit_ratio=dict(minimum=min(compared), median=float(np.median(compared)), maximum=max(compared))))
    # The separate warm-column files are repeated early points, not new crystals.
    duplicated = {}
    for axis in ('a', 'c'):
        all_rows = np.genfromtxt(io.BytesIO(members[axis+'-dimension-TF1991-raw.csv']), delimiter=',', skip_header=3)
        subset = np.loadtxt(io.BytesIO(members[axis+'-dim-5.3C-col.csv']), delimiter=',', skiprows=1)
        for row in subset:
            if not np.any(np.all(all_rows[:, :2] == row, axis=1)):
                raise ValueError('Separate column file no longer matches main digitization')
        duplicated[axis] = len(subset)
    return dict(cases=cases, duplicate_early_column_points=duplicated,
        compared_measurements=sum(row['compared_points'] for row in cases),
        evidence='Surrogate digitizations of the same 1991 experiment as the existing fits; not independent validation.',
        limits=['No pointwise errors or exact a/c pairing supplied.',
                'No extrapolation to digitized ages outside the existing selected fit interval.',
                'A size-fit comparison does not validate derivatives of that fit or reconstruct nucleation.'])


def run(output=RESULT):
    identity = {str(p.relative_to(ROOT)): shape.sha(p) for p in
                (Path(__file__), Path(lab.__file__), SOURCES)}
    identity.update(shape.producer())
    sources = read_sources(); members = shape.archive_members()
    for name, raw in members.items():
        if sources['archive_members'][name] != dict(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()):
            raise ValueError('Archive member identity differs from the source record: ' + name)
    result = dict(schema='terluna.atmosphere.hollow-growth-assessment/1', producer=identity,
        sources=sources, geometry_component=shape.audit(shape.ShapeLibrary()),
        chamber=chamber_controls(members), free_fall=free_fall_controls(members),
        evidence='Published geometry reproduced and measurement constraints assessed; no adoption of the 2026 growth model or a lunar prediction.',
        outstanding=['The full 2026 growth equations and implementation have not been obtained; publisher retrieval returned 403 and the accessible author preview supplied only the first page.',
                     'The archive supplies shape functions, not a complete time-stepping hollow-growth solver.',
                     'Rim evolution, facet step sources and thermal coupling need explicit tested equations.',
                     'No new free-fall or lunar aggregate validation follows from these substrate measurements.'])
    if any(shape.sha(ROOT/name) != digest for name,digest in identity.items()):
        raise RuntimeError('Producer changed during assessment')
    Path(output).write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULT)
    args=parser.parse_args()
    result=run(args.output)
    print(json.dumps(dict(output=str(args.output), geometry=result['geometry_component']['comparisons'],
                         step_response=result['chamber']['step_response']), indent=2))
