"""Check the frozen two-temperature calibration against newer growth evidence.

This is a bounded response audit, not another fit. The native executable and
all twelve controls are unchanged. Relative excess is measured between ice
and liquid saturation; it is distinct from ordinary ice supersaturation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

from atmosphere.snow import mcsnow_column_density as column
from atmosphere.snow import mcsnow_joint_acceptance as acceptance

HERE = Path(column.__file__).resolve().parent
CACHE = column.tuning.trace.CACHE / 'science_check'
PRODUCT = HERE / 'results/mcsnow_science_check.json'
BUILD = column.CACHE / 'build'
EXCESS = (.05, .1, .2, .4, 1.)
DIAMETERS_UM = (20., 60., 250.)
PROBE_DURATION = 1 / 32


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def constants():
    """Read the pinned native constants rather than introduce a second air law."""
    source = (BUILD / 'source/src/mo_atmo_types.f90').read_text()
    names = ('T_3', 'e_3', 'A_i', 'A_w', 'B_i', 'B_w', 'rhoi')
    result = {k: float(re.search(r'\b'+k+r'\s*=\s*([\d.eE+-]+)_wp', source).group(1)) for k in names}
    functions = (BUILD / 'source/src/mo_atmo.f90').read_text()
    for phase in ('i', 'w'):
        expected = f'psat{phase} = e_3 * EXP( A_{phase} * (T - T_3) / (T - B_{phase}) )'
        if expected not in functions:
            raise ValueError('Native saturation equation changed')
    return result


def forcing(temperature_C, relative_excess, values):
    if not math.isfinite(relative_excess) or not 0 <= relative_excess <= 1:
        raise ValueError('Relative excess must lie between ice and water saturation')
    temperature = values['T_3'] + temperature_C
    ice, water = [values['e_3'] * math.exp(values['A_'+p]*temperature_C / (temperature-values['B_'+p])) for p in ('i', 'w')]
    vapor = ice + relative_excess*(water-ice)
    return dict(temperature_K=temperature, water_saturation_ratio=vapor/water,
                ice_supersaturation=vapor/ice-1, relative_excess=relative_excess,
                ice_saturation_pressure_Pa=ice, water_saturation_pressure_Pa=water)


def native_cases(cases, controls, path):
    """Run one unchanged native worker; preserve the exact input and output."""
    path.mkdir(parents=True, exist_ok=True)
    lines = [str(len(cases))]
    for case in cases:
        if not 1 <= len(case['times_s']) <= 8:
            raise ValueError('The native driver accepts one to eight output times')
        lines.extend([' '.join(format(x, '.17g') for x in case['settings'])+' '+str(len(case['times_s'])),
                      ' '.join(format(x, '.17g') for x in case['times_s'])])
    case_path = path / 'cases.txt'
    case_path.write_text('\n'.join(lines)+'\n')
    control_text = ' '.join(format(x, '.17g') for x in controls)+'\n'
    (path / 'controls.txt').write_text(control_text)
    result = subprocess.run([str(BUILD / 'tuner'), str(case_path)], input=control_text,
                            cwd=BUILD, env=column.native.ENV, text=True, capture_output=True, timeout=120)
    (path / 'stdout.txt').write_text(result.stdout)
    (path / 'stderr.txt').write_text(result.stderr)
    result.check_returncode()
    before, marker, payload = result.stdout.partition('READY\n')
    output, done, trailing = payload.partition('DONE\n')
    if not marker or not done or trailing.strip():
        raise ValueError('Incomplete native response')
    rows = [list(map(float, line.split())) for line in output.splitlines() if line.strip()]
    if len(rows) != sum(len(c['times_s']) for c in cases):
        raise ValueError('Native row coverage changed')
    for row in rows:
        if len(row) != 8 or not all(map(math.isfinite, row)) or min(row[:3]) <= 0 or row[3] < 0:
            raise ValueError('Invalid native state')
        if abs(row[4]) > 1e-10 or abs(row[5]) > 1e-10 or not 0 < row[7] <= 1+1e-10:
            raise ValueError('Native mass, volume or density invariant failed')
    return rows


def incremental_rates(initial, final, seconds):
    dm = final[0]-initial[0]
    volume = lambda r: math.pi/6*r[1]**2*r[2]
    dv = volume(final)-volume(initial)
    da, dc = [math.log(final[i]/initial[i]) for i in (1, 2)]
    if min(dm, dv, da, dc, seconds) <= 0:
        raise ValueError('Expected positive growth in the selected finite-time probe')
    return dict(mass_rate_kg_s=dm/seconds, deposition_density_kg_m3=dm/dv,
                polar_over_equatorial_log_growth=dc/da)


def run():
    column.verify_build()
    sources = json.loads((HERE / 'mcsnow_science_check_sources.json').read_text())
    for record in sources['sources']:
        cached = record.get('cached_source')
        if cached and sha(HERE.parents[1] / cached['path']) != cached['sha256']:
            raise ValueError('Cached evidence identity changed: '+record['id'])
    baseline = acceptance.load_baseline()
    values = constants()
    records, parity, refinements = [], [], []
    for temperature in (-5.3, -14.4):
        fit = next(f for f in baseline['searches'] if f['name'] == baseline['preferred_search']) if temperature == -5.3 else baseline['previous_cold']
        controls = list(fit['physical_controls'])
        if temperature == -14.4:
            controls += [0., 1.]  # Neutral new warm-only density controls.
        tag = column.diagnosis.key(temperature)
        for step in (1/512, 1/1024):
            cases = []
            for diameter in DIAMETERS_UM:
                for sigma in EXCESS:
                    air = forcing(temperature, sigma, values)
                    d = diameter*1e-6
                    cases.append(dict(temperature_C=temperature, diameter_um=diameter, **air,
                        settings=[air['temperature_K'], 101000., air['water_saturation_ratio'], -1., d, d, step],
                        times_s=[0., PROBE_DURATION]))
            rows = native_cases(cases, controls, CACHE / f'probes_{tag}_{step:g}')
            for i, case in enumerate(cases):
                first, last = rows[2*i:2*i+2]
                records.append(dict(**case, step_s=step, initial=first, final=last,
                                    **incremental_rates(first, last, PROBE_DURATION)))
        # Reproduce supplied baseline values at their original timestep and starts.
        cases, targets = [], []
        for entry in column.tuning.cases():
            if entry['temperature_C'] != temperature or entry['initialization'] not in ('sphere', 'axes'):
                continue
            selected = sorted((r for r in fit['assessment']['rows'] if r['initialization'] == entry['initialization']
                               and r['age_min'] in (3., 5., 10., 15., 20., 30.)), key=lambda r: r['age_min'])
            c = entry['configuration']
            cases.append(dict(settings=[c['temperature_K'], c['pressure_Pa'], 1., c['initial_mass_kg'],
                               c['initial_equatorial_m'], c['initial_polar_m'], 1/128],
                              times_s=[(r['age_min']-entry['initial_age_min'])*60 for r in selected]))
            targets.extend(selected)
        actual = native_cases(cases, controls, CACHE / f'baseline_{tag}')
        difference = max(abs(v/t-1) for row, target in zip(actual, targets, strict=True) for v, t in zip(row[:4], target['values'], strict=True))
        parity.append(dict(temperature_C=temperature, states=len(targets), maximum_relative_difference=difference,
                           passed=difference <= 1e-12))
    for coarse in [r for r in records if r['step_s'] == 1/512]:
        fine = next(r for r in records if r['step_s'] == 1/1024 and
                    all(r[k] == coarse[k] for k in ('temperature_C', 'diameter_um', 'relative_excess')))
        difference = max(abs(fine[k]/coarse[k]-1) for k in ('mass_rate_kg_s', 'deposition_density_kg_m3', 'polar_over_equatorial_log_growth'))
        refinements.append(dict(temperature_C=coarse['temperature_C'], diameter_um=coarse['diameter_um'],
                                relative_excess=coarse['relative_excess'], maximum_relative_difference=difference,
                                passed=difference <= .01))
    checks = dict(baseline_parity=parity, timestep_refinement=refinements,
                  mass_balance_max=max(abs(r['final'][4]) for r in records),
                  volume_residual_max=max(abs(r['final'][5]) for r in records))
    if not all(c['passed'] for c in parity+refinements):
        raise ValueError('Numerical verification failed')
    raw_files = {str(p.relative_to(CACHE)): sha(p) for p in sorted(CACHE.glob('*/*')) if p.is_file()}
    source_files = ['mcsnow_science_check.py', 'mcsnow_science_check_sources.json', 'mcsnow_column_density.py',
                    'mcsnow_tuning.py', 'mcsnow_tune_driver.f90']
    result = dict(schema='terluna.atmosphere.mcsnow-science-check/1',
        evidence='Finite-time response probes of the unchanged calibrated component; not matched reconstructions of newer experiments.',
        baseline_commit=acceptance.BASELINE_COMMIT, baseline_sha256=acceptance.BASELINE_SHA256,
        producer_inputs={name: sha(HERE / name) for name in source_files},
        build_sha256=sha(BUILD / 'build.json'), executable_sha256=sha(BUILD / 'tuner'),
        source_sha256={name: sha(BUILD / 'source/src' / name) for name in
                      ('mo_depo.f90', 'mo_atmo.f90', 'mo_atmo_types.f90', 'mo_velocity.f90')},
        native_constants=values, probe_states=len(records), probe_records=records,
        checks=checks, raw_directory=str(CACHE), raw_files=raw_files,
        review_outcome='Humidity-dependent compact-to-complex morphology is absent; the proposed aggregation continuation condition is not met.',
        limits=['The inferred added-volume density is a model quantity, not a separately measured density.',
                'The published morphology threshold is an approximate range, not an exact 10-percent target.',
                'The new 2025 cold-chamber measurements lie outside the two fitted temperatures.',
                'No parameters, acceptance thresholds, earlier products or model source were changed.'])
    PRODUCT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(probe_states=len(records), parity=parity,
                         maximum_refinement=max(r['maximum_relative_difference'] for r in refinements),
                         product=str(PRODUCT))))


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    run()
