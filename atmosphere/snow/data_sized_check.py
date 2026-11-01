"""Check the frozen crystal model against each measured laboratory crystal.

The data-sized rule adopted on 24 October 2026 compares the model with every
crystal Takahashi et al. (1991) measured, at the crystal's own age and nominal
temperature.  On mass and on each measured axis a crystal passes when the model
lies within two pooled standard deviations of the scatter between repeat
crystals, in natural logarithms.  The archive holds no fall speeds crystal by
crystal, so fall speed is compared with the paper's Table 4 relation within its
stated 3 cm/s measurement accuracy.

Two frozen models are checked.  The improved McSnow component, which drives
every lunar column and storm calculation, is read from its retained native
traces at all ten laboratory temperatures.  The two-temperature calibration is
run unchanged through its frozen executable at -5.3 and -14.4 C.  Nothing here
fits or changes a model parameter.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from atmosphere.snow import takahashi_growth as lab
from atmosphere.snow import crystal_transport_replay as transport

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRODUCT = HERE / 'results/data_sized_check.json'
MEASUREMENTS = HERE / 'results/takahashi_measurements.json'
SHAPE_PRODUCT = HERE / 'results/mcsnow_shape_diagnosis.json'
CACHE = ROOT / 'research/runs/snow_agent/data_sized_check'
NOMINAL = tuple(sorted(lab.FITS))
WINDOW_C = .5           # crystals within this of a nominal temperature are its repeats
SIGMA_FACTOR = 2.       # about 95% of repeat crystals lie within two standard deviations
SPEED_TOLERANCE = .03   # m/s, the paper's stated fall-speed accuracy
QUANTITIES = ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m')
MEASURED = {'mass_kg': 'mass_kg', 'equatorial_dimension_m': 'a_full_mm', 'polar_dimension_m': 'c_full_mm'}
CALIBRATED = (-5.3, -14.4)
COMPONENT_DT = '0p0078125'
COARSE_DT = '0p03125'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(value):
    return format(value, '.8g').replace('-', 'm').replace('.', 'p')


def nominal(temperature_C):
    best = min(NOMINAL, key=lambda n: abs(n - temperature_C))
    return best if abs(best - temperature_C) <= WINDOW_C + 1e-9 else None


def crystals(product=None):
    """Measured crystals with SI values; excluded observables stay absent."""
    product = product or json.loads(MEASUREMENTS.read_text())
    out = []
    for row in product['rows']:
        excluded = set(row['calibration_excluded_observables'])
        values = {}
        for quantity, field in MEASURED.items():
            v = row.get(field)
            if v is None or field in excluded:
                continue
            values[quantity] = v if quantity == 'mass_kg' else v * 1e-3
        out.append(dict(excel_row=row['excel_row'], temperature_C=row['temperature_C'],
                        nominal_C=nominal(row['temperature_C']), age_min=row['age_min'], measured=values))
    return out


def scatter(items):
    """Pooled within-group standard deviation of ln values, by quantity.

    Groups are crystals of one nominal temperature and one age.  Deviations are
    taken from each group's own mean, so the variance estimate is unbiased with
    n - 1 degrees of freedom per group.
    """
    groups = collections.defaultdict(list)
    for c in items:
        if c['nominal_C'] is not None:
            groups[(c['nominal_C'], c['age_min'])].append(c)
    result = {}
    for q in QUANTITIES:
        total = [0., 0]
        per = collections.defaultdict(lambda: [0., 0])
        for (temperature, _age), members in groups.items():
            logs = [math.log(m['measured'][q]) for m in members if q in m['measured']]
            if len(logs) < 2:
                continue
            mean = sum(logs) / len(logs)
            ss = sum((x - mean) ** 2 for x in logs)
            total[0] += ss; total[1] += len(logs) - 1
            per[temperature][0] += ss; per[temperature][1] += len(logs) - 1
        result[q] = dict(pooled_sd_ln=math.sqrt(total[0] / total[1]), degrees_of_freedom=total[1],
                         by_temperature={format(t, 'g'): dict(sd_ln=math.sqrt(v[0] / v[1]), degrees_of_freedom=v[1])
                                         for t, v in sorted(per.items()) if v[1] > 0})
    result['replicate_groups'] = sorted([dict(nominal_C=t, age_min=a, crystals=len(m))
                                         for (t, a), m in groups.items() if len(m) > 1],
                                        key=lambda g: (g['nominal_C'], g['age_min']))
    return result


def verdict(model, measured, sd_ln):
    ratio = model / measured
    z = math.log(ratio) / sd_ln
    return dict(model=model, measured=measured, ratio=ratio, z=z, passed=abs(z) <= SIGMA_FACTOR)


def speed_verdict(model, temperature, age):
    reference = lab.fall_speed(temperature, age)
    difference = model - reference
    return dict(model=model, reference=reference, difference=difference, ratio=model / reference,
                passed=abs(difference) <= SPEED_TOLERANCE + 1e-12)


def speed_defined(temperature, age):
    return 3. <= age <= lab.FITS[temperature][0]


# --- the improved McSnow component: its retained native traces -------------------

def component_trace(temperature, dt_key=COMPONENT_DT, manifest=None):
    name = f'shape_trajectory_T{key(temperature)}_sphere_dt{dt_key}'
    manifest = manifest or {m['name']: m for m in json.loads(SHAPE_PRODUCT.read_text())['manifest']}
    entry = manifest[name]
    path = ROOT / entry['path'] / 'growth_trace.csv'
    digest = sha(path)
    if digest != entry['file_sha256']['growth_trace.csv']:
        raise ValueError('Retained trace changed: ' + name)
    with path.open() as f:
        rows = list(csv.DictReader(f))
    columns = {k: np.array([float(r[k]) for r in rows]) for k in
               ('state_elapsed_s', 'mass_kg', 'equatorial_dimension_m', 'polar_dimension_m', 'terminal_speed_m_s')}
    return dict(name=name, sha256=digest, **columns)


def at_age(trace, age_min):
    t = trace['state_elapsed_s']
    seconds = age_min * 60
    if not t[0] <= seconds <= t[-1] + 1e-9:
        return None
    return {q: float(np.interp(seconds, t, trace[q])) for q in (*QUANTITIES, 'terminal_speed_m_s')}


# --- the two-temperature calibration: its frozen executable ----------------------

def calibration_values(ages_by_temperature, dt=1 / 128, label='dt128'):
    """Run the unchanged calibrated executable from the 4 micrometre sphere."""
    from atmosphere.snow import mcsnow_column_density as column
    from atmosphere.snow import mcsnow_joint_acceptance as acceptance
    from atmosphere.snow import mcsnow_science_check as science
    column.verify_build()
    baseline = acceptance.load_baseline()
    out, parity = {}, []
    for temperature, ages in ages_by_temperature.items():
        fit = (next(f for f in baseline['searches'] if f['name'] == baseline['preferred_search'])
               if temperature == -5.3 else baseline['previous_cold'])
        controls = list(fit['physical_controls'])
        if temperature == -14.4:
            controls += [0., 1.]  # neutral warm-only density controls, as in the frozen baseline
        entry = next(e for e in column.tuning.cases()
                     if e['temperature_C'] == temperature and e['initialization'] == 'sphere')
        c = entry['configuration']
        settings = [c['temperature_K'], c['pressure_Pa'], 1., c['initial_mass_kg'],
                    c['initial_equatorial_m'], c['initial_polar_m'], dt]
        ages = sorted(ages)
        chunks = [ages[i:i + 8] for i in range(0, len(ages), 8)]
        cases = [dict(settings=settings, times_s=[a * 60 for a in chunk]) for chunk in chunks]
        rows = science.native_cases(cases, controls, CACHE / f'calibration_{key(temperature)}_{label}')
        out[temperature] = {a: dict(zip((*QUANTITIES, 'terminal_speed_m_s'), r[:4]))
                            for a, r in zip(ages, rows, strict=True)}
        if dt == 1 / 128:
            selected = sorted((r for r in fit['assessment']['rows'] if r['initialization'] == 'sphere'
                               and r['age_min'] in (3., 5., 10., 15., 20., 30.)), key=lambda r: r['age_min'])
            check = science.native_cases([dict(settings=settings, times_s=[r['age_min'] * 60 for r in selected])],
                                         controls, CACHE / f'calibration_{key(temperature)}_parity')
            difference = max(abs(v / t - 1) for row, target in zip(check, selected, strict=True)
                             for v, t in zip(row[:4], target['values'], strict=True))
            parity.append(dict(temperature_C=temperature, states=len(selected),
                               maximum_relative_difference=difference, passed=difference <= 1e-12))
    return out, parity, dict(executable_sha256=sha(column.CACHE / 'build/tuner'),
                             build_sha256=sha(column.CACHE / 'build/build.json'),
                             baseline_commit=acceptance.BASELINE_COMMIT, baseline_sha256=acceptance.BASELINE_SHA256)


# --- assessment -------------------------------------------------------------------

def assess(model_values, items, sd, temperatures, mapping='axes'):
    """Per-crystal verdicts; model_values[temperature][age] -> state or None."""
    records = []
    scale = 1. if mapping == 'axes' else transport.HEX_TO_SPHEROID ** (1 / 3)
    for c in items:
        t = c['nominal_C']
        if t not in temperatures:
            continue
        state = model_values[t].get(c['age_min'])
        if state is None:
            records.append(dict(excel_row=c['excel_row'], nominal_C=t, temperature_C=c['temperature_C'],
                                age_min=c['age_min'], compared=False, reason='age outside the model trajectory'))
            continue
        quantities = {}
        for q, measured in c['measured'].items():
            value = state[q] if q == 'mass_kg' else state[q] / scale
            quantities[q] = verdict(value, measured, sd[q])
        records.append(dict(excel_row=c['excel_row'], nominal_C=t, temperature_C=c['temperature_C'],
                            age_min=c['age_min'], compared=True, quantities=quantities))
    return records


def speed_records(model_values, items, temperatures):
    """One fall-speed comparison per temperature and age at which crystals exist."""
    pairs = sorted({(c['nominal_C'], c['age_min']) for c in items if c['nominal_C'] in temperatures})
    out = []
    for t, age in pairs:
        state = model_values[t].get(age)
        if state is None or not speed_defined(t, age):
            continue
        out.append(dict(nominal_C=t, age_min=age, **speed_verdict(state['terminal_speed_m_s'], t, age)))
    return out


def summarize(records, speeds, sd):
    expected = math.erfc(SIGMA_FACTOR / math.sqrt(2))
    rows = []
    for t in sorted({r['nominal_C'] for r in records}):
        for q in QUANTITIES:
            vs = [r['quantities'][q] for r in records if r['nominal_C'] == t and r['compared'] and q in r['quantities']]
            if not vs:
                continue
            logs = [math.log(v['ratio']) for v in vs]
            bias = sum(logs) / len(logs)
            rows.append(dict(nominal_C=t, quantity=q, crystals=len(vs),
                             outside=sum(not v['passed'] for v in vs),
                             expected_outside_by_chance=expected * len(vs),
                             mean_ln_ratio=bias, median_ratio=float(np.median([v['ratio'] for v in vs])),
                             bias_in_sd=bias / sd[q], bias_standard_error_in_sd=1 / math.sqrt(len(vs))))
        ss = [s for s in speeds if s['nominal_C'] == t]
        if ss:
            rows.append(dict(nominal_C=t, quantity='terminal_speed_m_s', ages=len(ss),
                             outside=sum(not s['passed'] for s in ss),
                             largest_difference_m_s=max(ss, key=lambda s: abs(s['difference']))['difference'],
                             largest_ratio_error=max(abs(s['ratio'] - 1) for s in ss)))
    return rows


def outside(records, speeds):
    lines = []
    for r in records:
        if not r['compared']:
            continue
        for q, v in r['quantities'].items():
            if not v['passed']:
                lines.append(dict(excel_row=r['excel_row'], temperature_C=r['temperature_C'], nominal_C=r['nominal_C'],
                                  age_min=r['age_min'], quantity=q, ratio=v['ratio'], z=v['z']))
    for s in speeds:
        if not s['passed']:
            lines.append(dict(nominal_C=s['nominal_C'], age_min=s['age_min'], quantity='terminal_speed_m_s',
                              difference_m_s=s['difference'], ratio=s['ratio']))
    return lines


def band(summary):
    """The component's systematic error: median model/measured ratio by quantity."""
    out = {}
    for q in (*QUANTITIES,):
        ratios = [(r['nominal_C'], r['median_ratio']) for r in summary if r['quantity'] == q]
        if ratios:
            lo = min(ratios, key=lambda x: x[1]); hi = max(ratios, key=lambda x: x[1])
            out[q] = dict(lowest=dict(nominal_C=lo[0], median_ratio=lo[1]), highest=dict(nominal_C=hi[0], median_ratio=hi[1]),
                          by_temperature={format(t, 'g'): v for t, v in ratios})
    return out


def run():
    measurement_product = json.loads(MEASUREMENTS.read_text())
    items = crystals(measurement_product)
    sd_info = scatter(items)
    sd = {q: sd_info[q]['pooled_sd_ln'] for q in QUANTITIES}
    assigned = [c for c in items if c['nominal_C'] is not None]
    ages = collections.defaultdict(set)
    for c in assigned:
        ages[c['nominal_C']].add(c['age_min'])

    # Improved component: fine and coarse retained traces.
    manifest = {m['name']: m for m in json.loads(SHAPE_PRODUCT.read_text())['manifest']}
    component, coarse, traces = {}, {}, {}
    for t in NOMINAL:
        fine = component_trace(t, COMPONENT_DT, manifest)
        rough = component_trace(t, COARSE_DT, manifest)
        traces[format(t, 'g')] = dict(fine=dict(name=fine['name'], sha256=fine['sha256']),
                                      coarse=dict(name=rough['name'], sha256=rough['sha256']),
                                      duration_min=float(fine['state_elapsed_s'][-1] / 60))
        component[t] = {a: at_age(fine, a) for a in ages[t]}
        coarse[t] = {a: at_age(rough, a) for a in ages[t]}
    component_refinement = max(abs(component[t][a][q] / coarse[t][a][q] - 1)
                               for t in NOMINAL for a in ages[t] if component[t][a] is not None
                               for q in (*QUANTITIES, 'terminal_speed_m_s'))

    # Frozen calibration: the unchanged executable at the measured ages.
    cal_ages = {t: ages[t] for t in CALIBRATED}
    calibration, parity, identity = calibration_values(cal_ages)
    finer, _, _ = calibration_values(cal_ages, dt=1 / 256, label='dt256')
    calibration_refinement = max(abs(finer[t][a][q] / calibration[t][a][q] - 1)
                                 for t in CALIBRATED for a in cal_ages[t] for q in (*QUANTITIES, 'terminal_speed_m_s'))
    if not all(p['passed'] for p in parity):
        raise ValueError('The frozen calibration did not reproduce its baseline')

    results = {}
    for name, values, temperatures in (('improved_component', component, NOMINAL),
                                       ('two_temperature_calibration', calibration, CALIBRATED)):
        records = assess(values, items, sd, temperatures)
        speeds = speed_records(values, items, temperatures)
        summary = summarize(records, speeds, sd)
        volume = assess(values, items, sd, temperatures, mapping='volume')
        results[name] = dict(records=records, speeds=speeds, summary=summary, outside=outside(records, speeds),
                             volume_mapping_outside=sum(not v['passed'] for r in volume if r['compared']
                                                        for k, v in r['quantities'].items() if k != 'mass_kg'),
                             axes_mapping_outside=sum(not v['passed'] for r in records if r['compared']
                                                      for k, v in r['quantities'].items() if k != 'mass_kg'))
        # Sensitivity: each temperature's own pooled scatter where it has at least 5 degrees of freedom.
        changed = 0
        for r in records:
            if not r['compared']:
                continue
            for q, v in r['quantities'].items():
                local = sd_info[q]['by_temperature'].get(format(r['nominal_C'], 'g'))
                if local and local['degrees_of_freedom'] >= 5:
                    changed += (abs(math.log(v['ratio'])) <= SIGMA_FACTOR * local['sd_ln']) != v['passed']
        results[name]['verdicts_changed_with_local_scatter'] = changed
    results['improved_component']['error_band'] = band(results['improved_component']['summary'])

    product = dict(
        schema='terluna.atmosphere.snow-data-sized-check/1',
        evidence=('Empirical comparison of two frozen models with individual laboratory crystals (mass and axes) '
                  'and with the published fall-speed relation; no parameter is fitted or changed.'),
        rule=dict(adopted='2026-10-24', window_C=WINDOW_C, sigma_factor=SIGMA_FACTOR,
                  speed_tolerance_m_s=SPEED_TOLERANCE, axes_mapping='model spheroid full axes against measured full axes',
                  speed_reference='Takahashi et al. (1991) Table 4, each line the average of about four sets; '
                                  'stated accuracy 3 cm/s; individual speeds are not archived',
                  expected_fraction_outside_by_chance=math.erfc(SIGMA_FACTOR / math.sqrt(2))),
        inputs=dict(measurements=str(MEASUREMENTS.relative_to(ROOT)), measurements_sha256=sha(MEASUREMENTS),
                    shape_product=str(SHAPE_PRODUCT.relative_to(ROOT)), shape_product_sha256=sha(SHAPE_PRODUCT),
                    component_traces=traces, calibration=identity),
        crystals=dict(total=len(items), assigned=len(assigned), outside_windows=len(items) - len(assigned),
                      by_nominal={format(t, 'g'): sum(c['nominal_C'] == t for c in assigned) for t in NOMINAL}),
        scatter=sd_info, results=results,
        numerical_checks=dict(calibration_parity=parity, calibration_dt_128_vs_256=calibration_refinement,
                              component_dt_32_vs_128=component_refinement),
        producer_sha256=sha(Path(__file__)),
        raw_directory=str(CACHE.relative_to(ROOT)),
        limits=['Repeat crystals within 0.5 C of a nominal temperature include that temperature spread.',
                'Pooled scatter assumes similar relative spread at every age and temperature; '
                'the per-temperature sensitivity is recorded.',
                'Thin planar thickness was inferred from mass and basal area by the authors and is absent per crystal.',
                'Fall speed is compared with fitted averages, so crystal-to-crystal speed scatter is not tested.',
                'The laboratory crystals grew at water saturation for at most 30 minutes at Earth gravity and pressure.'])
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def report(product=None):
    product = product or json.loads(PRODUCT.read_text())
    sc = product['scatter']
    print('scatter (ln sd, dof):', {q: (round(sc[q]['pooled_sd_ln'], 3), sc[q]['degrees_of_freedom']) for q in QUANTITIES})
    for name, r in product['results'].items():
        print('==', name, 'outside:', len(r['outside']), 'local-scatter verdict changes:', r['verdicts_changed_with_local_scatter'],
              'axes outside (axes/volume mapping):', r['axes_mapping_outside'], r['volume_mapping_outside'])
        for s in r['summary']:
            if s['quantity'] == 'terminal_speed_m_s':
                print(f"  {s['nominal_C']:6} speed ages={s['ages']} outside={s['outside']} "
                      f"largest diff={s['largest_difference_m_s']:+.3f} m/s")
            else:
                print(f"  {s['nominal_C']:6} {s['quantity'][:5]} n={s['crystals']:2} outside={s['outside']} "
                      f"(chance {s['expected_outside_by_chance']:.1f}) median ratio={s['median_ratio']:.3f} "
                      f"bias={s['bias_in_sd']:+.2f}sd")
    print('checks:', product['numerical_checks'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'report'))
    args = parser.parse_args()
    if args.action == 'run':
        report(run())
    else:
        report()


if __name__ == '__main__':
    main()
