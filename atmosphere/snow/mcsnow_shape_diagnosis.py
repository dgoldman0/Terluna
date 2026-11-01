"""Separate growth-rate and shape errors against identified laboratory fits.

Native traced McSnow is reused without another physical intervention. Equal-mass
comparisons are diagnostics of a rate-only repair, not empirical validation.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from atmosphere.snow import mcsnow_growth_trace as trace
from atmosphere.snow import takahashi_growth as lab
from atmosphere.snow import crystal_transport_replay as transport

native = trace.native
model = trace.model
HERE = Path(__file__).resolve().parent
PRODUCT = HERE / 'results/mcsnow_shape_diagnosis.json'
CACHE = trace.CACHE / 'shape_diagnosis'
MAPPINGS = ('axes', 'volume')
TIMESTEPS = (1 / 32, 1 / 128)
FIELDS = ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m', 'terminal_speed_m_s')


def key(value):
    return format(value, '.8g').replace('-', 'm').replace('.', 'p')


def exponent(segments, age):
    """The fit's selected smooth branch; joins are handled separately."""
    return next((p for upper, _, p in segments if age < upper), segments[-1][2])


def fit_state(temperature, age, mapping):
    obs = lab.observed(temperature, age)
    if mapping not in MAPPINGS:
        raise ValueError('Unknown observation-to-spheroid mapping')
    factor = 1. if mapping == 'axes' else transport.HEX_TO_SPHEROID
    scale = factor ** (1 / 3)
    a, c, mass = (obs[k] for k in ('equatorial_dimension_m', 'polar_dimension_m', 'mass_kg'))
    volume = factor * math.pi / 6 * a * a * c
    density = mass / volume
    return dict(**obs, mapping=mapping, volume_factor=factor, axis_scale=scale,
                model_a_m=scale*a, model_c_m=scale*c, volume_m3=volume,
                density_kg_m3=density, fall_speed_m_s=lab.fall_speed(temperature, age))


def fit_rates(temperature, age, mapping):
    state = fit_state(temperature, age, mapping)
    curves = lab.FITS[temperature][2:]
    joins = {float(s[0]) for curve in curves for s in curve[:-1]}
    if any(math.isclose(age, j, rel_tol=0, abs_tol=1e-10) for j in joins):
        raise ValueError('No derivative at a separately fitted join')
    b, a, c = [exponent(curve, age) for curve in curves]
    seconds = age * 60
    return dict(**state, mass_exponent=b, a_exponent=a, c_exponent=c,
                mass_rate_kg_s=state['mass_kg']*b/seconds,
                a_log_rate_s=a/seconds, c_log_rate_s=c/seconds,
                volume_log_rate_s=(2*a+c)/seconds,
                effective_growth_ratio=c/a,
                deposition_density_kg_m3=state['density_kg_m3']*b/(2*a+c))


def fit_joins(temperature):
    """Record discontinuities rather than differentiating or smoothing them."""
    rows = []
    for quantity, curve in zip(FIELDS[:3], lab.FITS[temperature][2:]):
        for left, right in zip(curve, curve[1:]):
            age = left[0]
            lo, hi = left[1]*age**left[2]*1e-3, right[1]*age**right[2]*1e-3
            rows.append(dict(quantity=quantity, age_min=age, left_limit=lo,
                             right_limit=hi, right_over_left=hi/lo))
    return rows


def sample_ages(temperature):
    end = lab.FITS[temperature][0]
    return sorted({3., float(end), *[float(t) for t in (5, 10, 15, 20, 30) if t <= end]})


def cases(ice_density):
    for temperature, (end, *_rest) in lab.FITS.items():
        for initialization in ('sphere', 'axes', 'volume'):
            state = None if initialization == 'sphere' else fit_state(temperature, 3., initialization)
            for dt in TIMESTEPS:
                config = model.configuration(273.15+temperature, dt=dt,
                                             duration=(end-(3 if state else 0))*60)
                config['output_interval_s'] = .25
                if state:
                    config.update(initial_mass_kg=state['mass_kg'],
                                  initial_equatorial_m=state['model_a_m'],
                                  initial_polar_m=state['model_c_m'])
                yield dict(name=f'shape_trajectory_T{key(temperature)}_{initialization}_dt{key(dt)}',
                           temperature_C=temperature, initialization=initialization,
                           dt_s=dt, initial_age_min=3. if state else 0.,
                           axis_scale=state['axis_scale'] if state else 1.,
                           admissible=state is None or state['density_kg_m3'] <= ice_density,
                           initial_state=state, configuration=config)


def probes(ice_density):
    for temperature in lab.FITS:
        end = lab.FITS[temperature][0]
        ages = sorted({3.5, 5.5, 8., 10., float(end)})
        for age, mapping, mode in itertools.product(ages, MAPPINGS, ('native', 'plate_base')):
            if age > end:
                continue
            state = fit_rates(temperature, age, mapping)
            config = model.configuration(273.15+temperature, dt=1/64, duration=1.)
            config.update(initial_mass_kg=state['mass_kg'],
                          initial_equatorial_m=state['model_a_m'],
                          initial_polar_m=state['model_c_m'], output_interval_s=1.)
            yield dict(name=f'shape_probe_T{key(temperature)}_age{key(age)}_{mapping}_{mode}',
                       temperature_C=temperature, age_min=age, mapping=mapping, mode=mode,
                       fit=state, admissible=state['density_kg_m3'] <= ice_density,
                       configuration=config)


def ice_density():
    path = trace.CACHE / 'traces/extended_plate_base_igf2_Tc-14.4_dt0.0078125'
    reference = trace.verify_run(path)
    return reference['air'][0]['native_ice_density_kg_m3']


def run_path(entry):
    return trace.CACHE / 'traces/shape_diagnosis' / entry['name']


def execute():
    rho = ice_density()
    for group, entries in (('trajectory', cases(rho)), ('probe', probes(rho))):
        for entry in entries:
            if not entry['admissible']:
                print(group, entry['name'], 'inadmissible fit-to-spheroid density', flush=True)
                continue
            mode = entry.get('mode', 'plate_base')
            result = trace.run(mode, 'shape_diagnosis/'+entry['name'], entry['configuration'])
            print(group, entry['name'], 'complete' if result['complete'] else 'FAILED', flush=True)


def snapshot_at(rows, elapsed):
    exact = next((r for r in rows if r['elapsed_s'] == elapsed), None)
    if exact is None:
        raise ValueError('Requested elapsed time is not a saved output')
    return exact


def ratios(snapshot, temperature, age, scale=1.):
    observed = lab.observed(temperature, age)
    observed['terminal_speed_m_s'] = lab.fall_speed(temperature, age)
    values = {k: snapshot[k] for k in FIELDS}
    for k in FIELDS[1:3]:
        values[k] /= scale
    return dict(values=values, observed=observed,
                ratio={k: values[k]/observed[k] for k in FIELDS})


def trajectory_record(entry, record):
    rows = record['snapshots']
    comparisons = []
    for age in sample_ages(entry['temperature_C']):
        elapsed = 60*(age-entry['initial_age_min'])
        if elapsed < 0 or elapsed > rows[-1]['elapsed_s']:
            continue
        same_age = snapshot_at(rows, elapsed)
        target_mass = lab.observed(entry['temperature_C'], age)['mass_kg']
        matched = trace.matched_mass(rows, target_mass)
        comparisons.append(dict(age_min=age,
            age_matched=ratios(same_age, entry['temperature_C'], age, entry['axis_scale']),
            mass_matched=(dict(**ratios(matched['state'], entry['temperature_C'], age, entry['axis_scale']),
                              model_age_min=entry['initial_age_min']+matched['state']['elapsed_s']/60,
                              interpolation_bracket_s=matched['bracket_elapsed_s']) if matched else None)))
    return dict(**{k:v for k,v in entry.items() if k != 'configuration'},
                complete=record['complete'], validation_errors=record['validation_errors'],
                comparisons=comparisons)


def probe_record(entry, record, rho):
    source = record['trace'][0]
    target = entry['fit']
    gamma = source['effective_growth_factor']
    rate = source['deposition_rate_kg_s']
    density = source['deposition_density_kg_m3']
    v_log_rate = rate/(density*target['volume_m3'])
    a_log_rate, c_log_rate = v_log_rate/(2+gamma), v_log_rate*gamma/(2+gamma)
    return dict(**{k:v for k,v in entry.items() if k != 'configuration'},
                complete=record['complete'], source_evaluation=source,
                fitted_deposition_density_exceeds_ice=target['deposition_density_kg_m3']>rho,
                ratio=dict(mass_rate=rate/target['mass_rate_kg_s'],
                           a_log_rate=a_log_rate/target['a_log_rate_s'],
                           c_log_rate=c_log_rate/target['c_log_rate_s'],
                           growth_ratio=gamma/target['effective_growth_ratio'],
                           deposition_density=density/target['deposition_density_kg_m3'],
                           speed=record['snapshots'][0]['terminal_speed_m_s']/target['fall_speed_m_s']),
                required_uptake_multiplier=dict(a=target['a_log_rate_s']/a_log_rate,
                                                c=target['c_log_rate_s']/c_log_rate,
                                                mass=target['mass_rate_kg_s']/rate))


def assemble():
    source = transport.verify_sources()
    rho = ice_density()
    runs, rejected, manifest, probe_results = [], [], [], []
    for group, entries in (('trajectory', cases(rho)), ('probe', probes(rho))):
        for entry in entries:
            if not entry['admissible']:
                rejected.append(dict(group=group, **entry,
                    reason='Fitted mass and mapped envelope imply density above native bulk ice; no clipping.'))
                continue
            path = run_path(entry)
            record = trace.verify_run(path)
            if record['configuration'] != entry['configuration'] or record['mode'] != entry.get('mode','plate_base'):
                raise ValueError('Existing run differs from the declared experiment')
            if not record['complete']:
                raise ValueError('Incomplete diagnostic retained at '+str(path))
            manifest.append(dict(name=entry['name'], group=group,
                                 path=str(path.relative_to(model.ROOT)),
                                 record_sha256=native.sha(path/'run.json'),
                                 build_record_sha256=record['build_record_sha256'],
                                 file_sha256=record['file_sha256']))
            if group == 'trajectory':
                runs.append(trajectory_record(entry, record))
            else:
                probe_results.append(probe_record(entry, record, rho))
    refinements = []
    for fine in [r for r in runs if r['dt_s']==min(TIMESTEPS)]:
        coarse = next(r for r in runs if r['temperature_C']==fine['temperature_C'] and
                      r['initialization']==fine['initialization'] and r['dt_s']==max(TIMESTEPS))
        rows = []
        for lo, hi in zip(coarse['comparisons'], fine['comparisons']):
            if lo['age_min'] != hi['age_min']:
                raise ValueError('Mismatched refinement ages')
            for comparison in ('age_matched', 'mass_matched'):
                if lo[comparison] is None or hi[comparison] is None:
                    continue
                rows.append(dict(age_min=hi['age_min'], comparison=comparison,
                                 relative_change={k:abs(lo[comparison]['values'][k]/hi[comparison]['values'][k]-1)
                                                  for k in FIELDS}))
        refinements.append(dict(temperature_C=fine['temperature_C'], initialization=fine['initialization'],
                                comparisons=rows))
    producer_paths = [Path(__file__), Path(trace.__file__), Path(model.__file__), Path(native.__file__),
                      Path(lab.__file__), Path(transport.__file__)]
    return dict(schema='terluna.atmosphere.mcsnow-shape-diagnosis/1',
                evidence='Native conditional calculations against laboratory fitted histories. No new growth law is fitted or adopted.',
                producer={str(p.relative_to(model.ROOT)):native.sha(p) for p in producer_paths},
                source_verification=source, native_ice_density_kg_m3=rho,
                observation_mapping=dict(
                    axes='Preserve measured full axes and mass; use enclosing spheroid volume. Reject densities above bulk ice.',
                    volume='Preserve laboratory hexagonal-prism envelope volume, mass and aspect ratio; multiply both model axes by (9 sqrt(3)/(4 pi))^(1/3). Invert that mapping only for reported axes. Drag uses actual model geometry.',
                    sphere='Native 4 micrometre sphere; unobserved nucleation history, retained as a declared sensitivity.'),
                fit_joins={str(t):fit_joins(t) for t in lab.FITS},
                trajectories=runs, fixed_state_probes=probe_results, excluded_states=rejected,
                numerical_refinements=refinements, manifest=manifest,
                limitations=[
                    'Laboratory fits summarize different crystals, lack pointwise uncertainty, and can jump at segment joins; rates apply within individual segments only.',
                    'Planar thickness was inferred from mass and basal area, so three fitted quantities are not independent observations.',
                    'Neither mapping proves that a spheroid represents the photographed crystal or its projected area.',
                    'Equal-mass interpolation changes model age. It diagnoses a rate-only repair and is not an age-matched prediction.',
                    'The IGF2+ shape-versus-mass argument assumes fixed temperature, air properties and geometry/drag/directional growth; it does not cover arbitrary changing weather.',
                    'Original author single-crystal generating configuration remains unrecovered. Laboratory fits and author numerical archives remain distinct targets.',
                    'No aggregate population, finite-water feedback, sublimation validation or lunar size distribution is established.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'assemble'))
    args = parser.parse_args()
    if args.action == 'run':
        execute()
    else:
        result = assemble()
        model.save(PRODUCT, result)
        print(len(result['trajectories']), 'trajectories;', len(result['fixed_state_probes']), 'fixed-state probes;',
              len(result['excluded_states']), 'explicit exclusions')


if __name__ == '__main__':
    main()
