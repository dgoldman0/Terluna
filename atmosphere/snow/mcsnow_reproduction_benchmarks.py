"""Reproducible candidate matrix and raw-output audit for the Welss benchmark.

The author archive is numerical reference output. Laboratory fits are evaluated
only inside their measured age intervals. This matrix does not adopt a repair.
"""
from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
import math
from pathlib import Path

import numpy as np
from atmosphere.snow import mcsnow_reproduction as model
from atmosphere.snow import takahashi_growth as lab

MODES = ('native', 'historical')
VARIANTS = (('IGF1', 1), ('IGF2plus', 2))
TEMPERATURES = (258.65, 262.65, 267.9)
REFINEMENTS = (.5, .25, .125, .0625, .03125, .015625)
PRODUCT = model.ROOT / 'atmosphere/snow/results/benchmark_reproduction.json'


def targets(cache):
    return {(v, t): model.author_data(cache, v, t)
            for v, times in [('HP', (600,)), ('IGF1', (180, 600)), ('IGF2plus', (180, 600))]
            for t in times}


def target_row(target, T):
    return next(r for r in target['rows'] if r['temperature_K'] == T)


def archived_configuration(author, prolate=True, offsets=(0, 0, 0)):
    """Reconstruct mass, V and phi; independently rounded axes are redundant."""
    values, errors = author['values'], author['rounding_half_width']
    mass, volume, phi = (values[k] + offset*errors[k] for k, offset in
                         zip(('mass_kg', 'ice_volume_m3', 'aspect_ratio'), offsets))
    equatorial = 2*(3*volume/(4*math.pi*phi))**(1/3)
    config = model.configuration(author['temperature_K'], duration=60, prolate=prolate)
    config.update(initial_mass_kg=mass, initial_equatorial_m=equatorial,
                  initial_polar_m=equatorial*phi)
    return config


def cases(cache=model.CACHE):
    """Declare every execution, including sensitivities that did not improve fit."""
    authors = targets(cache)
    entries = []

    def add(group, name, mode, config, **metadata):
        entries.append(dict(group=group, name=name, mode=mode, configuration=config, **metadata))

    for variant, growth in VARIANTS:
        for mode in MODES:
            for T in TEMPERATURES:
                for saturation in (True, False):
                    add('saturation', f'{variant}_{mode}_T{T}_sat{int(saturation)}', mode,
                        model.configuration(T, growth=growth, native_saturation=saturation), variant=variant)
                for dt in REFINEMENTS:
                    add('refinement', f'refine_{variant}_{mode}_T{T}_dt{dt}', mode,
                        model.configuration(T, dt=dt, growth=growth), variant=variant)
                add('thermal', f'heat_{variant}_{mode}_T{T}', mode+'_heat',
                    model.configuration(T, growth=growth), variant=variant)
            for row in authors[variant, 600]['rows']:
                T = row['temperature_K']
                add('grid', f'grid_{variant}_{mode}_T{T:.3f}', mode,
                    model.configuration(T, growth=growth), variant=variant)
            for Tc in lab.FITS:
                steps = [.125, .03125]
                if growth == 2 and Tc in (-14.4, -5.3):
                    steps.extend([.0078125, .001953125])
                for dt in steps:
                    add('extension', f'extended_{mode}_igf{growth}_Tc{Tc}_dt{dt}', mode,
                        model.configuration(Tc+273.15, dt=dt, duration=1800., growth=growth),
                        variant=variant, temperature_C=Tc)
    for variant in ('HP', 'IGF1', 'IGF2plus'):
        for T in (258.65, 267.9):
            author = target_row(authors[variant, 600], T)
            for prolate in (False, True):
                add('archived_state', f'archived_state_{variant}_T{T}_interp{int(prolate)}', 'native',
                    archived_configuration(author, prolate), variant=variant)
            if variant != 'HP' and T == 258.65:
                for i, offsets in enumerate(itertools.product((-1, 1), repeat=3)):
                    add('rounding', f'rounding_{variant}_T{T}_corner{i}', 'native',
                        archived_configuration(author, offsets=offsets), variant=variant,
                        rounding_offsets=list(offsets))
    baseline = model.ROOT / 'research/runs/snow_scale/mcsnow_growth/runs/2024_igf2_sphere_T14p4_dt0p05'
    config = dict(json.loads((baseline/'run.json').read_text())['configuration'], native_saturation=False)
    add('parity', 'parity_2024_igf2_sphere_T14p4_dt0p05', 'native', config,
        baseline=str(baseline.relative_to(model.ROOT)))
    if len({entry['name'] for entry in entries}) != len(entries):
        raise ValueError('Duplicate matrix execution')
    return entries, authors


def run_matrix(cache=model.CACHE):
    entries, _ = cases(cache)
    for mode in MODES:
        for heat in (False, True):
            model.build(cache / ('build_'+mode+('_heat' if heat else '')),
                        historical=mode == 'historical', leading_heat=heat)
    for i, entry in enumerate(entries):
        result = model.run(cache / ('build_'+entry['mode']), cache/'runs'/entry['name'], entry['configuration'])
        if not result['complete'] or (i+1) % 100 == 0:
            print(i+1, len(entries), entry['name'], 'complete', result['complete'], flush=True)


def endpoint(record, elapsed):
    return next((r for r in record['snapshots'] if r['elapsed_s'] == elapsed), None)


def endpoint_change(coarse, fine):
    return {k: (fine[k]-coarse[k])/fine[k] for k in model.QUANTITIES}


def assemble(cache=model.CACHE):
    entries, authors = cases(cache)
    build_records = {mode: model.verify_build(cache/('build_'+mode))
                     for mode in ('native', 'historical', 'native_heat', 'historical_heat')}
    build_hashes = {mode: model.native.sha(cache/('build_'+mode)/'build.json') for mode in build_records}
    records, manifest = {}, []
    for entry in entries:
        path = cache/'runs'/entry['name']
        record = model.verify_run(path)
        if record['configuration'] != entry['configuration']:
            raise ValueError('Matrix/input mismatch: '+entry['name'])
        if record['build_record_sha256'] != build_hashes[entry['mode']]:
            raise ValueError('Matrix/build mismatch: '+entry['name'])
        records[entry['name']] = record
        manifest.append(dict(**entry, path=str(path.relative_to(model.ROOT)),
                             run_record_sha256=model.native.sha(path/'run.json'),
                             complete=record['complete'], returncode=record['returncode'],
                             validation_errors=record['validation_errors']))
    comparisons, extensions, archived_states = [], [], []
    for entry in entries:
        record, config = records[entry['name']], entry['configuration']
        common = dict(case=entry['name'], mode=entry['mode'], group=entry['group'])
        if entry['group'] in ('grid', 'saturation', 'refinement', 'thermal'):
            for duration in (180, 600):
                row = endpoint(record, duration)
                target = authors[entry['variant'], duration]
                author = target_row(target, config['temperature_K'])
                comparisons.append(dict(**common, variant=entry['variant'], elapsed_s=duration,
                    temperature_K=config['temperature_K'], dt_s=config['dt_s'],
                    native_saturation=config['native_saturation'],
                    comparison=model.compare(row, author) if row else None))
        elif entry['group'] == 'extension':
            points = []
            for row in record['snapshots']:
                age, Tc = row['elapsed_s']/60, entry['temperature_C']
                observed = lab.observed(Tc, age) if 3 <= age <= lab.FITS[Tc][0] else None
                point = dict(age_min=age, state={k: row[k] for k in (*model.QUANTITIES, 'density_kg_m3')},
                             laboratory_fit=observed)
                if observed:
                    point['laboratory_fall_speed_m_s'] = lab.fall_speed(Tc, age)
                    point['fit_ratios'] = {k: row[k]/observed[k] for k in
                                          ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m', 'diameter_m')}
                    point['fit_ratios']['terminal_speed_m_s'] = row['terminal_speed_m_s']/lab.fall_speed(Tc, age)
                points.append(point)
            extensions.append(dict(**common, variant=entry['variant'], temperature_C=entry['temperature_C'],
                                   dt_s=config['dt_s'], complete=record['complete'], points=points))
        elif entry['group'] in ('archived_state', 'rounding'):
            author = target_row(authors[entry['variant'], 600], config['temperature_K'])
            row = endpoint(record, 0)
            archived_states.append(dict(**common, variant=entry['variant'], temperature_K=config['temperature_K'],
                prolate_interpolation=config['prolate_interpolation'],
                rounding_offsets=entry.get('rounding_offsets', [0, 0, 0]),
                author_speed_m_s=author['values']['terminal_speed_m_s'],
                author_speed_half_rounding=author['rounding_half_width']['terminal_speed_m_s'],
                reconstructed_speed_m_s=row['terminal_speed_m_s'] if row else None,
                speed_ratio=row['terminal_speed_m_s']/author['values']['terminal_speed_m_s'] if row else None))
    grid_summary = []
    for variant, _ in VARIANTS:
        for mode in MODES:
            for duration in (180, 600):
                selected = [r for r in comparisons if r['group']=='grid' and r['variant']==variant and
                            r['mode']==mode and r['elapsed_s']==duration and r['temperature_K'] < 273.15]
                valid = [r for r in selected if r['comparison'] is not None]
                grid_summary.append(dict(variant=variant, mode=mode, elapsed_s=duration,
                    growing_temperatures=len(selected), available_endpoints=len(valid),
                    matches_printed_precision=sum(r['comparison']['all_within_printed_rounding'] for r in valid),
                    absolute_relative_error={k: dict(median=float(np.median([abs(r['comparison']['relative_error'][k]) for r in valid])),
                        maximum=max(abs(r['comparison']['relative_error'][k]) for r in valid)) for k in model.QUANTITIES}))
    refinement = []
    for variant, _ in VARIANTS:
        for mode in MODES:
            for T in TEMPERATURES:
                names = [f'grid_{variant}_{mode}_T{T:.3f}'] + [f'refine_{variant}_{mode}_T{T}_dt{dt}' for dt in REFINEMENTS]
                last, previous = [endpoint(records[n], 600) for n in (names[-1], names[-2])]
                author = target_row(authors[variant, 600], T)
                residual = model.compare(last, author)
                refinement.append(dict(variant=variant, mode=mode, temperature_K=T, cases=names,
                    final_pair_relative_change=endpoint_change(previous, last),
                    final_comparison=residual,
                    outside_rounding_plus_last_change={k: abs(last[k]-author['values'][k]) >
                        author['rounding_half_width'][k]+abs(last[k]-previous[k]) for k in model.QUANTITIES},
                    criterion='Last-pair change is a numerical sensitivity, not a rigorous error bound.'))
    extension_refinement = []
    for variant, _ in VARIANTS:
        for mode in MODES:
            for Tc in lab.FITS:
                selected = sorted((r for r in extensions if r['variant']==variant and r['mode']==mode and
                                   r['temperature_C']==Tc), key=lambda r: -r['dt_s'])
                pairs = []
                for coarse, fine in zip(selected, selected[1:]):
                    a, b = endpoint(records[coarse['case']],1800), endpoint(records[fine['case']],1800)
                    pairs.append(dict(coarse_dt_s=coarse['dt_s'], fine_dt_s=fine['dt_s'],
                                      relative_change=endpoint_change(a, b) if a and b else None))
                extension_refinement.append(dict(variant=variant, mode=mode, temperature_C=Tc, pairs=pairs))
    parity_entry = next(e for e in entries if e['group']=='parity')
    parity_dir = cache/'runs'/parity_entry['name']
    baseline = model.ROOT/parity_entry['baseline']
    parity = {name: dict(baseline_sha256=model.native.sha(baseline/name),
                        candidate_sha256=model.native.sha(parity_dir/name),
                        byte_identical=(baseline/name).read_bytes()==(parity_dir/name).read_bytes())
              for name in ('trajectory.csv','air.csv','branch_diagnostic.csv','step_diagnostic.csv')}
    artifacts = ('archive_record.json','public_source_tree.json','deposition_history.json')
    sources = {name: dict(path=str((cache/name).relative_to(model.ROOT)), sha256=model.native.sha(cache/name))
               for name in artifacts}
    return dict(schema='terluna.atmosphere.mcsnow-benchmark-reproduction/1',
        producer={str(Path(p).relative_to(model.ROOT)): model.native.sha(p)
                  for p in (__file__, model.__file__, lab.__file__, model.native.__file__)},
        source_record='atmosphere/snow/benchmark_reproduction_sources.json',
        source_artifacts=sources, builds={k: dict(record_sha256=build_hashes[k], **v) for k,v in build_records.items()},
        author_targets={f'{v}_{t}': value for (v,t), value in authors.items()},
        policy=dict(author_match='All seven compared saved quantities inside their printed half-rounding interval.',
                    refinement='Compare last two steps separately; no asymptotic-convergence claim.',
                    extension='Unmatched candidate extension, not a reproduced published configuration.',
                    empirical='Takahashi Table 2 and 4 fits, no extrapolation beyond their recorded age window.'),
        counts=dict(executions=len(entries), complete=sum(r['complete'] for r in records.values()),
                    by_group=dict(Counter(e['group'] for e in entries))),
        parity=parity, grid_summary=grid_summary, refinements=refinement,
        extension_refinements=extension_refinement, comparisons=comparisons,
        archived_state_controls=archived_states, extensions=extensions, manifest=manifest,
        evidence='Native source controls and numerical archive comparison; exact paper reproduction remains unresolved; no adopted lunar closure.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run', 'assemble'])
    parser.add_argument('--cache', type=Path, default=model.CACHE)
    parser.add_argument('--output', type=Path, default=PRODUCT)
    args = parser.parse_args()
    if args.action == 'run':
        run_matrix(args.cache)
    product = assemble(args.cache)
    model.save(args.output, product)
    print(json.dumps(product['counts']))


if __name__ == '__main__':
    main()
