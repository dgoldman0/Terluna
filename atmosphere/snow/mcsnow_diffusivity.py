"""McSnow's vapour diffusivity in the lunar builds, with its 1 atm reference restored, and what the old factor moved.

McSnow diffuses vapour by Pruppacher and Klett's fit, Dv = 0.211e-4 (T / 273.15 K)^1.94 ps / p (mo_atmo.f90's
diffusivity_vapour, and the same line in mo_cond.f90), where ps is its constant surface pressure, 101325 Pa. The lunar
builds of mcsnow_aggregation.py set ps to the freezing level's pressure, 65925 Pa, to anchor their column, and that
also set the fit's reference: every lunar run diffused vapour at 0.651 of the fit.  This module gives the fit its own
constant, p_dv = 101325 Pa, in the lunar build only, checks that the sources then differ from the old lunar build in
those lines alone and that the Earth build is untouched, and reruns a minimal set of the growth-and-aggregation cases
with their recorded configurations, to compare with the old runs at the same statistics.

    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_diffusivity build
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_diffusivity run --case col_q3_f0.01_base
    OPENBLAS_NUM_THREADS=1 python -B -m atmosphere.snow.mcsnow_diffusivity compare   # results/mcsnow_diffusivity.json

The frozen component's physics is otherwise unchanged; runs go to research/runs/snow_agent_3/dv_fix.
"""
from __future__ import annotations

import argparse
import contextlib
import difflib
import hashlib
import json
from pathlib import Path

from atmosphere.snow import mcsnow_growth_aggregation as growagg

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CACHE = ROOT / 'research/runs/snow_agent_3/dv_fix'
OLD = growagg.CACHE
PRODUCT = HERE / 'results/mcsnow_diffusivity.json'
SCHEMA = 'terluna.atmosphere.snow-mcsnow-diffusivity/1'
REFERENCE_PA = 101325.
# The minimal set: water-limited columns at 3 and 0.3 mm/h with the reference crystal flux, 3 mm/h with a hundredth of
# it, one column held at saturation, one storm-path run in the archived humidity, and the deep light-snow column last.
CASES = ('col_q3_f0.01_base', 'col_sat_f0.01_base', 'path_archived_c0.1_base', 'col_q3_f1_base', 'col_q0.3_f1_base',
         'col_deep_q0.3_f1_base')
_original_sources = growagg.patched_sources


def _replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Expected one occurrence of: ' + old[:70])
    return text.replace(old, new)


def fixed_sources(world):
    """growagg's patched sources, with the diffusivity fit's reference pressure its own constant in the lunar build."""
    out = dict(_original_sources(world))
    if world != 'moon':
        return out
    types = out['mo_atmo_types.f90']
    ps_line = next(line for line in types.splitlines() if line.lstrip().startswith('ps    ='))
    out['mo_atmo_types.f90'] = _replace_once(
        types, ps_line + '\n', ps_line + '\n'
        f'    p_dv  = {REFERENCE_PA:.0f}._wp,     &  ! [Pa]    reference pressure of the vapour diffusivity fit (Terluna)\n')
    atmo = _replace_once(out['mo_atmo.f90'], 'B_w, ps, rhlr, e_3, T_3, Rd, tlr, Rv, T_3i, grav',
                         'B_w, ps, rhlr, e_3, T_3, Rd, tlr, Rv, T_3i, grav, p_dv')
    out['mo_atmo.f90'] = _replace_once(atmo, 'Dv = 0.211e-4_wp * (T * T_3i)**1.94_wp * ps / p',
                                       'Dv = 0.211e-4_wp * (T * T_3i)**1.94_wp * p_dv / p')
    cond = out.get('mo_cond.f90') or (growagg.SOURCE / 'src/mo_cond.f90').read_text()
    cond = _replace_once(cond, 'USE mo_atmo_types,   ONLY: T_3i, ps, rhol, rholi, t_atmo, Rv',
                         'USE mo_atmo_types,   ONLY: T_3i, ps, rhol, rholi, t_atmo, Rv, p_dv')
    out['mo_cond.f90'] = _replace_once(cond, 'Dv = 0.211e-4_wp * (atmo%T*T_3i)**1.94_wp * ps / atmo%p',
                                       'Dv = 0.211e-4_wp * (atmo%T*T_3i)**1.94_wp * p_dv / atmo%p')
    return out


@contextlib.contextmanager
def corrected():
    """growagg's build and run machinery on the corrected sources, in this module's cache."""
    saved = growagg.patched_sources, growagg.CACHE
    growagg.patched_sources, growagg.CACHE = fixed_sources, CACHE
    try:
        yield
    finally:
        growagg.patched_sources, growagg.CACHE = saved


def source_differences():
    """Lines in which the corrected lunar sources differ from the old lunar build's, and whether the Earth sources
    are the old ones."""
    old, new = _original_sources('moon'), fixed_sources('moon')
    changed = []
    for name in sorted(set(old) | set(new)):
        a = old.get(name, (growagg.SOURCE / 'src' / name).read_text())
        b = new[name]
        for line in difflib.unified_diff(a.splitlines(), b.splitlines(), lineterm='', n=0):
            if line.startswith(('+', '-')) and not line.startswith(('+++', '---')):
                changed.append(f'{name}: {line}')
    earth_same = fixed_sources('earth') == _original_sources('earth')
    allowed = ('p_dv',)
    only_dv = all(any(word in line for word in allowed) or 'ps / p' in line or 'ps / atmo%p' in line
                  or 'ONLY: T_3i, ps, rhol' in line or 'B_w, ps, rhlr' in line for line in changed)
    return dict(changed_lines=changed, earth_sources_unchanged=earth_same, lunar_changes_only_the_diffusivity=only_dv)


def build():
    with corrected():
        record = growagg.build('moon')
        path = growagg.build_path('moon')
    old = growagg.verify_build('moon')
    old_earth = growagg.verify_build('earth')
    return dict(corrected=dict(path=str(path.relative_to(ROOT)), executable_sha256=record['executable_sha256']),
                old_moon_executable_sha256=old['executable_sha256'],
                earth_executable_sha256=old_earth['executable_sha256'], **source_differences())


def run_case(name):
    """The old run's recorded configuration, with the corrected build."""
    old_dir = OLD / 'runs' / name
    record = json.loads((old_dir / 'run.json').read_text())
    link = old_dir / 'box_forcing.bin'
    forcing_link = link.resolve() if link.is_symlink() else None
    with corrected():
        out = growagg.run(name, 'moon', record['configuration'], seed=record['seed'], forcing_link=forcing_link)
    return dict(name=name, returncode=out['returncode'], wall_s=out['wall_s'], executable_sha256=out['executable_sha256'])


def _column_stats(name, cache):
    saved = growagg.CACHE
    growagg.CACHE = cache
    try:
        s = growagg.column_summary(name, pool_h=growagg.POOL_H)
    finally:
        growagg.CACHE = saved
    a = s['arriving']
    return dict(time_h=s['time_h'], onset_depth_km=s['onset_depth_km'],
                arriving={k: a.get(k) for k in ('water_mm_h', 'd50_mm', 'd90_mm', 'd99_mm', 'largest_mm', 'melted_d50_mm')},
                d90_by_depth_mm=[[b['depth_below_top_km'], round(b['d90_mm'], 3)] for b in s['bins']])


def _path_stats(name, cache):
    from atmosphere.snow import melting
    saved = growagg.CACHE
    growagg.CACHE = cache
    try:
        s = growagg.path_summary(name)
        budget = melting.path_budget(name)
    finally:
        growagg.CACHE = saved
    return dict(columns_with_arrivals=s['columns_with_arrivals'], arriving_d90_mm=s['arriving_d90_mm'],
                pooled={k: s['pooled'].get(k) for k in ('arrivals', 'd50_mm', 'd90_mm', 'd99_mm', 'largest_mm')},
                fraction_arriving_d90_over=s['fraction_arriving_d90_over'],
                budget={k: budget[k] for k in ('sublimated_share_of_grown', 'delivered_share_of_grown',
                                               'censored_share_of_grown', 'aloft_share_of_grown', 'grown_kg')})


MELTING_CASES = ('col_q0.3_f1_base', 'col_q3_f0.01_base', 'col_sat_f0.01_base', 'col_deep_q0.3_f1_base')


def _melting_stats(name, cache, tables):
    """The melting budget of a column's arriving flakes in the median composite columns (melting.py's saved response
    tables), for the old or the corrected run."""
    from atmosphere.snow import melting
    saved = growagg.CACHE
    growagg.CACHE = cache
    try:
        m, w, _ = melting.mcsnow_population(name)
    finally:
        growagg.CACHE = saved
    out = {}
    for box, (table, air, zf) in tables.items():
        pop = melting.population(table, m, w)
        lay = melting.layers(pop, air, zf)
        out[box] = dict(depths_ice_gone_m=lay['depths_ice_gone_m'], totals=lay['totals'],
                        ground=melting.water_by_ground(pop, zf))
    return out


def _melting_tables():
    import numpy as np
    from atmosphere.snow import melting, melting_environment as env
    tables = {}
    for box in ('land', 'sea'):
        data, _ = env.load(box)
        air, zf = melting.composite_air(data)
        with np.load(melting.RAW / f'table_{box}_median_central.npz') as z:
            stack = {k: z[k] for k in ('ice', 'water', 'sub', 'melt', 'evap', 'time')}
            spans, grid = z['spans'], z['grid']
        tables[box] = (dict(spans=spans, masses=melting.aggregate_mass(spans), grid=grid, stack=stack), air, zf)
    return tables


def compare():
    rows = {}
    tables = None
    for name in CASES:
        if not (CACHE / 'runs' / name / 'run.json').exists():
            continue
        stats = _path_stats if name.startswith('path_') else _column_stats
        old, new = stats(name, OLD), stats(name, CACHE)
        rows[name] = dict(old=old, corrected=new)
        if name in MELTING_CASES:
            tables = tables or _melting_tables()
            rows[name]['melting'] = dict(old=_melting_stats(name, OLD, tables),
                                         corrected=_melting_stats(name, CACHE, tables))
    product = dict(schema=SCHEMA,
                   evidence=('McSnow reruns of recorded growth-and-aggregation configurations with the vapour '
                             'diffusivity fit given its 1 atm reference in the lunar build, against the old runs at the '
                             'same statistics; numerical, no lunar validation.'),
                   reading_rule=('Columns: sizes by mass flux pooled over the last 12 hours, depths below the top where '
                                 'the 90th percentile passes 10 and 14 cm; paths: arrivals at 0 C and the column budgets '
                                 'of melting.path_budget.'),
                   reference_pressure_pa=REFERENCE_PA, old_factor=growagg.agg.lunar_column()['ps'] / REFERENCE_PA,
                   build=json.loads((CACHE / 'build.json').read_text()) if (CACHE / 'build.json').exists() else None,
                   cases=rows,
                   producer=dict(domain='atmosphere', files={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
                                                             for p in ('atmosphere/snow/mcsnow_diffusivity.py',
                                                                       'atmosphere/snow/mcsnow_growth_aggregation.py',
                                                                       'atmosphere/snow/mcsnow_aggregation.py')}))
    PRODUCT.write_text(json.dumps(product, indent=1) + '\n')
    return product


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('diff', 'build', 'run', 'compare'))
    parser.add_argument('--case')
    args = parser.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    if args.action == 'diff':
        print(json.dumps(source_differences(), indent=1))
    elif args.action == 'build':
        out = build()
        (CACHE / 'build.json').write_text(json.dumps(out, indent=1) + '\n')
        print(json.dumps(out, indent=1))
    elif args.action == 'run':
        print(json.dumps(run_case(args.case)), flush=True)
    elif args.action == 'compare':
        product = compare()
        print(len(product['cases']), 'cases compared')


if __name__ == '__main__':
    main()
