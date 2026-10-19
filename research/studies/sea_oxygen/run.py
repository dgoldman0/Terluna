"""The held fresh seas' oxygen and mixing in their basins at the resort's 301 K.

    python -m research.studies.sea_oxygen.seas                                    # the basins (about a minute)
    climate/gcm/.venv/bin/python -m research.studies.sea_oxygen.extract           # the GCM's weather (two minutes)
    python -m research.studies.sea_oxygen.run --columns                           # the columns (hours, two cores)
    python -m research.studies.sea_oxygen.run                                     # the product (columns cached)
    python -m research.studies.sea_oxygen.run --columns --cache-only              # the index and product from cache

First it proves the 301 K construction at 300 K: the lake kinds study's Setting301 gives climate_targets' own share
and lapse rate there, this study's extract holds lake_oxygen's waters and Sun exactly, and lake_oxygen's saturation
figures for its three targets, the resort's at 300 K among them, come out of this extract unchanged. Then, at 301 K,
the lake model (biosphere/ecology/lake.py) runs on each run basin of results/seas.json with lake_oxygen's physics()
and scenarios() (columns.py):
- the central columns: each basin's own hypsography to its deepest point, at its typical cell and at its coldest
  (setting.RULE), the deep mixing from the wind's central share, from the warm start through 80 recorded passes of the
  record and then lake_oxygen's settle test at the analysis step;
- the wind share's bracket at every typical cell, and at Procellarum and Ingenii the Lake number's diffusivity and the
  first round's fixed dissipation as sensitivity, with Ingenii's respiration products left out of the density and
  Procellarum's physics that of oligotrophic water; each from its basin's central typical column, 30 passes;
- the other form of the deep mixing at Procellarum's and Ingenii's typical cells: the wind's dissipation turned into a
  diffusivity by Shih et al.'s (2005) regimes of the buoyancy Reynolds number, as Bouffard & Boegman (2013) take them
  for lakes, in place of Osborn's rule (lake_oxygen/mixing_forms.py's swap of lake.py's diffusivity() inside the
  worker's process), from the warm start as the central columns run;
- the shelves and bays: lake_oxygen's 3 m productive pan at every column's cell and a 10 m pan at each typical cell,
  and a few of them again under the first round's fixed dissipation, the bracket of the wind's share in thin water;
- the slab check against the GCM's own water at the typical cells of the four largest basins and the two polar ones.
Each column is saved daily in the lake kinds study's format. It writes research/runs/sea_oxygen/columns.json, the
index of every saved column with its sha256s, and results/sea_oxygen.json (schema terluna.research.sea-oxygen/1).
With --cache-only it computes no column, slab or pan and stops at the first one not cached, or at an index that does
not hold every current case.
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import hashlib
import inspect
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_kinds import columns as KC
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import extract as X
from research.studies.lake_oxygen import mixing_forms as MF
from research.studies.lake_oxygen import run as R
from research.studies.sea_oxygen import analysis as A
from research.studies.sea_oxygen import columns as C
from research.studies.sea_oxygen import setting as T
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'sea_oxygen.json'
INDEX = T.RUNS / 'columns.json'
SCHEMA = 'terluna.research.sea-oxygen/1'
INDEX_SCHEMA = 'terluna.research.sea-oxygen-columns/1'
DEEP_SPM = R.STEPS_PER_MEAN_DEEP       # lake_oxygen's deep lakes step every 7150 s
PAN_SPM = R.STEPS_PER_MEAN             # and its shallows every 3575 s
CENTRAL_PHYSICS = R.COUPLED['deep']    # mesotrophic water, as lake_oxygen's deep lakes
PANS_M = (3.0, 10.0)
SENSITIVITY_BASINS = ('procellarum', 'ingenii')   # the Lake number's and the fixed dissipation's rows
SLAB_BASINS = ('procellarum', 'ingenii', 'nubium', 'humboldtianum', 'drygalski', 'rozhdestvenskiy')
FORM = 'shih_bouffard'                 # the deep mixing's other form (lake_oxygen/mixing_forms.py)
FORM_BASINS = ('procellarum', 'ingenii')
PAN_BRACKET = (('procellarum', 3.0), ('procellarum', 10.0), ('ingenii', 3.0), ('ingenii', 10.0), ('7n_302e', 3.0),
               ('drygalski', 3.0))     # typical cells' pans run again under the fixed dissipation
TROPHIC = tuple(R.TROPHIC)
CENTRAL = '/renewal/supply'


def _form_code():
    """sha256 of what the other form runs: mixing_forms.py's shih_bouffard() and mixing_form() and the constants they
    read, as this process loaded them, so that an edit elsewhere in that file leaves the form's columns standing."""
    h = hashlib.sha256()
    for name in ('shih_bouffard', 'mixing_form'):
        h.update(inspect.getsource(getattr(MF, name)).encode())
    h.update(json.dumps([MF.NU, MF.RE_MOLECULAR, MF.RE_ENERGETIC]).encode())
    return h.hexdigest()[:20]


FORM_CODE = _form_code()


def diffusivity_form(case):
    """The form that turns the dissipation into a diffusivity below the mixed layer in a case."""
    if case.get('form'):
        return case['form']
    return 'lake_number' if case.get('mixing') == 'lake_number' else 'osborn'


@contextlib.contextmanager
def deep_mixing_form(case):
    """lake.py's diffusivity() swapped for the case's form inside this process (none for Osborn's rule)."""
    if case.get('form'):
        with MF.mixing_form(case['form']):
            yield
    else:
        yield


def sha256(path) -> str:
    return KC.sha256(path)


def dumps(obj) -> str:
    text = json.dumps(obj, indent=1)
    return re.sub(r'\[\s+([-0-9.eE,\s]+?|null|true|false)\s+\]', lambda m: '[' + ' '.join(m.group(1).split()) + ']',
                  text)


# --- the setting ----------------------------------------------------------------------------------------------------

_SETTINGS = {}


def setting(ground_k=S.TARGET_K):
    if ground_k not in _SETTINGS:
        _SETTINGS[ground_k] = T.SeaSetting(ground_k)
    return _SETTINGS[ground_k]


def reproduction():
    """The construction at 300 K against the waters' oxygen study."""
    st300 = setting(S.CHECK_K)
    lo = R.Setting()                                   # lake_oxygen's own setting and its extract
    out = dict(construction=S.check_against_climate_targets(st300))
    out['waters_identical'] = {f'{run}__{k}': bool(np.array_equal(st300.f[f'{run}__{k}'], lo.f[f'{run}__{k}']))
                               for run in ('lock', 'hot') for k in ('waters_ts', 'waters_ps', 'waters_weight',
                                                                     'waters_lat', 'waters_surface_m')}
    out['sun_phase_identical'] = {run: st300.phase[run] == lo.phase[run] for run in ('lock', 'hot', 'eddies')}
    mine = json.loads(json.dumps(R.saturation(st300)))
    theirs = json.loads(R.OUT.read_text())['saturation']
    # the figures: every target's saturation, the waters' distribution and pressures (the checks' own labels may change
    # with lake_oxygen's code while its committed product keeps the labels it was written with)
    figures = ('targets', 'waters_pressure_atm', 'waters_temperature_check_c', 'pressure_atm', 'o2_fraction',
               'salinity_g_kg')
    out['saturation_identical'] = all(mine.get(k) == theirs.get(k) for k in figures)
    out['saturation_compared'] = list(figures)
    out['saturation_check_labels'] = dict(this_code=sorted(mine['checks']), lake_oxygen_product=sorted(theirs['checks']))
    out['resort_300k'] = dict(three_day_means_median_mg_l=mine['targets']['resort']['three_day_means']['median'],
                              seas_median_mg_l=mine['targets']['resort']['three_day_means']['seas_median'],
                              lake_oxygen=dict(three_day_means_median_mg_l=theirs['targets']['resort'][
                                  'three_day_means']['median'], seas_median_mg_l=theirs['targets']['resort'][
                                  'three_day_means']['seas_median']))
    out['proved'] = bool(out['construction']['share_equal'] and out['construction']['ladder_equal']
                         and out['construction']['lapse_equal'] and all(out['waters_identical'].values())
                         and all(out['sun_phase_identical'].values()) and out['saturation_identical'])
    return out


# --- the cases ------------------------------------------------------------------------------------------------------

def cases(st):
    """Every column: name -> case. The central columns, at each basin's typical and coldest cells, start warm and pass
    the record C.SPIN_PASSES times; the wind share's bracket at the typical cells, and the sensitivity rows at the two
    largest basins' typical cells, start from their basin's central typical column and pass C.FROM_PASSES times. The
    other form's columns run as the central ones do, with the form and the fingerprint of its code in the case (and so
    in the cache key); the pans' bracket takes the fixed dissipation of lake_oxygen's physics() and its cap."""
    chosen = T.choose_cells(st)
    out = {}
    deep = lambda key, cell, column, **extra: dict(dict(basin=key, cell=cell, column=column, coupled=CENTRAL_PHYSICS,
                                                        steps_per_mean=DEEP_SPM, passes=C.SPIN_PASSES), **extra)
    variant = lambda key, **extra: deep(key, chosen[key]['typical']['index'], 'typical', passes=C.FROM_PASSES,
                                        start=f'{key}_typical', **extra)
    for key, ch in chosen.items():
        out[f'{key}_typical'] = deep(key, ch['typical']['index'], 'typical')
        if not ch['one_column']:
            out[f'{key}_coldest'] = deep(key, ch['coldest']['index'], 'coldest')
    for key in chosen:
        for eta in R.ETA_BRACKET:
            out[f'{key}_typical_wind_share_{eta}'] = variant(key, eta=eta)
    for key in SENSITIVITY_BASINS:
        out[f'{key}_typical_lake_number'] = variant(key, mixing='lake_number')
        out[f'{key}_typical_fixed_dissipation'] = variant(key, mixing='dissipation')
    out['ingenii_typical_solute_none'] = variant('ingenii', solute='none')
    out['procellarum_typical_oligotrophic_physics'] = variant('procellarum', coupled='oligotrophic/renewal/supply')
    for key in FORM_BASINS:
        out[f'{key}_typical_{FORM}'] = deep(key, chosen[key]['typical']['index'], 'typical', form=FORM,
                                            form_code=FORM_CODE)
    for key, ch in chosen.items():
        for column in ('typical',) if ch['one_column'] else ('typical', 'coldest'):
            for pan in PANS_M if column == 'typical' else PANS_M[:1]:
                out[f'{key}_{column}_pan{pan:g}m'] = dict(basin=key, cell=ch[column]['index'], column=column,
                                                         coupled=R.COUPLED['shallow'], steps_per_mean=PAN_SPM,
                                                         pan_m=pan)
    for key, pan in PAN_BRACKET:
        out[f'{key}_typical_pan{pan:g}m_fixed_dissipation'] = dict(out[f'{key}_typical_pan{pan:g}m'],
                                                                   mixing='dissipation')
    return out, chosen


def cost(st, case):
    """Rough cost for ordering the work: layers times passes."""
    if case.get('pan_m'):
        return case['pan_m'] * 4.0
    return st.basins[case['basin']]['deepest_m'] / 10.0 * (case['passes'] + 4.0 * C.FINE_LIMIT)


def final_state(st, all_cases, name):
    path = C.CASES / f'{name}_{C.key(st, all_cases[name])}.json'
    return json.loads(path.read_text())['final_state']


def work_column(item):
    name, case, all_cases = item
    tic = time.time()
    st = setting()
    start = final_state(st, all_cases, case['start']) if case.get('start') else None
    with deep_mixing_form(case):
        res = C.run(st, name, case, start)
    return name, round(time.time() - tic, 1), res['case']['key']


def work_slab(item):
    i, run_name = item
    tic = time.time()
    res = C.slab(setting(), i, run_name)
    return f'slab {run_name} cell {i}', round(time.time() - tic, 1), res['key']


def run_columns(jobs):
    """The central columns, the pans and the slab checks first; then the variants, from the central columns."""
    st = setting()
    all_cases, chosen = cases(st)
    first = sorted([(n, c, all_cases) for n, c in all_cases.items() if not c.get('start')],
                   key=lambda x: -cost(st, x[1]))
    second = sorted([(n, c, all_cases) for n, c in all_cases.items() if c.get('start')], key=lambda x: -cost(st, x[1]))
    slabs = [(chosen[k]['typical']['index'], run) for k in SLAB_BASINS for run in ('lock', 'hot')]
    tic = time.time()
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for batch in (first, second):
            for name, seconds, k in pool.map(work_column, batch):
                print(f'  {name} ({k}): {seconds:.0f} s, {time.time() - tic:.0f} s in', flush=True)
        for label, seconds, k in pool.map(work_slab, slabs):
            print(f'  {label} ({k}): {seconds:.0f} s', flush=True)
    return write_index(st, all_cases, chosen)


def write_index(st, all_cases, chosen, write=True):
    """The index of every saved column, for the studies that read them (the lake methane study's tracer first);
    written to research/runs/sea_oxygen/columns.json unless write is False."""
    lake_py = Path(L.__file__)
    index = dict(schema=INDEX_SCHEMA, ground_k=st.ground_k, setting=T.describe(st),
                 format=('research/studies/lake_kinds/columns.py\'s arrays, one npz per case for its saved ten-year '
                         'record: t, oxygen, debt, age, solute, pool (days x layers, and x oxygen columns for oxygen, '
                         'debt and pool), z and area (the layers\' interfaces and their areas), mu, dt, pressure, '
                         'every and step, and the step-by-step series ' + ', '.join(C.SERIES + C.COLUMN_SERIES) +
                         ', with ' + ' and '.join(C.EXTRA_SERIES) + ' besides; the oxygen columns are lake_oxygen\'s '
                         'scenarios in the order given; a deep column\'s spin-up is in its case summary'),
                 scenarios=[s.name for s in R.scenarios()],
                 lake_py=dict(path=str(lake_py.resolve().relative_to(ROOT)), sha256=sha256(lake_py)),
                 producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p)
                                      for p in C.code_files() + (Path(__file__), HERE / 'seas.py', HERE / 'extract.py',
                                                                 Path(MF.__file__))},
                               inputs=st.inputs),
                 diffusivity_forms=dict(
                     osborn=('lake.py\'s diffusivity(): K = 0.2 eps / N^2 at every buoyancy Reynolds number, eps the '
                             'wind\'s (deep_mixing wind) or the fixed one (deep_mixing dissipation)'),
                     shih_bouffard=('lake_oxygen/mixing_forms.py\'s shih_bouffard() in place of lake.py\'s '
                                    'diffusivity(), eps the wind\'s: K = 2 nu Re_b^(1/2) above Re_b = 100, 0.2 nu Re_b '
                                    'from 7 to 100 and none below 7, with nu ' + f'{MF.NU:g}' + ' m2/s, the same cap '
                                    'and unstable rule'),
                     lake_number='lake.py\'s lake_number_diffusivity() (Yeates & Imberger 2003)',
                     form_code=FORM_CODE, reading='each case\'s diffusivity says which form its columns ran'),
                 cells=chosen, cases={})
    for name, case in all_cases.items():
        k = C.key(st, case)
        path = C.CASES / f'{name}_{k}.json'
        res = json.loads(path.read_text())
        b = st.basins[case['basin']]
        index['cases'][name] = dict(case=case, key=k, summary=path.name, summary_sha256=sha256(path),
                                    files=res['files'], kind=res['case']['kind'], physics=res['case']['physics'],
                                    diffusivity=diffusivity_form(case), started_from=case.get('start'),
                                    basin=dict(key=b['key'], area_km2=b['area_km2'], deepest_m=b['deepest_m'],
                                               mean_depth_m=b['mean_depth_m'], stands_for=b['stands_for']),
                                    lat=res['case']['lat'], lon=res['case']['lon'])
    planned, _ = cases(st)
    index['cases_planned'] = sorted(planned)
    index['complete'] = set(all_cases) == set(planned)
    if write:
        INDEX.write_text(json.dumps(index, indent=1) + '\n')
    return index


def index_complete(st, index):
    """Whether an index lists every case this code would run, each with its summary present and unchanged."""
    all_cases, _ = cases(st)
    if set(index.get('cases', {})) != set(all_cases):
        return False
    for name, case in all_cases.items():
        entry = index['cases'][name]
        path = C.CASES / entry['summary']
        if entry['key'] != C.key(st, case) or not path.exists() or sha256(path) != entry['summary_sha256']:
            return False
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=2, help='worker processes, at most 2 on the shared machine')
    ap.add_argument('--columns', action='store_true',
                    help='run the columns (hours on two cores; cached, resumable) and write their index')
    ap.add_argument('--cache-only', action='store_true', help='the index and product from the cache; stop at a miss')
    args = ap.parse_args(argv)
    if args.cache_only:
        os.environ[R.CACHE_ONLY] = '1'
    jobs = max(1, min(args.jobs, 2))
    tic = time.time()
    repro = reproduction()
    print('300 K construction', 'proved' if repro['proved'] else 'NOT PROVED', json.dumps(repro['resort_300k']),
          flush=True)
    if not repro['proved']:
        print(json.dumps(repro, indent=1))
        return 1
    st = setting()
    index = None
    if args.columns:
        index = run_columns(jobs)
        print('wrote', INDEX, f'with {len(index["cases"])} cases in {time.time() - tic:.0f} s', flush=True)
    elif INDEX.exists():
        index = json.loads(INDEX.read_text())
        if not index_complete(st, index):
            if R.cache_only():
                raise SystemExit(f'cache only: {INDEX} does not hold every current case; run --columns --cache-only')
            print(f'{INDEX} does not hold every current case; the product takes the setting alone', flush=True)
            index = None
    elif R.cache_only():
        raise SystemExit(f'cache only: {INDEX} is missing; run --columns --cache-only')
    product = analyse(st, index, repro) if index else setting_only(st, repro)
    product['seconds'] = round(time.time() - tic, 1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product) + '\n')
    print('wrote', OUT.relative_to(ROOT), 'with' if index else 'without', 'the columns',
          f'in {time.time() - tic:.0f} s')
    return 0


def producer(st, with_index=False):
    files = (Path(__file__), HERE / 'analysis.py', HERE / 'columns.py', HERE / 'setting.py', HERE / 'seas.py',
             HERE / 'extract.py', Path(S.__file__), Path(KC.__file__), Path(R.__file__), Path(X.__file__),
             Path(MF.__file__), Path(L.__file__))
    inputs = dict(st.inputs, **{str(R.OUT.relative_to(ROOT)): sha256(R.OUT)})
    if with_index:
        inputs[str(INDEX.relative_to(ROOT))] = sha256(INDEX)
    return dict(files={str(Path(f).resolve().relative_to(ROOT)): sha256(f) for f in files}, inputs=inputs,
                constants=dict(MOON_SURFACE_GRAVITY=MOON_SURFACE_GRAVITY, STANDARD_GRAVITY=STANDARD_GRAVITY,
                               SYNODIC_MONTH_DAYS=SYNODIC_MONTH_DAYS))


def common(st, repro, cells):
    """What the product holds whether or not the columns have run."""
    p = np.array([float(np.mean(st.f['lock__ps'][:, ch[c]['index']])) for ch in cells.values()
                  for c in ('typical', 'coldest')]) * 100.0 / 101325.0
    psl_note = ('the GCM\'s surface pressure at the cells, over T21\'s truncated orography; its sea-level pressure there '
                'is 1.199-1.200 atm')
    seas = st.seas
    return dict(
        setting=dict(T.describe(st), salinity_g_kg=st.salinity, o2_fraction=st.o2_fraction, gravity=st.gravity,
                     heat_flow_w_m2=st.heat_flow, design_pressure_pa=R.DESIGN_PRESSURE_PA, cells_rule=T.RULE,
                     central_physics=CENTRAL_PHYSICS, pans_m=list(PANS_M), spin_passes=C.SPIN_PASSES,
                     variant_passes=C.FROM_PASSES, analysis_step_limit=C.FINE_LIMIT,
                     settle_test=dict(deep_t_k=R.SETTLED_K, deepest_tenth_o2_g_m3=R.SETTLED_O2,
                                      bottom_age_share=R.SETTLED_AGE, bottom_solute_g_m3=R.SETTLED_SOLUTE),
                     deep_layer_max_m=C.DEEP_LAYER_MAX_M, wind_share=R.WIND_SHARE, cap_wind_m2_s=R.CAP_WIND,
                     dissipation_w_kg=R.DISSIPATION, fixed_cap_m2_s=R.CAP['central'], trophic=R.TROPHIC,
                     hypoxic_g_m3=R.HYPOXIC, anoxic_g_m3=R.ANOXIC,
                     shih_bouffard=dict(basins=list(FORM_BASINS), nu_m2_s=MF.NU,
                                        re_b_molecular=MF.RE_MOLECULAR, re_b_energetic=MF.RE_ENERGETIC,
                                        form_code=FORM_CODE),
                     pans_fixed_dissipation=[f'{k}_typical_pan{m:g}m' for k, m in PAN_BRACKET],
                     cell_pressure_atm=dict(min=round(float(p.min()), 4), median=round(float(np.median(p)), 4),
                                            max=round(float(p.max()), 4), note=psl_note)),
        reproduction=repro,
        seas=dict({k: seas['seas'][k] for k in seas['seas']}, connections=seas['connections'],
                  other_regions=seas['other_regions'], selection=seas['selection'],
                  run_basins={b['key']: {k: v for k, v in b.items() if k not in ('t21_cells', 'hypsography')}
                              for b in seas['run_basins']}),
        saturation=dict(at_301k=A.saturation_at(st, st.share['resort'])),
        cells=cells)


def setting_only(st, repro):
    """The product before the columns: the setting, the 300 K proof, the seas, their saturation at 301 K and the
    cells the columns take."""
    all_cases, chosen = cases(st)
    out = dict(schema=SCHEMA, status=STATUS_SETTING, producer=producer(st), evidence=EVIDENCE_SETTING,
               reading_rule=READING_RULE_SETTING)
    out.update(common(st, repro, chosen))
    out['columns'] = None
    out['cases_planned'] = {name: case for name, case in all_cases.items()}
    return out


def lakes_pans(st):
    """The lakes' 3 m productive pans at 301 K beside the seas' shallows: lake_oxygen's run_case() for its three band
    lakes in the lake kinds study's Setting301, taken from that study's cache when it holds them under the current
    code (its key includes the lake model's sha256), else run here and cached in this study's runs."""
    st301 = setting_lakes()
    found, files, source = {}, {}, set()
    for i, band in enumerate(R.SITE_BANDS):
        case = dict(kind='shallow', target='resort', site_i=i)
        theirs = ROOT / 'research' / 'runs' / 'lake_kinds' / 'cases' / f'lake_oxygen_shallow_{band}_{KC.key(st301, case)}.json'
        if theirs.exists():
            found[band] = json.loads(theirs.read_text())
            files[str(theirs.relative_to(ROOT))] = sha256(theirs)
            source.add('the lake kinds study\'s cache')
            continue
        mine = C.CASES / f'lakes_pan_{band}_{C.key(st, dict(lakes_pan=band, **case))}.json'
        if not mine.exists():
            if R.cache_only():
                R.uncached(theirs)
            mine.write_text(json.dumps(R.run_case(st301, 'shallow', 'resort', i)) + '\n')
        found[band] = json.loads(mine.read_text())
        files[str(mine.relative_to(ROOT))] = sha256(mine)
        source.add('run here')
    return dict(source=' and '.join(sorted(source)), ground_k=st301.ground_k, files=files, bands=A.lake_pans(found),
                forcing={b: v['case']['forcing'] for b, v in found.items()})


_LAKES = {}


def setting_lakes():
    if 'st' not in _LAKES:
        _LAKES['st'] = S.Setting301(S.TARGET_K)
    return _LAKES['st']


def heat_budget(st, case, res, path):
    """A column's heat budget over its saved decade as lake_oxygen's budget_errors() measures it: the error against
    the column's own heat (at least a kelvin's worth) and as a mean flux over its surface. columns.py saved the error
    against the decade's net heat, or a millionth of the column's heat where that is larger, and near equilibrium the
    net heat is almost nothing. The net heat is rebuilt here: the column's heat at the decade's end less its heat at
    the start, which is the first saved day's column (the column after the record's first step) less the heat that
    step took in; the step is run again from the column at the decade's end, which the next pass would start from as
    the saved decade started from the pass before it (a pan's spin-up passes the record's first year, so its start
    differs from its end by the weather; the step's heat then differs from the start's by about a hundredth of the
    net heat)."""
    rel = float(res['budget']['heat_relative_error'])
    bsn = C.basin(st, case)
    end = KC.state_from_dict(res['final_state'])
    h_end = L.heat_content(bsn, end)
    site = st.site(case['cell'], case['basin'])
    w, forcing, mu = R.forcing_for(st, 'resort', site, 'lock', None, case['steps_per_mean'])
    lake = L.Lake(bsn, C.physics(st, case), R.scenarios(), st.air)
    lake.calibrate_production(forcing)
    one = L.Forcing(**{f.name: (np.asarray(getattr(forcing, f.name))[:1] if np.ndim(getattr(forcing, f.name)) else
                                getattr(forcing, f.name)) for f in dataclasses.fields(L.Forcing)})
    with deep_mixing_form(case):
        first_step = float(lake.run(one, KC.state_from_dict(res['final_state']))['heat_in'])
    with np.load(path) as z:
        if not np.array_equal(np.asarray(z['z']), bsn.z):
            raise ValueError(f'{path.name}: the saved layers are not the case\'s basin')
        first_day = L.State(np.asarray(z['t'][0], float), end.oxygen, end.debt, end.age)
        first_day.ice = float(z['ice'][0])
        seconds = len(z['mu']) * float(z['dt'])
    h_first = L.heat_content(bsn, first_day)
    h0 = h_first - first_step
    net = h_end - h0
    err = rel * max(abs(net), abs(h0) * 1e-6)
    return dict(saved_heat_relative_error=rel, net_heat_share=rnd_sig(net / abs(h0)),
                heat_relative_error=rnd_sig(err / max(abs(h0), L.RHO_CP * float(bsn.volume.sum()))),
                heat_error_w_m2=rnd_sig(err / (float(bsn.area[0]) * seconds)),
                oxygen_relative_error=res['budget']['oxygen_relative_error'])


def rnd_sig(x, digits=3):
    return float(f'{x:.{digits}g}')


def analyse(st, index, repro):
    """The product from the saved columns."""
    cases_ = index['cases']
    load = lambda name: json.loads((C.CASES / cases_[name]['summary']).read_text())
    npz = lambda name: C.CASES / next(iter(cases_[name]['files']))
    basins = st.basins
    sill_of = {k: b['sills_m'] for k, b in basins.items()}
    deep, variants, pans = {}, {}, {}
    for name, entry in cases_.items():
        res = load(name)
        if entry['kind'] == 'shallow':
            pans[name] = A.pan_case(res)
            pans[name]['budget'] = heat_budget(st, entry['case'], res, npz(name))
            continue
        figures = A.deep_case(res, npz(name), sill_of[entry['case']['basin']])
        figures['deep_demand_g_m3_yr'] = A.demand_deepest_tenth(res, npz(name))
        figures['budget'] = heat_budget(st, entry['case'], res, npz(name))
        figures['diffusivity'] = entry.get('diffusivity', diffusivity_form(entry['case']))
        (deep if name.endswith(('_typical', '_coldest')) else variants)[name] = figures
    cells_all = T.cells_of(st)
    shares = {key: A.cold_end_shares(cs, index['cells'][key]['coldest']['min']) for key, cs in cells_all.items()}
    bounds = {}
    for key, ch in index['cells'].items():
        if ch['one_column']:
            continue
        bounds[key] = A.bounds(key, deep[f'{key}_typical'], deep[f'{key}_coldest'], shares[key])
    # each basin's coldest water at its sills (its coldest column where it has one) against its neighbours' deep
    # water in their typical columns
    pairs = []
    for key, b in basins.items():
        src = f'{key}_coldest' if f'{key}_coldest' in deep else f'{key}_typical'
        for other, depth in b['sills_m'].items():
            if f'{other}_typical' in deep and src in deep:
                pairs.append((src, f'{other}_typical', other, depth))
    connections = A.connected(deep, pairs)
    answers = A.answers(deep, variants, index['cells'], {k: b['stands_for'] for k, b in basins.items()})
    # the Moon's rotation against each basin's size, from its typical column's own stratification
    rotation = {key: A.rotation(A.internal_wave_speed(npz(f'{key}_typical'), st.gravity), basins[key])
                for key in index['cells']}
    slabs = [json.loads(p.read_text()) for p in sorted(C.CASES.glob('slab_*.json'))
             if p.name.rsplit('_', 1)[1].split('.')[0] == C.key(st, dict(slab=p.name.split('_')[1],
                                                                          cell=int(p.name.split('_')[2])))]
    out = dict(schema=SCHEMA, status=STATUS_COLUMNS, producer=producer(st, with_index=True), evidence=EVIDENCE,
               reading_rule=READING_RULE)
    out.update(common(st, repro, index['cells']))
    every = [*deep.values(), *variants.values(), *pans.values()]
    out.update(
        columns=dict(deep=len(deep), variants=len(variants), pans=len(pans),
                     diffusivity={form: sorted(n for n, e in cases_.items() if e.get('diffusivity') == form)
                                  for form in ('shih_bouffard', 'lake_number')}),
        answers=answers, seas_wide=A.seas_wide(answers), forms=A.forms(answers), deep=deep, variants=variants,
        pans=pans, pans_bracket=A.pan_bracket(pans), bounds=bounds, connections=connections, rotation=rotation,
        lakes_pans=lakes_pans(st),
        checks=dict(slab=slabs,
                    max_heat_relative_error=max(abs(c['budget']['heat_relative_error']) for c in every),
                    max_heat_error_w_m2=max(abs(c['budget']['heat_error_w_m2']) for c in every),
                    max_saved_heat_relative_error=max(abs(c['budget']['saved_heat_relative_error']) for c in every),
                    max_oxygen_relative_error=max(c['budget']['oxygen_relative_error'] for c in every),
                    ice_max_m=max(c['ice_max_m'] for c in every),
                    settled_at_the_analysis_step={k: v['spin_up']['at_the_analysis_step']['settled']
                                                  for k, v in {**deep, **variants}.items()}),
        index=dict(path=str(INDEX.relative_to(ROOT)), sha256=sha256(INDEX), cases=len(cases_)))
    return out


EVIDENCE = (
    'A one-dimensional lake model (biosphere/ecology/lake.py as corrected after the independent review of 18 October '
    '2026; checked by its own tests against analytic diffusion, convective and wind deepening, reaeration and '
    'conduction-limited ice) run with lake_oxygen\'s physics and scenarios on each basin\'s own hypsography from LOLA '
    'at 16 pixels a degree with the held seas\' 0.14 g/kg, forced by the GCM\'s own 3-day means at the basin\'s T21 '
    'water cells: the lock\'s weather moved to 301 K by the lake kinds study\'s Setting301, which interpolates between '
    'the design case and the hot world\'s test by the ground\'s temperature and so mixes two physics. The '
    'construction reproduces the waters\' oxygen study\'s 300 K figures exactly. Below the mixed layer the wind\'s '
    'share of its energy flux sets the dissipation (Earth lakes\' 0.3%, bracketed 0.2-0.4%), which Osborn\'s rule '
    'turns into a diffusivity in every column but two; at Procellarum\'s and Ingenii\'s typical cells Shih et '
    'al.\'s (2005) regimes of the buoyancy Reynolds number, as Bouffard & Boegman (2013) carry them to lakes, take its '
    'place (lake_oxygen/mixing_forms.py\'s swap inside the process), mixing the energetic deep water less. Oxygen '
    'demand is lake_oxygen\'s bounded scenarios from Earth\'s lakes. Every deep column passes the ten-year record 80 '
    'times from a warm start (variants 30 times from their basin\'s central column) in coarse steps, then passes it at '
    'the analysis step until lake_oxygen\'s settle test passes or ten passes, before the decade that is saved; eight '
    'of the 52 deep columns do not settle there, and analysis_step names the way each one\'s decade leans. A column '
    'has no horizontal circulation; these basins\' Lake numbers lie far under 1 and the large ones span tens of '
    'internal Rossby radii, so that the wind drives Ekman transport, upwelling along the shores and Kelvin and '
    'Poincare waves, which a column does not hold: each basin is run at its typical and its coldest cell as bounds, '
    'and the effect of exchange is a screen. No lunar sea has been measured.')
READING_RULE = (
    'Oxygen in g/m3 (mg/L), temperatures in C, depths in m. Quantile lists are the 10th, 50th and 90th percentiles. '
    'answers gives each basin\'s deep water: the typical column at the central wind share, the variants at its typical '
    'cell keyed typical_<variant> (the share\'s bracket wind_share_low and _high; at the two largest basins '
    'lake_number, fixed_dissipation, solute_none or oligotrophic_physics, and shih_bouffard, the other form of the deep '
    'mixing), and the coldest column; their oxygen is the saved decade\'s, after the passes at the analysis step, each '
    'trophic state\'s central column (surface-renewal transfer, demand fixed by supply, a 4-day labile pool) with its '
    'bracket across the other transfer and the pool\'s 2 and 7 days. analysis_step says whether those passes settled '
    'and, where they did not, which way the decade leans (reading; oxygen_is and age_is, upper or lower values). '
    'decade_overturns counts the saved decade\'s lunar cycles in which the mixed layer reached the bottom layer, '
    'decade_overturns_a_century the same as a rate; overturns_a_century and the long_ figures are the coarse passes\' '
    'last half, in 4-step means, which leave more oxygen in the deep water than the analysis step. seas_wide weighs '
    'each basin by the share of the seas it stands for, one column of each basin to a view, naming the unsettled '
    'columns it takes; forms sets the two forms of the deep mixing side by side where both ran, with their range. deep '
    'and variants hold each column\'s full figures, keyed basin_column[_variant]: spin_up gives the coarse passes\' '
    'long state and the passes at the analysis step; profiles gives figures from the daily columns: the deepest '
    'tenth\'s age since the water last touched the surface, the time-mean shares of the volume and the floor under '
    'hypoxic (2) and anoxic (0.5) water, the nights\' mixing and the temperature at the sills to the other basins of '
    'the sea; lake_number the Lake number by the day\'s warm layer and the month\'s thermocline; diffusivity the form '
    'of the deep mixing; budget the heat and oxygen errors in lake_oxygen\'s measures, the heat rebuilt from the saved '
    'decade beside the saved figure. pans are the shelves and bays, a 3 m and a 10 m pan with eutrophic physics, their '
    'nights from sunset to sunrise, beside lakes_pans, the lakes\' 3 m pans at 301 K; pans_bracket sets a few of them '
    'beside their runs under the fixed dissipation. bounds compare a basin\'s typical and coldest columns, with the '
    'screen of its deep water under exchange, read from the saved decade\'s overturns (screen) and the coarse passes\' '
    '(screen.coarse), the deepest tenth\'s reference demand and the factor the warm kinetics would raise it by. '
    'connections say whether a basin\'s water at its sill is ever colder than its neighbour\'s deep water. rotation '
    'gives each basin\'s internal Rossby radius from its typical column\'s stratification (wave in m/s) against its '
    'size. saturation.at_301k follows lake_oxygen\'s saturation() at the 301 K share. The saved columns are listed in '
    'research/runs/sea_oxygen/columns.json with their sha256s.')
STATUS_COLUMNS = 'the columns run on the lake model as its sha256 in producer.files gives it'
STATUS_SETTING = ('the sea columns wait for the lake model\'s revision after the independent review of 18 October '
                  '2026; this product holds what does not depend on the model\'s numbers')
EVIDENCE_SETTING = (
    'Hydrostatic geometry of the seas (results/seas.json, from LOLA at 16 pixels a degree above the GRAIL geoid), the '
    'GCM\'s own 3-day means at the seas\' T21 water cells moved to 301 K by the lake kinds study\'s Setting301, which '
    'interpolates between the design case and the hot world\'s test by the ground\'s temperature and so mixes two '
    'physics, and oxygen at saturation by the Benson and Krause fit (Garcia & Gordon 1992) with Benson & Krause\'s '
    'pressure factor, extended from its stated 1.1 atm to the Moon\'s 1.2. The construction reproduces the waters\' '
    'oxygen study\'s 300 K saturation figures exactly. No sea column has been run on the revised lake model.')
READING_RULE_SETTING = (
    'Oxygen in g/m3 (mg/L), temperatures in C. saturation.at_301k gives the seas at their median and band '
    'temperatures at 1.2 atm and the distribution of every water cell\'s 3-day means at its own pressure, at the 301 K '
    'share as lake_oxygen\'s saturation() reckons its targets. cells gives each run basin\'s typical and coldest T21 '
    'cell with their 301 K surface temperatures (the GCM\'s own, moved as the forcing is moved). cases_planned lists '
    'the columns run.py --columns would run; columns is null until they have run.')


if __name__ == '__main__':
    sys.exit(main())
