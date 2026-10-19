"""How the deep lakes' oxygen and methane rest on the form of the mixing below the mixed layer: a sensitivity study.

    python -m research.studies.lake_oxygen.mixing_forms                 # about 90 minutes on one core; resumes
    python -m research.studies.lake_oxygen.mixing_forms --jobs 2        # two cores
    python -m research.studies.lake_oxygen.mixing_forms --cache-only    # the product from the cached columns

The lake model (biosphere/ecology/lake.py) turns the dissipation below the mixed layer into a diffusivity by Osborn's
rule, K = 0.2 eps / N^2, at every buoyancy Reynolds number Re_b = eps / (nu N^2). In the near-neutral deep water the
wind's share of its energy keeps, Re_b runs to 10^3-10^5, where Shih et al. (2005) find the mixing efficiency falling:
K = 2 nu Re_b^(1/2) above Re_b = 100, Osborn's 0.2 nu Re_b between 7 and 100 and molecular diffusion alone below 7,
the three regimes Bouffard & Boegman (2013) carry to lakes. This script runs the lake model in its own process with
lake.py's diffusivity() replaced by those three regimes (lake.py's file unchanged; nu 1e-6 m2/s, the cap and the
unstable rule as lake.py's), on:
- the 352 m Mandel'shtam lake at 0-15 in the lock, the resort at 300 K and the hot world, as lake_oxygen runs it
  (mesotrophic physics carrying every trophic state's columns), beside lake_oxygen's product;
- the same lake at the resort's 301 K, each trophic state with its own physics, as the lake kinds' ladder runs it,
  and the 1,925 m giant at 301 K started from it carried down and oxygenated, as the ladder starts its giants;
- the bracket of the 3 m pans, whose share of the wind's energy extends a measurement made below deep lakes'
  thermoclines: lake_oxygen's shallows with the share at 0.15%, half its central value, and under the first round's
  fixed dissipation (1e-10 W/kg, cap 1e-4 m2/s), with Osborn's rule;
- one variant of Osborn's rule with the dissipation held at most at the wind's share spread over the water below 20 m,
  deeper than the ordinary nights mix the 352 m lakes, to show what its inverse scaling with the water left below the
  mixed layer does as an overturn nears.
Each deep column is saved daily in the lake kinds' format (lake_kinds/columns.py's capture: temperature, the oxygen
columns, debt, age, solutes and the labile pool, with the step-by-step series) to research/runs/lake_oxygen/
mixing_forms, indexed with each file's sha256 in mixing_forms_index.json, so the methane tracer can read them;
lake_anoxia's methane arithmetic runs on the 0-15 lakes of both forms, the lake standing for every lake.
results/mixing_forms.json (schema terluna.research.lake-oxygen-mixing-forms/1) holds the summaries beside the Osborn
runs they answer. Each column and pan is cached under a key of the code, what the runs read from lake_oxygen's setting
(its values_fingerprint()) and the case; --cache-only stops at the first one it would have to compute.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import run as A
from research.studies.lake_kinds import columns as C
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'lake_oxygen' / 'mixing_forms'
CASES = RUNS / 'cases'
INDEX = RUNS / 'mixing_forms_index.json'
OUT = HERE / 'results' / 'mixing_forms.json'
SCHEMA = 'terluna.research.lake-oxygen-mixing-forms/1'
NU = 1.0e-6                          # m2/s, the kinematic viscosity the energetic regime takes (0.77-1.0e-6 at 20-32 C)
RE_MOLECULAR, RE_ENERGETIC = 7.0, 100.0        # Shih et al. (2005): the diffusive, intermediate and energetic regimes
PANS_SHARE = 0.0015                  # the pans' low share of the wind's energy, half the central 0.3%
BOUND_M = 20.0                       # the bounded variant: the dissipation at most the share over the water below this
GIANT_M, GIANT_PASSES = 1925.0, 5    # the ladder's giant at 301 K and its fifty saved years
BAND = 0                             # the 0-15 lake, Mandel'shtam
TARGETS = ('lock', 'resort', 'hot')
FORMS = ('osborn', 'shih_bouffard')
_OSBORN = L.diffusivity
_WIND = L.wind_dissipation


def shih_bouffard(n2, physics, dissipation=None):
    """The diffusivity below the mixed layer in Shih et al.'s (2005) three regimes of the buoyancy Reynolds number
    Re_b = eps / (nu N^2), as Bouffard & Boegman (2013) take them for lakes: 2 nu Re_b^(1/2) above 100, Osborn's
    efficiency times nu Re_b from 7 to 100, and none beyond the molecular below 7; capped, with lake.py's rule across
    statically unstable interfaces."""
    eps = physics.dissipation if dissipation is None else float(dissipation)
    n2 = np.asarray(n2, dtype=float)
    reb = eps / (NU * np.maximum(n2, 1e-30))
    k = np.where(reb > RE_ENERGETIC, 2.0 * NU * np.sqrt(reb),
                 np.where(reb > RE_MOLECULAR, physics.mixing_efficiency * NU * reb, 0.0))
    unstable = max(physics.unstable_diffusivity, physics.max_diffusivity)
    return np.where(n2 > 0.0, np.minimum(k, physics.max_diffusivity), unstable)


@contextlib.contextmanager
def mixing_form(form, bound_kg_m2=None):
    """lake.py's diffusivity() and, for the bounded variant, its wind_dissipation() replaced inside this process."""
    L.diffusivity = shih_bouffard if form == 'shih_bouffard' else _OSBORN
    if bound_kg_m2:
        def bounded(p10, hypolimnion_kg_m2, physics):
            return physics.wind_share * p10 / max(hypolimnion_kg_m2, bound_kg_m2) if hypolimnion_kg_m2 > 0.0 else 0.0
        L.wind_dissipation = bounded
    try:
        yield
    finally:
        L.diffusivity, L.wind_dissipation = _OSBORN, _WIND


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def rel(p) -> str:
    """A path as the products give it: relative to the repository where it lies under it."""
    try:
        return str(Path(p).relative_to(ROOT))
    except ValueError:
        return str(p)


def code_files():
    return (Path(__file__), Path(L.__file__), Path(R.__file__), Path(C.__file__), Path(S.__file__))


_SETTINGS = {}


def setting(ground_k):
    """The resort at ground_k (lake_kinds' construction, which at 300 K is lake_oxygen's own); the lock and the hot
    world are lake_oxygen's whatever ground_k."""
    if ground_k not in _SETTINGS:
        _SETTINGS[ground_k] = S.Setting301(ground_k)
    return _SETTINGS[ground_k]


def key(st, case):
    """The cache key: the code that runs a column whole, what the runs read from the setting's inputs and the case."""
    h = hashlib.sha256()
    for p in code_files():
        h.update(sha256(p).encode())
    h.update(st.values_fingerprint().encode())
    h.update(json.dumps(dict(case, ground_k=st.ground_k), sort_keys=True).encode())
    return h.hexdigest()[:20]


def column(name, case, start=None):
    """One deep column, saved pass by pass as the lake kinds' columns.py saves them, and summarised as lake_oxygen
    summarises its runs. case: target, ground_k, depth_m, coupled, form, bound_m, spin ('deep' or 'long'), passes."""
    st = setting(case['ground_k'])
    CASES.mkdir(parents=True, exist_ok=True)
    k = key(st, case)
    final = CASES / f'{name}_{k}.json'
    if final.exists():
        return json.loads(final.read_text())
    if R.cache_only():
        R.uncached(final)
    tic = time.time()
    site = R.site_info(st, BAND)
    w, forcing, mu = R.forcing_for(st, case['target'], site, 'lock', None, R.STEPS_PER_MEAN_DEEP)
    bsn = R.basin('deep', site, case['depth_m'])
    scen = R.scenarios()
    names = [s.name for s in scen]
    ph = R.physics(st, 'deep', coupled=case['coupled'])
    lake = L.Lake(bsn, ph, scen, st.air)
    lake.calibrate_production(forcing)
    col = ph.solute_column
    bound = None
    if case.get('bound_m'):
        cum = np.cumsum(bsn.volume)
        above = float(np.interp(case['bound_m'], bsn.z[1:], cum))
        bound = L.RHO0 * (float(cum[-1]) - above) / float(bsn.area[0])
    every = max(1, int(round(C.PROFILE_S / forcing.dt)))
    with mixing_form(case['form'], bound):
        done = []
        while (CASES / f'{name}_{k}_pass{len(done):02d}.json').exists():
            done.append(CASES / f'{name}_{k}_pass{len(done):02d}.npz')
        if done:
            meta = json.loads(done[-1].with_suffix('.json').read_text())
            state, spun = C.state_from_dict(meta['state']), meta['spin_up']
        else:
            if start is not None:
                state = R.state_from(start, bsn, scen)
                state.oxygen[:] = C.saturated(st, float(state.t.mean()))       # oxygenated, as if just overturned
                state.debt[:], state.age[:], state.solute[:], state.pool[:] = 0.0, 0.0, 0.0, 0.0
                start_c = None
            else:
                start_c = float(forcing.air_k.mean()) - 273.15 + R.WARM_START_K
                state = L.initial_state(bsn, scen, start_c, C.saturated(st, start_c))
            which = 'deepest_long' if case['spin'] == 'long' else 'deep'
            spun = R.settle_passes(lake, R.coarse(forcing), state, col, bsn, R.MAX_PASSES[which])
            spun['at_the_analysis_step'] = R.settle_passes(lake, forcing, state, col, bsn, R.FINE_PASSES[which])
            spun['start_c'] = None if start_c is None else round(start_c, 2)
        for p in range(len(done), case['passes']):
            h0, o0 = L.heat_content(bsn, state), L.oxygen_content(bsn, state)
            rec, out = C.capture(lake, forcing, state, every)
            budget = dict(R.budget_errors(bsn, state, h0, o0, out, len(forcing) * forcing.dt),
                          produced=out['produced'].tolist(), respired=out['respired'].tolist(),
                          oxygen_in=out['oxygen_in'].tolist())
            series = {s: np.asarray(out[s], dtype=float) for s in C.SERIES + C.COLUMN_SERIES}
            path = CASES / f'{name}_{k}_pass{p:02d}.npz'
            np.savez_compressed(path, z=bsn.z, area=bsn.area, mu=mu, dt=forcing.dt,
                                pressure=forcing.pressure.astype(np.float32), every=every, **rec, **series)
            path.with_suffix('.json').write_text(json.dumps(dict(
                pass_index=p, budget=budget, spin_up=spun, interfaces={i: int(out[i]) for i in C.INTERFACES},
                state=C.state_dict(bsn, state), npz_sha256=sha256(path))) + '\n')
            del rec, out, series
    passes = [np.load(CASES / f'{name}_{k}_pass{p:02d}.npz') for p in range(case['passes'])]
    metas = [json.loads((CASES / f'{name}_{k}_pass{p:02d}.json').read_text()) for p in range(case['passes'])]
    out = {s: np.concatenate([z[s] for z in passes]).astype(float) for s in C.SERIES + C.COLUMN_SERIES}
    mu_all = np.concatenate([z['mu'] for z in passes]).astype(float)
    res = R.summarise('deep', out, mu_all, bsn, scen, forcing.dt)
    res['deep_mixing'] = dict(C.mixing_record(ph, metas, out['dissipation']), form=case['form'],
                              bound_kg_m2=bound)
    res['budget'] = {b: max(abs(m['budget'][b]) for m in metas)
                     for b in ('heat_relative_error', 'heat_error_w_m2', 'oxygen_relative_error',
                               'oxygen_error_g_m2_day')}
    res['settling'] = dict(spin_up=metas[0]['spin_up'])
    res['bottom_solute_g_m3'] = dict(mean=round(float(out['bottom_solute'].mean()), 3),
                                     max=round(float(out['bottom_solute'].max()), 3))
    days = len(mu_all) * forcing.dt / 86400.0
    res['case'] = dict(case, name=name, key=k, band=site['band'], depth_m=bsn.depth, layers=bsn.n,
                       area_km2=round(bsn.area[0] / 1e6, 1), step_s=forcing.dt, profile_every_steps=every,
                       years=round(days / 365.25, 2), scenarios=names, seconds=round(time.time() - tic, 1))
    res['files'] = {f'{name}_{k}_pass{p:02d}.npz': metas[p]['npz_sha256'] for p in range(case['passes'])}
    res['final_state'] = metas[-1]['state']
    for z in passes:
        z.close()
    final.write_text(json.dumps(res) + '\n')
    return res


def pan(name, case):
    """lake_oxygen's 3 m pan through its own run_case(), with the wind's share at the pans' low value or under the
    first round's fixed dissipation, cached here."""
    st = setting(300.0)
    CASES.mkdir(parents=True, exist_ok=True)
    path = CASES / f'{name}_{key(st, case)}.json'
    if path.exists():
        return json.loads(path.read_text())
    if R.cache_only():
        R.uncached(path)
    R.WIND_SHARE['pans_low'] = PANS_SHARE
    if case['bracket'] == 'share_low':
        res = R.run_case(st, 'shallow', case['target'], case['site_i'], eta='pans_low')
    else:
        res = R.run_case(st, 'shallow', case['target'], case['site_i'], mixing='dissipation', dissipation='central')
    res.pop('final_state', None)
    path.write_text(json.dumps(res) + '\n')
    return res


def cases():
    """name -> (case, the case it starts from)."""
    out = {}
    for trophic in R.TROPHIC:                                 # first, so the giants can start from them
        out[f'shih_301k_352m_0-15_{trophic}'] = (dict(target='resort', ground_k=301.0, depth_m=R.DEEP_M,
                                                      coupled=f'{trophic}/renewal/supply', form='shih_bouffard',
                                                      spin='deep', passes=1), None)
    for target in TARGETS:
        out[f'shih_{target}_352m_0-15'] = (dict(target=target, ground_k=300.0, depth_m=R.DEEP_M,
                                                coupled=R.COUPLED['deep'], form='shih_bouffard', spin='deep',
                                                passes=1), None)
    out['osborn_bounded_lock_352m_0-15'] = (dict(target='lock', ground_k=300.0, depth_m=R.DEEP_M,
                                                 coupled=R.COUPLED['deep'], form='osborn', bound_m=BOUND_M,
                                                 spin='deep', passes=1), None)
    for trophic in ('oligotrophic', 'mesotrophic', 'eutrophic', 'new_reservoir'):
        start = f'shih_301k_352m_0-15_{trophic}'
        out[f'shih_301k_giant_{GIANT_M:g}m_0-15_{trophic}'] = (dict(
            target='resort', ground_k=301.0, depth_m=GIANT_M, coupled=f'{trophic}/renewal/supply',
            form='shih_bouffard', spin='long', passes=GIANT_PASSES), start)
    return out


def pan_cases():
    return {f'pan_{bracket}_{target}_{R.SITE_BANDS[i]}': dict(bracket=bracket, target=target, site_i=i)
            for bracket in ('share_low', 'fixed_dissipation') for target in TARGETS for i in range(3)}


def _work(item):
    kind, name, case, start = item
    tic = time.time()
    res = pan(name, case) if kind == 'pan' else column(name, case, start)
    return kind, name, res, round(time.time() - tic, 1)


# --- the methane on both forms ------------------------------------------------------------------------------------

def methane(index):
    """lake_anoxia's methane arithmetic on the 0-15 lake of each target under both forms: Osborn's from lake_anoxia's
    own saved columns, Shih and Bouffard's from this study's, the lake standing for every lake Moon-wide."""
    class Columns(A.Setting):
        def load(self, name):
            path = CASES / name if name.endswith('.npz') else A.RUNS / f'{name}.npz'
            return dict(np.load(path))

    gas = Columns()
    classes = A.lake_classes(A.moon_lakes())
    out = {}
    for form in FORMS:
        deep = {}
        for target in TARGETS:
            if form == 'osborn':
                name = f'deep_{target}_0-15'
                case = dict(target=target, band='0-15', every=gas.meta['cases'][name]['every'])
            else:
                entry = index['cases'][f'shih_{target}_352m_0-15']
                name = list(entry['files'])[0]
                case = dict(target=target, band='0-15', every=int(np.load(CASES / name)['every']))
            deep[f'deep_{target}_0-15'] = A.deep_lake(gas, name, case)
        temps = {t: float(np.mean([v['bottom_t_c'] for v in deep.values() if v['target'] == t])) for t in TARGETS}
        b = A.budget(gas, deep, None, classes, temps)
        out[form] = dict(
            deep_lakes={k: dict(overturn_events=v['overturn_events'],
                                traits={t: dict(methane_made_mol_m2_yr=x['methane_made_mol_m2_yr'],
                                                released_at_overturns_mol_m2_yr={
                                                    r: sum(x['released_at_overturns_mol_m2'][r]) / v['years']
                                                    for r in A.RELEASE_RULES},
                                                anoxic_floor_days_per_year=x['anoxic_floor']['days_per_year'])
                                        for t, x in v['traits'].items()}) for k, v in deep.items()},
            moon_tg_yr={t: {tr: {r: b[t]['traits'][tr]['by_release_rule'][r]['moon_tg_yr'] for r in A.RELEASE_RULES}
                            for tr in b[t]['traits']} for t in b})
    return out


def compare(index):
    """Each form's figures beside the Osborn runs they answer: lake_oxygen's product for the 0-15 lake in each target,
    and the lake kinds' ladder for the 352 m lake and the giant at 301 K where its cases are saved."""
    product = json.loads(R.OUT.read_text())
    rows = {}
    for target in TARGETS:
        mine = index['summaries'][f'shih_{target}_352m_0-15']
        theirs = product['runs']['deep'][f'{target}/0-15']
        rows[target] = {form: dict(complete_overturns=s['mixing']['complete_overturns'], cycles=s['cycles'],
                                   bottom_age_years=s['mixing']['bottom_age_days']['mean'] / 365.25,
                                   deepest_tenth={t: dict(mean=s['oxygen'][f'{t}/renewal/supply']['deepest_tenth'][
                                       'mean'], anoxic_share=s['oxygen'][f'{t}/renewal/supply']['deepest_tenth'][
                                       'anoxic_share']) for t in R.TROPHIC})
                        for form, s in (('osborn', theirs), ('shih_bouffard', mine))}
    b = index['summaries']['osborn_bounded_lock_352m_0-15']
    rows['lock']['osborn_bounded'] = dict(complete_overturns=b['mixing']['complete_overturns'], cycles=b['cycles'],
                                          bottom_age_years=b['mixing']['bottom_age_days']['mean'] / 365.25,
                                          deepest_tenth={t: dict(mean=b['oxygen'][f'{t}/renewal/supply'][
                                              'deepest_tenth']['mean'], anoxic_share=b['oxygen'][
                                              f'{t}/renewal/supply']['deepest_tenth']['anoxic_share'])
                                              for t in R.TROPHIC})
    at301 = {}
    from research.studies.lake_kinds import ladder as LD
    st = setting(301.0)
    for trophic in R.TROPHIC:
        for depth, mine_name, ladder_name in ((R.DEEP_M, f'shih_301k_352m_0-15_{trophic}',
                                               LD.ladder_name(0, R.DEEP_M, trophic)),
                                              (GIANT_M, f'shih_301k_giant_{GIANT_M:g}m_0-15_{trophic}',
                                               LD.giant_name(0, GIANT_M, trophic))):
            s = index['summaries'][mine_name]
            entry = {'shih_bouffard': s}
            lcase = LD.cases().get(ladder_name)
            if lcase:
                path = C.CASES / f'{ladder_name}_{C.key(st, lcase[0])}.json'
                if path.exists():
                    entry['osborn'] = json.loads(path.read_text())
                elif R.cache_only():
                    R.uncached(path)
            at301[f'{depth:g}/{trophic}'] = {form: dict(
                complete_overturns=v['mixing']['complete_overturns'], cycles=v['cycles'],
                bottom_age_years=v['mixing']['bottom_age_days']['mean'] / 365.25,
                deepest_tenth=dict(mean=v['oxygen'][f'{trophic}/renewal/supply']['deepest_tenth']['mean'],
                                   anoxic_share=v['oxygen'][f'{trophic}/renewal/supply']['deepest_tenth'][
                                       'anoxic_share'])) for form, v in entry.items()}
    pans = {}
    for name, s in index['pans'].items():
        c = pan_cases()[name]
        theirs = product['runs']['shallow'][f"{c['target']}/{R.SITE_BANDS[c['site_i']]}"]
        pans[name] = {label: {t: dict(night_minimum=v['oxygen'][f'{t}/renewal/supply']['nights']['surface_minimum'],
                                      hypoxic_nights=v['oxygen'][f'{t}/renewal/supply']['nights']['hypoxic_share'],
                                      bottom_anoxic_share=v['oxygen'][f'{t}/renewal/supply']['bottom']['anoxic_share'],
                                      day_hours_bottom_anoxic=v['oxygen'][f'{t}/renewal/supply']['days'][
                                          'hours_bottom_anoxic'])
                              for t in R.TROPHIC}
                      for label, v in ((c['bracket'], s), ('central', theirs))}
    return dict(lake_oxygen_convention=rows, at_301k=at301, pans=pans)


def strip(x):
    return {k: v for k, v in x.items() if k not in ('final_state',)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=1, help='worker processes, at most 2 on the shared machine')
    ap.add_argument('--cache-only', action='store_true', help='the product from the cached columns; stop at a miss')
    args = ap.parse_args(argv)
    if args.cache_only:
        os.environ[R.CACHE_ONLY] = '1'
    jobs = max(1, min(args.jobs, 2))
    tic = time.time()
    all_cases = cases()
    first = [('column', n, c, None) for n, (c, s) in all_cases.items() if s is None]
    first += [('pan', n, c, None) for n, c in pan_cases().items()]
    summaries, pans = {}, {}
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for kind, name, res, seconds in pool.map(_work, first):
            (pans if kind == 'pan' else summaries)[name] = res
            print(f'  {name}: {seconds:.0f} s, {time.time() - tic:.0f} s in', flush=True)
        second = [('column', n, c, summaries[s]['final_state']) for n, (c, s) in all_cases.items() if s is not None]
        for kind, name, res, seconds in pool.map(_work, second):
            summaries[name] = res
            print(f'  {name}: {seconds:.0f} s, {time.time() - tic:.0f} s in', flush=True)
    st300, st301 = setting(300.0), setting(301.0)
    index = dict(schema=SCHEMA.replace('/1', '-index/1'), lake_py_sha256=sha256(L.__file__),
                 producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p) for p in code_files()},
                               inputs=st300.inputs),
                 cases={n: dict(case=all_cases[n][0], start=all_cases[n][1], key=s['case']['key'],
                                summary=f"{n}_{s['case']['key']}.json",
                                summary_sha256=sha256(CASES / f"{n}_{s['case']['key']}.json"), files=s['files'])
                        for n, s in summaries.items()},
                 summaries=summaries, pans=pans)
    INDEX.write_text(json.dumps({k: v for k, v in index.items() if k not in ('summaries', 'pans')}, indent=1) + '\n')
    gases = methane(index)
    cmp = compare(index)
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.resolve().relative_to(ROOT)): sha256(p) for p in code_files()},
                      inputs=dict(st300.inputs, **{rel(R.OUT): sha256(R.OUT), rel(INDEX): sha256(INDEX)})),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        setting=dict(nu_m2_s=NU, re_b_molecular=RE_MOLECULAR, re_b_energetic=RE_ENERGETIC,
                     wind_share=R.WIND_SHARE['central'], pans_share=PANS_SHARE, bound_m=BOUND_M,
                     cap_m2_s=R.CAP_WIND, fixed_dissipation_w_kg=R.DISSIPATION['central'],
                     fixed_cap_m2_s=R.CAP['central'],
                     share_of_the_way_301k=st301.share['resort'], lapse_k_per_km_301k=st301.lapse['resort'] * 1000.0),
        comparison=cmp, methane=gases,
        runs={n: strip(s) for n, s in summaries.items()}, pans=pans,
        index=rel(INDEX), seconds=round(time.time() - tic, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(A.dumps(A.rnd(out)) + '\n')
    print('wrote', rel(OUT), 'and', rel(INDEX), f'in {time.time() - tic:.0f} s')
    return 0


EVIDENCE = (
    'A sensitivity study on one lake model (biosphere/ecology/lake.py) through one GCM\'s climate at T21: the same '
    'runs as lake_oxygen\'s and the lake kinds\' ladder, with the diffusivity below the mixed layer from Shih et '
    'al.\'s (2005) three regimes of the buoyancy Reynolds number in place of Osborn\'s rule at every one, as Bouffard'
    ' & Boegman (2013) carry them to lakes; the pans\' share of the wind\'s energy at half its central value and '
    'under the first round\'s fixed dissipation; and Osborn\'s rule with the dissipation bounded near the floor. '
    'lake.py\'s file is unchanged: its diffusivity() and wind_dissipation() are replaced inside this process. The '
    'methane is lake_anoxia\'s arithmetic on the 0-15 lakes, standing for every lake. No lunar lake has been '
    'measured.')
READING_RULE = (
    'Oxygen in g/m3 (mg/L). comparison.lake_oxygen_convention gives the 0-15 lake in each target under Osborn\'s '
    'rule (lake_oxygen\'s product) and Shih and Bouffard\'s form (shih_<target>_352m_0-15, mesotrophic physics '
    'carrying every trophic state), with the lock\'s bounded variant; comparison.at_301k the 352 m lake and the '
    '1,925 m giant at 301 K, each trophic state with its own physics, under Shih and Bouffard\'s form beside the '
    'ladder\'s Osborn runs; comparison.pans the 3 m pans\' bracket beside lake_oxygen\'s central shallows. methane '
    'holds lake_anoxia\'s methane made in anoxic water, released at the overturns by its two rules and Moon-wide in '
    'Tg a year, from the 0-15 lake of each target under each form. runs and pans hold the summaries; the saved '
    'columns are listed with their sha256 in the index '
    '(research/runs/lake_oxygen/mixing_forms/mixing_forms_index.json), each pass a file in the lake kinds\' format.')


if __name__ == '__main__':
    sys.exit(main())
