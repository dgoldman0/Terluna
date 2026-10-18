"""The depth ladder at 301 K: the lake model on lakes whose floors lie at a ladder of depths, and the giants.

    python -m research.studies.lake_kinds.ladder              # about 2.5 hours on two cores; resumes where it stopped
    python -m research.studies.lake_kinds.ladder --jobs 1     # one core

First it proves the 301 K construction: lake_oxygen's own run_case() at the resort built here at 300 K must give the
waters' oxygen study's committed resort runs (its three 352 m lakes and the 0-15 pan), and this study's column runner
must give the same 352 m lake. Then, at 301 K, with the wind's deep mixing (lake_oxygen's central treatment, its share
of the wind's energy at 0.3%, with the bracket 0.2-0.4% where named):
- the ladder: each band lake's hypsography (Mandel'shtam, Gagarin, Rabbi Levi) scaled to floors at 5, 8, 12, 16, 20,
  30, 50, 100 and 352 m, ten years each after the spin-up in coarse steps, for the lakes of 100 m and more then at the
  analysis step until they settle there too (the shallower rungs overturn within months, so the coarse steps' deep
  water is renewed within the analysed decade), once for each trophic state with its own physics (its water's clarity
  and its respiration's products), and the 352 m lakes' oligotrophic and mesotrophic water at the wind share's
  bracket, the column saved daily;
- the pans: lake_oxygen's 3 m productive shallows, run by its own run_case();
- the giants: floors at 1,000, 1,925 and 3,500 m, started from the 352 m column of the same band, trophic state and
  wind share carried down and oxygenated, as lake_oxygen starts its 1,925 m lake, then passed over the record in
  coarse steps until they settle or for 300 years (lake_oxygen's long limit, thirty passes), and at the analysis
  step until they settle there or for 100 years more, so that they reach their own regime, and then fifty years
  (five passes) saved daily and checkpointed after each pass: at 0-15 the oligotrophic and mesotrophic states at all
  three depths and the eutrophic state at 1,925 m, at the other bands the oligotrophic and mesotrophic states at
  1,925 m, and the 1,925 m lake's oligotrophic and mesotrophic water at 0-15 at the bracket (a new reservoir is no
  giant);
- the checks: the 352 m lake with hour steps, the 30 m lake with two-hour steps and the 12 m lake with layers of
  0.25 m, at 0-15;
- the lakes held meromictic (meromictic.py): the 30 m and 50 m lakes at 0-15 with a chemocline at 20 m worth
  0.1 kg/m3 of density, and the 30 m lake's mesotrophic and eutrophic water with one worth 1 kg/m3, fifty years each
  from an oxygenated start, saved daily.
It writes research/runs/lake_kinds/ladder.json, the index run.py reads, with the sha256 of every saved file and of
the lake model.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from biosphere.ecology import lake as L
from research.studies.lake_kinds import columns as C
from research.studies.lake_kinds import meromictic as M
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INDEX = C.RUNS / 'ladder.json'
SCHEMA = 'terluna.research.lake-kinds-ladder/1'
LADDER_M = (5.0, 8.0, 12.0, 16.0, 20.0, 30.0, 50.0, 100.0, 352.0)
GIANTS_M = (1000.0, 1925.0, 3500.0)
GIANT_PASSES = 5
# the giants run: (depth, band) -> trophic states, at the central wind share; a new reservoir is no giant
GIANT_CASES = {(1000.0, 0): ('oligotrophic', 'mesotrophic'),
               (1925.0, 0): ('oligotrophic', 'mesotrophic', 'eutrophic'),
               (3500.0, 0): ('oligotrophic', 'mesotrophic'),
               (1925.0, 1): ('oligotrophic', 'mesotrophic'), (1925.0, 2): ('oligotrophic', 'mesotrophic')}
ETA_TROPHIC = ('oligotrophic', 'mesotrophic')        # the states run at the wind share's bracket
ETA_LADDER_M = (352.0,)                              # the ladder's lakes run at the bracket
ETA_GIANTS = ((1925.0, 0),)                          # and the giants
ETAS = ('low', 'high')
FINE_MIN_M = 100.0                      # lakes this deep settle at the analysis step after their coarse spin-up
MEROMICTIC = dict(band=0, step_depth_m=20.0, passes=5,                  # a chemocline at 20 m, five passes, and
                  steps=((30.0, 0.1, tuple(R.TROPHIC)), (50.0, 0.1, tuple(R.TROPHIC)),   # depth, step (kg/m3) and
                         (30.0, 1.0, ('mesotrophic', 'eutrophic'))))                     # the trophic states
FINE_MAX_M = 50.0                       # lakes to this depth step every hour, as lake_oxygen's shallows; deeper ones
FINE, DEEP = R.STEPS_PER_MEAN, R.STEPS_PER_MEAN_DEEP    # every two hours, as its deep lakes
TROPHIC = tuple(R.TROPHIC)
BANDS = (0, 1, 2)


def ladder_case(band, depth, trophic, eta='central', **extra):
    case = dict(band=band, depth_m=depth, coupled=f'{trophic}/renewal/supply',
                steps_per_mean=FINE if depth <= FINE_MAX_M else DEEP, passes=1, spin=True, fine=depth >= FINE_MIN_M)
    if eta != 'central':
        case['eta'] = eta
    return dict(case, **extra)


def ladder_name(band, depth, trophic, eta='central'):
    return f'ladder_{depth:g}m_{BAND_NAMES[band]}_{trophic}' + ('' if eta == 'central' else f'_eta_{eta}')


def giant_name(band, depth, trophic, eta='central'):
    return f'giant_{depth:g}m_{BAND_NAMES[band]}_{trophic}' + ('' if eta == 'central' else f'_eta_{eta}')


BAND_NAMES = ('0-15', '15-30', '30-45')


def cases():
    """Every case of the ladder, the giants and the checks: name -> (case, the case it starts from)."""
    out = {}
    for depth in sorted(LADDER_M, reverse=True):
        for band in BANDS:
            for trophic in TROPHIC:
                etas = ('central',) + (ETAS if depth in ETA_LADDER_M and trophic in ETA_TROPHIC else ())
                for eta in etas:
                    out[ladder_name(band, depth, trophic, eta)] = (ladder_case(band, depth, trophic, eta), None)
    for (depth, band), trophics in GIANT_CASES.items():
        for trophic in trophics:
            etas = ('central',) + (ETAS if (depth, band) in ETA_GIANTS and trophic in ETA_TROPHIC else ())
            for eta in etas:
                start = ladder_name(band, 352.0, trophic, eta)
                case = dict(band=band, depth_m=depth, coupled=f'{trophic}/renewal/supply', steps_per_mean=DEEP,
                            passes=GIANT_PASSES, spin='long', fine=True, start=start)
                if eta != 'central':
                    case['eta'] = eta
                out[giant_name(band, depth, trophic, eta)] = (case, start)
    out['check_352m_0-15_mesotrophic_hour_steps'] = (ladder_case(0, 352.0, 'mesotrophic', steps_per_mean=FINE), None)
    out['check_30m_0-15_mesotrophic_two_hour_steps'] = (ladder_case(0, 30.0, 'mesotrophic', steps_per_mean=DEEP),
                                                        None)
    out['check_12m_0-15_mesotrophic_quarter_metre_layers'] = (ladder_case(0, 12.0, 'mesotrophic', top_m=0.25), None)
    return out


def meromictic_cases():
    """The lakes held meromictic: name -> (case, the ladder case it starts from)."""
    out = {}
    b = MEROMICTIC['band']
    for depth, step, trophics in MEROMICTIC['steps']:
        for trophic in trophics:
            start = ladder_name(b, depth, trophic)
            case = dict(band=b, depth_m=depth, coupled=f'{trophic}/renewal/supply', steps_per_mean=DEEP,
                        passes=MEROMICTIC['passes'], start=start, step_kg_m3=step,
                        step_depth_m=MEROMICTIC['step_depth_m'])
            out[f'meromictic_{depth:g}m_{BAND_NAMES[b]}_{trophic}_step_{step:g}'] = (case, start)
    return out


# --- the workers ------------------------------------------------------------------------------------------------

_SETTINGS = {}


def setting(ground_k=S.TARGET_K):
    if ground_k not in _SETTINGS:
        _SETTINGS[ground_k] = S.Setting301(ground_k)
    return _SETTINGS[ground_k]


def final_state(name, case):
    st = setting()
    path = C.CASES / f'{name}_{C.key(st, case)}.json'
    return json.loads(path.read_text())['final_state']


def work_column(item):
    name, case, start_name, start_case = item
    tic = time.time()
    start = final_state(start_name, start_case) if start_name else None
    res = C.run(setting(), name, case, start)
    return name, round(time.time() - tic, 1), res['case']['key']


def work_meromictic(item):
    name, case, start_name, start_case = item
    tic = time.time()
    res = M.run(setting(), name, case, final_state(start_name, start_case))
    return name, round(time.time() - tic, 1), res['case']['key']


def lake_oxygen_case(st, kind, site_i):
    """lake_oxygen's run_case() for the resort as st builds it, cached in this study's runs."""
    case = dict(kind=kind, target='resort', site_i=site_i)
    C.CASES.mkdir(parents=True, exist_ok=True)
    path = C.CASES / f'lake_oxygen_{kind}_{BAND_NAMES[site_i]}_{C.key(st, case)}.json'
    if path.exists():
        return json.loads(path.read_text()), path
    res = R.run_case(st, kind, 'resort', site_i)
    path.write_text(json.dumps(res) + '\n')
    return res, path


def work_lake_oxygen(item):
    ground_k, kind, site_i = item
    tic = time.time()
    _, path = lake_oxygen_case(setting(ground_k), kind, site_i)
    return f'lake_oxygen {kind} {BAND_NAMES[site_i]} at {ground_k:g} K', round(time.time() - tic, 1), path.name


def strip(x, drop=('seconds', 'final_state')):
    if isinstance(x, dict):
        return {k: strip(v, drop) for k, v in x.items() if k not in drop}
    if isinstance(x, list):
        return [strip(v, drop) for v in x]
    return x


def differences(a, b, path=''):
    """Paths at which two summaries differ, with both values."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append((f'{path}/{k}', a.get(k), b.get(k)))
            else:
                out += differences(a[k], b[k], f'{path}/{k}')
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            out += differences(x, y, f'{path}[{i}]')
    elif a != b:
        out.append((path, a, b))
    return out


def reproduction():
    """The resort built here at 300 K against the waters' oxygen study's committed resort runs."""
    product = json.loads(R.OUT.read_text())
    st300 = setting(S.CHECK_K)
    out = dict(construction=S.check_against_climate_targets(st300), runs={})
    for kind, site_i in [('deep', i) for i in BANDS] + [('shallow', 0)]:
        mine, path = lake_oxygen_case(st300, kind, site_i)
        theirs = product['runs'][kind][f'resort/{BAND_NAMES[site_i]}']
        diff = differences(strip(theirs), strip(mine))
        out['runs'][f'{kind}/{BAND_NAMES[site_i]}'] = dict(file=path.name, identical=not diff,
                                                         differences=[list(d) for d in diff[:20]])
    # this study's own runner on the 0-15 deep lake at 300 K, against the same committed run
    case = ladder_case(0, R.DEEP_M, 'mesotrophic', steps_per_mean=DEEP)
    mine = C.run(st300, 'check_300k_352m_0-15_mesotrophic', case)
    theirs = product['runs']['deep']['resort/0-15']
    keys = ('mixing', 'oxygen', 'bottom_t_c', 'surface_t_c', 'cycles')
    diff = differences(strip({k: theirs[k] for k in keys}), strip({k: mine[k] for k in keys}))
    out['runner'] = dict(case='check_300k_352m_0-15_mesotrophic', key=mine['case']['key'], identical=not diff,
                         differences=[list(d) for d in diff[:20]], compared=list(keys))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=2, help='worker processes, at most 2 on the shared machine')
    args = ap.parse_args(argv)
    jobs = max(1, min(args.jobs, 2))
    tic = time.time()
    st = setting()
    all_cases = cases()
    first = [(n, c, None, None) for n, (c, s) in all_cases.items() if s is None]
    second = [(n, c, s, all_cases[s][0]) for n, (c, s) in all_cases.items() if s is not None]
    mero = meromictic_cases()
    third = [(n, c, s, all_cases[s][0]) for n, (c, s) in mero.items()]
    checks = [(S.CHECK_K, kind, i) for kind, i in [('deep', i) for i in BANDS] + [('shallow', 0)]] + \
             [(S.TARGET_K, 'shallow', i) for i in BANDS]
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for label, seconds, f in pool.map(work_lake_oxygen, checks):
            print(f'  {label}: {seconds:.0f} s', flush=True)
        for work, batch in ((work_column, first), (work_column, second), (work_meromictic, third)):
            for name, seconds, k in pool.map(work, batch):
                print(f'  {name} ({k}): {seconds:.0f} s, {time.time() - tic:.0f} s in', flush=True)
    repro = reproduction()
    index = dict(schema=SCHEMA, ground_k=st.ground_k, setting=S.describe(st), lake_py_sha256=C.sha256(L.__file__),
                 producer=dict(files={str(p.resolve().relative_to(ROOT)): C.sha256(p)
                                      for p in M.code_files() + (Path(__file__),)}, inputs=st.inputs),
                 reproduction=repro, cases={}, pans={})
    for name, (case, start) in list(all_cases.items()) + list(mero.items()):
        k = M.key(st, case) if name in mero else C.key(st, case)
        path = C.CASES / f'{name}_{k}.json'
        res = json.loads(path.read_text())
        index['cases'][name] = dict(case=case, key=k, summary=path.name, summary_sha256=C.sha256(path),
                                    files=res['files'], start=start)
    for i in BANDS:
        res, path = lake_oxygen_case(st, 'shallow', i)
        index['pans'][BAND_NAMES[i]] = dict(summary=path.name, summary_sha256=C.sha256(path))
    INDEX.write_text(json.dumps(index, indent=1) + '\n')
    print('wrote', INDEX, f'in {time.time() - tic:.0f} s')
    for k, v in repro['runs'].items():
        print('reproduced' if v['identical'] else 'DIFFERS', k, v['differences'][:3])
    print('runner', 'identical' if repro['runner']['identical'] else repro['runner']['differences'][:5])
    return 0


if __name__ == '__main__':
    sys.exit(main())
