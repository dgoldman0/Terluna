"""One lake column through the resort's record at 301 K, saved for the lake kinds and the studies after it.

The lake model (biosphere/ecology/lake.py) runs through lake_oxygen's forcing, basins, physics and scenarios, with
the resort moved to 301 K by setting.Setting301. A column starts warm and saturated, 1 K above the air's mean, and
passes over the whole record in coarse steps until its deep water's temperature, oxygen, age and solute all settle
(lake_oxygen's settle test) or the passes reach their limit, as lake_oxygen's deep lakes do; or it starts from
another column carried down and oxygenated, as lake_oxygen's 1,925 m lake starts from its 352 m lake, and either runs
several passes of the record from there or first passes over the record until it settles or reaches the long limit,
so that the recorded passes show its own regime. After the coarse passes it passes over the record again at the
analysis step until the settle test passes there too (lake_oxygen's FINE_PASSES), so the recorded passes start from
that step's own settled column. Each recorded pass is saved, with the column about every day
(lake_anoxia's profiles.py's spacing and arrays, with the solutes and the labile pool added) and the step-by-step
series lake_oxygen's summarise() reads, to research/runs/lake_kinds/cases; a stopped run resumes at its last saved
pass.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'lake_kinds'
CASES = RUNS / 'cases'
PROFILE_S = 12 * 7150.0              # lake_anoxia's spacing of the saved columns: 12 steps of 7150 s, about a day
SPIN_PASSES = dict(deep=R.MAX_PASSES['deep'], long=R.MAX_PASSES['deepest_long'])   # lake_oxygen's limits, in coarse
FINE_PASSES = dict(deep=R.FINE_PASSES['deep'], long=R.FINE_PASSES['deepest_long'])   # steps and at the analysis step
SERIES = ('surface_t', 'bottom_t', 'mixed_depth', 'ice', 'saturation', 'bottom_age', 'bottom_solute', 'dissipation')
INTERFACES = ('interfaces_below', 'interfaces_capped', 'interfaces_held_by_cap', 'interfaces_unstable')
COLUMN_SERIES = ('surface_o2', 'bottom_o2', 'deep_o2', 'mean_o2', 'bottom_debt')


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def code_files():
    return (Path(__file__), Path(S.__file__), Path(R.__file__), Path(L.__file__))


def key(st: S.Setting301, case: dict) -> str:
    """The cache key: the code that runs a column, lake_oxygen's inputs, the target and the case."""
    h = hashlib.sha256()
    for p in code_files():
        h.update(sha256(p).encode())
    h.update(json.dumps(st.inputs, sort_keys=True).encode())
    h.update(json.dumps(dict(ground_k=st.ground_k, **case), sort_keys=True).encode())
    return h.hexdigest()[:20]


def basin(site, depth, pan=False, top_m=None):
    """lake_oxygen's basins: the band lake's own hypsography scaled to the depth, or its 3 m pan; top_m refines the
    layers above 20 m (lake_oxygen's 0.5 m by default)."""
    if pan:
        return R.basin('shallow', site)
    if top_m is None:
        return R.basin('deep', site, float(depth))
    lake = site['lake']
    h = lake['hypsography']
    return L.Basin.from_hypsography(float(depth), lake['area_km2'] * 1e6, h['depth_fraction'], h['area_fraction'],
                                    top_m=top_m, max_m=8.0 if depth <= 500.0 else 25.0)


def saturated(st, t_c):
    return float(L.oxygen_saturation(t_c, R.DESIGN_PRESSURE_PA / 101325.0, st.salinity, st.o2_fraction))


def capture(lake, forcing, state, every):
    """One pass of the record with the column saved every `every` steps: lake_anoxia's profiles.py's arrays and the
    solutes, with the run's own record for lake_oxygen's summarise()."""
    at = list(range(0, len(forcing), every))
    out = lake.run(forcing, state, profiles_at=at)
    prof = out.pop('profiles')
    keys = sorted(prof)
    rec = dict(step=np.array(keys), t=np.array([prof[k]['t'] for k in keys], dtype=np.float32),
               oxygen=np.array([prof[k]['oxygen'] for k in keys], dtype=np.float32),
               debt=np.array([prof[k]['debt'] for k in keys], dtype=np.float32),
               age=np.array([prof[k]['age'] for k in keys], dtype=np.float32),
               solute=np.array([prof[k]['solute'] for k in keys], dtype=np.float32),
               pool=np.array([prof[k]['pool'] for k in keys], dtype=np.float32))
    del prof
    return rec, out


def spin_up(lake, forcing, state, col, limit, fine_limit):
    """lake_oxygen's spin-up of a deep lake: passes over the whole record in 4-step means until the deep water's mean
    temperature, the coupled column's deepest-tenth oxygen, the bottom water's age and its solute all settle
    (lake_oxygen's settled()), or the passes reach the limit, then passes at the analysis step until they settle
    there too, or fine_limit; each pass's overturn is recorded."""
    out = R.settle_passes(lake, R.coarse(forcing), state, col, lake.basin, limit)
    out['at_the_analysis_step'] = R.settle_passes(lake, forcing, state, col, lake.basin, fine_limit) \
        if fine_limit else None
    return out


def mixing_record(ph, metas, dissipation):
    """lake_oxygen's deep-mixing record over the saved passes: the treatment, the wind's share and the dissipation it
    gave, and the shares of the interfaces below the mixed layer at the cap, held back by it and unstable."""
    counts = {i: sum(m['interfaces'][i] for m in metas) for i in INTERFACES}
    return R.deep_mixing_record(ph, dict(counts, dissipation=dissipation))


def state_dict(bsn, state):
    return dict(centre=bsn.centre.tolist(), t=state.t.tolist(), oxygen=state.oxygen.tolist(),
                debt=state.debt.tolist(), age=state.age.tolist(), solute=state.solute.tolist(),
                pool=state.pool.tolist(), mixed_depth=state.mixed_depth, ice=state.ice,
                ice_surface_c=state.ice_surface_c)


def state_from_dict(saved):
    s = L.State(np.asarray(saved['t'], float), np.asarray(saved['oxygen'], float), np.asarray(saved['debt'], float),
                np.asarray(saved['age'], float), solute=np.asarray(saved['solute'], float),
                pool=np.asarray(saved['pool'], float) if 'pool' in saved else None)
    s.mixed_depth, s.ice, s.ice_surface_c = float(saved['mixed_depth']), float(saved['ice']), \
        float(saved['ice_surface_c'])
    return s


def run(st: S.Setting301, name: str, case: dict, start=None):
    """Run one case (or resume it) and return its summary. case: band (0, 1, 2), depth_m, coupled (the oxygen column
    whose trophic state sets the physics), steps_per_mean, passes, spin (True or 'long'), fine (whether the spin-up
    goes on at the analysis step), and optionally pan, top_m, eta (the wind's share), mixing, dissipation, cap and
    solute (lake_oxygen's names) and start, the name of the case whose final column, given as start, it begins from,
    carried down and oxygenated."""
    CASES.mkdir(parents=True, exist_ok=True)
    k = key(st, case)
    final = CASES / f'{name}_{k}.json'
    if final.exists():
        return json.loads(final.read_text())
    tic = time.time()
    site = R.site_info(st, case['band'])
    spm = case['steps_per_mean']
    w, forcing, mu = R.forcing_for(st, 'resort', site, 'lock', None, spm)
    bsn = basin(site, case['depth_m'], case.get('pan', False), case.get('top_m'))
    scen = R.scenarios()
    names = [s.name for s in scen]
    ph = R.physics(st, 'deep', mixing=case.get('mixing', R.DEEP_MIXING), eta=case.get('eta', 'central'),
                   dissipation=case.get('dissipation', 'central'), solute=case.get('solute', 'co2'),
                   coupled=case['coupled'], cap=case.get('cap'))
    lake = L.Lake(bsn, ph, scen, st.air)
    lake.calibrate_production(forcing)
    col = ph.solute_column
    every = max(1, int(round(PROFILE_S / forcing.dt)))
    done = []                                            # passes saved whole: the column and its record
    while (CASES / f'{name}_{k}_pass{len(done):02d}.json').exists():
        done.append(CASES / f'{name}_{k}_pass{len(done):02d}.npz')
    spun = None
    if done:                                             # resume after the last saved pass
        last = done[-1]
        state = state_from_dict(json.loads(last.with_suffix('.json').read_text())['state'])
        spun = json.loads(last.with_suffix('.json').read_text())['spin_up']
        first = len(done)
    else:
        first = 0
        if start is not None:
            state = R.state_from(start, bsn, scen)
            state.oxygen[:] = saturated(st, float(state.t.mean()))     # oxygenated, as if just overturned
            state.debt[:], state.age[:], state.solute[:], state.pool[:] = 0.0, 0.0, 0.0, 0.0
            start_c = None
        else:
            start_c = float(forcing.air_k.mean()) - 273.15 + R.WARM_START_K
            state = L.initial_state(bsn, scen, start_c, saturated(st, start_c))
        if case['spin']:
            which = 'long' if case['spin'] == 'long' else 'deep'
            spun = spin_up(lake, forcing, state, col, SPIN_PASSES[which],
                           FINE_PASSES[which] if case.get('fine', True) else 0)
        else:
            spun = dict(passes=0, settled=None, at_the_analysis_step=None)
        spun['start_c'] = None if start_c is None else round(start_c, 2)
    for p in range(first, case['passes']):
        h0, o0 = L.heat_content(bsn, state), L.oxygen_content(bsn, state)
        rec, out = capture(lake, forcing, state, every)
        budget = dict(R.budget_errors(bsn, state, h0, o0, out, len(forcing) * forcing.dt),
                      produced=out['produced'].tolist(), respired=out['respired'].tolist(),
                      oxygen_in=out['oxygen_in'].tolist())
        series = {s: np.asarray(out[s], dtype=float) for s in SERIES + COLUMN_SERIES}
        path = CASES / f'{name}_{k}_pass{p:02d}.npz'
        np.savez_compressed(path, z=bsn.z, area=bsn.area, mu=mu, dt=forcing.dt,
                            pressure=forcing.pressure.astype(np.float32), every=every, **rec, **series)
        path.with_suffix('.json').write_text(json.dumps(dict(pass_index=p, budget=budget, spin_up=spun,
                                                             interfaces={i: int(out[i]) for i in INTERFACES},
                                                             state=state_dict(bsn, state),
                                                             npz_sha256=sha256(path))) + '\n')
        del rec, out, series
    # the whole record, every pass in turn, summarised as lake_oxygen summarises its runs
    passes = [np.load(CASES / f'{name}_{k}_pass{p:02d}.npz') for p in range(case['passes'])]
    metas = [json.loads((CASES / f'{name}_{k}_pass{p:02d}.json').read_text()) for p in range(case['passes'])]
    out = {s: np.concatenate([z[s] for z in passes]).astype(float) for s in SERIES + COLUMN_SERIES}
    mu_all = np.concatenate([z['mu'] for z in passes]).astype(float)
    res = R.summarise('deep', out, mu_all, bsn, scen, forcing.dt)
    res['deep_mixing'] = mixing_record(ph, metas, out['dissipation'])
    days = len(mu_all) * forcing.dt / 86400.0
    produced = np.sum([m['budget']['produced'] for m in metas], 0)
    respired = np.sum([m['budget']['respired'] for m in metas], 0)
    from_air = np.sum([m['budget']['oxygen_in'] for m in metas], 0)
    res['budget'] = dict(heat_relative_error=max(abs(m['budget']['heat_relative_error']) for m in metas),
                         heat_error_w_m2=max(abs(m['budget']['heat_error_w_m2']) for m in metas),
                         oxygen_relative_error=max(m['budget']['oxygen_relative_error'] for m in metas),
                         oxygen_error_g_m2_day=max(m['budget']['oxygen_error_g_m2_day'] for m in metas),
                         production_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in produced],
                         respiration_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in respired],
                         from_the_air_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in from_air])
    year = int(round(365.25 * 86400.0 / forcing.dt))
    first2, last2 = slice(0, min(2 * year, len(mu_all))), slice(max(0, len(mu_all) - 2 * year), len(mu_all))
    res['settling'] = dict(spin_up=metas[0]['spin_up'],
                           bottom_t_first_two_years=round(float(out['bottom_t'][first2].mean()), 3),
                           bottom_t_last_two_years=round(float(out['bottom_t'][last2].mean()), 3),
                           deep_o2_first_two_years=[round(float(x), 3) for x in out['deep_o2'][first2].mean(0)],
                           deep_o2_last_two_years=[round(float(x), 3) for x in out['deep_o2'][last2].mean(0)])
    res['bottom_solute_g_m3'] = dict(mean=round(float(out['bottom_solute'].mean()), 3),
                                     max=round(float(out['bottom_solute'].max()), 3))
    res['case'] = dict(case, name=name, key=k, band_name=site['band'], lat=site['lat'], lon=site['lon'],
                       lake_level_above_sea_m=site['lake']['level_above_sea_m'], cell_elevation_m=site['elevation'],
                       depth_m=bsn.depth, layers=bsn.n, area_km2=round(bsn.area[0] / 1e6, 1),
                       mean_depth_m=round(float(bsn.volume.sum() / bsn.area[0]), 2), step_s=forcing.dt,
                       profile_every_steps=every, years=round(days / 365.25, 2), ground_k=st.ground_k,
                       share_of_the_way=st.share['resort'], lapse_k_per_km=round(st.lapse['resort'] * 1000.0, 4),
                       seconds=round(time.time() - tic, 1), scenarios=names,
                       forcing=dict(air_c=round(float(forcing.air_k.mean()) - 273.15, 2),
                                    sunlight_w_m2=round(float(forcing.sw_down.mean()), 1),
                                    sky_infrared_w_m2=round(float(forcing.lw_down.mean()), 1),
                                    wind10_m_s=round(float(forcing.wind10.mean()), 2),
                                    pressure_atm=round(float(forcing.pressure.mean()) / 101325.0, 4)))
    res['files'] = {f'{name}_{k}_pass{p:02d}.npz': metas[p]['npz_sha256'] for p in range(case['passes'])}
    res['final_state'] = metas[-1]['state']
    for z in passes:
        z.close()
    final.write_text(json.dumps(res) + '\n')
    return res
