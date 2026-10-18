"""One sea column through the resort's record at 301 K, saved as the lake kinds study saves its lakes.

The lake model (biosphere/ecology/lake.py) runs through lake_oxygen's forcing, physics and scenarios, on a run
basin's own hypsography (results/seas.json) with the held seas' salinity, forced by one of the basin's T21 water
cells in the 301 K resort (setting.SeaSetting). The physics is lake_oxygen's physics(): below the mixed layer the
wind's share of its energy flux at 10 m sets the dissipation (the author's central treatment), under a cap of
0.1 m2/s, with the Lake number's diffusivity and the first round's fixed dissipation as sensitivity rows.

A deep column starts warm and saturated, 1 K above the air's mean, as lake_oxygen's deep lakes do, and passes over
the whole record in its 4-step means SPIN_PASSES times, each pass recorded and the column saved after each so a stopped
run resumes; a variant starts from its basin's central column (lake_oxygen starts its variants from its base lakes)
and passes FROM_PASSES times. A sea's deep water takes centuries to find its regime, so the coarse passes are fixed
and recorded, and the long state is read from their last half. The column then passes over the record at the analysis
step until lake_oxygen's settle test passes there (the deep water's temperature, the coupled column's deepest-tenth
oxygen, the bottom water's age and its solute each steady from one pass to the next) or reaches lake_oxygen's limit;
a column cut off from the surface keeps ageing and does not pass. A pan (the shelves and bays) spins up over the
record's first year in hour steps, as lake_oxygen's productive shallows do. Then the record is run once more step by
step and saved with the column about every day in research/studies/lake_kinds/columns.py's arrays (its capture(),
with the labile pool), to research/runs/sea_oxygen/cases.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_kinds import columns as KC
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R
from research.studies.sea_oxygen import setting as T

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASES = T.RUNS / 'cases'
SPIN_PASSES = 80                     # a column from the warm start passes the ten-year record this many times first
FROM_PASSES = 30                     # and a variant started from its basin's central column this many
FINE_LIMIT = R.FINE_PASSES['deep']   # lake_oxygen's limit on the passes at the analysis step
EPISODE_DAYS = 15.0                  # complete overturns closer than this are one, as lake_anoxia groups them
DEEP_LAYER_MAX_M = 25.0              # lake_oxygen's layers below 20 m grow to 8 m, or 25 m in its 1,925 m lake
SERIES, COLUMN_SERIES = KC.SERIES, KC.COLUMN_SERIES
EXTRA_SERIES = ('u_star', 'dissipation')


def code_files():
    return (Path(__file__), Path(T.__file__), Path(KC.__file__), Path(S.__file__), Path(R.__file__), Path(L.__file__))


# What a column takes from the other studies' code, by name: the cache key holds their source and values, so that an
# edit elsewhere in those files (a new check, a docstring) leaves the cached columns standing, while any change to what
# the columns run re-keys them. This study's own files and the lake model enter whole.
USED = {
    R: dict(code=('kd', 'scenarios', 'window_humidity', 'windows', 'resort_windows', 'hourly', 'physics', 'nights',
                  'days', 'cycles', 'quantiles', 'summarise', 'coarse', 'forcing_for', 'lake_number', 'settled',
                  'deep_mixing_record', 'settle_passes', 'weighted_quantile', 'Setting'),
            values=('DESIGN_PRESSURE_PA', 'DESIGN_CO2_PPM', 'HELD_FRESH', 'TARGET_NAMES', 'SHALLOW_LAYER_M',
                    'STEPS_PER_MEAN', 'STEPS_PER_MEAN_DEEP', 'LOWEST_FLOOR_MS', 'SOLAR_HOUR_BINS', 'HYPOXIC', 'ANOXIC',
                    'G_C_TO_O2', 'COARSE', 'WARM_START_K', 'SETTLED_K', 'SETTLED_O2', 'SETTLED_AGE', 'SETTLED_SOLUTE',
                    'MAX_PASSES', 'FINE_PASSES', 'AFTER_OVERTURN_DAYS', 'TROPHIC', 'TRANSFER', 'KINETICS', 'SOURCE_C',
                    'LABILE_DAYS', 'METHANE_SHARE', 'IRON', 'WIND_SHARE', 'ETA_BRACKET', 'DEEP_MIXING', 'CAP_WIND',
                    'DISSIPATION', 'CAP', 'UNSTABLE_DIFFUSIVITY', 'LAKE_NUMBER_BELOW_M', 'SOLUTE', 'COUPLED')),
    KC: dict(code=('capture', 'state_dict', 'state_from_dict', 'saturated'),
             values=('PROFILE_S', 'SERIES', 'COLUMN_SERIES')),
    S: dict(code=('interpolation', 'Setting301'), values=('TARGET_K', 'LAPSE_LOW_KM', 'LAPSE_HIGH_KM')),
}


def _fingerprint():
    """sha256 of the source of every function and class the columns use from lake_oxygen and the lake kinds study,
    and of every constant they read there, as this process loaded them."""
    import inspect
    h = hashlib.sha256()
    for module, used in USED.items():
        h.update(module.__name__.encode())
        for name in used['code']:
            h.update(name.encode())
            h.update(inspect.getsource(getattr(module, name)).encode())
        for name in used['values']:
            h.update(name.encode())
            h.update(json.dumps(getattr(module, name), sort_keys=True, default=repr).encode())
    return h.hexdigest()


# Read once, as the modules are imported: inspect reads the source files from disk, so a later edit to one of them
# would otherwise be read against the code this process runs.
FINGERPRINT = _fingerprint()


def fingerprint():
    return FINGERPRINT


def key(st: T.SeaSetting, case: dict) -> str:
    """The cache key: this study's column code and the lake model whole, the fingerprint of what the columns use
    from lake_oxygen and the lake kinds study, the study's inputs, the target and the case."""
    h = hashlib.sha256()
    for p in (Path(__file__), Path(T.__file__), Path(L.__file__)):
        h.update(KC.sha256(p).encode())
    h.update(fingerprint().encode())
    h.update(json.dumps(st.inputs, sort_keys=True).encode())
    h.update(json.dumps(dict(ground_k=st.ground_k, **case), sort_keys=True).encode())
    return h.hexdigest()[:20]


def basin(st: T.SeaSetting, case: dict) -> L.Basin:
    """The run basin's hypsography to its deepest point, or a pan of the case's depth."""
    if case.get('pan_m'):
        return L.Basin.pan(float(case['pan_m']), R.SHALLOW_LAYER_M, 1e6)
    b = st.basins[case['basin']]
    h = b['hypsography']
    depth = float(b['deepest_m'])
    return L.Basin.from_hypsography(depth, b['area_km2'] * 1e6, h['depth_fraction'], h['area_fraction'],
                                    max_m=8.0 if depth <= 500.0 else DEEP_LAYER_MAX_M)


def physics(st, case):
    """lake_oxygen's physics(), with the case's deep mixing, wind share, fixed dissipation, solutes and coupled
    column, and its cap (0.1 m2/s under the wind and the Lake number, 1e-4 under the fixed dissipation)."""
    kind = 'shallow' if case.get('pan_m') else 'deep'
    return R.physics(st, kind, mixing=case.get('mixing', R.DEEP_MIXING), eta=case.get('eta', 'central'),
                     dissipation=case.get('dissipation', 'central'), solute=case.get('solute', 'co2'),
                     coupled=case.get('coupled'), cap=case.get('cap'))


def recorded_passes(lake, forcing, state, path, passes):
    """A fixed number of passes over the record in 4-step means, each recorded and the column saved after each to path
    (and resumed from it). A pass records the bottom water's mean temperature and solute and the deepest tenth's
    temperature at its end, the deepest tenth's oxygen in every column (its mean, its least, and the shares of the pass
    under 2 and 0.5 g/m3), the bottom water's mean age and the deepest tenth's at the end, the complete overturns (the
    mixed layer reaching the bottom layer, steps within EPISODE_DAYS one episode) and the deepest mixing."""
    spin = R.coarse(forcing)
    bsn = lake.basin
    v = bsn.volume
    cum = np.cumsum(v)
    deep = (cum - 0.5 * v) >= 0.9 * cum[-1]
    wd = np.where(deep, v, 0.0) / v[deep].sum()
    bottom_top = float(bsn.z[-2])
    gap = EPISODE_DAYS * 86400.0 / spin.dt
    history = []
    if path.exists():
        with np.load(path) as z:
            state.t[:], state.oxygen[:], state.debt[:] = z['t'], z['oxygen'], z['debt']
            state.age[:], state.solute[:], state.pool[:] = z['age'], z['solute'], z['pool']
            state.mixed_depth, state.ice, state.ice_surface_c = float(z['mixed_depth']), float(z['ice']), \
                float(z['ice_surface_c'])
            history = json.loads(str(z['history']))
    while len(history) < passes:
        o = lake.run(spin, state)
        d = o['deep_o2']
        full = np.flatnonzero(o['mixed_depth'] >= bottom_top - 1e-6)
        events = int(1 + np.count_nonzero(np.diff(full) > gap)) if full.size else 0
        history.append(dict(bottom_t=float(o['bottom_t'].mean()), deep_t_end=float(state.t @ wd),
                            bottom_solute=float(o['bottom_solute'].mean()),
                            deep_o2=[float(x) for x in d.mean(0)], deep_o2_min=[float(x) for x in d.min(0)],
                            deep_hypoxic=[float(x) for x in (d < R.HYPOXIC).mean(0)],
                            deep_anoxic=[float(x) for x in (d < R.ANOXIC).mean(0)],
                            bottom_age_days=float(o['bottom_age'].mean()), deep_age_days_end=float(state.age @ wd),
                            overturns=events, deepest_mixing_m=float(o['mixed_depth'].max())))
        np.savez(path, t=state.t, oxygen=state.oxygen, debt=state.debt, age=state.age, solute=state.solute,
                 pool=state.pool, mixed_depth=state.mixed_depth, ice=state.ice, ice_surface_c=state.ice_surface_c,
                 history=json.dumps(history))
    lists = {k: [h[k] for h in history] for k in history[0]} if history else {}
    return dict(passes=len(history), steps_per_pass=len(spin), step_s=spin.dt,
                by_pass={k: [_r(x) for x in vals] for k, vals in lists.items()})


def _r(x, nd=4):
    if isinstance(x, list):
        return [round(float(y), nd) for y in x]
    return round(float(x), nd) if isinstance(x, float) else x


def lake_numbers(st, bsn, rec, out, dt):
    """The Lake number of each saved column, as lake_oxygen reckons it (its lake_number()), with the day's mean u*^2,
    by the lunar month's thermocline and the day's warm layer."""
    day = int(round(86400.0 / dt))
    steps_at = np.asarray(rec['step'])
    u2 = [float((out['u_star'][k:k + day] ** 2).mean()) for k in steps_at]
    temps = [np.asarray(t, float) for t in rec['t']]
    res = {}
    for rule, below_m in R.LAKE_NUMBER_BELOW_M.items():
        ln = R.lake_number(bsn, temps, u2, st.gravity, below_m)
        res[rule] = dict(z_t_below_m=below_m, median=round(float(np.median(ln)), 4),
                         p10=round(float(np.quantile(ln, 0.1)), 4), p90=round(float(np.quantile(ln, 0.9)), 4),
                         share_under_1=round(float((ln < 1.0).mean()), 4), samples=int(ln.size))
    return res


def run(st: T.SeaSetting, name: str, case: dict, start=None):
    """Run one case (or resume it) and return its summary. case: basin (a run basin's key), cell (the site index),
    column ('typical' or 'coldest'), coupled (the oxygen column whose trophic state sets the physics),
    steps_per_mean, and optionally pan_m (a pan of that depth in place of the basin), mixing, eta, dissipation, cap
    and solute (lake_oxygen's names), passes (the recorded coarse passes) and start (the name of the case whose final
    column, given as start, it begins from)."""
    CASES.mkdir(parents=True, exist_ok=True)
    k = key(st, case)
    final = CASES / f'{name}_{k}.json'
    if final.exists():
        return json.loads(final.read_text())
    tic = time.time()
    site = st.site(case['cell'], case['basin'])
    spm = case['steps_per_mean']
    w, forcing, mu = R.forcing_for(st, 'resort', site, 'lock', None, spm)
    bsn = basin(st, case)
    scen = R.scenarios()
    names = [s.name for s in scen]
    ph = physics(st, case)
    lake = L.Lake(bsn, ph, scen, st.air)
    lake.calibrate_production(forcing)
    col = ph.solute_column
    every = max(1, int(round(KC.PROFILE_S / forcing.dt)))
    pass_path = CASES / f'{name}_{k}_pass00.npz'
    if pass_path.with_suffix('.json').exists():
        meta0 = json.loads(pass_path.with_suffix('.json').read_text())
    else:
        if start is not None:
            state = KC.state_from_dict(start)
            start_c = None
        else:
            start_c = float(forcing.air_k.mean()) - 273.15 + R.WARM_START_K
            state = L.initial_state(bsn, scen, start_c, KC.saturated(st, start_c))
        if case.get('pan_m'):
            spin = R.hourly(st, w, site, 0, min(120, w['air'].size), spm)[0]
            spun = R.settle_passes(lake, spin, state, col, bsn, R.MAX_PASSES['shallow'])
            spun['at_the_analysis_step'] = None
        else:
            spun = recorded_passes(lake, forcing, state, CASES / f'{name}_{k}_spin.npz', case['passes'])
            spun['at_the_analysis_step'] = R.settle_passes(lake, forcing, state, col, bsn, FINE_LIMIT)
        spun['start_c'] = None if start_c is None else round(start_c, 2)
        spun['started_from'] = case.get('start')
        h0, o0 = L.heat_content(bsn, state), L.oxygen_content(bsn, state)
        rec, out = KC.capture(lake, forcing, state, every)
        heat_err = (L.heat_content(bsn, state) - h0 - out['heat_in']) / max(abs(out['heat_in']), abs(h0) * 1e-6)
        ox_err = (L.oxygen_content(bsn, state) - o0 - out['oxygen_in'] - out['produced'] + out['respired'])
        budget = dict(heat_relative_error=float(heat_err),
                      oxygen_relative_error=float(np.max(np.abs(ox_err) / np.maximum(np.abs(o0), 1.0))),
                      produced=out['produced'].tolist(), respired=out['respired'].tolist(),
                      oxygen_in=out['oxygen_in'].tolist())
        extra = dict(deep_mixing=R.deep_mixing_record(ph, out),
                     lake_number=None if case.get('pan_m') else lake_numbers(st, bsn, rec, out, forcing.dt))
        series = {s: np.asarray(out[s], dtype=float) for s in SERIES + COLUMN_SERIES}
        series.update({s: np.asarray(out[s], dtype=np.float32) for s in EXTRA_SERIES})
        np.savez_compressed(pass_path, z=bsn.z, area=bsn.area, mu=mu, dt=forcing.dt,
                            pressure=forcing.pressure.astype(np.float32), every=every, **rec, **series)
        pass_path.with_suffix('.json').write_text(json.dumps(dict(pass_index=0, budget=budget, spin_up=spun,
                                                                  state=KC.state_dict(bsn, state), extra=extra,
                                                                  npz_sha256=KC.sha256(pass_path))) + '\n')
        del rec, out, series
        meta0 = json.loads(pass_path.with_suffix('.json').read_text())
    with np.load(pass_path) as z:
        out = {s: np.asarray(z[s], dtype=float) for s in SERIES + COLUMN_SERIES}
        mu_all = np.asarray(z['mu'], dtype=float)
    kind = 'shallow' if case.get('pan_m') else 'deep'
    res = R.summarise(kind, out, mu_all, bsn, scen, forcing.dt)
    days = len(mu_all) * forcing.dt / 86400.0
    b = meta0['budget']
    res['budget'] = dict(heat_relative_error=b['heat_relative_error'], oxygen_relative_error=b['oxygen_relative_error'],
                         production_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in b['produced']],
                         respiration_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in b['respired']],
                         from_the_air_g_m2_day=[round(float(x) / bsn.area[0] / days, 3) for x in b['oxygen_in']])
    year = int(round(365.25 * 86400.0 / forcing.dt))
    first2, last2 = slice(0, min(2 * year, len(mu_all))), slice(max(0, len(mu_all) - 2 * year), len(mu_all))
    res['settling'] = dict(spin_up=meta0['spin_up'],
                           bottom_t_first_two_years=round(float(out['bottom_t'][first2].mean()), 3),
                           bottom_t_last_two_years=round(float(out['bottom_t'][last2].mean()), 3),
                           deep_o2_first_two_years=[round(float(x), 3) for x in out['deep_o2'][first2].mean(0)],
                           deep_o2_last_two_years=[round(float(x), 3) for x in out['deep_o2'][last2].mean(0)])
    res['bottom_solute_g_m3'] = dict(mean=round(float(out['bottom_solute'].mean()), 3),
                                     max=round(float(out['bottom_solute'].max()), 3))
    res.update(meta0['extra'])
    res['case'] = dict(case, name=name, key=k, kind=kind, lat=site['lat'], lon=site['lon'], row=site['row'],
                       t21_column=site['column'], land_cell=site['land_cell'], eddy_factor=site['eddy_factor'],
                       depth_m=bsn.depth, layers=bsn.n, area_km2=round(bsn.area[0] / 1e6, 1),
                       mean_depth_m=round(float(bsn.volume.sum() / bsn.area[0]), 2), step_s=forcing.dt,
                       profile_every_steps=every, years=round(days / 365.25, 2), ground_k=st.ground_k,
                       share_of_the_way=st.share['resort'], seconds=round(time.time() - tic, 1), scenarios=names,
                       physics=dict(deep_mixing=ph.deep_mixing, wind_share=ph.wind_share, dissipation=ph.dissipation,
                                    cap=ph.max_diffusivity, unstable=ph.unstable_diffusivity,
                                    solute_density=ph.solute_density, debt_density=ph.debt_density,
                                    attenuation=ph.attenuation, coupled=names[ph.solute_column]),
                       pairing='lowest level' if forcing.wind_level is not None else '10 m',
                       forcing=dict(air_c=round(float(forcing.air_k.mean()) - 273.15, 2),
                                    sunlight_w_m2=round(float(forcing.sw_down.mean()), 1),
                                    sky_infrared_w_m2=round(float(forcing.lw_down.mean()), 1),
                                    wind10_m_s=round(float(forcing.wind10.mean()), 2),
                                    wind_level_m_s=round(float(np.mean(forcing.wind_level)), 2)
                                    if forcing.wind_level is not None else None,
                                    pressure_atm=round(float(forcing.pressure.mean()) / 101325.0, 4),
                                    gcm_surface_c=round(float(T.resort_surface_c(st, case['cell']).mean()), 2)))
    res['files'] = {pass_path.name: meta0['npz_sha256']}
    res['final_state'] = meta0['state']
    final.write_text(json.dumps(res) + '\n')
    return res


def slab(st: T.SeaSetting, i: int, run_name: str):
    """lake_oxygen's slab check (its run.slab_check) at a sea cell: a 50 m mixed slab forced as the columns are, from
    one run's own 3-day means, against the GCM's slab ocean there; with the exchange paired as PlaSim pairs it (the
    lowest level's wind and height, as the columns take it) and with the 10 m wind and coefficient, each with the free
    convection over water carried to lunar gravity, as the columns take it, and with Earth's constant, as PlaSim takes
    it. Cached."""
    CASES.mkdir(parents=True, exist_ok=True)
    k = key(st, dict(slab=run_name, cell=i))
    path = CASES / f'slab_{run_name}_{i}_{k}.json'
    if path.exists():
        return json.loads(path.read_text())
    site = st.site(i)
    w = R.windows(st, run_name, site, at_lake=False)
    bsn = L.Basin.pan(50.0, 50.0, 1.0)
    scen = [L.Scenario('none', 0.0, 0.5, 0.0, 0.0)]
    albedo = float(np.median(w['gcm_alb'][w['sw'] > 50.0])) if np.any(w['sw'] > 50.0) else 0.069
    gcm = w['gcm_ts'] - 273.15
    res = dict(cell=i, row=site['row'], column=site['column'], lat=site['lat'], lon=site['lon'], run=run_name,
               albedo=albedo, gcm_mean_c=round(float(gcm.mean()), 2), gcm_range_k=round(float(np.ptp(gcm)), 2))
    for pairing, paired in (('lowest_level', True), ('ten_metres', False)):
        forcing, _, _ = R.hourly(st, w, site, paired=paired)
        spin, _, _ = R.hourly(st, w, site, 0, 120, paired=paired)
        res[pairing] = {}
        for name, air in (('lunar_free_convection', st.air),
                          ('earth_free_convection', dataclasses.replace(st.air, free_convection_gravity=9.80665))):
            # lake_oxygen's slab physics, its solutes carried in the slab's one oxygen column (physics() points them
            # at the mesotrophic column of the full scenarios, which the slab does not hold)
            ph = dataclasses.replace(R.physics(st, 'deep', full_mixing=True, albedo=albedo, solute='none'),
                                     solute_column=0)
            lake = L.Lake(bsn, ph, scen, air)
            state = L.initial_state(bsn, scen, float(gcm[0]), 8.0)
            lake.run(spin, state)
            out = lake.run(forcing, state)
            model = out['surface_t'].reshape(-1, R.STEPS_PER_MEAN).mean(1)
            res[pairing][name] = dict(model_mean_c=round(float(model.mean()), 2),
                                      bias_k=round(float((model - gcm).mean()), 2),
                                      rms_k=round(float(np.sqrt(((model - gcm) ** 2).mean())), 2),
                                      correlation=round(float(np.corrcoef(model, gcm)[0, 1]), 3),
                                      model_range_k=round(float(np.ptp(model)), 2))
    res['key'] = k
    path.write_text(json.dumps(res) + '\n')
    return res
