"""The saved lake columns the methane tracer reads, with each run's surface transfer rebuilt from its own forcing.

The lake model (biosphere/ecology/lake.py) is not run here. Three kinds of saved run are read:
- lake_anoxia's columns (research/runs/lake_anoxia, with profiles.json): the waters' oxygen study's 352 m lakes in
  the lock, the resort at 300 K and the hot world, saved daily over ten years, and its 1,925 m lake saved every 30
  days over fifty; their physics follows the mesotrophic column and they carry no solute, which is rebuilt here as
  the coupled column's oxygen used since the water left the surface;
- the lake kinds' columns at 301 K (research/runs/lake_kinds, its ladder.json index and each case's passes), each
  trophic state with its own physics, with the solute saved;
- the seas' columns (research/runs/sea_oxygen), in the same format as the kinds'.

Each Column carries the basin, the model's physics and scenarios, the saved columns of the traced oxygen columns, the
mixed layer's depth at every step and the surface-renewal transfer velocity of O2 at every step. The velocity is the
model's own (lake.transfer_velocity with lake.bulk_fluxes, the same inputs in the same order): the bulk fluxes at the
surface temperature the step starts from, the infrared and the Schmidt number at the one it ends at, from the run's
forcing as lake_oxygen builds it (or the kinds' Setting301 at 301 K).
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.lake_oxygen import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'lake_methane'
ANOXIA_RUNS = ROOT / 'research' / 'runs' / 'lake_anoxia'
DEV_COLUMNS = RUNS / 'dev_columns'       # lake_anoxia's columns from the uncorrected model, copied for development
KINDS_RUNS = ROOT / 'research' / 'runs' / 'lake_kinds'
SEAS_RUNS = ROOT / 'research' / 'runs' / 'sea_oxygen'
BANDS = ('0-15', '15-30', '30-45')
CENTRAL = 'renewal/supply'
# The lake model as committed before the review of 18 October: the runs it made (lake_anoxia's columns, the copy the
# tracer was developed on, the kinds' first 301 K set) took a fixed euphotic respiration, no density from the debt and
# the bulk exchange at 10 m. Runs made since take the fixed model's own scenarios, physics and pairing.
PRE_FIX_LAKE_SHA256 = 'd99e5497498eee1eff6b4955778e785a8f6ffc5c0b6531843e6a3ef57b90d751'
PRE_FIX_SOLUTE = dict(none=0.0, co2=44.01 / 32.0 * 0.284e-3, carbonate=162.11 / 32.0 * 0.8e-3)
PRE_FIX_DISSIPATION = dict(low=1e-11, central=1e-10, high=1e-9)


def model_version(lake_sha256):
    return 'pre_fix' if lake_sha256 == PRE_FIX_LAKE_SHA256 else 'current'


def scenarios_for(names, version):
    """The scenarios a run was made with, in its own order: lake_oxygen's current ones by name, or for a pre-fix run
    the same trophic states, transfers and kinetics with the euphotic respiration at a fixed rate and the Arrhenius
    dependence from 20 C, as they were then."""
    if version == 'current':
        by = {x.name: x for x in R.scenarios()}
        return [by[n] for n in names]
    out = []
    for n in names:
        trophic, transfer, kinetics = n.split('/')[:3]
        p = R.TROPHIC[trophic]
        out.append(L.Scenario(n, gpp=p['gpp'], attenuation=R.kd(p['secchi_m']), sediment=p['sediment'],
                              deep_water=p['deep_water'], transfer=transfer, activation_ev=R.KINETICS[kinetics]))
    return out


def physics_for(st, coupled, version, names, mixing=None, eta='central', dissipation='central', solute='co2',
                cap=None):
    """A run's physics: lake_oxygen's current physics() with the run's own deep mixing, wind share, dissipation, solutes
    and cap (its defaults: the wind's 0.3% under a cap of 0.1 m2/s), or for a pre-fix run its physics then (the
    coupled state's clarity, the aerobic CO2 in the density and no more, the fixed dissipation under 1e-4 m2/s)."""
    if version == 'current':
        ph = R.physics(st, 'deep', mixing=mixing or R.DEEP_MIXING, eta=eta, dissipation=dissipation, solute=solute,
                       coupled=coupled, cap=cap)
    else:
        trophic = coupled.split('/')[0]
        ph = L.Physics(gravity=st.gravity, salinity=st.salinity, attenuation=R.kd(R.TROPHIC[trophic]['secchi_m']),
                       dissipation=PRE_FIX_DISSIPATION[dissipation], heat_flow=st.heat_flow, albedo=0.069,
                       solute_density=PRE_FIX_SOLUTE[solute], solute_column=0, deep_mixing='dissipation',
                       max_diffusivity=1e-4, unstable_diffusivity=1e-4)
    import dataclasses
    return dataclasses.replace(ph, solute_column=names.index(coupled))


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


@dataclass
class Column:
    """One saved lake run, read for the tracer. Snapshot s is the column after model step snap_step[s]; series are
    one value a step. oxygen and debt hold only the traced columns (names)."""
    name: str
    basin: L.Basin
    physics: L.Physics
    scenarios: list                      # the traced columns' scenarios, in the order of names
    names: list
    dt: float
    snap_step: np.ndarray
    t: np.ndarray                        # (S, n) C
    oxygen: np.ndarray                   # (S, n, J) g/m3
    debt: np.ndarray                     # (S, n, J) g O2/m3
    age: np.ndarray                      # (S, n) days
    solute: np.ndarray                   # (S, n) g/m3 of O2 respired in the coupled column
    solute_saved: bool
    mixed_depth: np.ndarray              # (N,) m
    surface_t: np.ndarray                # (N,) C
    pressure: np.ndarray                 # (N,) Pa
    k_o2: np.ndarray                     # (N,) m/s
    open_water: np.ndarray               # (N,) bool
    ml_o2: np.ndarray = None             # (N, J) the mixed layer's oxygen after each step, g/m3
    coupled_debt: np.ndarray = None      # (S, n) the debt of the column whose products set the density
    pool: np.ndarray = None              # (S, n, J) the labile pool, where the run saved it
    dissipation: np.ndarray = None       # (N,) the dissipation below the mixed layer each step, W/kg (wind mixing)
    u_star: np.ndarray = None            # (N,) the water's friction velocity each step, m/s, where the run saved it
    u_star_rebuilt: np.ndarray = None    # (N,) the same from the rebuilt forcing, for the check
    version: str = 'pre_fix'
    form: str = 'osborn'                 # the diffusivity's form below the mixed layer: lake.py's Osborn rule, or
                                         # 'shih_bouffard' (lake_oxygen's mixing_forms)
    meta: dict = field(default_factory=dict)

    @property
    def steps(self):
        return int(self.mixed_depth.size)

    def snapshot_index(self, n):
        """The snapshot that stands for model step n: the last one taken before it, or the first."""
        return max(int(np.searchsorted(self.snap_step, n, side='left')) - 1, 0)


def rebuilt_solute(st, t, oxygen_coupled, salinity, o2_fraction):
    """The coupled column's solute where it was not saved: the oxygen the water has used since it left the surface,
    saturation at its temperature less its oxygen. The model's solute grows by what respiration and the debt's
    reoxidation take and falls by what production gives, as the oxygen does the other way, and both relax at the
    surface at one rate, so below the surface their sum stays the oxygen the water left the surface with."""
    sat = L.oxygen_saturation(t, R.DESIGN_PRESSURE_PA / 101325.0, salinity, o2_fraction)
    return (sat - oxygen_coupled).astype(np.float32)


def transfer_velocity_series(forcing, surface_t, basin, physics, air, t_start=None, paired=True):
    """The model's surface-renewal transfer velocity of O2 (m/s) at every step, rebuilt from the forcing and the saved
    surface temperature: lake.run takes the bulk fluxes at the temperature the step starts from and the infrared,
    the expansion and the Schmidt number at the temperature after mixing, which the step ends with (nothing after the
    exchange changes it). paired: the bulk exchange at the GCM's lowest level, as the fixed model pairs it when the
    forcing gives that level; a pre-fix run took the 10 m wind and coefficient."""
    heat0 = float(L.light_shape(basin, physics)[0])
    sw = np.asarray(forcing.sw_down, float)
    lw = np.asarray(forcing.lw_down, float)
    ta, qa = np.asarray(forcing.air_k, float), np.asarray(forcing.q_air, float)
    u10, pa = np.asarray(forcing.wind10, float), np.asarray(forcing.pressure, float)
    level = paired and getattr(forcing, 'wind_level', None) is not None
    uz = np.asarray(forcing.wind_level, float) if level else u10
    zl = np.asarray(forcing.level_height, float) if level and forcing.level_height is not None else None
    ts = np.asarray(surface_t, float)
    g = physics.gravity
    k = np.empty(ts.size)
    u_stars = np.empty(ts.size)
    for n in range(ts.size):
        t0 = ts[n - 1] if n else (ts[0] if t_start is None else t_start)
        if zl is None:
            sens, lat_heat, _, tau = L.bulk_fluxes(t0, ta[n], qa[n], uz[n], pa[n], air)
        else:
            sens, lat_heat, _, tau = L.bulk_fluxes(t0, ta[n], qa[n], uz[n], pa[n], air, zl[n])
        u_star = math.sqrt(tau / L.RHO0)
        u_stars[n] = u_star
        t1 = float(ts[n])
        sw_water = sw[n] * (1.0 - physics.albedo)
        cooling = sens + lat_heat - L.WATER_EMISSIVITY * (lw[n] - L.STEFAN_BOLTZMANN * (t1 + 273.15) ** 4) \
            - sw_water * heat0
        buoyancy = g * L._expansion(t1) * max(cooling, 0.0) / L.RHO_CP
        k[n] = L.transfer_velocity('renewal', u10[n], u_star, buoyancy, t1)
    return np.stack([k, u_stars])


def wind_dissipation_series(forcing, basin, physics, air, mixed_depth, open_water):
    """The model's dissipation below the mixed layer each step under wind mixing (lake.run's own: the share of the
    wind's energy flux at 10 m over the mass of water below the mixed layer), where the run did not save it."""
    from research.studies.lake_methane.tracer import mixed_layers
    k_full, _ = mixed_layers(basin, mixed_depth)
    cum_v = np.cumsum(basin.volume)
    a0 = basin.area[0]
    u10, pa, ta = (np.asarray(x, float) for x in (forcing.wind10, forcing.pressure, forcing.air_k))
    p10 = L.wind_energy_flux(u10, pa / (air.gas_constant * ta))
    below = np.where(k_full < basin.n, L.RHO0 * (cum_v[-1] - cum_v[k_full - 1]) / a0, 0.0)
    eps = np.where(below > 0.0, physics.wind_share * p10 / np.maximum(below, 1e-300), 0.0)
    return np.where(np.asarray(open_water, bool), eps, 0.0)


def cached_transfer(key, build):
    """The transfer series cached in research/runs/lake_methane/transfer under a key of the saved file and the code;
    a cache-only run (lake_oxygen's cache_only()) builds none and stops at the first one not cached."""
    folder = RUNS / 'transfer'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f'{key}.npy'
    if path.exists():
        return np.load(path)
    if R.cache_only():
        R.uncached(path)
    k = build()
    np.save(path, k)
    return k


def code_key(st, *parts):
    """The transfer series' key: this module's, the lake model's and lake_oxygen's runner's code whole, what the
    series reads from the setting's inputs (its values_fingerprint()) and the parts, which give the saved files by their
    contents, so that a producer renaming its files (a re-key by hard links) leaves the key standing."""
    h = hashlib.sha256()
    for p in (Path(__file__), Path(L.__file__), Path(R.__file__)):
        h.update(sha256(p).encode())
    h.update(st.values_fingerprint().encode())
    h.update(json.dumps(parts, sort_keys=True, default=str).encode())
    return h.hexdigest()[:24]


# --- lake_anoxia's columns ----------------------------------------------------------------------------------------

def anoxia_meta(folder=ANOXIA_RUNS):
    meta = json.loads((Path(folder) / 'profiles.json').read_text())
    return meta


def anoxia_column(st, name, meta, traced=None, folder=ANOXIA_RUNS):
    """One of lake_anoxia's saved runs, from its folder (its own, or the development copy) with its profiles.json
    (meta). traced: the scenario names to keep (default the four central columns)."""
    case = meta['cases'][name]
    path = Path(folder) / f'{name}.npz'
    digest = sha256(path)
    if digest != case['npz_sha256']:
        raise SystemExit(f'{path} does not match its profiles.json')
    version = model_version(meta['producer']['files'].get('biosphere/ecology/lake.py'))
    z = np.load(path)
    names_all = list(meta['scenarios'])
    scen_all = scenarios_for(names_all, version)
    traced = traced or [f'{t}/{CENTRAL}' for t in R.TROPHIC]
    idx = [names_all.index(x) for x in traced]
    coupled = R.COUPLED['deep']
    physics = physics_for(st, coupled, version, names_all, eta=case.get('eta', 'central'))
    basin = L.Basin(z['z'], z['area'])
    site = R.site_info(st, BANDS.index(case['band']))
    _, forcing, _ = R.forcing_for(st, case['target'], site, 'lock', None, R.STEPS_PER_MEAN_DEEP)
    passes = case.get('passes', 1)
    if passes > 1:
        forcing = tile(forcing, passes)
    surface_t = z['surface_t'].astype(float)
    key = code_key(st, 'anoxia', name, digest, version)
    k, u_rebuilt = cached_transfer(key, lambda: transfer_velocity_series(forcing, surface_t, basin, physics, st.air,
                                                                         paired=version == 'current'))
    t = z['t']
    oxygen_all = z['oxygen']
    c = names_all.index(coupled)
    if 'solute' in z.files:
        solute, saved = z['solute'], True
    else:
        solute, saved = rebuilt_solute(st, t.astype(float), oxygen_all[:, :, c].astype(float), st.salinity,
                                       st.o2_fraction), False
    oxygen = np.ascontiguousarray(oxygen_all[:, :, idx])
    del oxygen_all
    debt_all = z['debt']
    mixed_depth = z['mixed_depth'].astype(float)
    snap_step = z['step'].astype(np.int64)
    if 'surface_o2' in z.files:
        ml_o2 = np.asarray(z['surface_o2'][:, idx], float)
    else:
        ml_o2 = mixed_layer_oxygen(basin, mixed_depth, z['mean_o2'][:, idx].astype(float), snap_step, oxygen)
    pool = np.ascontiguousarray(z['pool'][:, :, idx]) if 'pool' in z.files else None
    pressure = z['pressure'].astype(float)
    if passes > 1 and pressure.size * passes == surface_t.size:
        pressure = np.tile(pressure, passes)            # saved once, for the forcing each of the passes repeats
    return Column(name=name, basin=basin, physics=physics, scenarios=[scen_all[i] for i in idx], names=traced,
                  dt=float(z['dt']), snap_step=snap_step, t=t, oxygen=oxygen,
                  debt=np.ascontiguousarray(debt_all[:, :, idx]), age=z['age'], solute=solute, solute_saved=saved,
                  mixed_depth=mixed_depth, surface_t=surface_t, pressure=pressure, k_o2=k,
                  open_water=np.ones(surface_t.size, bool), ml_o2=ml_o2,
                  coupled_debt=np.ascontiguousarray(debt_all[:, :, c]), pool=pool, version=version,
                  dissipation=dissipation_of(z, forcing, basin, physics, st.air, mixed_depth,
                                             np.ones(surface_t.size, bool)),
                  u_star=np.asarray(z['u_star'], float) if 'u_star' in z.files else None, u_star_rebuilt=u_rebuilt,
                  meta=dict(source='lake_anoxia', file=str(path), sha256=digest, target=case['target'],
                            band=case['band'], kind=case['kind'], version=version))


def mixed_layer_oxygen(basin, mixed_depth, mean_o2, snap_step, oxygen):
    """The mixed layer's oxygen after each step where only the whole lake's mean was saved: that mean less what the
    water below the mixed layer holds in the last saved column, over the mixed layer's volume (exact when the whole
    lake is mixed)."""
    from research.studies.lake_methane.tracer import mixed_layers
    k_full, _ = mixed_layers(basin, mixed_depth)
    v = basin.volume
    cum_v = np.cumsum(v)
    total = cum_v[-1]
    s_of = np.maximum(np.searchsorted(snap_step, np.arange(mixed_depth.size), side='left') - 1, 0)
    below = np.cumsum((v[:, None, None] * oxygen.transpose(1, 0, 2))[::-1], axis=0)[::-1]  # (n, S, J) from k down
    below = np.concatenate([below, np.zeros((1,) + below.shape[1:])], axis=0)
    kk = np.maximum(k_full, 1)
    held_below = below[kk, s_of]                                                            # (N, J)
    out = (mean_o2 * total - held_below) / cum_v[kk - 1][:, None]
    return np.clip(out, 0.0, None)


def tile(forcing, k):
    return L.Forcing(dt=forcing.dt, sw_down=np.tile(forcing.sw_down, k), lw_down=np.tile(forcing.lw_down, k),
                     air_k=np.tile(forcing.air_k, k), q_air=np.tile(forcing.q_air, k),
                     wind10=np.tile(forcing.wind10, k), pressure=np.tile(forcing.pressure, k),
                     o2_fraction=forcing.o2_fraction,
                     **({} if getattr(forcing, 'wind_level', None) is None else
                        dict(wind_level=np.tile(forcing.wind_level, k),
                             level_height=np.tile(forcing.level_height, k))))


# --- the kinds' and the seas' columns -----------------------------------------------------------------------------

def dissipation_of(z, forcing, basin, physics, air, mixed_depth, open_water):
    """The dissipation each step: as the run saved it, or rebuilt where it mixed by the wind and did not save it."""
    if 'dissipation' in z.files:
        return np.asarray(z['dissipation'], float)
    if getattr(physics, 'deep_mixing', 'dissipation') == 'wind':
        return wind_dissipation_series(forcing, basin, physics, air, mixed_depth, open_water)
    return None


def saved_column(st, name, entry, folder, site, version, names_all, physics, target='resort', traced=None):
    """A case saved pass by pass (the lake kinds' columns.py and the seas' format), from its index entry (case, key,
    files with their sha256): the passes joined in time, the traced columns (default the case's own coupled column),
    the mixed layer's oxygen from the surface oxygen saved at every step, the coupled column's debt for the density,
    the labile pool, the friction velocity and the dissipation where the run saved them."""
    case = entry['case']
    scen_all = scenarios_for(names_all, version)
    coupled = names_all[physics.solute_column]
    traced = traced or [coupled]
    idx = [names_all.index(x) for x in traced]
    c = physics.solute_column
    files = sorted(entry['files'])
    keys = ('step', 't', 'oxygen', 'debt', 'coupled_debt', 'age', 'solute', 'mixed_depth', 'surface_t', 'pressure',
            'ice', 'surface_o2', 'pool', 'dissipation', 'u_star')
    parts = {k: [] for k in keys}
    offset, z0, digests = 0, None, {}
    for fname in files:
        f = Path(folder) / fname
        d = sha256(f)
        if d != entry['files'][fname]:
            raise SystemExit(f'{f} does not match its index')
        digests[fname] = d
        z = np.load(f)
        z0 = z0 or dict(z=z['z'], area=z['area'], dt=float(z['dt']))
        parts['step'].append(z['step'].astype(np.int64) + offset)
        for key in ('t', 'age', 'solute'):
            parts[key].append(z[key])
        debt = z['debt']
        parts['oxygen'].append(np.ascontiguousarray(z['oxygen'][:, :, idx]))
        parts['debt'].append(np.ascontiguousarray(debt[:, :, idx]))
        parts['coupled_debt'].append(np.ascontiguousarray(debt[:, :, c]))
        if 'pool' in z.files:
            parts['pool'].append(np.ascontiguousarray(z['pool'][:, :, idx]))
        for key in ('mixed_depth', 'surface_t', 'pressure', 'ice'):
            parts[key].append(np.asarray(z[key], float))
        for key in ('dissipation', 'u_star'):
            if key in z.files:
                parts[key].append(np.asarray(z[key], float))
        parts['surface_o2'].append(np.asarray(z['surface_o2'][:, idx], float))
        offset += z['mixed_depth'].size
        z.close()
        del debt
    joined = {k: (np.concatenate(v) if v else None) for k, v in parts.items()}
    basin = L.Basin(z0['z'], z0['area'])
    _, forcing, _ = R.forcing_for(st, target, site, 'lock', None, case['steps_per_mean'])
    if len(files) > 1:
        forcing = tile(forcing, len(files))
    key = code_key(st, 'saved', name, [digests[f] for f in files], version)      # the contents in their passes' order
    k, u_rebuilt = cached_transfer(key, lambda: transfer_velocity_series(forcing, joined['surface_t'], basin, physics,
                                                                         st.air, paired=version == 'current'))
    open_water = joined['ice'] <= 0.0
    eps = joined['dissipation']
    if eps is None and getattr(physics, 'deep_mixing', 'dissipation') == 'wind':
        eps = wind_dissipation_series(forcing, basin, physics, st.air, joined['mixed_depth'], open_water)
    return Column(name=name, basin=basin, physics=physics, scenarios=[scen_all[i] for i in idx], names=traced,
                  dt=z0['dt'], snap_step=joined['step'], t=joined['t'], oxygen=joined['oxygen'], debt=joined['debt'],
                  age=joined['age'], solute=joined['solute'], solute_saved=True,
                  mixed_depth=joined['mixed_depth'], surface_t=joined['surface_t'], pressure=joined['pressure'],
                  k_o2=k, open_water=open_water, ml_o2=joined['surface_o2'],
                  coupled_debt=joined['coupled_debt'], pool=joined['pool'], version=version, dissipation=eps,
                  u_star=joined['u_star'], u_star_rebuilt=u_rebuilt,
                  meta=dict(source=str(folder), files=digests, case=case, version=version))


def kinds_index(folder=KINDS_RUNS):
    return json.loads((Path(folder) / 'ladder.json').read_text())


def kinds_column(st301, name, index, folder=KINDS_RUNS, traced=None):
    """A lake kinds case at 301 K from the kinds' index: its forcing from the kinds' Setting301 (the resort moved to
    301 K), its scenarios in the order its summary records, its physics as the kinds' columns.py builds it (the case's
    own trophic state), its model version from the index's producer."""
    entry = index['cases'][name]
    case = entry['case']
    version = model_version(index['producer']['files'].get('biosphere/ecology/lake.py'))
    summary = json.loads((Path(folder) / 'cases' / entry['summary']).read_text())
    names_all = summary['case']['scenarios']
    site = R.site_info(st301, case['band'])
    physics = physics_for(st301, case['coupled'], version, names_all, mixing=case.get('mixing'),
                          eta=case.get('eta', 'central'), dissipation=case.get('dissipation', 'central'),
                          solute=case.get('solute', 'co2'), cap=case.get('cap'))
    return saved_column(st301, name, entry, Path(folder) / 'cases', site, version, names_all, physics, traced=traced)


# --- the seas ------------------------------------------------------------------------------------------------------

def seas_index(folder=SEAS_RUNS):
    """The seas' columns index (sea_oxygen's columns.json) and its sha256 as read."""
    path = Path(folder) / 'columns.json'
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def seas_column(st_sea, name, index, folder=SEAS_RUNS, traced=None):
    """A sea column at 301 K: sea_oxygen's setting (the seas' water cells in Setting301's resort), its site and its
    physics (sea_oxygen.columns.physics, the case's deep mixing, wind share, cap and solutes), the scenarios in the
    order its summary records, its model version from the index's lake.py."""
    from research.studies.sea_oxygen import columns as SC
    entry = index['cases'][name]
    case = entry['case']
    version = model_version(index['lake_py']['sha256'])
    summary = json.loads((Path(folder) / 'cases' / entry['summary']).read_text())
    names_all = summary['case']['scenarios']
    site = st_sea.site(case['cell'], case['basin'])
    physics = SC.physics(st_sea, case)
    if version != 'current':
        raise SystemExit(f'{name}: the seas were run on lake.py {index["lake_py"]["sha256"]}, not the current model')
    traced = traced or [f'{t}/{CENTRAL}' for t in R.TROPHIC]
    col = saved_column(st_sea, name, entry, Path(folder) / 'cases', site, version, names_all, physics, traced=traced)
    form = entry.get('diffusivity') or case.get('form') or 'osborn'
    # the index names the lake number's mixing as a diffusivity; the tracer takes it from the physics' deep mixing
    # (lake.lake_number_diffusivity), as the model did, over Osborn's form
    col.form = 'osborn' if form == 'lake_number' and getattr(physics, 'deep_mixing', None) == 'lake_number' else form
    col.meta['form'] = col.form
    return col


# --- the mixing forms' columns (lake_oxygen's Shih-Bouffard sensitivity) -------------------------------------------

MIXING_FORMS_RUNS = ROOT / 'research' / 'runs' / 'lake_oxygen' / 'mixing_forms'


def mixing_forms_index(folder=MIXING_FORMS_RUNS):
    path = Path(folder) / 'mixing_forms_index.json'
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def mixing_forms_column(name, index, folder=MIXING_FORMS_RUNS, traced=None):
    """A column of lake_oxygen's mixing-forms study, as its column() builds it: the resort at the case's ground
    temperature (the kinds' Setting301; the lock and the hot world are lake_oxygen's whatever it), the 0-15 lake's
    site, lake_oxygen's physics with the case's coupled state, and the diffusivity form the case ran under
    (Shih and Bouffard's three regimes, or Osborn's rule with the wind's dissipation bounded, which its saved
    dissipation carries). The 301 K rows trace their own coupled column, the three targets' rows the four central
    columns their mesotrophic physics carries."""
    from research.studies.lake_kinds import setting as KS
    entry = index['cases'][name]
    case = entry['case']
    version = model_version(index['lake_py_sha256'])
    if version != 'current':
        raise SystemExit(f'{name}: the mixing forms were run on lake.py {index["lake_py_sha256"]}')
    st = KS.Setting301(case['ground_k'])
    summary = json.loads((Path(folder) / 'cases' / entry['summary']).read_text())
    names_all = summary['case']['scenarios']
    site = R.site_info(st, 0)
    physics = physics_for(st, case['coupled'], version, names_all)
    one = case['coupled'] != R.COUPLED['deep'] or str(case.get('ground_k')) == '301.0'
    traced = traced or ([case['coupled']] if one else [f'{t}/{CENTRAL}' for t in R.TROPHIC])
    case_files = dict(entry, files=entry.get('files') or summary['files'],
                      case=dict(case, steps_per_mean=R.STEPS_PER_MEAN_DEEP))
    col = saved_column(st, name, case_files, Path(folder) / 'cases', site, version, names_all, physics,
                       target=case['target'], traced=traced)
    col.form = case.get('form', 'osborn')
    col.meta['form'] = col.form
    return col
