"""Why the Near-side Sea box's nights stay shallow, and why the GCM's simplified Betts-Miller scheme rains there anyway:
the box's clouds' entrainment measured from its own 3-D output, the height at which a plume mixing at that rate stops in
the night's and the afternoon's air, and the scheme's decision traced on the GCM's columns and its own snapshots.

    climate/gcm/.venv/bin/python -m climate.crm.box_entrainment    # ../results/crm/box_entrainment_near_side_sea.json

CM1. Every 3-D output of box_near_side_sea_day3's third lunar day (3-hourly) is sampled level by level, as Siebesma
and Cuijpers (1995, J. Atmos. Sci. 52) sampled their shallow cumulus: the cloud, points holding at least 0.01 g/kg of
cloud water and ice (box_onset's threshold); the cloudy updraft, cloud rising; and the core, cloudy updraft lighter than
its level's mean in density potential temperature, which counts every water species' weight. Each sample's environment
is every other point of its level. The outputs are pooled by quarter of the lunar day (night, after midnight; morning;
afternoon; evening), and the night-test boxes' 6-hourly outputs from hour 12 (their shallow cloud covers 29-34% of the
columns by then) are pooled by grid. Two variables that the cloud's own condensation leaves unchanged, the total
non-precipitating water qt = qv + qc + qi and the moist static energy h = cp T + g z + Lv qv - Lf qi, give the bulk
plume's fractional entrainment from each pooled profile, epsilon = -(d phi_c/dz)/(phi_c - phi_e), its detrainment as
epsilon less the mass flux's fractional growth. Rain forming in a sample takes qt away and freezing adds to h, so below
the freezing level the night's warm, barely raining cloud keeps both, and in the afternoon's deep cells neither holds
exactly. Resampling each quarter's outputs gives the range of the bulk rate, the mass-flux-weighted mean of epsilon
over the cloud below 16 km.

A plume then rises through each quarter's mean sounding: from the lowest level with the mean of the lowest 0.5 km (as
box_onset's parcel), undiluted to its condensation level, then mixing a share 1 - exp(-epsilon dz) of the surrounding
air's h and water between levels; its condensate stays in it as weight (reversible), or falls out as it forms. Its cap
is the first height above its level of free convection at which its buoyancy vanishes; a plume buoyant nowhere above
its condensation level has none. The caps are set against the columns' cloud tops, the top of each column's lowest
unbroken cloud below 40 km, which leaves out the thin ice that forms in place near 50 km.

GCM. On the rh 0.7 branch's site cells over years 10-14, two records of the scheme: its 3-day mean columns, on which
sbm.column (the tested Python copy) traces the decision at once, and the convection's heating at the last step of each
mean (dt_conv, one step, 18 degrees of hour angle after the mean's centre), whose highest heated layer is the top the
running scheme chose and whose column integral is the rain it formed. The scheme's columns are compared with CM1's
night (and with the Kuo GCM that forced the box) on PlaSim's levels, and the plume above, lifted from PlaSim's lowest
level through each mean column laid on CM1's levels and dropping its condensate as the scheme's parcel does, shows
where CM1's entrainment would stop it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from climate.crm import box_budget as bb
from climate.crm import box_onset as bo
from climate.crm import cm1_output as co
from climate.crm import ring_analysis as ra

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'crm'
OUT = RESULTS / 'box_entrainment_near_side_sea.json'
SCHEMA = 'terluna.climate.crm-box-entrainment/1'
DAY = 86400.0
LV, LS = 2.501e6, 2.834e6                  # CM1's latent heats of vaporization and sublimation (constants.F)
LF = LS - LV
EPSV = co.RD / co.RV
SOURCE = dict(case='box_near_side_sea_day3', first_day=59.0, last_day=88.5)
NIGHT_TESTS = {'6 km': 'box_near_side_sea_night_6km', '3 km': 'box_near_side_sea_night_3km',
               '2 km': 'box_near_side_sea_night_2km'}
NIGHT_TEST_FROM_H = 12.0
PERIODS = (('night', -180.0, -90.0), ('morning', -90.0, 0.0), ('afternoon', 0.0, 90.0), ('evening', 90.0, 180.0))
SAMPLES = ('cloud', 'updraft', 'core')
VARIABLES = ('qt', 'h')
MIN_POINTS = 50                            # sampled points at a level, pooled, for its profile to count
LOWER_M, RISING_M = 16000.0, 40000.0       # the lower cloud layer; cloud rising from below stays under 40 km
RESAMPLES, SEED = 300, 20261025
MIXED_M = 500.0
EPS_FIXED = (0.0, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3)    # per metre
HEIGHTS_M = (1000.0, 4000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 12000.0, 14000.0, 16000.0, 20000.0, 30000.0,
             40000.0, 50000.0)
LAYERS_M = (('0-6 km', 0.0, 6000.0), ('6-10 km', 6000.0, 10000.0), ('10-16 km', 10000.0, 16000.0),
            ('above 16 km', 16000.0, 1.0e6))
GCM_RUN, GCM_YEARS, GCM_RH = 'A28R_301k_cold_sbm_rh07', (10, 14), 0.7
KUO_RUN, KUO_YEARS = 'A28R_301k_cold_moist', (0, 5)       # the GCM that forced the boxes
FORMED_MM_DAY = 0.05                       # rain formed at a step that counts as the scheme raining (sbm_offline's)
SNAPSHOT_LAG_DEG = 18.0
LOW_TOP_LEVEL = 7                          # PlaSim's levels top first: 7 and 8 stand near 10 and 4 km at the sea
EFOLD_M = co.RD * 270.0 / 1.6242188765608128   # metres per e-folding of pressure at 270 K and lunar gravity, as box_onset


# ---------------------------------------------------------------------------------------------------- CM1's clouds

def thermo(d: dict, zh, g: float) -> dict:
    """A snapshot's fields (levels, columns) as the sampling takes them: temperature, pressure, vapour, cloud
    condensate (liquid and ice), total non-precipitating water qt, density potential temperature (every species'
    weight), moist static energy h (J/kg, frozen for the cloud ice), density and vertical wind."""
    th, p = np.asarray(d['th'], float), np.asarray(d['prs'], float)
    t = th * (p / co.P00) ** (co.RD / co.CP)
    qv, qc, qi = (np.asarray(d[k], float) for k in ('qv', 'qc', 'qi'))
    qp = sum(np.asarray(d[k], float) for k in ('qr', 'qs', 'qg') if k in d)
    qt = qv + qc + qi
    weight = (1.0 + qv / EPSV) / (1.0 + qt + qp)
    z = np.asarray(zh, float)[:, None]
    return dict(t=t, p=p, qv=qv, qcl=qc + qi, qt=qt, thrho=th * weight, h=co.CP * t + g * z + LV * qv - LF * qi,
                rho=p / (co.RD * t * weight), w=np.asarray(d['winterp'], float), th=th)


def masks(x: dict) -> dict:
    """The samples: cloud (at least ra.CLOUD_KG_KG of liquid and ice), cloudy updraft (cloud rising) and core (cloudy
    updraft lighter than its level's mean in density potential temperature)."""
    cloud = x['qcl'] >= ra.CLOUD_KG_KG
    updraft = cloud & (x['w'] > 0.0)
    core = updraft & (x['thrho'] > x['thrho'].mean(axis=1, keepdims=True))
    return dict(cloud=cloud, updraft=updraft, core=core)


def snapshot_sums(x: dict, m: dict) -> dict:
    """Level-by-level sums over one snapshot: for each sample its points and the sums of qt, h, w, rho w (its mass
    flux times the columns) and buoyancy (density potential temperature less the level's mean), and its environment's
    points, qt and h; the domain's sums of temperature, vapour, pressure, potential temperature, condensate and
    density, with the column count."""
    out = dict(columns=np.full(x['t'].shape[0], float(x['t'].shape[1])))
    excess = x['thrho'] - x['thrho'].mean(axis=1, keepdims=True)
    for name, mask in m.items():
        env = ~mask
        out[f'{name}_n'] = mask.sum(axis=1).astype(float)
        out[f'{name}_env_n'] = env.sum(axis=1).astype(float)
        for v in VARIABLES:
            out[f'{name}_{v}'] = (x[v] * mask).sum(axis=1)
            out[f'{name}_env_{v}'] = (x[v] * env).sum(axis=1)
        out[f'{name}_w'] = (x['w'] * mask).sum(axis=1)
        out[f'{name}_rhow'] = (x['rho'] * x['w'] * mask).sum(axis=1)
        out[f'{name}_b'] = (excess * mask).sum(axis=1)
    for k in ('t', 'qv', 'p', 'th', 'qcl', 'rho'):
        out[f'mean_{k}'] = x[k].sum(axis=1)
    return out


def cloud_tops(qcl, zh, limit_m: float = RISING_M) -> np.ndarray:
    """The top (m) of each column's lowest unbroken cloud whose base lies below limit_m, for the columns that have
    one: from its lowest cloudy level up to the last level before the first clear one."""
    cloudy = np.asarray(qcl) >= ra.CLOUD_KG_KG
    zh = np.asarray(zh, float)
    has = cloudy.any(axis=0)
    base = np.argmax(cloudy, axis=0)
    levels = np.arange(zh.size)[:, None]
    gap = (levels > base[None, :]) & ~cloudy
    first_gap = np.where(gap.any(axis=0), np.argmax(gap, axis=0), zh.size)
    keep = has & (zh[base] < limit_m)
    return zh[first_gap[keep] - 1]


def tops_summary(tops: list, columns: int) -> dict:
    """Pooled cloud tops (m) over a set of snapshots with `columns` columns in all."""
    t = np.concatenate(tops) if tops else np.zeros(0)
    if not t.size:
        return dict(columns=columns, cloudy_share=0.0)
    return dict(columns=columns, cloudy_share=float(t.size / columns), mean_m=float(t.mean()),
                median_m=float(np.median(t)), p90_m=float(np.percentile(t, 90)), p99_m=float(np.percentile(t, 99)),
                highest_m=float(t.max()), below_16km_share=float((t < LOWER_M).mean()))


def pooled(rows: list) -> dict:
    """The sums of several snapshots' snapshot_sums."""
    return {k: np.sum([r[k] for r in rows], axis=0) for k in rows[0]}


def profiles(s: dict, sample: str) -> dict:
    """A sample's pooled profiles: its points, area share, mean qt and h and its environment's, mean vertical wind,
    buoyancy (K) and mass flux (kg/m2/s), NaN where it has fewer than MIN_POINTS points."""
    n = s[f'{sample}_n']
    ok = n >= MIN_POINTS
    safe = np.where(ok, n, 1.0)
    out = dict(ok=ok, points=n, area=n / s['columns'], w=np.where(ok, s[f'{sample}_w'] / safe, np.nan),
               buoyancy_k=np.where(ok, s[f'{sample}_b'] / safe, np.nan), mass_flux=s[f'{sample}_rhow'] / s['columns'])
    for v in VARIABLES:
        out[v] = np.where(ok, s[f'{sample}_{v}'] / safe, np.nan)
        out[f'{v}_env'] = s[f'{sample}_env_{v}'] / np.maximum(s[f'{sample}_env_n'], 1.0)
    return out


def fractional(zh, phi_c, phi_e, ok) -> np.ndarray:
    """The bulk plume's fractional entrainment (per metre), -(d phi_c/dz)/(phi_c - phi_e), from a sample's profile
    and its environment's: the derivative centred where both neighbours count, one-sided where one does, NaN
    elsewhere."""
    zh, phi_c, phi_e, ok = (np.asarray(a) for a in (zh, phi_c, phi_e, ok))
    eps = np.full(zh.size, np.nan)
    for k in range(zh.size):
        if not ok[k]:
            continue
        below, above = k > 0 and ok[k - 1], k + 1 < zh.size and ok[k + 1]
        if below and above:
            d = (phi_c[k + 1] - phi_c[k - 1]) / (zh[k + 1] - zh[k - 1])
        elif above:
            d = (phi_c[k + 1] - phi_c[k]) / (zh[k + 1] - zh[k])
        elif below:
            d = (phi_c[k] - phi_c[k - 1]) / (zh[k] - zh[k - 1])
        else:
            continue
        diff = phi_c[k] - phi_e[k]
        if diff != 0.0:
            eps[k] = -d / diff
    return eps


def detrainment(zh, eps, mass_flux, ok) -> np.ndarray:
    """The fractional detrainment (per metre): entrainment less the mass flux's fractional growth with height."""
    m = np.where(np.asarray(ok) & (np.asarray(mass_flux) > 0.0), mass_flux, np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        growth = np.gradient(np.log(m), np.asarray(zh, float))
    return np.asarray(eps) - growth


def bulk(zh, eps, mass_flux, top_m: float = LOWER_M):
    """The mass-flux-weighted mean of epsilon over the levels below top_m where it is diagnosed and the flux rises;
    None where there are none."""
    use = np.isfinite(eps) & (np.asarray(zh) < top_m) & (np.asarray(mass_flux) > 0.0)
    if not use.any():
        return None
    return float(np.sum(eps[use] * mass_flux[use]) / np.sum(mass_flux[use]))


def entrainment(s: dict, zh) -> dict:
    """Each sample's pooled profiles with its entrainment from qt and from h and its detrainment, and the bulk rates
    below 16 km."""
    out = {}
    for sample in SAMPLES:
        p = profiles(s, sample)
        eps = {v: fractional(zh, p[v], p[f'{v}_env'], p['ok'] & (np.asarray(zh) < RISING_M)) for v in VARIABLES}
        out[sample] = dict(p, eps_qt=eps['qt'], eps_h=eps['h'],
                           delta_qt=detrainment(zh, eps['qt'], p['mass_flux'], p['ok']),
                           bulk_qt=bulk(zh, eps['qt'], p['mass_flux']), bulk_h=bulk(zh, eps['h'], p['mass_flux']))
    return out


def resampled(rows: list, zh, resamples: int = RESAMPLES, seed: int = SEED) -> dict:
    """The 5th and 95th percentiles of each sample's bulk rates over resamplings of the snapshots with replacement."""
    rng = np.random.default_rng(seed)
    got = {f'{s}_{v}': [] for s in SAMPLES for v in VARIABLES}
    for _ in range(resamples):
        e = entrainment(pooled([rows[i] for i in rng.integers(0, len(rows), len(rows))]), zh)
        for s in SAMPLES:
            for v in VARIABLES:
                if e[s][f'bulk_{v}'] is not None:
                    got[f'{s}_{v}'].append(e[s][f'bulk_{v}'])
    return {k: [float(np.percentile(v, 5)), float(np.percentile(v, 95))] if len(v) >= 0.9 * resamples else None
            for k, v in got.items()}


def sounding(s: dict, zh, zw) -> dict:
    """The pooled domain-mean sounding: temperature, vapour, pressure, potential temperature and condensate by level,
    with each level's mass per area from the mean density."""
    n = s['columns']
    mean = {k: s[f'mean_{k}'] / n for k in ('t', 'qv', 'p', 'th', 'qcl', 'rho')}
    mean['mass'] = mean['rho'] * np.diff(np.asarray(zw, float))
    return mean


# ------------------------------------------------------------------------------------------------- the entraining plume

def plume(z, p, t, q, t0, q0, eps, g, fallout: bool = False):
    """A plume from the lowest level with temperature t0 and vapour q0 through a sounding (levels the ground up):
    undiluted until it first saturates, then mixing a share 1 - exp(-eps dz) of the surrounding air's moist static
    energy and water between levels, eps per metre (a number, or a value for each level). Its condensate stays in it as weight,
    or with fallout falls out as it forms. Returns its buoyancy (m/s2) against the surrounding air's virtual
    temperature and whether it is saturated, level by level."""
    z, p, t, q = (np.asarray(a, float) for a in (z, p, t, q))
    rate = np.broadcast_to(np.asarray(eps, float), z.shape)
    h, qt = co.CP * t0 + g * z[0] + LV * q0, q0
    he = co.CP * t + g * z + LV * q
    b, wet = np.zeros(z.size), np.zeros(z.size, bool)
    for k in range(z.size):
        if k and wet[:k].any() and rate[k] > 0.0:
            mix = 1.0 - np.exp(-rate[k] * (z[k] - z[k - 1]))
            h, qt = (1.0 - mix) * h + mix * he[k], (1.0 - mix) * qt + mix * q[k]
        tk, qk, wet[k] = bo.saturate(h, qt, z[k], p[k], g)
        if fallout:
            qt = min(qt, qk)
        b[k] = g * (tk * (1.0 + qk / EPSV) / (1.0 + qt) - bo.virtual(t[k], q[k])) / bo.virtual(t[k], q[k])
    return b, wet


def cap(b, wet, z) -> dict:
    """A plume's condensation level, level of free convection (the first buoyant level at or above condensation and
    below 40 km, where cloud rising from below stays), cap
    (the height above that where its buoyancy first vanishes, between levels by the zero crossing) and highest buoyant
    level (m), and the CAPE of its first buoyant layer (J/kg); a plume buoyant nowhere above condensation has no lfc
    or cap, and its least negative buoyancy above condensation is given. cin_j_kg is the negative buoyancy it crosses
    below its level of free convection."""
    b, wet, z = np.asarray(b, float), np.asarray(wet, bool), np.asarray(z, float)
    if not wet.any():
        return dict(lcl_m=None, lfc_m=None, cap_m=None, el_m=None, cape_j_kg=0.0, cin_j_kg=None,
                    best_buoyancy_m_s2=None)
    lcl = int(np.argmax(wet))
    above = np.arange(z.size) >= lcl
    up = above & (b > 0.0) & (z < RISING_M)
    best = float(b[above].max())
    if not up.any():
        return dict(lcl_m=float(z[lcl]), lfc_m=None, cap_m=None, el_m=None, cape_j_kg=0.0, cin_j_kg=None,
                    best_buoyancy_m_s2=best)
    lfc = int(np.argmax(up))
    dz = np.gradient(z)
    cin = float(-np.sum(np.minimum(b[:lfc], 0.0) * dz[:lfc]))
    stop = next((k for k in range(lfc + 1, z.size) if b[k] <= 0.0), None)
    if stop is None:
        cap_m, last = float(z[-1]), z.size - 1
    else:
        cap_m, last = float(z[stop - 1] + (z[stop] - z[stop - 1]) * b[stop - 1] / (b[stop - 1] - b[stop])), stop - 1
    return dict(lcl_m=float(z[lcl]), lfc_m=float(z[lfc]), cap_m=cap_m, el_m=float(z[np.flatnonzero(b > 0.0).max()]),
                cape_j_kg=float(np.sum(b[lfc:last + 1] * dz[lfc:last + 1])), cin_j_kg=cin, best_buoyancy_m_s2=best)


def mixed_source(snd: dict, zw, depth: float = MIXED_M):
    """The lowest-level temperature and vapour of the lowest `depth` metres mixed by mass, as box_onset's parcel."""
    share = bo.layer_share(zw, depth)
    return bo.source(snd['p'], snd['t'], snd['qv'], snd['mass'] * share, depth)


def plume_cap(snd: dict, zh, zw, eps, g, fallout: bool = False) -> dict:
    t0, q0 = mixed_source(snd, zw)
    return cap(*plume(zh, snd['p'], snd['t'], snd['qv'], t0, q0, eps, g, fallout), zh)


def undiluted(snd: dict, zh, zw, g) -> dict:
    """The undiluted 0.5 km parcel of a mean sounding, reversible and pseudo-adiabatic: its buoyancy (m/s2) at the
    reporting heights and the CAPE of the pseudo-adiabatic parcel by layer (J/kg), its buoyant part only."""
    t0, q0 = mixed_source(snd, zw)
    zh = np.asarray(zh, float)
    out = {}
    for label, fallout in (('reversible', False), ('pseudo_adiabatic', True)):
        b, _ = plume(zh, snd['p'], snd['t'], snd['qv'], t0, q0, 0.0, g, fallout)
        out[f'{label}_buoyancy_m_s2'] = {f'{z / 1000:g} km': float(np.interp(z, zh, b)) for z in HEIGHTS_M}
        if fallout:
            dz = np.gradient(zh)
            out['pseudo_adiabatic_cape_by_layer_j_kg'] = {
                lab: float(np.sum(np.where((zh >= a) & (zh < c) & (b > 0.0), b * dz, 0.0)))
                for lab, a, c in (('below 10 km', 0.0, 10000.0), ('10-16 km', 10000.0, 16000.0),
                                  ('16-40 km', 16000.0, 40000.0), ('above 40 km', 40000.0, 1.0e6))}
    return out


# ------------------------------------------------------------------------------------------------- reading the boxes

def case_info(name: str) -> dict:
    case = ra.RUNS / name
    record = json.loads((case / 'case.json').read_text())
    geo = ra.case_geometry(case)
    return dict(case=case, record=record, g=record['build']['gravity_m_s2'], zh=geo['zh'], zw=geo['zw'],
                tap=record['configuration']['output_s'], dx_m=geo['dx'])


def period_of(hour_deg: float):
    return next((name for name, a, b in PERIODS if a <= hour_deg < b), None)


def read_source(name: str = SOURCE['case'], span=(SOURCE['first_day'], SOURCE['last_day'])) -> dict:
    """The source box's 3-D outputs over its third lunar day, sampled and grouped by quarter of the lunar day, with
    their cloud tops."""
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    info = case_info(name)
    cfg = info['record']['configuration']
    groups = {p[0]: dict(rows=[], tops=[], outputs=[], columns=0) for p in PERIODS}
    for n in sorted(int(f.name[8:14]) for f in info['case'].glob('cm1out_t*_s.dat')):
        t = (n - 1) * info['tap']
        if not span[0] * DAY + 1.0 < t <= span[1] * DAY + 1.0:
            continue
        x = thermo(ra.read_snapshot(info['case'], n), info['zh'], info['g'])
        group = groups[period_of(local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], t,
                                                  solar_day_s()))]
        group['rows'].append(snapshot_sums(x, masks(x)))
        group['tops'].append(cloud_tops(x['qcl'], info['zh']))
        group['outputs'].append(n)
        group['columns'] += x['t'].shape[1]
    return dict(info=info, groups=groups)


def read_night_test(name: str, from_h: float = NIGHT_TEST_FROM_H) -> dict:
    info = case_info(name)
    rows, tops, outputs, columns = [], [], [], 0
    for n in sorted(int(f.name[8:14]) for f in info['case'].glob('cm1out_t*_s.dat')):
        if (n - 1) * info['tap'] / 3600.0 < from_h - 1e-6:
            continue
        x = thermo(ra.read_snapshot(info['case'], n), info['zh'], info['g'])
        rows.append(snapshot_sums(x, masks(x)))
        tops.append(cloud_tops(x['qcl'], info['zh']))
        outputs.append(n)
        columns += x['t'].shape[1]
    return dict(info=info, rows=rows, tops=tops, outputs=outputs, columns=columns)


def forcing_by_period(name: str = SOURCE['case'], span=(SOURCE['first_day'], SOURCE['last_day'])) -> dict:
    """The prescribed forcing (the day-night table at the hour and the mean vertical wind on the box's own profile,
    box_onset.forcing) and the radiation as the air takes it, by quarter: layer totals in W/m2 and mm/day, and the
    heating (K/day) and moistening (g/kg/day) at the reporting heights."""
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    info = case_info(name)
    cfg = info['record']['configuration']
    zh, zw = info['zh'], info['zw']
    dz = np.diff(zw)
    layout = ra.grads_layout(info['case'] / 'cm1out_diag_s.ctl')[1]
    eq = co.eqtset(info['case'])
    acc = {p[0]: [] for p in PERIODS}
    for path in sorted(info['case'].glob('cm1out_diag_*_s.dat')):
        d = bb.read_diag(info['case'], path, layout)
        if not span[0] * DAY + 1.0 < d['mtime'] <= span[1] * DAY + 1.0:
            continue
        tt, qq = bo.forcing(info['case'], d, zh, zw)
        rad = co.applied_rate(np.asarray(d['swten']) + np.asarray(d['lwten']), eq)
        hour = local_hour_angle(cfg['start_hour_angle_deg'], cfg['site']['lon_deg'], d['mtime'], solar_day_s())
        acc[period_of(hour)].append((tt, qq, rad, d['rho'] * dz, d['hpbl'], d['prate'] * DAY))
    out = {}
    for name_, rows in acc.items():
        tt, qq, rad, mass = (np.mean([r[i] for r in rows], axis=0) for i in range(4))
        layer = lambda f, unit: {lab: float(np.sum((f * mass)[(zh >= a) & (zh < b)]) * unit) for lab, a, b in LAYERS_M}
        at = lambda f, unit: {f'{z / 1000:g} km': float(np.interp(z, zh, f) * unit) for z in HEIGHTS_M}
        out[name_] = dict(outputs=len(rows), boundary_layer_m=float(np.mean([r[4] for r in rows])),
                          rain_mm_day=float(np.mean([r[5] for r in rows])),
                          forcing_heating_w_m2=layer(tt, co.CP), forcing_moistening_mm_day=layer(qq, DAY),
                          radiation_w_m2=layer(rad, co.CP), forcing_k_day=at(tt, DAY), forcing_g_kg_day=at(qq, DAY * 1e3),
                          radiation_k_day=at(rad, DAY))
    return out


def environment(snd: dict, zh, g: float) -> dict:
    """A mean sounding's air at the reporting heights: temperature, potential temperature, vapour, relative humidity
    over liquid (Bolton's, as the parcels take it), moist static energy h and its saturated value h* (kJ/kg), and the
    layers' stability: the potential temperature's lapse (K/km), and the lowest h above the boundary layer."""
    zh = np.asarray(zh, float)
    qs = bo.saturation(snd['t'], snd['p'])
    h = co.CP * snd['t'] + g * zh + LV * snd['qv']
    hs = co.CP * snd['t'] + g * zh + LV * qs
    at = lambda f, scale=1.0: {f'{z / 1000:g} km': float(np.interp(z, zh, f) * scale) for z in HEIGHTS_M}
    lapse = {}
    for a, b in ((6000.0, 8000.0), (8000.0, 10000.0), (10000.0, 16000.0), (16000.0, 30000.0)):
        th = np.interp([a, b], zh, snd['th'])
        lapse[f'{a / 1000:g}-{b / 1000:g} km'] = float((th[1] - th[0]) / (b - a) * 1000.0)
    free = (zh > 6000.0) & (zh < 60000.0)
    k = int(np.flatnonzero(free)[np.argmin(h[free])])
    return dict(temperature_k=at(snd['t']), potential_temperature_k=at(snd['th']), vapour_g_kg=at(snd['qv'], 1e3),
                relative_humidity=at(snd['qv'] / qs), h_kj_kg=at(h, 1e-3), h_saturated_kj_kg=at(hs, 1e-3),
                theta_lapse_k_km=lapse, h_lowest_above_6km=dict(z_m=float(zh[k]), h_kj_kg=float(h[k] * 1e-3)),
                h_at_ground_kj_kg=float(h[0] * 1e-3))


def level_table(e: dict, zh, keep) -> list:
    """The levels of a sample's profile that count, as rows for the record."""
    rows = []
    for k in np.flatnonzero(keep):
        r = lambda a, scale=1.0: None if not np.isfinite(a[k]) else float(a[k] * scale)
        rows.append(dict(z_m=float(zh[k]), points=int(e['points'][k]), area=float(e['area'][k]), w_m_s=r(e['w']),
                         buoyancy_k=r(e['buoyancy_k']), mass_flux_kg_m2_s=float(e['mass_flux'][k]),
                         qt_minus_env_g_kg=r(e['qt'] - e['qt_env'], 1e3), h_minus_env_j_kg=r(e['h'] - e['h_env']),
                         eps_qt_per_m=r(e['eps_qt']), eps_h_per_m=r(e['eps_h']), delta_qt_per_m=r(e['delta_qt'])))
    return rows


def sample_record(e: dict, zh, ranges=None) -> dict:
    out = {}
    for sample in SAMPLES:
        s = e[sample]
        keep = s['ok'] & (np.asarray(zh) < RISING_M)
        out[sample] = dict(bulk_eps_qt_per_m=s['bulk_qt'], bulk_eps_h_per_m=s['bulk_h'],
                           levels=level_table(s, zh, keep))
        if ranges is not None:
            out[sample]['bulk_eps_qt_5_95'] = ranges[f'{sample}_qt']
            out[sample]['bulk_eps_h_5_95'] = ranges[f'{sample}_h']
    return out


def own_rates(e: dict) -> dict:
    """The rates a plume takes from one pooled set: the bulk core and cloudy-updraft rates from qt and h."""
    return {f'{s} {v}': e[s][f'bulk_{v}'] for s in ('core', 'updraft') for v in VARIABLES
            if e[s][f'bulk_{v}'] is not None}


def analyse_cm1() -> dict:
    src = read_source()
    info = src['info']
    zh, zw, g = info['zh'], info['zw'], info['g']
    forcing = forcing_by_period()
    periods, soundings, rates = {}, {}, {}
    for name, grp in src['groups'].items():
        s = pooled(grp['rows'])
        e = entrainment(s, zh)
        soundings[name] = sounding(s, zh, zw)
        rates[name] = own_rates(e)
        periods[name] = dict(snapshots=len(grp['rows']), outputs=[grp['outputs'][0], grp['outputs'][-1]],
                             samples=sample_record(e, zh, resampled(grp['rows'], zh)),
                             cloud_tops=tops_summary(grp['tops'], grp['columns']),
                             environment=environment(soundings[name], zh, g), forcing=forcing[name],
                             undiluted_parcel=undiluted(soundings[name], zh, zw, g))
    plumes = {}
    for env in ('night', 'morning', 'afternoon', 'evening'):
        snd = soundings[env]
        choices = {'undiluted': 0.0}
        for source in ('night', 'afternoon', 'evening'):
            choices.update({f'{source} {k}': v for k, v in rates[source].items()})
        choices.update({f'fixed {x:g}/m': x for x in EPS_FIXED[1:]})
        plumes[env] = {label: dict(eps_per_m=float(x), reversible=plume_cap(snd, zh, zw, x, g),
                                   pseudo_adiabatic=plume_cap(snd, zh, zw, x, g, fallout=True))
                       for label, x in choices.items()}
    resolution = {}
    for grid, name in NIGHT_TESTS.items():
        nt = read_night_test(name)
        s = pooled(nt['rows'])
        e = entrainment(s, nt['info']['zh'])
        snd = sounding(s, nt['info']['zh'], nt['info']['zw'])
        rates_t = own_rates(e)
        second = [r for r, n in zip(nt['rows'], nt['outputs']) if (n - 1) * nt['info']['tap'] / 3600.0 > 24.0 + 1e-6]
        e2 = entrainment(pooled(second), nt['info']['zh'])
        resolution[grid] = dict(
            case=name, dx_m=nt['info']['dx_m'], snapshots=len(nt['rows']), outputs=nt['outputs'],
            columns_per_snapshot=int(nt['columns'] // len(nt['rows'])),
            samples=sample_record(e, nt['info']['zh']),
            second_day=dict(snapshots=len(second), bulk={f'{s} {v}': e2[s][f'bulk_{v}'] for s in SAMPLES
                                                         for v in VARIABLES}),
            cloud_tops=tops_summary(nt['tops'], nt['columns']),
            plumes={label: plume_cap(snd, nt['info']['zh'], nt['info']['zw'], x, nt['info']['g'])
                    for label, x in dict(undiluted=0.0, **rates_t).items()})
    return dict(source=dict(case=SOURCE['case'], span_days=[SOURCE['first_day'], SOURCE['last_day']],
                            dx_m=info['dx_m'], columns_per_snapshot=int(info['record']['grid']['nx'] *
                                                                         info['record']['grid']['ny']),
                            output_s=info['tap'], gravity_m_s2=g),
                periods=periods, plumes=plumes, resolution=resolution,
                eps_per_efolding_at_270k=float(EFOLD_M))


# ------------------------------------------------------------------------------------------------- the GCM's scheme

def quarter(hour_deg) -> str:
    """The quarter of the lunar day (sbm_offline's names) of an hour angle."""
    from climate.gcm import sbm_offline as so
    h = (float(hour_deg) + 180.0) % 360.0 - 180.0
    return next(name for name, a, b in so.QUARTERS if a <= h < b)


def gcm_records(run: str, years, box: str = 'box_near_side_sea') -> dict:
    """The GCM's 3-day means in the cells of a box's site, with the convection's heating at each mean's last step:
    temperature (K), specific humidity, surface pressure (Pa), cos-latitude weight, the site's hour angle at the
    mean's centre (degrees), year, convective rain at the ground (mm/day) and dt_conv (K/s), levels top first."""
    import netCDF4
    from climate.crm import gcm_site as gs
    record = json.loads((ra.RUNS / box / 'case.json').read_text())
    cells = [tuple(c) for c in record['gcm_site']['cells']]
    lon_site = record['configuration']['site']['lon_deg']
    out = {k: [] for k in ('t', 'q', 'ps', 'weight', 'hour', 'year', 'prc', 'dtc')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(HERE.parent / 'gcm' / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            idx = [(int(np.argmin(abs(lat - a))), int(np.argmin(abs(lon - b)))) for a, b in cells]
            f = {k: np.asarray(d[k][:], float) for k in ('ta', 'hus', 'ps', 'prc', 'dt_conv')}
            hour = (lon_site - gs.subsolar_longitude(np.asarray(d['czen'][:], float), lat, lon) + 180.0) % 360.0 - 180.0
        for n in range(f['ta'].shape[0]):
            for j, i in idx:
                out['t'].append(f['ta'][n, :, j, i])
                out['q'].append(f['hus'][n, :, j, i])
                out['ps'].append(f['ps'][n, j, i] * 100.0)
                out['weight'].append(np.cos(np.radians(lat[j])))
                out['hour'].append(hour[n])
                out['year'].append(year)
                out['prc'].append(f['prc'][n, j, i] * 1000.0 * DAY)
                out['dtc'].append(f['dt_conv'][n, :, j, i])
    return dict({k: np.asarray(v) for k, v in out.items()}, sigma=sigma, edges=gs.plasim_edges(sigma), cells=cells)


def formed_rain(dtc, q, ps, dsigma, c) -> np.ndarray:
    """The rain (mm/day) the scheme formed at a step: the column integral of cp times its heating over L, which its
    moist-enthalpy correction makes the rain's latent heat; the shallow scheme's mixing adds nothing to it."""
    m = np.asarray(ps)[..., None] * np.asarray(dsigma) / c['ga']
    return (c['acpd'] * (1.0 + c['ADV'] * np.asarray(q)) * np.asarray(dtc) * m).sum(axis=-1) / c['ALV'] * DAY


def heated_top(dtc) -> int:
    """The highest layer (index, top first) the convection heats or cools at a step, -1 for none."""
    hot = np.flatnonzero(np.abs(np.asarray(dtc)) > 1e-12)
    return int(hot.min()) if hot.size else -1


def vapour_above_reference(t, q, sigma, dsigma, ps, tp, top, rh, c) -> np.ndarray:
    """Each convecting layer's vapour above the scheme's reference, rh of saturation at the parcel's temperature
    (kg/m2; the scheme rains their sum over the relaxation time when it is positive), top first, nil above the top."""
    out = np.zeros(len(t))
    if top is None:
        return out
    from climate.gcm import sbm
    for k in range(top, len(t)):
        out[k] = (q[k] - rh * sbm.saturation(tp[k], ps * sigma[k], c)) * ps * dsigma[k] / c['ga']
    return out


def fine_column(t, q, ps, sigma, zh, gravity):
    """A GCM column laid on CM1's levels for a plume: its full levels' heights from their virtual temperatures, then
    temperature, vapour (mixing ratio) and log pressure linear in height between them, from the lowest full level to the
    highest; returns the heights, pressure, temperature and mixing ratio there."""
    from climate.crm import gcm_site as gs
    z = gs.level_heights(np.asarray(t, float), np.asarray(q, float), np.asarray(sigma, float), gravity)[::-1]
    p = (ps * np.asarray(sigma, float))[::-1]
    use = (zh >= z[0]) & (zh <= z[-1])
    zz = np.concatenate([[z[0]], np.asarray(zh, float)[use & (zh > z[0])]])
    r = np.asarray(q, float)[::-1] / (1.0 - np.asarray(q, float)[::-1])
    return (zz, np.exp(np.interp(zz, z, np.log(p))), np.interp(zz, z, np.asarray(t, float)[::-1]),
            np.interp(zz, z, r))


def analyse_gcm(cm1_rates: dict, cm1_night: dict) -> dict:
    from climate.crm import gcm_site as gs
    from climate.gcm import sbm, sbm_offline as so
    c = so.run_constants(GCM_RUN)
    g = gcm_records(GCM_RUN, GCM_YEARS)
    sigma, edges = g['sigma'], g['edges']
    dsigma = np.diff(edges)
    z_mean = gs.level_heights(g['t'], g['q'], sigma, c['ga'])
    names = [q[0] for q in so.QUARTERS]
    q_mean = np.array([quarter(h) for h in g['hour']])
    q_snap = np.array([quarter(h + SNAPSHOT_LAG_DEG) for h in g['hour']])
    w = g['weight']
    # the running scheme at each mean's last step
    formed = formed_rain(g['dtc'], g['q'], g['ps'], dsigma, c)
    raining = formed > FORMED_MM_DAY
    tops = np.array([heated_top(x) if r else -1 for x, r in zip(g['dtc'], raining)])
    level_z = {int(k): float(np.average(z_mean[:, k], weights=w)) for k in range(sigma.size)}
    steps = {}
    for name in names:
        s = q_snap == name
        rain = float(np.sum(formed[s] * w[s]) / np.sum(w[s]))
        by_top = {}
        for k in range(sigma.size):
            m = s & raining & (tops == k)
            if m.any():
                by_top[f'level {k}'] = dict(sigma=float(sigma[k]), height_m=level_z[k], steps=int(m.sum()),
                                            rain_share=float(np.sum(formed[m] * w[m]) / np.sum(formed[s] * w[s])),
                                            mean_rate_mm_day=float(formed[m].mean()))
        low = s & raining & (tops >= LOW_TOP_LEVEL)
        def share(i):
            total = np.sum((formed * w)[i])
            return float(np.sum((formed * w)[i] * (tops[i] >= LOW_TOP_LEVEL) * raining[i]) / total) if total else np.nan
        idx = np.flatnonzero(s)
        rng = np.random.default_rng(SEED)
        spread = [share(rng.choice(idx, idx.size)) for _ in range(RESAMPLES)] if rain > 0 else None
        steps[name] = dict(steps=int(s.sum()), raining=int((s & raining).sum()), formed_mm_day=rain,
                           by_top=by_top, topped_at_or_below_10km_rain_share=float(
                               np.sum(formed[low] * w[low]) / np.sum(formed[s] * w[s])) if rain > 0 else None,
                           topped_at_or_below_10km_5_95=[float(np.nanpercentile(spread, 5)),
                                                         float(np.nanpercentile(spread, 95))] if rain > 0 else None)
    # the decision at once on each mean column
    decisions = []
    for t, q, ps in zip(g['t'], g['q'], g['ps']):
        o = sbm.column(t, q, sigma, dsigma, ps, tau=so.TAU_S, rh=GCM_RH, eps=0.0, step=so.STEP_S, akap=so.KAPPA, c=c)
        excess = vapour_above_reference(t, q, sigma, dsigma, ps, o['reference_t'], o['top'], GCM_RH, c)
        decisions.append(dict(top=o['top'], cin=o['cin'], deep=o['top'] is not None and not o['shallow'],
                              rain=float(o['rain'].sum()) * DAY, excess=excess))
    deep = np.array([d['deep'] for d in decisions])
    buoyant = np.array([d['top'] is not None for d in decisions])
    top_index = np.array([d['top'] if d['top'] is not None else -1 for d in decisions])
    offline = {}
    for name in names:
        s = q_mean == name
        sd = s & deep
        ex = np.array([d['excess'] for d, k in zip(decisions, sd) if k])
        offline[name] = dict(
            means=int(s.sum()), buoyant=int((s & buoyant).sum()), deep_and_raining=int(sd.sum()),
            shallow=int((s & buoyant & ~deep).sum()),
            buoyant_tops={f'level {k}': int(np.sum(s & buoyant & (top_index == k))) for k in range(sigma.size)},
            deep_tops={f'level {k}': int(np.sum(sd & (top_index == k))) for k in range(sigma.size)},
            parcel_inhibition_j_kg=float(np.median([d['cin'] for d, k in zip(decisions, s & buoyant) if k]))
            if (s & buoyant).any() else None,
            vapour_above_reference_kg_m2=dict(
                total_median=float(np.median(ex.sum(axis=1))) if sd.any() else None,
                lowest_two_layers_median=float(np.median(ex[:, -2:].sum(axis=1))) if sd.any() else None,
                by_level_mean=[float(v) for v in ex.mean(axis=0)] if sd.any() else None),
            rain_at_once_mm_day=float(np.sum(np.array([d['rain'] for d in decisions])[s] * w[s]) / np.sum(w[s])))
    # the columns against CM1's night, on PlaSim's levels
    sat = lambda t, p: np.array([sbm.saturation(a, b, c) for a, b in zip(t, p)])

    def level_means(mask):
        t = np.average(g['t'][mask], axis=0, weights=w[mask])
        q = np.average(g['q'][mask], axis=0, weights=w[mask])
        ps = float(np.average(g['ps'][mask], weights=w[mask]))
        return dict(t_k=t.tolist(), q_g_kg=(q * 1e3).tolist(), rh=(q / sat(t, ps * sigma)).tolist(), ps_pa=ps,
                    height_m=gs.level_heights(t, q, sigma, c['ga']).tolist())
    night = q_mean == 'midnight-sunrise'
    compare = dict(sigma=sigma.tolist(), gcm_night_mean=level_means(night),
                   gcm_night_deep_means=level_means(night & deep) if (night & deep).any() else None,
                   gcm_afternoon_mean=level_means(q_mean == 'noon-sunset'))
    for mapping in ('level', 'layer'):
        t, q, _, _ = so.on_levels(cm1_night, mapping, sigma, edges, cm1_night['zh'], cm1_night['zw'])
        compare[f'cm1_night_{mapping}'] = dict(t_k=t.tolist(), q_g_kg=(q * 1e3).tolist(),
                                              rh=(q / sat(t, cm1_night['ps'] * sigma)).tolist(),
                                              ps_pa=float(cm1_night['ps']))
    kuo = gcm_records(KUO_RUN, KUO_YEARS)
    kuo_night = np.array([quarter(h) for h in kuo['hour']]) == 'midnight-sunrise'
    kt = np.average(kuo['t'][kuo_night], axis=0, weights=kuo['weight'][kuo_night])
    kq = np.average(kuo['q'][kuo_night], axis=0, weights=kuo['weight'][kuo_night])
    kps = float(np.average(kuo['ps'][kuo_night], weights=kuo['weight'][kuo_night]))
    compare['kuo_forcing_gcm_night_mean'] = dict(run=KUO_RUN, years=list(KUO_YEARS), t_k=kt.tolist(),
                                                 q_g_kg=(kq * 1e3).tolist(), rh=(kq / sat(kt, kps * sigma)).tolist(),
                                                 ps_pa=kps)
    # an entraining plume from PlaSim's lowest level through each mean column on CM1's levels
    zh = cm1_night['zh']
    rates = dict(undiluted=0.0, **cm1_rates)
    plumes = {}
    for name in names:
        rows = {}
        for label, sel in (('deep and raining at once', (q_mean == name) & deep),
                           ('buoyant', (q_mean == name) & buoyant)):
            idx = np.flatnonzero(sel)
            if not idx.size:
                continue
            caps = {label_: [] for label_ in rates}
            scheme_top = []
            for i in idx:
                zz, pp, tt, rr = fine_column(g['t'][i], g['q'][i], g['ps'][i], sigma, zh, c['ga'])
                for label_, x in rates.items():
                    cp_ = cap(*plume(zz, pp, tt, rr, tt[0], rr[0], x, c['ga'], fallout=True), zz)
                    caps[label_].append(cp_['cap_m'] if cp_['cap_m'] is not None else np.nan)
                scheme_top.append(float(z_mean[i, decisions[i]['top']]))
            st = np.array(scheme_top)

            def summary(x, v):
                return dict(eps_per_m=float(x), eps_per_efolding=float(x * EFOLD_M),
                            cap_median_m=float(np.nanmedian(v)) if np.isfinite(v).any() else None,
                            never_buoyant_share=float(np.mean(~np.isfinite(v))),
                            stops_below_scheme_top_share=float(np.mean(~np.isfinite(v) | (v < st))))
            rows[label] = dict(columns=int(idx.size), scheme_top_median_m=float(np.median(st)),
                               plumes={label_: summary(x, np.array(caps[label_], float)) for label_, x in rates.items()})
        plumes[name] = rows
    return dict(run=GCM_RUN, rh=GCM_RH, years=list(GCM_YEARS), cells=[list(x) for x in g['cells']],
                cell_means=int(len(g['hour'])), level_heights_m=level_z, sigma=sigma.tolist(),
                running_scheme_steps=steps, decision_on_means=offline, against_cm1=compare, plumes=plumes)


def cm1_night_column() -> dict:
    """CM1's night: the source box's domain-mean column over its post-midnight outputs (sbm_offline's column form),
    with its levels."""
    from climate.gcm import sbm_offline as so
    b = so.cm1_columns(SOURCE['case'], (SOURCE['first_day'], SOURCE['last_day']))
    cols = [x for x in b['columns'] if quarter(x['hour']) == 'midnight-sunrise']
    return dict(t=np.mean([x['t'] for x in cols], axis=0), q=np.mean([x['q'] for x in cols], axis=0),
                p=np.mean([x['p'] for x in cols], axis=0), ps=float(np.mean([x['ps'] for x in cols])),
                zh=b['zh'], zw=b['zw'], outputs=len(cols))


def analyse() -> dict:
    cm1 = analyse_cm1()
    night = cm1['periods']['night']['samples']
    rates = {f'CM1 night {s} {v}': night[s][f'bulk_eps_{v}_per_m'] for s in ('core', 'updraft') for v in VARIABLES
             if night[s][f'bulk_eps_{v}_per_m'] is not None}
    after = cm1['periods']['afternoon']['samples']
    rates.update({f'CM1 afternoon {s} qt': after[s]['bulk_eps_qt_per_m'] for s in ('core', 'updraft')
                  if after[s]['bulk_eps_qt_per_m'] is not None})
    column = cm1_night_column()
    gcm = analyse_gcm(rates, column)
    gcm['against_cm1']['cm1_night_outputs'] = column['outputs']
    return dict(schema=SCHEMA, cm1=cm1, gcm=gcm, evidence=EVIDENCE, reading_rule=READING_RULE)


EVIDENCE = ('Diagnosis from existing output, no new runs. CM1: one realization of the sea box\'s third lunar day on its '
            '6 km grid (3-D output every 3 hours, 59 snapshots a quarter) and the three night-test boxes at 6, 3 and '
            '2 km (6-hourly, 7 snapshots each from hour 12), every one with CM1\'s boundary-layer scheme and 111 '
            'levels about 0.5-0.7 km apart in the night\'s cloud layer; the entrainment is the resolved clouds\' bulk '
            'rate on those grids, a measurement of the model, and the resampling range treats 3-hourly snapshots as '
            'independent. The plume is a bulk-plume calculation on mean soundings. GCM: 3-day means and one-step '
            'snapshots of the convection\'s heating at the end of each mean, five years of the rh 0.7 branch at its six '
            'sea cells; the scheme\'s decision on the means is the tested Python copy run at once, and the running '
            'scheme\'s tops come from its heating at single steps.')
READING_RULE = ('cm1.periods: by quarter of the lunar day at the site (hour angle from noon: night -180 to -90, morning '
                '-90 to 0, afternoon 0 to 90, evening 90 to 180). samples: cloud (cloud water and ice at least 0.01 '
                'g/kg), updraft (cloud rising), core (cloudy updraft lighter than its level\'s mean in density '
                'potential temperature), each against every other point of its level; levels with at least 50 pooled '
                'points below 40 km: area share, w, buoyancy (K), mass flux (kg/m2/s), qt and h less the environment\'s, '
                'entrainment from qt and h and detrainment from qt (per metre); bulk: the mass-flux-weighted mean '
                'below 16 km, with the 5-95% range over 300 resamplings of the snapshots. cloud_tops: the top of each '
                'column\'s lowest unbroken cloud below 40 km, pooled. plumes: by environment (the quarter\'s mean '
                'sounding), a plume from the lowest 0.5 km mixed, undiluted to its condensation level, then mixing '
                'eps_per_m; reversible keeps its condensate as weight, pseudo_adiabatic drops it; cap is the first '
                'height above its free-convection level where its buoyancy vanishes, None where it is buoyant nowhere '
                'above condensation. eps per e-folding of pressure is eps per metre times eps_per_efolding_at_270k (m). '
                'gcm.running_scheme_steps: the convection\'s heating at the last step of each mean, by quarter of that '
                'step; raining where its column integral forms more than 0.05 mm/day; by_top, the highest heated '
                'layer (index top first, its mean height). decision_on_means: sbm.column at rh 0.7 without mixing on '
                'each 3-day mean, by quarter of the mean\'s centre; vapour_above_reference by level top first (kg/m2). '
                'against_cm1: means on PlaSim\'s full levels, top first, rh by PlaSim\'s saturation; CM1\'s night '
                'interpolated in log pressure (level) and as layer means by mass (layer). plumes: the plume above, '
                'pseudo-adiabatic, from PlaSim\'s lowest level through each mean column laid on CM1\'s levels, against '
                'the scheme\'s top on the same column.')


def compact(x, digits: int = 5):
    """The record with every float to `digits` significant figures, NaN as None."""
    if isinstance(x, dict):
        return {k: compact(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [compact(v, digits) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(x) else float(f'{float(x):.{digits}g}')
    if isinstance(x, np.integer):
        return int(x)
    return x


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args(argv)
    out = compact(analyse())
    out['producer'] = dict(domain='climate', files={
        f'{p.parent.name}/{p.name}': hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        for p in (Path(__file__).resolve(), HERE / 'box_onset.py', HERE.parent / 'gcm' / 'sbm.py',
                  HERE.parent / 'gcm' / 'sbm_offline.py')})
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(HERE.parents[1])}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
