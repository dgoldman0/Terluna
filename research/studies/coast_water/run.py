"""How well the fresh seas and lakes water their coasts and islands at the resort's 301 K: the groundwater the shores
hold within plants' reach under the deserts, and the dew and fog of the long night against the rain.

    climate/gcm/.venv/bin/python -m research.studies.coast_water.extract   # once, about a minute; needs netCDF4
    python -m research.studies.coast_water.shores                          # once, about twenty seconds
    python -m research.studies.coast_water.run                             # about seven minutes, 2.4 GB

It binds the GCM extract and the 16 px/deg shores (both ignored, on the research drive), climate_targets' product for
the resort's basis, geography's water table with its 2 px/deg grids, the four local coasts at 100 m, the CM1 rings'
and the Jura box's products and the Jura Island screen, and writes results/coast_water.json. No model is run beyond
the steady belts and the surface balances themselves.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

from research.studies.climate_targets import extract as CT
from research.studies.coast_water import belts as B
from research.studies.coast_water import dew as D
from research.studies.coast_water import extract as EX
from research.studies.coast_water import shores as S
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'coast_water.json'
SCHEMA = 'terluna.research.coast-water/1'
CODE = ('research/studies/coast_water/run.py', 'research/studies/coast_water/belts.py',
        'research/studies/coast_water/dew.py', 'research/studies/coast_water/shores.py',
        'research/studies/coast_water/extract.py', 'research/studies/climate_targets/extract.py')
TARGETS = 'research/studies/climate_targets/results/climate_targets.json'
TABLE = 'geography/results/groundwater_table.json'
RINGS = {'45n_corrected': 'climate/results/crm/ring_ring_45n_corrected_twilight.json',
         '45s_corrected': 'climate/results/crm/ring_ring_45s_corrected_twilight.json',
         '45n_wet': 'climate/results/crm/ring_ring_45n_lsw.json',
         '45s_wet': 'climate/results/crm/ring_ring_45s_lsw.json',
         'equator': 'climate/results/crm/ring_ring_equator.json'}
JURA_BOX = 'climate/results/crm/jura_box_jura.json'
JURA_SHORES = 'research/studies/jura_island/results/shores.json'
INPUTS = {TARGETS: 'terluna.research.climate-targets/5', TABLE: 'terluna.geography.groundwater-table/1',
          JURA_BOX: 'terluna.climate.crm-jura/1', JURA_SHORES: 'terluna.research.jura-shores/1',
          **{p: 'terluna.climate.crm-ring/1' for p in RINGS.values()}}
GCM = 'research/runs/coast_water/gcm_extract.npz'
GCM_SIDECAR = 'research/runs/coast_water/gcm_extract.json'
SHORES = 'research/runs/coast_water/shores_16ppd.npz'
TABLE_GRIDS = 'research/runs/resources/groundwater/groundwater_table_2ppd.npz'
COASTS = {name: f'research/runs/sea_appearance/terrain/{name}.npz'
          for name in ('w_procellarum', 'ingenii_coast', 's_nubium', 'smythii_headland')}
IGNORED = (GCM, GCM_SIDECAR, SHORES, TABLE_GRIDS) + tuple(COASTS.values())
TARGET_K = 301.0
DESERT_MM_DAY = CT.DESERT_MM_DAY
BANDS = CT.RAIN_BANDS
REACH_M = (15.0, 25.0)                 # height above the water table from which plants still draw (geography)
SANDY_FACTOR = 10.0                    # sandy ground: ten times the central transmissivity (geography's screen)
DAY_S = JULIAN_DAY
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
LUNAR_DAY_DAYS = SYNODIC_MONTH_DAYS
LUNAR_NIGHT_H = 0.5 * SYNODIC_MONTH_DAYS * 24.0
TEN_TO_TWO_M = 4.87 / math.log(67.8 * 10.0 - 5.42)    # FAO-56 equation 47, a wind at 10 m carried to 2 m
SHORE_CLASSES = ('sea_coast', 'island_coast', 'common_level_lake_shore', 'open_lake_shore', 'closed_lake_shore')
LAND_ALBEDO, WATER_ALBEDO = 0.20, 0.069          # the GCM's own surface albedos (constant in both runs)
LAND_EMISSIVITY, WATER_EMISSIVITY = 1.0, 0.98    # PlaSim's radmod.f90: zeps = dls + 0.98 (1 - dls)
NIGHT_GROUND_OFFSET_K = 0.3                      # the rings' night air at 2 m above their ground, design ground
FOG_LWC_KG_M3 = (1e-5, 5e-4)                     # radiation fog, 0.01-0.5 g/m3 (the food study's sources)
FOG_DEPOSITION_M_S = (5e-4, 0.021)               # lunar settling of a 10 um drop to Earth's short vegetation's least
FOG_DEPOSITION_HIGH_M_S = 0.08                   # Earth's short vegetation's greatest (Katata 2014), winds 1.1-9 m/s
CM1_RUNS = 'climate/crm/runs'
LOCAL_MARGIN_KM = S.NEAR_KM                      # the local coasts' edges left out of their check
EVIDENCE = (
    'A runnable calculation on GCM output and geography\'s products; no model is run beyond the steady belts and '
    'the surface balances. The resort at 301 K is interpolated between the design case and the hot world\'s test by '
    'the ground\'s temperature, as climate_targets builds it, each cell\'s rain and each per-cell result taken '
    'across the two runs, so it mixes their two physics; the rain is taken linearly and by a constant factor a '
    'kelvin, which the two runs cannot choose between. The belts are steady Dupuit flow from each shore\'s water '
    'with the plants\' uptake as the sink, at three transmissivities (stage 2\'s slow and central crust at the '
    'ground and ten times the central, sandy ground) and reaches of 15 and 25 m, on LOLA\'s 16 px/deg heights with '
    'the drainage product\'s lakes and the atlas\'s seas; their sub-pixel widths are a screen, which four coasts at '
    '100 m, none in a desert, put at 1.0-2.7 times the belts resolved there, and the 100 m belts converge on one of '
    'the four. The plants\' demand is FAO-56\'s Penman-Monteith '
    'reference surface on the GCM\'s own radiation, air, humidity and wind, the wind carried down from the lowest '
    'level by the CM1 rings\' ratios. The night\'s dew and fog are an energy balance on the GCM\'s 3-day means, '
    'checked against its 6.5-hour means and against Earth\'s measured dew; ExoPlaSim itself forms no dew, since its '
    'surface flux can only raise the lowest level\'s humidity (fluxmod.f90, line 664). The 15 km water table is '
    'geography\'s on the earlier design climate\'s recharge. T21 resolves no coast: each cell\'s air is all land or '
    'all water. The slice\'s cost is scaled from the CM1 runs\' own records; no slice is run.')
READING_RULE = (
    'Lengths in km, areas in km2, widths in m, heights in m above the water\'s level. Deserts are T21 cells whose '
    'rain at 301 K is under 0.5 mm a day, the rain taken linearly (rain_mask) or by a constant factor a kelvin '
    '(log_rain_mask), or the linear rain under the threshold that gives climate_targets\' interpolated desert share '
    '(share_matched); every_cell counts all the land in such cells at 16 px/deg, '
    'model_land only the cells the GCM runs as land. Shares are of that desert land. Belt cases are named by '
    'transmissivity (slow, central, sandy) and reach; uptake is the plants\' demand less the linear rain where it '
    'falls. within_reach_corrected divides each case\'s share by its own range of 16 px/deg over 100 m on the '
    'local coasts. '
    'Demand and dew are mm a day as time means over the whole lunar day (dew falls only at night), mm per lunar day '
    'for a lunar cycle of 29.53 days; a 301 K value is the per-cell interpolation of the two runs\' own values, '
    'weighted by land area. Islands are 8-connected land bodies other than the main land, assigned to the desert by '
    'their centroid\'s cell. Pairs and ranges are [low, high].')


def digest(path) -> str:
    h = hashlib.sha256()
    with open(ROOT / path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def load():
    data = {p: json.loads((ROOT / p).read_text()) for p in INPUTS}
    for p, schema in INPUTS.items():
        if data[p].get('schema') != schema:
            raise ValueError(f'{p}: schema {data[p].get("schema")} is not {schema}')
    sidecar = json.loads((ROOT / GCM_SIDECAR).read_text())
    if sidecar.get('schema') != EX.SCHEMA or sidecar['npz']['sha256'] != digest(GCM):
        raise ValueError('the GCM extract does not match its sidecar')
    return data, sidecar


# ------------------------------------------------------------------------------------------------- the 301 K basis
def basis(targets):
    lo, hi = targets['runs']['design']['ground_k'], targets['runs']['test']['ground_k']
    f = (TARGET_K - lo) / (hi - lo)
    row = next(r for r in targets['resort_band']['rows'] if abs(r['ground_k'] - TARGET_K) < 1e-9)
    if abs(row['share_of_the_way'] - f) > 1e-12:
        raise ValueError('the share of the way differs from climate_targets\'')
    return dict(target_k=TARGET_K, design=dict(run=targets['runs']['design']['run'], ground_k=lo),
                test=dict(run=targets['runs']['test']['run'], ground_k=hi), share_of_the_way=f,
                climate_targets=dict(land_rain_mm_day=row['rain']['land_rain_mm_day'],
                                     land_desert=row['rain']['land_desert'],
                                     model_land_desert=row['rain']['model_land_desert']),
                climate_targets_runs=dict(every_cell=[targets['runs'][k]['land_desert'] for k in ('design', 'test')],
                                          model_land=[targets['runs'][k]['model_land_desert']
                                                      for k in ('design', 'test')]))


def interpolate(a, b, f):
    return a + f * (b - a)


def t21(g):
    lat = g['lat']
    lsm = g['design_lsm']
    wl = g['land_fraction'] * g['cell_area_fraction']
    alat = np.abs(lat)[:, None] * np.ones((1, 64))
    coastal = lsm & CT.beside(~lsm)
    return dict(lat=lat, lsm=lsm, wl=wl, alat=alat, coastal=coastal, inland=lsm & ~coastal,
                water_with_land=~lsm & (g['land_fraction'] > 0), area=g['cell_area_fraction'])


def band_of(alat, lo, hi):
    return (alat >= lo) & ((alat < hi) if hi < 90 else (alat <= hi))


def share_matched(rain, weights, share):
    """The rain threshold (mm a day) under which the land's share first reaches the given one, and that share: the
    threshold lies midway between the last cell it takes in and the next."""
    order = np.argsort(rain.ravel(), kind='stable')
    r = rain.ravel()[order]
    c = np.cumsum(weights.ravel()[order]) / weights.sum()
    k = int(np.searchsorted(c, share))
    return float(0.5 * (r[k] + r[min(k + 1, r.size - 1)])), float(c[k])


FORMS = dict(additive='each cell\'s rain moved linearly in the ground\'s temperature',
             log_linear='each cell\'s rain moved by a constant factor for each kelvin')
SLOPE_K = (301.0, 302.0)              # where the forms' desert curves are compared, per kelvin


def rain_at(rd, rt, f, form):
    """Each cell's rain (mm a day) a share f of the way from the design case to the test."""
    return interpolate(rd, rt, f) if form == 'additive' else rd ** (1.0 - f) * rt ** f


def deserts(g, grid, f, base, shore_cells, targets_product):
    """The deserts at 301 K: each cell's rain interpolated in each form and the land under 0.5 mm a day, both counts;
    the mask whose threshold gives climate_targets' interpolated share; the curves from 295 to 307 K; and by band how
    far each form's mask stands from the interpolated shares."""
    rd, rt = g['design_pr_mean'] * DAY_S * 1000.0, g['test_pr_mean'] * DAY_S * 1000.0
    rain = {form: rain_at(rd, rt, f, form) for form in FORMS}
    wl, lsm, alat = grid['wl'], grid['lsm'], grid['alat']
    counts = dict(every_cell=wl, model_land=wl * lsm)
    targets = dict(every_cell=base['climate_targets']['land_desert'],
                   model_land=base['climate_targets']['model_land_desert'])
    masks, out = {}, dict(rain_mm_day=dict(design=rd, test=rt), forms=FORMS, curve=[])

    def entry(m, w, count, cut):
        land16 = shore_cells['land'] * (m if count == 'every_cell' else m & lsm)
        return dict(threshold_mm_day=cut, share=float((w * m).sum() / w.sum()), land_km2_16ppd=float(land16.sum()),
                    land_km2_boundary=float((w * m).sum() * 4.0 * math.pi * (MOON_RADIUS / 1000.0) ** 2),
                    bands={f'{lo}-{hi}': float((w * m * band_of(alat, lo, hi)).sum()
                                               / max((w * band_of(alat, lo, hi)).sum(), 1e-30)) for lo, hi in BANDS})

    for count, w in counts.items():
        for form, name in (('additive', 'rain_mask'), ('log_linear', 'log_rain_mask')):
            m = rain[form] < DESERT_MM_DAY
            masks[(name, count)] = m
            out.setdefault(name, {})[count] = entry(m, w, count, DESERT_MM_DAY)
        threshold, got = share_matched(rain['additive'], np.where(w > 0, w, 0.0), targets[count])
        m = rain['additive'] < threshold
        masks[('share_matched', count)] = m
        out.setdefault('share_matched', {})[count] = dict(entry(m, w, count, threshold),
                                                          climate_targets_share=targets[count])
    lo, hi = base['design']['ground_k'], base['test']['ground_k']
    for k in np.arange(295.0, 307.0, 0.5):
        ff = (min(max(float(k), lo), hi) - lo) / (hi - lo)
        row = dict(ground_k=float(k))
        for form in FORMS:
            r = rain_at(rd, rt, ff, form)
            row.update({f'{form}_{c}': float((w * (r < DESERT_MM_DAY)).sum() / w.sum()) for c, w in counts.items()})
        row['interpolated_shares_every_cell'] = (base['climate_targets_runs']['every_cell'][0]
                                                 + ff * (base['climate_targets_runs']['every_cell'][1]
                                                         - base['climate_targets_runs']['every_cell'][0]))
        out['curve'].append(row)
    by_k = {row['ground_k']: row for row in out['curve']}
    a, b = by_k[SLOPE_K[0]], by_k[SLOPE_K[1]]
    out['points_per_kelvin'] = dict(
        between_k=list(SLOPE_K),
        **{form: 100.0 * (b[f'{form}_every_cell'] - a[f'{form}_every_cell']) for form in FORMS},
        interpolated_shares=100.0 * (b['interpolated_shares_every_cell'] - a['interpolated_shares_every_cell']))
    # by band, each form's mask against climate_targets' interpolated band shares, in points of all the land
    runs = targets_product['runs']
    out['against_interpolated_shares'] = {}
    for form, name in (('additive', 'rain_mask'), ('log_linear', 'log_rain_mask')):
        rows = {}
        for lo_, hi_ in BANDS:
            key = f'{lo_}-{hi_}'
            dd = runs['design']['land']['every_cell']['bands'][key]['desert']
            dt = runs['test']['land']['every_cell']['bands'][key]['desert']
            share = (wl * band_of(alat, lo_, hi_)).sum() / wl.sum()
            mask = out[name]['every_cell']['bands'][key]
            band_rain = {k: runs[k]['land']['every_cell']['bands'][key]['rain_mm_day'] for k in ('design', 'test')}
            rows[key] = dict(land_share=float(share), interpolated=dd + f * (dt - dd), mask=mask,
                             points=float(100.0 * share * (dd + f * (dt - dd) - mask)), rain_mm_day=band_rain,
                             design_desert=dd)
        out['against_interpolated_shares'][form] = dict(bands=rows,
                                                        points=float(sum(v['points'] for v in rows.values())))
    return rain['additive'], masks, out


# ------------------------------------------------------------------------------------------------ the GCM's surface
GARDNER = dict(conductivity_m_s=1e-5, alpha_per_m=(2.0, 1.0))     # geography's sandy loam and its finer soil
ROOTS_M = (3.0, 5.0)


def capillary_reach(uptake_m_s):
    """Geography's screen of the reach at lunar gravity: the height above the water table from which Gardner's
    exponential soil still lifts the plants' uptake, ln(1 + eK/q)/(alpha e) with e the ratio of the gravities, for
    its sandy loam and its finer soil, and with 3-5 m of roots added."""
    from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
    e = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    rise = [math.log(1.0 + e * GARDNER['conductivity_m_s'] / uptake_m_s) / (a * e) for a in GARDNER['alpha_per_m']]
    return dict(uptake_mm_day=uptake_m_s * 1000.0 * DAY_S, capillary_rise_m=rise,
                with_roots_m=[rise[0] + ROOTS_M[0], rise[1] + ROOTS_M[1]])


def ring_ratios(rings):
    """The CM1 rings' wind at 10 m over that at 1 km, mean speeds by day (within 90 degrees of noon) and by night,
    over land and water: the two rings on the corrected design ground for the central value, all five for the range."""
    out = {}
    for key, ring in rings.items():
        b = ring['by_hour_angle']
        ha = np.abs(np.array(b['hour_angle_deg'], float))
        for part in ('land', 'water'):
            u10 = np.array(b[part]['wind_10m_m_s'], float)
            u1 = np.array(b[part]['wind_by_height_m_s'], float)[0]
            for name, m in (('day', ha <= 90.0), ('night', ha > 90.0)):
                out.setdefault(part, {}).setdefault(name, {})[key] = float(np.nanmean(u10[m]) / np.nanmean(u1[m]))
    return {part: {name: dict(by_ring=v, central=float(np.mean([v['45n_corrected'], v['45s_corrected']])),
                              range=[float(min(v.values())), float(max(v.values()))])
                   for name, v in parts.items()} for part, parts in out.items()}


def vector_speed_factor(g, window, grid):
    """How far the 3-day means' mean vector understates the mean of the 6.5-hour speeds within them, over land and
    water, by day and by night, from the design case's two years written every 6.5 hours."""
    ua, va, czen = g['eddies_ua'], g['eddies_va'], g['eddies_czen']
    n = ua.shape[0] // 11 * 11
    u, v = ua[:n].reshape(-1, 11, 32, 64), va[:n].reshape(-1, 11, 32, 64)
    six = np.hypot(u, v).mean(axis=1)
    three = np.hypot(u.mean(axis=1), v.mean(axis=1))
    cut = CT.half_sunlit(window) * np.cos(np.radians(grid['lat']))[:, None]
    day = czen[:n].reshape(-1, 11, 32, 64).mean(axis=1) > cut
    out = {}
    for part, m in (('land', grid['lsm']), ('water', ~grid['lsm'])):
        for name, dm in (('day', day), ('night', ~day)):
            sel = dm & m[None]
            w = np.broadcast_to(grid['area'], six.shape)[sel]
            out.setdefault(part, {})[name] = float((six[sel] * w).sum() / (three[sel] * w).sum())
    return out


def surface(g, key, sidecar, grid, ratios, vector, night_ratio=None):
    """A run's means at the surface: day mask, air at 2 m, the ground, humidity, pressure, the sky's longwave, the
    sunlight reaching the ground and the wind at 2 m for the reference surface and the plants."""
    q = CT.window_humidity(g[f'{key}_hus'].astype(np.float64))
    n = q.shape[0]
    take = lambda k: g[f'{key}_{k}'][:n].astype(np.float64)
    lat = grid['lat']
    window = sidecar['runs'][key]['sun_window_rad']
    day = take('czen') > CT.half_sunlit(window) * np.cos(np.radians(lat))[:, None]
    tas, ts, ps = take('tas'), take('ts'), take('ps') * 100.0
    land = np.broadcast_to(grid['lsm'], tas.shape)
    emissivity = np.where(land, LAND_EMISSIVITY, WATER_EMISSIVITY)
    sky = take('rls') / emissivity + STEFAN_BOLTZMANN * ts ** 4
    sun = take('rss') / (1.0 - np.where(land, LAND_ALBEDO, WATER_ALBEDO))
    speed = np.hypot(take('ua'), take('va'))
    factor = np.ones_like(speed)
    for part, m in (('land', land), ('water', ~land)):
        for name, dm in (('day', day), ('night', ~day)):
            r = ratios[part][name]['central']
            if night_ratio is not None and part == 'land' and name == 'night':
                r = night_ratio
            vf = vector[part][name] if key != 'eddies' else 1.0
            factor = np.where(m & dm, vf * r, factor)
    u2 = speed * factor * TEN_TO_TWO_M
    return dict(day=day, tas=tas, ts=ts, ps=ps, q=np.minimum(q, D.saturation_humidity(tas, ps)), sky=sky, sun=sun,
                u2=u2, land=land, n=n)


def cell_mean(values, mask=None):
    v = values if mask is None else np.where(mask, values, 0.0)
    return v.mean(axis=0)


def demand_by_cell(s):
    """The reference surface's demand (mm a day), each cell's time mean, with its day and night parts."""
    e = D.demand((1.0 - D.REFERENCE_ALBEDO) * s['sun'], s['sky'], s['tas'], s['q'], s['ps'], s['u2'], s['day'])
    mm = e * DAY_S
    return dict(all=cell_mean(mm), day=cell_mean(mm, s['day']), night=cell_mean(mm, ~s['day']))


def night_by_cell(s, cool=False, width=D.LEAF_WIDTHS_M[1]):
    """Each cell's night: dew on the canopy, on an exposed leaf and on bare land, fog, the night's humidity and the
    sky's deficit, all as time means over the whole lunar day (mm a day) or means over the night's means."""
    night = ~s['day']
    idx = np.nonzero(night)
    t = s['tas'][idx]
    if cool:
        t = np.where(s['land'][idx], np.minimum(t, s['ts'][idx] + NIGHT_GROUND_OFFSET_K), t)
    q, p, sky, u2, ts = (s[k][idx] for k in ('q', 'ps', 'sky', 'u2', 'ts'))
    q = np.minimum(q, D.saturation_humidity(t, p))
    canopy, tc = D.canopy_dew(t, q, p, sky, u2, width=width)
    leaf, _ = D.leaf_dew(t, q, p, sky, u2, width=width)
    land = s['land'][idx]
    soil = np.where(land, D.soil_dew(ts, q, p, u2, t), 0.0)
    fog = np.where(land, D.saturated_at_ground(ts, q, p), q >= D.saturation_humidity(t, p))
    rh = q / D.saturation_humidity(t, p)
    depression = t - D.dew_point(q, p)
    deficit = STEFAN_BOLTZMANN * t ** 4 - sky

    def spread(values):
        out = np.zeros(night.shape)
        out[idx] = values
        return out

    nights = night.sum(axis=0)
    per_night = lambda v: spread(v).sum(axis=0) / np.maximum(nights, 1)
    return dict(canopy_mm_day=cell_mean(spread(canopy)) * DAY_S, leaf_mm_day=cell_mean(spread(leaf)) * DAY_S,
                soil_mm_day=cell_mean(spread(soil)) * DAY_S, fog_share=per_night(fog.astype(float)),
                rh=per_night(rh), dew_point_depression_k=per_night(depression), sky_deficit_w_m2=per_night(deficit),
                canopy_below_air_k=per_night(t - tc), dew_share=per_night((canopy > 0).astype(float)),
                air_k=per_night(t), wind_2m_m_s=per_night(u2), pressure_pa=per_night(p))


def dew_onset(state, rhs=(0.90, 0.95, 0.99)):
    """The night humidity at which a canopy begins to gather dew, and what it gathers a lunar night at higher
    humidity, for a mean night: its air (K), sky deficit (W/m2), wind at 2 m and pressure."""
    t, deficit, u2, p = (np.array([state[k]]) for k in ('air_k', 'sky_deficit_w_m2', 'wind_2m_m_s', 'pressure_pa'))
    sky = STEFAN_BOLTZMANN * t ** 4 - deficit

    def gathered(rh):
        e = rh * D.saturation_pressure(t)
        q = D.EPSILON * e / (p - (1.0 - D.EPSILON) * e)
        return float(D.canopy_dew(t, q, p, sky, u2)[0][0])

    lo, hi = 0.5, 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if gathered(mid) <= 0.0 else (lo, mid)
    return dict(onset_rh=0.5 * (lo + hi),
                canopy_mm_lunar_night={f'{rh:g}': gathered(rh) * LUNAR_NIGHT_H * 3600.0 for rh in rhs})


# --------------------------------------------------------------------------------------------------------- belts
def transmissivities(table):
    t = table['setting']['transmissivity_m2_s_at_depth']
    central = t['central']['0 m']
    return dict(slow=t['slow']['0 m'], central=central, sandy=SANDY_FACTOR * central)


def cases(table):
    return {f'{name}_{int(r)}m': dict(transmissivity_m2_s=tv, reach_m=r)
            for name, tv in transmissivities(table).items() for r in REACH_M}


def shore_pixels(sh):
    """The 16 px/deg shores' arrays with each pixel's T21 cell, row and pixel area."""
    row_area = sh['row_area_km2']
    t21_row, t21_col = sh['t21_row'].astype(np.int64), sh['t21_col'].astype(np.int64)
    cidx = sh['coast_index'].astype(np.int64)
    pidx = sh['pixel_index'].astype(np.int64)
    coast_row, coast_col = np.divmod(cidx, S.COLS)
    pix_row, pix_col = np.divmod(pidx, S.COLS)
    coastal = np.zeros(S.ROWS * S.COLS, bool)
    coastal[cidx] = True
    inland = ~coastal[pidx]
    return dict(
        coast=dict(cell=t21_row[coast_row] * 64 + t21_col[coast_col], row=coast_row, col=coast_col,
                   area=row_area[coast_row], edge=sh['coast_edge_km'].astype(np.float64),
                   across=sh['coast_across_km'].astype(np.float64), rise=sh['coast_rise_m'].astype(np.float64),
                   island=sh['coast_island'], lake=sh['coast_lake'], body=sh['coast_island_body']),
        inland=dict(cell=(t21_row[pix_row] * 64 + t21_col[pix_col])[inland], row=pix_row[inland], col=pix_col[inland],
                    area=row_area[pix_row][inland], distance=sh['pixel_distance_km'][inland].astype(np.float64),
                    step_rows=sh['pixel_step_rows'][inland], step_cols=sh['pixel_step_cols'][inland],
                    height=sh['pixel_height_m'][inland].astype(np.float64), water=sh['pixel_water_class'][inland],
                    island=sh['pixel_island'][inland], lake=sh['pixel_lake'][inland],
                    body=sh['pixel_island_body'][inland]))


def classify(water, island):
    """Shore class index (0..4 in SHORE_CLASSES) from the water class (1..4) and the island flag."""
    w = np.asarray(water).astype(np.int64)
    return np.where(w == S.SEA, np.where(island, 1, 0), w)


def belts_for(px, uptake_cell, transmissivity, reach):
    """Belt areas, lengths and widths by shore class and T21 cell for one case. uptake_cell is m/s by T21 cell."""
    c_cell = uptake_cell.ravel() / (2.0 * transmissivity)
    co = px['coast']
    c = c_cell[co['cell']][:, None]
    areas, widths, lengths = B.first_ring(co['edge'], co['across'], co['rise'], c, reach)
    total = areas.sum(axis=1)
    scale = np.where(total > co['area'], co['area'] / np.maximum(total, 1e-30), 1.0)
    areas = areas * scale[:, None]
    cls = np.stack([np.where(co['island'], 1, 0), np.full(total.size, 2), np.full(total.size, 3),
                    np.full(total.size, 4)], axis=1)
    il = px['inland']
    share = B.inland(il['distance'], il['step_rows'], il['step_cols'], il['height'], c_cell[il['cell']], reach)
    il_area = share * il['area']
    il_cls = classify(il['water'], il['island'])
    belt = np.zeros((5, 2048))
    length = np.zeros((5, 2048))
    for k in range(4):
        np.add.at(belt, (cls[:, k], co['cell']), areas[:, k])
        np.add.at(length, (cls[:, k], co['cell']), lengths[:, k])
    np.add.at(belt, (il_cls, il['cell']), il_area)
    strips = dict(width=widths, length=lengths, cls=cls, cell=co['cell'])
    return dict(belt=belt, length=length, strips=strips, inland_area=il_area, inland_cls=il_cls,
                coast_areas=areas)


def weighted_percentiles(x, w, qs=(0.1, 0.5, 0.9)):
    x, w = np.asarray(x, float).ravel(), np.asarray(w, float).ravel()
    keep = w > 0
    if not keep.any():
        return [None] * len(qs)
    return [CT.wquantile(x[keep], w[keep], q) for q in qs]


def belt_summary(runs, uptake_cell, masks, grid, shore_cells, cs):
    out = dict(cases=cs, by_class={}, deserts={})
    cell_lsm = grid['lsm'].ravel()
    alat = grid['alat'].ravel()
    up = uptake_cell.ravel()
    for name, case in cs.items():
        r = runs[name]
        st = r['strips']
        regions = {'all': np.ones(2048, bool)}
        for (mask_name, count), m in masks.items():
            regions[f'{mask_name}.{count}'] = m.ravel() & (cell_lsm if count == 'model_land' else True)
        for region, cells_in in regions.items():
            row = {}
            for k, cls in enumerate(SHORE_CLASSES):
                sel = (st['cls'] == k) & cells_in[st['cell']][:, None] & (st['length'] > 0)
                p10, p50, p90 = weighted_percentiles(st['width'][sel], st['length'][sel])
                length = float(r['length'][k][cells_in].sum())
                area = float(r['belt'][k][cells_in].sum())
                uptake = float((r['belt'][k][cells_in] * 1e6 * up[cells_in]).sum())       # m3/s
                row[cls] = dict(shore_km=length, belt_km2=area, width_m=dict(p10=p10, median=p50, p90=p90),
                                first_ring_km2=float((r['coast_areas'][sel]).sum()),
                                uptake_km3_yr=uptake * YEAR_S / 1e9,
                                supply_m3_yr_per_km=(uptake * YEAR_S / length) if length > 0 else None)
            land16 = float(shore_cells['land'].ravel()[cells_in].sum())
            total_belt = float(r['belt'][:, cells_in].sum())
            island16 = float(shore_cells['island_land'].ravel()[cells_in].sum())
            row['island_coast']['island_land_km2'] = island16
            row['island_coast']['share_of_island_land'] = (row['island_coast']['belt_km2'] / island16
                                                           if island16 > 0 else None)
            entry = dict(by_class=row, land_km2=land16, belt_km2=total_belt,
                         share_within_reach=(total_belt / land16) if land16 > 0 else None,
                         uptake_km3_yr=float(sum(v['uptake_km3_yr'] for v in row.values())),
                         shore_km=float(sum(v['shore_km'] for v in row.values())))
            if region != 'all':
                entry['bands'] = {}
                for lo, hi in BANDS:
                    b = cells_in & band_of(alat, lo, hi)
                    l16 = float(shore_cells['land'].ravel()[b].sum())
                    bb = float(r['belt'][:, b].sum())
                    entry['bands'][f'{lo}-{hi}'] = dict(land_km2=l16, belt_km2=bb,
                                                        share_within_reach=(bb / l16) if l16 > 0 else None)
                entry['land_per_shore_km'] = land16 / entry['shore_km'] if entry['shore_km'] > 0 else None
            out['deserts' if region != 'all' else 'by_class'].setdefault(region, {})[name] = entry
    return out


# --------------------------------------------------------------------------------------- the check at 100 m
def aeqd(lat0, lon0, lat, lon):
    """Metres east and north of (lat0, lon0) by the azimuthal equidistant projection geography.coast_terrain uses."""
    p0, l0, p, l = map(np.radians, (lat0, lon0, lat, lon))
    c = np.arccos(np.clip(np.sin(p0) * np.sin(p) + np.cos(p0) * np.cos(p) * np.cos(l - l0), -1.0, 1.0))
    az = np.arctan2(np.sin(l - l0) * np.cos(p), np.cos(p0) * np.sin(p) - np.sin(p0) * np.cos(p) * np.cos(l - l0))
    return MOON_RADIUS * c * np.sin(az), MOON_RADIUS * c * np.cos(az)


def fine_belts(height, spacing_m, c, reach):
    """Belt area (km2) and shore length (km) on a local grid at fine spacing, with the same rules as at 16 px/deg:
    common-level water where the height is below 0."""
    water = height < 0.0
    dist, (ii, jj) = ndimage.distance_transform_edt(~water, sampling=spacing_m, return_indices=True)
    rows, cols = np.indices(height.shape)
    steps_r, steps_c = ii - rows, jj - cols
    land = ~water
    coastal = np.zeros_like(land)
    edge = np.zeros(height.shape)
    for dr, dc in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        nb = np.zeros_like(water)
        src = water[max(dr, 0):height.shape[0] + min(dr, 0), max(dc, 0):height.shape[1] + min(dc, 0)]
        nb[max(-dr, 0):height.shape[0] + min(-dr, 0), max(-dc, 0):height.shape[1] + min(-dc, 0)] = src
        facing = land & nb
        coastal |= facing
        edge += facing * spacing_m / 1000.0
    pixel_area = (spacing_m / 1000.0) ** 2
    across = np.full(height.shape, spacing_m / 1000.0)
    areas, widths, lengths = B.first_ring(edge, across, height, c, reach)
    first = np.minimum(areas, pixel_area)
    inner = land & ~coastal
    inland_area = np.zeros(height.shape)
    inland_area[inner] = B.inland(dist[inner] / 1000.0, steps_r[inner], steps_c[inner], height[inner], c,
                                  reach) * pixel_area
    return np.where(coastal, first, 0.0) + inland_area, np.where(coastal, lengths, 0.0), np.where(coastal, widths, 0.0)


RESOLUTIONS = (1, 2, 4)                 # the local coasts at 100 m and in blocks of 2 and 4 (200 and 400 m)


def subset(px, coast, inland):
    """The 16 px/deg shores' pixels kept by two masks, one over the pixels beside the water and one over the rest."""
    return dict(coast={k: v[coast] for k, v in px['coast'].items()},
                inland={k: v[inland] for k, v in px['inland'].items()})


def blocks(height, n):
    """A local grid's heights in blocks of n by n pixels, and the block centres' offsets (m) from the grid's centre
    in units of the fine spacing."""
    m = (height.shape[0] // n) * n
    h = height[:m, :m].reshape(m // n, n, m // n, n).mean(axis=(1, 3))
    offsets = (np.arange(m // n) * n + 0.5 * (n - 1)) - 0.5 * (height.shape[0] - 1)
    return h, offsets


def local_check(px, sh, uptake_cell, cs, grid, data_coasts):
    """The belts on the four local coasts at 100 m, and in blocks of 200 and 400 m, against the 16 px/deg pixels in
    the same squares, the common-level water's shores only and one uptake for both, the tile centre's cell's."""
    out = {}
    t21_lat = grid['lat']
    row_lat = sh['row_lat']
    for name, z in data_coasts.items():
        meta = json.loads(str(z['metadata']))
        lon0, lat0 = meta['centre_lon_lat_deg']
        spacing = meta['spacing_m']
        inner = meta['half_width_m'] - LOCAL_MARGIN_KM * 1000.0
        h = z['height_m'].astype(np.float64)
        row = int(np.argmin(np.abs(t21_lat - lat0)))
        col = int(np.floor(((lon0 % 360.0) + 180.0 / 64.0) / (360.0 / 64.0))) % 64
        e = float(uptake_cell[row, col])
        co, il = px['coast'], px['inland']

        def inside(rows, cols):
            lat_ = row_lat[rows]
            near = np.abs(lat_ - lat0) < 3.0
            x_, y_ = aeqd(lat0, lon0, lat_[near], (cols[near] + 0.5) / S.PPD)
            keep = np.zeros(lat_.size, bool)
            keep[np.flatnonzero(near)] = (np.abs(x_) <= inner) & (np.abs(y_) <= inner)
            return keep

        in_c = inside(co['row'], co['col'])
        in_i = inside(il['row'], il['col']) & np.isin(il['water'], (S.SEA, S.LAKE_COMMON))
        tile = subset(px, in_c, in_i)
        constant = np.full((32, 64), e)
        xx, yy = np.meshgrid(z['x_m'], z['y_m'], indexing='xy')
        square = (np.abs(xx) <= inner) & (np.abs(yy) <= inner)
        entry = dict(centre_lat_lon_deg=[lat0, lon0], cell=[row, col], inner_half_width_km=inner / 1000.0,
                     uptake_mm_day=e * DAY_S * 1000, land_share_100m=float((h[square] >= 0).mean()), cases={})
        for cname, case in cs.items():
            c = e / (2.0 * case['transmissivity_m2_s'])
            r = belts_for(tile, constant, case['transmissivity_m2_s'], case['reach_m'])
            coarse = float(r['coast_areas'][:, :2].sum() + r['inland_area'].sum())
            v = dict(belt_km2_16ppd=coarse, shore_km_16ppd=float(r['strips']['length'][:, :2].sum()))
            for n in RESOLUTIONS:
                hb, off = blocks(h, n) if n > 1 else (h, np.arange(h.shape[0]) - 0.5 * (h.shape[0] - 1))
                xb = off * spacing
                keep = (np.abs(xb)[None, :] <= inner) & (np.abs(xb)[:, None] <= inner)
                fine, lengths, widths = fine_belts(hb, spacing * n, c, case['reach_m'])
                tag = f'{int(spacing * n)}m'
                v[f'belt_km2_{tag}'] = float(fine[keep].sum())
                v[f'shore_km_{tag}'] = float(lengths[keep].sum())
                if n == 1:
                    sel = keep & (lengths > 0)
                    v['width_m_100m'] = dict(zip(('p10', 'median', 'p90'),
                                                 weighted_percentiles(widths[sel], lengths[sel])))
            v['ratio_16ppd_over_100m'] = v['belt_km2_16ppd'] / v['belt_km2_100m'] if v['belt_km2_100m'] > 0 else None
            v['ratio_200m_over_100m'] = v['belt_km2_200m'] / v['belt_km2_100m'] if v['belt_km2_100m'] > 0 else None
            v['ratio_400m_over_100m'] = v['belt_km2_400m'] / v['belt_km2_100m'] if v['belt_km2_100m'] > 0 else None
            entry['cases'][cname] = v
        out[name] = entry
    ratios = {c: [out[n]['cases'][c]['ratio_16ppd_over_100m'] for n in out] for c in cs}
    halving = {c: [out[n]['cases'][c]['ratio_200m_over_100m'] for n in out] for c in cs}
    return dict(coasts=out, ratio_range={c: [min(v), max(v)] for c, v in ratios.items()},
                ratio_200m_over_100m_range={c: [min(v), max(v)] for c, v in halving.items()})


def corrected_shares(deserts_out, check):
    """Each desert's share within reach divided, case by case, by that case's own range of 16 px/deg over 100 m on
    the local coasts: the lower end by its largest ratio and the upper end by its smallest; over every case and over
    the central crust's two reaches."""
    out = {}
    for region, cases_ in deserts_out.items():
        rows = {name: [e['share_within_reach'] / check['ratio_range'][name][1],
                       e['share_within_reach'] / check['ratio_range'][name][0]] for name, e in cases_.items()}
        central = [v for k, v in rows.items() if k.startswith('central')]
        out[region] = dict(by_case=rows, all_cases=[min(v[0] for v in rows.values()), max(v[1] for v in rows.values())],
                           central=[min(v[0] for v in central), max(v[1] for v in central)])
    return out


# ------------------------------------------------------------------------------------------- the 15 km water table
def water_table(grids, masks, grid):
    d = grids['depth_m_mature'].astype(np.float64)
    seep = grids['seepage_mm_yr_mature'].astype(np.float64)
    lat, lon = grids['lat_deg'].astype(np.float64), grids['lon_deg'].astype(np.float64)
    r = MOON_RADIUS / 1000.0
    area = (r * r * np.radians(0.5) * (np.sin(np.radians(lat + 0.25)) - np.sin(np.radians(lat - 0.25))))[:, None] \
        * np.ones((1, lon.size))
    land = np.isfinite(d)
    _, weights = np.polynomial.legendre.leggauss(32)
    edges = np.concatenate([[-1.0], -1.0 + np.cumsum(weights)])
    k = np.clip(np.searchsorted(edges, np.sin(np.radians(lat))) - 1, 0, 31)
    south_first = np.sort(grid['lat'])
    rows = np.array([int(np.argmin(np.abs(grid['lat'] - south_first[i]))) for i in k])
    cols = (np.floor((lon + 180.0 / 64.0) / (360.0 / 64.0)).astype(int) % 64)
    cell = rows[:, None] * 64 + cols[None, :]
    lsm = grid['lsm'].ravel()
    alat = np.abs(lat)[:, None] * np.ones((1, lon.size))
    tiled = np.concatenate([~land, ~land, ~land], axis=1)       # cells from the nearest water cell, round the seam
    steps = ndimage.distance_transform_cdt(~tiled, metric='chessboard')[:, lon.size:2 * lon.size]
    out = dict(check=dict(land_km2=float(area[land].sum()), median_depth_m=float(np.median(d[land])),
                          seeping_share=float(area[land & (d <= 0.0)].sum() / area[land].sum())))
    for (mask_name, count), m in masks.items():
        cells_in = m.ravel() & (lsm if count == 'model_land' else True)
        sel = land & cells_in[cell]
        a = area[sel].sum()

        def shares(s):
            aa = area[s].sum()
            return dict(land_km2=float(aa), at_ground=float(area[s & (d <= 0.0)].sum() / aa) if aa else None,
                        within_15_m=float(area[s & (d <= 15.0)].sum() / aa) if aa else None,
                        within_25_m=float(area[s & (d <= 25.0)].sum() / aa) if aa else None,
                        median_depth_m=float(np.median(d[s])) if s.any() else None)
        entry = shares(sel)
        entry['seepage_km3_yr'] = float((seep[sel] / 1000.0 * area[sel] * 1e6).sum() / 1e9)
        shallow = sel & (d <= 25.0)
        sa = area[shallow].sum()
        entry['within_25_m_beside_water'] = dict(
            one_cell=float(area[shallow & (steps <= 1)].sum() / sa) if sa else None,
            two_cells=float(area[shallow & (steps <= 2)].sum() / sa) if sa else None)
        for state in ('young_1000_yr', 'young_10000_yr'):
            dy_ = grids[f'depth_m_{state}'].astype(np.float64)
            entry[state] = dict(at_ground=float(area[sel & (dy_ <= 0.0)].sum() / a) if a else None,
                                within_25_m=float(area[sel & (dy_ <= 25.0)].sum() / a) if a else None,
                                median_depth_m=float(np.median(dy_[sel])) if sel.any() else None)
        entry['bands'] = {f'{lo}-{hi}': shares(sel & band_of(alat, lo, hi)) for lo, hi in BANDS}
        out.setdefault(mask_name, {})[count] = entry
    return out


# -------------------------------------------------------------------------------------------------- dew and fog
def night_classes(grid, masks, per_cell, rain, shore_lengths, mask_name='rain_mask'):
    """The night's figures at 301 K by T21 class (the shore land in the model's water cells, coastal land cells,
    inland land cells), in and out of the rain mask's deserts, weighted by land area, with the rain beside them."""
    wl = grid['wl']
    desert = masks[(mask_name, 'every_cell')]
    classes = dict(shore_land_in_water_cells=grid['water_with_land'], coastal_land_cells=grid['coastal'],
                   inland_land_cells=grid['inland'])
    out = {}
    for cname, cm in classes.items():
        for dname, dm in (('desert', desert), ('not_desert', ~desert)):
            m = cm & dm
            w = wl * m
            if w.sum() <= 0:
                continue
            row = dict(land_km2=float(w.sum() * 4.0 * math.pi * (MOON_RADIUS / 1000.0) ** 2), cells=int(m.sum()),
                       rain_mm_lunar_day=float((w * rain).sum() / w.sum() * LUNAR_DAY_DAYS))
            for k, v in per_cell.items():
                mean = float((w * v).sum() / w.sum())
                if k.endswith('_mm_day'):
                    row[k.replace('_mm_day', '_mm_lunar_day')] = mean * LUNAR_DAY_DAYS
                else:
                    row[k] = mean
            rain_ = row['rain_mm_lunar_day']
            row['canopy_over_rain'] = row['canopy_mm_lunar_day'] / rain_ if rain_ else None
            out.setdefault(cname, {})[dname] = row
    # the 16 px/deg shores' classes, by the length of each in each desert cell
    by_shore = {}
    for k, cls in enumerate(SHORE_CLASSES):
        w = shore_lengths[k].reshape(32, 64) * desert
        if w.sum() <= 0:
            continue
        by_shore[cls] = dict(shore_km=float(w.sum()),
                             canopy_mm_lunar_day=float((w * per_cell['canopy_mm_day']).sum() / w.sum()
                                                       * LUNAR_DAY_DAYS),
                             rh=float((w * per_cell['rh']).sum() / w.sum()),
                             rain_mm_lunar_day=float((w * rain).sum() / w.sum() * LUNAR_DAY_DAYS))
    out['desert_shores_by_class'] = by_shore
    return out


def earth_check():
    """The dew balance on an Earth night: 10 hours under a clear sky's 87 W/m2 deficit (the radiative-convective
    column's Earth surface) and under a cloudier 40 W/m2, at 90 and 96% humidity and 1 m/s."""
    out = []
    for deficit in (87.0, 40.0):
        for rh in (0.90, 0.96):
            c, l = D.earth_night(rh, deficit, 10.0)
            out.append(dict(deficit_w_m2=deficit, rh=rh, canopy_mm=float(c), leaf_mm=float(l)))
    return out


def ring_pressures(rings):
    """Each ring's surface pressure (Pa) from its own case file in the CM1 runs, with the files' sha256."""
    out, files = {}, {}
    for key, ring in rings.items():
        path = f'{CM1_RUNS}/{ring["case"]}/case.json'
        out[key] = float(json.loads((ROOT / path).read_text())['build']['air']['surface_pa'])
        files[path] = digest(path)
    return out, files


def rings_and_jura(rings, jura, gcm_45_deficit, pressures):
    """The rings' land nights and the Jura box's air beside the GCM's: night humidity, the ground below the 2 m air,
    the wind, fog, and the canopy's dew on each ring's mean night under the GCM's 45-degree night deficit, at each
    ring's own surface pressure."""
    out = {}
    for key, ring in rings.items():
        b = ring['by_hour_angle']
        ha = np.abs(np.array(b['hour_angle_deg'], float))
        night = ha > 90.0
        L = b['land']
        t = np.nanmean(np.array(L['air_2m_c'], float)[night]) + 273.15
        rh = np.nanmean(np.array(L['humidity_2m'], float)[night])
        ground = np.nanmean(np.array(L['ground_c'], float)[night]) + 273.15
        u10 = np.nanmean(np.array(L['wind_10m_m_s'], float)[night])
        p = pressures[key]
        e = rh * D.saturation_pressure(t)
        q = D.EPSILON * e / (p - (1.0 - D.EPSILON) * e)
        sky = STEFAN_BOLTZMANN * t ** 4 - gcm_45_deficit
        canopy, tc = D.canopy_dew(np.array([t]), np.array([q]), np.array([p]), np.array([sky]),
                                  np.array([u10 * TEN_TO_TWO_M]))
        out[key] = dict(night_air_c=float(t - 273.15), night_rh=float(rh), ground_below_air_k=float(t - ground),
                        night_wind_10m_m_s=float(u10), dew_point_depression_k=float(t - D.dew_point(q, p)),
                        fog_land=float(ring['fog']['land']),
                        land_evaporation_mm_day=float(ring['water_budget']['evaporation_mm_day']['land']),
                        canopy_mm_lunar_night=float(canopy[0] * LUNAR_NIGHT_H * 3600.0),
                        canopy_below_air_k=float(t - tc[0]), surface_pa=p)
    dist = {g['group']: dict(dew_point_day_c=g['dewpoint_c']['day'], dew_point_night_c=g['dewpoint_c']['night'],
                             wind_10m_median_m_s=g['wind_10m_m_s']['median'])
            for g in jura['by_distance'] + jura['waters']}
    onshore = jura['through_the_day']['onshore_wind_m_s']
    ha = np.array(jura['through_the_day']['hour_angle_deg'], float)
    jb = dict(by_distance=dist,
              onshore_wind_m_s={k: dict(day_max=float(max(v)), night_mean=float(np.mean(np.array(v)[np.abs(ha) > 90])))
                                for k, v in onshore.items()},
              land_dew_point_below_sea_k=float(np.mean([dist['the bay']['dew_point_night_c'],
                                                        dist['the open sea']['dew_point_night_c']])
                                               - dist['0-7 km from the sea']['dew_point_night_c']))
    return dict(rings=out, jura_box=jb, gcm_45_night_deficit_w_m2=gcm_45_deficit)


def eddies_check(g, sidecar, grid, ratios, vector, desert):
    """Dew in the design case's 6.5-hour means against the 3-day means their elevens make, over the same two years:
    how much the 3-day means smooth the night."""
    s6 = surface(g, 'eddies', sidecar, grid, ratios, vector)
    n6 = night_by_cell(s6)
    # 3-day means of the 6.5-hour fields: the humidity from each 6.5-hour mean's start and end, averaged
    k = s6['n'] // 11 * 11
    agg = {}
    for key in ('tas', 'ts', 'ps', 'q', 'sky', 'sun'):
        agg[key] = s6[key][:k].reshape(-1, 11, 32, 64).mean(axis=1)
    czen = g['eddies_czen'][:k].astype(np.float64).reshape(-1, 11, 32, 64).mean(axis=1)
    window3 = sidecar['runs']['design']['sun_window_rad']
    agg['day'] = czen > CT.half_sunlit(window3) * np.cos(np.radians(grid['lat']))[:, None]
    ua = g['eddies_ua'][:k].astype(np.float64).reshape(-1, 11, 32, 64).mean(axis=1)
    va = g['eddies_va'][:k].astype(np.float64).reshape(-1, 11, 32, 64).mean(axis=1)
    land = np.broadcast_to(grid['lsm'], agg['tas'].shape)
    factor = np.ones_like(ua)
    for part, m in (('land', land), ('water', ~land)):
        for name, dm in (('day', agg['day']), ('night', ~agg['day'])):
            factor = np.where(m & dm, vector[part][name] * ratios[part][name]['central'], factor)
    agg['u2'] = np.hypot(ua, va) * factor * TEN_TO_TWO_M
    agg['land'] = land
    agg['q'] = np.minimum(agg['q'], D.saturation_humidity(agg['tas'], agg['ps']))
    agg['n'] = agg['tas'].shape[0]
    n3 = night_by_cell(agg)
    wl = grid['wl'] * grid['lsm']
    out = {}
    for name, m in (('land', np.ones_like(desert)), ('desert', desert)):
        w = wl * m
        out[name] = {f'{k}_{res}': float((w * v[k]).sum() / w.sum())
                     for res, v in (('6_5_h', n6), ('3_day', n3)) for k in ('canopy_mm_day', 'soil_mm_day', 'rh',
                                                                            'fog_share', 'dew_share')}
    pet6, pet3 = demand_by_cell(s6), demand_by_cell(agg)
    out['demand_land_mm_day'] = dict(six_5_h=float((wl * pet6['all']).sum() / wl.sum()),
                                     three_day=float((wl * pet3['all']).sum() / wl.sum()))
    return out


# ---------------------------------------------------------------------------------------------- the slice's cost
SLICE_RUNS = ('ring_45n_corrected_twilight', 'ring_45s_corrected_twilight', 'ring_45n_corrected', 'box_jura')
SLICE_LOG_SEGMENT = 30
SLICE_UPPER_WIND_M_S = 25.0              # the rings' fastest air aloft, near 70 km, for the horizontal Courant limit
SLICE_DESIGNS = dict(uniform_2_km=(2000.0, 300), stretched_500_m=(500.0, 212), uniform_1_km=(1000.0, 600),
                     stretched_250_m=(250.0, 263))
SLICE_MODEL_DAYS = 59.0


def slice_log(path):
    """A CM1 segment log's mean step (s), mean acoustic sub-steps and the share of its time in the acoustic steps."""
    steps, sound, total = [], [], None
    for line in (ROOT / path).read_text().splitlines():
        if line.startswith(' cflmax,ksmax,dt,nsound:'):
            parts = line.split(':')[1].split()
            steps.append((float(parts[2]), float(parts[3])))
        elif line.strip().startswith('sound    :'):
            sound.append(float(line.split('%')[0].split()[-1]) / 100.0)
    a = np.array(steps)
    return dict(mean_step_s=float(a[:, 0].mean()), mean_sub_steps=float(a[:, 1].mean()), sound_share=sound[-1],
                steps=int(len(a)))


def slice_cost(runs_dir='climate/crm/runs'):
    """The 2-D slice's wall time scaled from the CM1 runs' own records (progress.json and a mid-run segment's log):
    wall time goes as the columns, the steps a model day and the acoustic sub-steps' share. The step stays at the
    runs' own where the horizontal Courant number at the finest spacing allows, with the fastest air aloft, and the
    sub-steps keep the runs' acoustic Courant number. Two lunar days, the first a spin-up, as the rings ran."""
    measured, files = {}, {}
    for run in SLICE_RUNS:
        prog = f'{runs_dir}/{run}/progress.json'
        log = f'{runs_dir}/{run}/cm1_segment_{SLICE_LOG_SEGMENT:03d}.log'
        p = json.loads((ROOT / prog).read_text())
        cfg = json.loads((ROOT / f'{runs_dir}/{run}/case.json').read_text())
        walls = [x['wall_s'] for x in p['segments']]
        case = f'{runs_dir}/{run}/case.json'
        files.update({prog: digest(prog), log: digest(log), case: digest(case)})
        g = cfg['grid']
        measured[run] = dict(columns=int(g['nx'] * g.get('ny', 1)), spacing_m=float(g['dx_m']), levels=int(g['nz']),
                             threads=sorted({x['threads'] for x in p['segments']}),
                             median_wall_s_per_model_day=float(np.median(walls)), wall_hours=float(sum(walls) / 3600.0),
                             **slice_log(log))
    rings = [measured[r] for r in SLICE_RUNS if r.startswith('ring')]
    per_column = [r['median_wall_s_per_model_day'] / r['columns'] for r in rings]
    step = float(np.mean([r['mean_step_s'] for r in rings]))
    sub = float(np.mean([r['mean_sub_steps'] for r in rings]))
    share = float(np.mean([r['sound_share'] for r in rings]))
    spacing = rings[0]['spacing_m']
    designs = {}
    for name, (fine, columns) in SLICE_DESIGNS.items():
        dt = min(step, fine / SLICE_UPPER_WIND_M_S)
        subs = sub * (dt / step) * (spacing / fine)
        factor = (1.0 - share) + share * subs / sub
        wall = [c * columns * (step / dt) * factor for c in (min(per_column), max(per_column))]
        designs[name] = dict(finest_spacing_m=fine, columns=columns, step_s=dt, sub_steps=subs,
                             wall_s_per_model_day=wall, wall_hours=[w * SLICE_MODEL_DAYS / 3600.0 for w in wall])
    return dict(measured=measured, designs=designs, model_days=SLICE_MODEL_DAYS,
                upper_wind_m_s=SLICE_UPPER_WIND_M_S, threads=6), files


# -------------------------------------------------------------------------------------------------------- the run
def results():
    data, sidecar = load()
    targets, table = data[TARGETS], data[TABLE]
    base = basis(targets)
    f = base['share_of_the_way']
    g = np.load(ROOT / GCM)
    sh = np.load(ROOT / SHORES)
    sh_meta = json.loads(str(sh['metadata']))
    if sh_meta.get('schema') != S.SCHEMA:
        raise ValueError('the shores product has another schema')
    grid = t21(g)
    shore_cells = {k: sh[f'cell_{k}'] for k in ('land', 'island_land') + tuple(S.WATER_CLASSES.values())}
    rain, masks, desert_out = deserts(g, grid, f, base, shore_cells, targets)
    rings = {k: data[p] for k, p in RINGS.items()}
    ratios = ring_ratios(rings)
    vector = vector_speed_factor(g, sidecar['runs']['eddies']['sun_window_rad'], grid)
    # the plants' demand and the night in each run, then at 301 K cell by cell
    per_run, nights, demands = {}, {}, {}
    for key in ('design', 'test'):
        s = surface(g, key, sidecar, grid, ratios, vector)
        demands[key] = demand_by_cell(s)
        nights[key] = dict(central=night_by_cell(s), cool=night_by_cell(s, cool=True),
                           **{f'leaf_{int(w * 100)}_cm': night_by_cell(s, width=w) for w in D.LEAF_WIDTHS_M
                              if w != D.LEAF_WIDTHS_M[1]})
        for end, r in zip(('low', 'high'), ratios['land']['night']['range']):
            nights[key][f'night_wind_{end}'] = night_by_cell(surface(g, key, sidecar, grid, ratios, vector, r))
        per_run[key] = s['n']
        del s
    pet = {k: interpolate(demands['design'][k], demands['test'][k], f) for k in ('all', 'day', 'night')}
    uptake = np.maximum(pet['all'] - rain, 0.0) / 1000.0 / DAY_S           # m/s by T21 cell
    night301 = {variant: {k: interpolate(nights['design'][variant][k], nights['test'][variant][k], f)
                          for k in nights['design'][variant]} for variant in nights['design']}
    # belts
    px = shore_pixels(sh)
    cs = cases(table)
    runs = {name: belts_for(px, uptake, c['transmissivity_m2_s'], c['reach_m']) for name, c in cs.items()}
    belt = belt_summary(runs, uptake, masks, grid, shore_cells, cs)
    coast_npz = {n: np.load(ROOT / p) for n, p in COASTS.items()}
    belt['resolution_check'] = local_check(px, sh, uptake, cs, grid, coast_npz)
    belt['within_reach_corrected'] = corrected_shares(belt['deserts'], belt['resolution_check'])
    jura = data[JURA_SHORES]
    belt['jura_island_screen'] = dict(belts_km2={k: v['area_km2'] for k, v in jura['belts'].items()},
                                      share={k: v['share'] for k, v in jura['belts'].items()})
    desert_cells = masks[('rain_mask', 'every_cell')]
    wl = grid['wl']
    dm = wl * desert_cells
    demand_out = dict(
        wind=dict(ring_ratios_10m_over_1km=ratios, vector_speed_factor=vector, ten_to_two_m=TEN_TO_TWO_M),
        land_mm_day={k: float((wl * grid['lsm'] * v).sum() / (wl * grid['lsm']).sum()) for k, v in pet.items()},
        desert_mm_day={k: float((dm * v).sum() / dm.sum()) for k, v in pet.items()},
        desert_uptake_mm_day=float((dm * uptake * 1000.0 * DAY_S).sum() / dm.sum()),
        desert_rain_mm_day=float((dm * rain).sum() / dm.sum()),
        desert_uptake_percentiles_mm_day=weighted_percentiles(uptake * 1000.0 * DAY_S, dm),
        by_run={k: dict(land_mm_day=float((wl * grid['lsm'] * v['all']).sum() / (wl * grid['lsm']).sum()),
                        desert_mm_day=float((dm * v['all']).sum() / dm.sum())) for k, v in demands.items()})
    med_uptake = demand_out['desert_uptake_percentiles_mm_day'][1] / 1000.0 / DAY_S
    demand_out['capillary_reach_m'] = capillary_reach(demand_out['desert_uptake_mm_day'] / 1000.0 / DAY_S)
    rain_every = belt['deserts']['rain_mask.every_cell']['central_15m']
    demand_out['hinterland_m3_yr_per_km_of_shore'] = (rain_every['land_per_shore_km'] * 1e6
                                                      * demand_out['desert_mm_day']['all'] / 1000.0
                                                      * JULIAN_YEAR_DAYS)
    belt['flat_shore'] = {name: dict(width_m=float(B.flat_width(c['transmissivity_m2_s'], med_uptake, c['reach_m'])),
                                     supply_m3_yr_per_km=float(math.sqrt(2.0 * c['transmissivity_m2_s'] * med_uptake
                                                                         * c['reach_m']) * 1000.0 * YEAR_S))
                          for name, c in cs.items()}
    belt['closed_lakes'] = closed_lakes(px, sh, uptake, cs, masks, grid, g, f, runs)
    belt['islands'] = islands(px, sh, cs, masks, grid, runs, uptake)
    grids = np.load(ROOT / TABLE_GRIDS)
    wt = water_table(grids, masks, grid)
    # the night
    lengths = runs['central_15m']['length']
    dew_out = dict(earth_check=earth_check(),
                   at_301_k={v: night_classes(grid, masks, night301[v], rain, lengths) for v in night301},
                   by_run={k: night_classes(grid, masks, nights[k]['central'], desert_out['rain_mm_day'][k], lengths)
                           for k in nights})
    lat45 = band_of(grid['alat'], 40, 50) & grid['lsm']
    deficit45 = float((wl * lat45 * nights['design']['central']['sky_deficit_w_m2']).sum() / (wl * lat45).sum())
    pressures, ring_cases = ring_pressures(rings)
    dew_out['cm1'] = rings_and_jura(rings, data[JURA_BOX], deficit45, pressures)
    dew_out['at_301_k_log_mask'] = night_classes(grid, masks, night301['central'],
                                                 rain_at(desert_out['rain_mm_day']['design'],
                                                         desert_out['rain_mm_day']['test'], f, 'log_linear'),
                                                 lengths, 'log_rain_mask')
    dew_out['onset'] = {cls: dew_onset(dew_out['at_301_k']['central'][cls]['desert'])
                        for cls in ('shore_land_in_water_cells', 'coastal_land_cells', 'inland_land_cells')}
    dew_out['eddies_check'] = eddies_check(g, sidecar, grid, ratios, vector, desert_cells)
    eddy = dew_out['eddies_check']['desert']
    factor = eddy['fog_share_6_5_h'] / eddy['fog_share_3_day']
    fog301 = [dew_out['at_301_k']['central'][c]['desert']['fog_share'] for c in ('coastal_land_cells',
                                                                                 'inland_land_cells')]
    night_s = LUNAR_NIGHT_H * 3600.0
    desert_rain = [dew_out['at_301_k']['central'][c]['desert']['rain_mm_lunar_day']
                   for c in ('shore_land_in_water_cells', 'coastal_land_cells', 'inland_land_cells')]
    most = max(fog301) * factor * night_s * FOG_LWC_KG_M3[1]
    dew_out['desert_fog_at_301_k'] = dict(
        share_of_night_3_day=fog301, six_hour_factor=factor,
        water_mm_lunar_night=[min(fog301) * factor * night_s * FOG_LWC_KG_M3[0] * FOG_DEPOSITION_M_S[0],
                              most * FOG_DEPOSITION_M_S[1]],
        water_mm_lunar_night_at_earths_most=most * FOG_DEPOSITION_HIGH_M_S,
        share_of_desert_rain_at_earths_most=most * FOG_DEPOSITION_HIGH_M_S / float(np.mean(desert_rain)))
    dew_out['fog_water_mm_lunar_night_per_percent'] = dict(
        lwc_kg_m3=list(FOG_LWC_KG_M3), deposition_m_s=list(FOG_DEPOSITION_M_S),
        range=[0.01 * LUNAR_NIGHT_H * 3600.0 * FOG_LWC_KG_M3[0] * FOG_DEPOSITION_M_S[0],
               0.01 * LUNAR_NIGHT_H * 3600.0 * FOG_LWC_KG_M3[1] * FOG_DEPOSITION_M_S[1]])
    cost, cost_files = slice_cost()
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, basis=base,
               deserts={k: v for k, v in desert_out.items() if k != 'rain_mm_day'}, demand=demand_out, belts=belt,
               water_table=wt, dew=dew_out, slice_plan=dict(cost=cost),
               means=dict(per_run, shores=sh_meta['seas'], islands_km2=sh_meta['island_km2']))
    out['producer'] = dict(study='research/studies/coast_water', files={p: digest(p) for p in CODE},
                           inputs={p: digest(p) for p in INPUTS}, ignored_inputs={p: digest(p) for p in IGNORED},
                           shores_inputs=sh_meta['inputs'], gcm_files=sidecar['files'],
                           cm1_records=dict(cost_files, **ring_cases))
    return out


ISLAND_SIZES_KM2 = (0.0, 10.0, 100.0, 1000.0, 10000.0, math.inf)
WOODED = 0.95                       # an island wholly within reach: its belts cover this much of its land


def pixel_area_at(sh, lat):
    """The 16 px/deg pixel's area (km2) at each latitude."""
    rows = np.clip(((90.0 - np.asarray(lat)) * S.PPD).astype(int), 0, S.ROWS - 1)
    return sh['row_area_km2'][rows]


def islands(px, sh, cs, masks, grid, runs, uptake):
    """Each island's land within reach of fresh water against its whole land: how many islands, and how much island
    land, lie wholly within reach (the strips cover WOODED of their land), by size and in each mask's deserts (an
    island's cell is its centroid's). Beside it, the disc rule: an island of the same area as a disc, wooded
    throughout where its highest pixel and the fall at its middle under uptake everywhere, E a^2 / 4T, stay within
    the reach."""
    area = sh['island_area_km2']
    ids = np.arange(area.size)
    real = (ids > 0) & (area > 0)
    lat, lon = sh['island_lat'], sh['island_lon']
    row = np.array([int(np.argmin(np.abs(grid['lat'] - x))) for x in lat])
    col = (np.floor((lon + 180.0 / 64.0) / (360.0 / 64.0)).astype(int) % 64)
    e = uptake[row, col]
    radius = np.sqrt(area * 1e6 / np.pi)
    top = sh['island_top_m'].astype(np.float64)
    co, il = px['coast'], px['inland']
    small = real & (area < ISLAND_SIZES_KM2[1])
    pixels = area / pixel_area_at(sh, lat)
    out = dict(count=int(real.sum()), land_km2=float(area[real].sum()), wooded_share=WOODED,
               small=dict(count=int(small.sum()), mean_km2=float(area[small].mean()),
                          pixels=dict(zip(('p10', 'median', 'p90'),
                                          [float(np.percentile(pixels[small], q)) for q in (10, 50, 90)]))),
               cases={})
    deserts_ = {'all': real, 'desert': real & masks[('rain_mask', 'every_cell')][row, col],
                'log_desert': real & masks[('log_rain_mask', 'every_cell')][row, col]}
    for name, case in cs.items():
        r = runs[name]
        belt = np.bincount(co['body'], weights=r['coast_areas'].sum(axis=1), minlength=area.size)
        belt += np.bincount(il['body'], weights=r['inland_area'], minlength=area.size)
        share = np.where(real, belt / np.maximum(area, 1e-30), 0.0)
        disc = real & (np.maximum(top, 0.0) + e * radius ** 2 / (4.0 * case['transmissivity_m2_s'])
                       <= case['reach_m'])
        rows = {}
        for region, m in deserts_.items():
            wooded = m & (share >= WOODED)
            sizes = {}
            for lo, hi in zip(ISLAND_SIZES_KM2[:-1], ISLAND_SIZES_KM2[1:]):
                s_ = m & (area >= lo) & (area < hi)
                sizes[f'{lo:g}-{hi:g}_km2'] = dict(count=int(s_.sum()), land_km2=float(area[s_].sum()),
                                                  wooded_count=int((s_ & wooded).sum()),
                                                  disc_rule_count=int((s_ & disc).sum()),
                                                  both_count=int((s_ & wooded & disc).sum()),
                                                  within_reach_share=float(belt[s_].sum() / area[s_].sum())
                                                  if s_.any() else None)
            rows[region] = dict(count=int(m.sum()), land_km2=float(area[m].sum()),
                                within_reach_share=float(belt[m].sum() / area[m].sum()) if m.any() else None,
                                wooded_count=int(wooded.sum()), wooded_land_km2=float(area[wooded].sum()),
                                largest_wooded_km2=float(area[wooded].max()) if wooded.any() else None,
                                disc_rule_count=int((m & disc).sum()), by_size=sizes)
        out['cases'][name] = rows
    return out


def closed_lakes(px, sh, uptake, cs, masks, grid, g, f, runs):
    """The rain-fed closed lakes' shores in each mask's deserts: the belts' draw against the lake's own evaporation,
    its area times the GCM's open-water evaporation at 301 K along its row, which is the row's sea cells' under
    their marine air; with the lakes' sizes and the rows' evaporation."""
    evap = interpolate(-g['design_evap_mean'], -g['test_evap_mean'], f) * DAY_S * 1000.0       # mm a day
    water = ~grid['lsm']
    row_evap = np.array([float((grid['area'][r] * water[r] * evap[r]).sum()
                               / max((grid['area'][r] * water[r]).sum(), 1e-30)) if water[r].any() else np.nan
                         for r in range(32)])
    for r in range(32):
        if np.isnan(row_evap[r]):
            good = np.flatnonzero(~np.isnan(row_evap))
            row_evap[r] = row_evap[good[np.argmin(np.abs(good - r))]]
    lake_lat = sh['lake_lat']
    lake_row = np.array([int(np.argmin(np.abs(grid['lat'] - x))) for x in lake_lat])
    lake_evap_m3_s = sh['lake_area_km2'] * 1e6 * row_evap[lake_row] / 1000.0 / DAY_S
    lake_pixels = sh['lake_area_km2'] / pixel_area_at(sh, lake_lat)
    out = {}
    co, il = px['coast'], px['inland']
    up = uptake.ravel()
    for mask_name in ('rain_mask', 'log_rain_mask'):
        desert = masks[(mask_name, 'every_cell')].ravel()
        rows = {}
        for name in cs:
            r = runs[name]
            draw = np.zeros(sh['lake_area_km2'].size)
            sel = desert[co['cell']]
            np.add.at(draw, co['lake'][sel], r['coast_areas'][sel][:, 3] * 1e6 * up[co['cell'][sel]])
            isel = desert[il['cell']] & (il['water'] == S.LAKE_CLOSED)
            np.add.at(draw, il['lake'][isel], r['inland_area'][isel] * 1e6 * up[il['cell'][isel]])
            lakes = np.flatnonzero(draw > 0)
            lakes = lakes[lakes > 0]
            ratio = draw[lakes] / np.maximum(lake_evap_m3_s[lakes], 1e-30)
            rows[name] = dict(
                lakes=int(lakes.size), draw_km3_yr=float(draw[lakes].sum() * YEAR_S / 1e9),
                evaporation_km3_yr=float(lake_evap_m3_s[lakes].sum() * YEAR_S / 1e9),
                draw_over_evaporation=dict(zip(('p10', 'median', 'p90'),
                                               [float(np.percentile(ratio, q)) for q in (10, 50, 90)])),
                lake_km2=dict(zip(('p10', 'median', 'p90'),
                                  [float(np.percentile(sh['lake_area_km2'][lakes], q)) for q in (10, 50, 90)])),
                lake_pixels_median=float(np.median(lake_pixels[lakes])),
                row_evaporation_mm_day=[float(row_evap[lake_row[lakes]].min()), float(row_evap[lake_row[lakes]].max())],
                lake_latitudes_deg=[float(np.abs(lake_lat[lakes]).min()), float(np.abs(lake_lat[lakes]).max())])
        out[mask_name] = rows
    return out


def rounded(x, digits=6):
    if isinstance(x, dict):
        return {k: rounded(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rounded(v, digits) for v in x]
    if isinstance(x, (float, np.floating)):
        if not math.isfinite(float(x)):
            return None
        return float(f'{float(x):.{digits}g}')
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def main(argv=None) -> int:
    product = rounded(results())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    d = product['deserts']
    print('deserts at 301 K:', {k: {c: d[k][c]['share'] for c in ('every_cell', 'model_land')}
                                for k in ('rain_mask', 'log_rain_mask', 'share_matched')})
    return 0


if __name__ == '__main__':
    sys.exit(main())
