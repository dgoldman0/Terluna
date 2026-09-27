"""What a ring shows, on the equator or a circle of latitude: weather through the lunar day, storms, winds and
clouds, from CM1 output.

    climate/gcm/.venv/bin/python -m climate.crm.ring_analysis ring --from-day 10

writes the summary to ../results/crm/ring_<case>.json, the fields a page or figure needs to
../results/crm/ring_<case>_fields.npz, and the full arrays to products/ring_<case>.npz.

CM1 writes a snapshot every three model hours (GrADS files: two-dimensional surface fields along x and
x-z sections). Every column of the ring sits at its own local time, the hour angle growing eastward
360 degrees around the ring, so each snapshot holds the whole lunar day at once and composites by local
time use every column. The analysis starts after a spin-up (--from-day) and reports, for land and water
separately: the air at 2 m and the ground by local time, rain, cloud cover and cloud heights, winds at
10 m and at heights that matter for flight, the depth of the mixed layer, the cloud base, CAPE; storms
as contiguous stretches of rain; and the same quantities from the GCM design case for comparison.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE / 'runs'
PRODUCTS = HERE / 'products'                                           # full arrays, regenerated, not in Git
RESULTS = HERE.parent / 'results' / 'crm'                              # compact summaries, kept in Git
SCHEMA = 'terluna.climate.crm-ring/1'
HOUR_BINS = 36                        # 10 degrees of hour angle, 0.82 Earth days
CLOUD_KG_KG = 1.0e-5                  # condensate (liquid + ice) that counts as cloud
RAIN_MM_H = 1.0                       # rain rate that counts as a storm's footprint
WIND_HEIGHTS_M = (1000.0, 5000.0, 10000.0, 20000.0, 30000.0, 40000.0, 55000.0, 70000.0)
FLIGHT_BAND_M = (35000.0, 45000.0)    # where low-gravity flight is expected to concentrate
TRACK_REACH = 5                       # columns (30 km) a storm may shift between snapshots and stay the same storm
STORM_HALF_WIDTH = 120                # columns each side of the heaviest storm kept as a section (1,440 km in all)
FIELDS_SCHEMA = 'terluna.climate.crm-ring-fields/1'
FIELDS_READING_RULE = ('Sections are means over every column by local time (10-degree bins, 0 = noon) and height, to '
                       '90 km; eastward wind is toward later local times. The storm is the snapshot with the heaviest '
                       'rain of the span: 240 columns (about 1,450 km) centred on its core, heights to 100 km, instantaneous. '
                       'The rain map keeps the heaviest rain in each block of four columns (24 km), every 3 hours, with '
                       'the position of local noon.')
EVIDENCE = ('CM1 r22.0 in two dimensions along the equator at lunar gravity (6-km columns, 111 levels to 150 km), '
            'with the Morrison microphysics given lunar fall speeds, RRTMG radiation for the design air, the Sun '
            'crossing the ring once a lunar day, the 28% scenario\'s seas and lakes laid flat at sea level with a '
            'placeholder land surface, and the domain-mean temperature, vapour and eastward wind above 8-16 km held '
            'to the GCM design case (A28_dim5, the 5% dimmer shield). It resolves storms and the day-night '
            'circulation along the equator; it has no north-south dimension, so air cannot converge on the storms '
            'from the north and south, and two-dimensional convection organises into lines more readily than '
            'real storms do.')


def evidence(latitude: float) -> str:
    """What a ring on the given circle of latitude is: EVIDENCE on the equator."""
    if latitude == 0.0:
        return EVIDENCE
    return (f'CM1 r22.0 in two dimensions along the circle of latitude {abs(latitude):.0f} {"N" if latitude > 0 else "S"} at '
            'lunar gravity (6-km columns, 111 levels to 150 km), with the Morrison microphysics given lunar fall speeds, '
            'RRTMG radiation for the design air, the Sun crossing the ring once a lunar day and rising at most '
            f'{90 - abs(latitude):.0f} degrees (equinox sunlight), the Coriolis force of the latitude on departures from '
            'the reference wind, the 28% scenario\'s seas and lakes laid flat at sea level with a placeholder land '
            'surface that evaporates as readily as the GCM\'s land on that circle, and the domain-mean temperature, '
            'vapour and winds above 8-16 km held to the GCM design case\'s row nearest the latitude (A28_dim5, the 5% '
            'dimmer shield). It has no north-south dimension, so air cannot cross the pole or converge on the storms '
            'from the north and south; the land\'s wetness is prescribed, so rain does not wet it; two-dimensional '
            'convection organises into lines more readily than real storms do.')


READING_RULE = ('Composites are by local hour angle (0 = noon, negative = morning) in 10-degree bins; heights are '
                'above the flat ground at sea level. Snapshots are instantaneous, every 3 model hours. Storm '
                'widths are contiguous stretches of surface rain above 1 mm/h in one snapshot.')


def grads_layout(ctl: Path):
    lines = ctl.read_text().splitlines()
    nx = int(next(l for l in lines if l.startswith('xdef')).split()[1])
    i = next(k for k, l in enumerate(lines) if l.lower().startswith('vars'))
    return nx, [(l.split()[0], max(1, int(l.split()[1]))) for l in lines[i + 1:i + 1 + int(lines[i].split()[1])]]


def read_snapshot(case: Path, n: int, kind: str = 's') -> dict:
    """One output time: surface fields as (nx,) arrays, sections as (nz, nx)."""
    nx, layout = grads_layout(case / f'cm1out_{kind}.ctl')
    raw = np.fromfile(case / f'cm1out_t{n:06d}_{kind}.dat', dtype='<f4')
    out, offset = {}, 0
    for name, levels in layout:
        size = levels * nx
        field = raw[offset:offset + size].reshape(levels, nx)
        out[name] = field[0] if levels == 1 else field
        offset += size
    if offset != raw.size:
        raise ValueError(f'{kind} output {n}: {raw.size} values, layout expects {offset}')
    return out


def case_geometry(case: Path) -> dict:
    record = json.loads((case / 'case.json').read_text())
    grid = record['grid']
    nx, dx = grid['nx'], grid['dx_m']
    x = (np.arange(nx) + 0.5) * dx
    zw = np.loadtxt(case / 'input_grid_z')
    land = np.zeros(nx, bool)
    rows = np.loadtxt(case / 'terluna_surface.txt', skiprows=1, ndmin=2)
    for x0, x1, xland, *_ in rows:
        land[(x >= x0) & (x < x1)] = xland == 1
    sun = dict(h0=record['configuration']['start_hour_angle_deg'], length=grid['length_m'])
    return dict(record=record, x=x, dx=dx, zh=0.5 * (zw[1:] + zw[:-1]), zw=zw, land=land, sun=sun)


def hour_angle(geo: dict, t_s: float, day_s: float) -> np.ndarray:
    """Local hour angle (degrees, -180..180, 0 = noon) of every column at model time t."""
    h = geo['sun']['h0'] + 360.0 * t_s / day_s + 360.0 * geo['x'] / geo['sun']['length']
    return (h + 180.0) % 360.0 - 180.0


def level_index(zh, z):
    return int(np.argmin(np.abs(zh - z)))


def saturation_pa(t):
    """Saturation vapour pressure over water (Pa) at temperature t (K), Bolton (1980)."""
    return 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))


def vapour_pa(r, p):
    """Vapour pressure (Pa) of air with water-vapour mixing ratio r (kg/kg) at pressure p, as CM1 writes it."""
    return r * p / (0.622 + r)


def wet_bulb(t, e, p):
    """Isobaric wet-bulb temperature (K) of air at t (K) with vapour pressure e (Pa) at pressure p, from the
    psychrometric equation by bisection: the same method as the GCM's wet-bulb figures in ../gcm/README.md."""
    t, e, p = (np.asarray(a, dtype=float) for a in (t, e, p))
    gamma = 1004.0 / (0.622 * 2.45e6)
    lo, hi = np.full_like(t, 200.0), t.copy()
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        above = saturation_pa(mid) - gamma * p * (t - mid) - e > 0
        hi, lo = np.where(above, mid, hi), np.where(above, lo, mid)
    return 0.5 * (lo + hi)


def storms(prate_mm_h, threshold=RAIN_MM_H):
    """Contiguous stretches of rain above threshold on the periodic ring: (start index, length in columns)."""
    wet = prate_mm_h >= threshold
    if wet.all():
        return [(0, wet.size)]
    shift = int(np.argmin(wet))                                        # start the scan in a dry column
    rolled = np.roll(wet, -shift)
    found, start = [], None
    for i, w in enumerate(np.append(rolled, False)):
        if w and start is None:
            start = i
        elif not w and start is not None:
            found.append(((start + shift) % wet.size, i - start))
            start = None
    return found


def analyse(name: str, from_day: float) -> dict:
    case = RUNS / name
    geo = case_geometry(case)
    record = geo['record']
    tap = record['configuration']['output_s']
    from climate.crm.cm1_run import solar_day_s
    day_s = solar_day_s()
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    use = [n for n in outputs if (n - 1) * tap >= from_day * 86400.0]
    if not use:
        raise RuntimeError(f'no output after day {from_day}')
    zh, land = geo['zh'], geo['land']
    latitude = record.get('latitude_deg', 0.0)
    has_v = any(key == 'vinterp' for key, _ in grads_layout(case / 'cm1out_s.ctl')[1])   # rings off the equator
    bins = np.linspace(-180.0, 180.0, HOUR_BINS + 1)
    surface_names = ('t2', 'tsk', 'rh2', 'tw2', 'prate', 'evap', 's10', 'u10', 'hpbl', 'lcl', 'cape', 'lwp', 'cwp', 'pwat')
    comp = {k: np.zeros((2, HOUR_BINS)) for k in (*surface_names, 'cloud', 'cloud_high', 'cloud_low', 'cloud_flight',
                                                    'cloud_top', 'cloud_base', 'wind_levels')}
    comp['wind_levels'] = np.zeros((2, HOUR_BINS, len(WIND_HEIGHTS_M)))
    count = np.zeros((2, HOUR_BINS))
    cloud_top_count = np.zeros((2, HOUR_BINS))
    cloud_profile = np.zeros((2, 2, zh.size))                          # land/water, day/night, height
    section = {k: np.zeros((HOUR_BINS, zh.size)) for k in ('u', 'w', 'cloud', 'th_anomaly', 'qv_anomaly', *(('v',) if has_v else ()))}
    section_count = np.zeros(HOUR_BINS)                                # all columns, by local time and height
    profile_count = np.zeros((2, 2))
    wind_samples = {z: [] for z in WIND_HEIGHTS_M}
    vertical_samples = {z: [] for z in WIND_HEIGHTS_M}
    snapshot_storms = []
    wind10 = {'land': [], 'water': []}
    t2_all = {'land': [], 'water': []}
    tw2_all = {'land': [], 'water': []}
    rh2_all = {'land': [], 'water': []}
    rain_all = {'land': [], 'water': []}
    pw = []
    heaviest = dict(rain=-1.0)
    storm_rows = []
    hov = dict(time_day=[], prate=[], cloud=[], t2=[], u10=[])
    ilev = {z: level_index(zh, z) for z in WIND_HEIGHTS_M}
    high = zh >= 25000.0                                              # above the freezing level
    flight = (zh >= FLIGHT_BAND_M[0]) & (zh <= FLIGHT_BAND_M[1])
    for n in use:
        t = (n - 1) * tap
        d = read_snapshot(case, n)
        e2 = vapour_pa(d['q2'], d['psfc'])
        d['rh2'], d['tw2'] = e2 / saturation_pa(d['t2']), wet_bulb(d['t2'], e2, d['psfc'])
        d['evap'] = d['qvflux'] * air_density(d) * 3600.0              # mm/h
        h = hour_angle(geo, t, day_s)
        k = np.clip(np.digitize(h, bins) - 1, 0, HOUR_BINS - 1)
        condensate = d['qc'] + d['qi']
        cloudy = condensate >= CLOUD_KG_KG                              # (nz, nx)
        any_cloud = cloudy.any(axis=0)
        top = np.where(any_cloud, zh[np.where(cloudy, np.arange(zh.size)[:, None], -1).max(axis=0)], np.nan)
        base_i = np.where(cloudy, np.arange(zh.size)[:, None], zh.size).min(axis=0)
        base = np.where(any_cloud, zh[np.minimum(base_i, zh.size - 1)], np.nan)
        prate = d['prate'] * 3600.0                                    # mm/h
        pw.append((t / 86400.0, float(d['pwat'].mean() * 1000.0)))
        if prate.max() > heaviest['rain']:
            i0 = int(np.argmax(prate))
            cols = (i0 + np.arange(-STORM_HALF_WIDTH, STORM_HALF_WIDTH)) % prate.size
            heaviest = dict(rain=float(prate.max()), day=t / 86400.0, hour_angle=float(h[i0]),
                            condensate=(d['qc'] + d['qi'] + d['qs'] + d['qg'])[:, cols] * 1000.0,
                            rain_water=d['qr'][:, cols] * 1000.0, updraft=d['winterp'][:, cols],
                            rain_mm_h=prate[cols], wind_10m=d['s10'][cols], air_2m=d['t2'][cols] - 273.15, land=land[cols])
        u = d['uinterp']
        speed = np.hypot(u, d['vinterp']) if has_v else np.abs(u)          # horizontal wind speed
        th_anomaly = d['th'] - d['th'].mean(axis=1, keepdims=True)
        qv_anomaly = d['qv'] - d['qv'].mean(axis=1, keepdims=True)
        for b in range(HOUR_BINS):
            m = k == b
            if m.any():
                section_count[b] += m.sum()
                for key, field in (('u', u), ('w', d['winterp']), ('cloud', cloudy), ('th_anomaly', th_anomaly),
                                   ('qv_anomaly', qv_anomaly), *((('v', d['vinterp']),) if has_v else ())):
                    section[key][b] += field[:, m].sum(axis=1)
        for s, sel in ((0, land), (1, ~land)):
            for b in range(HOUR_BINS):
                m = sel & (k == b)
                if not m.any():
                    continue
                count[s, b] += m.sum()
                for key in surface_names:
                    value = prate if key == 'prate' else d[key]
                    comp[key][s, b] += value[m].sum()
                comp['cloud'][s, b] += any_cloud[m].sum()
                comp['cloud_high'][s, b] += cloudy[high][:, m].any(axis=0).sum()
                comp['cloud_low'][s, b] += cloudy[~high][:, m].any(axis=0).sum()
                comp['cloud_flight'][s, b] += cloudy[flight][:, m].any(axis=0).sum()
                ok = m & any_cloud
                cloud_top_count[s, b] += ok.sum()
                comp['cloud_top'][s, b] += np.nansum(top[ok])
                comp['cloud_base'][s, b] += np.nansum(base[ok])
                for j, z in enumerate(WIND_HEIGHTS_M):
                    comp['wind_levels'][s, b, j] += speed[ilev[z], m].sum()
            for dn, sun in ((0, np.cos(np.radians(h)) > 0.0), (1, np.cos(np.radians(h)) <= 0.0)):
                m = sel & sun
                if m.any():
                    cloud_profile[s, dn] += cloudy[:, m].sum(axis=1)
                    profile_count[s, dn] += m.sum()
        for z in WIND_HEIGHTS_M:
            wind_samples[z].append(speed[ilev[z]])
            vertical_samples[z].append(d['winterp'][ilev[z]])
        for key, sel in (('land', land), ('water', ~land)):
            wind10[key].append(d['s10'][sel])
            t2_all[key].append(d['t2'][sel])
            tw2_all[key].append(d['tw2'][sel])
            rh2_all[key].append(d['rh2'][sel])
            rain_all[key].append(prate[sel])
        found = storms(prate)
        snapshot_storms.append(found)
        for start, length in found:
            idx = (start + np.arange(length)) % prate.size
            near = (start + np.arange(-17, length + 17)) % prate.size          # within about 100 km
            col_top = np.nanmax(top[idx]) if np.isfinite(top[idx]).any() else np.nan
            storm_rows.append(dict(time_day=t / 86400.0, width_km=length * geo['dx'] / 1000.0,
                                   rain_max_mm_h=float(prate[idx].max()), rain_mean_mm_h=float(prate[idx].mean()),
                                   cloud_top_km=float(col_top / 1000.0) if np.isfinite(col_top) else None,
                                   updraft_max_m_s=float(d['winterp'][:, idx].max()),
                                   gust_10m_max_m_s=float(d['s10'][near].max()),
                                   t2_drop_k=float(np.median(d['t2'][near]) - d['t2'][near].min()),
                                   land_share=float(land[idx].mean()), hour_angle_deg=float(h[idx[length // 2]])))
        hov['time_day'].append(t / 86400.0)
        hov['prate'].append(prate.astype(np.float32))
        hov['cloud'].append(any_cloud.astype(np.int8))
        hov['t2'].append(d['t2'].astype(np.float32))
        hov['u10'].append(d['u10'].astype(np.float32))
    safe = np.maximum(count, 1)
    composites = {key: (comp[key] / safe) for key in (*surface_names, 'cloud', 'cloud_high', 'cloud_low', 'cloud_flight')}
    composites['cloud_top_m'] = comp['cloud_top'] / np.maximum(cloud_top_count, 1)
    composites['cloud_base_m'] = comp['cloud_base'] / np.maximum(cloud_top_count, 1)
    composites['wind_levels'] = comp['wind_levels'] / safe[:, :, None]
    profiles = cloud_profile / np.maximum(profile_count, 1)[:, :, None]
    section = {key: v / np.maximum(section_count, 1)[:, None] for key, v in section.items()}
    pct = lambda a, q: float(np.percentile(np.concatenate(a), q))
    summary = dict(
        schema=SCHEMA, case=name, evidence=evidence(latitude), reading_rule=READING_RULE, latitude_deg=latitude,
        span_days=[(use[0] - 1) * tap / 86400.0, (use[-1] - 1) * tap / 86400.0], snapshots=len(use),
        land_share=float(land.mean()),
        air_2m_c={key: dict(p1=pct(v, 1) - 273.15, median=pct(v, 50) - 273.15, p99=pct(v, 99) - 273.15)
                  for key, v in t2_all.items()},
        wet_bulb_2m_c={key: dict(median=pct(v, 50) - 273.15, p99=pct(v, 99) - 273.15, max=float(np.concatenate(v).max()) - 273.15)
                       for key, v in tw2_all.items()},
        humidity_2m={key: dict(p10=pct(v, 10), median=pct(v, 50), p90=pct(v, 90)) for key, v in rh2_all.items()},
        wind_10m_m_s={key: dict(median=pct(v, 50), p90=pct(v, 90), p99=pct(v, 99), max=float(np.concatenate(v).max()))
                      for key, v in wind10.items()},
        wind_aloft_m_s={f'{z / 1000:.0f}_km': dict(median=pct(v, 50), p90=pct(v, 90), p99=pct(v, 99),
                                                  max=float(np.concatenate(v).max())) for z, v in wind_samples.items()},
        rain_mm_day={key: float(np.concatenate(v).mean() * 24.0) for key, v in rain_all.items()},
        water_budget=water_budget(composites, count, land, pw),
        rain_rate_mm_h={key: dict(wet_share=float((np.concatenate(v) > 0.1).mean()), p99=pct(v, 99),
                                  max=float(np.concatenate(v).max())) for key, v in rain_all.items()},
        cloud_cover={key: float(composites['cloud'][s].mean()) for s, key in ((0, 'land'), (1, 'water'))},
        vertical_wind_m_s={f'{z / 1000:.0f}_km': dict(abs_median=float(np.median(np.abs(np.concatenate(v)))),
                                                     abs_p99=float(np.percentile(np.abs(np.concatenate(v)), 99)),
                                                     up_max=float(np.concatenate(v).max()), down_max=float(-np.concatenate(v).min()))
                           for z, v in vertical_samples.items()},
        flight_band_km=[FLIGHT_BAND_M[0] / 1000.0, FLIGHT_BAND_M[1] / 1000.0],
        cloud_in_flight_band={key: float(composites['cloud_flight'][s].mean()) for s, key in ((0, 'land'), (1, 'water'))},
        storms=storm_summary(storm_rows, len(use)),
        storm_tracks=track_summary(track_storms(snapshot_storms, land.size), tap, geo['dx'], land.size),
        rain_spells=rain_spells(np.array(hov['prate']), land, tap),
        by_hour_angle=hour_angle_table(0.5 * (bins[1:] + bins[:-1]), composites),
        circulation=circulation_table(0.5 * (bins[1:] + bins[:-1]), zh, section),
        gcm=gcm_equator_composites(band=record['reference']['band_deg'], latitude=latitude))
    arrays = dict(hour_angle_deg=0.5 * (bins[1:] + bins[:-1]), heights_m=zh, wind_heights_m=np.array(WIND_HEIGHTS_M),
                  cloud_profile_land_day=profiles[0, 0], cloud_profile_land_night=profiles[0, 1],
                  cloud_profile_water_day=profiles[1, 0], cloud_profile_water_night=profiles[1, 1],
                  hov_time_day=np.array(hov['time_day']), hov_prate_mm_h=np.array(hov['prate']),
                  hov_cloud=np.array(hov['cloud']), hov_t2_k=np.array(hov['t2']), hov_u10_m_s=np.array(hov['u10']),
                  x_km=geo['x'] / 1000.0, land=land,
                  hov_noon_x_km=np.array([noon_x_km(geo, day * 86400.0, day_s) for day in hov['time_day']]),
                  storm_day=heaviest['day'], storm_hour_angle_deg=heaviest['hour_angle'], storm_rain_max_mm_h=heaviest['rain'],
                  storm_condensate_g_kg=heaviest['condensate'], storm_rain_water_g_kg=heaviest['rain_water'],
                  storm_updraft_m_s=heaviest['updraft'], storm_rain_mm_h=heaviest['rain_mm_h'],
                  storm_wind_10m_m_s=heaviest['wind_10m'], storm_air_2m_c=heaviest['air_2m'], storm_land=heaviest['land'],
                  **{f'section_{key}': v for key, v in section.items()},
                  **{f'{key}_land': v[0] for key, v in composites.items()},
                  **{f'{key}_water': v[1] for key, v in composites.items()})
    return summary, arrays


def noon_x_km(geo: dict, t_s: float, day_s: float) -> float:
    """Where on the ring it is local noon at model time t (km east of the ring's origin)."""
    return float((-(geo['sun']['h0'] + 360.0 * t_s / day_s) / 360.0 * geo['sun']['length']) % geo['sun']['length'] / 1000.0)


def write_fields(path: Path, summary: dict, arrays: dict) -> None:
    """The fields a page or figure needs, small enough to keep in Git: the circulation by local time and
    height, the heaviest storm's section, and the rain map pooled to 24-km blocks."""
    zh = arrays['heights_m'] / 1000.0
    sec, low = zh <= 90.0, zh <= 100.0
    rain = arrays['hov_prate_mm_h']
    blocks = rain.shape[1] // 4
    x = arrays['x_km'][:blocks * 4].reshape(blocks, 4).mean(axis=1)
    meta = dict(schema=FIELDS_SCHEMA, case=summary['case'], span_days=summary['span_days'], evidence=summary['evidence'],
                reading_rule=FIELDS_READING_RULE, producer=summary.get('producer'))
    np.savez_compressed(
        path, metadata=json.dumps(meta),
        section_hour_angle_deg=arrays['hour_angle_deg'].astype(np.float32), section_height_km=zh[sec].astype(np.float32),
        section_eastward_m_s=arrays['section_u'][:, sec].astype(np.float32),
        section_upward_m_s=arrays['section_w'][:, sec].astype(np.float32),
        section_cloud_fraction=arrays['section_cloud'][:, sec].astype(np.float32),
        **({'section_northward_m_s': arrays['section_v'][:, sec].astype(np.float32)} if 'section_v' in arrays else {}),
        storm_day=np.float32(arrays['storm_day']), storm_hour_angle_deg=np.float32(arrays['storm_hour_angle_deg']),
        storm_rain_max_mm_h=np.float32(arrays['storm_rain_max_mm_h']), storm_height_km=zh[low].astype(np.float32),
        storm_x_km=(np.arange(2 * STORM_HALF_WIDTH) * (arrays['x_km'][1] - arrays['x_km'][0])).astype(np.float32),
        **{k: arrays[k][low].astype(np.float16) for k in ('storm_condensate_g_kg', 'storm_rain_water_g_kg', 'storm_updraft_m_s')},
        **{k: arrays[k].astype(np.float32) for k in ('storm_rain_mm_h', 'storm_wind_10m_m_s', 'storm_air_2m_c')},
        storm_land=arrays['storm_land'].astype(bool),
        rain_map_time_day=arrays['hov_time_day'].astype(np.float32), rain_map_x_km=x.astype(np.float32),
        rain_map_land_share=arrays['land'][:blocks * 4].reshape(blocks, 4).mean(axis=1).astype(np.float32),
        rain_map_max_mm_h=rain[:, :blocks * 4].reshape(rain.shape[0], blocks, 4).max(axis=2).astype(np.float16),
        rain_map_noon_x_km=arrays['hov_noon_x_km'].astype(np.float32))


def air_density(d):
    """Air density (kg/m3) in the lowest model level of a snapshot."""
    t = d['th'][0] * (d['prs'][0] / 1.0e5) ** (287.04 / 1005.7)
    return d['prs'][0] / (287.04 * t * (1.0 + 0.61 * d['qv'][0]))


def water_budget(composites, count, land, pw) -> dict:
    """Evaporation against rain over the analysed span, and how fast the air's water was still changing.
    Rain below evaporation with water still building up means the lower air had not yet come into balance."""
    weight = count / count.sum()
    mean = lambda key: float((composites[key] * weight).sum() * 24.0)
    by_surface = lambda key, s: float((composites[key][s] * count[s]).sum() / count[s].sum() * 24.0)
    (t0, w0), (t1, w1) = pw[0], pw[-1]
    return dict(evaporation_mm_day=dict(all=mean('evap'), land=by_surface('evap', 0), water=by_surface('evap', 1)),
                rain_mm_day=dict(all=mean('prate'), land=by_surface('prate', 0), water=by_surface('prate', 1)),
                precipitable_water_mm=dict(start=w0, end=w1),
                storage_mm_day=(w1 - w0) / (t1 - t0) if t1 > t0 else 0.0,
                residence_days=w1 / max(mean('evap'), 1e-9))


def circulation_table(centres, zh, section) -> dict:
    """The day-night circulation at the wind heights, all columns, by local time: the mean eastward wind
    (positive toward later local times, since places to the east are further into their day) and the mean
    vertical wind."""
    rows = [level_index(zh, z) for z in WIND_HEIGHTS_M]
    out = dict(hour_angle_deg=[round(float(c), 1) for c in centres], heights_km=[round(z / 1000.0) for z in WIND_HEIGHTS_M],
               eastward_m_s=[[round(float(v), 2) for v in section['u'][:, k]] for k in rows],
               upward_cm_s=[[round(float(v) * 100.0, 2) for v in section['w'][:, k]] for k in rows],
               cloud_fraction=[[round(float(v), 3) for v in section['cloud'][:, k]] for k in rows])
    if 'v' in section:                                                 # off the equator, the wind turns
        out['northward_m_s'] = [[round(float(v), 2) for v in section['v'][:, k]] for k in rows]
    return out


def hour_angle_table(centres, composites) -> dict:
    """The composites by local time as rounded lists, land and water, for the summary kept in Git."""
    fields = (('air_2m_c', 't2', 1.0, -273.15, 1), ('ground_c', 'tsk', 1.0, -273.15, 1), ('humidity_2m', 'rh2', 1.0, 0.0, 2),
              ('wet_bulb_2m_c', 'tw2', 1.0, -273.15, 1), ('rain_mm_h', 'prate', 1.0, 0.0, 3),
              ('evaporation_mm_h', 'evap', 1.0, 0.0, 3), ('wind_10m_m_s', 's10', 1.0, 0.0, 1),
              ('mixed_layer_km', 'hpbl', 1e-3, 0.0, 1),
              ('cloud_base_lcl_km', 'lcl', 1e-3, 0.0, 1), ('cape_j_kg', 'cape', 1.0, 0.0, 0),
              ('cloud_cover', 'cloud', 1.0, 0.0, 2), ('cloud_above_25_km', 'cloud_high', 1.0, 0.0, 2),
              ('cloud_in_flight_band', 'cloud_flight', 1.0, 0.0, 2),
              ('cloud_top_km', 'cloud_top_m', 1e-3, 0.0, 1), ('precipitable_water_mm', 'pwat', 1e3, 0.0, 1))
    out = dict(hour_angle_deg=[round(float(c), 1) for c in centres],
               wind_heights_km=[round(z / 1000.0) for z in WIND_HEIGHTS_M])
    for s, surface in ((0, 'land'), (1, 'water')):
        rows = {name: [round(float(v) * scale + offset, digits) for v in composites[key][s]]
                for name, key, scale, offset, digits in fields}
        rows['wind_by_height_m_s'] = [[round(float(v), 1) for v in level] for level in composites['wind_levels'][s].T]
        out[surface] = rows
    return out


def track_storms(snapshot_storms, nx, reach=TRACK_REACH):
    """Follow storms through consecutive snapshots. A storm continues the track whose latest footprint most
    overlaps its own footprint widened by `reach` columns each side; when a storm splits, the rest start new
    tracks. Returns tracks as lists of (snapshot index, start column, width in columns)."""
    tracks, current = [], []
    for ti, found in enumerate(snapshot_storms):
        following, taken = [], set()
        for start, length in found:
            reach_cols = set(((start + np.arange(-reach, length + reach)) % nx).tolist())
            best, best_overlap = None, 0
            for j, track in enumerate(current):
                if j in taken:
                    continue
                _, s0, l0 = track[-1]
                overlap = len(reach_cols & set(((s0 + np.arange(l0)) % nx).tolist()))
                if overlap > best_overlap:
                    best, best_overlap = j, overlap
            if best is None:
                track = [(ti, start, length)]
                tracks.append(track)
            else:
                track = current[best]
                track.append((ti, start, length))
                taken.add(best)
            following.append(track)
        current = following
    return tracks


def track_summary(tracks, output_s, dx, nx) -> dict:
    """Lifetimes (the snapshots a storm was seen in, times the output interval) and drift (+ eastward)."""
    if not tracks:
        return dict(count=0)
    seen = np.array([len(t) for t in tracks])
    drift = []
    for t in tracks:
        if len(t) < 2:
            continue
        centre = lambda item: item[1] + 0.5 * item[2]
        shift = (centre(t[-1]) - centre(t[0]) + 0.5 * nx) % nx - 0.5 * nx
        drift.append(shift * dx / ((t[-1][0] - t[0][0]) * output_s))
    drift = np.array(drift)
    hours = seen * output_s / 3600.0
    return dict(count=len(tracks), output_interval_h=output_s / 3600.0,
                seen_once_share=float((seen == 1).mean()),
                lifetime_h=dict(median=float(np.median(hours)), p90=float(np.percentile(hours, 90)), max=float(hours.max())),
                drift_m_s=dict(p10=float(np.percentile(drift, 10)), median=float(np.median(drift)),
                               p90=float(np.percentile(drift, 90))) if drift.size else None)


def rain_spells(prate_mm_h, land, output_s, threshold=RAIN_MM_H) -> dict:
    """Rain as a place meets it: spells of consecutive snapshots with rain above threshold at one column,
    how many per lunar day and how long, over land and water. The count per lunar day needs a span of at
    least one lunar day, since only then has every place passed through every local time."""
    from climate.crm.cm1_run import solar_day_s
    nt = prate_mm_h.shape[0]
    lunar_days = nt * output_s / solar_day_s()
    out = dict(span_lunar_days=lunar_days)
    for key, sel in (('land', land), ('water', ~land)):
        lengths = []
        for column in prate_mm_h[:, sel].T >= threshold:
            run = 0
            for wet in np.append(column, False):
                if wet:
                    run += 1
                elif run:
                    lengths.append(run)
                    run = 0
        lengths = np.array(lengths) * output_s / 3600.0
        per_day = len(lengths) / sel.sum() / lunar_days
        out[key] = dict(spells_per_lunar_day=float(per_day) if lunar_days >= 0.95 else None,
                        duration_h=dict(median=float(np.median(lengths)), p90=float(np.percentile(lengths, 90)),
                                        max=float(lengths.max())) if lengths.size else None)
    return out


def storm_summary(rows, snapshots):
    if not rows:
        return dict(count=0)
    get = lambda k: np.array([r[k] for r in rows if r[k] is not None], dtype=float)
    stats = lambda a: dict(median=float(np.median(a)), p90=float(np.percentile(a, 90)), max=float(a.max())) if a.size else None
    return dict(count=len(rows), per_snapshot=len(rows) / snapshots, width_km=stats(get('width_km')),
                rain_max_mm_h=stats(get('rain_max_mm_h')),
                cloud_top_km=stats(get('cloud_top_km')), updraft_max_m_s=stats(get('updraft_max_m_s')),
                gust_10m_max_m_s=stats(get('gust_10m_max_m_s')), t2_drop_k=stats(get('t2_drop_k')),
                land_share=float(np.mean([r['land_share'] for r in rows])),
                hour_angle_deg=stats(get('hour_angle_deg')))


def gcm_equator_composites(run='A28_dim5', years=(15, 24), band=10.0, bins=10, latitude=0.0):
    """The GCM design case on the equator (or the rows within `band` degrees of another latitude) by local
    time, for comparison: air near the ground, ground, rain and cloud cover of sea-level seas and of land
    (3-day means, ten 36-degree bins)."""
    import netCDF4
    from climate.gcm.climatology import subsolar_longitudes
    from climate.crm.cm1_run import planet
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    get = {k: [] for k in ('tas', 'ts', 'pr', 'clt', 'czen')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            for k in get:
                get[k].append(np.asarray(d[k][:], dtype=float))
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
    f = {k: np.concatenate(v) for k, v in get.items()}
    height = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(lsm.shape) / planet()['gravity_m_s2']
    rows = (np.abs(lat - latitude) < band)[:, None]
    groups = dict(land=rows & lsm, water=rows & ~lsm & (height < 1.0))
    sub = subsolar_longitudes(f['czen'], lat, lon)
    hour = (lon[None, :] - sub[:, None] + 180.0) % 360.0 - 180.0                       # (t, lon)
    k = np.clip(((hour + 180.0) / (360.0 / bins)).astype(int), 0, bins - 1)
    out = dict(run=run, years=list(years), band_deg=band, latitude_deg=latitude, hour_angle_deg=list((np.arange(bins) + 0.5) * 360.0 / bins - 180.0))
    for name, sel in groups.items():
        rows_out = {}
        for key, scale, offset in (('tas', 1.0, -273.15), ('ts', 1.0, -273.15), ('pr', 86400e3 / 24.0, 0.0), ('clt', 1.0, 0.0)):
            vals = []
            for b in range(bins):
                m = sel[None] & (k[:, None, :] == b)
                vals.append(float(f[key][m].mean() * scale + offset))
            rows_out[{'tas': 'air_c', 'ts': 'ground_c', 'pr': 'rain_mm_h', 'clt': 'cloud_cover'}[key]] = vals
        out[name] = rows_out
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=10.0)
    args = parser.parse_args(argv)
    summary, arrays = analyse(args.case, args.from_day)
    summary['producer'] = dict(domain='climate', files={'crm/ring_analysis.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    PRODUCTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f'ring_{args.case}.json').write_text(json.dumps(summary, indent=1) + '\n')
    write_fields(RESULTS / f'ring_{args.case}_fields.npz', summary, arrays)
    np.savez_compressed(PRODUCTS / f'ring_{args.case}.npz', metadata=json.dumps({k: summary[k] for k in ('schema', 'case', 'evidence', 'reading_rule', 'span_days')}),
                        **arrays)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('evidence', 'reading_rule')}, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
