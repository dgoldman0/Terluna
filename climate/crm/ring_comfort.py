"""Comfort to be out in on the flat rings' land, by latitude and by the ground's real height, beside the GCM.

    climate/gcm/.venv/bin/python -m climate.crm.ring_comfort ring_equator ring_70_45e ring_70_135e

writes the comparison to ../results/crm/comfort_flat_rings.json; with --gcm RUN:FIRST-LAST it compares another GCM
run instead of the design case and writes comfort_flat_rings_<RUN>.json.

The rings lay their land at sea level. The author's plan of 2026-09-29 (research/decisions.md) reads them for
low-lying land and takes high ground from them corrected for height, with what the terrain twin showed before its
gales: ring_70_45e_terrain, the same ring with its ground at its real heights, minus ring_70_45e, in the air and
dewpoint at 2 m over land through model days 1-11, by day and by night, in classes of ground height. The correction
is interpolated in height between the class means, is nothing at sea level and holds the highest class's value above
it; every land column of every ring takes it at its own height from the atlas (cm1_run.path_heights, as the terrain
ring did). The GCM's own land cells are the check: how their air and lowest-layer dewpoint change with the cell's
height within bands of latitude. For each ring the analysis also reports how settled its second lunar day is, its
comfort at sea level by latitude with what keeps the rest uncomfortable, and its land evaporation against the GCM's
nearest land cell, with the share of its land that CM1's driest land class makes wetter than the GCM's ground.
"""
from __future__ import annotations
import argparse
import functools
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from climate.crm import ring_analysis as ra
from climate.crm import ring_crossings as rc

RESULTS = ra.RESULTS
SCHEMA = 'terluna.climate.crm-comfort/1'
TWIN = ('ring_70_45e', 'ring_70_45e_terrain')          # the flat ring and the same ring with its ground heights
TWIN_DAYS = (1.0, 11.0)                                 # after the first day's adjustment, before the gales
HEIGHT_CLASSES_M = ((0, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 6000), (6000, 9000))
LATITUDE_BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))
REPORT_BANDS = ((0, 30), (30, 60), (60, 90))
REPORT_HEIGHTS_M = ((0, 500), (500, 2000), (2000, 4000), (4000, 9000))
SHARES_H = (100, 200, 350)                              # comfortable hours per lunar day the land shares are counted at
EVIDENCE_WITH = ('CM1 r22.0 rings along great circles at lunar gravity, their land laid flat at sea level (each ring\'s own '
            'summary says what a ring is), read over their second lunar day after a lunar day of spin-up. The height '
            'correction comes from one ring\'s path over its first eleven days, in a two-dimensional run whose air, '
            'forced over its mountains, later broke into gales along the whole ring; it assumes high ground everywhere '
            'cools and dries with height as that path\'s did then, and it carries none of the highlands\' own weather '
            '(slope winds, rain on the slopes facing the wind, cold air pooling in valleys). The GCM run '
            '({run}, years {first}-{last}, 3-day means, cells about 170 km wide) is the other view.')
EVIDENCE = EVIDENCE_WITH.format(run=rc.GCM_RUN, first=rc.GCM_YEARS[0], last=rc.GCM_YEARS[1])
READING_RULE = ('Hours per lunar day are the share of the 3-hourly snapshots of the second lunar day (model days 29.53 to '
                '59.0) in which a column\'s air at 2 m and dewpoint fall in the band, times the lunar day, averaged over '
                'the columns: strict means 18-26 C with a dewpoint of at most 15 C, loose 17-27 C and 16 C. Corrected '
                'values add the twin\'s difference at the column\'s height, the day\'s or the night\'s as the Sun is up '
                'or down at the column; flat values are the ring as run. The GCM\'s hours come from its 3-day means at '
                'the nearest land cell, its air carried to the ground along a dry adiabat and its dewpoint that of its '
                'lowest layer, about 0.86 km up; the cell stands at its own height, given beside. Settling is the mean '
                'over the last third of the second lunar day minus that over the same third of the first, the same '
                'local times at every column. Evaporation is over the second lunar day; the GCM\'s is its mean over '
                'its years at the nearest land cell, whose wetness the ring took.')


def surface_reader(case: Path, sigma: float | None = None):
    """A reader of one snapshot of case: every column's air and dewpoint at 2 m (C), whether the Sun is up and the
    evaporation (mm/day), and with sigma the dewpoint where the pressure is sigma times the surface's, reading only
    those fields and the sections' lowest levels (to 2.5 km)."""
    nx, layout = ra.grads_layout(case / 'cm1out_s.ctl')
    offsets, offset = {}, 0
    for name, levels in layout:
        offsets[name] = offset
        offset += levels * nx
    levels = int(np.searchsorted(ra.case_geometry(case)['zh'], 2500.0)) + 1 if sigma else 1

    def read(n: int) -> dict:
        raw = np.memmap(case / f'cm1out_t{n:06d}_s.dat', dtype='<f4', mode='r')
        f = {k: np.array(raw[offsets[k]:offsets[k] + nx], dtype=float) for k in ('t2', 'q2', 'psfc', 'coszen', 'qvflux', 'prate')}
        section = {k: np.array(raw[offsets[k]:offsets[k] + levels * nx], dtype=float).reshape(levels, nx)
                   for k in ('th', 'prs', 'qv')}                               # a section's lowest level first
        out = dict(air_c=f['t2'] - 273.15, dew_c=rc.dewpoint_c(ra.vapour_pa(f['q2'], f['psfc'])),
                   day=f['coszen'] > 0.0, evap_mm_day=f['qvflux'] * ra.air_density(section) * 86400.0,
                   rain_mm_day=f['prate'] * 86400.0)
        if sigma:
            layer = rc.at_sigma(section['qv'], section['prs'], f['psfc'], sigma)
            out['layer_dew_c'] = rc.dewpoint_c(ra.vapour_pa(layer, sigma * f['psfc']))
        return out
    return read


def column_wetness(case: Path, x) -> np.ndarray:
    """The moisture availability each column of a case was given: the SLMO of its land-use row in the case's own
    table (the classes in force when the case was set up), NaN on water."""
    lu = np.zeros(len(x), int)
    for x0, x1, land, row, *_ in np.loadtxt(case / 'terluna_surface.txt', skiprows=1, ndmin=2):
        lu[(x >= x0) & (x < x1)] = int(row) if land == 1 else 0
    lines = (case / 'LANDUSE.TBL').read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == 'SUMMER') + 1
    slmo = {}
    for line in lines[start:]:
        if line.strip() == 'WINTER':
            break
        fields = line.split(',')
        slmo[int(fields[0])] = float(fields[2])                         # index, ALBD, SLMO
    return np.array([slmo.get(k, np.nan) for k in lu])


def snapshots(case: Path, output_s: float, from_s: float, to_s: float) -> list:
    """The output numbers of case whose model times lie from from_s to to_s (seconds)."""
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    return [n for n in outputs if from_s <= (n - 1) * output_s <= to_s]


def class_table(heights, land, d_air, d_dew, day, classes=HEIGHT_CLASSES_M) -> list:
    """Differences (time, column) over land in classes of ground height, by day and by night: for each class with
    land, its bounds, its column count, its mean height (m) and the mean differences in air and dewpoint (C)."""
    rows = []
    for lo, hi in classes:
        m = land & (heights >= lo) & (heights < hi)
        if not m.any():
            continue
        row = dict(from_m=lo, to_m=hi, columns=int(m.sum()), mean_height_m=float(heights[m].mean()))
        for part, sun in (('day', day), ('night', ~day)):
            sel = sun & m[None, :]
            row[f'air_{part}_c'] = float(d_air[sel].mean()) if sel.any() else 0.0
            row[f'dewpoint_{part}_c'] = float(d_dew[sel].mean()) if sel.any() else 0.0
        rows.append(row)
    return rows


def per_km(table: list, key: str, above_m: float = 500.0) -> float:
    """The change per kilometre of a class table's column key, fitted through zero over the classes whose mean
    height is at least above_m."""
    rows = [r for r in table if r['mean_height_m'] >= above_m]
    h = np.array([r['mean_height_m'] for r in rows])
    v = np.array([r[key] for r in rows])
    return float(np.sum(h * v) / np.sum(h * h) * 1000.0)


def correction(table: list, height_m, day) -> tuple:
    """The differences in air and dewpoint (C) at heights height_m, the day's where day is true and the night's
    elsewhere: interpolated between the class means of table, nothing at sea level, and the highest class's value
    above it."""
    x = np.array([0.0] + [r['mean_height_m'] for r in table])
    out = []
    for key in ('air', 'dewpoint'):
        by_day = np.interp(height_m, x, np.array([0.0] + [r[f'{key}_day_c'] for r in table]))
        by_night = np.interp(height_m, x, np.array([0.0] + [r[f'{key}_night_c'] for r in table]))
        out.append(np.where(day, by_day, by_night))
    return tuple(out)


def twin_differences(flat: str = TWIN[0], rough: str = TWIN[1], days=TWIN_DAYS) -> dict:
    """What the terrain twin showed over land by height class, terrain minus flat, through the given model days."""
    flat_case, rough_case = ra.RUNS / flat, ra.RUNS / rough
    geo = ra.case_geometry(rough_case)
    heights = np.fromfile(rough_case / 'perts.dat', dtype='<f4').astype(float)
    output_s = geo['record']['configuration']['output_s']
    read_flat, read_rough = surface_reader(flat_case), surface_reader(rough_case)
    d_air, d_dew, day = [], [], []
    for n in snapshots(rough_case, output_s, days[0] * 86400.0, days[1] * 86400.0):
        a, b = read_flat(n), read_rough(n)
        d_air.append(b['air_c'] - a['air_c'])
        d_dew.append(b['dew_c'] - a['dew_c'])
        day.append(a['day'])
    d_air, d_dew, day = np.array(d_air), np.array(d_dew), np.array(day)
    table = class_table(heights, geo['land'], d_air, d_dew, day)
    return dict(rings=[flat, rough], days=list(days), snapshots=len(d_air), classes=table,
                per_km_above_500_m={k: per_km(table, k) for k in ('air_day_c', 'dewpoint_day_c', 'air_night_c',
                                                                   'dewpoint_night_c')})


def gcm_land(run: str = rc.GCM_RUN, years=rc.GCM_YEARS) -> dict:
    """The GCM design case's 3-day means over its years at every cell: air at 2 m and the lowest layer's dewpoint (C),
    with its mean evaporation (mm/day), grid, land mask and ground heights (m)."""
    import netCDF4
    from climate.crm.cm1_run import planet
    folder = ra.HERE.parent / 'gcm' / 'runs' / run / 'model'
    get = {k: [] for k in ('tas', 'hus', 'ps', 'evap', 'pr')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            for k in get:
                get[k].append(np.asarray(d[k][:, -1] if k == 'hus' else d[k][:], dtype=float))   # the lowest layer's vapour
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            land = np.asarray(d['lsm'][:], float)[0] > 0.5
            sigma = float(d['lev'][-1])
    f = {k: np.concatenate(v) for k, v in get.items()}
    height = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(land.shape) / planet()['gravity_m_s2']
    q, p = f['hus'], f['ps'] * 100.0
    return dict(run=run, years=list(years), sigma=sigma, lat=lat, lon=lon, land=land, height_m=height, air_c=f['tas'] - 273.15,
                dew_c=rc.dewpoint_c(q * p / (0.622 + 0.378 * q)),
                evap_mm_day=-f['evap'].mean(axis=0) * 8.64e7,             # PlaSim counts evaporation as negative
                rain_mm_day=f['pr'].mean(axis=0) * 8.64e7)


def height_rates(gcm: dict, bands=((0, 30), (30, 60), (60, 90))) -> list:
    """Within each band of latitude, the least-squares change per kilometre of cell height of the GCM land cells'
    mean air and median lowest-layer dewpoint, and their values at sea level on that line."""
    alat = np.abs(gcm['lat'])[:, None] * np.ones(gcm['lon'].size)[None, :]
    air, dew = gcm['air_c'].mean(axis=0), np.median(gcm['dew_c'], axis=0)
    rows = []
    for lo, hi in bands:
        m = gcm['land'] & (alat >= lo) & (alat < hi)
        h = gcm['height_m'][m] / 1000.0
        a, d = np.polyfit(h, air[m], 1), np.polyfit(h, dew[m], 1)
        rows.append(dict(from_deg=lo, to_deg=hi, cells=int(m.sum()), heights_km=[float(h.min()), float(h.max())],
                         air_per_km_c=float(a[0]), air_at_sea_level_c=float(a[1]),
                         dewpoint_per_km_c=float(d[0]), dewpoint_at_sea_level_c=float(d[1])))
    return rows


def nearest_land_cells(lat_deg, lon_deg, grid_lat, grid_lon, land) -> tuple:
    """Row and column indices of the land cell nearest each point, along great circles."""
    la, lo = np.radians(np.meshgrid(grid_lat, grid_lon, indexing='ij'))
    rows, cols = [], []
    for a, b in zip(np.radians(lat_deg), np.radians(lon_deg)):
        cosd = np.where(land, np.sin(la) * np.sin(a) + np.cos(la) * np.cos(a) * np.cos(lo - b), -2.0)
        i, j = np.unravel_index(int(np.argmax(cosd)), cosd.shape)
        rows.append(i)
        cols.append(j)
    return np.array(rows), np.array(cols)


@functools.cache
def lunar_day_h() -> float:
    from climate.crm.cm1_run import solar_day_s
    return solar_day_s() / 3600.0


def hours(share):
    """A share of the lunar day in hours (a float for a single share)."""
    out = np.asarray(share, dtype=float) * lunar_day_h()
    return float(out) if out.ndim == 0 else out


def ring_columns(name: str, table: list, gcm: dict) -> dict:
    """Column by column over the second lunar day of one flat ring: latitude, ground height, land, comfortable hours
    per lunar day as run and corrected for height, what keeps the rest uncomfortable, the mean corrected air and
    dewpoint, evaporation, and the GCM's nearest land cell; with the ring's settling and its air by day and night."""
    from climate.crm.cm1_run import path_heights, solar_day_s
    case = ra.RUNS / name
    geo = ra.case_geometry(case)
    record = geo['record']
    lat, lon = np.asarray(record['path']['lat']), np.asarray(record['path']['lon'])
    land, output_s, day_s = geo['land'], record['configuration']['output_s'], solar_day_s()
    height = np.where(land, path_heights(lat, lon), 0.0)
    read = surface_reader(case, gcm['sigma'])
    numbers = snapshots(case, output_s, day_s, 2.0 * day_s)
    second = [read(n) for n in numbers]
    air, dew, layer_dew, sun = (np.array([s[k] for s in second]) for k in ('air_c', 'dew_c', 'layer_dew_c', 'day'))
    evap = np.mean([s['evap_mm_day'] for s in second], axis=0)
    rain = np.mean([s['rain_mm_day'] for s in second], axis=0)
    d_air, d_dew = correction(table, height[None, :], sun)
    air2, dew2 = air + d_air, dew + d_dew
    late = [s for n, s in zip(numbers, second) if (n - 1) * output_s >= 5.0 * day_s / 3.0]
    early = [read(n) for n in snapshots(case, output_s, 2.0 * day_s / 3.0, day_s - 1.0)]
    settling = {surface: {f'{key}_c': float(np.mean([s[f'{key}_c'][sel] for s in late])
                                              - np.mean([s[f'{key}_c'][sel] for s in early])) for key in ('air', 'dew')}
                for surface, sel in (('land', land), ('water', ~land))}
    rows, cols = nearest_land_cells(lat, lon, gcm['lat'], gcm['lon'], gcm['land'])
    g_air, g_dew = gcm['air_c'][:, rows, cols], gcm['dew_c'][:, rows, cols]
    driest = float(np.nanmin(column_wetness(case, geo['x'])[land]))       # the driest class the ring's land was given
    rows_any, cols_any = nearest_land_cells(lat, lon, gcm['lat'], gcm['lon'], np.ones_like(gcm['land']))
    cell_land, cell_height = gcm['land'][rows_any, cols_any], gcm['height_m'][rows_any, cols_any]
    against = {}
    for surface, sel in (('sea', ~land & ~cell_land & (cell_height < 1.0)), ('low_land', land & cell_land & (cell_height < 1000.0))):
        if not sel.any():
            against[surface] = dict(columns=0)
            continue
        c_air, c_dew = gcm['air_c'][:, rows_any, cols_any][:, sel], gcm['dew_c'][:, rows_any, cols_any][:, sel]
        against[surface] = dict(
            columns=int(sel.sum()), gcm_cell_height_m=float(cell_height[sel].mean()),
            cm1=dict(air_c=float(air[:, sel].mean()), dewpoint_2m_c=float(np.median(dew[:, sel])),
                     dewpoint_gcm_layer_c=float(np.median(layer_dew[:, sel])), rain_mm_day=float(rain[sel].mean()),
                     strict_hours=hours(rc.comfortable(air[:, sel], dew[:, sel], 'strict').mean())),
            gcm=dict(air_c=float(c_air.mean()), dewpoint_gcm_layer_c=float(np.median(c_dew)),
                     rain_mm_day=float(gcm['rain_mm_day'][rows_any, cols_any][sel].mean()),
                     strict_hours=hours(rc.comfortable(c_air, c_dew, 'strict').mean())))
    against['all_land_rain_mm_day'] = dict(cm1=float(rain[land].mean()), gcm_nearest_land_cell=float(gcm['rain_mm_day'][rows, cols][land].mean()))
    return dict(
        name=name, lat=np.abs(lat), height=height, land=land, wetness=np.asarray(record['reference']['wetness']), driest=driest,
        strict=rc.comfortable(air2, dew2, 'strict').mean(axis=0), loose=rc.comfortable(air2, dew2, 'loose').mean(axis=0),
        strict_flat=rc.comfortable(air, dew, 'strict').mean(axis=0), loose_flat=rc.comfortable(air, dew, 'loose').mean(axis=0),
        too_hot=(air > 26.0).mean(axis=0), too_cold=(air < 18.0).mean(axis=0),
        too_humid=((air >= 18.0) & (air <= 26.0) & (dew > 15.0)).mean(axis=0),
        air=air2.mean(axis=0), dew=dew2.mean(axis=0), evap=evap, flat_air=air, flat_dew=dew, sun=sun,
        gcm_strict=rc.comfortable(g_air, g_dew, 'strict').mean(axis=0), gcm_loose=rc.comfortable(g_air, g_dew, 'loose').mean(axis=0),
        gcm_height=gcm['height_m'][rows, cols], gcm_air=g_air.mean(axis=0), gcm_dew=np.median(g_dew, axis=0),
        gcm_evap=gcm['evap_mm_day'][rows, cols], settling=settling, against_gcm=against, snapshots=len(second),
        span_days=[(numbers[0] - 1) * output_s / 86400.0, (numbers[-1] - 1) * output_s / 86400.0])


def sea_level_bands(r: dict, bands=LATITUDE_BANDS) -> list:
    """One ring's land at sea level (as run) by band of latitude."""
    out = []
    for lo, hi in bands:
        m = r['land'] & (r['lat'] >= lo) & (r['lat'] < hi)
        if not m.any():
            continue
        day, night = r['sun'] & m[None, :], ~r['sun'] & m[None, :]
        out.append(dict(
            from_deg=lo, to_deg=hi, columns=int(m.sum()),
            comfortable_hours={'strict': hours(r['strict_flat'][m].mean()), 'loose': hours(r['loose_flat'][m].mean())},
            uncomfortable_hours={'too_hot': hours(r['too_hot'][m].mean()), 'too_cold': hours(r['too_cold'][m].mean()),
                                 'too_humid': hours(r['too_humid'][m].mean())},
            air_day_c=dict(median=float(np.median(r['flat_air'][day])), p95=float(np.percentile(r['flat_air'][day], 95))),
            air_night_c=dict(p5=float(np.percentile(r['flat_air'][night], 5)), median=float(np.median(r['flat_air'][night]))),
            dewpoint_c=dict(day=float(np.median(r['flat_dew'][day])), night=float(np.median(r['flat_dew'][night]))),
            evaporation_mm_day=dict(cm1=float(r['evap'][m].mean()), gcm=float(r['gcm_evap'][m].mean())),
            gcm_wetness_median=float(np.median(r['wetness'][m])),
            share_given_the_driest_class=float(np.mean(r['wetness'][m] < r['driest'])),
            gcm=dict(comfortable_hours={'strict': hours(r['gcm_strict'][m].mean()), 'loose': hours(r['gcm_loose'][m].mean())},
                     cell_height_km=float(r['gcm_height'][m].mean() / 1000.0))))
    return out


def by_latitude_and_height(rings: list, bands=REPORT_BANDS, heights=REPORT_HEIGHTS_M) -> list:
    """Land of all the rings together, corrected to its real height, by band of latitude and class of height."""
    cat = lambda key: np.concatenate([r[key] for r in rings])
    land, lat, height = cat('land'), cat('lat'), cat('height')
    out = []
    for lo, hi in bands:
        for h0, h1 in heights:
            m = land & (lat >= lo) & (lat < hi) & (height >= h0) & (height < h1)
            if not m.any():
                continue
            out.append(dict(
                from_deg=lo, to_deg=hi, from_m=h0, to_m=h1, columns=int(m.sum()), mean_height_km=float(height[m].mean() / 1000.0),
                comfortable_hours={'strict': hours(cat('strict')[m].mean()), 'loose': hours(cat('loose')[m].mean()),
                                   'strict_flat': hours(cat('strict_flat')[m].mean())},
                air_c=float(cat('air')[m].mean()), dewpoint_c=float(cat('dew')[m].mean()),
                gcm=dict(comfortable_hours={'strict': hours(cat('gcm_strict')[m].mean()), 'loose': hours(cat('gcm_loose')[m].mean())},
                         cell_height_km=float(cat('gcm_height')[m].mean() / 1000.0), air_c=float(cat('gcm_air')[m].mean()),
                         dewpoint_c=float(cat('gcm_dew')[m].mean()))))
    return out


def land_shares(rings: list, thresholds=SHARES_H) -> dict:
    """The share of all the rings' land with at least each number of strict comfortable hours per lunar day,
    corrected for height and in the GCM, and the mean hours."""
    cat = lambda key: np.concatenate([r[key][r['land']] for r in rings])
    cm1, gcm = hours(cat('strict')), hours(cat('gcm_strict'))
    return dict(columns=int(cm1.size), mean_hours=dict(corrected=float(cm1.mean()), gcm=float(gcm.mean()),
                                                       corrected_loose=hours(cat('loose').mean()), gcm_loose=hours(cat('gcm_loose').mean())),
                at_least={str(t): dict(corrected=float(np.mean(cm1 >= t)), gcm=float(np.mean(gcm >= t))) for t in thresholds})


def analyse(names, gcm_run: str = rc.GCM_RUN, gcm_years=rc.GCM_YEARS) -> dict:
    twin = twin_differences()
    gcm = gcm_land(gcm_run, gcm_years)
    rings = [ring_columns(name, twin['classes'], gcm) for name in names]
    return dict(schema=SCHEMA, rings=list(names), gcm=dict(run=gcm_run, years=list(gcm_years)),
                evidence=EVIDENCE_WITH.format(run=gcm_run, first=gcm_years[0], last=gcm_years[1]), reading_rule=READING_RULE,
                span_days={r['name']: r['span_days'] for r in rings}, snapshots={r['name']: r['snapshots'] for r in rings},
                twin=twin, gcm_height_rates=height_rates(gcm),
                settling={r['name']: r['settling'] for r in rings},
                against_gcm={r['name']: r['against_gcm'] for r in rings},
                sea_level={r['name']: sea_level_bands(r) for r in rings},
                corrected_for_height=by_latitude_and_height(rings), land_shares=land_shares(rings))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('rings', nargs='+')
    parser.add_argument('--gcm', metavar='RUN:FIRST-LAST', help='another GCM run and its years (default: the design case)')
    parser.add_argument('--output', type=Path, help='where to write the comparison (default: the results folder)')
    args = parser.parse_args(argv)
    gcm_run, gcm_years = rc.GCM_RUN, rc.GCM_YEARS
    if args.gcm:
        gcm_run, span = args.gcm.split(':')
        gcm_years = tuple(int(v) for v in span.split('-'))
    result = analyse(args.rings, gcm_run, gcm_years)
    result['producer'] = dict(domain='climate', files={'crm/ring_comfort.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    path = args.output or RESULTS / ('comfort_flat_rings.json' if gcm_run == rc.GCM_RUN and tuple(gcm_years) == tuple(rc.GCM_YEARS)
                                     else f'comfort_flat_rings_{gcm_run}.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rc.rounded(result), indent=1) + '\n')
    shares = result['land_shares']
    print(f'{path}: {shares["columns"]} land columns; strict hours per lunar day corrected for height '
          f'{shares["mean_hours"]["corrected"]:.0f}, GCM {shares["mean_hours"]["gcm"]:.0f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
