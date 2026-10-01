"""What the highland box shows: the weather of high ground through the lunar day, by the ground's height and its place
in the terrain, beside the flat ring corrected for height and the GCM.

    climate/gcm/.venv/bin/python -m climate.crm.highland_analysis box_highland --from-day 29.5

writes ../results/crm/highland_<case>.json.

Every column of the box stands at its own height. The analysis reads the span's snapshots (every 3 model hours) of the
columns beyond the box's blended edges and groups them by height, in bands of a kilometre, and by their place in the
terrain: valleys and basins (at least 150 m below the mean ground within 30 km), ridges and summits (at least 150 m
above it) and the slopes between. For each group it reports the air and dewpoint at 2 m by day and by night, the wind
at 10 m, the rain, cloud and fog (cloud within 200 m of the ground), and the hours comfortable to be out in; how the
air and dewpoint change with height by day and by night; the wind blowing up or down the slopes through the lunar
day; how much warmer or colder valleys run than ridges against the box's change with height, by day and by night; and
the rain on slopes where the wind at 10 m blows uphill against those where it blows downhill. Beside them stand the
flat ring (the ring's land within the site's radius, at sea level) corrected to each band's height with the terrain
twin's differences (ring_comfort), and the GCM's land cells near the site.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from climate.crm import ring_analysis as ra
from climate.crm import ring_crossings as rc
from climate.crm import ring_comfort as rcf

RESULTS = ra.RESULTS
SCHEMA = 'terluna.climate.crm-highland/1'
BANDS_M = tuple((lo, lo + 1000) for lo in range(0, 7000, 1000))
POSITION_RADIUS_M = 30.0e3              # the neighbourhood whose mean ground sets a column's place in the terrain
POSITION_M = 150.0                       # below it a valley, above it a ridge
STEEP = 0.05                             # slopes at least this steep count for slope winds and windward rain
FOG_TOP_M = 200.0
HOUR_BINS = 12                           # 30 degrees of local time, about 2.5 Earth days
EVIDENCE = ('CM1 r22.0 in three dimensions at lunar gravity over the atlas\'s ground (smoothed to slopes of at most '
            'about 0.2 and blended to a common height at the edges, where the box wraps round), the upper air held to '
            'the corrected GCM\'s along the ring, the GCM\'s mean vertical wind and the ring\'s own day-night '
            'circulation prescribed, the ground\'s deep temperature carried from sea level at the GCM\'s lapse rate. It '
            'resolves slope winds, cold air in valleys and storms over the terrain; it has no land and water around it '
            'and takes the planet-wide circulation from the ring. Like every CM1 case here it runs warmer and more humid '
            'near the ground than the GCM, for reasons not yet found, so its absolute warmth and humidity are CM1\'s end '
            'of the range; the contrasts within the box rest on the terrain and the Sun.')
READING_RULE = ('Columns within 50 km of the box\'s edge, where its ground is blended, are left out. Heights are the '
                'ground\'s above sea level. Day and night follow the Sun at the site; hours per lunar day are shares of '
                'the snapshots times the lunar day. Comfortable means air at 2 m of 18-26 C with a dewpoint of at most '
                '15 C (strict), or 17-27 C and 16 C (loose). The change with height is a least-squares line over the '
                'columns; a valley\'s or ridge\'s warmth is its mean departure from that line. The upslope wind is the '
                'wind at 10 m along the ground\'s steepest ascent on slopes of at least 0.05. Windward slopes are those '
                'where the wind at 10 m blows uphill at that moment. The corrected flat ring adds the terrain twin\'s '
                'difference at each band\'s mean height to the ring\'s land near the site over the same span; the GCM\'s '
                'values are its 3-day means on its land cells within 150 km, at their own heights.')


def periodic_mean(field, radius_cells: int):
    """The mean of a periodic field (ny, nx) over a disk of radius_cells about each point."""
    total, count = np.zeros_like(field, dtype=float), 0
    for dy in range(-radius_cells, radius_cells + 1):
        for dx in range(-radius_cells, radius_cells + 1):
            if dx * dx + dy * dy <= radius_cells * radius_cells:
                total += np.roll(np.roll(field, dy, 0), dx, 1)
                count += 1
    return total / count


def terrain_place(heights, dx: float, radius_m=POSITION_RADIUS_M, margin_m=POSITION_M):
    """-1 for valleys and basins, +1 for ridges and summits, 0 for the slopes between."""
    rel = heights - periodic_mean(heights, int(round(radius_m / dx)))
    return np.where(rel <= -margin_m, -1, np.where(rel >= margin_m, 1, 0))


def slopes(heights, dx: float):
    """The ground's gradient (m/m) along x and y by centred differences on the periodic box."""
    gx = (np.roll(heights, -1, 1) - np.roll(heights, 1, 1)) / (2.0 * dx)
    gy = (np.roll(heights, -1, 0) - np.roll(heights, 1, 0)) / (2.0 * dx)
    return gx, gy


def interior(shape, dx: float, margin_m: float):
    """Columns farther than margin_m from every edge."""
    ny, nx = shape
    m = int(np.ceil(margin_m / dx))
    keep = np.zeros(shape, bool)
    keep[m:ny - m, m:nx - m] = True
    return keep


def line(x, y):
    """Slope and intercept of the least-squares line."""
    slope, icpt = np.polyfit(np.asarray(x, float), np.asarray(y, float), 1)
    return float(slope), float(icpt)


def reader(case: Path):
    """A reader of one snapshot: the 2-D fields the analysis needs and whether any level, or one within 200 m of the
    ground, holds cloud."""
    nx, layout = ra.grads_layout(case / 'cm1out_s.ctl')
    offsets, offset = {}, 0
    for name, levels in layout:
        offsets[name] = (offset, levels)
        offset += levels * nx
    zh = ra.case_geometry(case)['zh']
    low = int(np.searchsorted(zh, FOG_TOP_M))

    def read(n: int) -> dict:
        raw = np.memmap(case / f'cm1out_t{n:06d}_s.dat', dtype='<f4', mode='r')
        flat = {k: np.array(raw[offsets[k][0]:offsets[k][0] + nx], float) for k in
                ('t2', 'q2', 'psfc', 'prate', 'u10', 'v10', 's10', 'coszen')}
        cond = sum(np.array(raw[offsets[k][0]:offsets[k][0] + offsets[k][1] * nx], float).reshape(offsets[k][1], nx)
                   for k in ('qc', 'qi'))
        cloudy = cond >= ra.CLOUD_KG_KG
        flat['cloud'] = cloudy.any(axis=0)
        flat['fog'] = cloudy[:max(low, 1)].any(axis=0)
        return flat
    return read


def stat(fn, values, *args):
    """fn of the values, or None where there are none (a span without daylight, say)."""
    values = np.asarray(values)
    return float(fn(values, *args)) if values.size else None


def comfortable_hours(air, dew, band='strict') -> float:
    return rcf.hours(rc.comfortable(np.asarray(air), np.asarray(dew), band).mean())


def analyse(name: str, from_day: float) -> dict:
    from climate.crm.cm1_run import CASES, solar_day_s
    case = ra.RUNS / name
    cfg = CASES[name]
    record = json.loads((case / 'case.json').read_text())
    grid, site = record['grid'], cfg['site']
    ny, nx, dx = grid['ny'], grid['nx'], grid['dx_m']
    heights = np.fromfile(case / 'perts.dat', dtype='<f4').astype(float).reshape(ny, nx)
    keep = interior(heights.shape, dx, record['terrain']['taper_km'] * 1000.0)
    place = terrain_place(heights, dx)
    gx, gy = slopes(heights, dx)
    steep = keep & (np.hypot(gx, gy) >= STEEP)
    day_s, output_s = solar_day_s(), record['configuration']['output_s']
    numbers = rcf.snapshots(case, output_s, from_day * 86400.0, np.inf)
    read = reader(case)
    shape2 = (len(numbers), ny, nx)
    f = {k: np.zeros(shape2, np.float32) for k in ('air', 'dew', 'rain', 'wind', 'up', 'cloud', 'fog')}
    sun, hour = np.zeros(len(numbers), bool), np.zeros(len(numbers))
    start = record['configuration']['start_hour_angle_deg'] + site['lon_deg']
    for i, n in enumerate(numbers):
        d = read(n)
        f['air'][i] = (d['t2'] - 273.15).reshape(ny, nx)
        f['dew'][i] = rc.dewpoint_c(ra.vapour_pa(d['q2'], d['psfc'])).reshape(ny, nx)
        f['rain'][i] = (d['prate'] * 86400.0).reshape(ny, nx)                       # mm/day
        f['wind'][i] = d['s10'].reshape(ny, nx)
        u, v = d['u10'].reshape(ny, nx), d['v10'].reshape(ny, nx)
        f['up'][i] = np.where(steep, (u * gx + v * gy) / np.maximum(np.hypot(gx, gy), 1e-9), 0.0)
        f['cloud'][i], f['fog'][i] = d['cloud'].reshape(ny, nx), d['fog'].reshape(ny, nx)
        sun[i] = d['coszen'].mean() > 0.0
        hour[i] = (start + 360.0 * (n - 1) * output_s / day_s + 180.0) % 360.0 - 180.0
    day, night = sun[:, None, None] & keep[None], ~sun[:, None, None] & keep[None]
    h = np.broadcast_to(heights, shape2)

    def group(mask2d, label):
        m = keep & mask2d
        if m.sum() < 5:
            return None
        md, mn = sun[:, None, None] & m[None], ~sun[:, None, None] & m[None]
        all_ = np.broadcast_to(m, shape2)
        return dict(group=label, columns=int(m.sum()), mean_height_m=float(heights[m].mean()),
                    air_day_c=dict(median=stat(np.median, f['air'][md]), p95=stat(np.percentile, f['air'][md], 95)),
                    air_night_c=dict(p5=stat(np.percentile, f['air'][mn], 5), median=stat(np.median, f['air'][mn])),
                    dewpoint_c=dict(day=stat(np.median, f['dew'][md]), night=stat(np.median, f['dew'][mn])),
                    wind_10m_m_s=dict(median=float(np.median(f['wind'][all_])), p95=float(np.percentile(f['wind'][all_], 95))),
                    rain_mm_per_lunar_day=float(f['rain'][all_].mean() * day_s / 86400.0),
                    cloud=float(f['cloud'][all_].mean()), fog_night=stat(np.mean, f['fog'][mn]),
                    comfortable_hours=dict(strict=comfortable_hours(f['air'][all_], f['dew'][all_]),
                                           loose=comfortable_hours(f['air'][all_], f['dew'][all_], 'loose')),
                    uncomfortable_hours=dict(too_hot=rcf.hours((f['air'][all_] > 26.0).mean()),
                                             too_cold=rcf.hours((f['air'][all_] < 18.0).mean()),
                                             too_humid=rcf.hours(((f['air'][all_] >= 18.0) & (f['air'][all_] <= 26.0)
                                                                  & (f['dew'][all_] > 15.0)).mean())))

    by_height = [g for lo, hi in BANDS_M if (g := group((heights >= lo) & (heights < hi), f'{lo / 1000:.0f}-{hi / 1000:.0f} km'))]
    by_place = [g for p, label in ((-1, 'valleys and basins'), (0, 'slopes'), (1, 'ridges and summits'))
                if (g := group(place == p, label))]
    change, warmth = {}, {}
    for part, sel in (('day', day), ('night', night)):
        if not sel.any():
            continue
        for key in ('air', 'dew'):
            change[f'{key}_{part}_per_km'] = line(h[sel] / 1000.0, f[key][sel])[0]
        slope, icpt = line(h[sel], f['air'][sel])
        resid = f['air'] - (icpt + slope * h)
        for p, label in ((-1, 'valleys'), (1, 'ridges')):
            m = sel & (place == p)[None]
            warmth[f'{label}_{part}_c'] = float(resid[m].mean())
    bins = np.linspace(-180.0, 180.0, HOUR_BINS + 1)
    k = np.clip(np.digitize(hour, bins) - 1, 0, HOUR_BINS - 1)
    upslope = [dict(hour_angle_deg=float(0.5 * (bins[b] + bins[b + 1])),
                    upslope_m_s=float(f['up'][k == b][:, steep].mean()) if (k == b).any() else None,
                    share_uphill=float((f['up'][k == b][:, steep] > 0).mean()) if (k == b).any() else None)
               for b in range(HOUR_BINS)]
    windward = (f['up'] > 0) & steep[None]
    leeward = (f['up'] < 0) & steep[None]
    rain_slopes = dict(windward_mm_day=float(f['rain'][windward].mean()), leeward_mm_day=float(f['rain'][leeward].mean()),
                       windward_share_of_slope_time=float(windward.sum() / max(1, (windward | leeward).sum())))
    return dict(schema=SCHEMA, case=name, evidence=EVIDENCE, reading_rule=READING_RULE,
                span_days=[(numbers[0] - 1) * output_s / 86400.0, (numbers[-1] - 1) * output_s / 86400.0],
                snapshots=len(numbers), columns=int(keep.sum()), site=dict(site, heading_deg=record['terrain']['heading_deg']),
                heights_m=dict(min=float(heights[keep].min()), max=float(heights[keep].max()), mean=float(heights[keep].mean())),
                places=dict(valleys=int((keep & (place == -1)).sum()), slopes=int((keep & (place == 0)).sum()),
                            ridges=int((keep & (place == 1)).sum()), steep=int(steep.sum())),
                box=group(np.ones_like(keep), 'all'), by_height=by_height, by_place=by_place, change_with_height=change,
                valleys_and_ridges=warmth, upslope_wind_by_local_time=upslope, rain_on_slopes=rain_slopes,
                flat_ring=flat_ring_corrected(cfg, record, from_day, by_height),
                gcm=gcm_near_site(site))


def flat_ring_corrected(cfg, record, from_day: float, bands: list) -> dict:
    """The ring's land within the site's radius over the same span, at sea level and corrected to each band's mean
    height with the terrain twin's differences, by day and by night: air, dewpoint and comfortable hours."""
    from climate.crm.cm1_run import solar_day_s
    ring = ra.RUNS / cfg['site']['ring']
    geo = ra.case_geometry(ring)
    cols = np.asarray(record['site']['ring_columns'])
    cols = cols[geo['land'][cols]]
    read = rcf.surface_reader(ring)
    output_s = geo['record']['configuration']['output_s']
    snaps = [read(n) for n in rcf.snapshots(ring, output_s, from_day * 86400.0, np.inf)]
    air = np.array([s['air_c'][cols] for s in snaps])
    dew = np.array([s['dew_c'][cols] for s in snaps])
    day = np.array([s['day'][cols] for s in snaps])
    table = rcf.twin_differences()['classes']
    out = dict(ring=cfg['site']['ring'], columns=int(cols.size), sea_level=dict(
        air_c=float(air.mean()), dewpoint_c=float(np.median(dew)), comfortable_hours=comfortable_hours(air, dew)), bands=[])
    for b in bands:
        d_air, d_dew = rcf.correction(table, np.full(air.shape, b['mean_height_m']), day)
        a2, w2 = air + d_air, dew + d_dew
        out['bands'].append(dict(group=b['group'], mean_height_m=b['mean_height_m'], air_c=float(a2.mean()),
                                 dewpoint_c=float(np.median(w2)), comfortable_hours=comfortable_hours(a2, w2)))
    return out


def gcm_near_site(site: dict) -> dict:
    """The GCM's land cells within 150 km of the site: their heights, air, lowest-layer dewpoint and comfortable hours."""
    near = rc.gcm_near([dict(lat_deg=site['lat_deg'], lon_deg=site['lon_deg'])])[0]
    land = near.get('land')
    if not land:
        return dict(near)
    g = rcf.gcm_land()
    cells = [(int(np.argmin(np.abs(g['lat'] - la))), int(np.argmin(np.abs((g['lon'] - lo + 180.0) % 360.0 - 180.0))))
             for la, lo in land['cells']]
    return dict(run=near['run'], years=near['years'], cells=len(cells),
                heights_m=[float(g['height_m'][i, j]) for i, j in cells],
                air_c=land['air_c']['mean'], dewpoint_c=land['dewpoint_c']['mean'],
                comfortable_hours=land['comfortable_hours_per_lunar_day'], rain_mm_per_lunar_day=land['rain_mm_per_lunar_day'])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--output', type=Path, help='where to write the summary (default: the results folder)')
    args = parser.parse_args(argv)
    result = analyse(args.case, args.from_day)
    result['producer'] = dict(domain='climate', files={f'crm/{Path(p).name}': hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
                                                       for p in (__file__, rcf.__file__)})
    path = args.output or RESULTS / f'highland_{args.case}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rc.rounded(result), indent=1) + '\n')
    b = result['box']
    print(f'{path}: {result["snapshots"]} snapshots, {result["columns"]} columns; strict comfortable hours '
          f'{b["comfortable_hours"]["strict"]:.0f} per lunar day')
    return 0


if __name__ == '__main__':
    sys.exit(main())
