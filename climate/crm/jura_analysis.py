"""What the Jura Island box shows: the island's weather through the lunar day by the ground's height, its distance
from the sea and the water it faces, the bay's air against the open sea's, the sea breeze and the arc's cloud, beside
the flat ring corrected for height and the GCM.

    climate/gcm/.venv/bin/python -m climate.crm.jura_analysis box_jura --from-day 29.5

writes ../results/crm/jura_<case>.json.

Water columns are those the setup made sea (land use 16). Land faces the bay when its nearest water lies within Sinus
Iridum's circle (the IAU gazetteer's centre and diameter) and the open sea otherwise. The analysis reads the span's
snapshots (every 3 model hours) of the columns beyond the box's blended edges. Each column keeps its own day and night
and its own local time, since the box spans about 28 degrees of longitude. For each group it reports, as
highland_analysis does, the air and dewpoint at 2 m by day and by night, the wind at 10 m, the rain, cloud and fog
(cloud within 200 m of the ground) and the hours comfortable to be out in; then the onshore wind on the coasts and the
cloud over the land by height, through the lunar day.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from climate.crm import highland_analysis as ha
from climate.crm import ring_analysis as ra
from climate.crm import ring_crossings as rc
from climate.crm import ring_comfort as rcf

RESULTS = ra.RESULTS
SCHEMA = 'terluna.climate.crm-jura/1'
BANDS_M = ((0, 100), (100, 250), (250, 500), (500, 1000), (1000, 1500), (1500, 2500))
DISTANCES_KM = ((0, 7), (7, 13), (13, 25), (25, 1000))         # about one, two and four columns from the water
LOW_M = 500.0                                                   # the low ground, where the gardens of the layout lie
HOUR_BINS = 12                                                  # 30 degrees of local time, about 2.5 Earth days
WATER_LU = 16
EVIDENCE = ('CM1 r22.0 in three dimensions at lunar gravity over Jura Island at the atlas\'s heights (smoothed to the '
            'box\'s 6 km and blended to the sea at the edges, where the box wraps round) with the atlas\'s sea and '
            'islets round it at one temperature; the upper air held to the corrected GCM\'s along the 45 N row, the '
            'GCM\'s mean vertical wind and the 45 N twilight ring\'s day-night circulation prescribed, the curved '
            'air\'s twilight on. The land is the placeholder land, bare and at the GCM\'s dryness there, so the '
            'island\'s shore belts and gardens are missing. CM1 runs warmer near the ground than the GCM, for reasons '
            'still under test, so its absolute warmth is CM1\'s end of the range; the contrasts within the island '
            'rest on the terrain, the sea and the Sun.')
READING_RULE = ('Columns within 50 km of the box\'s edge are left out. Heights are the ground\'s above sea level and '
                'distances the nearest water column\'s, between column centres. Day and night are each column\'s own '
                '(the Sun above its horizon or not); hours per lunar day are shares of the snapshots times the lunar '
                'day. Comfortable means air at 2 m of 18-26 C with a dewpoint of at most 15 C (strict), or 17-27 C and '
                '16 C (loose). The onshore wind is the wind at 10 m along the line from a coastal column\'s nearest '
                'water to it, positive from the water. The flat ring and the GCM are as highland_analysis gives them.')


def sinus_iridum():
    """The bay's circle: centre (degrees north, east) and radius (km), from the IAU gazetteer."""
    from geography import nomenclature as nm
    f = nm.by_name(nm.features())['Sinus Iridum']
    return float(f['lat']), float(f['lon']) % 360.0, 0.5 * float(f['diameter_km'])


def great_circle_km(lat1, lon1, lat2, lon2, radius_km):
    a, b = np.radians(lat1), np.radians(lat2)
    d = np.sin(0.5 * (b - a)) ** 2 + np.cos(a) * np.cos(b) * np.sin(0.5 * np.radians(lon2 - lon1)) ** 2
    return 2.0 * radius_km * np.arcsin(np.sqrt(d))


def coast_geometry(water, dx: float):
    """For every column: the distance (km) to the nearest water column and the unit vector (x, y) from that water to
    the column, (0, 0) on water."""
    from scipy.ndimage import distance_transform_edt
    distance, (jw, iw) = distance_transform_edt(~water, return_indices=True)
    j, i = np.indices(water.shape)
    ex, ey = (i - iw).astype(float), (j - jw).astype(float)
    norm = np.maximum(np.hypot(ex, ey), 1e-9)
    return distance * dx / 1000.0, ex / norm, ey / norm, jw, iw


def gcm_near_site(site: dict, run: str, years) -> dict:
    """The GCM's land cells within 150 km of the site in the run that forced the box: their heights, air, lowest-layer
    dewpoint, comfortable hours and rain (highland_analysis.gcm_near_site for another run)."""
    near = rc.gcm_near([dict(lat_deg=site['lat_deg'], lon_deg=site['lon_deg'])], run=run, years=years)[0]
    land = near.get('land')
    if not land:
        return dict(near)
    g = rcf.gcm_land(run, years)
    cells = [(int(np.argmin(np.abs(g['lat'] - la))), int(np.argmin(np.abs((g['lon'] - lo + 180.0) % 360.0 - 180.0))))
             for la, lo in land['cells']]
    return dict(run=near['run'], years=near['years'], cells=len(cells),
                heights_m=[float(g['height_m'][i, j]) for i, j in cells],
                air_c=land['air_c']['mean'], dewpoint_c=land['dewpoint_c']['mean'],
                comfortable_hours=land['comfortable_hours_per_lunar_day'],
                rain_mm_per_lunar_day=land['rain_mm_per_lunar_day'])


def analyse(name: str, from_day: float) -> dict:
    from climate.crm.cm1_run import CASES, box_columns, planet, solar_day_s
    case = ra.RUNS / name
    cfg = CASES[name]
    record = json.loads((case / 'case.json').read_text())
    grid, site = record['grid'], cfg['site']
    ny, nx, dx = grid['ny'], grid['nx'], grid['dx_m']
    heights = np.fromfile(case / 'perts.dat', dtype='<f4').astype(float).reshape(ny, nx)
    lu = np.loadtxt(case / 'terluna_surface2d.txt', skiprows=1)[:, 2].reshape(ny, nx)
    water = lu == WATER_LU
    keep = ha.interior(heights.shape, dx, record['terrain']['taper_km'] * 1000.0)
    land = keep & ~water
    lat, lon = box_columns(site['lat_deg'], site['lon_deg'], record['terrain']['heading_deg'], nx, ny, dx)
    blat, blon, bradius = sinus_iridum()
    in_bay = great_circle_km(lat, lon, blat, blon, planet()['radius_m'] / 1000.0) <= bradius
    distance, ex, ey, jw, iw = coast_geometry(water, dx)
    faces_bay = in_bay[jw, iw]
    day_s, output_s = solar_day_s(), record['configuration']['output_s']
    numbers = rcf.snapshots(case, output_s, from_day * 86400.0, np.inf)
    read = ha.reader(case)
    shape2 = (len(numbers), ny, nx)
    f = {k: np.zeros(shape2, np.float32) for k in ('air', 'dew', 'rain', 'wind', 'onshore', 'cloud', 'fog', 'sun')}
    start = record['configuration']['start_hour_angle_deg'] + site['lon_deg']
    offset = (lon - site['lon_deg'] + 180.0) % 360.0 - 180.0                       # each column's longitude east
    hour = np.zeros(shape2, np.float32)
    for i, n in enumerate(numbers):
        d = read(n)
        f['air'][i] = (d['t2'] - 273.15).reshape(ny, nx)
        f['dew'][i] = rc.dewpoint_c(ra.vapour_pa(d['q2'], d['psfc'])).reshape(ny, nx)
        f['rain'][i] = (d['prate'] * 86400.0).reshape(ny, nx)                       # mm/day
        f['wind'][i] = d['s10'].reshape(ny, nx)
        f['onshore'][i] = d['u10'].reshape(ny, nx) * ex + d['v10'].reshape(ny, nx) * ey
        f['cloud'][i], f['fog'][i] = d['cloud'].reshape(ny, nx), d['fog'].reshape(ny, nx)
        f['sun'][i] = d['coszen'].reshape(ny, nx) > 0.0
        hour[i] = (start + 360.0 * (n - 1) * output_s / day_s + offset + 180.0) % 360.0 - 180.0
    sun = f['sun'].astype(bool)

    def group(mask2d, label):
        m = keep & mask2d
        if m.sum() < 5:
            return None
        sel = np.broadcast_to(m, shape2)
        md, mn = sel & sun, sel & ~sun
        a, w = f['air'][sel], f['dew'][sel]
        return dict(group=label, columns=int(m.sum()), mean_height_m=float(heights[m].mean()),
                    air_day_c=dict(median=ha.stat(np.median, f['air'][md]),
                                   p95=ha.stat(np.percentile, f['air'][md], 95)),
                    air_night_c=dict(p5=ha.stat(np.percentile, f['air'][mn], 5),
                                     median=ha.stat(np.median, f['air'][mn])),
                    dewpoint_c=dict(day=ha.stat(np.median, f['dew'][md]), night=ha.stat(np.median, f['dew'][mn])),
                    wind_10m_m_s=dict(median=float(np.median(f['wind'][sel])),
                                      p95=float(np.percentile(f['wind'][sel], 95))),
                    rain_mm_per_lunar_day=float(f['rain'][sel].mean() * day_s / 86400.0),
                    cloud=float(f['cloud'][sel].mean()), fog_night=ha.stat(np.mean, f['fog'][mn]),
                    comfortable_hours=dict(strict=ha.comfortable_hours(a, w),
                                           loose=ha.comfortable_hours(a, w, 'loose')),
                    uncomfortable_hours=dict(too_hot=rcf.hours((a > 26.0).mean()),
                                             too_cold=rcf.hours((a < 18.0).mean()),
                                             too_humid=rcf.hours(((a >= 18.0) & (a <= 26.0) & (w > 15.0)).mean())))

    by_height = [g for lo, hi in BANDS_M if (g := group(~water & (heights >= lo) & (heights < hi), f'{lo}-{hi} m'))]
    by_distance = [g for lo, hi in DISTANCES_KM
                   if (g := group(~water & (distance > lo) & (distance <= hi),
                                  f'{lo}-{hi} km from the sea' if hi < 1000 else f'over {lo} km from the sea'))]
    by_side = [g for mask, label in ((faces_bay, 'facing the bay'), (~faces_bay, 'facing the open sea'),
                                     (faces_bay & (heights < LOW_M), 'facing the bay, below 500 m'),
                                     (~faces_bay & (heights < LOW_M), 'facing the open sea, below 500 m'))
               if (g := group(~water & mask, label))]
    waters = [g for mask, label in ((water & in_bay, 'the bay'), (water & ~in_bay, 'the open sea'))
              if (g := group(mask, label))]
    bins = np.linspace(-180.0, 180.0, HOUR_BINS + 1)
    k = np.clip(np.digitize(hour, bins) - 1, 0, HOUR_BINS - 1)
    coastal = land & (distance <= 13.0)

    def by_time(field, mask2d):
        sel = np.broadcast_to(mask2d, shape2)
        return [float(field[sel & (k == b)].mean()) if (sel & (k == b)).any() else None for b in range(HOUR_BINS)]

    through_the_day = dict(
        hour_angle_deg=[float(0.5 * (bins[b] + bins[b + 1])) for b in range(HOUR_BINS)],
        onshore_wind_m_s={label: by_time(f['onshore'], coastal & mask) for label, mask in
                          (('coasts facing the bay', faces_bay), ('coasts facing the open sea', ~faces_bay))},
        cloud={label: by_time(f['cloud'], mask) for label, mask in
               (('land below 250 m', land & (heights < 250.0)),
                ('land 250-1,000 m', land & (heights >= 250.0) & (heights < 1000.0)),
                ('land above 1,000 m', land & (heights >= 1000.0)), ('the bay', keep & water & in_bay),
                ('the open sea', keep & water & ~in_bay))},
        air_c={label: by_time(f['air'], mask) for label, mask in
               (('land below 250 m', land & (heights < 250.0)), ('land above 1,000 m', land & (heights >= 1000.0)),
                ('the bay', keep & water & in_bay), ('the open sea', keep & water & ~in_bay))})
    return dict(schema=SCHEMA, case=name, purpose=record['configuration'].get('purpose'),
                evidence=EVIDENCE, reading_rule=READING_RULE,
                span_days=[(numbers[0] - 1) * output_s / 86400.0, (numbers[-1] - 1) * output_s / 86400.0],
                snapshots=len(numbers), columns=dict(kept=int(keep.sum()), land=int(land.sum()),
                                                     bay=int((keep & water & in_bay).sum()),
                                                     open_sea=int((keep & water & ~in_bay).sum()),
                                                     land_facing_bay=int((land & faces_bay).sum())),
                bay=dict(lat_deg=blat, lon_deg=blon, radius_km=bradius),
                heights_m=dict(max=float(heights[land].max()), mean=float(heights[land].mean())),
                island=group(~water, 'the island'), by_height=by_height, by_distance=by_distance, by_side=by_side,
                waters=waters, through_the_day=through_the_day,
                flat_ring=ha.flat_ring_corrected(cfg, record, from_day, by_height),
                gcm=gcm_near_site(site, *cfg.get('gcm', (rc.GCM_RUN, rc.GCM_YEARS))))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    args = parser.parse_args(argv)
    result = analyse(args.case, args.from_day)
    files = (__file__, ha.__file__, rcf.__file__)
    result['producer'] = dict(domain='climate', files={f'crm/{Path(p).name}': hashlib.sha256(Path(p).read_bytes())
                                                       .hexdigest()[:16] for p in files})
    path = RESULTS / f'jura_{args.case}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rc.rounded(result), indent=1) + '\n')
    i = result['island']
    print(f'{path}: {result["snapshots"]} snapshots, {result["columns"]["land"]} land columns; strict comfortable '
          f'hours {i["comfortable_hours"]["strict"]:.0f} per lunar day, rain {i["rain_mm_per_lunar_day"]:.2f} mm')
    return 0


if __name__ == '__main__':
    sys.exit(main())
