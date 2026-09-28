"""What a three-dimensional box shows against the rings it extends and the GCM, at its site by local time: the air
near the ground and as the GCM measures it, humidity, cloud, fog, rain and wind.

    climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e --from-day 29.5

writes ../results/crm/box_<case>.json.

A box is a ring's patch given a second horizontal dimension, forced as the ring was there. The analysis takes every
column of the box and each ring's land columns within 100 km of the site, in every 3-hour snapshot of the analysed
span, and compares them snapshot by snapshot at the same local time (ring_crossings.compare): each one's composites
by local time and means over the day, the night and the span, and the differences, box minus ring, with their
moving-block bootstrap ranges. Two measures join the crossing comparison's: the air as the GCM measures it, and the
share of the ground under rain. Cloud by height, by day and by night, shows how the cloud is built.
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

RUNS, RESULTS = ra.RUNS, ra.RESULTS
SCHEMA = 'terluna.climate.crm-box/1'
KEYS = (*rc.VARIABLES, 'air_gcm_layer_c', 'raining')
PROFILE_TOP_M = 90000.0
EVIDENCE = ('CM1 r22.0 in three dimensions at lunar gravity: a box 385 km square of 64 x 64 columns 6.0 km wide on 111 '
            'levels to 150 km, wrapping round on all sides, at 0 N 0 E, all land (the placeholder surface at wetness '
            '0.90), with the Morrison microphysics given lunar fall speeds, RRTMG radiation for the design air and the '
            'Sun of the site, the same across the box; the upper air above 8-16 km held to ring A\'s reference; the GCM '
            'design case\'s mean vertical wind there imposed as large-scale vertical advection; and the heating and '
            'moistening of the rings\' own day-night circulation there prescribed by local time and height. It '
            'resolves storms in three dimensions. It cannot make the planet-wide circulation, which it takes from the '
            'rings, and has no land and water around it.')
READING_RULE = ('The box takes all its columns; each ring, its land columns within 100 km of the site along the ring. '
                'Values are means over those columns of snapshots every 3 model hours, composited by the local hour angle '
                'at the site (0 = noon, negative = morning) in 36-degree bins; day is when the Sun is up. Differences are '
                'the box minus the ring, snapshot by snapshot at the same local time; p5 and p95 bound 90% of the means '
                'of series rebuilt from one-day blocks drawn at random (a moving-block bootstrap, 2,000 draws). '
                'air_gcm_layer_c is the mass-weighted mean temperature of the GCM\'s lowest layer (the bottom 3.3% of '
                'the air\'s mass, about 1.8 km deep) carried to the ground along a dry adiabat, as the GCM\'s own air near '
                'the ground is; raining is the share of columns with rain of at least 1 mm/h. Comfortable means air at '
                '2 m of 18-26 C with a dewpoint of at most 15 C (strict), or 17-27 C and 16 C (loose). The GCM values '
                'are its cells whose centres lie within 150 km, as in ring_crossings. Cloud profiles are the share of '
                'columns with condensate of at least 0.01 g/kg at each height, by day and by night.')


def gcm_layer_air_c(d: dict, sigma: float) -> np.ndarray:
    """The air as the GCM measures it, column by column: the mass-weighted mean temperature of the GCM's lowest
    layer (from the ground up to 2 sigma - 1 of the surface pressure, the full level sigma midway), carried to the
    ground along a dry adiabat."""
    kappa = 287.04 / 1004.64
    t = d['th'] * (d['prs'] / 1.0e5) ** kappa
    p, ps = d['prs'], d['psfc']
    lower = np.vstack([ps[None, :], 0.5 * (p[1:] + p[:-1])])            # each level's lower and upper face
    upper = np.vstack([lower[1:], p[-1:]])
    mass = np.clip(np.minimum(lower, ps[None, :]) - np.maximum(upper, (2.0 * sigma - 1.0) * ps[None, :]), 0.0, None)
    return (t * mass).sum(axis=0) / mass.sum(axis=0) * sigma ** -kappa - 273.15


def fields(d: dict, zh, sigma: float) -> dict:
    """The crossing comparison's variables column by column, with the air as the GCM measures it and rain cover."""
    out = rc.snapshot_fields(d, zh, sigma)
    out.update(air_gcm_layer_c=gcm_layer_air_c(d, sigma), raining=d['prate'] * 3600.0 >= ra.RAIN_MM_H)
    return out


def cloud_profiles(case: Path, columns, times, start_deg: float, lon_deg: float, day_s: float) -> dict:
    """The share of the columns with cloud at each height up to PROFILE_TOP_M, by day and by night at the site,
    over the snapshots at the given model times (s)."""
    geo = ra.case_geometry(case)
    tap = geo['record']['configuration']['output_s']
    keep = geo['zh'] <= PROFILE_TOP_M
    total, count = np.zeros((2, int(keep.sum()))), np.zeros(2)
    for n in sorted(int(round(t / tap)) + 1 for t in times):
        d = ra.read_snapshot(case, n)
        cloudy = ((d['qc'] + d['qi'])[:, columns] >= ra.CLOUD_KG_KG)[keep]
        night = int(np.cos(np.radians(rc.local_hour_angle(start_deg, lon_deg, (n - 1) * tap, day_s))) <= 0.0)
        total[night] += cloudy.mean(axis=1)
        count[night] += 1
    share = total / np.maximum(count, 1)[:, None]
    return dict(height_km=[round(float(z) / 1000.0, 2) for z in geo['zh'][keep]],
                day=[round(float(v), 4) for v in share[0]], night=[round(float(v), 4) for v in share[1]])


def analyse(name: str, from_day: float) -> dict:
    from climate.crm.cm1_run import CASES, solar_day_s
    cfg = CASES[name]
    site = cfg['site']
    rings = list(cfg['day_night']['rings']) if cfg.get('day_night') else [site['ring']]
    sigma = rc.gcm_lowest_sigma()
    day_s = solar_day_s()
    case = RUNS / name
    record = json.loads((case / 'case.json').read_text())
    start = record['configuration']['start_hour_angle_deg']
    everything = np.arange(record['grid']['nx'] * record['grid'].get('ny', 1))
    series = {name: rc.patch_series(case, [{'land': everything}], from_day, None, sigma, fields, KEYS)}
    columns = {name: everything}
    for ring in rings:
        geo = ra.case_geometry(RUNS / ring)
        cols = rc.patch_columns(geo['record']['grid']['nx'], rc.along_ring(rc.ring_path(geo['record']), site['lat_deg'],
                                                                           site['lon_deg']), site['radius_m'], geo['dx'])
        columns[ring] = cols[geo['land'][cols]]
        series[ring] = rc.patch_series(RUNS / ring, [{'land': columns[ring]}], from_day, None, sigma, fields, KEYS)
    times = sorted(set.intersection(*(set(s['time_s']) for s in series.values())))
    if not times:
        raise RuntimeError(f'no output common to {name} and {", ".join(rings)} after day {from_day}')
    pick = lambda s: {key: np.asarray(s['patches'][0]['land'][key])[[s['time_s'].index(t) for t in times]] for key in KEYS}
    aligned = {k: pick(s) for k, s in series.items()}
    hour = [rc.local_hour_angle(start, site['lon_deg'], t, day_s) for t in times]
    against = {ring: rc.compare(hour, {name: aligned[name], ring: aligned[ring]}, day_s / 3600.0, KEYS) for ring in rings}
    return dict(schema=SCHEMA, case=name, site=dict(site, land_columns={k: int(v.size) for k, v in columns.items()}),
                rings=rings, evidence=EVIDENCE, reading_rule=READING_RULE, gcm_lowest_sigma=sigma,
                span_days=[times[0] / 86400.0, times[-1] / 86400.0], snapshots=len(times), against=against,
                gcm=rc.gcm_near([dict(lat_deg=site['lat_deg'], lon_deg=site['lon_deg'])])[0],
                cloud_profiles={k: cloud_profiles(RUNS / k, columns[k], times, start, site['lon_deg'], day_s)
                                for k in (name, *rings)})


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--output', type=Path, help='where to write the comparison (default: the results folder)')
    args = parser.parse_args(argv)
    result = analyse(args.case, args.from_day)
    result['producer'] = dict(domain='climate', files={f'crm/{Path(f).name}': hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]
                                                       for f in (__file__, rc.__file__)})
    path = args.output or RESULTS / f'box_{args.case}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rc.rounded(result), indent=1) + '\n')
    print(f'{path}: {result["snapshots"]} snapshots, days {result["span_days"][0]:.1f}-{result["span_days"][1]:.1f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
