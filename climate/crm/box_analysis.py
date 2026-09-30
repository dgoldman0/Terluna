"""What a three-dimensional box shows against the rings it extends and the GCM, at its site by local time: the air
near the ground and as the GCM measures it, humidity, cloud, fog, rain and wind.

    climate/gcm/.venv/bin/python -m climate.crm.box_analysis box_0e --from-day 29.5

writes ../results/crm/box_<case>.json.

A box is a ring's patch given a second horizontal dimension, forced as the ring was there. The analysis takes every
column of the box and each ring's land columns within 100 km of the site, in every 3-hour snapshot of the analysed
span, and compares them snapshot by snapshot at the same local time (ring_crossings.compare): each one's composites
by local time and means over the day, the night and the span, and the differences, box minus ring, with their
moving-block bootstrap ranges. Two measures join the crossing comparison's: the air as the GCM measures it, and the
share of the ground under rain. Cloud by height, by day and by night, shows how the cloud is built, and the
energy at the ground (sunlight and infrared where the case writes them, the heat carried into the air and the heat
of evaporation) shows where the warmth comes from, against the GCM's own surface budget.
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
            'rings, and has no land and water around it. The boxes built on it (size, land, fall speeds) say how they '
            'depart from it in their purpose, given with the results.')
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
                'columns with condensate of at least 0.01 g/kg at each height, by day and by night. Surface energy '
                'terms are in W/m2: net sunlight and net infrared at the ground positive downward, written only by the '
                'boxes that write the radiation at the ground; the heat carried into the air and the heat of evaporation '
                'positive upward; the GCM\'s are 3-day means on its land cells within 150 km.')


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


def surface_energy(case: Path, columns, times, start_deg: float, lon_deg: float, day_s: float) -> dict:
    """Where the energy at the ground goes, over the columns, by local time in the comparison's bins (W/m2): the
    heat carried into the air and the heat of evaporation, and where the case writes it, the net sunlight and net
    infrared at the ground (positive down)."""
    geo = ra.case_geometry(case)
    tap = geo['record']['configuration']['output_s']
    edges = np.linspace(-180.0, 180.0, rc.BINS + 1)
    sums, count = {}, np.zeros(rc.BINS)
    for n in sorted(int(round(t / tap)) + 1 for t in times):
        d = ra.read_snapshot(case, n)
        rho = ra.air_density(d)[columns]
        terms = dict(sensible=rho * 1005.7 * d['thflux'][columns], latent=rho * 2.501e6 * d['qvflux'][columns])
        if 'swdnb' in d:
            terms.update(net_sunlight=(d['swdnb'] - d['swupb'])[columns], net_infrared=(d['lwdnb'] - d['lwupb'])[columns])
        b = min(int(np.digitize(rc.local_hour_angle(start_deg, lon_deg, (n - 1) * tap, day_s), edges)) - 1, rc.BINS - 1)
        for key, value in terms.items():
            sums.setdefault(key, np.zeros(rc.BINS))[b] += float(value.mean())
        count[b] += 1
    return {key: dict(by_local_time=[float(v) for v in total / np.maximum(count, 1)],
                      mean=float(total.sum() / count.sum())) for key, total in sums.items()}


def gcm_surface_energy(lat_deg: float, lon_deg: float, radius_m=rc.GCM_RADIUS_M, run=rc.GCM_RUN, years=rc.GCM_YEARS) -> dict:
    """The GCM's energy at the ground on its land cells near a site, by local time in the same bins (W/m2, 3-day
    means): net sunlight and net infrared (positive down), and the heat carried into the air and of evaporation."""
    import netCDF4
    from climate.gcm.climatology import subsolar_longitudes
    from climate.crm.cm1_run import planet
    folder = ra.HERE.parent / 'gcm' / 'runs' / run / 'model'
    get = {k: [] for k in ('rss', 'rls', 'hfss', 'hfls', 'czen')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            for k in get:
                get[k].append(np.asarray(d[k][:], dtype=float))
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
    f = {k: np.concatenate(v) for k, v in get.items()}
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    lat2, lon2 = np.radians(np.meshgrid(lat, lon, indexing='ij'))
    near = np.arccos(np.clip(np.sin(lat2) * np.sin(la) + np.cos(lat2) * np.cos(la) * np.cos(lon2 - lo), -1.0, 1.0))
    cells = (near * planet()['radius_m'] <= radius_m) & lsm
    hour = (lon[None, :] - subsolar_longitudes(f['czen'], lat, lon)[:, None] + 180.0) % 360.0 - 180.0
    k = np.clip(((hour + 180.0) / (360.0 / rc.BINS)).astype(int), 0, rc.BINS - 1)
    sel = np.broadcast_to(cells, f['rss'].shape)
    kk = np.broadcast_to(k[:, None, :], f['rss'].shape)
    out = dict(cells=int(cells.sum()))
    for key, v in (('net_sunlight', f['rss']), ('net_infrared', f['rls']), ('sensible', -f['hfss']), ('latent', -f['hfls'])):
        out[key] = dict(by_local_time=[float(v[sel & (kk == b)].mean()) for b in range(rc.BINS)], mean=float(v[sel].mean()))
    return out


def microphysics_budget(case: Path, from_day: float) -> dict:
    """What the microphysics did with the box's water from from_day to the end, from CM1's running totals in
    cm1out_stats (kg): water condensed, evaporated back from cloud, evaporated from rain, snow and graupel on the way
    down, and reaching the ground, in mm over the box, with the share of the falling water that evaporated."""
    ctl = (case / 'cm1out_stats.ctl').read_text().splitlines()
    i = next(k for k, line in enumerate(ctl) if line.lower().startswith('vars'))
    names = [line.split()[0] for line in ctl[i + 1:i + 1 + int(ctl[i].split()[1])]]
    totals = dict(zip(names, np.fromfile(case / 'cm1out_stats.dat', dtype='<f4').reshape(-1, len(names)).T))
    grid = json.loads((case / 'case.json').read_text())['grid']
    area = grid['nx'] * grid.get('ny', 1) * grid['dx_m'] ** 2
    day = totals['mtime'] / 86400.0
    mm = lambda key: float((totals[key][-1] - np.interp(from_day, day, totals[key])) / area)      # kg/m2 = mm
    out = dict(span_days=[from_day, float(day[-1])], condensed_mm=mm('tcond'), evaporated_from_cloud_mm=mm('tevac'),
               evaporated_on_the_way_down_mm=mm('tevar'), reached_the_ground_mm=mm('train'))
    falling = out['evaporated_on_the_way_down_mm'] + out['reached_the_ground_mm']
    out['share_of_falling_water_evaporated'] = out['evaporated_on_the_way_down_mm'] / falling if falling > 0 else None
    return out


def analyse(name: str, from_day: float, reference: str | None = None) -> dict:
    """The box against the rings it extends, the GCM and, if given, a reference box at the same site."""
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
    if reference:                                                      # another box at the site, all its columns
        ref = json.loads((RUNS / reference / 'case.json').read_text())['grid']
        columns[reference] = np.arange(ref['nx'] * ref.get('ny', 1))
        series[reference] = rc.patch_series(RUNS / reference, [{'land': columns[reference]}], from_day, None, sigma,
                                            fields, KEYS)
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
    others = [*rings, *([reference] if reference else [])]
    against = {k: rc.compare(hour, {name: aligned[name], k: aligned[k]}, day_s / 3600.0, KEYS) for k in others}
    return dict(schema=SCHEMA, case=name, purpose=cfg.get('purpose'), build=record['build']['label'],
                site=dict(site, land_columns={k: int(v.size) for k, v in columns.items()}),
                rings=rings, evidence=EVIDENCE, reading_rule=READING_RULE, gcm_lowest_sigma=sigma,
                span_days=[times[0] / 86400.0, times[-1] / 86400.0], snapshots=len(times), against=against,
                gcm=rc.gcm_near([dict(lat_deg=site['lat_deg'], lon_deg=site['lon_deg'])])[0],
                cloud_profiles={k: cloud_profiles(RUNS / k, columns[k], times, start, site['lon_deg'], day_s)
                                for k in (name, *others)},
                surface_energy=dict({k: surface_energy(RUNS / k, columns[k], times, start, site['lon_deg'], day_s)
                                     for k in (name, *others)}, gcm=gcm_surface_energy(site['lat_deg'], site['lon_deg'])),
                microphysics={k: microphysics_budget(RUNS / k, from_day) for k in (name, *([reference] if reference else []))
                              if (RUNS / k / 'cm1out_stats.dat').exists()},
                **({'reference': reference} if reference else {}))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case')
    parser.add_argument('--from-day', type=float, default=29.5)
    parser.add_argument('--reference', help='another box at the same site to compare with, column for column')
    parser.add_argument('--output', type=Path, help='where to write the comparison (default: the results folder)')
    args = parser.parse_args(argv)
    result = analyse(args.case, args.from_day, args.reference)
    result['producer'] = dict(domain='climate', files={f'crm/{Path(f).name}': hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]
                                                       for f in (__file__, rc.__file__)})
    path = args.output or RESULTS / f'box_{args.case}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rc.rounded(result), indent=1) + '\n')
    print(f'{path}: {result["snapshots"]} snapshots, days {result["span_days"][0]:.1f}-{result["span_days"][1]:.1f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
