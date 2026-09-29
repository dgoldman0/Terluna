"""The GCM's clear-sky sunlight against the repository's line-by-line model, column by column.

    climate/gcm/.venv/bin/python -m climate.gcm.radiation_check A28_dim5 --years 15 24

writes ../results/gcm/radiation_check_<run>.json. Needs the line-by-line model's spectroscopic inputs
(python -m atmosphere.radiative_convective.fetch_inputs --download) and takes about half an hour on 8 cores.

The GCM's radiation was calibrated by matching its clear-sky planetary albedo to line-by-line on the global-mean
column (README, "Radiation, checked against line-by-line and calibrated"). That leaves open how it divides the
rest of the sunlight between the air and the ground, and how either varies from column to column. Here each chosen
column's own sunlight, clear sky, is split three ways over the lunar day: reflected to space, absorbed in the air
and absorbed by the ground, as the GCM's own clear-sky diagnostics have it and as line-by-line has it on the GCM's
mean profile of that column, its surface albedo and sunlight through the titania stack. The columns span the design
case's land from its driest to its most humid and from lowland to highland, with seas from the equator to 60
degrees.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-radiation-check/1'
EVIDENCE = ('Clear-sky solar fluxes only. The GCM\'s are its own clear-sky diagnostics (3-day means over the years). '
            'The line-by-line model (HITRAN lines, MT_CKD water continuum, collision-induced absorption, Rayleigh '
            'scattering, a delta-scaled two-stream solver with a pseudo-spherical direct beam) runs on each column\'s '
            'time-mean profile on the GCM\'s ten levels, extended isothermal and dry above its top, with the column\'s '
            'mean surface albedo and the Sun on the equator (the Moon\'s 1.5 degree tilt left out), so it leaves out how '
            'the fluxes vary with the profile from day to day and season to season.')
READING_RULE = ('Shares are of the sunlight arriving at the top of each column over the lunar day, clear sky. The GCM\'s '
                'arriving sunlight is its own time mean, which its uneven Sun makes up to 12% less or more than it '
                'should be depending on longitude; the shares compare how each model divides whatever arrives. '
                'difference is the GCM minus line-by-line, in points of the arriving sunlight.')


def mean_fields(run: str, years) -> dict:
    """Time means over the years of what the check needs, on the GCM grid."""
    import netCDF4
    from climate.crm.cm1_run import planet
    folder = HERE / 'runs' / run / 'model'
    keys = ('ta', 'hus', 'ps', 'tas', 'prw', 'alb', 'rst', 'rsut', 'rstcs', 'rsscs')
    sums, count = {}, 0
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            for k in keys:
                sums[k] = sums.get(k, 0.0) + np.asarray(d[k][:], float).sum(axis=0)
            count += d['tas'].shape[0]
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
            sigma = np.asarray(d['lev'][:], float)
    f = {k: v / count for k, v in sums.items()}
    height = np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel().reshape(lsm.shape)
    f.update(lat=lat, lon=lon, lsm=lsm, sigma=sigma, height_m=height / planet()['gravity_m_s2'],
             arriving=f['rst'] - f['rsut'])                                       # rsut is negative, upward
    return f


def choose_columns(f: dict, site=(0.0, 0.0)) -> list:
    """Land from the driest to the most humid (water vapour at its 5th, 25th, 50th, 75th and 95th percentiles), a
    tropical highland and the land nearest a site, within 60 degrees of the equator; sea-level seas nearest the
    equator, 25, 45 and 60 degrees, each the one of median water vapour in its band."""
    lat2 = np.abs(f['lat'])[:, None] * np.ones(f['lon'].size)[None, :]
    land = f['lsm'] & (lat2 <= 60.0)
    sea = ~f['lsm'] & (f['height_m'] < 1.0)
    picks, taken = [], set()

    def add(label, j, i):
        if (j, i) not in taken:
            taken.add((j, i))
            picks.append(dict(label=label, j=int(j), i=int(i)))

    j0 = int(np.argmin(np.abs(f['lat'] - site[0])))
    i0 = int(np.argmin(np.abs((f['lon'] - site[1] + 180.0) % 360.0 - 180.0)))
    add('site land', j0, i0) if f['lsm'][j0, i0] else None
    for q in (5, 25, 50, 75, 95):
        target = np.percentile(f['prw'][land], q)
        j, i = np.unravel_index(np.argmin(np.where(land, np.abs(f['prw'] - target), np.inf)), land.shape)
        add(f'land, water vapour at its {q}th percentile', j, i)
    tropical = land & (lat2 <= 30.0)
    j, i = np.unravel_index(np.argmax(np.where(tropical, f['height_m'], -np.inf)), land.shape)
    add('tropical highland', j, i)
    for band in (0.0, 25.0, 45.0, 60.0):
        rows = np.abs(lat2 - band) <= 5.0 if band else lat2 <= 5.0
        cells = sea & rows
        if cells.any():
            target = np.median(f['prw'][cells])
            j, i = np.unravel_index(np.argmin(np.where(cells, np.abs(f['prw'] - target), np.inf)), cells.shape)
            add(f'sea near {band:.0f} degrees', j, i)
    return picks


def column_profile(f: dict, j: int, i: int) -> dict:
    """A column's time-mean levels from the ground up, extended above the GCM's top."""
    ps = float(f['ps'][j, i]) * 100.0                                             # hPa in the output
    above = np.exp(np.linspace(np.log(f['sigma'][0] * ps * 0.8), np.log(50.0), 12))
    t = f['ta'][:, j, i]
    return dict(p=np.concatenate([f['sigma'] * ps, above]), t=np.concatenate([t, np.full(above.size, t[0])]),
                qv=np.concatenate([f['hus'][:, j, i], np.full(above.size, 1e-7)]), ps=ps, t2=float(f['tas'][j, i]))


def check(run: str, years, processes: int = 8, site=(0.0, 0.0), only: str | None = None) -> dict:
    from climate.crm.radiation_check import line_by_line_column, solar_parts, daily_mean
    f = mean_fields(run, years)
    nodes, weights = np.polynomial.legendre.leggauss(6)
    hours, hour_weights = 0.25 * np.pi * (nodes + 1.0), 0.25 * np.pi * weights
    rows = []
    for pick in choose_columns(f, site):
        if only and pick['label'] != only:
            continue
        j, i = pick['j'], pick['i']
        lat, albedo = float(f['lat'][j]), float(f['alb'][j, i])
        prof = column_profile(f, j, i)
        col, water = line_by_line_column(prof['p'], prof['t'], prof['qv'], prof['ps'], prof['t2'])
        parts = solar_parts(col, np.cos(np.radians(lat)) * np.cos(hours), processes, albedo, shields=('titania',))
        lbl = daily_mean(parts['titania'], hours, hour_weights, 1.0, lat)
        total = sum(lbl.values())
        lbl = {k: v / total for k, v in lbl.items()}
        arriving = float(f['arriving'][j, i])
        own = dict(reflected=arriving - float(f['rstcs'][j, i]), air=float(f['rstcs'][j, i] - f['rsscs'][j, i]),
                   ground=float(f['rsscs'][j, i]))
        own = {k: v / arriving for k, v in own.items()}
        rows.append(dict(pick, lat_deg=lat, lon_deg=float(f['lon'][i]), land=bool(f['lsm'][j, i]),
                         height_m=float(f['height_m'][j, i]), surface_albedo=albedo, water_kg_m2=float(f['prw'][j, i]),
                         line_by_line_water_kg_m2=water, arriving_w_m2=arriving, gcm=own, line_by_line=lbl,
                         difference_points={k: 100.0 * (own[k] - lbl[k]) for k in own}))
        print(f"{pick['label']:40s} {lat:6.1f} {f['lon'][i]:6.1f}  water {f['prw'][j, i]:4.0f}  albedo {albedo:.3f}  "
              f"GCM {own['reflected']:.3f} {own['air']:.3f} {own['ground']:.3f}  line-by-line {lbl['reflected']:.3f} "
              f"{lbl['air']:.3f} {lbl['ground']:.3f}", flush=True)
    mean = lambda sel, k: float(np.mean([r['difference_points'][k] for r in rows if sel(r)])) if any(sel(r) for r in rows) else None
    summary = {surface: {k: mean(sel, k) for k in ('reflected', 'air', 'ground')}
               for surface, sel in (('land', lambda r: r['land']), ('sea', lambda r: not r['land']))}
    return dict(schema=SCHEMA, run=run, years=list(years), evidence=EVIDENCE, reading_rule=READING_RULE,
                columns=rows, mean_difference_points=summary)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('--years', type=int, nargs=2, default=(15, 24))
    parser.add_argument('--processes', type=int, default=8)
    parser.add_argument('--only', help="one column's label, e.g. 'site land'; the result is named after it")
    args = parser.parse_args(argv)
    result = check(args.run, tuple(args.years), args.processes, only=args.only)
    result['producer'] = dict(domain='climate', files={f'{Path(p).parent.name}/{Path(p).name}': hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
                                                       for p in (__file__, HERE.parent / 'crm' / 'radiation_check.py')})
    from climate.crm.ring_crossings import rounded
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"radiation_check_{args.run}{'_' + args.only.split()[0] if args.only else ''}.json"
    path.write_text(json.dumps(rounded(result, 4), indent=1) + '\n')
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
