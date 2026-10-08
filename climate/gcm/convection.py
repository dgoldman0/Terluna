"""Convective rain of an ExoPlaSim run over the whole Moon, written as a data product for other lanes.

    climate/gcm/.venv/bin/python -m climate.gcm.convection A28_dim5_moon:20-29
    # -> results/gcm/convection_A28_dim5_moon.json

For the chosen model years it maps, cell by cell, the mean of the convection scheme's rain (prc) and of the total
rain (pr), the means of prc to the powers 1.5 and 2 (so a consumer can weight a nonlinear process by the same
output), the share of outputs in which the scheme rains at all, and the land share from the model's land-sea mask.
Each value is a mean over the 3-day outputs at T21, so single storms are averaged out; the moments carry how the
3-day means are spread. Rain is in mm per day of water.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from climate.gcm.site_winds import RESULTS

SCHEMA = 'terluna.climate.gcm-convection/1'
RAINING_MM_DAY = 0.1          # a 3-day mean above this counts as convective rain falling
EVIDENCE = ('ExoPlaSim 3-day output means at T21 (about 170 km cells) over the whole Moon. prc is the convection '
            'scheme\'s rain and pr all rain; the 3-day means average single storms out, so the maps describe where '
            'and how much the scheme convects, which a storm model must translate into storms.')
READING_RULE = ('Maps run [lat][lon] over the model grid (lat_deg north to south, lon_deg east from 0). Rain in mm '
                'per day of liquid water. prc_p15_mean and prc_p2_mean are the means of prc**1.5 and prc**2 with prc '
                'in mm/day. raining_share is the share of 3-day outputs with prc above raining_mm_day. land_share is '
                'the mean of the model\'s land-sea mask (1 land). cell_area_share sums to 1 over the Moon.')


def area_shares(lat_deg, nlon):
    """Each cell's share of the sphere, from the Gaussian latitudes' midpoints."""
    lat = np.radians(np.asarray(lat_deg, dtype=float))
    edges = np.concatenate([[np.pi / 2], 0.5 * (lat[:-1] + lat[1:]), [-np.pi / 2]])
    band = 0.5 * (np.sin(edges[:-1]) - np.sin(edges[1:]))
    return np.repeat(band[:, None] / nlon, nlon, axis=1)


def convection(folder: str, first: int, last: int) -> tuple[dict, dict]:
    import netCDF4
    from climate.gcm.exoplasim_run import RUNS
    sums, n, inputs, lat, lon = {}, 0, {}, None, None
    for year in range(first, last + 1):
        path = RUNS / folder / 'model' / f'MOST.{year:05d}.nc'
        inputs[f'climate/gcm/runs/{folder}/model/{path.name}'] = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
        with netCDF4.Dataset(path) as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            prc = np.asarray(d['prc'][:], float) * 86400e3
            pr = np.asarray(d['pr'][:], float) * 86400e3
            lsm = np.asarray(d['lsm'][:], float)
        prc = np.clip(prc, 0.0, None)
        for key, value in (('prc', prc), ('pr', np.clip(pr, 0.0, None)), ('prc_p15', prc ** 1.5),
                           ('prc_p2', prc ** 2), ('raining', (prc > RAINING_MM_DAY).astype(float)), ('land', lsm)):
            sums[key] = sums.get(key, 0.0) + value.sum(axis=0)
        n += prc.shape[0]
    mean = {k: v / n for k, v in sums.items()}
    share = area_shares(lat, lon.size)
    r = lambda a, d=4: np.round(np.asarray(a, float), d).tolist()
    land = mean['land'] > 0.5
    out = dict(lat_deg=r(lat, 3), lon_deg=r(lon, 3), outputs=n, raining_mm_day=RAINING_MM_DAY,
               cell_area_share=r(share, 7), prc_mean=r(mean['prc']), pr_mean=r(mean['pr']),
               prc_p15_mean=r(mean['prc_p15']), prc_p2_mean=r(mean['prc_p2']), raining_share=r(mean['raining'], 3),
               land_share=r(mean['land'], 3),
               moon_mean=dict(prc=round(float((share * mean['prc']).sum()), 4),
                              pr=round(float((share * mean['pr']).sum()), 4),
                              prc_land=round(float((share * mean['prc'] * land).sum() / (share * land).sum()), 4),
                              prc_sea=round(float((share * mean['prc'] * ~land).sum() / (share * ~land).sum()), 4),
                              land_area_share=round(float((share * land).sum()), 4)))
    return out, inputs


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    if len(args) != 1 or ':' not in args[0]:
        raise SystemExit('Give a run and year range, for example A28_dim5_moon:20-29')
    folder, span = args[0].split(':')
    first, last = (int(v) for v in span.split('-'))
    from climate.gcm.exoplasim_run import RUNS
    progress = json.loads((RUNS / folder / 'progress.json').read_text())
    maps, inputs = convection(folder, first, last)
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    product = dict(schema=SCHEMA, run=folder, years=[first, last], configuration=progress['configuration'],
                   producer=dict(domain='climate', files={'gcm/convection.py': digest(__file__)}, inputs=inputs),
                   evidence=EVIDENCE, reading_rule=READING_RULE, **maps)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / f'convection_{folder}.json'
    out.write_text(json.dumps(product, indent=1) + '\n')
    m = maps['moon_mean']
    print(f"{folder} years {first}-{last}: convective rain {m['prc']} mm/day over the Moon (land {m['prc_land']}, "
          f"sea {m['prc_sea']}), all rain {m['pr']}; land {m['land_area_share']:.0%} of the area")
    print(out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
