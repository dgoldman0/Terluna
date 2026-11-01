"""What the storm-water study reads from the design years: the drylands' cells of the design run
A28R_301k_design_l15_snow_s934, years 10-19, with every 3-day mean of the fields their storms, soils and pools need.

    /home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.storm_water.extract   # about a minute

The cells are classed as session-tools/design_run/design_years.py classes them: each cell's mean over the ten tenths of
the lunar day (each tenth the mean of the 3-day means whose hour angle falls in it, from climate.crm.gcm_site's
subsolar longitude), its rain over Priestley and Taylor's demand on that mean's net radiation, air and pressure, in
UNEP's aridity classes. The drylands are the arid, semi-arid and dry sub-humid cells, aridity 0.05-0.65. For each of
them the extract keeps all 1,200 means of the rain and its convective part, the evaporation, the bucket's soil water,
the net sunlight and longwave at the ground, the air at 2 m, the ground, the surface pressure, the wind at 10 m and the
sensible and latent heat, with the tenth each mean falls in. It writes research/runs/storm_water/design_years.npz
(ignored) and a sidecar design_years.json with the sha256 and size of every year file read. A change of design years
reruns it in about a minute.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from climate.crm import gcm_site as GS
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = Path('/media/projectspace/terluna-research/gcm-runs')
RUN, YEARS = 'A28R_301k_design_l15_snow_s934', (10, 19)
OUT = ROOT / 'research' / 'runs' / 'storm_water' / 'design_years.npz'
SIDECAR = OUT.with_suffix('.json')
SCHEMA = 'terluna.research.storm-water-design-years/1'
SERIES = ('pr', 'prc', 'evap', 'mrso', 'rss', 'rls', 'tas', 'ts', 'ps', 'uas', 'vas', 'hfss', 'hfls')
# design_years.py's demand: Priestley and Taylor's alpha, water's latent heat (J/kg) and air's heat capacity (J/kg/K)
PT_ALPHA, LATENT, CP = 1.26, 2.5e6, 1004.0
UNEP = (('hyper-arid', 0.0, 0.05), ('arid', 0.05, 0.2), ('semi-arid', 0.2, 0.5), ('dry sub-humid', 0.5, 0.65),
        ('humid', 0.65, np.inf))
DRYLANDS = ('arid', 'semi-arid', 'dry sub-humid')
MEAN_S = 86400.0                                     # PlaSim's rates are m/s of water


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def tenths(czen, lat, lon):
    """The tenth of the lunar day (0-9, from midnight through noon at 5) each 3-day mean falls in, by longitude."""
    hour = (lon[None, :] - GS.subsolar_longitude(czen, lat, lon)[:, None] + 180.0) % 360.0 - 180.0
    return (np.floor((hour + 180.0) / 36.0).astype(int) % 10).astype(np.int8)          # (time, lon)


def demand(tas, ps, net):
    """Priestley and Taylor's demand (m/s of water) on mean net radiation (W/m2), air at 2 m (K) and pressure (hPa)."""
    t = tas - 273.15
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
    s = 4098.0 * es / (t + 237.3) ** 2
    gamma = CP * ps * 100.0 / 1000.0 / (0.622 * LATENT)
    return np.maximum(PT_ALPHA * s / (s + gamma) * net / LATENT, 1e-3 / 86400.0) / 1000.0


def cell_areas(lat, nlon):
    """Each grid cell's area (km2): the band between the midpoints of the Gaussian latitudes, poles at the ends."""
    edges = np.concatenate(([90.0], 0.5 * (lat[:-1] + lat[1:]), [-90.0]))
    band = np.abs(np.sin(np.radians(edges[:-1])) - np.sin(np.radians(edges[1:])))
    return (MOON_RADIUS / 1e3) ** 2 * (2.0 * np.pi / nlon) * band[:, None] * np.ones((1, nlon))


def read():
    import netCDF4
    files, series, tenth, lat, lon, land = {}, {k: [] for k in SERIES}, [], None, None, None
    for year in range(YEARS[0], YEARS[1] + 1):
        path = RUNS / RUN / 'model' / f'MOST.{year:05d}.nc'
        files[str(path)] = dict(sha256=sha256(path), bytes=path.stat().st_size)
        with netCDF4.Dataset(path) as d:
            la, lo = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            mask = np.asarray(d['lsm'][:], float) > 0.5
            if lat is not None and not (np.array_equal(la, lat) and np.array_equal(mask[0], land)):
                raise ValueError(f'{path}: the grid or land mask changes')
            if not (mask == mask[0]).all():
                raise ValueError(f'{path}: the land mask changes within the year')
            lat, lon, land = la, lo, mask[0]
            for k in SERIES:
                series[k].append(np.asarray(d[k][:], np.float64))
            tenth.append(tenths(np.asarray(d['czen'][:], float), lat, lon))
    return files, {k: np.concatenate(v) for k, v in series.items()}, np.concatenate(tenth), lat, lon, land


def classify(f, tenth, land):
    """Each cell's yearly means as design_years.py forms them, its demand, aridity and UNEP class."""
    by = np.zeros((10,) + land.shape + (len(SERIES),))
    cnt = np.zeros((10,) + land.shape)
    stack = np.stack([f[k] for k in SERIES], axis=-1)                                # (time, lat, lon, field)
    for k in range(10):
        sel = (tenth == k)[:, None, :]                                               # (time, 1, lon)
        by[k] = np.where(sel[..., None], stack, 0.0).sum(axis=0)
        cnt[k] = np.broadcast_to(sel, sel.shape[:1] + land.shape).sum(axis=0)
    by /= np.maximum(cnt, 1)[..., None]
    mean = {k: by[..., n].mean(axis=0) for n, k in enumerate(SERIES)}
    pet = demand(mean['tas'], mean['ps'], mean['rss'] + mean['rls'])
    ai = mean['pr'] / pet
    cls = np.full(land.shape, '', dtype=object)
    for name, lo, hi in UNEP:
        cls[land & (ai >= lo) & (ai < hi)] = name
    return mean, pet, ai, cls


def extract():
    files, f, tenth, lat, lon, land = read()
    mean, pet, ai, cls = classify(f, tenth, land)
    area = cell_areas(lat, lon.size)
    jj, ii = np.nonzero(np.isin(cls, DRYLANDS))
    arrays = {k: f[k][:, jj, ii].astype(np.float32) for k in SERIES}
    arrays.update(tenth=tenth[:, ii], lat=lat[jj], lon=lon[ii], j=jj, i=ii, area_km2=area[jj, ii],
                  aridity=ai[jj, ii], demand_m_s=pet[jj, ii], cls=np.array([str(c) for c in cls[jj, ii]]),
                  grid_lat=lat, grid_lon=lon, land=land)
    wl = area * land
    shares = {name: float(wl[cls == name].sum() / wl.sum()) for name, _, _ in UNEP}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT, **arrays)
    meta = dict(schema=SCHEMA, run=RUN, years=list(YEARS), means=int(f['pr'].shape[0]), cells=int(jj.size),
                land_shares=shares, files=files, npz=dict(sha256=sha256(OUT), bytes=OUT.stat().st_size),
                producer=dict(study='research/studies/storm_water', file='research/studies/storm_water/extract.py',
                              sha256=sha256(Path(__file__))),
                rule='UNEP classes of rain over Priestley-Taylor demand (alpha 1.26) on each cell\'s mean over the '
                     'ten tenths of the lunar day, as session-tools/design_run/design_years.py reads them')
    SIDECAR.write_text(json.dumps(meta, indent=1) + '\n')
    return meta


def composite(x):
    """Each dryland cell's means by tenth of the lunar day and over the ten tenths, as design_years.py composites
    them, for every field the extract keeps: ({field: (tenth, cell)}, {field: (cell,)}), in the year files' units."""
    tenth = x['tenth']
    count = np.stack([(tenth == k).sum(axis=0) for k in range(10)])
    by = {k: np.stack([(x[k].astype(float) * (tenth == t)).sum(axis=0) for t in range(10)]) / np.maximum(count, 1)
          for k in SERIES}
    return by, {k: v.mean(axis=0) for k, v in by.items()}


def load():
    """The extract and its sidecar, checked against each other."""
    meta = json.loads(SIDECAR.read_text())
    if meta.get('schema') != SCHEMA or meta['npz']['sha256'] != sha256(OUT):
        raise ValueError('the design-years extract does not match its sidecar; rerun extract.py')
    with np.load(OUT, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}, meta


if __name__ == '__main__':
    m = extract()
    print(f"{m['cells']} dryland cells over {m['means']} means; land shares " +
          ', '.join(f'{k} {v:.3f}' for k, v in m['land_shares'].items()), file=sys.stderr)
