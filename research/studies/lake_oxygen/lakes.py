"""The study's lakes: the rain-fed lakes' depths, and a crater lake near their median depth in each latitude band.

    python -m research.studies.lake_oxygen.lakes        # a few seconds; reads LOLA and the drainage product

Geography's drainage product holds the rain-fed lakes on LOLA's 16 pixel/degree grid, each lake cell with its surface
level above the geoid ([geography](../../../geography/README.md#rivers-and-rain-fed-lakes-a-first-estimate)). Their
depths are those levels less LOLA's heights above the GRAIL geoid, as the drainage summary takes them, and the lakes
are its connected bodies (eight neighbours). This writes results/lakes.json: how the lakes' area divides by local
depth, by latitude and by each lake's deepest point; and for each band of 0-15, 15-30 and 30-45 degrees the crater
lake whose deepest point lies nearest the median local depth, with its hypsography, its level and the IAU crater it
fills. The lake model takes each chosen lake's hypsography scaled to the median depth.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DRAINAGE = ROOT / 'geography' / 'products' / 'drainage_28pct_16ppd.npz'
DRAINAGE_SUMMARY = ROOT / 'geography' / 'results' / 'drainage.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
OUT = HERE / 'results' / 'lakes.json'
SCHEMA = 'terluna.research.lake-oxygen-lakes/1'
PPD = 16
BANDS = ((0, 15), (15, 30), (30, 45))
DEPTH_WINDOW = (0.85, 1.15)         # a chosen lake's deepest point lies within 15% of the median local depth
MIN_AREA_KM2 = 300.0                # and it is large enough to carry its own wind and fetch
HYPSOGRAPHY_POINTS = 41


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def lake_depths():
    """Every lake cell's depth (m), area (m2), latitude, longitude and level above the geoid (m)."""
    from scipy.interpolate import RegularGridInterpolator
    from geography import topography as tp
    z = np.load(DRAINAGE)
    cells = z['lake_cell'].astype(np.int64)
    level = z['lake_level_m'].astype(float)
    atlas = json.loads(ATLAS.read_text())
    radius_h, grid = tp.read_ldem(PPD)
    _, coarse = tp.read_ldem(4)
    n4 = tp.geoid(coarse, atlas['grid']['geoid_degree'])
    # the geoid on the 0.25-degree grid interpolated to the lake cells, as topography.height_above_geoid does
    lat_c = np.concatenate([[90.0], coarse.lat_deg, [-90.0]])
    pad = np.vstack([n4[:1].mean(1, keepdims=True).repeat(n4.shape[1], 1), n4,
                     n4[-1:].mean(1, keepdims=True).repeat(n4.shape[1], 1)])
    lon_c = np.concatenate([coarse.lon_deg - 360.0, coarse.lon_deg, coarse.lon_deg + 360.0])
    geoid = RegularGridInterpolator((lat_c[::-1], lon_c), np.hstack([pad, pad, pad])[::-1])
    cols = 360 * PPD
    lat, lon = grid.lat_deg[cells // cols], grid.lon_deg[cells % cols]
    height = radius_h.ravel()[cells] - geoid(np.stack([lat, lon], 1))
    area = grid.cell_area_m2.ravel()[cells]
    return dict(cells=cells, depth=level - height, area=area, lat=lat, lon=lon, level=level,
                sea_level=float(atlas['sea_level_m']))


def bodies(d):
    """Connected lake bodies (eight neighbours) and each body's area, volume, deepest point and centre."""
    from scipy import ndimage
    rows, cols = 180 * PPD, 360 * PPD
    mask = np.zeros(rows * cols, bool)
    mask[d['cells']] = True
    labels, n = ndimage.label(mask.reshape(rows, cols), structure=np.ones((3, 3), bool))
    body = labels.ravel()[d['cells']]
    a = d['area']
    area = np.bincount(body, a, n + 1)
    volume = np.bincount(body, a * d['depth'], n + 1)
    deepest = np.zeros(n + 1)
    np.maximum.at(deepest, body, d['depth'])
    lat = np.bincount(body, a * d['lat'], n + 1) / np.maximum(area, 1.0)
    x = np.bincount(body, a * np.cos(np.radians(d['lon'])), n + 1)
    y = np.bincount(body, a * np.sin(np.radians(d['lon'])), n + 1)
    level = np.full(n + 1, -np.inf)
    np.maximum.at(level, body, d['level'])
    return dict(body=body, n=int(n), area=area, volume=volume, deepest=deepest, lat=lat,
                lon=np.degrees(np.arctan2(y, x)) % 360.0, level=level)


def weighted_quantile(x, w, q):
    i = np.argsort(x, kind='stable')
    c = np.cumsum(w[i])
    return float(np.interp(q * c[-1], c - 0.5 * w[i], x[i]))


def crater_of(lat, lon, craters):
    """The smallest named crater whose rim (half its diameter) contains the point, or None."""
    best = None
    for c in craters:
        dlat = np.radians(lat - c['lat'])
        dlon = np.radians(((lon - c['lon'] + 180.0) % 360.0) - 180.0)
        h = np.sin(dlat / 2) ** 2 + np.cos(np.radians(lat)) * np.cos(np.radians(c['lat'])) * np.sin(dlon / 2) ** 2
        km = 2.0 * 1737.4 * np.arcsin(np.sqrt(min(1.0, h)))
        if km <= 0.5 * c['diameter_km'] and (best is None or c['diameter_km'] < best['diameter_km']):
            best = dict(name=c['name'], diameter_km=round(c['diameter_km'], 1), centre_km=round(float(km), 1))
    return best


def hypsography(depth, area, deepest):
    """The lake's area below each depth, as fractions of its surface area at fractions of its deepest point."""
    f = np.linspace(0.0, 1.0, HYPSOGRAPHY_POINTS)
    order = np.argsort(depth)
    d, a = depth[order], area[order]
    below = a.sum() - np.concatenate([[0.0], np.cumsum(a)])      # area deeper than each sorted depth
    frac = np.interp(f * deepest, np.concatenate([[0.0], d]), below) / a.sum()
    frac[0], frac[-1] = 1.0, 0.0
    return [round(float(x), 4) for x in f], [round(float(x), 5) for x in np.minimum.accumulate(frac)]


def select(d, b, median_depth, craters):
    chosen = []
    for lo, hi in BANDS:
        ok = ((np.abs(b['lat']) >= lo) & (np.abs(b['lat']) < hi) & (b['area'] / 1e6 >= MIN_AREA_KM2)
              & (b['deepest'] >= DEPTH_WINDOW[0] * median_depth) & (b['deepest'] <= DEPTH_WINDOW[1] * median_depth))
        ok[0] = False
        candidates = np.flatnonzero(ok)
        named = [(k, crater_of(b['lat'][k], b['lon'][k], craters)) for k in candidates]
        pool = [(k, c) for k, c in named if c] or [(k, None) for k in candidates]
        k, crater = min(pool, key=lambda kc: abs(b['deepest'][kc[0]] - median_depth))
        cells = b['body'] == k
        f, area_fraction = hypsography(d['depth'][cells], d['area'][cells], b['deepest'][k])
        chosen.append(dict(band=f'{lo}-{hi}', candidates=int(candidates.size),
                           named_candidates=sum(1 for _, c in named if c),
                           lat=round(float(b['lat'][k]), 3), lon_east=round(float(b['lon'][k]), 3),
                           area_km2=round(float(b['area'][k]) / 1e6, 1),
                           mean_depth_m=round(float(b['volume'][k] / b['area'][k]), 1),
                           deepest_m=round(float(b['deepest'][k]), 1),
                           level_above_sea_m=round(float(b['level'][k] - d['sea_level']), 1),
                           crater=crater, hypsography=dict(depth_fraction=f, area_fraction=area_fraction)))
    return chosen


def summary(d, b):
    a, depth = d['area'], d['depth']
    median = weighted_quantile(depth, np.ones_like(depth), 0.5)        # by cell, as the drainage summary counts it
    lake_area = b['area'][1:] / 1e6
    deepest = b['deepest'][1:]
    alat = np.abs(b['lat'][1:])
    out = dict(lake_cells=int(depth.size), bodies=b['n'], median_local_depth_m=round(median, 1),
               median_local_depth_by_area_m=round(weighted_quantile(depth, a, 0.5), 1),
               local_depth_quantiles_by_area_m={str(q): round(weighted_quantile(depth, a, q), 1)
                                                for q in (0.1, 0.25, 0.75, 0.9)},
               area_share_deeper_than_median=round(float(a[depth > median].sum() / a.sum()), 4),
               area_share_in_lakes_deepest_beyond_median=round(float(lake_area[deepest > median].sum()
                                                                    / lake_area.sum()), 4),
               lake_area_km2=round(float(lake_area.sum())),
               bands={})
    for lo, hi in BANDS + ((45, 90),):
        m = (alat >= lo) & (alat < hi)
        out['bands'][f'{lo}-{hi}'] = dict(bodies=int(m.sum()),
                                          area_share=round(float(lake_area[m].sum() / lake_area.sum()), 4))
    for lo, hi in ((0, 1), (1, 10), (10, 100), (100, 1000), (1000, 1e4), (1e4, 1e7)):
        m = (lake_area >= lo) & (lake_area < hi)
        out.setdefault('by_size_km2', {})[f'{lo:g}-{hi:g}'] = dict(
            bodies=int(m.sum()), area_share=round(float(lake_area[m].sum() / lake_area.sum()), 4),
            median_deepest_m=round(float(np.median(deepest[m])), 1) if m.any() else None)
    return out, median


def main(argv=None) -> int:
    from geography import nomenclature
    from geography import fetch_inputs
    d = lake_depths()
    b = bodies(d)
    stats, median = summary(d, b)
    reference = json.loads(DRAINAGE_SUMMARY.read_text())['lakes']
    if abs(stats['median_local_depth_m'] - reference['median_depth_m']) > 1.0 or stats['bodies'] != reference['bodies']:
        raise SystemExit(f"Reproduced {stats['median_local_depth_m']} m and {stats['bodies']} bodies; the drainage "
                         f"summary has {reference['median_depth_m']} m and {reference['bodies']}")
    craters = [f for f in nomenclature.features() if f['code'] == 'AA']
    lakes = select(d, b, median, craters)
    inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (DRAINAGE, DRAINAGE_SUMMARY, ATLAS)}
    for name in ('ldem_16.img', 'ldem_4.img', 'jggrx_0420a_sha.tab', nomenclature.ARCHIVE):
        inputs[f'geography/inputs/{name}'] = sha256(fetch_inputs.path(name))
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(Path(__file__).resolve().relative_to(ROOT)): sha256(__file__)}, inputs=inputs),
        evidence=('Hydrostatic geometry: the drainage product\'s rain-fed lakes (a first estimate from the climate run '
                  'A28_dim5_moon, before the GCM\'s corrections; the lakes are regenerated once the design case '
                  'settles) on LOLA at 16 pixels a degree above the GRAIL geoid. Depths are lake levels less ground '
                  'heights at 1.9 km cells; a lake\'s deepest point is its deepest cell.'),
        reading_rule=('Depths in m below each lake\'s surface. The median local depth counts lake cells as the '
                      'drainage summary does; the by-area figures weight cells by area. A chosen lake\'s hypsography '
                      'gives the share of its surface area lying deeper than each fraction of its deepest point. '
                      'Levels are above the seas\' surface, the atlas\'s sea level.'),
        lakes=stats, selection=dict(bands=[f'{lo}-{hi}' for lo, hi in BANDS], depth_window=list(DEPTH_WINDOW),
                                    min_area_km2=MIN_AREA_KM2,
                                    rule=('in each band, among lakes of at least min_area_km2 whose deepest point '
                                          'lies within the window of the median local depth, the one inside a named '
                                          'IAU crater whose deepest point lies nearest the median')),
        chosen=lakes)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    for k in lakes:
        print(k['band'], k['lat'], k['lon_east'], k['area_km2'], 'km2', k['mean_depth_m'], k['deepest_m'], 'm',
              k['level_above_sea_m'], 'm up', k['crater'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
