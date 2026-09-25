"""The Open Moon atlas at the selected standing-water share, written as a data product.

    python -m geography.atlas             # results/atlas.json and products/atlas_<pct>pct_4ppd.npz
    python -m geography.atlas --share 0.3 # another share, for comparison

One common water level above the geoid covers the scenario's share of the Moon
(shared/scenarios/water.json). For that level the atlas lists the separate seas and
lakes with their depths and the IAU water names inside them, the islands with their
summits, how much of each named mare floods, and a comparison of candidate shares.
The filling is hydrostatic, as in water_inventory.py; rainfall, runoff, evaporation,
groundwater and crustal loading then move each shoreline.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from shared.constants import MOON_RADIUS
from geography import topography as tp, basins as bs, fetch_inputs, nomenclature as nm

HERE = Path(__file__).resolve().parent
SCENARIO = HERE.parent / 'shared' / 'scenarios' / 'water.json'
SCHEMA = 'terluna.geography.atlas/1'
GRID_SCHEMA = 'terluna.geography.atlas-grid/1'
EVIDENCE = ('Hydrostatic storage geometry: LOLA topography (4 pixels per degree, polar caps at 16) above the GRAIL '
            'geoid (degree 200, with static Earth tide and rotation), filled to one common water level covering the '
            'scenario share. Every closed depression below that level counts as water. Shorelines and lake levels '
            'follow once rainfall, runoff, evaporation, groundwater and crustal loading are added. Names come from '
            'the IAU gazetteer; names for whole seas and islands are working labels.')
READING_RULE = ('Heights and levels are metres above the geoid (a zero-mean equipotential). Longitudes are east. '
                'Areas are exact spherical sums; depths are level minus ground. A body or island is one 8-connected '
                'region, joined across the 0/360 seam and the poles.')
COMPARE_SHARES = (0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.32, 0.35)
NEAR_SIDE_MARIA = ('Mare Imbrium', 'Mare Serenitatis', 'Mare Tranquillitatis', 'Mare Crisium', 'Mare Fecunditatis',
                   'Mare Nectaris', 'Mare Nubium', 'Mare Humorum', 'Mare Cognitum', 'Mare Insularum', 'Mare Vaporum',
                   'Oceanus Procellarum')
BODY_MIN_SHARE = 5e-4        # smallest sea or lake listed by name
ISLAND_MIN_KM2 = 400.0       # smallest island listed
ISLAND_NAME_RADIUS_KM = 60.0
MARE_CAP = 0.6               # inner share of each mare's radius used for its flooding fraction


def scenario_share(path: Path = SCENARIO) -> float:
    return float(json.loads(Path(path).read_text())['surface_water_share'])


def producer():
    files = ['topography.py', 'basins.py', 'nomenclature.py', 'atlas.py']
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    return dict(domain='geography', files={f: digest(HERE / f) for f in files},
                inputs_manifest=digest(HERE / 'inputs.json'), scenario=digest(SCENARIO))


def great_circle_km(lat0, lon0, lat, lon):
    p0, l0, p, l = map(np.radians, (lat0, lon0, lat, lon))
    c = np.sin(p0) * np.sin(p) + np.cos(p0) * np.cos(p) * np.cos(l - l0)
    return MOON_RADIUS / 1e3 * np.arccos(np.clip(c, -1.0, 1.0))


def cell_index(grid, lat, lon):
    i = np.clip(np.round((90.0 - np.asarray(lat)) * grid.pixels_per_degree - 0.5).astype(int), 0, grid.lat_deg.size - 1)
    j = np.round((np.asarray(lon) % 360.0) * grid.pixels_per_degree - 0.5).astype(int) % grid.lon_deg.size
    return i, j


def _centroids(labels, area, grid, n):
    la, lo = np.meshgrid(np.radians(grid.lat_deg), np.radians(grid.lon_deg), indexing='ij')
    lab, w = labels.ravel(), area.ravel()
    x = np.bincount(lab, w * (np.cos(la) * np.cos(lo)).ravel(), n)
    y = np.bincount(lab, w * (np.cos(la) * np.sin(lo)).ravel(), n)
    z = np.bincount(lab, w * np.sin(la).ravel(), n)
    return np.degrees(np.arctan2(z, np.hypot(x, y))), np.degrees(np.arctan2(y, x))


def bodies(height, area, level, grid, features=()):
    """Seas and lakes below a common level, largest first, with depth statistics and IAU water names."""
    _, lab = bs.seas(height, area, level, min_area_fraction=0.0)
    n = int(lab.max()) + 1
    depth = np.where(lab > 0, level - height, 0.0)
    a = np.bincount(lab.ravel(), area.ravel(), n)
    vol = np.bincount(lab.ravel(), (area * depth).ravel(), n)
    maxd = np.zeros(n)
    np.maximum.at(maxd, lab.ravel(), depth.ravel())
    clat, clon = _centroids(lab, area, grid, n)
    names = {}
    for f in features:
        if f['code'] in nm.WATER_CODES:
            i, j = cell_index(grid, f['lat'], f['lon'])
            if lab[i, j]:
                names.setdefault(int(lab[i, j]), []).append(f)
    out = []
    for c in np.argsort(-a):
        if c == 0 or a[c] < BODY_MIN_SHARE * bs.AREA:
            continue
        named = sorted(names.get(int(c), []), key=lambda f: -f['diameter_km'])
        out.append(dict(label=int(c), share=round(float(a[c] / bs.AREA), 5), area_km2=round(float(a[c] / 1e6)),
                        mean_depth_m=round(float(vol[c] / a[c])), max_depth_m=round(float(maxd[c])),
                        volume_km3=round(float(vol[c] / 1e9)), centroid=[round(float(clat[c]), 2), round(float(clon[c]), 2)],
                        water_names=[f['name'] for f in named]))
    wet_share = float(a[1:].sum() / bs.AREA)
    listed = sum(b['share'] for b in out)
    return out, lab, dict(water_share=round(wet_share, 4), small_lakes_share=round(wet_share - listed, 4))


def islands(height, area, level, grid, features=()):
    """Land regions above a common level: the main land mass and the islands, with their summits."""
    _, lab = bs.seas(-height, area, -level, min_area_fraction=0.0)
    n = int(lab.max()) + 1
    a = np.bincount(lab.ravel(), area.ravel(), n)
    a[0] = 0.0
    main = int(np.argmax(a))
    named = [f for f in features if f['code'] in ('MO', 'AA') and f['diameter_km'] >= 8.0]
    flat, flon = (np.array([f[k] for f in named]) for k in ('lat', 'lon'))
    out = []
    for c in np.argsort(-a):
        if c in (0, main) or a[c] < ISLAND_MIN_KM2 * 1e6:
            continue
        top = np.argmax(np.where(lab == c, height, -np.inf))
        i, j = np.unravel_index(top, height.shape)
        slat, slon = float(grid.lat_deg[i]), float(grid.lon_deg[j])
        slon = slon - 360.0 if slon > 180.0 else slon
        near = []
        if named:
            d = great_circle_km(slat, slon, flat, flon)
            fi, fj = cell_index(grid, flat, flon)
            on = (lab[fi, fj] == c) | (d < ISLAND_NAME_RADIUS_KM)
            near = [named[k]['name'] for k in sorted(np.flatnonzero(on), key=lambda k: -named[k]['diameter_km'])][:3]
        out.append(dict(label=int(c), area_km2=round(float(a[c] / 1e6)), summit_m=round(float(height[i, j] - level)),
                        summit=[round(slat, 2), round(slon, 2)], features=near))
    land_share = float(a.sum() / bs.AREA)
    return out, lab, dict(main_land_share=round(float(a[main] / bs.AREA), 4),
                          island_share=round(land_share - float(a[main] / bs.AREA), 4))


def mare_flooding(height, level, grid, features):
    """Share of each large mare's inner cap below the level, and its median height relative to the level."""
    la, lo = np.meshgrid(grid.lat_deg, grid.lon_deg, indexing='ij')
    w = np.cos(np.radians(la))
    rows = []
    for f in sorted((f for f in features if f['code'] in ('ME', 'OC') and f['diameter_km'] > 200.0),
                    key=lambda f: -f['diameter_km']):
        cap = great_circle_km(f['lat'], f['lon'], la, lo) <= MARE_CAP * f['diameter_km'] / 2
        rows.append(dict(name=f['name'], flooded=round(float((w * cap * (height < level)).sum() / (w * cap).sum()), 3),
                         median_height_m=round(float(np.median(height[cap]) - level))))
    return rows


def compare_shares(height, area, grid, features, shares=COMPARE_SHARES):
    by = nm.by_name(features)
    la, lo = np.meshgrid(np.radians(grid.lat_deg), np.radians(grid.lon_deg), indexing='ij')
    mu = np.cos(la) * np.cos(lo)
    near = mu > 0
    rows = []
    for s in shares:
        level = bs.level_for_area(height, area, s)
        wet = height < level
        sizes, lab = bs.seas(height, area, level, min_area_fraction=1e-3)
        _, vol = bs.level_curve(height, area, [level])
        maria = [by[m] for m in NEAR_SIDE_MARIA if m in by]
        home = lab[cell_index(grid, by['Mare Imbrium']['lat'], by['Mare Imbrium']['lon'])]
        joined = [m['name'] for m in maria if lab[cell_index(grid, m['lat'], m['lon'])] == home and home]
        flood = [r['flooded'] for r in mare_flooding(height, level, grid, maria)]
        rows.append(dict(share=s, level_m=round(level), global_layer_m=round(float(vol[0] / bs.AREA)),
                         water_mass_kg=float(vol[0] * bs.WATER_DENSITY),
                         near_side_water=round(float((area * wet)[near].sum() / area[near].sum()), 3),
                         earth_facing_disk_water=round(float((area * wet * mu)[near].sum() / (area * mu)[near].sum()), 3),
                         bodies_over_0_1pct=len(sizes), largest_body_share=round(float(sizes[0]), 4),
                         near_side_maria_joined=len(joined), near_side_maria_listed=len(maria),
                         near_side_mare_interiors_flooded=round(float(np.mean(flood)), 3)))
    return rows


def point_heights(lat, lon, ppd=16, degree=200):
    """Heights above the geoid at points: LOLA at ppd, geoid from the 4-pixel/degree synthesis (smooth there).
    Also returns the height above the 1737.4 km sphere, for comparison with published landing heights."""
    radius, grid = tp.read_ldem(ppd)
    _, g4 = tp.read_ldem(4)
    geoid = tp.geoid(g4, degree)
    i, j = cell_index(grid, lat, lon)
    i4, j4 = cell_index(g4, lat, lon)
    sphere = radius[i, j]
    return sphere - geoid[i4, j4], sphere


def polar_caps(min_abs_lat=79.5, ppd=16, degree=200):
    """Heights above the geoid poleward of min_abs_lat at ppd, for polar maps."""
    radius, grid = tp.read_ldem(ppd)
    _, g4 = tp.read_ldem(4)
    geoid = tp.geoid(g4, degree)
    step = ppd // 4
    caps = {}
    for name, rows in (('south', np.flatnonzero(grid.lat_deg <= -min_abs_lat)),
                       ('north', np.flatnonzero(grid.lat_deg >= min_abs_lat))):
        caps[name] = (radius[rows] - geoid[rows // step][:, np.arange(radius.shape[1]) // step]).astype(np.float32)
        caps[name + '_lat'] = grid.lat_deg[rows]
    return caps


def compute(share=None, ppd=4, degree=200, gazetteer=None):
    share = scenario_share() if share is None else float(share)
    height, _, grid = tp.height_above_geoid(ppd, degree)
    area = grid.cell_area_m2
    features = nm.features(gazetteer)
    level = bs.level_for_area(height, area, share)
    water, wlab, wstats = bodies(height, area, level, grid, features)
    land, llab, lstats = islands(height, area, level, grid, features)
    _, vol = bs.level_curve(height, area, [level])
    result = dict(share=share, sea_level_m=round(level, 1), global_layer_m=round(float(vol[0] / bs.AREA), 1),
                  water_mass_kg=float(vol[0] * bs.WATER_DENSITY), **wstats, **lstats,
                  bodies=water, islands=land, mare_flooding=mare_flooding(height, level, grid, features),
                  share_comparison=compare_shares(height, area, grid, features),
                  grid=dict(pixels_per_degree=ppd, geoid_degree=degree))
    arrays = dict(height_m=height.astype(np.float32), water_label=wlab.astype(np.int32), land_label=llab.astype(np.int32),
                  lat_deg=grid.lat_deg, lon_deg=grid.lon_deg)
    return result, arrays


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--share', type=float, help='standing-water share (default: shared/scenarios/water.json)')
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--products', type=Path, default=HERE / 'products')
    args = parser.parse_args(argv)
    if fetch_inputs.restore()['status'] != 'PASS':
        raise SystemExit('BLOCKED: geography inputs missing; run python -m geography.fetch_inputs --download')
    result, arrays = compute(args.share)
    pct = round(100 * result['share'])
    meta = dict(schema=SCHEMA, producer=producer(), evidence=EVIDENCE, reading_rule=READING_RULE,
                units=dict(heights='m above the geoid', areas='km2', shares='fraction of the whole Moon',
                           water_mass='kg at 1000 kg/m3'), **result)
    args.out.mkdir(parents=True, exist_ok=True)
    name = 'atlas.json' if args.share is None else f'atlas_{pct}pct.json'
    (args.out / name).write_text(json.dumps(meta, indent=1) + '\n')
    args.products.mkdir(parents=True, exist_ok=True)
    grid_meta = dict(schema=GRID_SCHEMA, producer=producer(), evidence=EVIDENCE, reading_rule=READING_RULE,
                     share=result['share'], sea_level_m=result['sea_level_m'], atlas=name)
    np.savez_compressed(args.products / f'atlas_{pct}pct_4ppd.npz', metadata=json.dumps(grid_meta), **arrays, **polar_caps())
    print(f"{result['share']:.0%} water: sea level {result['sea_level_m']:+.0f} m, layer {result['global_layer_m']:.0f} m; "
          f"main land {result['main_land_share']:.1%}, islands {result['island_share']:.1%}")
    for b in result['bodies']:
        print(f"  {b['share']:6.2%} {b['area_km2']:>9,} km2  mean {b['mean_depth_m']:5d} m  max {b['max_depth_m']:5d} m  "
              f"{', '.join(b['water_names'][:4])}")


if __name__ == '__main__':
    main()
