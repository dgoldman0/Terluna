"""The globe's close views at 64 pixels per degree (about 470 m at the equator): packs of tiles of the design surface's
finer ground, in today's frame, which the page loads only for the ground in view when the camera comes close.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 visualization/atlas/fine_tiles.py \\
        --surface /media/projectspace/terluna-research/atlas-round/globe_surface/globe_surface_16ppd.npz \\
        --surface-record research/studies/atlas_round/results/globe_surface.json \\
        --moisture research/studies/atlas_round/results/moisture_map_final.json \\
        --falls research/studies/atlas_round/results/falls.json     # about a quarter of an hour on three workers

The finer ground is research/studies/atlas_round/fine_surface.py's product: the design heights at 64 px/deg, each
1.9 km cell's mean kept at the design surface's, and every sea and lake at its own level on them. Each tile, 8 degrees a
side (512 x 512 pixels), from 80 N to 80 S, holds three layers:

- colour: on land, the colour the page's 1.9 km true colour gives the land before the water is laid over it
  (appearance.true_colour; rivers and white water included, at 1.9 km), read at each fine pixel's place in the settled
  frame; on the water, the designed waters' colours as the page sets them, by the depth the fine heights give under the
  water's level (the shelves' productive-coast colour in the shallows) and the 1.9 km plumes;
- shading: the relief's tangent-space normals from the fine heights, at the page's exaggeration (as
  appearance.relief_normals makes them);
- water and sky: the fine water (red), the share of an even sky each place sees past the fine relief (green, as
  appearance.sky_view finds it) and a flag where the tile holds data (blue).

A pack is 4 x 5 tiles, 32 by 40 degrees, as three images: the colour and the shading as WebP, the water and sky as
lossless WebP. Tiles of open water deeper than the shelves are left empty, and packs holding no tile are left out. The
page reads the tiles in today's frame and turns its lookup, and the shading's north and east, by the settled frame's
rotation. Poleward of 80 degrees the page keeps its 1.9 km ground. The ground's colour stays the 1.9 km cover's: what
is finer is the relief, its shading and sky, the shores and the water's depth.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import appearance  # noqa: E402
from research.studies.atlas_round import fine_surface as FS  # noqa: E402

OUT = HERE / 'out' / 'fine'
SCHEMA = 'terluna.visualization.globe-fine-tiles/1'
PPD, TILE, TILE_DEG = FS.PPD, 512, 8
LAT_TOP = 80                                  # the tiles span 80 N to 80 S
TILE_ROWS, TILE_COLS = 2 * LAT_TOP // TILE_DEG, 360 // TILE_DEG
PACK_ROWS, PACK_COLS = 4, 5                   # tiles in a pack
MARGIN_KM = appearance.SKY_REACH_KM + 2.0     # the sky search's reach, with the gradients' few pixels
QUALITY = dict(colour=85, normal=80)
WORKERS = 3

_G = {}                                       # the arrays the workers read, set before the pool forks


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


# -------------------------------------------------------------------------------------------- per tile

def normals(h, lat, ppd):
    """appearance.relief_normals for a window of rows whose latitudes are given."""
    dy = appearance.MOON_KM * 1000.0 * np.radians(1.0 / ppd)
    gy, gx = appearance.gradients(h, dy, lat)
    n = np.stack([-gx * appearance.NORMAL_STRENGTH, gy * appearance.NORMAL_STRENGTH, np.ones_like(h)], axis=-1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return np.round((n + 1) / 2 * 255).astype(np.uint8)


def sky(h, lat, ppd, water):
    """appearance.sky_view for a window of rows whose latitudes are given, held at the window's edges (its margins
    reach past the search)."""
    import cv2
    rows, cols = h.shape
    km = appearance.MOON_KM * np.radians(1.0 / ppd)
    per_col = (1.0 / (km * np.maximum(np.cos(np.radians(lat)), 0.02))).astype(np.float32)[:, None]
    yy, xx = np.mgrid[0:rows, 0:cols].astype(np.float32)
    seen = np.zeros_like(h)
    for k in range(appearance.SKY_AZIMUTHS):
        north, east = np.cos(2.0 * np.pi * k / appearance.SKY_AZIMUTHS), np.sin(2.0 * np.pi * k / appearance.SKY_AZIMUTHS)
        steepest = np.zeros_like(h)
        for d in np.geomspace(0.5 * km, appearance.SKY_REACH_KM, appearance.SKY_STEPS):
            ahead = cv2.remap(h, xx + np.float32(d * east) * per_col, yy - np.float32(d * north / km),
                              cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
            np.maximum(steepest, ((ahead - h) / 1000.0 - d * d / (2.0 * appearance.MOON_KM)) / d, out=steepest)
        seen += 1.0 / (1.0 + steepest * steepest)
    seen /= appearance.SKY_AZIMUTHS
    seen[water] = 1.0
    return seen


def sample16(field, lat_t, lon_t):
    """A 16 px/deg field in the settled frame (180 W at the left edge) read bilinear at settled latitudes and
    longitudes, wrapping in longitude."""
    import cv2
    x = ((lon_t + 180.0) * 16.0 - 0.5).astype(np.float32)
    y = ((90.0 - lat_t) * 16.0 - 0.5).astype(np.float32)
    return cv2.remap(field, x, y, cv2.INTER_LINEAR, borderMode=cv2.BORDER_WRAP)


def tile(tr, tc):
    """One tile's three layers (uint8 arrays), or None for open water deeper than the shallows."""
    g = _G
    r0, c0 = (90 - LAT_TOP + TILE_DEG * tr) * PPD, TILE_DEG * tc * PPD
    code = g['water'][r0:r0 + TILE, c0:c0 + TILE]
    if code.all():
        depth = (g['level'][r0:r0 + TILE, c0:c0 + TILE].astype(np.float32)
                 - g['height'][r0:r0 + TILE, c0:c0 + TILE].astype(np.float32)) * 0.5
        if depth.min() > appearance.SHELF_M[0]:
            return None
    lat_edge = LAT_TOP - TILE_DEG * tr
    worst = max(abs(lat_edge), abs(lat_edge - TILE_DEG)) + 2.0
    km = appearance.MOON_KM * np.radians(1.0 / PPD)
    mr = int(np.ceil(MARGIN_KM / km))
    mc = int(min(np.ceil(MARGIN_KM / (km * np.cos(np.radians(worst)))), 6 * TILE))
    rows = np.clip(np.arange(r0 - mr, r0 + TILE + mr), 0, FS.ROWS - 1)
    cols = np.arange(c0 - mc, c0 + TILE + mc) % FS.COLS
    h = g['height'][rows][:, cols].astype(np.float32) * 0.5
    lev = g['level'][rows][:, cols].astype(np.float32) * 0.5
    wet = g['water'][rows][:, cols] > 0
    lat = 90.0 - (rows + 0.5) / PPD
    surf = np.where(wet, lev, h)
    inner = (slice(mr, mr + TILE), slice(mc, mc + TILE))
    normal = normals(h, lat, PPD)[inner]
    seen = sky(surf, lat, PPD, wet)[inner]
    h, lev, wet = h[inner], lev[inner], wet[inner]
    aux = np.stack([wet.astype(np.uint8) * 255, np.round(np.clip(seen, 0, 1) * 255).astype(np.uint8),
                    np.full(wet.shape, 255, np.uint8)], axis=-1)
    return dict(colour=colour_layer(r0, c0, h, lev, wet), normal=normal, aux=aux)


def colour_layer(r0, c0, h, lev, wet):
    """A tile's colour (uint8, sRGB) from its fine heights, water levels and water: on land the page's 1.9 km land
    colour read at each fine pixel's place in the settled frame; on the water the designed waters' colours by the
    depth under the water's level (appearance.SHELF_M) and the 1.9 km plumes."""
    from climate.gcm.turned import to_turned
    g = _G
    la = np.repeat((90.0 - (np.arange(r0, r0 + TILE) + 0.5) / PPD)[:, None], TILE, axis=1)
    lo = np.repeat(((np.arange(c0, c0 + TILE) + 0.5) / PPD)[None, :], TILE, axis=0)
    lat_t, lon_t = to_turned(la, lo, g['rotation'])
    land = np.stack([sample16(g['land'][..., k], lat_t, lon_t) for k in range(3)], axis=-1)
    plume = np.clip(sample16(g['plume'], lat_t, lon_t), 0.0, 1.0)[..., None]
    shelf = appearance.smoothstep(appearance.SHELF_M[0], appearance.SHELF_M[1], np.where(wet, lev - h, 0.0))[..., None]
    w = g['waters']
    sea = (w['open_sea'] * (1 - shelf) + w['productive_coast'] * shelf) * (1 - plume) + w['river_mouth'] * plume
    colour = np.where(wet[..., None], sea, land)
    return np.round(appearance.srgb_encode(colour) * 255).astype(np.uint8)


def water_tile(tr, tc):
    """One tile's colour and water alone, as tile() makes them, without the shading and the sky search."""
    g = _G
    r0, c0 = (90 - LAT_TOP + TILE_DEG * tr) * PPD, TILE_DEG * tc * PPD
    part = (slice(r0, r0 + TILE), slice(c0, c0 + TILE))
    h = g['height'][part].astype(np.float32) * 0.5
    lev = g['level'][part].astype(np.float32) * 0.5
    wet = g['water'][part] > 0
    return dict(colour=colour_layer(r0, c0, h, lev, wet), wet=wet)


def _tile(args):
    return args, tile(*args)


def _water_tile(args):
    return args, water_tile(*args)


# -------------------------------------------------------------------------------------------- build

def colour_parts(surface_path, surface_record, moisture_path, falls_path, log=print):
    """The page's 1.9 km true colour of the land (before the water is laid over it) and the plumes' share, in the
    settled frame with 180 W at the left edge, computed as render.write_globe_surface and appearance.build compute it."""
    import render
    grid = render.surface_grid(surface_path)
    record = json.loads(Path(surface_record).read_text())
    meta = grid['meta']
    level0 = float(meta['common_level_m'])
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    turned = appearance.turn_texture(grid['height'].astype(np.float32), rotation, left_lon=0.0, nearest=True)
    h = np.roll(turned, turned.shape[1] // 2, axis=1)
    routing = render.routing_product(record)
    drainage = appearance.surface_grids(h.shape, surface_path, rotation, routing=routing)
    moisture = appearance.moisture_record(moisture_path)
    rock = appearance.rock_fields(h.shape, 16, rotation)
    waters = appearance.water_albedos()
    falls = None
    if falls_path:
        falls_record = json.loads(Path(falls_path).read_text())
        falls = appearance.falls_grids(h.shape, falls_record, Path(falls_record['product']['path']), rotation)
    parts = {}
    appearance.true_colour(h, level0, 16, drainage, moisture, rock, waters, falls, parts=parts)
    log('the 1.9 km land colour and plumes computed')
    return parts, waters, rotation, routing


def build(surface, surface_record, moisture, falls, out=OUT, log=print, workers=WORKERS):
    from PIL import Image
    started = time.time()
    fine_record = json.loads(FS.OUT.read_text())
    if sha256(FS.PRODUCT) != fine_record['product']['sha256']:
        raise SystemExit(f'{FS.PRODUCT} no longer matches its record {FS.OUT}')
    parts, waters, rotation, routing = colour_parts(surface, surface_record, moisture, falls, log)
    with np.load(FS.PRODUCT, allow_pickle=False) as z:
        _G.update(height=z['height'], water=z['water'], level=z['level'])
    _G.update(land=parts['land'], plume=parts['plume'].astype(np.float32), waters=waters, rotation=rotation)
    log(f'fine surface read, {time.time() - started:.0f} s')
    jobs = [(tr, tc) for tr in range(TILE_ROWS) for tc in range(TILE_COLS)]
    tiles = {}
    with Pool(workers) as pool:
        for k, (key, t) in enumerate(pool.imap_unordered(_tile, jobs, chunksize=4)):
            tiles[key] = t
            if (k + 1) % 100 == 0:
                log(f'{k + 1} of {len(jobs)} tiles, {time.time() - started:.0f} s')
    out.mkdir(parents=True, exist_ok=True)
    packs, total = [], 0
    for pr in range(TILE_ROWS // PACK_ROWS):
        for pc in range(TILE_COLS // PACK_COLS):
            keys = [(pr * PACK_ROWS + i, pc * PACK_COLS + j) for i in range(PACK_ROWS) for j in range(PACK_COLS)]
            if not any(tiles[k] is not None for k in keys):
                continue
            layers = {name: np.zeros((PACK_ROWS * TILE, PACK_COLS * TILE, 3), np.uint8) for name in ('colour', 'normal', 'aux')}
            valid = []
            for (tr, tc) in keys:
                t = tiles[(tr, tc)]
                valid.append(t is not None)
                if t is None:
                    continue
                i, j = (tr - pr * PACK_ROWS) * TILE, (tc - pc * PACK_COLS) * TILE
                for name in layers:
                    layers[name][i:i + TILE, j:j + TILE] = t[name]
            files = {}
            for name, a in layers.items():
                path = out / f'p{pr}_{pc}_{name}.webp'
                img = Image.fromarray(a, 'RGB')
                if name == 'aux':
                    img.save(path, 'WEBP', lossless=True, quality=60, method=4)
                else:
                    img.save(path, 'WEBP', quality=QUALITY[name], method=4)
                files[name] = path.name
                total += path.stat().st_size
            packs.append(dict(r=pr, c=pc, lat_top=LAT_TOP - pr * PACK_ROWS * TILE_DEG, lon_left=pc * PACK_COLS * TILE_DEG,
                              files=files, valid=valid))
    manifest = dict(
        schema=SCHEMA, ppd=PPD, tile_px=TILE, tile_deg=TILE_DEG, lat_top=LAT_TOP, pack_tiles=[PACK_ROWS, PACK_COLS],
        frame='today', rotation_today_to_turned=np.asarray(rotation).tolist(), packs=packs, bytes=total,
        tiles=sum(1 for t in tiles.values() if t is not None),
        sources={str(FS.OUT.relative_to(ROOT)): sha256(FS.OUT), str(FS.PRODUCT): fine_record['product']['sha256'],
                 str(surface): sha256(surface), str(moisture): sha256(moisture), str(routing): sha256(routing)},
        credit=('Finer heights: SLDEM2015 (Barker et al. 2016; NASA\'s LOLA Team and JAXA\'s SELENE/Kaguya Team) within '
                '60 degrees and LOLA\'s LDEM_64 (LOLA Science Team, NASA GSFC; PDS Geosciences Node) beyond, carried onto '
                'the design surface.'),
        seconds=round(time.time() - started))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    log(f"{manifest['tiles']} tiles in {len(packs)} packs, {total / 1e6:.1f} MB, {manifest['seconds']} s")
    return manifest


def refresh_water(surface, surface_record, moisture, falls, out=OUT, base=None, log=print, workers=WORKERS):
    """The packs' colour and water made again by the present rules on the present fine surface, each tile keeping its
    shading and its sky: for a change to the waters' colours, the plumes or the fine water, none of which the slow
    sky search depends on beyond the few cells where a shore has moved (water sees the whole sky; ground the water
    has left keeps the water's). In `out` the tiles of fine_look.py's trial windows are left as they stand; in `base`
    (the folder holding the packs as they stood before the trial, for the packs kept there) every tile is made again,
    so that the trial can be cut again on them."""
    from PIL import Image
    started = time.time()
    out = Path(out)
    manifest = json.loads((out / 'manifest.json').read_text())
    fine_record = json.loads(FS.OUT.read_text())
    if sha256(FS.PRODUCT) != fine_record['product']['sha256']:
        raise SystemExit(f'{FS.PRODUCT} no longer matches its record {FS.OUT}')
    parts, waters, rotation, routing = colour_parts(surface, surface_record, moisture, falls, log)
    with np.load(FS.PRODUCT, allow_pickle=False) as z:
        _G.update(height=z['height'], water=z['water'], level=z['level'])
    _G.update(land=parts['land'], plume=parts['plume'].astype(np.float32), waters=waters, rotation=rotation)
    log(f'fine surface read, {time.time() - started:.0f} s')
    trial = {tuple(t) for w in manifest.get('look', {}).get('windows', []) for t in w['tiles']}
    folders = [(out, trial)] + ([(Path(base), set())] if base else [])
    done = 0
    with Pool(workers) as pool:
        for n, pack in enumerate(manifest['packs']):
            pr, pc = pack['r'], pack['c']
            keys = [(pr * PACK_ROWS + i, pc * PACK_COLS + j) for i in range(PACK_ROWS) for j in range(PACK_COLS)]
            tiles = dict(pool.imap_unordered(_water_tile, [k for k, ok in zip(keys, pack['valid']) if ok]))
            for folder, kept in folders:
                paths = {name: folder / pack['files'][name] for name in ('colour', 'aux')}
                if not all(p.exists() for p in paths.values()):
                    continue
                colour = np.array(Image.open(paths['colour']).convert('RGB'))
                aux = np.array(Image.open(paths['aux']).convert('RGB'))
                for (tr, tc), t in tiles.items():
                    if (tr, tc) in kept:
                        continue
                    here = (slice((tr - pr * PACK_ROWS) * TILE, (tr - pr * PACK_ROWS + 1) * TILE),
                            slice((tc - pc * PACK_COLS) * TILE, (tc - pc * PACK_COLS + 1) * TILE))
                    colour[here] = t['colour']
                    aux[here + (0,)] = t['wet'] * np.uint8(255)
                    aux[here + (1,)] = np.where(t['wet'], 255, aux[here + (1,)])
                    done += folder is out
                Image.fromarray(colour).save(paths['colour'], 'WEBP', quality=QUALITY['colour'], method=4)
                Image.fromarray(aux).save(paths['aux'], 'WEBP', lossless=True, quality=60, method=4)
            log(f"{n + 1} of {len(manifest['packs'])} packs, {time.time() - started:.0f} s")
    manifest['bytes'] = sum((out / f).stat().st_size for p in manifest['packs'] for f in p['files'].values())
    manifest['water'] = dict(
        note=('The colour and the water of every tile outside the trial windows were made again after the first cut, '
              'by the rules and on the fine surface named here; the shading and the sky are the first cut\'s.'),
        tiles=done, shelf_m=list(appearance.SHELF_M),
        sources={str(FS.OUT.relative_to(ROOT)): sha256(FS.OUT), str(FS.PRODUCT): fine_record['product']['sha256'],
                 str(surface): sha256(surface), str(moisture): sha256(moisture), str(routing): sha256(routing),
                 'visualization/atlas/appearance.py': sha256(HERE / 'appearance.py')},
        seconds=round(time.time() - started))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    log(f"{done} tiles' colour and water made again in {len(manifest['packs'])} packs, {manifest['water']['seconds']} s")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--surface', type=Path, required=True)
    ap.add_argument('--surface-record', type=Path, required=True)
    ap.add_argument('--moisture', type=Path, required=True)
    ap.add_argument('--falls', type=Path)
    ap.add_argument('--workers', type=int, default=WORKERS)
    ap.add_argument('--water-only', action='store_true',
                    help='make the colour and the water of the packs in place again, keeping the shading and the sky')
    ap.add_argument('--base', type=Path, help='with --water-only: the folder of the packs as they stood before the trial look')
    a = ap.parse_args(argv)
    args = (a.surface.absolute(), a.surface_record.absolute(), a.moisture.absolute(), a.falls.absolute() if a.falls else None)
    if a.water_only:
        refresh_water(*args, base=a.base, workers=a.workers)
    else:
        build(*args, workers=a.workers)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
