"""A trial of the land's close-up look, in a few windows: the ground's colour and texture at 474 m as a guesstimate,
with the rivers drawn as water at their computed width.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python3 visualization/atlas/fine_look.py \\
        --surface /media/projectspace/terluna-research/atlas-round/globe_surface/globe_surface_16ppd.npz \\
        --surface-record research/studies/atlas_round/results/globe_surface.json \\
        --moisture research/studies/atlas_round/results/moisture_map_final.json \\
        --window Kohulim:14,2 --window "Karadala to Hamada Aswad:12,319" ...     # a minute or two for five windows

fine_tiles.py's tiles carry the page's 1.9 km colour onto the finer ground. Here, in each window (2 x 2 tiles, 16
degrees a side, round a point on today's Moon), the colour is made again on the 474 m ground by the page's own cover
rules (appearance.ground_colour and soil_colour: plants over weathered soil by the moisture map's dryness), with the
terms that follow the ground evaluated on the finer heights, and a grain for each kind of land:

- wetter, so greener, in hollows (the smoothed ground's curvature at about 1.4 and 5.7 km), along the lines water drains
  down (the cells draining through each place on the fine heights), beside rivers (a belt about a kilometre wide) and on
  fresh shores under 20 m above the water; drier, and rockier, where the ground steepens;
- rivers as water at the width the page's rule gives their flow (appearance.RIVER_WIDTH), on courses through the lowest
  ground of each routed 1.9 km cell and its receiver, for rivers of LOOK['river_m3s'] or more on the design routing;
- closed canopy mottled at about 1.6 km, grassland patchy at about 8 km, and on flat, bare, dry ground fields of dune
  crests a few kilometres apart and several times as long, with the drainage lines paling to dry washes where little
  grows;
- water in the open water's colour, taking the productive coast's green only in its shallows (LOOK['littoral_m'],
  from 20 m, fully at 3 m). The trial set this rule for its windows, where the page's former one (from 150 m) painted
  a whole shallow lake or bay the green of the forest beside it; the page has since taken it for all its water
  (appearance.SHELF_M), and the windows read it from there.

All of it is a guesstimate, labelled so on the page: it makes the moisture map's dryness look like land, and says
nothing of which plants grow. The dunes' crests run along the parallels as an illustration: the design climate's wind
directions are not yet a product. The windows' tiles replace fine_tiles.py's in their packs; the rest stay.
"""
from __future__ import annotations

import argparse
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
import fine_tiles as FT  # noqa: E402
from research.studies.atlas_round import design_routing as DRt  # noqa: E402
from research.studies.atlas_round import fine_surface as FS  # noqa: E402

LOOK = dict(
    hollow_ln=(0.22, 0.15),     # wetter in hollows: times the standardized curvature at 1.4 and 5.7 km, clipped at 2.5
    slope_ln=0.5,               # drier where the ground steepens from 0.15 to 0.5
    thread_ln=0.6,              # wetter along drainage lines, from 60 to 2,000 fine cells draining through
    river_ln=1.6,               # the green belt beside a river, about a kilometre wide
    shore_ln=1.0,               # and on fresh shores under 20 m above the water, within about 2 km of it
    patch_ln=0.22,              # patchiness, in log aridity
    canopy=0.20, grass=0.07,    # brightness mottling of closed canopy and of grassland
    dune=0.09, dune_km=(3.5, 20.0),   # dune crests on flat, bare, dry ground: brightness; spacing and length (km)
    littoral_m=appearance.SHELF_M,    # water shallower than these takes the productive coast's colour: the page's rule
    wash=0.35,                  # how far a drainage line pales toward the bare soil where little grows
    river_m3s=30.0)             # rivers drawn from this flow up
SALT = 3                        # the water code of the salt sea, whose shores keep no fresh belt
PPD, TILE, K = FT.PPD, FT.TILE, FS.K
_G = FT._G


def sphere_points(lat, lon):
    """Unit vectors of latitudes and longitudes (deg), for noise that has no seam anywhere on the sphere."""
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def fbm(p, km, seed, octaves=3):
    """Value noise summed over octaves from a wavelength of `km`, between 0 and 1 with a mean near 0.5."""
    total, norm, gain = 0.0, 0.0, 1.0
    for k in range(octaves):
        total = total + gain * appearance.value_noise(p * (appearance.MOON_KM / km * 2 ** k), seed + 17 * k)
        norm += gain
        gain *= 0.5
    return (total / norm).astype(np.float32)


def crests(p, across_km, along_km, seed):
    """Crest lines (1 on a crest, falling to -1 between) of noise stretched along the parallels: `across_km` from crest
    to crest going north, about `along_km` long. The stretch is made by scaling the sphere's polar axis, so the lines
    have no seam; toward the poles they shorten."""
    scale = np.array([appearance.MOON_KM / along_km, appearance.MOON_KM / along_km, appearance.MOON_KM / (2.0 * across_km)],
                     np.float32)
    n = 0.7 * appearance.value_noise(p * scale, seed) + 0.3 * appearance.value_noise(p * scale * 2.3, seed + 5)
    return (1.0 - 4.0 * np.abs(n - 0.5)).clip(-1.0, 1.0).astype(np.float32)


def drainage(surf, wet):
    """The number of fine cells draining through each cell, down the steepest of its eight neighbours; water takes
    what reaches it and passes nothing on."""
    rows, cols = surf.shape
    pad = np.pad(surf, 1, mode='edge')
    best = np.zeros(surf.shape, np.float32)
    rcv = np.full(surf.shape, -1, np.int64)
    rr, cc = np.arange(rows)[:, None], np.arange(cols)[None, :]
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            drop = (surf - pad[1 + dr:1 + dr + rows, 1 + dc:1 + dc + cols]) / (1.414 if dr and dc else 1.0)
            better = drop > best
            target = np.clip(rr + dr, 0, rows - 1) * cols + np.clip(cc + dc, 0, cols - 1)
            rcv = np.where(better, target, rcv)
            best = np.where(better, drop, best)
    rcv[wet] = -1
    order = np.argsort(surf, axis=None)[::-1].tolist()
    r, acc = rcv.ravel().tolist(), [1.0] * (rows * cols)
    for i in order:
        j = r[i]
        if j >= 0:
            acc[j] += acc[i]
    return np.asarray(acc, np.float32).reshape(rows, cols)


def river_shares(h, rows, cols, lat):
    """The share of each fine pixel a river's water covers, and a strength (0 to 1) for its green belt, over a padded
    window of fine heights whose global rows and columns are given: each routed river cell joined to its receiver
    through the lowest fine cell of each one's 4 x 4 block."""
    import cv2
    g = _G
    cols16 = FS.COLS // K
    share = np.zeros(h.shape, np.float32)
    belt = np.zeros(h.shape, np.float32)
    r_first, c_first = int(rows[0]), int(cols[0])
    if rows[-1] - rows[0] != len(rows) - 1:              # a window clamped at the grid's edge: no rivers drawn
        return share, belt
    i0, i1 = (r_first + K - 1) // K, (int(rows[-1]) + 1) // K
    span = len(cols)
    j_first = (c_first + K - 1) // K
    levels = {}

    def low_point(i16, j16):
        lr = i16 * K - r_first
        lc = (j16 * K - c_first) % FS.COLS
        if lr < 0 or lr + K > h.shape[0] or lc + K > span:
            return None
        block = h[lr:lr + K, lc:lc + K]
        y, x = divmod(int(np.argmin(block)), K)
        return lc + x, lr + y

    for i16 in range(i0, i1):
        row_q = g['q16'][i16 * cols16:(i16 + 1) * cols16]
        for dj in range((span - K) // K):
            j16 = (j_first + dj) % cols16
            q = float(row_q[j16])
            if q < LOOK['river_m3s']:
                continue
            cell = i16 * cols16 + j16
            if g['water16'][cell]:
                continue
            rcv = int(g['rcv16'][cell])
            if rcv < 0 or rcv == cell:
                continue
            a, b = low_point(i16, j16), low_point(rcv // cols16, rcv % cols16)
            if a is None or b is None or abs(a[0] - b[0]) > 3 * K:
                continue
            width = appearance.RIVER_WIDTH * np.sqrt(q) / (appearance.MOON_KM * 1000.0 * np.radians(1.0 / PPD))
            key = (max(1, int(round(width))), round(min(width, 1.0), 1) or 0.1, round(min(1.0, np.sqrt(q / 300.0)), 1))
            levels.setdefault(key, []).append((a, b))
    for (thick, value, strength), segments in levels.items():
        canvas = np.zeros(h.shape, np.float32)
        for a, b in segments:
            cv2.line(canvas, a, b, 1.0, thick, cv2.LINE_AA)
        np.maximum(share, canvas * value, out=share)
        np.maximum(belt, canvas * strength, out=belt)
    belt = np.clip(cv2.GaussianBlur(belt, (0, 0), 2.0) * 2.5, 0.0, 1.0)
    return share, belt


def look_tile(tr, tc):
    """One tile's three layers with the trial look."""
    import cv2
    from climate.gcm.turned import to_turned
    g = _G
    r0, c0 = (90 - FT.LAT_TOP + FT.TILE_DEG * tr) * PPD, FT.TILE_DEG * tc * PPD
    lat_edge = FT.LAT_TOP - FT.TILE_DEG * tr
    worst = max(abs(lat_edge), abs(lat_edge - FT.TILE_DEG)) + 2.0
    km = appearance.MOON_KM * np.radians(1.0 / PPD)
    mr = int(np.ceil(FT.MARGIN_KM / km))
    mc = int(min(np.ceil(FT.MARGIN_KM / (km * np.cos(np.radians(worst)))), 6 * TILE))
    rows = np.clip(np.arange(r0 - mr, r0 + TILE + mr), 0, FS.ROWS - 1)
    cols = np.arange(c0 - mc, c0 + TILE + mc) % FS.COLS
    h = g['height'][rows][:, cols].astype(np.float32) * 0.5
    lev = g['level'][rows][:, cols].astype(np.float32) * 0.5
    code = g['water'][rows][:, cols]
    wet = code > 0
    lat = 90.0 - (rows + 0.5) / PPD
    surf = np.where(wet, lev, h)
    inner = (slice(mr, mr + TILE), slice(mc, mc + TILE))
    normal = FT.normals(h, lat, PPD)[inner]
    seen = FT.sky(surf, lat, PPD, wet)[inner]

    # what follows the ground, on the padded window
    gy, gx = appearance.gradients(surf, km * 1000.0, lat)
    slope = np.hypot(gx, gy)
    ln = np.zeros(h.shape, np.float32)
    for sigma, weight in zip((3.0, 12.0), LOOK['hollow_ln']):
        curve = cv2.Laplacian(cv2.GaussianBlur(surf, (0, 0), sigma), cv2.CV_32F)
        ln += weight * np.clip(curve / (curve[~wet].std() + 1e-6), -2.5, 2.5)
    acc = drainage(surf, wet)
    thread = appearance.smoothstep(np.log10(60.0), np.log10(2000.0), np.log10(acc))
    share, belt = river_shares(h, rows, cols, lat)
    fresh = wet & (code != SALT)
    near_level = cv2.dilate(np.where(fresh, lev, -1e9).astype(np.float32),
                            cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9)))
    away = cv2.distanceTransform((~fresh).astype(np.uint8), cv2.DIST_L2, 3)        # fine cells from fresh water
    shore = np.where((near_level > -1e8) & ~wet,
                     appearance.smoothstep(20.0, 3.0, h - near_level) * appearance.smoothstep(4.5, 1.0, away), 0.0)
    ln += LOOK['thread_ln'] * thread + LOOK['river_ln'] * belt + LOOK['shore_ln'] * shore
    ln -= LOOK['slope_ln'] * appearance.smoothstep(0.15, 0.5, slope)
    ln, slope, thread, share = ln[inner], slope[inner], thread[inner], share[inner]
    h, lev, wet = h[inner], lev[inner], wet[inner]

    # the climate's dryness and the rock, read in the settled frame; the grain, on the sphere
    la = np.repeat((90.0 - (np.arange(r0, r0 + TILE) + 0.5) / PPD)[:, None], TILE, axis=1)
    lo = np.repeat(((np.arange(c0, c0 + TILE) + 0.5) / PPD)[None, :], TILE, axis=0)
    lat_t, lon_t = to_turned(la, lo, g['rotation'])
    p = sphere_points(la, lo)
    ln += FT.sample16(g['ln16'], lat_t, lon_t) + LOOK['patch_ln'] * (fbm(p, 6.0, 11) - 0.5) / 0.146
    aridity = np.exp(ln)
    soil = appearance.soil_colour(aridity, FT.sample16(g['feo'], lat_t, lon_t), FT.sample16(g['tio2'], lat_t, lon_t))
    rocky = 0.6 * appearance.smoothstep(0.35, 0.7, slope)
    colour = appearance.ground_colour(aridity, soil, FT.sample16(g['mafic'], lat_t, lon_t), rocky)
    cover = appearance.plant_cover(aridity)
    trees = appearance.share_of(aridity, appearance.TREES)
    colour *= (1.0 + LOOK['canopy'] * 2.0 * (fbm(p, 1.6, 23) - 0.5) * trees * cover)[..., None]
    colour *= (1.0 + LOOK['grass'] * 2.0 * (fbm(p, 8.0, 31) - 0.5) * (1.0 - trees) * cover)[..., None]
    bare = (1.0 - cover) * appearance.smoothstep(0.12, 0.03, slope)
    fields = appearance.smoothstep(0.45, 0.6, fbm(p, 60.0, 41))
    colour *= (1.0 + LOOK['dune'] * (crests(p, *LOOK['dune_km'], 53) - 0.5) * bare * fields)[..., None]
    pale = (LOOK['wash'] * thread * (1.0 - cover))[..., None]
    colour = colour * (1.0 - pale) + np.clip(soil * 1.25, 0.0, 1.0) * pale
    w = g['waters']
    colour = colour * (1.0 - share[..., None]) + w['river_mouth'] * share[..., None]
    # the plumes' 1.9 km share broken into streaks, as mud mixing into clear water is
    plume = np.clip(FT.sample16(g['plume'], lat_t, lon_t) * (0.4 + 1.2 * fbm(p, 4.0, 71)), 0.0, 0.8)[..., None]
    depth = np.where(wet, lev - h, 0.0)
    shelf = appearance.smoothstep(*LOOK['littoral_m'], depth)[..., None]
    sea = (w['open_sea'] * (1 - shelf) + w['productive_coast'] * shelf) * (1 - plume) + w['river_mouth'] * plume
    colour = np.where(wet[..., None], sea, colour)
    colour = np.round(appearance.srgb_encode(np.clip(colour, 0.0, 1.0)) * 255).astype(np.uint8)
    water = np.maximum(wet.astype(np.float32), share)
    aux = np.stack([np.round(water * 255).astype(np.uint8), np.round(np.clip(seen, 0, 1) * 255).astype(np.uint8),
                    np.full(wet.shape, 255, np.uint8)], axis=-1)
    return dict(colour=colour, normal=normal, aux=aux)


def _look(args):
    return args, look_tile(*args)


def window_tiles(lat, lon):
    """The 2 x 2 tiles round a point on today's Moon."""
    tr = int(np.clip(round((FT.LAT_TOP - lat) / FT.TILE_DEG) - 1, 0, FT.TILE_ROWS - 2))
    tc = int(round((lon % 360.0) / FT.TILE_DEG) - 1)
    return [(tr + i, (tc + j) % FT.TILE_COLS) for i in (0, 1) for j in (0, 1)]


def build(surface, surface_record, moisture, falls, windows, out=FT.OUT, log=print, workers=FT.WORKERS):
    import render
    from PIL import Image
    started = time.time()
    manifest = json.loads((out / 'manifest.json').read_text())
    grid = render.surface_grid(surface)
    meta = grid['meta']
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    # the 1.9 km fields the look reads, kept on the research drive while their inputs and appearance.py stand, so that
    # a change to the look's settings costs the tiles alone
    key = FT.sha256(surface)[:16] + FT.sha256(moisture)[:16] + FT.sha256(HERE / 'appearance.py')[:16]
    cache = FS.BULK / 'look_fields.npz'
    fields = None
    if cache.exists():
        with np.load(cache, allow_pickle=False) as z:
            if str(z['key']) == key:
                fields = {k: z[k] for k in ('plume', 'ln16', 'feo', 'tio2', 'mafic')}
    if fields is None:
        parts, _, _, _ = FT.colour_parts(surface, surface_record, moisture, falls, log)
        shape16 = (grid['rows'], grid['cols'])
        rock = appearance.rock_fields(shape16, 16, rotation)
        fields = dict(plume=parts['plume'].astype(np.float32),
                      ln16=appearance.moisture_field(appearance.moisture_record(moisture), shape16).astype(np.float32),
                      feo=rock['feo'], tio2=rock['tio2'], mafic=rock['mafic'])
        np.savez(cache, key=key, **fields)
    waters = appearance.water_albedos()
    routing = json.loads(DRt.OUT.read_text())['product']
    if FT.sha256(routing['path']) != routing['sha256']:
        raise SystemExit(f"{routing['path']} no longer matches its record")
    with np.load(routing['path'], allow_pickle=False) as z:
        q16 = z['discharge_m3s'].astype(np.float32)
        rcv16 = z['receiver'].astype(np.int32)
        water16 = z['routing_lake'] | (z['sea'] >= 0)
        water16[z['lake_cell']] = True
    with np.load(FS.PRODUCT, allow_pickle=False) as z:
        _G.update(height=z['height'], water=z['water'], level=z['level'])
    _G.update(waters=waters, rotation=rotation, q16=q16, rcv16=rcv16, water16=water16, **fields)
    log(f'fields read, {time.time() - started:.0f} s')
    jobs = sorted({t for _, lat, lon in windows for t in window_tiles(lat, lon)})
    tiles = {}
    with Pool(workers) as pool:
        for k, (key, t) in enumerate(pool.imap_unordered(_look, jobs, chunksize=1)):
            tiles[key] = t
            log(f'{k + 1} of {len(jobs)} tiles, {time.time() - started:.0f} s')
    packs = {(p['r'], p['c']): p for p in manifest['packs']}
    touched = sorted({(tr // FT.PACK_ROWS, tc // FT.PACK_COLS) for tr, tc in tiles})
    for pr, pc in touched:
        pack = packs[(pr, pc)]
        for name in ('colour', 'normal', 'aux'):
            path = out / pack['files'][name]
            a = np.array(Image.open(path).convert('RGB'))
            for (tr, tc), t in tiles.items():
                if (tr // FT.PACK_ROWS, tc // FT.PACK_COLS) != (pr, pc):
                    continue
                i, j = (tr - pr * FT.PACK_ROWS) * TILE, (tc - pc * FT.PACK_COLS) * TILE
                a[i:i + TILE, j:j + TILE] = t[name]
            img = Image.fromarray(a)
            if name == 'aux':
                img.save(path, 'WEBP', lossless=True, quality=60, method=4)
            else:
                img.save(path, 'WEBP', quality=FT.QUALITY[name], method=4)
    manifest['bytes'] = sum((out / f).stat().st_size for p in manifest['packs'] for f in p['files'].values())
    manifest['look'] = dict(
        windows=[dict(name=name, today=[lat, lon], tiles=window_tiles(lat, lon)) for name, lat, lon in windows],
        packs=[packs[k]['files'] for k in touched], settings=LOOK,
        note=('In these trial windows the land\'s colour and texture at 474 m are a guesstimate: the page\'s cover rules '
              'on the finer ground, greener in hollows, along drainage lines, beside rivers and on fresh shores, barer on '
              'steep ground; rivers of 30 m³/s or more as water at their computed width; canopy mottling, grassland '
              'patches and, on flat bare ground, dune crests whose direction is an illustration; water green only in '
              'its shallows.'))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    log(f"{len(tiles)} tiles in {len(touched)} packs rewritten, {time.time() - started:.0f} s")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--surface', type=Path, required=True)
    ap.add_argument('--surface-record', type=Path, required=True)
    ap.add_argument('--moisture', type=Path, required=True)
    ap.add_argument('--falls', type=Path)
    ap.add_argument('--window', action='append', required=True, help='name:latitude,longitude on today\'s Moon')
    ap.add_argument('--workers', type=int, default=FT.WORKERS)
    a = ap.parse_args(argv)
    windows = []
    for w in a.window:
        name, at = w.rsplit(':', 1)
        lat, lon = (float(v) for v in at.split(','))
        windows.append((name, lat, lon))
    build(a.surface.absolute(), a.surface_record.absolute(), a.moisture.absolute(),
          a.falls.absolute() if a.falls else None, windows, workers=a.workers)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
