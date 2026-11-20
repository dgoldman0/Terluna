"""The close views' colour and water cut again from the ground's look (ground_look.py), for chosen packs or all of them.

    ulimit -v 12000000; OPENBLAS_NUM_THREADS=1 python3 visualization/atlas/look_tiles.py \\
        --surface <globe_surface_16ppd.npz> --moisture research/studies/atlas_round/results/moisture_map_final.json \\
        --base <folder of the packs to start from> --out <folder> [--sand <sand_16ppd.npz>] [--shading <folder>] \\
        [--packs 1_2,3_0] [--workers 3]

fine_tiles.py cuts each tile's three layers; ground_look.py makes the land's colour again from the ground at 474 m.
This module joins them. Every tile in the packs of `base` keeps its sky; its colour becomes the look's, and the red of
its third layer the look's share of each pixel under water (rivers included). A tile keeps its shading too, except
where the sand product holds a made sea: there the dunes' slopes are added to the ground's own shading, which is read
from `--shading` (packs no cut has changed) so that a cut made twice gives the same packs. The ground is worked in
blocks of 2 x 2 tiles, each with ground_look's margin round it, so that water gathers across the blocks' edges and
neighbouring blocks agree. A tile fine_tiles.py left empty (open water deeper than the shelves) stays empty.

In a folder beside `out`, named <out>_look, it writes look_average_16ppd.npz: the look averaged to the far view's
1.9 km on today's grid (linear colour, the share under water and the cells covered), from which a far texture can be
made that agrees with the tiles (far_from_look.py); what the look keeps between runs lies there too, so that `out`
holds only what a page serves. `out` must be a folder other than the page's own sources; the packs are copied into it
and changed there.
"""
from __future__ import annotations

import argparse
import json
import shutil
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
import ground_look as GL  # noqa: E402

PPD, TILE = FT.PPD, FT.TILE
BLOCK_ROWS, BLOCK_COLS = 2, 2       # the tiles worked as one piece of ground: 16 degrees a side, about 2.5 GB a worker
FAR = 4                             # fine pixels to a far cell's side (64 to 16 a degree)
NOTE = ("the land's colour at 474 m is made from the ground itself. Each pixel mixes open water, ground kept wet, "
        'rain-fed plants and open ground of several kinds (steppe and desert soils, stony ground, pale pans, dark clay on '
        'basalt, silt and sand), placed by display rules on the fine heights, the water, the rock maps and the climate. '
        'Rivers of 30 m³/s or more lie in their valleys at their computed widths. It is a guess at how the design\'s '
        'dryness, rock and water look as land, and says nothing of which plants grow. A made sand sea is drawn with dunes '
        'of the form, spacing and height the sand study gives it; where each dune stands is made up.')

_G = {}                             # what the workers read, set before the pool forks


def blocks(pack):
    """A pack's valid tiles as blocks of at most BLOCK_ROWS by BLOCK_COLS tiles: each (first tile row, first tile
    column, tile rows, tile columns). A block with no valid tile is left out."""
    pr, pc = pack['r'], pack['c']
    valid = np.asarray(pack['valid'], bool).reshape(FT.PACK_ROWS, FT.PACK_COLS)
    out = []
    for i in range(0, FT.PACK_ROWS, BLOCK_ROWS):
        for j in range(0, FT.PACK_COLS, BLOCK_COLS):
            part = valid[i:i + BLOCK_ROWS, j:j + BLOCK_COLS]
            if part.any():
                out.append((pr * FT.PACK_ROWS + i, pc * FT.PACK_COLS + j, part.shape[0], part.shape[1]))
    return out


def block_look(block):
    """A block's look: its colour (uint8 sRGB), the share of each pixel under water, and the slopes of a made sea's
    dunes where the block holds any (eastward and northward, float16), else None."""
    tr, tc, nr, nc = block
    origin = ((90 - FT.LAT_TOP + FT.TILE_DEG * tr) * PPD, FT.TILE_DEG * tc * PPD)
    out = GL.look_window(tr, tc, _G, size=(nr * TILE, nc * TILE), origin=origin)
    dunes = None if out['dunes'] is None else tuple(a.astype(np.float16) for a in out['dunes'])
    return block, out['colour'], np.round(np.clip(out['water'], 0.0, 1.0) * 255).astype(np.uint8), dunes


def with_dunes(normal, east, north):
    """A shading layer (the relief's tangent-space normals, uint8, as appearance.relief_normals bakes them) with the
    slopes of a made sea's dunes added to the ground's own."""
    n = normal.astype(np.float32) / 255.0 * 2.0 - 1.0
    up = np.maximum(n[..., 2], 0.05)
    tilt = np.stack([n[..., 0] / up - east.astype(np.float32) * appearance.NORMAL_STRENGTH,
                     n[..., 1] / up - north.astype(np.float32) * appearance.NORMAL_STRENGTH, np.ones_like(up)], axis=-1)
    tilt /= np.linalg.norm(tilt, axis=-1, keepdims=True)
    return np.round((tilt + 1.0) / 2.0 * 255.0).astype(np.uint8)


def far_mean(a, n=FAR):
    """The mean of each n x n block of pixels along the first two axes."""
    r, c = a.shape[0] // n, a.shape[1] // n
    return a[:r * n, :c * n].reshape(r, n, c, n, *a.shape[2:]).mean(axis=(1, 3))


def cut(surface, moisture, base, out, sand=None, packs=None, workers=3, log=print, shading=None):
    from PIL import Image
    started = time.time()
    base, out = Path(base), Path(out)
    if out.resolve() == (HERE / 'out' / 'fine').resolve():
        raise SystemExit('cut the look into a folder other than the page\'s own sources')
    manifest = json.loads((base / 'manifest.json').read_text())
    chosen = [p for p in manifest['packs'] if packs is None or (p['r'], p['c']) in packs]
    out.mkdir(parents=True, exist_ok=True)
    if out.resolve() != base.resolve():
        for p in manifest['packs']:
            for name in p['files'].values():
                if not (out / name).exists():
                    shutil.copy2(base / name, out / name)
    side = out.parent / (out.name + '_look')
    side.mkdir(exist_ok=True)
    g = GL.read_inputs(Path(surface), Path(moisture), cache=side, log=log, sand=Path(sand) if sand is not None else None)
    _G.update(g)
    rows16, cols16 = 180 * 16, 360 * 16
    far = dict(colour=np.zeros((rows16, cols16, 3), np.float32), water=np.zeros((rows16, cols16), np.float32),
               covered=np.zeros((rows16, cols16), bool))
    kept = side / 'look_average_16ppd.npz'
    if kept.exists() and packs is not None:              # a cut of chosen packs adds to the average kept from before
        with np.load(kept, allow_pickle=False) as z:
            far = {k: z[k].astype(far[k].dtype) for k in far}
    shading = Path(shading) if shading is not None else base     # where the ground's own shading is read: packs no cut has changed
    shaded = []
    jobs = [(p, b) for p in chosen for b in blocks(p)]
    log(f'{len(chosen)} packs in {len(jobs)} blocks, {workers} workers')
    done = {(p['r'], p['c']): [] for p in chosen}
    need = {(p['r'], p['c']): len(blocks(p)) for p in chosen}
    by_key = {(p['r'], p['c']): p for p in chosen}
    tiles = 0
    with Pool(workers) as pool:
        for n, (block, colour, water, dunes) in enumerate(pool.imap_unordered(block_look, [b for _, b in jobs])):
            tr, tc, nr, nc = block
            key = (tr // FT.PACK_ROWS, tc // FT.PACK_COLS)
            done[key].append((block, colour, water, dunes))
            r16, c16 = (90 - FT.LAT_TOP + FT.TILE_DEG * tr) * 16, FT.TILE_DEG * tc * 16
            part = (slice(r16, r16 + nr * TILE // FAR), slice(c16, c16 + nc * TILE // FAR))
            far['colour'][part] = far_mean(decode(colour))
            far['water'][part] = far_mean(water.astype(np.float32) / 255.0)
            far['covered'][part] = True
            if len(done[key]) == need[key]:              # the pack is whole: write its colour and water
                pack = by_key[key]
                paths = {name: out / pack['files'][name] for name in ('colour', 'aux')}
                img = {name: np.array(Image.open(paths[name]).convert('RGB')) for name in paths}
                valid = np.asarray(pack['valid'], bool).reshape(FT.PACK_ROWS, FT.PACK_COLS)
                ground = shading / pack['files']['normal']               # the ground's own shading, kept apart from the packs cut here
                relief = None
                if any(d is not None for _, _, _, d in done[key]):
                    relief = np.array(Image.open(ground).convert('RGB'))
                for (btr, btc, bnr, bnc), col, wat, dun in done[key]:
                    for i in range(bnr):
                        for j in range(bnc):
                            ti, tj = btr + i - key[0] * FT.PACK_ROWS, btc + j - key[1] * FT.PACK_COLS
                            if not valid[ti, tj]:
                                continue
                            here = (slice(ti * TILE, (ti + 1) * TILE), slice(tj * TILE, (tj + 1) * TILE))
                            src = (slice(i * TILE, (i + 1) * TILE), slice(j * TILE, (j + 1) * TILE))
                            img['colour'][here] = col[src]
                            img['aux'][here + (0,)] = wat[src]
                            if dun is not None:
                                relief[here] = with_dunes(relief[here], dun[0][src], dun[1][src])
                            tiles += 1
                Image.fromarray(img['colour']).save(paths['colour'], 'WEBP', quality=FT.QUALITY['colour'], method=4)
                Image.fromarray(img['aux']).save(paths['aux'], 'WEBP', lossless=True, quality=60, method=4)
                if relief is not None:                           # dunes stand in this pack: its shading is the ground's with theirs added
                    Image.fromarray(relief).save(out / pack['files']['normal'], 'WEBP', quality=FT.QUALITY['normal'], method=4)
                    shaded.append(pack['files']['normal'])
                elif ground.resolve() != (out / pack['files']['normal']).resolve() and FT.sha256(ground) != FT.sha256(out / pack['files']['normal']):
                    shutil.copy2(ground, out / pack['files']['normal'])   # no dunes here now: the ground's own shading back
                done[key] = []
            log(f'{n + 1} of {len(jobs)} blocks, {time.time() - started:.0f} s')
    np.savez_compressed(kept, colour=far['colour'].astype(np.float16), water=far['water'].astype(np.float16), covered=far['covered'])
    manifest['bytes'] = sum((out / f).stat().st_size for p in manifest['packs'] for f in p['files'].values())
    whole = packs is None
    manifest['look'] = dict(note=NOTE, windows=[], whole=whole, packs=[[p['r'], p['c']] for p in chosen], tiles=tiles, dunes_shaded=sorted(shaded),
                            settings={k: v for k, v in GL.LOOK.items()},
                            sources=dict(g['sources'], **{f'visualization/atlas/{name}': FT.sha256(HERE / name)
                                                          for name in ('ground_look.py', 'ground_kinds.py')}),
                            seconds=round(time.time() - started))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1, default=list) + '\n')
    log(f"{tiles} tiles' colour and water cut from the look in {len(chosen)} packs, {manifest['look']['seconds']} s")
    return manifest


def decode(colour):
    """sRGB bytes as linear values."""
    v = np.asarray(colour, np.float32) / 255.0
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4).astype(np.float32)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--surface', type=Path, required=True)
    ap.add_argument('--moisture', type=Path, required=True)
    ap.add_argument('--base', type=Path, required=True, help='the folder of packs to start from (manifest.json beside them)')
    ap.add_argument('--out', type=Path, required=True, help='where the changed packs go; never the page\'s own sources')
    ap.add_argument('--sand', type=Path, help="the sand study's grid product (sand_16ppd.npz), or the made seas' (sand_seas_16ppd.npz)")
    ap.add_argument('--shading', type=Path, help="a folder of packs whose shading no cut has changed (the page's own sources will do; only "
                                                 "read): where a made sea's dunes stand, the pack's shading is made from it")
    ap.add_argument('--packs', help='packs to cut, as row_column pairs joined by commas (all of them without it)')
    ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args(argv)
    packs = {tuple(int(v) for v in p.split('_')) for p in a.packs.split(',')} if a.packs else None
    cut(a.surface.absolute(), a.moisture.absolute(), a.base.absolute(), a.out.absolute(),
        sand=a.sand.absolute() if a.sand else None, packs=packs, workers=a.workers, log=lambda s: print(s, flush=True),
        shading=a.shading.absolute() if a.shading else None)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
