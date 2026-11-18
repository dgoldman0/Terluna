"""The design surface at 64 pixels per degree (about 470 m at the equator), for the globe's close views: heights that
keep every 1.9 km cell's mean at the design surface's height, and every sea and lake at its own level on them.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.fine_surface   # a few minutes

The heights' shape comes from the finest seamless grids tied to LOLA: within 60 degrees of the equator SLDEM2015
(LOLA's altimetry merged with stereo heights from Kaguya's Terrain Camera, 128 px/deg averaged to 64), which fills the
gaps between LOLA's tracks that LOLA's own grids interpolate across, and poleward LOLA's LDEM_64, the two joined across
57-60 degrees. The design surface's 16 px/deg heights (globe_surface.py: LOLA less the geoid plus the settled water's
load) differ from that grid's 4 x 4 block means by the geoid, the load and the small differences between LOLA's
separately gridded products; that difference is carried to the fine grid as a smooth correction whose block means
equal it (bilinear enlargement, its residuals enlarged again twelve times, then the last centimetres put in exactly), and
added. Every 1.9 km cell's mean is therefore the height the routing, the lakes' levels and spills, the water's share
and the page's checks were computed on, and inside each cell the ground keeps the fine grid's own shape.

The water: a fine cell under one of the design surface's water cells is water where its height lies below that cell's
water level, so islets and points the fine heights raise above the water show as land; and a fine cell under a land
cell beside the water is water where it lies below the level of the nearest water cell and joins water across fine
cells, so coves and inlets narrower than 1.9 km show as water. Each sea and lake keeps the level the design gives it:
a notch narrower than 1.9 km that would lower a lake's spill is not looked for, and a crater pond smaller than a 1.9 km
cell that the routing never filled stays dry.

The product, on the research drive, holds height (int16 half-metres above the design surface's datum), water (uint8:
0 dry ground, 1-7 the design surface's seas by their body number, 8 a lake) and level (int16 half-metres, the water's
level where there is water), rows from 90 N and columns from 0 E, 64 per degree, on today's Moon.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import globe_surface as GS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'fine_surface.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/fine_surface')
PRODUCT = BULK / 'fine_surface_64ppd.npz'
INPUTS = Path('/media/projectspace/terluna-research/atlas-round/fine_ground/inputs')
LDEM64 = INPUTS / 'ldem_64.img'
SLDEM = INPUTS / 'sldem2015' / 'sldem2015_128_60s_60n_000_360_int16.img'
SLDEM_JP2 = INPUTS / 'sldem2015' / 'sldem2015_128_60s_60n_000_360.jp2'
SCHEMA = 'terluna.research.atlas-round-fine-surface/1'
PRODUCT_SCHEMA = 'terluna.research.atlas-round-fine-surface-grid/1'
CODE = ('research/studies/atlas_round/fine_surface.py',)

PPD, K = 64, 4                       # the fine grid, and fine cells along a 16 px/deg cell
ROWS, COLS = 180 * PPD, 360 * PPD
JOIN_DEG = (57.0, 60.0)              # SLDEM2015 gives way to LDEM_64 across these latitudes
LIFT_ROUNDS = 12                     # the smooth correction's rounds; each cuts the block means' residual by half or more
MOON_KM = 1737.4
LAKE = 8                             # the water code of a lake off the seas


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def rel(path) -> str:
    p = Path(path)
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)


# ------------------------------------------------------------------------------------------- grids

def block_mean(a, k=K):
    """Mean over k x k blocks."""
    r, c = a.shape[0] // k, a.shape[1] // k
    return a[:r * k, :c * k].reshape(r, k, c, k).mean(axis=(1, 3), dtype=np.float64).astype(np.float32)


def expand(a, k=K):
    """Each cell repeated over its k x k block."""
    return np.repeat(np.repeat(a, k, axis=0), k, axis=1)


def enlarge(a, k=K):
    """A coarse field enlarged k times, bilinear between cell centres, wrapping in longitude (columns) and held at the
    first and last rows."""
    import cv2
    pad = np.concatenate([a[:, -1:], a, a[:, :1]], axis=1).astype(np.float32)
    big = cv2.resize(pad, (pad.shape[1] * k, a.shape[0] * k), interpolation=cv2.INTER_LINEAR)
    return big[:, k:k + a.shape[1] * k]


def lift_to_means(coarse, k=K, rounds=LIFT_ROUNDS):
    """A smooth fine field whose k x k block means equal the coarse field: its bilinear enlargement, with the block
    means' residuals enlarged and added back `rounds` times, and what remains (a small fraction of a metre) added to
    each block whole. Returns the field and the largest residual before that last step."""
    fine = enlarge(coarse, k)
    for _ in range(rounds):
        fine += enlarge(coarse - block_mean(fine, k), k)
    rest = coarse - block_mean(fine, k)
    fine += expand(rest, k)
    return fine, float(np.abs(rest).max())


# ------------------------------------------------------------------------------------------- heights

def ldem64():
    """LOLA's LDEM_64: metres above the 1,737.4 km sphere, rows from 90 N, columns from 0 E."""
    return np.fromfile(LDEM64, dtype='<i2').reshape(ROWS, COLS).astype(np.float32) * 0.5


def sldem64():
    """SLDEM2015 at 128 px/deg averaged to 64: metres above the 1,737.4 km sphere for the rows from 60 N to 60 S,
    columns turned to start at 0 E (the product starts at 180 W). Cells the product leaves empty are NaN."""
    raw = np.memmap(SLDEM, dtype='<i2', mode='r', shape=(120 * 128, 360 * 128))
    out = np.empty((120 * PPD, COLS), np.float32)
    for r0 in range(0, raw.shape[0], 1024):
        band = np.asarray(raw[r0:r0 + 1024]).astype(np.float32)
        band[band <= -32767] = np.nan
        out[r0 // 2:(r0 + band.shape[0]) // 2] = band.reshape(band.shape[0] // 2, 2, -1, 2).mean(axis=(1, 3)) * 0.5
    return np.roll(out, -180 * PPD, axis=1)


def merged_heights(log=print):
    """The detail's source on the fine grid: SLDEM2015 within 57 degrees, LDEM_64 beyond 60, a linear blend between;
    LDEM_64 wherever SLDEM2015 is empty."""
    fine = ldem64()
    s = sldem64()
    lat = 90.0 - (np.arange(ROWS) + 0.5) / PPD
    r0 = 30 * PPD                                       # the fine grid's row at 60 N
    w = np.clip((JOIN_DEG[1] - np.abs(lat[r0:r0 + s.shape[0]])) / (JOIN_DEG[1] - JOIN_DEG[0]), 0.0, 1.0)[:, None]
    part = fine[r0:r0 + s.shape[0]]
    gap = ~np.isfinite(s)
    s[gap] = part[gap]
    part += w * (s - part)
    log(f'merged heights: SLDEM2015 empty in {gap.mean():.2e} of its cells, filled from LDEM_64')
    return fine


def fine_heights(design16, merged):
    """The design heights at 64 px/deg: the merged grid with a smooth correction added whose 4 x 4 block means are
    the design surface's heights less the merged grid's block means (the geoid, the water's load and the small
    differences between LOLA's gridded products). Each block's mean becomes the design height, and inside a block the
    ground keeps the merged grid's own shape with only the correction's gentle slope added. Returns the heights and
    the lift's last residual (m)."""
    correction, residual = lift_to_means(design16 - block_mean(merged))
    merged += correction
    return merged, residual


# ------------------------------------------------------------------------------------------- water

def fine_water(height, water16, level16, body16):
    """The water on the fine heights: under a water cell, the fine cells below its level; under a land cell beside the
    water, the fine cells below the nearest water cell's level that join water across fine cells. Returns the water
    code (uint8) and the level (float32, NaN on dry ground)."""
    from scipy import ndimage
    code16 = np.where(water16, np.where(body16 == 0, LAKE, body16), 0).astype(np.uint8)
    # each land cell beside the water takes its nearest water cell's level and code
    _, (ii, jj) = ndimage.distance_transform_edt(~water16, return_indices=True)
    near_level, near_code = level16[ii, jj], code16[ii, jj]
    beside = ~water16 & ndimage.binary_dilation(water16, structure=np.ones((3, 3), bool))
    beside[:, 0] |= ~water16[:, 0] & water16[:, -1]                 # across the date line
    beside[:, -1] |= ~water16[:, -1] & water16[:, 0]
    lev16 = np.where(water16, level16, np.where(beside, near_level, np.nan)).astype(np.float32)
    cod16 = np.where(water16, code16, np.where(beside, near_code, 0)).astype(np.uint8)
    level = expand(lev16)
    under = expand(water16)
    wet = height < level                                            # NaN levels compare false
    candidate = wet & ~under
    core = wet & under
    labels, n = ndimage.label(wet, structure=np.ones((3, 3), bool))
    keep = np.bincount(labels[core].ravel(), minlength=n + 1) > 0
    keep[0] = False
    joined = keep[labels]
    del labels
    water = (core | (candidate & joined))
    code = np.where(water, expand(cod16), 0).astype(np.uint8)
    level = np.where(water, level, np.nan).astype(np.float32)
    return code, level


# ------------------------------------------------------------------------------------------- build

def area_weights(rows, ppd):
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    return np.cos(np.radians(lat))


def build(log=print):
    started = time.time()
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        design16 = z['height_m'].astype(np.float32)
        cells, levels, bodies = z['water_cell'], z['water_level_m'], z['water_body']
    rows16, cols16 = design16.shape
    water16 = np.zeros(rows16 * cols16, bool)
    water16[cells] = True
    level16 = np.full(rows16 * cols16, np.nan, np.float32)
    level16[cells] = levels
    body16 = np.zeros(rows16 * cols16, np.int16)
    body16[cells] = bodies
    water16, level16, body16 = (a.reshape(rows16, cols16) for a in (water16, level16, body16))
    log(f'design surface read, {time.time() - started:.0f} s')
    merged = merged_heights(log)
    log(f'merged heights, {time.time() - started:.0f} s')
    height, residual = fine_heights(design16, merged)
    del merged
    off = block_mean(np.round(height * 2.0) * 0.5) - design16      # as stored, in half-metres
    log(f'fine heights, {time.time() - started:.0f} s; block means off the design by at most {np.abs(off).max():.3f} m')
    code, level = fine_water(height, water16, level16, body16)
    log(f'fine water, {time.time() - started:.0f} s')

    # checks: the water's share at each grid, and where the fine water departs from the 1.9 km cells
    w16 = area_weights(rows16, 16)[:, None] * np.ones((1, cols16), np.float32)
    share16 = float((water16 * w16).sum() / w16.sum())
    wf = area_weights(ROWS, PPD)
    fine_wet = code > 0
    share = float((fine_wet.sum(axis=1) * wf).sum() / (wf.sum() * COLS))
    under = expand(water16)
    raised = float(((under & ~fine_wet).sum(axis=1) * wf).sum() / (wf.sum() * COLS))
    flooded = float(((~under & fine_wet).sum(axis=1) * wf).sum() / (wf.sum() * COLS))
    lakes_share = float(((code == LAKE).sum(axis=1) * wf).sum() / (wf.sum() * COLS))
    del under, fine_wet
    log(f'water: {share16:.4f} of the Moon at 1.9 km, {share:.4f} at the fine grid ({raised:.4f} raised to land, '
        f'{flooded:.4f} flooded beside the water)')

    BULK.mkdir(parents=True, exist_ok=True)
    pmeta = dict(schema=PRODUCT_SCHEMA, grid=dict(pixels_per_degree=PPD, rows=ROWS, cols=COLS, first_row='90 N',
                                                  first_col='0 E', frame='today'),
                 bodies=meta['bodies'], codes={'0': 'dry ground', **{str(b): meta['bodies'][b] for b in range(1, 8)},
                                               str(LAKE): 'a lake off the seas'},
                 reading_rule=('Rows from 90 N and columns from 0 E, 64 per degree, on today\'s Moon. height and level are '
                               'int16 half-metres above the design surface\'s datum (its height_m); level is 0 on dry '
                               'ground. water is 0 on dry ground, 1-7 a sea by the design surface\'s body number, 8 a '
                               'lake off the seas.'))
    half = lambda a: np.clip(np.round(np.nan_to_num(a, nan=0.0) * 2.0), -32768, 32767).astype(np.int16)
    np.savez_compressed(PRODUCT, metadata=json.dumps(pmeta), height=half(height), water=code, level=half(level))
    log(f'product written, {time.time() - started:.0f} s')
    inputs = {rel(GS.OUT): sha256(GS.OUT)}
    ignored = {str(GS.PRODUCT): sha256(GS.PRODUCT), str(LDEM64): sha256(LDEM64), str(SLDEM_JP2): sha256(SLDEM_JP2)}
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}), inputs=inputs, ignored_inputs=ignored,
        evidence=('A calculation on products: the design surface (globe_surface.py) carried to 64 px/deg with its 1.9 km '
                  'cells\' means kept, the detail from SLDEM2015 (Barker et al. 2016; LOLA with Kaguya Terrain Camera '
                  'stereo; public domain, credit NASA\'s LOLA Team and JAXA\'s SELENE/Kaguya Team) within 60 degrees and '
                  'LOLA\'s LDEM_64 beyond, and the design\'s seas and lakes at their levels on those heights. The detail '
                  'is today\'s ground shape; the young Moon keeps its landforms as LOLA shows them, barring erosion.'),
        reading_rule=('product holds the fine grid (see its metadata). checks compares the fine grid with the design '
                      'surface: the block means\' largest departure (m, from the int16 rounding), and the water\'s share '
                      'of the Moon at both grids, with the share the fine heights raise above the water inside the 1.9 km '
                      'water cells and the share they flood beside it.'),
        parameters=dict(ppd=PPD, join_deg=list(JOIN_DEG), lift_rounds=LIFT_ROUNDS),
        checks=dict(block_mean_largest_off_m=round(float(np.abs(off).max()), 3), lift_last_residual_m=round(residual, 4),
                    water_share_16ppd=round(share16, 5), water_share_fine=round(share, 5),
                    raised_to_land_share=round(raised, 5), flooded_beside_share=round(flooded, 5),
                    lake_share_fine=round(lakes_share, 5)),
        product=dict(path=str(PRODUCT), sha256=sha256(PRODUCT)), seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')
    print(f"wrote {OUT}: {json.dumps(record['checks'])}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
