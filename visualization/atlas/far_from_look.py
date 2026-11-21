"""A far view that agrees with the close look, for a trial copy of the page: its 1.9 km colour and water textures taken
from the average of the close look.

    ulimit -v 8000000; python3 visualization/atlas/far_from_look.py --page <trial page folder> \\
        --average <fine folder>_look/look_average_16ppd.npz [--before <folder to keep the textures as they were>]

look_tiles.py leaves the look averaged to 1.9 km on today's grid. The page's far textures are in the settled frame,
180 W at the left edge. Here each pixel of textures/albedo.jpg and textures/water.png is read from that average at its
place on today's Moon, wherever the look was cut; within two degrees of the tiles' edge at 80 degrees it runs into the
texture as it was, and poleward of that and over open water the tiles leave empty the texture stays. With
--match-caps the land poleward of the tiles takes the colour of the look's last ring before the edge at its own
longitude (cap_ring) times the earlier texture's detail (land_detail), since the look has not been worked there.
The average holds linear colour and the texture sRGB. It changes a copy of the page only: it refuses the page's own sources. Built into
render.py this would be a part of the page's colour (appearance.true_colour); as a trial it is a step after the build.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import appearance  # noqa: E402
import fine_tiles as FT  # noqa: E402
from climate.gcm.turned import to_today  # noqa: E402

EDGE_DEG = 2.0                      # the look runs into the texture as it was over this many degrees at the tiles' edge
RING = (76.0, 80.0)                 # the look's last ring before the tiles' edge: its land's colour, longitude by longitude,
RING_SMOOTH_DEG = 12.0              # smoothed round the circle over this, is the colour the caps take
DETAIL = (0.75, 1.3)                # and the caps keep the earlier texture's detail: its land's brightness over its own mean
DETAIL_PX = 20.0                    # within this many pixels, held between these bounds
FAR_NOTE = (' The far view takes the same look averaged to 1.9 km; poleward of 80° it keeps the earlier ground\'s '
            'detail in the colours of the look\'s last ring.')
QUALITY = 88                        # render.py's for the colour texture


def decode(colour):
    """sRGB bytes as linear values."""
    v = np.asarray(colour, np.float32) / 255.0
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4).astype(np.float32)


def _chunks(albedo, average, rotation, chunk):
    """The texture in bands of rows: each band's rows, its pixels' latitudes and longitudes on today's Moon, and the
    look's colour, water share and weight read at them."""
    import cv2
    rows, cols = albedo.shape[:2]
    colour = np.ascontiguousarray(average['colour'], np.float32)
    share = np.ascontiguousarray(average['water'], np.float32)
    covered = np.ascontiguousarray(average['covered'], np.float32)
    tex, avg = cols / 360.0, colour.shape[1] / 360.0       # pixels a degree: the texture's, the average's
    lon_t = -180.0 + (np.arange(cols) + 0.5) / tex
    wrap = dict(interpolation=cv2.INTER_LINEAR, borderMode=cv2.BORDER_WRAP)
    for r0 in range(0, rows, chunk):
        lat_t = 90.0 - (np.arange(r0, min(r0 + chunk, rows)) + 0.5) / tex
        la, lo = to_today(np.repeat(lat_t[:, None], cols, axis=1), np.repeat(lon_t[None, :], lat_t.size, axis=0), rotation)
        x = (np.mod(lo, 360.0) * avg - 0.5).astype(np.float32)
        y = np.clip((90.0 - la) * avg - 0.5, 0, colour.shape[0] - 1).astype(np.float32)
        # every one of the four cells a pixel is read from must hold the look, or an empty cell's black bleeds in
        k = (cv2.remap(covered, x, y, **wrap) > 0.999) * appearance.smoothstep(FT.LAT_TOP, FT.LAT_TOP - EDGE_DEG, np.abs(la))
        yield slice(r0, r0 + lat_t.size), la, lo, cv2.remap(colour, x, y, **wrap), cv2.remap(share, x, y, **wrap), k.astype(np.float32)


def cap_ring(average):
    """For the north (1) and the south (-1): the look's mean land colour (linear) in RING at each column of the
    average's grid, smoothed round the circle; None where the ring holds no land the look was cut on."""
    from scipy.ndimage import gaussian_filter1d
    colour, share, covered = (np.asarray(average[k], np.float32) for k in ('colour', 'water', 'covered'))
    ppd = colour.shape[0] / 180.0
    out = {}
    for s in (1, -1):
        r0, r1 = sorted((int(round((90.0 - s * RING[0]) * ppd)), int(round((90.0 - s * RING[1]) * ppd))))
        land = (covered[r0:r1] > 0.5) & (share[r0:r1] < 0.1)
        weight = gaussian_filter1d(land.sum(axis=0).astype(np.float64), RING_SMOOTH_DEG * ppd, mode='wrap')
        summed = gaussian_filter1d((colour[r0:r1] * land[..., None]).sum(axis=0).astype(np.float64), RING_SMOOTH_DEG * ppd,
                                   axis=0, mode='wrap')
        out[s] = (summed / np.maximum(weight, 1e-9)[:, None]).astype(np.float32) if weight.max() > 1e-6 else None
    return out


def land_detail(albedo, water):
    """The earlier texture's detail on land: each pixel's brightness over the mean brightness of the land within
    DETAIL_PX of it, held within DETAIL."""
    import cv2
    lum = decode(albedo) @ np.array([0.2126, 0.7152, 0.0722], np.float32)
    land = (np.asarray(water) < 128).astype(np.float32)
    mean = cv2.GaussianBlur(lum * land, (0, 0), DETAIL_PX) / np.maximum(cv2.GaussianBlur(land, (0, 0), DETAIL_PX), 1e-3)
    return np.clip(lum / np.maximum(mean, 1e-4), DETAIL[0], DETAIL[1]).astype(np.float32)


def compose(albedo, water, average, rotation, chunk=240, caps=False):
    """The far colour (uint8 sRGB) and water (uint8) with the look's average laid in, and the share of the texture's
    pixels that took it. `average` holds colour (linear), water (share) and covered on today's grid, rows from 90 N,
    columns from 0 E; the textures are in the settled frame, rows from 90 N, columns from 180 W. With `caps` the land
    the look does not reach poleward of the tiles takes the colour of the look's last ring at its longitude (cap_ring)
    times the earlier texture's own detail (land_detail)."""
    rows, cols = albedo.shape[:2]
    out_c, out_w = albedo.copy(), water.copy()
    taken = 0.0
    ring = cap_ring(average) if caps else None
    detail = land_detail(albedo, water) if caps else None
    per_degree = np.asarray(average['colour']).shape[1] / 360.0
    for part, la, lo, c, w, k in _chunks(albedo, average, rotation, chunk):
        old = decode(albedo[part])
        if caps:
            column = (np.mod(lo, 360.0) * per_degree).astype(np.int64) % np.asarray(average['colour']).shape[1]
            for s in (1, -1):
                if ring[s] is None:
                    continue
                cap = (s * la > FT.LAT_TOP - EDGE_DEG) & (water[part] < 128)
                old = np.where(cap[..., None], ring[s][column] * detail[part][..., None], old)
        mixed = old * (1.0 - k[..., None]) + c * k[..., None]
        out_c[part] = np.round(appearance.srgb_encode(np.clip(mixed, 0.0, 1.0)) * 255).astype(np.uint8)
        out_w[part] = np.round(water[part] * (1.0 - k) + 255.0 * np.clip(w, 0.0, 1.0) * k).astype(np.uint8)
        taken += float(k.sum())
    return out_c, out_w, taken / (rows * cols)


def main(argv=None) -> int:
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--page', type=Path, required=True, help='a copy of the built page: the folder holding textures/ and fine/')
    ap.add_argument('--average', type=Path, required=True, help="look_tiles.py's look_average_16ppd.npz")
    ap.add_argument('--before', type=Path, help='a folder that keeps the two textures as they were (made once, never overwritten)')
    ap.add_argument('--match-caps', action='store_true',
                    help="give the land poleward of the tiles the colour of the look's last ring at its longitude, keeping its own detail")
    a = ap.parse_args(argv)
    page = a.page.absolute()
    if page.resolve() == (HERE / 'out').resolve():
        raise SystemExit('this changes a copy of the page, never the page\'s own sources')
    textures = page / 'textures'
    before = a.before.absolute() if a.before else page.parent / (page.name + '_far_before')
    before.mkdir(parents=True, exist_ok=True)
    for name in ('albedo.jpg', 'water.png'):
        if not (before / name).exists():
            shutil.copy2(textures / name, before / name)
    albedo = np.array(Image.open(before / 'albedo.jpg').convert('RGB'))
    water = np.array(Image.open(before / 'water.png').convert('L'))
    rotation = np.asarray(json.loads((page / 'fine' / 'manifest.json').read_text())['rotation_today_to_turned'], float)
    with np.load(a.average, allow_pickle=False) as z:
        average = {k: z[k].astype(np.float32) for k in ('colour', 'water', 'covered')}
    colour, wet, share = compose(albedo, water, average, rotation, caps=a.match_caps)
    was, now = water.astype(np.float32) / 255.0, wet.astype(np.float32) / 255.0
    Image.fromarray(colour).save(textures / 'albedo.jpg', format='JPEG', quality=QUALITY, subsampling=0, optimize=True)
    Image.fromarray(wet).save(textures / 'water.png', format='PNG', optimize=True)
    manifest_path = page / 'fine' / 'manifest.json'                 # the page's note on the look says what the far view is
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('look') and FAR_NOTE.strip() not in manifest['look'].get('note', ''):
        manifest['look']['note'] = manifest['look'].get('note', '') + FAR_NOTE
        manifest_path.write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps(dict(share_of_texture_from_the_look=round(share, 4),
                          water_share_before=round(float(was.mean()), 4), water_share_after=round(float(now.mean()), 4),
                          water_agreement=round(float(1.0 - np.abs(was - now).mean()), 4), kept_before=str(before),
                          caps='the look\'s last ring by longitude, with the earlier detail' if a.match_caps else 'as they were')))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
