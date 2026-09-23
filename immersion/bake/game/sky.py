#!/usr/bin/env python3
"""Export the sky for the Unreal game: engine-sky parameters, the clear-sky atlas as
images, the site ephemeris, the stars, Earth's imagery and the cloud columns.

Usage: python bake/game/sky.py --out <dir>

Every file comes from a domain product or the experience's baked copy of one;
manifest.json records which, with hashes. Nothing here computes new physics.

The atlas images hold, for each of 174 Sun elevations, the clear sky's radiance
(cd/m2, photopic-weighted linear sRGB) over 33 view azimuths (0 = toward the Sun,
32 = away, symmetric) and 41 view elevations (tile row r from the top: 0 = straight
up, 20 = the horizon, 40 = straight down; elevation = sign(t) * t^2 * 90 degrees with
t = (20 - r) / 20). Frame f sits in tile (f % 14, f // 14) from the top left of a
14 x 13 grid.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
IMMERSION = HERE.parent.parent
REPO = IMMERSION.parent
sys.path.insert(0, str(REPO))

from illumination.sky.engine_atmosphere import product as engine_atmosphere  # noqa: E402

COLUMNS, ROWS = 14, 13


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unpack(world: dict) -> np.ndarray:
    """The atlas frames of one world as float32 cd/m2, shape (frames, 41, 33, 3)."""
    raw = gzip.decompress(base64.b64decode(world['data']))
    halves = np.frombuffer(raw, dtype='<f2').astype(np.float32)
    w, h, n = world['width'], world['height'], len(world['suns'])
    frames = halves.reshape(n, h, w, 3)
    return frames * np.asarray(world['scale'], dtype=np.float32)[:, None, None, None]


def write_hdr(path: Path, rgb: np.ndarray) -> None:
    """Radiance RGBE (.hdr), rows from the top."""
    h, w, _ = rgb.shape
    peak = rgb.max(axis=2)
    mantissa, exponent = np.frexp(peak)
    scale = np.where(peak > 1e-32, mantissa * 256.0 / np.maximum(peak, 1e-32), 0.0)
    rgbe = np.zeros((h, w, 4), dtype=np.uint8)
    rgbe[..., :3] = np.clip(np.floor(rgb * scale[..., None]), 0, 255).astype(np.uint8)
    rgbe[..., 3] = np.where(peak > 1e-32, exponent + 128, 0).astype(np.uint8)
    with open(path, 'wb') as f:
        f.write(b'#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n')
        f.write(f'-Y {h} +X {w}\n'.encode())
        f.write(rgbe.tobytes())


def atlas_image(frames: np.ndarray) -> np.ndarray:
    n, h, w, _ = frames.shape
    image = np.zeros((ROWS * h, COLUMNS * w, 3), dtype=np.float32)
    for f in range(n):
        r, c = divmod(f, COLUMNS)
        image[r * h:(r + 1) * h, c * w:(c + 1) * w] = frames[f][::-1]  # straight up at the top
    return image


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', required=True, type=Path)
    out = ap.parse_args().out.resolve() / 'sky'
    out.mkdir(parents=True, exist_ok=True)
    sources = []

    def source(path: Path) -> None:
        sources.append({'file': str(path.relative_to(REPO)), 'sha256': sha256(path)})

    # Engine-sky parameters from the illumination domain.
    (out / 'engine_atmosphere.json').write_text(json.dumps(engine_atmosphere(), indent=1))
    source(REPO / 'illumination/sky/atmospheres.py')

    # The clear-sky atlas, per world, as images and exact floats plus its tables.
    atlas_path = IMMERSION / 'assets/sky/atmosphere.json'
    atlas = json.loads(atlas_path.read_text())
    source(atlas_path)
    worlds = {}
    for name, world in atlas['worlds'].items():
        frames = unpack(world)
        if frames.shape[0] > COLUMNS * ROWS:
            raise ValueError('More atlas frames than the tile grid holds')
        image = atlas_image(frames)
        write_hdr(out / f'atlas_{name}.hdr', image)
        (out / f'atlas_{name}.f32').write_bytes(image.astype('<f4').tobytes())
        worlds[name] = {
            'image': f'atlas_{name}.hdr',
            'floats': f'atlas_{name}.f32 (float32 little-endian RGB, rows from the top)',
            'size_px': [image.shape[1], image.shape[0]],
            'tile_px': [world['width'], world['height']],
            'grid': [COLUMNS, ROWS],
            'sun_elevations_deg': world['suns'],
            'radius_m': world['radius'],
            # Irradiance at the ground (lux) and cloud-altitude light, per Sun elevation.
            **{k: world[k] for k in ('direct', 'direct_horizontal', 'diffuse', 'cloud_direct', 'cloud_diffuse')},
        }
    (out / 'atlas.json').write_text(json.dumps({
        'units': atlas['units'],
        'layout': __doc__.split('The atlas images hold, ')[1].split('\n"""')[0].replace('\n', ' '),
        'provenance': atlas['provenance'],
        'worlds': worlds,
    }, indent=1))

    # The site ephemeris, stars, Earth imagery and cloud columns, as published.
    for rel in ['assets/sky/site-sky.json', 'assets/sky/stars.json', 'assets/earth/manifest.json',
                'assets/earth/calibration.json', 'assets/earth/land_shallow_topo_2048.jpg',
                'assets/earth/cloud_combined_2048.jpg']:
        path = IMMERSION / rel
        shutil.copyfile(path, out / path.name)
        source(path)
    columns = REPO / 'atmosphere/column/products/columns.json'
    if not columns.exists():
        raise FileNotFoundError('Run node atmosphere/column/export_columns.cjs first')
    shutil.copyfile(columns, out / 'columns.json')
    source(columns)

    files = [{'file': p.name, 'sha256': sha256(p)} for p in sorted(out.iterdir()) if p.name != 'manifest.json']
    (out / 'manifest.json').write_text(json.dumps({
        'schema': 'terluna.game.sky/1',
        'exporter': 'immersion/bake/game/sky.py',
        'sources': sources,
        'files': files,
    }, indent=1))
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
