#!/usr/bin/env python3
"""Export the whole Moon at the scenario water share for the Unreal game: heights, water depth and the
land-sea mask as equirectangular images, with the sea level and projection in moon.json.

Usage: python bake/game/moon.py --out <dir> [--ppd 4|16]

Heights and water are computed: they come from the geography atlas product (LOLA above the GRAIL geoid,
filled to one common sea level at the share in shared/scenarios/water.json). Surface cover (soil, rock,
vegetation, colour) is a placeholder for the game until the biosphere and climate supply it.
manifest.json records the product hashes.

Images put 180 degrees west at the left edge and 90 degrees north at the top, one pixel per 1/ppd degree:
  moon_height.png       16-bit; metres above the geoid = value * height_scale + height_offset
  moon_water_depth.png  16-bit; metres of water = value * depth_scale, zero on land
  moon_land_sea.png     8-bit; 255 water, 0 land
At --ppd 16 the heights come from the 16 pixel/degree LOLA grid through the geography domain.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
IMMERSION = HERE.parent.parent
REPO = IMMERSION.parent
sys.path.insert(0, str(REPO))

from geography import atlas as ga  # noqa: E402
from geography import topography as tp  # noqa: E402
from shared.constants import MOON_RADIUS  # noqa: E402

SCHEMA = 'terluna.immersion.game-moon/1'


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_png(path: Path, array: np.ndarray) -> None:
    from PIL import Image
    Image.fromarray(array).save(path, format='PNG', optimize=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--ppd', type=int, choices=(4, 16), default=4)
    args = parser.parse_args(argv)

    atlas_path = REPO / 'geography' / 'results' / 'atlas.json'
    atlas = json.loads(atlas_path.read_text())
    if atlas.get('schema') != ga.SCHEMA:
        raise SystemExit(f'Expected {ga.SCHEMA} in {atlas_path}; run python -m geography.atlas')
    grid_path = REPO / 'geography' / 'products' / f"atlas_{round(100 * atlas['share'])}pct_4ppd.npz"
    if not grid_path.is_file():
        raise SystemExit(f'{grid_path.name} is missing; run python -m geography.atlas')
    level = float(atlas['sea_level_m'])
    if args.ppd == 4:
        height = np.load(grid_path)['height_m'].astype(np.float64)
    else:
        height, _, _ = tp.height_above_geoid(16, atlas['grid']['geoid_degree'])

    height = np.roll(height, height.shape[1] // 2, axis=1)
    depth = np.maximum(level - height, 0.0)
    water = depth > 0

    offset = float(np.floor(height.min()))
    height_scale = float((height.max() - offset) / 65535.0)
    depth_scale = float(max(depth.max(), 1.0) / 65535.0)
    args.out.mkdir(parents=True, exist_ok=True)
    files = {
        'moon_height.png': np.round((height - offset) / height_scale).astype(np.uint16),
        'moon_water_depth.png': np.round(depth / depth_scale).astype(np.uint16),
        'moon_land_sea.png': np.where(water, 255, 0).astype(np.uint8),
    }
    for name, array in files.items():
        write_png(args.out / name, array)

    moon = dict(
        schema=SCHEMA,
        share=atlas['share'], sea_level_m=level, moon_radius_m=MOON_RADIUS,
        projection=('equirectangular; column 0 starts at 180 degrees west, row 0 at 90 degrees north; '
                    f'{args.ppd} pixels per degree; longitudes east'),
        size=[int(height.shape[1]), int(height.shape[0])],
        height=dict(file='moon_height.png', scale_m=height_scale, offset_m=offset,
                    rule='metres above the geoid = value * scale_m + offset_m'),
        water_depth=dict(file='moon_water_depth.png', scale_m=depth_scale, rule='metres of water = value * scale_m'),
        land_sea=dict(file='moon_land_sea.png', rule='255 water, 0 land'),
        kinds=dict(heights='computed', water='computed',
                   surface_cover='placeholder until the biosphere and climate supply it'),
        evidence=atlas['evidence'],
    )
    (args.out / 'moon.json').write_text(json.dumps(moon, indent=2) + '\n')
    manifest = dict(
        exporter=str(Path(__file__).relative_to(REPO)), exporter_sha256=sha256(Path(__file__)),
        sources={str(atlas_path.relative_to(REPO)): sha256(atlas_path),
                 str(grid_path.relative_to(REPO)): sha256(grid_path)},
        outputs={name: sha256(args.out / name) for name in [*files, 'moon.json']},
    )
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"{atlas['share']:.0%} Moon at {args.ppd} px/deg: {height.shape[1]}x{height.shape[0]}, "
          f"sea level {level:,.0f} m, water {water.mean():.1%} of pixels -> {args.out}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
