"""Supplement the 24 coastal frames with southern Nubium at local midnight.

Consumes the same sky, tide, sea-slope and SWAN products as scenes.py. The
single-process runner temporarily supplies a portrait camera to that existing
calculation; it does not alter or regenerate the original 24-scene product.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from research.studies.sea_appearance import lighting, scenes


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build():
    inputs = scenes.Inputs()
    site = next(s for s in inputs.calendar['sites'] if s['short'] == 'S Nubium')
    series = site['series']
    # Lowest Sun in this month at the study station, on the hourly grid.
    i = int(np.argmin(series['sun_elevation_deg']))
    hour = float(inputs.hours[i])
    record = dict(coast='S Nubium', moment='midnight', hour=hour,
                  light_on_ground_lux=float(series['ground_lux_from_sun'][i]
                                            + series['ground_lux_from_earth'][i]
                                            + lighting.FLOOR_LUX))
    old = (scenes.WIDTH, scenes.HEIGHT, scenes.FOV_DEG)
    scenes.WIDTH, scenes.HEIGHT, scenes.FOV_DEG = 1200, 1800, 65.0
    scenes.VIEWS[('S Nubium', 'midnight')] = (22.0, 16.0)
    try:
        frame = scenes.scene(inputs, 0, record)
    finally:
        scenes.WIDTH, scenes.HEIGHT, scenes.FOV_DEG = old
        del scenes.VIEWS[('S Nubium', 'midnight')]
    source_paths = [Path(__file__), Path(scenes.__file__), Path(lighting.__file__),
                    scenes.CALENDAR, scenes.SLOPE_SERIES, scenes.SOLUTION,
                    scenes.TIDES, scenes.WATERS, scenes.ROOT/frame['waves']['file'],
                    scenes.TERRAIN/'s_nubium.npz', scenes.ROOT/'shared/constants.json',
                    scenes.ROOT/'research/studies/sea_appearance/regimes.py',
                    scenes.ROOT/'illumination/earthlight/model.py',
                    scenes.ROOT/'illumination/water_surface/reflection.py',
                    scenes.ROOT/'illumination/water_column/model.py',
                    scenes.ROOT/'geography/lunar_ephemeris.py']
    return dict(
        schema='terluna.research.nubium-midnight/1',
        producer={'model':'research/studies/sea_appearance/nubium_midnight.py',
                  'sha256':digest(__file__)},
        evidence='Existing sea-appearance calculation at the hourly minimum of solar elevation, with the same offshore camera position and a new portrait view toward Earth. No new atmospheric or wave simulation.',
        reading_rule='A supplemental scene, not SN-6 (the minimum-illuminance frame). Calendar ground illuminance is at the study station; directional colours use the offshore camera. The nearest available SWAN restart is identified in waves. Clear molecular atmosphere, illustrative productive-coast composition and bare-soil land placeholders; disk lighting still uses a uniform Earth. Texture and exposure blending belong to the illustration, not this radiance product.',
        inputs={str(p.relative_to(scenes.ROOT)):digest(p) for p in source_paths},
        jd_tt=inputs.jd(hour),
        selection={'criterion':'minimum solar elevation on the hourly lighting calendar',
                   'hour':hour, 'station_sun_lux':series['ground_lux_from_sun'][i],
                   'station_earth_lux':series['ground_lux_from_earth'][i]},
        scene=frame)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=Path(__file__).with_name('results')/'nubium-midnight.json')
    args = p.parse_args()
    product = build()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(product, indent=2)+'\n')
    s = product['scene']
    print(json.dumps({k:s[k] for k in ['hour','camera','light','earth','waves','land','stars']}, indent=2))


if __name__ == '__main__':
    main()
