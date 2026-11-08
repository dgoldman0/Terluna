"""A turned run's climatology carried back onto today's latitude and longitude, for the lanes that work on today's
ground, such as the geography domain's drainage routing (7 November 2026).

    python -m climate.gcm.turned_climatology <orientation product> products/climatology_<run>.npz \
        --source <branch@commit:path>          # writes products/climatology_<run>_today.npz

A run on the turned Moon (climate.gcm.turned) keeps its fields in the turned frame. The ground and its water turn with
the crust, so a field that belongs to the ground, such as the rain that falls on a slope, belongs at the ground's place
in today's frame as the atlas, the terrain and the drainage hold it. Each cell of today's T21 grid takes the turned
field at the cell's turned place, bilinear in latitude and longitude on the turned grid (longitude periodic, latitudes
beyond the outermost Gaussian rows taking those rows). For a smooth field that keeps it within 0.25% of the field's
range away from the poles and within 2% in the two outermost rows of each side. The Sun composites stay in the turned
frame's latitude and are left out. The product keeps the climatology's schema, so any reader of a climatology reads it, and its metadata name
the frame, the turned climatology and the orientation product it was carried with.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from climate.gcm import turned

SCHEMA = 'terluna.climate.gcm-climatology/1'
LATLON = ('clt_mean', 'low_mean', 'high_mean', 'pr_mm_day', 'evap_mm_day', 'mrso_m', 'ts_k', 'prw_kg_m2', 'sic', 'snd_m',
          'lsm')


def sample(field, lat_src, lon_src, lat_pts, lon_pts):
    """Bilinear values of a field on a latitude-longitude grid (latitudes descending, longitudes evenly spaced from
    lon_src[0] and periodic) at the given points; latitudes beyond the outermost rows take those rows."""
    lat_src, lon_src = np.asarray(lat_src, float), np.asarray(lon_src, float)
    rows = np.interp(-np.asarray(lat_pts, float), -lat_src, np.arange(lat_src.size))
    cols = ((np.asarray(lon_pts, float) - lon_src[0]) / (lon_src[1] - lon_src[0])) % lon_src.size
    r0 = np.clip(np.floor(rows).astype(int), 0, lat_src.size - 2)
    fr = np.clip(rows - r0, 0.0, 1.0)
    c0 = np.floor(cols).astype(int) % lon_src.size
    c1 = (c0 + 1) % lon_src.size
    fc = cols - np.floor(cols)
    f = np.asarray(field, float)
    top = f[r0, c0] * (1 - fc) + f[r0, c1] * fc
    bottom = f[r0 + 1, c0] * (1 - fc) + f[r0 + 1, c1] * fc
    return top * (1 - fr) + bottom * fr


def to_today(clim, rotation):
    """The latitude-longitude fields of a turned climatology on the same grid in today's frame."""
    lat, lon = np.asarray(clim['lat'], float), np.asarray(clim['lon'], float)
    grid_lat, grid_lon = np.meshgrid(lat, lon, indexing='ij')
    lat_t, lon_t = turned.to_turned(grid_lat, grid_lon, rotation)
    return {k: sample(clim[k], lat, lon, lat_t, lon_t).astype(np.float32) for k in LATLON}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('orientation', type=Path, help="domain/geology's orientation product (moon_orientation.json)")
    parser.add_argument('climatology', type=Path, help="a turned run's climatology (climatology_<run>.npz)")
    parser.add_argument('--source', required=True, help='where the orientation product was read, as branch@commit:path')
    args = parser.parse_args(argv)
    frame = turned.read_frame(args.orientation)
    z = np.load(args.climatology)
    meta = json.loads(str(z['metadata']))
    if meta.get('schema') != SCHEMA:
        raise SystemExit(f'Expected {SCHEMA} in {args.climatology}')
    clim = {k: z[k] for k in z.files if k != 'metadata'}
    fields = to_today(clim, frame['rotation'])
    meta.update(frame='today', carried=dict(
        method=('each cell of today\'s grid takes the turned field at its turned place, bilinear on the turned grid; '
                'the Sun composites are left out'),
        climatology=dict(file=args.climatology.name, sha256=digest(args.climatology)),
        orientation=dict(source=args.source, sha256=digest(args.orientation)),
        producer={'climate/gcm/turned_climatology.py': digest(Path(__file__))[:16],
                  'climate/gcm/turned.py': digest(Path(turned.__file__))[:16]}))
    out = args.climatology.with_name(args.climatology.stem + '_today.npz')
    np.savez_compressed(out, metadata=json.dumps(meta), lat=clim['lat'], lon=clim['lon'], **fields)
    weights = np.cos(np.radians(np.asarray(clim['lat'], float)))[:, None] * np.ones_like(fields['pr_mm_day'])
    mean = lambda f: float((f * weights).sum() / weights.sum())  # noqa: E731
    print(f"{out.name}: rain {mean(fields['pr_mm_day']):.3f} mm/day against {mean(clim['pr_mm_day']):.3f} in the turned "
          f"frame; land {mean(fields['lsm']):.3f} against {mean(clim['lsm']):.3f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
