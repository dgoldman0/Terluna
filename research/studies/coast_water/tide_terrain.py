"""LOLA's 256 px/deg rows across the deserts' sea coasts at 45-60 degrees, for the tide study's profiles at 100 m.

    python -m research.studies.coast_water.tide_terrain --download     # about 95 MB by HTTP byte range

The coast-water study's four local coasts lie at 2-37 degrees, outside the deserts. These slices take one degree of
latitude across half the Moon's longitudes each, where the linear rain mask's desert sea coasts at 16 px/deg are
densest at 45-60 degrees (north 48-49 and 60-61 degrees in the eastern half, 50-51 degrees north in the western,
57-58 degrees south in the eastern), from the archived tiles geography's coast inputs came from. Each slice is fetched
by an exact byte range with the tile's etag, its size checked and its sha256 recorded in tide_inputs.json; a failed
fetch or a changed hash drops that slice. The bytes go to research/runs/coast_water/tide (ignored, on the research
drive). Heights above the atlas's sea level on local grids follow geography.coast_terrain.local_grid: bilinear
between LOLA's pixel centres, less the degree-200 GRAIL geoid on a 0.05-degree grid.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RECORD = HERE / 'tide_inputs.json'
STORE = 'research/runs/coast_water/tide'
BASE = 'https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/'
TILES = {'00n_90n_000_180': dict(line_projection_offset=23039.5, sample_projection_offset=46079.5),
         '00n_90n_180_360': dict(line_projection_offset=23039.5, sample_projection_offset=-0.5),
         '90s_00s_000_180': dict(line_projection_offset=-0.5, sample_projection_offset=46079.5)}
SLICES = (('n48_e', '00n_90n_000_180', 48.0, 49.0), ('n60_e', '00n_90n_000_180', 60.0, 61.0),
          ('n50_w', '00n_90n_180_360', 50.0, 51.0), ('s57_e', '90s_00s_000_180', -58.0, -57.0))
PPD, SAMPLES, LINES, ROW_BYTES, TILE_BYTES = 256, 46080, 23040, 92160, 2123366400
CREDIT = ('NASA LRO LOLA Science Team, NASA Goddard Space Flight Center; archived by the PDS Geosciences Node at '
          'Washington University in St. Louis (Smith et al. 2010, Space Sci. Rev. 150, 209).')


def digest(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 22), b''):
            h.update(block)
    return h.hexdigest()


def entry(name, tile, south, north):
    """The slice's rows (whole degrees with a row's margin) and its byte range in the tile."""
    lo = TILES[tile]['line_projection_offset']
    first = int(math.floor(lo - north * PPD)) - 1
    last = int(math.ceil(lo - south * PPD)) + 1
    first, last = max(first, 0), min(last, LINES - 1)
    count = last - first + 1
    return dict(name=name, tile=tile, source_url=f'{BASE}ldem_256_{tile}.img', latitude_deg=[south, north],
                first_row=first, row_count=count, row_bytes=ROW_BYTES, samples=SAMPLES, pixels_per_degree=PPD,
                selected_byte_range=[first * ROW_BYTES, (first + count) * ROW_BYTES - 1],
                selected_bytes=count * ROW_BYTES, sample_type='<i2', scale_m=0.5, center_longitude_deg=180.0,
                **TILES[tile], local_image=f'{STORE}/ldem_256_{tile}_rows_{first}_{first + count - 1}.img')


def fetch(e):
    """One slice by byte range; returns its record, or the reason it was dropped."""
    path = ROOT / e['local_image']
    path.parent.mkdir(parents=True, exist_ok=True)
    head = urlopen(Request(e['source_url'], method='HEAD'), timeout=60)
    etag, size = head.headers.get('ETag'), int(head.headers.get('Content-Length', -1))
    if size != TILE_BYTES:
        return dict(e, dropped=f'the tile is {size} bytes, not {TILE_BYTES}')
    first, last = e['selected_byte_range']
    partial = path.with_suffix('.partial')
    with urlopen(Request(e['source_url'], headers={'Range': f'bytes={first}-{last}', 'If-Range': etag}),
                 timeout=120) as r:
        if r.status != 206 or r.headers.get('Content-Range') != f'bytes {first}-{last}/{TILE_BYTES}':
            return dict(e, dropped='the server did not return the exact byte range')
        with open(partial, 'wb') as f:
            remaining = e['selected_bytes']
            while remaining:
                chunk = r.read(min(1 << 20, remaining))
                if not chunk:
                    break
                f.write(chunk)
                remaining -= len(chunk)
    if partial.stat().st_size != e['selected_bytes']:
        partial.unlink()
        return dict(e, dropped='incomplete download')
    partial.rename(path)
    return dict(e, source_etag=etag, source_last_modified=head.headers.get('Last-Modified'),
                selected_sha256=digest(path), retrieved='2026-11-04')


def load_record():
    return json.loads(RECORD.read_text()) if RECORD.exists() else None


def restore():
    """The recorded slices, each checked against its size and sha256; a mismatch is an error."""
    record = load_record()
    if record is None:
        raise FileNotFoundError('Fetch the slices with --download first')
    out = []
    for e in record['slices']:
        if 'dropped' in e:
            continue
        path = ROOT / e['local_image']
        if not path.exists() or path.stat().st_size != e['selected_bytes'] or digest(path) != e['selected_sha256']:
            raise ValueError(f'{path} differs from its record')
        out.append(e)
    return out


def radius_heights(e, latitude, longitude):
    """LOLA heights above 1737.4 km, bilinear between pixel centres (geography.coast_terrain.radius_heights)."""
    raw = np.memmap(ROOT / e['local_image'], dtype='<i2', mode='r', shape=(e['row_count'], e['samples']))
    row = e['line_projection_offset'] - np.asarray(latitude) * PPD - e['first_row']
    col = e['sample_projection_offset'] + (np.asarray(longitude) % 360 - e['center_longitude_deg']) * PPD
    if row.min() < 0 or row.max() > e['row_count'] - 1 or col.min() < 0 or col.max() > e['samples'] - 1:
        raise ValueError(f'The local grid leaves the rows of {e["name"]}')
    r0, r1 = int(np.floor(row.min())), int(np.ceil(row.max())) + 1
    c0, c1 = int(np.floor(col.min())), int(np.ceil(col.max())) + 1
    block = np.asarray(raw[r0:r1, c0:c1], dtype=float) * e['scale_m']
    return map_coordinates(block, [row - r0, col - c0], order=1, mode='nearest')


def local_grid(e, lon0, lat0, half_width_m, spacing_m=100.0):
    """Heights (m) above the atlas's sea level on a local grid round (lon0, lat0), as geography.coast_terrain
    builds them; x east and y north (m) by the azimuthal equidistant projection."""
    from geography import coast_terrain as CT
    from geography.topography import Grid, geoid
    axis = np.arange(-half_width_m, half_width_m + spacing_m / 2, spacing_m)
    x, y = np.meshgrid(axis, axis, indexing='xy')
    lat, lon = CT.destination(lon0, lat0, x, y)
    radius = radius_heights(e, lat, lon)
    step = 0.05
    glat = np.arange(lat.max() + step, lat.min() - step, -step)
    glon = np.arange(lon.min() - step, lon.max() + step, step)
    datum = geoid(Grid(glat, glon, int(round(1 / step))), max_degree=200)
    geoid_height = map_coordinates(datum, [(glat[0] - lat) / step, (lon - glon[0]) / step], order=1, mode='nearest')
    with np.load(CT.ATLAS, allow_pickle=False) as atlas:
        level = float(json.loads(atlas['metadata'].item())['sea_level_m'])
    return axis, (radius - geoid_height - level).astype(np.float32), level


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args(argv)
    if not args.download:
        print(json.dumps([e['name'] for e in restore()]))
        return 0
    slices = []
    for name, tile, south, north in SLICES:
        e = entry(name, tile, south, north)
        try:
            rec = fetch(e)
        except OSError as err:                                   # dropped, never substituted
            rec = dict(e, dropped=f'fetch failed: {err}')
        slices.append(rec)
        print(name, rec.get('selected_sha256') or rec['dropped'])
    RECORD.write_text(json.dumps(dict(
        schema='terluna.research.coast-tide-inputs/1',
        status='external inputs; fetched by HTTP byte range and hash-checked by tide_terrain.py; ignored by Git',
        credit=CREDIT, rights='NASA scientific data; retain source attribution. The bytes stay in ignored storage.',
        evidence=('Archived altimetry-derived raster (LOLA GDR, 256 pixels per degree, about 118 m at the equator '
                  'and 59-69 m east-west at 45-60 degrees), with interpolation for incomplete coverage.'),
        slices=slices), indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
