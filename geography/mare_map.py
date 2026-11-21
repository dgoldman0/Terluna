"""Where the maria's basalt lies: the LROC team's mapped mare boundaries on a latitude-longitude grid.

    python -m geography.mare_map            # the basalt's share of the Moon and its largest maria

The boundaries are LROC_GLOBAL_MARE_180 (Nelson et al. 2014, LPSC 45, abstract 2861): 644 polygons digitized at
1:200,000 to 1:500,000 on the LROC WAC 643 nm mosaic at 100 m per pixel, a WAC band-ratio mosaic that picks out the
high-titanium basalts and a Clementine colour ratio that separates the old highlands from the younger maria. They
reach from 66.5 S to 64.8 N; no mare is mapped poleward of them. fetch_surface_inputs.py restores the archive,
hash-checked, and this module reads the shapefile from it. A pixel counts as basalt where its centre lies inside a
mapped mare, in today's frame (the IAU Moon 2000 coordinates the polygons and LOLA share).
"""
from __future__ import annotations

import io
import zipfile

import numpy as np

from geography import fetch_surface_inputs as fetch_inputs

ARCHIVE = 'LROC_GLOBAL_MARE_180.ZIP'
LAYER = 'LROC_GLOBAL_MARE_180'
REFERENCE = ('Nelson, D. M., Koeber, S. D., Daud, K., Robinson, M. S., Watters, T. R., Banks, M. E. and Williams, N. R. '
             '(2014), Mapping lunar maria extents and lobate scarps using LROC image products, Lunar and Planetary '
             'Science Conference 45, abstract 2861')


def outlines() -> list[tuple[dict, str]]:
    """Each mapped mare's polygon (GeoJSON-like, longitudes -180 to 180) and its name."""
    import shapefile
    with zipfile.ZipFile(fetch_inputs.path(ARCHIVE)) as z:
        part = lambda ext: io.BytesIO(z.read(f'{LAYER}.{ext}'))
        reader = shapefile.Reader(shp=part('SHP'), shx=part('SHX'), dbf=part('DBF'))
        names = [r['MARE_NAME'] for r in reader.records()]
        return [(s.__geo_interface__, n) for s, n in zip(reader.shapes(), names)]


def basalt(ppd: int, left_lon: float = -180.0) -> np.ndarray:
    """True where a mapped mare covers the pixel's centre: rows from 90 N, columns eastward from left_lon."""
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    grid = rasterize(((shape, 1) for shape, _ in outlines()), out_shape=(180 * ppd, 360 * ppd),
                     transform=from_origin(-180.0, 90.0, 1.0 / ppd, 1.0 / ppd), fill=0, dtype='uint8')
    return np.roll(grid.astype(bool), -int(round((left_lon + 180.0) * ppd)), axis=1)


def area_share(mask: np.ndarray) -> float:
    """A mask's share of the sphere, rows from 90 N."""
    rows = mask.shape[0]
    weight = np.cos(np.radians(90.0 - (np.arange(rows) + 0.5) * 180.0 / rows))
    return float((mask.sum(axis=1) * weight).sum() / (weight.sum() * mask.shape[1]))


def main() -> int:
    grid = basalt(16)
    print(f'basalt covers {area_share(grid):.1%} of the Moon at 16 px/deg ({len(outlines())} polygons)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
