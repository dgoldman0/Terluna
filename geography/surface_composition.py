"""The surface's iron and titanium from Lunar Prospector's gamma rays: FeO and TiO2 (weight per cent) on a grid.

    python -m geography.surface_composition      # the oxides at a few known places and their spread

Prettyman et al. (2006, JGR 111, E12007) unmixed Lunar Prospector's gamma-ray spectra into the major oxides'
weight fractions. The table read here holds them on 11,306 equal-area pixels 2 degrees across from the high-altitude
phase (about 100 km), where the instrument saw about 5 degrees of arc, so the 2-degree pixels sample a smoother field:
the maria's 15-21% FeO, the feldspathic highlands' 3-5% and the South Pole-Aitken basin's floor between them, and the
high-titanium maria (Tranquillitatis near 9% TiO2). fetch_surface_inputs.py restores the table, hash-checked.
grid() puts each pixel's values on a latitude-longitude grid and smooths them over the instrument's footprint, in
today's frame.
"""
from __future__ import annotations

import numpy as np

from geography import fetch_surface_inputs as fetch_inputs

TABLE = 'lpgrs_high1_elem_abundance_2deg.tab'
REFERENCE = ('Prettyman, T. H., Hagerty, J. J., Elphic, R. C., Feldman, W. C., Lawrence, D. J., McKinney, G. W. and '
             'Vaniman, D. T. (2006), Elemental composition of the lunar surface: analysis of gamma ray spectroscopy '
             'data from Lunar Prospector, Journal of Geophysical Research 111, E12007')
FOOTPRINT_DEG = 2.0       # Gaussian smoothing (standard deviation) for the 5-degree footprint the 2-degree pixels sample


def pixels() -> dict:
    """The table's pixels: bounds (degrees, east longitude -180 to 180) and FeO and TiO2 in weight per cent."""
    rows = np.loadtxt(fetch_inputs.path(TABLE), usecols=(1, 2, 3, 4, 11, 12))
    return dict(lat0=rows[:, 0], lat1=rows[:, 1], lon0=rows[:, 2], lon1=rows[:, 3],
                tio2=100.0 * rows[:, 4], feo=100.0 * rows[:, 5])


def grid(ppd: int = 1, left_lon: float = -180.0, smooth_deg: float = FOOTPRINT_DEG, exclude=None) -> dict:
    """FeO and TiO2 (weight per cent) at each cell's centre: rows from 90 N, columns eastward from left_lon, each cell
    taking the pixel that holds its centre, then smoothed over smooth_deg (round on the ground, wrapping in longitude).
    exclude: cells (True), such as the maria, whose oxides the result leaves out: each cell then counts by the share of
    its surroundings outside them, smoothed over the same footprint, to the fourth power, since the instrument mixed
    a pixel's surroundings into it; cells whose footprint saw mostly maria count for little, and the rest fill in."""
    from scipy.ndimage import gaussian_filter1d
    p = pixels()
    rows, cols = 180 * ppd, 360 * ppd
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    lon = ((left_lon + (np.arange(cols) + 0.5) / ppd) + 180.0) % 360.0 - 180.0
    index = np.full((rows, cols), -1, np.int64)
    for low in np.unique(p['lat0']):
        band = np.flatnonzero(p['lat0'] == low)
        band = band[np.argsort(p['lon0'][band])]
        inside = (lat >= low) & (lat < p['lat1'][band[0]])
        if not inside.any():
            continue
        at = np.clip(np.searchsorted(p['lon0'][band], lon, side='right') - 1, 0, band.size - 1)
        index[inside] = band[at][None, :]
    if (index < 0).any():
        raise ValueError(f'{int((index < 0).sum())} cells lie in no pixel of {TABLE}')
    across = smooth_deg * ppd / np.maximum(np.cos(np.radians(lat)), 0.05)
    response = np.exp(-2.0 * np.pi ** 2 * (across[:, None] * np.fft.rfftfreq(cols)[None, :]) ** 2)

    def smooth(f):
        f = gaussian_filter1d(f, smooth_deg * ppd, axis=0, mode='nearest')
        return np.fft.irfft(np.fft.rfft(f, axis=1) * response, n=cols, axis=1)

    weight = (np.ones((rows, cols)) if exclude is None
              else np.clip(1.0 - smooth(np.asarray(exclude, float)), 0.0, 1.0) ** 4 * ~np.asarray(exclude, bool))
    total = np.maximum(smooth(weight), 1e-6)
    return {name: (smooth(p[name][index] * weight) / total).astype(np.float32) for name in ('feo', 'tio2')}


def main() -> int:
    g = grid(4, left_lon=0.0)
    at = lambda name, la, lo: float(g[name][int((90 - la) * 4), int((lo % 360) * 4)])
    for place, la, lo in (('Mare Imbrium', 32.8, -15.6), ('Mare Tranquillitatis', 8.5, 31.4),
                          ('South Pole-Aitken floor', -50.0, -170.0), ('far-side highlands', 10.0, 175.0)):
        print(f"{place}: FeO {at('feo', la, lo):.1f}%, TiO2 {at('tio2', la, lo):.1f}%")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
