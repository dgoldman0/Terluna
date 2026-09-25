"""Surface climatology of an ExoPlaSim run, written as a data product for display and other lanes.

    climate/gcm/.venv/bin/python -m climate.gcm.climatology A:30-39     # products/climatology_A.npz

For the chosen model years it averages the 3-day output means of cloud cover, precipitation,
evaporation, soil moisture, surface temperature, water vapour, sea ice, snow and the land mask, on the
model grid. It splits cloud into a low deck (model layers with sigma above 0.5, about the lowest 34 km) and a
high deck (the layers above), each the largest layer cover in its column (maximum overlap), and gives
the mean cover of every layer.
It also composites total, low and high cloud against the Sun, in latitude and local hour angle (degrees
east of the subsolar point), so that a display can move clouds with the Sun as the slowly rotating Moon
does. Each output interval is a 3-day mean, a tenth of a lunar day, so the hour angle is resolved to
about 36 degrees.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.climate.gcm-climatology/1'
FIELDS = ('clt', 'pr', 'evap', 'mrso', 'ts', 'prw', 'sic', 'snd', 'lsm', 'czen')
HOUR_BINS = 36
DECK_SIGMA = 0.5   # layers above this sigma form the low deck, layers below it the high deck
EVIDENCE = ('Means of ExoPlaSim 3-day output means over the stated model years at T21 (32 x 64). PlaSim diagnoses '
            'cloud cover from humidity with Earth-tuned schemes, and its cloud radiative effect is the weak point of the '
            'run (see README, "Clouds"); cover and precipitation serve here as a climatology for display and for '
            'first estimates of surface cover.')
READING_RULE = ('Latitudes north first as in the model; longitudes east, 0-360. Precipitation and evaporation in mm per '
                'day of liquid water, evaporation positive upward. clt_sun[lat, k] is mean cloud cover at hour angle '
                'hour_angle_deg[k] (degrees east of the subsolar point, -180..180); low_sun and high_sun likewise for '
                'the two decks. cl_layer_mean[k] is the area-weighted mean cover of the layer at sigma[k].')


def decks(cl, sigma, split=DECK_SIGMA):
    """Low and high cloud decks from layer cover (time, layer, lat, lon): the largest cover among the layers
    below and above the split, as maximum overlap within a deck."""
    return cl[:, sigma > split].max(axis=1), cl[:, sigma <= split].max(axis=1)


def read(folder: str, first: int, last: int):
    import netCDF4
    from climate.gcm.exoplasim_run import RUNS
    fields = {k: [] for k in (*FIELDS, 'low', 'high', 'cl_layer')}
    for year in range(first, last + 1):
        with netCDF4.Dataset(RUNS / folder / 'model' / f'MOST.{year:05d}.nc') as d:
            for k in FIELDS:
                fields[k].append(np.asarray(d[k][:], dtype=float))
            sigma = np.asarray(d['lev'][:], dtype=float)
            cl = np.asarray(d['cl'][:], dtype=float)
            low, high = decks(cl, sigma)
            fields['low'].append(low)
            fields['high'].append(high)
            fields['cl_layer'].append(cl.mean(axis=(0, 3))[None])
            lat = np.asarray(d['lat'][:], dtype=float)
            lon = np.asarray(d['lon'][:], dtype=float)
    fields = {k: np.concatenate(v) for k, v in fields.items()}
    fields['sigma'] = sigma
    return fields, lat, lon


def subsolar_longitudes(czen, lat, lon):
    """Longitude of greatest mean Sun height near the equator, for each output interval."""
    rows = np.abs(lat) < 20.0
    return lon[np.argmax(czen[:, rows, :].mean(axis=1), axis=1)]


def sun_composite(clt, czen, lat, lon, bins=HOUR_BINS):
    sub = subsolar_longitudes(czen, lat, lon)
    hour = (lon[None, :] - sub[:, None] + 180.0) % 360.0 - 180.0
    index = np.clip(((hour + 180.0) / 360.0 * bins).astype(int), 0, bins - 1)
    total = np.zeros((lat.size, bins))
    count = np.zeros((lat.size, bins))
    for t in range(clt.shape[0]):
        for j in range(lon.size):
            total[:, index[t, j]] += clt[t, :, j]
            count[:, index[t, j]] += 1
    centres = (np.arange(bins) + 0.5) / bins * 360.0 - 180.0
    return total / np.maximum(count, 1), centres


def climatology(folder: str, first: int, last: int) -> dict:
    f, lat, lon = read(folder, first, last)
    mm_day = 86400.0 * 1000.0
    clt_sun, hour = sun_composite(f['clt'], f['czen'], lat, lon)
    low_sun, _ = sun_composite(f['low'], f['czen'], lat, lon)
    high_sun, _ = sun_composite(f['high'], f['czen'], lat, lon)
    evap = f['evap'].mean(axis=0) * mm_day
    evap = -evap if evap.mean() < 0 else evap
    weight = np.cos(np.radians(lat))
    layer = (f['cl_layer'].mean(axis=0) * weight[None, :]).sum(axis=1) / weight.sum()
    return dict(lat=lat, lon=lon, clt_mean=f['clt'].mean(axis=0), low_mean=f['low'].mean(axis=0),
                high_mean=f['high'].mean(axis=0), pr_mm_day=f['pr'].mean(axis=0) * mm_day,
                evap_mm_day=evap, mrso_m=f['mrso'].mean(axis=0), ts_k=f['ts'].mean(axis=0),
                prw_kg_m2=f['prw'].mean(axis=0), sic=f['sic'].mean(axis=0), snd_m=f['snd'].mean(axis=0),
                lsm=f['lsm'].mean(axis=0), clt_sun=clt_sun, low_sun=low_sun, high_sun=high_sun,
                hour_angle_deg=hour, sigma=f['sigma'], cl_layer_mean=layer)


def main(argv=None) -> int:
    args = list(argv or sys.argv[1:])
    if len(args) != 1 or ':' not in args[0]:
        raise SystemExit('Give one run and year range, for example A:30-39')
    folder, span = args[0].split(':')
    first, last = (int(v) for v in span.split('-'))
    data = climatology(folder, first, last)
    from climate.gcm.exoplasim_run import RUNS
    progress = json.loads((RUNS / folder / 'progress.json').read_text())
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    meta = dict(schema=SCHEMA, run=folder, years=[first, last], configuration=progress['configuration'],
                deck_sigma=DECK_SIGMA,
                producer=dict(domain='climate', files={'gcm/climatology.py': digest(__file__)}),
                evidence=EVIDENCE, reading_rule=READING_RULE)
    out = HERE / 'products' / f'climatology_{folder}.npz'
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, metadata=json.dumps(meta), **{k: v.astype(np.float32) for k, v in data.items()})
    w = np.cos(np.radians(data['lat']))[:, None]
    area = lambda f: float((f * w).sum() / (w.sum() * f.shape[1]))
    print(f"{out}: area-weighted cloud cover {area(data['clt_mean']):.3f} (low deck {area(data['low_mean']):.3f}, high "
          f"{area(data['high_mean']):.3f}), precipitation {area(data['pr_mm_day']):.2f} mm/day, subsolar-noon cloud "
          f"{area(data['clt_sun'][:, HOUR_BINS // 2:HOUR_BINS // 2 + 1]):.3f}, midnight {area(data['clt_sun'][:, :1]):.3f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
