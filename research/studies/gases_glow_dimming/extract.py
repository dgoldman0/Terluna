"""What the study reads from the GCM's own output beyond the climatologies: the air's temperature, humidity and cloud
by layer and the rain's convective part over the compared years, and the first 3-day means of each branch, before its
climate responds.

    climate/gcm/.venv/bin/python -m research.studies.gases_glow_dimming.extract     # seconds; needs netCDF4

It reads the ExoPlaSim output of the runs in climate/gcm/runs (ignored by Git) and writes
research/runs/gases_glow_dimming/extract.json (ignored), with the sha256 of every file it read; run.py binds it.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'research' / 'runs' / 'gases_glow_dimming' / 'extract.json'
SCHEMA = 'terluna.research.gases-glow-dimming-extract/1'
RUNS_DIR = ROOT / 'climate' / 'gcm' / 'runs'
# The compared years of each run, as run.py compares them.
LAYERS = {'A28_dim5_moon': (20, 29), 'A28_gases': (10, 19), 'A28_gases_glow': (10, 19),
          'A28_gases_glow_held': (10, 19)}
# Branches from the design case's year-29 restart whose first means are read: the calibration pair, the control
# year flown on the patched build, and the three runs.
FIRST = ('A28_cal_control', 'A28_cal_co2_1022', 'A28_cal_control_patched', 'A28_gases', 'A28_gases_glow',
         'A28_gases_glow_held')
FIRST_KEYS = ('ntr', 'rst', 'rlut', 'rss', 'rls', 'nbr', 'hfls', 'hfss', 'ts', 'pr')
FIRST_COUNTS = (1, 2, 4, 10)


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def output(run, year):
    return RUNS_DIR / run / 'model' / f'MOST.{year:05d}.nc'


def weights(lat, shape):
    return np.cos(np.radians(lat))[:, None] * np.ones(shape)


def layers(run, years):
    """Area-weighted means by layer of the air temperature (K), specific humidity (kg/kg) and cloud cover, and the
    global means of the rain and its convective part (m/s)."""
    import netCDF4
    acc, rain, files = None, np.zeros(2), {}
    for year in range(years[0], years[1] + 1):
        path = output(run, year)
        with netCDF4.Dataset(path) as d:
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            fields = [np.asarray(d[k][:], float).mean(axis=0) for k in ('ta', 'hus', 'cl')]   # (layer, lat, lon)
            surface = [np.asarray(d[k][:], float).mean(axis=0) for k in ('pr', 'prc')]
        w = weights(lat, fields[0].shape[-1])
        row = np.array([(f * w).sum(axis=(-2, -1)) / w.sum() for f in fields])
        acc = row if acc is None else acc + row
        rain += [(f * w).sum() / w.sum() for f in surface]
        files[str(path.relative_to(ROOT))] = sha256(path)
    n = years[1] - years[0] + 1
    acc /= n
    rain /= n
    return dict(years=list(years), sigma=sigma.tolist(), ta_k=acc[0].tolist(), hus=acc[1].tolist(),
                cl=acc[2].tolist(), rain_m_s=float(rain[0]), convective_rain_m_s=float(rain[1])), files


def first_means(run):
    """Global means of the first 1, 2, 4 and 10 output means of the branch's year 0. Fluxes W/m2 as PlaSim writes
    them (ntr, rst, rlut at the top of the air; rss, rls, nbr, hfls, hfss at the surface, positive downward),
    ts K, pr m/s."""
    import netCDF4
    path = output(run, 0)
    with netCDF4.Dataset(path) as d:
        lat = np.asarray(d['lat'][:], float)
        v = {k: np.asarray(d[k][:], float) for k in FIRST_KEYS}
    w = weights(lat, v['ts'].shape[-1])
    series = {k: (x * w).sum(axis=(-2, -1)) / w.sum() for k, x in v.items()}
    return {str(n): {k: float(np.mean(s[:n])) for k, s in series.items()} for n in FIRST_COUNTS}, \
        {str(path.relative_to(ROOT)): sha256(path)}


def extract():
    files, out = {}, dict(schema=SCHEMA, layers={}, first_means={})
    for run, years in LAYERS.items():
        out['layers'][run], f = layers(run, years)
        files.update(f)
    for run in FIRST:
        out['first_means'][run], f = first_means(run)
        files.update(f)
    out['producer'] = dict(files={str(Path(__file__).resolve().relative_to(ROOT)): sha256(__file__)}, read=files)
    return out


def main(argv=None) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(extract(), indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
