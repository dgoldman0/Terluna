"""The seas' weather at 301 K: the lake kinds study's resort setting with the GCM's own water cells over the seas.

The resort at 301 K is research/studies/lake_kinds/setting.py's Setting301: lake_oxygen's setting with the share of
the way and the lapse rate moved to 301 K on climate_targets' ground-temperature basis. SeaSetting keeps all of it
(the air, gravity, the held seas' salinity, the floor's heat flow, the share and the lapse rate) and replaces only
the weather it reads, lake_oxygen's lake cells, with this study's extract of the seas' water cells
(research/runs/sea_oxygen/forcing.npz), from which it recomputes the Sun's phase and the wind's eddy factor exactly as
lake_oxygen's Setting computes them. lake_oxygen's windows(), resort_windows(), hourly() and forcing_for() then run
unchanged on these cells: a sea cell sits at the seas' level (elevation 0 above the common level, the lake level 0),
so its air, sky and pressure are the cell's own; a basin with no water cell would take its largest cell, carried
down to the seas' level along the target's lapse rate as lake_oxygen carries its lakes.

Each run basin takes two cells (RULE): the typical cell, whose mean surface temperature in the 301 K weather lies
nearest the area-weighted median of the basin's cells, and the coldest cell, whose coldest 3-day mean is the lowest
among the cells holding at least COLD_MIN_KM2 of the basin. A basin whose two cells' coldest 3-day means lie within
SAME_K of each other runs one column.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R
from shared.constants import SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'sea_oxygen'
FORCING = RUNS / 'forcing.npz'
FORCING_META = RUNS / 'forcing.json'
SEAS = HERE / 'results' / 'seas.json'
COLD_MIN_KM2 = 1000.0        # the coldest cell holds at least this much of the basin's water
SAME_K = 0.5                 # typical and coldest cells this close in their coldest 3-day means run one column
RULE = ('typical: the water cell whose mean 301 K surface temperature lies nearest the area-weighted median of the '
        'basin\'s cells; coldest: the cell, holding at least 1,000 km2 of the basin, whose coldest 3-day mean is the '
        'lowest; one column where the two coldest 3-day means lie within 0.5 K')


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def weighted_quantile(x, w, q):
    i = np.argsort(x, kind='stable')
    c = np.cumsum(w[i])
    return float(np.interp(q * c[-1], c - 0.5 * w[i], x[i]))


class SeaSetting(S.Setting301):
    """Setting301 reading the seas' water cells."""

    def __init__(self, ground_k=S.TARGET_K):
        super().__init__(ground_k=ground_k)
        self.lake_oxygen_inputs = dict(self.inputs)
        self.seas = json.loads(SEAS.read_text())
        meta = json.loads(FORCING_META.read_text())
        if meta['npz_sha256'] != sha256(FORCING):
            raise SystemExit('forcing.npz does not match forcing.json; run the extract again')
        if meta['producer']['read']['research/studies/sea_oxygen/results/seas.json'] != sha256(SEAS):
            raise SystemExit('the forcing was extracted for another seas.json; run the extract again')
        with np.load(FORCING) as z:
            self.f = {k: np.asarray(z[k], dtype=np.float64) for k in z.files}
        self.meta = meta
        self.phase = {run: self._phase(run) for run in ('lock', 'hot', 'eddies')}
        self.eddy_factor = self._eddy_factor()
        self.inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (SEAS, FORCING, FORCING_META)}
        for p in (R.TARGETS, R.SEA_CHEMISTRY, R.UPTAKE):
            self.inputs[str(p.relative_to(ROOT))] = self.lake_oxygen_inputs[str(p.relative_to(ROOT))]
        self.basins = {b['key']: b for b in self.seas['run_basins']}

    def site(self, i, basin_key=None):
        """lake_oxygen's site for cell i: at the seas' level, with the run basin's own hypsography as its lake."""
        s = self.meta['sites'][i]
        b = self.basins[basin_key] if basin_key else None
        lake = dict(level_above_sea_m=0.0)
        if b:
            lake.update(area_km2=b['area_km2'], hypsography=b['hypsography'], key=b['key'])
        return dict(index=i, band=basin_key or f'cell {s["row"]},{s["column"]}', lat=s['lat'], lon=s['lon'],
                    elevation=s['elevation_m'] if s['land_cell'] else 0.0, lake=lake,
                    eddy_factor=float(self.eddy_factor[i]), row=s['row'], column=s['column'],
                    land_cell=s['land_cell'])


def solar_bins(hour):
    return np.clip(((hour + 180.0) / 360.0 * R.SOLAR_HOUR_BINS).astype(int), 0, R.SOLAR_HOUR_BINS - 1)


def resort_surface_c(st: SeaSetting, i):
    """The GCM's surface temperature at cell i in the resort's weather (C): the lock's 3-day means moved by the
    share of the way of the two runs' mean cycles by solar hour, as resort_windows() moves the forcing."""
    site = st.site(i)
    lock = R.windows(st, 'lock', site, target='resort')
    hot = R.windows(st, 'hot', site, target='resort')

    def composite(w):
        b = solar_bins(w['hour'])
        return np.array([w['gcm_ts'][b == k].mean() for k in range(R.SOLAR_HOUR_BINS)])
    shifted = lock['gcm_ts'] + st.share['resort'] * (composite(hot) - composite(lock))[solar_bins(lock['hour'])]
    return shifted - 273.15


def cells_of(st: SeaSetting):
    """Each run basin's cells with the basin's water in each (km2) and their 301 K surface temperatures."""
    out = {k: [] for k in st.basins}
    cache = {}
    for i, s in enumerate(st.meta['sites']):
        for key, km2 in s['basins'].items():
            if i not in cache:
                t = resort_surface_c(st, i)
                cache[i] = dict(mean=float(t.mean()), min=float(t.min()), p01=float(np.quantile(t, 0.01)),
                                lock_min=float(st.f['lock__ts'][:, i].min() - 273.15))
            out[key].append(dict(index=i, row=s['row'], column=s['column'], lat=s['lat'], lon=s['lon'],
                                 area_km2=km2, land_cell=s['land_cell'], **cache[i]))
    return out


def choose_cells(st: SeaSetting, cells=None):
    """The typical and the coldest cell of each run basin (RULE), with the figures that chose them."""
    cells = cells or cells_of(st)
    out = {}
    for key, cs in cells.items():
        means = np.array([c['mean'] for c in cs])
        w = np.array([c['area_km2'] for c in cs])
        median = weighted_quantile(means, w, 0.5)
        typical = cs[int(np.argmin(np.abs(means - median)))]
        pool = [c for c in cs if c['area_km2'] >= COLD_MIN_KM2] or cs
        coldest = min(pool, key=lambda c: c['min'])
        same = coldest['index'] == typical['index'] or abs(coldest['min'] - typical['min']) < SAME_K
        total = w.sum()
        within = [c for c in cs if c['min'] <= coldest['min'] + SAME_K]
        out[key] = dict(typical=typical, coldest=coldest, one_column=bool(same), cells=len(cs),
                        area_km2_in_cells=round(float(total), 1), median_mean_c=round(median, 3),
                        mean_c_range=[round(float(means.min()), 3), round(float(means.max()), 3)],
                        coldest_spell_c_range=[round(min(c['min'] for c in cs), 3), round(max(c['min'] for c in cs), 3)],
                        area_share_within_same_k_of_coldest=round(sum(c['area_km2'] for c in within) / total, 4))
    return out


def describe(st: SeaSetting):
    d = S.describe(st)
    d.update(sites=len(st.meta['sites']), synodic_month_days=SYNODIC_MONTH_DAYS,
             eddy_wind_factor_range=[round(float(st.eddy_factor.min()), 3), round(float(st.eddy_factor.max()), 3)])
    return d
