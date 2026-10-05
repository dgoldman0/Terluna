"""What the project's own models give the aerosol study: the regions' areas and their weather by day and night.

The regions follow the climate runs and the ecology register's landscapes (research/studies/lunar_cycle_ecology,
F1-F3): the seas; the rain-fed lakes; the wet land, where the GCM rains 0.5 mm/day or more, mostly within 45 degrees
of the equator, for the forests of deep canopies and their wet margins; the fog desert, the land between 30 and 60
degrees with less rain, for the plains of storage-organ plants, grasses and crusts; and the polar dry land poleward of
60 degrees. Areas come from the GCM climatology (A28_dim5_moon, years 20-29) and geography's lakes. The weather near
the ground comes from CM1's rings, which resolve the air near the ground that the GCM's lowest layer, 2.5 km deep in
the lunar air, averages over: humidity, fog, rain, wind at 10 m and the mixed layer's depth, by day and by night.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CLIMATOLOGY = ROOT / 'climate' / 'gcm' / 'products' / 'climatology_A28_dim5_moon.npz'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
RINGS = ROOT / 'climate' / 'results' / 'crm'
DRY_MM_DAY = 0.5                                  # the GCM README's dry land
SUNSET_DEG = 95.0                                 # hour angle of sunset and dawn over the rings
# the rings whose land or water stands for each region's weather
REGIME_RINGS = {
    'seas': (('ring_a', 'water'), ('ring_a_prime', 'water'), ('ring_45n_lsw', 'water'), ('ring_45s_lsw', 'water')),
    'lakes': (('ring_a', 'water'), ('ring_a_prime', 'water')),
    'wet_land': (('ring_a', 'land'), ('ring_a_prime', 'land')),
    'fog_desert': (('ring_45n_lsw', 'land'), ('ring_45s_lsw', 'land'), ('ring_70_45e', 'land'), ('ring_70_135e', 'land')),
    'polar_dry_land': (('ring_80n_lsw', 'land'), ('ring_80s_lsw', 'land')),
}


def areas() -> dict:
    """Each region's share of the Moon's surface: the seas from the atlas and the lakes from the drainage, the land
    between them split by the GCM's rain and latitude (the GCM's T21 coastline, which held the earlier, larger lakes as
    water, scaled to the land left by the atlas and the drainage)."""
    d = np.load(CLIMATOLOGY)
    lat, lsm, pr = d['lat'], d['lsm'], d['pr_mm_day']
    w = np.cos(np.radians(lat))[:, None] * np.ones_like(lsm)
    land = lsm > 0.5
    alat = np.abs(lat)[:, None] * np.ones_like(lsm)
    seas = json.loads(ATLAS.read_text())['water_share']
    lakes = json.loads(DRAINAGE.read_text())['lakes']['area_share']
    scale = (1.0 - seas - lakes) / float((w * land).sum() / w.sum())
    share = lambda mask: scale * float((w * mask).sum() / w.sum())
    return dict(seas=seas, lakes=lakes, wet_land=share(land & (pr >= DRY_MM_DAY)),
                fog_desert=share(land & (pr < DRY_MM_DAY) & (alat < 60.0)),
                polar_dry_land=share(land & (pr < DRY_MM_DAY) & (alat >= 60.0)))


def ring(case: str) -> dict:
    return json.loads((RINGS / f'ring_{case}.json').read_text())


def regime(region: str) -> dict:
    """A region's weather near the ground by day and by night (hour angle within or beyond 95 degrees of noon),
    averaged over its rings: relative humidity and temperature at 2 m, fog, rain, wind at 10 m and the mixed layer's
    depth; the wind's 99th percentile and largest value over the lunar day, and the storms' gusts."""
    rows = []
    for case, surface in REGIME_RINGS[region]:
        d = ring(case)
        hour = np.array(d['by_hour_angle']['hour_angle_deg'])
        night = np.abs(hour) > SUNSET_DEG
        by = {k: np.array(v, float) for k, v in d['by_hour_angle'][surface].items() if len(v) == hour.size}
        row = {}
        for phase, sel in (('day', ~night), ('night', night)):
            row[phase] = dict(relative_humidity=np.nanmean(by['humidity_2m'][sel]), air_c=np.nanmean(by['air_2m_c'][sel]),
                              fog=np.nanmean(by['fog'][sel]), rain_mm_day=24.0 * np.nanmean(by['rain_mm_h'][sel]),
                              wind_m_s=np.nanmean(by['wind_10m_m_s'][sel]),
                              mixed_layer_km=np.nanmean(by['mixed_layer_km'][sel]))
        row['wind_p99_m_s'] = d['wind_10m_m_s'][surface]['p99']
        row['wind_max_m_s'] = d['wind_10m_m_s'][surface]['max']
        row['day_hours'] = float(24.0 * 29.53 * (~night).mean())
        rows.append(row)
    out = {}
    for phase in ('day', 'night'):
        out[phase] = {k: float(np.mean([r[phase][k] for r in rows])) for k in rows[0][phase]}
    for k in ('wind_p99_m_s', 'wind_max_m_s', 'day_hours'):
        out[k] = float(np.mean([r[k] for r in rows]))
    out['rings'] = [f'{c} ({s})' for c, s in REGIME_RINGS[region]]
    return out
