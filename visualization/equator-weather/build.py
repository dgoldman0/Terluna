"""Build the equator weather page from the climate domain's cloud-resolving results.

    python3 visualization/equator-weather/build.py          # writes Open_Moon_Equator_Weather.html here

It reads two products that climate/crm/ring_analysis.py writes and Git keeps:
climate/results/crm/ring_ring.json (terluna.climate.crm-ring/1: composites by local time, storm, wind
and water statistics, the GCM's equatorial composites) and ring_ring_fields.npz
(terluna.climate.crm-ring-fields/1: the circulation by local time and height, the heaviest storm's
section and the rain map). The charts draw those arrays, and every number in the page's text is
filled from them, so a rerun of the model changes the page with it. The page carries the products'
evidence statements, reading rules and sha256.
"""
from __future__ import annotations
import base64
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / 'climate' / 'results' / 'crm'
SUMMARY = RESULTS / 'ring_ring.json'
FIELDS = RESULTS / 'ring_ring_fields.npz'
TEMPLATE = HERE / 'equator_weather.template.html'
OUTPUT = HERE / 'Open_Moon_Equator_Weather.html'
SCHEMAS = ('terluna.climate.crm-ring/1', 'terluna.climate.crm-ring-fields/1')
RAIN_FLOOR, RAIN_TOP = 0.1, 50.0          # mm/h: the rain map's colour scale, coded into one byte per cell
DOCUMENT = '<!doctype html>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n'


def load():
    summary = json.loads(SUMMARY.read_text())
    fields = dict(np.load(FIELDS))
    meta = json.loads(str(fields.pop('metadata')))
    if summary.get('schema') != SCHEMAS[0] or meta.get('schema') != SCHEMAS[1]:
        raise ValueError(f'expected {SCHEMAS}, found {summary.get("schema")} and {meta.get("schema")}')
    return summary, fields, meta


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def code_rain(rain_mm_h: np.ndarray) -> np.ndarray:
    """One byte per cell: 0 below RAIN_FLOOR, then 1-255 logarithmic up to RAIN_TOP and above."""
    lo, hi = np.log10(RAIN_FLOOR), np.log10(RAIN_TOP)
    scaled = 1 + 254 * (np.log10(np.maximum(rain_mm_h, RAIN_FLOOR)) - lo) / (hi - lo)
    return np.where(rain_mm_h < RAIN_FLOOR, 0, np.clip(np.round(scaled), 1, 255)).astype(np.uint8)


def decode_rain(code: np.ndarray) -> np.ndarray:
    """The lower edge of each byte's rain class (mm/h), zero for 0: the inverse the page uses."""
    lo, hi = np.log10(RAIN_FLOOR), np.log10(RAIN_TOP)
    return np.where(code == 0, 0.0, 10 ** (lo + (code.astype(float) - 1) / 254 * (hi - lo)))


def page_data(summary: dict, fields: dict, lunar_day_days: float) -> dict:
    r = lambda a, d: [round(float(v), d) for v in np.asarray(a).ravel()]
    rows = lambda a, d: [r(row, d) for row in np.asarray(a, dtype=float)]
    b = summary['by_hour_angle']
    keys = ('air_2m_c', 'ground_c', 'humidity_2m', 'wet_bulb_2m_c', 'rain_mm_h', 'cloud_cover', 'wind_10m_m_s',
            'mixed_layer_km', 'cloud_top_km')
    rain = fields['rain_map_max_mm_h'].astype(float)
    return dict(
        lunar_day_days=lunar_day_days,
        by_hour=dict(h=b['hour_angle_deg'], **{s: {k: b[s][k] for k in keys} for s in ('land', 'water')}),
        sections=dict(h=r(fields['section_hour_angle_deg'], 1), z=r(fields['section_height_km'], 2),
                      u=rows(fields['section_eastward_m_s'], 2), w=rows(fields['section_upward_m_s'] * 100.0, 2),
                      cloud=rows(fields['section_cloud_fraction'], 3)),
        storm=dict(day=round(float(fields['storm_day']), 2), hour_angle=round(float(fields['storm_hour_angle_deg']), 1),
                   rain_max_mm_h=round(float(fields['storm_rain_max_mm_h']), 1),
                   dx_km=round(float(fields['storm_x_km'][1] - fields['storm_x_km'][0]), 3), z=r(fields['storm_height_km'], 2),
                   cond_log10=rows(np.log10(np.maximum(fields['storm_condensate_g_kg'].astype(float), 1e-4)), 2),
                   rainwater_g_kg=rows(fields['storm_rain_water_g_kg'], 2), updraft_m_s=rows(fields['storm_updraft_m_s'], 1),
                   rain_mm_h=r(fields['storm_rain_mm_h'], 1), wind10_m_s=r(fields['storm_wind_10m_m_s'], 1),
                   air2m_c=r(fields['storm_air_2m_c'], 1), land=[int(v) for v in fields['storm_land']]),
        rain_map=dict(time_day=r(fields['rain_map_time_day'], 3), x_km=r(fields['rain_map_x_km'], 1),
                      land=[int(v >= 0.5) for v in fields['rain_map_land_share']], nt=int(rain.shape[0]), nx=int(rain.shape[1]),
                      code_b64=base64.b64encode(code_rain(rain).tobytes()).decode(), floor=RAIN_FLOOR, top=RAIN_TOP,
                      noon_x_km=r(fields['rain_map_noon_x_km'], 1), length_km=round(float(fields['rain_map_x_km'][-1] + fields['rain_map_x_km'][0]), 1)),
        **{k: summary[k] for k in ('wind_aloft_m_s', 'wind_10m_m_s', 'storms', 'storm_tracks', 'rain_spells', 'water_budget')})


def crossing(h, u, rising: bool) -> float:
    """Local time (degrees) where u changes sign, going up (-to+) or down (+to-), interpolated, on the ring."""
    h, u = np.asarray(h, float), np.asarray(u, float)
    for i in range(len(h)):
        a, b = u[i], u[(i + 1) % len(h)]
        if (a < 0 <= b) if rising else (a > 0 >= b):
            h0, h1 = h[i], h[(i + 1) % len(h)] + (360.0 if i == len(h) - 1 else 0.0)
            x = h0 + (h1 - h0) * (-a) / (b - a)
            return (x + 180.0) % 360.0 - 180.0
    return float('nan')


def text_values(summary: dict, lunar_day_days: float) -> dict:
    """Every number the page's text quotes, from the products."""
    b, st, tr, wb = summary['by_hour_angle'], summary['storms'], summary['storm_tracks'], summary['water_budget']
    h = b['hour_angle_deg']
    at = lambda surface, key, angle: b[surface][key][h.index(angle)]
    days = lambda deg: deg / 360.0 * lunar_day_days
    wet = [a for a, v in zip(h, b['land']['rain_mm_h']) if v > 0.05]
    c = summary['circulation']
    u1 = c['eastward_m_s'][c['heights_km'].index(1)]
    w5 = np.array(c['upward_cm_s'][c['heights_km'].index(5)])
    quiet = np.array([not (0.0 <= a <= 90.0) for a in c['hour_angle_deg']])
    sink_cm_s = -float(w5[quiet].mean())
    parts, meets = crossing(c['hour_angle_deg'], u1, True), crossing(c['hour_angle_deg'], u1, False)
    flight = [v for v in b['land']['cloud_in_flight_band'] if v > 0]
    wind40 = summary['wind_aloft_m_s']['40_km']
    gcm_land = summary['gcm']['land']
    f = lambda v, d=0: f'{v:,.{d}f}'
    return dict(
        lunar_day=f(lunar_day_days, 1), week=f(lunar_day_days / 4, 1),
        water_mm=f(np.mean([wb['precipitable_water_mm']['start'], wb['precipitable_water_mm']['end']]), 0),
        residence_days=f(wb['residence_days'], 0),
        top_max_km=f(st['cloud_top_km']['max'], 0), top_median_km=f(st['cloud_top_km']['median'], 0),
        top_p90_km=f(st['cloud_top_km']['p90'], 0), rain_max_mm_h=f(st['rain_max_mm_h']['max'], 0),
        sunset_air=f(at('land', 'air_2m_c', 95), 0), dawn_air=f(at('land', 'air_2m_c', -95), 1),
        fog_pct=f(100 * at('land', 'cloud_cover', -145), 0), fog_top_m=f(round(1000 * at('land', 'cloud_top_km', -145), -2), 0),
        noon_air=f(at('land', 'air_2m_c', 5), 0), noon_ground=f(at('land', 'ground_c', 5), 0),
        noon_rh=f(100 * at('land', 'humidity_2m', 5), 0), noon_wet_bulb=f(at('land', 'wet_bulb_2m_c', 5), 0),
        noon_tower_km=f(at('land', 'cloud_top_km', 5), 0), late_air=f(at('land', 'air_2m_c', 75), 0),
        storms_from_days=f(days(min(wet) - 5), 0), storms_to_days=f(days(max(wet) + 5), 0),
        sea_surface=f(np.mean(b['water']['ground_c']), 1), sea_air_min=f(min(b['water']['air_2m_c']), 0),
        sea_air_max=f(max(b['water']['air_2m_c']), 0),
        storms_at_once=f(st['per_snapshot'], 1), storm_land_pct=f(100 * st['land_share'], 0),
        width_median=f(st['width_km']['median'], 0), width_p90=f(st['width_km']['p90'], 0), width_max=f(st['width_km']['max'], 0),
        life_median_h=f(tr['lifetime_h']['median'], 0), life_max_h=f(tr['lifetime_h']['max'], 0),
        parts_before_dawn=f(-90 - parts, 0), meets_after_noon=f(meets, 0),
        morning_flow=f(max(u1), 1), evening_flow=f(-min(u1), 1),
        sink_cm_s=f(sink_cm_s, 1), sink_km_day=f(sink_cm_s / 100 * 86400 / 1000, 1),
        wind10_land=f(summary['wind_10m_m_s']['land']['median'], 1), wind10_sea=f(summary['wind_10m_m_s']['water']['median'], 1),
        gust_max=f(max(summary['wind_10m_m_s']['land']['max'], summary['wind_10m_m_s']['water']['max']), 0),
        band_median=f(wind40['median'], 0), band_p99=f(wind40['p99'], 0),
        band_w_p99=f(summary['vertical_wind_m_s']['40_km']['abs_p99'], 1),
        band_cloud_lo=f(100 * min(flight), 0) if flight else '0', band_cloud_hi=f(100 * max(flight), 0) if flight else '0',
        band_updraft=f(summary['vertical_wind_m_s']['40_km']['up_max'], 0),
        evap_all=f(wb['evaporation_mm_day']['all'], 1), evap_land=f(wb['evaporation_mm_day']['land'], 1),
        evap_sea=f(wb['evaporation_mm_day']['water'], 1), rain_all=f(wb['rain_mm_day']['all'], 1),
        spells=f(summary['rain_spells']['land']['spells_per_lunar_day'], 1),
        rain_per_lunar_day=f(wb['rain_mm_day']['land'] * lunar_day_days, 0), rain_per_year=f(round(wb['rain_mm_day']['land'] * 365.25, -1), 0),
        gcm_rain_land=f(np.mean(gcm_land['rain_mm_h']) * 24, 1), gcm_cloud_land=f(100 * np.mean(gcm_land['cloud_cover']), 0),
        cloud_land=f(100 * summary['cloud_cover']['land'], 0),
        span_from=f(summary['span_days'][0], 1), span_to=f(summary['span_days'][1], 1), snapshots=f(summary['snapshots'], 0))


def fill(template: str, values: dict) -> str:
    missing = sorted(set(re.findall(r'\{\{(\w+)\}\}', template)) - set(values))
    if missing:
        raise KeyError(f'the template quotes values the products do not give: {missing}')
    return re.sub(r'\{\{(\w+)\}\}', lambda m: values[m.group(1)], template)


def build() -> str:
    summary, fields, meta = load()
    constants = json.loads((ROOT / 'shared' / 'constants.json').read_text())
    lunar_day_days = constants['moon']['synodic_month_days']
    data = page_data(summary, fields, lunar_day_days)
    values = text_values(summary, lunar_day_days)
    values['source'] = (f'Products: <code>{SUMMARY.relative_to(ROOT)}</code> ({summary["schema"]}, sha256 {sha256(SUMMARY)[:16]}) and '
                        f'<code>{FIELDS.relative_to(ROOT)}</code> ({meta["schema"]}, sha256 {sha256(FIELDS)[:16]}), '
                        f'producer <code>climate/crm/ring_analysis.py</code> {summary["producer"]["files"]["crm/ring_analysis.py"]}. '
                        f'Evidence: {summary["evidence"]} Reading rule: {summary["reading_rule"]} {meta["reading_rule"]}')
    page = fill(TEMPLATE.read_text(), values)
    payload = json.dumps(data, separators=(',', ':'))
    if '</' in payload:
        raise ValueError('the data would close its script element')
    return page.replace('__DATA__', payload)


def main(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--fragment', type=Path, help='also write the page without its document head, for hosts that add their own')
    args = parser.parse_args(argv)
    page = build()
    OUTPUT.write_text(DOCUMENT + page)
    print(f'wrote {OUTPUT.relative_to(ROOT)} ({OUTPUT.stat().st_size / 1e3:.0f} kB)')
    if args.fragment:
        args.fragment.write_text(page)
    return 0


if __name__ == '__main__':
    sys.exit(main())
