"""The gravity pair: the same air over the same sea at lunar and at Earth gravity, level by level.

    climate/gcm/.venv/bin/python -m climate.crm.pair_analysis     # ../results/crm/gravity_pair.json, products/gravity_pair.npz

The Earth case's lengths and times are the Moon's times g_moon/g_earth, so its levels sit at the same
pressures as the Moon case's and the two can be set side by side level by level. If lunar convection
were Earth's stretched in size and duration, every profile here would match; the differences measure
what does not stretch: the fall speeds of drops and ice and the rates at which rain forms.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

from climate.crm.ring_analysis import read_snapshot, RUNS, PRODUCTS, RESULTS, CLOUD_KG_KG

SCHEMA = 'terluna.climate.crm-gravity-pair/1'
EVIDENCE = ('Two CM1 r22.0 runs in two dimensions over a sea held at 298.9 K, without radiation, starting from the '
            'equatorial ring\'s own air over its seas, with the domain-mean air above 8-16 km (lunar) relaxed over '
            '3 days to the GCM design case\'s equatorial profile (A28_dim5): one at lunar gravity (6-km columns, '
            '150-km top, 20 days), one at Earth gravity with every '
            'length and time scaled by g_moon/g_earth (0.99-km columns, 24.8-km top, 3.3 days). The Morrison '
            'microphysics falls at each gravity\'s speeds. The pair isolates what gravity changes in the clouds '
            'at a fixed air mass; radiation, land and the day-night cycle are left out.')
READING_RULE = ('Levels are matched by index: the same fraction of the surface pressure in both runs. Heights are '
                'given in lunar kilometres (the Earth case\'s times 6.04). Rates are per Earth day. Means are over '
                'the second half of each run. If lunar convection were Earth\'s stretched, cover, rates, updraft '
                'speeds and pressure-weighted mixing ratios would match (ratio 1) and paths per square metre would '
                'differ by g_earth/g_moon; departures come from what does not stretch.')


def rho_surface(d):
    """Air density (kg/m3) at the lowest model level, from its pressure and potential temperature."""
    p, th = d['prs'][0], d['th'][0]
    return p / (287.04 * th * (p / 1e5) ** (287.04 / 1005.7))


def run_means(case: Path, first_fraction=0.5):
    record = json.loads((case / 'case.json').read_text())
    scale = record['time_scale']
    tap = record['configuration']['output_s'] * scale
    outputs = sorted(int(p.name[8:14]) for p in case.glob('cm1out_t*_s.dat'))
    last = outputs[-1]
    use = [n for n in outputs if n - 1 >= first_fraction * (last - 1)]
    acc, count = {}, 0
    updrafts = []
    for n in use:
        d = read_snapshot(case, n)
        cloudy = (d['qc'] + d['qi']) >= CLOUD_KG_KG
        fields = dict(cloud_fraction=cloudy.mean(axis=1), qc=d['qc'].mean(axis=1), qr=d['qr'].mean(axis=1),
                      qi=d['qi'].mean(axis=1), qs=d['qs'].mean(axis=1), qg=d['qg'].mean(axis=1),
                      prs=d['prs'].mean(axis=1), rh_proxy=d['qv'].mean(axis=1),
                      cover=np.array([cloudy.any(axis=0).mean()]), rain_mm_day=np.array([d['prate'].mean() * 86400.0]),
                      evap_mm_day=np.array([(d['qvflux'] * rho_surface(d)).mean() * 86400.0]))
        up = d['winterp'][cloudy & (d['winterp'] > 1.0)]
        updrafts.append(up)
        for k, v in fields.items():
            acc[k] = acc.get(k, 0) + v
        count += 1
    means = {k: v / count for k, v in acc.items()}
    ups = np.concatenate(updrafts) if updrafts else np.array([0.0])
    means['updraft_p50_m_s'] = float(np.percentile(ups, 50)) if ups.size else 0.0
    means['updraft_p99_m_s'] = float(np.percentile(ups, 99)) if ups.size else 0.0
    zw = np.loadtxt(case / 'input_grid_z')
    means['z_lunar_km'] = 0.5 * (zw[1:] + zw[:-1]) / scale / 1000.0
    means['snapshots'] = len(use)
    means['span_model_days'] = [(use[0] - 1) * tap / 86400.0, (use[-1] - 1) * tap / 86400.0]
    return means


def main(argv=None) -> int:
    moon, earth = run_means(RUNS / 'pair_moon'), run_means(RUNS / 'pair_earth')
    path = lambda m, k, scale: float((m[k] * np.gradient(m['prs']) * -1 / scale).sum())     # kg per m2 of condensate
    g_moon = json.loads((RUNS / 'pair_moon' / 'case.json').read_text())['build']['gravity_m_s2']
    g_earth = json.loads((RUNS / 'pair_earth' / 'case.json').read_text())['build']['gravity_m_s2']
    weighted = lambda m, k: float((m[k] * np.gradient(m['prs']) * -1).sum() / (m['prs'][0] - m['prs'][-1]) * 1000.0)
    z = moon['z_lunar_km']
    cloudy_top = lambda m: float(z[np.where(m['cloud_fraction'] > 0.01)[0].max()]) if (m['cloud_fraction'] > 0.01).any() else None
    table = {}
    for name, m, g in (('moon', moon, g_moon), ('earth_g', earth, g_earth)):
        table[name] = dict(cloud_cover=float(m['cover'][0]), rain_mm_day=float(m['rain_mm_day'][0]),
                           evaporation_mm_day=float(m['evap_mm_day'][0]),
                           cloud_liquid_path_kg_m2=path(m, 'qc', g), ice_path_kg_m2=path(m, 'qi', g),
                           snow_path_kg_m2=path(m, 'qs', g), graupel_path_kg_m2=path(m, 'qg', g),
                           rain_path_kg_m2=path(m, 'qr', g),
                           mean_mixing_ratio_g_kg={k: weighted(m, k) for k in ('qc', 'qr', 'qi', 'qs', 'qg')},
                           updraft_p50_m_s=m['updraft_p50_m_s'], updraft_p99_m_s=m['updraft_p99_m_s'],
                           max_cloud_fraction=float(m['cloud_fraction'].max()),
                           most_cloud_at_lunar_km=float(z[int(np.argmax(m['cloud_fraction']))]),
                           cloud_top_lunar_km=cloudy_top(m), snapshots=m['snapshots'], span_model_days=m['span_model_days'])
    ratio = lambda a, b: a / b if b else None
    ratios = {k: ratio(table['moon'][k], table['earth_g'][k]) for k in
              ('cloud_cover', 'rain_mm_day', 'evaporation_mm_day', 'updraft_p50_m_s', 'updraft_p99_m_s', 'max_cloud_fraction')}
    ratios['mean_mixing_ratio'] = {k: ratio(table['moon']['mean_mixing_ratio_g_kg'][k], table['earth_g']['mean_mixing_ratio_g_kg'][k])
                                   for k in ('qc', 'qr', 'qi', 'qs', 'qg')}
    table['moon_over_earth_g'] = ratios
    summary = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, runs=table,
                   producer=dict(domain='climate', files={'crm/pair_analysis.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}))
    PRODUCTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'gravity_pair.json').write_text(json.dumps(summary, indent=1) + '\n')
    np.savez_compressed(PRODUCTS / 'gravity_pair.npz', metadata=json.dumps({k: summary[k] for k in ('schema', 'evidence', 'reading_rule')}),
                        z_lunar_km=moon['z_lunar_km'],
                        **{f'{k}_moon': moon[k] for k in ('cloud_fraction', 'qc', 'qr', 'qi', 'qs', 'qg', 'prs')},
                        **{f'{k}_earth_g': earth[k] for k in ('cloud_fraction', 'qc', 'qr', 'qi', 'qs', 'qg', 'prs')})
    print(json.dumps(table, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
