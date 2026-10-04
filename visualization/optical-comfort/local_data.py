"""Plot native-sky, dusk-cloud and resolved water-slope products."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / 'research/runs/optical_comfort/figures'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    files = [ROOT / 'research/studies/optical_comfort/results/native_sky.json',
             ROOT / 'research/studies/cloud_twilight/results/local_clouds.json',
             ROOT / 'research/studies/optical_comfort/results/water_slopes.json']
    sky, cloud, water = [json.loads(p.read_text()) for p in files]
    expected = ('terluna.research.optical-comfort-native-sky/1', 'terluna.research.cloud-twilight-local/1',
                'terluna.research.optical-water-slopes/1')
    if any(p['schema'] != e for p, e in zip((sky, cloud, water), expected)):
        raise ValueError('Unsupported local optical products')
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    blue, teal, amber = '#3978ae', '#26968e', '#c98c38'
    a = axes[0, 0]
    values = np.array([sky['checks'][key] for key in ('maximum_ambient_relative_change', 'maximum_scene_total_relative_change')]) * 100
    a.barh(['Gaze ambient light', 'Finite-scene light'], values, color=blue)
    for i, v in enumerate(values):
        a.text(v + .004, i, f'{v:.3f}%', va='center')
    a.set(xlim=(0, .28), xlabel='Largest absolute change (%)', title='Native sky compared with packed sky')
    a = axes[0, 1]
    bands = cloud['height_bands']
    labels = [f"{b['height_km'][0]:g}–{b['height_km'][1]:g}" for b in bands]
    values = np.array([b['cloud_column_fraction'] for b in bands]) * 100
    a.bar(labels, values, color=teal)
    for i, v in enumerate(values):
        a.text(i, v + .3, f'{v:.2f}%', ha='center')
    a.set(ylim=(0, 16), xlabel='Altitude band (km)', ylabel='Share of all dusk columns (%)', title='Cloud occurrence in the equatorial ring')
    sites = list(water['summaries'])
    labels = ['West face', 'East face', 'South shore', 'Offshore']
    a = axes[1, 0]
    q = [water['summaries'][k]['rms_gradient'] for k in sites]
    mid = np.array([v['p50'] for v in q]); lo = np.array([v['p10'] for v in q]); hi = np.array([v['p90'] for v in q])
    a.errorbar(np.arange(4), mid, yerr=[mid - lo, hi - mid], fmt='o', color=blue, capsize=5)
    a.set(xticks=np.arange(4), xticklabels=labels, ylim=(0, .16), ylabel='RMS surface gradient',
          title='Resolved slopes: hourly 10th–90th percentiles')
    a = axes[1, 1]
    for key, label, offset, color in [('elevation_variance_fraction_above_025hz', 'Surface-elevation variance', -.18, blue),
                                      ('slope_fraction_above_025hz', 'Slope variance', .18, amber)]:
        values = [water['summaries'][k][key]['p50'] * 100 for k in sites]
        a.bar(np.arange(4) + offset, values, width=.35, color=color, label=label)
    a.set(xticks=np.arange(4), xticklabels=labels, ylim=(0, 100), ylabel='Median share above 0.25 Hz (%)',
          title='Shorter resolved waves carry substantial slope')
    a.legend(frameon=False, fontsize=9)
    fig.suptitle('Optical-comfort study · local-data follow-up', fontsize=16)
    fig.text(.5, .025, f"Clouds: {cloud['source']['snapshots']} second-cycle snapshots, {cloud['source']['columns']:,} dusk columns; height bands can overlap.\n"
             'Water: second-cycle coastal spectra, 0.005–0.497 Hz; shorter-wave roughness remains open.\n'
             'Sky: existing atmospheric-shape transfer. These are physical light and surface diagnostics.',
             ha='center', fontsize=9, color='#444444')
    fig.tight_layout(rect=(0, .105, 1, .95))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for suffix in ('png', 'svg'):
        fig.savefig(OUTPUT / f'local_data.{suffix}', dpi=160)
    (OUTPUT / 'manifest.json').write_text(json.dumps(dict(
        schema='terluna.visualization.optical-local-data/1',
        producer=dict(file=str(Path(__file__).relative_to(ROOT)), sha256=sha(Path(__file__))),
        inputs={str(p.relative_to(ROOT)): sha(p) for p in files},
        reading_rule='All bars, medians and ranges read the stored study products; no scene renderer or exposure map.',
        outputs={f'local_data.{s}': sha(OUTPUT / f'local_data.{s}') for s in ('png', 'svg')}), indent=2) + '\n')
    print(OUTPUT / 'local_data.png')


if __name__ == '__main__':
    main()
