"""Plot the committed cloud-twilight product; no new light transport here."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / 'research/studies/cloud_twilight/results/cloud_twilight.json'
OUT = Path(__file__).resolve().parent / 'figures'
COLOURS = {20.: '#2878a0', 40.: '#16857d', 60.: '#bd8b2b', 80.: '#a24e69'}


def main():
    data = json.loads(PRODUCT.read_text())
    if data['schema'] != 'terluna.research.cloud-twilight/1':
        raise ValueError('Unsupported cloud-twilight product')
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(14.8, 12))
    a, b, c, d = axes.ravel()

    for height, colour in COLOURS.items():
        rows = [r for r in data['overhead'] if r['model'] == 'current_design' and r['height_km'] == height]
        a.semilogy([r['hours_after_surface_solar_centre_set'] for r in rows],
                   [r['cloud_direct_normal_lux'] or np.nan for r in rows],
                   color=colour, lw=2.3, label=f'{height:g} km')
    a.set(xlabel='Hours after geometric surface sunset', ylabel='Direct normal illumination at cloud (lux)',
          xlim=(0, 36), ylim=(10, 70000))
    a.legend(frameon=False, ncol=2, loc='upper right')
    a.text(.04, .08, 'Current column + design shield\nOverhead cloud; molecular extinction only',
           transform=a.transAxes, fontsize=9.5, color='#555555', linespacing=1.6)

    ground = data['ground_context']['rows']
    b.semilogy([r['hours'] for r in ground], [r['horizontal_total_lux'] for r in ground],
               color='#525f72', lw=2.5)
    for dep in (5., 10., 15., 20.):
        row = next(r for r in ground if r['depression_deg'] == dep)
        b.scatter(row['hours'], row['horizontal_total_lux'], color='#525f72', s=28)
        b.annotate(f"{row['horizontal_total_lux']/1000:.2g} klux", (row['hours'], row['horizontal_total_lux']),
                   xytext=(0, 11), textcoords='offset points', fontsize=10, ha='center')
    b.set(xlabel='Hours after geometric surface sunset', ylabel='Horizontal illumination at ground (lux)',
          xlim=(0, 50), ylim=(80, 10000))
    b.text(.04, .08, 'Inherited unfiltered sky proxy, albedo 0.1\nContext: different atmosphere from panel A',
           transform=b.transAxes, fontsize=9.5, color='#555555', linespacing=1.6)

    for dep, colour in [(5., '#16857d'), (10., '#bd8b2b'), (15., '#a24e69')]:
        rows = [r for r in data['views'] if r['model'] == 'current_design' and r['height_km'] == 40
                and r['observer_sun_depression_deg'] == dep and 0 <= r['westward_surface_range_km'] <= 300]
        hours = rows[0]['hours_after_surface_solar_centre_set']
        c.semilogy([r['westward_surface_range_km'] for r in rows],
                   [100 * r['two_path_source_fraction'] or np.nan for r in rows],
                   color=colour, lw=2.1, marker='o', ms=4, label=f'{hours:.1f} hours; Sun {dep:g}° down')
    c.set(xlabel='Surface distance west, toward sunset (km)',
          ylabel='Sun → cloud → eye source factor (%)', xlim=(0, 310), ylim=(.03, 8))
    c.legend(frameon=False, loc='lower right', fontsize=9.5)
    c.text(.03, .95, '40 km feature; current column + shield\nCloud scattering and path glow not calculated',
           transform=c.transAxes, va='top', fontsize=9.5, color='#555555', linespacing=1.6)

    section = data['occurrence']
    hour, height = np.asarray(section['hour_angle_deg']), np.asarray(section['height_km'])
    mask = (hour >= 65) & (hour <= 125)
    times = (hour[mask] - 90) * data['assumptions']['hours_per_degree']
    cloud = np.asarray(section['fraction'])[mask].T * 100
    colour = d.pcolormesh(times, height, np.maximum(cloud, .01), shading='nearest', cmap='YlGnBu',
                         norm=LogNorm(vmin=.1, vmax=10), rasterized=True)
    d.axvline(0, color='#cf695b', ls='--', lw=1.4)
    d.set(xlabel='Hours from local geometric sunset', ylabel='Height above flat ground (km)',
          ylim=(0, 85), xlim=(-59, 79))
    d.text(.04, .94, 'Corrected 2-D equatorial ring\nFraction of samples cloudy at each height',
           transform=d.transAxes, va='top', fontsize=9.5, color='#333333', linespacing=1.6,
           bbox=dict(facecolor='white', edgecolor='none', alpha=.75, pad=4))
    cb = fig.colorbar(colour, ax=d, pad=.025, fraction=.055, ticks=[.1, 1, 10], extend='both')
    cb.ax.set_yticklabels(['0.1%', '1%', '10%'])
    cb.set_label('Cloud occurrence per cell')

    titles = ['A   Elevated cloud material keeps a fading direct beam',
              'B   The twilight background stays substantial',
              'C   Illumination and viewing paths favour different ranges',
              'D   The highest clouds are sparse near dusk']
    for ax, title in zip(axes.ravel(), titles):
        ax.set_title(title, loc='left', fontsize=12, weight='bold', pad=14)
        ax.grid(axis='y', alpha=.15)
        ax.set_axisbelow(True)
    fig.text(.07, .959, 'High clouds in Terluna’s long twilight', fontsize=23, weight='bold')
    fig.text(.07, .923, 'A first mechanism study from committed climate and optical data', fontsize=12.5, color='#555555')
    fig.text(.07, .04, 'Equatorial equinox; smooth spherical ground; observer 1.6 m above sea level. '
                        'Zero time is ground-level solar-centre sunset.\n'
                        'Illumination, atmospheric transmission and cloud occurrence are separate diagnostics. '
                        'These panels do not predict cloud luminance, colour or visibility.',
             fontsize=10, color='#555555', linespacing=1.7)
    fig.subplots_adjust(left=.07, right=.96, top=.86, bottom=.15, hspace=.43, wspace=.27)
    fig.savefig(OUT / 'cloud_twilight.png', dpi=160, facecolor='white')
    manifest = dict(product=str(PRODUCT.relative_to(ROOT)),
                    product_sha256=hashlib.sha256(PRODUCT.read_bytes()).hexdigest(),
                    plot_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    panels=titles, evidence=data['evidence'])
    (OUT / 'cloud_twilight_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(OUT / 'cloud_twilight.png')


if __name__ == '__main__':
    main()
