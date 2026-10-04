"""Display the optical-comfort product; no physical model lives in this plot."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / 'research/studies/optical_comfort/results/optical_comfort.json'
OUT = Path(__file__).resolve().parent / 'figures'
COLOURS = {'moon': '#087f83', 'earth': '#6a617a'}
LABELS = {'moon': 'Terluna', 'earth': 'Earth control'}


def main():
    data = json.loads(PRODUCT.read_text())
    if data['schema'] != 'terluna.research.optical-comfort/1':
        raise ValueError('Unsupported optical-comfort schema')
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titleweight': 'bold', 'axes.labelcolor': '#333333'})
    fig, grid = plt.subplots(2, 2, figsize=(13, 10))
    axes = grid.ravel()
    for world in ('moon', 'earth'):
        rows = sorted([r for r in data['surface_scenes'] if r['world'] == world
                       and r['sun_deg'] == 90 and r['landscape_albedo'] in (0.1, 0.3, 0.5, 0.8)],
                      key=lambda r: r['landscape_albedo'])
        x = [r['landscape_albedo'] for r in rows]
        colour = COLOURS[world]
        axes[0].plot(x, [r['ground_luminance_cd_m2'] / 1000 for r in rows],
                     marker='o', color=colour, label=LABELS[world], lw=2.5)
        axes[1].plot(x, [r['vertical_ambient_lux'] / 1000 for r in rows],
                     marker='o', color=colour, lw=2.5)
        axes[1].fill_between(x, [r['vertical_zenith_bright_lux'] / 1000 for r in rows],
                             [r['vertical_horizon_bright_lux'] / 1000 for r in rows], color=colour, alpha=0.13)
        axes[2].plot(x, [r['disk_only_shadow_fraction'] * 100 for r in rows],
                     marker='o', color=colour, lw=2.5)
    titles = ['Horizontal ground is dimmer', 'Level-gaze light under imposed skies',
              'Solar-only shadows retain more light']
    labels = ['Ground luminance (thousand cd/m²)', 'Vertical ambient illuminance (klux)',
              'Shadow / sunlit illumination (%)']
    limits = [(0, 38), (0, 100), (0, 65)]
    for ax, title, label, ylim in zip(axes, titles, labels, limits):
        ax.set_title(title, loc='left', fontsize=13, pad=17)
        ax.set_xlabel('Landscape albedo')
        ax.set_ylabel(label)
        ax.set_xticks([0.1, 0.3, 0.5, 0.8])
        ax.set_ylim(*ylim)
        ax.grid(axis='y', alpha=0.18)
        ax.set_axisbelow(True)
    axes[0].legend(loc='upper left', frameon=False)
    axes[1].text(0.03, 0.95, 'Line: uniform sky\nBand: three illustrative sky shapes',
                 transform=axes[1].transAxes, va='top', fontsize=9, color='#555555')
    axes[2].text(0.03, 0.95, 'Only the solar disk is screened.\nSurroundings remain illuminated.',
                 transform=axes[2].transAxes, va='top', fontsize=9, color='#555555')
    for albedo, colour in ((0.1, '#454b50'), (0.8, '#b57928')):
        selected = [r for r in data['gaze_scenes'] if r['sun_deg'] == 90
                    and r['landscape_albedo'] == albedo and r['azimuth_from_sun_deg'] == 0 and r['gaze_deg'] <= 0]
        gazes = sorted({r['gaze_deg'] for r in selected})
        lookup = {(r['world'], r['gaze_deg']): r['total_lux'] for r in selected}
        ratios = [lookup['moon', e] / lookup['earth', e] for e in gazes]
        axes[3].plot(gazes, ratios, marker='o', color=colour, lw=2.5, label=f'Albedo {albedo}')
    axes[3].axhline(1, color='#777777', lw=1, linestyle=':')
    axes[3].set(title='Looking down can reverse the comparison', xlabel='Gaze elevation (degrees)',
                ylabel='Moon / Earth eye-plane illuminance', ylim=(0.5, 2.4))
    axes[3].set_xticks([-90, -60, -30, 0], ['−90\nStraight down', '−60', '−30', '0\nLevel'])
    axes[3].text(0.03, 0.95, 'Uniform sky; solar disk outside these views', transform=axes[3].transAxes,
                 va='top', fontsize=9, color='#555555')
    axes[3].legend(frameon=False, loc='upper left', bbox_to_anchor=(0, 0.87))
    axes[3].grid(axis='y', alpha=0.18)
    axes[3].set_axisbelow(True)
    fig.suptitle('Terluna daylight: surfaces, shadows and gaze', x=0.08, y=0.97,
                 ha='left', fontsize=20, fontweight='bold')
    fig.text(0.08, 0.93, 'Clear sky • Sun overhead • neutral matte landscapes • 1.2-atm design Moon',
             fontsize=12, color='#555555')
    fig.text(0.08, 0.025, 'Conditional first pass. Eye values depend on assumed angular fields and gaze.\n'
                        'Sky-shape bands are scenarios, not bounds. No human discomfort threshold is inferred.',
             fontsize=10, color='#555555', linespacing=1.6)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.84, bottom=0.16, wspace=0.30, hspace=0.63)
    for extension in ('png', 'svg'):
        fig.savefig(OUT / f'optical_comfort.{extension}', dpi=160, facecolor='white')
    manifest = dict(product=str(PRODUCT.relative_to(ROOT)),
                    product_sha256=hashlib.sha256(PRODUCT.read_bytes()).hexdigest(),
                    plot_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    reading='Noon surface rows for albedos 0.1, 0.3, 0.5, 0.8; downward/level gaze ratios '
                            'for albedos 0.1 and 0.8. Lines join samples for display only.')
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(OUT / 'optical_comfort.png')


if __name__ == '__main__':
    main()
