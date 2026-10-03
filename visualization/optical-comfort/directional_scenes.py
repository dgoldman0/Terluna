"""Display stored directional-sky and finite-scene results without modelling."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / 'research/studies/optical_comfort/results/directional_scenes.json'
OUT = Path(__file__).resolve().parent / 'figures'
COLOURS = {'moon': '#087f83', 'earth': '#796484'}
LABELS = {'moon': 'Terluna', 'earth': 'Earth control'}
NAMES = ['open_dark', 'open_pale_patch', 'pale_courtyard', 'canopy_dark',
         'canopy_pale', 'screened_canopy_pale']
TICKS = ['Dark\nopen', 'Pale\npatch', 'Pale\ncourtyard', 'Canopy\ndark floor',
         'Canopy\npale floor', 'Screened\npale floor']


def main():
    data = json.loads(PRODUCT.read_text())
    if data['schema'] != 'terluna.research.optical-comfort-directional-scenes/1':
        raise ValueError('Unsupported scene product')
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titleweight': 'bold'})
    fig, axes = plt.subplots(2, 2, figsize=(15, 11.5))
    a, b, c, d = axes.ravel()
    views = {(r['world'], r['gaze_deg']): r for r in data['angular_views']
             if r['sun_deg'] == 90 and r['azimuth_from_sun_deg'] == 90 and r['gaze_deg'] <= 0}
    gazes = sorted(g for w, g in views if w == 'moon')
    for key, label, style, colour in [('ambient_lux', 'Atlas angular patterns', '-', '#087f83'),
                                      ('uniform_ambient_lux', 'Uniform sky', '--', '#92989d')]:
        values = [views['moon', g][key] / views['earth', g][key] for g in gazes]
        a.plot(gazes, values, style, color=colour, marker='o', lw=2.3, label=label)
    a.axhline(1, color='#444444', lw=.8)
    a.set(title='Gaze changes the Moon / Earth comparison', xlabel='Gaze elevation (degrees; negative looks down)',
          ylabel='Ambient eye illuminance, Moon / Earth', ylim=(.6, 2.3))
    a.set_xticks([-90, -60, -30, 0])
    a.legend(frameon=False, loc='upper left')
    a.text(.03, .57, 'Noon, open dark ground\nLevel gaze: 2.11 → 1.83', transform=a.transAxes,
           color='#555555', fontsize=10, va='top')
    rows = {(r['world'], r['sun_deg'], r['scene']): r for r in data['scenes']}
    x = np.arange(len(NAMES))
    for world, offset in [('moon', -.18), ('earth', .18)]:
        eye = [rows[world, 90., n]['probes']['eye_level']['ambient_lux'] / 1000 for n in NAMES]
        b.bar(x + offset, eye, .34, color=COLOURS[world], label=LABELS[world])
        contrast = [rows[world, 90., n]['tasks']['step_signed_michelson'] for n in NAMES]
        d.bar(x + offset, contrast, .34, color=COLOURS[world])
    b.set(title='Local surfaces and enclosure shape the view', ylabel='Level-eye ambient illumination at noon (klux)',
          ylim=(0, 51))
    b.legend(frameon=False, loc='upper left', ncol=2)
    for h, offset, colour in [(90., -.18, '#087f83'), (30., .18, '#cf9b45')]:
        base = rows['moon', h, 'open_dark']['probes']['workplane']['total_lux']
        values = [100 * rows['moon', h, n]['probes']['workplane']['total_lux'] / base for n in NAMES]
        c.bar(x + offset, values, .34, color=colour, label=f'Sun {h:g}°')
    c.set(title='A roof can lose its shade as the Sun lowers', ylabel='Lunar workplane light / open workplane (%)',
          ylim=(0, 119))
    c.legend(frameon=False, loc='upper left', ncol=2)
    d.set(title='Step contrast depends on the setting', ylabel='Signed tread–riser Michelson contrast at noon',
          ylim=(-.12, .97))
    d.axhline(0, color='#444444', lw=.8)
    d.text(.02, .98, 'Same 0.2 reflectance on both faces.\nNegative: riser brighter than tread.',
           transform=d.transAxes, va='top', fontsize=9, color='#555555')
    for ax in (b, c, d):
        ax.set_xticks(x, TICKS)
        ax.tick_params(axis='x', labelsize=9)
    for ax in axes.ravel():
        ax.set_title(ax.get_title(), loc='left', fontsize=13, pad=14)
        ax.set_title('')  # remove the centred copy; retain the left-hand title
        ax.grid(axis='y', alpha=.17)
        ax.set_axisbelow(True)
    fig.text(.07, .96, 'Optical comfort: sky direction, shelter and surface contrast', fontsize=21, weight='bold')
    fig.text(.07, .92, 'Current light totals × transferred clear-sky angular patterns • six prescribed local scenes',
             fontsize=12, color='#555555')
    fig.text(.07, .033, 'All atmospheric boundaries use landscape albedo 0.1. Pale patches are finite, 8 × 8 m. '
                        'Canopies are 6 × 6 m at 3 m height.\n'
                        'Reflected light and cast shadows are included. Geometry and reflectance are scenarios; '
                        'visual comfort and detection thresholds are unvalidated.',
             fontsize=9.7, color='#555555', linespacing=1.7)
    fig.subplots_adjust(left=.07, right=.97, top=.85, bottom=.14, hspace=.45, wspace=.28)
    for extension in ('png', 'svg'):
        fig.savefig(OUT / f'directional_scenes.{extension}', dpi=160, facecolor='white')
    manifest = dict(product=str(PRODUCT.relative_to(ROOT)),
                    product_sha256=hashlib.sha256(PRODUCT.read_bytes()).hexdigest(),
                    plot_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    reading='Noon open-ground gaze ratios; six noon scene eye values; lunar workplane '
                            'fractions at solar elevations 90 and 30 degrees; noon signed step contrasts.')
    (OUT / 'directional_scenes_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(OUT / 'directional_scenes.png')


if __name__ == '__main__':
    main()
