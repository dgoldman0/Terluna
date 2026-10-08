"""Render saved finite-fleet products; the study owns all geometry and dynamics."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PRODUCTS = ROOT/'research/studies/solar_shield_array/results'
OUTPUT = Path(__file__).with_name('results')/'fleet.png'


def main():
    names = ['fleet_dense_double.json', 'fleet_refinement.json', 'fleet_handover.json']
    products = {name: json.loads((PRODUCTS/name).read_text()) for name in names}
    doubled = products[names[0]]['cases']['dense_double_1000km']
    refined = products[names[1]]['cases']['dense_1000km']
    minute = products[names[2]]
    variant = refined['variants']['4R_nominal']
    teal, orange, red = '#007f79', '#b06c00', '#ad4150'
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.titleweight': 'bold'})
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.subplots_adjust(top=.84, bottom=.14, left=.08, right=.95, hspace=.43, wspace=.32)
    fig.suptitle('Finite-Sun fleet coverage: gaps remain', x=.08, y=.975,
                 ha='left', fontsize=21, fontweight='bold')
    fig.text(.08, .933, '4 lunar radii · 50 g/m² · common-epoch DE440 histories · 24 orbital planes', fontsize=11)
    fig.text(.08, .901, '1,000 km square geometry probes. Nominal populations violate collision and mutual-shadow guards.',
             fontsize=10, color='#555555')

    ax = axes[0, 0]
    t = refined['times_days']
    for key, label, color in [('ray_coverage', 'Mean over aperture and solar disk', teal),
                              ('all_sampled_sun_coverage', 'Area covered for all sampled Sun directions', orange)]:
        ax.plot(t, [100*r[key] for r in variant['records']], color=color, lw=1.1, label=label)
    strict = refined['variants']['4R_collision_and_shadow_screened']
    ax.plot(t, [100*r['ray_coverage'] for r in strict['records']], color=red, lw=1.1,
            label='17-member screened subset: ray interception')
    ax.set(xlabel='Days after 2026-10-04 TDB', ylabel='Coverage (%)', ylim=(0, 100), xlim=(0, 30),
           title='1,536 squares · 3-hour / 32-Sun samples')
    ax.legend(fontsize=7.5, loc='upper right', frameon=False)
    ax.grid(axis='y', alpha=.2)

    ax = axes[0, 1]
    records = doubled['selections']['nominal']['records']; t = doubled['coverage_times_days']
    actual = np.array([r['ray_coverage'] for r in records])*100
    summed = actual+np.array([r['excess_intersection_area_fraction'] for r in records])*100
    ax.plot(t, summed, color='#77828a', lw=1.1, label='Sum of intercepted areas (overlap counted)')
    ax.plot(t, actual, color=teal, lw=1.3, label='Union: actual ray interception')
    ax.axhline(100, color=red, ls='--', lw=.8)
    ax.set(xlabel='Days after epoch', ylabel='Fraction of aperture (%)', xlim=(0, 30),
           title='3,072 squares · summed area hides gaps')
    ax.legend(fontsize=7.5, loc='upper right', frameon=False)
    ax.grid(axis='y', alpha=.2)

    ax = axes[1, 0]
    spatial = variant['spatial_map_at_worst_date']; axis = np.array(spatial['axis_km'])/1000
    field = np.ma.masked_less(np.array(spatial['intercepted_solar_fraction']), 0)
    extent = (axis[0], axis[-1], axis[0], axis[-1])
    im = ax.imshow(field, extent=extent, origin='lower', vmin=0, vmax=1, cmap='viridis', interpolation='nearest')
    ax.set(xlabel='Aperture coordinate (1,000 km)', ylabel='Aperture coordinate (1,000 km)',
           title=f"1,536 squares · spatial gaps on day {variant['worst_date_days']:.3f}")
    colorbar = fig.colorbar(im, ax=ax, fraction=.042, pad=.025)
    colorbar.set_label('Fraction of 64 solar directions intercepted', fontsize=8)
    ax.text(.02, .02, f"{100*spatial['points_with_zero_filtering_fraction']:.1f}% of sampled positions entirely unfiltered",
            transform=ax.transAxes, color='white', fontsize=8,
            bbox={'facecolor': '#333333', 'alpha': .8, 'edgecolor': 'none'})

    ax = axes[1, 1]
    hours = (np.array(minute['times_days'])-minute['times_days'][0])*24
    ax.plot(hours, [100*r['ray_coverage'] for r in minute['records']], color=teal, lw=1.4, label='Ray interception')
    ax.plot(hours, [100*r['all_sampled_sun_coverage'] for r in minute['records']], color=orange,
            lw=1.2, label='All sampled Sun directions')
    ax.set(xlabel='Hours after day 26.0', ylabel='Coverage (%)', ylim=(0, 100), xlim=(0, 12),
           title='1,536 squares · local 1-minute / 32-Sun audit')
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis='y', alpha=.2)

    fig.text(.08, .080, 'Lines join computed samples; they do not certify coverage between them. The map is a spatial sample of the same finite-ray model.', fontsize=9)
    fig.text(.08, .057, 'Centre-force approximation: large-square finite-extent gravity is comparable to available photon thrust. Structural dynamics are omitted.', fontsize=9)
    fig.text(.08, .034, 'Independent 10 km phase replay differs by 898 m (50 m budget). No operational fleet, collector/habitat load or delivered power is demonstrated.', fontsize=9)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=170, facecolor='white')
    plt.close(fig)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    metadata = dict(products={str((PRODUCTS/n).relative_to(ROOT)): digest(PRODUCTS/n) for n in names},
                    renderer_sha256=digest(Path(__file__)), image_sha256=digest(OUTPUT),
                    reading_rule='Displays saved products only; finite sampling, collision exclusions and centre-force limitations remain attached.')
    OUTPUT.with_suffix('.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(OUTPUT)


if __name__ == '__main__':
    main()
