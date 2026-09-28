"""Figures for the sky-fleet study: how the flyers share the air above the summit metropolis, as a section to scale.

    python3 visualization/sky-fleet/figures.py        # build/heights.png

Read from the products: the ground along an east-west line through the summit (the geography atlas, 7.6 km cells),
the tower's height and width profile (summit tower study), each class's height band and the port's docks (sky-fleet
study), and the daytime mixed layer and storm tops of the cloud-resolving ring (climate). The ring is flat, at sea
level and on the equator, so its mixed layer is drawn over the metropolis's ground as a reading of its depth, not a
model of the summit's.
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FLEET = ROOT / 'research' / 'studies' / 'sky_fleet' / 'results' / 'sky_fleet.json'
TOWER = ROOT / 'research' / 'studies' / 'summit_tower' / 'results' / 'summit_tower.json'
ATLAS = ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz'
RING = ROOT / 'climate' / 'results' / 'crm' / 'ring_ring.json'
BUILD = HERE / 'build'
INK, MUTED, RULE = '#1f2a33', '#5d6b76', '#c3ccd3'
GROUND, TOWER_C = '#b9ad98', '#8c959c'
BANDS = {'city': '#3f7fb3', 'wings': '#d08a3a', 'soar': '#e8c26a', 'regional': '#6aa37a', 'band': '#5a6fb0',
         'fast': '#9a6ab0', 'platform': '#b0607a', 'storm': '#9aa6b0'}


def fonts():
    names = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams.update({'font.family': 'Inter' if 'Inter' in names else 'DejaVu Sans', 'font.size': 9.5,
                         'axes.edgecolor': RULE, 'axes.labelcolor': MUTED, 'xtick.color': MUTED,
                         'ytick.color': MUTED, 'text.color': INK})


def section(summit_lat, summit_lon, half_km=60.0):
    """The ground along the east-west line through the summit, in km above sea level, against distance in km."""
    a = np.load(ATLAS, allow_pickle=False)
    meta = json.loads(str(a['metadata']))
    lat, lon, h = a['lat_deg'], a['lon_deg'], a['height_m']
    i = int(np.argmin(np.abs(lat - summit_lat)))
    km_per_deg = np.pi / 180.0 * 1737.4 * np.cos(np.radians(summit_lat))
    x = ((lon - summit_lon + 180.0) % 360.0 - 180.0) * km_per_deg
    order = np.argsort(x)
    x, z = x[order], (h[i, order] - meta['sea_level_m']) / 1e3
    keep = np.abs(x) <= half_km + 10.0
    return x[keep], z[keep]


def main():
    fonts()
    fleet = json.loads(FLEET.read_text())
    tower = json.loads(TOWER.read_text())
    ring = json.loads(RING.read_text())
    site = tower['site']
    design = tower['port']['design']
    H = design['height_km']
    top_w = fleet['long_haul']['crown']['top_width_m'] / 1e3
    base_w, p = design['base_width_m'] / 1e3, design['flare_exponent']
    ground0 = site['ground_above_sea_m'] / 1e3
    lon0 = site['lon_deg']
    x, zg = section(site['lat_deg'], lon0)

    fig, ax = plt.subplots(figsize=(11.0, 8.2), dpi=200)
    half = 60.0
    ax.set_xlim(-half, half)
    ax.set_ylim(0, 80)
    ax.set_xlabel('Distance from the tower along an east–west line (km)')
    ax.set_ylabel('Height above sea level (km)')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    # Heights of the air's features over the whole column
    storms = ring['storms']['cloud_top_km']
    ax.axhspan(35.0, 45.0, color=BANDS['band'], alpha=0.10, lw=0)
    ax.axhline(storms['median'], color=BANDS['storm'], lw=1.0, ls=(0, (5, 3)))
    ax.axhline(storms['p90'], color=BANDS['storm'], lw=0.8, ls=(0, (2, 3)))
    ax.text(half - 1, storms['median'] + 0.6, f"storm tops, median {storms['median']:.0f} km", ha='right',
            va='bottom', fontsize=8.5, color=MUTED)
    ax.text(-half + 1, storms['p90'] + 0.6, f"storm tops, tallest tenth {storms['p90']:.0f} km and more", ha='left',
            va='bottom', fontsize=8.5, color=MUTED)

    # The daytime mixed layer (the ring's, over land at its deepest) over the metropolis's ground
    soar = fleet['people_on_wings']['soaring']
    city_r = fleet['metropolis']['city']['radius_km']
    inside = np.abs(x) <= city_r
    zi = np.interp(np.linspace(-city_r, city_r, 200), x, zg)
    xs = np.linspace(-city_r, city_r, 200)
    ax.fill_between(xs, zi, zi + soar['mixed_layer_max_km'], color=BANDS['soar'], alpha=0.22, lw=0)
    ax.text(-city_r + 1.5, np.interp(-city_r + 1.5, xs, zi) + soar['mixed_layer_max_km'] - 1.5,
            f"thermals by day: the mixed layer\nup to {soar['mixed_layer_max_km']:.0f} km deep (the ring's land)",
            fontsize=8, color='#8a6a1f', va='top')

    # The city's flyers: bands above the ground, drawn thick enough to read
    for name, colour, key in (('wings, gliders, canopies', BANDS['wings'], 'wings, gliders and canopies'),
                              ('sky boats and air taxis', BANDS['city'], 'sky boats and air taxis'),
                              ('sky ferries', '#2c5f87', 'sky ferries')):
        b = fleet['metropolis']['bands'][key]
        ax.fill_between(xs, zi + b['bottom_m'] / 1e3, zi + b['top_m'] / 1e3, color=colour, alpha=0.85, lw=0)

    # Ground
    ax.fill_between(x, 0, zg, color=GROUND, lw=0)
    ax.plot(x, zg, color='#8c8170', lw=0.8)
    for edge in (-city_r, city_r):
        ze = float(np.interp(edge, x, zg))
        ax.plot([edge, edge], [ze - 2.2, ze - 0.4], color='#5f5646', lw=1.0)
    ax.annotate('', xy=(city_r, 3.0), xytext=(-city_r, 3.0),
                arrowprops=dict(arrowstyle='<|-|>', color='#5f5646', lw=0.9, shrinkA=0, shrinkB=0))
    ax.text(0, 3.6, f'the metropolis, within {city_r:.0f} km of the tower', ha='center', fontsize=8.5,
            color='#3d372c')

    # The tower, to scale
    zz = np.linspace(0, H, 200)
    w = top_w + (base_w - top_w) * (1 - zz / H) ** p
    ax.fill_betweenx(ground0 + zz, -w / 2, w / 2, color=TOWER_C, lw=0)
    reg = fleet['long_haul']['crown']
    ax.plot([0.9, 7.0], [ground0 + 12.5, ground0 + 12.5], color=BANDS['regional'], lw=1.0)
    ax.text(7.3, ground0 + 12.5, 'regional docks, 10–15 km above the summit', va='center', fontsize=8.5,
            color='#3f7a50')
    ax.plot([0.3, 7.0], [35.3, 35.3], color=BANDS['band'], lw=1.0)
    ax.text(7.3, 35.3, 'the crown: long-haul berths in the flight band', va='center', fontsize=8.5, color='#3d4f8c')
    ax.text(-half + 1, 44.3, 'flight band, 35–45 km: long-haul liners and freighters', fontsize=8.5,
            color='#3d4f8c', va='top')
    fast = fleet['air']['fast_cruise']['height_km']
    ax.axhline(fast, xmin=0.02, xmax=0.45, color=BANDS['fast'], lw=2.2)
    ax.text(-half + 1, fast + 0.8, f'winged liners cruise near {fast:.0f} km', fontsize=8.5, color='#6d4585')
    plat = fleet['air']['platforms']['height_km']
    ax.axhline(plat, xmin=0.02, xmax=0.30, color=BANDS['platform'], lw=2.2)
    ax.text(-half + 1, plat + 0.8, f'high platforms near {plat:.0f} km', fontsize=8.5, color='#8a4059')

    # Legend for the city's bands
    y0 = 62.0
    for i, (label, colour, key) in enumerate((('wings, gliders, canopies', BANDS['wings'], 'wings, gliders and canopies'),
                                              ('sky boats and air taxis', BANDS['city'], 'sky boats and air taxis'),
                                              ('sky ferries', '#2c5f87', 'sky ferries'))):
        b = fleet['metropolis']['bands'][key]
        ax.add_patch(plt.Rectangle((18, y0 - 2.4 * i), 3.2, 1.3, color=colour, lw=0))
        ax.text(22.2, y0 - 2.4 * i + 0.65, f"{label}: {b['bottom_m'] / 1e3:.2g}–{b['top_m'] / 1e3:.2g} km above the ground",
                va='center', fontsize=8.5)
    ax.text(18, y0 + 3.0, 'Over the metropolis', fontsize=9.5, fontweight='bold')
    ax.set_title('How the flyers share the air above the summit metropolis (to scale)', loc='left', fontsize=12,
                 fontweight='bold', color=INK, pad=12)
    fig.text(0.01, 0.005, ('Ground along 5.4° N from the geography atlas (7.6 km cells); the tower from the summit tower '
                           'study; bands from the sky-fleet study; mixed layer and storm tops from the cloud-resolving '
                           'ring (flat, equatorial, at sea level).'), fontsize=7.5, color=MUTED)
    BUILD.mkdir(parents=True, exist_ok=True)
    out = BUILD / 'heights.png'
    fig.savefig(out, bbox_inches='tight', facecolor='white')
    print(out)


if __name__ == '__main__':
    main()
