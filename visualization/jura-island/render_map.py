#!/usr/bin/env python3
"""A map of Jura Island's topography and water depth at the common water level, from the flyover's terrain.

    python3 visualization/jura-island/render_map.py    # out/jura-island-topography.png and map-record.json

Writes to visualization/jura-island/out/ (kept out of Git), or to --out; the record holds the soundings and the
sha256 of everything read. The main map is a local equirectangular projection true to scale at 45 deg N. The 3-D
inset exaggerates relief 20 times and the profile has independent axis scales, and the sheet says so beside each.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, LightSource, LinearSegmentedColormap, ListedColormap, Normalize

import film
import terrain as terrain_data
from water_depth import DEPTH_COLOURS, DEPTH_LEVELS_M, colours

OUT = HERE / 'out'
WINDOW = dict(south=34.0, north=56.0, west=-52.0, east=-19.0)
PROFILE = dict(latitude=44.5, west=-42.0, east=-23.0)
INK, INK2 = '#213c46', '#55717b'
CODE = ('film.py', 'terrain.py', 'water_depth.py', 'render_map.py')


def draw(data, path):
    xscale, yscale = film.R_KM * np.pi / 180 * np.cos(np.radians(film.LAT0)), film.R_KM * np.pi / 180

    def px(lon):
        return (np.asarray(lon) - film.LON0) * xscale

    def py(lat):
        return (np.asarray(lat) - film.LAT0) * yscale

    rows = (data.lat >= WINDOW['south']) & (data.lat <= WINDOW['north'])
    cols = (data.lon >= WINDOW['west']) & (data.lon <= WINDOW['east'])
    lat, lon = data.lat[rows], data.lon[cols]
    z, island = data.z[np.ix_(rows, cols)], data.island[np.ix_(rows, cols)]
    X, Y = np.meshgrid(px(lon), py(lat))
    ppd = data.pixels_per_degree
    land = LinearSegmentedColormap.from_list('land', ['#e1d9bd', '#c5ba92', '#a89a70', '#8d7c5b', '#b3a590', '#f0eadd'])
    land_norm = Normalize(0, 3000)
    water = ListedColormap(DEPTH_COLOURS)
    water_norm = BoundaryNorm(DEPTH_LEVELS_M, len(DEPTH_COLOURS), clip=True)
    rgb = np.where((z > 0)[..., None], land(land_norm(z))[..., :3], colours(np.maximum(-z, 0)))
    shade = LightSource(azdeg=315, altdeg=38).hillshade(z, vert_exag=5, dx=xscale / ppd * 1000, dy=yscale / ppd * 1000)
    strength = np.where(z > 0, 0.37, 0.20)[..., None]
    rgb = np.clip(rgb * (1 - strength + strength * (0.4 + shade[..., None] * 0.9)), 0, 1)
    captions = film.captions(data.sea_level_m, data.water_share, ppd)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'text.color': INK, 'axes.labelcolor': INK,
                         'xtick.color': '#46616b', 'ytick.color': '#46616b', 'axes.edgecolor': '#9aafb3'})
    fig = plt.figure(figsize=(18, 12), facecolor='#f6f4ee')
    fig.text(.04, .955, 'JURA ISLAND', fontsize=29, fontweight='bold')
    fig.text(.04, .919, 'Topography and water depth at the common water level  ·  Open Moon', fontsize=17, color=INK2)
    fig.text(.955, .956, captions['share'], ha='right', fontsize=12, fontweight='bold')
    fig.text(.955, .930, captions['level'], ha='right', fontsize=10, color=INK2)

    ax = fig.add_axes([.04, .16, .608, .714])
    half_x, half_y = xscale * 0.5 / ppd, yscale * 0.5 / ppd
    ax.imshow(rgb, extent=[X.min() - half_x, X.max() + half_x, Y.min() - half_y, Y.max() + half_y], origin='upper',
              interpolation='bilinear')
    depth_lines = ax.contour(X, Y, np.where(z < 0, -z, np.nan), levels=DEPTH_LEVELS_M[1:-1], colors='#cae3e3',
                             linewidths=.5, alpha=.7)
    # long contours carry labels, which keeps the small craters' rings readable
    labels = ax.clabel(depth_lines, levels=[500, 750, 1000, 1250, 1500, 2000], fmt=lambda d: f'{int(d):,} m',
                       fontsize=6.8, inline=True, inline_spacing=3)
    for text in labels:
        text.set_path_effects([patheffects.withStroke(linewidth=1.8, foreground='#24495d')])
    height_lines = ax.contour(X, Y, np.where(z > 0, z, np.nan), levels=[500, 1000, 1500, 2000, 2500], colors='#675e45',
                              linewidths=.35, alpha=.48)
    ax.clabel(height_lines, levels=[1000, 2000], fmt='%d m', fontsize=6.2, inline=True)
    ax.contour(X, Y, z, levels=[0], colors='#e4eddb', linewidths=.45)
    coast = ax.contour(X, Y, island.astype(float), levels=[.5], colors='#fff1b8', linewidths=1.4)
    coast.set_path_effects([patheffects.Stroke(linewidth=2.4, foreground='#665a35'), patheffects.Normal()])
    ax.set_xlim(px(WINDOW['west']), px(WINDOW['east']))
    ax.set_ylim(py(WINDOW['south']), py(WINDOW['north']))
    ax.set_aspect('equal')
    xt, yt = np.arange(-50, -19, 5), np.arange(35, 56, 5)
    ax.set_xticks(px(xt), [f'{abs(i)}°W' for i in xt])
    ax.set_yticks(py(yt), [f'{i}°N' for i in yt])
    ax.tick_params(length=0, pad=6, labelsize=9)
    ax.grid(color='white', alpha=.16, linewidth=.5)
    halo = [patheffects.withStroke(linewidth=3, foreground='#f4efd9', alpha=.95)]
    sea_halo = [patheffects.withStroke(linewidth=2, foreground='#22425b')]

    def label(lon_, lat_, text, **kwargs):
        return ax.text(px(lon_), py(lat_), text, ha='center', va='center', **kwargs)

    label(-41.2, 43.2, 'JURA\nISLAND', fontsize=19, fontweight='bold', linespacing=1.2, path_effects=halo,
          color='#3b3d2e')
    label(-29.8, 45.8, 'SINUS IRIDUM', fontsize=11.5, fontweight='bold', color='#e5efeb', path_effects=sea_halo)
    label(-29.8, 45.3, 'Bay of Rainbows', fontsize=10, fontstyle='italic', color='#d3e9e8')
    label(-25.2, 36.5, 'MARE IMBRIUM', fontsize=14, fontweight='bold', color='#e5efeb', path_effects=sea_halo)
    label(-25.2, 35.9, 'Sea of Rains', fontsize=10, fontstyle='italic', color='#d3e9e8')
    label(-35.3, 48.2, 'MONTES JURA', fontsize=10, fontweight='bold', rotation=20, path_effects=halo, color='#4f4c37')

    soundings = terrain_data.soundings(data)
    for s in soundings:
        x, y = px(s['longitude_deg']), py(s['latitude_deg'])
        ax.plot(x, y, marker='o', ms=2.8, color='#fff7d2', mec='#244b61', mew=.4, zorder=12)
        offset = (6, -13) if s['name'] == 'Bay floor' else (6, 4)
        ax.annotate(f'{s["display_depth_m"]:,} m', xy=(x, y), xytext=offset, textcoords='offset points',
                    fontsize=8.8, fontweight='bold', color='#fff9dc', zorder=13,
                    bbox=dict(boxstyle='round,pad=.22', fc='#173d53', ec='none', alpha=.87))

    # the profile runs along 44.5 N from the western uplands into Imbrium
    section_lon = np.linspace(PROFILE['west'], PROFILE['east'], 650)
    section_lat = np.full_like(section_lon, PROFILE['latitude'])
    section_z = data.heights(section_lon, section_lat)
    section_km = (section_lon - section_lon[0]) * film.R_KM * np.pi / 180 * np.cos(np.radians(PROFILE['latitude']))
    ax.plot(px(section_lon), py(section_lat), '--', lw=1.3, color='#ffe4a0', alpha=.85, zorder=10)
    for name, lon_ in (('A', PROFILE['west']), ('B', PROFILE['east'])):
        ax.text(px(lon_), py(PROFILE['latitude']) - 13, name, ha='center', va='top', fontsize=10, fontweight='bold',
                color='#ffe4a0', zorder=12, path_effects=[patheffects.withStroke(linewidth=2, foreground='#264c5b')])
    summit_lon, summit_lat, summit_m = data.summit
    ax.scatter([px(summit_lon)], [py(summit_lat)], marker='^', s=47, c='#793e28', edgecolors='#fff4df',
               linewidths=.8, zorder=11)
    ax.annotate(f'Highest sampled ground\n{summit_m / 1000:.2f} km above the water',
                xy=(px(summit_lon), py(summit_lat)), xytext=(px(-46.2), py(48.6)), ha='center', fontsize=9,
                color='#4b422e', bbox=dict(boxstyle='round,pad=.4', fc='#f6f1dd', ec='none', alpha=.95),
                arrowprops=dict(arrowstyle='-', lw=1, color='#4b422e'))
    sx, sy = px(-50), py(35.3)
    ax.plot([sx, sx + 100], [sy, sy], color='#edf4df', lw=2.7)
    for xx in (sx, sx + 50, sx + 100):
        ax.plot([xx, xx], [sy - 4, sy + 4], color='#edf4df', lw=1)
    ax.text(sx, sy + 9, '0', ha='center', fontsize=8, color='#edf4df')
    ax.text(sx + 100, sy + 9, '100 km', ha='center', fontsize=8, color='#edf4df')
    ax.annotate('N', xy=(px(-50), py(54.5)), xytext=(px(-50), py(53.2)), ha='center', fontsize=12, fontweight='bold',
                arrowprops=dict(arrowstyle='-|>', color='#e1f2ed', lw=1.3), color='#e1f2ed')

    # the seafloor at its measured height, with the inset's exaggeration stated
    fig.text(.689, .862, 'LAND AND SEAFLOOR', fontsize=12, fontweight='bold')
    fig.text(.689, .840, 'Water surface removed  ·  vertical relief ×20', fontsize=9.5, color=INK2)
    ax3 = fig.add_axes([.66, .596, .33, .235], projection='3d', computed_zorder=False)
    ax3.set_facecolor('#f6f4ee')
    step = 2
    ax3.plot_surface(X[::step, ::step], Y[::step, ::step], z[::step, ::step] / 1000, facecolors=rgb[::step, ::step],
                     rcount=z[::step, ::step].shape[0], ccount=z[::step, ::step].shape[1], linewidth=0,
                     antialiased=False, shade=False)
    ax3.set_xlim(X.min(), X.max())
    ax3.set_ylim(Y.min(), Y.max())
    ax3.set_zlim(-3, 3.3)
    ax3.set_box_aspect([np.ptp(X), np.ptp(Y), 6.3 * 20])
    ax3.view_init(elev=48, azim=-63)
    ax3.set_axis_off()

    fig.text(.689, .590, 'A–B  FROM JURA INTO THE BAY', fontsize=11.5, fontweight='bold')
    axp = fig.add_axes([.717, .385, .238, .167])
    axp.set_facecolor('#f6f4ee')
    axp.fill_between(section_km, section_z, -2200, color='#b7a885')
    axp.fill_between(section_km, section_z, 0, where=section_z < 0, color='#80bcc8')
    axp.plot(section_km, section_z, color='#5a5b48', lw=1.1)
    axp.axhline(0, color='#28738c', lw=1)
    axp.text(section_km.max() - 3, 85, 'Water level', ha='right', fontsize=8, color='#28738c')
    axp.set_ylim(-1600, 3000)
    axp.set_xlim(0, section_km.max())
    axp.set_yticks([-1500, -1000, 0, 1000, 2000, 3000],
                   ['−1,500', '−1,000', '0', '1,000', '2,000', '3,000'])
    axp.set_xticks([0, 100, 200, 300, 400])
    axp.tick_params(labelsize=7.8, length=2)
    axp.set_ylabel('Metres from the water level', fontsize=8, labelpad=3)
    axp.set_xlabel('Distance from A (km)', fontsize=8, labelpad=3)
    axp.spines[['top', 'right']].set_visible(False)
    axp.grid(axis='y', alpha=.2, lw=.4)
    axp.text(4, 2700, 'A', fontsize=9, fontweight='bold')
    axp.text(section_km.max() - 4, 2700, 'B', ha='right', fontsize=9, fontweight='bold')
    fig.text(.689, .343, 'The profile\'s vertical and horizontal scales differ.', fontsize=8.5, color=INK2)

    fig.text(.689, .306, 'WATER DEPTH BELOW THE SURFACE', fontsize=10.5, fontweight='bold')
    cax = fig.add_axes([.689, .267, .266, .018])
    bar = fig.colorbar(ScalarMappable(norm=water_norm, cmap=water), cax=cax, orientation='horizontal',
                       ticks=DEPTH_LEVELS_M, spacing='uniform')
    bar.ax.set_xticklabels([f'{int(v):,}' for v in DEPTH_LEVELS_M[:-1]] + [f'{int(DEPTH_LEVELS_M[-1]):,}+'])
    bar.ax.tick_params(labelsize=6.7, length=2)
    bar.outline.set_visible(False)
    fig.text(.689, .229, 'Depths in metres  ·  bright dots mark the sampled depths.', fontsize=8.5, color=INK2)
    fig.text(.689, .198, 'LAND HEIGHT ABOVE THE SURFACE', fontsize=10.5, fontweight='bold')
    cax = fig.add_axes([.689, .163, .266, .015])
    bar = fig.colorbar(ScalarMappable(norm=land_norm, cmap=land), cax=cax, orientation='horizontal',
                       ticks=[0, 1000, 2000, 3000])
    bar.ax.set_xticklabels(['0', '1,000', '2,000', '3,000 m'])
    bar.ax.tick_params(labelsize=8, length=2)
    bar.outline.set_visible(False)

    notes = (f'Terrain: NASA LRO LOLA at {ppd} pixels per degree, about 1.9 km north–south. Water datum: the '
             'degree-200 GRAIL GL0420A geoid with the static tide and rotation.',
             'Depth is the vertical distance from the common water surface to the sampled terrain. Point labels are '
             'interpolated and rounded to the nearest 10 m.',
             'Common-water-level view: all ground below the level is flooded, enclosed crater floors included. Lake '
             'levels, tides, erosion and crustal loading are left out; the hillshade is illustrative.',
             'Pale outline: the island\'s coast. Dashed A–B line: the profile. North up; local equirectangular '
             'projection, true to scale at 45°N. Land contours every 500 m.')
    for i, note in enumerate(notes):
        fig.text(.04, .108 - .024 * i, note, fontsize=9, color='#526972')
    fig.savefig(path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    return soundings, dict(latitude_deg=PROFILE['latitude'], start_longitude_deg=PROFILE['west'],
                           end_longitude_deg=PROFILE['east'], distance_km=float(section_km[-1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = terrain_data.build()
    soundings, profile = draw(data, args.out / 'jura-island-topography.png')
    files = [HERE / f for f in CODE] + [ROOT / 'geography/topography.py', ROOT / 'geography/fetch_inputs.py']
    record = dict(schema='terluna.visualization.jura-map-render/1',
                  code={terrain_data.relative(p): terrain_data.digest(p) for p in files}, read=data.hashes,
                  sea_level_m=data.sea_level_m, water_share=data.water_share,
                  depth_levels_m=DEPTH_LEVELS_M.tolist(), soundings=soundings, profile=profile,
                  written=['jura-island-topography.png'])
    (args.out / 'map-record.json').write_text(json.dumps(record, indent=1) + '\n')
    print(args.out / 'jura-island-topography.png', flush=True)


if __name__ == '__main__':
    main()
