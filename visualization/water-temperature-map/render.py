#!/usr/bin/env python3
"""The Open Moon's waters' warmth as one sheet: the open water's body on average and in its hottest hundredth of
hours, still water 4 m and 2 m deep in theirs, and the shares of the seas, the sea coast and the rain-fed lakes in each
warmth class.

    python3 visualization/water-temperature-map/render.py [--out DIR]   # PNG and a manifest; default the research drive

Reads research/studies/littoral_community/results/water_map.json (schema terluna.research.littoral-water-map/1) and
draws its cells as they stand, in the run's own frame (the turned Moon's): each T21 cell a box of its class's colour,
no smoothing or downscaling, as the moisture map's sheet draws its cells. The climate model's water cells are outlined:
inside the outline a cell shows the sea and its shores, outside it the lakes and ponds of a land cell at its mean
ground level. The sheet's footnote gives the product's own check on how to read a mixed layer: the lake model's calm
column at five sites, whose surface lies above the body and between the two still layers. The manifest records the
product's hash and a faithfulness check, the seas' and the land's class shares counted from the drawn grids against
the same shares counted from the product's cells.

Colours: one ordered ramp for every panel, eight classes of 2 K from under 28 C to 40 C and more, a single orange hue
family falling in lightness (validated as an ordered ramp: lightness monotone, neighbouring steps at least 0.06 apart
in OKLCH lightness, hue spread 34 degrees). Its lightest step recedes toward the page, as a map's lowest class may;
the map's edge and the outline carry the boundary. The images stay out of Git.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import textwrap
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
import matplotlib.patheffects as patheffects                     # noqa: E402
from matplotlib.collections import LineCollection                # noqa: E402
from matplotlib.patches import Patch                             # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / 'research' / 'studies' / 'littoral_community' / 'results' / 'water_map.json'
SCHEMA = 'terluna.research.littoral-water-map/1'
OUT = Path('/media/projectspace/terluna-research/session-tools/water_temperature_map')
RAMP = ['#ffeedf', '#fbd2b3', '#f6b387', '#ef925c', '#e06e36', '#c44e20', '#9c351c', '#6e231a']
INK, INK2, INK3 = '#1f1f1d', '#52514e', '#8a8984'
OUTLINE = '#1f1f1d'
WIDTH_IN, HEIGHT_IN = 13.5, 14.6
# the product's share, its statistic, the row's label, what the share is of, and the row's group
BARS = (('seas_surface', 'mean_c', "The seas' upper 50 m, on average", "of the seas' area", 0),
        ('seas_surface', 'p99_c', "The seas' upper 50 m, hottest 1% of hours", "of the seas' area", 0),
        ('seas_still_4m', 'p99_c', 'On the seas: a still layer 4 m deep, hottest 1%', "of the seas' area", 0),
        ('seas_still_2m', 'p99_c', 'On the seas: a still layer 2 m deep, hottest 1%', "of the seas' area", 0),
        ('sea_coast_shore_2m_x1', 'p99_c', "Sea coast: 2 m shore water changing all of its water a day with the sea's body, hottest 1%",
         "of the coast's length", 1),
        ('sea_coast_shore_2m_x0.3', 'p99_c', 'Sea coast: 2 m shore water changing 30% a day, hottest 1%', "of the coast's length", 1),
        ('sea_coast_still_2m', 'p99_c', 'Sea coast: still shallows 2 m deep, hottest 1%', "of the coast's length", 1),
        ('lakes_mixed_10m', 'p99_c', 'Rain-fed lakes: the upper 10 m as one mixed layer, hottest 1%', "of the lakes' area", 2),
        ('lakes_still_4m', 'p99_c', 'Rain-fed lakes: still water 4 m deep, hottest 1%', "of the lakes' area", 2),
        ('lakes_still_2m', 'p99_c', 'Rain-fed lakes: still water 2 m deep, hottest 1%', "of the lakes' area", 2))
GROUP_GAP = 0.35


def grid(product, values):
    g = product['grid']
    return np.asarray(values, float).reshape(len(g['lat']), len(g['lon']))


def shifted(product, a):
    """Columns reordered to run from 180 W to 180 E, the dateline column repeated at both ends, with the edges."""
    lon = np.array(product['grid']['lon'])
    lon180 = (lon + 180.0) % 360.0 - 180.0
    order = np.argsort(lon180)
    step = 360.0 / lon.size
    centres = np.append(lon180[order], lon180[order][0] + 360.0)
    a2 = np.concatenate([a[:, order], a[:, order][:, :1]], axis=1)
    edges = np.append(centres - step / 2.0, centres[-1] + step / 2.0)
    return a2, edges


def classes_of(product, t_c):
    return np.searchsorted(np.array(product['classes']['edges_c']), np.asarray(t_c), side='right')


def draw_map(ax, product, values, sea, title):
    """One field's classes on an equal-area (Mollweide) map, the near side at the centre and the limbs dashed, the
    climate model's water cells outlined. Returns the classes as drawn on the product's own grid."""
    drawn = classes_of(product, values)
    a2, xedges = shifted(product, drawn.astype(float))
    s2, _ = shifted(product, sea.astype(float))
    yedges = np.array(product['grid']['lat_edges'])
    xr, yr = np.radians(np.clip(xedges, -180.0, 180.0)), np.radians(yedges)
    cmap = matplotlib.colors.ListedColormap(RAMP)
    ax.pcolormesh(xr, yr, a2, cmap=cmap, vmin=-0.5, vmax=len(RAMP) - 0.5, shading='flat', edgecolors='none')
    segments = []
    water = s2 > 0.5
    for i, j in zip(*np.nonzero(water[:, 1:] != water[:, :-1])):          # between neighbours along a row
        segments.append([(xr[j + 1], yr[i]), (xr[j + 1], yr[i + 1])])
    for i, j in zip(*np.nonzero(water[1:, :] != water[:-1, :])):          # between neighbours along a column
        if xr[j + 1] > xr[j]:
            segments.append([(xr[j], yr[i + 1]), (xr[j + 1], yr[i + 1])])
    ax.add_collection(LineCollection(segments, colors=OUTLINE, linewidths=0.45, alpha=0.9, zorder=4))
    for x in (-90.0, 90.0):
        ax.plot(np.radians([x, x]), np.radians([-89.9, 89.9]), color=INK2, lw=0.7, ls=(0, (4, 3)), zorder=5)
    ax.grid(True, color='white', lw=0.35, alpha=0.6)
    ax.set_xticklabels([])
    for lon, name in ((-90.0, '90° W'), (0.0, '0°'), (90.0, '90° E')):
        ax.text(np.radians(lon), np.radians(-3.0), name, ha='center', va='top', fontsize=5.8, color=INK, zorder=6,
                path_effects=[patheffects.withStroke(linewidth=1.8, foreground='white')])
    ax.tick_params(axis='y', labelsize=6, colors=INK2)
    ax.set_title(title, loc='left', fontsize=9.2, color=INK, pad=6)
    return drawn


def area_shares(product, drawn, mask):
    """Each class's share of the masked cells' area, from the drawn classes and the grid's own cell areas."""
    edges = np.radians(np.array(product['grid']['lat_edges']))
    band = (np.sin(edges[:-1]) - np.sin(edges[1:]))[:, None] * np.ones(drawn.shape)
    total = band[mask].sum()
    return [float(band[mask & (drawn == k)].sum() / total) for k in range(len(RAMP))]


def whole(values):
    """A range of temperatures to the whole degree, as '31–37'."""
    return f'{min(values):.0f}–{max(values):.0f}'


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(OUT))
    args = ap.parse_args(argv)
    p = json.loads(PRODUCT.read_text())
    if p['schema'] != SCHEMA:
        raise ValueError(f"expected {SCHEMA}, found {p['schema']}")
    c = p['cells']
    sea = grid(p, c['sea']) > 0.5
    names = p['classes']['names']
    if len(names) != len(RAMP):
        raise ValueError('the ramp and the classes differ in number')
    body_mean = np.where(sea, grid(p, c['sea_surface']['mean_c']), grid(p, c['mixed_10m']['mean_c']))
    body_p99 = np.where(sea, grid(p, c['sea_surface']['p99_c']), grid(p, c['mixed_10m']['p99_c']))
    still4_p99, still2_p99 = grid(p, c['still_4m']['p99_c']), grid(p, c['still_2m']['p99_c'])
    apart = max(float(np.abs(grid(p, c[k]['mean_c']) - body_mean).max()) for k in ('still_4m', 'still_2m'))

    top = lambda inches: 1.0 - inches / HEIGHT_IN                 # a height measured down from the sheet's top edge
    fig = plt.figure(figsize=(WIDTH_IN, HEIGHT_IN), dpi=150, facecolor='white')
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.06, left=0.07, right=0.985, top=top(1.37), bottom=top(8.34))
    gs_bars = fig.add_gridspec(1, 1, left=0.07, right=0.985, top=top(8.99), bottom=1.12 / HEIGHT_IN)
    s = p['settings']
    fig.text(0.02, top(0.43), "How warm the Open Moon's waters run at the design climate", fontsize=14, color=INK, weight='bold')
    heading = textwrap.wrap(
        f"{s['run']}, years {s['years'][0]}–{s['years'][1]}. T21 cells about 170 km across, drawn as computed on an "
        f"equal-area map in the turned Moon's frame, the near side centred and the limbs dashed. Outlined: the climate "
        f"model's water cells, the seas and large lakes; the other cells show the lakes and ponds of the land at each "
        f"cell's mean ground level. Each value is a mixed layer's temperature under its cell's sunlight, sky, air and "
        f"wind; still water exchanges nothing with other water. On average the still layers of C and D run within "
        f"{math.ceil(apart * 10.0) / 10.0:.1f} K of A.", width=194)
    if len(heading) > 3:
        raise ValueError('the heading no longer fits its three lines')
    fig.text(0.02, top(0.58), '\n'.join(heading), fontsize=8.6, color=INK2, linespacing=1.5, va='top')

    drawn = {}
    for (row, col), key, values, title in (
            ((0, 0), 'body_mean', body_mean, "A. The open water's body, on average: the sea's upper 50 m; on land, a lake's upper 10 m"),
            ((0, 1), 'body_p99', body_p99, "B. The open water's body in its hottest 1% of hours"),
            ((1, 0), 'still4_p99', still4_p99, 'C. Still water 4 m deep in its hottest 1% of hours'),
            ((1, 1), 'still2_p99', still2_p99, 'D. Still water 2 m deep in its hottest 1% of hours')):
        ax = fig.add_subplot(gs[row, col], projection='mollweide')
        drawn[key] = draw_map(ax, p, values, sea, title)
    handles = [Patch(facecolor=col, edgecolor=INK3, linewidth=0.5, label=n) for col, n in zip(RAMP, names)]
    handles.append(Patch(facecolor='white', edgecolor=OUTLINE, linewidth=0.9, label="outline: the climate model's water cells"))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.53, top(4.44)), ncol=5, fontsize=7.6, frameon=False,
               handlelength=1.6, columnspacing=1.4)

    ax = fig.add_subplot(gs_bars[0, 0])
    y = [-(n + GROUP_GAP * row[4]) for n, row in enumerate(BARS)]
    for yi, (key, stat, label, of, _) in zip(y, BARS):
        shares = np.array(p['shares'][stat][key]) * 100.0
        left = 0.0
        for k, v in enumerate(shares):
            if v > 0:
                ax.barh(yi, v, left=left, height=0.46, color=RAMP[k], edgecolor='white', linewidth=1.2)
                if v >= 7.0:
                    ax.text(left + v / 2.0, yi, f'{v:.0f}%', ha='center', va='center', fontsize=7.0,
                            color='white' if k >= 4 else INK)
                left += v
        hot = float(shares[5:].sum())
        ax.text(0.0, yi + 0.30, label, ha='left', va='bottom', fontsize=7.8, color=INK)
        ax.text(101.0, yi, f'{hot:.0f}% {of} at 36 °C and more', va='center', fontsize=7.4, color=INK)
    ax.set_yticks([])
    ax.set_xlim(0, 132)
    ax.set_ylim(min(y) - 0.45, 0.7)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("share, %  (the seas by their area, the sea coast by its length, the rain-fed lakes by their area)",
                  fontsize=7.6, color=INK2, loc='left')
    ax.tick_params(axis='x', labelsize=7, colors=INK2)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(INK3)
    ax.spines['bottom'].set_bounds(0, 100)
    sm = p['summary']
    ax.set_title(f"E. The waters by warmth class. The seas and large lakes cover {sm['seas']['area_km2'] / 1e6:.1f} million km², "
                 f"their upper 50 m averaging {sm['seas']['surface_mean_c']:.1f} °C; {sm['coast']['length_km']:,.0f} km of sea "
                 f"coast; {sm['lakes']['area_km2'] / 1e6:.2f} million km² of rain-fed lakes on the land.",
                 loc='left', fontsize=9.4, color=INK, pad=8)

    # how to read a mixed layer, from the product's own check: the lake model's calm column at the forcing's sites
    sites = list(p['checks']['sites'].values())
    column = [v['calm_column']['surface_p99_c'] for v in sites]
    above = [v['calm_column']['surface_p99_c'] - v['body_p99_c'] for v in sites]
    note = textwrap.wrap(
        f"A mixed layer has one temperature through its depth. Late in the lunar day the top few metres of calm deep water "
        f"run warmer than the body below: at the study's five sites the project's lake model puts the surface of a still "
        f"{sites[0]['calm_column']['depth_m']:.0f} m water at {whole(column)} °C in its hottest 1% of hours, "
        f"{min(above):.0f}–{max(above):.0f} K above the body there and between the still 4 m layer "
        f"({whole([v['still_4m_p99_c'] for v in sites])} °C) and the still 2 m layer ({whole([v['still_2m_p99_c'] for v in sites])} °C). "
        f"The body's figures in A, B and the shore rows hold for a surface that waves, currents and tides keep mixed into it; "
        f"where the water lies calm, read C and D.", width=205)
    if len(note) > 3:
        raise ValueError('the footnote no longer fits its three lines')
    fig.text(0.02, 0.60 / HEIGHT_IN, '\n'.join(note), fontsize=8.2, color=INK2, linespacing=1.5, va='top')

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    png = out / 'water_temperature_map.png'
    fig.savefig(png, dpi=150, facecolor='white')
    # faithfulness: the class shares counted from the drawn grids against the same shares from the product's cells
    area = grid(p, c['area_km2'])
    def from_cells(values, mask):
        k = classes_of(p, values)
        return [float(area[mask & (k == n)].sum() / area[mask].sum()) for n in range(len(RAMP))]
    worst = 0.0
    for key, values in (('body_mean', body_mean), ('body_p99', body_p99), ('still4_p99', still4_p99), ('still2_p99', still2_p99)):
        for mask in (sea, ~sea):
            worst = max(worst, max(abs(a - b) for a, b in zip(area_shares(p, drawn[key], mask), from_cells(values, mask))))
    against = {'body_mean': ('mean_c', 'seas_surface'), 'body_p99': ('p99_c', 'seas_surface'),
               'still4_p99': ('p99_c', 'seas_still_4m'), 'still2_p99': ('p99_c', 'seas_still_2m')}
    worst_product = max(abs(a - b) for key, (stat, share) in against.items()
                        for a, b in zip(area_shares(p, drawn[key], sea), p['shares'][stat][share]))
    manifest = dict(product='research/studies/littoral_community/results/water_map.json',
                    product_sha256=hashlib.sha256(PRODUCT.read_bytes()).hexdigest(), schema=p['schema'],
                    renderer='visualization/water-temperature-map/render.py',
                    renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    image=png.name, image_sha256=hashlib.sha256(png.read_bytes()).hexdigest(),
                    ramp=RAMP, classes=names, panels=list(drawn),
                    faithfulness=dict(worst_share_difference_drawn_against_cells=worst,
                                      worst_difference_seas_drawn_against_product_shares=worst_product))
    (out / 'water_temperature_map_figure.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(png, json.dumps(manifest['faithfulness']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
