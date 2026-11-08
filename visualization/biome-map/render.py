#!/usr/bin/env python3
"""The Open Moon's first biome map as one sheet: the biomes as run with their uncertain cells hatched, where the
uncertainty comes from, the night's light, and the land's shares with their robust parts.

    python3 visualization/biome-map/render.py [--out DIR]     # PNG and a manifest; default the research drive's folder

Reads biosphere/ecology/results/biome_map.json (schema terluna.biosphere.ecology-biome-map/1) and draws its cells as
they stand: each T21 land cell a box of its biome's colour, sea cells grey-blue, no smoothing or downscaling. The
manifest records the product's hash and a faithfulness check, the land shares counted from the drawn grid against the
product's. Colours: the biomes' ramp runs from sand to dark green in falling lightness; the uncertainty sources
take the first three slots of the reference categorical palette; the night ramp runs from the dark belt to the bright
nights. The images stay out of Git.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
import matplotlib.patheffects as patheffects                     # noqa: E402
from matplotlib.collections import PatchCollection               # noqa: E402
from matplotlib.patches import Patch, Rectangle                  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / 'biosphere' / 'ecology' / 'results' / 'biome_map.json'
SCHEMA = 'terluna.biosphere.ecology-biome-map/1'
OUT = Path('/media/projectspace/terluna-research/session-tools/agents/biome_map')
BIOME_COLOURS = ['#faf5ea', '#ebc98c', '#aeb041', '#6a9c50', '#2a6a4a', '#0c3528']   # neighbours 15.9-22 apart (CIEDE2000)
POLAR_DOTS_DEG = 83.0              # the outermost rows' cells crowd into a line on the projection: counted there, not dotted
SEA = '#cfdbe3'
UNCERTAIN = '#d8d6d0'
SOURCE_COLOURS = {'robust': '#e4e3df', 'the gap alone': '#eb6834', 'the years alone': '#2a78d6', 'both together': '#1baf7a'}
NIGHT_COLOURS = ['#1d2740', '#56678f', '#a9c3ea', '#e3b264']
INK, INK2, INK3 = '#1f1f1d', '#52514e', '#8a8984'


def grid_array(product, values, fill=np.nan):
    g = product['grid']
    a = np.full((len(g['lat']), len(g['lon'])), fill, dtype=float)
    a[np.array(product['cells']['row']), np.array(product['cells']['col'])] = values
    return a


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


def draw_map(ax, product, values, colours, title, hatch=None, dots=None):
    """One field on an equal-area (Mollweide) map, the near side at the centre and the limbs dashed."""
    a = grid_array(product, values)
    a2, xedges = shifted(product, a)
    yedges = np.array(product['grid']['lat_edges'])
    xr, yr = np.radians(np.clip(xedges, -180.0, 180.0)), np.radians(yedges)
    land = np.isfinite(a2)
    ax.pcolormesh(xr, yr, np.where(land, np.nan, 0.0), cmap=matplotlib.colors.ListedColormap([SEA]), shading='flat',
                  edgecolors='none')
    cmap = matplotlib.colors.ListedColormap(colours)
    ax.pcolormesh(xr, yr, np.ma.masked_invalid(a2), cmap=cmap, vmin=-0.5, vmax=len(colours) - 0.5, shading='flat',
                  edgecolors='none')
    if hatch is not None:
        h2, _ = shifted(product, grid_array(product, hatch.astype(float), 0.0))
        boxes = [Rectangle((xr[j], yr[i + 1]), xr[j + 1] - xr[j], yr[i] - yr[i + 1]) for i, j in zip(*np.nonzero(h2 > 0.5))]
        ax.add_collection(PatchCollection(boxes, facecolor='none', edgecolor=INK, hatch='//////', linewidth=0.0,
                                          alpha=0.6))
    if dots is not None:
        d2, _ = shifted(product, grid_array(product, dots.astype(float), 0.0))
        i, j = np.nonzero(d2 > 0.5)
        y = 0.5 * (yr[i] + yr[i + 1])
        keep = np.abs(y) < np.radians(POLAR_DOTS_DEG)
        ax.scatter(0.5 * (xr[j[keep]] + xr[j[keep] + 1]), y[keep], s=5, c=INK, edgecolors='white', linewidths=0.35,
                   zorder=5)
        flat, lat = grid_array(product, dots.astype(float), 0.0) > 0.5, np.array(product['grid']['lat'])
        north, south = int(flat[lat >= POLAR_DOTS_DEG].sum()), int(flat[lat <= -POLAR_DOTS_DEG].sum())
        if north or south:
            ax.text(1.0, 1.0, f'poleward of {POLAR_DOTS_DEG:.0f}°, not dotted: {north} N, {south} S',
                    transform=ax.transAxes, ha='right', va='top', fontsize=6, color=INK2)
    for x in (-90.0, 90.0):
        ax.plot(np.radians([x, x]), np.radians([-89.9, 89.9]), color=INK2, lw=0.7, ls=(0, (4, 3)))
    ax.grid(True, color='white', lw=0.35, alpha=0.7)
    ax.set_xticklabels([])
    for lon, name in ((-90.0, '90° W'), (0.0, '0°'), (90.0, '90° E')):
        ax.text(np.radians(lon), np.radians(-3.0), name, ha='center', va='top', fontsize=5.8, color=INK, zorder=6,
                path_effects=[patheffects.withStroke(linewidth=1.8, foreground='white')])
    ax.tick_params(axis='y', labelsize=6, colors=INK2)
    ax.set_title(title, loc='left', fontsize=9.2, color=INK, pad=6)
    return a


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(OUT))
    args = ap.parse_args(argv)
    p = json.loads(PRODUCT.read_text())
    if p['schema'] != SCHEMA:
        raise ValueError(f"expected {SCHEMA}, found {p['schema']}")
    c = {k: np.array(v) for k, v in p['cells'].items() if not isinstance(v[0], list)}
    names = [x['name'] for x in p['biomes']]
    robust = c['robust'].astype(bool)
    store_only = robust & ~c['robust_store'].astype(bool)
    soda = (c['closed_lake_fraction'].astype(float) > 0) & (c['biome'] <= 3)

    fig = plt.figure(figsize=(13.5, 12.6), dpi=150, facecolor='white')
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.82], hspace=0.30, wspace=0.06,
                          left=0.155, right=0.985, top=0.905, bottom=0.05)
    run = p['producer']
    fig.text(0.02, 0.968, 'The Open Moon\'s biomes on the design climate: a first map', fontsize=14,
             color=INK, weight='bold')
    fig.text(0.02, 0.928, f"{run['run']}, years {run['years'][0]}–{run['years'][1]}; T21 cells about 170 km across, "
             f"drawn as computed on an equal-area map, the near side at the centre and the limbs dashed.\nProvisional: "
             f"remade if the convection work changes the design climate. The rain scaled by {p['gap']['scale']:.2f} toward CM1 is a "
             f"sensitivity of the classes from three boxes, never a corrected climate.", fontsize=8.6, color=INK2,
             linespacing=1.5)
    biome_handles = [Patch(facecolor=col, edgecolor=INK3, linewidth=0.4, label=n)
                         for col, n in zip(BIOME_COLOURS, names)]
    sea_handle = Patch(facecolor=SEA, edgecolor='none', label='seas and large lakes (GCM water cells)')

    ax = fig.add_subplot(gs[0, 0], projection='mollweide')
    drawn = draw_map(ax, p, c['biome'].astype(float), BIOME_COLOURS,
                     'A. Biomes as run; dots: closed lakes in dry cells', dots=soda)
    axb = fig.add_subplot(gs[0, 1], projection='mollweide')
    core = np.where(robust, c['biome'].astype(float), 6.0)
    draw_map(axb, p, core, BIOME_COLOURS + [UNCERTAIN], 'B. The robust core: grey where the class is uncertain')
    fig.legend(handles=biome_handles + [Patch(facecolor=UNCERTAIN, edgecolor=INK3, linewidth=0.4,
                                                  label='uncertain (B)'), sea_handle],
               loc='upper center', bbox_to_anchor=(0.57, 0.648), ncol=4, fontsize=7.6, frameon=False,
               handlelength=1.6, columnspacing=1.4)

    ax2 = fig.add_subplot(gs[1, 0], projection='mollweide')
    src = np.where(robust, 0, np.where(~c['robust_gap'].astype(bool) & ~c['robust_years'].astype(bool), 3,
                   np.where(~c['robust_gap'].astype(bool), 1, np.where(~c['robust_years'].astype(bool), 2, 3))))
    draw_map(ax2, p, src.astype(float), list(SOURCE_COLOURS.values()),
             'C. What makes a cell uncertain; dots: robust cells a hollows\' store run moves', dots=store_only)
    ax2.legend(handles=[Patch(facecolor=v, edgecolor=INK3, linewidth=0.4, label=k) for k, v in SOURCE_COLOURS.items()],
               loc='upper center', bbox_to_anchor=(0.5, -0.03), ncol=2, fontsize=7.4, frameon=False)
    ax3 = fig.add_subplot(gs[1, 1], projection='mollweide')
    draw_map(ax3, p, c['night_regime'].astype(float), NIGHT_COLOURS,
             'D. The night\'s darkest natural light (clear sky, mean geometry)')
    labels = ['dark, under 0.01 lux', 'dim, 0.01–1 lux', 'Earthlit, 1 lux and more', 'twilit, 1 lux and more']
    ax3.legend(handles=[Patch(facecolor=v, edgecolor=INK3, linewidth=0.4, label=k) for v, k in zip(NIGHT_COLOURS, labels)],
               loc='upper center', bbox_to_anchor=(0.5, -0.03), ncol=2, fontsize=7.4, frameon=False)

    ax4 = fig.add_subplot(gs[2, :])
    rows = [p['shares']['by_biome'][x['key']] for x in p['biomes']]
    y = np.arange(len(rows))[::-1]
    for yi, row, col in zip(y, rows, BIOME_COLOURS):
        ax4.barh(yi, row['robust_share'] * 100, height=0.62, color=col, edgecolor=INK3, linewidth=0.4)
        ax4.barh(yi, (row['share'] - row['robust_share']) * 100, left=row['robust_share'] * 100, height=0.62,
                 color=col, edgecolor=INK3, linewidth=0.4, hatch='////', alpha=0.85)
        lo, hi = row['share_years_95']
        ax4.errorbar(row['share'] * 100, yi, xerr=[[row['share'] * 100 - lo * 100], [hi * 100 - row['share'] * 100]],
                     fmt='none', ecolor=INK, elinewidth=0.9, capsize=2.5)
        ax4.plot([row['share_at_gap_scale'] * 100] * 2, [yi - 0.36, yi + 0.36], color='#eb6834', lw=2.0)
        ax4.text(max(row['share'], row['share_at_gap_scale'], hi) * 100 + 0.8, yi,
                 f"{row['share'] * 100:.1f}% as run (years' 95%: {lo * 100:.1f}–{hi * 100:.1f}), "
                 f"{row['robust_share'] * 100:.1f}% robust; {row['share_at_gap_scale'] * 100:.1f}% with the rain × "
                 f"{p['gap']['scale']:.2f}", va='center', fontsize=7.6, color=INK)
    ax4.set_yticks(y)
    ax4.set_yticklabels([n.replace(' over ', '\nover ') for n in names], fontsize=8, color=INK, linespacing=1.1)
    ax4.set_xlabel('share of the land, %   (solid: robust cells; hatched: uncertain cells; black bar: the years\' 95% '
                   'range of the share; orange tick: the share with the rain scaled toward CM1)', fontsize=7.6, color=INK2)
    ax4.tick_params(axis='x', labelsize=7, colors=INK2)
    ax4.tick_params(axis='y', length=0)
    for side in ('top', 'right'):
        ax4.spines[side].set_visible(False)
    ax4.spines['left'].set_color(INK3)
    ax4.spines['bottom'].set_color(INK3)
    ax4.set_xlim(0, max(r['share'] for r in rows) * 100 * 2.1)
    ax4.grid(axis='x', color='#ecebe7', lw=0.6)
    ax4.set_axisbelow(True)
    sh = p['shares']
    ax4.set_title(f"E. Shares of the land. Robust {sh['robust'] * 100:.0f}% (the gap alone {sh['robust_gap'] * 100:.0f}%, "
                  f"the years alone {sh['robust_years'] * 100:.0f}%; with the hollows' store {sh['robust_all'] * 100:.0f}%). "
                  f"Drylands {p['mix']['drylands_share'] * 100:.0f}% of the land, Earth's {p['mix']['earth_drylands']['drylands'] * 100:.0f}%.",
                  loc='left', fontsize=9.4, color=INK, pad=5)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    png = out / 'biome_map.png'
    fig.savefig(png, dpi=150, facecolor='white')
    # faithfulness: the shares counted from the drawn grid against the product's
    lat = np.array(p['grid']['lat'])
    w = np.cos(np.radians(lat))[:, None] * np.ones(drawn.shape)
    land = np.isfinite(drawn)
    drawn_shares = {n: float((w * (drawn == k)).sum() / (w * land).sum()) for k, n in enumerate(names)}
    worst = max(abs(drawn_shares[n] - r['share']) for n, r in zip(names, rows))
    manifest = dict(product='biosphere/ecology/results/biome_map.json',
                    product_sha256=hashlib.sha256(PRODUCT.read_bytes()).hexdigest(),
                    renderer='visualization/biome-map/render.py',
                    renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    image=str(png), image_sha256=hashlib.sha256(png.read_bytes()).hexdigest(),
                    faithfulness=dict(drawn_land_cells=int(land.sum()), product_land_cells=p['land']['cells'],
                                      worst_share_difference=worst, passes=bool(worst < 1e-3 and land.sum() == p['land']['cells'])))
    (out / 'biome_map_figure.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps(manifest['faithfulness']), png)
    return 0 if manifest['faithfulness']['passes'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
