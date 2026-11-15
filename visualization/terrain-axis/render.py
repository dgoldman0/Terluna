#!/usr/bin/env python3
"""The terrain beside the moisture map, as three sheets on the turned Moon: humidity × terrain, the six terrain settings
with regional close-ups, and the night's darkness with the rain's timing.

    python3 visualization/terrain-axis/render.py [--out DIR]     # PNGs and a manifest; default the research drive's folder

Reads biosphere/ecology/results/terrain_axis.json (schema terluna.biosphere.ecology-terrain-axis/1) with its per-pixel
fields, and biosphere/ecology/results/moisture_map.json for the darkness and the rain's timing. The global maps are equal-
area Mollweide rasters in the turned Moon's frame, the near side at the centre and the limbs dashed; each raster pixel
takes the 16 px/deg pixel of today's ground under it (carried with the orientation product's rotation), or the T21 cell
holding it, with no smoothing. The layout, palettes and Mollweide raster follow the Codex mockups of 7 November. The
manifest records the products' hashes and a faithfulness check: the land shares counted from the drawn rasters against
the product's. The images stay out of Git.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                  # noqa: E402
import matplotlib.patheffects as pe                              # noqa: E402
from matplotlib.colors import to_rgb                             # noqa: E402
from matplotlib.patches import Ellipse, Patch, Rectangle        # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
import sys                                                       # noqa: E402
sys.path.insert(0, str(ROOT))
from climate.gcm.turned import read_frame, to_today, to_turned  # noqa: E402
from biosphere.ecology.moisture_map import cell_index            # noqa: E402

PRODUCT = ROOT / 'biosphere' / 'ecology' / 'results' / 'terrain_axis.json'
MOISTURE = ROOT / 'biosphere' / 'ecology' / 'results' / 'moisture_map.json'
SCHEMA = 'terluna.biosphere.ecology-terrain-axis/1'
OUT = Path('/media/projectspace/terluna-research/session-tools/terrain_axis')
INK, INK2, INK3 = '#20221f', '#62665e', '#8a8984'
SEA = '#cfdbe3'
TERRAIN_COLOURS = ['#8c6999', '#d8d4b7', '#c1844f']
BASE = ['#c48a74', '#c3aa58', '#a1ab5d', '#60a16e', '#439a9e', '#4d4f9e']       # hue: moisture class; wet in indigo, clear of the lakes
PALETTE = np.array([v for h in BASE for v in (np.array(to_rgb(h)) * .65, np.array(to_rgb(h)), .52 * np.array(to_rgb(h)) + .48)])
CUE_COLOURS = ['#eee8dc', '#43865c', '#9164a8', '#153e38']
SETTING_COLOURS = ['#995b47', '#d48e6a', '#e2c29b', '#e8e2ce', '#a3b69e', '#bcaec8']
DARK_COLOURS = ['#f6f0d8', '#e0b96e', '#a0637f', '#3f2a56']        # kept clear of the seas' blue-grey
RAIN_COLOURS = ['#f4efe3', '#a9d2b5', '#3f9c8f', '#1d5a6e']
DRY_DARK_COLOURS = ['#e9e5da', '#d3a35a', '#7d86b8', '#5a3b6e']
CRATERS = (('Copernicus', 9.62, -20.08), ('Tycho', -43.31, -11.36), ('Aristarchus', 23.73, -47.49))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mollweide(lon, lat):
    phi = np.radians(np.asarray(lat, dtype=float))
    theta = phi.copy()
    for _ in range(12):
        den = 2 + 2 * np.cos(2 * theta)
        theta -= np.divide(2 * theta + np.sin(2 * theta) - np.pi * np.sin(phi), den, out=np.zeros_like(theta),
                           where=abs(den) > 1e-10)
    return 2 * np.sqrt(2) / np.pi * np.radians(lon) * np.cos(theta), np.sqrt(2) * np.sin(theta)


def raster_points(width):
    """The turned frame's latitude and longitude at each equal-area raster pixel, and which pixels lie on the globe."""
    height = width // 2
    xx = ((np.arange(width) + .5) / width * 4 * np.sqrt(2) - 2 * np.sqrt(2))[None, :]
    yy = (np.sqrt(2) - (np.arange(height) + .5) / height * 2 * np.sqrt(2))[:, None]
    theta = np.arcsin(yy / np.sqrt(2))
    phi = np.arcsin(np.clip((2 * theta + np.sin(2 * theta)) / np.pi, -1, 1))
    lam = np.pi * xx / (2 * np.sqrt(2) * np.cos(theta))
    inside = np.abs(lam) <= np.pi
    return np.degrees(np.broadcast_to(phi, inside.shape))[inside], np.degrees(lam)[inside], inside


def ground_index(lat_t, lon_t, rotation, ppd, rows, cols):
    lat0, lon0 = to_today(lat_t, lon_t, rotation)
    r = np.clip(np.floor((90.0 - lat0) * ppd).astype(int), 0, rows - 1)
    c = np.floor((lon0 % 360.0) * ppd).astype(int) % cols
    return r, c


def image(values, colours, inside, missing=None):
    rgb = np.array([to_rgb(c) if isinstance(c, str) else c for c in colours])
    out = np.ones(inside.shape + (4,))
    out[..., 3] = 0.0
    px = np.where(values >= 0, values, 0)
    col = rgb[px]
    col[values < 0] = to_rgb(SEA)
    if missing is not None:
        col[missing] = to_rgb('#ecebe7')
    out[inside, :3] = col
    out[inside, 3] = 1.0
    return out


def continuous(values, cmap, lo, hi, inside, water):
    out = np.ones(inside.shape + (4,))
    out[..., 3] = 0.0
    col = plt.get_cmap(cmap)(np.clip((values - lo) / (hi - lo), 0, 1))[:, :3]
    col[water] = to_rgb(SEA)
    out[inside, :3] = col
    out[inside, 3] = 1.0
    return out


def draw_globe(fig, box, rgba, title, marks=()):
    ax = fig.add_axes(box)
    ax.imshow(rgba, extent=(-2 * np.sqrt(2), 2 * np.sqrt(2), -np.sqrt(2), np.sqrt(2)), origin='upper',
              interpolation='nearest', rasterized=True, zorder=1)
    lat = np.linspace(-89.99, 89.99, 300)
    for lo in range(-150, 180, 30):
        ax.plot(*mollweide(lo, lat), color='white', alpha=.55, lw=.3, zorder=2)
    for la in range(-75, 90, 15):
        lon = np.linspace(-180, 180, 360)
        ax.plot(*mollweide(lon, np.full(lon.shape, la)), color='white', alpha=.55, lw=.3, zorder=2)
        xl, yl = mollweide(-180, np.array(float(la)))
        ax.text(xl - .07, yl, f'{la}°', ha='right', va='center', fontsize=6.2, color=INK2)
    for lo in (-90, 90):
        ax.plot(*mollweide(lo, lat), color=INK2, alpha=.85, lw=.65, ls=(0, (4, 3)), zorder=3)
    for lo, label in ((-90, '90° W'), (0, '0°'), (90, '90° E')):
        x, y = mollweide(lo, np.array(-2.0))
        ax.text(x, y, label, ha='center', va='top', fontsize=6.2, color=INK, zorder=4,
                path_effects=[pe.withStroke(linewidth=1.7, foreground='white')])
    for name, la, lo in marks:
        x, y = mollweide(lo, np.array(la))
        ax.add_patch(Rectangle((x - .09, y - .09), .18, .18, fill=False, edgecolor=INK, lw=.9, zorder=6))
        ax.text(x + .12, y + .1, name, fontsize=6.6, color=INK, zorder=6,
                path_effects=[pe.withStroke(linewidth=1.7, foreground='white')])
    ax.add_patch(Ellipse((0, 0), 4 * np.sqrt(2), 2 * np.sqrt(2), fill=False, edgecolor=INK, lw=.85, zorder=5))
    ax.set_xlim(-3.07, 2.86)
    ax.set_ylim(-1.46, 1.46)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=9.6, color=INK, loc='left', pad=8)
    return ax


def drawn_shares(values, inside_count, n):
    """Shares of the land counted from an equal-area raster (each pixel the same area)."""
    land = values >= 0
    return np.bincount(values[land], minlength=n) / land.sum()


def legend(fig, x, y, colours, names, ncol, size=7.4):
    fig.legend(handles=[Patch(facecolor=c, edgecolor='#797b73', lw=.4, label=n) for c, n in zip(colours, names)],
               loc='upper left', bbox_to_anchor=(x, y), ncol=ncol, fontsize=size, frameon=False, handlelength=1.5,
               columnspacing=1.2)


def note(fig, x, y, text, size=7.6, width=None):
    fig.text(x, y, text, fontsize=size, color=INK2, va='top', linespacing=1.45)


def pick_windows(p, mp, f, rotation, classes=(1, 2, 5), half=6.0):
    """For each class, a window in the turned frame centred on a robust land cell with at least six of its eight
    neighbours land cells of the same class, within 55 degrees of the equator, choosing the most rugged by the share of
    local lows and highs."""
    g, c = mp['grid'], mp['cells']
    glat, glon = np.array(g['lat']), np.array(g['lon'])
    cls = np.full((glat.size, glon.size), -1)
    rob = np.zeros(cls.shape, bool)
    cls[np.array(c['row']), np.array(c['col'])] = c['moisture']
    rob[np.array(c['row']), np.array(c['col'])] = np.array(c['robust'], bool)
    picks = []
    for k in classes:
        best = None
        for i in range(1, glat.size - 1):
            if abs(glat[i]) > 55:
                continue
            for j in range(glon.size):
                if cls[i, j] != k or not rob[i, j]:
                    continue
                block = cls[i - 1:i + 2, [(j - 1) % glon.size, j, (j + 1) % glon.size]]
                if (block == k).sum() - 1 < 6:
                    continue
                la, lo = float(glat[i]), float((glon[j] + 180) % 360 - 180)
                LA, LO = np.meshgrid(la + np.linspace(-half, half, 48), lo + np.linspace(-half, half, 48) / max(math.cos(math.radians(la)), .3), indexing='ij')
                r, cc = ground_index(LA.ravel(), LO.ravel(), rotation, 16, *f['position'].shape)
                pos = f['position'][r, cc]
                score = np.mean((pos == 0) | (pos == 2)) * np.mean(pos >= 0)
                if best is None or score > best[0]:
                    best = (score, la, lo)
        if best:
            picks.append((p['classes']['moisture'][k], best[1], best[2]))
    return picks


def window_image(f, rotation, field, colours, centre, half, ppd_out=16):
    la, lo = centre
    coslat = max(math.cos(math.radians(la)), .3)
    n = int(2 * half * ppd_out)
    lat_t = la + half - (np.arange(n) + .5) / ppd_out
    lon_t = lo - half / coslat + (np.arange(n) + .5) / ppd_out / coslat
    LA, LO = np.meshgrid(lat_t, lon_t, indexing='ij')
    r, c = ground_index(LA.ravel(), LO.ravel(), rotation, 16, *field.shape)
    v = field[r, c].reshape(LA.shape)
    rgb = np.array([to_rgb(x) if isinstance(x, str) else x for x in colours])
    img = rgb[np.maximum(v, 0)]
    img[v < 0] = to_rgb(SEA)
    return img, (lon_t[0], lon_t[-1], lat_t[-1], lat_t[0]), 1 / coslat


def humidity_sheet(p, mp, f, idx, inside, rotation, out):
    names = p['classes']['moisture']
    m, pos, cue = f['moisture'][idx], f['position'][idx], f['cue'][idx]
    water = m < 0
    joint = np.where(water, -1, m.astype(int) * 3 + pos.astype(int))
    fig = plt.figure(figsize=(14.4, 17.0), dpi=160, facecolor='white')
    run = p['producer']
    fig.text(.028, .978, "The Open Moon: humidity × terrain on the turned Moon's design climate", fontsize=17, color=INK,
             weight='bold', va='top')
    fig.text(.028, .958, f"The moisture map's classes ({run['run']}, years {run['years'][0]}–{run['years'][1]}, T21 cells "
             f"about 170 km across) under LOLA's terrain at 16 px/deg (1.9 km), with the lakes routed on the same forty "
             f"years.\nEqual-area maps in the turned Moon's frame, the near side at the centre and the limbs dashed. "
             f"{p['land']['spread_share']:.1%} of the land lies in cells the run treats as sea and takes the nearest land "
             f"cell's class.", fontsize=8.8, color=INK2, linespacing=1.5, va='top')
    wins = pick_windows(p, mp, f, rotation)
    marks = [(f'{"GHI"[k]}', la, lo) for k, (_, la, lo) in enumerate(wins)]
    draw_globe(fig, [.02, .675, .47, .245], image(joint, PALETTE, inside),
               'A. Moisture class (hue) and terrain position (shade); boxes: G–I', marks=marks)
    draw_globe(fig, [.51, .675, .47, .245], image(cue, CUE_COLOURS, inside),
               'B. Wetland cues: humid local lows and low ground near water')
    legend(fig, .53, .683, CUE_COLOURS, p['classes']['cues'], 4)
    note(fig, .53, .662, 'Humid local lows: humid or wet class, local low, slope under 2°. Low ground near water:\n'
         'within 15 km of a sea or lake, 0–100 m above its surface, under 2°, not a local high.\nThe cues mark candidate '
         'settings; no flooding or groundwater model enters.', size=7.4)
    draw_globe(fig, [.02, .385, .47, .245], image(np.where(water, -1, pos), TERRAIN_COLOURS, inside),
               'C. Terrain position against the land within 30 km')
    legend(fig, .04, .393, TERRAIN_COLOURS, ['local lows', 'intermediate', 'local highs'], 3)
    note(fig, .04, .372, 'Low or high when the height differs from the area-weighted mean of the land within 30 km\n'
         'by more than max(100 m, half the standard deviation there). Position is relative to the\nsurroundings, apart '
         'from height above sea level.', size=7.4)
    tw = f['wetness'][idx].astype(float)
    lo, hi = p['wetness']['land_p10'], p['wetness']['land_p90'] + 2.0
    draw_globe(fig, [.51, .385, .47, .245], continuous(tw, 'YlGnBu', lo, hi, inside, water),
               'D. Topographic wetness index: where water gathers')
    cax = fig.add_axes([.555, .388, .3, .007])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap='YlGnBu', norm=plt.Normalize(lo, hi)), cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=6.6, colors=INK2)
    cb.outline.set_visible(False)
    note(fig, .53, .370, f"ln(a / tan b): a, the area drained through a pixel per unit width in the forty years' "
         f"routing;\nb, its slope (3.8 km baseline, floored at 0.1°). Land median {p['wetness']['land_median']:.1f}; the "
         f"top tenth, above {p['wetness']['land_p90']:.1f},\nmarks valley floors, footslopes and shores wetter than their "
         f"climate.", size=7.4)
    for k, (name, la, lo_) in enumerate(wins):
        ax = fig.add_axes([.05 + k * .315, .175, .27, .15])
        img, ext, asp = window_image(f, rotation, np.where(f['moisture'] < 0, -1, f['moisture'].astype(int) * 3 + np.maximum(f['position'], 0)), PALETTE, (la, lo_), 6.0)
        ax.imshow(img, extent=ext, interpolation='nearest', aspect=asp)
        ax.tick_params(labelsize=6.2, colors=INK2)
        for s_ in ax.spines.values():
            s_.set_color(INK3)
        ax.set_title(f'{"GHI"[k]}. {name} ground at {la:.0f}°, {lo_:.0f}° (turned), 12° across', fontsize=8.4, color=INK,
                     loc='left')

    shares = np.array(p['shares']['moisture_by_position']) * 100
    fig.text(.028, .135, 'E. Colour key and shares of the land', fontsize=9.6, color=INK)
    ax = fig.add_axes([.085, .018, .28, .108])
    ax.imshow(PALETTE.reshape(6, 3, 3), aspect='auto', interpolation='nearest', origin='upper')
    for h in range(6):
        for t_ in range(3):
            lum = PALETTE[h * 3 + t_] @ np.array([.299, .587, .114])
            ax.text(t_, h, f'{shares[h, t_]:.1f}%', ha='center', va='center', color='white' if lum < .53 else INK, fontsize=8.2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(['local low', 'intermediate', 'local high'], fontsize=7.4, color=INK)
    ax.set_yticks(range(6))
    ax.set_yticklabels(names, fontsize=7.4, color=INK)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-.5, 3), minor=True)
    ax.set_yticks(np.arange(-.5, 6), minor=True)
    ax.grid(which='minor', color='white', lw=1.2)
    ax.tick_params(which='minor', length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ws = p['wetness']['by_moisture']
    rows = '\n'.join(f"{n}: {100 * sum(p['shares']['moisture_by_position'][k]):.1f}% of the land; top-tenth wetness on "
                     f"{100 * ws[k]['share_above_land_p90']:.1f}% of it" for k, n in enumerate(names) if ws[k])
    fig.text(.43, .135, 'F. By moisture class', fontsize=9.6, color=INK)
    note(fig, .43, .122, rows, size=7.6)
    note(fig, .70, .122, f"Robust moisture class (the moisture map's flag)\non {p['land']['robust_share']:.0%} of the land. "
         f"Thresholds are proposals from\nthe 7 November mockups; the register has\nadopted none of them as a biome "
         f"boundary.\nLand {p['land']['area_km2'] / 1e6:.2f} million km², {p['land']['share_of_moon']:.1%} of the Moon.",
         size=7.6)
    fig.text(.028, .004, 'Sources: biosphere/ecology/results/terrain_axis.json and moisture_map.json; LOLA above the GRAIL '
             'geoid; the drainage of geography/turned_corrected_y10-49.', fontsize=7.0, color=INK2)
    path = out / 'humidity_terrain.png'
    fig.savefig(path, dpi=160, facecolor='white')
    plt.close(fig)
    return path, dict(joint=drawn_shares(joint, inside.sum(), 18).tolist(), windows=wins)


def closeup(ax, f, rotation, centre, half_deg, ppd_out=24):
    """A box in the turned frame at about the source resolution, each pixel the 16 px/deg pixel of today's ground under
    it."""
    la, lo = centre
    n = int(2 * half_deg * ppd_out)
    lat_t = la + half_deg - (np.arange(n) + .5) / ppd_out
    coslat = max(math.cos(math.radians(la)), .2)
    lon_t = lo - half_deg / coslat + (np.arange(int(n / 1)) + .5) / ppd_out / coslat
    LA, LO = np.meshgrid(lat_t, lon_t, indexing='ij')
    r, c = ground_index(LA.ravel(), LO.ravel(), rotation, 16, *f['setting'].shape)
    s = f['setting'][r, c].reshape(LA.shape)
    rgb = np.array([to_rgb(x) for x in SETTING_COLOURS])
    img = rgb[np.maximum(s, 0)]
    img[s < 0] = to_rgb(SEA)
    ax.imshow(img, extent=(lon_t[0], lon_t[-1], lat_t[-1], lat_t[0]), interpolation='nearest', aspect=1 / coslat)
    ax.tick_params(labelsize=6.4, colors=INK2)
    for s_ in ax.spines.values():
        s_.set_color(INK3)
    return s


def settings_sheet(p, f, idx, inside, rotation, out):
    names = [s['name'] for s in p['classes']['settings']]
    s = f['setting'][idx]
    fig = plt.figure(figsize=(14.4, 12.4), dpi=160, facecolor='white')
    fig.text(.028, .978, 'The Open Moon: terrain settings on the turned Moon', fontsize=17, color=INK, weight='bold',
             va='top')
    fig.text(.028, .950, "The landscape mockup's six settings, computed on LOLA's terrain at 16 px/deg (1.9 km) in place of "
             "an invented crater. Each pixel takes the first that applies:\nbasin floor (at least 10 m below the spill "
             "level of the closed depression that holds it, from the forty years' routing), steep (over 10°), crest "
             "(local high),\nhollow (local low), plain (under 2°), gentle slope. Slopes on a 3.8 km baseline; steeper "
             "ground within a pixel is unresolved.", fontsize=8.6, color=INK2, linespacing=1.5, va='top')
    marks = []
    for name, la, lo in CRATERS:
        tl, tn = to_turned(np.array(la), np.array(lo), rotation)
        marks.append((name, float(tl), float(tn)))
    draw_globe(fig, [.02, .515, .62, .37], image(np.where(s >= 0, s, -1), SETTING_COLOURS, inside),
               'A. Six terrain settings; boxes: the close-ups', marks=marks)
    legend(fig, .66, .86, SETTING_COLOURS, names, 1, size=8.2)
    sh = p['shares']['settings']
    note(fig, .66, .72, '\n'.join(f'{n}: {100 * x:.1f}% of the land' for n, x in zip(names, sh)), size=8.2)
    for k, (name, la, lo) in enumerate(marks):
        ax = fig.add_axes([.05 + k * .315, .265, .27, .2])
        closeup(ax, f, rotation, (la, lo), 1.6)
        ax.set_title(f'{"BCD"[k]}. {name}, at {la:.1f}°, {lo:.1f}° on the turned Moon', fontsize=8.6,
                     color=INK, loc='left')
    ax = fig.add_axes([.17, .05, .62, .155])
    within = np.array(p['shares']['settings_within_moisture'])
    mnames = p['classes']['moisture']
    left = np.zeros(6)
    for j in range(6):
        ax.barh(np.arange(6)[::-1], within[:, j] * 100, left=left * 100, color=SETTING_COLOURS[j], edgecolor='white', lw=.5)
        left += within[:, j]
    ax.set_yticks(np.arange(6)[::-1])
    ax.set_yticklabels(mnames, fontsize=7.8, color=INK)
    ax.set_xlim(0, 100)
    ax.tick_params(axis='x', labelsize=7, colors=INK2)
    ax.tick_params(axis='y', length=0)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_xlabel('share of each moisture class\'s land, %', fontsize=7.6, color=INK2)
    ax.set_title('E. The settings within each moisture class', fontsize=9.6, color=INK, loc='left')
    path = out / 'terrain_settings.png'
    fig.savefig(path, dpi=160, facecolor='white')
    plt.close(fig)
    return path, dict(settings=drawn_shares(np.where(s >= 0, s, -1), inside.sum(), 6).tolist())


def night_sheet(mp, lat_t, lon_t, inside, out):
    g, c = mp['grid'], mp['cells']
    glat, glon = np.array(g['lat']), np.array(g['lon'])
    grid = np.full((glat.size, glon.size), -1)
    cr, cc = np.array(c['row']), np.array(c['col'])
    grid[cr, cc] = np.arange(cr.size)
    r, k = cell_index(lat_t, lon_t, glat, glon)
    cell = grid[r, k]
    sea = cell < 0
    dark = np.array(c['plant_dark_hours'], float)
    rain = np.array(c['night_rain_share'], float)
    ai = np.array(c['aridity'], float)
    area = np.array(c['weight'], float)                         # the moisture map's land weights

    def classed(values, edges):
        v = np.digitize(values, edges)
        out_ = np.where(sea, -1, v[np.maximum(cell, 0)])
        return out_
    dark_cls = classed(dark, [160, 220, 240])
    rain_cls = classed(rain, [0.10, 0.25, 0.40])
    dd = (ai < 0.65).astype(int) + 2 * (dark > 220).astype(int)
    dd_px = np.where(sea, -1, dd[np.maximum(cell, 0)])
    fig = plt.figure(figsize=(14.4, 10.0), dpi=160, facecolor='white')
    run = mp['producer']
    fig.text(.028, .975, 'The Open Moon: the night\'s darkness and the rain\'s timing on the turned Moon', fontsize=17,
             color=INK, weight='bold', va='top')
    fig.text(.028, .938, f"The two panels of the ecological-strategy mockup that the data support, from the moisture map "
             f"({run['run']}, years {run['years'][0]}–{run['years'][1]}), drawn as T21 cells, about 170 km across.\n"
             "The plant-strategy classes themselves need a vegetation model the project does not have.", fontsize=8.8,
             color=INK2, linespacing=1.5, va='top')
    draw_globe(fig, [.02, .515, .47, .37], image(dark_cls, DARK_COLOURS, inside),
               'A. Hours per lunar day with too little light for photosynthesis')
    legend(fig, .04, .522, DARK_COLOURS, ['under 160 h', '160–220 h', '220–240 h', '240 h and more'], 4)
    draw_globe(fig, [.51, .515, .47, .37], image(rain_cls, RAIN_COLOURS, inside),
               'B. Share of the rain that falls at night')
    legend(fig, .53, .522, RAIN_COLOURS, ['under 10%', '10–25%', '25–40%', '40% and more'], 4)
    draw_globe(fig, [.02, .085, .47, .37], image(dd_px, DRY_DARK_COLOURS, inside),
               'C. Drylands and long darkness together')
    legend(fig, .04, .092, DRY_DARK_COLOURS, ['neither', 'dryland (aridity under 0.65)', 'dark over 220 h',
                                              'both'], 2)
    w = area / area.sum()
    lines = [f"Darkness: under 160 h on {100 * w[dark < 160].sum():.0f}% of the land, 220 h and more on "
             f"{100 * w[dark >= 220].sum():.0f}%.",
             f"Night rain: 25% or more of the rain at night on {100 * w[rain >= .25].sum():.0f}% of the land, "
             f"40% or more on {100 * w[rain >= .40].sum():.0f}%.",
             f"Drylands under long darkness (C, both): {100 * w[(ai < .65) & (dark > 220)].sum():.0f}% of the land.",
             '',
             'The darkness is the moisture map\'s clear-sky plant-dark hours, which depend on latitude alone: the plant',
             'model\'s bands at 0°, 30° and 60°, interpolated, and poleward of 60° the polar light study\'s hours,',
             'which end at 66° N and 67° S.',
             'Night rain falls while plants cannot use it, so the share that falls at night is the',
             'lunar counterpart of a dry season. Both are descriptors; plant stress is not diagnosed.']
    fig.text(.53, .40, 'D. Shares of the land (the moisture map\'s land weights)', fontsize=9.6, color=INK)
    note(fig, .53, .38, '\n'.join(lines), size=8.2)
    path = out / 'darkness_rain_timing.png'
    fig.savefig(path, dpi=160, facecolor='white')
    plt.close(fig)
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(OUT))
    ap.add_argument('--width', type=int, default=1800)
    args = ap.parse_args(argv)
    p = json.loads(PRODUCT.read_text())
    if p['schema'] != SCHEMA:
        raise ValueError(f"expected {SCHEMA}, found {p['schema']}")
    if digest(p['fields']['path']) != p['fields']['sha256']:
        raise ValueError('the per-pixel fields differ from the product\'s record')
    mp = json.loads(MOISTURE.read_text())
    rotation = read_frame(ROOT / p['producer']['frame']['path'])['rotation']
    f = dict(np.load(p['fields']['path']))
    lat_t, lon_t, inside = raster_points(args.width)
    rows, cols = f['setting'].shape
    r, c = ground_index(lat_t, lon_t, rotation, 16, rows, cols)
    idx = (r, c)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    a, sa = humidity_sheet(p, mp, f, idx, inside, rotation, out)
    b, sb = settings_sheet(p, f, idx, inside, rotation, out)
    n = night_sheet(mp, lat_t, lon_t, inside, out)
    joint = np.array(p['shares']['moisture_by_position']).ravel()
    faithful = dict(joint_max_difference=float(np.max(np.abs(np.array(sa['joint']) - joint))),
                    settings_max_difference=float(np.max(np.abs(np.array(sb['settings']) - np.array(p['shares']['settings'])))))
    manifest = dict(product=str(PRODUCT.relative_to(ROOT)), product_sha256=digest(PRODUCT),
                    moisture_sha256=digest(MOISTURE), fields_sha256=p['fields']['sha256'],
                    images={x.name: digest(x) for x in (a, b, n)}, raster_width=args.width,
                    humidity_windows=sa['windows'],
                    faithfulness=dict(faithful, note='land shares counted from the equal-area rasters against the '
                                      'product\'s; the raster samples about one source pixel in fifty'))
    (out / 'terrain_axis_figures.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps(manifest['faithfulness']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
