#!/usr/bin/env python3
"""Map sheets and a 3D globe of the Open Moon atlas: near and far side, a global sheet, the polar regions,
and an interactive globe page with an appearance mode and a map mode.

    python visualization/atlas/render.py      # writes visualization/atlas/out/ (kept out of Git)

Reads the geography atlas product (geography/results/atlas.json and its grid in
geography/products) and the conservation study's register and targets. Water is
drawn in depth bands of one blue hue, land as neutral shaded relief. The globe's
appearance mode also reads the illumination domain's engine atmosphere and sky atlas
and the climate domain's GCM climatology (see appearance.py). The manifest records
product hashes and faithfulness checks: the water share of the rendered near-side
disk against the product's Earth-facing disk share, and the water share of each
globe texture against the product's.
"""
from __future__ import annotations
import base64
import csv
import hashlib
import io
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from geography import atlas as ga
from shared.constants import MOON_RADIUS
import appearance

HERE = Path(__file__).resolve().parent
OUT = HERE / 'out'
STUDY = ROOT / 'research' / 'studies' / 'conservation' / 'results'
SURF, INK, INK2 = '#fcfcfb', '#0b0b0b', '#52514e'
HERITAGE, SCIENCE = '#eb6834', INK
BANDS = [(250, '#b7d3f6'), (500, '#86b6ef'), (1000, '#5598e7'), (2000, '#2a78d6'), (4000, '#1c5cab'), (1e9, '#104281')]
BAND_LABELS = ['0–250 m', '250–500 m', '500–1,000 m', '1–2 km', '2–4 km', 'over 4 km']
HALO = [pe.withStroke(linewidth=2.6, foreground='white')]
HALO_DARK = [pe.withStroke(linewidth=2.2, foreground='#0d366b')]
SEA_LABELS = [('Oceanus\nProcellarum', 14, -52), ('Mare\nImbrium', 36, -17), ('Mare\nSerenitatis', 30.5, 21),
              ('Mare\nCrisium', 16, 55), ('Mare\nNectaris', -17.5, 38), ('Mare\nFecunditatis', -6, 52),
              ('Mare\nNubium', -21, -16.6), ('Mare\nHumorum', -24.4, -38.6), ('Mare Frigoris', 57, -5),
              ('Tranquillitatis\narm', 3.5, 26.5), ('Mare Cognitum', -11, -23), ('South Pole–Aitken Sea', -50, -175),
              ('Mare\nIngenii', -33.5, 165), ('Mare\nMoscoviense', 27, 148), ('Mare\nOrientale', -19.5, -94.5),
              ('Smythii–Marginis\nSea', 6, 86), ('Mare\nHumboldtianum', 57, 83)]
LAND_LABELS = [('Tranquillitatis\nplateau', 11, 34), ('Vaporum', 13.5, 4.5), ('Insularum', 7, -29.5)]
SITE_LABELS = {'apollo-11': ('Apollo 11', (5, -9)), 'apollo-12': ('Apollo 12', (-38, 4)), 'apollo-14': ('Apollo 14', (5, 3)),
               'apollo-15': ('Apollo 15', (-8, -12)), 'apollo-16': ('Apollo 16', (5, 3)), 'apollo-17': ('Apollo 17', (5, 3)),
               'luna-2': ('Luna 2', (-30, 3)), 'luna-9': ('Luna 9', (4, 3)), 'surveyor-1': ('Surveyor 1', (-44, -3)),
               'lunokhod-1': ('Lunokhod 1', (5, 3)), 'lunokhod-2': ('Lunokhod 2', (5, 3)), 'luna-16': ('Luna 16', (-31, -3)),
               'luna-24': ('Luna 24', (4, -9)), 'blue-ghost-1': ('Blue Ghost 1', (5, -8)), 'beresheet': ('Beresheet', (5, 3)),
               'change-3': ("Chang'e 3", (5, 3)), 'change-4': ("Chang'e 4", (-44, 5)), 'change-5': ("Chang'e 5", (5, 3)),
               'change-6': ("Chang'e 6", (5, 3)), 'chandrayaan-3': ('Chandrayaan-3', (5, 3)), 'slim': ('SLIM', (-24, -8)),
               'im-1': ('IM-1', (5, 3)), 'lunar-prospector': ("Shoemaker's ashes", (12, -12)), 'lcross': ('LCROSS impact', (-20, 8)),
               'hakuto-r-m1': ('Hakuto-R', (5, 3)), 'luna-25': ('Luna 25', (5, 3)), 'resilience': ('Resilience', (5, 3))}
GLOBE_LABEL_MAX_LAT = 78.0   # sites nearer the poles are labelled on the polar sheet
TARGET_LABELS = {'Ina', 'Reiner Gamma', 'Aristarchus', 'Sulpicius Gallus', 'Mons Rumker', 'Von Karman', 'Schrodinger',
                 'Tsiolkovskiy', 'Copernicus', 'Tycho'}
POLAR_TARGET_OFFSETS = {'Shackleton': (6, -12), 'de Gerlache': (-100, 4), 'Sverdrup': (-50, -12), 'Haworth': (-46, 6),
                        'Shoemaker': (10, 10), 'Faustini': (8, -13), 'Nobile': (6, -4), 'Cabeus': (-112, -12),
                        'Hermite': (-42, 4), 'Whipple': (6, 4), 'Peary': (6, -4), 'Rozhdestvenskiy': (6, 4)}
DISPLAY = {'Rumker': 'Rümker', 'Karman': 'Kármán', 'Schrodinger': 'Schrödinger'}
LAND_PLACEHOLDER = '#b9b4aa'   # neutral land until the biosphere and climate supply surface cover
NORMAL_STRENGTH = 10.0         # relief is baked into the normal map at this exaggeration; the page rescales it


def hexrgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def shaded(height, lat, ppd):
    """Neutral hillshade for land (sun from the north-west, 40 degrees up) and depth bands for water."""
    dy = 1737.4e3 * np.radians(1 / ppd)
    dx = dy * np.maximum(np.cos(np.radians(lat)), 1e-3)[:, None]
    gy, gx = np.gradient(height)
    sx, sy = gx / dx, -gy / dy
    slope, aspect = np.arctan(np.hypot(sx, sy)), np.arctan2(-sx, sy)
    shade = np.clip(np.sin(np.radians(40)) * np.cos(slope) + np.cos(np.radians(40)) * np.sin(slope) * np.cos(np.radians(315) - aspect), 0, 1)
    return shade


def composite(height, level, shade):
    rgb = np.repeat((0.74 + 0.22 * shade)[..., None], 3, axis=2)
    depth = level - height
    lo = 0.0
    for hi, col in BANDS:
        m = (depth > 0) & (depth >= lo) & (depth < hi)
        rgb[m] = hexrgb(col) * (0.93 + 0.07 * shade[m, None])
        lo = hi
    return rgb


def sample_index(lat, lon, ppd, shape):
    i = np.clip(np.round((90 - lat) * ppd - 0.5).astype(int), 0, shape[0] - 1)
    j = np.round((lon % 360) * ppd - 0.5).astype(int) % shape[1]
    return i, j


def ortho(rgb, water, lon0, ppd, n=900):
    x, y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(1, -1, n))
    r2 = x ** 2 + y ** 2
    inside = r2 <= 1
    z = np.sqrt(np.clip(1 - r2, 0, 1))
    la = np.degrees(np.arcsin(np.clip(y, -1, 1)))
    lo = lon0 + np.degrees(np.arctan2(x, z))
    img = np.ones((n, n, 3)) * hexrgb(SURF)
    i, j = sample_index(la[inside], lo[inside], ppd, rgb.shape)
    img[inside] = rgb[i, j]
    img[inside & (r2 > 0.994)] = hexrgb('#8f8e8a')
    return img, float(water[i, j].mean())


def project(lon0, la, lo):
    p, l, l0 = np.radians(la), np.radians(lo), np.radians(lon0)
    return np.cos(p) * np.sin(l - l0), np.sin(p), np.cos(p) * np.cos(l - l0) > 0.05


def markers(ax, fx, sites, targets, small=False):
    s0 = 26 if small else 34
    for st in sites:
        x, y, vis = fx(st['lat'], st['lon'])
        if not vis:
            continue
        ax.scatter(x, y, s=s0, marker='o' if st['status'] == 'submerged' else '^', c=HERITAGE, edgecolors='white', linewidths=0.9, zorder=6)
        if st['id'] in SITE_LABELS and abs(st['lat']) < GLOBE_LABEL_MAX_LAT:
            text, off = SITE_LABELS[st['id']]
            ax.annotate(text, (x, y), xytext=off, textcoords='offset points', fontsize=6.3, color=INK, path_effects=HALO, zorder=7)
    for t in targets:
        x, y, vis = fx(float(t['lat']), float(t['lon']))
        if not vis:
            continue
        ax.scatter(x, y, s=s0 * 0.55, marker='s', c=SCIENCE, edgecolors='white', linewidths=0.8, zorder=5)
        if t['name'] in TARGET_LABELS:
            name = t['name']
            for k, v in DISPLAY.items():
                name = name.replace(k, v)
            off = (-50, -9) if t['name'] == 'Von Karman' else (4, -8)
            ax.annotate(name, (x, y), xytext=off, textcoords='offset points', fontsize=6, color=INK2, style='italic',
                        path_effects=HALO, zorder=7)


def legend(fig, y, science_label='Scientific sampling target'):
    handles = [Patch(facecolor=c, edgecolor='none', label=l) for (_, c), l in zip(BANDS, BAND_LABELS)]
    handles += [Line2D([], [], marker='o', ls='', color=HERITAGE, mec='white', ms=7, label='Heritage site, submerged'),
                Line2D([], [], marker='^', ls='', color=HERITAGE, mec='white', ms=7, label='Heritage site on land'),
                Line2D([], [], marker='s', ls='', color=SCIENCE, mec='white', ms=5.5, label=science_label),
                Patch(facecolor='#c9c8c4', edgecolor='none', label='Land (shaded relief)')]
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=7.2, frameon=False, bbox_to_anchor=(0.5, y),
               title='Water depth, and map symbols', title_fontsize=7.5)


def island_labels(atlas, min_km2=4000.0):
    return [(i['features'][0].replace('Montes ', '') + ' I.', i['summit'][0], i['summit'][1])
            for i in atlas['islands'] if i['area_km2'] >= min_km2 and i['features']]



def png_data_uri(array):
    from PIL import Image
    buffer = io.BytesIO()
    Image.fromarray(array).save(buffer, format='PNG', optimize=True)
    return 'data:image/png;base64,' + base64.b64encode(buffer.getvalue()).decode('ascii')


def globe_textures(height, level, lat, ppd):
    """Equirectangular textures, 180 W at the left edge: colour by depth band with neutral land, the
    sea-clamped surface as 8-bit displacement, and a tangent-space normal map of that surface."""
    h = np.roll(height, height.shape[1] // 2, axis=1)
    depth = level - h
    water = depth > 0
    colour = np.empty(h.shape + (3,), dtype=np.uint8)
    colour[:] = np.round(hexrgb(LAND_PLACEHOLDER) * 255).astype(np.uint8)
    lo = 0.0
    for hi, col in BANDS:
        colour[water & (depth >= lo) & (depth < hi)] = np.round(hexrgb(col) * 255).astype(np.uint8)
        lo = hi
    surf = np.maximum(h, level)
    top = float(surf.max())
    displacement = np.round((surf - level) / (top - level) * 255).astype(np.uint8)
    dy = MOON_RADIUS * np.radians(1 / ppd)
    dx = dy * np.maximum(np.cos(np.radians(lat)), 1e-3)[:, None]
    gy, gx = np.gradient(surf)
    n = np.stack([-(gx / dx) * NORMAL_STRENGTH, (gy / dy) * NORMAL_STRENGTH, np.ones_like(surf)], axis=-1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    normal = np.round((n + 1) / 2 * 255).astype(np.uint8)
    return colour, displacement, normal, (top - level) / 1000.0, water


def write_globe(atlas, height, lat, ppd, sites, targets, hashes):
    level = atlas['sea_level_m']
    colour, displacement, normal, relief_km, water = globe_textures(height, level, lat, ppd)
    looks, params, sources = appearance.build(height, level, lat, ppd)
    area = np.cos(np.radians(lat))[:, None] * np.ones_like(water, dtype=float)
    texture_water = float((area * water).sum() / area.sum())
    appearance_water = float((area * (looks['cloud_noise'][..., 2] > 127)).sum() / area.sum())
    bodies = atlas['bodies']
    comparison = next(r for r in atlas['share_comparison'] if abs(r['share'] - atlas['share']) < 1e-9)
    minus = lambda v: f'{v:,.0f}'.replace('-', '\u2212')
    data = dict(
        share=atlas['share'], levelText=minus(level) + ' m', landColour=LAND_PLACEHOLDER,
        heightRangeKm=relief_km, normalStrength=NORMAL_STRENGTH,
        bands=[[label, col] for (_, col), label in zip(BANDS, BAND_LABELS)],
        figures=[['Water', f"{atlas['share']:.0%} of the surface"], ['Sea level (geoid)', minus(level) + ' m'],
                 ['Near side under water', f"{comparison['near_side_water']:.0%}"],
                 ['Main land mass', f"{atlas['main_land_share']:.1%}"], ['Islands', f"{atlas['island_share']:.1%}"],
                 ['Near-side Sea', f"{bodies[0]['share']:.1%} \u00b7 mean {bodies[0]['mean_depth_m']:,} m"],
                 ['South Pole\u2013Aitken Sea', f"{bodies[1]['share']:.1%} \u00b7 mean {bodies[1]['mean_depth_m']:,} m"],
                 ['Deepest water', f"{max(b['max_depth_m'] for b in bodies):,} m"]],
        textures=dict(colour=png_data_uri(colour), height=png_data_uri(displacement), normal=png_data_uri(normal),
                      albedo=png_data_uri(looks['albedo']), cloudNoise=png_data_uri(looks['cloud_noise']),
                      cloudSun=png_data_uri(looks['cloud_sun']), cloudGeo=png_data_uri(looks['cloud_geo'])),
        appearance=params,
        labels=dict(
            seas=[[name.replace('\n', ' '), la, lo] for name, la, lo in SEA_LABELS],
            islands=[[name, la, lo, i['summit_m']] for (name, la, lo), i in
                     zip(island_labels(atlas), [i for i in atlas['islands'] if i['area_km2'] >= 4000.0 and i['features']])],
            sites=[dict(id=s['id'], name=s['name'], year=s['year'], lat=s['lat'], lon=s['lon'], status=s['status'],
                        short=SITE_LABELS.get(s['id'], ('', None))[0], depth=s.get('depth_m', 0),
                        pressure=s.get('pressure_atm', 0.0), nearestLand=s.get('nearest_land_km', 0),
                        nearestWater=s.get('nearest_water_km', 0), aboveSea=max(s.get('above_sea_m', 0), 0),
                        marks=s['marks'], treatment=s['treatment_type'].capitalize(),
                        evaluation=s['evaluation_status']) for s in sites],
            targets=[dict(id=t['name'], name=next((t['name'].replace(k, v) for k, v in DISPLAY.items() if k in t['name']), t['name']),
                          lat=float(t['lat']), lon=float(t['lon']), aboveSea=int(t['height_relative_to_sea_m']),
                          zone=t['zone'], record=t['record'][0].upper() + t['record'][1:] + '.', sample=t['sample'])
                     for t in targets]),
        provenance=('Geography atlas product <code>terluna.geography.atlas/1</code> (atlas.json '
                    f"{hashes['atlas_json'][:12]}, grid {hashes['atlas_grid'][:12]}); the conservation register "
                    f"({hashes['register'][:12]}); the illumination domain's engine atmosphere "
                    f"({digest(sources['engine'])[:12]}) and Open Moon sky atlas ({digest(sources['sky_atlas'])[:12]}); "
                    f"the climate domain's <code>{params['climatology']['schema']}</code> for run "
                    f"{params['climatology']['run']} ({digest(sources['climatology'])[:12]}). "
                    f"Textures {colour.shape[1]}\u00d7{colour.shape[0]}, 4 pixels per degree."))
    template = (HERE / 'globe.template.html').read_text()
    page = template.replace('__ATLAS_DATA__', json.dumps(data, ensure_ascii=False))
    (OUT / 'globe.html').write_text(page)
    products = {str(Path(v).relative_to(ROOT)): digest(v) for v in sources.values()}
    passed = all(abs(v - atlas['water_share']) <= 0.002 for v in (texture_water, appearance_water))
    return dict(texture_water_share=round(texture_water, 4), appearance_water_share=round(appearance_water, 4),
                product_water_share=atlas['water_share'], tolerance=0.002, passed=passed,
                page_bytes=len(page.encode()), appearance_products=products)


def main():
    atlas_path = ROOT / 'geography' / 'results' / 'atlas.json'
    atlas = json.loads(atlas_path.read_text())
    if atlas.get('schema') != ga.SCHEMA:
        raise SystemExit(f'Expected {ga.SCHEMA} in {atlas_path}')
    pct = round(100 * atlas['share'])
    grid_path = ROOT / 'geography' / 'products' / f'atlas_{pct}pct_4ppd.npz'
    z = np.load(grid_path)
    height, lat, lon = z['height_m'].astype(float), z['lat_deg'], z['lon_deg']
    level, ppd = atlas['sea_level_m'], atlas['grid']['pixels_per_degree']
    reg_path, tgt_path = STUDY / 'heritage_register.json', STUDY / 'science_targets.csv'
    sites = json.loads(reg_path.read_text())['sites']
    with open(tgt_path) as handle:
        targets = list(csv.DictReader(handle))
    rgb = composite(height, level, shaded(height, lat, ppd))
    water = height < level
    note = (f"Approximate atlas. {atlas['share']:.0%} of the surface under water; sea level {level:,.0f} m above the GRAIL geoid, "
            'filling every depression below it. LOLA topography; names from the IAU gazetteer; island and sea names are working labels. '
            'Heritage and targets from the conservation study.')
    OUT.mkdir(parents=True, exist_ok=True)
    islands = island_labels(atlas)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6.9), facecolor=SURF)
    disk_water = None
    for ax, (title, lon0) in zip(axes, [('Near side (facing Earth)', 0.0), ('Far side', 180.0)]):
        img, share = ortho(rgb, water, lon0, ppd)
        disk_water = share if lon0 == 0.0 else disk_water
        ax.imshow(img, extent=(-1, 1, -1, 1), interpolation='bilinear')
        fx = lambda la, lo, lon0=lon0: project(lon0, la, lo)
        for name, la, lo in SEA_LABELS:
            x, y, vis = fx(la, lo)
            if vis:
                ax.text(x, y, name, ha='center', va='center', fontsize=7, color='white', style='italic',
                        path_effects=HALO_DARK, zorder=4, linespacing=0.95)
        for name, la, lo in LAND_LABELS + islands:
            x, y, vis = fx(la, lo)
            if vis:
                ax.text(x, y, name, ha='center', va='center', fontsize=6.3, color=INK2, path_effects=HALO, zorder=4, linespacing=0.95)
        markers(ax, fx, sites, targets)
        ax.set_title(title, fontsize=10, color=INK, loc='left')
        ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect('equal'); ax.axis('off')
    fig.suptitle(f"The Open Moon at {atlas['share']:.0%} water: an approximate atlas", fontsize=13, color=INK, x=0.02, ha='left', y=0.985)
    fig.text(0.02, 0.935, note, fontsize=7, color=INK2, ha='left', va='top', wrap=True)
    legend(fig, 0.035)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.14, wspace=0.03)
    fig.savefig(OUT / f'atlas_{pct}pct_near_far.png', dpi=200, facecolor=SURF)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 7.4), facecolor=SURF)
    ax.imshow(np.roll(rgb, rgb.shape[1] // 2, axis=1), extent=(-180, 180, -90, 90), interpolation='bilinear', aspect='auto')
    for v in range(-150, 181, 30):
        ax.axvline(v, color='white', lw=0.4, alpha=0.6, zorder=2)
    for v in range(-60, 61, 30):
        ax.axhline(v, color='white', lw=0.4, alpha=0.6, zorder=2)
    fg = lambda la, lo: (((lo + 180) % 360) - 180, la, True)
    for name, la, lo in SEA_LABELS:
        ax.text(((lo + 180) % 360) - 180, la, name, ha='center', va='center', fontsize=6.8, color='white', style='italic',
                path_effects=HALO_DARK, zorder=4)
    for name, la, lo in LAND_LABELS + islands:
        ax.text(lo, la, name.replace('\n', ' '), ha='center', va='center', fontsize=6, color=INK2, path_effects=HALO, zorder=4)
    markers(ax, fg, sites, targets, small=True)
    ax.set_xticks(range(-180, 181, 30)); ax.set_yticks(range(-90, 91, 30))
    ax.tick_params(labelsize=7, colors=INK2)
    ax.set_xlabel('Longitude (°E)', fontsize=8, color=INK2); ax.set_ylabel('Latitude (°N)', fontsize=8, color=INK2)
    for s in ax.spines.values():
        s.set_color('#c9c8c4')
    fig.suptitle(f"The Open Moon at {atlas['share']:.0%} water: global sheet", fontsize=13, color=INK, x=0.02, ha='left', y=0.985)
    fig.text(0.02, 0.94, note, fontsize=7, color=INK2, ha='left', va='top', wrap=True)
    legend(fig, 0.01)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.87, bottom=0.17)
    fig.savefig(OUT / f'atlas_{pct}pct_global.png', dpi=200, facecolor=SURF)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.6), facecolor=SURF)
    for ax, south in zip(axes, (True, False)):
        cap, clat = z['south' if south else 'north'].astype(float), z['south_lat' if south else 'north_lat']
        cap_rgb = composite(cap, level, shaded(cap, clat, 16))
        n, rmax = 1000, np.cos(np.radians(80))
        x, y = np.meshgrid(np.linspace(-rmax, rmax, n), np.linspace(rmax, -rmax, n))
        inside = x ** 2 + y ** 2 <= rmax ** 2
        la = (-1 if south else 1) * np.degrees(np.arcsin(np.sqrt(np.clip(1 - x ** 2 - y ** 2, 0, 1))))
        lo = np.degrees(np.arctan2(x, y if south else -y)) % 360
        i = np.clip(np.round(np.abs((clat[0] - la) * 16)).astype(int), 0, cap_rgb.shape[0] - 1)
        j = np.round(lo * 16 - 0.5).astype(int) % cap_rgb.shape[1]
        img = np.ones((n, n, 3)) * hexrgb(SURF)
        img[inside] = cap_rgb[i[inside], j[inside]]
        ax.imshow(img, extent=(-rmax, rmax, -rmax, rmax), interpolation='bilinear')

        def fp(la_, lo_):
            if (south and la_ > -80) or (not south and la_ < 80):
                return 0, 0, False
            c, s_, c_ = np.cos(np.radians(la_)), np.sin(np.radians(lo_)), np.cos(np.radians(lo_))
            return (c * s_, c * c_, True) if south else (c * s_, -c * c_, True)
        for t in targets:
            xx, yy, vis = fp(float(t['lat']), float(t['lon']))
            if vis and t['name'] in POLAR_TARGET_OFFSETS:
                ax.scatter(xx, yy, s=22, marker='s', c=SCIENCE, edgecolors='white', lw=0.8, zorder=5)
                rel = int(t['height_relative_to_sea_m'])
                label = f"{t['name']} ({t['becomes']}, {rel:+,d} m)".replace('-', '−')
                ax.annotate(label, (xx, yy), xytext=POLAR_TARGET_OFFSETS[t['name']], textcoords='offset points', fontsize=6.6,
                            color=INK2, style='italic', path_effects=HALO, zorder=7)
        for st in sites:
            xx, yy, vis = fp(st['lat'], st['lon'])
            if vis:
                ax.scatter(xx, yy, s=40, marker='o' if st['status'] == 'submerged' else '^', c=HERITAGE, edgecolors='white', lw=0.9, zorder=6)
                text, off = SITE_LABELS.get(st['id'], (st['name'], (6, 4)))
                ax.annotate(text, (xx, yy), xytext=off, textcoords='offset points', fontsize=6.8, color=INK, path_effects=HALO, zorder=7)
        for deg in range(0, 360, 90):
            c, s_, c_ = np.cos(np.radians(80)), np.sin(np.radians(deg)), np.cos(np.radians(deg))
            xx, yy = (c * s_ * 1.07, c * c_ * 1.07) if south else (c * s_ * 1.07, -c * c_ * 1.07)
            ax.text(xx, yy, f'{deg}°E', ha='center', va='center', fontsize=6.5, color=INK2)
        ax.set_title('South polar region (poleward of 80°S)' if south else 'North polar region (poleward of 80°N)',
                     fontsize=10, color=INK, loc='left')
        ax.set_xlim(-rmax * 1.13, rmax * 1.13); ax.set_ylim(-rmax * 1.13, rmax * 1.13); ax.set_aspect('equal'); ax.axis('off')
    fig.suptitle(f"Polar cold traps at {atlas['share']:.0%} water", fontsize=13, color=INK, x=0.02, ha='left', y=0.985)
    fig.text(0.02, 0.94, ('Squares mark permanently shadowed craters with floor height relative to sea level. Depressions below sea level '
                          'are drawn filled; lake levels follow groundwater and each catchment\'s rainfall balance. LOLA 16 px/deg cylindrical '
                          'grid over the GRAIL geoid; polar stereographic grids are the next refinement. 0°E faces Earth.'),
             fontsize=7, color=INK2, ha='left', va='top', wrap=True)
    legend(fig, 0.01, 'Permanently shadowed crater (sampling target)')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.14, wspace=0.04)
    fig.savefig(OUT / f'atlas_{pct}pct_poles.png', dpi=200, facecolor=SURF)
    plt.close(fig)

    product_disk = next(r['earth_facing_disk_water'] for r in atlas['share_comparison'] if abs(r['share'] - atlas['share']) < 1e-9)
    products = {str(p.relative_to(ROOT)): digest(p) for p in (atlas_path, grid_path, reg_path, tgt_path)}
    globe = write_globe(atlas, height, lat, ppd, sites, targets,
                        dict(atlas_json=products[str(atlas_path.relative_to(ROOT))],
                             atlas_grid=products[str(grid_path.relative_to(ROOT))],
                             register=products[str(reg_path.relative_to(ROOT))]))
    manifest = dict(products=products, renderer=digest(Path(__file__)), template=digest(HERE / 'globe.template.html'),
                    faithfulness=dict(rendered_near_side_disk_water=round(disk_water, 4), product_disk_water=product_disk,
                                      difference=round(disk_water - product_disk, 4), tolerance=0.01,
                                      passed=abs(disk_water - product_disk) <= 0.01),
                    globe=globe)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps(dict(sheets=manifest['faithfulness'], globe=globe)))
    return 0 if manifest['faithfulness']['passed'] and globe['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
