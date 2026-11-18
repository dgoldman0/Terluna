#!/usr/bin/env python3
"""Map sheets and a 3D globe of the Open Moon atlas: near and far side, a global sheet, the polar regions,
and an interactive globe page with an appearance mode and a topography mode.

    python visualization/atlas/render.py      # writes visualization/atlas/out/ (kept out of Git)
    python visualization/atlas/render.py --frame <moon_orientation.json> --source <branch@commit:path> \
        --climatology <a turned run's climatology> --drainage geography/turned       # the globe on the turned Moon

With --frame the globe is drawn on the turned Moon (7 November 2026): its heights, rivers, lakes and labels are carried
into the frame of domain/geology's orientation product with climate.gcm.turned's rotation, its clouds and ground
cover come from a run on the turned Moon, and its rivers and lakes from the drainage routed on that run's climate
carried back to today's ground (climate.gcm.turned_climatology). The sheets stay on today's grid.

Reads the geography atlas product (geography/results/atlas.json and its grid in
geography/products) and the conservation study's register and targets. Water is
drawn in depth bands of one blue hue, land as neutral shaded relief. The globe's
appearance mode also reads the illumination domain's engine atmosphere and sky atlas,
the climate domain's GCM climatology and the geography domain's drainage product
(see appearance.py). The topography, as a mode and as an overlay on the appearance, reads
one relief texture, the atlas grid's heights relative to sea level in whole metres. The manifest records
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
from geography import topography as tp
from shared.constants import MOON_RADIUS
from climate.gcm.turned import read_frame, to_turned
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
# Topography: metres below and above sea level, interpolated in sRGB; the two stops at 0 m make the coast an edge.
# Below sea level each of the sheets' depth bands keeps its colour at the band's middle depth (the last band's at
# 5.5 km), between a lighter shore and darker deepest water, so that the globe and the sheets read alike.
BAND_MIDDLES = [125, 375, 750, 1500, 3000, 5500]
TOPO = ([(-7000, '#0a2a5a')] + [(-d, col) for d, (_, col) in reversed(list(zip(BAND_MIDDLES, BANDS)))]
        + [(0, '#d6e8fb'), (0, '#4e8f54'), (300, '#78a85c'), (800, '#a9bd6a'), (1500, '#d5cf86'), (2500, '#d9b47a'),
           (3500, '#c5916a'), (5000, '#a87664'), (7000, '#9a8278'), (9000, '#c4bcb5'), (12000, '#f4f1ed')])
SURFACE_PPD = 16               # the appearance mode's surface and water mask; relief at half this


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


def jpeg_data_uri(array, quality=88):
    from PIL import Image
    buffer = io.BytesIO()
    Image.fromarray(array).save(buffer, format='JPEG', quality=quality, subsampling=0, optimize=True)
    return 'data:image/jpeg;base64,' + base64.b64encode(buffer.getvalue()).decode('ascii')


TEXTURES = 'textures'           # the page's textures as files beside it (the design surface's page), under OUT


def texture_file(name, array, kind='png', quality=88):
    """A texture written as a file beside the page (OUT/textures/<name>.<kind>); returns the path the page reads."""
    from PIL import Image
    (OUT / TEXTURES).mkdir(parents=True, exist_ok=True)
    path = OUT / TEXTURES / f'{name}.{"jpg" if kind == "jpeg" else "png"}'
    image = Image.fromarray(array)
    if kind == 'jpeg':
        image.save(path, format='JPEG', quality=quality, subsampling=0, optimize=True)
    else:
        image.save(path, format='PNG', optimize=True)
    return f'{TEXTURES}/{path.name}'


STOPS_SCHEMA = 'terluna.research.drainage-network-globe-stops/1'
SEA_NAMES = dict(near_side='the Near-side Sea', south_pole_aitken='the South Pole–Aitken Sea',
                 smythii_marginis='Smythii–Marginis', humboldtianum='Humboldtianum', orientale='Orientale',
                 moscoviense='Moscoviense')


def stops_textures(path, rotation, ppd):
    """The drainage network study's stops (globe_stops.json and its maps) on the globe's surface grid, 180 W at the
    left edge, carried into the turned frame when a rotation is given: a bit per stop where one of the network's
    bodies covers the pixel, with the first and second body (one texture), the overflow river's bit per stop and share
    of the pixel with the stops at which the second body is the surface (another), and quarter-resolution copies of
    the first and of those stops for the page's readout."""
    record = json.loads(Path(path).read_text())
    if record.get('schema') != STOPS_SCHEMA:
        raise SystemExit(f'Expected {STOPS_SCHEMA} in {path}')
    z = np.load(ROOT / record['maps'])
    meta = json.loads(str(z['metadata']))
    if meta['grid']['pixels_per_degree'] != ppd or meta['grid']['left_lon'] != 0.0:
        raise SystemExit(f"{record['maps']} is not on the {ppd} px/deg grid with 0 E at the left edge")

    def carried(a, channel=False):
        if rotation is not None:
            a = (appearance.turn_channels(a, rotation, left_lon=0.0) if channel
                 else appearance.turn_texture(a, rotation, left_lon=0.0, nearest=True))
        return np.roll(a, a.shape[1] // 2, axis=1)

    bodies, first, second = (carried(z[k]) for k in ('body_bits', 'first', 'second'))
    second_bits = carried(z['second_bits'])
    channel_q = carried(z['channel_q'].astype(np.float32), channel=True)
    channel_bits = carried(z['channel_bits'], channel=True)
    km = appearance.MOON_KM * np.radians(1.0 / ppd)
    share = np.clip(appearance.RIVER_WIDTH * np.sqrt(channel_q) / (km * 1000.0), 0.0, 1.0)
    stops = np.dstack([bodies, first, second]).astype(np.uint8)
    channels = np.dstack([channel_bits, np.round(share * 255), second_bits]).astype(np.uint8)
    readout = (np.ascontiguousarray(stops[2::4, 2::4]), np.ascontiguousarray(second_bits[2::4, 2::4]).astype(np.uint8))
    return record, stops, channels, readout


def stops_data(record, stops, channels, readout, share_of, colours, atlas_water):
    """The page's years control: each stop's year, levels and seas, the textures, and the faithfulness of the drawn
    water to the product's share at each stop. The water drawn in a year is the bodies' that year with the atlas's own
    water that no body covered at year 0, which stays as drawn."""
    names = [''] + [SEA_NAMES.get(b['name'], b['name']) if b['kind'] == 'sea'
                    else f"a hollow on {SEA_NAMES.get(b['fed_by'], b['fed_by'])}'s overflow" for b in record['bodies']]
    z = np.load(ROOT / record['maps'])
    levels = [[-1.0e5 if not np.isfinite(v) else round(float(v), 1) for v in row] for row in z['levels']]
    held = atlas_water & ((stops[..., 0] & 1) == 0)
    drawn = [share_of((((stops[..., 0] >> k) & 1) == 1) | held) for k in range(len(record['stops']))]
    checks = [dict(year=s['year'], drawn=round(d, 4), product=s['water_share'], difference=round(d - s['water_share'], 4))
              for s, d in zip(record['stops'], drawn)]
    texts = [dict(year=s['year'], share=s['water_share'], full=[SEA_NAMES[k] for k, v in s['states'].items() if v == 'full'],
                  levels={SEA_NAMES[k]: v for k, v in s['levels_m'].items()}) for s in record['stops']]
    data = dict(years=[s['year'] for s in record['stops']], levels=levels, names=names, stops=texts,
                textures=dict(stops=png_data_uri(stops), channels=png_data_uri(channels),
                              readout=png_data_uri(readout[0]), readoutSecond=png_data_uri(readout[1])),
                colours={k: [round(float(c), 4) for c in v] for k, v in colours.items()},
                producer=record['producer'], case=record['case'])
    return data, checks


def relief_texture(height, level):
    """The map and topography modes' texture, 180 W at the left edge: the height relative to sea level in whole
    metres above the lowest cell, high byte in red and low byte in green. The page adds back half a metre, so
    that every cell stays on its side of sea level. Returns the texture, the lowest whole metre and the decoded
    heights."""
    h = np.roll(height, height.shape[1] // 2, axis=1) - level
    low = int(np.floor(h.min()))
    v = np.floor(h - low).astype(np.int64)
    if v.max() > 0xFFFF:
        raise ValueError(f'relief of {v.max()} m does not fit 16 bits')
    rgb = np.zeros(h.shape + (3,), np.uint8)
    rgb[..., 0], rgb[..., 1] = v >> 8, v & 0xFF
    return rgb, low, decode_relief(rgb, low)


def decode_relief(rgb, low):
    """Heights relative to sea level from the relief texture, as the page reads them."""
    return rgb[..., 0].astype(np.int64) * 256 + rgb[..., 1] + low + 0.5


# What the page says of each turned run's climate over the years the globe draws: the ground, the rain and the drylands
# from the run's land readings (climate/gcm/README.md on domain/climate).
TURNED_CLIMATES = {
    'A28R_301k_design_l15_snow_s924_turned_refit':
        'the design run on the turned Moon at the target, {label}: the ground at 301.3 K, rain of 1.98 mm/day over the '
        'globe and 2.19 over land, 59.5% of the land dry (an aridity below 0.65) and the drylands moved south by the turn',
    'A28R_301k_design_l15_snow_s938_turned_refit_precession':
        "the design climate on the turned Moon, {label}, with PlaSim's night-side twilight counted once and the equinox "
        'regressing with the node: the ground at 301.0 K, rain of 1.95 mm/day over the globe and 2.15 over land, and '
        '60.7% of the land dry (an aridity below 0.65), 63.2% in the north and 58.4% in the south',
    'A28S_301k_design_l15_snow_s930_settled':
        'the design climate on the final surface, {label}, at sunlight 0.930: the ground at 300.3 K, rain of 1.93 '
        'mm/day over the globe and 2.13 over land, and 57.1% of the land dry (an aridity below 0.65)',
    # an estimate (climate.gcm.pattern_scaling): its figures are read from the climatology and moisture map drawn
    'A28D_301k_design_l15_snow_s936_settled':
        'the design climate on the design surface, an estimate: {label}, at sunlight 0.936, carried by pattern scaling '
        'from its 301.0 K to the 300.3 K at which the run at 0.930 on the final surface settled, which leaves the ground '
        'at {ground:.1f} K, rain of {rain:.2f} mm/day over the globe and {rain_land:.2f} over land, and {dry:.1%} of the '
        'land dry (an aridity below 0.65)',
}


def climate_figures(climatology, moisture=None):
    """The figures a climate sentence reads when they are not written into it: the ground's mean temperature and the
    rain over the globe and over the land (cos-latitude weighted) from the climatology, and the dry share of the run's
    land cells (the four classes below an aridity of 0.65) from its moisture map."""
    with np.load(climatology) as z:
        w = np.cos(np.radians(z['lat']))[:, None] * np.ones((1, z['lon'].size))
        land = z['lsm'] > 0.5
        out = dict(ground=float((z['ts_k'] * w).sum() / w.sum()), rain=float((z['pr_mm_day'] * w).sum() / w.sum()),
                   rain_land=float((z['pr_mm_day'] * w * land).sum() / (w * land).sum()), dry=float('nan'))
    if moisture:
        m = json.loads(Path(moisture).read_text())
        out['dry'] = sum(m['shares']['by_class'][c['key']]['share'] for c in m['classes'][:4])
    return out


def turned_note(turn, drained, run, years):
    """The page's paragraph on the turned Moon: the turn, the run whose clouds and rain it shows, and the lakes."""
    if run not in TURNED_CLIMATES:
        raise KeyError(f'no climate readings for {run} in TURNED_CLIMATES; add them before drawing the globe from it')
    climate = TURNED_CLIMATES[run].format(label=f'{run}, model years {years[0]}–{years[-1]}')
    lakes, rivers = drained['lakes'], drained['rivers']
    return ('<b>The turned Moon.</b> The weight of the seas and lakes turns the Moon about '
            f"{turn['pole_moved_deg']:.1f}° into a new stable orientation, the orientation study's elastic state, so the "
            "globe's poles, equator, grid and the point beneath Earth are the turned Moon's, and every place keeps its "
            'ground. In today\'s coordinates the turned north pole stands at 81.95° N, 167.14° W and the point beneath '
            'Earth at 8.00° N, 6.09° E: near-side places move about 8° south and far-side places north. The water level, '
            "the lakes' levels and every height stay on the dry Moon's geoid until the atlas folds in the water's own "
            f"gravity and the crust's sag. The clouds and the ground cover follow {climate}. The rivers and rain-fed lakes "
            f"are routed on that run's rain and evaporation carried back to the ground: {lakes['area_share']:.1%} of the "
            f"Moon in {lakes['bodies']:,} lakes holding {lakes['volume_km3'] / 1e6:.1f} million km³, with "
            f"{rivers['discharge_to_sea_m3s']:,} m³/s reaching the seas. The run itself was set with the earlier "
            "estimate's lakes, routed on September's climate.")


def write_globe(atlas, height, ppd, sites, targets, hashes, turn=None, stops=None):
    level = atlas['sea_level_m']
    rotation = turn['rotation'] if turn else None
    if turn:
        height = appearance.turn_texture(height.astype(np.float32), rotation, left_lon=0.0, nearest=True).astype(float)
    place = (lambda la, lo: tuple(round(float(v), 3) for v in to_turned(la, lo, rotation))) if turn \
        else (lambda la, lo: (la, lo))
    relief, low, decoded = relief_texture(height, level)
    water = decoded < 0
    exact = np.roll(height, height.shape[1] // 2, axis=1) - level
    relief_error = float(np.abs(decoded - exact).max())
    relief_off_side = int(((decoded < 0) != (exact < 0)).sum())
    highest, lowest = round(float(exact.max())), round(float(exact.min()))
    del exact
    fine, _, _ = tp.height_above_geoid(SURFACE_PPD, atlas['grid']['geoid_degree'])
    fine = fine.astype(np.float32)
    if turn:
        fine = appearance.turn_texture(fine, rotation, left_lon=0.0, nearest=True)   # keeps the atlas's water share
    looks, params, sources, masks = appearance.build(
        fine, level, SURFACE_PPD, climatology=turn['climatology'] if turn else appearance.CLIMATOLOGY,
        drainage_path=(turn['drainage'] if turn else ROOT / 'geography') / 'products' / appearance.DRAINAGE.name,
        rotation=rotation)
    del fine
    share = lambda mask, ppd_: float((mask * np.cos(np.radians(90 - (np.arange(mask.shape[0]) + 0.5) / ppd_))[:, None]).sum()
                                     / (np.cos(np.radians(90 - (np.arange(mask.shape[0]) + 0.5) / ppd_)).sum() * mask.shape[1]))
    texture_water = share(water, ppd)
    surface_water = share(masks['sea'], SURFACE_PPD)
    lake_water = share(masks['water'], SURFACE_PPD) - surface_water
    mask_water = share(looks['water'] / 255.0, SURFACE_PPD)
    drained = json.loads(((turn['drainage'] if turn else ROOT / 'geography') / 'results' / 'drainage.json').read_text()) \
        if 'drainage' in sources else None
    bodies = atlas['bodies']
    comparison = next(r for r in atlas['share_comparison'] if abs(r['share'] - atlas['share']) < 1e-9)
    minus = lambda v: f'{v:,.0f}'.replace('-', '\u2212')
    data = dict(
        share=atlas['share'], levelText=minus(level) + ' m', normalStrength=appearance.NORMAL_STRENGTH,
        relief=dict(low=low, highest=highest, lowest=lowest),
        topo=[[h, col] for h, col in TOPO],
        figures=[['Water', f"{atlas['share']:.0%} of the surface"], ['Sea level (geoid)', minus(level) + ' m'],
                 ['Near side under water', f"{comparison['near_side_water']:.0%}"],
                 ['Main land mass', f"{atlas['main_land_share']:.1%}"], ['Islands', f"{atlas['island_share']:.1%}"],
                 ['Near-side Sea', f"{bodies[0]['share']:.1%} \u00b7 mean {bodies[0]['mean_depth_m']:,} m"],
                 ['South Pole\u2013Aitken Sea', f"{bodies[1]['share']:.1%} \u00b7 mean {bodies[1]['mean_depth_m']:,} m"],
                 ['Deepest water', f"{max(b['max_depth_m'] for b in bodies):,} m"]]
                + ([['Rain-fed lakes above sea level', f"{drained['lakes']['area_share']:.1%} \u00b7 "
                     f"{drained['lakes']['bodies']:,} lakes"],
                    ['Water in those lakes', f"{drained['lakes']['volume_km3'] / 1e6:.1f} million km\u00b3"],
                    ['Rivers into the seas', f"{drained['rivers']['discharge_to_sea_m3s']:,} m\u00b3/s"]] if drained else [])
                + ([['Turn of the pole', f"{turn['pole_moved_deg']:.2f}\u00b0"]] if turn else []),
        textures=dict(relief=png_data_uri(relief), normal=jpeg_data_uri(looks['normal']), albedo=jpeg_data_uri(looks['albedo']),
                      water=png_data_uri(looks['water']), cloudSun=png_data_uri(looks['cloud_sun']),
                      cloudGeo=png_data_uri(looks['cloud_geo'])),
        appearance=params,
        labels=dict(
            seas=[[name.replace('\n', ' '), *place(la, lo)] for name, la, lo in SEA_LABELS],
            islands=[[name, *place(la, lo), i['summit_m']] for (name, la, lo), i in
                     zip(island_labels(atlas), [i for i in atlas['islands'] if i['area_km2'] >= 4000.0 and i['features']])],
            sites=[dict(id=s['id'], name=s['name'], year=s['year'], **dict(zip(('lat', 'lon'), place(s['lat'], s['lon']))),
                        status=s['status'],
                        short=SITE_LABELS.get(s['id'], ('', None))[0], depth=s.get('depth_m', 0),
                        pressure=s.get('pressure_atm', 0.0), nearestLand=s.get('nearest_land_km', 0),
                        nearestWater=s.get('nearest_water_km', 0), aboveSea=max(s.get('above_sea_m', 0), 0),
                        marks=s['marks'], treatment=s['treatment_type'].capitalize(),
                        evaluation=s['evaluation_status']) for s in sites],
            targets=[dict(id=t['name'], name=next((t['name'].replace(k, v) for k, v in DISPLAY.items() if k in t['name']), t['name']),
                          **dict(zip(('lat', 'lon'), place(float(t['lat']), float(t['lon'])))),
                          aboveSea=int(t['height_relative_to_sea_m']),
                          zone=t['zone'], record=t['record'][0].upper() + t['record'][1:] + '.', sample=t['sample'])
                     for t in targets]),
        provenance=('Geography atlas product <code>terluna.geography.atlas/1</code> (atlas.json '
                    f"{hashes['atlas_json'][:12]}, grid {hashes['atlas_grid'][:12]}); the conservation register "
                    f"({hashes['register'][:12]}); the illumination domain's engine atmosphere "
                    f"({digest(sources['engine'])[:12]}) and Open Moon sky atlas ({digest(sources['sky_atlas'])[:12]}); "
                    f"the climate domain's <code>{params['climatology']['schema']}</code> for run "
                    f"{params['climatology']['run']} ({digest(sources['climatology'])[:12]})"
                    + (f"; the geography domain's drainage product ({digest(sources['drainage'])[:12]})" if drained else '')
                    + (f"; the turned Moon's frame, domain/geology's orientation product ({turn['source']}, "
                       f"{digest(turn['orientation'])[:12]})" if turn else '')
                    + '. '
                    f"Surface {looks['albedo'].shape[1]}\u00d7{looks['albedo'].shape[0]} ({SURFACE_PPD} pixels per degree), "
                    f"relief at {SURFACE_PPD // 2}, the map at {ppd}."))
    if turn:
        years = params['climatology']['years']
        data['frame'] = dict(
            note=turned_note(turn, drained, params['climatology']['run'], years),
            lede=' It shows the Moon as its seas turn it, with the clouds and rain of the design run on the turned Moon.')
    stop_checks = None
    if stops:
        record, stop_tex, channel_tex, readout = stops_textures(stops, rotation, SURFACE_PPD)
        colours = dict(deep=appearance.DEEP_WATER, mid=appearance.MID_WATER, turbid=appearance.TURBID_WATER,
                       river=appearance.RIVER_WATER, bared=appearance.BASALT_SOIL)
        data['stops'], checks = stops_data(record, stop_tex, channel_tex, readout,
                                           lambda m: share(m.astype(np.float32), SURFACE_PPD), colours, masks['sea'])
        unlike = int(((stop_tex[..., 0] & 1).astype(bool) & ~masks['sea']).sum())     # bodies at year 0 off the sea
        stop_checks = dict(stops=checks, year_0_body_cells_off_the_surface_sea=unlike, tolerance=0.002,
                           passed=all(abs(c['difference']) <= 0.002 for c in checks) and unlike <= 0.0005 * masks['sea'].size)
        data['provenance'] += (f" The years: the drainage network study's <code>{STOPS_SCHEMA}</code> "
                               f"({digest(Path(stops))[:12]}; its code {record['producer']['code_sha256'][:12]}).")
        params['legend'].append(['Sea floor laid bare as a sea falls (placeholder)',
                                 '#' + ''.join(f'{int(round(v * 255)):02x}'
                                               for v in appearance.srgb_encode(np.asarray(appearance.BASALT_SOIL)))])
    template = (HERE / 'globe.template.html').read_text()
    page = template.replace('__ATLAS_DATA__', json.dumps(data, ensure_ascii=False))
    (OUT / 'globe.html').write_text(page)
    products = {str(Path(v).relative_to(ROOT)): digest(v) for v in sources.values()}
    lakes = drained['lakes']['area_share'] if drained else 0.0
    passed = (all(abs(v - atlas['water_share']) <= 0.002 for v in (texture_water, surface_water))
              and abs(lake_water - lakes) <= 0.002 and 0.0 <= mask_water - atlas['water_share'] - lakes <= 0.004
              and relief_error <= 0.5 and relief_off_side == 0 and (stop_checks is None or stop_checks['passed']))
    if stops:
        products[str(Path(stops).relative_to(ROOT))] = digest(Path(stops))
    return dict(relief_max_error_m=round(relief_error, 3), relief_cells_off_side=relief_off_side,
                map_texture_water_share=round(texture_water, 4), surface_sea_share=round(surface_water, 4),
                surface_lake_share=round(lake_water, 4), product_water_share=atlas['water_share'],
                product_lake_share=lakes, water_mask_share_with_rivers=round(mask_water, 4), tolerance=0.002,
                years=stop_checks, passed=passed, page_bytes=len(page.encode()), appearance_products=products)


# ---------------------------------------------------------------- the globe on the atlas round's design surface

SURFACE_PRESSURE_ATM = 1.2               # the air at the surface, as the conservation register reckons a site's water
WATER_KG_M3, LUNAR_G = 1000.0, 1.62
NEAREST_WINDOW_DEG = 8.0


def surface_grid(path):
    """The design surface (the atlas round's globe_surface product): heights, each cell's water level (NaN on dry
    ground) and body, on today's grid from 90 N and 0 E."""
    z = np.load(path)
    meta = json.loads(str(z['metadata']))
    rows, cols = meta['grid']['rows'], meta['grid']['cols']
    level = np.full(rows * cols, np.nan, np.float32)
    level[z['water_cell']] = z['water_level_m']
    body = np.zeros(rows * cols, np.int8)
    body[z['water_cell']] = np.where(z['water_body'] > 0, z['water_body'], len(meta['bodies']))
    return dict(meta=meta, height=z['height_m'].astype(np.float32).reshape(rows, cols), level=level.reshape(rows, cols),
                body=body.reshape(rows, cols), rows=rows, cols=cols, ppd=meta['grid']['pixels_per_degree'])


def block(a, k, how):
    r, c = a.shape[0] // k, a.shape[1] // k
    b = a[:r * k, :c * k].reshape(r, k, c, k)
    return getattr(np, how)(b, axis=(1, 3))


def surface_map(grid, k):
    """The map grid's heights and water from the 16 px/deg surface, k cells to a side: each cell's mean ground, and its
    water where most of it lies under water, as the mean depth of that water and the body most of it is."""
    h = block(grid['height'].astype(np.float64), k, 'mean')
    wet16 = np.isfinite(grid['level'])
    share = block(wet16.astype(np.float64), k, 'mean')
    lev = block(np.where(wet16, grid['level'] - grid['height'], 0.0).astype(np.float64), k, 'sum') / np.maximum(block(wet16.astype(np.float64), k, 'sum'), 1)
    bodies = grid['body']
    counts = np.stack([block((bodies == b).astype(np.float64), k, 'sum') for b in range(int(bodies.max()) + 1)])
    counts[0] = -1
    body = np.argmax(counts, axis=0).astype(np.uint8)
    wet = share >= 0.5
    return h, np.where(wet, lev, np.nan), np.where(wet, body, 0).astype(np.uint8)


def surface_texture(depth, body):
    """The water for the page, 180 W at the left edge: its depth in whole metres in two bytes (red and green) and the
    body over the place in blue (0 on dry ground)."""
    lv = np.roll(depth, depth.shape[1] // 2, axis=1)
    bd = np.roll(body, body.shape[1] // 2, axis=1)
    v = np.where(np.isfinite(lv), np.clip(np.round(np.nan_to_num(lv, nan=0.0)), 0, 0xFFFF), 0).astype(np.int64)
    rgb = np.zeros(lv.shape + (3,), np.uint8)
    rgb[..., 0], rgb[..., 1], rgb[..., 2] = v >> 8, v & 0xFF, np.where(np.isfinite(lv), bd, 0)
    return rgb


def nearest_km(mask, i, j, want, ppd, window=NEAREST_WINDOW_DEG):
    """The distance (km) from cell (i, j) to the nearest cell where mask equals `want`, within a window, on the sphere."""
    from scipy import ndimage
    rows, cols = mask.shape
    k = int(window * ppd)
    r_ = np.arange(max(i - k, 0), min(i + k + 1, rows))
    c_ = np.arange(j - k, j + k + 1) % cols
    win = mask[np.ix_(r_, c_)] == want
    if not win.any():
        return round(MOON_RADIUS / 1000.0 * np.radians(window))
    d, (ii, jj) = ndimage.distance_transform_edt(~win, return_indices=True)
    a, b = i - r_[0], k if c_.size > 2 * k else int(np.flatnonzero(c_ == j)[0])
    la1, lo1 = 90.0 - (i + 0.5) / ppd, (j + 0.5) / ppd
    la2, lo2 = 90.0 - (r_[ii[a, b]] + 0.5) / ppd, (c_[jj[a, b]] + 0.5) / ppd
    dlon = (lo2 - lo1 + 180.0) % 360.0 - 180.0
    return round(MOON_RADIUS / 1000.0 * np.radians(np.hypot(la2 - la1, dlon * np.cos(np.radians(0.5 * (la1 + la2))))))


def surface_places(sites, targets, grid, level0, names):
    """The heritage sites and sampling targets on the design surface: each one's ground against the common level, the
    water over it and its depth, the pressure there as the register reckons it, and its nearest land or water."""
    ppd, wet = grid['ppd'], np.isfinite(grid['level'])
    cell = lambda la, lo: (min(int((90.0 - la) * ppd), grid['rows'] - 1), int((lo % 360.0) * ppd) % grid['cols'])
    changed = []
    out_sites = []
    for st in sites:
        i, j = cell(st['lat'], st['lon'])
        ground = float(grid['height'][i, j] - level0)
        under = bool(wet[i, j])
        x = dict(st)
        x['height'] = round(ground)
        x['water_name'] = names[int(grid['body'][i, j])] if under else ''

        if under:
            depth = float(grid['level'][i, j] - grid['height'][i, j])
            x.update(status='submerged', depth_m=round(depth),
                     pressure_atm=round(SURFACE_PRESSURE_ATM + WATER_KG_M3 * LUNAR_G * depth / 101325.0, 1),
                     nearest_land_km=nearest_km(wet, i, j, False, ppd), above_sea_m=0)
        else:
            x.update(status='land', depth_m=0, pressure_atm=0.0, above_sea_m=round(ground),
                     nearest_water_km=nearest_km(wet, i, j, True, ppd))
        if (st['status'] == 'submerged') != under:
            x['status_note'] = ('The register\'s treatment was written for the atlas\'s water, which '
                                + ('left it on land.' if under else 'covered it.'))
            changed.append(st['id'])
        out_sites.append(x)
    out_targets = []
    for t in targets:
        i, j = cell(float(t['lat']), float(t['lon']))
        x = dict(t)
        x['height_relative_to_sea_m'] = str(round(float(grid['height'][i, j] - level0)))
        x['under_m'] = round(float(grid['level'][i, j] - grid['height'][i, j])) if wet[i, j] else -1
        out_targets.append(x)
    return out_sites, out_targets, changed


def surface_note(turn, summary, climate):
    """The page's paragraph on the design surface: the turn, the seas at their levels, the parting and the lakes."""
    seas = summary['seas']
    m = lambda v: f'{abs(v):,.0f} m ' + ('above' if v >= 0 else 'below')
    return ('<b>The design surface.</b> The globe shows the Moon as the atlas round\'s design gives it. The floor is '
            "LOLA's as the water's weight redraws it at the settled state, and the Moon is turned "
            f"{turn['pole_moved_deg']:.2f}° into the orientation that water settles it in, so the globe's poles, "
            'equator, grid and the point beneath Earth are the turned Moon\'s and every place keeps its ground. Heights '
            'are measured from the common level, the level 28% of the surface lies below, '
            f"{abs(summary['common_level_above_the_dry_geoid_m']):,.0f} m below the dry Moon's geoid. Each sea stands "
            f"at its own level: Moscoviense, Orientale and Smythii–Marginis full at their rims, "
            f"{seas['Mare Moscoviense']['level_m']:,.0f}, {seas['Mare Orientale']['level_m']:,.0f} and "
            f"{seas['Smythii–Marginis Sea']['level_m']:,.0f} m above the common level; the South Pole–Aitken Sea and "
            "Mare Humboldtianum at the common level, which fixed weirs on Moscoviense's overflow hold; and the "
            "near-side sea parted at its neck in western Procellarum, its main body "
            f"{m(seas['near-side sea']['level_m'])} the common level and Cognitum, Humorum, Nubium and southern "
            f"Procellarum {m(seas['southern near-side seas']['level_m'])} it, joined by a river of about "
            f"{round(summary['neck_spill_m3s'], -3):,.0f} m³/s over the neck. The clouds and the ground cover follow {climate}. The "
            "rivers and lakes are routed on that run's rain and evaporation carried back to the ground, its open "
            f"water's evaporation raised to close the water budget: {summary['lake_share']:.1%} of the Moon in "
            f"{summary['lakes']:,} lakes holding {summary['lake_volume_km3'] / 1e6:.1f} million km³, with "
            f"{summary['rivers_into_the_seas_m3s']:,} m³/s reaching the seas.")


def routing_product(record):
    """The routing product the design surface was drawn from: named by the routing record among the surface record's
    inputs, and read only while its bytes match that record."""
    name = next(k for k in record['inputs'] if k.endswith('design_routing.json'))
    product = json.loads((ROOT / name).read_text())['product']
    if digest(product['path']) != product['sha256']:
        raise SystemExit(f"{product['path']} no longer matches its record {name}")
    return Path(product['path'])


def cuts_product(record):
    """The works' trenches the design surface was cut with: named by the cuts record among the surface record's inputs,
    and read only while its bytes match that record; None for a surface without cuts."""
    name = next((k for k in record['inputs'] if k.endswith('design_cuts.json')), None)
    if name is None:
        return None
    product = json.loads((ROOT / name).read_text())['product']
    if digest(product['path']) != product['sha256']:
        raise SystemExit(f"{product['path']} no longer matches its record {name}")
    return Path(product['path'])


def work_lines(path, rotation):
    """Each stretch of the works' cuts as the page's canal line: its cells' centres in the order the water takes them,
    as latitude and longitude in the globe's frame. Empty for a cuts product without courses."""
    if path is None:
        return []
    with np.load(path, allow_pickle=False) as z:
        if 'course_cell' not in z.files:
            return []
        cells, stretch = z['course_cell'].astype(np.int64), z['course_stretch']
        cols = json.loads(str(z['metadata']))['grid']['cols']
    ppd = cols / 360.0
    lat = 90.0 - (cells // cols + 0.5) / ppd
    lon = ((cells % cols + 0.5) / ppd + 180.0) % 360.0 - 180.0
    return [[[round(float(v), 4) for v in to_turned(la, lo, rotation)] for la, lo in zip(lat[stretch == s], lon[stretch == s])]
            for s in range(int(stretch.max()) + 1)] if cells.size else []


WATER_REACH_KM = 100.0          # the sea within this of land is drawn as its nearest region's, for borders across it


def water_regions(land_region, sea, reach_km=WATER_REACH_KM):
    """For drawing only: each sea cell within reach_km of a shore takes the region of the nearest shore cell, and 0
    farther out, so that a land border runs on across straits, inlets and island chains and stops at the open sea.
    land_region holds each land cell's page region; the grid runs from 90 N and from 0 E, as the districts' grid."""
    from scipy.spatial import cKDTree
    rows, cols = sea.shape
    ppd = cols / 360.0
    lat = np.radians(90.0 - (np.arange(rows) + 0.5) / ppd)
    lon = np.radians((np.arange(cols) + 0.5) / ppd)
    unit = lambda r, c: np.stack((np.cos(lat[r]) * np.cos(lon[c]), np.cos(lat[r]) * np.sin(lon[c]), np.sin(lat[r])), -1)  # noqa: E731
    beside = np.roll(sea, 1, axis=1) | np.roll(sea, -1, axis=1)
    beside[1:] |= sea[:-1]
    beside[:-1] |= sea[1:]
    sr, sc = np.nonzero(~sea & beside & (land_region > 0))
    out = np.zeros(sea.shape, np.int64)
    if not sr.size:
        return out
    wr, wc = np.nonzero(sea)
    chord = 2.0 * np.sin(reach_km / (MOON_RADIUS / 1000.0) / 2.0)
    dist, k = cKDTree(unit(sr, sc)).query(unit(wr, wc), distance_upper_bound=chord)
    near = np.isfinite(dist)
    out[wr[near], wc[near]] = land_region[sr[k[near]], sc[k[near]]]
    return out


def water_cells(surface_path, shape):
    """Each cell of a coarser grid that is mostly water, a sea's or a lake's, on the design surface's 16 px/deg grid."""
    with np.load(surface_path, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        rows, cols = meta['grid']['rows'], meta['grid']['cols']
        fine = np.zeros(rows * cols, bool)
        fine[z['water_cell']] = True
    k = rows // shape[0]
    return fine.reshape(shape[0], k, shape[1], k).mean(axis=(1, 3)) >= 0.5


def districts_layer(record_path, rotation, surface_path=None):
    """The atlas round's draft regions and districts for the page: a texture in the globe's frame (180 W at the left)
    holding at each 8 px/deg cell the district's place in the page's list (red's low seven bits and green, 0 on open
    water), whether the cell is water (red's top bit, a lake or a sea), and the region
    the page draws there (blue): a land cell's own, and a sea cell within WATER_REACH_KM of a shore its nearest
    region's, 0 farther out; and the lists with each unit's name, code, tongue, label point, area and short
    description; for a coastal sea, the land districts it faces and that feed it, and for a land district its coastal
    seas."""
    record = json.loads(Path(record_path).read_text())
    product = Path(record['product']['path'])
    if digest(product) != record['product']['sha256']:
        raise SystemExit(f'{product} no longer matches its record {record_path}')
    with np.load(product, allow_pickle=False) as z:
        district, region = z['district'].astype(np.int64), z['region'].astype(np.int64)
    regions = record['regions']
    page_region = np.zeros(int(region.max()) + 1, np.int64)
    for k, r in enumerate(regions):
        page_region[int(r['code'][1:]) if r['kind'] == 'land' else 100 + int(r['code'][1:])] = k + 1
    units = record['districts'] + record['coastal_seas']
    if district.max() > 65535 or len(regions) > 255:
        raise SystemExit('too many districts or regions for the page\'s texture')
    if district.max() > 32767:
        raise SystemExit('too many districts for the page\'s texture with its water flag')
    sea = region >= 100
    shown = np.where(sea, 0, page_region[region])
    shown[sea] = water_regions(shown, sea)[sea]
    wet = sea | (water_cells(surface_path, district.shape) if surface_path else np.zeros(sea.shape, bool))
    tex = np.zeros(district.shape + (3,), np.uint8)
    tex[..., 0], tex[..., 1], tex[..., 2] = (district >> 8) | (wet.astype(np.int64) << 7), district & 0xFF, shown
    tex = np.roll(tex, tex.shape[1] // 2, axis=1)
    turned = appearance.turn_texture(tex, rotation, left_lon=-180.0, nearest=True)
    return turned, districts_data(record), Path(record_path), product


def districts_data(record):
    """The page's lists of regions, districts and coastal seas from a districts record: each unit's name, its meaning,
    code, tongue, label point, area, short description and links, indexed as the districts texture holds them."""
    regions = record['regions']
    units = record['districts'] + record['coastal_seas']
    index = {u['code']: k + 1 for k, u in enumerate(units)}
    tongue = {k: t['label'] for k, t in record['tongues'].items()}
    tongue_where = {t['label']: t.get('where', '') for t in record['tongues'].values()}
    region_name = {r['code']: r['name'] for r in regions}
    unit = lambda u, kind: dict(
        index=index.get(u['code'], 0), code=u['code'], name=u['name'], kind=kind, tongue=tongue[u['tongue']],
        lat=u['at']['turned'][0], lon=u['at']['turned'][1], area=u['area_km2'], description=u['description'],
        region=region_name.get(u.get('region'), ''), meaning=u.get('name_meaning', ''),
        links=[index[c] for c in (u.get('faces', []) + u.get('fed_by', []) + u.get('coasts', [])) if c in index][:7])
    data = dict(regions=[unit(r, 'region') | dict(index=0, sea=r['kind'] == 'sea', meaning=r.get('name_meaning', ''))
                         for r in regions],
                districts=[unit(d, 'district') for d in record['districts']]
                + [unit(c, 'coast') for c in record['coastal_seas']],
                summary=record['summary'], tongues={l: f'{l} for {w}' if w else l for l, w in tongue_where.items()},
                coined=record.get('terlunan_coined', {}))
    return data


DATA_OPEN, DATA_CLOSE = '<script id="atlas-data" type="application/json">', '</script>'


def rename_page(page_path, record_path):
    """The built page's regions, districts and coastal seas renamed from a districts record, its textures kept: the
    record must hold the units the page drew, in the same order, so that the districts texture's indices still hold.
    A rename then needs no rebuild, only the page republished."""
    page = Path(page_path).read_text()
    start = page.index(DATA_OPEN) + len(DATA_OPEN)
    end = page.index(DATA_CLOSE, start)
    data = json.loads(page[start:end])
    if 'districts' not in data:
        raise SystemExit(f'{page_path} draws no districts to rename')
    record = json.loads(Path(record_path).read_text())
    new = districts_data(record)
    for key in ('regions', 'districts'):
        if [u['code'] for u in new[key]] != [u['code'] for u in data['districts'][key]]:
            raise SystemExit(f'{record_path} holds other {key} than the page drew; build the page again')
    if record['product']['sha256'] != (data['districts'].get('product_sha256') or record['product']['sha256']):
        raise SystemExit(f'{record_path} names another grid than the page drew; build the page again')
    data['districts'] = dict(new, product_sha256=record['product']['sha256'])
    marker = "Regions and districts: the atlas round's first draft (districts.json "
    at = data['provenance'].find(marker)
    if at >= 0:
        at += len(marker)
        data['provenance'] = data['provenance'][:at] + digest(record_path)[:12] + data['provenance'][at + 12:]
    page = page[:start] + json.dumps(data, ensure_ascii=False) + page[end:]
    Path(page_path).write_text(page)
    return dict(renamed=len(new['regions']) + len(new['districts']), page_bytes=len(page.encode()))


def work_marks(summary, grid, level0, rotation):
    """The page's marks for the works' cuts: each piece of trench at its deepest point, with its card's lines."""
    pieces = summary.get('cuts', {}).get('pieces', [])
    out = []
    for i, x in enumerate(pieces):
        lat, lon = x['at']['turned']
        first = i == 0
        out.append(dict(lat=lat, lon=lon, height=round(x['level_m']), short=f"cut {x['deepest_m']:,.0f} m",
                        title=("Moscoviense's northern cut" if first else 'A cut on the northern branch'),
                        meta=(f"{x['deepest_m']:,.0f} m deep at its deepest · about {x['km']} km of trench · water "
                              f"{x['level_m']:,.0f} m above the common level · {x['flow_m3s']:,} m³/s"),
                        note=('The works of 16 November: ' + ('the cut through Moscoviense\'s northern ridge carries its '
                              'surplus north, Humboldtianum\'s share and Smythii–Marginis\'s. ' if first else
                              'a further cut on the way north. ') + 'Drawn one 1.9 km cell wide, about as wide as a '
                              'channel for its flow needs; the works have not designed its width. The amber line follows '
                              'the course of each cut, from the water it leaves to the water it falls into, a few pixels '
                              'wide at any zoom, as atlases draw canals: far wider than the trench at most zooms.')))
    return out


FALL_MARKS = 12                 # the largest falls the page marks


def fall_marks(record, grid, level0):
    """The page's marks for the largest falls of the falls study: place in the globe's frame, the ground's height there
    above the common level, a short label and the card's lines."""
    heights = grid['height'].ravel()
    great = record['summary'].get('moscoviense_river')
    out = []
    for x in record['reaches'][:FALL_MARKS]:
        lat, lon = x['at']['turned']
        near = ('beside ' if x['beside']['beside'] else 'near ') + x['beside']['name']
        river = (f"On the river that carries Mare Moscoviense's overflow over its natural pass to the "
                 f"{great['into']}: {great['length_km']:,} km, falling {great['falls_m']:,} m in all. "
                 if great and x['on_moscoviense_river'] else '')
        out.append(dict(lat=lat, lon=lon, height=round(float(heights[x['start']] - level0)), short=f"{x['drop_m']:,} m",
                        title=f"{x['drop_m']:,} m drop {near}",
                        meta=(f"{x['flow_m3s']:,} m³/s · {x['drop_m']:,} m in {x['length_km']} km · steepest "
                              f"{x['steepest_deg']}° over a 1.9 km step · {x['power_gw']:.0f} GW of falling water"),
                        note=(river + f"It ends in {x['into']}. LOLA's 59 m heights show the six largest drops as long "
                              'cascades down crater walls, their steepest 60 m 16–28° and no step steeper than 45°: white '
                              'water rather than a free fall, though cliffs of a few tens of metres would not show.')))
    return out


def sea_text(sea, sources):
    """The page's sentences on the water's slopes, from appearance.build's sea parameters."""
    tilt = lambda mss: f"{np.degrees(np.arctan(np.sqrt(mss))):.0f}°"
    return ('The water mirrors the Sun and the sky through the slopes of its waves, which follow the design climate\'s '
            f"daytime wind: the surface wind of the same run and years ({digest(sources['wind'])[:12]}) gives each "
            'climate cell\'s friction velocity over water, on the open seas the median of its stress\'s, '
            f"{sea['openWaterMedianUstar']:.3f} m/s, and over the run's land cells, which hold the lakes and the seas' "
            'edges, the wind 100 m up, where land and water share it, brought down over the water; the illumination '
            "domain's model of a wind-roughened sea at lunar "
            'gravity turns it into slopes: the short waves of Elfouhaily et al. (1997) above 1 m and a fully developed '
            f"sea's longer waves below, a mean square slope of {sea['openWaterMedianMss']:.3f} over the open seas "
            f"({sea['openWaterMss'][0]:.3f}–{sea['openWaterMss'][1]:.3f} across them), an rms tilt of "
            f"{tilt(sea['openWaterMedianMss'])}. Facets with Gaussian slopes, as Cox and Munk measured Earth's, with "
            'Fresnel\'s law and Smith\'s shadowing, give the Sun\'s glitter and the sky the water reflects, which toward '
            'the limb falls to about half what level water would mirror. The longer waves of a fully developed sea are an '
            'upper estimate in winds this light, and the run\'s 3-day means leave out gusts and calm hours, which would '
            'break the glitter into patches. ')


def cover_text(cover, sources, sea=None):
    """The page's paragraph on the ground in true colour, from appearance.build's cover parameters (and its sea
    parameters, when the winds set the water's slopes)."""
    m = cover['moisture']
    at = lambda a: f"{appearance.plant_cover(a):.0%}"
    return ('<b>The ground: informed and guesstimate.</b> The ground\'s colour follows the design climate\'s dryness: the '
            f"aridity of the moisture map (rain over the water the air could draw off; run {m['run']}, model years "
            f"{m['years'][0]}–{m['years'][-1]}; {digest(sources['moisture'])[:12]}) in each of the run's 170 km land cells, "
            'carried between the cells\' centres and smoothed as the clouds\' fields are, the run\'s water cells taking '
            'their nearest land cell\'s. Earth\'s land under the same dryness gives the colours, a guesstimate: plants '
            f"cover {at(0.02)} of the ground at an aridity of 0.02, {at(0.1)} at 0.1, {at(0.2)} at 0.2, {at(0.5)} at 0.5 "
            f"and {at(1.0)} at 1, as Earth's drylands lose their cover fast below 0.2 (Berdugo et al. 2020); trees take a "
            'growing share of it above 0.3, and the grass is straw where it is dry and green where it is wet. Bare '
            'ground is weathered soil. Under the air space weathering ends and oxygen and water rust the rock\'s iron, '
            'so the soils take the dry colours of Earth\'s soils on like rock under like dryness, at the moderate '
            'weathering the terrain\'s upkeep keeps: the rock dust it carries back uphill renews the top metre every '
            '2,400 to 100,000 years. On the feldspathic highland rock, whose 4–6% iron matches Earth\'s continental '
            'crust, they run from light brown in the deserts (Munsell 7.5YR 6/4, with lime) through brown (7.5YR 5/4) to '
            'yellowish red where it is humid (5YR 5/6); on basalt from reddish brown (5YR 4/4) through dark reddish brown '
            '(5YR 3/4) to the dark red of terra roxa (2.5YR 3/6), darker with more titanium. The rock\'s iron and titanium '
            f"come from Lunar Prospector's gamma rays (Prettyman et al. 2006; {digest(sources['composition'])[:12]}), the "
            f"maria's outlines from the LROC team's map (Nelson et al. 2014; {digest(sources['mare'])[:12]}): basalt "
            f"holds {cover['basaltLandShare']:.1%} of the land. Rock shows on steep slopes. The drylands' soils and their "
            'living crusts are chosen with the biome design, so these are the soils the rock and climate give unaided. '
            'Within a cell, ground along larger rivers and in hollows is wetter and steep slopes are drier, and noise '
            'draws patches. The colours read the climate model\'s rain: where it rains at night, about a quarter of the '
            'land, the cloud-resolving boxes\' dry end would leave the ground drier, and lake districts err dry, the '
            'moisture map\'s own limits. The waters take the colours the sea appearance study calculates for the '
            f"designed waters ({digest(sources['water_colours'])[:12]}), the light leaving them under sunlight: deep blue "
            'offshore, the productive coast\'s green over water shallower than about 150 m and in shallow lakes, and the '
            'river mouth\'s ochre in rivers and where they carry sediment into a lake or sea, a plume about 0.3 km wide '
            'per (m³/s)^½ of the runoff gathered since the last lake or sea, since lakes and seas trap their rivers\' '
            'sediment and send clear water on' + ('. ' + sea_text(sea, sources) if sea else '; the water also mirrors the sky. ')
            + ('Rivers turn to white water where they drop steeply, wholly where a 1.9 km step falls 200 m a km, and '
               'spray whitens the ground and water up to about 5 km round the falls of 2 GW or more, display '
               'guesstimates over the falls study\'s steep reaches; the largest twelve falls are marked. '
               if cover.get('falls') else '')
            + 'The relief shades the ground from '
            'its heights at 16 pixels per degree, and the sky\'s diffuse light reaches each place in the share of the '
            'sky its surroundings leave open, found from the heights out to 50 km in 16 directions with the Moon\'s '
            'curvature: under this air the sky gives 43% of the green light on level ground with the Sun 60° up, so '
            'hollows and crater walls darken where they see less of it. River width follows '
            'discharge, by Earth\'s relation widened for the slower flow of lunar gravity; most channels are narrower '
            'than a pixel. Cloud shapes are noise the page draws in 3D, thresholded to the run\'s cover, and the decks\' '
            'opacity is a guess (the low deck holds the run\'s cloud water, the high deck little).')


def write_globe_surface(surface_path, record_path, climatology, sites, targets, hashes, moisture=None, falls=None,
                        wind=None, districts=None):
    """The globe page on the design surface: the relief and the water at each body's own level, the appearance from
    the same water (the ground in true colour from the moisture map when one is given), the sites and targets placed
    on it, and the figures from the surface's record."""
    grid = surface_grid(surface_path)
    record = json.loads(Path(record_path).read_text())
    summary = record['summary']
    meta = grid['meta']
    level0 = float(meta['common_level_m'])
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    pole = float(np.degrees(np.arccos(np.clip(rotation[2, 2], -1, 1))))
    turn = dict(rotation=rotation, pole_moved_deg=pole)
    ppd_map = SURFACE_PPD // 4
    h_map, lev_map, body_map = surface_map(grid, grid['ppd'] // ppd_map)
    turned = lambda a: appearance.turn_texture(a.astype(np.float32), rotation, left_lon=0.0, nearest=True)
    h_map_t = turned(h_map).astype(float)
    lev_map_t = turned(lev_map).astype(float)
    body_map_t = turned(body_map.astype(np.float32)).astype(np.uint8)
    relief, low, decoded = relief_texture(h_map_t, level0)
    exact = np.roll(h_map_t, h_map_t.shape[1] // 2, axis=1) - level0
    relief_error = float(np.abs(decoded - exact).max())
    highest, lowest = round(float(exact.max())), round(float(exact.min()))
    del exact
    surface_tex = surface_texture(lev_map_t, body_map_t)
    fine = turned(grid['height'])
    routing = routing_product(record) if moisture else None
    falls_record = json.loads(Path(falls).read_text()) if falls else None
    if falls_record and digest(falls_record['product']['path']) != falls_record['product']['sha256']:
        raise SystemExit(f"{falls_record['product']['path']} no longer matches its record {falls}")
    looks, params, sources, masks = appearance.build(
        fine, level0, SURFACE_PPD, climatology=climatology, rotation=rotation, surface_path=surface_path,
        moisture_path=moisture, routing_path=routing,
        falls_paths=(Path(falls), Path(falls_record['product']['path'])) if falls_record and moisture else None,
        wind_path=wind)
    del fine
    truly = 'cover' in params
    share = lambda mask, ppd_: float((mask * np.cos(np.radians(90 - (np.arange(mask.shape[0]) + 0.5) / ppd_))[:, None]).sum()
                                     / (np.cos(np.radians(90 - (np.arange(mask.shape[0]) + 0.5) / ppd_)).sum() * mask.shape[1]))
    texture_water = share(surface_tex[..., 2] > 0, ppd_map)
    surface_water = share(masks['sea'], SURFACE_PPD)
    mask_water = share(looks['water'] / 255.0, SURFACE_PPD)
    body_names = ['dry ground'] + ['the ' + b if b.startswith(('near-side', 'southern')) else b for b in meta['bodies'][1:]] + ['a lake']
    places_sites, places_targets, changed = surface_places(sites, targets, grid, level0, body_names)
    place = lambda la, lo: tuple(round(float(v), 3) for v in to_turned(la, lo, rotation))
    seas = summary['seas']
    minus = lambda v: f'{v:,.0f}'.replace('-', '−')
    sea_fig = lambda name, key: [name, f"{seas[key]['share']:.1%} · mean {seas[key]['mean_depth_m']:,} m · "
                                       f"{minus(seas[key]['level_m'])} m"]
    climate = TURNED_CLIMATES[params['climatology']['run']].format(
        label=f"{params['climatology']['run']}, model years {params['climatology']['years'][0]}–{params['climatology']['years'][-1]}",
        **climate_figures(climatology, moisture))
    data = dict(
        share=summary['water_share'], levelText=minus(level0) + ' m', normalStrength=appearance.NORMAL_STRENGTH,
        relief=dict(low=low, highest=highest, lowest=lowest),
        topo=[[h, col] for h, col in TOPO],
        reference='the common level',
        headline='The Open Moon on its design surface',
        lede=('An estimate of how the Open Moon looks from space as the atlas round designs it: its air, clouds and '
              + ('its ground in true colour, plants and soils as the design climate\'s dryness leaves them, '
                 if truly else 'a guessed cover of plants ')
              + 'over the floor its water redraws, with every sea and lake at its own level and the near-side sea parted '
              'at its neck. Drag to turn it, scroll or pinch to zoom, click for the height of a place and the water over '
              'it, and move the Sun to see the halo, the twilight and the view from Earth.'),
        coverText=cover_text(params['cover'], sources, params.get('sea')) if truly else None,
        figures=[['Water', f"{summary['water_share']:.0%} of the surface"],
                 ['Common level (geoid)', minus(level0) + ' m'],
                 sea_fig('Near-side sea', 'near-side sea'), sea_fig('Southern near-side seas', 'southern near-side seas'),
                 sea_fig('South Pole–Aitken Sea', 'South Pole–Aitken Sea'), sea_fig('Mare Humboldtianum', 'Mare Humboldtianum'),
                 ['Upper seas full', f"{sum(seas[k]['share'] for k in ('Mare Moscoviense', 'Mare Orientale', 'Smythii–Marginis Sea')):.1%}"
                                     ' · Moscoviense, Orientale, Smythii–Marginis'],
                 ['Lakes', f"{summary['lake_share']:.1%} · {summary['lakes']:,} lakes"],
                 ['Water in those lakes', f"{summary['lake_volume_km3'] / 1e6:.1f} million km³"],
                 ['Rivers into the seas', f"{summary['rivers_into_the_seas_m3s']:,} m³/s"],
                 ['Over the neck', f"{summary['neck_spill_m3s']:,} m³/s"],
                 ['Deepest water', f"{summary['deepest_water_m']:,} m"],
                 ['Turn of the pole', f"{pole:.2f}°"]],
        lakeFigures=(f"the lakes cover a further {summary['lake_share']:.1%} of the Moon and hold "
                     f"{summary['lake_volume_km3'] / 1e6:.1f} million km³"),
        waterSource=('Land, water and water depth come from LOLA topography at 16 pixels per degree (1.9 km) as the '
                     "water's weight redraws it at the settled state, with every sea and lake at its own level as the atlas "
                     "round's design surface sets it"
                     + (f". The works' cuts from Moscoviense towards Humboldtianum and Smythii–Marginis are cut into the "
                        f"ground down to the water they carry, {summary['cuts']['deepest_m']:,.0f} m at the deepest, and "
                        f"its {summary['cuts']['flow_m3s']:,} m³/s runs through them as a river"
                        if summary.get('cuts', {}).get('trench_cells') else '')
                     + ". Rivers and lakes come from the round's routing of the design climate's runoff over that floor:"),
        topoText=('The ground\'s height above or below the common level, from the design surface\'s heights at 4 pixels '
                  'per degree (7.6 km), with the water taken away. A sea takes the colour of the depth band for its '
                  'depth under its own surface, averaged over each 7.6 km cell most of which it covers, so a sea above '
                  'the common level and the parted near-side sea\'s lower body read as seas like any other. The lakes '
                  'off the seas are left out of these colours, so their beds show; the appearance draws them as water, '
                  'and a click reads them. Ground below the common level that no sea covers, such as the floor the '
                  'parted sea leaves and the hollows that dry on the design rain, takes the lowland\'s colour, and above '
                  'it the colours run from lowland green to the near-white of the highest ground. Shading comes from '
                  f'the full relief at {SURFACE_PPD if truly else SURFACE_PPD // 2} pixels per degree under light that every place has alike, 45° up in its own '
                  'north-west, and contour lines mark the chosen interval, with the common level darker. In the '
                  'topography mode the slider exaggerates the relief down into the sea floors as well as up into the '
                  'highlands. Over the appearance the same colours and contours lie on the rendered globe at the chosen '
                  'opacity, unlit by the Sun. Heritage sites and sampling targets come from the conservation study\'s '
                  'register, each placed on the design surface with the water over it there. In both modes, a click on '
                  'the globe reads the 0.25° cell under it.'),
        surface=dict(names=body_names, levels=[None] + [round(seas[b]['level_m'], 1) for b in meta['bodies'][1:]] + [None]),
        textures=(dict(relief=texture_file('relief', relief), surface=texture_file('surface', surface_tex),
                       normal=texture_file('normal', looks['normal'], 'jpeg'),
                       albedo=texture_file('albedo', looks['albedo'], 'jpeg'), water=texture_file('water', looks['water']),
                       sky=texture_file('sky', looks['sky'], 'jpeg', 92), cloudSun=texture_file('cloud_sun', looks['cloud_sun']),
                       cloudGeo=texture_file('cloud_geo', looks['cloud_geo']),
                       **({'seaSlope': texture_file('sea_slope', looks['sea_slope'])} if 'sea_slope' in looks else {}))
                  if truly else
                  dict(relief=png_data_uri(relief), surface=png_data_uri(surface_tex), normal=jpeg_data_uri(looks['normal']),
                       albedo=jpeg_data_uri(looks['albedo']), water=png_data_uri(looks['water']),
                       cloudSun=png_data_uri(looks['cloud_sun']), cloudGeo=png_data_uri(looks['cloud_geo']))),
        appearance=params,
        labels=dict(
            seas=[[name.replace('\n', ' '), *place(la, lo)] for name, la, lo in SEA_LABELS],
            islands=[],
            falls=fall_marks(falls_record, grid, level0) if falls_record else [],
            works=work_marks(summary, grid, level0, rotation),
            workLines=work_lines(cuts_product(record), rotation),
            sites=[dict(id=s['id'], name=s['name'], year=s['year'], **dict(zip(('lat', 'lon'), place(s['lat'], s['lon']))),
                        status=s['status'], short=SITE_LABELS.get(s['id'], ('', None))[0], depth=s.get('depth_m', 0),
                        height=s['height'], pressure=s.get('pressure_atm', 0.0), nearestLand=s.get('nearest_land_km', 0),
                        nearestWater=s.get('nearest_water_km', 0), aboveSea=max(s.get('above_sea_m', 0), 0),
                        marks=s['marks'], treatment=s['treatment_type'].capitalize(), evaluation=s['evaluation_status'],
                        statusNote=s.get('status_note', ''), waterName=s.get('water_name', '')) for s in places_sites],
            targets=[dict(id=t['name'], name=next((t['name'].replace(k, v) for k, v in DISPLAY.items() if k in t['name']), t['name']),
                          **dict(zip(('lat', 'lon'), place(float(t['lat']), float(t['lon'])))),
                          aboveSea=int(t['height_relative_to_sea_m']), under=t['under_m'],
                          zone=t['zone'], record=t['record'][0].upper() + t['record'][1:] + '.', sample=t['sample'])
                     for t in places_targets]),
        provenance=("The atlas round's design surface <code>" + meta['schema'] + '</code> (globe_surface_16ppd.npz '
                    f"{digest(surface_path)[:12]}; its record {digest(record_path)[:12]}, the settled frame "
                    f"{meta['frame']['sha256'][:12]}); the conservation register ({hashes['register'][:12]}); the "
                    f"illumination domain's engine atmosphere ({digest(sources['engine'])[:12]}) and Open Moon sky atlas "
                    f"({digest(sources['sky_atlas'])[:12]}); the climate domain's <code>{params['climatology']['schema']}"
                    f"</code> for run {params['climatology']['run']} ({digest(sources['climatology'])[:12]})"
                    + (f"; the moisture map <code>{params['cover']['moisture']['schema']}</code> of the same run and years "
                       f"({digest(sources['moisture'])[:12]}), the LROC mare boundaries ({digest(sources['mare'])[:12]}) and "
                       f"the sea appearance study's water colours ({digest(sources['water_colours'])[:12]})" if truly else '')
                    + (f"; the climate domain's <code>{params['sea']['schema']}</code> of the same run and years "
                       f"({digest(sources['wind'])[:12]}), which sets the sea's slopes" if 'wind' in sources else '')
                    + f". Surface {looks['albedo'].shape[1]}×{looks['albedo'].shape[0]} ({SURFACE_PPD} pixels per "
                    f"degree), relief at {SURFACE_PPD if truly else SURFACE_PPD // 2}, the map at {ppd_map}."),
        frame=dict(note=surface_note(turn, summary, climate),
                   lede=''))
    if districts and truly:
        texture, layer, rec_path, product = districts_layer(districts, rotation, surface_path)
        data['textures']['districts'] = texture_file('districts', texture)
        data['districts'] = dict(layer, product_sha256=digest(product))
        data['provenance'] += (f" Regions and districts: the atlas round's first draft (districts.json "
                               f"{digest(rec_path)[:12]}, its grid {digest(product)[:12]}).")
        sources.update(districts=rec_path, districts_grid=product)
    template = (HERE / 'globe.template.html').read_text()
    page = template.replace('__ATLAS_DATA__', json.dumps(data, ensure_ascii=False))
    (OUT / 'globe.html').write_text(page)
    products = {str(Path(v)): digest(v) for v in sources.values()}
    passed = (abs(surface_water - summary['water_share']) <= 0.002 and abs(texture_water - summary['water_share']) <= 0.005
              and mask_water >= surface_water - 1e-4 and relief_error <= 0.5)
    return dict(relief_max_error_m=round(relief_error, 3), map_texture_water_share=round(texture_water, 4),
                surface_water_share=round(surface_water, 4), product_water_share=summary['water_share'],
                water_mask_share_with_rivers=round(mask_water, 4), tolerances=dict(surface=0.002, map=0.005),
                sites_whose_water_changed=changed, passed=passed, page_bytes=len(page.encode()),
                appearance_products=products)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--frame', type=Path, help="domain/geology's orientation product: the globe on the turned Moon")
    parser.add_argument('--source', help='where the orientation product was read, as branch@commit:path')
    parser.add_argument('--climatology', type=Path, default=appearance.CLIMATOLOGY,
                        help="the clouds and ground cover: a climatology in the globe's frame")
    parser.add_argument('--drainage', type=Path, default=ROOT / 'geography',
                        help="the folder holding the drainage product's results/ and products/")
    parser.add_argument('--stops', type=Path, help="the drainage network study's globe_stops.json: the years control")
    parser.add_argument('--surface', type=Path, help="the atlas round's design surface (globe_surface_16ppd.npz): "
                                                     'the globe alone, every water body at its own level')
    parser.add_argument('--surface-record', type=Path, help="the design surface's record (results/globe_surface.json)")
    parser.add_argument('--falls', type=Path, help="with --moisture, the falls study's record (results/falls.json): "
                                                    'the largest falls marked, white water on steep reaches')
    parser.add_argument('--moisture', type=Path, help="with --surface, the moisture map of the climatology's run and years "
                                                      '(results/moisture_map_final.json): the ground in true colour')
    parser.add_argument('--districts', type=Path, help="with --moisture, the atlas round's draft regions and districts "
                                                       '(research/studies/atlas_round/results/districts.json): regions when '
                                                       'zoomed out, districts and their names when zoomed in')
    parser.add_argument('--wind', type=Path, help="with --moisture, the climate domain's surface wind of the climatology's "
                                                  'run and years (products/surface_wind_<run>_y<years>.npz): the sea\'s '
                                                  'slopes, and so its glint and its reflected sky, from the winds')
    parser.add_argument('--names-only', action='store_true',
                        help="with --districts: rename the built page's regions, districts and coastal seas from the "
                             'record, its textures kept (out/globe.html)')
    args = parser.parse_args(argv)
    if args.names_only:
        if not args.districts:
            raise SystemExit('--names-only needs --districts')
        print(json.dumps(rename_page(OUT / 'globe.html', args.districts.absolute())))
        return 0
    if args.surface:
        if not args.surface_record:
            raise SystemExit('--surface needs --surface-record')
        reg_path, tgt_path = STUDY / 'heritage_register.json', STUDY / 'science_targets.csv'
        sites = json.loads(reg_path.read_text())['sites']
        with open(tgt_path) as handle:
            targets = list(csv.DictReader(handle))
        OUT.mkdir(parents=True, exist_ok=True)
        globe = write_globe_surface(args.surface.absolute(), args.surface_record.absolute(), args.climatology.absolute(),
                                    sites, targets, dict(register=digest(reg_path)),
                                    moisture=args.moisture.absolute() if args.moisture else None,
                                    falls=args.falls.absolute() if args.falls else None,
                                    wind=args.wind.absolute() if args.wind else None,
                                    districts=args.districts.absolute() if args.districts else None)
        products = {str(p.relative_to(ROOT)): digest(p) for p in (reg_path, tgt_path)}
        products[str(args.surface)] = digest(args.surface)
        products[str(args.surface_record)] = digest(args.surface_record)
        manifest = dict(products=products, renderer=digest(Path(__file__)), template=digest(HERE / 'globe.template.html'),
                        mode='design surface', globe=globe)
        (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
        print(json.dumps(dict(globe={k: v for k, v in globe.items() if k != 'appearance_products'})))
        return 0 if globe['passed'] else 1
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
    turn = None
    if args.frame:
        if not args.source:
            raise SystemExit('--frame needs --source, where the orientation product was read')
        frame = read_frame(args.frame)
        turn = dict(rotation=frame['rotation'], orientation=args.frame, source=args.source,
                    climatology=args.climatology.absolute(), drainage=args.drainage.absolute(),
                    pole_moved_deg=float(np.degrees(np.arccos(np.clip(frame['rotation'][2, 2], -1, 1)))))
    globe = write_globe(atlas, height, ppd, sites, targets,
                        dict(atlas_json=products[str(atlas_path.relative_to(ROOT))],
                             atlas_grid=products[str(grid_path.relative_to(ROOT))],
                             register=products[str(reg_path.relative_to(ROOT))]), turn,
                        args.stops.absolute() if args.stops else None)
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
