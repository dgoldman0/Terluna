"""Fact sheets for the atlas's land regions, and the texts written from them.

    ulimit -v 8000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.region_profiles   # about a minute
    python -m research.studies.atlas_round.region_profiles --texts-only      # seconds: after a rename or new texts

Each land region of the districts draft (districts.py) gets a fact sheet from the atlas round's products: where it lies
and how high Earth stands in its sky; its land, lakes, rivers and steep reaches; its coastal waters and the seas they lie
on; the heritage sites on its ground and under its waters; its climate on the moisture map of the design climate; and
its neighbours. Each figure has one display form, the form a text shows it in.

The texts are written by hand from the fact sheets and kept in region_texts.json, keyed by each region's anchor, so that
a rename or a redraw leaves them in place. A text writes a figure as {key}, filled from the current fact sheet, and
another place as [name], tied to that place's anchor when the text was written and shown under the place's current
name. Each text records the facts its wording rests on. When a redraw or a new climate moves one of them past its
tolerance the text is listed as stale, and the tests fail until it is rewritten.
"""

import argparse
import json
import math
import re
import time
from pathlib import Path

import numpy as np

from climate.gcm.turned import to_turned
from research.studies.atlas_round import district_names as DN
from research.studies.atlas_round import districts as D
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import moisture_final as MF

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'region_profiles.json'
TEXTS = HERE / 'region_texts.json'
SCHEMA = 'terluna.research.atlas-round-region-profiles/1'
CODE = ('research/studies/atlas_round/region_profiles.py',)

MOON_KM = D.MOON_KM
SMOOTH_DEG = D.SMOOTH_DEG             # the climate fields are averaged as the moisture map's classes are
NEIGHBOUR_MIN_KM = 20.0               # a shorter shared border makes no neighbour
SNOW_M = 0.01                         # land with a year-mean snow depth of at least this counts as snow-covered
SHARE_TOLERANCE = 0.07                # a share a text rests on may move this much before the text is stale
RELATIVE_TOLERANCE = 0.2              # any other figure, this fraction of itself
SIDE_SHARE = 0.95                     # near side: Earth above the horizon from this share of the land; far side: below

# Each figure a text can show, and the kind that sets its display form.
FIGURES = {
    'area': 'km2', 'waters': 'km2', 'lake_area': 'km2',
    'districts': 'int', 'land_pieces': 'int', 'falls': 'int', 'crater_floors': 'int', 'heritage_sites': 'int',
    'lakes': 'count',
    'lake_share': 'pct', 'basalt_share': 'pct', 'salt_share': 'pct', 'cloud': 'pct', 'snow_share': 'pct',
    'night_rain': 'pct', 'earth_share': 'pct', 'ice_share': 'pct', 'drains_main_share': 'pct',
    'river': 'flow',
    'relief': 'height', 'highest': 'height', 'mean_height': 'height',
    'coast_depth': 'depth',
    'earth': 'deg', 'earth_low': 'deg', 'earth_high': 'deg',
    'air': 'temp', 'air_cool': 'temp', 'air_warm': 'temp', 'air_low': 'temp', 'air_high': 'temp',
    'swing': 'diff',
    'rain': 'mm', 'rain_low': 'mm', 'rain_high': 'mm',
    'falls_gw': 'gw', 'falls_total_gw': 'gw',
    'tide': 'tide',
    'largest_crater': 'km',
    'centre': 'latlon',
}
FIGURE_RE = re.compile(r'\{([a-z_]+)\}')
NB = '\u00a0'                          # a figure's number and unit stay on one line
PLACE_RE = re.compile(r'\[([^\[\]]+)\]')


def sha256(path) -> str:
    return D.sha256(path)


def rel(path) -> str:
    return D.rel(path)


# ------------------------------------------------------------------------------------------- display forms

def _step(v, step):
    return step * math.floor(v / step + 0.5)


def _sig(v, digits):
    if v == 0:
        return 0.0
    e = math.floor(math.log10(abs(v))) - digits + 1
    return _step(v, 10.0 ** e)


def _int(v):
    return f'{int(round(v)):,}'.replace('-', '−')


def display(key, value):
    """A figure in the form a text shows it, by its kind in FIGURES."""
    kind = FIGURES[key]
    if value is None:
        return None
    if kind == 'latlon':
        la, lo = value
        lo = (lo + 180.0) % 360.0 - 180.0
        part = lambda v, pos, neg: f'{abs(v):.0f}°' + ('' if round(abs(v)) == 0 else ' ' + (pos if v >= 0 else neg))
        return part(la, 'N', 'S') + ', ' + part(lo, 'E', 'W')
    v = float(value)
    if kind == 'km2':
        return _int(_step(v, 100) if v < 10_000 else _sig(v, 3)) + NB + 'km²'
    if kind == 'int':
        return _int(v)
    if kind == 'count':
        return _int(v if v < 20 else _step(v, 5) if v < 200 else _step(v, 10) if v < 2000 else _step(v, 100))
    if kind == 'pct':
        p = v * 100.0
        if 0 < p < 1:
            return 'under 1%'
        q = int(_step(p, 1))
        if q == 100 and v < 1:
            q = 99
        return f'{q}%'
    if kind == 'flow':
        return _int(_sig(v, 2) if v >= 1000 else _step(v, 10) if v >= 100 else _step(v, 1)) + NB + 'm³/s'
    if kind == 'height':
        return _int(_step(v, 100) if abs(v) >= 1000 else _step(v, 10)) + NB + 'm'
    if kind == 'depth':
        return _int(_step(v, 10)) + NB + 'm'
    if kind == 'deg':
        return _int(v) + '°'
    if kind == 'temp':
        return _int(v) + NB + '°C'
    if kind == 'diff':
        return f'{v:.1f}' + NB + '°C'
    if kind == 'mm':
        return _int(_step(v, 10) if v >= 100 else _step(v, 5)) + NB + 'mm'
    if kind == 'gw':
        return (_int(v) if v >= 10 else f'{v:.1f}') + NB + 'GW'
    if kind == 'tide':
        return f'{v:.1f}' + NB + 'm'
    if kind == 'km':
        return _int(v) + NB + 'km'
    raise KeyError(kind)


# ------------------------------------------------------------------------------------------------ names

def gloss(region):
    """What a region's name means: the names table's or a calque's meaning when the record holds one, otherwise the
    generator's words, found by composing every head and modifier of its tongue until one gives the name."""
    if region.get('name_meaning'):
        return region['name_meaning']
    t = DN.TONGUES[region['tongue']]
    for head, h in t['heads'].items():
        for mod, m in t['mods'].items():
            if DN.compose(region['tongue'], head, mod) == region['name']:
                return f'{mod} {head}'.title() + ' (' + _parts(t['join'], head, h, mod, m) + ')'
    return ''


def _parts(join, head, h, mod, m):
    h, m = h.lower(), m.lower()
    if join in ('fused', 'rendaku'):
        return f'{m} {mod} + {h} {head}'
    if join == 'olivo':
        return f"h' the + {h} {head} + {m} {mod}"
    if join == 'tagalog':
        return f"{h} {head} + {'-ng' if h[-1] in 'aeiou' else 'na'} + {m} {mod}"
    return f'{h} {head} + {m} {mod}'


def names_by_anchor(record):
    """Every named unit's current name by its anchor: regions, seas, districts and coastal seas."""
    out = {}
    for x in record['regions'] + record['districts'] + record['coastal_seas']:
        if x.get('anchor'):
            out[x['anchor']] = x['name']
    return out


# ------------------------------------------------------------------------------------------------ fields

def _weighted_quantile(v, w, q):
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, q * c[-1])])


def _filled_grid(lat_g, lon_g, rows, cols, values):
    """A field of the moisture map's land cells on the climate run's whole grid, each water cell taking its nearest
    land cell's value, as the map assigns them."""
    from scipy.spatial import cKDTree
    g = np.full((lat_g.size, lon_g.size), np.nan)
    g[rows, cols] = values
    pts = D.unit(*np.meshgrid(lat_g, lon_g, indexing='ij'))
    land = np.isfinite(g)
    _, nearest = cKDTree(pts[land]).query(pts[~land])
    g[~land] = g[land][nearest]
    return g


class Smoother:
    """Gaussian weights over the 12 nearest cells of a climate grid within SMOOTH_DEG, at points in the turned frame:
    the moisture map's classes are averaged so, and its other fields are read the same way here."""

    def __init__(self, lat_g, lon_g, lat_t, lon_t, chunk=400_000):
        from scipy.spatial import cKDTree
        tree = cKDTree(D.unit(*np.meshgrid(lat_g, lon_g, indexing='ij')).reshape(-1, 3))
        sigma = np.radians(SMOOTH_DEG)
        self.index = np.empty((lat_t.size, 12), np.int32)
        self.weight = np.empty((lat_t.size, 12), np.float32)
        self.nearest = np.empty(lat_t.size, np.int32)
        for a in range(0, lat_t.size, chunk):
            d, i = tree.query(D.unit(lat_t[a:a + chunk], lon_t[a:a + chunk]), k=12)
            w = np.exp(-0.5 * (d / sigma) ** 2)
            self.index[a:a + chunk] = i
            self.weight[a:a + chunk] = w / w.sum(axis=1, keepdims=True)
            self.nearest[a:a + chunk] = i[:, 0]

    def mean(self, grid):
        v = np.asarray(grid, float).ravel()
        return (self.weight * v[self.index]).sum(axis=1)

    def nearest_of(self, grid):
        return np.asarray(grid).ravel()[self.nearest]


def earth_elevation(lat, lon, sub):
    """Earth's height above the horizon in degrees, from the frame's point beneath Earth (today's coordinates)."""
    la, lo = np.radians(lat), np.radians(lon)
    sla, slo = np.radians(sub['lat_deg']), np.radians(sub['lon_deg'])
    return np.degrees(np.arcsin(np.clip(np.cos(la) * np.cos(sla) * np.cos(lo - slo) + np.sin(la) * np.sin(sla), -1, 1)))


def borders(region, lat, ppd):
    """Shared border lengths in km between the region grid's values, both ways, from its 4-neighbour cell edges."""
    cell = MOON_KM * np.radians(1.0 / ppd)
    out = {}
    a, b = region[:-1, :], region[1:, :]                       # edges between rows: a cell's width at the boundary
    edge = cell * np.cos(np.radians(lat[:-1] - 0.5 / ppd))[:, None] * np.ones((1, region.shape[1]))
    pairs = [(a, b, edge)]
    a2, b2 = region, np.roll(region, -1, axis=1)               # edges between columns: a cell's height
    pairs.append((a2, b2, np.full(region.shape, cell)))
    for x, y, e in pairs:
        m = x != y
        key = x[m].astype(np.int64) * 1000 + y[m].astype(np.int64)
        u, inv = np.unique(key, return_inverse=True)
        s = np.bincount(inv, e[m])
        for k, km in zip(u, s):
            p, q = divmod(int(k), 1000)
            out[(p, q)] = out.get((p, q), 0.0) + float(km)
            out[(q, p)] = out.get((q, p), 0.0) + float(km)
    return out


# ------------------------------------------------------------------------------------------------ build

def build(log=print):
    from scipy import ndimage
    from scipy.spatial import cKDTree
    started = time.time()
    record = json.loads(D.OUT.read_text())
    if sha256(D.PRODUCT) != record['product']['sha256']:
        raise ValueError('the districts grid does not match the districts record; run districts.py again')
    lat, lon, area = D.grid()
    with np.load(D.PRODUCT) as g:
        region, district = g['region'].astype(np.int32), g['district'].astype(np.int32)
    nland = len(record['districts'])
    land8 = (district >= 1) & (district <= nland)
    waters8 = district > nland
    regions = [r for r in record['regions'] if r['kind'] == 'land']
    seas = {r['code']: r for r in record['regions'] if r['kind'] == 'sea'}
    names = names_by_anchor(record)

    # the design surface: heights, the water and its levels, the frame
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        height16 = z['height_m'].astype(np.float32)
        wcell, wlevel, wbody = z['water_cell'], z['water_level_m'], z['water_body']
    rows16, cols16 = height16.shape
    is_water16 = np.zeros(rows16 * cols16, bool)
    is_water16[wcell] = True
    body16 = np.full(rows16 * cols16, -1, np.int8)
    body16[wcell] = wbody
    level16 = np.full(rows16 * cols16, np.nan, np.float32)
    level16[wcell] = wlevel
    bodies = meta['bodies']
    body_level = {b: float(np.median(wlevel[wbody == b])) for b in range(len(bodies)) if np.any(wbody == b)}
    sea_body = {code: bodies.index(s['anchor'].split(':', 1)[1]) for code, s in seas.items()}
    rotation = np.asarray(meta['frame']['rotation_today_to_turned'], float)
    sub = meta['frame']['sub_earth_point']
    region16 = np.repeat(np.repeat(region, 2, axis=0), 2, axis=1).ravel()
    district16 = np.repeat(np.repeat(district, 2, axis=0), 2, axis=1).ravel()
    lat16 = 90.0 - (np.arange(rows16) + 0.5) / 16
    area16 = np.repeat((MOON_KM * np.radians(1 / 16)) ** 2 * np.cos(np.radians(lat16)), cols16)
    dry16 = (district16 >= 1) & (district16 <= nland) & ~is_water16
    lake16 = (district16 >= 1) & (district16 <= nland) & is_water16 & (body16 == 0)
    log(f'surface read, {time.time() - started:.0f} s')

    # the climate on the moisture map, read at every land and coastal-water cell in the run's turned frame
    m = json.loads(MF.OUT.read_text())
    cells = m['cells']
    lat_g, lon_g = np.asarray(m['grid']['lat'], float), np.asarray(m['grid']['lon'], float)
    rr, cc = np.asarray(cells['row']), np.asarray(cells['col'])
    r8, c8 = np.nonzero(land8 | waters8)
    lat_t, lon_t = to_turned(lat[r8], lon[c8], rotation)
    sm = Smoother(lat_g, lon_g, np.asarray(lat_t, float), np.asarray(lon_t, float))
    field = {}
    for f in ('air_c', 'coldest_3day_air_c', 'warmest_3day_air_c', 'lunar_day_swing_k', 'rain_mm_yr',
              'night_rain_share'):
        field[f] = sm.mean(_filled_grid(lat_g, lon_g, rr, cc, np.asarray(cells[f], float)))
    field['night_regime'] = sm.nearest_of(_filled_grid(lat_g, lon_g, rr, cc, np.asarray(cells['night_regime'], float)))
    with np.load(ROOT / FD.TURNED_CLIMATOLOGY) as z:
        if z['lat'].size != lat_g.size or not np.allclose(z['lat'], lat_g, atol=1e-3) or \
                z['lon'].size != lon_g.size or not np.allclose(z['lon'], lon_g, atol=1e-3):
            raise ValueError('the climatology and the moisture map are on different grids')
        field['cloud'] = sm.mean(z['clt_mean'])
        field['snow_m'] = sm.mean(z['snd_m'])
        field['sic'] = sm.mean(z['sic'])
    # the tide on a coast: each coastal-water cell takes the tide of the nearest of the map's cells beside a sea
    tide_m = np.asarray([np.nan if t is None else t for t in cells['tide_m']], float)
    has = np.isfinite(tide_m)
    tide_tree = cKDTree(D.unit(np.asarray(cells['lat'])[has], np.asarray(cells['lon'])[has]))
    _, ti = tide_tree.query(D.unit(np.asarray(lat_t, float), np.asarray(lon_t, float)))
    field['tide'] = tide_m[has][ti]
    # how far each place lies from the run's nearest land cell, whose climate the map gives it, and whether the run
    # holds the place itself as water: an island the run cannot see takes a shore's climate
    land_tree = cKDTree(D.unit(lat_g[rr], lon_g[cc]))
    chord, _ = land_tree.query(D.unit(np.asarray(lat_t, float), np.asarray(lon_t, float)))
    field['reach_km'] = 2 * np.arcsin(np.clip(chord / 2, 0, 1)) * MOON_KM
    run_land = np.zeros(lat_g.size * lon_g.size, bool)
    run_land[rr * lon_g.size + cc] = True
    field['in_run_water'] = ~run_land[sm.nearest]
    log(f'climate read at {r8.size:,} cells, {time.time() - started:.0f} s')

    w8 = area[r8]
    reg8 = region[r8, c8]
    is_land8 = land8[r8, c8]
    earth8 = earth_elevation(lat[r8], lon[c8], sub)
    falls = json.loads(D.FALLS.read_text())['reaches'] if D.FALLS.exists() else []
    heritage = json.loads(D.HERITAGE.read_text())['sites']
    edges = borders(region, lat, D.PPD)
    regimes = m['night_regimes']
    tide_classes = m['tide_classes']

    profiles = []
    for r in regions:
        n = int(r['code'][1:])
        L = is_land8 & (reg8 == n)
        W = ~is_land8 & (reg8 == n)
        wl = w8[L]
        mean = lambda f, mask=L, ww=None: float(np.average(field[f][mask], weights=w8[mask] if ww is None else ww))
        ds = [d for d in record['districts'] if d['region'] == r['code']]
        cs = [c for c in record['coastal_seas'] if c['region'] == r['code']]

        # where it lies
        e = earth8[L]
        earth_share = float(np.average(e > 0, weights=wl))
        side = 'near' if earth_share >= SIDE_SHARE else 'far' if earth_share <= 1 - SIDE_SHARE else 'limb'
        lat_c = r['at']['turned'][0]
        band = 'equatorial' if abs(lat_c) < 15 else 'middle' if abs(lat_c) < 45 else 'high' if abs(lat_c) < 65 \
            else 'polar'

        # land and its heights above the sea its waters lie on
        sea_area = {}
        for c in cs:
            sea_area[c['sea']] = sea_area.get(c['sea'], 0.0) + c['area_km2']
        main_sea = max(sea_area, key=sea_area.get) if sea_area else None
        base = body_level[sea_body[main_sea]] if main_sea else meta['common_level_m']
        d16 = dry16 & (region16 == n)
        h = height16.ravel()[d16]
        top = int(np.flatnonzero(d16)[np.argmax(h)])
        top_r, top_c = divmod(top, cols16)
        lake_area = float(area16[lake16 & (region16 == n)].sum())
        pieces = int(ndimage.label(land8 & (region == n), structure=np.ones((3, 3)))[1])

        # rivers, lakes, drainage and steep reaches
        drains = {}
        for d in ds:
            drains[d['drains_to']] = drains.get(d['drains_to'], 0.0) + d['area_km2']
        tot = sum(drains.values())
        drains = {k: v / tot for k, v in sorted(drains.items(), key=lambda kv: -kv[1])}
        drains_main = next(iter(drains))
        river_d = max(ds, key=lambda d: d['largest_flow_m3s'])
        here_falls = []
        for x in falls:
            rf, cf = divmod(int(x['start']), cols16)
            if region[rf // 2, cf // 2] == n:
                here_falls.append(x)
        here_falls.sort(key=lambda x: -x['power_gw'])
        craters = [d for d in ds if d.get('centre') == 'crater' and d.get('centre_km')]

        # heritage on its ground and under its waters
        sites = []
        for s in heritage:
            r16 = min(int((90.0 - s['lat']) * 16), rows16 - 1)
            c16 = int((s['lon'] % 360.0) * 16) % cols16
            k = r16 * cols16 + c16
            if region[r16 // 2, c16 // 2] != n:
                continue
            item = dict(name=s['name'], year=s.get('year'), kind=s.get('kind'))
            dk = int(district16[k])
            if is_water16[k]:
                b = int(body16[k])
                item.update(where='water', depth_m=round(float(level16[k] - height16.ravel()[k])),
                            water=('a lake' if b == 0 else next((sv['name'] for code, sv in seas.items()
                                                                  if sea_body[code] == b), bodies[b])))
            else:
                item.update(where='land', height_m=round(float(height16.ravel()[k] - base)))
            if 1 <= dk <= nland + len(record['coastal_seas']):
                unit = (record['districts'] + record['coastal_seas'])[dk - 1]
                item['unit'] = unit['name']
            sites.append(item)

        # neighbours: land regions sharing a border, and the seas whose open waters touch it
        nb = sorted(((q, km) for (p, q), km in edges.items() if p == n and 1 <= q <= len(regions) and
                     km >= NEIGHBOUR_MIN_KM), key=lambda t: -t[1])
        beyond = sorted({q - 100 for (p, q), km in edges.items() if p == n and q >= 100 and km >= NEIGHBOUR_MIN_KM})
        by_number = {int(x['code'][1:]): x for x in regions}
        sea_by_body = {sea_body[code]: code for code in seas}

        regime = np.bincount(field['night_regime'][L].astype(int), wl, len(regimes)) / wl.sum()
        tide = float(np.average(field['tide'][W], weights=w8[W])) if W.any() else None
        tide_class = None
        if tide is not None:
            for tc in tide_classes:
                lo_, hi_ = tc['range_m']
                if tide >= lo_ and (hi_ is None or tide < hi_):
                    tide_class = tc['name']
        salt_area = sum(c['area_km2'] for c in cs if c['salt'])
        cs_area = sum(c['area_km2'] for c in cs)
        rain = field['rain_mm_yr'][L]
        air = field['air_c'][L]

        facts = dict(
            area=r['area_km2'], waters=r['waters_km2'], lake_area=lake_area, districts=r['districts'],
            land_pieces=pieces, falls=len(here_falls), crater_floors=len(craters), heritage_sites=len(sites),
            lakes=sum(d['lakes'] for d in ds), lake_share=r['lake_share'], basalt_share=r['basalt_share'],
            salt_share=salt_area / cs_area if cs_area else 0.0, cloud=mean('cloud'),
            snow_share=float(np.average(field['snow_m'][L] >= SNOW_M, weights=wl)),
            night_rain=mean('night_rain_share'), earth_share=earth_share,
            ice_share=float(np.average(field['sic'][W], weights=w8[W])) if W.any() else 0.0,
            drains_main_share=drains[drains_main], river=river_d['largest_flow_m3s'], relief=r['relief_m'],
            highest=float(h.max() - base), mean_height=float(np.average(h, weights=area16[d16]) - base),
            coast_depth=(sum(c['mean_depth_m'] * c['area_km2'] for c in cs) / cs_area) if cs_area else None,
            earth=float(np.average(e, weights=wl)), earth_low=float(e.min()), earth_high=float(e.max()),
            air=mean('air_c'), air_cool=mean('coldest_3day_air_c'), air_warm=mean('warmest_3day_air_c'),
            air_low=_weighted_quantile(air, wl, 0.1), air_high=_weighted_quantile(air, wl, 0.9),
            swing=mean('lunar_day_swing_k'), rain=mean('rain_mm_yr'), rain_low=_weighted_quantile(rain, wl, 0.1),
            rain_high=_weighted_quantile(rain, wl, 0.9),
            falls_gw=here_falls[0]['power_gw'] if here_falls else 0.0,
            falls_total_gw=sum(x['power_gw'] for x in here_falls), tide=tide,
            largest_crater=max((d['centre_km'] for d in craters), default=0), centre=r['at']['turned'],
            # the categories a text may rest on
            side=side, hemisphere='north' if lat_c >= 0 else 'south', band=band,
            main_zone=max(r['zones'], key=r['zones'].get), main_sea=seas[main_sea]['anchor'] if main_sea else None,
            drains_main=seas[drains_main]['anchor'] if drains_main in seas else drains_main,
            districts_draining_elsewhere=sum(1 for d in ds if d['drains_to'] != drains_main),
            climate_reach_km=float(np.average(field['reach_km'][L], weights=wl)),
            climate_in_water_share=float(np.average(field['in_run_water'][L], weights=wl)),
            climate_cells=int(np.unique(sm.nearest[L]).size),
            night=regimes[int(np.argmax(regime))], island_group=r['island_group'], tongue=r['tongue'],
            tide_class=tide_class, salt=salt_area > 0.5 * cs_area if cs_area else False,
            **{'zone_' + k.replace('-', '_').replace(' ', '_'): v for k, v in r['zones'].items()})
        figures = {k: display(k, facts[k]) for k in FIGURES if facts.get(k) is not None}
        places = {r['anchor']: r['name']}
        detail = dict(
            zones=r['zones'],
            night_regimes={k: round(float(v), 3) for k, v in zip(regimes, regime)},
            drains={(seas[k]['anchor'] if k in seas else k): round(v, 3) for k, v in drains.items()},
            seas={seas[k]['anchor']: round(v / cs_area, 3) for k, v in sea_area.items()} if cs_area else {},
            highest_at=dict(today=[round(float(90.0 - (top_r + 0.5) / 16), 2), round(float((top_c + 0.5) / 16), 2)]),
            heights_above=('the level of ' + seas[main_sea]['name']) if main_sea else 'the common level',
            heights_above_anchor=seas[main_sea]['anchor'] if main_sea else None,
            largest_river=dict(flow_m3s=river_d['largest_flow_m3s'], unit=river_d['name']),
            falls=[dict(power_gw=x['power_gw'], drop_m=x['drop_m'], flow_m3s=x['flow_m3s'], into=x['into'])
                   for x in here_falls[:3]],
            heritage=sites,
            neighbours=[dict(name=by_number[q]['name'], anchor=by_number[q]['anchor'], border_km=round(km))
                        for q, km in nb],
            open_waters=[seas[sea_by_body[b]]['name'] for b in beyond if b in sea_by_body],
            districts=[dict(name=d['name'], anchor=d['anchor'], area_km2=d['area_km2'], zone=d['zone'],
                            lakes=d['lakes'], largest_flow_m3s=d['largest_flow_m3s'],
                            feature=(d['feature'] or {}).get('name'), meaning=d.get('name_meaning') or None,
                            crater_km=d.get('centre_km') if d.get('centre') == 'crater' else None,
                            description=d['description']) for d in sorted(ds, key=lambda d: -d['area_km2'])],
            coastal_seas=[dict(name=c['name'], anchor=c['anchor'], sea=seas[c['sea']]['name'], area_km2=c['area_km2'],
                               mean_depth_m=c['mean_depth_m'], salt=c['salt'], enclosed=c['enclosed'],
                               strait=c['strait']) for c in sorted(cs, key=lambda c: -c['area_km2'])])
        for x in detail['neighbours'] + detail['districts'] + detail['coastal_seas']:
            places[x['anchor']] = x['name']
        for k in list(detail['seas']) + [a for a in detail['drains'] if isinstance(a, str) and a.startswith('sea:')]:
            places[k] = names[k]
        profiles.append(dict(code=r['code'], anchor=r['anchor'], name=r['name'], tongue=r['tongue'],
                             tongue_label=DN.TONGUES[r['tongue']]['label'], tongue_where=DN.WHERE[r['tongue']],
                             gloss=gloss(r), island_group=r['island_group'],
                             facts={k: (round(v, 4) if isinstance(v, float) else v) for k, v in facts.items()},
                             figures=figures, detail=detail, places=places))
    log(f'{len(profiles)} fact sheets, {time.time() - started:.0f} s')
    tracked = {rel(p): sha256(p) for p in (D.OUT, GS.OUT, MF.OUT, D.FALLS, D.HERITAGE)}
    ignored = {str(p): sha256(p) for p in (D.PRODUCT, GS.PRODUCT)}
    ignored[FD.TURNED_CLIMATOLOGY] = sha256(ROOT / FD.TURNED_CLIMATOLOGY)
    return dict(
        schema=SCHEMA, producer=dict(files={f: sha256(ROOT / f) for f in CODE}), inputs=tracked, ignored_inputs=ignored,
        evidence=('Computed facts of the atlas round\'s draft regions on the design surface and the estimated design '
                  'climate: areas, heights, water and drainage from the districts record and the design surface; '
                  'climate from the moisture map and the climatology of the climate run, read on its cells of about '
                  '170 km; Earth\'s height in the sky from the settled frame\'s point beneath Earth. The texts are '
                  'written by hand from these facts and describe the land, water, climate and sky alone.'),
        reading_rule=('One fact sheet per land region. facts holds raw values (areas in km2, heights in metres above '
                      'the level of the sea that holds most of the region\'s coastal waters, temperatures in deg C of '
                      'the air at 2 m, rain in mm a year, shares from 0 to 1, Earth\'s height in degrees above the '
                      'horizon; air_cool and air_warm are the land\'s means of each place\'s coldest and warmest '
                      'three-day means, air_low and air_high the 10th and 90th percentiles of its mean air across the '
                      'land, swing the mean difference between day and night over a lunar day; climate_reach_km is the '
                      'land\'s mean distance to the climate run\'s nearest land cell, whose climate the moisture map '
                      'gives it, and climate_in_water_share the share of its land the run holds as water, where the '
                      'climate is a shore\'s); figures holds the '
                      'same in the form a text shows them. text is the region\'s written text with its figures and '
                      'names filled in, or null where none is written; text_stale lists the facts that have moved '
                      'past their tolerance since it was written.'),
        parameters=dict(smooth_deg=SMOOTH_DEG, neighbour_min_km=NEIGHBOUR_MIN_KM, snow_m=SNOW_M,
                        share_tolerance=SHARE_TOLERANCE, relative_tolerance=RELATIVE_TOLERANCE, side_share=SIDE_SHARE),
        common_level_seas=[s['anchor'] for code, s in seas.items()
                           if abs(body_level[sea_body[code]] - meta['common_level_m']) < 1.0],
        districts_product=record['product']['sha256'], regions=profiles, seconds=round(time.time() - started))


# ------------------------------------------------------------------------------------------------ texts

def load_texts(path=TEXTS):
    return json.loads(Path(path).read_text())['texts'] if Path(path).exists() else {}


def fill(entry, figures, names):
    """A text with its figures and place names filled in; KeyError for a figure or place the entry does not hold."""
    text = entry['text']
    for key in FIGURE_RE.findall(text):
        if key not in FIGURES:
            raise KeyError(f'no figure {key!r}')
        if figures.get(key) is None:
            raise KeyError(f'figure {key!r} has no value here')
    for name in PLACE_RE.findall(text):
        if name not in entry.get('places', {}):
            raise KeyError(f'place [{name}] is not tied to an anchor')
    text = FIGURE_RE.sub(lambda mt: figures[mt.group(1)], text)
    return PLACE_RE.sub(lambda mt: names.get(entry['places'][mt.group(1)], mt.group(1)), text)


def stale(entry, facts, names):
    """The facts a text rests on that have moved past their tolerance, and the places it names that are gone."""
    out = []
    for k, then in entry.get('rests_on', {}).items():
        now = facts.get(k)
        if now is None and then is not None:
            out.append(f'{k}: gone')
        elif isinstance(then, (bool, str)) or then is None:
            if now != then:
                out.append(f'{k}: {then} -> {now}')
        elif isinstance(then, (int, float)):
            is_share = FIGURES.get(k) == 'pct' or k.startswith('zone_') or k.endswith('_share')
            off = abs(now - then) > (SHARE_TOLERANCE if is_share else RELATIVE_TOLERANCE * max(abs(then), 1e-9))
            if off:
                out.append(f'{k}: {then} -> {now}')
    for name, anchor in entry.get('places', {}).items():
        if anchor not in names:
            out.append(f'[{name}]: {anchor} gone')
    return out


def apply_texts(product, record, texts=None, log=print):
    """Fill in each region's text from the table and list the stale ones; the names come from the record."""
    texts = load_texts() if texts is None else texts
    names = names_by_anchor(record)
    written = stale_n = 0
    for p in product['regions']:
        p['name'] = names.get(p['anchor'], p['name'])
        p['places'] = {a: names.get(a, n) for a, n in p['places'].items()}
        entry = texts.get(p['anchor'])
        if not entry:
            p['text'], p['text_stale'] = None, []
            continue
        written += 1
        p['text'] = fill(entry, p['figures'], names)
        p['text_stale'] = stale(entry, p['facts'], names)
        p['text_written'] = entry.get('written')
        if p['text_stale']:
            stale_n += 1
            log(f"stale: {p['name']} ({p['anchor']}): {'; '.join(p['text_stale'])}")
    product['texts'] = dict(written=written, missing=len(product['regions']) - written, stale=stale_n,
                            table=rel(TEXTS), table_sha256=sha256(TEXTS) if TEXTS.exists() else None)
    return product


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--texts-only', action='store_true',
                    help='keep the fact sheets, take the names from the districts record and fill in the texts again')
    args = ap.parse_args(argv)
    record = json.loads(D.OUT.read_text())
    if args.texts_only:
        product = json.loads(OUT.read_text())
        if product['districts_product'] != record['product']['sha256']:
            raise SystemExit('the districts were drawn again since the fact sheets: run without --texts-only')
    else:
        product = build()
    apply_texts(product, record)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, ensure_ascii=False) + '\n')
    t = product['texts']
    print(f"wrote {OUT}: {len(product['regions'])} regions, {t['written']} texts, {t['missing']} missing, "
          f"{t['stale']} stale")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
