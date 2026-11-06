"""Where a seawater nursery lake can go on the turned Moon: the atlas's closed hollows below the sea level near the
large fresh seas, screened for a native migration (adults down a one-way channel from a fresh sea to spawn in salt
water, their young climbing back).

For every water body of geography's atlas grid (28%, 4 px/deg) other than the large ones the atlas lists, within
MAX_DISTANCE_KM of one of them:
- its lake at a level MARGIN_M below the common sea level: floor, greatest and mean depth, area, volume and salt;
- the distance between the two waters' nearest pixel centres, and the pass, the lowest ground a channel between them
  has to cross (the minimax height along 8-connected paths within a local window), above the sea level;
- its net evaporation on the turned design years with the refit, A28R_301k_design_l15_snow_s924_turned_refit, years
  10-19, read from the biome map's cached sums of that run (session-tools/turned_design/peek/cache) as
  biosphere/ecology/biome_map.py reads them: rain, and Priestley-Taylor's demand on the mean net radiation. The place
  is followed with the turn by the frame record's rotation, as climate/gcm/turned.py does. The lake evaporates as
  the storm-water study's pools did on the design years' semi-arid air, scaled by the nearest land cell's demand over
  that air's (714 mm a year), and its net evaporation is that less the cell's rain.

GCM output and geography's products only; seconds of computation. The ecology and climate modules are not imported:
the reader reproduces their arithmetic, and the tests check it against the turned years' own land reading.
"""
from __future__ import annotations

import hashlib
import heapq
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
ATLAS_GRID = ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz'
FRAME = Path('/home/kir/Documents/Projects/terluna-climate/climate/results/gcm/turned_surface_8d71fd350ca4c6cd.json')
DRIVE = Path('/media/projectspace/terluna-research/session-tools/turned_design')
CACHE = DRIVE / 'peek' / 'cache'
LAND_READING = DRIVE / 'turned_years_land_A28R_301k_design_l15_snow_s924_turned_refit.json'
RUN, YEARS = 'A28R_301k_design_l15_snow_s924_turned_refit', (10, 19)

KEYS = ('pr', 'rss', 'rls', 'tas', 'ps')
LATENT, CP_AIR = 2.5e6, 1004.0          # as biome_map.py and the design run's land readings take them
MM_DAY = 86400e3                        # m/s of water to mm/day
R_KM = MOON_RADIUS / 1e3

MARGIN_M = 1.0                          # the lake's level below the common sea level (the brief: a metre or more)
MAX_DISTANCE_KM = 50.0                  # design choice: hollows screened within this distance of a large water
SEMI_ARID_DEMAND_MM_YR = 714.0          # climate/gcm README: the design years' semi-arid land, where the pools ran
POOL_EVAPORATION_M_YR = (0.804, 0.932)  # this study: the storm-water pools, 3 m deep to 0.3 m with Lake Mead's factor
NAMED = {433: 'Near-side Sea', 1417: 'South Pole-Aitken Sea', 874: 'Smythii-Marginis Sea', 213: 'Humboldtianum Sea',
         1296: 'Orientale', 714: 'Moscoviense'}
# Nursery, stream, route and cost thresholds (design guesses; their basis is in the README)
NURSERY_AREA_KM2 = (100.0, 300.0)       # adequate, good: a larval cohort's spread over months (larval_spread_km)
NURSERY_MEAN_DEPTH_M = 20.0             # open water above the floor for drifting larvae
STREAM_M3_S = 1.0                       # the inflow, the fishway's attraction flow, as a real stream
ROUTE_KM = 20.0                         # to one of the six named seas
SALT_CAP_MT = 10000.0                   # about a decade of the crust top-up's chloride (24 km3 of seawater a year)
LARVAL_DAYS = (120.0, 270.0)            # wild leptocephali's months (recalled) to the eel dossier's 270-day larval cycle
LAKE_DIFFUSIVITY_M2_S = (1.0, 10.0)     # recalled: horizontal eddy diffusivity in lakes at kilometre scales
SALT_MT_PER_KM3 = 35.94                 # this study: seawater at 29.5 C
OUT = HERE / 'results' / 'sites.json'
SCHEMA = 'terluna.research.seawater-lake-sites/1'


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# --- The turned frame and the run ------------------------------------------------------------------------------------
def unit(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) * np.ones_like(lon)], axis=-1)


def to_turned(lat_deg, lon_deg, rotation):
    """Today's latitude and longitude to the turned Moon's (climate/gcm/turned.py's to_turned), longitude in (-180, 180]."""
    p = unit(np.asarray(lat_deg, float), np.asarray(lon_deg, float)) @ np.asarray(rotation, float).T
    return np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0]))


def demand(mean):
    """Priestley-Taylor's demand (mm/day) on the mean net radiation at the ground, as biome_map.demand computes it."""
    t = mean['tas'] - 273.15
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
    s = 4098.0 * es / (t + 237.3) ** 2
    gamma = CP_AIR * mean['ps'] * 100.0 / 1000.0 / (0.622 * LATENT)
    return np.maximum(1.26 * s / (s + gamma) * (mean['rss'] + mean['rls']) / LATENT * 86400.0, 1e-3)


def read_run(cache=CACHE, run=RUN, years=YEARS):
    """The run's means over the ten tenths of the lunar day, from the biome map's cached yearly sums."""
    acc, cnt, inputs = None, None, {}
    for y in range(years[0], years[1] + 1):
        path = Path(cache) / f'{run}.{y:05d}.0.npz'
        z = np.load(path)
        a = {k: np.asarray(z[f'acc_{k}'], float) for k in KEYS}
        if acc is None:
            acc, cnt = a, np.asarray(z['cnt'], float)
            lat, lon, land = np.asarray(z['lat'], float), np.asarray(z['lon'], float), np.asarray(z['land']).astype(bool)
        else:
            acc = {k: acc[k] + a[k] for k in KEYS}
            cnt = cnt + np.asarray(z['cnt'], float)
        inputs[str(path)] = sha256(path)
    mean = {k: (v / np.maximum(cnt, 1)).mean(axis=0) for k, v in acc.items()}
    rain = mean['pr'] * MM_DAY
    pet = demand(mean)
    return dict(lat=lat, lon=lon, land=land, rain=rain, pet=pet, ai=rain / pet, tas=mean['tas']), inputs


def nearest_cell(c, lat_deg, lon_deg, want_land=None):
    """The cell whose centre lies nearest (largest dot product), among land or water cells if asked, and its distance."""
    la, lo = np.meshgrid(c['lat'], c['lon'], indexing='ij')
    dot = (unit(la, lo) * unit(lat_deg, lon_deg)).sum(axis=-1)
    if want_land is not None:
        dot = np.where(c['land'] == want_land, dot, -2.0)
    j, i = np.unravel_index(np.argmax(dot), dot.shape)
    return int(j), int(i), float(np.degrees(np.arccos(min(1.0, dot[j, i]))))


def semi_arid_demand(c):
    """The run's mean demand over its semi-arid land (UNEP aridity 0.2-0.5), area-weighted, mm a year."""
    w = np.cos(np.radians(c['lat']))[:, None] * c['land'] * ((c['ai'] >= 0.2) & (c['ai'] < 0.5))
    return float((c['pet'] * w).sum() / w.sum() * 365.0)


# --- The atlas's hollows ---------------------------------------------------------------------------------------------
def pixel_areas_km2(lat_deg, step_deg):
    lat = np.radians(np.asarray(lat_deg, float))
    h = math.radians(step_deg) / 2.0
    return R_KM ** 2 * math.radians(step_deg) * (np.sin(lat + h) - np.sin(lat - h))


def great_circle_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2.0 * R_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def lake_geometry(heights, areas, level):
    """Floor, depths, area and volume of the water in a hollow's pixels below a level."""
    wet = heights < level
    if not wet.any():
        return None
    depth = (level - heights)[wet]
    area = float(areas[wet].sum())
    volume = float((depth * areas[wet]).sum()) * 1e-3            # km2 m -> km3
    return dict(floor_m=float(heights.min()), level_m=float(level), greatest_depth_m=float(depth.max()),
                area_km2=area, mean_depth_m=volume * 1e3 / area, volume_km3=volume,
                salt_mt=volume * SALT_MT_PER_KM3)


def minimax_pass(h, start, goal):
    """The lowest height a path must reach between two pixel sets of a window: Dijkstra on the largest height met
    along 8-connected paths, the start pixels' own heights included."""
    ny, nx = h.shape
    best = np.full(h.shape, np.inf)
    heap = []
    for j, i in zip(*np.nonzero(start)):
        best[j, i] = h[j, i]
        heap.append((float(h[j, i]), int(j), int(i)))
    heapq.heapify(heap)
    while heap:
        cost, j, i = heapq.heappop(heap)
        if cost > best[j, i]:
            continue
        if goal[j, i]:
            return cost
        for dj in (-1, 0, 1):
            for di in (-1, 0, 1):
                jj, ii = j + dj, i + di
                if (dj or di) and 0 <= jj < ny and 0 <= ii < nx:
                    c = max(cost, float(h[jj, ii]))
                    if c < best[jj, ii]:
                        best[jj, ii] = c
                        heapq.heappush(heap, (c, jj, ii))
    return math.inf


def rank_key(c):
    """Most criteria met first (nursery, stream, route, affordable salt), then a good nursery, the shorter route, the
    lower pass and the stronger stream."""
    m = c['meets']
    return (-sum(m[k] for k in ('nursery', 'stream', 'route', 'affordable')), not m['good_nursery'], c['distance_km'],
            c['pass_above_sea_m'], -min(c['inflow_m3_s']))


def larval_spread_km(diffusivity_m2_s, days):
    """One standard deviation of a larval cohort's spread by horizontal eddy diffusion, sqrt(2 K t)."""
    return math.sqrt(2.0 * diffusivity_m2_s * days * 86400.0) * 1e-3


def screen(margin_m=MARGIN_M, max_km=MAX_DISTANCE_KM):
    """The candidate hollows near the large waters, with their lakes, routes and turned climate."""
    from scipy import ndimage

    inputs = {str(p.relative_to(ROOT)): sha256(p) for p in (ATLAS, ATLAS_GRID)}
    atlas = json.loads(ATLAS.read_text())
    grid = np.load(ATLAS_GRID)
    h, lab = np.asarray(grid['height_m'], float), np.asarray(grid['water_label'])
    lat, lon = np.asarray(grid['lat_deg'], float), np.asarray(grid['lon_deg'], float)
    step = abs(float(lat[1] - lat[0]))      # the grid runs north to south
    sea_level = float(atlas['sea_level_m'])
    large = {b['label']: NAMED.get(b['label'], f"unnamed common-level lake (atlas body {b['label']})")
             for b in atlas['bodies']}
    frame = json.loads(FRAME.read_text())
    inputs[str(FRAME)] = sha256(FRAME)
    rot = np.asarray(frame['frame']['rotation_today_to_turned'], float)
    clim, run_inputs = read_run()
    inputs.update(run_inputs)
    ref_turned = semi_arid_demand(clim)

    areas = pixel_areas_km2(lat, step)[:, None] * np.ones((1, lon.size))
    ny, nx = h.shape
    w_lat = int(math.ceil(max_km / (R_KM * math.radians(step)))) + 1
    level = sea_level - margin_m
    out = []
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        if sl is None or k in large:
            continue
        r0, r1 = sl[0].start, sl[0].stop
        c0, c1 = sl[1].start, sl[1].stop
        lat_max = max(abs(lat[r0]), abs(lat[r1 - 1]))
        if lat_max > 79.5:                       # the cylindrical grid misstates small polar craters (atlas README)
            continue
        w_lon = int(math.ceil(w_lat / max(math.cos(math.radians(min(lat_max + w_lat * step, 89.0))), 0.05)))
        rows = np.arange(max(r0 - w_lat, 0), min(r1 + w_lat, ny))
        if c1 - c0 > nx // 2:
            cols = np.arange(nx)
        else:
            cols = np.arange(c0 - w_lon, c1 + w_lon) % nx
        wl, wh = lab[np.ix_(rows, cols)], h[np.ix_(rows, cols)]
        mine = wl == k
        seas_here = [s for s in large if (wl == s).any()]
        if not seas_here:
            continue
        geo = lake_geometry(wh[mine], areas[np.ix_(rows, cols)][mine], level)
        if geo is None:
            continue
        jy, ix = np.nonzero(mine)
        plat, plon = lat[rows][jy], lon[cols][ix]
        best = None
        for s in seas_here:
            sy, sx = np.nonzero(wl == s)
            d = great_circle_km(plat[:, None], plon[:, None], lat[rows][sy][None, :], lon[cols][sx][None, :]).min()
            if best is None or d < best[1]:
                best = (s, float(d))
        if best[1] > max_km:
            continue
        pass_m = minimax_pass(wh, mine, wl == best[0]) - sea_level
        v = (unit(plat, plon) * areas[np.ix_(rows, cols)][mine][:, None]).sum(axis=0)
        clat = float(np.degrees(np.arcsin(v[2] / np.linalg.norm(v))))
        clon = float(np.degrees(np.arctan2(v[1], v[0])))
        tlat, tlon = (float(x) for x in to_turned(clat, clon, rot))
        jl, il, dl = nearest_cell(clim, tlat, tlon, True)
        ja, ia, da = nearest_cell(clim, tlat, tlon, None)
        rain = float(clim['rain'][jl, il]) * 365.0 * 1e-3
        pet = float(clim['pet'][jl, il]) * 365.0
        e_lake = tuple(x * pet / SEMI_ARID_DEMAND_MM_YR for x in POOL_EVAPORATION_M_YR)
        net = tuple(x - rain for x in e_lake)
        out.append(dict(
            label=k, sea=large[best[0]], sea_label=best[0], distance_km=best[1], pass_above_sea_m=pass_m,
            today=[round(clat, 3), round(((clon + 180.0) % 360.0) - 180.0, 3)], turned=[round(tlat, 3), round(tlon, 3)],
            **geo,
            cell=dict(lat=float(clim['lat'][jl]), lon_e=float(clim['lon'][il]), distance_deg=dl,
                      nearest_cell_is_land=bool(clim['land'][ja, ia]), rain_m_yr=rain, demand_mm_yr=pet,
                      aridity=float(clim['ai'][jl, il]), air_c=float(clim['tas'][jl, il] - 273.15)),
            lake_evaporation_m_yr=e_lake, net_evaporation_m_yr=net,
            inflow_m3_s=tuple(x * geo['area_km2'] * 1e6 / (365.25 * 86400.0) for x in net)))
    for c in out:
        c['meets'] = dict(nursery=c['area_km2'] >= NURSERY_AREA_KM2[0] and c['mean_depth_m'] >= NURSERY_MEAN_DEPTH_M,
                          good_nursery=c['area_km2'] >= NURSERY_AREA_KM2[1] and c['mean_depth_m'] >= NURSERY_MEAN_DEPTH_M,
                          stream=min(c['inflow_m3_s']) >= STREAM_M3_S,
                          route=c['distance_km'] <= ROUTE_KM and c['sea_label'] in NAMED,
                          affordable=c['salt_mt'] <= SALT_CAP_MT)
    out.sort(key=rank_key)
    for n, c in enumerate(out, start=1):
        c['rank'] = n
    t21 = dict(cell_deg=float(abs(clim['lat'][0] - clim['lat'][1])), cell_km_north_south=R_KM * math.radians(360.0 / clim['lon'].size),
               equatorial_cell_km2=R_KM ** 2 * math.radians(360.0 / clim['lon'].size) * math.radians(360.0 / clim['lon'].size))
    spread = {f'{K:g}_m2_s_{t:g}_days': larval_spread_km(K, t) for K in LAKE_DIFFUSIVITY_M2_S for t in LARVAL_DAYS}
    return dict(run=f'{RUN}:{YEARS[0]}-{YEARS[1]}', sea_level_m=sea_level, margin_m=margin_m, max_distance_km=max_km,
                semi_arid_demand_mm_yr=dict(design_years_pools=SEMI_ARID_DEMAND_MM_YR, turned_refit=ref_turned),
                thresholds=dict(nursery_area_km2=NURSERY_AREA_KM2, nursery_mean_depth_m=NURSERY_MEAN_DEPTH_M,
                                stream_m3_s=STREAM_M3_S, route_km=ROUTE_KM, salt_cap_mt=SALT_CAP_MT),
                larval_spread_km=spread, t21=t21, large_waters={str(k): v for k, v in large.items()},
                counts=dict(candidates=len(out), meeting_all=sum(all(c['meets'][k] for k in
                                                                     ('nursery', 'stream', 'route', 'affordable'))
                                                                 for c in out)),
                candidates=out), inputs


def rounded(x, digits=4):
    if isinstance(x, float):
        return float(f'{x:.{digits}g}') if math.isfinite(x) else None
    if isinstance(x, (list, tuple)):
        return [rounded(v, digits) for v in x]
    if isinstance(x, dict):
        return {k: rounded(v, digits) for k, v in x.items()}
    return x


def build() -> dict:
    result, inputs = screen()
    return dict(
        schema=SCHEMA,
        producer=dict(file='research/studies/seawater_lakes/sites.py', sha256=sha256(Path(__file__))),
        inputs=inputs,
        evidence=('A screen on geography\'s atlas (4 px/deg hydrostatic fill at 28%) and the turned design years with '
                  'the refit (GCM output at T21, read from the biome map\'s cached sums). Lake evaporation carries the '
                  'storm-water pools\' energy-balance result to each cell by its Priestley-Taylor demand; thresholds '
                  'for nursery, stream, route and salt are design guesses. Distances and passes are good to about a '
                  'pixel, 7.6 km and the pixel\'s relief; each hollow takes the climate of a cell of about 29,000 km2.'),
        reading_rule=('Heights in metres above the geoid; the lake stands margin_m below the common sea level. today '
                      'and turned are [lat, lon east in (-180, 180]]; turned follows the frame record\'s rotation. '
                      'Ranges are (low, high) across the pools\' evaporation. pass_above_sea_m is the lowest ground a '
                      'channel crosses to the named sea, from the 4 px/deg grid. Candidates are ranked by rank_key.'),
        **rounded(result))


def main(argv=None) -> int:
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}: {record["counts"]}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
