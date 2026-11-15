"""Where the rivers end: the large seas' river inflow on the design years, under three rules for the waters below
the common sea level.

    ulimit -v 9000000; python -m research.studies.south_sea_ocean.inflows   # results/inflows.json; 3.5 minutes, 5 GB

geography's drainage counts every cell below the common sea level as sea, so a river that ends in a crater hollow
below that level counts as reaching the seas, and the hollow stands at the common level whatever its own balance. A
sea's own balance needs to know which rivers reach it. This runner repeats geography's routing of the design years'
runoff (geography.drainage.route, on the climatology geography/turned_corrected names) three times, changing only
which waters take rivers in without spilling:

- every_water: every water below the common level, as geography's drainage has it. The run must reproduce that
  product's discharge to the seas.
- six_seas: the six large seas only (Near-side, South Pole-Aitken, Smythii-Marginis, Humboldtianum, Orientale,
  Moscoviense). Every other hollow below the common level is a depression like those above it: it keeps the lake its
  own balance of inflow, rain and evaporation allows, and spills over its lowest pass once full.
- two_seas: the Near-side and South Pole-Aitken seas only, so the four smaller seas spill too once full.

The rules bracket the answer; none is a prediction of the levels the waters would find. The lakes' open-water balance
and the runoff stay the design years' in every run, although a hollow that shrinks or a sea that spills would change
the climate that made them. The mouths of each run go to research/runs/south_sea_ocean (ignored by Git).
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'inflows.json'
RUNS = ROOT / 'research' / 'runs' / 'south_sea_ocean'
SCHEMA = 'terluna.research.south-sea-inflows/1'
ATLAS = 'geography/results/atlas.json'
DRAINAGE = 'geography/turned_corrected/results/drainage.json'
FRAME = 'climate/results/gcm/turned_surface_8d71fd350ca4c6cd.json'
RUN = 'A28R_301k_design_l15_snow_s938_turned_refit_precession'
CLIMATOLOGY = f'climate/gcm/products/climatology_{RUN}_today.npz'
YEARS = [10, 19]
CODE = ('research/studies/south_sea_ocean/inflows.py', 'geography/drainage.py', 'geography/topography.py',
        'geography/atlas.py', 'climate/gcm/turned.py')
# The six large seas by the atlas's labels; each is matched to the 16 pixel/degree water body nearest its centroid.
SEAS = {433: 'near_side', 1417: 'south_pole_aitken', 874: 'smythii_marginis', 213: 'humboldtianum',
        1296: 'orientale', 714: 'moscoviense'}
RULES = dict(
    every_water='every water below the common level takes its rivers in (geography\'s drainage)',
    six_seas='the six large seas take rivers in; every other hollow below the common level fills by its own balance '
             'and spills when full',
    two_seas='the Near-side and South Pole-Aitken seas take rivers in; every other water fills by its own balance '
             'and spills when full')
TWO = ('near_side', 'south_pole_aitken')
TOP = 8                                 # mouths listed for each sea
HOLLOW_M3S = 5000.0                     # hollows listed: those taking this much river water at the common level
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY
KM3_YR = YEAR_S / 1e9                   # m3/s to km3 a year


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def label_waters(below):
    """Connected waters (8-neighbour) on a latitude-longitude grid whose longitude is periodic. Returns the label of
    each cell (0 on land) with the labels across the seam merged."""
    from scipy import ndimage
    body, n = ndimage.label(below, structure=np.ones((3, 3), bool))
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    rows = below.shape[0]
    for dr in (-1, 0, 1):               # a cell on the last column touches three cells of the first
        a = body[max(0, -dr):rows - max(0, dr), -1]
        b = body[max(0, dr):rows - max(0, -dr), 0]
        for u, v in zip(a, b):
            if u and v:
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[max(ru, rv)] = min(ru, rv)
    root = np.array([find(i) for i in range(n + 1)])
    return root[body]


def centroid(lat_deg, lon_deg, weight):
    """The area-weighted mean direction of a set of cells, as latitude and longitude in degrees."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    x = float((weight * np.cos(lat) * np.cos(lon)).sum())
    y = float((weight * np.cos(lat) * np.sin(lon)).sum())
    z = float((weight * np.sin(lat)).sum())
    return float(np.degrees(np.arctan2(z, np.hypot(x, y)))), float(np.degrees(np.arctan2(y, x)))


def separation_deg(a, b):
    """Great-circle angle between two (latitude, longitude) places in degrees."""
    la, lo, lb, lp = map(np.radians, (a[0], a[1], b[0], b[1]))
    c = np.sin(la) * np.sin(lb) + np.cos(la) * np.cos(lb) * np.cos(lo - lp)
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def name_seas(body, cell_area, lat, lon, atlas_bodies):
    """Match the atlas's six large seas to the water bodies of the routing grid: the largest bodies, each named by
    the atlas sea whose centroid lies nearest its own."""
    rows, cols = body.shape
    flat = body.ravel()
    area = np.bincount(flat, weights=np.repeat(cell_area, cols), minlength=int(flat.max()) + 1)
    area[0] = 0.0
    largest = np.argsort(area)[::-1][:len(SEAS)]
    wanted = {SEAS[b['label']]: b for b in atlas_bodies if b['label'] in SEAS}
    named = {}
    for k in largest:
        r, c = np.nonzero(body == k)
        here = centroid(lat[r], lon[c], cell_area[r])
        name = min(wanted, key=lambda s: separation_deg(here, wanted[s]['centroid']))
        if name in named:
            raise SystemExit(f'Two water bodies match the atlas sea {name}')
        named[name] = dict(label=int(k), area_km2=float(area[k]) / 1e6, centroid=here,
                           atlas_area_km2=float(wanted[name]['area_km2']))
    return named, area


def filling(r, cells, heights, sea_level):
    """How a water that is not an outlet fills in a routing. Every pit on one basin's floor shares the basin's pass,
    so the level of the depressions under the water's cells is the level at which the basin spills: it is given as
    the area-weighted median with its lowest and highest values, with the area that level floods in those
    depressions' catchments and the water that stands between the common level and the spill level there."""
    node, area = r['label'], r['area']
    nodes = np.unique(node[cells])
    above = r['level'][node[cells]] - sea_level
    order = np.argsort(above)
    cumulative = np.cumsum(area[cells][order])
    median = float(above[order][np.searchsorted(cumulative, 0.5 * cumulative[-1])])
    inside = np.flatnonzero(np.isin(node, nodes))
    flooded = inside[heights[inside] < r['level'][node[inside]]]
    rise = r['level'][node[flooded]] - np.maximum(heights[flooded], sea_level)
    return dict(spill_level_above_common_m=round(median, 1), lowest_m=round(float(above.min()), 1),
                highest_m=round(float(above.max()), 1), area_at_the_common_level_km2=round(float(area[cells].sum()) / 1e6),
                area_flooded_at_the_spill_level_km2=round(float(area[flooded].sum()) / 1e6),
                volume_above_the_common_level_km3=round(float((area[flooded] * rise).sum()) / 1e9))


def build():
    from geography import drainage as dr
    from geography import topography as tp
    from climate.gcm.turned import to_turned

    started = time.time()
    atlas = json.loads((ROOT / ATLAS).read_text())
    committed = json.loads((ROOT / DRAINAGE).read_text())
    frame = json.loads((ROOT / FRAME).read_text())
    rotation = np.array(frame['frame']['rotation_today_to_turned'])
    clim = dict(np.load(ROOT / CLIMATOLOGY))
    meta = json.loads(str(clim.pop('metadata')))
    if meta.get('run') != RUN or list(meta.get('years')) != YEARS:
        raise SystemExit(f'{CLIMATOLOGY} is not {RUN} years {YEARS}')
    if committed['climatology']['run'] != RUN or list(committed['climatology']['years']) != YEARS:
        raise SystemExit(f'{DRAINAGE} was not routed on {RUN} years {YEARS}')

    ppd = dr.PPD
    level = float(atlas['sea_level_m'])
    h, _, _ = tp.height_above_geoid(ppd, atlas['grid']['geoid_degree'])
    h = h.astype(np.float32)
    rows, cols = h.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / ppd
    lon = ((np.arange(cols) + 0.5) / ppd + 180.0) % 360.0 - 180.0
    cell_km = MOON_RADIUS / 1000.0 * np.radians(1.0 / ppd)
    cell_area = (cell_km * 1e3) ** 2 * np.cos(np.radians(lat))        # m2, by row
    below = h < level
    body = label_waters(below)
    seas, body_area = name_seas(body, cell_area, lat, lon, atlas['bodies'])
    flat_body = body.ravel()
    sea_of = {v['label']: k for k, v in seas.items()}
    hollows_km2 = float(body_area.sum() - sum(body_area[v['label']] for v in seas.values())) / 1e6

    runoff = dr.runoff_field(clim, rows, cols, None)
    precip = dr.climate_on_grid(clim, 'pr_mm_day', rows, cols)
    lake_net = precip - dr.lake_evaporation(clim, lat)[:, None].astype(np.float32)
    del precip

    outlets = dict(every_water=below,
                   six_seas=np.isin(body, [v['label'] for v in seas.values()]),
                   two_seas=np.isin(body, [seas[s]['label'] for s in TWO]))
    RUNS.mkdir(parents=True, exist_ok=True)
    rules = {}
    for rule, outlet in outlets.items():
        r = dr.route(h, outlet, lat, cell_km, runoff, lake_net)
        q, out = r['discharge'], r['outlet']
        mouths = np.flatnonzero(out & (q > 0))
        basin = r['drained_m2'][mouths] / 1e6
        total = float(q[out].sum())
        per_sea = {}
        for name in seas:
            if rule == 'two_seas' and name not in TWO:
                continue
            here = np.flatnonzero(flat_body[mouths] == seas[name]['label'])
            qs = q[mouths[here]]
            listed = []
            for j in here[np.argsort(qs)[::-1][:TOP]]:
                c = int(mouths[j])
                place = (float(lat[c // cols]), float(lon[c % cols]))
                t_lat, t_lon = to_turned(np.array([place[0]]), np.array([place[1]]), rotation)
                listed.append(dict(today=[round(place[0], 2), round(place[1], 2)],
                                   turned=[round(float(t_lat[0]), 2), round(float((t_lon[0] + 180.0) % 360.0 - 180.0), 2)],
                                   discharge_m3s=round(float(q[c]), 1), basin_km2=round(float(basin[j]))))
            per_sea[name] = dict(inflow_m3s=round(float(qs.sum()), 1), inflow_km3_yr=round(float(qs.sum()) * KM3_YR, 1),
                                 mouths_of_1_m3s=int((qs >= 1.0).sum()), mouths_of_100_m3s=int((qs >= 100.0).sum()),
                                 largest=listed)
        elsewhere = total - sum(v['inflow_m3s'] for v in per_sea.values())
        rules[rule] = dict(rule=RULES[rule], to_the_outlets_m3s=round(total, 1), to_other_waters_m3s=round(elsewhere, 1),
                           seas=per_sea)
        if rule == 'every_water':
            # The hollows, other than the six seas, that take a large river in at the common level
            taken = np.bincount(flat_body[mouths], weights=q[mouths], minlength=body_area.size)
            taken[[v['label'] for v in seas.values()]] = 0.0
            large_hollows = [int(k) for k in np.argsort(taken)[::-1] if taken[k] >= HOLLOW_M3S]
            hollow_inflow = {k: float(taken[k]) for k in large_hollows}
        else:
            heights = h.ravel()
            hollows = []
            for k in large_hollows:
                cells = np.flatnonzero(flat_body == k)
                place = centroid(lat[cells // cols], lon[cells % cols], r['area'][cells])
                hollows.append(dict(centroid_today=[round(place[0], 2), round(place[1], 2)],
                                    inflow_at_the_common_level_m3s=round(hollow_inflow[k], 1),
                                    **filling(r, cells, heights, level)))
            rules[rule]['hollows_filling'] = hollows
        if rule == 'two_seas':
            rules[rule]['smaller_seas_filling'] = {
                name: filling(r, np.flatnonzero(flat_body == seas[name]['label']), h.ravel(), level)
                for name in seas if name not in TWO}
        np.savez_compressed(RUNS / f'mouths_{rule}.npz', mouth_cell=mouths.astype(np.int32),
                            mouth_discharge_m3s=q[mouths].astype(np.float32),
                            mouth_water=flat_body[mouths].astype(np.int32), mouth_basin_km2=basin.astype(np.float32))
        del r
        print(f'{rule}: {total:,.0f} m3/s to the outlets, {time.time() - started:.0f} s', flush=True)

    reproduced = rules['every_water']['to_the_outlets_m3s']
    if round(reproduced) != committed['rivers']['discharge_to_sea_m3s']:
        raise SystemExit(f'The routing gives {reproduced} m3/s to the seas; {DRAINAGE} holds '
                         f"{committed['rivers']['discharge_to_sea_m3s']}")
    return dict(
        schema=SCHEMA,
        producer=dict(files={p: sha256(ROOT / p) for p in CODE}),
        inputs={ATLAS: sha256(ROOT / ATLAS), DRAINAGE: sha256(ROOT / DRAINAGE), FRAME: sha256(ROOT / FRAME)},
        ignored_inputs={CLIMATOLOGY: sha256(ROOT / CLIMATOLOGY)},
        climatology=dict(run=RUN, years=YEARS),
        evidence=('A calculation: geography\'s steady annual routing of the design years\' runoff (ExoPlaSim at T21, '
                  'rain less evaporation over land) over LOLA\'s 16 pixel/degree heights, repeated with three choices '
                  'of the waters that take rivers in. The first reproduces geography\'s drainage product. The rules '
                  'bracket where the rivers end; the climate, the lakes\' open-water balance and the sea level stay '
                  'the design years\' in all three, so none is a prediction of the levels the waters would find. The '
                  'model\'s rain differs from the cloud-resolving boxes\' by a factor the climate lane carries as a '
                  'range.'),
        reading_rule=('Discharges in m3/s and km3 a year (a Julian year). A mouth is a cell of the receiving water '
                      'that a channel enters; its basin is the land draining through it. Places are [latitude, '
                      'longitude east] in degrees on today\'s Moon and on the turned Moon of the frame record. '
                      'to_other_waters_m3s is the discharge reaching outlet waters other than the seas listed. '
                      'hollows_filling lists the hollows, other than the six seas, that take 5,000 m3/s or more of '
                      'river water at the common level, and smaller_seas_filling the four smaller seas in the two_seas '
                      'run: the height above the common level at which each basin spills (the level of the depressions '
                      'under its cells, as the area-weighted median with its lowest and highest values), its area at '
                      'the common level, the area that spill level floods and the volume of water between the two '
                      'levels.'),
        sea_level_m=level,
        seas={k: dict(area_km2=round(v['area_km2']), atlas_area_km2=round(v['atlas_area_km2']),
                      centroid_today=[round(v['centroid'][0], 2), round(v['centroid'][1], 2)]) for k, v in seas.items()},
        hollows_below_the_common_level_km2=round(hollows_km2),
        reproduces=dict(product=DRAINAGE, discharge_to_sea_m3s=committed['rivers']['discharge_to_sea_m3s']),
        rules=rules)


def main(argv=None) -> int:
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
