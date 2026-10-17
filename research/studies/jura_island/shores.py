"""Jura Island's shores: how much of the island lies near the sea, how steep its shores rise, and a screen of the belts
where roots reach the fresh sea's water table under the shore; python -m research.studies.jura_island.shores.

Heights come from LOLA at its native 16 px/deg less the degree-200 GRAIL geoid, as geography.coastal_grid takes them,
against the 28% atlas's sea level. The belt screen sizes the belts from assumed reaches and supply limits; it maps none.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, label

from geography import fetch_inputs
from geography.topography import Grid, geoid
from shared.constants import MOON_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.jura-shores/1'
OUT = HERE / 'results/shores.json'
INPUTS = ('geography/inputs/ldem_16.img', 'geography/inputs/jggrx_0420a_sha.tab', 'geography/results/atlas.json')
PPD = 16
FRAME = dict(south=36.0, north=54.0, west=309.0, east=339.0)    # degrees: box_jura's domain with a margin
CENTRE = (45.0, 323.875)                                         # box_jura's centre (degrees north, east)
BANDS_M = (0, 10, 25, 50, 100, 250, 500, 1000, 5000)
WITHIN_KM = (2, 5, 10, 20)
RINGS_KM = ((0.0, 2.5), (2.5, 4.5), (4.5, 7.0))                 # the first three pixel rings inland
# The belt screen's cases: the height above the sea from which roots and capillary rise still draw on the water table
# (3-5 m of roots and 12-23 m of capillary reach at lunar gravity, geography's shore screen on domain/resources), and
# the widest belt the ground can feed on a flat shore (the central transmissivity at a metre a year, to sandy ground)
BELTS = dict(low=dict(reach_m=15.0, supply_km=0.56), high=dict(reach_m=25.0, supply_km=2.5))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def terrain():
    """Heights above the atlas's sea level (m) on LOLA's native pixel centres over the frame, and the sea level."""
    manifest = fetch_inputs.manifest()
    for name in ('ldem_16.img', 'jggrx_0420a_sha.tab'):
        entry = next(f for f in manifest['files'] if f['name'] == name)
        if digest(fetch_inputs.path(name)) != entry['sha256']:
            raise ValueError(f'Terrain source bytes differ: {name}')
    lat = FRAME['south'] + (np.arange(round((FRAME['north'] - FRAME['south']) * PPD)) + 0.5) / PPD
    lon = FRAME['west'] + (np.arange(round((FRAME['east'] - FRAME['west']) * PPD)) + 0.5) / PPD
    rows = ((90 - lat) * PPD - 0.5).round().astype(int)
    cols = (lon * PPD - 0.5).round().astype(int)
    raw = np.memmap(fetch_inputs.path('ldem_16.img'), dtype='<i2', mode='r', shape=(180 * PPD, 360 * PPD))
    ground = np.asarray(raw[np.ix_(rows, cols)], dtype=float) * 0.5 - geoid(Grid(lat, lon, PPD), max_degree=200)
    level = float(json.loads((ROOT / 'geography/results/atlas.json').read_text())['sea_level_m'])
    return lat, lon, ground - level, level


def measure(lat, lon, height, ppd=PPD, centre=CENTRE):
    """The island's numbers from heights above the sea (m, rows by latitude, columns by longitude, in degrees): the
    island is the connected land holding the highest ground within 4 degrees of latitude and 6 of longitude of the
    centre."""
    land = height > 0
    sea = ~land
    parts, count = label(land, structure=np.ones((3, 3)))
    near = (np.abs(lat[:, None] - centre[0]) < 4) & (np.abs(lon[None, :] - centre[1]) < 6)
    top = np.unravel_index(np.argmax(np.where(near & land, height, -np.inf)), height.shape)
    island = parts == parts[top]
    dy = 2 * np.pi * MOON_RADIUS / 1000 / 360 / ppd                       # km per pixel north-south
    dx = dy * np.cos(np.radians(lat))[:, None] * np.ones_like(height)    # km per pixel east-west
    area = dx * dy
    total = float(area[island].sum())
    h = height
    # distance from the sea, with the east-west spacing taken at the centre's latitude
    distance = distance_transform_edt(land, sampling=(dy, dy * np.cos(np.radians(centre[0]))))
    bands = {f'{lo}-{hi} m': float(area[island & (h >= lo) & (h < hi)].sum()) for lo, hi in zip(BANDS_M, BANDS_M[1:])}
    within = []
    for km in WITHIN_KM:
        m = island & (distance <= km)
        a = float(area[m].sum())
        within.append(dict(within_km=km, area_km2=round(a), share=round(a / total, 4),
                           median_height_m=round(float(np.median(h[m])))))
    rings = []
    for lo, hi in RINGS_KM:
        q = np.percentile(h[island & (distance > lo) & (distance <= hi)], [10, 25, 50, 75, 90])
        rings.append(dict(from_km=lo, to_km=hi, **dict(zip(('p10', 'p25', 'median', 'p75', 'p90'), map(round, q)))))
    # coast: island pixel edges facing the sea; an east-west neighbour shares a north-south edge, and the reverse
    ew = np.zeros_like(h)
    ns = np.zeros_like(h)
    ew[:, 1:] += island[:, 1:] & sea[:, :-1]
    ew[:, :-1] += island[:, :-1] & sea[:, 1:]
    ns[1:, :] += island[1:, :] & sea[:-1, :]
    ns[:-1, :] += island[:-1, :] & sea[1:, :]
    edges_km = ew * dy + ns * dx
    staircase = float(edges_km.sum())
    # belt screen: each coastal pixel's shore rises its mean height over half its width; the belt reaches up to the
    # case's height, no wider than the ground can feed
    coastal = edges_km > 0
    half = 0.25 * (dx + dy)
    slope = np.maximum(h, 0.5) / 1000.0 / half
    belts = {}
    for name, case in BELTS.items():
        width = np.minimum(case['reach_m'] / 1000.0 / slope, case['supply_km'])
        belt = float((edges_km * width)[coastal].sum() * np.pi / 4)
        belts[name] = dict(case, area_km2=round(belt), share=round(belt / total, 4))
    others = sorted(((float(area[parts == k].sum()), float(h[parts == k].max())) for k in range(1, count + 1)
                     if k != parts[top]), reverse=True)
    rows, cols = lat[island.any(1)], lon[island.any(0)]
    return dict(
        island=dict(area_km2=round(total), summit_m=round(float(h[top])), summit_lat_deg=round(float(lat[top[0]]), 3),
                    summit_lon_deg=round(float(lon[top[1]]), 3),
                    lat_range_deg=[round(float(rows.min()), 2), round(float(rows.max()), 2)],
                    lon_range_deg=[round(float(cols.min()), 2), round(float(cols.max()), 2)],
                    mean_height_m=round(float((h * area)[island].sum() / total))),
        heights_km2={k: round(v) for k, v in bands.items()},
        heights_share={k: round(v / total, 4) for k, v in bands.items()},
        below_reach_km2={f'{case["reach_m"]:g} m': round(float(area[island & (h < case['reach_m'])].sum()))
                         for case in BELTS.values()},
        near_sea=within, inland_heights_m=rings,
        coast_staircase_km=round(staircase), coast_km=round(staircase * np.pi / 4),
        belts=belts,
        islets=dict(count=len(others), largest=[dict(area_km2=round(a), top_m=round(t)) for a, t in others[:12]]))


def run():
    lat, lon, height, level = terrain()
    files = [Path(__file__), ROOT / 'geography/topography.py', ROOT / 'geography/fetch_inputs.py']
    product = dict(
        schema=SCHEMA,
        producer=dict(study='research/studies/jura_island', files={str(p.relative_to(ROOT)): digest(p) for p in files},
                      inputs={p: digest(ROOT / p) for p in INPUTS}, constants=constants_used(files)),
        evidence=('Heights from LOLA at 16 px/deg, about 1.9 km north-south by 1.3 km east-west at 45 N, less the '
                  'degree-200 GRAIL geoid and the 28% atlas\'s sea level; the island is the connected land holding '
                  'the highest ground near box_jura\'s centre. The belt screen takes as assumed ranges how high above '
                  'the sea roots and capillary rise reach the water table and how wide a belt the ground can feed on '
                  'a flat shore, with each coastal pixel\'s mean height over half its width as its shore slope. It '
                  'sizes the belts and maps none; the shore study at tens of metres would.'),
        reading_rule=('Areas in km2 and heights in m above the atlas\'s sea level. Distances in km from the nearest '
                      'sea pixel, the east-west spacing taken at 45 N. The coast length counts pixel edges, a '
                      'staircase; coast_km is pi/4 of it, the mean for a coast of random direction. Shares are of the '
                      'island.'),
        frame_deg=FRAME, pixels_per_degree=PPD, sea_level_m=level,
        **measure(lat, lon, height))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    return product


if __name__ == '__main__':
    p = run()
    print(json.dumps({k: p[k] for k in ('island', 'heights_share', 'near_sea', 'coast_km', 'belts')}, indent=1))
