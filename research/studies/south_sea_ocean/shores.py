"""The south sea's shores: how long they are, what a fresh sea feeds along them and what a salt sea would wet.

    ulimit -v 9000000; python -m research.studies.south_sea_ocean.shores    # results/shores.json; a minute, 2 GB

The coast-water study (research/studies/coast_water) found that a fresh sea waters a belt of ground a few hundred
metres wide along its desert shores: plants draw the sea's water through the ground wherever the water table lies
within their reach. Around a salt sea that water is sea water. This module takes the South Pole-Aitken Sea's own shore
off the same 16 pixel/degree shores (coast_water.shores) and counts, with that study's rules:

- the shore's length, on the mainland and on the sea's islands, by the design years' rain behind it;
- the belts a fresh sea feeds along its desert shores (coast_water.belts: steady flow from the water's edge under
  the plants' draw), their area and the water they draw; around a salt sea the same draw would bring sea water under
  the roots, so the salt that water carries is given with it;
- the ground the monthly high water covers (half the typical monthly range of geography's tide over each shore
  pixel's slope), which a salt sea wets with salt water each month;
- the small waters below the common sea level within 12 km of the sea, which the atlas holds at the sea's level and
  which stand in the same ground as the sea.

The rain is the design years' (the corrected run's years 10-19, carried to today's frame), read at the model's land
cell nearest each shore pixel and, for comparison, at the cell over it. Desert is land under 0.5 mm of rain a day, as
the climate lane and the coast-water study have it. The plants' draw in the deserts is the coast-water study's at
301 K, which came from the earlier design case on the unturned Moon; it is carried here as a transfer. The tide is
geography's on the unturned Moon's geoid.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from climate.gcm.turned import to_turned
from research.studies.coast_water import belts as B
from research.studies.coast_water import load as L
from research.studies.coast_water import shores as S
from research.studies.south_sea_ocean import inflows as I
from research.studies.south_sea_ocean import run as R
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'shores.json'
SCHEMA = 'terluna.research.south-sea-shores/1'
COAST_WATER = 'research/studies/coast_water/results/coast_water.json'
TIDE_GRID = 'geography/products/tides_28pct_4ppd.npz'
CLIMATOLOGY = I.CLIMATOLOGY                  # the design years carried to today's frame
CODE = ('research/studies/south_sea_ocean/shores.py', 'research/studies/south_sea_ocean/inflows.py',
        'research/studies/south_sea_ocean/run.py', 'research/studies/coast_water/shores.py',
        'research/studies/coast_water/belts.py', 'research/studies/coast_water/load.py', 'climate/gcm/turned.py')
SEA = 'south_pole_aitken'
DESERT_MM_DAY = 0.5                          # land under this much rain a day is desert (climate/gcm README)
RAIN_CLASSES = (('desert', 0.0, 0.5), ('half_dry', 0.5, 1.0), ('moist', 1.0, 2.0), ('wet', 2.0, np.inf))
ROCK_SALT_KG_M3 = 2165.0                     # halite
NEAR_KM = 12.0                               # small waters counted: within this distance of the sea
BANDS = R.BANDS
DAY_S = JULIAN_DAY
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)


def facing(index, target, cols):
    """Which of the land pixels at the flat indices have a pixel of the target water beside them (north, south, east
    or west; the east and west sides join across the seam, and nothing lies beyond a pole)."""
    rows = target.shape[0]
    r, c = np.divmod(np.asarray(index, np.int64), cols)
    out = np.zeros(r.size, bool)
    for dr, dc in ((0, 1), (0, -1), (-1, 0), (1, 0)):
        rr = r + dr
        inside = (rr >= 0) & (rr < rows)
        out |= inside & target[np.clip(rr, 0, rows - 1), (c + dc) % cols]
    return out


def rain_at(lat_deg, lon_deg, clim):
    """The design years' rain (mm a day) at each place: at the model cell over it and at the nearest land cell."""
    glat, glon = clim['lat'].astype(float), clim['lon'].astype(float)
    pr, land = clim['pr_mm_day'].astype(float), clim['lsm'] > 0.5
    i = np.abs(np.asarray(lat_deg, float)[:, None] - glat[None, :]).argmin(axis=1)
    step = 360.0 / glon.size
    j = np.rint(np.mod(np.asarray(lon_deg, float) - glon[0], 360.0) / step).astype(int) % glon.size
    li, lj = np.nonzero(land)
    _, k = cKDTree(unit(glat[li], glon[lj])).query(unit(lat_deg, lon_deg))
    return pr[i, j], pr[li[k], lj[k]], land[i, j]


def rain_class(rain_mm_day):
    """Each rain's class index in RAIN_CLASSES."""
    edges = np.array([c[1] for c in RAIN_CLASSES[1:]])
    return np.searchsorted(edges, np.asarray(rain_mm_day, float), side='right')


def high_water_strip(amplitude_m, rise_m, across_km):
    """The width (m) of ground the monthly high water covers on a shore pixel: the tide's amplitude over the slope
    from the shoreline to the pixel's height at half its depth, and no wider than the pixel."""
    across = np.asarray(across_km, float) * 1000.0
    rise = np.maximum(np.nan_to_num(np.asarray(rise_m, float), nan=0.0), 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = np.where(across > 0, rise / (0.5 * across), 0.0)
        w = np.where(slope > 0, np.asarray(amplitude_m, float) / slope, across)
    return np.minimum(w, across)


def tide_amplitudes(lat_deg, lon_deg):
    """Half the typical monthly range (m) at the nearest sea node of geography's tide grid."""
    with np.load(ROOT / TIDE_GRID) as z:
        lat, lon = z['lat_deg'][z['rows']], z['lon_deg'][z['cols']]
        node_range = z['monthly_range_median_m'].astype(float)
    _, node = cKDTree(unit(lat, lon)).query(unit(lat_deg, lon_deg))
    return 0.5 * node_range[node]


def percentiles(x, w, qs=(0.1, 0.5, 0.9)):
    x, w = np.asarray(x, float), np.asarray(w, float)
    keep = w > 0
    if not keep.any():
        return [None] * len(qs)
    order = np.argsort(x[keep])
    cum = np.cumsum(w[keep][order])
    return [float(x[keep][order][min(np.searchsorted(cum, q * cum[-1]), order.size - 1)]) for q in qs]


def belts_for(coast, inland, uptake_m_s, transmissivity, reach):
    """Belt area (km2) of each shore pixel and of each pixel further inland under one case of ground and reach, as
    the coast-water study counts them: a strip along each edge facing the sea, capped at the pixel's own area, and
    the share of each inland pixel the belt on the steady slope from the shoreline covers."""
    c = uptake_m_s / (2.0 * transmissivity)
    area, width, length = B.first_ring(coast['edge'], coast['across'], coast['rise'], c, reach)
    area = np.minimum(area, coast['area'])
    share = B.inland(inland['distance'], inland['step_rows'], inland['step_cols'], inland['height'], c, reach)
    return area, width, length, share * inland['area']


def build() -> dict:
    atlas = json.loads((ROOT / I.ATLAS).read_text())
    cw = json.loads((ROOT / COAST_WATER).read_text())
    frame = json.loads((ROOT / I.FRAME).read_text())['frame']
    rotation = np.array(frame['rotation_today_to_turned'], float)
    clim = np.load(ROOT / CLIMATOLOGY, allow_pickle=True)
    sh = S.build()
    rows, cols = S.ROWS, S.COLS
    lat = sh['lat']
    lon = (np.arange(cols) + 0.5) / S.PPD
    wclass = sh['wclass'].reshape(rows, cols)
    body = I.label_waters(wclass == S.SEA)
    named, _ = I.name_seas(body, sh['row_area'] * 1e6, lat, lon, atlas['bodies'])
    sea = body == named[SEA]['label']
    sea_km2 = float((sea.sum(axis=1) * sh['row_area']).sum())

    # the shore pixels beside the sea, and the pixels further inland whose nearest water is the sea
    co, px = sh['coast'], sh['pixels']
    beside = facing(co['index'], sea, cols) & (co['edge_km'][:, 0] > 0)
    two_seas = int((beside & facing(co['index'], (wclass == S.SEA) & ~sea, cols)).sum())
    cr, cc = np.divmod(co['index'][beside].astype(np.int64), cols)
    coast = dict(row=cr, col=cc, edge=co['edge_km'][beside, 0].astype(float),
                 across=co['across_km'][beside, 0].astype(float), rise=co['rise_m'][beside, 0].astype(float),
                 island=co['island'][beside], area=sh['row_area'][cr])
    coastal = np.zeros(rows * cols, bool)
    coastal[co['index']] = True
    pr_, pc_ = np.divmod(px['index'].astype(np.int64), cols)
    nearest_r, nearest_c = pr_ + px['step_rows'], (pc_ + px['step_cols']) % cols
    far = ~coastal[px['index']] & (px['water_class'] == S.SEA) & sea[nearest_r, nearest_c]
    inland = dict(row=pr_[far], col=pc_[far], distance=px['distance_km'][far].astype(float),
                  step_rows=px['step_rows'][far], step_cols=px['step_cols'][far],
                  height=px['height_m'][far].astype(float), island=px['island'][far], area=sh['row_area'][pr_[far]])

    length = B.STAIRCASE * coast['edge']
    c_lat, c_lon = lat[coast['row']], lon[coast['col']]
    rain_over, rain_land, over_land = rain_at(c_lat, c_lon, clim)
    cls = rain_class(rain_land)
    t_lat, _ = to_turned(c_lat, c_lon, rotation)
    i_rain = rain_at(lat[inland['row']], lon[inland['col']], clim)[1]
    desert, i_desert = rain_land < DESERT_MM_DAY, i_rain < DESERT_MM_DAY

    def by_class(values):
        return {name: float(values[cls == k].sum()) for k, (name, _, _) in enumerate(RAIN_CLASSES)}

    shore = dict(
        length_km=float(length.sum()), mainland_km=float(length[~coast['island']].sum()),
        island_km=float(length[coast['island']].sum()), pixels=int(length.size), pixels_beside_two_seas=two_seas,
        by_rain_km=by_class(length),
        by_rain_km_mainland=by_class(np.where(coast['island'], 0.0, length)),
        by_rain_at_the_cell_over_it_km={name: float(length[rain_class(rain_over) == k].sum())
                                        for k, (name, _, _) in enumerate(RAIN_CLASSES)},
        share_under_the_models_water_cells=float(length[~over_land].sum() / length.sum()),
        by_turned_latitude_km={f'{a}-{b}S': dict(
            all=float(length[(-t_lat >= a) & (-t_lat < b)].sum()),
            desert=float(length[(-t_lat >= a) & (-t_lat < b) & desert].sum())) for a, b in BANDS},
        moon_km=dict(seas=float(B.STAIRCASE * co['edge_km'][:, 0].astype(float).sum()),
                     every_water=float(B.STAIRCASE * co['edge_km'].astype(float).sum())))

    # the belts a fresh sea feeds along the desert shores, case by case as the coast-water study has them
    uptake_mm = cw['demand']['desert_uptake_mm_day']
    uptake = uptake_mm / 1000.0 / DAY_S
    cases = cw['belts']['cases']
    salt_kg_m3 = R.salt(1.0, 1.0)['salt_kg'] / 1e9          # the salt in a cubic metre of sea water
    salt_density = float(R.seawater_density(R.practical_from_reference(R.sc.SEAWATER_G_KG), R.T_REF_C))
    fresh_density = float(R.seawater_density(0.0, R.T_REF_C))
    belts = {}
    for name, case in cases.items():
        area, width, strip_length, far_area = belts_for(coast, inland, uptake, case['transmissivity_m2_s'],
                                                        case['reach_m'])
        km2 = float(area[desert].sum() + far_area[i_desert].sum())
        draw_m3 = km2 * 1e6 * uptake * YEAR_S
        belts[name] = dict(
            **case, belt_km2=km2, beside_the_water_km2=float(area[desert].sum()),
            further_inland_km2=float(far_area[i_desert].sum()),
            on_the_islands_km2=float(area[desert & coast['island']].sum()
                                     + far_area[i_desert & inland['island']].sum()),
            belt_per_km_of_shore_m=km2 / max(float(length[desert].sum()), 1e-9) * 1000.0,
            strip_width_m_p10_p50_p90=percentiles(width[desert], length[desert]),
            draw_km3_yr=draw_m3 / 1e9, draw_per_km_of_shore_m3_yr=draw_m3 / max(float(length[desert].sum()), 1e-9),
            salt_in_that_draw_mt_yr=draw_m3 * salt_kg_m3 / 1e9,
            salt_per_m2_of_belt_kg_yr=uptake * YEAR_S * salt_kg_m3)
    central = sorted(n for n in belts if n.startswith('central'))
    # what that draw does to a root zone, by the coast-water study's rule (soil water at 4 dS/m is saline)
    ions = R.ion_grams_per_kg()
    chloride_mol_m3 = salt_kg_m3 * ions['cl'] / sum(ions.values()) / R.sc.M_CL * 1000.0
    draw_m_yr = uptake * YEAR_S
    pore_m = L.CEMENT_POROSITY * L.ROOT_ZONE_M
    root_zone = dict(
        draw_m_yr=draw_m_yr, chloride_mol_per_m3_of_sea_water=chloride_mol_m3,
        saline_at_mol_m3=L.SALINE_MOL_M3, soil_water_m=pore_m,
        days_to_a_saline_root_zone=L.SALINE_MOL_M3 * pore_m / (draw_m_yr * chloride_mol_m3) * JULIAN_YEAR_DAYS,
        water_to_flush_below_the_roots_m_yr=draw_m_yr * chloride_mol_m3 / L.SALINE_MOL_M3,
        rock_salt_mm_yr=draw_m_yr * salt_kg_m3 / ROCK_SALT_KG_M3 * 1000.0,
        fresh_water_reaches_below_the_sea_times_its_height_above=fresh_density / (salt_density - fresh_density))
    # the coast-water study's own desert sea and island coasts, Moon-wide, on its 301 K mask (unturned Moon)
    theirs = cw['belts']['deserts']['rain_mask.every_cell']
    moon_wide = {n: dict(shore_km=sum(theirs[n]['by_class'][k]['shore_km'] for k in ('sea_coast', 'island_coast')),
                         belt_km2=sum(theirs[n]['by_class'][k]['belt_km2'] for k in ('sea_coast', 'island_coast')))
                 for n in central}

    # the monthly tide on the shore and the ground its high water covers
    amp = tide_amplitudes(c_lat, c_lon)
    strip = high_water_strip(amp, coast['rise'], coast['across'])
    flooded = length * strip / 1000.0
    tide = dict(
        typical_monthly_range_m_p10_p50_p90=percentiles(2.0 * amp, length),
        high_water_strip_m_p10_p50_p90=percentiles(strip, length),
        covered_at_high_water_km2=float(flooded.sum()), by_rain_km2=by_class(flooded),
        on_the_islands_km2=float(flooded[coast['island']].sum()),
        share_of_shore_capped_at_the_pixel=float(length[strip >= coast['across'] * 1000.0 - 1e-6].sum()
                                                 / length.sum()))

    # the small waters below the common level near the sea
    small = wclass == S.LAKE_COMMON
    dist, _ = S.distances(sea, lat, sh['dy'])
    near = small & (dist <= NEAR_KM)
    lakes = S.merge_labels(ndimage.label(small, structure=np.ones((3, 3), bool))[0])
    touched = np.unique(lakes[near])
    touched = touched[touched > 0]
    lake_area = S.body_areas(lakes, sh['row_area'])
    nr, nc = np.nonzero(near)
    near_rain = rain_at(lat[nr], lon[nc], clim)[1]
    near_area = sh['row_area'][nr]
    small_waters = dict(
        within_km=NEAR_KM, bodies=int(touched.size), area_of_those_bodies_km2=float(lake_area[touched].sum()),
        area_within_reach_km2=float(near_area.sum()),
        area_within_reach_in_desert_km2=float(near_area[near_rain < DESERT_MM_DAY].sum()),
        largest_km2=sorted((float(a) for a in lake_area[touched]), reverse=True)[:5])

    inputs = (I.ATLAS, COAST_WATER, I.FRAME, CLIMATOLOGY, TIDE_GRID, S.DRAINAGE) + S.LOLA
    return dict(
        schema=SCHEMA,
        evidence=('Counts on the coast-water study\'s 16 pixel/degree shores (LOLA above the GRAIL geoid at the '
                  '28% atlas\'s sea level) with that study\'s belt rules, the design years\' rain at T21 and '
                  'geography\'s equilibrium monthly tide. Nothing is modelled beyond those products. Shore lengths '
                  'carry the staircase factor pi/4; at 100 m the coast-water study found belts 1.0-2.7 times '
                  'narrower than this resolution gives, and the same caution holds for the ground the high water '
                  'covers. The plants\' draw in the deserts is the coast-water study\'s at 301 K on the earlier '
                  'design case; the tide is on the unturned Moon\'s geoid.'),
        reading_rule=('Lengths in km, areas in km2, widths in m. by_rain_* split the shore by the design years\' '
                      'rain at the model\'s land cell nearest each shore pixel: desert under 0.5 mm a day, half_dry '
                      '0.5-1, moist 1-2, wet 2 or more. fresh_sea_belts are the belts a fresh sea feeds along the '
                      'desert shore, one entry for each of the coast-water study\'s cases of ground and reach; '
                      'salt_in_that_draw is what the same draw would carry from a sea at seawater\'s strength. '
                      'tide.covered_at_high_water_km2 is the ground between the mean shoreline and the typical '
                      'month\'s high water.'),
        sea=dict(name=SEA, area_km2=sea_km2, atlas_area_km2=named[SEA]['atlas_area_km2'],
                 sea_level_m=sh['level'], pixels_per_degree=S.PPD),
        climate=dict(run=I.RUN, years=I.YEARS, frame='today', desert_mm_day=DESERT_MM_DAY),
        shore=shore,
        fresh_sea_belts=dict(desert_shore_km=float(length[desert].sum()), uptake_mm_day=uptake_mm,
                             uptake_source='coast_water.json demand.desert_uptake_mm_day (301 K, earlier design case)',
                             salt_kg_per_m3_of_sea_water=salt_kg_m3, root_zone=root_zone, cases=belts,
                             central_km2=[belts[n]['belt_km2'] for n in central],
                             coast_water_desert_sea_and_island_coasts=moon_wide),
        tide=tide, small_waters_below_the_common_level=small_waters,
        inputs={p: sha256(ROOT / p) for p in inputs},
        producer=dict(files={p: sha256(ROOT / p) for p in CODE}))


def main(argv=None) -> int:
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    s, b, t = out['shore'], out['fresh_sea_belts'], out['tide']
    print(f"shore {s['length_km']:,.0f} km ({s['island_km']:,.0f} on islands); by rain {s['by_rain_km']}")
    print(f"desert shore {b['desert_shore_km']:,.0f} km; belts {b['central_km2']} km2 (central crust, 15 and 25 m)")
    print(f"high water covers {t['covered_at_high_water_km2']:,.0f} km2; range {t['typical_monthly_range_m_p10_p50_p90']}")
    print(out['small_waters_below_the_common_level'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
