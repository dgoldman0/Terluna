"""Cities on the Open Moon's land: the lava tubes and the highlands.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.land_cities.run   # -> results/land_cities.json

A coupled screen on committed products. Part one sets the lunar lava-tube candidates of the literature against the
28% atlas and the water the Open Moon brings, solves the water table about each land candidate in two dimensions,
and sizes what keeping a tube dry takes, what it holds, and its roles before the air and after it. Part two reads the
highlands from the atlas and the climate products, sizes how tall buildings get at a sixth of Earth's gravity with the
engineering lattice model, and gives the land, water and travel the highland people take. Grid sums, closed-form
physics, one small groundwater solve per land candidate and arithmetic on cited values.

Main's corrected per-place product is read in-tree, its hash checked (BOUND).
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from geography import topography as tp
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY
from shared.provenance import constants_used

from research.studies.land_cities import grid, highlands as hl, tubes as tb

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.land-cities/1'
OUT = HERE / 'results' / 'land_cities.json'

INPUTS = {
    'geography/results/atlas.json': 'terluna.geography.atlas/1',
    'geography/results/drainage.json': 'terluna.geography.drainage/1',
    'geography/results/groundwater.json': 'terluna.geography.groundwater/1',
    'climate/results/crm/comfort_flat_rings.json': 'terluna.climate.crm-comfort/1',
    'climate/results/crm/highland_box_highland_own_height.json': 'terluna.climate.crm-highland/1',
    'climate/results/gcm/global_winds_A28_dim5_moon.json': 'terluna.climate.gcm-global-winds/1',
    'atmosphere/electricity/results/occurrence_A28_dim5_moon.json': 'terluna.atmosphere.lightning-occurrence/1',
    'illumination/aerosol/results/haze.json': 'terluna.illumination.haze/1',
    'research/studies/sea_appearance/results/lighting_calendar.json': 'terluna.research.sea-lighting-calendar/1',
    'research/studies/summit_tower/results/summit_tower.json': 'terluna.research.summit-tower/1',
    'research/studies/summit_tower/results/summit_tower_form.json': 'terluna.research.summit-tower-form/1',
    'research/studies/sky_fleet/results/sky_fleet.json': 'terluna.research.sky-fleet/1',
    'research/studies/inhabited_volume/results/inhabited_volume.json': 'terluna.research.inhabited-volume/1',
    'research/studies/ways_of_living/results/ways_of_living.json': 'terluna.research.ways-of-living/1',
    'research/studies/daily_life/results/daily_life.json': 'terluna.research.daily-life/1',
}
# Grid products and external inputs kept out of Git (regenerated or fetched by their domains).
IGNORED = {
    'geography/products/atlas_28pct_4ppd.npz': 'terluna.geography.atlas-grid/1',
    'geography/products/drainage_28pct_16ppd.npz': 'terluna.geography.drainage/1',
    'climate/gcm/products/climatology_A28_dim5_moon.npz': 'terluna.climate.gcm-climatology/1',
    'geography/inputs/ldem_16.img': None,
    'geography/inputs/jggrx_0420a_sha.tab': None,
    'geography/inputs/MOON_nomenclature_center_pts.zip': None,
}
BOUND = {
    'places': dict(path='research/studies/joint_synthesis/results/places.json',
                   schema='terluna.research.joint-synthesis-places/1', branch='main',
                   commit=None,
                   sha256='9707e4d1c7332ba352b90dbf77c66fde080ec5e0fe24e57b54067c61ecdcfc8d',
                   note="Main's corrected per-place product (each place's own sky, lightning and haze)."),
}
FILES = ('research/studies/land_cities/run.py', 'research/studies/land_cities/tubes.py',
         'research/studies/land_cities/groundwater.py', 'research/studies/land_cities/highlands.py',
         'research/studies/land_cities/grid.py', 'engineering/towers/lattice.py', 'geography/topography.py',
         'geography/nomenclature.py')
GEOID_DEGREE = 200          # the atlas's geoid


class MissingInput(RuntimeError):
    """An ignored grid product, an external input or a bound input is not available here."""


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(path) -> str:
    return sha256(path)[:16]


def load_json(rel, schema):
    data = json.loads((ROOT / rel).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema in {rel}: {data.get("schema")}')
    return data


def load_npz(rel, schema):
    path = ROOT / rel
    if not path.exists():
        raise MissingInput(f'{rel} is not present; regenerate it in its domain')
    with np.load(path, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files}
    meta = json.loads(str(data.pop('metadata')))
    if meta.get('schema') != schema:
        raise ValueError(f'Unexpected schema in {rel}: {meta.get("schema")}')
    return data, meta


def bound_bytes(key) -> bytes:
    """The bound input's bytes: the in-tree file when its hash matches, as it does since the joint integration;
    else the Git object at a pinned commit."""
    b = BOUND[key]
    local = ROOT / b['path']
    raw = local.read_bytes() if local.exists() else b''
    if hashlib.sha256(raw).hexdigest() != b['sha256']:
        try:
            raw = subprocess.run(['git', '-C', str(ROOT), 'show', f"{b['commit']}:{b['path']}"], capture_output=True,
                                 check=True).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            raise MissingInput(f"{b['path']} at {b['branch']} {b['commit'][:7]} is not readable here") from error
    if hashlib.sha256(raw).hexdigest() != b['sha256']:
        raise ValueError(f"{b['path']} at {b['commit'][:7]} does not match its pinned hash")
    return raw


def load_bound(key):
    data = json.loads(bound_bytes(key))
    if data.get('schema') != BOUND[key]['schema']:
        raise ValueError(f'Unexpected schema in bound input {key}')
    return data


def missing_inputs():
    missing = [p for p in IGNORED if not (ROOT / p).exists()]
    for key in BOUND:
        try:
            bound_bytes(key)
        except MissingInput as error:
            missing.append(str(error))
    return missing


def r(x, n=3):
    return tb.r(x, n)


# ---------------------------------------------------------------------------------------------------------------

def lola_points(cands):
    """LOLA's 16 px/deg radius grid and the degree-200 geoid (with Earth's static tide and rotation, as the atlas)
    at the candidates."""
    lola, _ = tp.read_ldem(16)
    lat = np.array([c['lat'] for c in cands], float)
    lon = np.array([c['lon'] % 360.0 for c in cands], float)
    geo = np.diag(tp.geoid(tp.Grid(lat, lon, 16), GEOID_DEGREE))
    return lola, geo


def results() -> dict:
    data = {rel: load_json(rel, schema) for rel, schema in INPUTS.items()}
    atlas_grid, atlas_meta = load_npz('geography/products/atlas_28pct_4ppd.npz',
                                      IGNORED['geography/products/atlas_28pct_4ppd.npz'])
    drainage, _ = load_npz('geography/products/drainage_28pct_16ppd.npz',
                           IGNORED['geography/products/drainage_28pct_16ppd.npz'])
    clim, _ = load_npz('climate/gcm/products/climatology_A28_dim5_moon.npz',
                       IGNORED['climate/gcm/products/climatology_A28_dim5_moon.npz'])
    for rel in ('geography/inputs/ldem_16.img', 'geography/inputs/jggrx_0420a_sha.tab',
                'geography/inputs/MOON_nomenclature_center_pts.zip'):
        if not (ROOT / rel).exists():
            raise MissingInput(f'{rel} is not present; fetch it with geography/fetch_inputs.py')
    places = load_bound('places')['places']
    at = grid.build_atlas(atlas_grid, atlas_meta, drainage)
    if abs(at.sea_level_m - data['geography/results/atlas.json']['sea_level_m']) > 0.05:
        raise ValueError('The atlas grid and atlas product disagree on the sea level')
    rain = grid.regrid(clim['pr_mm_day'], clim['lat'], clim['lon'], at.lat, at.lon)
    evap = grid.regrid(clim['evap_mm_day'], clim['lat'], clim['lon'], at.lat, at.lon)
    lake16, level16 = grid.lake_mask16(drainage)
    return dict(tubes=tube_part(data, at, rain, evap, lake16, level16),
                highlands=highland_part(data, at, rain, evap, places))


def highland_part(data, at, rain, evap, places):
    occ = data['atmosphere/electricity/results/occurrence_A28_dim5_moon.json']
    flashes = grid.regrid(np.asarray(occ['central_flash_density'], float), occ['lat_deg'], occ['lon_deg'],
                          at.lat, at.lon)
    tower = data['research/studies/summit_tower/results/summit_tower.json']
    wol = data['research/studies/ways_of_living/results/ways_of_living.json']
    settings, per_billion = wol['settings'], wol['per_billion_common']
    iv = data['research/studies/inhabited_volume/results/inhabited_volume.json']
    efficiency = iv['assumptions']['pump_efficiency']
    towns, districts = settings['highland_towns'], settings['elevated_districts']
    home = dict(km2=towns['home_dry_land_km2'] + districts['home_dry_land_km2'],
                far_side_share=towns['physical']['earth_up_shares']['never'])
    land = hl.land(at, rain, evap, flashes, home)
    summit = dict(lat_deg=tower['site']['lat_deg'], lon_deg=tower['site']['lon_deg'])
    rise = hl.summit_rise(at, summit, rain, flashes)
    heat = [per_billion['heat_on_moon_tw'][k] for k in ('efficient', 'european', 'high_income')]
    design = settings['elevated_districts']['design']
    return dict(
        land=land, summit_rise=rise, wettest_and_stormiest_ground=hl.wettest_and_stormiest(at, summit, rain, flashes),
        comfort=hl.comfort(data['climate/results/crm/comfort_flat_rings.json'],
                           data['climate/results/crm/highland_box_highland_own_height.json'], land),
        lightning_at_places={k: places[k]['lightning'] for k in ('summit port', 'highland box', 'polar plain')},
        towers=hl.towers(data['climate/results/gcm/global_winds_A28_dim5_moon.json'], tower, heat),
        ground_people=hl.ground_people(land, rise),
        water=hl.water(at, land, per_billion['domestic_water_gross_km3_year'], efficiency),
        access=hl.access(data['research/studies/sky_fleet/results/sky_fleet.json']),
        scattered=hl.scattered(settings),
        tower_district=dict(people_per_km2=design['tower_district_people_per_km2'],
                            building_t_per_person=design['building_t_per_person'],
                            summit_tower_steel_t_per_person=design['summit_tower_steel_t_per_person_at_70_m2']))


NAMED_RILLES = ('Vallis Schroteri', 'Rima Sharp', 'Rima Marius', 'Rima Hadley', 'Rimae Prinz',
                'Rimae Sulpicius Gallus')
HEADS_M = (15.0, 50.0, 100.0, 500.0, 1000.0, 2000.0)
ROOFS_M = (5.0, 10.0, 20.0, 40.0, 100.0)
LAG_DEPTHS_M = (20.0, 50.0, 100.0, 150.0)
DAYLIGHT_LUX = (250.0, 500.0)           # Brown et al. 2022's daytime floor of 250 lux, and 500 lux for desk work
DELIVERED_SHARE = (0.3, 0.5)            # heliostat and light-guide train, mirror to floor (a design guess)
LATITUDES = (0.0, 15.0, 30.0, 45.0)
OUTDOOR_AIR_L_S = (5.0, 10.0)           # per person: ASHRAE 62.1's office composite to Approved Document F's 10
FAN_PRESSURE_PA = (300.0, 1000.0)       # through shafts, filters and ducts (a design range)
FAN_EFFICIENCY = 0.6
CITIES = (1e6, 7e6)                     # people in a tube city
STACK = dict(height_m=100.0, delta_t_k=5.0)
ISA_SEA_LEVEL = dict(density_kg_m3=1.225, temperature_k=288.15)
FUNNELS_M = (150.0, 300.0)              # pit funnels across, from the atlas's mouths
BREACH_HEADS_M = (50.0, 500.0, 2000.0)
# The author's estimate, under a roof like the 20-25 m lava sheet over the Marius Hills pit (Haruyama et al. 2009).
REFERENCE_TUBE = dict(width_m=1000.0, height_m=50.0, length_km=50.0, roof_m=25.0)


def gazetteer_rilles():
    from geography import nomenclature as nm
    found = {f['name']: f for f in nm.features() if f['name'] in NAMED_RILLES}
    return [dict(key=n.lower().replace(' ', '_'), name=n, lat=round(found[n]['lat'], 3),
                 lon=round(found[n]['lon'] % 360.0, 3), length_km=round(found[n]['diameter_km'], 1),
                 notes=tb.NAMED_RILLE_NOTES.get(n)) for n in NAMED_RILLES if n in found]


def candidates():
    """Every candidate to place on the atlas, with its type."""
    out = [dict(p, type='mare pit') for p in tb.PITS]
    out += [dict(key=x['key'], name=x['name'], lat=x['lat'], lon=x['lon'], type='roofed rille segments')
            for x in tb.RILLES]
    out += [dict(x, type='collapse chain') for x in tb.CHAINS]
    out += [dict(x, type='gravity candidate') for x in tb.GRAVITY]
    out += [dict(x, type='radar void') for x in tb.RADAR]
    out += [dict(x, type='named rille (centre)') for x in gazetteer_rilles()]
    return out


def place_candidates(data, at, rain, evap, lake16, level16):
    """The candidates on the atlas, with each land candidate's water table and each one's floor against it."""
    cands = candidates()
    lola, geo = lola_points(cands)
    water_km, nearest = at.km_to_water()
    runoff = rain - evap
    located = tb.locate(cands, at, lola, geo, lake16, level16, data['geography/results/atlas.json']['bodies'],
                        rain, runoff, water_km, nearest)
    by_key = {c['key']: c for c in cands}
    sites = {}
    for row in located:
        c = by_key[row['key']]
        raw = row.pop('raw')
        row['type'] = c['type']
        floor = tb.floor_depth_m(c)
        row['floor_below_ground_m'] = tb.rr(floor) if floor else None
        if row['class'] in LAND:
            if floor:
                row['floor_above_sea_m'] = tb.rr([raw['ground_m'] - at.sea_level_m - f for f in floor[::-1]], 4)
            place = (round(c['lat'], 4), round(c['lon'] % 360.0, 4))
            if place not in sites:
                sites[place] = tb.groundwater_site(at, c['lat'], c['lon'], rain, runoff)
            row['water_table'], depth = tb.water_table(sites[place], raw['ground_m'], floor)
            if floor:
                row['floor_below_water_table'] = {
                    case: ('yes' if depth[case] < floor[0] else 'partly' if depth[case] < floor[1] else 'no')
                    for case in ('low', 'central', 'high')}
        else:
            row['head_on_floor_m'] = tb.rr([raw['water_depth_m'] + f for f in floor]) if floor else None
    classes = {}
    for row in located:
        classes.setdefault(row['type'], {}).setdefault(row['class'], []).append(row['name'])
    counts = {t: {k: len(v) for k, v in d.items()} for t, d in classes.items()}
    return located, by_key, classes, counts


LAND, UNDER_WATER = ('dry', 'coastal'), ('under the sea', 'under a lake')
STRADDLE_SEA_SHARE = (0.1, 0.9)     # a search box this far under the sea, and no further, straddles the coast


def add(pair, other):
    return [pair[0] + other[0], pair[1] + other[1]]


def capacity_totals(located, by_key):
    """What the sized candidates hold, each candidate's sum over its segments, the sums by kind, and the totals by
    class with each region's tube counted once at its largest estimate, without and with the GRAIL candidates. A
    region whose candidates lie both on land and under water, or whose search box lies 10-90% under the sea,
    straddles the coast."""
    rows = {row['key']: row for row in located}
    sized = []
    for row in located:
        c = by_key[row['key']]
        for case in tb.sizes(c):
            sized.append(dict(key=row['key'], name=row['name'], type=row['type'], cls=row['class'],
                              region=c.get('region') or tb.REGIONS.get(row['key'], row['key']), **case))
    by_candidate = {}
    for s_ in sized:
        e = by_candidate.setdefault(s_['key'], dict(name=s_['name'], type=s_['type'], cls=s_['cls'],
                                                    region=s_['region'], small=[0.0, 0.0], large=[0.0, 0.0]))
        for size in ('small', 'large'):
            e[size] = add(e[size], s_[size]['people_million_at_70_m2'])
    by_type = {}
    for e in by_candidate.values():
        t = by_type.setdefault(e['type'], dict(small=[0.0, 0.0], large=[0.0, 0.0], candidates=0))
        t['small'], t['large'] = add(t['small'], e['small']), add(t['large'], e['large'])
        t['candidates'] += 1
    regions = {}
    for k, e in by_candidate.items():
        regions.setdefault(e['region'], {})[k] = e
    totals = {}
    for label, keep in (('mapped', lambda t: t != 'gravity candidate'), ('with_gravity', lambda t: True)):
        by_class = {}
        for name, reg in regions.items():
            options = {k: v for k, v in reg.items() if keep(v['type'])}
            if not options:
                continue
            k, v = max(options.items(), key=lambda kv: kv[1]['large'][1])
            # Every candidate of the region in this total, sized or not, tells which side of the coast it lies.
            members = [x for x in rows if (by_key[x].get('region') or tb.REGIONS.get(x, x)) == name
                       and keep(rows[x]['type'])]
            sides = {rows[x]['class'] for x in members}
            seas = [rows[x]['box']['sea_share'] for x in members
                    if rows[x].get('box') and rows[x]['box']['sea_share'] is not None]
            straddles = (sides & set(LAND) and sides & set(UNDER_WATER)) or any(
                STRADDLE_SEA_SHARE[0] <= x <= STRADDLE_SEA_SHARE[1] for x in seas)
            cls = 'straddles the coast' if straddles else v['cls']
            t = by_class.setdefault(cls, dict(small=[0.0, 0.0], large=[0.0, 0.0], regions=[]))
            t['small'], t['large'] = add(t['small'], v['small']), add(t['large'], v['large'])
            t['regions'].append(f'{name} ({k})')
        all_small = [sum(v['small'][i] for v in by_class.values()) for i in (0, 1)]
        all_large = [sum(v['large'][i] for v in by_class.values()) for i in (0, 1)]
        totals[label] = dict(by_class={c: dict(small=tb.rr(v['small']), large=tb.rr(v['large']), regions=v['regions'])
                                       for c, v in by_class.items()},
                             all=dict(small=tb.rr(all_small), large=tb.rr(all_large)))
    candidates_out = {k: dict(e, small=tb.rr(e['small']), large=tb.rr(e['large'])) for k, e in by_candidate.items()}
    types_out = {k: dict(e, small=tb.rr(e['small']), large=tb.rr(e['large'])) for k, e in by_type.items()}
    return sized, candidates_out, types_out, totals


def sealing_settings(located):
    """Seepage and pumping for the author's tube above the water table, with the water table at the ground, and
    under 500 m and 2 km of sea."""
    ref = REFERENCE_TUBE
    mid, crown = ref['roof_m'] + ref['height_m'] / 2.0, ref['roof_m'] + ref['height_m']
    out = []
    # Above the water table only the rain that percolates through the roof arrives: at most the recharge over the
    # tube's footprint, at the land candidates' rain range and the recharge shares.
    land_rain = [r_['rain_mm_day'] for r_ in located if r_['class'] in ('dry', 'coastal')]
    rain = [min(land_rain), max(land_rain)]
    footprint = ref['width_m'] * ref['length_km'] * 1e3
    q = [tb.RECHARGE_SHARE[0] * rain[0] / 1e3 / JULIAN_DAY * footprint,
         tb.RECHARGE_SHARE[2] * rain[1] / 1e3 / JULIAN_DAY * footprint]
    out.append(dict(setting='above the water table, percolation only', head_m=0.0, rain_mm_day=tb.rr(rain),
                    percolation_m3_s=tb.rr(q), pumping_mw=tb.rr([tb.pumping_mw(x, crown) for x in q])))
    for label, sea_depth in (('the water table at the ground', 0.0), ('under 500 m of sea', 500.0),
                             ('under 2,000 m of sea', 2000.0)):
        sea = tb.sealing(sea_depth + mid, ref['width_m'], ref['height_m'], ref['length_km'], sea_depth + crown)
        area_m2 = sea['envelope_km2'] * 1e6
        flux = {k: tb.r(v['seepage_m3_s'] / area_m2 * 1e3 * JULIAN_DAY) for k, v in sea['by_ring'].items()}
        per_100m = {k: tb.r(v['seepage_m3_s'] / (ref['length_km'] * 10.0) * 1e3 * 60.0)
                    for k, v in sea['by_ring'].items()}
        out.append(dict(setting=label, **sea, seepage_l_m2_day=flux, seepage_l_min_per_100_m=per_100m))
    return out


def tube_part(data, at, rain, evap, lake16, level16):
    located, by_key, classes, counts = place_candidates(data, at, rain, evap, lake16, level16)
    sized, by_candidate, by_type, totals = capacity_totals(located, by_key)
    sealing = sealing_settings(located)
    ref = REFERENCE_TUBE
    energy = data['research/studies/ways_of_living/results/ways_of_living.json']['per_billion_common']['energy_tw']
    vault = tb.check_estimate()['vault']
    residents = [vault[f'{tb.STOREY_M[1]:g}_m_storeys']['people_million']['70_m2'][0],
                 vault[f'{tb.STOREY_M[0]:g}_m_storeys']['people_million']['70_m2'][1]]
    residents_energy = dict(people_million=residents, kw_per_person=[energy['efficient'], energy['high_income']],
                            gw=tb.rr([residents[0] * energy['efficient'], residents[1] * energy['high_income']]),
                            basis='the vault\'s residents at the ways-of-living energy a person, efficient to '
                                  'high-income (TW a billion is kW a person)')
    haze = data['illumination/aerosol/results/haze.json']['regions']['seas']['cases']['clean']['day']['daylight']
    calendar = data['research/studies/sea_appearance/results/lighting_calendar.json']
    above_air = calendar['sky']['shielded_sunlight_above_air_lux']
    lit = data['research/studies/daily_life/results/daily_life.json']['work']['indoor_daylight_through_the_night']
    winds = data['climate/results/gcm/global_winds_A28_dim5_moon.json']
    lowest_level_m3 = tb.level_widths(ref['width_m'], ref['height_m'], tb.STOREY_M[0])[0] * tb.STOREY_M[0] \
        * ref['length_km'] * 1e3
    first = dict(
        roof_shielding=tb.roof_shielding(ROOFS_M), air_matching_earth=tb.air_matching_earth_column(),
        roof_to_hold_air=[tb.roof_to_hold_air(p) for p in (101.325, 50.0)],
        warming_lag_years={f'{d:g}_m': tb.warming_lag_years(d) for d in LAG_DEPTHS_M},
        temperatures=tb.TEMPERATURES, dose=tb.DOSE,
        limit_msv_year=tb.r(tb.DOSE['open_moon_limit_msv_day'] * JULIAN_YEAR_DAYS))
    moon_air = dict(density_kg_m3=winds['density_mean_kg_m3'][0], temperature_k=winds['temperature_k'][0])
    lasting = dict(
        daylight=[tb.heliostat_share(lat, above_air, haze, lux, DELIVERED_SHARE)
                  for lat in LATITUDES for lux in DAYLIGHT_LUX],
        lamps_kwh_per_person_half_cycle=lit['kwh_per_person'],
        ventilation=[tb.ventilation(n, q, FAN_PRESSURE_PA, FAN_EFFICIENCY) for n in CITIES for q in OUTDOOR_AIR_L_S],
        stack_pa=dict(moon=tb.r(tb.stack_pressure_pa(STACK['height_m'], STACK['delta_t_k'], moon_air['density_kg_m3'],
                                                     moon_air['temperature_k'])),
                      earth=tb.r(tb.stack_pressure_pa(STACK['height_m'], STACK['delta_t_k'],
                                                      ISA_SEA_LEVEL['density_kg_m3'], ISA_SEA_LEVEL['temperature_k'],
                                                      gravity=STANDARD_GRAVITY)),
                      moon_air=moon_air, **STACK),
        floods=dict(
            breach_m3_s={f'{h:g}_m': tb.r(tb.breach_inflow_m3_s(h, 1.0)) for h in BREACH_HEADS_M},
            lowest_level_m3=tb.r(lowest_level_m3),
            lowest_level_fill_years_by_central_seepage={
                s_['setting']: tb.r(tb.fill_years(lowest_level_m3, s_['by_ring']['central']['seepage_m3_s']))
                for s_ in sealing if 'by_ring' in s_},
            storm_rain_mm_h=tb.STORM_RAIN_MM_H,
            funnel_inflow_m3_s=[tb.r(tb.STORM_RAIN_MM_H[1] / 1e3 / 3600.0 * math.pi * (d / 2.0) ** 2)
                                for d in FUNNELS_M]),
        quakes=tb.QUAKES)
    return dict(
        inventory=dict(pits=tb.PITS, rilles=tb.rille_summary(), chains=tb.chain_summary(), gravity=tb.GRAVITY,
                       radar=tb.RADAR, named_rilles=gazetteer_rilles()),
        stability=dict(blair=tb.STABILITY, theinat=tb.STABILITY_THEINAT, sources=tb.STABILITY_SOURCES,
                       checks=tb.stability_checks(),
                       one_km=dict(blair_lithostatic_roof_m=2, theinat_roof_m=tb.theinat_min_roof_m(1000.0))),
        on_the_atlas=dict(rows=located, classes=classes, counts=counts),
        water_table_method=dict(
            solve='steady Dupuit mound in two dimensions on the atlas within about 400 km of each land candidate, the '
                  'seas and lakes at their levels, the ground capping the mound; rise by the Boussinesq equation',
            window_half_km=tb.gw.WINDOW_KM, moon_over_earth_mound=tb.mound_ratio_limits()),
        pressures=[dict(head_m=h, pressure_kpa=tb.r(float(tb.pressure_pa(h)) / 1e3),
                        earth_head_same_pressure_m=tb.r(float(tb.earth_head_m(h)))) for h in HEADS_M],
        ita_membrane_head_moon_m=tb.r(tb.ITA_MEMBRANE_HEAD_EARTH_M / tb.GRAVITY_RATIO),
        membrane_class_heads=dict(earth_m=list(tb.MEMBRANE_CLASS_HEADS_EARTH_M),
                                  moon_m=tb.rr([h / tb.GRAVITY_RATIO for h in tb.MEMBRANE_CLASS_HEADS_EARTH_M])),
        leakage_allowed_l_m2_day=tb.LEAKAGE_ALLOWED_L_M2_DAY,
        norwegian_inflow_l_min_100m=tb.NORWEGIAN_INFLOW_L_MIN_100M,
        sealing=dict(reference_tube=ref, settings=sealing, residents_energy=residents_energy),
        capacity=dict(estimate_check=tb.check_estimate(), candidates=sized, by_candidate=by_candidate,
                      by_type=by_type, totals=totals, grail_size=tb.GRAIL_SIZE,
                      grail_section=dict(levels=tb.GRAIL_SIZE['levels'], height_to_width=list(tb.HEIGHT_TO_WIDTH),
                                         basis='the lowest ten levels in a semi-elliptical vault')),
        first_phase=first, lasting=lasting)


def product() -> dict:
    res = results()
    files = {p: digest(ROOT / p) for p in FILES}
    inputs = {p: digest(ROOT / p) for p in INPUTS}
    ignored = {p: digest(ROOT / p) for p in IGNORED if (ROOT / p).exists()}
    bound = {k: dict(path_at_integration=b['path'], branch=b['branch'], commit=b['commit'], sha256=b['sha256'],
                     schema=b['schema'], note=b['note']) for k, b in BOUND.items()}
    return dict(
        schema=SCHEMA,
        producer=dict(study='land_cities', files=files, inputs=inputs, ignored_inputs=ignored, bound_inputs=bound,
                      constants=constants_used(FILES)),
        evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=assumptions(), **res)


EVIDENCE = ('A literature inventory of lunar lava-tube candidates set against the 28% atlas (LOLA 16 px/deg less the '
            'degree-200 geoid); a two-dimensional Dupuit groundwater solve on the atlas about each land candidate, '
            'checked against the analytic mound; closed-form seepage, conduction and view-factor physics at lunar '
            'gravity with Earth basalt and tunnelling values carried over; and grid sums of the atlas, the design '
            'climatology and the lightning map over the highlands, with the engineering lattice model sized for '
            'building towers. Tube sizes rest on remote sensing (pit imaging, rille morphology, gravity and radar) '
            'and are bounds and estimates; no tube interior beyond tens of metres from a pit has been seen. The '
            'groundwater solve takes uniform rock over a flat base at the atlas\'s 7.6 km cells; nothing here is a '
            'structural or radiation-transport model.')
READING_RULE = ('Heights in m or km above the 28% sea level unless named; longitudes east. Areas km2, people in '
                'millions or billions as named, power MW or GW as named, water km3 a year. Ranges are [low, high] '
                'or low/central/high cases as labelled. Water tables are the atlas cell\'s steady level from the '
                'two-dimensional solve, capped by the cell\'s mean ground, set against the candidate\'s own ground; '
                'heads are metres of water at lunar gravity, with the Earth head that gives the same pressure beside '
                'them.')


def assumptions():
    return dict(
        water_table=dict(recharge_share_of_rain=tb.RECHARGE_SHARE, k_earth_m_s=tb.K_EARTH,
                         basalt_thickness_m=tb.BASALT_THICKNESS_M, specific_yield=tb.SPECIFIC_YIELD,
                         rule='low: least recharge, most permeable and thickest basalt; high: the reverse',
                         base='a flat base the basalt thickness below the window\'s lowest water; Gong et al. measure '
                              'the thickness down from the mare surface, so under high ground the basalt may end '
                              'above that base and the saturated column then runs into the crust below',
                         rise_start='the water table at the sea\'s level under the land when the rain begins'),
        sealing=dict(grouted_k_earth_m_s=tb.GROUTED_K_EARTH, grouted_ring_m=tb.GROUTED_RING_M,
                     pump_efficiency=tb.PUMP_EFFICIENCY),
        capacity=dict(storey_m=tb.STOREY_M, usable_share=tb.USABLE_SHARE, floor_m2_per_person=tb.FLOOR_M2_PER_PERSON,
                      section='a semi-elliptical vault over a flat floor; a level counts where a full storey of '
                              'headroom remains'),
        rock=dict(basalt_bulk_density_kg_m3=tb.BASALT_BULK_DENSITY, rock_diffusivity_m2_s=tb.ROCK_DIFFUSIVITY),
        towers=dict(storey_m=hl.STOREY_M, plate_m=hl.PLATE_M, floor_mass_kg_m2=hl.FLOOR_MASS_KG_M2,
                    moon_ground_km=hl.TOWER_GROUND_KM, earth_design_gust_m_s=hl.EARTH_DESIGN_GUST),
        people=dict(densities_per_km2=hl.DENSITIES, ground_billion=hl.GROUND_BILLIONS,
                    highland_shares=hl.HIGHLAND_SHARES, kept_shares=hl.KEPT_SHARES),
        water_recovery=hl.RECOVERY)


def dumps(obj) -> str:
    return json.dumps(obj, indent=1, ensure_ascii=False) + '\n'


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product()))
    print(OUT)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
