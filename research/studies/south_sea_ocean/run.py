"""The South Pole-Aitken Sea as an Earth-like ocean: its salt, its water balance on the design years, the rivers
against it, what its floor takes and how salt changes its mixing.

    python -m research.studies.south_sea_ocean.run     # results/south_sea_ocean.json; seconds

On 15 November 2026 the author asked to test the viability and the cost and benefit of making the south sea an
Earth-like ocean, the more so if it connects to brackish and fresh waterways, which would also give eels and other
fish that hatch in salt water their sea. This module answers the parts that the project's products already settle:

- the salt: what the sea holds at seawater's make-up and at two weaker strengths, ion by ion;
- the water: rain and evaporation over the sea's own cells in the design years, and so the fresh water it needs to
  hold its level;
- the rivers: what reaches it under the three rules of inflows.py, against that need;
- the floor: the water the crust takes from under the sea (geography's stage 1) and the salt that goes with it;
- the mixing: how fast evaporation salts a surface layer, set beside the cooling that sinks the fresh seas' water.

Where the salt comes from is resources' question and is not answered here. The sea's circulation, its brackish
reaches and its oxygen are reasoned from these numbers in the README and need the ocean models before they are results.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from climate.gcm.turned import to_turned
from research.studies.sea_chemistry import run as sc
from research.studies.seawater_lakes.run import practical_from_reference, seawater_density
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'south_sea_ocean.json'
SCHEMA = 'terluna.research.south-sea-ocean/1'
ATLAS = 'geography/results/atlas.json'
ATLAS_GRID = 'geography/products/atlas_28pct_4ppd.npz'
FRAME = 'climate/results/gcm/turned_surface_8d71fd350ca4c6cd.json'
RUN = 'A28R_301k_design_l15_snow_s938_turned_refit_precession'
CLIMATOLOGY = f'climate/gcm/products/climatology_{RUN}.npz'
YEARS = [10, 19]
INFLOWS = 'research/studies/south_sea_ocean/results/inflows.json'
UPTAKE = 'geography/results/groundwater_uptake.json'
SEA_CHEMISTRY = 'research/studies/sea_chemistry/results/sea_chemistry.json'
IMPORTS = ('research/studies/sea_chemistry/run.py', 'research/studies/seawater_lakes/run.py', 'climate/gcm/turned.py',
           'shared/constants.json')

SEA_LABEL, SEA_NAME = 1417, 'south_pole_aitken'
SEAS = {433: 'near_side', 1417: 'south_pole_aitken', 874: 'smythii_marginis', 213: 'humboldtianum',
        1296: 'orientale', 714: 'moscoviense'}
T_REF_C = 29.5                          # the seawater-lakes study's sea-surface temperature, kept so the salt agrees
STRENGTHS = dict(ocean=1.0, half=0.5, fifth=0.2)   # shares of seawater's salt: 35.2, 17.6 and 7.0 g/kg
# climate/gcm README, "The gap carried as a range": against the model's rain over the turned sea box, the
# cloud-resolving box rains 0.75 of it under the model's imposed motion and twice it making its own. Applied here to
# the rain on the sea as a sensitivity; the rivers are not rescaled.
RAIN_FACTORS = (0.75, 2.0)
LAYERS_M = (50.0, 100.0, 200.0)         # surface layers for the salting rate (the GCM's slab is 50 m)
UPTAKE_YEARS = ('100', '1000', '10000', '100000', '1e+06')
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90))
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY
MM_DAY_TO_M_YR = JULIAN_YEAR_DAYS / 1000.0
R_M = MOON_RADIUS
# Grams of each ion in a mole (the carbon cycle's and the sea chemistry's values).
ION_G_MOL = dict(na=sc.cc.M_NA, k=sc.cc.M_K, mg=sc.cc.M_MG, ca=sc.cc.M_CA, cl=sc.M_CL, so4=sc.M_SO4,
                 alkalinity=sc.M_HCO3)


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# --- The salt ---------------------------------------------------------------------------------------------------------
def ion_grams_per_kg() -> dict:
    """Seawater's major ions in g per kg of seawater (Millero et al. 2008 as the sea chemistry carries it), the
    alkalinity as bicarbonate, and the rest of the 35.165 g/kg (bromide, strontium, boron, fluoride) as 'minor'."""
    grams = {ion: sc.SEAWATER[ion] * ION_G_MOL[ion] for ion in sc.SEAWATER}
    grams['minor'] = sc.SEAWATER_G_KG - sum(grams.values())
    return grams


def salt(volume_km3: float, share: float, t_c: float = T_REF_C) -> dict:
    """The salt in a volume of water at a share of seawater's strength: its salinity, density, mass and ions (kg)."""
    salinity = share * sc.SEAWATER_G_KG
    density = seawater_density(practical_from_reference(salinity), t_c)
    water_kg = volume_km3 * 1e9 * density
    return dict(salinity_g_kg=salinity, density_kg_m3=density, sea_mass_kg=water_kg,
                salt_kg=water_kg * salinity * 1e-3,
                ions_kg={ion: water_kg * share * g * 1e-3 for ion, g in ion_grams_per_kg().items()})


# --- The sea on the model's grid --------------------------------------------------------------------------------------
def cell_areas_km2(lat_deg, lon_deg) -> np.ndarray:
    """Areas of a regular latitude-longitude grid's cells, km2, by row."""
    dlat = np.radians(abs(lat_deg[1] - lat_deg[0]))
    dlon = np.radians(abs(lon_deg[1] - lon_deg[0]))
    return R_M ** 2 * np.cos(np.radians(lat_deg)) * dlat * dlon / 1e6


def on_model_grid(sea_mask, lat_deg, lon_deg, rotation, glat, glon) -> np.ndarray:
    """A sea's area (km2) in each cell of the model's grid on the turned Moon: every atlas cell of the sea is turned
    with the frame's rotation and given to the model cell whose centre is nearest."""
    rows, cols = np.nonzero(sea_mask)
    t_lat, t_lon = to_turned(lat_deg[rows], lon_deg[cols], rotation)
    i = np.abs(t_lat[:, None] - np.asarray(glat, float)[None, :]).argmin(axis=1)
    step = 360.0 / len(glon)
    j = np.rint(np.mod(t_lon - float(glon[0]), 360.0) / step).astype(int) % len(glon)
    area = np.zeros((len(glat), len(glon)))
    np.add.at(area, (i, j), cell_areas_km2(lat_deg, lon_deg)[rows])
    return area


def water_balance(area_km2, clim, sea_area_km2: float) -> dict:
    """Rain and evaporation over a sea in the model's own water cells, area-weighted, and the net evaporation over
    the atlas's area of the sea."""
    wet = area_km2 * (clim['lsm'] < 0.5)
    w = wet.sum()
    rain = float((clim['pr_mm_day'] * wet).sum() / w)
    evaporation = float((clim['evap_mm_day'] * wet).sum() / w)
    surface_c = float((clim['ts_k'] * wet).sum() / w) - 273.15
    rows = []
    for r in np.flatnonzero(wet.sum(axis=1) > 0):
        a = wet[r].sum()
        rows.append(dict(latitude_turned=round(float(clim['lat'][r]), 2), area_km2=round(float(a)),
                         rain_mm_day=round(float((clim['pr_mm_day'][r] * wet[r]).sum() / a), 3),
                         evaporation_mm_day=round(float((clim['evap_mm_day'][r] * wet[r]).sum() / a), 3),
                         surface_c=round(float((clim['ts_k'][r] * wet[r]).sum() / a) - 273.15, 2)))
    net = (evaporation - rain) * MM_DAY_TO_M_YR
    net_range = tuple(sorted((evaporation - f * rain) * MM_DAY_TO_M_YR for f in RAIN_FACTORS))
    return dict(share_of_area_in_model_water_cells=float(w / area_km2.sum()), model_cells=int((area_km2 > 0).sum()),
                model_water_cells=int((wet > 0).sum()), rain_mm_day=rain, evaporation_mm_day=evaporation,
                surface_c=surface_c, net_evaporation_m_yr=net, net_evaporation_km3_yr=net * sea_area_km2 * 1e-3,
                net_evaporation_m_yr_over_the_rain_range=net_range,
                net_evaporation_km3_yr_over_the_rain_range=tuple(x * sea_area_km2 * 1e-3 for x in net_range),
                rows=rows)


# --- The rivers against the need --------------------------------------------------------------------------------------
def rivers_against_need(inflow_km3_yr: float, need_km3_yr: float, net_m_yr: float, area_km2: float) -> dict:
    """A sea's river inflow beside its net evaporation: the share met, the surplus (negative for a deficit) and the
    area whose net evaporation the inflow would balance if the sea found its own size at the same loss per area."""
    surplus = inflow_km3_yr - need_km3_yr
    return dict(inflow_km3_yr=inflow_km3_yr, share_of_need=inflow_km3_yr / need_km3_yr, surplus_km3_yr=surplus,
                surplus_m3s=surplus * 1e9 / YEAR_S, level_change_m_yr=surplus / (area_km2 * 1e-3),
                balanced_area_km2=inflow_km3_yr / (net_m_yr * 1e-3),
                balanced_area_share_of_sea=inflow_km3_yr / (net_m_yr * 1e-3) / area_km2)


# --- The floor ---------------------------------------------------------------------------------------------------------
def band_areas_km2(sea_mask, lat_deg, lon_deg) -> dict:
    """The sea's area in each band of latitude on today's Moon (both hemispheres together, as stage 1's bands are)."""
    row_area = cell_areas_km2(lat_deg, lon_deg) * sea_mask.sum(axis=1)
    return {f'{a}-{b}': float(row_area[(np.abs(lat_deg) >= a) & (np.abs(lat_deg) < b)].sum()) for a, b in BANDS}


def floor_uptake(bands_km2: dict, uptake: dict, ocean: dict, volume_km3: float) -> dict:
    """Water the crust holds from under the sea's floor by each year of stage 1, and the salt that goes with it when
    the sea is at seawater's strength. Stage 1's columns are fresh water under each band's mean sea depth."""
    out = {}
    for case in ('slow', 'central', 'fast'):
        rows, before, year_before = {}, 0.0, 0.0
        for year in UPTAKE_YEARS:
            km3 = sum(km2 * uptake['cases'][case]['bands'][band]['columns']['sea_floor']['by_year'][year]['held_m'] * 1e-3
                      for band, km2 in bands_km2.items())
            salt_kg = km3 * 1e9 * ocean['density_kg_m3'] * ocean['salinity_g_kg'] * 1e-3
            rows[year] = dict(water_km3=km3, share_of_the_sea=km3 / volume_km3, salt_kg=salt_kg,
                              water_km3_yr_since_last=(km3 - before) / (float(year) - year_before),
                              salt_kg_yr_since_last=(km3 - before) / (float(year) - year_before) * 1e9
                              * ocean['density_kg_m3'] * ocean['salinity_g_kg'] * 1e-3)
            before, year_before = km3, float(year)
        out[case] = rows
    return out


# --- The mixing --------------------------------------------------------------------------------------------------------
def density_slopes(salinity_g_kg: float, t_c: float) -> dict:
    """Seawater's density change per kelvin of cooling and per g/kg of salt (EOS-80, centred differences)."""
    s = practical_from_reference(salinity_g_kg)
    per_k = (seawater_density(s, t_c - 0.5) - seawater_density(s, t_c + 0.5))
    per_g = (seawater_density(practical_from_reference(salinity_g_kg + 0.5), t_c)
             - seawater_density(practical_from_reference(salinity_g_kg - 0.5), t_c))
    return dict(kg_m3_per_k_of_cooling=per_k, kg_m3_per_g_kg=per_g, k_of_cooling_per_g_kg=per_g / per_k)


def salting(net_m_yr: float, salinity_g_kg: float, slopes: dict) -> dict:
    """How fast net evaporation salts a surface layer that exchanges nothing, and the cooling of equal weight."""
    out = {}
    for h in LAYERS_M:
        rate = salinity_g_kg * net_m_yr / h
        out[f'{h:g}_m'] = dict(g_kg_a_year=rate, kg_m3_a_year=rate * slopes['kg_m3_per_g_kg'],
                               equal_to_cooling_k_a_year=rate * slopes['k_of_cooling_per_g_kg'])
    return out


def build() -> dict:
    atlas = json.loads((ROOT / ATLAS).read_text())
    grid = np.load(ROOT / ATLAS_GRID)
    frame = json.loads((ROOT / FRAME).read_text())
    rotation = np.array(frame['frame']['rotation_today_to_turned'])
    clim = dict(np.load(ROOT / CLIMATOLOGY))
    meta = json.loads(str(clim.pop('metadata')))
    if meta.get('run') != RUN or list(meta.get('years')) != YEARS:
        raise SystemExit(f'{CLIMATOLOGY} is not {RUN} years {YEARS}')
    inflows = json.loads((ROOT / INFLOWS).read_text())
    if inflows['climatology'] != dict(run=RUN, years=YEARS):
        raise SystemExit(f'{INFLOWS} was not routed on {RUN} years {YEARS}')
    uptake = json.loads((ROOT / UPTAKE).read_text())
    sea_chem = json.loads((ROOT / SEA_CHEMISTRY).read_text())
    held = [r for r in sea_chem['runs']
            if (r['corner'], r['delivery'], r['sinks'], r['sediment']) == ('central', 'fresh', 'managed', 'kept_open')]
    fresh_g_kg = float(held[0]['rows'][-1]['salinity_g_kg'])

    labels, lat, lon = grid['water_label'], grid['lat_deg'], grid['lon_deg']
    bodies = {b['label']: b for b in atlas['bodies']}
    body = bodies[SEA_LABEL]
    area, volume = float(body['area_km2']), float(body['volume_km3'])
    mask = labels == SEA_LABEL
    t_lat, _ = to_turned(lat[np.nonzero(mask)[0]], lon[np.nonzero(mask)[1]], rotation)
    sea = dict(atlas_label=SEA_LABEL, area_km2=area, share_of_the_moon=float(body['share']), volume_km3=volume,
               mean_depth_m=float(body['mean_depth_m']), greatest_depth_m=float(body['max_depth_m']),
               latitudes_today=[float(lat[mask.any(axis=1)].min()), float(lat[mask.any(axis=1)].max())],
               latitudes_turned=[float(t_lat.min()), float(t_lat.max())])

    # The salt
    strengths = {name: salt(volume, share) for name, share in STRENGTHS.items()}
    ocean = strengths['ocean']
    all_seas = sum(float(bodies[k]['volume_km3']) for k in SEAS)
    salt_block = dict(ion_g_per_kg=ion_grams_per_kg(), temperature_c=T_REF_C, strengths=strengths,
                      salt_kg_per_km3=ocean['salt_kg'] / volume,
                      share_of_the_six_seas_volume=volume / all_seas)

    # The water: every large sea on the same reading, the south sea with its rows
    balances = {}
    for label, name in SEAS.items():
        a = on_model_grid(labels == label, lat, lon, rotation, clim['lat'], clim['lon'])
        b = water_balance(a, clim, float(bodies[label]['area_km2']))
        if name != SEA_NAME:
            b.pop('rows')
        balances[name] = b
    own = balances[SEA_NAME]
    need = own['net_evaporation_km3_yr']
    need_range = own['net_evaporation_km3_yr_over_the_rain_range']

    # The rivers against the need, for each rule and each sea the rule lists
    rivers = {}
    for rule, r in inflows['rules'].items():
        rivers[rule] = {name: rivers_against_need(s['inflow_km3_yr'], balances[name]['net_evaporation_km3_yr'],
                                                  balances[name]['net_evaporation_m_yr'],
                                                  float(bodies[{v: k for k, v in SEAS.items()}[name]]['area_km2']))
                        for name, s in r['seas'].items()}
    for rule in rivers:
        rivers[rule][SEA_NAME]['share_of_need_over_the_rain_range'] = tuple(
            sorted(inflows['rules'][rule]['seas'][SEA_NAME]['inflow_km3_yr'] / n for n in need_range))

    # How soon the waters that take rivers in at the common level would fill to their own spill levels
    km3_yr = YEAR_S / 1e9
    filling = dict(
        hollows=[dict(centroid_today=x['centroid_today'], inflow_km3_yr=x['inflow_at_the_common_level_m3s'] * km3_yr,
                      spill_level_above_common_m=x['spill_level_above_common_m'],
                      volume_above_the_common_level_km3=x['volume_above_the_common_level_km3'],
                      years_to_fill=x['volume_above_the_common_level_km3']
                      / (x['inflow_at_the_common_level_m3s'] * km3_yr))
                 for x in inflows['rules']['six_seas']['hollows_filling']],
        smaller_seas={name: dict(spill_level_above_common_m=x['spill_level_above_common_m'],
                                 volume_above_the_common_level_km3=x['volume_above_the_common_level_km3'],
                                 surplus_km3_yr=rivers['six_seas'][name]['surplus_km3_yr'],
                                 years_to_fill=(x['volume_above_the_common_level_km3']
                                                / rivers['six_seas'][name]['surplus_km3_yr']
                                                if rivers['six_seas'][name]['surplus_km3_yr'] > 0 else None))
                      for name, x in inflows['rules']['two_seas']['smaller_seas_filling'].items()})

    # Held at its level with its rivers matched to its loss: the water's turnover and the load the rivers bring
    held_sea = dict(fresh_inflow_km3_yr=need, fresh_inflow_m3s=need * 1e9 / YEAR_S,
                    fresh_inflow_km3_yr_over_the_rain_range=need_range,
                    fresh_inflow_m3s_over_the_rain_range=tuple(x * 1e9 / YEAR_S for x in need_range),
                    years_to_evaporate_its_volume=volume / need,
                    inflow_salinity_g_kg=fresh_g_kg,
                    inflow_load_kg_yr=need * 1e9 * seawater_density(practical_from_reference(fresh_g_kg), T_REF_C)
                    * fresh_g_kg * 1e-3)
    held_sea['inflow_load_share_of_salt_a_year'] = held_sea['inflow_load_kg_yr'] / ocean['salt_kg']

    # The floor
    bands = band_areas_km2(mask, lat, lon)
    floor = dict(sea_area_by_band_km2=bands, by_case=floor_uptake(bands, uptake, ocean, volume))

    # The mixing
    slopes = density_slopes(ocean['salinity_g_kg'], T_REF_C)
    rows = own['rows']
    warm = max(rows, key=lambda r: r['surface_c'])
    cold = min(rows, key=lambda r: r['surface_c'])
    mixing = dict(density_slopes=slopes, salting_of_a_closed_surface_layer=salting(own['net_evaporation_m_yr'],
                                                                                   ocean['salinity_g_kg'], slopes),
                  warmest_row=dict(latitude_turned=warm['latitude_turned'], surface_c=warm['surface_c']),
                  coldest_row=dict(latitude_turned=cold['latitude_turned'], surface_c=cold['surface_c']),
                  cold_end_denser_kg_m3=(warm['surface_c'] - cold['surface_c']) * slopes['kg_m3_per_k_of_cooling'],
                  salt_that_would_offset_it_g_kg=(warm['surface_c'] - cold['surface_c']) / slopes['k_of_cooling_per_g_kg'],
                  fresh_water_lighter_than_the_sea_kg_m3=ocean['density_kg_m3']
                  - seawater_density(practical_from_reference(fresh_g_kg), T_REF_C))

    return dict(
        schema=SCHEMA,
        producer=dict(file='research/studies/south_sea_ocean/run.py', sha256=sha256(Path(__file__)),
                      imports={p: sha256(ROOT / p) for p in IMPORTS}),
        inputs={p: sha256(ROOT / p) for p in (ATLAS, FRAME, INFLOWS, UPTAKE, SEA_CHEMISTRY)},
        ignored_inputs={p: sha256(ROOT / p) for p in (ATLAS_GRID, CLIMATOLOGY)},
        climatology=dict(run=RUN, years=YEARS),
        evidence=('Calculations on project products; no model was run here. The sea is the atlas\'s body at the 28% '
                  'water share. Rain, evaporation and surface temperature are ten-year means of ExoPlaSim at T21 with '
                  'a prescribed 50 m slab of fresh water, read at the model\'s water cells under the sea; the model\'s '
                  'rain carries the range the climate lane states against the cloud-resolving boxes, applied to the '
                  'rain on the sea alone. The rivers are inflows.py\'s three routings of the same years. The floor\'s '
                  'uptake is geography\'s stage 1, computed for fresh water. Densities are UNESCO\'s EOS-80 at one '
                  'atmosphere. Salt lowers the vapour pressure of water by about 2%, which this reading leaves out. '
                  'Nothing here is validated against a lunar observation.'),
        reading_rule=('Masses in kg, volumes in km3, rates per Julian year. Ranges are (low, high). A positive net '
                      'evaporation is water the sea loses to the air beyond its rain; the fresh inflow that holds its '
                      'level equals it. In rivers, surplus is inflow less net evaporation and balanced_area is the '
                      'sea area whose net evaporation the inflow would meet at the same loss per area. In filling, '
                      'years_to_fill is the water between the common level and a basin\'s spill level over the '
                      'river water it takes (a hollow) or over its surplus under the six_seas rule (a smaller sea; '
                      'null where it has none), without the evaporation of the growing lake. Salinities '
                      'are reference salinities in g per kg of sea water. Latitudes of rows are the model\'s, on the '
                      'turned Moon.'),
        sea=sea, salt=salt_block, water=balances, rivers=rivers, filling=filling, held=held_sea, floor=floor,
        mixing=mixing)


def main(argv=None) -> int:
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
