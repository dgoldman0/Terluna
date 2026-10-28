"""Food systems for the Open Moon: what each costs per billion residents, how each fits the lunar cycle, and mixes.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.food_systems.run   # -> results/food_systems.json

A coupled screen. Each food system's costs per billion lunar residents on the lunar design's plant-rich diet come from
the project's products (the crop model's yields and night stores, the design climate's rain, fog and wind, the atlas,
the heat ledger) and from primary literature values listed in parameters.py with their sources in sources.json. The
systems are set together in mixes for 4, 6 and 8 billion lunar residents. Grid sums, closed-form physics and
arithmetic on cited values only: no crop, climate, hydrology or ecosystem model runs here.

Seven inputs came from other branches; since the joint integration they are read in-tree, their hashes checked (BOUND).
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

from research.studies.food_systems import mixes as mixes_mod
from research.studies.food_systems import parameters as par
from research.studies.food_systems import systems as sys_mod
from research.studies.food_systems import zones as zones_mod
from research.studies.food_systems.zones import MissingInput

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.food-systems/1'
OUT = HERE / 'results' / 'food_systems.json'
CODE = ['run.py', 'zones.py', 'parameters.py', 'systems.py', 'mixes.py']

# Committed products read here: path -> schema.
INPUTS = {
    'research/studies/ways_of_living/results/ways_of_living.json': 'terluna.research.ways-of-living/1',
    'provisioning/results/population.json': 'terluna.provisioning.population/1',
    'research/studies/daily_life/results/daily_life.json': 'terluna.research.daily-life/1',
    'biosphere/canopy/results/canopy.json': 'terluna.biosphere.canopy-photosynthesis/1',
    'biosphere/canopy/results/plant.json': 'terluna.biosphere.plant-carbon-cycle/1',
    'biosphere/canopy/results/fruit.json': 'terluna.biosphere.fruit-carbon/2',
    'research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json': 'terluna.research.lunar-cycle-ecology/1',
    'research/studies/open_moon_aerosol/results/open_moon_aerosol.json': 'terluna.research.open-moon-aerosol/1',
    'climate/results/crm/ring_ring_45n_lsw.json': 'terluna.climate.crm-ring/1',
    'climate/results/crm/ring_ring_45s_lsw.json': 'terluna.climate.crm-ring/1',
    'climate/results/crm/ring_ring_70_45e.json': 'terluna.climate.crm-ring/1',
    'climate/results/crm/ring_ring_70_135e.json': 'terluna.climate.crm-ring/1',
    'geography/results/drainage.json': 'terluna.geography.drainage/1',
    'climate/results/gcm/radiation_check_A28_dim5_moon.json': 'terluna.climate.gcm-radiation-check/1',
    'research/studies/summit_tower/results/summit_tower.json': 'terluna.research.summit-tower/1',
    'research/studies/summit_tower/results/summit_tower_form.json': 'terluna.research.summit-tower-form/1',
    'research/studies/inhabited_volume/results/inhabited_volume.json': 'terluna.research.inhabited-volume/1',
}
# The fog desert's weather: CM1's land rings that the aerosol study reads for it.
FOG_RINGS = ['climate/results/crm/ring_ring_45n_lsw.json', 'climate/results/crm/ring_ring_45s_lsw.json',
             'climate/results/crm/ring_ring_70_45e.json', 'climate/results/crm/ring_ring_70_135e.json']
FOG_DESERT_DEG = (zones_mod.RAIN_BELT_DEG, zones_mod.POLAR_DEG)   # the fog desert's latitudes; the tilted
                                                                  # rings' land spans every latitude
WIND_ROUNDING_M_S = 0.1           # the ring products give bin-mean winds to one decimal (climate/crm/ring_analysis.py)
RAYLEIGH_SHAPE = 2.0              # the Weibull shape of the Rayleigh distribution (definition)
# Inputs from the former domain branches, read in-tree since the joint integration, their hashes checked.
BOUND = {
    'ecology_people': dict(
        path='biosphere/ecology/results/people.json', schema='terluna.biosphere.ecology-people/1',
        branch='main', worktree='terluna', commit=None,
        blob='e76556670bf282fa37635f5a4f3976a1f9adf0f1',
        sha256='758a8dcf68c6eb58082e5749e22c35e2da7ffd530f7483622e6c1b3011779ee8',
        note="The lunar design's diets, crop land, crop water, rain-fed bands, nitrogen, phosphorus, N2O and methane "
             'per person (people_and_land.md).'),
    'ecology_first_screen': dict(
        path='biosphere/ecology/results/first_screen.json', schema='terluna.biosphere.ecology-first-screen/1',
        branch='main', worktree='terluna', commit=None,
        blob='08a4e7b179f8176f656b580657a9db7f2423270a',
        sha256='d2a42b6e51f60250f671ca7cf57aa948d1291e2d26d31f84064f038c3829969d',
        note="Methane's soil sink and lifetime and N2O's accumulation per Tg N in the design air."),
    'ecology_people_and_land': dict(
        path='biosphere/ecology/people_and_land.md', schema=None,
        branch='main', worktree='terluna', commit=None,
        blob='39751c8d07e0b74deb42009d31d533d5517c6db3',
        sha256='fbd83a2137a12eb59c2cc4a8ffaed6e37a5aecbedb48c32bd47d15a3c02e38da',
        note="The write-up of the ecology people product; bound for its definitions of the lunar design."),
    'resources_phosphorus': dict(
        path='resources/results/phosphorus.json', schema='terluna.resources.phosphorus/1',
        branch='main', worktree='terluna', commit=None,
        blob='d924ac5ff6ca8c75bc490ee67deccebe9a838c72',
        sha256='adc90d7349332c31bb1c94af3c8061f022a5e8639e676534187c9f6b09e4421a',
        note='The aerial phosphorus return account, whose complete-diet P demand provisioning supplies; its lift '
             'energy per kg to 10 km.'),
    'resources_phosphorus_py': dict(
        path='resources/phosphorus.py', schema=None,
        branch='main', worktree='terluna', commit=None,
        blob='0b8568bc8af9c5ffc09e6baf9287bc86404cdf09',
        sha256='48f3ef951dc6d781f42de706df7652cbd53d945a67fbd8a6e857c253c9e81e43',
        note='The producer of the phosphorus account.'),
    'resources_phosphorus_md': dict(
        path='resources/phosphorus.md', schema=None,
        branch='main', worktree='terluna', commit=None,
        blob='030a5a80e4b6edab39c6a908d5bb519e7b31f721',
        sha256='afa8db4d5e8b698c6e4d8bb4b6d5e20f424f5edabd9144ccab7baf4c0a7af138',
        note='The write-up of the phosphorus account; the author\'s direction of 2026-10-11 to draw phosphorus from '
             'the seas is recorded in the decisions register, Phosphorus from the seas.'),
    'resources_sea_phosphorus': dict(
        path='resources/results/sea_phosphorus.json', schema='terluna.resources.sea-phosphorus/1',
        branch='main', worktree='terluna', commit=None,
        blob='c81e61c5527b67fa1e0950a983ce59aff5dcf4b2',
        sha256='07fea2bba81ab5472853ea9d35a7cfd8155737e52cb1711d560d9d6fc0bbf0a1',
        note='The land\'s phosphorus through the rain-fed lakes to the seas: the dissolved input the common-level '
             'waters receive, and what the waters\' life can spare to seaweed, fish and sorption.'),
}
HOURS_PER_YEAR = JULIAN_YEAR_DAYS * 24.0
CYCLE_HOURS = SYNODIC_MONTH_DAYS * 24.0
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(path) -> str:
    return sha256(path)[:16]


def load_json(rel, schema):
    data = json.loads((ROOT / rel).read_text())
    if data.get('schema') != schema:
        raise ValueError(f'Unexpected schema in {rel}: {data.get("schema")}')
    return data


def bound_bytes(key) -> bytes:
    """A bound input's bytes: the in-tree file, which since the joint integration of 12 October 2026 holds every
    former branch's inputs; its hash is checked against the pin."""
    b = BOUND[key]
    local = ROOT / b['path']
    if not local.exists():
        raise MissingInput(f"{b['path']} is not in this tree")
    raw = local.read_bytes()
    if hashlib.sha256(raw).hexdigest() != b['sha256']:
        raise ValueError(f"{b['path']} does not match its pinned hash")
    return raw


def load_bound(key):
    data = json.loads(bound_bytes(key))
    if data.get('schema') != BOUND[key]['schema']:
        raise ValueError(f'Unexpected schema in bound input {key}')
    return data


def r(x, n=3):
    """Round a float to n significant figures for the product."""
    if x is None or isinstance(x, (bool, int, str)):
        return x
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    return float(f'{x:.{n}g}')


def rr(seq, n=3):
    return [r(x, n) for x in seq]


def rd(obj, n=3):
    """Round every float in a nested structure."""
    if isinstance(obj, dict):
        return {k: rd(v, n) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [rd(v, n) for v in obj]
    return r(obj, n)


def dumps(product):
    """JSON with one key per line and lists of plain values kept on one line."""
    text = json.dumps(product, indent=1)
    return re.sub(r'\[([^\[\]{}]*)\]', lambda m: '[' + re.sub(r'\s+', ' ', m.group(1)).strip() + ']', text)


# ---------------------------------------------------------------------------------------------------------------
# The basis every system shares.

def interp(x, xs, ys):
    """Linear interpolation, held flat beyond the ends."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


def light_basis(data):
    """The crop model's yield and light by latitude, and the sunlight the ground absorbs."""
    fruit = data['biosphere/canopy/results/fruit.json']['results']['near']
    lats = [0.0, 30.0, 60.0]

    def harvest(band, strategy):
        return fruit[band][strategy]['day_fruit']['best']['fed']['shares']['0.6']['harvest_kg']
    ratios = [sum(harvest(b, s) / harvest('0', s) for s in ('earth_like', 'idle_quarter')) / 2.0
              for b in ('0', '30', '60')]
    canopy = data['biosphere/canopy/results/canopy.json']['summaries']['moon']
    par = [canopy[f'{lat:.1f}']['par_mean_umol_m2_s'] for lat in lats]
    cycle_s = CYCLE_HOURS * 3600.0
    plant = data['biosphere/canopy/results/plant.json']['results']['near']
    dark = [plant[b]['earth_like']['dark_hours'] for b in ('0', '30', '60')]
    rad = data['climate/results/gcm/radiation_check_A28_dim5_moon.json']['columns']
    land = sorted((abs(c['lat_deg']), c['arriving_w_m2'] * c['line_by_line']['ground'],
                   c['arriving_w_m2'] * c['gcm']['ground']) for c in rad if c['land'] and c['label'] != 'tropical highland')
    absorbed_lat = [x[0] for x in land]
    absorbed_lo = [min(x[1], x[2]) for x in land]
    absorbed_hi = [max(x[1], x[2]) for x in land]
    return dict(
        fruit_ratio_by_lat=dict(lat_deg=lats, ratio=ratios),
        lat_yield=lambda lat: interp(lat, lats, ratios),
        par_mean_umol_m2_s=dict(lat_deg=lats, par=par),
        par_mol_m2_cycle=lambda lat: interp(lat, lats, par) * 1e-6 * cycle_s,
        dark_hours_by_lat=dict(lat_deg=lats, hours=dark),
        dark_hours=lambda lat: interp(lat, lats, dark),
        absorbed_clear_w_m2=dict(lat_deg=absorbed_lat, low=absorbed_lo, high=absorbed_hi),
        absorbed_clear=lambda lat: [interp(lat, absorbed_lat, absorbed_lo), interp(lat, absorbed_lat, absorbed_hi)])


def rain_basis(eco, zones, gpp_equator):
    """Ecology's rain-fed bands: the rain's share of the crop model's water need, the bands' production against the
    equator's (the plant model's gross production at 0 degrees, as ecology divides by), and the rain-fed land on the
    atlas grid."""
    rows = eco['land']['by_latitude']
    gpp_eq = gpp_equator
    keys = [f"{int(r['abs_lat_deg'][0])}-{int(r['abs_lat_deg'][1])}" for r in rows]
    unscaled = {k: list(r['rain_supports_share_of_model']) for k, r in zip(keys, rows)}
    gpp = {k: r['gross_production_g_c_m2_day'] / gpp_eq for k, r in zip(keys, rows)}
    scaled = {k: [x * gpp[k] for x in unscaled[k]] for k in keys}
    by_key = {k: b for k, b in zip(keys, zones['bands'])}
    cap = [sum(by_key[k]['dry_land_km2'] * scaled[k][i] for k in keys) for i in (0, 1)]
    wet = {k: r['wet_share'] for k, r in zip(keys, rows)}
    return dict(support_unscaled=unscaled, gpp_factor=gpp, support_scaled=scaled, rain_fed_eq_km2=cap,
                rain_fed_eq_km2_ecology=list(eco['land']['potential_equivalent_share_km2']), wet_share=wet), by_key


def water_capped_rain(rain, by_key, need_eco, need, gpp_equator, wue_ecology):
    """Ecology's rain shares at the crop water need the absorbed sunlight allows. The crop model's yields at the
    capped need imply a water-use efficiency raised by ecology's need over the capped one, and the rain's share of
    the model's growth rises in the same proportion. Ecology's shares and rain-fed land are kept for comparison."""
    factor = [need_eco[1] / need[1], need_eco[0] / need[0]]      # pairs run [low share, high share]
    unscaled = {k: [min(1.0, v[0] * factor[0]), min(1.0, v[1] * factor[1])]
                for k, v in rain['support_unscaled'].items()}
    scaled = {k: [x * rain['gpp_factor'][k] for x in unscaled[k]] for k in unscaled}
    cap = [sum(by_key[k]['dry_land_km2'] * scaled[k][i] for k in scaled) for i in (0, 1)]
    return dict(rain, support_unscaled=unscaled, support_scaled=scaled, rain_fed_eq_km2=cap,
                support_unscaled_at_ecology_need=rain['support_unscaled'],
                support_scaled_at_ecology_need=rain['support_scaled'],
                rain_fed_eq_km2_at_ecology_need=rain['rain_fed_eq_km2'], share_factor=factor,
                water_use_efficiency_g_c_per_kg=dict(ecology=list(wue_ecology),
                                                     needed_at_capped_need=gpp_equator / need[1]))


def diet_basis(eco):
    """The lunar design's diets in the ecology screen's terms, and the reference diet's food groups with their
    protein, water and food energy from the food composition table."""
    food = eco['food']
    kcal = [v for k, v in food['kcal_m2_day'].items() if k.split('/')[1] == '0']
    eb = par.value('ecology_basis')
    groups = par.value('eat_lancet_g_day')
    comp = par.value('composition')
    match = par.value('eat_lancet_composition')
    protein = {g: groups[g] * comp[match[g]]['protein'] / 100.0 for g in groups if match.get(g)}
    water_g = sum(groups[g] * comp[match[g]]['water'] / 100.0 for g in groups if match.get(g))
    animal = sum(protein[g] for g in par.value('eat_lancet_animal_groups'))
    produce_km3 = {}
    for diet, factor in (('plant_based', [1.0, 1.0]), ('plant_rich', food['diets']['plant_rich']['crop_factor'])):
        produce_km3[diet] = [water_g * 1e-3 * loss * f * JULIAN_YEAR_DAYS * 1e9 / 1e12
                             for loss, f in zip(eb['loss_factor'], factor)]
    vegfruit_kcal = sum(par.value('eat_lancet_kcal_day')[g] for g in ('vegetables', 'fruits'))
    return dict(
        intake_kcal_day=eb['intake_kcal_day'], loss_factor=eb['loss_factor'], variety_factor=eb['variety_factor'],
        cloud_share=eb['cloud_share'],
        crop_factor={'plant_based': [1.0, 1.0], 'plant_rich': list(food['diets']['plant_rich']['crop_factor'])},
        crop_m2={'plant_based': list(food['diets']['plant_based']['crop_m2']),
                 'plant_rich': list(food['diets']['plant_rich']['crop_m2'])},
        grazing_m2=list(food['diets']['plant_rich']['grazing_m2']),
        crop_kcal_m2_day=[min(kcal), max(kcal)],
        settlement_m2=list(eco['needs']['practices']['lunar_design']['settlement_m2']),
        reference_diet_protein_g_day=sum(protein.values()), animal_source_protein_g_day=animal,
        reference_diet_water_g_day=water_g,
        vegetables_g_day=groups['vegetables'], fruit_g_day=groups['fruits'],
        vegfruit_energy_share=vegfruit_kcal / eb['intake_kcal_day'],
        # the lamps grow the vegetables as lettuce and the fruit as tomato, whose food energy is less than the
        # reference diet's mixed produce
        lamp_produce_kcal_day=(groups['vegetables'] * comp['lettuce']['kcal'] +
                               groups['fruits'] * comp['tomato']['kcal']) / 100.0,
        protein_kcal_per_g=[min(par.value('protein_calorie_factor_kcal_g').values()),
                            max(par.value('protein_calorie_factor_kcal_g').values())],
        produce_water_km3_per_billion=produce_km3)


def fog_rings(data):
    """CM1's fog-desert land rings by hour of the lunar day: fog, the 10 m wind, the Weibull shape of each ring's
    instantaneous land winds over the whole cycle (from their median and 90th percentile), and the night's
    evaporation."""
    rings, night_evap = [], []
    for path in FOG_RINGS:
        ring = data[path]
        land = ring['by_hour_angle']['land']
        hours = ring['by_hour_angle']['hour_angle_deg']
        spread = ring['wind_10m_m_s']['land']
        shape = math.log(math.log(10.0) / math.log(2.0)) / math.log(spread['p90'] / spread['median'])
        rings.append(dict(name=ring['case'], fog=list(land['fog']), wind=list(land['wind_10m_m_s']),
                          shape_whole_cycle=shape, wind_median_p90_m_s=[spread['median'], spread['p90']],
                          in_fog_desert_latitudes=FOG_DESERT_DEG[0] <= abs(ring['latitude_deg']) <= FOG_DESERT_DEG[1]))
        night_evap += [e for h, e in zip(hours, land['evaporation_mm_h']) if abs(h) > 90.0]
    return rings, [min(night_evap), max(night_evap)]


def seas_phosphorus(sea_p):
    """What the waters can give sea farms, from resources' sea-phosphorus study: the dissolved phosphorus reaching the
    common-level waters past the rain-fed lakes, and the share of each water's dissolved input the waters' life can
    spare to seaweed, fish and sorption, in the common-level waters and in all waters, lakes included."""
    routes, easy = sea_p['routes'], sea_p['statements']['easy_route']
    kt = 1e6                                                                            # kg in a kt
    return dict(
        lakes_keep_share=list(easy['kept_by_lakes_share']),
        common_level_dissolved_input_kt_year=[x / kt for x in routes['bioavailable_input_kg_yr']['seas']['base']],
        common_level_ceiling_kt_year=[x / kt for x in routes['seaweed']['ceiling_kg_yr']['seas']['base']],
        all_waters_ceiling_kt_year=[x / kt for x in routes['waters_life_ceiling_kg_yr']['base']],
        all_waters_ceiling_with_particle_release_kt_year=[
            x / kt for x in routes['waters_life_ceiling_kg_yr']['with_particle_release']],
        ceiling_share_of_dissolved_input=routes['ceiling_share'],
        reading=('Sea farms draw on what the waters\' life can spare: 30% of each water\'s dissolved input, a design '
                 'guess of resources\' study for the largest harvest that leaves the food webs whole.'))


def build_basis(data, eco, screen, zones, derived, sea_p):
    """Everything the systems share, from the products, the bound inputs and the literature values."""
    pop = data['provisioning/results/population.json']
    tower = data['research/studies/summit_tower/results/summit_tower.json']
    form = data['research/studies/summit_tower/results/summit_tower_form.json']
    plant = data['biosphere/canopy/results/plant.json']['results']['near']['0']
    lce = data['research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json']
    light = light_basis(data)
    gpp_equator = data['biosphere/canopy/results/plant.json']['results']['near']['0']['earth_like']['gpp'] / \
        SYNODIC_MONTH_DAYS
    rain, bands_by_key = rain_basis(eco, zones, gpp_equator)
    zones = dict(zones, bands_by_key=bands_by_key)
    lunar = eco['needs']['practices']['lunar_design']
    eb = par.value('ecology_basis')
    latent = par.value('latent_heat_mj_kg') * 1e6
    absorbed_eq = light['absorbed_clear'](0.0)
    cap = [x * eb['cloud_share'] / latent * JULIAN_DAY for x in absorbed_eq]          # mm a day
    need_eco = list(eco['land']['crop_water_need_mm_day_equator'])
    need = [need_eco[0], min(need_eco[1], max(cap))]
    rain = water_capped_rain(rain, bands_by_key, need_eco, need, gpp_equator, eco['land']['wue_g_c_per_kg'])
    # the rain-fed land the mixes' fill can reach: the wet bands, the 30-45 degree band held to its wet lowland
    # (the wet middle latitudes) and its own highlands
    middle_wet = zones['zones']['middle_wet']['dry_land_km2']

    def fillable_dry(k):
        b = bands_by_key[k]
        return b['dry_land_km2'] if k != '30-45' else min(b['dry_land_km2'], middle_wet + b['highland_dry_land_km2'])
    for key, shares in (('rain_fed_eq_km2_fillable', rain['support_scaled']),
                        ('rain_fed_eq_km2_fillable_at_ecology_need', rain['support_scaled_at_ecology_need'])):
        rain[key] = [sum(fillable_dry(k) * shares[k][i] for k in mixes_mod.WET_BANDS) for i in (0, 1)]

    def mm_day(w_m2):
        return w_m2 / latent * JULIAN_DAY
    # The absorbed sunlight caps what a watered field transpires, and so what it grows at the crops' efficiency:
    # [high-cost end, low-cost end] in crop-model equivalent per km2, never above the band's production.
    rain['evaporation_cap_mm_day_by_band'] = {
        k: [mm_day(x * eb['cloud_share']) for x in light['absorbed_clear'](0.5 * sum(b['abs_lat_deg']))]
        for k, b in bands_by_key.items()}
    rain['watered_yield_cap'] = {
        k: [min(rain['gpp_factor'][k], c[0] / need[1]), min(rain['gpp_factor'][k], c[1] / need[0])]
        for k, c in rain['evaporation_cap_mm_day_by_band'].items()}
    gt = par.value('greenhouse_transmission')
    absorbed_inside = {}
    for zone, info in zones['zones'].items():
        a = light['absorbed_clear'](info['abs_lat_deg_mean'])
        absorbed_inside[zone] = [a[0] * eb['cloud_share'] * gt[0], a[1] * eb['cloud_share'] * gt[1]]
    light['absorbed_inside_w_m2'] = absorbed_inside
    light['evaporation_cap_inside_mm_day'] = {z: [mm_day(x) for x in v] for z, v in absorbed_inside.items()}
    vm = par.value('venlo_materials')
    lamp_p = par.value('led_efficacy_umol_j')
    nasa = par.value('nasa_edible_g_per_mol')
    comp = par.value('composition')
    sets = par.value('lamp_set_points')

    def kcal_dry(food):
        c = comp[food]
        return c['kcal'] / (100.0 - c['water'])
    rings, night_evap = fog_rings(data)
    pump = data['research/studies/inhabited_volume/results/inhabited_volume.json']['assumptions']['pump_efficiency']
    aq = par.value('earth_aqueducts')
    acre_foot_m3 = 1233.48                       # m3 in an acre-foot (definition: 43,560 cubic feet)
    basis = dict(
        diet=diet_basis(eco),
        nitrogen=dict(new_kg=list(lunar['n_new_kg']), n2o_kg_n=list(lunar['n2o_kg_n']),
                      fixation_w=list(lunar['fixation_w']), eaten_kg=lunar['n_eaten_kg'],
                      applied_kg=[lunar['n_eaten_kg'] / eb['chain_nue'][1], lunar['n_eaten_kg'] / eb['chain_nue'][0]],
                      chain_nue=eb['chain_nue'], field_n2o_share=eb['field_n2o_share']),
        phosphorus=dict(eaten_kg=list(lunar['p_eaten_kg']), lost_kg=list(lunar['p_lost_kg'])),
        methane=dict(kg_person=list(lunar['ch4_kg'])),
        gases=dict(n2o_tg_n_per_ppb=screen['trace_gases']['nitrous_oxide']['tg_n_per_ppb'],
                   ch4_soil_sink_tg_yr_per_ppm=[min(screen['trace_gases']['methane']['soil_sink_tg_yr_per_ppm']),
                                                max(screen['trace_gases']['methane']['soil_sink_tg_yr_per_ppm'])],
                   ch4_lifetime_years=[min(screen['trace_gases']['methane']['lifetime_years']),
                                       max(screen['trace_gases']['methane']['lifetime_years'])],
                   natural_ch4_tg=eco['capacity']['natural_ch4_tg'],
                   natural_n2o_tg_n=list(eco['capacity']['natural_n2o_tg_n']),
                   land_p_loss_tg=list(eco['capacity']['land_p_loss_tg'])),
        heat=dict(moon_w_m2_per_tw=pop['heat']['moon_w_m2_per_tw'],
                  dimming_percent_per_tw=list(pop['heat']['dimming_percent_per_tw']),
                  heat_per_w=dict(array=pop['earth_service']['heat_on_earth_per_w_delivered'],
                                  fusion=[1.0 / eb['fusion_conversion'][1], 1.0 / eb['fusion_conversion'][0]]),
                  residents_heat_tw_per_billion=[pop['affluence']['band_primary_kw']['efficient'] *
                                                 pop['heat']['placement_factor']['industry_in_orbit_thermal'],
                                                 pop['affluence']['band_primary_kw']['high_income'] *
                                                 pop['heat']['placement_factor']['industry_in_orbit_thermal']]),
        supply_land_km2_per_tw=dict(
            array_rectenna=pop['earth_service']['rectenna_km2_per_tw'],
            ground_solar_with_night_storage=1e12 / tower['sunlight']['horizontal_panel_w_m2'] / 1e6 *
            (0.5 + 0.5 / tower['night_storage']['round_trip'])),
        water=dict(crop_need_mm_day=need, crop_need_ecology_mm_day=need_eco, evaporation_cap_mm_day=cap,
                   pump_efficiency=[pump, max(par.value('design_pipe')['pump_efficiency'])],
                   application_efficiency=par.value('application_efficiency')['range'],
                   irrigation_pressure_kpa=par.value('irrigation_pressure_kpa'),
                   river_km3_year=data['geography/results/drainage.json']['rivers']['discharge_to_sea_m3s'] *
                   SECONDS_PER_YEAR / 1e9,
                   air_extraction_l_per_kwh=par.value('air_extraction_l_per_kwh'),
                   earth_aqueducts_kwh_m3=dict(
                       state_water_project=aq['state_water_project_kwh_af'] / acre_foot_m3,
                       colorado_tijuana=aq['colorado_tijuana_kwh'] / aq['colorado_tijuana_m3'],
                       colorado_tijuana_lift_m=aq['colorado_tijuana_lift_m'],
                       colorado_tijuana_km=aq['colorado_tijuana_km'])),
        rain=rain, zones=zones, light=light,
        glass=dict(glass_kg_m2=vm['glass_kg'] / vm['area_m2'], steel_kg_m2=vm['steel_kg'] / vm['area_m2'],
                   aluminium_kg_m2=vm['aluminium_kg'] / vm['area_m2'],
                   concrete_m3_m2=vm['concrete_m3'] / vm['area_m2'], life_years=vm['life_years'],
                   transmission=gt, light_response=par.value('light_response_fruit_vegetables'),
                   co2_gain_vented=derived['co2_gain']['vented'], co2_gain_closed=derived['co2_gain']['closed'],
                   co2_water_use_factor_vented=derived['co2_water_use_factor']['vented'],
                   co2_water_use_factor_closed=derived['co2_water_use_factor']['closed'],
                   operations_kwh_m2_yr=par.value('glasshouse_electricity_kwh_m2_yr'),
                   cooling_cop=par.value('heat_pump_cooling_cop'), roof_capture=par.value('roof_capture')),
        soilless=dict(n2o_share_of_applied=par.value('soilless_n2o_share'),
                      p_discharge_share=par.value('soilless_p_discharge_share')),
        lamps=dict(efficacy_umol_j=[lamp_p['best_tested_2020'][0],
                                    lamp_p['potential'][1] * lamp_p['humidity_protection_factor']],
                   lighting_share=par.value('lighting_share_of_load'),
                   floor_kg_m2=form['basis']['masses_kg_m2']['floor'],
                   staples=dict(potato=dict(edible_g_per_mol=[nasa['potato'][0], nasa['potato_best_transplanted']],
                                            kcal_per_g_dry=kcal_dry('potato'), ppfd=sets['potato']['ppfd'],
                                            photoperiod_h=sets['potato']['photoperiod_h']),
                                wheat=dict(edible_g_per_mol=[nasa['wheat'][0], derived['asseng_dry_grain_g_per_mol']],
                                           kcal_per_g_dry=kcal_dry('wheat'), ppfd=sets['wheat']['ppfd'],
                                           photoperiod_h=sets['wheat']['photoperiod_h'])),
                   layers=par.value('vertical_farm_layers'),
                   leafy_fresh_g_per_mol=derived['lamp_produce_fresh_g_per_mol']['leafy'],
                   fruit_fresh_g_per_mol=derived['lamp_produce_fresh_g_per_mol']['fruit'],
                   night_ppfd=par.value('night_lamp_ppfd')),
        stores=dict(earth_like=plant['earth_like']['store'], idle_quarter=plant['idle_quarter']['store'],
                    night_upkeep_w_m2=lce['night']['night_upkeep_w_m2'], dark_hours=light['dark_hours_by_lat']),
        night=dict(hours=pop['night']['night_hours'], european_final_kwh=pop['night']['final_kwh_per_person']['european']),
        fog=dict(rings=rings, night_land_evaporation_mm_h=night_evap, lwc_g_m3=par.value('fog_lwc_g_m3')['range_used'],
                 wind_rounding_m_s=WIND_ROUNDING_M_S, rayleigh_shape=RAYLEIGH_SHAPE,
                 net_efficiency=[par.value('fog_net_efficiency')['overall_large_collector'],
                                 par.value('fog_net_efficiency')['upper_bound']],
                 stops_below_m_s=par.value('fog_collection_stops_m_s'),
                 earth_operational_l_m2_day=par.value('fog_operational_l_m2_day')),
        dew=dict(earth_l_m2_night=par.value('dew_l_m2_night')),
        microbial=derived['microbial'], fermentation=derived['fermentation'], electro_sugar=derived['electro_sugar'],
        aquaponics=derived['fish'], sea=derived['sea'],
        kept_shares=sorted(eco['capacity']['kept_shares'].values()),
        land_rules=sys_mod.land_rules(eco['capacity']['kept_shares'].values()), farm_yield=par.value('yield_gap'),
        array_intercepted_pw=list(pop['array']['intercepted_pw']),
        seas_phosphorus=seas_phosphorus(sea_p))
    basis['fog_net_l_m2_year'] = \
        sys_mod.fog_desert_water(basis)['fog_nets']['yields_low_to_high']['l_per_m2_year_wind_above_cutoff']
    return basis


# ---------------------------------------------------------------------------------------------------------------
# Derived literature values the systems read.

def co2_gain(ppm_from, ppm_to):
    """Production gain from raising CO2 from ppm_from to ppm_to (ppm at Earth pressure) by Nederhoff's rule,
    integrated: d ln P / dC = coefficient / 100 / 100 x (reference / C)^2."""
    rule = par.value('co2_rule')
    k = rule['coefficient'] / 1e4 * rule['reference_ppm'] ** 2
    return math.exp(k * (1.0 / ppm_from - 1.0 / ppm_to))


def derived_parameters(data, eco):
    """Values computed from the literature entries: the CO2 gain the glass can hold at the design air's CO2, the
    microbial protein's electricity across three accounts, the lamps' produce yields, fish and seas."""
    canopy = data['biosphere/canopy/results/canopy.json']['settings']['co2_pa']
    # the canopy model holds both worlds at 400 ppm (biosphere/canopy/model.py), so the design air's partial
    # pressure is that share of Earth's at Earth's pressure
    ambient = 400.0 * canopy['moon'] / canopy['earth']
    enrich = par.value('co2_enrichment_ppm_earth_equivalent')
    gains = {k: sorted([co2_gain(ambient, max(v[0], ambient)), co2_gain(ambient, max(v[1], ambient))])
             for k, v in enrich.items()}
    # the crops' water per unit growth falls as the CO2 raises production (the gain) and closes the stomata
    # (Nederhoff's 3-4% per 100 ppm): [high-cost end, low-cost end]
    stomata = par.value('co2_stomatal_reduction_per_100_ppm')
    water_use = {k: [gains[k][0] / (1.0 - stomata[0] * (max(v[1], ambient) - ambient) / 100.0),
                     gains[k][1] / (1.0 - stomata[1] * (max(v[0], ambient) - ambient) / 100.0)]
                 for k, v in enrich.items()}
    # microbial protein: electricity per kg protein, three accounts
    aw = par.value('atomic_weights')
    nh3 = aw['N'] + 3 * aw['H']
    n_kwh = par.value('ammonia_electricity_kwh_per_kg_n')
    j = par.value('microbial_jarvio')
    j_kwh = sum(j['electricity_kwh_per_kg_product']) / j['protein_share']
    j_nh3 = j['ammonia_kg_per_kg_product'] * aw['N'] / nh3 * n_kwh / j['protein_share']
    g = par.value('microbial_garcia_martinez')
    g_kwh = [g['kwh_per_kg_dry'] / s for s in g['protein_share'][::-1]]
    leg = par.value('microbial_leger')
    leg_mj = [leg['biomass_mj_kg'] / (leg['electrolysis_efficiency'][1] * leg['biomass_efficiency'][1]) +
              leg['dac_mj_kg'][0] + leg['nutrients_mj_kg'][0] + leg['bioreactor_mj_kg'][0] +
              leg['downstream_food_mj_kg'][0],
              leg['biomass_mj_kg'] / (leg['electrolysis_efficiency'][0] * leg['biomass_efficiency'][0]) +
              leg['dac_mj_kg'][1] + leg['nutrients_mj_kg'][1] + leg['bioreactor_mj_kg'][1] +
              leg['downstream_food_mj_kg'][1]]
    leg_kwh = [leg_mj[0] / 3.6 / leg['protein_share'][1], leg_mj[1] / 3.6 / leg['protein_share'][0]]
    accounts = dict(jarvio=j_kwh + j_nh3, garcia_martinez=[x + 1.0 / 6.25 * n_kwh for x in g_kwh], leger=leg_kwh)
    lo = min(accounts['jarvio'], accounts['garcia_martinez'][0], accounts['leger'][0])
    hi = max(accounts['jarvio'], accounts['garcia_martinez'][1], accounts['leger'][1])
    p_leger = [x * aw['P'] / (3 * aw['H'] + aw['P'] + 4 * aw['O']) for x in leg['h3po4_kg']]
    kh2po4 = aw['K'] + 2 * aw['H'] + aw['P'] + 4 * aw['O']
    na2hpo4 = 2 * aw['Na'] + aw['H'] + aw['P'] + 4 * aw['O'] + 2 * (2 * aw['H'] + aw['O'])
    p_jarvio = (j['kh2po4_kg'] * aw['P'] / kh2po4 + j['na2hpo4_2h2o_kg'] * aw['P'] / na2hpo4) * j['p_use_share'] / \
        (1.0 - j['moisture'])
    na = par.value('single_cell_protein_nucleic_acid')
    limit = par.value('nucleic_acid_limit_g_day')
    microbial = dict(kwh_per_kg_protein=[lo, hi], accounts=accounts, ammonia_kwh_per_kg_protein=[0.0, 0.0],
                     protein_share=[min(g['protein_share'][0], leg['protein_share'][0]),
                                    max(g['protein_share'][1], leg['protein_share'][1])],
                     p_share_of_dw=[min(p_leger), max(max(p_leger), p_jarvio)],
                     co2_kg_per_kg_dw=leg['co2_kg'], protein_per_n=6.25,
                     protein_g_day_at_nucleic_acid_limit=dict(
                         raw=[limit / na['raw'][1] * g['protein_share'][0], limit / na['raw'][0] * g['protein_share'][1]],
                         processed=[limit / na['processed'] * g['protein_share'][0],
                                    limit / na['processed'] * g['protein_share'][1]]))
    pf = par.value('precision_fermentation')
    es = par.value('electro_sugar')
    fermentation = dict(kwh_per_kg_protein=sum(pf['electricity_kwh_per_kg_product']) / pf['protein_share'],
                        heat_kwh_per_kg_protein=pf['heat_kwh_per_kg_product'] / pf['protein_share'],
                        glucose_kg_per_kg_protein=pf['glucose_kg_per_kg_product'] / pf['protein_share'],
                        glucose_kcal_per_g=par.value('glucose_combustion_kj_per_mol') /
                        (6 * aw['C'] + 12 * aw['H'] + 6 * aw['O']) / (sys_mod.KCAL_J / 1e3))
    comp = par.value('composition')
    fish = dict(fillet_protein_share=comp['tilapia']['protein'] / 100.0,
                fillet_yield=sorted([par.value('fish_edible_and_fcr')['tilapia_edible'][1],
                                     par.value('ras_fish')['fillet_yield']['tilapia']]),
                fcr=[par.value('ras_fish')['fcr'], par.value('aquaponics_nutrients')['tilapia_fcr_trial'][1]],
                feed_kcal_per_g=comp['wheat']['kcal'] / 100.0,
                feed_protein_share=par.value('cage_feed')['protein_share'],
                feed_p_share=par.value('cage_feed')['p_share'], protein_per_n=6.25,
                kwh_per_kg_fish=par.value('ras_fish')['electricity_kwh_per_t_live']['tilapia'] / 1e3,
                water_l_per_kg_fish=par.value('ras_fish')['water_m3_per_t_live']['tilapia'],
                n2o_share_of_n=par.value('aquaponics_nutrients')['n2o_share_of_n'],
                p_share_sludge=sorted(1.0 - x for x in par.value('aquaponics_nutrients')['p_kept_in_fish_and_plants']),
                feed_per_kg_crops=par.value('aquaponics_units')['love_crops']['feed_kg'],
                feed_g_m2_bed_day=par.value('aquaponics_units')['delaide_feed_g_m2_day'])
    sw, kelp = par.value('seaweed_china'), par.value('kelp_elite_t_air_dry_ha_yr')
    spro = par.value('seaweed_protein')
    mus, zeb = par.value('mussel_farm'), par.value('zebra_mussel_farm')
    meth = par.value('fresh_water_farm_methane')
    sea = dict(seaweed_t_dw_ha_yr=[sw['t_dw_km2_yr'] / 100.0, kelp[1]],
               seaweed_protein_share=[spro['mean_brown'], sw['n_share'] * spro['n_to_protein']],
               seaweed_p_share=sw['p_share'], seaweed_n_share=sw['n_share'],
               mussel_t_live_ha_yr=[zeb['t_live_ha_yr'], mus['t_live_ha_yr']],
               mussel_meat_share=par.value('mussel_meat_share'),
               mussel_meat_protein_share=comp['mussel']['protein'] / 100.0,
               mussel_p_kg_per_t=[1e3 * x for x in par.value('mussel_whole_np_share')['p']],
               mussel_n_kg_per_t=[1e3 * x for x in par.value('mussel_whole_np_share')['n']],
               farm_sediment_ch4_g_m2_yr=[0.0, (meth['farm_g_c_m2_yr'] - meth['reference_g_c_m2_yr']) *
                                          (aw['C'] + 4 * aw['H']) / aw['C']])
    nasa = par.value('nasa_edible_g_per_mol')
    lettuce_water = comp['lettuce']['water'] / 100.0
    tomato_water = comp['tomato']['water'] / 100.0
    produce = dict(leafy=[nasa['lettuce'][0] / (1.0 - lettuce_water), par.value('lettuce_lamp_fresh_g_per_mol')['carotti_best']],
                   fruit=[nasa['tomato'][0] / (1.0 - tomato_water), par.value('tomato_lamp_and_sun_fresh_g_per_mol')])
    asseng = par.value('wheat_vertical_farm_asseng')
    asseng_mol = asseng['ppfd'] * 1e-6 * 3600.0 * asseng['photoperiod_h'] * asseng['lit_share'] * JULIAN_YEAR_DAYS
    asseng_g_mol = asseng['grain_kg_m2_yr'] * 1e3 * (1.0 - asseng['grain_moisture']) / asseng_mol
    acetate = par.value('acetate_food_hann')
    electro_sugar = dict(power_to_food=es['power_to_food'],
                         acetate_route_power_to_food=dict(
                             yeast=acetate['acetate_efficiency'] * acetate['yeast_efficiency'],
                             algae=acetate['acetate_efficiency'] * acetate['algae_efficiency']))
    return dict(co2_gain=gains, co2_water_use_factor=water_use, ambient_ppm_earth_equivalent=ambient,
                microbial=microbial, fermentation=fermentation,
                electro_sugar=electro_sugar, fish=fish, sea=sea,
                lamp_produce_fresh_g_per_mol=produce, asseng_dry_grain_g_per_mol=asseng_g_mol)


# ---------------------------------------------------------------------------------------------------------------
# The product.

def basis_record(basis):
    """The basis as plain data for the product (the interpolating functions dropped)."""
    out = {}
    for k, v in basis.items():
        if k in ('zones',):
            continue
        if isinstance(v, dict):
            out[k] = {kk: vv for kk, vv in v.items() if not callable(vv)}
        elif k == 'land_rules':
            out[k] = {sys_mod.kept_key(r): dict(closed=float(r), hybrid=r.hybrid) for r in v}
        else:
            out[k] = v
    out['zones'] = {k: v for k, v in basis['zones'].items() if k not in ('bands_by_key',)}
    return out


def ecology_check(basis, eco):
    """The ecology basis reproduces the bound product's crop land and N2O."""
    d = basis['diet']
    kcal = d['crop_kcal_m2_day']
    pb = [d['intake_kcal_day'] / kcal[1] * d['loss_factor'][0] * d['variety_factor'][0] / d['cloud_share'],
          d['intake_kcal_day'] / kcal[0] * d['loss_factor'][1] * d['variety_factor'][1] / d['cloud_share']]
    eb = par.value('ecology_basis')
    n2o = [basis['nitrogen']['applied_kg'][0] * eb['field_n2o_share'][0],
           basis['nitrogen']['applied_kg'][1] * eb['field_n2o_share'][1]]
    return dict(plant_based_crop_m2_recomputed=pb, plant_based_crop_m2_bound=d['crop_m2']['plant_based'],
                n2o_kg_n_recomputed=n2o, n2o_kg_n_bound=basis['nitrogen']['n2o_kg_n'])


def cycle_fit(basis, systems):
    """How each system meets the 354-hour night, with the products' numbers."""
    s, light = basis['stores'], basis['light']
    vf = systems['vertical_farms']['potato']
    mp = systems['microbial_protein_animal_part']
    lamps_glass = systems['glass_night_lamps']
    return {
        'rain_fed_and_watered_fields': dict(
            how=('Evergreen crops idle at a quarter of their daytime upkeep on their own carbon balance; the day fruit '
                 'is set at sunrise and ripe at sunset, so the harvest comes at dusk and towns hold a cycle of food.'),
            store_g_c_m2=dict(earth_like=s['earth_like'], idle_quarter=s['idle_quarter']),
            night_upkeep_w_m2=s['night_upkeep_w_m2'], dark_hours=s['dark_hours'], night_energy_kwh_per_resident=0.0),
        'sunlit_glasshouses': dict(
            how=('As the fields under the glass, on the crops\' own stores; or lamps through the dark hours, which grow '
                 'that share of the food at the vertical farms\' electricity per kcal.'),
            night_lamps_fog_desert=lamps_glass['fog_desert'], night_lamps_highlands=lamps_glass['highlands']),
        'vertical_farms': dict(
            how='Lamps on their own clock; the night\'s energy must come from storage, a day-side grid, fusion or the array.',
            night_energy_kwh_per_resident=[x * basis['night']['hours'] for x in vf['power_kw_per_resident']],
            over_european_night=vf['night_energy_over_european_night']),
        'microbial_protein': dict(
            how='Reactors on firm power, or on day-side power with hydrogen stored through the night.',
            night_energy_kwh_per_resident=list(mp['night_energy_twh'])),
        'aquaponics': dict(how=('Fish fed on a 24-hour routine from stored feed; the vegetables on the fish water idle '
                                'through the night like the glasshouse crops.')),
        'sea_farms': dict(how=('Macroalgae idle through the night on their stores: Arctic kelp keeps its photosynthesis '
                               'through 4 months of polar night and survived 7 months of darkness, against the lunar '
                               'night\'s 15 days; shellfish feed on plankton that idles with them.'),
                          kelp_darkness=par.value('kelp_darkness')),
    }


B_NUM = 1e9


def all_systems(basis):
    """Each system's costs per billion residents."""
    d = basis['diet']
    protein = d['animal_source_protein_g_day']
    out = dict(
        rain_fed_fields=sys_mod.rain_fed_fields(basis, 'plant_rich'),
        rain_fed_fields_plant_based=sys_mod.rain_fed_fields(basis, 'plant_based'),
        irrigated_fields=sys_mod.irrigated_fields(basis, 'plant_rich'),
        fog_desert_glass_vented=sys_mod.glasshouse(basis, 'fog_desert', 'plant_rich', closed=False),
        fog_desert_glass_closed=sys_mod.glasshouse(basis, 'fog_desert', 'plant_rich', closed=True),
        highland_glass_vented=sys_mod.glasshouse(basis, 'highlands', 'plant_rich', closed=False),
        rain_belt_glass_vented=sys_mod.glasshouse(basis, 'rain_belt', 'plant_rich', closed=False),
        glass_night_lamps=dict(fog_desert=sys_mod.night_lamps_for_glass(basis, 'fog_desert'),
                               highlands=sys_mod.night_lamps_for_glass(basis, 'highlands'),
                               rain_belt=sys_mod.night_lamps_for_glass(basis, 'rain_belt')),
        fog_desert_water=sys_mod.fog_desert_water(basis),
        vertical_farms=sys_mod.vertical_farms(basis, 'plant_rich'),
        vertical_farms_vegetables_and_fruit=sys_mod.both(lambda end: dict(
            power_kw_per_resident=sys_mod.lamp_produce(basis, d['vegetables_g_day'], basis['lamps']['leafy_fresh_g_per_mol'], end) +
            sys_mod.lamp_produce(basis, d['fruit_g_day'], basis['lamps']['fruit_fresh_g_per_mol'], end))),
        microbial_protein_animal_part=sys_mod.microbial_protein(basis, protein),
        microbial_protein_half_protein=sys_mod.microbial_protein(basis, d['reference_diet_protein_g_day'] / 2.0),
        precision_fermentation_animal_part=sys_mod.precision_fermentation(basis, protein),
        power_made_food=sys_mod.power_made_food(basis),
        aquaponics_animal_part=sys_mod.aquaponics(basis, protein),
        sea_farms_animal_part=sys_mod.sea_farms(basis, protein),
    )
    return out


def results():
    data = {p: load_json(p, s) for p, s in INPUTS.items()}
    eco = load_bound('ecology_people')
    screen = load_bound('ecology_first_screen')
    phosphorus = load_bound('resources_phosphorus')
    sea_p = load_bound('resources_sea_phosphorus')
    for key in ('ecology_people_and_land', 'resources_phosphorus_py', 'resources_phosphorus_md'):
        bound_bytes(key)
    zones = zones_mod.zone_screen(ROOT)
    derived = derived_parameters(data, eco)
    basis = build_basis(data, eco, screen, zones, derived, sea_p)
    systems = all_systems(basis)
    mixes = mixes_mod.all_mixes(basis)
    record = basis_record(basis)
    phos_lift = next(t for t in phosphorus['transport'] if t['altitude_m'] == 10000 and t['p_mass_fraction'] == 0.1)
    files = [HERE / f for f in CODE]
    product = dict(
        schema=SCHEMA,
        producer=dict(study='food_systems', files={str(f.relative_to(ROOT)): digest(f) for f in files},
                      constants=constants_used(files),
                      inputs={p: digest(ROOT / p) for p in INPUTS},
                      sources=digest(HERE / 'sources.json'),
                      ignored_inputs={p: digest(ROOT / p) for p in zones_mod.GRID_SCHEMAS},
                      bound_inputs={k: dict(path_at_integration=v['path'], branch=v['branch'], commit=v['commit'],
                                            blob=v['blob'], sha256=v['sha256'], worktree=v['worktree'],
                                            read_from='in-tree',
                                            note=v['note']) for k, v in BOUND.items()}),
        evidence=EVIDENCE, reading_rule=READING_RULE, marks=MARKS,
        assumptions=dict(design_guesses={k: dict(value=v['value'], unit=v['unit'], note=v.get('note'))
                                         for k, v in par.P.items() if 'design_guess' in v['source']},
                         mixes={k: v['rule'] for k, v in mixes_mod.MIXES.items()},
                         aloft_share=mixes_mod.ALOFT_SHARE, populations_billion=list(mixes_mod.POPULATIONS),
                         wet_bands_filled=list(mixes_mod.WET_BANDS),
                         hybrid_rules=[dict(r, closed=r['wild'] + r['hybrid'])
                                       for r in sys_mod.HYBRID_RULES]),
        literature={k: dict(value=v['value'], unit=v['unit'], source=v['source']) for k, v in par.P.items()},
        derived=rd(derived),
        basis=rd(record),
        checks=rd(dict(ecology_basis=ecology_check(basis, eco),
                       food_lift_to_10_km=food_lift(phos_lift),
                       crop_model_against_dutch_glass=dutch_check(basis))),
        systems=rd(systems),
        cycle_fit=rd(cycle_fit(basis, systems)),
        mixes=rd(mixes))
    return product


def food_lift(row):
    """Power to lift the food of the residents aloft to 10 km, at resources' energy per kg of material lifted
    (ideal, and at its pump efficiency); the phosphorus rides in that food."""
    aloft = max(mixes_mod.POPULATIONS) * mixes_mod.ALOFT_SHARE                      # billion
    kg_day = sum(par.value('eat_lancet_g_day').values()) / 1e3                       # food eaten a person a day
    ideal_gw = aloft * 1e9 * kg_day / JULIAN_DAY * row['ideal_j_kg'] / 1e9
    return dict(residents_aloft_billion=aloft, food_kg_per_person_day=kg_day, ideal_j_kg=row['ideal_j_kg'],
                efficiency=row['efficiency'], ideal_gw=ideal_gw, gw_at_efficiency=ideal_gw / row['efficiency'])


def dutch_check(basis):
    """The crop model's energy-only fruit against Dutch glasshouse tomatoes in food energy per m2 a year: unlit
    houses, the sunlit house's match, from the norm to the energy-saving crop, and the lit hybrid crop."""
    y = par.value('dutch_yields_kg_m2_yr')
    tomato_kcal_kg = par.value('composition')['tomato']['kcal'] * 10.0
    return dict(crop_model_kcal_m2_yr=[x * JULIAN_YEAR_DAYS for x in basis['diet']['crop_kcal_m2_day']],
                dutch_unlit_tomato_kcal_m2_yr=[y['tomato_unlit_norm'] * tomato_kcal_kg,
                                               y['tomato_unlit_energy_saving'] * tomato_kcal_kg],
                dutch_lit_hybrid_tomato_kcal_m2_yr=y['tomato_lit_hybrid'] * tomato_kcal_kg)


EVIDENCE = ('A screen on committed products (the crop model\'s yields, stores and light, the design climate\'s rain, '
            'fog and wind, the clear-sky radiation check, the atlas and drainage grids, the heat and power ledgers) '
            'and the products of two other domains (ecology\'s diets and nutrients; resources\' phosphorus accounts '
            'for the sky and the waters), with Earth literature values for each food system. Yields, efficiencies and emission '
            'factors are measured on Earth or modelled there and serve as reference points for the Moon\'s engineered '
            'organisms. The study runs no crop, hydrology, climate or ecosystem model and chooses no mix.')
READING_RULE = ('Per billion lunar residents unless named: areas km2, water km3 a year, power and heat TW (TW per '
                'billion equals kW a person), phosphorus and nitrogen Mt a year, N2O kt N a year and ppb, methane Mt a '
                'year and ppm at steady state, structure Gt. Results computed at the two ends are pairs [low-cost end, '
                'high-cost end]: every uncertain input at its favourable value, then at its unfavourable one, so each '
                'end is one consistent case. Inputs, derived values and capacities are ranges, smallest first: the '
                'basis and derived blocks, billions_the_rain_feeds, billions_the_waters_supply, residents_per_km2, the '
                'rain on glass and the fog nets\' yields (yields_low_to_high and each ring\'s). Heat is given for two '
                'supplies: power from the array (its receivers\' losses included) and fusion on the Moon. Land rules: '
                'kept_30 and kept_50 keep 30% and half of each zone wild; wild_30_hybrid_20 is the author\'s '
                'compromise, 30% wild and another 20% hybrid, and wild_50_hybrid_20 half wild and another 20% hybrid '
                'for comparison. The hybrid land is heavily manicured megaparks and great gardens, so fields, glass, '
                'towns and grazing keep outside a rule\'s wild and hybrid shares (the kept share), and the wild shares '
                'count against the land outside the kept share, which differs by rule, so dry_land_share_used '
                'compares the rules directly. Mixes: '
                'cases[population][land rule]; zone_use_share is each zone\'s use against its land outside the kept '
                'share (grazing on the fog desert included), sea farms against the seas shallower than 200 m, and '
                'their phosphorus against what the waters\' life can spare; feasible says which ends fit, '
                'feasible_to the fraction of the way to the high-cost end, every input moved together, at which the '
                'first zone fills, and feasible_range the case there beside the low-cost end.')
MARKS = dict(screen=('Numbers computed here: the zones, the fog-net and water screen, every system\'s costs per '
                     'billion, the farm-yield sensitivity, every mix and its feasibility.'),
             derived=('Arithmetic on cited values: the derived block (CO2 gain and water-use factor in glass, '
                      'microbial protein electricity across three accounts, the acetate route\'s electricity to food, '
                      'glucose\'s food energy, fish and sea values, produce yields under lamps) and the food-lift '
                      'check.'))


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()) + '\n')
    print(OUT.relative_to(ROOT))


if __name__ == '__main__':
    main()
