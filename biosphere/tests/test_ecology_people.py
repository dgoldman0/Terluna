"""The people screen binds its inputs and constants, and its sky, land, needs and capacities keep their arithmetic."""
import hashlib
import json
import math
from pathlib import Path

from biosphere.ecology import people
from shared.constants import MOON_RADIUS, STANDARD_GRAVITY, MOON_SURFACE_GRAVITY, SYNODIC_MONTH_DAYS
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(people.OUT.read_text())
AREA = 4 * math.pi * MOON_RADIUS ** 2


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def close(a, b, rel=0.02):
    return abs(a - b) <= rel * max(abs(a), abs(b))


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == people.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert set(PRODUCT['producer']['inputs']) == {str(p.relative_to(ROOT)) for p, _ in people.INPUTS.values()}
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_inputs_carry_the_schemas_the_screen_expects():
    for path, schema in people.INPUTS.values():
        if schema is not None:
            assert json.loads(path.read_text())['schema'] == schema, path


def test_the_air_bands_add_up_to_the_column():
    sky = PRODUCT['sky']
    bands = list(sky['bands'].values())
    for lower, upper in zip(bands, bands[1:]):
        assert lower['top_km'] == upper['bottom_km']
    mass = sum(b['mass_kg'] for b in bands) + sky['above_90_km_mass_share'] * sky['air_mass_kg']
    assert close(mass, sky['air_mass_kg'], 0.01)
    inventory = json.loads(people.INPUTS['reference_air'][0].read_text())['inventory']['mass_kg']
    assert close(sky['air_mass_kg'], inventory, 0.06)            # the GCM's profile against the design inventory
    assert close(sum(b['mass_share'] for b in bands) + sky['above_90_km_mass_share'], 1.0, 0.01)
    # each square metre of ground carries about 7 Earth columns; pressure over gravity at the ground
    assert close(sky['column_vs_earth'], sky['column_ground_t_m2'] * 1e3 / (people.EARTH_SEA_LEVEL_PA / STANDARD_GRAVITY))


def test_light_aloft_follows_the_geometry_and_the_column():
    sky = PRODUCT['sky']['bands']
    ground = sky['low_sky']['light_at_bottom']
    assert ground['horizon_dip_deg'] == 0.0
    assert close(ground['direct_sun_overhead'], math.exp(-people.RAYLEIGH_TAU_550), 0.005)
    per_hour = 360.0 / (SYNODIC_MONTH_DAYS * 24.0)
    previous = -1.0
    for b in sky.values():
        top = b['light_at_top']
        dip = math.degrees(math.acos(MOON_RADIUS / (MOON_RADIUS + b['top_km'] * 1e3)))
        assert close(top['horizon_dip_deg'], dip, 0.005)
        assert close(top['sun_stays_up_longer_h_each_end_equator'], dip / per_hour, 0.01)
        assert top['direct_sun_overhead'] > previous and top['direct_sun_on_horizon'] < top['direct_sun_overhead']
        previous = top['direct_sun_overhead']


def test_cosmic_ray_ionization_stays_below_earths_sea_level_to_90_km():
    ratios = [b['ionization_mid_vs_earth_sea_level'] for b in PRODUCT['sky']['bands'].values()]
    assert all(r < 1.0 for r in ratios) and ratios == sorted(ratios)
    assert PRODUCT['sky']['ionization_ground_vs_earth_sea_level'] < ratios[0]


def test_lift_scales_by_the_rule_of_similarity():
    lift = PRODUCT['lift']
    k = STANDARD_GRAVITY / MOON_SURFACE_GRAVITY
    assert close(lift['similarity_factor'], k, 0.002)
    for name, (span, mass) in people.EARTH_FLYERS.items():
        row = lift['flyers'][name]
        assert all(close(m, s * k, 0.003) for m, s in zip(row['moon_span_m'], span))
        assert all(close(m * 1e3, x * k ** 3, 0.003) for m, x in zip(row['moon_mass_t'], mass))
    assert 8.0 < lift['earth_sea_level_density_height_km'] < 10.0     # the sky-ship study's 9 km


def test_near_passes_follow_density_speed_and_radius():
    for band in PRODUCT['traffic']['busy_hour'].values():
        expected = band['per_km3'] / 1e9 * band['speed_m_s'] * math.pi * people.NEAR_PASS_M ** 2 * 3600.0
        assert close(band['near_passes_per_hour'], expected, 0.005)


def test_land_shares_are_consistent():
    land = PRODUCT['land']
    assert close(land['moon_km2'], AREA / 1e6, 0.001)
    assert close(land['dry_land'], land['land'] - land['lakes'], 0.001)
    assert close(sum(land['regions'].values()), land['dry_land'], 0.002)
    assert close(sum(b['dry_land_share'] for b in land['by_latitude']), land['dry_land'], 0.005)
    lo, hi = land['potential_equivalent_share']
    assert 0 < lo < hi < land['dry_land']
    for b in land['by_latitude']:
        assert 0.0 < b['rain_supports_share_of_model'][0] <= b['rain_supports_share_of_model'][1] <= 1.0


def test_food_land_follows_the_day_fruit_and_the_factors():
    food = PRODUCT['food']
    equator = [v for k, v in food['kcal_m2_day'].items() if k.split('/')[1] == '0']
    assert close(food['potential_m2_person'][0], people.INTAKE_KCAL_DAY / max(equator), 0.005)
    assert close(food['potential_m2_person'][1], people.INTAKE_KCAL_DAY / min(equator), 0.005)
    lo = food['potential_m2_person'][0] / people.CLOUD_SHARE * people.LOSS_FACTOR[0] * people.VARIETY_FACTOR[0]
    hi = food['potential_m2_person'][1] / people.CLOUD_SHARE * people.LOSS_FACTOR[1] * people.VARIETY_FACTOR[1]
    assert close(food['plant_based_m2_person'][0], lo, 0.01) and close(food['plant_based_m2_person'][1], hi, 0.01)
    for name, diet in people.DIETS.items():
        assert close(food['diets'][name]['crop_m2'][0], lo * diet['crop'][0], 0.01)
        assert close(food['diets'][name]['crop_m2'][1], hi * diet['crop'][1], 0.01)
    # Peters et al. 2016: the US baseline's 0.34 ha of cropland against 0.13 ha for the vegan diet
    assert close(food['diets']['earth_affluent']['crop_factor'][0], 0.34 / 0.13, 0.005)


def test_lit_sky_per_person_follows_falchi_shares_over_land_and_people():
    lit = PRODUCT['needs']['lit_sky_m2_by_region']
    km2, n = people.LAND_AND_PEOPLE_2015['world']
    assert close(lit['world'][0], 0.23 * km2 * 1e6 / n, 0.005)
    design = PRODUCT['needs']['practices']['lunar_design']['lit_sky_m2']
    earth = PRODUCT['needs']['practices']['earth_practice']['lit_sky_m2']
    assert close(design[1], earth[1] * people.LIT_DESIGN_FACTOR, 0.005)


def test_ethylene_settles_where_plant_sources_meet_soil_uptake():
    eth = PRODUCT['ethylene']
    lo, hi = eth['steady_ppb']['earth_soils']
    e_lo, e_hi = eth['steady_ppb']['engineered_soils']
    assert close(lo / e_lo, 10.0, 0.01) and close(hi / e_hi, 10.0, 0.01)       # uptake ten times faster
    assert close(hi / lo, people.ETHYLENE_EARTH_LAND_TG_YR[1] / people.ETHYLENE_EARTH_LAND_TG_YR[0], 0.01)
    assert close(eth['lifetime_years']['earth_soils'], 10 * eth['lifetime_years']['engineered_soils'], 0.01)
    # the steady state is the source over the uptake: burden / lifetime = source
    air_mol = json.loads(people.INPUTS['screen'][0].read_text())['trace_gases']['air_mol']
    burden_tg = air_mol * lo * 1e-9 * 28.05 / 1e12
    assert close(burden_tg / eth['lifetime_years']['earth_soils'], eth['source_tg_yr'][0], 0.02)


def test_capacity_rows_keep_their_arithmetic():
    cap = PRODUCT['capacity']
    screen = json.loads(people.INPUTS['screen'][0].read_text())['trace_gases']
    sink = screen['methane']['soil_sink_tg_yr_per_ppm']
    per_ppb = screen['nitrous_oxide']['tg_n_per_ppb']
    dry_m2 = PRODUCT['land']['dry_land'] * AREA
    fog_m2 = PRODUCT['land']['regions']['fog_desert'] * AREA
    for practice, block in cap['practices'].items():
        need = PRODUCT['needs']['practices'][practice]
        for row in block['rows']:
            n = row['people']
            assert close(row['ch4_tg'][1], n * need['ch4_kg'][1] / 1e9, 0.01)
            assert close(row['ch4_added_steady_ppm'][1], row['ch4_tg'][1] / sink[0], 0.01)
            assert close(row['ch4_added_steady_ppm'][0], row['ch4_tg'][0] / sink[2], 0.01)
            assert close(row['n2o_ppb_per_year'][1], row['n2o_tg_n'][1] / per_ppb, 0.01)
            assert close(row['n2o_ppm_after_years']['1000'][1], row['n2o_ppb_per_year'][1], 0.01)
            assert close(row['heat_w_m2'][1], row['heat_tw'][1] * 1e12 / AREA, 0.01)
            assert close(row['dimming_percent_to_answer'][1], row['warming_k_unanswered'][1] / people.K_PER_PERCENT_DIMMING, 0.01)
            assert close(row['lit_sky_share_of_dry_land'][1], n * need['lit_sky_m2'][1] / dry_m2, 0.01)
        limits = block['limits_people']
        # at each limit the population times its need meets the threshold
        assert close(limits['methane_parity'][0] * need['ch4_kg'][1] / 1e9, cap['natural_ch4_tg'], 0.01)
        assert close(limits['heat_1k_unanswered'][0] * need['heat_kw'][1] * 1e3 / 1e12, cap['tw_per_kelvin'][0], 0.01)
        assert close(limits['lit_sky_half'][0] * need['lit_sky_m2'][1], 0.5 * dry_m2, 0.01)
        if need['grazing_m2'][1]:
            assert close(limits['grazing_half'][0] * need['grazing_m2'][1], 0.5 * fog_m2, 0.01)
        for kept in ('half', '30x30'):
            assert limits[f'land_{kept}'][0] <= limits[f'crops_{kept}'][0] * 1.001
        assert block['order_by_geometric_mean'][0] in limits
    # a terawatt used on the Moon adds about 0.026 W/m2; a kelvin takes 38-61 TW (the joint ledger, provisioning)
    assert close(cap['tw_per_kelvin'][0], people.W_M2_PER_K[0] * AREA / 1e12, 0.01)
    assert 37 < cap['tw_per_kelvin'][0] < 39 and 60 < cap['tw_per_kelvin'][1] < 62


def test_the_lunar_design_needs_less_per_person_at_the_same_nitrogen_eaten():
    earth, lunar = PRODUCT['needs']['practices']['earth_practice'], PRODUCT['needs']['practices']['lunar_design']
    for key in ('food_m2', 'grazing_m2', 'settlement_m2', 'n_new_kg', 'p_lost_kg', 'ch4_kg', 'n2o_kg_n', 'heat_kw',
                'lit_sky_m2'):
        assert lunar[key][1] <= earth[key][1] and lunar[key][0] <= earth[key][0], key
    assert lunar['n_eaten_kg'] == earth['n_eaten_kg'] == people.N_EATEN_KG_YR
    # new nitrogen is what the fields need for the food eaten, less what returns from sewage
    nue, back = people.DESIGN_CHAIN_NUE, people.DESIGN_SEWAGE_N_RETURN
    assert close(lunar['n_new_kg'][1], people.N_EATEN_KG_YR / nue[0] - back[0] * people.N_EATEN_KG_YR, 0.01)
    assert close(lunar['n2o_kg_n'][1], people.N_EATEN_KG_YR / nue[0] * people.N2O_EF[1], 0.01)
    assert earth['n_new_kg'] == list(people.FOOD_N_FOOTPRINT_KG_YR)


def test_horizons_approach_steady_methane_and_grow_n2o_linearly():
    h = PRODUCT['horizons']['by_years']
    shares = [h[k]['methane_share_of_steady_state'][0] for k in ('500', '1000', '5000')]
    assert shares == sorted(shares) and shares[-1] < 1.0
    assert close(h['5000']['natural_n2o_ppm'][1], 10 * h['500']['natural_n2o_ppm'][1], 0.01)
    growth = PRODUCT['horizons']['years_from_0_1_billion']
    assert close(growth['0.01']['10'], math.log(100) / 0.01, 0.005)
