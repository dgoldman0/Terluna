"""Joint ledgers of the three domains at the author's choices (the joint integration of 12 October 2026).

    python -m research.studies.domain_synthesis.ledgers
    # -> results/ledgers.json

The resources, ecology and provisioning branches merged on 12 October 2026. Each study was built on the others'
products as they then stood; these ledgers recombine them for what the author has since chosen: distribution C'
(8 billion on the Moon, 4 on the ground and 4 aloft, the land 30% wild and 20% hybrid) fed by the food blend
(power-heavy 40%, power-made 50% and fully lamp-lit 10% of the diet) with heritage meat from the homesteads. The
distributions study placed C' with the power-heavy mix; here the blend's food takes its place in the totals, and
placing the blend's fields band by band, as the distributions study does, is left to that study.

- heat: C''s heat by component, with fusion and with power from the array, the blend's food in place of the
  power-heavy mix's, close computing beside it, the W/m2 and the added dimming that answers it, and the ring
  fleet's infrared glow, its dimming taken as provisioning takes it (the glow over the design sunlight);
- land: the dry land under the rule: C''s settlements, the blend's fields, the heritage meat's pasture and the wild
  share left, beside C' with the power-heavy mix;
- gases: methane and N2O from the wild waters and soils, the blend's farms and the heritage cattle, against the
  soils' methane sink and N2O's lack of one on any human timescale: the steady states, the removal that holds
  each gas at a reference level, the share of the atmosphere that removal from open air would treat each year at
  full capture, and the least power that complete separation takes;
- phosphorus: resources' ledger at 8 billion, with the blend's losses at both ends, and the sky's need beside it;
- heights: the sky's layers and what uses them, from the products, with the storms' tops, the cloud at the
  network's height and the sky towns' passes of sky reefs, by the distributions study's encounter rule;
- array power: with power from the array, what the sky towns' tops and the ground's receivers take, at the 1978
  reference rectenna's intensities, and how wide a beam from the ring lands.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path

from shared.constants import GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, SPEED_OF_LIGHT

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'ledgers.json'
SCHEMA = 'terluna.research.domain-synthesis-ledgers/1'
INPUTS = {
    'research/studies/distributions/results/distributions.json': 'terluna.research.distributions/1',
    'research/studies/food_systems/results/food_blend.json': 'terluna.research.food-blend/1',
    'provisioning/results/population.json': None,
    'biosphere/ecology/results/first_screen.json': None,
    'biosphere/ecology/results/people.json': None,
    'resources/results/sea_phosphorus.json': None,
    'research/studies/sky_ecology/results/sky_ecology.json': 'terluna.research.sky-ecology/1',
    'research/studies/inhabited_volume/results/inhabited_volume.json': None,
    'research/studies/sky_ships/results/sky_ships.json': None,
    'research/studies/summit_tower/results/summit_tower_form.json': None,
    'research/studies/sky_fleet/results/sky_fleet.json': None,
    'research/studies/solar_shield_array/results/ring_layout.json': None,
    'research/studies/solar_shield_array/results/integrated_ledger.json': None,
}
RULE = 'wild_30_hybrid_20'              # the author's compromise: 30% of the land wild and another 20% hybrid
KEPT, HYBRID = 0.3, 0.2
MOON_BILLION = 8
ALOFT_BILLION = 4
HERITAGE_KG = ('0.5', '1', '2')         # heritage meat a person a year, the blend's three cases
# Reference levels for holding the gases: Earth's present air, as the ecology screen and brief give it (the screen's
# methane value; machines.md's 0.34 ppm of N2O).
N2O_REFERENCE_PPM = 0.34
# Complete denitrification's cut to N2O, from one N2O-respiring bacterium on Earth's farm soils (Hiis et al. 2024,
# Nature 630, 421-428: 50-95% depending on the soil); applied here to every source, the wild waters included, as a
# design bracket.
DENITRIFICATION_CUT = (0.5, 0.95)
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY
# Sky reefs: a kilometre or two across (sky ecology README; aerophytes README), half the aerophyte cover (sky
# ecology's design guess, in its product). Towns meet them by the distributions study's encounter rule.
REEF_WIDTH_KM = (1.0, 2.0)
PASS_MARGIN_KM = 0.05                    # the distributions study's margin on each side of a pass
TOWN_DISTRICT_T = '50.0'                 # the distributions study's town at 50 t of district a resident
# Beamed power: the two microwave bands of the 1978 reference system and the shield study's links; that system's
# rectenna peak, 23 mW/cm2 at its centre (US DOE & NASA 1978), beside the mean provisioning reads; and ICNIRP's (2020)
# reference level for the general public at 2-300 GHz, whole body.
BEAM_GHZ = (2.45, 5.8)
RECTENNA_PEAK_W_M2 = 230.0
PUBLIC_LIMIT_W_M2 = 10.0
EVIDENCE = ("Arithmetic on committed screen products of the three domains and main's studies: C' and the food blend "
            "recombined in totals, the trace gases' sources and sinks summed, resources' phosphorus ledger read at the "
            "blend's losses, and the sky's layers read from their products. Nothing is simulated here; each figure "
            "keeps its source's evidence state, and the blend's fields are not placed band by band.")
READING_RULE = ("Pairs are [low-cost end, high-cost end] as in the food and distributions studies unless labelled; "
                "areas km2, power and heat TW, W/m2 over the whole Moon, dimming in per cent of sunlight, methane Tg "
                "and ppm, N2O Tg N and ppb or ppm, phosphorus Tg P (Mt) a year; heights km above sea level except the "
                "local air traffic, above the metropolis's ground, and the ring fleet, by radius. Gas cases "
                "pair the sources' low end with the soils' fastest sink and the high end with the slowest; central "
                "takes the screen's central sources and sink with the farms at both ends and the cattle at 1 kg a "
                "person (the larger of the two carcass yields' figures). The share of the atmosphere "
                "treated is the air a full-capture removal from open air would pass each year, over the whole air.")


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, schema=None):
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected schema {data.get("schema")}')
    return data


def rnd(x, sig=3):
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def close(a, b, rel=5e-3):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-9)


def pair(f):
    return [f(0), f(1)]


def heat_ledger(cprime, blend, power_heavy, population):
    """C''s heat with the blend's food in place of the power-heavy mix's, for each supply."""
    h = population['heat']
    dim, w_m2 = h['dimming_percent_per_tw'], h['moon_w_m2_per_tw']
    out = {}
    for supply in ('fusion', 'array'):
        energy = cprime['energy'][supply]
        parts = {k: list(v) for k, v in energy['heat_by_component_tw'].items()}
        placed = power_heavy['heat'][supply]['heat_tw']
        if not all(close(parts['food'][e], placed[e]) for e in (0, 1)):
            raise ValueError(f"C''s food heat ({supply}) is not the power-heavy mix's")
        parts['food'] = list(blend['heat'][supply]['heat_tw'])
        total = pair(lambda e: sum(v[e] for v in parts.values()))
        computing = pair(lambda e: energy['close_computing_dimming_percent'][e] / dim[e])
        with_computing = pair(lambda e: total[e] + computing[e])
        out[supply] = dict(
            heat_by_component_tw=parts,
            heat_tw=total,
            heat_w_m2=pair(lambda e: total[e] * w_m2),
            added_dimming_percent=pair(lambda e: total[e] * dim[e]),
            food_share=pair(lambda e: parts['food'][e] / total[e]),
            close_computing_heat_tw=computing,
            with_close_computing=dict(heat_tw=with_computing, heat_w_m2=pair(lambda e: with_computing[e] * w_m2),
                                      added_dimming_percent=pair(lambda e: with_computing[e] * dim[e])),
            as_placed_with_power_heavy=dict(
                heat_tw=list(energy['heat_tw']), added_dimming_percent=list(energy['added_dimming_percent']),
                with_close_computing=dict(heat_tw=pair(lambda e: energy['heat_tw'][e] + computing[e]),
                                          added_dimming_percent=pair(lambda e: (energy['heat_tw'][e] + computing[e])
                                                                     * dim[e]))))
    glow = h['fleet_glow_w_m2']
    out['fleet_glow'] = dict(w_m2=list(glow), as_tw_used=list(h['fleet_glow_as_tw_used']),
                             added_dimming_percent=[g / h['design_sunlight_w_m2'] * 100.0 for g in glow],
                             basis='the glow over the design sunlight, as provisioning gives it')
    for supply in ('fusion', 'array'):
        f = out[supply]['with_close_computing']
        out[f'{supply}_with_glow'] = dict(w_m2=pair(lambda e: f['heat_w_m2'][e] + glow[e]),
                                          added_dimming_percent=pair(lambda e: f['added_dimming_percent'][e]
                                                                     + out['fleet_glow']['added_dimming_percent'][e]))
    out['factors'] = dict(dimming_percent_per_tw=list(dim), moon_w_m2_per_tw=w_m2)
    return out


def land_ledger(cprime, blend, heritage, people):
    """The dry land's uses under the rule, with the blend's fields and the heritage pasture."""
    dry = people['capacity']['dry_land_km2']
    kept, parks = KEPT * dry, HYBRID * dry
    if not close(parks, cprime['land']['parks_km2'][0]):
        raise ValueError("C''s parks are not the rule's hybrid fifth")
    outside = dry - kept - parks
    settled = list(cprime['land']['settled_km2'])
    towns = blend['towns_km2']
    zones = blend['land_km2']
    fields = list(blend['fields_km2'])
    if not all(close(fields[e], zones['total'][e] - towns[e]) for e in (0, 1)):
        raise ValueError("The blend's fields are not its land less its towns")
    fields_by_zone = {z: v for z, v in zones.items() if z not in ('total', 'glass') and max(v) > 0}
    fields_by_zone['highlands'] = pair(lambda e: zones['highlands'][e] - towns[e])
    cases = {}
    for kg in HERITAGE_KG:
        pasture = heritage['cases'][f'{MOON_BILLION}_billion_{kg}_kg']['land_km2']
        used = pair(lambda e: settled[e] + fields[e] + pasture)
        cases[f'{kg}_kg'] = dict(heritage_pasture_km2=pasture, used_km2=used,
                                 used_share_of_dry_land=pair(lambda e: used[e] / dry),
                                 wild_share_outside_kept=pair(lambda e: (outside - used[e]) / outside),
                                 wild_share_of_dry_land=pair(lambda e: (kept + outside - used[e]) / dry))
    return dict(dry_land_km2=dry, kept_wild_km2=kept, parks_km2=parks, outside_kept_and_parks_km2=outside,
                settled_km2=settled, settled_by_setting_km2=cprime['land']['settled_by_setting_km2'],
                blend_fields_km2=fields, blend_fields_by_zone_km2=fields_by_zone,
                middle_wet_fields_km2=list(blend['middle_wet_fields_km2']),
                homestead_region_km2=heritage['homestead_region_km2'],
                with_heritage=cases,
                as_placed_with_power_heavy=dict(food_km2=list(cprime['land']['food_km2']),
                                                used_km2=list(cprime['land']['used_km2']),
                                                wild_share_outside_kept=sorted(cprime['land']['wild_share_outside_kept'],
                                                                               reverse=True)))


def removal(rate_tg, tg_per_unit, level_units, air_mol, temperature_k, per_unit):
    """For a removal rate in Tg (or Tg N) a year: the share of the atmosphere a full-capture removal from open air
    passes each year at the level held, and the least power that complete separation takes, RT(ln(1/x) +
    (1-x)/x ln(1/(1-x))) a mole of the gas, near RT(ln(1/x) + 1) at these levels."""
    units_per_year = rate_tg / tg_per_unit                   # ppm or ppb a year
    share = units_per_year / level_units
    mol = units_per_year * per_unit * air_mol
    x = level_units * per_unit
    per_mol = math.log(1.0 / x) - (1.0 - x) / x * math.log1p(-x)
    watts = mol * GAS_CONSTANT * temperature_k * per_mol / SECONDS_PER_YEAR
    return dict(atmosphere_share_per_year=share, least_power_mw=watts / 1e6)


def gas_ledger(screen, people, blend, heritage, temperature_k):
    """Methane and N2O: sources, sinks, steady states and the removal that holds each at a reference level."""
    tg = screen['trace_gases']
    air_mol = tg['air_mol']
    m = tg['methane']
    sink = m['soil_sink_tg_yr_per_ppm']                    # the soils' uptake, slowest to fastest
    fast, central_sink, slow = max(sink), sorted(sink)[1], min(sink)
    tg_per_ppm, ref = m['tg_per_ppm'], m['earth_ppm']
    farms = blend['methane']['mt_year']
    cattle = {kg: heritage['cases'][f'{MOON_BILLION}_billion_{kg}_kg']['cattle_ch4_mt_year'] for kg in HERITAGE_KG}
    wild_c = m['central']['sources_tg_yr']
    if not close(wild_c, people['capacity']['natural_ch4_tg']):
        raise ValueError('The screen and the people product disagree on the wild methane')

    def methane_case(wild, farm, meat, k):
        s = wild + farm + meat
        hold = max(0.0, s - k * ref)
        return dict(sources_tg_yr=dict(wild=wild, farms=farm, heritage_cattle=meat, total=s),
                    people_share=(farm + meat) / s, soil_sink_tg_yr_per_ppm=k, steady_ppm=s / k,
                    hold_at_reference_tg_yr=hold, people_part_tg_yr=farm + meat,
                    open_air_hold=removal(hold, tg_per_ppm, ref, air_mol, temperature_k, 1e-6))
    methane = dict(
        reference_ppm=ref, tg_per_ppm=tg_per_ppm,
        wild_sources_tg_yr=dict(lakes=m['sources_tg_yr']['lakes'], wetlands=m['sources_tg_yr']['wetlands'],
                                total=m['sources_tg_yr']['total'], central=wild_c),
        low=methane_case(m['sources_tg_yr']['total'][0], farms[0], min(cattle['0.5']), fast),
        central=[methane_case(wild_c, farms[e], max(cattle['1']), central_sink) for e in (0, 1)],
        high=methane_case(m['sources_tg_yr']['total'][1], farms[1], max(cattle['2']), slow))
    n = tg['nitrous_oxide']
    wild_n = n['sources_tg_n_yr']
    farm_n = [x / 1e3 for x in blend['n2o']['kt_n_year']]
    per_ppb = n['tg_n_per_ppb']

    def n2o_case(wild, farm):
        s = wild + farm
        rise = s / per_ppb
        left = [s * (1.0 - c) for c in DENITRIFICATION_CUT[::-1]]
        return dict(sources_tg_n_yr=dict(wild=wild, farms=farm, total=s), people_share=farm / s,
                    rise_ppb_per_year=rise, ppm_after_1000_years=rise, years_to_reference=N2O_REFERENCE_PPM * 1e3 / rise,
                    hold_tg_n_yr=s, left_after_complete_denitrification_tg_n_yr=left,
                    open_air_hold=removal(s, per_ppb, N2O_REFERENCE_PPM * 1e3, air_mol, temperature_k, 1e-9))
    nitrous = dict(reference_ppm=N2O_REFERENCE_PPM, tg_n_per_ppb=per_ppb,
                   lifetime_years_from_gap_light=n['lifetime_years_from_gap_light'],
                   low=n2o_case(wild_n[0], farm_n[0]), high=n2o_case(wild_n[1], farm_n[1]),
                   denitrification_cut=list(DENITRIFICATION_CUT),
                   heritage_manure='not counted: the blend prices the cattle\'s enteric methane only')
    return dict(methane=methane, n2o=nitrous, ground_temperature_k=temperature_k)


def phosphorus_ledger(sea, blend, cprime, sky):
    """Resources' ledger at 8 billion, with the blend's people's losses at both ends."""
    lost = blend['phosphorus_mt_year']['lost']
    if not all(close(lost[e] * 1e3, cprime['phosphorus']['need_kt'][e]) for e in (0, 1)):
        raise ValueError("The blend's phosphorus losses are not C''s need")
    corners = sea['ledger']['design_target']['corners']
    out = {}
    for land in ('low_loss', 'high_loss'):
        c = corners[f'{land}_{MOON_BILLION}_billion']
        if not close(c['people_losses'] / 1e9, lost[1]):
            raise ValueError("Resources' 8-billion corner is not the blend's high-cost loss")
        kept_by_lakes = c['to_lakes'] / c['into_waters']
        rows = []
        for e in (0, 1):
            people_tg = lost[e]
            into = c['land_loss'] / 1e9 + people_tg
            to_lakes = into * kept_by_lakes
            rows.append(dict(people_losses=people_tg, into_waters=into, to_lakes=to_lakes,
                             dredged_to_farms=people_tg, dredged_share_of_lake_input=people_tg / to_lakes))
        out[land] = dict(land_loss=c['land_loss'] / 1e9, carriers_to_wild_land=c['carriers_to_wild_land'] / 1e9,
                         kept_by_lakes_share=kept_by_lakes, by_food_end=rows,
                         resources_corner=dict(to_lakes=c['to_lakes'] / 1e9, lake_burial=c['lake_burial'] / 1e9,
                                               sea_burial=c['sea_burial'] / 1e9,
                                               dredged_share_of_lake_input=c['dredged_share_of_lake_input']))
    p = blend['phosphorus_mt_year']
    need = sky['headlines']['sky_need_t_year_at_1pct']
    return dict(by_land_loss=out, eaten=list(p['eaten']), lost=list(lost),
                in_microbial_biomass=list(p['in_microbial_biomass']),
                sewage_recoverable_power_heavy=[x / 1e3 for x in cprime['phosphorus']['sewage_recoverable_kt']],
                sky_need_at_one_percent_cover={k: v / 1e6 for k, v in need.items()},
                reef_p_subsidy_g_m2_year=sky['headlines']['reef_p_subsidy_g_m2_year'])


def reef_passes(cover, town_km, drift, share, widths=REEF_WIDTH_KM):
    """Bodies a town passes a day at an aerophyte cover by the distributions study's rule, density x swept width x
    relative speed, when bodies of the given widths make up the share of the cover: [fewest, most]."""
    rates = []
    for d in widths:
        density = cover * share / (math.pi * (d / 2.0) ** 2)
        swept = town_km + d + 2.0 * PASS_MARGIN_KM
        rates += [density * swept * v * JULIAN_DAY / 1e3 for v in drift]
    return [min(rates), max(rates)]


def heights_ledger(people, volume, ships, tower, dist, sky, fleet):
    """The sky's layers and what uses them, with the sky towns' passes of sky reefs."""
    bands = {k: dict(bottom_km=v['bottom_km'], top_km=v['top_km'], mass_share=v['mass_share'],
                     cloud_fraction_ring=v['cloud_fraction_ring']) for k, v in people['sky']['bands'].items()}
    ground_km = tower['basis']['ground_above_sea_m'] / 1e3
    local = {k: dict(above_ground_km=[v['bottom_m'] / 1e3, v['top_m'] / 1e3],
                     above_sea_level_at_the_summit_km=[ground_km + v['bottom_m'] / 1e3, ground_km + v['top_m'] / 1e3],
                     aloft_per_million_people=v['aloft_per_million_people'])
             for k, v in people['traffic']['busy_hour'].items()}
    storms = ships['storms']
    crown_km = (tower['basis']['ground_above_sea_m'] + tower['width_profile']['height_m']) / 1e3
    traffic = dist['moon']['cprime'][RULE]['traffic'][TOWN_DISTRICT_T]
    town_km = 2.0 * dist['basis']['sky']['envelope_radius_km']
    drift = dist['assumptions']['design_guesses']['relative_drift_m_s']['value']
    small = traffic['aerophyte_passes_per_town_day']
    for cover, published in small.items():          # the rule reproduces the study's 100 m bodies
        mine = reef_passes(float(cover), town_km, drift, 1.0, widths=(0.1,))
        if not all(close(a, b, 1e-2) for a, b in zip(mine, published)):
            raise ValueError('The encounter rule does not reproduce the distributions study')
    share = sky['phosphorus']['design_guesses']['class_shares_of_cover']['sky_reef']
    return dict(bands=bands, summit_ground_km=ground_km, local_air_traffic_above_the_metropolis=local,
                sky_towns=dict(reference_height_km=volume['assumptions']['reference_height_km'],
                               district_t_per_resident=float(TOWN_DISTRICT_T), towns=traffic['towns'],
                               people_per_town=ALOFT_BILLION * 1e9 / traffic['towns'],
                               spacing_even_km=traffic['spacing_even_km'],
                               spacing_busiest_band_km=traffic['spacing_busiest_band_km'],
                               encounters_per_town_day_mean=traffic['encounters_per_town_day_mean'],
                               passes_per_town_day_if_100_m_bodies_made_the_cover=small,
                               reef_passes_per_town_day={c: reef_passes(float(c), town_km, drift, share) for c in small},
                               reef_share_of_cover=share, reef_width_km=list(REEF_WIDTH_KM),
                               relative_drift_m_s=list(drift)),
                storms=dict(cloud_top_km=storms['cloud_top_km'], width_km=storms['width_km'],
                            per_snapshot_on_ring=storms['per_snapshot'], drift_m_s=storms['drift_m_s']),
                flight_band_storm_cloud_share=ships['band']['storm_cloud_share'],
                cloud_at_30_km=dict(day=bands['cool_air']['cloud_fraction_ring']['day'][0],
                                    night=bands['cool_air']['cloud_fraction_ring']['night'][0])
                if bands['cool_air']['cloud_fraction_ring']['heights_km'] == [30.0] else None,
                regional_docks_km=fleet['air']['regional_docks']['height_km'],
                summit_crown_km=crown_km)


def array_power(heat, heights, dist, population, ring, ledger, blend):
    """Power from the array, divided between the sky towns' tops and the ground's receivers. Everything received
    ends as heat on the Moon, so the intake is the heat ledger's array parts; the aloft share of the residents goes
    to the towns. Beams follow the 1978 system's capture parameter tau (At Ar = (tau lambda d)^2)."""
    parts = heat['array']['heat_by_component_tw']
    aloft = ALOFT_BILLION / MOON_BILLION
    sky_keys = ('sky_control', 'sky_hydrogen_remade', 'sky_trips_and_food_lift')
    sky = pair(lambda e: sum(parts[k][e] for k in sky_keys) + aloft * parts['residents'][e])
    ground = pair(lambda e: parts['food'][e] + parts['water'][e] + (1.0 - aloft) * parts['residents'][e]
                  + heat['array']['close_computing_heat_tw'][e])
    towns = heights['sky_towns']['towns']
    radius_m = dist['basis']['sky']['envelope_radius_km'] * 1e3
    tops_m2 = towns * math.pi * radius_m ** 2
    service = population['earth_service']
    tau, mean = service['tau_1978'], service['rectenna_mean_w_m2']
    r = ring['design']['middle_radius_km'] * 1e3
    distances = [r - MOON_RADIUS, math.sqrt(r ** 2 - MOON_RADIUS ** 2)]          # overhead and at the horizon
    tile_m2 = ledger['requirements']['tile_side_m'] ** 2
    beams = {}
    for ghz in BEAM_GHZ:
        lam = SPEED_OF_LIGHT / (ghz * 1e9)
        receiver = [(tau * lam * d) ** 2 / tile_m2 for d in distances]
        transmitter = [(tau * lam * d) ** 2 / (math.pi * radius_m ** 2) for d in distances]
        beams[f'{ghz:g}_ghz'] = dict(receiver_m_from_a_tile=[2.0 * math.sqrt(a / math.pi) for a in receiver],
                                     transmitter_side_m_for_a_town_top=[math.sqrt(a) for a in transmitter])
    fields = blend['fields_km2']
    return dict(
        sky_towns=dict(intake_tw=sky, tops_km2=tops_m2 / 1e6, w_m2=pair(lambda e: sky[e] * 1e12 / tops_m2)),
        ground=dict(intake_tw=ground,
                    receiver_km2_at_1978_mean=pair(lambda e: ground[e] * 1e12 / mean / 1e6),
                    receiver_km2_at_1978_peak=pair(lambda e: ground[e] * 1e12 / RECTENNA_PEAK_W_M2 / 1e6),
                    w_m2_over_the_blends_fields=pair(lambda e: ground[e] * 1e12 / (fields[e] * 1e6))),
        intensities_w_m2=dict(rectenna_1978_mean=mean, rectenna_1978_peak=RECTENNA_PEAK_W_M2,
                              public_reference_level=PUBLIC_LIMIT_W_M2),
        beams=beams, ring_middle_radius_km=r / 1e3, distance_km=[d / 1e3 for d in distances], tau=tau)


def results():
    data = {p: load(p, s) for p, s in INPUTS.items()}
    dist = data['research/studies/distributions/results/distributions.json']
    food = data['research/studies/food_systems/results/food_blend.json']
    if food['criterion']['planned_land_rule'] != RULE or food['criterion']['planned_population_billion'] != MOON_BILLION:
        raise ValueError("The blend's planned case is not C''s")
    cprime = dist['moon']['cprime'][RULE]['full']
    blend = food['blend'][str(MOON_BILLION)][RULE]
    heritage = food['heritage_meat']
    power_heavy = food['reference_planned_case']['power_heavy_alone']
    people = data['biosphere/ecology/results/people.json']
    volume = data['research/studies/inhabited_volume/results/inhabited_volume.json']
    ground = next(r for r in volume['air'] if r['height_km'] == 0.0)
    heat = heat_ledger(cprime, blend, power_heavy, data['provisioning/results/population.json'])
    heights = heights_ledger(people, volume, data['research/studies/sky_ships/results/sky_ships.json'],
                             data['research/studies/summit_tower/results/summit_tower_form.json'], dist,
                             data['research/studies/sky_ecology/results/sky_ecology.json'],
                             data['research/studies/sky_fleet/results/sky_fleet.json'])
    return dict(
        schema=SCHEMA,
        producer=dict(files={'research/studies/domain_synthesis/ledgers.py': sha256(__file__)},
                      inputs={p: sha256(ROOT / p) for p in INPUTS}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        choices=dict(distribution="C'", moon_billion=MOON_BILLION, aloft_billion=ALOFT_BILLION, land_rule=RULE,
                     blend_weights=food['weights'], heritage_kg_per_person_year=[float(k) for k in HERITAGE_KG]),
        heat=rnd(heat),
        land=rnd(land_ledger(cprime, blend, heritage, people)),
        gases=rnd(gas_ledger(data['biosphere/ecology/results/first_screen.json'], people, blend, heritage,
                             ground['temperature_c'] + 273.15)),
        phosphorus=rnd(phosphorus_ledger(data['resources/results/sea_phosphorus.json'], blend, cprime,
                                         data['research/studies/sky_ecology/results/sky_ecology.json'])),
        heights=rnd(heights),
        array_power=rnd(array_power(heat, heights, dist, data['provisioning/results/population.json'],
                                    data['research/studies/solar_shield_array/results/ring_layout.json'],
                                    data['research/studies/solar_shield_array/results/integrated_ledger.json'],
                                    blend)))


def dumps(product) -> str:
    return json.dumps(product, indent=2, allow_nan=False) + '\n'


def main(argv=None) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
