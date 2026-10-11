"""Mixes of food systems for 4, 6 and 8 billion lunar residents.

Each mix is a rule that assigns the diet's parts to systems:
- the plant foods: the plant-based diet's crops, at the crop model's yields (ecology);
- the animal part: what the plant-rich diet adds to them, either as ecology's feed crops and grazing or as the
  same protein from fish, shellfish, macroalgae or bacteria.
Fields fill the rain-fed bands wettest first outside the kept share and the towns; what the rain cannot grow is
grown by watering those fields to what the absorbed sunlight lets them grow, and beyond that by watered fields in the
fog desert. Half the residents live aloft (2, 3 and 4 billion, within the author's 2-4 billion), the rest in towns on
the highlands at the lunar design's 25-100 m2 a person. Each case is evaluated at the low-cost and high-cost ends of
every input and for a kept 30% and a kept half, and checked against its zones: a case that overflows one cannot be
built, and the feasible range runs to the case where the first zone fills.
"""
from __future__ import annotations

from research.studies.food_systems import systems as S
from research.studies.food_systems.systems import B, KM2, both, pick, yearly, zip_ends

POPULATIONS = (4.0, 6.0, 8.0)
ALOFT_SHARE = 0.5
WET_BANDS = ('0-15', '15-30', '30-45')          # bands where most land is wet (ecology), as ways of living fills
BAND_ZONE = {'0-15': 'rain_belt', '15-30': 'rain_belt', '30-45': 'middle_wet'}
FIT_TOLERANCE = 1e-9              # rounding room when a fill exactly reaches a zone's edge

MIXES = {
    'rain_fed_first': dict(
        rule=('The plant-rich diet on fields: rain-fed first, wettest land first; the rest by watering those fields '
              'from the seas, then watered fields in the fog desert; ecology\'s grazing on the fog desert.'),
        plant='fields', animal='ecology'),
    'greenhouse_heavy': dict(
        rule=('Half the plant foods in vented sunlit glasshouses on the highlands beside the towns, half on rain-fed '
              'fields; the animal part as aquaponic fish fed from the fields.'),
        plant='half_glass_highlands', animal='aquaponics'),
    'fog_desert_glass': dict(
        rule=('All plant foods and the fish feed in vented glasshouses in the fog desert, watered by rain on the glass '
              'and pipes from the seas; the animal part as aquaponic fish; the rain belt left to the wild.'),
        plant='glass_fog_desert', animal='aquaponics_under_glass'),
    'power_heavy': dict(
        rule=('Vegetables and fruit in lamp-lit vertical farms, the animal part as microbial protein on power, the rest '
              'of the plant foods on rain-fed fields.'),
        plant='fields_less_vegfruit', animal='microbial'),
    'fully_lamp_lit': dict(
        rule=('Every plant food in lamp-lit vertical farms on a potato-like staple, the animal part as microbial '
              'protein: no farmland.'),
        plant='lamps', animal='microbial'),
    'power_made': dict(
        rule=('No farmland: microbial protein for half the reference diet\'s protein and fungal protein fermented on '
              'sugar made from power for the other half, sugar made by electrochemistry for the food energy, '
              'vegetables and fruit from lamp-lit farms.'),
        plant='power_made', animal='none'),
    'sea_supported': dict(
        rule=('The animal part from the fresh seas, half as macroalgae and half as filter-feeding shellfish; the plant '
              'foods on rain-fed fields.'),
        plant='fields', animal='sea'),
}


def field_fill(basis, eq_need, end, kept, reserved):
    """Fields for eq_need km2 of crop-model cropland: rain-fed in the wet bands, wettest first, on land outside the
    kept share and the reserved land; the shortfall grown by watering the filled fields to what the sunlight lets
    them grow; what still remains on watered fields in the fog desert."""
    rain, bands = basis['rain'], basis['zones']['bands_by_key']
    support = {b: pick(rain['support_scaled'][b], end, low_is_good=False) for b in WET_BANDS}
    order = sorted(WET_BANDS, key=lambda b: -support[b])
    remaining = eq_need
    phys, eq_rain, by_zone = {}, {}, dict(rain_belt=0.0, middle_wet=0.0, highlands=0.0)
    middle_wet = basis['zones']['zones']['middle_wet']['dry_land_km2']
    for b in order:
        # the 30-45 degree band holds its wet lowland, the wet middle latitudes, and its own highlands
        dry_b = (bands[b]['dry_land_km2'] if b != '30-45' else
                 min(bands[b]['dry_land_km2'], middle_wet + bands[b]['highland_dry_land_km2']))
        avail = max(dry_b * (1.0 - kept) - reserved.get(b, 0.0), 0.0)
        take = min(remaining / support[b], avail) if support[b] > 0 else 0.0
        phys[b], eq_rain[b] = take, take * support[b]
        remaining -= eq_rain[b]
        high = max(bands[b]['highland_dry_land_km2'] * (1.0 - kept) - reserved.get(b, 0.0), 0.0)
        hs = high / avail if avail > 0 else 0.0
        by_zone['highlands'] += take * hs
        by_zone[BAND_ZONE[b]] += take * (1.0 - hs)
    shortfall = max(remaining, 0.0)
    grows = {b: pick(rain['watered_yield_cap'][b], end, low_is_good=False) for b in (*WET_BANDS, '45-60')}
    headroom = sum(phys[b] * max(grows[b] - support[b], 0.0) for b in WET_BANDS)
    watered_eq = min(shortfall, headroom)
    fog_eq = shortfall - watered_eq
    w = basis['water']
    need = S.crop_water_mm(basis, end)
    eff = pick(w['application_efficiency'], end, low_is_good=False)
    belt = basis['zones']['zones']['rain_belt']
    fog = basis['zones']['zones']['fog_desert']
    s_fog = pick(rain['support_scaled']['45-60'], end, low_is_good=False) / grows['45-60']
    km3_belt = yearly(need, watered_eq) / eff
    km3_fog = yearly(need, fog_eq) * (1.0 - s_fog) / eff
    kwh_belt = S.pumping_kwh_m3(basis, belt['height_above_sea_km_mean'], belt['km_to_sea_or_great_lake_mean'], end,
                                emitters=True)
    kwh_fog = S.pumping_kwh_m3(basis, fog['height_above_sea_km_mean'], fog['km_to_sea_or_great_lake_mean'], end,
                               emitters=True)
    fog_km2 = fog_eq / grows['45-60']
    # the other way: water the whole need on the wettest band, each field to what it can grow
    grows0 = pick(rain['watered_yield_cap']['0-15'], end, low_is_good=False)
    s0 = pick(rain['support_scaled']['0-15'], end, low_is_good=False) / grows0
    alt_km3 = yearly(need, eq_need) * (1.0 - s0) / eff
    return dict(
        all_watered_alternative=dict(field_km2=eq_need / grows0, withdrawal_km3_year=alt_km3,
                                     lift_tw=S.tw_of(alt_km3, kwh_belt)),
        eq_need_km2=eq_need, rain_fed_eq_km2=sum(eq_rain.values()), field_km2=sum(phys.values()),
        field_km2_by_zone=by_zone, watered_eq_km2=watered_eq, fog_desert_watered_km2=fog_km2,
        rain_used_km3_year=yearly(need, sum(eq_rain.values())),
        withdrawal_km3_year=km3_belt + km3_fog,
        lift_tw=S.tw_of(km3_belt, kwh_belt) + S.tw_of(km3_fog, kwh_fog))


def soilless_n2o_kt(basis, people_billion, end):
    """N2O from soilless culture growing people_billion residents' share of the diet: the lunar design's applied
    nitrogen at the soilless emission factor."""
    return (people_billion * pick(basis['nitrogen']['applied_kg'], end) *
            pick(basis['soilless']['n2o_share_of_applied'], end) * 1e3)


def reserve_by_band(basis, highland_km2):
    """Spread a highland use over the bands in proportion to their highland dry land."""
    bands = basis['zones']['bands_by_key']
    total = sum(b['highland_dry_land_km2'] for b in bands.values())
    return {k: highland_km2 * b['highland_dry_land_km2'] / total for k, b in bands.items()}


def mix_case(basis, mix, n, end, kept):
    """One mix for n billion lunar residents at one end and one kept share."""
    rule = MIXES[mix]
    d = basis['diet']
    aloft = n * ALOFT_SHARE
    surface = n - aloft
    towns = surface * pick(d['settlement_m2'], end) * B / KM2
    plant_eq = n * S.eq_km2(basis, 'plant_based', end)
    animal_protein = d['animal_source_protein_g_day']
    land = dict(rain_belt=0.0, middle_wet=0.0, fog_desert=0.0, highlands=towns, polar=0.0)
    glass_km2 = 0.0
    glass_zone = None
    # new nitrogen for the diet's crops; bacterial protein carries its own ammonia in its electricity accounts
    microbial_g = dict(microbial=animal_protein, power_made=d['reference_diet_protein_g_day'] / 2.0)
    bacterial = microbial_g.get(rule['animal'], microbial_g.get(rule['plant'], 0.0))
    crop_protein_share = 1.0 - bacterial / d['reference_diet_protein_g_day']
    power = dict(irrigation=0.0, glasshouse_operations=0.0, glasshouse_water=0.0, lamps=0.0, microbial=0.0,
                 aquaponics=0.0, power_made_food=0.0,
                 nitrogen_fixation=n * pick(basis['nitrogen']['fixation_w'], end) * B / 1e12 * crop_protein_share)
    water = dict(rain_used=0.0, withdrawn_from_seas=0.0, rain_on_glass=0.0)
    gases = dict(n2o_kt_n=0.0, ch4_mt=n * pick(basis['methane']['kg_person'], end))
    phosphorus = dict(eaten=n * pick(basis['phosphorus']['eaten_kg'], end),
                      lost=n * pick(basis['phosphorus']['lost_kg'], end), drawn_from_seas=0.0, glass_discharge=0.0,
                      in_microbial_biomass=0.0, to_fish_sludge=0.0)
    structure = dict(glass_gt=0.0, steel_gt=0.0, aluminium_gt=0.0, building_gt=0.0)
    night = dict(twh=0.0)
    sea = dict(sea_km2=0.0)
    field_eq = 0.0
    glass_eq = 0.0
    grazing = 0.0
    feed_on_glass = rule['animal'] == 'aquaponics_under_glass'
    feed = dict(eq=0.0, applied_n_mt=0.0, p_mt=0.0)          # fish feed: its crops' area, nitrogen and phosphorus

    # ---- the animal part
    if rule['animal'] == 'ecology':
        field_eq += n * (S.eq_km2(basis, 'plant_rich', end) - S.eq_km2(basis, 'plant_based', end))
        grazing = n * pick(d['grazing_m2'], end) * B / KM2
    elif rule['animal'] in ('aquaponics', 'aquaponics_under_glass'):
        aq = S.aquaponics(basis, animal_protein, end=end)
        feed.update(eq=n * aq['feed_crop_equivalent_km2'], p_mt=n * aq['feed_phosphorus_mt_year'],
                    applied_n_mt=n * aq['feed_nitrogen_mt_year'] /
                    pick(basis['nitrogen']['chain_nue'], end, low_is_good=False))
        if feed_on_glass:
            glass_eq += feed['eq']
        else:
            field_eq += feed['eq']
        power['aquaponics'] += n * aq['power_tw']
        water['withdrawn_from_seas'] += n * aq['water_km3_year']
        phosphorus['to_fish_sludge'] += n * aq['phosphorus_to_sludge_mt_year']
        gases['n2o_kt_n'] += n * aq['n2o']['kt_n_year']
        night['twh'] += n * aq['power_tw'] * basis['night']['hours']
    elif rule['animal'] == 'microbial':
        mp = S.microbial_protein(basis, animal_protein, end=end)
        power['microbial'] += n * mp['power_tw']
        phosphorus['in_microbial_biomass'] += n * mp['phosphorus_mt_year']
        night['twh'] += n * mp['night_energy_twh']
    elif rule['animal'] == 'sea':
        sf = S.sea_farms(basis, animal_protein / 2.0, end=end)
        sea['sea_km2'] = n * (sf['seaweed']['sea_km2'] + sf['shellfish']['sea_km2'])
        phosphorus['drawn_from_seas'] = n * (sf['seaweed']['phosphorus_drawn_mt_year'] +
                                             sf['shellfish']['phosphorus_drawn_mt_year'])
        gases['ch4_mt'] += sea['sea_km2'] * KM2 * pick(basis['sea']['farm_sediment_ch4_g_m2_yr'], end) / 1e12

    # ---- the plant foods
    if rule['plant'] == 'fields':
        field_eq += plant_eq
    elif rule['plant'] == 'half_glass_highlands':
        field_eq += plant_eq / 2.0
        glass_eq += plant_eq / 2.0
        glass_zone = 'highlands'
    elif rule['plant'] == 'glass_fog_desert':
        glass_eq += plant_eq
        glass_zone = 'fog_desert'
    elif rule['plant'] == 'fields_less_vegfruit':
        field_eq += plant_eq * (1.0 - d['vegfruit_energy_share'])
        kw = (S.lamp_produce(basis, d['vegetables_g_day'], basis['lamps']['leafy_fresh_g_per_mol'], end) +
              S.lamp_produce(basis, d['fruit_g_day'], basis['lamps']['fruit_fresh_g_per_mol'], end))
        power['lamps'] += n * kw
        night['twh'] += n * kw * basis['night']['hours']
        gases['n2o_kt_n'] += soilless_n2o_kt(basis, n * d['vegfruit_energy_share'], end)
    elif rule['plant'] == 'power_made':
        pm = S.power_made_food(basis, end=end)
        power['power_made_food'] = n * pm['power_tw']
        night['twh'] += n * pm['night_energy_twh']
        phosphorus['in_microbial_biomass'] = n * pm['phosphorus_in_microbial_biomass_mt_year']
        gases['n2o_kt_n'] += soilless_n2o_kt(basis, n * d['vegfruit_energy_share'], end)
    elif rule['plant'] == 'lamps':
        vf = S.vertical_farms(basis, 'plant_based', end=end)['potato']
        power['lamps'] += n * vf['power_tw']
        structure['building_gt'] += n * vf['building_gt']
        night['twh'] += n * vf['night_energy_twh']
        gases['n2o_kt_n'] += soilless_n2o_kt(basis, n, end)

    # ---- glass
    if glass_eq > 0.0:
        zone = glass_zone or 'fog_desert'
        g = S.glasshouse(basis, zone, 'plant_based', closed=False, end=end)
        scale = glass_eq / (S.eq_km2(basis, 'plant_based', end))      # billions' worth of plant-based crops
        feed_scale = feed['eq'] / S.eq_km2(basis, 'plant_based', end) if feed_on_glass else 0.0
        glass_km2 = scale * g['glass_km2']
        land[zone] += glass_km2
        power['glasshouse_operations'] += scale * g['power_tw']['operations']
        power['glasshouse_water'] += scale * g['power_tw']['water']
        water['withdrawn_from_seas'] += scale * g['from_seas_km3_year']
        water['rain_on_glass'] += scale * min(g['rain_on_glass_km3_year'], g['consumed_km3_year'])
        for key in ('glass_gt', 'steel_gt', 'aluminium_gt'):
            structure[key] += scale * g[key]
        # the plant foods at the diet's nitrogen and phosphorus; the fish feed at its own
        gases['n2o_kt_n'] += (scale - feed_scale) * g['nitrogen']['n2o']['kt_n_year']
        phosphorus['glass_discharge'] += (scale - feed_scale) * g['phosphorus_mt_year']['discharged']
        if feed_on_glass:
            gases['n2o_kt_n'] += feed['applied_n_mt'] * pick(basis['soilless']['n2o_share_of_applied'], end) * 1e3
            phosphorus['glass_discharge'] += feed['p_mt'] * pick(basis['soilless']['p_discharge_share'], end)

    # ---- fields
    fill = None
    if field_eq > 0.0:
        reserved = reserve_by_band(basis, towns + (glass_km2 if glass_zone == 'highlands' else 0.0))
        fill = field_fill(basis, field_eq, end, kept, reserved)
        for zname, km2 in fill['field_km2_by_zone'].items():
            land[zname] += km2
        land['fog_desert'] += fill['fog_desert_watered_km2']
        power['irrigation'] += fill['lift_tw']
        water['rain_used'] += fill['rain_used_km3_year']
        water['withdrawn_from_seas'] += fill['withdrawal_km3_year']
        feed_on_fields = 0.0 if feed_on_glass else feed['eq']
        share = (field_eq - feed_on_fields) / (n * S.eq_km2(basis, 'plant_rich', end))
        gases['n2o_kt_n'] += n * pick(basis['nitrogen']['n2o_kg_n'], end) * 1e3 * share
        if feed_on_fields > 0.0:
            gases['n2o_kt_n'] += feed['applied_n_mt'] * pick(basis['nitrogen']['field_n2o_share'], end) * 1e3

    total_power = sum(power.values())
    seas_p = basis['seas_phosphorus']
    zones = basis['zones']['zones']
    zone_area = {z: zones[z]['dry_land_km2'] for z in land}
    used = sum(land.values())
    dry = basis['zones']['dry_land_km2']
    out = dict(
        people_billion=dict(lunar=n, aloft=aloft, surface=surface),
        land_km2=dict(land, total=used, glass=glass_km2),
        dry_land_share_used=used / dry,
        dry_land_wild_share_outside_kept=1.0 - (used + grazing) / (dry * (1.0 - kept)),
        zone_wild_share_outside_kept={z: 1.0 - (land[z] + (grazing if z == 'fog_desert' else 0.0)) /
                                      (zone_area[z] * (1.0 - kept)) for z in land},
        sea_farms_km2=sea['sea_km2'],
        grazing_km2=grazing, grazing_share_of_fog_desert_outside_kept=grazing / (zone_area['fog_desert'] * (1.0 - kept)),
        glass_share_of_zone=(glass_km2 / zone_area[glass_zone]) if glass_zone else 0.0,
        residents_heat_tw=n * S.pick(basis['heat']['residents_heat_tw_per_billion'], end),
        power_tw=dict(power, total=total_power),
        heat=S.heat(total_power, basis, end),
        water_km3_year=water,
        phosphorus_mt_year=dict(phosphorus, net_new_input=phosphorus['lost'] + phosphorus['glass_discharge'] -
                                phosphorus['drawn_from_seas'],
                                drawn_over_common_level_dissolved_input=phosphorus['drawn_from_seas'] * 1e3 /
                                pick(seas_p['common_level_dissolved_input_kt_year'], end, low_is_good=False),
                                drawn_over_all_waters_ceiling=phosphorus['drawn_from_seas'] * 1e3 /
                                pick(seas_p['all_waters_ceiling_kt_year'], end, low_is_good=False)),
        n2o=S.n2o(gases['n2o_kt_n'], basis),
        methane=S.ch4(gases['ch4_mt'], basis, end),
        structure_gt=structure,
        night_energy_twh=night['twh'],
        night_energy_over_european_night=night['twh'] / (n * basis['night']['european_final_kwh']))
    if fill is not None:
        rf_cap = pick(basis['rain']['rain_fed_eq_km2_fillable'], end, low_is_good=False)
        out['fields'] = dict(fill, rain_fed_land_share_used=fill['rain_fed_eq_km2'] / rf_cap,
                             rain_fed_wild_share_outside_kept=1.0 - fill['rain_fed_eq_km2'] / (rf_cap * (1.0 - kept)))
    # Each zone's use, grazing on the fog desert included, against its land outside the kept share; sea farms
    # against the seas shallower than 200 m, and their phosphorus against what the waters' life can spare in all
    # waters, lakes included. A case that overflows any of them cannot be built.
    use = dict(land)
    use['fog_desert'] += grazing
    shares = {z: use[z] / (zone_area[z] * (1.0 - kept)) for z in use}
    shares['sea_farms_in_seas_shallower_than_200_m'] = sea['sea_km2'] / basis['zones']['sea_km2_by_depth_m']['200']
    shares['sea_farms_phosphorus_in_the_waters'] = out['phosphorus_mt_year']['drawn_over_all_waters_ceiling']
    out['zone_use_share'] = shares
    out['sea_farms_share_of_seas_shallower_than_50_m'] = sea['sea_km2'] / basis['zones']['sea_km2_by_depth_m']['50']
    out['overflowing'] = sorted(z for z, v in shares.items() if v > 1.0 + FIT_TOLERANCE)
    out['feasible'] = not out['overflowing']
    return out


def feasible_limit(basis, mix, n, kept, tol=1e-4):
    """The largest fraction t of the way from the low-cost end to the high-cost end, every input moved together, at
    which the mix still fits its zones: 1 when the high-cost end fits, None when the low-cost end does not."""
    if mix_case(basis, mix, n, 'high', kept)['feasible']:
        return 1.0
    if not mix_case(basis, mix, n, 'low', kept)['feasible']:
        return None
    lo, hi = 0.0, 1.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mix_case(basis, mix, n, mid, kept)['feasible']:
            lo = mid
        else:
            hi = mid
    return lo


def all_mixes(basis):
    """Every mix at every population and kept share: both ends, which of them fit, and the feasible range, which runs
    from the low-cost end to the case where the first zone fills."""
    out = {}
    for mix, rule in MIXES.items():
        cases = {}
        for n in POPULATIONS:
            cases[f'{n:g}'] = {}
            for k in basis['kept_shares']:
                case = both(lambda end, n=n, k=k: mix_case(basis, mix, n, end, k))
                t = feasible_limit(basis, mix, n, k)
                case['feasible_to'] = t
                case['feasible_range'] = (zip_ends(mix_case(basis, mix, n, 'low', k), mix_case(basis, mix, n, t, k))
                                          if t is not None and t < 1.0 else None)
                cases[f'{n:g}'][f'kept_{int(k * 100)}'] = case
        out[mix] = dict(rule=rule['rule'], cases=cases)
    return out
