"""Each food system's costs per billion lunar residents, its yields, and how it meets the 354-hour night.

Every system is evaluated twice: at its low-cost end, with every uncertain input at its favourable value, and at its
high-cost end, with every input at its unfavourable value. A range [low, high] in the product is those two ends; each
end is one consistent set of inputs. The basis (diet, heat, zones and the products' numbers) is assembled in run.py;
literature values come from parameters.py.
"""
from __future__ import annotations

import math

from scipy.special import gamma, gammaincc

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, SYNODIC_MONTH_DAYS

from research.studies.food_systems.parameters import P, value

B = 1e9                                   # residents in a billion
KM2 = 1e6                                 # m2 in a km2
SECONDS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY
CYCLE_HOURS = SYNODIC_MONTH_DAYS * 24.0
CYCLES_PER_YEAR = JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS
KWH_J = 3.6e6                             # J in a kWh (definition)
KCAL_J = 4184.0                           # J in a thermochemical kilocalorie (definition)
WATER_KG_M3 = 1000.0                      # fresh water, as ways of living takes it
ENDS = ('low', 'high')


def pick(pair, end, low_is_good=True):
    """The favourable value of an uncertain input at the low-cost end, the unfavourable one at the high-cost end. A
    number t from 0 to 1 gives the case with every input moved that fraction of the way from its favourable value to
    its unfavourable one (the mixes' feasibility scan)."""
    lo, hi = (pair if isinstance(pair, (list, tuple)) else (pair, pair))
    good, bad = (lo, hi) if low_is_good else (hi, lo)
    if end == 'low':
        return good
    if end == 'high':
        return bad
    return good + float(end) * (bad - good)


def both(fn):
    """Run fn at both ends and zip its (nested dict) outputs into [low, high] pairs."""
    return zip_ends(fn('low'), fn('high'))


def zip_ends(a, b):
    if isinstance(a, dict):
        return {k: zip_ends(a[k], b[k]) for k in a}
    if isinstance(a, str) or a is None:
        return a
    return [a, b]


def yearly(mm_day, area_km2):
    """km3 a year of mm_day over area_km2."""
    return mm_day * 1e-3 * area_km2 * KM2 * JULIAN_YEAR_DAYS / 1e9


def tw_of(km3_year, kwh_m3):
    """Mean power (TW) to move km3_year at kwh_m3."""
    return km3_year * 1e9 * kwh_m3 * KWH_J / SECONDS_PER_YEAR / 1e12


def lift_kwh_m3(height_km, pump_efficiency):
    """Electricity to lift a cubic metre of fresh water height_km at the Moon's surface gravity."""
    return WATER_KG_M3 * MOON_SURFACE_GRAVITY * height_km * 1e3 / pump_efficiency / KWH_J


def pipe_kwh_m3(distance_km, end):
    """Electricity to push a cubic metre of water distance_km along a large pipe (Darcy-Weisbach); the pressure loss
    does not depend on gravity."""
    pipe = P['design_pipe']['value']
    dp_per_m = pick(pipe['friction_factor'], end) / pipe['diameter_m'] * WATER_KG_M3 * pipe['speed_m_s'] ** 2 / 2.0
    return dp_per_m * distance_km * 1e3 / pick(pipe['pump_efficiency'], end, low_is_good=False) / KWH_J


def pumping_kwh_m3(basis, height_km, distance_km, end, emitters=False):
    """Electricity to lift a cubic metre from the seas to height_km and pipe it distance_km, and for irrigation to
    press it through the emitters, with one pump efficiency at each end."""
    w = basis['water']
    eff = pick(w['pump_efficiency'], end, low_is_good=False)
    kwh = lift_kwh_m3(height_km, eff) + pipe_kwh_m3(distance_km, end)
    if emitters:
        kwh += pick(w['irrigation_pressure_kpa'], end) * 1e3 / eff / KWH_J
    return kwh


def heat(power_tw, basis, end):
    """Heat on the Moon from power used on it, for each supply, with the W/m2 and the added dimming that answers
    it. Power from the array releases its receivers' losses as well; fusion its plants' waste heat."""
    h = basis['heat']
    out = {}
    for supply, factor in h['heat_per_w'].items():
        q = power_tw * pick(factor, end)
        out[supply] = dict(heat_tw=q, w_m2=q * h['moon_w_m2_per_tw'],
                           added_dimming_percent=q * pick(h['dimming_percent_per_tw'], end))
    return out


def n2o(kt_n_year, basis):
    """N2O, which the design air cannot destroy: what a year adds and what a thousand years add."""
    per_ppb = basis['gases']['n2o_tg_n_per_ppb'] * 1e3
    return dict(kt_n_year=kt_n_year, ppb_per_year=kt_n_year / per_ppb, ppb_after_1000_years=1e3 * kt_n_year / per_ppb)


def ch4(mt_year, basis, end):
    """Methane, which soils alone take up: its steady addition at the soil sink's slow or fast end."""
    sink = pick(basis['gases']['ch4_soil_sink_tg_yr_per_ppm'], end, low_is_good=False)
    return dict(mt_year=mt_year, steady_ppm=mt_year / sink)


def diet_kcal(basis, diet, end):
    """Food energy grown per resident per day on a diet, in the ecology screen's terms: intake, the loss factor, the
    variety factor against an energy-only crop and the diet's crop factor."""
    d = basis['diet']
    return (d['intake_kcal_day'] * pick(d['loss_factor'], end) * pick(d['variety_factor'], end) *
            pick(d['crop_factor'][diet], end))


def eq_km2(basis, diet, end):
    """Cropland at the crop model's equatorial yields for a billion residents."""
    return pick(basis['diet']['crop_m2'][diet], end) * B / KM2


def crop_water_mm(basis, end):
    return pick(basis['water']['crop_need_mm_day'], end)


# ---------------------------------------------------------------------------------------------------------------
# Field crops.

def field_nutrients(basis, end, people_billion=1.0):
    """The lunar design's nitrogen, phosphorus, N2O and methane per person (ecology), for people_billion."""
    n, p, m = basis['nitrogen'], basis['phosphorus'], basis['methane']
    eaten = pick(p['eaten_kg'], end) * people_billion
    return dict(
        phosphorus_mt_year=dict(eaten=eaten, in_harvest=eaten * pick(basis['diet']['loss_factor'], end),
                                lost=pick(p['lost_kg'], end) * people_billion),
        nitrogen=dict(new_mt_n_year=pick(n['new_kg'], end) * people_billion,
                      applied_mt_n_year=pick(n['applied_kg'], end) * people_billion,
                      n2o=n2o(pick(n['n2o_kg_n'], end) * people_billion * 1e3, basis)),
        methane=ch4(pick(m['kg_person'], end) * people_billion, basis, end))


def rain_fed_fields(basis, diet='plant_rich'):
    """The lunar design's field crops on the rain alone (ecology's crop model yields and rain-fed bands), against the
    rain-fed land the mixes can fill: the bands where most land is wet, 0-45 degrees."""
    rain = basis['rain']

    def one(end):
        eq = eq_km2(basis, diet, end)
        best = pick(rain['support_scaled']['0-15'], end, low_is_good=False)
        cap = pick(rain['rain_fed_eq_km2_fillable'], end, low_is_good=False)
        power = pick(basis['nitrogen']['fixation_w'], end) * B / 1e12
        out = dict(crop_model_equivalent_km2=eq, land_km2_on_the_wettest_band=eq / best,
                   share_of_rain_fed_land=eq / cap,
                   billions_the_rain_feeds={f'kept_{int(k * 100)}': cap * (1 - k) / eq for k in basis['kept_shares']},
                   rain_water_used_km3_year=yearly(crop_water_mm(basis, end), eq),
                   power_tw=dict(nitrogen_fixation=power), heat=heat(power, basis, end))
        if diet == 'plant_rich':
            graze = pick(basis['diet']['grazing_m2'], end) * B / KM2
            out.update(grazing_km2=graze,
                       grazing_share_of_fog_desert=graze / basis['zones']['zones']['fog_desert']['dry_land_km2'])
        out.update(field_nutrients(basis, end))
        out['billions_the_rain_feeds'] = {k: v for k, v in out['billions_the_rain_feeds'].items()}
        return out
    res = both(one)
    # the rain feeds most people at the low-cost end; keep [fewest, most] order
    res['billions_the_rain_feeds'] = {k: sorted(v) for k, v in res['billions_the_rain_feeds'].items()}
    res['yield'] = dict(kcal_m2_day_energy_only_crop=basis['diet']['crop_kcal_m2_day'],
                        residents_per_km2=[KM2 / x for x in basis['diet']['crop_m2'][diet][::-1]])
    res['farm_yield_sensitivity'] = farm_yield_sensitivity(basis, diet)
    res['night'] = night_fields(basis)
    return res


def farm_yield_sensitivity(basis, diet='plant_rich'):
    """The crop model's yields are potential yields. Fields whose farms reach a share of that potential need land in
    proportion, and the rain then feeds fewer people."""
    rain, fy = basis['rain'], basis['farm_yield']
    out = {}
    for name, share in (('practical_ceiling', fy['irrigated_share_of_potential']),
                        ('earth_rain_fed_common', fy['rain_fed_share_of_potential'])):
        def one(end, share=share):
            s = pick(share, end, low_is_good=False)
            eq = eq_km2(basis, diet, end) / s
            cap = pick(rain['rain_fed_eq_km2_fillable'], end, low_is_good=False)
            return dict(share_of_potential=s, land_factor=1.0 / s, crop_model_equivalent_km2=eq,
                        share_of_rain_fed_land=eq / cap,
                        billions_the_rain_feeds={f'kept_{int(k * 100)}': cap * (1 - k) / eq
                                                 for k in basis['kept_shares']})
        res = both(one)
        res['billions_the_rain_feeds'] = {k: sorted(v) for k, v in res['billions_the_rain_feeds'].items()}
        out[name] = res
    return out


def night_fields(basis):
    s = basis['stores']
    return dict(
        strategy=('Evergreen crops store carbon by day and idle at a quarter of their daytime upkeep once their own '
                  'carbon balance turns negative; the day fruit is set at sunrise and ripe at sunset, so the harvest '
                  'comes once a cycle at dusk and the towns hold a cycle of food.'),
        store_g_c_m2=dict(earth_like=s['earth_like'], idle_quarter=s['idle_quarter']),
        night_upkeep_w_m2=s['night_upkeep_w_m2'], dark_hours=s['dark_hours'], night_energy_twh_per_billion=0.0)


def irrigated_fields(basis, diet='plant_rich'):
    """Field crops watered with water lifted from the fresh seas and lakes, to the crop model's yield or to what the
    absorbed sunlight lets them transpire, whichever is less: topping up the rain belt's rain-fed fields, or watering
    fields in the fog desert."""
    rain, w, z = basis['rain'], basis['water'], basis['zones']['zones']
    out = dict(what=('Field crops watered to the crop model\'s yield; the water is lifted from the seas and lakes to '
                     'the zone\'s mean height and piped its mean distance.'))
    for name, zone, band in (('rain_belt_top_up', 'rain_belt', '0-15'), ('fog_desert', 'fog_desert', '45-60')):
        zz = z[zone]

        def one(end, zz=zz, band=band):
            eq = eq_km2(basis, diet, end)
            grows = pick(rain['watered_yield_cap'][band], end, low_is_good=False)          # eq per km2 watered
            s = pick(rain['support_scaled'][band], end, low_is_good=False) / grows          # rain's share of it
            net = yearly(crop_water_mm(basis, end), eq) * (1.0 - s)
            km3 = net / pick(w['application_efficiency'], end, low_is_good=False)
            kwh = pumping_kwh_m3(basis, zz['height_above_sea_km_mean'], zz['km_to_sea_or_great_lake_mean'], end,
                                 emitters=True)
            power = tw_of(km3, kwh)
            return dict(land_km2=eq / grows, yield_against_equator=grows,
                        energy_limited_share=grows / rain['gpp_factor'][band], rain_supplies=s,
                        net_irrigation_km3_year=net,
                        withdrawal_km3_year=km3, share_of_river_flow=km3 / w['river_km3_year'], kwh_m3=kwh,
                        power_tw=power, heat=heat(power, basis, end))
        out[name] = both(one)
    out.update(both(lambda end: field_nutrients(basis, end)))
    out['night'] = night_fields(basis)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Glass.

def glass_yield_factor(basis, end, closed):
    """Yield per m2 of glass against the open field at the same latitude: the house's light transmission raised to
    the crop's light response, times the CO2 enrichment the house can hold."""
    g = basis['glass']
    tau = pick(g['transmission'], end, low_is_good=False)
    k = pick(g['light_response'], end)
    gain = pick(g['co2_gain_closed'] if closed else g['co2_gain_vented'], end, low_is_good=False)
    return tau ** k * gain


def glasshouse(basis, zone, diet='plant_rich', closed=False, end=None):
    """Sunlit hydroponic glasshouses in a zone growing the diet at the crop model's yields, with Venlo-type
    structure. Vented houses lose the crop's transpiration and replace it from rain on the glass and the seas;
    closed houses condense and reuse it and pump out the absorbed sunlight with heat pumps."""
    g, w = basis['glass'], basis['water']
    zz = basis['zones']['zones'][zone]
    lat_factor = basis['light']['lat_yield'](zz['abs_lat_deg_mean'])

    def one(end):
        eq = eq_km2(basis, diet, end)
        gf = glass_yield_factor(basis, end, closed)
        # the CO2 held in the house raises the crops' water-use efficiency; the sunlight absorbed inside caps what
        # they transpire, and with it what they grow
        need = crop_water_mm(basis, end) / pick(g['co2_water_use_factor_closed' if closed else
                                                  'co2_water_use_factor_vented'], end, low_is_good=False)
        grows = min(lat_factor * gf,
                    pick(basis['light']['evaporation_cap_inside_mm_day'][zone], end, low_is_good=False) / need)
        area = eq / grows
        m2 = area * KM2
        transp = yearly(need, eq)
        roof = yearly(zz['rain_mm_day'] * pick(g['roof_capture'], end, low_is_good=False), area)
        produce = pick(basis['diet']['produce_water_km3_per_billion'][diet], end)
        consumed = produce if closed else transp
        from_seas = max(consumed - roof, 0.0)
        kwh = pumping_kwh_m3(basis, zz['height_above_sea_km_mean'], zz['km_to_sea_or_great_lake_mean'], end)
        operations = m2 * pick(g['operations_kwh_m2_yr'], end) * KWH_J / SECONDS_PER_YEAR / 1e12
        cooling = (m2 * pick(basis['light']['absorbed_inside_w_m2'][zone], end, low_is_good=False) /
                   pick(g['cooling_cop'], end, low_is_good=False) / 1e12) if closed else 0.0
        water = tw_of(from_seas, kwh)
        power = operations + cooling + water
        air_kwh = 1e3 / pick(w['air_extraction_l_per_kwh'], end, low_is_good=False)
        nets_l_m2 = pick(basis['fog_net_l_m2_year'], end, low_is_good=False)
        out = dict(glass_km2=area, share_of_zone=area / zz['dry_land_km2'], yield_against_open_field=gf,
                   yield_against_equator_field=grows, energy_limited_share=grows / (lat_factor * gf),
                   transpiration_km3_year=transp, consumed_km3_year=consumed, rain_on_glass_km3_year=roof,
                   from_seas_km3_year=from_seas,
                   power_tw=dict(operations=operations, cooling=cooling, water=water, total=power),
                   heat=heat(power, basis, end),
                   glass_gt=m2 * g['glass_kg_m2'] / 1e12, steel_gt=m2 * g['steel_kg_m2'] / 1e12,
                   aluminium_gt=m2 * g['aluminium_kg_m2'] / 1e12,
                   glass_renewed_mt_year=m2 * g['glass_kg_m2'] / g['life_years'] / 1e9,
                   if_from_seas_were_drawn_from_air_tw=tw_of(from_seas, air_kwh),
                   fog_net_km2_for_the_sea_water=(from_seas * 1e12 / nets_l_m2 / KM2) if nets_l_m2 > 0 else None)
        out.update(soilless_nutrients(basis, end))
        return out
    if end is not None:
        return one(end)
    res = both(one)
    res['latitude_yield_factor'] = lat_factor
    res['zone'] = zone
    return res


def soilless_nutrients(basis, end, people_billion=1.0):
    """Nitrogen and phosphorus of recirculating soilless culture growing the diet: the lunar design's applied
    nitrogen per person with the soilless emission factor for N2O, the people's phosphorus at the lunar design's
    losses, and the share of the nutrient solution a house discharges."""
    s, p, n = basis['soilless'], basis['phosphorus'], basis['nitrogen']
    applied = pick(n['applied_kg'], end) * people_billion                            # Mt N a year
    eaten = pick(p['eaten_kg'], end) * people_billion
    harvest = eaten * pick(basis['diet']['loss_factor'], end)
    return dict(
        phosphorus_mt_year=dict(eaten=eaten, in_harvest=harvest, lost=pick(p['lost_kg'], end) * people_billion,
                                discharged=harvest * pick(s['p_discharge_share'], end)),
        nitrogen=dict(applied_mt_n_year=applied, n2o=n2o(applied * pick(s['n2o_share_of_applied'], end) * 1e3, basis)),
        methane=ch4(pick(basis['methane']['kg_person'], end) * people_billion, basis, end))


def night_lamps_for_glass(basis, zone, diet='plant_rich'):
    """Lighting a glasshouse crop through its dark hours at a growth level: the power over the cycle per m2, the
    photons it adds against the cycle's sunlight under the glass, and the cost per billion if that share of the diet
    grows under the lamps."""
    g, lamp = basis['glass'], basis['lamps']
    zz = basis['zones']['zones'][zone]
    dark = basis['light']['dark_hours'](zz['abs_lat_deg_mean'])
    sun_mol = basis['light']['par_mol_m2_cycle'](zz['abs_lat_deg_mean'])

    def one(end):
        ppfd = pick(lamp['night_ppfd'], end, low_is_good=False)
        eff = pick(lamp['efficacy_umol_j'], end, low_is_good=False)
        w_m2_on = ppfd / eff
        w_m2_cycle = w_m2_on * dark / CYCLE_HOURS
        tau = pick(g['transmission'], end, low_is_good=False)
        added = ppfd * 1e-6 * dark * 3600.0 / (sun_mol * basis['diet']['cloud_share'] * tau)
        glass = glasshouse(basis, zone, diet, end=end)['glass_km2']
        lit_glass = glass / (1.0 + added)                 # the lamps' photons grow food as the sun's do
        power = lit_glass * KM2 * w_m2_cycle / 1e12
        return dict(dark_hours=dark, w_m2_while_lit=w_m2_on, w_m2_over_the_cycle=w_m2_cycle,
                    photons_added_against_sunlight_under_glass=added, glass_km2_with_lamps=lit_glass,
                    glass_saved_km2=glass - lit_glass, power_tw=power, heat=heat(power, basis, end))
    return both(one)


# ---------------------------------------------------------------------------------------------------------------
# Lamps.

def vertical_farms(basis, diet='plant_rich', end=None):
    """Indoor farms lit only by lamps, growing the whole diet: photons for the diet's food energy at a staple crop's
    edible light-use efficiency under lamps, electricity at the fixtures' photon efficacy and lighting's share of the
    load, shelves at the lamps' light, floor at the stated layers and building mass at the summit tower study's floor
    load."""
    lamp = basis['lamps']

    def one(end, info):
        kcal_mol = pick(info['edible_g_per_mol'], end, low_is_good=False) * info['kcal_per_g_dry']
        mol_year = diet_kcal(basis, diet, end) * JULIAN_YEAR_DAYS / kcal_mol                # per resident
        eff = pick(lamp['efficacy_umol_j'], end, low_is_good=False)
        share = pick(lamp['lighting_share'], end, low_is_good=False)
        kw = mol_year / (eff * 1e-6) / SECONDS_PER_YEAR / share / 1e3
        shelf = mol_year / (info['ppfd'] * 1e-6 * info['photoperiod_h'] * 3600.0 * JULIAN_YEAR_DAYS)
        floor = shelf / pick(lamp['layers'], end, low_is_good=False)
        power = kw * B / 1e9
        return dict(power_kw_per_resident=kw, power_tw=power, heat=heat(power, basis, end),
                    photons_mol_per_resident_year=mol_year, kcal_per_mol=kcal_mol,
                    kwh_per_kcal=1.0 / (kcal_mol * eff * 1e-6) / KWH_J / share,
                    shelf_km2=shelf * B / KM2, floor_km2=floor * B / KM2,
                    building_gt=floor * B * lamp['floor_kg_m2'] / 1e12,
                    night_energy_twh=power * 1e3 * basis['night']['hours'] / 1e3,
                    night_energy_over_european_night=kw * basis['night']['hours'] /
                    basis['night']['european_final_kwh'],
                    supply_land_km2={k: power * pick(f, end) for k, f in basis['supply_land_km2_per_tw'].items()})
    if end is not None:
        return {staple: one(end, info) for staple, info in lamp['staples'].items()}
    return {staple: both(lambda e, info=info: one(e, info)) for staple, info in lamp['staples'].items()}


def lamp_produce(basis, grams_day, fresh_g_per_mol, end):
    """kW a resident to grow grams_day of fresh produce (before losses) under lamps at fresh_g_per_mol."""
    lamp = basis['lamps']
    grams = grams_day * pick(basis['diet']['loss_factor'], end)
    mol_day = grams / pick(fresh_g_per_mol, end, low_is_good=False)
    eff = pick(lamp['efficacy_umol_j'], end, low_is_good=False)
    share = pick(lamp['lighting_share'], end, low_is_good=False)
    return mol_day / 86400.0 / (eff * 1e-6) / share / 1e3


# ---------------------------------------------------------------------------------------------------------------
# Power-made protein.

def microbial_protein(basis, protein_g_day, end=None):
    """Hydrogen-oxidising bacteria on electrolytic hydrogen, CO2 and ammonia supplying protein_g_day of protein eaten
    a resident, grown with the diet's losses: electricity per kg protein across three process accounts, the ammonia's
    electricity where an account leaves it out, and the biomass's phosphorus and nitrogen."""
    mp = basis['microbial']

    def one(end):
        kg = protein_g_day * pick(basis['diet']['loss_factor'], end) * 1e-3 * JULIAN_YEAR_DAYS * B   # kg grown
        kwh = pick(mp['kwh_per_kg_protein'], end) + pick(mp['ammonia_kwh_per_kg_protein'], end)
        power = kg * kwh * KWH_J / SECONDS_PER_YEAR / 1e12
        biomass = kg / pick(mp['protein_share'], end, low_is_good=False)
        return dict(kwh_per_kg_protein=kwh, power_tw=power, heat=heat(power, basis, end),
                    biomass_mt_year=biomass / 1e9, co2_mt_year=biomass * mp['co2_kg_per_kg_dw'] / 1e9,
                    phosphorus_mt_year=biomass * pick(mp['p_share_of_dw'], end) / 1e9,
                    nitrogen_mt_year=kg / mp['protein_per_n'] / 1e9,
                    night_energy_twh=power * basis['night']['hours'])
    if end is not None:
        return one(end)
    res = both(one)
    res['protein_g_day'] = protein_g_day
    return res


# ---------------------------------------------------------------------------------------------------------------
# Fish and the seas.

def aquaponics(basis, protein_g_day, end=None):
    """Fish in recirculating tanks fed on crops, their water feeding vegetables under glass: the fish for
    protein_g_day of fillet protein eaten a resident, grown with the diet's losses, their feed and its crop land,
    power, water, and where the feed's nitrogen and phosphorus go."""
    a = basis['aquaponics']

    def one(end):
        fillet = protein_g_day * pick(basis['diet']['loss_factor'], end) / a['fillet_protein_share']   # g a day
        live = fillet / pick(a['fillet_yield'], end, low_is_good=False)
        feed = live * pick(a['fcr'], end)
        kg_live = live * 1e-3 * JULIAN_YEAR_DAYS * B
        kg_feed = feed * 1e-3 * JULIAN_YEAR_DAYS * B
        feed_kcal_day = feed * a['feed_kcal_per_g']
        feed_eq = feed_kcal_day * pick(basis['diet']['variety_factor'], end) / \
            pick(basis['diet']['crop_kcal_m2_day'], end, low_is_good=False) * B / KM2
        power = kg_live * pick(a['kwh_per_kg_fish'], end) * KWH_J / SECONDS_PER_YEAR / 1e12
        feed_n = kg_feed * pick(a['feed_protein_share'], end) / a['protein_per_n']
        feed_p = kg_feed * pick(a['feed_p_share'], end)
        return dict(fish_live_mt_year=kg_live / 1e9, feed_mt_year=kg_feed / 1e9,
                    feed_crop_equivalent_km2=feed_eq, power_tw=power, heat=heat(power, basis, end),
                    water_km3_year=kg_live * pick(a['water_l_per_kg_fish'], end) * 1e-3 / 1e9,
                    feed_nitrogen_mt_year=feed_n / 1e9,
                    n2o=n2o(feed_n * pick(a['n2o_share_of_n'], end) / 1e6, basis),
                    feed_phosphorus_mt_year=feed_p / 1e9,
                    phosphorus_to_sludge_mt_year=feed_p * pick(a['p_share_sludge'], end) / 1e9)
    if end is not None:
        return one(end)
    res = both(one)
    res['protein_g_day'] = protein_g_day
    return res


def sea_farms(basis, protein_g_day, end=None):
    """Farms in the fresh seas supplying protein_g_day of protein eaten a resident, grown with the diet's losses:
    macroalgae or filter-feeding shellfish, with the sea area they take and the nitrogen and phosphorus their harvest
    draws from the seas."""
    s = basis['sea']
    out = {}

    def grown_kg(end):
        return protein_g_day * pick(basis['diet']['loss_factor'], end) * 1e-3 * JULIAN_YEAR_DAYS * B

    def seaweed(end):
        kg_protein = grown_kg(end)
        protein_share = pick(s['seaweed_protein_share'], end, low_is_good=False)
        dw = kg_protein / protein_share
        area = dw / (pick(s['seaweed_t_dw_ha_yr'], end, low_is_good=False) * 1e3) / 100.0       # km2
        return dict(dry_mass_mt_year=dw / 1e9, sea_km2=area,
                    share_of_seas_shallower_than_50_m=area / basis['zones']['sea_km2_by_depth_m']['50'],
                    share_of_seas_shallower_than_200_m=area / basis['zones']['sea_km2_by_depth_m']['200'],
                    phosphorus_drawn_mt_year=dw * pick(s['seaweed_p_share'], end) / 1e9,
                    nitrogen_drawn_mt_year=dw * pick(s['seaweed_n_share'], end) / 1e9)

    def shellfish(end):
        kg_protein = grown_kg(end)
        protein_per_t = (pick(s['mussel_meat_share'], end, low_is_good=False) *
                         pick(s['mussel_meat_protein_share'], end, low_is_good=False) * 1e3)      # kg per t live
        t_live = kg_protein / protein_per_t
        yield_t_ha = pick(s['mussel_t_live_ha_yr'], end, low_is_good=False)
        area = t_live / yield_t_ha / 100.0
        return dict(live_mt_year=t_live / 1e6, sea_km2=area,
                    share_of_seas_shallower_than_50_m=area / basis['zones']['sea_km2_by_depth_m']['50'],
                    share_of_seas_shallower_than_200_m=area / basis['zones']['sea_km2_by_depth_m']['200'],
                    phosphorus_drawn_mt_year=t_live * pick(s['mussel_p_kg_per_t'], end) / 1e9,
                    nitrogen_drawn_mt_year=t_live * pick(s['mussel_n_kg_per_t'], end) / 1e9)
    def against_waters(farm):
        """The draw a billion residents' farm makes against what the waters' life can spare, and the billions the
        waters could supply: [low-cost end, high-cost end] for the draw, smallest first for the billions."""
        def one(end):
            draw_kt = farm(end)['phosphorus_drawn_mt_year'] * 1e3
            ceiling = pick(basis['seas_phosphorus']['all_waters_ceiling_kt_year'], end, low_is_good=False)
            common = pick(basis['seas_phosphorus']['common_level_ceiling_kt_year'], end, low_is_good=False)
            return dict(draw_over_all_waters_ceiling=draw_kt / ceiling,
                        draw_over_common_level_ceiling=draw_kt / common,
                        billions_the_waters_supply=ceiling / draw_kt)
        res = both(one)
        res['billions_the_waters_supply'] = sorted(res['billions_the_waters_supply'])
        return res

    if end is not None:
        return dict(seaweed=seaweed(end), shellfish=shellfish(end))
    out['seaweed'] = both(seaweed)
    out['shellfish'] = both(shellfish)
    out['against_the_waters_phosphorus'] = dict(seaweed=against_waters(seaweed), shellfish=against_waters(shellfish))
    out['protein_g_day'] = protein_g_day
    return out


# ---------------------------------------------------------------------------------------------------------------
# Water for the fog desert.

def weibull_above(mean, shape, cut):
    """Mean of the wind speed counted only above cut, E[U; U > cut], for Weibull-distributed speeds of the given mean
    and shape."""
    if mean <= 0.0:
        return 0.0
    a = 1.0 + 1.0 / shape
    scale = mean / gamma(a)
    return mean * gammaincc(a, (cut / scale) ** shape)


def bin_fog_flux(fog_share, mean, shape, cut, pairing):
    """Mean of fog x wind speed over one bin of hour angle, counting wind above cut only, with the bin's instantaneous
    speeds spread as a Weibull distribution about its mean. 'independent' pairs fog and wind at random within the
    bin; 'windiest' puts the fog in the windiest share of the bin, the most fog-net water the bin's means allow."""
    if fog_share <= 0.0 or mean <= 0.0:
        return 0.0
    if pairing == 'independent':
        return fog_share * weibull_above(mean, shape, cut)
    scale = mean / gamma(1.0 + 1.0 / shape)
    start = scale * (-math.log(fog_share)) ** (1.0 / shape) if fog_share < 1.0 else 0.0
    return weibull_above(mean, shape, max(cut, start))


FOG_NET_METHOD = ('The ring products give fog and the mean 10 m wind speed in 10-degree bins of hour angle over every '
                  'land column and snapshot, the speeds to 0.1 m/s, and hold no pairing of fog and wind snapshot by '
                  'snapshot. Each bin\'s instantaneous speeds are spread as a Weibull distribution about its mean. The '
                  'low-yield end takes a Rayleigh spread (shape 2), fog and wind independent within the bin, the bin '
                  'means rounded down and the poorer ring; the high-yield end takes each ring\'s whole-cycle spread '
                  '(shape fitted to its median and 90th-percentile wind, day-night swing included), the fog in the '
                  'windiest share of each bin, the means rounded up and the better ring. The bracket takes the two '
                  'rings on the fog desert\'s latitudes, 45 N and 45 S; the tilted rings are shown beside them.')


def fog_desert_water(basis):
    """Where a fog-desert farm's water can come from, per m2 of collector or per m3: rain on glass, fog nets in the
    design climate's fogs and winds (CM1's land rings), dew, extraction from the air, and pipes from the seas."""
    w, fog = basis['water'], basis['fog']
    zz = basis['zones']['zones']['fog_desert']
    dt = CYCLE_HOURS / len(fog['rings'][0]['fog']) * 3600.0
    half_step = fog['wind_rounding_m_s'] / 2.0

    def flux(ring, end, cut):
        """m of foggy air through a metre of net a cycle, counting wind above cut."""
        if end == 'low':
            shape, pairing, shift = ring['shape_whole_cycle'], 'windiest', half_step
        else:
            shape, pairing, shift = fog['rayleigh_shape'], 'independent', -half_step
        return dt * sum(bin_fog_flux(f, max(u + shift, 0.0), shape, cut, pairing)
                        for f, u in zip(ring['fog'], ring['wind']))

    def litres(m_per_cycle, end):
        lwc = pick(fog['lwc_g_m3'], end, low_is_good=False)
        eff = pick(fog['net_efficiency'], end, low_is_good=False)
        return m_per_cycle * lwc * eff * CYCLES_PER_YEAR / 1e3                           # kg (L) per m2 a year

    rings = {}
    for ring in fog['rings']:
        f_sum = sum(ring['fog'])
        info = dict(in_fog_desert_latitudes=ring.get('in_fog_desert_latitudes', True),
                    fog_hours_per_cycle=f_sum * dt / 3600.0,
                    mean_wind_in_fog_m_s=(sum(f * u for f, u in zip(ring['fog'], ring['wind'])) / f_sum) if f_sum else 0.0,
                    weibull_shape_whole_cycle=ring['shape_whole_cycle'])
        for key, cut in (('all_wind', 0.0), ('wind_above_cutoff', fog['stops_below_m_s'])):
            m = [flux(ring, 'low', cut), flux(ring, 'high', cut)]
            info[f'foggy_air_m_per_cycle_{key}'] = sorted(m)
            info[f'l_per_m2_year_{key}'] = sorted([litres(m[0], 'low'), litres(m[1], 'high')])
        rings[ring['name']] = info
    in_desert = [rings[r['name']] for r in fog['rings'] if r.get('in_fog_desert_latitudes', True)]
    yields = {f'l_per_m2_year_{key}': [min(r[f'l_per_m2_year_{key}'][0] for r in in_desert),
                                       max(r[f'l_per_m2_year_{key}'][1] for r in in_desert)]
              for key in ('all_wind', 'wind_above_cutoff')}
    out = dict(fog_nets=dict(method=FOG_NET_METHOD, rings=rings, yields_low_to_high=yields,
                             earth_operational_l_m2_year=[x * JULIAN_YEAR_DAYS for x in fog['earth_operational_l_m2_day']]))
    out['rain'] = dict(mm_day=zz['rain_mm_day'], l_per_m2_glass_year=[zz['rain_mm_day'] * JULIAN_YEAR_DAYS * c for c in
                                                                      basis['glass']['roof_capture']])
    out['dew'] = dict(earth_clear_night_l_m2=basis['dew']['earth_l_m2_night'],
                      night_land_evaporation_mm_h=fog['night_land_evaporation_mm_h'])
    out['from_air'] = dict(l_per_kwh=w['air_extraction_l_per_kwh'],
                           kwh_m3=[1e3 / w['air_extraction_l_per_kwh'][1], 1e3 / w['air_extraction_l_per_kwh'][0]])
    out['from_seas'] = both(lambda end: dict(
        lift_kwh_m3=lift_kwh_m3(zz['height_above_sea_km_mean'], pick(w['pump_efficiency'], end, low_is_good=False)),
        pipe_kwh_m3=pipe_kwh_m3(zz['km_to_sea_or_great_lake_mean'], end),
        total_kwh_m3=pumping_kwh_m3(basis, zz['height_above_sea_km_mean'], zz['km_to_sea_or_great_lake_mean'], end)))
    out['from_seas']['height_km'] = zz['height_above_sea_km_mean']
    out['from_seas']['distance_km'] = zz['km_to_sea_or_great_lake_mean']
    out['from_seas']['earth_aqueducts_kwh_m3'] = w['earth_aqueducts_kwh_m3']
    return out


def precision_fermentation(basis, protein_g_day, end=None):
    """Proteins made by a fungus on sugar (the egg-white protein case), supplying protein_g_day eaten a resident,
    grown with the diet's losses: the process electricity and heat, and the sugar from crops at the crop model's yield
    or from electrochemistry."""
    f, es = basis['fermentation'], basis['electro_sugar']

    def one(end):
        kg = protein_g_day * pick(basis['diet']['loss_factor'], end) * 1e-3 * JULIAN_YEAR_DAYS * B
        process_tw = kg * (f['kwh_per_kg_protein'] + f['heat_kwh_per_kg_protein']) * KWH_J / SECONDS_PER_YEAR / 1e12
        sugar_kg = kg * f['glucose_kg_per_kg_protein']
        sugar_kcal = sugar_kg * 1e3 * f['glucose_kcal_per_g']
        crop_km2 = sugar_kcal / (pick(basis['diet']['crop_kcal_m2_day'], end, low_is_good=False) * JULIAN_YEAR_DAYS) / KM2
        electro_tw = sugar_kcal * KCAL_J / pick(es['power_to_food'], end, low_is_good=False) / SECONDS_PER_YEAR / 1e12
        return dict(process_power_tw=process_tw, sugar_mt_year=sugar_kg / 1e9, sugar_crop_equivalent_km2=crop_km2,
                    electro_sugar_power_tw=electro_tw, total_power_with_electro_sugar_tw=process_tw + electro_tw)
    if end is not None:
        return one(end)
    res = both(one)
    res['protein_g_day'] = protein_g_day
    return res


def power_made_food(basis, end=None):
    """A diet made on power without farmland: microbial protein for half the reference diet's protein (the most its
    nucleic acids allow), fungal protein from precision fermentation on sugar made from power for the other half,
    sugar made by electrochemistry for the food energy, and vegetables and fruit from lamp-lit farms."""
    d, es = basis['diet'], basis['electro_sugar']
    half = d['reference_diet_protein_g_day'] / 2.0

    def one(end):
        microbes = microbial_protein(basis, half, end=end)
        microbial_kw = microbes['power_tw']                                               # TW a billion = kW each
        fermented_kw = precision_fermentation(basis, half, end=end)['total_power_with_electro_sugar_tw']
        protein_kcal = d['reference_diet_protein_g_day'] * pick(d['protein_kcal_per_g'], end, low_is_good=False)
        sugar_kcal = d['intake_kcal_day'] - d['lamp_produce_kcal_day'] - protein_kcal
        food_w = sugar_kcal * pick(d['loss_factor'], end) * KCAL_J / 86400.0
        sugar_kw = food_w / pick(es['power_to_food'], end, low_is_good=False) / 1e3
        produce_kw = (lamp_produce(basis, d['vegetables_g_day'], basis['lamps']['leafy_fresh_g_per_mol'], end) +
                      lamp_produce(basis, d['fruit_g_day'], basis['lamps']['fruit_fresh_g_per_mol'], end))
        kw = microbial_kw + fermented_kw + sugar_kw + produce_kw
        power = kw * B / 1e9
        return dict(microbial_protein_kw_per_resident=microbial_kw, fermented_protein_kw_per_resident=fermented_kw,
                    protein_kw_per_resident=microbial_kw + fermented_kw, food_energy_kw_per_resident=sugar_kw,
                    vegetables_and_fruit_kw_per_resident=produce_kw, total_kw_per_resident=kw, power_tw=power,
                    heat=heat(power, basis, end), night_energy_twh=power * basis['night']['hours'],
                    sugar_kcal_day=sugar_kcal, protein_kcal_day=protein_kcal,
                    phosphorus_in_microbial_biomass_mt_year=microbes['phosphorus_mt_year'])
    if end is not None:
        return one(end)
    res = both(one)
    res['note'] = ('Sugar supplies the food energy the protein and the produce leave; the fats and the many plant '
                   'compounds of a varied diet have no measured electrochemical route.')
    return res
