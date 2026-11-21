"""The settlements from their designs: what a resident of a crater town or of a lava-tube city takes to build and to run,
worked out for the year it is built in from the design studies' own quantities and that year's technology.

Two design products are read, not their prices in U:
  research/studies/vacuum_crater_buildout/results/first_hundred.json   the first hundred crater towns, 500 residents each
  research/studies/lava_tube_cities/results/buildout.json              four lava-tube cities of 50,000 and those hundred
                                                                        towns, with the region's farms and light industry

A resident's build is carried in the streams the designs' own schedules use: regolith handled, finished bulk products,
finished metal, pressure-system components, equipment (the designs' imports), and the first air and water. The city
halls' streams are the founding system's less the hundred towns'; they are shared over the halls' floor, so that a
resident takes the halls of a city home and, wherever they live, the region's farms and light industry. Towns built
after the founding system put those regional halls on the surface under three metres of regolith, since no cave is
assumed for them; the founding towns use the cities' cave halls, as their design does.
"""
from __future__ import annotations
from functools import lru_cache
import json
from pathlib import Path

from . import technology as tech

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TOWNS = ROOT / 'research/studies/vacuum_crater_buildout/results/first_hundred.json'
CITIES = ROOT / 'research/studies/lava_tube_cities/results/buildout.json'
STREAMS = ('handling_kg', 'bulk_kg', 'metal_kg', 'pressure_kg', 'equipment_kg', 'volatiles_kg')


@lru_cache(maxsize=None)
def inputs():
    return json.loads((HERE / 'period_inputs.json').read_text())


def v(entry):
    return entry['value'] if isinstance(entry, dict) else entry


@lru_cache(maxsize=None)
def base(town_case='500'):
    """The designs' quantities per resident and per square metre of hall floor, at their own engineering (2041). The
    city halls are always the founding system's less its hundred 500-resident towns; `town_case` names the first
    hundred's case whose towns are carried (500, 750 or 1000 residents)."""
    town = json.loads(TOWNS.read_text())['cases'][town_case]['quantities']
    founding_towns = json.loads(TOWNS.read_text())['cases']['500']['quantities']
    system = json.loads(CITIES.read_text())['cases']['reference']
    n_town = town['population']
    n_all = system['population']
    n_city = system['city_population']
    sd = system['schedule']['production_demand']

    def streams(q):
        d = q['production_demand']
        out = dict(handling_kg=d['handling_kg'], bulk_kg=d['bulk_finished_kg'], metal_kg=d['metal_finished_kg'],
                   pressure_kg=d['pressure_system_kg'], equipment_kg=d['special_import_kg'], volatiles_kg=d['volatile_import_kg'])
        s = q['special_imports']
        out['equipment_kg'] -= s['power_hardware_kg'] + s['reserve_hardware_kg'] + s['construction_plant_kg']
        return out
    town_streams = streams(town)
    founding_streams = streams(founding_towns)
    system_streams = dict(handling_kg=sd['handling'], bulk_kg=sd['local_bulk'], metal_kg=sd['local_metal'],
                          pressure_kg=sd['pressure_components'], equipment_kg=sd['other_imports'], volatiles_kg=sd['volatile_imports'])
    # The model sizes the power supply and the industry's plant itself, so the designs' power and storage hardware and
    # the towns' construction plant leave their equipment: the towns' own power, reserve and construction-plant imports,
    # and the city system's integrated power and storage (their yearly renewal imports times their lives).
    special = founding_towns['special_imports']
    renew, lives = system['resources']['renewal_import_kg_year'], system['maintenance']['life_years']
    system_streams['equipment_kg'] -= renew['new_power'] * lives['new_power'] + renew['new_storage'] * lives['new_storage']
    system_streams['equipment_kg'] -= special['power_hardware_kg'] + special['reserve_hardware_kg'] + special['construction_plant_kg']
    floors = system['geometry']['floor_demand_m2']
    hall_floor = system['geometry']['required_floor_m2']
    halls = {k: (system_streams[k] - founding_streams[k]) / hall_floor for k in STREAMS}
    living = floors['ordinary'] + floors['parks'] + floors['municipal']
    food = system['food']
    loads = system['utilities']['mean_load_components_W']
    resources = system['resources']
    ops = resources['operating_import_kg_year']
    s = inputs()['settlements']
    # The town's shielding is handled and placed but neither made into products nor renewed; the city halls' grading
    # likewise. Their shares of the handling stay fixed when structures grow lighter.
    town_fixed = town['shielding_kg'] / town_streams['handling_kg']
    hall_fixed = system['geometry']['grading_kg'] / (system_streams['handling_kg'] - founding_streams['handling_kg'])
    return dict(
        town=dict(per_resident={k: town_streams[k] / n_town for k in STREAMS}, fixed_handling_share=town_fixed,
                  floor_m2=v(s['town_floor_m2_resident']), gardens_W=town['services']['garden_average_W'] / n_town),
        halls=dict(per_m2=halls, fixed_handling_share=hall_fixed),
        city_living_m2=living / n_city,
        shared_m2=dict(agriculture=floors['agriculture'] / n_all, aquaculture=floors['aquaculture'] / n_all,
                       fermentation=floors['fermentation'] / n_all, food_handling=floors['food_handling'] / n_all,
                       light_industry=floors['light_industry'] / n_all),
        city_gardens_W=(loads['civic_gardens'] - town['services']['garden_average_W']) / n_city,
        loads_W=dict(food=loads['food'] / n_all, basic=loads['basic_services'] / n_all, industry=loads['everyday_industry'] / n_all),
        cooling_share=loads['cooling'] / (loads['food'] + loads['basic_services'] + loads['civic_gardens'] + loads['everyday_industry']),
        firm_over_bus=system['utilities']['occupied_firm_W'] / system['utilities']['mean_electric_bus_W'],
        water_makeup_kg_year=ops['water'] / n_all,
        landscape_water_kg_year=v(s['landscape_water_kg_day']) * 365.25,
        air_makeup_kg_year=ops['air'] / n_all,
        earth_supplies_kg_year=(ops['biosphere'] + ops['food_nutrients_and_carbon'] + ops['everyday_carbon_feed'] + ops['specialist_consumables']) / n_all,
        local_goods_kg_year=resources['local_goods_makeup_kg_year'] / n_all,
        operators=(resources['food_operators'] + resources['industry_operators']) / n_all,
        diet_kcal_day=v(s['diet_kcal_day']),
        microbial_share_design=food['microbial_counted_kcal_person_day'] / v(s['diet_kcal_day']),
        population=dict(towns=n_town, cities=n_city, all=n_all),
        shielding_t_m2=v(s['regional_hall_shielding_t_per_m2_floor']))


def food_W(case, year):
    """Food's average electrical load per resident: crops under lamps for the share of the diet they supply and
    electro-microbial food for the rest, scaled from the founding design's food load at its lamps and crops."""
    b = base()
    t = tech.Year(case, year)
    t0 = tech.Year('static', year)
    kcal = b['diet_kcal_day']
    # The design's food load less its electro-microbial food at the 2041 efficiency is the crops' share.
    microbial0_W = kcal * t0.microbial_diet_share * 4184.0 / t0.microbial_efficiency / 86400.0
    crop_J_kcal0 = (b['loads_W']['food'] - microbial0_W) * 86400.0 / (kcal * (1 - t0.microbial_diet_share))
    crop_J_kcal = crop_J_kcal0 * (t0.led_umol_J / t.led_umol_J) * (t0.crop_kcal_per_mol / t.crop_kcal_per_mol)
    microbial_J_kcal = 4184.0 / t.microbial_efficiency
    share = t.microbial_diet_share
    return (kcal * (1 - share) * crop_J_kcal + kcal * share * microbial_J_kcal) / 86400.0


def shared_floor_m2(case, year):
    """The region's farms and light industry per resident: crops follow the share of the diet grown and the food a
    square metre yields at the design's light; fermentation follows the microbial share."""
    b = base()
    t, t0 = tech.Year(case, year), tech.Year('static', year)
    s = b['shared_m2']
    crops = s['agriculture'] * (1 - t.microbial_diet_share) / (1 - t0.microbial_diet_share) * t0.crop_kcal_per_mol / t.crop_kcal_per_mol
    ferment = s['fermentation'] * t.microbial_diet_share / t0.microbial_diet_share
    return dict(agriculture=crops, aquaculture=s['aquaculture'], fermentation=ferment, food_handling=s['food_handling'],
                light_industry=s['light_industry'])


def build(kind, case, year, shared_in_cave):
    """What one resident of a kind ('town' or 'city') takes to build in a calendar year, in kilograms by stream, with
    the floor. `shared_in_cave` puts the regional halls in a cave (the founding system and cities, whose towns are the
    first hundred's 500-resident towns) or on the surface under shielding (towns after the founding system, of the size
    period_inputs.json sets)."""
    b = base('500' if (shared_in_cave or kind == 'city') else str(v(inputs()['settlements']['town_residents'])))
    t = tech.Year(case, year)
    lighter = t.built_mass_multiple

    def scaled(streams, fixed_share):
        out = dict(streams)
        out['handling_kg'] = streams['handling_kg'] * (fixed_share + (1 - fixed_share) * lighter)
        for k in ('bulk_kg', 'metal_kg', 'pressure_kg'):
            out[k] = streams[k] * lighter
        return out

    halls = b['halls']
    shared = sum(shared_floor_m2(case, year).values())
    hall_floor = shared + (b['city_living_m2'] if kind == 'city' else 0.0)
    out = scaled({k: halls['per_m2'][k] * hall_floor for k in STREAMS}, halls['fixed_handling_share'])
    if not shared_in_cave:
        out['handling_kg'] += shared * b['shielding_t_m2'] * 1000.0
    if kind == 'town':
        town = scaled(b['town']['per_resident'], b['town']['fixed_handling_share'])
        out = {k: out[k] + town[k] for k in STREAMS}
    out['floor_m2'] = hall_floor + (b['town']['floor_m2'] if kind == 'town' else 0.0)
    out['hall_floor_m2'] = hall_floor
    return out


def needs(kind, case, year):
    """What one resident needs each year: average power by end use, the electric bus and the firm power to be
    supplied; water, air and supplies; local everyday goods; the settlement's own operators."""
    b = base()
    t, t0 = tech.Year(case, year), tech.Year('static', year)
    gardens = (b['town']['gardens_W'] if kind == 'town' else b['city_gardens_W']) * t0.led_umol_J / t.led_umol_J
    end = dict(food=food_W(case, year), basic=t.basic_kW_resident * 1000.0, gardens=gardens, industry=b['loads_W']['industry'])
    end_use = sum(end.values())
    bus = end_use * (1 + b['cooling_share'])
    losses = (1 - t.water_recovery) / (1 - t0.water_recovery)
    water = b['landscape_water_kg_year'] + (b['water_makeup_kg_year'] - b['landscape_water_kg_year']) * losses
    return dict(end_use_W=end, mean_W=bus, firm_W=bus * b['firm_over_bus'], water_kg_year=water,
                air_kg_year=b['air_makeup_kg_year'], earth_supplies_kg_year=b['earth_supplies_kg_year'],
                local_goods_kg_year=b['local_goods_kg_year'], operators=b['operators'] / t.output_per_job)
