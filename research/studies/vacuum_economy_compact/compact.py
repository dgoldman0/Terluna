"""A compact model of the vacuum buildout's economy: what a resident needs and what the settlement then needs, year by
year, with the technology of each year.

usage: python3 -m research.studies.vacuum_economy_compact.compact [--ramp NAME] [--technology CASE]   (repository root)

The settlement's needs per resident start from the founding design (research/studies/lava_tube_cities, its saved
build-out for 250,000 residents in four lava-tube cities and a hundred crater towns) and move with the technology paths
of technology.json: food and basic services set the power, water recovery the water drawn, built mass the construction.
The supply of firm power is chosen each year as the cheapest of the options the technology offers, and transport from
Earth follows the launch price. The population follows a ramp of the programme definition. Every quantity is stated
per resident and for the whole, beside the checks that say whether the whole is sensible.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path

from research.studies.vacuum_economy import programme as definition

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / 'research/studies/lava_tube_cities/results/buildout.json'
YEAR_ZERO = 2041
HOURS = 8766.0
DAY = 86400.0
KCAL = 4184.0
CASES = ('static', 'slow', 'central', 'fast')


@lru_cache(maxsize=None)
def technology():
    return json.loads((HERE / 'technology.json').read_text())


def path(name, case, year):
    """The value of one technology quantity in a calendar year: the 2041 value in the static case, interpolated
    between the stated years otherwise."""
    t = technology()
    q = t['quantities'][name]
    years = t['years']
    values = q['slow' if case == 'static' else case]
    if case == 'static' or year <= years[0]:
        return values[0]
    if year >= years[-1]:
        return values[-1]
    i = next(k for k in range(1, len(years)) if year <= years[k])
    f = (year - years[i - 1]) / (years[i] - years[i - 1])
    a, b = values[i - 1], values[i]
    if q.get('interpolate') == 'geometric' and a > 0 and b > 0:
        return a * (b / a) ** f
    return a + (b - a) * f


def fusion(case, year):
    """Fusion's specific mass, price and life once it has arrived on the Moon in this case, None before."""
    q = technology()['quantities']['fusion']
    arrival = None if case == 'static' else q['arrives'][case]
    if arrival is None or year < arrival:
        return None
    return dict(kg_kWe=q['kg_kWe'][case], USD2024_kWe=q['USD2024_kWe'][case], life_years=q['life_years'][case])


@lru_cache(maxsize=None)
def design():
    """The founding design per resident: its reference case for 250,000 residents."""
    d = json.loads(DESIGN.read_text())['cases']['reference']
    n = d['population']
    u, r, g, f = d['utilities'], d['resources'], d['geometry'], d['food']
    loads = u['mean_load_components_W']
    return dict(
        population=n,
        food_W=loads['food'] / n, basic_W=loads['basic_services'] / n, gardens_W=loads['civic_gardens'] / n,
        industry_W=loads['everyday_industry'] / n, cooling_W=loads['cooling'] / n,
        bus_over_end_use=u['mean_electric_bus_W'] / u['mean_end_use_W'],
        firm_over_bus=u['occupied_firm_W'] / u['mean_electric_bus_W'],
        local_built_kg=r['local_finished_kg'] / n, imported_built_kg=r['special_imports_kg'] / n,
        initial_air_water_kg=r['initial_air_water_kg'] / n, construction_kWh=r['construction_electricity_kWh'] / n,
        water_makeup_kg_year=r['operating_import_kg_year']['water'] / n,
        other_operating_imports_kg_year=sum(v for k, v in r['operating_import_kg_year'].items() if k != 'water') / n,
        renewal_imports_kg_year=sum(r['renewal_import_kg_year'].values()) / n,
        floor_m2=g['required_floor_m2'] / n, farm_floor_m2=f['farm_floor_m2'] / n,
        operators=(r['food_operators'] + r['industry_operators']) / n,
        diet_kcal_day=2500.0)


def food_W(case, year):
    """Food's average electrical load per resident: crops under lamps for the share of the diet they supply, and
    electro-microbial food for the rest, scaled from the founding design's food load."""
    y = year
    base = design()
    # The founding design's crop electricity per kcal eaten, with its losses and auxiliaries, at 3.5 umol/J and 1.25 kcal/mol.
    share0 = path('microbial_diet_share', 'static', y)
    crop_J_kcal0 = base['food_W'] * DAY / (base['diet_kcal_day'] * (1 - share0))
    share = path('microbial_diet_share', case, y)
    crop_J_kcal = crop_J_kcal0 * (3.5 / path('led_umol_J', case, y)) * (1.25 / path('crop_kcal_per_mol', case, y))
    microbial_J_kcal = KCAL / path('microbial_efficiency', case, y)
    kcal = base['diet_kcal_day']
    return (kcal * (1 - share) * crop_J_kcal + kcal * share * microbial_J_kcal) / DAY


def needs(case, year):
    """What one resident needs in a calendar year: average and firm power, built mass made on the Moon and bought from
    Earth, the water and other imports of each year, the floor and the jobs of the settlement's own operation."""
    b = design()
    food = food_W(case, year)
    basic = path('basic_kW_resident', case, year) * 1000
    end_use = food + basic + b['gardens_W'] + b['industry_W']
    # Cooling follows the heat of the end uses, as in the design, where the electric bus is end use plus cooling.
    cooling = b['cooling_W'] * end_use / (b['food_W'] + b['basic_W'] + b['gardens_W'] + b['industry_W'])
    bus = end_use + cooling
    firm = bus * b['firm_over_bus']
    built = path('built_mass_multiple', case, year)
    losses = (1 - path('water_recovery', case, year)) / (1 - path('water_recovery', 'static', year))
    # Farm floor follows the crop food grown and the food a square metre yields at the design's light.
    share0 = path('microbial_diet_share', 'static', year)
    farm = b['farm_floor_m2'] * (1 - path('microbial_diet_share', case, year)) / (1 - share0) * 1.25 / path('crop_kcal_per_mol', case, year)
    return dict(food_W=food, basic_W=basic, mean_W=bus, firm_W=firm,
                local_built_kg=b['local_built_kg'] * built, imported_built_kg=b['imported_built_kg'] * built,
                water_kg_year=b['water_makeup_kg_year'] * losses, other_imports_kg_year=b['other_operating_imports_kg_year'],
                renewal_imports_kg_year=b['renewal_imports_kg_year'] * path('equipment_import_share', case, year) / 0.05,
                floor_m2=b['floor_m2'] - b['farm_floor_m2'] + farm,
                operators=b['operators'] / path('output_per_job', case, year))


def supply(case, year, rate=0.03):
    """The cost a year of one watt of firm power on the Moon, for each option the technology offers in that year, in
    dollars of 2026: its capital at the price of the equipment and of landing what is bought from Earth, annuitized over
    its life at the rate, beside the mass bought from Earth per watt."""
    launch = path('launch_USD2026_kg', case, year) * path('landed_multiple', case, year)
    equipment = 1000.0 * path('equipment_price_multiple', case, year)       # USD2026 a kg made on Earth, the study's value
    local = 50.0 / path('output_per_job', case, year) + 5.0                  # USD2026 a kg made on the Moon, the study's value
    life = path('equipment_life_years', case, year)

    def annual(capital, years):
        return capital * rate / (1 - (1 + rate) ** -years)

    out = {}
    # Solar with storage through the 354-hour night: arrays carry the day's load and charge the night's.
    rt = path('storage_round_trip', case, year)
    array_W = 1 + 1 / rt
    storage_kg = 354.0 / path('storage_Wh_kg', case, year)
    array_kg = array_W / path('array_W_kg', case, year)
    share = path('array_local_share', case, year)
    imported = array_kg * (1 - share) + storage_kg * 0.2
    made = array_kg * share + storage_kg * 0.8
    out['solar_with_storage'] = dict(USD2026_W_year=annual(imported * (equipment + launch) + made * local, life), imported_kg_W=imported)
    # Beamed from stations in lunar orbit, built on Earth, with receivers made on the Moon. A delivered watt needs the
    # arrays' watts divided by the beam's efficiency and the end-of-life factor, times the share of the constellation that
    # is not over the settlement: the founding design's ten stations with two carrying the load, 5. Its arrays are the
    # same arrays as on the surface; at 2041 this is 11.5 W of array, 0.12 kg, per delivered watt.
    beam_kg = path('array_W_kg', case, year) ** -1 * 5.0 / (path('beam_efficiency', case, year) * 0.85)
    out['beamed'] = dict(USD2026_W_year=annual(beam_kg * (equipment + launch) + 0.02 * local, 20), imported_kg_W=beam_kg)
    # Fission plant built on Earth at its own price, landed whole.
    fission_kg = path('fission_kg_kWe', case, year) / 1000
    price = path('fission_USD2024_kWe', case, year) / 1000
    out['fission'] = dict(USD2026_W_year=annual(price + fission_kg * launch, path('fission_life_years', case, year)), imported_kg_W=fission_kg)
    f = fusion(case, year)
    if f:
        out['fusion'] = dict(USD2026_W_year=annual(f['USD2024_kWe'] / 1000 + f['kg_kWe'] / 1000 * launch, f['life_years']), imported_kg_W=f['kg_kWe'] / 1000)
    best = min(out, key=lambda k: out[k]['USD2026_W_year'])
    return dict(options=out, chosen=best, launch_landed_USD2026_kg=launch)


# The share of the world's people who could visit in a year, the study's own bound (no source gives one): slow, central
# and fast technology carry 0.01%, 0.1% and 1%; the static case the slow bound.
VISITOR_CAP_SHARE = dict(static=1e-4, slow=1e-4, central=1e-3, fast=1e-2)
LAUNCH_PAYLOAD_KG = 100000.0               # a fully reusable heavy launcher's payload [E4]
LAUNCHES_PER_PAD_YEAR = 76.0               # the largest authorised rate found for a heavy-lift complex [T5]
USD2026_TO_2041 = 1.557967416600765        # the staged studies' normalization (provisioning/vacuum_industry)
# Unit values the staged model carried without a source, USD2026: the study's starting values, moving with the
# technology paths. Local bulk and metal follow the price of equipment as industry learns; imports land at the launch
# price times the landing multiple.
STAGED = dict(bulk_kg=1.0, metal_kg=10.0, imported_equipment_kg=500.0, special_kg=150.0, water_at_depot_kg=5.0,
              wage_year=200000.0)
METAL_SHARE_OF_LOCAL = 0.125                # the founding demand's 27.4 t of metal in 219 t of local stock a resident
LOCAL_LIFE_YEARS = 80.0                     # the design's ordinary durable stock


def unit_values(case, year):
    """Unit values of a calendar year in USD2026: a kilogram made on the Moon, bought from Earth and landed, and a watt of
    firm power for a year by the cheapest supply."""
    s = supply(case, year)
    price = path('equipment_price_multiple', case, year)
    landed = s['launch_landed_USD2026_kg']
    return dict(bulk_kg=STAGED['bulk_kg'] * price, metal_kg=STAGED['metal_kg'] * price,
                imported_kg=STAGED['imported_equipment_kg'] * price + landed, special_kg=STAGED['special_kg'] + landed,
                water_kg=STAGED['water_at_depot_kg'] + landed * 0.1, landed_kg=landed,
                firm_W_year=s['options'][s['chosen']]['USD2026_W_year'], supply=s['chosen'],
                wage_year=STAGED['wage_year'] * path('output_per_job', case, year))


def resident_costs(case, year, rate=0.03):
    """What a new resident's share of the settlement costs to build, and a resident's share to run for a year, in
    USD2041: local and imported construction, firm power, the first air and water; then renewal, power, water and other
    imports, and the settlement's own operators at the wage of the year."""
    n, u = needs(case, year), unit_values(case, year)
    share = path('equipment_import_share', case, year) / 0.05
    imported_built = n['imported_built_kg'] * share
    local = n['local_built_kg']
    local_value = local * (METAL_SHARE_OF_LOCAL * u['metal_kg'] + (1 - METAL_SHARE_OF_LOCAL) * u['bulk_kg'])
    power_capital = u['firm_W_year'] * (1 - (1 + rate) ** -20) / rate * n['firm_W']
    build = dict(local=local_value, imported=imported_built * u['imported_kg'], power=power_capital,
                 air_water=design()['initial_air_water_kg'] * u['water_kg'])
    run_year = dict(renewal=local_value / LOCAL_LIFE_YEARS + n['renewal_imports_kg_year'] * u['imported_kg'],
                    power=n['firm_W'] * u['firm_W_year'], water=n['water_kg_year'] * u['water_kg'],
                    other_imports=n['other_imports_kg_year'] * u['special_kg'], operators=n['operators'] * u['wage_year'])
    return dict(build_USD2041={k: v * USD2026_TO_2041 for k, v in build.items()},
                run_USD2041_year={k: v * USD2026_TO_2041 for k, v in run_year.items()},
                imported_built_kg=imported_built, supply=u['supply'], landed_USD2026_kg=u['landed_kg'])


def sales(scenario, case, y, residents):
    """What the system sells in programme year y under a scenario of demand, priced from the compact model's costs:
    receipts and the costs each line itself carries, in USD2041 a year, by line."""
    from research.studies.vacuum_economy import demand
    v = demand.value
    inputs = demand.load()
    sc, lines = inputs['scenarios'][scenario], inputs['lines']
    world = demand.World(inputs, sc, YEAR_ZERO, 8, USD2026_TO_2041)
    date, age = YEAR_ZERO + y, y - 8
    u = unit_values(case, date)
    out = {}
    # Sport and entertainment, and hosted work: receipts the scenario states, with the shares others keep.
    culture = sc['culture']
    per_follower = world.dollars(v(culture['receipts_per_follower_USD2024']), 2024) * (world.income(date) / world.income(2024)) ** v(lines['culture']['income_elasticity'])
    cultural = min(v(culture['followers_share_of_world']) * world.population(date) * per_follower,
                   v(lines['culture']['opening_receipts_USD2041_year']) * (1 + v(culture['growth_of_receipts'])) ** age)
    out['culture'] = (cultural, cultural * v(culture['share_others_keep']))
    hosted = v(sc['institutions']['share_of_world_output']) * world.output(date)
    out['institutions'] = (hosted, hosted * v(sc['institutions']['share_others_keep']))
    # Computing at the price of the same service on Earth; the system pays its power, cooled, and the Earth's share.
    watts = v(sc['computing']['capacity_W_2049']) * world.output(date) / world.output(world.opening)
    price = world.dollars(v(sc['computing']['price_USD2024_W_year']), 2024)
    cost = watts * (u['firm_W_year'] * 1.3 * USD2026_TO_2041) + watts * price * v(lines['computing']['Earth_service_share_of_receipts'])
    out['computing'] = (watts * price, cost) if price * watts > cost else (0.0, 0.0)
    # Propellant and material at the depot: the price ceiling is the share of the same kilogram launched from Earth to
    # the depot; the system's cost is water at the depot and its electrolysis.
    ceiling = v(lines['space']['share_of_Earth_delivered_price']) * path('launch_USD2026_kg', case, date) * 2.8 * USD2026_TO_2041
    anchor = world.dollars(v(lines['space']['anchor_price_USD2026_kg']), 2026)
    e = v(sc['space']['price_elasticity'])
    kg_cost = (STAGED['water_at_depot_kg'] * 1.1 + 65 / 9 * u['firm_W_year'] / HOURS * 1000 / 1000) * USD2026_TO_2041
    sell = min(ceiling, kg_cost * e / (e - 1)) if e > 1 else ceiling
    activity = world.space_activity(date)
    for line, mass in (('propellant', 1.0), ('materials', v(lines['space']['material_per_propellant']))):
        kg = v(sc['space']['propellant_kg_year_2049']) * mass * activity * (sell / anchor) ** -e
        out[line] = (kg * sell, kg * kg_cost) if sell > kg_cost else (0.0, 0.0)
    # Visits: the passage both ways and a stay of the stated nights at a resident's daily cost and the hotel's staff.
    visits = lines['visits']
    passage = (v(visits['launched_kg_per_visitor']) + v(visits['cabin_kg_per_visitor'])) * u['landed_kg'] * 1.3
    stay = v(visits['stay_nights']) * (sum(resident_costs(case, date)['run_USD2041_year'].values()) / 365.25 / USD2026_TO_2041 + 250000 * 0.25 / 365.25 / v(visits['occupancy']) / USD2026_TO_2041)
    visit_cost = (passage + stay) * USD2026_TO_2041
    ev = v(sc['visits']['price_elasticity'])
    anchor_visit = world.dollars(visits['anchor']['price_USD2022'], 2022)
    visit_price = visit_cost * ev / (ev - 1) if ev > 1 else anchor_visit
    wanted = visits['anchor']['visitors_year'] * world.wealthy_buyers(date, v(sc['visits']['income_elasticity'])) * (visit_price / anchor_visit) ** -ev
    # The scenarios' demand answers price with one elasticity at every price; as a passage falls from millions of dollars
    # to tens of thousands it would carry more people than live on Earth. Visits are held to a share of the world's people
    # a year, the study's own bound, saturating smoothly.
    cap = VISITOR_CAP_SHARE[case] * world.population(date)
    visitors = wanted / (1 + wanted / cap)
    keep = v(visits['Earth_service_share_of_receipts'])
    out['visits'] = (visitors * visit_price, visitors * (visit_cost + visit_price * keep)) if visit_price > visit_cost else (0.0, 0.0)
    detail = {k: dict(receipts_USD2041=r, cost_USD2041=c) for k, (r, c) in out.items()}
    detail['visits'].update(visitors=visitors if out['visits'][0] else 0.0, price_USD2041=visit_price, cost_per_visit_USD2041=visit_cost)
    detail['visits']['launched_kg'] = detail['visits']['visitors'] * v(visits['launched_kg_per_visitor'])
    return detail


def account(ramp='constant_growth', scenario='central', case='central', rate=0.03):
    """The whole programme's money by year: the founding build before the opening, each year's build for new residents
    and running for all, the lines' receipts and own costs; with the funding needed, the net present value at the rate
    and the year the cumulative account turns positive."""
    programme = definition.load()
    rows, cumulative, worst, npv, payback = [], 0.0, 0.0, 0.0, None
    opening = definition.population(programme, ramp, 8)
    founding = opening * sum(resident_costs(case, YEAR_ZERO + 4)['build_USD2041'].values())
    for y in range(0, 221):
        if y < 8:
            build, running, receipts, line_costs = founding / 8, 0.0, 0.0, 0.0
            people, checks = 0.0, {}
        else:
            people = definition.population(programme, ramp, y)
            new = max(0.0, definition.population(programme, ramp, y + 1) - people)
            c = resident_costs(case, YEAR_ZERO + y)
            build, running = new * sum(c['build_USD2041'].values()), people * sum(c['run_USD2041_year'].values())
            s = sales(scenario, case, y, people)
            receipts = sum(x['receipts_USD2041'] for x in s.values())
            line_costs = sum(x['cost_USD2041'] for x in s.values())
            n = needs(case, YEAR_ZERO + y)
            imported = people * (n['other_imports_kg_year'] + n['renewal_imports_kg_year']) + new * c['imported_built_kg']
            launches = (imported + s['visits']['launched_kg']) / LAUNCH_PAYLOAD_KG
            checks = dict(visitors=s['visits']['visitors'], imported_t=imported / 1000, launches=launches, pads=launches / LAUNCHES_PER_PAD_YEAR,
                          mean_GW=people * n['mean_W'] / 1e9, receipts_by_line={k: x['receipts_USD2041'] for k, x in s.items()})
        net = receipts - line_costs - build - running
        cumulative += net
        worst = min(worst, cumulative)
        npv += net / (1 + rate) ** y
        if payback is None and y >= 8 and cumulative >= 0:
            payback = y
        rows.append(dict(year=y, residents=people, receipts=receipts, line_costs=line_costs, build=build, running=running, net=net, cumulative=cumulative, checks=checks))
    return dict(ramp=ramp, scenario=scenario, technology=case, funding_needed_USD2041=-worst, NPV_3pct_USD2041=npv,
                payback_year=payback, rows=rows)


def run(ramp='constant_growth', case='central', first=8, last=220):
    """Year by year: residents, their needs, the power supply chosen and the whole's power, imports and launches."""
    programme = definition.load()
    rows = []
    for y in range(first, last + 1):
        year = YEAR_ZERO + y
        people = definition.population(programme, ramp, y)
        new = max(0.0, definition.population(programme, ramp, y + 1) - people)
        n = needs(case, year)
        s = supply(case, year)
        chosen = s['options'][s['chosen']]
        imports = people * (n['other_imports_kg_year'] + n['renewal_imports_kg_year']) + new * n['imported_built_kg'] \
            + new * n['firm_W'] * chosen['imported_kg_W']
        rows.append(dict(year=y, calendar_year=year, residents=people, new_residents=new,
                         mean_W_resident=n['mean_W'], firm_W_resident=n['firm_W'], food_W_resident=n['food_W'],
                         power_GW_mean=people * n['mean_W'] / 1e9, supply=s['chosen'],
                         power_USD2026_W_year=chosen['USD2026_W_year'],
                         imports_t_year=imports / 1000, launches_year=imports / 100000.0,
                         water_t_year=people * n['water_kg_year'] / 1000,
                         operators_resident=n['operators'], local_built_t_year=new * n['local_built_kg'] / 1000))
    return rows


SCHEMA = 'terluna.research.vacuum-economy-compact/1'
HORIZONS = (100, 200, 220)
NEEDS_YEARS = (2049, 2070, 2100, 2150, 2241)


def digest(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dependencies():
    from research.studies.vacuum_economy import demand
    files = (HERE / 'compact.py', HERE / 'technology.json', DESIGN, ROOT / 'research/studies/vacuum_economy/scenarios.json',
             Path(demand.__file__), ROOT / 'research/studies/vacuum_economy/programme.json', Path(definition.__file__))
    return {str(Path(p).resolve().relative_to(ROOT)): digest(p) for p in files}


def product():
    """Every ramp, scenario of demand and technology case: the funding needed, the net present value and the year the
    account turns positive, and at each horizon the account and the checks of the whole; with what a resident needs and
    costs in chosen years under each technology case."""
    accounts = []
    for ramp in ('constant_growth', 'front_loaded'):
        for scenario in ('low', 'central', 'high'):
            for case in CASES:
                a = account(ramp, scenario, case)
                at = {str(h): {k: a['rows'][h][k] for k in ('residents', 'receipts', 'line_costs', 'build', 'running', 'net', 'cumulative', 'checks')}
                      for h in HORIZONS}
                accounts.append(dict(ramp=ramp, scenario=scenario, technology=case, funding_needed_USD2041=a['funding_needed_USD2041'],
                                     NPV_3pct_USD2041=a['NPV_3pct_USD2041'], payback_year=a['payback_year'], at_horizon=at))
    residents = {case: {str(year): dict(needs=needs(case, year), costs=resident_costs(case, year), supply=supply(case, year)['chosen'])
                        for year in NEEDS_YEARS} for case in CASES}
    return dict(schema=SCHEMA, producer=dict(path=str((HERE / 'compact.py').relative_to(ROOT)), sha256=digest(HERE / 'compact.py')),
                dependencies=dependencies(),
                evidence='A compact model on stated paths of technology. The settlement\'s needs per resident are the founding design\'s, '
                         'moving with the paths; the paths beyond their sourced starting values, the unit values of local work and of '
                         'imported equipment, the bound on visits and the transport multiples are the study\'s own. Results are '
                         'conditional on all of them and show how far the technology carries the outcome.',
                units='Years of the programme (zero is 2041), residents, kilograms, watts, constant USD2041.',
                reading_rule='Each account is one ramp of residents, one scenario of demand and one technology case. Funding needed is '
                             'the deepest cumulative shortfall from year zero; the net present value discounts every year\'s net at 3% '
                             'to year zero; the payback year is the first after the opening in which the cumulative account is not '
                             'negative. Checks give the visitors, the mass imported, the launches from Earth at 100 t and the pads at '
                             '76 launches a year, the settlement\'s mean power and the receipts by line in the horizon\'s year.',
                accounts=accounts, residents=residents)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--ramp', default='constant_growth')
    parser.add_argument('--technology', default='central', choices=CASES)
    parser.add_argument('--write', action='store_true', help='write results/compact.json')
    args = parser.parse_args()
    if args.write:
        out = HERE / 'results/compact.json'
        out.parent.mkdir(exist_ok=True)
        out.write_text(json.dumps(product(), indent=1) + '\n')
        print(out)
        return
    rows = run(args.ramp, args.technology)
    print('%5s %11s %9s %8s %8s %16s %10s %12s %10s %10s' % ('year', 'residents', 'food kW', 'mean kW', 'GW', 'supply', '$/W-yr', 'imports t/yr', 'launches', 'water kt'))
    for r in rows:
        if r['year'] in (8, 25, 50, 75, 100, 125, 150, 175, 200, 220):
            print('%5d %11.0f %9.1f %8.1f %8.2f %16s %10.3f %12.0f %10.0f %10.1f' % (
                r['calendar_year'], r['residents'], r['food_W_resident'] / 1e3, r['mean_W_resident'] / 1e3, r['power_GW_mean'],
                r['supply'], r['power_USD2026_W_year'], r['imports_t_year'], r['launches_year'], r['water_t_year'] / 1e3))


if __name__ == '__main__':
    main()
