"""What the settlement sells in a year: the coupled study's scenarios of demand judged at their own date
(research/studies/vacuum_economy/scenarios.json, read through its demand.py), priced from this model's own costs.

Each line returns its receipts, what it pays to Earth (the shares partners and services on Earth keep, hardware bought
there, passages), its cost for deciding whether to serve it (what it pays to Earth with its hands at the wage and its
local work at the lines' costs), the hands it takes, and what it asks of the settlement's industry and transport.
The seller's price follows the coupled study's rule: where buyers answer a price cut with more than a matching rise in
quantity (elasticity e above one) the price is the cost of a unit times e/(e-1); otherwise the ceiling, or for visits
the price the travel model is anchored at.
"""
from __future__ import annotations
from functools import lru_cache
import math

from research.studies.vacuum_economy import demand
from . import technology as tech
from .designs import inputs, v

LINES = ('culture', 'institutions', 'computing', 'propellant', 'materials', 'metals', 'visits')


@lru_cache(maxsize=None)
def scenarios():
    return demand.load()


@lru_cache(maxsize=None)
def world(scenario):
    s = scenarios()
    to2041 = v(inputs()['money']['USD2026_to_USD2041'])
    return demand.World(s, s['scenarios'][scenario], 2041, 8, to2041)


def seller(cost, elasticity, ceiling):
    if elasticity > 1:
        price = cost * elasticity / (elasticity - 1)
        return price if ceiling is None else min(price, ceiling)
    return ceiling


def year(scenario, tcase, y, prices, depot_water_kg_cost, stay_cost_day, feed_kg, metals_keep_USD=0.0):
    """The lines of one programme year. `prices` is the year's engineering.Prices; `depot_water_kg_cost` the cost of a
    kilogram of water at the depot; `stay_cost_day` what a visitor's day on the Moon costs beside the hotel's staff;
    `feed_kg` the feed mined at the sources in the year, which bounds the metals recovered; `metals_keep_USD` the value
    of metals a reserve of U would keep from the refined output instead of selling it."""
    s = scenarios()
    sc, lines = s['scenarios'][scenario], s['lines']
    w = world(scenario)
    date = 2041 + y
    age = y - 8
    p, t = prices, prices.t
    wage = p.wage
    out = {}

    # Sport and entertainment: receipts grow from the opening portfolio toward the mature following.
    cu = sc['culture']
    per_follower = w.dollars(demand.value(cu['receipts_per_follower_USD2024']), 2024) * (w.income(date) / w.income(2024)) ** demand.value(lines['culture']['income_elasticity'])
    mature = demand.value(cu['followers_share_of_world']) * w.population(date) * per_follower
    grown = demand.value(lines['culture']['opening_receipts_USD2041_year']) * (1 + demand.value(cu['growth_of_receipts'])) ** age
    receipts = min(mature, grown)
    jobs = (lines['culture']['resident_jobs']['fixed'] + lines['culture']['resident_jobs']['per_billion_USD2041'] * receipts / 1e9) / t.output_per_job
    power = lines['culture']['power_W']['fixed'] + lines['culture']['power_W']['per_billion_USD2041'] * receipts / 1e9
    earth = receipts * demand.value(cu['share_others_keep'])
    out['culture'] = dict(receipts=receipts, to_Earth=earth, jobs=jobs, power_W=power, cost=earth + jobs * wage + power * p.kWh * 8.766)

    # Hosted research and professional work: a share of world output, staffed by residents.
    ins = sc['institutions']
    receipts = demand.value(ins['share_of_world_output']) * w.output(date)
    kept = 1 - demand.value(ins['share_others_keep'])
    jobs = receipts * kept / (demand.value(ins['retained_USD2041_per_job']) * t.output_per_job)
    earth = receipts * demand.value(ins['share_others_keep'])
    out['institutions'] = dict(receipts=receipts, to_Earth=earth, jobs=jobs, power_W=0.0, cost=earth + jobs * wage)

    # Computing at the price of the same service on Earth. Its hardware is bought from Earth until chips are made on
    # the Moon, then made locally; its power is the year's firm supply.
    co, cl = sc['computing'], inputs()['sales']
    watts = demand.value(co['capacity_W_2049']) * w.output(date) / w.output(w.opening)
    price = w.dollars(demand.value(co['price_USD2024_W_year']), 2024)
    hw_kg = watts * v(cl['computing_hardware_kg_W'])
    life = v(cl['computing_hardware_life_years'])
    hw_price = (p.t.equipment_price_multiple * v(inputs()['money']['computing_hardware_USD2026_kg']) * v(inputs()['money']['USD2026_to_USD2041']) + p.landed)
    local_hw = t.arrived('chips_on_moon')
    hw_year = hw_kg / life * (p.line_kg['manufactured'] if local_hw else hw_price)
    load_W = watts * v(cl['computing_overhead'])
    earth = watts * price * demand.value(lines['computing']['Earth_service_share_of_receipts']) + (0.0 if local_hw else hw_year)
    cost = earth + (hw_year if local_hw else 0.0) + load_W / 1000.0 * p.supply_here['USD_kW_year']
    out['computing'] = dict(receipts=watts * price, to_Earth=earth, jobs=0.0, power_W=load_W, cost=cost,
                            hardware_kg_year=hw_kg / life, hardware_local=local_hw, watts=watts)

    # Propellant and material at the depot, against the same kilogram launched from Earth to the depot.
    sp = sc['space']
    ceiling = demand.value(lines['space']['share_of_Earth_delivered_price']) * p.launch * (1 + 0.6 * (t.landed_multiple - 1))
    anchor = w.dollars(demand.value(lines['space']['anchor_price_USD2026_kg']), 2026)
    activity = w.space_activity(date)
    e = demand.value(sp['price_elasticity'])
    prop_cost = depot_water_kg_cost + p.line_kg['propellant']
    lift = p.lift()
    mat_cost = p.line_kg['metal'] + lift['propellant_kg'] * prop_cost + lift['kWh'] * p.kWh + lift['launcher_kg_per_kg_year'] * p.line_kg['manufactured'] * 0.1
    bound = inputs()['sales']['depot_sales_cap_multiple'][tcase]
    for line, mass, unit_cost in (('propellant', 1.0, prop_cost), ('materials', demand.value(lines['space']['material_per_propellant']), mat_cost)):
        price = seller(unit_cost, e, ceiling)
        base = demand.value(sp['propellant_kg_year_2049']) * mass * activity
        wanted = base * (price / anchor) ** -e
        kg = wanted / (1 + wanted / (bound * base))
        out[line] = dict(receipts=kg * price, to_Earth=0.0, jobs=0.0, power_W=0.0, cost=kg * unit_cost, kg=kg, price=price)
    out['materials']['lift'] = lift['how']

    # Precious metals from the feed already mined; the platinum-group price falls as the supply to Earth grows. A share
    # kept for a reserve of U reaches no buyer on Earth, so the rest sells at a higher price.
    me, ml = sc['metals'], lines['metals']
    world_supply = ml['PGM_world_production_kg_year'] * w.output(date) / w.output(2025)
    wanted = demand.value(me['PGM_supply_share_of_world']) * world_supply
    refined = min(wanted, feed_kg * demand.value(me['refined_PGM_per_kg_feed']))
    gold = refined * demand.value(ml['gold_per_PGM'])
    pgm_price = lambda sold: w.dollars(ml['PGM_price_USD2025_kg'], 2025) * (1 + sold / world_supply) ** (-1 / demand.value(ml['price_elasticity']))
    gold_price = w.dollars(ml['gold_price_USD2025_kg'], 2025)
    kept = 0.0
    if metals_keep_USD > 0 and refined > 0:
        kept = min(1.0, metals_keep_USD / (refined * pgm_price(refined) + gold * gold_price))
    metal_price = pgm_price(refined * (1 - kept))
    receipts = refined * (1 - kept) * metal_price + gold * (1 - kept) * gold_price
    retained = (refined * metal_price + gold * gold_price) * kept
    treated = refined / demand.value(me['refined_PGM_per_kg_feed']) if refined else 0.0
    treatment = treated * demand.value(me['treatment_USD2041_per_kg_feed'])
    out['metals'] = dict(receipts=receipts, to_Earth=0.0, jobs=0.0, power_W=0.0, cost=treatment, kg=refined + gold, treated_feed_kg=treated,
                         retained_value=retained)

    # Visits: the passage both ways and a stay; the hotel's staff are residents' hands.
    vi, vl = sc['visits'], lines['visits']
    nights = demand.value(vl['stay_nights'])
    staff_per_bed = v(inputs()['sales']['hotel_staff_per_bed'])
    occupancy = demand.value(vl['occupancy'])
    passage = 2 * p.fare()
    stay = nights * (stay_cost_day + staff_per_bed * wage / (365.25 * occupancy))
    unit = passage + stay
    ev = demand.value(vi['price_elasticity'])
    anchor_visit = w.dollars(vl['anchor']['price_USD2022'], 2022)
    keep = demand.value(vl['Earth_service_share_of_receipts'])
    price = unit / (1 - keep) * ev / (ev - 1) if ev > 1 else anchor_visit
    wanted = vl['anchor']['visitors_year'] * w.wealthy_buyers(date, demand.value(vi['income_elasticity'])) * (price / anchor_visit) ** -ev
    cap = inputs()['sales']['visitor_cap_share'][tcase] * w.population(date)
    visitors = wanted / (1 + wanted / cap)
    beds = visitors * nights / (365.25 * occupancy)
    out['visits'] = dict(receipts=visitors * price, to_Earth=visitors * (passage + price * keep), jobs=beds * staff_per_bed, power_W=0.0,
                         cost=visitors * (unit + price * keep), visitors=visitors, price=price, unit_cost=unit,
                         launched_kg=visitors * 2 * v(inputs()['people']['launched_kg_per_passenger']))
    for line in out.values():
        worth = line['receipts'] + line.get('retained_value', 0.0)
        line['served'] = worth > line['cost'] and worth > 0
    return out
