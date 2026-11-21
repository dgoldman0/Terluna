"""What the customers of the vacuum buildout buy in each year, judged at that year's date.

A scenario of demand states, for seven lines of sale, how many customers there are in a year, what they pay and how
the quantity answers the price. The world the customers live in moves with the calendar: its population follows the
United Nations' projection, its income per person grows, and activity in space grows with the space economy. No
quantity is a multiple of the founding system. The lines are

  culture       sport and entertainment for followers on Earth
  institutions  research, operations and professional work that outside institutions pay the settlement to host
  computing     computing capacity sold at the cost of the same service on Earth
  propellant    propellant handed to customers at the depot
  materials     structural material handed to customers at the depot
  metals        platinum-group metals and gold recovered from the feed of the sources, returned to Earth
  visits        visitors who pay for the passage and a stay

(propellant and materials share one description, `space`, in scenarios.json). `book` returns the orders year by year
in the form the industry model and the financing model read. Its callers may state the price of the lines whose
price the system sets, and may scale any line down in any year, as when the residents cannot staff it.
"""
from __future__ import annotations
from bisect import bisect_right
import json
import math
from pathlib import Path

from provisioning.vacuum_industry import model
from shared.constants import JULIAN_YEAR_DAYS
from . import programme as definition

HERE = Path(__file__).resolve().parent
LINES = ('culture', 'institutions', 'computing', 'propellant', 'materials', 'metals', 'visits')
# The lines whose price the system sets from its own cost, and the unit each is priced in.
PRICED = {'propellant': 'kg', 'materials': 'kg', 'visits': 'visit'}
SECTORS = ('materials', 'source_products', 'propellant', 'mission_support', 'computing', 'destination', 'sports',
           'entertainment', 'licenses', 'research', 'professional')


def load(path=HERE / 'scenarios.json'):
    return json.loads(Path(path).read_text())


def value(entry):
    return entry['value'] if isinstance(entry, dict) else entry


def interpolate(points, x):
    """A straight line between the listed points, held at the first and the last value outside them."""
    xs = [p[0] for p in points]
    if x <= xs[0]:
        return points[0][1]
    if x >= xs[-1]:
        return points[-1][1]
    i = bisect_right(xs, x)
    (x0, y0), (x1, y1) = points[i - 1], points[i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


class World:
    """The customers' world in a calendar year, under one scenario."""

    def __init__(self, inputs, scenario, year_zero, first_productive_year, to_2041):
        self.inputs, self.scenario = inputs, scenario
        self.year_zero, self.opening = year_zero, year_zero + first_productive_year
        self.to_2041 = to_2041                                  # dollars of 2026 to constant dollars of 2041
        self.growth = value(scenario['income_per_person_growth_after_2049'])

    def dollars(self, amount, year):
        """An amount in the dollars of a stated year, in constant dollars of 2041."""
        factor = 1.0 if year == 2026 else self.inputs['money']['USD%d_to_USD2026' % year]
        return amount * factor * self.to_2041

    def population(self, calendar):
        return interpolate(self.inputs['world']['population']['by_calendar_year'], calendar)

    def income(self, calendar):
        """Income per person, as a multiple of 2024's."""
        before = value(self.inputs['world']['income_per_person_growth_to_2049'])
        return (1 + before) ** (min(calendar, 2049) - 2024) * (1 + self.growth) ** max(0, calendar - 2049)

    def output(self, calendar):
        """World output in constant dollars of 2041."""
        stated = self.inputs['world']['output']
        at = stated['calendar_year']
        scale = self.population(calendar) * self.income(calendar) / (self.population(at) * self.income(at))
        return self.dollars(stated['USD2015'], 2015) * scale

    def space_activity(self, calendar):
        """Activity in space, as a multiple of the opening year's: its own growth to a stated year, world output after."""
        rate = value(self.scenario['space']['activity_growth'])
        until = self.inputs['lines']['space']['activity_growth_until']
        own = (1 + rate) ** (min(calendar, until) - self.opening)
        return own * (self.output(calendar) / self.output(until) if calendar > until else 1.0)

    def wealthy_buyers(self, calendar, income_elasticity):
        """People able to buy a visit, as a multiple of the travel model's year."""
        lines = self.inputs['lines']['visits']
        start = lines['anchor']['calendar_year']
        rate = value(self.inputs['world']['wealthy_buyers_growth_to_2049'])
        to_opening = (1 + rate) ** (self.opening - start)
        people = self.population(calendar) / self.population(self.opening)
        richer = (self.income(calendar) / self.income(self.opening)) ** income_elasticity
        return to_opening * people * richer


def earth_delivered_price(coefficients, year, to_2041):
    """What a kilogram launched from Earth costs at the depot in a year, in constant dollars of 2041.

    The kilogram and the propellant that carries it from low Earth orbit are both launched, as in the industry
    model's arrangement `earth`, at that year's launch price; the propellant is bought on Earth.
    """
    transport = coefficients['industrial']['transport']
    data = dict(coefficients, industrial=dict(coefficients['industrial'], transport=dict(transport, earth_leg='earth')))
    route = model.routes(data, year)['Earth']
    burned = model.ferry(1e6, 0.0, route)['burned'] / 1e6
    launch = model.at(transport['LEO_USD2026_kg'], year)
    return ((1 + burned) * launch + burned * transport['earth_propellant_USD2026_kg']) * to_2041


PRICING = ('seller', 'cost_plus')


def seller_price(inputs, name, line, cost, to_2041, pricing='seller'):
    """The system's price for a line whose price it sets, given the line's cost per unit.

    Under `seller` it is the price that brings the seller most. Where a price cut brings more than a matching rise
    in quantity (an elasticity above one), that price is the cost times e/(e-1). Otherwise no price below the ceiling
    pays: propellant and materials are sold at the ceiling the demand states (returned as None), and visits at the
    price at which the travel model is anchored. Under `cost_plus` it is the cost with the stated margin, which
    carries more travellers and sells more propellant for less profit.
    """
    scenario, lines = inputs['scenarios'][name], inputs['lines']
    if pricing not in PRICING:
        raise ValueError('Unknown rule of pricing: ' + pricing)
    if pricing == 'cost_plus':
        return cost * (1 + value(lines['cost_plus_margin']))
    elasticity = value(scenario['visits' if line == 'visits' else 'space']['price_elasticity'])
    if elasticity > 1:
        return cost * elasticity / (elasticity - 1)
    if line != 'visits':
        return None
    return lines['visits']['anchor']['price_USD2022'] * inputs['money']['USD2022_to_USD2026'] * to_2041


def book(programme, ramp, inputs, name, coefficients, prices=None, scale=None, feed=None):
    """The orders of every year of the programme for one scenario on one ramp.

    prices: for a line in PRICED, the system's own price by year in dollars of 2041; where none is given the line
    is sold at its ceiling (propellant and materials) or its first price (visits).
    scale: for any line, the share of its demand that is served, by year (default all of it).
    feed: the feed mined at the sources by year, in kilograms, which limits the metals recovered from it; where
    none is given no metals are sold.
    """
    if name not in inputs['scenarios']:
        raise ValueError('Unknown scenario of demand: ' + name)
    scenario, lines = inputs['scenarios'][name], inputs['lines']
    calendar = programme['calendar']
    first, last = calendar['first_productive_year'], calendar['last_year']
    to_2041 = coefficients['normalization']['USD2026_to_USD2041']
    world = World(inputs, scenario, calendar['year_zero'], first, to_2041)
    prices, scale = prices or {}, scale or {}
    for line in list(prices) + list(scale):
        if line not in LINES or (line in prices and line not in PRICED):
            raise ValueError('No such line, or a line whose price the system does not set: ' + line)
    wage = value(lines['service_wage_USD2041_year'])
    culture, institutions, computing = scenario['culture'], scenario['institutions'], scenario['computing']
    space, metals, visits = scenario['space'], scenario['metals'], scenario['visits']
    rows = []

    def served(line, year):
        share = scale.get(line, {}).get(year, 1.0)
        if not 0.0 <= share <= 1.0:
            raise ValueError('A line is served in a share between none and all of its demand')
        return share

    for year in range(first, last + 1):
        date = calendar['year_zero'] + year
        age = year - first
        residents = definition.population(programme, ramp, year)
        # The work a job does grows with labour productivity, and the wage rises in step.
        output_per_job = model.labour_index(coefficients, year)
        detail = {}

        # Sport and entertainment: receipts grow from the opening portfolio toward the mature following.
        per_follower = (world.dollars(value(culture['receipts_per_follower_USD2024']), 2024)
                        * (world.income(date) / world.income(2024)) ** value(lines['culture']['income_elasticity']))
        mature = value(culture['followers_share_of_world']) * world.population(date) * per_follower
        grown = value(lines['culture']['opening_receipts_USD2041_year']) * (1 + value(culture['growth_of_receipts'])) ** age
        cultural = min(mature, grown) * served('culture', year)
        cultural_jobs = (lines['culture']['resident_jobs']['fixed'] + lines['culture']['resident_jobs']['per_billion_USD2041'] * cultural / 1e9) * (1.0 if cultural else 0.0) / output_per_job
        cultural_W = (lines['culture']['power_W']['fixed'] + lines['culture']['power_W']['per_billion_USD2041'] * cultural / 1e9) * (1.0 if cultural else 0.0)
        detail['culture'] = dict(receipts_USD2041=cultural, quantity=cultural / per_follower, unit='followers', price_USD2041=per_follower)

        # Hosted research and professional work: a share of world output, staffed by residents.
        hosted = value(institutions['share_of_world_output']) * world.output(date) * served('institutions', year)
        kept = 1 - value(institutions['share_others_keep'])
        retained = value(institutions['retained_USD2041_per_job']) * output_per_job
        hosted_jobs = hosted * kept / retained
        detail['institutions'] = dict(receipts_USD2041=hosted, quantity=hosted_jobs, unit='resident jobs', price_USD2041=retained / kept)

        # Computing at the cost of the same service on Earth.
        watts = value(computing['capacity_W_2049']) * world.output(date) / world.output(world.opening) * served('computing', year)
        watt_price = world.dollars(value(computing['price_USD2024_W_year']), 2024)
        detail['computing'] = dict(receipts_USD2041=watts * watt_price, quantity=watts, unit='W', price_USD2041=watt_price)

        # Propellant and structural material at the depot, against the same kilogram launched from Earth.
        ceiling = value(lines['space']['share_of_Earth_delivered_price']) * earth_delivered_price(coefficients, year, to_2041)
        anchor = world.dollars(value(lines['space']['anchor_price_USD2026_kg']), 2026)
        activity = world.space_activity(date)
        for line, mass in (('propellant', 1.0), ('materials', value(lines['space']['material_per_propellant']))):
            own = prices.get(line, {}).get(year)
            price = ceiling if own is None else min(ceiling, own)
            kg = (value(space['propellant_kg_year_2049']) * mass * activity * (price / anchor) ** -value(space['price_elasticity'])
                  * served(line, year))
            detail[line] = dict(receipts_USD2041=kg * price, quantity=kg, unit='kg', price_USD2041=price, ceiling_USD2041=ceiling)

        # Precious metals from the feed already mined. The platinum-group price falls as the supply to Earth grows.
        world_supply = lines['metals']['PGM_world_production_kg_year'] * world.output(date) / world.output(2025)
        wanted = value(metals['PGM_supply_share_of_world']) * world_supply
        from_feed = (feed or {}).get(year, 0.0) * value(metals['refined_PGM_per_kg_feed'])
        refined = min(wanted, from_feed) * served('metals', year)
        gold = refined * value(lines['metals']['gold_per_PGM'])
        metal_price = (world.dollars(lines['metals']['PGM_price_USD2025_kg'], 2025)
                       * (1 + refined / world_supply) ** (-1 / value(lines['metals']['price_elasticity'])))
        metal_receipts = refined * metal_price + gold * world.dollars(lines['metals']['gold_price_USD2025_kg'], 2025)
        treated = refined / value(metals['refined_PGM_per_kg_feed'])
        treatment = treated * value(metals['treatment_USD2041_per_kg_feed'])
        detail['metals'] = dict(receipts_USD2041=metal_receipts, quantity=refined + gold, unit='kg', price_USD2041=metal_receipts / (refined + gold) if refined else metal_price,
                                treated_feed_kg=treated)

        # Visits: the price covers the passage and the stay.
        first_price = value(lines['visits']['first_price_USD2041'])
        visit_price = prices.get('visits', {}).get(year, first_price)
        anchor_visit = world.dollars(lines['visits']['anchor']['price_USD2022'], 2022)
        visitors = (lines['visits']['anchor']['visitors_year'] * world.wealthy_buyers(date, value(visits['income_elasticity']))
                    * (visit_price / anchor_visit) ** -value(visits['price_elasticity']) * served('visits', year))
        nights = visitors * value(lines['visits']['stay_nights'])
        beds = math.ceil(nights / (JULIAN_YEAR_DAYS * value(lines['visits']['occupancy']))) if visitors else 0
        detail['visits'] = dict(receipts_USD2041=visitors * visit_price, quantity=visitors, unit='visitors', price_USD2041=visit_price)

        jobs = cultural_jobs + hosted_jobs
        payroll = jobs * wage * output_per_job
        outside = (cultural * value(culture['share_others_keep']) + hosted * value(institutions['share_others_keep'])
                   + detail['computing']['receipts_USD2041'] * value(lines['computing']['Earth_service_share_of_receipts'])
                   + detail['visits']['receipts_USD2041'] * value(lines['visits']['Earth_service_share_of_receipts']) + treatment)
        sport = value(lines['culture']['sport_share'])
        receipts = dict.fromkeys(SECTORS, 0.0)
        receipts.update(sports=cultural * sport, entertainment=cultural * (1 - sport), research=hosted,
                        computing=detail['computing']['receipts_USD2041'], propellant=detail['propellant']['receipts_USD2041'],
                        materials=detail['materials']['receipts_USD2041'], source_products=metal_receipts,
                        destination=detail['visits']['receipts_USD2041'])
        own_cost = dict(culture=cultural * value(culture['share_others_keep']) + cultural_jobs * wage * output_per_job,
                        institutions=hosted * value(institutions['share_others_keep']) + hosted_jobs * wage * output_per_job,
                        computing=detail['computing']['receipts_USD2041'] * value(lines['computing']['Earth_service_share_of_receipts']),
                        propellant=0.0, materials=0.0, metals=treatment,
                        visits=detail['visits']['receipts_USD2041'] * value(lines['visits']['Earth_service_share_of_receipts']))
        for line in LINES:
            detail[line]['own_cost_USD2041'] = own_cost[line]
        rows.append(dict(
            year=year, population=residents,
            commissioned_people=0.0 if year == first else residents - definition.population(programme, ramp, year - 1),
            new_residential_people=definition.population(programme, ramp, year + 2) - definition.population(programme, ramp, year + 1),
            goods=0.0, source_goods=refined + gold, space_goods=0.0, space_source_goods=detail['materials']['quantity'],
            propellant_sales=detail['propellant']['quantity'], compute=watts,
            visitor_beds=beds, visitor_nights_year=nights, visitors_year=visitors,
            passenger_launch_kg_year=visitors * value(lines['visits']['launched_kg_per_visitor']),
            passenger_cabin_kg_year=visitors * value(lines['visits']['cabin_kg_per_visitor']),
            service_jobs=jobs, premium_jobs=0.0, logistics_jobs=0.0, service_kit_count=0,
            service_power_W=cultural_W, premium_power_W=0.0,
            service_installed_capital_USD2041=0.0, service_new_capital_USD2041=0.0, service_replacement_allowance_USD2041=0.0,
            commercial_operating_USD2041=outside + payroll, commercial_payroll_USD2041=payroll,
            gross_USD2041=math.fsum(receipts.values()), revenue_by_sector_USD2041=receipts, lines=detail))
    return dict(name='%s:%s' % (name, ramp), scenario=name, ramp=ramp, portfolio=name, market=name, market_growth=world.growth,
                market_growth_name=name, annual=rows)
