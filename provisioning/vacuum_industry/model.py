"""The industry of the vacuum buildout: its stocks of equipment by cohort, its dated orders and its physical flows.

For a schedule of residents and customers' orders, year by year, the model sizes the equipment that houses the
residents and fills the orders, orders each component two years before it enters service, retires it by age, and
closes the flows of material, freight, propellant and power that the equipment and the orders require. It reads the
engineering coefficients it is given and takes its physical constants from shared/constants.json.

Three points in the calculation accept measures that change it: after the ledger of cohorts, after the physical
flows of a year, and after the equipment a year requires. A measures object with the methods after_ledger,
after_physical and after_requirements is called at each; recovery.py holds the measures of material recovery and
equipment refurbishment. Without measures the model is the plain calculation.
"""
from __future__ import annotations
from copy import deepcopy
from functools import lru_cache
import math
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY, STEFAN_BOLTZMANN
HOURS = JULIAN_DAY * JULIAN_YEAR_DAYS / 3600
# The share of the residents who can hold the industry's and the settlement's jobs.
WORKFORCE_SHARE = 0.4
DESTS = ('surface', 'source', 'depot', 'leo')

def at(value, year):
    """A coefficient's value in a year: the value itself, or read from a schedule {'by_year': [[year, value], ...]}.

    A number is interpolated in a straight line between the listed years and held outside them. A name takes the
    entry of the latest listed year that is not after the given one.
    """
    if not isinstance(value, dict):
        return value
    points = value['by_year']
    if year <= points[0][0]:
        return points[0][1]
    for (y0, v0), (y1, v1) in zip(points, points[1:]):
        if year < y1:
            return v0 if isinstance(v0, str) else v0 + (v1 - v0) * (year - y0) / (y1 - y0)
    return points[-1][1]


def burns(cargo, back, dry, dv, isp, margin, mode, retention=1.0):
    """Propellant for one round trip: `out` for the laden leg, `back` for the return.

    In mode `earth_chemical` both burns' propellant is loaded where the trip begins. In `split_chemical` each burn's
    propellant is loaded where that burn begins. In `space_chemical` the far end supplies everything: the return
    carries the propellant of the next outward burn, at a handling loss.
    """
    k = math.expm1(dv / (isp * STANDARD_GRAVITY))
    if mode == 'split_chemical':
        a = (cargo + dry) * k * margin
        b = (dry + back) * k * margin
        return dict(out=a, back=b, tank=max(a, b), burned=a + b, space_supply=b / retention)
    if mode == 'earth_chemical':
        b = (dry + back) * k
        a = (cargo + dry + b) * k
        a *= margin
        b *= margin
        return dict(out=a, back=b, tank=a + b, burned=a + b, space_supply=a + b)
    a = (cargo + dry) * k * margin
    carried = a / retention
    b = (dry + carried + back) * k * margin
    return dict(out=a, back=b, tank=carried + b, burned=a + b, space_supply=(carried + b) / retention)

def ferry(cargo, back, route):
    """Algebraic equivalent of the inherited affine tank-limited trip solver."""
    args = {k: route[k] for k in ('dry', 'dv', 'isp', 'margin', 'mode', 'retention')}
    empty = burns(0, 0, **args)
    unit_args = dict(args, dry=0.0)
    one = burns(1, 0, **unit_args)
    ret = burns(0, 1, **unit_args)
    room = route['tank_limit'] - empty['tank']
    if room <= 0:
        raise ValueError('Vehicle cannot carry its own return propellant')
    if route['mode'] == 'split_chemical':
        # Each leg carries its own burn alone, so the tanks limit a leg exactly as the payload does.
        n = max(cargo, back) / route['payload']
    else:
        n = max(cargo / route['payload'], back / route['payload'], (one['tank'] * cargo + ret['tank'] * back) / room)
    if n <= 0:
        return dict(trips=0.0, out=0.0, back=0.0, tank=0.0, burned=0.0, space_supply=0.0)
    b = burns(cargo / n, back / n, **args)
    return dict(trips=n, **{k: v * n for k, v in b.items()})

EARTH_LEGS = {'space': 'space_chemical', 'split': 'split_chemical', 'earth': 'earth_chemical'}


def earth_route(data, year):
    """The ferry between low Earth orbit and the depot, under the arrangement transport.earth_leg names for the year.

    `space` supplies both burns from the depot, the returning ferry carrying the next outward burn's propellant
    down to low Earth orbit. `earth` buys both burns' propellant on Earth and launches it. `split` launches the
    outward burn's propellant and fuels the return at the depot. The founding study's ferry is sized for `space`;
    for the others the ferry keeps its fixed mass and carries tanks for what it holds, at the same mass for each
    kilogram stored.
    """
    t = data['industrial']['transport']
    v = data['original_architecture']['peak_flows']['Earth_vehicle']
    leg = at(t.get('earth_leg', 'space'), year)
    if leg not in EARTH_LEGS:
        raise ValueError('Unknown arrangement of the Earth leg: ' + str(leg))
    route = dict(payload=t['earth_payload_kg'], dry=v['dry_kg'], dv=v['delta_v_oneway_m_s'], isp=v['Isp_s'], margin=t['propellant_margin'], mode='space_chemical', retention=data['transport']['chemical']['handling_retention'], tank_limit=v['burns']['peak_tank_contents_kg'], cycles=JULIAN_YEAR_DAYS / v['cycle_days'])
    if leg == 'space':
        return dict(route, leg=leg, dry_scale=1.0)
    fixed = data['transport']['chemical']['fixed_dry_to_payload'] * route['payload']
    per_stored = (v['dry_kg'] - fixed) / v['burns']['peak_tank_contents_kg']
    args = {k: route[k] for k in ('dv', 'isp', 'margin', 'retention')}
    dry = fixed
    for _ in range(200):
        tank = burns(route['payload'], 0.0, dry, mode=EARTH_LEGS[leg], **args)['tank']
        closed = fixed + per_stored * tank
        if abs(closed - dry) <= 1e-12 * closed:
            break
        dry = closed
    else:
        raise ValueError('The Earth-leg ferry does not close on its tanks')
    return dict(route, dry=dry, mode=EARTH_LEGS[leg], tank_limit=tank, leg=leg, dry_scale=dry / v['dry_kg'])


def routes(data, year):
    s = data['industrial']
    t = s['transport']
    v = data['original_architecture']['peak_flows']['Earth_vehicle']
    out = {'Earth': earth_route(data, year)}
    for name, prefix in [('Moon', 'lander'), ('source', 'nea')]:
        payload = t[prefix + '_payload_kg']
        dry = payload * t[prefix + '_dry_to_payload']
        dv = t[prefix + '_dv_oneway_m_s']
        args = dict(dry=dry, dv=dv, isp=t['Isp_s'], margin=t['propellant_margin'], mode='earth_chemical', retention=1.0)
        tank = burns(payload, 0, **args)['tank']
        if name == 'source':
            tank = payload * data['original_architecture']['peak_flows']['nea_round_trip']['total'] * t['propellant_margin']
        out[name] = dict(payload=payload, **args, tank_limit=tank, cycles=t['lander_cycles_year'] if name == 'Moon' else 1 / t['nea_cycle_years'])
    return out

def catalogue(data):
    """One mass unit per asset; active and structural lifetimes stay distinct."""
    out = deepcopy(data['templates'])
    q = data['economy']
    s = data['industrial']
    ex = q['exports']
    rv = data['revenue']
    for a in out.values():
        a.update(material='metal', civil=False)

    def add(name, owner, price, life, share=0.0, structlife=None, dest='surface', material='metal', civil=False):
        out[name] = dict(name=name, owner=owner, price=price, life=life, structural_share=share, structural_life=structlife or life, destination=dest, material=material, civil=civil, mass=0.0)
    add('export_factory', 'commerce', s['prices']['lunar_hardware_USD2026_kg'], ex['factory_equipment_life_years'], 0.6, 30)
    add('source_computer', 'commerce', ex['compute_IT_USD2026_W'] / ex['compute_IT_kg_W'], ex['compute_life_years'])
    out['source_computer']['destination'] = 'source'
    add('source_communications', 'commerce', ex['compute_communications_USD2026_W'] / ex['compute_communications_kg_W'], ex['compute_life_years'], dest='source')
    sv = rv['services']
    add('service_equipment', 'service', s['prices']['lunar_hardware_USD2026_kg'], sv['research_hardware_life_years'], 0.6, 30)
    add('visitor_equipment', 'hospitality', s['prices']['lunar_hardware_USD2026_kg'], rv['visitors']['equipment_life_years'], 0.6, 30)
    rr = rv['return']
    add('return_capsule_buses', 'commerce', s['prices']['spacecraft_USD2026_kg'], rr['durable_calendar_life_years'], dest='leo')
    for name, owner, mat, life in [('residential_bulk', 'settlement', 'bulk', s['lunar']['bulk_life_years']), ('residential_roads', 'settlement', 'bulk', s['lunar']['road_life_years']), ('residential_metal', 'settlement', 'metal', s['lunar']['metal_life_years']), ('industrial_civil_bulk', 'industry', 'bulk', s['lunar']['bulk_life_years']), ('industrial_civil_metal', 'industry', 'metal', s['lunar']['metal_life_years']), ('commercial_civil_bulk', 'commerce', 'bulk', s['lunar']['bulk_life_years']), ('commercial_civil_metal', 'commerce', 'metal', s['lunar']['metal_life_years']), ('service_civil_bulk', 'service', 'bulk', s['lunar']['bulk_life_years']), ('service_civil_metal', 'service', 'metal', s['lunar']['metal_life_years']), ('visitor_civil_bulk', 'hospitality', 'bulk', s['lunar']['bulk_life_years']), ('visitor_civil_metal', 'hospitality', 'metal', s['lunar']['metal_life_years'])]:
        add(name, owner, q['production']['local_bulk_value_USD2026_kg'] if mat == 'bulk' else q['production']['structural_metal_value_USD2026_kg'], life, 1.0, life, material=mat, civil=True)
    return out

def commercial_targets(data, row):
    ex = data['economy']['exports']
    rr = data['revenue']['return']
    sv = data['revenue']['services']
    s = data['industrial']
    goods = row.get('goods_capacity', row.get('goods', 0.0) + row.get('space_goods', 0.0))
    compute = row.get('compute_capacity', row.get('compute', 0.0))
    plate = goods / (ex['factory_yield'] * ex['factory_utilization'])
    jobs = row.get('service_jobs', 0.0)
    floor = plate * ex['factory_floor_m2_per_nameplate_kg_year']
    service_floor = jobs * sv['research_floor_m2_staff']
    beds = row.get('visitor_beds', 0.0)
    v = data['revenue']['visitors']
    flights = (row.get('goods', 0.0) + row.get('source_goods', 0.0)) * (1 + data['revenue']['manufacturing']['transport_packaging_mass_fraction']) / rr['capsule_payload_kg']
    capsules = math.ceil(flights * s['economics']['capacity_margin'] / (rr['cycles_year'] * s['transport']['fleet_availability']))
    bus = rr['capsule_payload_kg'] * rr['capsule_dry_to_payload'] * (1 - rr['shield_fraction_of_dry'])
    return dict(export_factory=plate * ex['factory_equipment_kg_per_nameplate_kg_year'], source_computer=compute * ex['compute_IT_kg_W'], source_communications=compute * ex['compute_communications_kg_W'], service_equipment=jobs * sv['research_hardware_kg_staff'], return_capsule_buses=capsules * bus, commercial_civil_bulk=floor * ex['factory_bulk_kg_m2'], commercial_civil_metal=floor * ex['factory_metal_kg_m2'], service_civil_bulk=service_floor * sv['bulk_kg_m2'], service_civil_metal=service_floor * sv['metal_kg_m2'], visitor_equipment=beds * v['surface_hardware_kg_bed'], visitor_civil_bulk=beds * v['local_bulk_kg_bed'], visitor_civil_metal=beds * v['local_metal_kg_bed'])

def residential_targets(data, population):
    f = population / data['normalization']['population']
    d = data['demand']
    out = {n: a['mass'] * f for n, a in data['templates'].items() if a['owner'] == 'settlement'}
    out.update(residential_bulk=(d['bulk'] - d['roads']) * f, residential_roads=d['roads'] * f, residential_metal=d['metal'] * f)
    return out

def local_working_parts(data):
    """How far the industry makes its own working parts: the share by order year, its labour, and the assets left out.

    `share` is the share the industry is able to make in an order year, and `staffed` the part of that share its
    residents can staff; their product is the share made. A share of 0 buys every working part from Earth, as the
    staged model does. The share made is fabricated from local metal like structure, and costs resident labour at
    the stated fraction of the part's import price. All three may be schedules by year.
    """
    stated = data['economy']['production'].get('local_working_parts', {})
    return dict(share=stated.get('share', 0.0), staffed=stated.get('staffed', 1.0),
                labour_fraction_of_import_price=stated.get('labour_fraction_of_import_price', 0.0),
                excluded_assets=tuple(stated.get('excluded_assets', ())))


def technology(data):
    """How the industry's technology improves after the founding: the fall in the price of machinery's working parts
    for each doubling of the working parts acquired since the founding deployment, and the yearly growth of labour
    productivity from the stated year, the wage rising with it. Both are nought where not stated, which keeps the
    founding design's technology throughout."""
    stated = data['economy'].get('technology', {})
    out = dict(learning_per_doubling=stated.get('learning_per_doubling', 0.0),
               labour_productivity_growth=stated.get('labour_productivity_growth', 0.0), from_year=stated.get('from_year', 8))
    if not (0.0 <= out['learning_per_doubling'] < 1.0 and out['labour_productivity_growth'] > -1.0):
        raise ValueError('Learning lowers a price by less than all of it per doubling, and productivity grows by more than -100%')
    return out


def labour_index(data, year):
    """Output per job in a year as a multiple of the founding design's; the wage is the same multiple of its own."""
    tech = technology(data)
    return (1.0 + tech['labour_productivity_growth']) ** max(0, year - tech['from_year'])


def learning_index(data, seed, events):
    """The price of working parts ordered in each year as a share of the founding design's, by order year.

    The founding deployment's working parts are the first unit of the curve, and each later order adds its working
    parts to what has been acquired, whether bought from Earth or made in place; orders before the stated year keep the
    founding price. None where no learning is stated."""
    tech = technology(data)
    rate = tech['learning_per_doubling']
    if rate <= 0.0:
        return None
    cats = catalogue(data)
    founding = sum(mass * (1 - cats[n]['structural_share']) for n, mass in seed.items() if n in cats)
    if founding <= 0.0:
        raise ValueError('Learning needs the founding deployment\'s working parts as its first unit')
    ordered = {}
    for ev in events:
        if ev['component'] == 'active':
            ordered[ev['order_year']] = ordered.get(ev['order_year'], 0.0) + ev['mass_kg']
    exponent = math.log2(1 - rate)
    index, total = {}, founding
    for y in sorted(ordered):
        total += ordered[y]
        index[y] = 1.0 if y < tech['from_year'] else (total / founding) ** exponent
    return index


def blank_flow():
    return dict(imports={k: 0.0 for k in DESTS}, local={k: 0.0 for k in DESTS}, recycled={k: 0.0 for k in DESTS}, bulk=0.0, bulk_recycled=0.0, external_capital=0.0, external_replacement=0.0, local_capital=0.0, local_replacement=0.0, service_capital=0.0, service_replacement=0.0, retired_mass=0.0, new_floor=0.0, capsule_purchases_kg=0.0, new_visitor_beds=0.0, industry_replacement_import_kg=0.0, service_external_capital=0.0, service_external_replacement=0.0, hospitality_external_capital=0.0, hospitality_external_replacement=0.0, capital_imports={k: 0.0 for k in DESTS}, capital_local={k: 0.0 for k in DESTS}, capital_bulk=0.0, hospitality_capital=0.0, hospitality_replacement=0.0, local_working_capital=0.0, local_working_replacement=0.0, local_working_labour_if_staffed=0.0)

@lru_cache(maxsize=None)
def lifetime_scale(life, shape):
    # Annual cohorts use beginning-of-year surviving stock. Calibrate the
    # discrete survival sum to the inherited mean; otherwise sampling a
    # continuous Weibull adds approximately half a service year for free.
    lo = 0.0
    hi = 2 * life / math.gamma(1 + 1 / shape)
    for _ in range(60):
        mid = (lo + hi) / 2
        mean = sum((math.exp(-(age / mid) ** shape) for age in range(math.ceil(15 * life) + 1)))
        if mean < life:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

@lru_cache(maxsize=None)
def survival(age, life, shape=3.0):
    if age < 0:
        return 1.0
    if shape is None:
        return float(age < life)
    return math.exp(-(age / lifetime_scale(life, shape)) ** shape)

def cohort_ledger(data, rows, targets, seed, measures=None, price_index=None):
    """Annual components expire at their recorded age; surplus is not replaced.

    Orders and their manufacture/freight occur two years before service. Seed
    hardware has eight deployment cohorts (short-lived units use their life),
    corresponding to the explicitly assumed preproduction deployment window.
    A price index by order year, where given, scales the price of working parts.
    """
    cats = catalogue(data)
    start = rows[0]['year']
    end = rows[-1]['year']
    lead = 2
    flows = {y: blank_flow() for y in range(start - lead, end + 1)}
    stocks = {y: {} for y in range(start, end + 1)}
    events = []
    working = local_working_parts(data)
    shape = data.get('lifetime_shape', 3.0)
    for name, a in cats.items():
        for part, share, life in [('structure', a['structural_share'], a['structural_life']), ('active', 1 - a['structural_share'], a['life'])]:
            if share <= 0:
                continue
            phases = max(1, min(8, math.ceil(life)))
            denominator = sum((survival(j + 1, life, shape) for j in range(phases)))
            if denominator <= 0:
                raise ValueError('Seed deployment ages exceed every component lifetime')
            cohorts = [(start - j - 1, seed.get(name, 0.0) * share / denominator) for j in range(phases)]
            curve = [survival(age, life, shape) for age in range(end - start + phases + 2)]
            for row in rows:
                y = row['year']
                expired = alive = 0.0
                surviving = []
                for birth, mass in cohorts:
                    age = y - birth
                    fraction = curve[age]
                    if y > start:
                        expired += mass * (curve[max(0, age - 1)] - fraction)
                    if fraction > 1e-12:
                        surviving.append((birth, mass))
                        alive += mass * fraction
                cohorts = surviving
                required = targets[y].get(name, 0.0) * share
                added = max(0.0, required - alive)
                if added > 1e-06:
                    cohorts.append((y, added))
                    order = y - lead
                    replacement = min(added, expired)
                    growth = added - replacement
                    able = 0.0 if part != 'active' or name in working['excluded_assets'] else at(working['share'], order)
                    made = able * at(working['staffed'], order)
                    labour = at(working['labour_fraction_of_import_price'], order)
                    if not (0 <= made <= 1 and labour >= 0):
                        raise ValueError('The locally made share of working parts lies between 0 and 1, and its labour is not negative')
                    learned = price_index.get(order, 1.0) if price_index and part == 'active' else 1.0
                    for workyear in (order, order + 1):
                        f = flows[workyear]
                        amount = added / 2
                        growth_work = growth / 2
                        replace_work = replacement / 2
                        unit = a['price'] * learned if part == 'active' or a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
                        if part == 'active':
                            # Working parts: the share bought from Earth, and the share the industry makes itself.
                            bought = 1 - made
                            f['imports'][a['destination']] += amount * bought
                            f['external_capital'] += growth_work * bought * unit
                            f['external_replacement'] += replace_work * bought * unit
                            f['capital_imports'][a['destination']] += growth_work * bought
                            if a['owner'] == 'industry':
                                f['industry_replacement_import_kg'] += replace_work * bought
                            if a['owner'] in ('service', 'hospitality'):
                                f[a['owner'] + '_external_capital'] += growth_work * bought * unit
                                f[a['owner'] + '_external_replacement'] += replace_work * bought * unit
                            if name == 'return_capsule_buses':
                                f['capsule_purchases_kg'] += amount
                            f['local'][a['destination']] += amount * made
                            f['capital_local'][a['destination']] += growth_work * made
                            f['local_working_capital'] += growth_work * made * unit * labour
                            f['local_working_replacement'] += replace_work * made * unit * labour
                            f['local_working_labour_if_staffed'] += amount * able * unit * labour
                        elif a['material'] == 'bulk':
                            f['bulk'] += amount
                            f['capital_bulk'] += growth_work
                        else:
                            f['local'][a['destination']] += amount
                            f['capital_local'][a['destination']] += growth_work
                        if part == 'structure':
                            f['local_capital'] += growth_work * unit
                            f['local_replacement'] += replace_work * unit
                        if a['owner'] == 'service':
                            f['service_capital'] += growth_work * unit
                            f['service_replacement'] += replace_work * unit
                        if a['owner'] == 'hospitality':
                            f['hospitality_capital'] += growth_work * unit
                            f['hospitality_replacement'] += replace_work * unit
                        if name in ('commercial_civil_bulk', 'service_civil_bulk'):
                            f['new_floor'] += growth_work / data['economy']['exports']['factory_bulk_kg_m2']
                    events.append(dict(asset=name, component=part, order_year=order, commission_year=y, mass_kg=added, replacement_kg=replacement, locally_made_share=made, local_labour_fraction=labour, price_index=learned))
                if expired:
                    f = flows[y]
                    f['retired_mass'] += expired
                    if part == 'structure':
                        if a['material'] == 'bulk':
                            f['bulk_recycled'] += expired * data['industrial']['lunar']['recovery']
                        else:
                            f['recycled'][a['destination']] += expired * (data['industrial']['lunar']['recovery'] if a['destination'] == 'surface' else data['economy']['production']['source_recovery'])
                available = alive + added
                old = stocks[y].get(name, math.inf)
                stocks[y][name] = min(old, available / share)
    previous_beds = rows[0].get('visitor_beds', 0.0)
    for row in rows:
        beds = row.get('visitor_beds', 0.0)
        added = max(0.0, beds - previous_beds)
        previous_beds = beds
        for y in (row['year'] - 2, row['year'] - 1):
            flows[y]['new_visitor_beds'] += added / 2
    if measures is not None:
        measures.after_ledger(data, rows, seed, flows, events)
    return (flows, stocks, events)

def physical(data, row, f, bank=None, measures=None):
    """Mass balance followed by the inherited chemical freight and process laws."""
    s = data['industrial']
    d = data['demand']
    z = s['source']
    l = s['lunar']
    q = data['economy']
    ex = q['exports']
    rv = data['revenue']
    rr = rv['return']
    t = s['transport']
    p = q['production']
    ps = row['population'] / data['normalization']['population']
    gs = row.get('new_residential_people', 0.0) / data['normalization']['population']
    bank = deepcopy(bank or {**{k: 0.0 for k in DESTS}, 'bulk': 0.0})
    imports = dict(f['imports'])
    local = dict(f['local'])
    virgin = {}
    for dest in DESTS:
        bank[dest] += f['recycled'][dest]
        used = min(local[dest], bank[dest])
        bank[dest] -= used
        virgin[dest] = local[dest] - used
    bank['bulk'] += f['bulk_recycled']
    rebulk = min(f['bulk'], bank['bulk'])
    bank['bulk'] -= rebulk
    goods = row.get('goods', 0.0)
    source_goods = row.get('source_goods', 0.0)
    fuel = row.get('propellant_sales', 0.0)
    compute = row.get('compute', 0.0)
    # Goods handed to customers at the depot are made as the goods for Earth are and travel no further.
    space_goods = row.get('space_goods', 0.0)
    space_source_goods = row.get('space_source_goods', 0.0)
    # Travellers: the mass launched for them from Earth, and the cabins that carry them each way beyond low orbit.
    travellers_launched = row.get('passenger_launch_kg_year', 0.0)
    cabins = row.get('passenger_cabin_kg_year', 0.0)
    resident_consumables = sum((v for k, v in data['settlement']['operating_import_kg_year'].items() if k not in ('water', 'air'))) * ps
    industrial_consumables = p['process_consumables_fraction_of_imported_renewal'] * f.get('industry_replacement_import_kg', 0.0)
    guest_years = row.get('visitor_nights_year', 0.0) / JULIAN_YEAR_DAYS
    visitor_consumables = resident_consumables / row['population'] * guest_years if row['population'] else 0.0
    research_consumables = row.get('service_jobs', 0.0) * rv['services']['research_consumables_kg_staff_year']
    imports['surface'] += resident_consumables + visitor_consumables + research_consumables
    imports['source'] += industrial_consumables
    factory_feed = (goods + space_goods) / ex['factory_yield']
    earthfeed = factory_feed * ex['factory_Earth_input_fraction']
    imports['surface'] += earthfeed
    factory_local = factory_feed - earthfeed + goods * rv['manufacturing']['transport_packaging_mass_fraction']
    source_product = (source_goods + space_source_goods) / ex['factory_yield'] + source_goods * rv['manufacturing']['transport_packaging_mass_fraction']
    surface_alloy = virgin['surface'] * p['alloy_import_fraction_of_virgin_local_metal']
    source_alloy = (sum((virgin[k] for k in ('source', 'depot', 'leo'))) + source_product) * p['alloy_import_fraction_of_virgin_local_metal']
    imports['surface'] += surface_alloy
    imports['source'] += source_alloy
    freshbulk = f['bulk'] - rebulk + d['goods_virgin'] * ps
    freshmetal = virgin['surface'] - surface_alloy + factory_local
    bulk = f['bulk'] + d['goods_virgin'] * ps
    metal = local['surface'] + factory_local
    raw = freshbulk / l['bulk_yield'] + freshmetal / l['metal_yield']
    handling = raw + d['loose_and_grading'] * gs
    lunar_W = (bulk * l['bulk_kWh_kg'] + metal * l['metal_kWh_kg'] + handling * l['handling_kWh_kg'] + local['surface'] * p['local_fabrication_finish_kWh_kg'] + (goods + space_goods) * ex['factory_finish_kWh_kg']) * 1000 / HOURS + row.get('service_power_W', 0.0) + row.get('premium_power_W', 0.0) + row.get('visitor_beds', 0.0) * rv['visitors']['mean_W_bed']
    water = d['yearly_water'] * (ps + guest_years / data['normalization']['population']) + d['water'] * gs + f.get('new_visitor_beds', 0.0) * rv['visitors']['initial_water_kg_bed']
    air = d['yearly_air'] * (ps + guest_years / data['normalization']['population']) + d['air'] * gs + f.get('new_visitor_beds', 0.0) * rv['visitors']['initial_air_kg_bed']
    # New commercial enclosure air follows the explicitly priced floor proxy.
    newfloor = f['new_floor']
    air += newfloor * rv['manufacturing']['air_kg_m2']
    packed = (goods + source_goods) * (1 + rv['manufacturing']['transport_packaging_mass_fraction'])
    source_packed = source_goods * (1 + rv['manufacturing']['transport_packaging_mass_fraction'])
    flights = packed / rr['capsule_payload_kg']
    bus = rr['capsule_payload_kg'] * rr['capsule_dry_to_payload'] * (1 - rr['shield_fraction_of_dry'])
    shield = rr['capsule_payload_kg'] * rr['capsule_dry_to_payload'] * rr['shield_fraction_of_dry']
    deorbit = (rr['capsule_payload_kg'] + bus + shield) * math.expm1(rr['deorbit_dv_m_s'] / (rr['deorbit_Isp_s'] * STANDARD_GRAVITY))
    # Relaunch used buses once per sortie; newly bought buses are already within
    # this stream and must not receive a second launch charge in their purchase year.
    imports['leo'] += flights * (bus + shield + deorbit) - f.get('capsule_purchases_kg', 0.0)
    recycle_back = f['recycled']['depot'] + f['recycled']['leo']
    leo_recycle = f['recycled']['leo']
    r = routes(data, row['year'])
    earth = ferry(imports['surface'] + imports['depot'] + imports['source'] + leo_recycle + cabins, packed + local['leo'] + cabins, r['Earth'])
    moon = ferry(imports['surface'] + water + air + cabins, packed - source_packed + space_goods + cabins, r['Moon'])
    nitrogen = air * z['air_nitrogen_mass_fraction']
    oxygen = air - nitrogen
    # The Earth leg's propellant: what the depot supplies from asteroid water, and what is bought on Earth and launched.
    leg = r['Earth']['leg']
    earth_space = 0.0 if leg == 'earth' else earth['space_supply']
    earth_launched = earth['burned'] if leg == 'earth' else earth['out'] if leg == 'split' else 0.0
    depot_water = moon['burned'] + oxygen + earth_space + fuel
    source_cargo = water + nitrogen + depot_water + local['depot'] + local['leo'] + source_packed + space_source_goods
    source = ferry(source_cargo, imports['source'] + recycle_back, r['source'])
    source_water = water + depot_water + source['burned']
    source_metal = sum((local[k] for k in ('source', 'depot', 'leo'))) + source_product
    source_virgin = sum((virgin[k] for k in ('source', 'depot', 'leo'))) + source_product - source_alloy
    feed = max(source_water / (z['water_grade'] * z['water_recovery']), nitrogen / (z['nitrogen_grade'] * z['nitrogen_recovery']), source_virgin / (z['metal_grade'] * z['metal_recovery']))
    process_kWh = feed * z['ore_kWh_kg'] + nitrogen * z['nitrogen_kWh_kg'] + source_metal * z['metal_refining_and_fabrication_kWh_kg']
    chem = source['burned'] / z['water_to_hydrogen_mass_ratio'] * z['electrolysis_liquefaction_kWh_kg_H2']
    source_W = (process_kWh + chem) * 1000 / HOURS + compute
    heat = (process_kWh + chem * (1 - z['hydrogen_HHV_kWh_kg'] / z['electrolysis_liquefaction_kWh_kg_H2'])) * 1000 / HOURS + compute
    depot_W = depot_water / z['water_to_hydrogen_mass_ratio'] * z['electrolysis_liquefaction_kWh_kg_H2'] * 1000 / HOURS
    ground = d['occupied_W'] * ps + lunar_W * s['economics']['capacity_margin'] / l['process_availability']
    operating_material = (resident_consumables + visitor_consumables + research_consumables + industrial_consumables + surface_alloy + source_alloy + earthfeed) * s['prices']['other_special_USD2026_kg']
    capsule_op = flights * (shield * s['prices']['spacecraft_USD2026_kg'] + deorbit * t['earth_propellant_USD2026_kg']) + packed * rr['Earth_ground_service_USD2026_per_payload_kg']
    LEO = sum(imports.values()) + earth_launched + travellers_launched
    out = dict(lunar_raw_kg_year=raw, lunar_bulk_kg_year=bulk, lunar_metal_kg_year=metal, lunar_recycled_bulk_kg_year=rebulk, lunar_virgin_metal_kg_year=freshmetal, lunar_handling_kg_year=handling, lunar_process_W=lunar_W, source_feed_kg_year=feed, source_water_kg_year=source_water, source_nitrogen_kg_year=nitrogen, source_virgin_metal_kg_year=source_virgin, source_metal_kg_year=source_metal, source_process_W=source_W, source_heat_W=heat, source_feed_process_W=feed * z['ore_kWh_kg'] * 1000 / HOURS, source_product_refining_W=(nitrogen * z['nitrogen_kWh_kg'] + source_metal * z['metal_refining_and_fabrication_kWh_kg']) * 1000 / HOURS, depot_W=depot_W, depot_water_kg_year=depot_water, ground_firm_W=ground, Earth=earth, Moon=moon, source=source, source_cargo_kg_year=source_cargo, source_cargo_breakdown_kg_year=dict(resident_and_guest_water=water, nitrogen=nitrogen, depot_water=depot_water, local_offsite_structures=local['depot'] + local['leo'], source_exports=source_packed + space_source_goods), depot_water_breakdown_kg_year=dict(lunar_ferry_fuel=moon['burned'], oxygen_makeup=oxygen, Earth_ferry_fuel=earth_space, external_fuel_sales=fuel), source_water_breakdown_kg_year=dict(resident_and_guest_water=water, depot_water=depot_water, source_ferry_fuel=source['burned']), source_retained_tails_kg_year=feed - source_water - nitrogen - source_virgin, LEO_kg_year=LEO, imports_by_destination_kg_year=imports, local_fabrication_by_destination_kg_year=local, recycled_inventory_kg=bank, factory_floor_new_m2=newfloor, construction_floor_m2=data['settlement']['production_demand']['floor'] * gs + newfloor, operating_material_USD2026=operating_material, capsule_operations_USD2026=capsule_op, Earth_launch_USD2026=LEO * at(t['LEO_USD2026_kg'], row['year']) + earth_launched * t['earth_propellant_USD2026_kg'], Earth_leg_propellant_launched_kg_year=earth_launched, resource_provision_USD2026=raw * s['economics']['local_restoration_USD2026_kg_feed'] + feed * rv['finance']['source_transition_USD2026_per_kg_feed'])
    if measures is not None:
        measures.after_physical(data, row, f, out)
    return out

def requirements(data, row, f, measures=None):
    """Rated machinery required by actual flows, with inherited spare margins."""
    s = data['industrial']
    l = s['lunar']
    z = s['source']
    p = s['power']
    t = s['transport']
    d = data['demand']
    margin = s['economics']['capacity_margin']
    original = data['original_architecture']
    out = residential_targets(data, row['population'])
    out.update(commercial_targets(data, row))

    def machine(name, flow, rate, unit, availability):
        out[name] = math.ceil(flow * margin / (rate * HOURS * availability)) * unit
    machine('excavators', f['lunar_handling_kg_year'], l['excavator_kg_hr'], l['excavator_kg'], l['mine_availability'])
    machine('haulers', f['lunar_handling_kg_year'], l['hauler_payload_kg'] / l['haul_roundtrip_hr'], l['hauler_kg'], l['mine_availability'])
    for prefix, key in [('', 'lunar_metal_kg_year'), ('source_', 'source_metal_kg_year')]:
        machine(prefix + 'metal_trains', f[key], l['metal_train_kg_hr'], l['metal_train_hardware_kg'], l['process_availability'])
        machine(prefix + 'fabrication_cells', f[key], l['fabrication_cell_kg_hr'], l['fabrication_cell_hardware_kg'], l['process_availability'])
    machine('bulk_trains', f['lunar_bulk_kg_year'], l['bulk_train_kg_hr'], l['bulk_train_hardware_kg'], l['process_availability'])
    machine('source_mine_trains', f['source_feed_kg_year'], z['mine_module_kg_hr'], z['mine_module_hardware_kg'], z['mine_availability'])
    out['construction_and_service_robots'] = l['construction_plant_kg'] * f['construction_floor_m2'] / data['settlement']['production']['floor_m2_year']
    ground = f['ground_firm_W']
    industrialfirm = f['lunar_process_W'] * margin / l['process_availability']
    out['additional_ground_receivers'] = industrialfirm * d['power_hardware_kg_W']
    area = industrialfirm * d['old_radiator_m2'] / d['old_total_W']
    out['additional_ground_radiator_controls'] = area * d['radiator_controls_kg_m2']
    sourceDC = ground / (p['RF_efficiency'] * p['beam_capture'] * p['receiver_efficiency']) * p['stations'] / 2 / p['end_of_life_factor']
    out['lunar_beam_constellation'] = sourceDC / p['space_system_W_kg']
    radflux = p['emissivity'] * STEFAN_BOLTZMANN * p['radiator_K'] ** 4
    for pre, W, heat in [('source', f['source_process_W'], f['source_heat_W']), ('depot', f['depot_W'], f['depot_W'] * (1 - z['hydrogen_HHV_kWh_kg'] / z['electrolysis_liquefaction_kWh_kg_H2']))]:
        chemical_W = max(0.0, W - row.get('compute', 0.0)) if pre == 'source' else W
        out[pre + '_chemistry'] = chemical_W * margin * z['chemistry_hardware_kg_W']
        out[pre + '_solar'] = W * margin / p['end_of_life_factor'] / (p['source_solar_flux_fraction'] if pre == 'source' else 1) / p['space_system_W_kg']
        out[pre + '_radiators'] = heat * margin / radflux * p['radiator_kg_m2']
    r = routes(data, row['year'])
    counts = {name: math.ceil(f[name]['trips'] * margin / (a['cycles'] * t['fleet_availability'])) for name, a in r.items()}
    out['lunar_landers'] = counts['Moon'] * r['Moon']['dry']
    out['asteroid_tugs'] = counts['source'] * r['source']['dry']
    for name, a in original['assets'].items():
        if name.startswith('Earth_'):
            out[name] = counts['Earth'] * a['mass_kg'] / original['engineering']['counts']['Earth_route_ships'] * r['Earth']['dry_scale']
    out['depots_and_orbital_workshops'] = (f['source_cargo_kg_year'] + f['Moon']['burned']) * t['depot_inventory_years'] * t['tank_dry_to_contents']
    # Propellant held in low Earth orbit and at the depot for the Earth leg, by where each burn is fuelled.
    leg = r['Earth']['leg']
    held = dict(LEO=f['Earth']['burned'] if leg == 'earth' else f['Earth']['out'], depot=0.0 if leg == 'earth' else f['Earth']['space_supply'])
    for name, place in [('LEO_route_tanks', 'LEO'), ('LEO_route_storage_controls', 'LEO'), ('depot_route_storage_controls', 'depot')]:
        peak = original['peak_flows']
        oldstock = peak['Earth_route_LEO_inventory_kg'] if place == 'LEO' else peak['Earth_route_depot_inventory_kg']
        out[name] = original['assets'][name]['mass_kg'] * held[place] * data['transport']['depots']['inventory_years'] / oldstock
    surface = sum((v for n, v in out.items() if n in data['templates'] and data['templates'][n]['owner'] == 'industry' and (data['templates'][n]['destination'] == 'surface')))
    out['industrial_civil_bulk'] = surface * l['infrastructure_bulk_per_hardware_kg']
    out['industrial_civil_metal'] = surface * l['infrastructure_metal_per_hardware_kg'] + area * d['radiator_metal_kg_m2']
    if measures is not None:
        measures.after_requirements(data, row, f, out)
    return out

def seed_stock(data, first, bootstrap_years=8, measures=None):
    """Size one seed against construction throughput, then retain its assets.

    Complete seed industrial equipment is procured from outside. This avoids
    assuming an unbuilt local industry manufactures its own first active plant.
    The eight-year average is a throughput screen, not a qualified deployment.
    """
    cats = catalogue(data)
    base = residential_targets(data, first['population'])
    base.update(commercial_targets(data, first))
    seed = dict(base)
    # The founding deployment's average year stands at the middle of its years, for coefficients that follow a schedule.
    row = dict(first, new_residential_people=first['population'] / bootstrap_years, population=first['population'] / 2, year=first['year'] - bootstrap_years // 2)

    def demands(stock):
        f = blank_flow()
        for n, mass in stock.items():
            a = cats[n]
            share = a['structural_share']
            amount = mass / bootstrap_years

            def gross_factor(life):
                phases = max(1, min(8, math.ceil(life)))
                average = sum((survival(j + 1, life, data.get('lifetime_shape', 3.0)) for j in range(phases))) / phases
                return 1 / average
            structural_factor = gross_factor(a['structural_life'])
            active_factor = gross_factor(a['life'])
            if a['owner'] == 'industry' and (not a['civil']):
                amount *= share * structural_factor + (1 - share) * active_factor
                f['imports'][a['destination']] += amount
                f['external_capital'] += amount * a['price']
            else:
                imported = amount * (1 - share) * active_factor
                made = amount * share * structural_factor
                f['imports'][a['destination']] += imported
                f['external_capital'] += imported * a['price']
                if a['material'] == 'bulk':
                    f['bulk'] += made
                else:
                    f['local'][a['destination']] += made
                unit = a['price'] if a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
                f['local_capital'] += made * unit
                if a['owner'] == 'service':
                    f['service_capital'] += imported * a['price'] + made * unit
                    f['service_external_capital'] += imported * a['price']
        ex = data['economy']['exports']
        sv = data['revenue']['services']
        f['new_floor'] = (first.get('goods_capacity', first.get('goods', 0.0)) / (ex['factory_yield'] * ex['factory_utilization']) * ex['factory_floor_m2_per_nameplate_kg_year'] + first.get('service_jobs', 0.0) * sv['research_floor_m2_staff']) / bootstrap_years
        f['new_visitor_beds'] = first.get('visitor_beds', 0.0) / bootstrap_years
        return f
    history = []
    for iteration in range(300):
        flow = demands(seed)
        phys = physical(data, row, flow, measures=measures)
        needed = requirements(data, row, phys, measures)
        needed.update(base)
        # Seed also supports opening operation and the first construction order.
        opening = requirements(data, first, physical(data, first, blank_flow(), measures=measures), measures)
        needed = {k: max(v, opening.get(k, 0.0)) for k, v in needed.items()}
        error = max((max(0.0, v - seed.get(k, 0.0)) / max(1.0, v) for k, v in needed.items()), default=0.0)
        history.append(error)
        seed = {k: max(seed.get(k, 0.0), v) for k, v in needed.items()}
        if error < 1e-09:
            break
    else:
        raise ValueError('Seed stock/flow requirements did not converge')
    flow = demands(seed)
    phys = physical(data, row, flow, measures=measures)
    I = data['normalization']['USD2026_to_USD2041']
    original_value = sum((a['mass'] * a['price'] for a in cats.values() if a['owner'] == 'industry' and (not a['civil'])))
    seed_value = sum((seed.get(n, 0.0) * a['price'] for n, a in cats.items() if a['owner'] == 'industry' and (not a['civil'])))
    ratio = seed_value / original_value
    econ = data['industrial']['economics']
    staff = econ['staff'] * ratio
    resident_jobs = staff / 4 + (data['settlement']['food_operators'] + data['settlement']['industry_operators']) * first['population'] / data['normalization']['population'] / 2
    resident_payroll = resident_jobs * econ['loaded_staff_USD2026_year'] * bootstrap_years * I
    external_payroll = staff / 4 * econ['loaded_staff_USD2026_year'] * bootstrap_years * I
    cost = dict(external_hardware_USD2041=flow['external_capital'] * bootstrap_years * I, resident_payroll_USD2041=resident_payroll, external_payroll_USD2041=external_payroll, external_engineering_support_USD2041=data['economy']['production']['external_engineering_contract_USD2026_year'] * ratio / 2 * bootstrap_years * I, delivery_and_operating_material_USD2041=(phys['Earth_launch_USD2026'] + phys['operating_material_USD2026'] + phys['capsule_operations_USD2026']) * bootstrap_years * I, resource_provision_USD2041=phys['resource_provision_USD2026'] * bootstrap_years * I, common_programme_USD2041=data['common_programme_USD2041'], identified_service_hardware_USD2041=flow['service_external_capital'] * bootstrap_years * I)
    return (seed, dict(stock_mass_kg=seed, bootstrap_years=bootstrap_years, cost=cost, resident_payroll_USD2041=resident_payroll, local_fabrication_account_value_USD2041=flow['local_capital'] * bootstrap_years * I, lunar_raw_kg=phys['lunar_raw_kg_year'] * bootstrap_years, source_feed_kg=phys['source_feed_kg_year'] * bootstrap_years, physical_annual_average=phys, iterations=iteration + 1, maximum_relative_requirement_residual=history[-1], commission_years=list(range(first['year'] - 8, first['year'])), evidence='Seed cohorts commissioned in the eight years before opening; gross procurement includes survival losses to opening. Complete industrial equipment imported once; local enclosure/civil construction and average half-occupied seed operations included. Detailed preproduction staging, staff accommodation and source qualification remain unsized.'))

def solve(data, annual, *, measures=None, max_iterations=100, tolerance=1e-07):
    """Close dated component procurements and demand-sized productive capacity.

    Capacity plans increase monotonically within the fixed-point solution.
    Each year's target is independent of the population multiplier and can
    decline after construction demand subsides. Surviving surplus stays on hand.
    """
    rows = deepcopy(annual)
    if not rows or any((b['year'] != a['year'] + 1 for a, b in zip(rows, rows[1:]))):
        raise ValueError('Consecutive annual request rows required')
    if any((row.get(k, 0) < 0 for row in rows for k in ('population', 'goods', 'source_goods', 'space_goods', 'space_source_goods', 'compute', 'propellant_sales', 'new_residential_people', 'passenger_launch_kg_year', 'passenger_cabin_kg_year'))):
        raise ValueError('Negative physical demand')
    # Two forward years prevent final observations receiving free expired plant.
    end = rows[-1]['year']
    for y in range(end + 1, end + 41):
        rows.append(dict(rows[-1], year=y, new_residential_people=0.0, commissioned_people=0.0))
    seed, seed_record = seed_stock(data, rows[0], measures=measures)
    targets = {r['year']: requirements(data, r, physical(data, r, blank_flow(), measures=measures), measures) for r in rows}
    residuals = []
    for iteration in range(max_iterations):
        flows, stocks, events = cohort_ledger(data, rows, targets, seed, measures)
        bank = None
        maximum = 0.0
        next_targets = {}
        for index, row in enumerate(rows):
            workrow = dict(row, new_residential_people=(row.get('new_residential_people', 0.0) + (rows[index - 1].get('new_residential_people', 0.0) if index else 0.0)) / 2)
            y = row['year']
            f = physical(data, workrow, flows[y], bank, measures)
            bank = f['recycled_inventory_kg']
            need = requirements(data, row, f, measures)
            previous = targets[y]
            error = max((max(0.0, v - previous.get(k, 0.0)) / max(1.0, v) for k, v in need.items()), default=0.0)
            maximum = max(maximum, error)
            next_targets[y] = {k: max(v, previous.get(k, 0.0)) for k, v in need.items()}
        residuals.append(maximum)
        targets = next_targets
        if maximum < tolerance:
            break
    converged = residuals[-1] < tolerance
    flows, stocks, events = cohort_ledger(data, rows, targets, seed, measures)
    # The quantities do not depend on prices, so learning prices the same orders once they are known.
    learned = learning_index(data, seed, events)
    if learned:
        flows, stocks, events = cohort_ledger(data, rows, targets, seed, measures, learned)
    bank = None
    result = []
    I = data['normalization']['USD2026_to_USD2041']
    cats = catalogue(data)
    original_value = sum((a['mass'] * a['price'] for a in cats.values() if a['owner'] == 'industry' and (not a['civil'])))
    cumulative_lunar = cumulative_source = 0.0
    for index, row in enumerate(rows):
        y = row['year']
        f = flows[y]
        opening_bank = deepcopy(bank)
        workrow = dict(row, new_residential_people=(row.get('new_residential_people', 0.0) + (rows[index - 1].get('new_residential_people', 0.0) if index else 0.0)) / 2)
        physical_row = physical(data, workrow, f, bank, measures)
        bank = physical_row['recycled_inventory_kg']
        need = requirements(data, row, physical_row, measures)
        ratios = {n: need[n] / max(1e-30, stocks[y].get(n, 0.0)) for n in need if need[n] > 1e-08}
        # Capacity comparison uses required machinery including the inherited25%
        # margin. A ratio>1 is an unsatisfied requirement, never added capacity.
        maximum = max(ratios.values(), default=0.0)
        feasible = converged and maximum <= 1 + 1e-06
        industrial_value = sum((min(stocks[y].get(n, 0.0), need.get(n, 0.0)) * a['price'] for n, a in cats.items() if a['owner'] == 'industry' and (not a['civil'])))
        # Labour productivity: each job does more as the years pass, and its wage rises in step.
        output_per_job = labour_index(data, y)
        wage = data['industrial']['economics']['loaded_staff_USD2026_year'] * output_per_job
        staffing = data['industrial']['economics']['staff'] * industrial_value / original_value / output_per_job
        ex = data['economy']['exports']
        rv = data['revenue']
        P0 = data['normalization']['population']
        factory_jobs = row.get('goods_capacity', row.get('goods', 0.0) + row.get('space_goods', 0.0)) / (ex['factory_yield'] * ex['factory_utilization']) / (P0 * rv['manufacturing']['nameplate_finished_kg_per_resident_year']) * rv['manufacturing']['staff'] / output_per_job
        compute_jobs = row.get('compute_capacity', row.get('compute', 0.0)) / 1000000.0 * rv['computing']['staff_per_MW'] / output_per_job
        base_jobs = (data['settlement']['food_operators'] + data['settlement']['industry_operators']) * row['population'] / P0 / output_per_job
        working_parts_jobs = (f['local_working_capital'] + f['local_working_replacement']) / wage
        core_resident_jobs = staffing / 2 + base_jobs + factory_jobs + compute_jobs + working_parts_jobs
        commercial_jobs = row.get('service_jobs', 0.0) + row.get('premium_jobs', 0.0) + row.get('logistics_jobs', 0.0) + row.get('visitor_beds', 0.0) * rv['visitors']['staff_per_bed'] / output_per_job
        total_resident_jobs = core_resident_jobs + commercial_jobs
        ratios['resident_workforce'] = total_resident_jobs / max(1.0, row['population'] * WORKFORCE_SHARE)
        maximum = max(ratios.values(), default=0.0)
        feasible = converged and maximum <= 1 + 1e-06
        payroll = core_resident_jobs * wage
        external_staff = staffing / 2 * wage + data['economy']['production']['external_engineering_contract_USD2026_year'] * industrial_value / original_value
        # Evaluate the same physical system without new capital manufacture.
        # Replacement, occupied residents, sold outputs and recycling stay in it.
        opflow = deepcopy(f)
        for dest in DESTS:
            opflow['imports'][dest] -= f['capital_imports'][dest]
            opflow['local'][dest] -= f['capital_local'][dest]
        opflow['bulk'] -= f['capital_bulk']
        opflow['new_floor'] = 0.0
        opflow['new_visitor_beds'] = 0.0
        # Capsule purchase masses travel in the sortie stream; subtract only
        # replacement purchases when forming the operating counterfactual.
        capbus = sum((ev['mass_kg'] - ev['replacement_kg'] for ev in events if ev['asset'] == 'return_capsule_buses' and ev['order_year'] in (y, y - 1))) / 2
        opflow['capsule_purchases_kg'] = max(0.0, f['capsule_purchases_kg'] - capbus)
        op_physical = physical(data, dict(workrow, new_residential_people=0.0), opflow, opening_bank, measures)
        variable = ('Earth_launch_USD2026', 'operating_material_USD2026', 'capsule_operations_USD2026')
        construction_logistics = sum((physical_row[k] - op_physical[k] for k in variable))
        cap_resource = physical_row['resource_provision_USD2026'] - op_physical['resource_provision_USD2026']
        production_value = (row.get('goods', 0.0) + row.get('source_goods', 0.0) + row.get('space_goods', 0.0) + row.get('space_source_goods', 0.0)) * data['economy']['production']['structural_metal_value_USD2026_kg']
        construction_fraction = f['local_capital'] / max(1.0, f['local_capital'] + f['local_replacement'] + production_value)
        construction_payroll = staffing / 2 * wage * construction_fraction
        costs = dict(capital_external_USD2041=(f['external_capital'] + construction_logistics) * I, capital_local_USD2041=(construction_payroll + f['local_working_capital']) * I, local_capital_account_value_USD2041=f['local_capital'] * I, local_replacement_account_value_USD2041=f['local_replacement'] * I, capital_resource_provision_USD2041=cap_resource * I, operating_external_replacement_USD2041=f['external_replacement'] * I, operating_external_logistics_materials_USD2041=sum((op_physical[k] for k in variable)) * I, operating_external_services_USD2041=external_staff * I, operating_local_replacement_USD2041=f['local_working_replacement'] * I, operating_local_payroll_USD2041=(payroll - construction_payroll - f['local_working_capital'] - f['local_working_replacement']) * I, resident_payroll_USD2041=payroll * I, operating_resource_provision_USD2041=op_physical['resource_provision_USD2026'] * I, resource_provision_USD2041=physical_row['resource_provision_USD2026'] * I, identified_service_capital_USD2041=(f['service_external_capital'] + construction_payroll * max(0.0, f['service_capital'] - f['service_external_capital']) / max(1.0, f['local_capital'])) * I, identified_service_replacement_USD2041=f['service_external_replacement'] * I, identified_hospitality_capital_USD2041=f['hospitality_external_capital'] * I, identified_hospitality_replacement_USD2041=f['hospitality_external_replacement'] * I, Earth_launch_USD2041=physical_row['Earth_launch_USD2026'] * I, operating_material_USD2041=physical_row['operating_material_USD2026'] * I, capsule_operations_USD2041=physical_row['capsule_operations_USD2026'] * I)
        service_stock = service_mean_renewal = 0.0
        for n, mass in commercial_targets(data, row).items():
            a = cats[n]
            if a['owner'] != 'service':
                continue
            share = a['structural_share']
            localprice = a['price'] if a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
            service_stock += mass * ((1 - share) * a['price'] + share * localprice)
            service_mean_renewal += mass * ((1 - share) * a['price'] / a['life'] + share * localprice / a['structural_life'])
        costs['identified_service_stock_value_USD2041'] = service_stock * I
        costs['identified_service_mean_renewal_allowance_USD2041'] = service_mean_renewal * I
        costs['operating_total_USD2041'] = sum((costs[k] for k in ('operating_external_replacement_USD2041', 'operating_external_logistics_materials_USD2041', 'operating_external_services_USD2041', 'operating_local_replacement_USD2041', 'operating_local_payroll_USD2041', 'operating_resource_provision_USD2041')))
        costs['capital_total_USD2041'] = costs['capital_external_USD2041'] + costs['capital_local_USD2041'] + costs['capital_resource_provision_USD2041']
        direct_total = (f['external_capital'] + f['external_replacement'] + sum((physical_row[k] for k in variable)) + external_staff + payroll + physical_row['resource_provision_USD2026']) * I
        costs['accounting_residual_USD2041'] = costs['capital_total_USD2041'] + costs['operating_total_USD2041'] - direct_total
        yearly_events = [ev for ev in events if ev['order_year'] in (y, y - 1)]
        drivers = {}
        for ev in yearly_events:
            n = ev['asset']
            a = cats[n]
            unit = a['price'] if ev['component'] == 'active' or a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
            if ev['component'] == 'active':
                drivers[n] = drivers.get(n, 0.0) + ev['mass_kg'] * (1 - ev['locally_made_share']) / 2 * unit * ev.get('price_index', 1.0) * I
        tug_import = sum((ev['mass_kg'] * (1 - ev['locally_made_share']) / 2 for ev in yearly_events if ev['asset'] == 'asteroid_tugs' and ev['component'] == 'active'))
        physical_row['source_tug_component_import_kg_year'] = tug_import
        physical_row['source_tug_count'] = stocks[y].get('asteroid_tugs', 0.0) / routes(data, y)['source']['dry']
        cumulative_lunar += physical_row['lunar_raw_kg_year']
        cumulative_source += physical_row['source_feed_kg_year']
        if y <= end:
            result.append(dict(year=y, population=row['population'], request=row, physical=physical_row, cash=costs, industrial_staff_total=staffing, technical_resident_jobs=total_resident_jobs, core_resident_jobs=core_resident_jobs, working_parts_jobs=working_parts_jobs, commercial_resident_jobs=commercial_jobs, top_imported_component_purchases_USD2041=dict(sorted(drivers.items(), key=lambda kv: kv[1], reverse=True)[:8]), required_asset_mass_kg=need, available_asset_mass_kg=stocks[y], capacity_use=ratios, maximum_capacity_use=maximum, feasible=feasible, actual_sales_capacity_fraction=1.0 if feasible else 0.0, shortfalls={k: v - 1 for k, v in ratios.items() if v > 1 + 1e-06}, cumulative_lunar_raw_kg=cumulative_lunar, cumulative_source_feed_kg=cumulative_source))
    # Procurements before opening belong to the seed once, not year8 operations.
    late = blank_flow()
    for y in range(rows[0]['year'] - 2, rows[0]['year']):
        f = flows[y]
        for dest in DESTS:
            late['imports'][dest] += f['imports'][dest]
            late['local'][dest] += f['local'][dest]
        for k in ('bulk', 'external_capital', 'external_replacement', 'local_capital', 'local_replacement', 'service_capital', 'service_replacement', 'new_floor', 'new_visitor_beds', 'capsule_purchases_kg', 'local_working_capital', 'local_working_replacement'):
            late[k] += f[k]
    late_phys = physical(data, dict(rows[0], population=0.0, new_residential_people=0.0, goods=0.0, source_goods=0.0, propellant_sales=0.0, compute=0.0, service_jobs=0.0, visitor_beds=0.0, visitor_nights_year=0.0, service_power_W=0.0, premium_power_W=0.0), late, measures=measures)
    seed_record['additional_preproduction_cost_USD2041'] = (late['external_capital'] + late['external_replacement'] + late['local_working_capital'] + late['local_working_replacement'] + late_phys['Earth_launch_USD2026'] + late_phys['operating_material_USD2026'] + late_phys['resource_provision_USD2026']) * I
    seed_record['lunar_raw_kg'] += late_phys['lunar_raw_kg_year']
    seed_record['source_feed_kg'] += late_phys['source_feed_kg_year']
    seed_record['cost']['additional_preproduction_cost_USD2041'] = seed_record['additional_preproduction_cost_USD2041']
    seed_record['total_resource_provision_USD2041'] = seed_record['cost']['resource_provision_USD2041'] + late_phys['resource_provision_USD2026'] * I
    seed_record['cost']['total_USD2041'] = sum((v for k, v in seed_record['cost'].items() if k != 'identified_service_hardware_USD2041'))
    return dict(seed=seed_record, annual=result, component_orders=[r for r in events if r['order_year'] <= end], convergence=dict(converged=converged, iterations=iteration + 1, relative_requirement_residual=residuals[-1], history=residuals), first_capacity_shortfall=next((r['year'] for r in result if not r['feasible']), None))
