"""One run of the vacuum buildout's economy: each case's orders, the industry that fills them and its financing.

usage: python3 -m research.studies.vacuum_economy.run [--case NAME ...] [--jobs N] [--output PATH]      (in the repository root)

A case names a scenario of demand, a population ramp and the recovery measures the industry applies. For each case
the run draws up the customers' orders year by year, sizes the industry that fills them with the cohort model of
provisioning/vacuum_industry, finances that industry with the financing model of the returns study, and writes one
compact record: the funding, the receipts and the returns at each horizon, the receipts and costs of each period with
the share of the costs paid to Earth, the material drawn from the Moon and the near-Earth sources, and what the mined
feed makes available beside what the system uses and sells of it.

The scenarios `low`, `central` and `high` are the demand of scenarios.json, judged at each year's date, each with the
launch price and the local manufacture that scenario states. For them mix.py decides from the system's own costs
which lines of sale are served, at what price and how far, and the record adds the account of each line.

The scenario `first_pass` is the demand of the staged studies of 12 to 14 November 2026, with their prices, audiences
and buyers' budgets. On the ramp of the same name it reproduces the figures of those studies, so that every later
change can be measured against them; on the programme's two ramps it shows what the schedule of residents alone changes.

A scenario case on a solved ramp (`resident_led`, `resident_led_beyond_plan`) takes its residents from the schedule
residents.py found for that scenario and saved in results/residents.json; the record adds how that schedule was found.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import importlib.util
import json
import math
import multiprocessing
from pathlib import Path
import time

from provisioning.vacuum_industry import industry, model
from research.studies.enclosed_returns import finance, markets, revised
from . import demand as scenarios, mix, programme as definition
from .inputs import HERE, INDUSTRY, RETURNS, ROOT, digest, engineering, sources
from . import inputs

SCHEMA = 'terluna.research.vacuum-economy/1'
RESIDENTS = HERE / 'results/residents.json'
FINANCING = ('equity', 'funded_credit', 'U_credit')
EARTH = ('capital_external_USD2041', 'operating_external_replacement_USD2041',
         'operating_external_logistics_materials_USD2041', 'operating_external_services_USD2041')


def dependencies():
    """The files the cases are computed from: the study's inputs, the list of cases and the solved schedules."""
    return inputs.dependencies((HERE / 'cases.json', RESIDENTS))


@lru_cache(maxsize=None)
def solved():
    """The schedules residents.py found, refused where any file it read has changed since."""
    product = json.loads(RESIDENTS.read_text())
    stale = [path for path, value in product['dependencies'].items() if not (ROOT / path).exists() or digest(ROOT / path) != value]
    if stale:
        raise ValueError('The solved schedules are out of date (run residents.py again): ' + ', '.join(stale))
    return product


def programme_for(case):
    """The programme a case runs on: the definition, with a solved ramp's residents where the case names one."""
    programme = sources()['programme']
    if programme['ramps'][case['ramp']]['kind'] != 'solved':
        return programme
    if case['scenario'] == 'first_pass':
        raise ValueError('A solved ramp belongs to a scenario of demand, not to the first pass')
    found = solved()['scenarios'][case['scenario']]['schedules'][case['ramp']]
    return definition.solved(programme, case['ramp'], found)


@lru_cache(maxsize=None)
def staged_demand(ramp):
    """The staged studies' order-book module, drawing its residents from one of the programme's ramps.

    That module takes its residents from a function of its own. For any ramp but the staged one, a private copy of
    the module is loaded and that one function replaced; the staged module and its figures stay as they are.
    """
    programme = sources()['programme']
    if ramp not in programme['ramps']:
        raise ValueError('Unknown population ramp: ' + ramp)
    if ramp == 'first_pass':
        return markets
    spec = importlib.util.spec_from_file_location('research.studies.enclosed_returns._markets_on_' + ramp, RETURNS / 'markets.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.population = lambda year, rules: definition.population(programme, ramp, year)
    return module


def orders(case):
    """The customers' orders and the residents of each year, for one case."""
    s = sources()
    if case['scenario'] != 'first_pass':
        raise ValueError('Unknown scenario of demand: ' + case['scenario'])
    rules, strategies, staged, services, _, _ = s['demand']
    book = staged_demand(case['ramp']).case(rules, strategies, staged, services, case['portfolio'], case['market'],
                                            case['growth'], rules['market_growth'][case['growth']])
    for row in book['annual']:
        planned = definition.population(s['programme'], case['ramp'], row['year'])
        if not math.isclose(row['population'], planned, rel_tol=1e-12):
            raise ValueError('The order book and the programme definition disagree on the residents of year %d' % row['year'])
    return book


def traffic(physical, coefficients, horizons):
    """The flights of the year of each horizon: launches from Earth and the pads they imply, and the trips of each route."""
    transport = coefficients['industrial']['transport']
    rows = {row['year']: row['physical'] for row in physical['annual']}
    out = {}
    for h in horizons:
        year = rows[h]
        launches = year['LEO_kg_year'] / transport['launch_payload_kg']
        out[str(h)] = dict(launched_kg_year=year['LEO_kg_year'], launches_year=launches,
                           launch_pads=launches / transport['launches_per_pad_year'],
                           Earth_leg_propellant_launched_kg_year=year['Earth_leg_propellant_launched_kg_year'],
                           ferry_trips_year=year['Earth']['trips'], lunar_landings_year=year['Moon']['trips'],
                           source_tug_trips_year=year['source']['trips'])
    return out


def workforce(physical, loops, coefficients, horizons):
    """The residents' jobs in the year of each horizon, and the working parts of the machinery they make."""
    working = model.local_working_parts(coefficients)
    staffed = loops.get('local_working_parts', {}).get('staffed_by_order_year')
    rows = {row['year']: row for row in physical['annual']}
    out = {}
    for h in horizons:
        part = staffed[str(h)] if staffed else model.at(working['staffed'], h)
        out[str(h)] = dict(jobs=rows[h]['technical_resident_jobs'], working_parts_jobs=rows[h]['working_parts_jobs'],
                           share_of_workforce_in_use=rows[h]['capacity_summary']['resident_workforce'],
                           working_parts_made_share=model.at(working['share'], h) * part)
    return out


def periods_account(book, physical, path, financed, money, periods):
    """Receipts by sector and costs by category in each period, undiscounted, with the share of costs paid to Earth."""
    account = revised.stratification(book, physical, path, financed, money, periods)['periods']
    out = []
    for period in account:
        values = period['values']['0.0']
        costs, receipts = values['costs_by_category'], values['revenue_by_sector']
        total = math.fsum(costs.values())
        # The founding deployment is one amount, with its purchases from Earth and its local work together.
        itemised = total > 0 and not costs['seed_deployment_and_allowance_USD2041']
        earth = math.fsum(costs[k] for k in EARTH) if itemised else None
        out.append(dict(first_year=period['first_year'], last_year=period['last_year'],
                        receipts_USD2041=math.fsum(receipts.values()), receipts_by_sector_USD2041=receipts,
                        costs_USD2041=total, costs_by_category_USD2041=costs,
                        paid_to_Earth_USD2041=earth, paid_to_Earth_share=earth / total if itemised else None,
                        distributions_USD2041=values['distributions'], additional_equity_USD2041=values['additional_equity']))
    return out


def decide(case):
    """The orders of a scenario case: prices from the system's own costs, and the lines its residents can staff."""
    s = sources()
    hotel = s['money']['hotel_staff_per_bed'] * s['money']['hotel_staff_USD2041_year']
    decision = mix.Decision(programme_for(case), case['ramp'], s['scenarios'], case['scenario'], engineering(case), s['recovery'],
                            case['recovery'], s['periods'], hotel, log=lambda *a: print('  ', case['name'], *a, flush=True),
                            pricing=case.get('pricing', 'seller'), allowance=s['money']['capital_cost_allowance_fraction'])
    return decision.decide()


def lines_record(decided, horizons):
    """Each line of sale as served: its worth, its account by period, and its price and quantity at each horizon."""
    rows = {row['year']: row for row in decided['book']['annual']}
    out = {}
    for line in scenarios.LINES:
        account = decided['final'].get(line)
        at = {}
        for h in horizons:
            detail = rows[h]['lines'][line]
            at[str(h)] = dict(quantity=detail['quantity'], price_USD2041=detail['price_USD2041'], receipts_USD2041=detail['receipts_USD2041'],
                              share_of_demand_served=decided['scale'][line][h])
        out[line] = dict(unit=rows[horizons[0]]['lines'][line]['unit'], served=account is not None,
                         dropped_as_unprofitable=line in decided['dropped'],
                         first_year_served=(decided['starts'].get(line, min(decided['scale'][line])) if account is not None else None),
                         worth_USD2041=account['worth_USD2041'] if account else None,
                         worth_at_first_prices_USD2041=decided['first'][line]['worth_USD2041'] if line in decided['first'] else None,
                         periods=account['periods'] if account else [], at_horizon=at)
    return dict(comparison_rate=mix.RATE, order_of_worth_per_job=decided['order'], lines=out, industry_solutions=decided['solves'])


def resources(book, physical, loops, coefficients, horizons):
    """What the mined feed makes available to each horizon, what the system uses and sells of it, and what lies idle."""
    source = coefficients['industrial']['source']
    metals = sources()['scenarios']['lines']['metals']
    out = {}
    for h in horizons:
        drawn = revised.cumulative_demand(physical, loops, h)
        mined = drawn['feed_kg']

        def sold(key):
            return math.fsum(row.get(key, 0.0) for row in book['annual'] if row['year'] <= h)

        def stream(name, used):
            recoverable = mined * source[name + '_grade'] * source[name + '_recovery']
            return dict(recoverable_kg=recoverable, used_kg=used, idle_kg=max(0.0, recoverable - used))

        # In a scenario the goods the sources return to Earth are precious metals; in the first pass they are metal products.
        precious = sold('source_goods') if 'lines' in book['annual'][0] else 0.0
        products = sold('space_source_goods') + sold('source_goods') - precious
        contained = mined * (metals['contained_PGM_per_kg_feed'] + metals['contained_gold_per_kg_feed'])
        out[str(h)] = dict(asteroid_feed_kg=mined, lunar_raw_kg=drawn['lunar_raw_kg'],
                           water=dict(stream('water', drawn['water_kg_bounds'][0]), sold_as_propellant_kg=sold('propellant_sales')),
                           nitrogen=stream('nitrogen', drawn['nitrogen_kg']),
                           metal=dict(stream('metal', drawn['virgin_metal_kg']), sold_as_products_kg=products),
                           precious_metals=dict(contained_kg=contained, sold_kg=precious, left_in_tailings_kg=contained - precious),
                           held_back_kg=0.0)
    return out


def residents_record(case):
    """How the schedule of a solved ramp was found: its residents at each horizon, the rate and the binding condition
    of each period, the year it reaches the planned population, and whether the opening settlement is sustained."""
    entry = solved()['scenarios'][case['scenario']]
    programme = programme_for(case)
    trial = entry['limited_trial'] if case['ramp'] == 'resident_led' and entry['limited_trial'] else entry['final_trial']
    return dict(residents_at={str(h): definition.population(programme, case['ramp'], h) for h in programme['horizons']},
                planned_population_reached_year=entry['planned_population_reached_year'],
                opening_settlement_sustained=entry['opening_settlement_sustained'],
                periods=[dict(first_year=p['first_year'], last_year=p['last_year'], growth_rate=p['rate'],
                              lowest_unsustained_rate=p.get('lowest_unsustained_rate'), binding=p['binding']) for p in entry['periods']],
                search=dict(trials=entry['trials'], margins_of_the_schedule=trial.get('margins'), sustained=trial['sustained']))


def checkpoint(measure):
    """The headline measures of one financing arrangement at one horizon."""
    return dict(horizon=measure['horizon'], financing=measure['mode'], population=measure['population'],
                receipts_USD2041_year=measure['annual_gross'], operating_cash_USD2041_year=measure['annual_operating_available'],
                cumulative_receipts_USD2041=measure['cumulative_gross_receipts'], payouts_USD2041=measure['owner_payouts'],
                first_payout_year=measure['first_payout_year'], initial_funding_USD2041=measure['initial_full_funding'],
                additional_equity_USD2041=measure['additional_equity_required'],
                all_contributions_USD2041=measure['all_funder_contributions'], real_IRR=measure['real_IRR_roundoff_screened'],
                funder_cash_NPV_USD2041={rate: v['sponsor_distributed_cash_NPV'] for rate, v in measure['discounted'].items()},
                system_cash_NPV_USD2041={rate: v['unlevered_cash_NPV'] for rate, v in measure['discounted'].items()},
                socioeconomic_NPV_USD2041={rate: v['full_price_socioeconomic_NPV'] for rate, v in measure['discounted'].items()})


def solve(case):
    """Orders, industry and financing for one case, as a compact record."""
    s = sources()
    money = dict(s['money'], horizons=s['programme']['horizons'])
    coefficients = engineering(case)
    if case['scenario'] == 'first_pass':
        decided = None
        book = orders(case)
        physical, loops = industry.build(coefficients, s['recovery'], book, case['recovery'])
    else:
        decided = decide(case)
        book, physical, loops = decided['book'], decided['physical'], decided['loops']
    record = dict(case, convergence=physical['convergence'], first_capacity_shortfall=physical['first_capacity_shortfall'])
    if not revised.valid(physical):
        return dict(record, status='physical_requirements_unresolved')
    path = finance.plan(book, physical, money)
    equity = finance.finance(path, money, 'equity')
    financed = {'equity': equity, **{mode: finance.finance(path, money, mode, equity) for mode in FINANCING[1:]}}
    drawn = {}
    for horizon in money['horizons']:
        demand = revised.cumulative_demand(physical, loops, horizon)
        drawn[str(horizon)] = dict(asteroid_feed_kg=demand['feed_kg'], lunar_raw_kg=demand['lunar_raw_kg'],
                                   water_kg=demand['water_kg_bounds'][0], nitrogen_kg=demand['nitrogen_kg'],
                                   virgin_metal_kg=demand['virgin_metal_kg'])
    if decided:
        record = dict(record, sales=lines_record(decided, money['horizons']))
    if s['programme']['ramps'][case['ramp']]['kind'] == 'solved':
        record = dict(record, residents=residents_record(case))
    return dict(record, status='conditional_account', resources=resources(book, physical, loops, coefficients, money['horizons']),
                water_against_pool=inputs.water_drawn(physical, loops, coefficients, money['horizons']),
                opening=dict(founding_cost_USD2041=path['seed_hardware_and_carry'], working_cash_USD2041=path['initial_working'],
                             reserve_USD2041=equity['initial_reserve']),
                checkpoints=[checkpoint(revised.measures(path, financed[mode], money, horizon))
                             for mode in FINANCING for horizon in money['horizons']],
                periods=periods_account(book, physical, path, financed['U_credit'], money, s['periods']),
                material_drawn=drawn, traffic=traffic(physical, coefficients, money['horizons']),
                workforce=workforce(physical, loops, coefficients, money['horizons']))


def cases(names=None):
    listed = json.loads((HERE / 'cases.json').read_text())['cases']
    known = {c['name']: c for c in listed}
    if len(known) != len(listed):
        raise ValueError('Two cases share a name')
    if names:
        missing = [n for n in names if n not in known]
        if missing:
            raise ValueError('Unknown case: ' + ', '.join(missing))
        return [known[n] for n in names]
    return listed


def timed(case):
    start = time.monotonic()
    record = solve(case)
    print(case['name'], record['status'], round(time.monotonic() - start, 1), 'seconds', flush=True)
    return record


def product(selected, jobs=1):
    if jobs > 1:
        with multiprocessing.Pool(jobs) as pool:
            records = pool.map(timed, selected, chunksize=1)
    else:
        records = [timed(case) for case in selected]
    return dict(schema=SCHEMA, producer=dict(path=str(Path(__file__).resolve().relative_to(ROOT)), sha256=digest(__file__)),
                dependencies=dependencies(),
                evidence='Conditional accounts. The industry is a numerically closed sizing of equipment, freight and '
                         'power for the stated orders; the orders are scenarios of demand, and the financing follows '
                         'stated rules. Sources, plant, markets and prices await qualification.',
                units='Years of the programme, residents, kilograms, constant USD2041.',
                reading_rule='Each case stands alone: one scenario of demand on one population ramp with one set of '
                             'recovery measures. A checkpoint gives one financing arrangement at one horizon. The periods '
                             'use the financing in U credit and are undiscounted; costs paid to Earth are imported '
                             'machinery and its replacement, launch and imported materials, and outside services; '
                             'the period of the founding deployment gives none, because that deployment is one amount. '
                             'Material drawn is cumulative from year zero. Traffic is for the year of each horizon: the mass launched '
                             'from Earth, cargo and any propellant for the Earth leg together, the launches that carry it at the '
                             'stated payload and the pads those launches need at the stated launches a pad a year, and the round '
                             'trips of the ferries, the lunar landers and the tugs to the sources. Workforce gives the jobs of the same '
                             'year, the share of the available workforce they fill, the jobs that make working parts of machinery '
                             'and the share of those parts made in place for orders of that year. Resources give, cumulatively to each '
                             'horizon, the feed mined at the sources, what its water, nitrogen and metal yield at the stated grades '
                             'and recoveries, what the system uses of each, what it sells, and what lies idle; nothing is held back '
                             'while the protection gates reserve no quantity. A scenario case adds its sales: for each line whether '
                             'it is served and from which year or was dropped as costing more than it brings, its worth at the '
                             'comparison rate, its receipts, cost, quantity and mean jobs by period against the same orders without '
                             'the line, and its quantity, price, receipts and the share of its demand served in the year of each '
                             'horizon. A line\'s cost counts the financing\'s allowance on capital and, where residents cannot staff '
                             'all the working parts the industry is able to make, the parts its jobs would have made. Water against '
                             'the pool gives the water contained in the feed mined to each horizon (the feed times its stated grade, '
                             'before recovery) beside the water of the source survey\'s two pools, the wider near-Earth screen '
                             'within which the programme\'s sources lie and the broader near-Earth population that transport '
                             'reaching beyond that screen would open, at their 5th, 50th and 95th percentiles, with the year the '
                             'draw passes each. A case on a solved ramp adds how its '
                             'schedule was found: the residents at each horizon, the growth rate of each period with the lowest '
                             'rate found unsustained and the condition that failed there, the year the planned population is '
                             'reached, and the margins of the schedule as residents.py measured them while the surplus was kept '
                             'for growth; its checkpoints and periods finance the same schedule under the ordinary rules.',
                cases=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', action='append', help='a case named in cases.json; repeat for several (default: all)')
    parser.add_argument('--output', type=Path, default=HERE / 'results/economy.json')
    parser.add_argument('--jobs', type=int, default=1, help='cases solved at once (default: one)')
    args = parser.parse_args()
    result = product(cases(args.case), args.jobs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1, sort_keys=True, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
