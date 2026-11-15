"""One run of the vacuum buildout's economy: each case's orders, the industry that fills them and its financing.

usage: python3 -m research.studies.vacuum_economy.run [--case NAME ...] [--output PATH]      (in the repository root)

A case names a scenario of demand, a population ramp and the recovery measures the industry applies. For each case
the run draws up the customers' orders year by year, sizes the industry that fills them with the cohort model of
provisioning/vacuum_industry, finances that industry with the financing model of the returns study, and writes one
compact record: the funding, the receipts and the returns at each horizon, the receipts and costs of each period with
the share of the costs paid to Earth, and the material drawn from the Moon and the near-Earth sources.

The scenario `first_pass` is the demand of the staged studies of 12 to 14 November 2026, with their prices, audiences
and buyers' budgets. On the ramp of the same name it reproduces the figures of those studies, so that every later
change can be measured against them; on the programme's two ramps it shows what the schedule of residents alone changes.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

from provisioning.vacuum_industry import industry
from research.studies.enclosed_returns import finance, markets, revised
from . import programme as definition

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.vacuum-economy/1'
INDUSTRY = ROOT / 'provisioning/vacuum_industry'
RETURNS = ROOT / 'research/studies/enclosed_returns'
FINANCING = ('equity', 'funded_credit', 'U_credit')
EARTH = ('capital_external_USD2041', 'operating_external_replacement_USD2041',
         'operating_external_logistics_materials_USD2041', 'operating_external_services_USD2041')
# The files a first-pass case is computed from, beside those the order book and the financing name themselves.
READ = (HERE / 'programme.json', HERE / 'programme.py', HERE / 'cases.json',
        INDUSTRY / 'model.py', INDUSTRY / 'recovery.py', INDUSTRY / 'industry.py',
        INDUSTRY / 'engineering.json', INDUSTRY / 'recovery.json',
        RETURNS / 'markets.py', RETURNS / 'enclosed_inputs.json', RETURNS / 'finance.py', RETURNS / 'finance_inputs.json',
        RETURNS / 'revised.py', RETURNS / 'revised_inputs.json', RETURNS / 'bottlenecks.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def sources():
    """The programme, the engineering coefficients, the recovery measures, the demand inputs and the financing rules."""
    programme = definition.load()
    demand = markets.inputs()                                  # checks the digests of the inputs it inherits
    money = json.loads((RETURNS / 'finance_inputs.json').read_text())
    for name, pinned in money['inherited_dependencies'].items():
        if digest(RETURNS / name) != pinned:
            raise ValueError('Changed financing dependency: ' + name)
    periods = json.loads((RETURNS / 'revised_inputs.json').read_text())['periods']
    calendar, people, rules = programme['calendar'], programme['population'], demand[0]
    stated = dict(first_productive_year=calendar['first_productive_year'], completion_year=calendar['completion_year'],
                  last_year=calendar['last_year'], seed_population=people['opening'],
                  founding_population=people['founding_system'], cap_population=people['planned'])
    for key, value in stated.items():
        if rules[key] != value:
            raise ValueError('The demand inputs and the programme definition disagree on ' + key)
    if money['horizons'] != programme['horizons'] or money['rates'] != programme['discount_rates']:
        raise ValueError('The financing rules and the programme definition disagree on horizons or discount rates')
    if periods[0][0] != 0 or periods[-1][1] != calendar['last_year']:
        raise ValueError('The reporting periods do not span the programme')
    return dict(programme=programme, demand=demand, money=money, periods=periods,
                engineering=json.loads((INDUSTRY / 'engineering.json').read_text()),
                recovery=json.loads((INDUSTRY / 'recovery.json').read_text()))


def dependencies():
    s = sources()
    paths = {Path(p).resolve() for p in READ}
    paths |= {(RETURNS / name).resolve() for name in s['demand'][0]['inherited_inputs']}
    paths |= {(RETURNS / name).resolve() for name in s['money']['inherited_dependencies']}
    return {str(p.relative_to(ROOT)): digest(p) for p in sorted(paths)}


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
    book = orders(case)
    physical, loops = industry.build(s['engineering'], s['recovery'], book, case['recovery'])
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
    return dict(record, status='conditional_account',
                opening=dict(founding_cost_USD2041=path['seed_hardware_and_carry'], working_cash_USD2041=path['initial_working'],
                             reserve_USD2041=equity['initial_reserve']),
                checkpoints=[checkpoint(revised.measures(path, financed[mode], money, horizon))
                             for mode in FINANCING for horizon in money['horizons']],
                periods=periods_account(book, physical, path, financed['U_credit'], money, s['periods']),
                material_drawn=drawn)


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


def product(selected):
    records = []
    for case in selected:
        start = time.monotonic()
        records.append(solve(case))
        print(case['name'], records[-1]['status'], round(time.monotonic() - start, 1), 'seconds', flush=True)
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
                             'Material drawn is cumulative from year zero.',
                cases=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', action='append', help='a case named in cases.json; repeat for several (default: all)')
    parser.add_argument('--output', type=Path, default=HERE / 'results/economy.json')
    args = parser.parse_args()
    result = product(cases(args.case))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=1, sort_keys=True, allow_nan=False) + '\n')
    print(args.output)


if __name__ == '__main__':
    main()
