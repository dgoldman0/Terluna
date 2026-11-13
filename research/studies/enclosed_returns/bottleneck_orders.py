"""Controlled demand and commercial-rollout experiments on the saved baseline.

Orders remain conditional requests. The provisioning owner must solve every
changed request before the financial consumer admits its revenues.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from . import markets

HERE = Path(__file__).resolve().parent


def demand_variant(q, b, market, multiplier):
    """Expand budgets, funded quantities and participation; keep audience fixed."""
    q, b = deepcopy(q), deepcopy(b)
    for key in (
        'industrial_budget_USD2041_year', 'mission_budget_USD2041_year',
        'household_budget_USD2041_year', 'institutional_services_budget_USD2041_year',
        'premium_demand_kg_year', 'compute_contracted_capacity_W',
        'funded_compatible_missions_year',
    ):
        q['markets'][market][key] *= multiplier
    scenario = b['scenarios'][market]
    for key in ('paid_fraction', 'licenses', 'research_contracts', 'professional_clients',
                'advertising_and_rights_USD2041_fan_year'):
        scenario[key] *= multiplier
    if scenario['paid_fraction'] > 1:
        raise ValueError('Paid participation exceeds the fixed audience')
    return q, b


def scheduled_case(cfg, q, d, b, name, market, growth_name, rate, max_orders):
    """Retain baseline scheduling, allowing a declared annual package-order limit."""
    kits = 1
    pending = []
    paid = licenses = last_stock = 0.
    rows = []
    for year in range(cfg['first_productive_year'], cfg['last_year'] + 1):
        age = year - cfg['first_productive_year']
        kits += sum(n for due, n in pending if due == year)
        pending = [(y, n) for y, n in pending if y > year]
        maturation = min(1., b['service_demand_initial_fraction'] *
                         (1 + b['service_demand_annual_growth'])**min(age, 100))
        factor = (1 + rate)**min(age, cfg['completion_year'] - cfg['first_productive_year'])
        scenario = b['scenarios'][market]
        wanted = math.ceil(max(scenario['research_contracts'] / b['service_kit']['research_contracts_year'],
                               scenario['professional_clients'] / b['service_kit']['professional_clients_year']) * maturation * factor)
        staff_limit = max(1, math.floor((cfg['service_staff_fraction'] * markets.population(year, cfg)
                                       - b['cultural_platform']['staff']) / b['service_kit']['staff']))
        gap = max(0, min(wanted, staff_limit) - kits - sum(n for _, n in pending))
        order = min(max_orders, gap) if year <= cfg['completion_year'] - cfg['commissioning_years'] else 0
        if order:
            pending.append((year + cfg['commissioning_years'], order))
        row = markets.annual_market(cfg, q, d, b, name, market, rate, year, kits, paid, licenses)
        paid, licenses = row['paid_accounts'], row['licenses']
        stock = row['service_installed_capital_USD2041']
        row.update(year=year, population=markets.population(year, cfg),
                   commissioned_people=0. if year == cfg['first_productive_year'] else markets.population(year, cfg) - markets.population(year-1, cfg),
                   new_residential_people=markets.population(year+2, cfg) - markets.population(year+1, cfg),
                   service_new_capital_USD2041=max(0., stock-last_stock), service_kit_orders=order)
        last_stock = stock
        rows.append(row)
    return dict(name=f'{name}:{market}:{growth_name}', portfolio=name, market=market,
                market_growth=rate, market_growth_name=growth_name, annual=rows)


def product():
    cfg, q, d, b, _, _ = markets.inputs()
    spec = json.loads((HERE/'bottleneck_inputs.json').read_text())
    for path, identity in spec['baseline_dependencies'].items():
        if markets.digest(HERE/path) != identity:
            raise ValueError('Changed diagnostic baseline: ' + path)
    cases = []
    for portfolio in spec['experimental_portfolios']:
        for variant, changes in spec['experiments'].items():
            cc, qq, bb = deepcopy(cfg), deepcopy(q), deepcopy(b)
            cc.update(changes.get('market_config', {}))
            qq, bb = demand_variant(qq, bb, spec['market'], changes.get('demand_multiplier', 1.))
            case = scheduled_case(cc, qq, d, bb, portfolio, spec['market'], spec['growth'],
                                  cfg['market_growth'][spec['growth']], changes.get('service_orders_per_year', 1))
            base_name = case['name']
            case['name'] += ':' + variant
            case['experiment'] = dict(baseline_name=base_name, variant=variant,
                                      engineering_overrides=changes.get('engineering_overrides', {}))
            cases.append(case)
    return dict(schema='terluna.humanities.enclosed-bottleneck-orders/1',
                producer=dict(path='research/studies/enclosed_returns/bottleneck_orders.py', sha256=markets.digest(__file__)),
                dependencies={p: markets.digest(HERE/p) for p in ('bottleneck_inputs.json', 'markets.py', 'enclosed_inputs.json', *cfg['inherited_inputs'])},
                evidence=spec['evidence'], rules=spec['rules'], experiments=spec['experiments'], cases=cases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE.parents[1]/'runs/enclosed_returns/bottleneck_orders.json')
    args = parser.parse_args()
    result = product()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'), allow_nan=False) + '\n')
    print('Diagnostic order cases:', len(result['cases']))
    print(args.output)


if __name__ == '__main__':
    main()
