"""Conditional U price indices; no currency conversion or inflation forecast.

Each experiment holds a dated portfolio's fixed-price financial path as the
counterfactual. It changes expenditure directed at a synthetic domestic final-use
basket and solves an explicit supply response. The index is a level difference,
not an annual rate, and does not feed back into the financing simulation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from . import strategies

HERE = strategies.HERE


def price_response(demand_shift, supply_elasticity, extra_supply=0.0):
    """Solve P*Q=1+d and Q=(1+a)*P**epsilon, with base P=Q=1."""
    if (not all(math.isfinite(x) for x in
                (demand_shift, supply_elasticity, extra_supply))
            or demand_shift <= -1 or extra_supply <= -1
            or supply_elasticity < 0):
        raise ValueError('Invalid expenditure or supply experiment')
    price = ((1 + demand_shift) / (1 + extra_supply)) ** (1 / (1 + supply_elasticity))
    quantity = (1 + extra_supply) * price ** supply_elasticity
    return dict(price_index=100 * price, quantity_index=100 * quantity,
                expenditure_index=100 * (1 + demand_shift),
                price_change_percent=100 * (price - 1))


def window(rows, horizon, years, failure_year=None):
    if failure_year is not None and failure_year <= horizon:
        return None
    chosen = [row for row in rows if horizon - years < row['year'] <= horizon]
    expected = list(range(horizon - years + 1, horizon + 1))
    return chosen if [row['year'] for row in chosen] == expected else None


def experiment(distributions, investment, basket, case):
    if basket <= 0 or distributions < 0 or investment < 0:
        raise ValueError('Positive basket and nonnegative funded expenditure required')
    if not 0 <= case['distribution_share'] <= 1 or not 0 <= case['investment_local_share'] <= 1:
        raise ValueError('Spending shares must be in [0,1]')
    if case['basket_multiplier'] <= 0:
        raise ValueError('Positive valuation multiplier required')
    distribution_demand = case['distribution_share'] * distributions
    investment_demand = case['investment_local_share'] * investment
    shift = (distribution_demand + investment_demand) / (basket * case['basket_multiplier'])
    return dict(demand_shift=shift,
                distribution_demand_USD2041_equivalent=distribution_demand,
                investment_demand_USD2041_equivalent=investment_demand,
                **price_response(shift, case['supply_elasticity'], case['extra_supply']))


def compare_window(rows, profile, config):
    domestic = profile['domestic']['USD2041']
    # Producer consumption includes civic-service and housing imputations.
    # Domestic capital formation excludes purchased Earth capital; the profile
    # has no expansion. No wage or intermediate-input total is added again.
    basket_per_scale = domestic['resident_consumption'] + domestic['domestic_capital_formation']
    baskets = [basket_per_scale * row['scale'] for row in rows]
    dividends = [row['dividend_USD2041'] for row in rows]
    investments = [row['industrial_orders_USD2041'] + row['service_orders_USD2041'] for row in rows]
    n = len(rows)
    avg = lambda values: math.fsum(values) / n
    basket, dividend, investment = avg(baskets), avg(dividends), avg(investments)
    cases = {}
    for name, case in config['cases'].items():
        result = experiment(dividend, investment, basket, case)
        peaks = [experiment(f, i, b, case) for f, i, b in zip(dividends, investments, baskets)]
        result['maximum_annual_price_index'] = max(row['price_index'] for row in peaks)
        result['maximum_annual_demand_shift'] = max(row['demand_shift'] for row in peaks)
        cases[name] = result
    # These thresholds concern the distribution-only channel, with no extra
    # supply shift. Values >1 mean even recycling all distributions is below
    # the selected index limit in this experiment.
    thresholds = {str(eps): basket * (config['index_limit'] ** (1 + eps) - 1) / dividend
                  if dividend > 0 else None for eps in (0, 1, 3)}
    imports = avg([row['operation']['core_imports_USD2041'] for row in rows])
    cash = avg([row['operation']['available_USD2041'] for row in rows])
    return dict(window_start=rows[0]['year'], window_end=rows[-1]['year'],
                mean_scale=avg([row['scale'] for row in rows]),
                domestic_final_use_basket_USD2041_year=basket,
                distributions_USD2041_year=dividend, investment_orders_USD2041_year=investment,
                distributions_over_basket=dividend / basket,
                investment_orders_over_basket=investment / basket,
                cases=cases, distribution_share_at_index_limit=thresholds,
                extra_supply_for_no_price_rise_at_one_percent=0.01 * dividend / basket,
                mean_available_cash_USD2041_year=cash,
                mean_identified_core_imports_USD2041_year=imports,
                cash_after_import_price_shock_USD2041_year=cash - config['identified_core_import_price_shock'] * imports,
                core_import_price_rise_at_zero_cash=cash / imports if imports > 0 else None)


def results():
    config = json.loads((HERE / 'prices_inputs.json').read_text())
    saved = json.loads((HERE / 'results/strategies.json').read_text())
    q, d, b, e, t = strategies.inputs()
    if saved['producer_sha256'] != strategies.run.digest(HERE / 'strategies.py'):
        raise ValueError('Strategy source changed without regenerated results')
    for path, digest in saved['dependencies'].items():
        if digest != strategies.run.digest(HERE / path):
            raise ValueError('Strategy dependency changed: ' + path)
    expected = {(row['strategy'], row['market'], row['horizon_elapsed_years']): row for row in saved['horizons']}
    output = []
    for name in q['strategies']:
        profile = strategies.profile(q, d, name)
        if profile['engineering']['growth_people_year'] != 0:
            raise ValueError('Final-use baseline already includes expansion')
        for market in q['markets']:
            simulation = strategies.simulate(q, d, b, e, t, name, market, max(config['horizons']))
            if simulation['summary'] != expected[name, market, max(config['horizons'])]:
                raise ValueError('Underlying financial path failed exact reproduction')
            for horizon in config['horizons']:
                failure = simulation['summary']['failure']
                selected = window(simulation['annual'], horizon, config['averaging_years'],
                                  failure['year'] if failure else None)
                row = dict(strategy=name, market=market, horizon=horizon)
                if selected is None:
                    row.update(status='unfunded_path_no_sustained_price_estimate', failure=failure)
                else:
                    row.update(status='conditional_price_level_experiment', **compare_window(selected, profile, config))
                output.append(row)
    dependencies = ('prices_inputs.json', 'strategies.py', 'strategy_inputs.json',
                    'results/strategies.json', 'inputs/strategy_profiles.json',
                    'staged.py', 'staged_inputs.json', 'run.py', 'returns_inputs.json')
    return dict(schema='terluna.research.price-pressure/1', evidence=config['evidence'],
                unit=config['unit'], boundary=config['boundary'],
                producer_sha256=strategies.run.digest(HERE / 'prices.py'),
                dependencies={p: strategies.run.digest(HERE / p) for p in dependencies},
                baseline_valuation='Producer domestic final-use imputation, including its illustrative 30 USD2041/U convention; this is not a market exchange rate or observed expenditure. Basket multipliers test 0.5 and 2 times that valuation.',
                reading_rule='Indices compare each dated portfolio with its own no-additional-spending/supply counterfactual. They are not annual inflation rates or an index path from year zero. Read assumptions beside every result. Full financial paths reproduce the saved strategy summaries exactly; pricing does not feed back into them.',
                cases=config['cases'], rows=output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=strategies.run.ROOT / 'research/runs/complex_returns/prices.json')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results(), indent=2, allow_nan=False) + '\n')
    print(args.output)
