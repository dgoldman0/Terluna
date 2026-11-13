"""Conditional U money-stock transmission, domestic prices and liquidity funding.

The external quote is controlled, not forecast. Each dated ten-year window is a
separate policy experiment. Fixed-price financial paths remain the counterfactual.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from research.studies.complex_returns import strategies

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def reserve_allocation(available, distribution, share):
    if not 0 <= share <= 1 or distribution < 0:
        raise ValueError('Invalid reserve allocation')
    contribution = min(share * max(available, 0), distribution)
    return contribution, distribution - contribution


def money_price(basket, planned_spending, extra_outside_balances, case, extra_supply):
    """Stock identities at reference quotes and fixed domestic turnover."""
    if (not all(math.isfinite(v) for v in (basket, planned_spending, extra_outside_balances,
                                          extra_supply, *case.values()))
            or basket <= 0 or planned_spending <= 0 or extra_supply <= -1
            or case['domestic_holding_days'] <= 0 or case['supply_elasticity'] < 0
            or case['issuance_accommodation'] < 0
            or not 0 <= case['dealer_inventory_absorption'] <= 1):
        raise ValueError('Invalid monetary experiment')
    k = case['domestic_holding_days'] / 365
    required_local = k * planned_spending
    dealer_transfer = case['dealer_inventory_absorption'] * extra_outside_balances
    outside = extra_outside_balances - dealer_transfer
    issued = case['issuance_accommodation'] * outside
    local_change = issued - outside
    local_stock = required_local + local_change
    result = dict(planned_spending_value=planned_spending, required_local_balances_value=required_local,
                  extra_customer_balances_value=extra_outside_balances,
                  outside_dealer_inventory_reduction_value=dealer_transfer,
                  extra_net_outside_balances_value=outside, issued_for_outside_balances_value=issued,
                  local_balance_change_value=local_change, remaining_local_balances_value=local_stock,
                  monetary_balance_gap_fraction=local_change / required_local,
                  outside_acquisition_not_accommodated_by_issuance_value=outside-issued)
    if local_stock <= 0:
        return dict(result, status='infeasible_domestic_balance_exhaustion')
    spending = local_stock / k
    price = (spending / basket / (1 + extra_supply)) ** (1 / (1 + case['supply_elasticity']))
    quantity = (1 + extra_supply) * price ** case['supply_elasticity']
    return dict(result, status='conditional_reference_quote', price_level=price,
                quantity_level=quantity, actual_spending_value=spending)


def run_window(annual, profile, horizon, config, case):
    n = config['window_years']
    rows = [row for row in annual if horizon-n <= row['year'] <= horizon]
    if [row['year'] for row in rows] != list(range(horizon-n, horizon+1)):
        raise ValueError('Complete opening observation and transition required')
    domestic = profile['domestic']['USD2041']
    base = (domestic['resident_consumption'] + domestic['domestic_capital_formation']) * case['basket_multiplier']
    holding = case['customer_holding_days'] / 365
    opening_foreign = config['opening_U_settlement_share'] * rows[0]['operation']['gross_USD2041'] * holding
    dealer_opening = case['opening_dealer_inventory_multiple_of_customer_balances'] * opening_foreign
    trajectory = []
    previous = 1.
    contributed = 0.
    available_sum = 0.
    distribution_sum = 0.
    for j, row in enumerate(rows):
        fraction = j / n
        op = row['operation']
        reserve, payout = reserve_allocation(op['available_USD2041'], row['dividend_USD2041'], case['reserve_cash_share']) if j else (0., row['dividend_USD2041'])
        if j:
            contributed += reserve
            available_sum += op['available_USD2041']
            distribution_sum += row['dividend_USD2041']
        basket = base * row['scale']
        local_payout = fraction * case['distribution_spending_share'] * payout
        orders = row['industrial_orders_USD2041'] + row['service_orders_USD2041']
        local_orders = fraction * case['investment_local_share'] * orders
        share = config['opening_U_settlement_share'] + fraction * (config['closing_U_settlement_share']-config['opening_U_settlement_share'])
        foreign = share * op['gross_USD2041'] * holding
        result = money_price(basket, basket+local_payout+local_orders, foreign-opening_foreign,
                             case, fraction*case['extra_supply'])
        dealer_remaining = dealer_opening-result['outside_dealer_inventory_reduction_value']
        result['outside_dealer_remaining_inventory_value'] = dealer_remaining
        if dealer_remaining < -1e-8:
            result['status'] = 'infeasible_outside_dealer_inventory'
            for key in ('price_level','quantity_level','actual_spending_value'):
                result.pop(key,None)
        result.update(year=row['year'], basket_value=basket, gross_receipts_value=op['gross_USD2041'],
                      reserve_contribution_value=reserve, payout_after_reserves_value=payout,
                      local_payout_spending_value=local_payout, local_order_spending_value=local_orders,
                      customer_balances_value=foreign, settlement_share=share)
        if 'price_level' in result:
            result['annual_price_change_percent'] = 100*(result['price_level']/previous-1) if previous is not None else None
            previous = result['price_level']
        else:
            previous = None
        trajectory.append(result)
    end = trajectory[-1]
    valid = all('price_level' in r for r in trajectory)
    summary = dict(status='conditional_reference_quote' if valid else 'infeasible_balance_transfer',
                   opening_outside_dealer_inventory_value=dealer_opening,
                   closing=end, reserve_contributions_value=contributed,
                   original_distributions_value=distribution_sum,
                   reduced_distributions_value=distribution_sum-contributed,
                   mean_available_cash_value=available_sum/n)
    if valid:
        changes = [r['annual_price_change_percent'] for r in trajectory[1:]]
        summary.update(price_level_change_percent=100*(end['price_level']-1),
                       geometric_annual_change_percent=100*(end['price_level']**(1/n)-1),
                       smallest_annual_change_percent=min(changes), largest_annual_change_percent=max(changes))
    return summary, trajectory


def import_cost_stresses(config):
    rows = []
    for e in config['U_external_value_factors']:
        for weight in config['basket_external_cost_shares']:
            for beta in config['pass_through_fractions']:
                rows.append(dict(U_external_value_factor=e, basket_external_cost_share=weight,
                                 pass_through_fraction=beta, imported_cost_change_percent=100*(1/e-1),
                                 domestic_basket_change_percent=100*weight*beta*(1/e-1)))
    return rows


def results():
    config = json.loads((HERE/'pressure_inputs.json').read_text())
    saved_path = strategies.HERE/'results/strategies.json'
    digest = strategies.run.digest
    if digest(saved_path) != config['strategy_result_sha256']:
        raise ValueError('Pinned financial product changed')
    saved = json.loads(saved_path.read_text())
    if saved['producer_sha256'] != digest(strategies.HERE/'strategies.py'):
        raise ValueError('Financial producer changed')
    for path, sha in saved['dependencies'].items():
        if digest(strategies.HERE/path) != sha:
            raise ValueError('Financial dependency changed: '+path)
    expected = {(r['strategy'],r['market'],r['horizon_elapsed_years']):r for r in saved['horizons']}
    q,d,b,e,t = strategies.inputs()
    output = []
    for name in q['strategies']:
        profile = strategies.profile(q,d,name)
        for market in q['markets']:
            simulation = strategies.simulate(q,d,b,e,t,name,market,max(config['horizons']))
            summary = simulation['summary']
            if summary != expected[name,market,max(config['horizons'])]:
                raise ValueError('Financial path failed exact reproduction')
            for horizon in config['horizons']:
                row = dict(strategy=name, market=market, horizon=horizon)
                failure = summary['failure']
                if failure and failure['year'] <= horizon:
                    output.append(dict(row,status='unfunded_path',failure=failure))
                    continue
                cases, paths = {}, {}
                for key, case in config['cases'].items():
                    cases[key], path = run_window(simulation['annual'],profile,horizon,config,case)
                    if key in ('supply_expansion','issuance_shortfall','issuance_excess'):
                        paths[key] = path
                final = next(r for r in simulation['annual'] if r['year']==horizon)
                op = final['operation']
                core = op['core_imports_USD2041']
                services = sum(op['revenue_by_sector_USD2041'][k] for k in ('sports','entertainment','licenses','research','professional'))
                initial_reserve = config['initial_reserve_fraction_of_original_seed']*summary['external_seed_USD2041']
                closing = cases['supply_expansion']['closing']
                stress = config['reversal_stress']
                retained = stress['export_receipt_fraction']*stress['new_customer_holding_days']/config['cases']['supply_expansion']['customer_holding_days']
                reversal = closing['customer_balances_value']*(1-retained)
                shock_loss = (1-stress['export_receipt_fraction'])*op['gross_USD2041']
                shock_cash = op['available_USD2041']-shock_loss
                shock_payout = max(0.,final['dividend_USD2041']-shock_loss)
                shock_reserve, shock_payout_after = reserve_allocation(shock_cash,shock_payout,
                    config['cases']['supply_expansion']['reserve_cash_share'])
                funding_gap = max(0.,shock_loss-final['dividend_USD2041'])
                reversal_cases = {}
                for retirement in stress['retirement_fractions']:
                    case = dict(config['cases']['supply_expansion'],issuance_accommodation=retirement)
                    if funding_gap > 0:
                        reversal_cases[str(retirement)] = dict(status='operating_and_order_funding_gap',gap_value=funding_gap)
                        continue
                    target = closing['basket_value']+case['distribution_spending_share']*shock_payout_after
                    r = money_price(closing['basket_value'],target,-reversal,case,case['extra_supply'])
                    r['price_change_from_pre_shock_percent'] = 100*(r['price_level']/closing['price_level']-1)
                    r['issuer_external_asset_use_value'] = retirement*reversal
                    r['private_external_asset_use_value'] = (1-retirement)*reversal
                    r['automatic_local_balance_supply_change_value'] = case['domestic_holding_days']/365*(target-closing['planned_spending_value'])
                    reversal_cases[str(retirement)] = r
                output.append(dict(row,status='conditional_monetary_and_price_experiment',
                    original_seed_value=summary['external_seed_USD2041'], additional_initial_reserve_capital_value=initial_reserve,
                    annual_core_external_purchases_value=core, annual_total_operating_cost_value=op['total_cost_USD2041'],
                    worldwide_service_receipt_fraction=services/op['gross_USD2041'],
                    cases=cases, paths=paths,
                    reversal=dict(customer_balance_release_value=reversal,
                        available_cash_after_receipt_shock_before_cost_adjustment_value=shock_cash,
                        post_shock_reserve_contribution_value=shock_reserve, post_shock_payout_value=shock_payout_after,
                        funding_gap_value=funding_gap, retirement_cases=reversal_cases,
                        initial_reserve_only_surplus_after_isolated_release_value=initial_reserve-reversal,
                        initial_plus_window_contributions_value=initial_reserve+cases['supply_expansion']['reserve_contributions_value'],
                        reserve_surplus_after_isolated_release_value=initial_reserve+cases['supply_expansion']['reserve_contributions_value']-reversal,
                        warning='This tests only the specified customer-balance release at the reference quote. It is not total liability coverage, an adopted convertibility promise or a rescue of a failed financial path.')))
    dependencies = [HERE/'pressure_inputs.json', HERE/'trade_transition_inputs.json',
                    HERE/'results/trade_transition.json', saved_path, strategies.HERE/'strategies.py',
                    strategies.HERE/'prices.py', *[strategies.HERE/p for p in saved['dependencies']]]
    return dict(schema='terluna.research.u-pressure/1', evidence=config['evidence'], boundary=config['boundary'],
        valuation='All value fields are source USD2041 equivalents at an unspecified common reference quote, never nominal U or a selected U/USD conversion.',
        reading_rule='Each horizon is a separate ten-year transition experiment, not repeated adoption along one history. Prices assume the controlled reference quote and fixed turnover. Annual changes are relative to a constant-price counterfactual, not unconditional forecasts. Floating FX response is directional only. Reserve construction reallocates payouts and leaves baseline orders unchanged; feedback into orders, credit and exports is outside this experiment.',
        producer_sha256=digest(HERE/'pressure.py'), dependencies={str(p.relative_to(ROOT)):digest(p) for p in dependencies},
        cases=config['cases'], import_cost_stresses=import_cost_stresses(config['import_stress']), rows=output)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/monetary_system/pressure.json')
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(args.output)
