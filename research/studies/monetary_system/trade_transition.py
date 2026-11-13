"""Fixed-stock U settlement, constituent baskets and funded liquidity.

The existing L/S ledger examples remain unchanged. This consumer uses the
generic ledger with new actors and U_nom contracts, and reads financial
portfolios through their pinned data product. It does not solve an FX market.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from .ledger import Ledger, claim

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs(path=HERE / 'trade_transition_inputs.json'):
    config = json.loads(Path(path).read_text())
    if config['schema'] != 'terluna.research.trade-transition-inputs/1':
        raise ValueError('Unexpected trade-transition input schema')
    return config


def portfolio_product(config, root=HERE):
    source = config['portfolio_source']
    path = Path(root) / source['path']
    if digest(path) != source['sha256']:
        raise ValueError('Unreviewed portfolio product')
    product = json.loads(path.read_text())
    if product['schema'] != source['schema']:
        raise ValueError('Unexpected portfolio schema')
    return product


def capture(book):
    """Report the unchanged base and deposits held by each transacting actor."""
    snap = book.snapshot()
    base = book.amount(claim('zone_issuer', 'zone_bank', 'base_reserve', 'U_nom'))
    customer = book.amount(claim('zone_bank', 'earth_customer', 'deposit', 'U_nom'))
    exporter = book.amount(claim('zone_bank', 'zone_exporter', 'deposit', 'U_nom'))
    dealer = book.amount(claim('zone_bank', 'earth_dealer', 'deposit', 'U_nom'))
    return dict(
        base_U_nom=base, bank_deposits_U_nom=customer + exporter + dealer,
        customer_U_nom=customer, exporter_U_nom=exporter, dealer_U_nom=dealer,
        dealer_external_assets_S=book.amount(claim('external_issuer', 'earth_dealer', 'settlement', 'S')),
        exporter_external_assets_S=book.amount(claim('external_issuer', 'zone_exporter', 'settlement', 'S')),
        foreign_U_liability_constant_U=(customer + dealer) * book.prices_U['U_nom'],
        zone_net_foreign_assets_constant_U=snap['net_foreign_assets_U']['Zone'],
        actor_net_worth_constant_U={name: row['net_worth_U'] for name, row in snap['balance_sheets'].items()},
        global_financial_residual_constant_U=snap['global_financial_residual_U'],
        aggregate_net_worth_constant_U=snap['aggregate_net_worth_U'],
        net_worth_minus_real_assets_constant_U=snap['net_worth_minus_real_assets_U'])


def settlement_loop(config, name):
    """Keep exports/imports fixed; vary who retains U after settlement."""
    spec = config['settlement_loops']
    case = spec['cases'][name]
    actors = {'external_issuer': 'Earth', 'earth_customer': 'Earth', 'earth_dealer': 'Earth',
              'earth_supplier': 'Earth', 'zone_issuer': 'Zone',
              'zone_bank': 'Zone', 'zone_exporter': 'Zone'}
    initial = spec['customer_opening_S']
    fixed_stock = spec['fixed_U_stock']
    prices = config['valuation_prices_constant_U']
    # The opening 120-U stock already belongs to the dealer. The issuer's
    # synthetic domestic real asset closes its opening balance sheet; it is
    # not foreign collateral or a redemption promise. A matching external
    # loan finances the customer's separate S deposit.
    book = Ledger(actors, config['valuation_prices_constant_U'], {
        claim('external_issuer', 'earth_customer', 'settlement', 'S'): initial,
        claim('earth_customer', 'external_issuer', 'opening_loan', 'S'): initial,
        claim('zone_issuer', 'zone_bank', 'base_reserve', 'U_nom'): fixed_stock,
        claim('zone_bank', 'earth_dealer', 'deposit', 'U_nom'): fixed_stock},
        {'zone_issuer': fixed_stock * prices['U_nom']})
    price_U = book.prices_U['U_nom']
    price_S = book.prices_U['S']
    rate_S_per_U = price_U / price_S
    stages = {'opening': capture(book)}
    purchase = case['customer_purchase_U_nom']
    purchase_S = purchase * rate_S_per_U
    book.post('Customer buys existing U from the FX dealer at the assumed quote', {
        claim('external_issuer', 'earth_customer', 'settlement', 'S'): -purchase_S,
        claim('external_issuer', 'earth_dealer', 'settlement', 'S'): purchase_S,
        claim('zone_bank', 'earth_dealer', 'deposit', 'U_nom'): -purchase,
        claim('zone_bank', 'earth_customer', 'deposit', 'U_nom'): purchase})
    stages['customer_acquired_U'] = capture(book)
    exports = spec['service_exports_U_nom']
    book.post('Customer consumes the exported service and pays the exporter', {
        claim('zone_bank', 'earth_customer', 'deposit', 'U_nom'): -exports,
        claim('zone_bank', 'zone_exporter', 'deposit', 'U_nom'): exports},
        net_worth_changes_U={'earth_customer': -exports * price_U,
                             'zone_exporter': exports * price_U},
        flow_class='current_service_export')
    stages['export_paid'] = capture(book)
    sale = case['exporter_sale_U_nom']
    sale_S = sale * rate_S_per_U
    book.post('Exporter sells U back to the dealer for existing S', {
        claim('zone_bank', 'zone_exporter', 'deposit', 'U_nom'): -sale,
        claim('zone_bank', 'earth_dealer', 'deposit', 'U_nom'): sale,
        claim('external_issuer', 'earth_dealer', 'settlement', 'S'): -sale_S,
        claim('external_issuer', 'zone_exporter', 'settlement', 'S'): sale_S})
    stages['exporter_reconverted'] = capture(book)
    imports = spec['service_imports_S']
    book.post('Exporter consumes an Earth service and pays in S', {
        claim('external_issuer', 'zone_exporter', 'settlement', 'S'): -imports,
        claim('external_issuer', 'earth_supplier', 'settlement', 'S'): imports},
        net_worth_changes_U={'zone_exporter': -imports * price_S,
                             'earth_supplier': imports * price_S},
        flow_class='current_service_import')
    stages['imports_paid'] = capture(book)
    final = stages['imports_paid']
    trade = exports * price_U - imports * price_S
    nfa_change = final['zone_net_foreign_assets_constant_U'] - stages['opening']['zone_net_foreign_assets_constant_U']
    if not math.isclose(nfa_change, trade, abs_tol=1e-9):
        raise ValueError('Trade and net external financial position do not reconcile')
    if any(stage['base_U_nom'] != fixed_stock or stage['bank_deposits_U_nom'] != fixed_stock
           for stage in stages.values()):
        raise ValueError('FX settlement changed the fixed U stock')
    return dict(
        stages=stages, final_claims=book.snapshot()['claims'], journal=book.journal,
        exports_constant_U=exports * price_U, imports_constant_U=imports * price_S,
        net_service_receipts_constant_U=trade,
        change_in_zone_net_foreign_assets_constant_U=nfa_change,
        interpretation=f'Each case starts with the same {fixed_stock:g}-U base and deposit stock held initially by an Earth FX dealer. Customers and exporters trade this existing inventory at an exogenous quote; no issuer creates or redeems U and no foreign reserve backing is assumed. The opening domestic real asset is a synthetic balance-sheet initialization, not new wealth from trade or a conversion guarantee. Export receipts and imported services have explicit income/expenditure counterparts. The ledger supplies neither a complete production-cost/GVA account nor an endogenous exchange rate.')


def basket_value(constituents, quotes):
    """Price fixed quantities of U-related currencies in external unit S."""
    if (not constituents or set(constituents) != set(quotes) or
            any(not math.isfinite(x) or x < 0 for x in constituents.values()) or
            not any(constituents.values()) or
            any(not math.isfinite(x) or x <= 0 for x in quotes.values())):
        raise ValueError('Invalid constituent quantities or external quotes')
    return math.fsum(constituents[name] * quotes[name] for name in constituents)


def currency_basket(spec):
    quantities = spec['constituent_units_per_U']
    reference = basket_value(quantities, spec['reference_quotes_S_per_constituent'])
    scenario = basket_value(quantities, spec['scenario_quotes_S_per_constituent'])
    return dict(
        constituent_units_per_U=quantities,
        reference_quotes_S_per_constituent=spec['reference_quotes_S_per_constituent'],
        scenario_quotes_S_per_constituent=spec['scenario_quotes_S_per_constituent'],
        reference_value_S_per_U=reference, scenario_value_S_per_U=scenario,
        external_value_ratio=scenario / reference,
        external_value_change_percent=100 * (scenario / reference - 1),
        rule=spec['rule'])


def partial_external_liquidity(spec):
    """Allocate existing funding once, then stress its external liquid assets."""
    initial = spec['initial_investment_S_value']
    initial_share = spec['initial_reserve_allocation_fraction']
    annual_cash = spec['annual_available_cash_S_value']
    retained_share = spec['earned_cash_reserve_allocation_fraction']
    liability = spec['convertible_liability_value_at_stress_quote_S']
    if (not math.isfinite(initial) or initial < 0 or
            not math.isfinite(liability) or liability <= 0 or
            not 0 <= initial_share <= 1 or not 0 <= retained_share <= 1 or
            any(not math.isfinite(x) or x < 0 for x in annual_cash)):
        raise ValueError('Invalid reserve funding inputs')
    initial_reserves = initial * initial_share
    earned_allocations = [cash * retained_share for cash in annual_cash]
    productive_capital = initial - initial_reserves
    other_cash = math.fsum(annual_cash) - math.fsum(earned_allocations)
    reserves = initial_reserves + math.fsum(earned_allocations)
    funding = initial + math.fsum(annual_cash)
    residual = funding - productive_capital - reserves - other_cash
    if not math.isclose(residual, 0, abs_tol=1e-9):
        raise ValueError('Reserve assets duplicate funding uses')
    stresses = []
    for haircut in spec['reserve_asset_haircuts']:
        if not math.isfinite(haircut) or not 0 <= haircut <= 1:
            raise ValueError('Invalid reserve asset haircut')
        assets = reserves * (1 - haircut)
        for withdrawal in spec['withdrawal_fractions']:
            if not math.isfinite(withdrawal) or not 0 <= withdrawal <= 1:
                raise ValueError('Invalid withdrawal fraction')
            request = liability * withdrawal
            payment = min(assets, request)
            stresses.append(dict(
                reserve_asset_haircut=haircut, withdrawal_fraction=withdrawal,
                reserve_asset_value_after_haircut_S=assets,
                reserve_valuation_loss_S=reserves - assets,
                coverage_of_liabilities=assets / liability,
                requested_external_payment_S=request, available_payment_S=payment,
                remaining_reserves_S=assets - payment,
                unmet_external_payment_S=max(0., request - assets),
                payment_fully_covered=assets >= request))
    return dict(
        initial_investment_S_value=initial,
        initial_reserve_allocation_S=initial_reserves,
        annual_available_cash_S_value=annual_cash,
        annual_earned_cash_reserve_allocation_S=earned_allocations,
        total_funding_S_value=funding, productive_capital_S_value=productive_capital,
        accumulated_external_reserves_S=reserves,
        remaining_cash_for_other_uses_S=other_cash, funding_residual_S=residual,
        convertible_liability_value_at_stress_quote_S=liability,
        coverage_of_liabilities=reserves / liability,
        stresses=stresses, rule=spec['rule'])


def transaction_balance(annual_flow, settlement_share, holding_days, year_days):
    values = (annual_flow, settlement_share, holding_days, year_days)
    if (not all(math.isfinite(v) for v in values) or annual_flow < 0 or
            not 0 <= settlement_share <= 1 or holding_days < 0 or year_days <= 0):
        raise ValueError('Invalid transaction-balance scenario')
    return annual_flow * settlement_share * holding_days / year_days


def portfolio_balances(product, spec):
    summaries = {(r['strategy'], r['market'], r['horizon_elapsed_years']): r
                 for r in product['horizons']}
    checkpoints = {(r['strategy'], r['market']): {x['year']: x for x in r['rows']}
                   for r in product['checkpoints']}
    rows = []
    for (strategy, market), annual in checkpoints.items():
        for horizon in spec['horizons']:
            summary = summaries[strategy, market, horizon]
            row = dict(strategy=strategy, market=market, horizon=horizon)
            failure = summary['failure']
            if failure and failure['year'] <= horizon:
                row.update(status='unfunded_path_no_future_trade_estimate', failure=failure)
            else:
                if horizon not in annual:
                    raise ValueError('Missing funded horizon checkpoint')
                flow = annual[horizon]['operation']['gross_USD2041']
                if flow != summary['final_gross_USD2041']:
                    raise ValueError('Checkpoint and horizon receipts disagree')
                row.update(
                    status='conditional_transaction_balance',
                    annual_gross_receipts_source_value=flow,
                    balance_source_value_by_share_and_holding_days=[[
                        transaction_balance(flow, share, days, spec['settlement_year_days'])
                        for days in spec['holding_days']] for share in spec['settlement_shares']])
            rows.append(row)
    return dict(
        source_value_unit='constant USD2041 from the financial source; valuation equivalent only, not a nominal U quantity or a USD/U exchange rate',
        matrix_rows_settlement_share=spec['settlement_shares'],
        matrix_columns_holding_days=spec['holding_days'],
        parameter_combinations_per_funded_row=len(spec['settlement_shares']) * len(spec['holding_days']),
        settlement_year_days=spec['settlement_year_days'],
        rule=spec['rule'], rows=rows)


def results(input_path=HERE / 'trade_transition_inputs.json'):
    config = inputs(input_path)
    product = portfolio_product(config)
    return dict(
        schema='terluna.research.trade-transition/1', evidence=config['evidence'],
        boundary=config['boundary'], units=config['units'],
        producer_sha256=digest(HERE / 'trade_transition.py'),
        input_sha256=digest(input_path), ledger_sha256=digest(HERE / 'ledger.py'),
        portfolio_source={**config['portfolio_source'],
                          'producer_sha256': product['producer_sha256'],
                          'input_sha256': product['input_sha256'],
                          'evidence': product['evidence']},
        reading_rule='U invoicing, actual settlement and sustained balances are distinct. Customer purchases can be offset by recipient reconversion or supplied from dealer inventory without changing the U stock. Claims always have matching liabilities, and foreign holdings add no unpaired wealth. The basket contains U-related currencies, with S only an external measurement unit. A separate partial external-liquidity example funds reserves once from specified investment and earned-cash allocations; those assets cannot also finance other uses. These scenarios add no revenue, return or capital-financing credit to the existing portfolios. No adoption probability, issuer, peg, unlimited convertibility or endogenous exchange-rate response is estimated.',
        valuation_prices_constant_U=config['valuation_prices_constant_U'],
        settlement_loops={name: settlement_loop(config, name) for name in config['settlement_loops']['cases']},
        currency_basket=currency_basket(config['currency_basket']),
        partial_external_liquidity=partial_external_liquidity(config['partial_external_liquidity']),
        portfolio_transaction_balances=portfolio_balances(product, config['transaction_balance']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'research/runs/monetary_system/trade_transition.json')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results(), indent=2, allow_nan=False) + '\n')
    print(args.output)
