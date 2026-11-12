"""Annual cultural-export accounts; real U, conditional prices and costs.

The boundary consolidates participating lunar producers. Supplier purchases
are intermediate use; wages are value added. Debt principal is financing.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def nonnegative(value):
    if isinstance(value, bool) or not math.isfinite(value) or value < 0:
        raise ValueError("Expected a finite nonnegative quantity")
    return value


def fraction(value):
    if nonnegative(value) > 1:
        raise ValueError("Expected a fraction in [0, 1]")
    return value


def positive(value):
    if nonnegative(value) == 0:
        raise ValueError("Expected a positive quantity")
    return value


def count(value):
    if positive(value) != int(value):
        raise ValueError("Expected a positive integer count")
    return int(value)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_bound_inputs(path, schema):
    data = json.loads(Path(path).read_text())
    if data['schema'] != schema:
        raise ValueError('Unexpected input schema')
    for name, identity in data.get('products', {}).items():
        product = ROOT / name
        if digest(product) != identity['sha256']:
            raise ValueError(f'Input changed; review before rebinding: {name}')
        if json.loads(product.read_text())['schema'] != identity['schema']:
            raise ValueError(f'Unexpected product schema: {name}')
    return data


def producer(*paths):
    return {'domain': 'humanities', 'files': {p: digest(ROOT / p) for p in paths},
            'constants': constants_used(paths)}


def payment(principal, rate, years):
    nonnegative(principal)
    nonnegative(rate)
    count(years)
    return principal / years if rate == 0 else principal * rate / -math.expm1(-years * math.log1p(rate))


def account(receipts, wages, purchases, capital, finance, tax_fraction,
            supplier_import_fraction, nonresident_wages=0):
    """No inventory change, product taxes or capitalized own-account work.

    Receipts are earned revenue net of offworld distributors' retained fees.
    Each purchase has a direct import fraction; the remaining domestic supply
    has a separately assumed total upstream import content. Capital renewal
    equals depreciation. Employees work on the Moon; nonresident compensation
    is separately identified as primary income paid outside the territory.
    """
    for row in receipts.values():
        nonnegative(row['U_year'])
        fraction(row['export_fraction'])
    for row in purchases.values():
        nonnegative(row['U_year'])
        fraction(row['import_fraction'])
    nonnegative(wages)
    if nonnegative(nonresident_wages) > wages:
        raise ValueError('Nonresident compensation exceeds total compensation')
    fraction(tax_fraction)
    fraction(supplier_import_fraction)
    installed = nonnegative(capital['installed_U'])
    renewal = installed / positive(capital['life_years'])
    renewal_imports = renewal * fraction(capital['import_fraction'])
    revenue = sum(r['U_year'] for r in receipts.values())
    exports = sum(r['U_year'] * r['export_fraction'] for r in receipts.values())
    intermediate = sum(r['U_year'] for r in purchases.values())
    direct_imports = sum(r['U_year'] * r['import_fraction'] for r in purchases.values())
    domestic_inputs = intermediate - direct_imports
    upstream_imports = domestic_inputs * supplier_import_fraction
    gva = revenue - intermediate
    operating_surplus = gva - wages
    net_surplus = operating_surplus - renewal
    tax = max(0, net_surplus) * tax_fraction
    cash = net_surplus - tax
    debt = nonnegative(finance['debt_U'])
    rate = nonnegative(finance['real_rate'])
    service = payment(debt, rate, finance['years'])
    interest = debt * rate
    foreign_debt = fraction(finance['foreign_debt_fraction'])
    dividends = max(0, cash - service) * fraction(finance['distribution_fraction'])
    foreign_dividends = dividends * fraction(finance['foreign_equity_fraction'])
    # Equity income includes retained earnings (and the modeled share of losses).
    # A principal payment changes distributions, not the year's accrued profit.
    foreign_equity_earnings = (net_surplus-tax-interest)*finance['foreign_equity_fraction']
    trade = exports - direct_imports - upstream_imports - renewal_imports
    external_cash = trade - nonresident_wages - service * foreign_debt - foreign_dividends
    return dict(
        receipts=receipts, purchases=purchases, revenue_U_year=revenue,
        exports_U_year=exports, wages_U_year=wages,
        resident_wages_U_year=wages-nonresident_wages,
        nonresident_wages_U_year=nonresident_wages,
        intermediate_U_year=intermediate, direct_imports_U_year=direct_imports,
        domestic_purchases_U_year=domestic_inputs,
        upstream_imports_U_year=upstream_imports,
        gross_value_added_U_year=gva, gross_operating_surplus_U_year=operating_surplus,
        installed_capital_U=installed, replacement_U_year=renewal, replacement_imports_U_year=renewal_imports,
        net_value_added_U_year=gva-renewal, net_operating_surplus_U_year=net_surplus,
        tax_allowance_U_year=tax, cash_before_debt_U_year=cash,
        debt_service_U_year=service, interest_U_year=interest,
        principal_repayment_U_year=service-interest,
        debt_coverage=cash/service if service else None,
        cash_after_debt_U_year=cash-service,
        funding_gap_U_year=max(0, service-cash),
        foreign_dividends_U_year=foreign_dividends,
        foreign_equity_earnings_U_year=foreign_equity_earnings,
        direct_net_resident_income_U_year=gva-renewal-nonresident_wages-interest*foreign_debt-foreign_equity_earnings,
        external_trade_after_replacement_U_year=trade,
        external_cash_after_financing_U_year=external_cash,
        value_added_identity_residual_U=gva-wages-operating_surplus,
    )
