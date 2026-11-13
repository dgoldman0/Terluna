import json
from copy import deepcopy

import pytest

from research.studies.complex_returns import policy as p


@pytest.fixture
def config():
    return json.loads((p.HERE / 'policy_inputs.json').read_text())


@pytest.fixture
def snap():
    # Synthetic annual accounts: a solvent operation with no opening usable cash.
    return dict(revenue=120_000., operating_cost=60_000., core_imports=12_000.,
                cash=0., debt_service=1_200., external_debt_service=600.,
                reserves=100_000., existing_U_loans=10_000.,
                existing_convertible_U=6_000., basket=120_000., price=1.,
                relative_price=1., holding_years=.25, elasticity=1.,
                cumulative_limit=.05, seed=1_000_000.,
                baseline_annual_payout=10_000., baseline_annual_new_orders=20_000.)


def test_rate_relief_excludes_external_and_supplier_interest_and_is_a_transfer():
    sim = {'annual': [
        dict(year=10, U_interest=300., supplier_interest=150., external_interest=900.),
        dict(year=15, failure={'type': 'payment_shortfall'}),
        dict(year=20, U_interest=120., supplier_interest=120., external_interest=900.),
        dict(year=30, U_interest=60., supplier_interest=0., external_interest=900.),
    ]}
    original = deepcopy(sim)
    result = p.rate_relief(sim, 20, .005)
    # Only 150 of the first row's interest belongs to local U lending.
    assert result['realized_borrower_interest_saving'] == pytest.approx(25.)
    assert result['borrower_interest_saving_NPV_3_percent'] == pytest.approx(25 / 1.03**10)
    assert result['bank_gross_income_reduction_NPV_3_percent'] == pytest.approx(
        result['borrower_interest_saving_NPV_3_percent'])
    assert result['consolidated_interest_transfer_NPV'] == 0
    assert result['operating_less_capital_NPV_change'] == 0
    assert result['affected_portion_coverage_capacity_increase'] == pytest.approx(1 / 15)
    assert sim == original
    for row in sim['annual']:
        if 'failure' not in row:
            row['external_interest'] *= 100
            row['U_interest'] += 1_000
            row['supplier_interest'] += 1_000
    assert p.rate_relief(sim, 20, .005) == result
    assert p.rate_relief(original, 30, .005)['realized_borrower_interest_saving'] == pytest.approx(35.)


def test_snapshot_excludes_restricted_and_owner_assets():
    row = dict(gross=200_000., available=60_000., core_imports=20_000.,
               scale=2., cash=8_000., FX_restricted=6_000., U_escrow=90_000.,
               supplier_prepay=10_000., custody=60_000.,
               interest=800., principal=1_000., external_interest=500., U_principal=300.,
               reserves=10_000., U_debt=5_000., outside_U=2_000., price=1.01,
               relative_price=1.002, cash_distribution_to_owners=4_000.,
               industrial_order=3_000., service_order=2_000.)
    profile = {'domestic': {'USD2041': {
        'resident_consumption': 30_000., 'domestic_capital_formation': 10_000.}}}
    mon = dict(domestic_holding_days=90., supply_elasticity=1.,
               cumulative_relative_price_limit=.05)
    before = deepcopy(row)
    result = p.snapshot(row, profile, mon, 1_000_000.)
    assert result['cash'] == 2_000
    assert result['reserves'] == 10_000
    assert result['operating_cost'] == 140_000
    assert result['basket'] == 80_000
    assert result['debt_service'] == 1_800
    assert result['external_debt_service'] == 1_200
    assert result['existing_convertible_U'] == 3_000
    assert result['baseline_annual_new_orders'] == 5_000
    assert row == before
    row['FX_restricted'] = 9_000.
    assert p.snapshot(row, profile, mon, 1_000_000.)['cash'] == 0


@pytest.mark.parametrize('months', [6, 12])
def test_delayed_invoices_conserve_receipts_and_bridge_repayment(config, snap, months):
    cfg = config['liquidity']
    result = p.bridge(snap, cfg, config['policies']['U_term_25'], months, trace=True)
    assert result['status'] == 'paid_and_repaid'
    rows = result['rows']
    assert len(rows) == 2 * months + cfg['collection_tail_months']
    monthly_invoice = snap['revenue'] / 12
    cumulative_collections = 0.
    for row in rows:
        month = row['month']
        expected = 0. if month <= months else (2 * monthly_invoice if month <= 2 * months else monthly_invoice)
        assert row['receipts'] == pytest.approx(expected)
        cumulative_collections += row['receipts']
        assert cumulative_collections + row['receivables'] == pytest.approx(month * monthly_invoice)
        assert row['loan'] == pytest.approx(sum(row[k] for k in ('FX', 'outside_term', 'local', 'local_term')))
        assert row['accessible_FX'] == pytest.approx(snap['reserves'] - row['FX'])
        assert row['loan'] <= (1 - cfg['collateral_haircut']) * row['receivables'] + 1e-7
        assert row['loan'] <= cfg['facility_fraction_of_seed'] * snap['seed'] + 1e-7
        assert cfg['capital_fraction_of_loans'] * (snap['existing_U_loans'] + row['loan']) <= snap['reserves'] + 1e-7
    assert cumulative_collections == pytest.approx(len(rows) * monthly_invoice)
    assert result['terminal_receivables'] == pytest.approx(0., abs=1e-7)
    assert result['terminal_bridge'] == pytest.approx(0., abs=1e-7)
    assert result['gross_draws'] == pytest.approx(result['principal_repaid'])
    assert result['gross_draws'] > 0
    assert result['first_draw_month'] == 1
    assert months < result['final_repayment_month'] <= len(rows)
    assert rows[-1]['accessible_FX'] == pytest.approx(snap['reserves'])
    expected_cash = snap['cash'] + cumulative_collections - len(rows) * (
        snap['operating_cost'] + snap['debt_service']) / 12 - result['borrower_financing_cost']
    assert result['closing_cash'] == pytest.approx(expected_cash)
    assert result['max_accounting_residual'] < 1e-7


def test_term_holding_allows_payments_without_spendable_balance_breach(config, snap):
    cfg = config['liquidity']
    blocked = p.bridge(snap, cfg, config['policies']['U_unabsorbed'], 6)
    held = p.bridge(snap, cfg, config['policies']['U_term_25'], 6, trace=True)
    assert blocked['status'] == 'payment_shortfall'
    assert 'price_balance_allowance' in blocked['failure']['reasons']
    assert blocked['gross_draws'] == 0
    assert held['status'] == 'paid_and_repaid'
    assert held['peak_local_term_claims'] > 0
    assert held['peak_outside_term_claims'] > 0
    assert held['term_holder_interest'] > 0
    assert held['term_principal_redeemed_into_external_savings'] > 0
    assert held['intermediary_margin_before_costs'] == pytest.approx(
        held['borrower_financing_cost'] - held['term_holder_interest'])
    for row in held['rows']:
        implied_price = (snap['price'] ** (1 + snap['elasticity']) +
                         row['local'] / (snap['holding_years'] * snap['basket'])) ** (1 / (1 + snap['elasticity']))
        assert implied_price / snap['price'] <= 1 + cfg['maximum_incremental_price_premium'] + 1e-10


def test_balance_allowance_uses_remaining_cumulative_headroom(config, snap):
    cfg = config['liquidity']
    ordinary = p.balance_allowance(snap, cfg)
    snap.update(relative_price=1.048, price=1.1)
    allowance = p.balance_allowance(snap, cfg)
    implied_price = (snap['price'] ** 2 + allowance / (snap['holding_years'] * snap['basket'])) ** .5
    assert implied_price / snap['price'] * snap['relative_price'] == pytest.approx(1.05)
    assert allowance < ordinary
    snap['relative_price'] = 1.05
    assert p.balance_allowance(snap, cfg) == 0


def test_frozen_FX_access_blocks_external_payments_without_destroying_assets(config, snap):
    normal = config['policies']['U_term_50']
    result = p.bridge(snap, config['liquidity'], dict(normal, FX_access=0.), 6)
    assert result['status'] == 'payment_shortfall'
    assert 'external_liquidity' in result['failure']['reasons']
    assert result['failure']['month'] == 1
    assert result['opening_reserve_assets'] == snap['reserves']
    assert result['minimum_accessible_FX'] == 0
    assert result['gross_draws'] == 0
    assert p.bridge(snap, config['liquidity'], normal, 6)['status'] == 'paid_and_repaid'


@pytest.mark.parametrize('bound', ['collateral', 'facility_limit', 'capital'])
def test_new_credit_cannot_cross_coverage_or_capital_bounds(config, snap, bound):
    cfg = deepcopy(config['liquidity'])
    if bound == 'collateral':
        cfg['collateral_haircut'] = .9
    elif bound == 'facility_limit':
        cfg['facility_fraction_of_seed'] = .001
    else:
        # The old bank meets its capital rule; this advance would breach it.
        snap.update(reserves=10_000., existing_U_loans=65_000.,
                    existing_convertible_U=0., core_imports=0., external_debt_service=0.)
        cfg['unclassified_operating_external_fraction'] = 0.
        assert cfg['capital_fraction_of_loans'] * snap['existing_U_loans'] < snap['reserves']
    result = p.bridge(snap, cfg, config['policies']['U_term_25'], 6)
    assert result['status'] == 'payment_shortfall'
    assert result['failure']['reasons'] == [bound]
    assert result['failure']['month'] == 1
    assert result['gross_draws'] == result['peak_bridge'] == 0


def test_supplier_acceptance_does_not_refinance_external_debt(config, snap):
    # The first month has no cash. Every operating supplier accepts U, but
    # existing external debt service still needs foreign settlement assets.
    spec = dict(config['policies']['U_term_50'], outside_U_share=1.)
    result = p.bridge(snap, config['liquidity'], spec, 6, trace=True)
    assert result['status'] == 'paid_and_repaid'
    first = result['rows'][0]
    assert first['FX'] == pytest.approx(snap['external_debt_service'] / 12)
    expected_suppliers = (snap['core_imports'] + config['liquidity']['unclassified_operating_external_fraction'] *
                          (snap['operating_cost'] - snap['core_imports'])) / 12
    assert first['outside_term'] == pytest.approx(expected_suppliers)
    snap['cash'] = snap['external_debt_service'] / 12
    prepaid_debt = p.bridge(snap, config['liquidity'], spec, 6, trace=True)
    assert prepaid_debt['rows'][0]['FX'] == pytest.approx(0.)
    assert prepaid_debt['rows'][0]['outside_term'] == pytest.approx(expected_suppliers)


def test_permanent_receivable_loss_does_not_cancel_created_claims(config, snap):
    cfg = dict(config['liquidity'], loss_fractions=[0., 1.])
    result = p.bridge(snap, cfg, config['policies']['U_term_25'], 6)
    good, bad = result['collateral_loss_stress']
    assert good['loan_loss'] == 0
    assert bad['loan_loss'] == pytest.approx(result['peak_bridge'])
    assert bad['U_fund_loan_loss'] == bad['loan_loss']
    assert bad['external_lender_loan_loss'] == 0
    assert bad['fund_equity_after_loss'] == pytest.approx(snap['reserves'] - result['peak_bridge'])
    assert bad['issued_U_claims_before_resolution'] == good['issued_U_claims_before_resolution'] > 0


def test_external_lender_loss_does_not_reduce_U_fund_equity(config, snap):
    cfg = dict(config['liquidity'], loss_fractions=[1.])
    result = p.bridge(snap, cfg, config['policies']['external_bridge'], 6)
    assert result['status'] == 'paid_and_repaid'
    loss = result['collateral_loss_stress'][0]
    assert loss['external_lender_loan_loss'] == pytest.approx(result['peak_bridge'])
    assert loss['U_fund_loan_loss'] == 0
    assert loss['fund_equity_after_loss'] == snap['reserves']
    assert loss['issued_U_claims_before_resolution'] == 0
    assert result['minimum_accessible_FX'] == snap['reserves']
    assert result['term_holder_interest'] == 0


def test_liquidity_facility_does_not_roll_structural_operating_losses(config, snap):
    snap['revenue'] = snap['operating_cost'] + snap['debt_service'] - 1_000
    result = p.bridge(snap, config['liquidity'], config['policies']['U_term_25'], 6)
    assert result['status'] == 'payment_shortfall'
    assert 'nonpositive_operating_coverage' in result['failure']['reasons']
    assert result['gross_draws'] == 0


def test_saved_policy_product_reproduces_from_pinned_sources():
    saved = json.loads((p.HERE / 'results/policy.json').read_text())
    assert saved['producer_sha256'] == p.c.s.run.digest(p.HERE / 'policy.py')
    for name, digest in saved['dependencies'].items():
        assert p.c.s.run.digest(p.HERE / name) == digest, name
    assert p.product() == saved
    assert len(saved['funding']) == 9 * 3 * 4
    available = [r for r in saved['liquidity'] if 'policy' in r]
    assert available
    assert any(r['peak_bridge'] > 0 for r in available)
    for row in available:
        assert row['max_accounting_residual'] < 1
        if row['status'] == 'paid_and_repaid':
            assert row['terminal_bridge'] <= 1
            assert abs(row['terminal_receivables']) < 1
            assert row['gross_draws'] == pytest.approx(row['principal_repaid'], abs=1)
        for stress in row['collateral_loss_stress']:
            assert stress['U_fund_loan_loss'] + stress['external_lender_loan_loss'] == pytest.approx(stress['loan_loss'])
    for row in saved['funding']:
        for relief in row['rate_relief'] + row['wider_rate_relief']:
            assert relief['consolidated_interest_transfer_NPV'] == 0
            assert relief['operating_less_capital_NPV_change'] == 0
