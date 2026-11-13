import json
import math
from copy import deepcopy

import pytest

from research.studies.complex_returns import consolidated as m


@pytest.fixture
def sim():
    rows = []
    for year, payout, ui, si, ei, principal, up, draw, sp, origin in (
        (3, 120., 30., 12., 8., 40., 25., 60., 10., 20.),
        (4, 200., 18., 6., 6., 30., 18., 0., 8., 0.),
        (5, 400., 12., 3., 4., 20., 10., 0., 6., 0.),
    ):
        rows.append(dict(year=year, cash_distribution_to_owners=payout,
                         U_interest=ui, supplier_interest=si, external_interest=ei,
                         principal=principal, U_principal=up, external_draw=draw,
                         supplier_principal=sp, supplier_credit_origination=origin,
                         cash=80., U_escrow=30., supplier_prepay=20., debt=50.,
                         custody=15., reserves=70.))
    return dict(initial_hardware_and_carry=600., initial_working_cash=90.,
                initial_reserve=30., external_seed=690., annual=rows)


def test_seed_working_cash_and_monetary_reserve_keep_their_dates(sim):
    flows = m.cashflows(sim, 5, 3, .005)
    assert flows['initial'] == [200., 200., 290., 0., 0., 0.]
    assert flows['dividends'] == [0., 0., 0., 120., 200., 400.]
    assert flows['relief'] == pytest.approx([0., 0., 0., 3., 2., 1.5])
    assert flows['equity'] == pytest.approx([-200., -200., -290., 123., 202., 401.5])
    assert flows['sponsor'] == pytest.approx([-230., -200., -290., 123., 202., 401.5])
    assert sum(flows['initial']) == sim['external_seed']
    result = m.financial_metrics(sim, 5, 3, [.0, .1], .005)
    discounted = result['discounted']['0.1']
    expected_funding = 230 + 200 / 1.1 + 290 / 1.1**2
    assert discounted['initial_sponsor_funding_PV'] == pytest.approx(expected_funding)
    assert discounted['sponsor_cash_NPV_with_saving_release'] == pytest.approx(
        -expected_funding + 123 / 1.1**3 + 202 / 1.1**4 + 401.5 / 1.1**5)


def test_local_interest_release_has_equal_bank_cost_and_excludes_other_claims(sim):
    neutral = m.cashflows(sim, 5, 3)
    changed = m.cashflows(sim, 5, 3, .005)
    assert changed['bank_margin'] == pytest.approx([0., 0., 0., 15., 10., 7.5])
    assert changed['external_lender'] == neutral['external_lender'] == [0., 0., 0., -37., 18., 14.]
    assert changed['supplier_claims'] == neutral['supplier_claims'] == [0., 0., 0., 2., 14., 9.]
    for year in range(6):
        assert changed['sponsor'][year] + changed['bank_margin'][year] == pytest.approx(
            neutral['sponsor'][year] + neutral['bank_margin'][year])
    altered = deepcopy(sim)
    for row in altered['annual']:
        row['external_interest'] *= 100
        row['U_interest'] += 1_000
        row['supplier_interest'] += 1_000
    other_claims = m.cashflows(altered, 5, 3, .005)
    assert other_claims['relief'] == pytest.approx(changed['relief'])
    assert other_claims['bank_margin'] == pytest.approx(changed['bank_margin'])
    metric = m.financial_metrics(sim, 5, 3, [.03], .005)['discounted']['0.03']
    assert metric['U_bank_margin_reduction_PV'] == metric['conditional_saving_release_PV']
    assert metric['sponsor_cash_NPV_with_saving_release'] - metric['sponsor_contractual_cash_NPV'] == pytest.approx(
        metric['conditional_saving_release_PV'])


def test_terminal_book_claims_are_separate_from_paid_cash(sim):
    result = m.financial_metrics(sim, 5, 3, [.0, .1])
    terminal = result['terminal_book_claims']
    assert terminal['net_project_liquid_assets'] == 40
    assert terminal['owner_custody'] == 15
    assert terminal['monetary_reserve_assets'] == 70
    assert terminal['supplier_prepayments_excluded'] == 20
    assert result['cumulative_owner_payout'] == 720
    assert result['sponsor_nondiscounted_cash_return'] == 0
    at10 = result['discounted']['0.1']
    assert at10['terminal_book_claim_PV'] == pytest.approx(125 / 1.1**5)
    assert at10['sponsor_NPV_including_terminal_book_claims'] - at10['sponsor_contractual_cash_NPV'] == pytest.approx(125 / 1.1**5)


@pytest.mark.parametrize('flows,expected', [
    ([-100., 110.], .1),
    ([-100., 0., 121.], .1),
    ([-100., 100.], 0.),
    ([-100., 99.], None),
    ([-100., 230., -132.], None),
    ([100., 10.], None),
    ([0., 0.], None),
])
def test_IRR_requires_a_nonnegative_unambiguous_root(flows, expected):
    result = m.nonnegative_irr(flows)
    if expected is None:
        assert result is None
    else:
        assert result == pytest.approx(expected, abs=1e-12)
        assert m.pv(flows, result) == pytest.approx(0., abs=1e-9)


def test_durable_payback_survives_later_deficits_and_discounting():
    assert m.durable_recovery([-100., 120., -50., 40.]) == 3
    assert m.durable_recovery([-100., 120., -50., 20.]) is None
    assert m.durable_recovery([-100., 110.]) == 1
    assert m.durable_recovery([-100., 110.], .2) is None
    assert m.durable_recovery([-100., 60., 60.], .1) == 2


@pytest.fixture
def economic_inputs():
    q = dict(strategies={'example': {'premium_share': .5}},
             premium={'staff_per_founding_scale': 4.},
             logistics={'support_staff_per_founding_scale': 2.})
    b = dict(cultural_platform={'staff': 3.}, service_kit={'staff': 5.},
             service_staff_USD2041_year=10., seed_scale=.5,
             bootstrap_carrying_fraction=.25, marginal_capital_multiplier=1.5)
    profile = dict(repricing={'propellant_kg_year': 1.},
                   commissioning={'additional_resource_provision_USD2041': 20.},
                   external={'USD2041': {'resource_provision': 10.}},
                   remaining_balanced_growth={
                       'fraction_of_residents_year': .1,
                       'at_limit': {'external': {'USD2041': {'resource_provision': 14.}}}},
                   domestic={'USD2041': {'resident_consumption': 200.}},
                   engineering={'population': 100.})
    old = dict(resident_work_fraction=.6, alternative_USD2041_worker_year=50.)
    return q, b, profile, old


def test_known_resident_payroll_uses_actual_kits_scale_and_activities(economic_inputs):
    q, b, profile, _ = economic_inputs
    row = dict(scale=.5, service_kits=2)
    # Platform3 + kits10 + premium1 + logistics1, each paid10.
    assert m.embedded_resident_payroll(q, b, profile, 'example', row) == 150
    profile['repricing']['propellant_kg_year'] = 0.
    assert m.embedded_resident_payroll(q, b, profile, 'example', row) == 140
    q['strategies']['example']['premium_share'] = 0.
    assert m.embedded_resident_payroll(q, b, profile, 'example', row) == 130


def test_resource_reserves_change_at_founding_boundary(economic_inputs):
    _, b, profile, _ = economic_inputs
    assert m.capital_reserve(profile, b, .4, .2) == pytest.approx(6.)
    assert m.capital_reserve(profile, b, .9, .2) == pytest.approx(12.)
    assert m.capital_reserve(profile, b, 1.2, .2) == pytest.approx(18.)
    assert m.capital_reserve(profile, b, 1.2, 0.) == 0


def test_economic_costs_remove_identified_reserves_and_payroll_once(sim, economic_inputs):
    q, b, profile, old = economic_inputs
    sim['annual'] = [dict(year=year, scale=.5, service_kits=2, starts=starts,
                          industrial_order=100., service_order=50.,
                          operation=dict(gross_USD2041=900., total_cost_USD2041=500.,
                                         resource_provision_USD2041=5.))
                     for year, starts in ((3, .6), (4, .1))]
    result = m.economic_components(sim, q, b, profile, 'example', 5, 3, old)
    # Initial reserve removal: .5*20 commissioning + .5*10*.25*3 carry.
    assert sum(result['procurement_cost'][:3]) == pytest.approx(586.25)
    assert sum(result['removed_resource_provision'][:3]) == pytest.approx(13.75)
    # The second order is post-founding even while its capacity is uncommissioned.
    assert result['removed_resource_provision'][3:5] == pytest.approx([29., 14.])
    assert result['removed_resident_payroll'][3:5] == [150., 150.]
    assert result['procurement_cost'][3:5] == pytest.approx([471., 486.])
    assert result['exports'][3:5] == [900., 900.]
    assert result['resident_services'][3:5] == [100., 100.]
    assert result['resident_opportunity_cost'][3:5] == [1_500., 1_500.]
    # Liquidity capital is neither procurement nor resident opportunity cost.
    sim['initial_working_cash'] += 1_000_000
    sim['initial_reserve'] += 1_000_000
    assert m.economic_components(sim, q, b, profile, 'example', 5, 3, old) == result


def test_economic_valuation_and_switching_values_use_net_new_benefits():
    components = dict(exports=[0., 100., 100.], resident_services=[0., 20., 20.],
                      procurement_cost=[100., 10., 10.], resident_opportunity_cost=[0., 5., 5.],
                      removed_resource_provision=[1e9, 1e9, 1e9],
                      removed_resident_payroll=[1e9, 1e9, 1e9])
    cfg = dict(economic_export_value_fractions=[0., .5, 1.],
               resident_service_value_fractions=[0., .5, 1.],
               resident_opportunity_cost_multipliers=[.5, 1., 2.])
    result = m.economic_metrics(components, 0., [1, 2], cfg)
    grid = {(v['export_value_fraction'], v['resident_value_fraction']): v
            for v in result['valuation_grid']}
    assert grid[1., 1.]['net_value'] == 110
    assert grid[1., 1.]['benefit_cost_ratio'] == pytest.approx(240 / 130)
    assert grid[0., 0.]['additional_constant_net_benefit_to_close_gap'] == 65
    assert grid[.5, 0.]['additional_constant_net_benefit_to_close_gap'] == 15
    assert result['export_value_break_even_with_full_resident_value'] == pytest.approx(.45)
    assert [v['full_price_proxy_net_value'] for v in result['labour_sensitivity']] == [115., 110., 100.]
    discounted = m.economic_metrics(components, .1, [1, 2], cfg)
    assert discounted['constant_annual_benefit_factor'] == pytest.approx(1 / 1.1 + 1 / 1.1**2)
    for value in discounted['valuation_grid']:
        gap = value['additional_constant_net_benefit_to_close_gap']
        if value['net_value'] < 0:
            assert value['net_value'] + gap * (1 / 1.1 + 1 / 1.1**2) == pytest.approx(0., abs=1e-10)
        else:
            assert gap == 0


def test_failed_paths_get_no_future_cash_or_hypothetical_terminal_claim(sim):
    sim['annual'] = [sim['annual'][0], dict(year=4, failure={'type': 'payment_shortfall'})]
    result = m.financial_metrics(sim, 8, 3, [.0, .03], .005)
    flows = m.cashflows(sim, 8, 3, .005)
    assert result['status'] == 'failed_path_cash_to_failure'
    assert result['terminal_book_claims'] is None
    assert result['sponsor_nonnegative_IRR'] is None
    assert result['cumulative_owner_payout'] == 120
    assert result['cumulative_saving_release'] == pytest.approx(3.)
    assert flows['dividends'][4:] == [0.] * 5
    assert flows['relief'][4:] == [0.] * 5
    for row in result['discounted'].values():
        assert row['terminal_book_claim_PV'] is None
        assert row['sponsor_NPV_including_terminal_book_claims'] is None
        assert row['additional_annual_sponsor_cash_to_close_gap'] is None


def test_saved_consolidated_product_reproduces_and_reconciles_sources():
    saved = json.loads((m.HERE / 'results/consolidated.json').read_text())
    assert saved['producer_sha256'] == m.c.s.run.digest(m.HERE / 'consolidated.py')
    for path, digest in saved['dependencies'].items():
        assert m.c.s.run.digest(m.HERE / path) == digest, path
    assert m.product() == saved
    assert len(saved['financial']) == 9 * 3 * 5 * 4
    assert len(saved['economic']) == len(saved['reconciliation']) == 9 * 3 * 4
    for row in saved['reconciliation']:
        steps = ('original_distribution_NPV', 'initial_monetary_reserve',
                 'matched_controller_change', 'commercial_U25_change',
                 'wider_supplier_change', 'conditional_saving_release_change')
        assert sum(row[k] for k in steps) == pytest.approx(row['latest_sponsor_cash_NPV'], abs=.01)
    financial = {(r['strategy'], r['market'], r['horizon'], r['mode']): r for r in saved['financial']}
    for row in saved['financial']:
        if row['mode'] != 'latest_U_package':
            continue
        baseline = financial[row['strategy'], row['market'], row['horizon'], 'wider_U']
        assert row['final_scale'] == baseline['final_scale']
        assert row['cumulative_capital_orders'] == baseline['cumulative_capital_orders']
        assert row['incremental_operating_less_capital_NPV'] == baseline['incremental_operating_less_capital_NPV']
        assert row['price_breaches'] == baseline['price_breaches'] == []
        assert row['max_price_premium_percent'] <= 5 + 1e-7
        for rate, account in row['discounted'].items():
            other = baseline['discounted'][rate]
            # Undiscounted800-year totals reach hundreds of trillions; one
            # representable floating-point step can exceed a cent there.
            expected = other['sponsor_contractual_cash_NPV']
            assert account['sponsor_contractual_cash_NPV'] == pytest.approx(
                expected, rel=0., abs=max(.01, 4*math.ulp(expected)))
            expected_transfer = other['sponsor_cash_NPV_with_saving_release'] + other['U_bank_gross_margin_PV']
            assert account['sponsor_cash_NPV_with_saving_release'] + account['U_bank_gross_margin_PV'] == pytest.approx(
                expected_transfer, rel=0., abs=max(.01, 4*math.ulp(expected_transfer)))
    for row in saved['economic']:
        finance = financial[row['strategy'], row['market'], row['horizon'], 'latest_U_package']
        if finance['status'] == 'failed_path_cash_to_failure':
            assert row['status'] == 'unavailable_after_baseline_failure'
            assert 'rates' not in row


def test_latest_physical_path_satisfies_the_actual_annual_price_gate():
    _, _, _, _, base = m.inputs()
    q, d, b, e, t = m.c.s.inputs()
    control = m.c.simulate(q, d, b, e, t, 'diversified', 'breakthrough', 100,
                           base, 'retained_earnings')
    wider = deepcopy(base)
    wider['financing']['international_U_mix']['outside_supplier_U_fraction'] = .5
    sim = m.c.simulate(q, d, b, e, t, 'diversified', 'breakthrough', 100,
                      wider, 'international_U_mix', control)
    previous = 1.
    assert sim['failure'] is None
    assert sum(r['U_draw'] for r in sim['annual']) > 0
    for row in sim['annual']:
        assert row['pressure_comparison_available']
        assert row['relative_price'] <= 1 + base['monetary']['cumulative_relative_price_limit'] + 1e-8
        assert row['relative_price'] / max(1., previous) - 1 <= base['monetary']['annual_relative_price_limit'] + 1e-8
        previous = row['relative_price']
