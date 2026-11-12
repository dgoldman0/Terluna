"""Accounting identities, economic boundaries and audience break-even controls."""
from copy import deepcopy
import json

import pytest

from humanities.economics import accounts, sports


@pytest.fixture
def spec():
    return accounts.read_bound_inputs(sports.INPUT, 'terluna.humanities.sports-economics-inputs/1')


def test_sports_product_reproduces_and_checks_all_bound_inputs():
    saved = json.loads((accounts.HERE/'results/sports.json').read_text())
    assert sports.results() == saved
    for path, expected in saved['producer']['files'].items():
        assert accounts.digest(accounts.ROOT/path) == expected
    assert saved['legacy_reference'] == dict(revenue_U_year=39000000, cash_before_debt_U_year=7120000)


def test_unreviewed_product_bytes_are_rejected(spec, tmp_path):
    next(iter(spec['products'].values()))['sha256'] = '0'*64
    p = tmp_path/'inputs.json'
    p.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match='Input changed'):
        accounts.read_bound_inputs(p, spec['schema'])


def test_schedule_counts_games_and_reserves_real_arena_hours(spec):
    out = sports.schedule(spec['season'])
    assert out['matches_year'] == 51
    assert out['athletes'] == 32 and out['workers'] == 56
    assert out['occupied_arena_hours'] == 1924
    bad = dict(spec['season'], arena_days=200)
    with pytest.raises(ValueError, match='Arena use'):
        sports.schedule(bad)
    with pytest.raises(ValueError, match='population'):
        sports.schedule(dict(spec['season'], host_residents=40))


def test_wages_remain_in_value_added_and_club_distributions_cancel(spec):
    case = spec['cases']['reference']
    a = sports.evaluate(spec, case)
    assert a['revenue_U_year'] == pytest.approx(25100000)
    assert a['wages_U_year'] == pytest.approx(7840000)
    assert a['intermediate_U_year'] == pytest.approx(13640000)
    assert a['gross_value_added_U_year'] == pytest.approx(11460000)
    assert a['cash_before_debt_U_year'] == pytest.approx(496000)
    clubs = a['clubs']
    assert clubs['aggregate_surplus_before_tax_U_year'] + clubs['central_surplus_before_tax_U_year'] == pytest.approx(a['net_operating_surplus_U_year'])
    spec['club_revenue_fraction'] = 0.8
    b = sports.evaluate(spec, case)
    for field in ('revenue_U_year','gross_value_added_U_year','cash_before_debt_U_year'):
        assert a[field] == b[field]
    assert b['clubs']['central_surplus_before_tax_U_year'] < 0


@pytest.mark.parametrize('coverage', [0, 1, 1.5])
def test_analytic_audience_threshold_against_independent_bisection(spec, coverage):
    ref = spec['cases']['reference']
    lo, hi = 0, 1e7
    for _ in range(60):
        mid = (lo+hi)/2
        row = sports.evaluate(spec, dict(ref, paying_audience_year=mid))
        if row['cash_before_debt_U_year'] < coverage*row['debt_service_U_year']:
            lo = mid
        else:
            hi = mid
    assert sports.audience_threshold(spec, ref, coverage) == pytest.approx(hi, rel=1e-10)


def test_acquisition_can_remove_positive_marginal_returns(spec):
    ref = spec['cases']['reference']
    assert sports.audience_threshold(spec, dict(ref, retention=0, acquisition_U_per_new_customer=30), 1) is None
    low = sports.audience_threshold(spec, dict(ref, retention=0.4), 0)
    assert low > sports.audience_threshold(spec, ref, 0)


def test_debt_and_ownership_change_income_and_cash_without_changing_production(spec):
    case = spec['cases']['global_niche']
    a = sports.evaluate(spec, case)
    spec['finance']['foreign_equity_fraction'] = 1
    spec['finance']['foreign_debt_fraction'] = 0
    spec['finance']['debt_U'] = 0
    b = sports.evaluate(spec, case)
    assert a['gross_value_added_U_year'] == b['gross_value_added_U_year']
    assert a['net_value_added_U_year'] == b['net_value_added_U_year']
    assert a['direct_net_resident_income_U_year'] != b['direct_net_resident_income_U_year']
    assert a['cash_after_debt_U_year'] < b['cash_after_debt_U_year']


def test_import_content_changes_trade_but_not_direct_producer_value_added(spec):
    ref = spec['cases']['reference']
    a = sports.evaluate(spec, ref)
    spec['supplier_import_fraction'] = 0
    b = sports.evaluate(spec, ref)
    assert a['gross_value_added_U_year'] == b['gross_value_added_U_year']
    assert b['external_trade_after_replacement_U_year']-a['external_trade_after_replacement_U_year'] == pytest.approx(a['upstream_imports_U_year'])


def test_retained_foreign_earnings_count_as_income_even_without_dividends(spec):
    case = spec['cases']['global_niche']
    a = sports.evaluate(spec, case)
    spec['finance']['distribution_fraction'] = 0
    b = sports.evaluate(spec, case)
    assert b['foreign_dividends_U_year'] == 0
    assert b['foreign_equity_earnings_U_year'] > 0
    assert a['direct_net_resident_income_U_year'] == b['direct_net_resident_income_U_year']
    assert b['external_cash_after_financing_U_year'] > a['external_cash_after_financing_U_year']


def test_annuity_amortizes_principal_including_zero_rate():
    for rate in (0, 0.05, 0.1):
        balance = 20000000
        pay = accounts.payment(balance, rate, 10)
        for _ in range(10):
            balance = balance*(1+rate)-pay
        assert abs(balance) < 1e-6


@pytest.mark.parametrize('bad', [-1, float('nan'), float('inf'), True])
def test_invalid_audience_rejected(spec, bad):
    with pytest.raises(ValueError):
        sports.evaluate(spec, dict(spec['cases']['reference'], paying_audience_year=bad))


def test_failed_remote_business_retains_payroll_and_host_costs(spec):
    a = sports.evaluate(spec, spec['cases']['remote_receipts_lost'])
    assert a['exports_U_year'] == 0
    assert a['wages_U_year'] > 0 and a['funding_gap_U_year'] > 0
    assert a['external_trade_after_replacement_U_year'] < 0


def test_large_fanbase_scales_without_population_ceiling_or_duplicate_members(spec):
    row = spec['fan_support_cases']['hundred_million_fans']
    a = sports.fan_support(spec, row)
    assert a['exports_U_year'] == pytest.approx(750000000)
    assert a['wages_U_year'] == pytest.approx(150000000)
    assert a['schedule']['workers'] == 56
    assert a['support']['membership_buyers'] < a['support']['paid_viewers']
    assert a['support']['earned_export_U_per_fan'] == 7.5


def test_fan_gifts_and_investment_change_funding_without_becoming_revenue(spec):
    row = spec['fan_support_cases']['hundred_million_fans']
    a = sports.fan_support(spec, row)
    b = sports.fan_support(spec, dict(row, patron_share=0, founding_investor_share=0))
    assert a['revenue_U_year'] == b['revenue_U_year']
    assert a['gross_value_added_U_year'] == b['gross_value_added_U_year']
    assert a['support']['patron_transfers_U_year'] == 20000000
    assert a['support']['one_time_fan_investment_U'] == 50000000
    assert a['support']['cash_after_debt_and_patron_transfers_U_year'] > b['support']['cash_after_debt_and_patron_transfers_U_year']


def test_many_clubs_can_visit_without_being_permanent_residents(spec):
    out = sports.touring_circuit(spec)
    c = out['capacity']
    assert (c['competitions'],c['clubs'],c['registered_athletes']) == (4,32,256)
    assert (c['resident_workers'],c['visiting_athletes']) == (72,224)
    assert c['visitor_arrivals_year'] == 512
    assert c['guest_nights_year'] == 28672
    assert c['peak_guests'] == 128
    assert c['occupied_arena_hours'] == 4264
    assert out['nonresident_wages_U_year'] == 80000000
    assert out['audience']['unique_fans'] == 80000000
    assert out['audience']['sum_of_sport_followings'] == 110000000
    assert out['exports_U_year'] == 600000000


@pytest.mark.parametrize('change', [dict(guest_beds=100),dict(arrival_capacity_year=500),dict(arena_hours_per_day=6)])
def test_touring_capacity_failure_is_not_silently_clipped(spec, change):
    with pytest.raises(ValueError, match='capacity'):
        sports.touring_circuit(spec,dict(spec['touring_circuit'],**change))


def test_visiting_employee_residence_changes_income_without_relocating_output(spec):
    a = sports.touring_circuit(spec)
    b = accounts.account(a['receipts'],a['wages_U_year'],a['purchases'],
        dict(spec['capital'],installed_U=a['installed_capital_U']),spec['finance'],
        spec['tax_fraction'],spec['supplier_import_fraction'],nonresident_wages=0)
    assert b['gross_value_added_U_year'] == a['gross_value_added_U_year']
    assert b['direct_net_resident_income_U_year']-a['direct_net_resident_income_U_year'] == pytest.approx(80000000)
    assert b['external_cash_after_financing_U_year']-a['external_cash_after_financing_U_year'] == pytest.approx(80000000)
