import json
import math

import pytest

from research.studies.monetary_system import pressure as m


CONFIG = json.loads((m.HERE/'pressure_inputs.json').read_text())
BASE = CONFIG['cases']['supply_expansion']


def test_issuance_meeting_new_outside_holdings_preserves_local_prices():
    for outside in (0, 10, 10000):
        r = m.money_price(100, 110, outside, BASE, .05)
        assert r['price_level'] == pytest.approx(math.sqrt(1.1/1.05))
        assert r['local_balance_change_value'] == 0
        assert r['issued_for_outside_balances_value'] == outside
        assert r['price_level']*r['quantity_level'] == pytest.approx(1.1)


def test_stock_conservation_and_persistent_shortfall():
    local = 100*90/365
    for nu in (.5,1,1.5):
        case = dict(BASE,issuance_accommodation=nu)
        r = m.money_price(100,100,.2*local,case,0)
        assert r['remaining_local_balances_value']+r['extra_net_outside_balances_value'] == pytest.approx(local+r['issued_for_outside_balances_value'])
        assert r['price_level'] == pytest.approx(math.sqrt(1+(nu-1)*.2))
        # A stock shortfall persists without a further flow shock.
        assert m.money_price(100,100,.2*local,case,0)['price_level'] == r['price_level']


def test_outside_dealer_transfers_do_not_drain_domestic_balances():
    case = dict(BASE,issuance_accommodation=0.,dealer_inventory_absorption=1.)
    r = m.money_price(100,100,1000,case,0)
    assert r['extra_net_outside_balances_value'] == 0
    assert r['outside_dealer_inventory_reduction_value'] == 1000
    assert r['price_level'] == 1


def test_outside_demand_reversal_requires_retirement():
    fixed = m.money_price(100,100,-5,dict(BASE,issuance_accommodation=0),0)
    retired = m.money_price(100,100,-5,BASE,0)
    assert fixed['price_level'] > 1
    assert retired['price_level'] == 1
    assert retired['issued_for_outside_balances_value'] == -5


def test_exhausted_money_stock_is_rejected_as_a_price_result():
    r = m.money_price(100,100,100,dict(BASE,issuance_accommodation=0),0)
    assert r['status'] == 'infeasible_domestic_balance_exhaustion'
    assert 'price_level' not in r


@pytest.mark.parametrize('cash,payout,expected',[(100,50,20),(100,5,5),(-100,5,0),(100,0,0)])
def test_reserves_use_distributions_without_spending_cash_twice(cash,payout,expected):
    reserve, remaining = m.reserve_allocation(cash,payout,.2)
    assert reserve == expected
    assert reserve+remaining == payout


def test_unit_scaling_is_not_a_currency_conversion():
    a = m.money_price(100,110,10,BASE,.05)
    b = m.money_price(3700,4070,370,BASE,.05)
    assert a['price_level'] == pytest.approx(b['price_level'])
    assert a['monetary_balance_gap_fraction'] == pytest.approx(b['monetary_balance_gap_fraction'])


def test_import_quote_direction_and_asymmetry():
    rows = m.import_cost_stresses(CONFIG['import_stress'])
    assert all(r['domestic_basket_change_percent']<0 for r in rows if r['U_external_value_factor']>1)
    assert all(r['domestic_basket_change_percent']>0 for r in rows if r['U_external_value_factor']<1)
    assert max(r['domestic_basket_change_percent'] for r in rows) == pytest.approx(2.5)


def test_saved_product_paths_and_price_changes():
    saved = json.loads((m.HERE/'results/pressure.json').read_text())
    assert m.results() == saved
    assert len(saved['rows']) == 108
    funded = [r for r in saved['rows'] if 'cases' in r]
    assert len(funded) == 73
    for row in saved['rows']:
        if 'cases' not in row:
            assert row['failure']['year'] <= row['horizon']
            continue
        cases = row['cases']
        assert row['additional_initial_reserve_capital_value'] == pytest.approx(.05*row['original_seed_value'])
        for c in cases.values():
            assert c['reduced_distributions_value'] >= 0
            assert c['reduced_distributions_value']+c['reserve_contributions_value'] == pytest.approx(c['original_distributions_value'])
        for key,path in row['paths'].items():
            assert len(path) == 11
            assert path[0]['price_level'] == pytest.approx(1.)
            if cases[key]['status'] != 'conditional_reference_quote':
                continue
            annual_product = math.prod(1+r['annual_price_change_percent']/100 for r in path[1:])
            assert annual_product == pytest.approx(path[-1]['price_level'])
            assert (1+cases[key]['geometric_annual_change_percent']/100)**10 == pytest.approx(annual_product)
        # Extra reserve funding reduces spending; it cannot reduce supply costs.
        assert cases['reserves_already_funded']['closing']['price_level'] >= cases['supply_expansion']['closing']['price_level']
