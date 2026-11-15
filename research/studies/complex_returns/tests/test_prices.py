import json
import math

import pytest

from research.studies.complex_returns import prices as m


def test_market_clearing_and_fixed_supply_limit():
    # A 21% larger budget supports price and quantity increases of 10%
    # when supply is proportional to price; with fixed supply all is price.
    response = m.price_response(.21, 1)
    assert response['price_index'] == pytest.approx(110)
    assert response['quantity_index'] == pytest.approx(110)
    assert m.price_response(.21, 0)['price_index'] == pytest.approx(121)
    for demand, elasticity, extra in [(.1, 3, .05), (-.2, 1, 0), (.3, 0, .1)]:
        row = m.price_response(demand, elasticity, extra)
        price, quantity = row['price_index'] / 100, row['quantity_index'] / 100
        assert price * quantity == pytest.approx(1 + demand)
        assert quantity == pytest.approx((1 + extra) * price ** elasticity)


def test_no_spending_recycling_and_supply_driven_deflation():
    assert m.price_response(0, 1)['price_index'] == 100
    assert m.price_response(0, 1, .05)['price_index'] == pytest.approx(100 / math.sqrt(1.05))
    assert m.price_response(.05, 1, .05)['price_index'] == 100


def test_unit_rescaling_and_basket_valuation_are_different():
    case = dict(distribution_share=.01, investment_local_share=.1,
                supply_elasticity=1, extra_supply=0, basket_multiplier=1)
    baseline = m.experiment(1000, 10, 100, case)
    redenominated = m.experiment(7000, 70, 700, case)
    assert baseline['price_index'] == redenominated['price_index']
    larger = m.experiment(1000, 10, 100, dict(case, basket_multiplier=2))
    assert larger['demand_shift'] == baseline['demand_shift'] / 2
    assert larger['price_index'] < baseline['price_index']


def test_only_allocated_distributions_and_expansion_enter_demand():
    case = dict(distribution_share=.01, investment_local_share=.1,
                supply_elasticity=0, extra_supply=0, basket_multiplier=1)
    row = m.experiment(1000, 100, 100, case)
    assert row['distribution_demand_USD2041_equivalent'] == 10
    assert row['investment_demand_USD2041_equivalent'] == 10
    assert row['price_index'] == pytest.approx(120)
    assert m.experiment(1000, 100, 100, dict(case, distribution_share=0,
                        investment_local_share=0))['price_index'] == 100


def test_incomplete_or_failed_window_is_not_a_continuing_economy():
    rows = [dict(year=y) for y in range(8, 31)]
    assert [r['year'] for r in m.window(rows, 30, 10)] == list(range(21, 31))
    assert m.window(rows[:-1], 30, 10) is None
    assert m.window(rows, 100, 10) is None
    assert m.window(rows, 30, 10, failure_year=30) is None
    assert m.window(rows, 30, 10, failure_year=31) is not None


@pytest.mark.parametrize('args', [(-1, 1, 0), (0, -1, 0), (0, 1, -1), (math.nan, 1, 0)])
def test_invalid_price_experiments(args):
    with pytest.raises(ValueError):
        m.price_response(*args)


def test_saved_product_and_no_double_counted_final_uses():
    saved = json.loads((m.HERE / 'results/prices.json').read_text())
    assert m.results() == saved
    assert len(saved['rows']) == 9 * 3 * 4
    for row in saved['rows']:
        if row['status'] != 'conditional_price_level_experiment':
            assert row['failure'] is not None
            assert 'cases' not in row
            continue
        cases = row['cases']
        assert cases['no_recycling']['price_index'] == 100
        assert cases['extra_supply_only']['price_index'] < 100
        assert cases['one_percent_fixed_supply']['price_index'] >= cases['one_percent_responsive_supply']['price_index'] >= cases['one_percent_elastic_supply']['price_index']
        assert cases['one_percent_plus_construction']['price_index'] >= cases['one_percent_responsive_supply']['price_index']
        assert row['cash_after_import_price_shock_USD2041_year'] <= row['mean_available_cash_USD2041_year']
        for epsilon, share in row['distribution_share_at_index_limit'].items():
            if share is not None:
                demand = share * row['distributions_over_basket']
                assert m.price_response(demand, float(epsilon))['price_index'] == pytest.approx(102)
    # All economic values must be included once: domestic final uses exclude
    # the intermediates, gross output, exports and imported capital.
    q, d, *_ = m.strategies.inputs()
    for row in saved['rows']:
        if 'mean_scale' not in row:
            continue
        p = m.strategies.profile(q, d, row['strategy'])['domestic']['USD2041']
        expected = row['mean_scale'] * (p['resident_consumption'] + p['domestic_capital_formation'])
        assert row['domestic_final_use_basket_USD2041_year'] == pytest.approx(expected)
