"""Portfolio overlap, inherited capacity and incremental production controls."""
from copy import deepcopy
import json

import pytest

from humanities.economics import accounts, entertainment


@pytest.fixture
def data():
    spec = accounts.read_bound_inputs(entertainment.INPUT,'terluna.humanities.entertainment-inputs/1')
    sport = json.loads((accounts.ROOT/entertainment.SPORTS_PRODUCT).read_text())
    return spec,sport


def test_product_reproduces_with_pinned_sports_input():
    actual = entertainment.results()
    assert actual == json.loads((accounts.HERE/'results/entertainment.json').read_text())
    assert actual['sports_input_sha256'] == accounts.digest(accounts.ROOT/entertainment.SPORTS_PRODUCT)


def test_removing_new_programmes_exactly_recovers_sports_accounts(data):
    spec,sport = data
    out = entertainment.evaluate(spec,sport,spec['cases']['sports_control'])
    for k in ('revenue_U_year','exports_U_year','gross_value_added_U_year','resident_wages_U_year',
              'nonresident_wages_U_year','intermediate_U_year','cash_before_debt_U_year',
              'external_cash_after_financing_U_year','direct_net_resident_income_U_year'):
        assert out[k] == pytest.approx(sport['touring_circuit'][k]),k


def test_bundle_counts_people_once_and_does_not_add_sports_again(data):
    spec,sport = data
    out = entertainment.evaluate(spec,sport,spec['cases']['expanded'])
    assert out['audience']['unique_fans'] == 140000000
    assert out['audience']['sum_genre_followings'] == 200000000
    assert out['receipts']['portfolio_viewing']['U_year'] == 700000000
    assert out['exports_U_year'] == 1202000000
    assert out['increment_against_sport']['exports_U_year'] == 602000000


def test_extra_programmes_need_extra_demand_or_price_to_improve_cash(data):
    spec,sport = data
    out = entertainment.evaluate(spec,sport,spec['cases']['shared_audience_same_price'])
    assert out['exports_U_year'] > sport['touring_circuit']['exports_U_year']
    assert out['increment_against_sport']['cash_before_debt_U_year'] < 0


def test_added_cast_and_programmes_fit_in_remaining_shared_capacity(data):
    spec,sport = data
    out = entertainment.evaluate(spec,sport,spec['cases']['expanded'])
    cap = out['capacity']
    assert (cap['resident_workers'],cap['arrivals_year'],cap['peak_guests']) == (104,592,168)
    assert cap['guest_nights_year'] == 30912
    assert cap['arena_hours_year'] == 4552
    assert cap['studio_hours_year'] == 900
    assert cap['additional_crew_hours'] == 17568


@pytest.mark.parametrize('field,value', [('visiting_workers',50),('resident_workers',1),('studio_available_hours_year',500)])
def test_capacity_shortfall_is_reported(data,field,value):
    spec,sport = data
    spec['production'][field] = value
    with pytest.raises(ValueError,match='capacity'):
        entertainment.evaluate(spec,sport,spec['cases']['expanded'])


def test_common_revenue_shock_does_not_assume_instant_payroll_reductions():
    out = entertainment.results()['cases']
    assert out['shared_revenue_shock']['wages_U_year'] == out['expanded']['wages_U_year']
    assert out['shared_revenue_shock']['exports_U_year'] == out['expanded']['exports_U_year']/2
    assert out['shared_revenue_shock']['cash_before_debt_U_year'] < out['expanded']['cash_before_debt_U_year']/2


def test_creator_patron_transfers_are_outside_value_added(data):
    spec,sport = data
    a = entertainment.evaluate(spec,sport,spec['cases']['expanded'])
    spec['patron_share'] = 0
    b = entertainment.evaluate(spec,sport,spec['cases']['expanded'])
    assert a['gross_value_added_U_year'] == b['gross_value_added_U_year']
    assert a['patron_transfers_U_year'] == 28000000
    assert b['patron_transfers_U_year'] == 0
