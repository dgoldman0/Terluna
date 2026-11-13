"""Settlement, external claims and transaction stocks have distinct identities."""
import copy
import json

import pytest

from research.studies.monetary_system import trade_transition as m


def final(config, name):
    return m.settlement_loop(config, name)['stages']['imports_paid']


def test_same_real_trade_and_paired_claims_in_all_loops():
    config = m.inputs()
    for name in config['settlement_loops']['cases']:
        loop = m.settlement_loop(config, name)
        assert loop['exports_constant_U'] == 100
        assert loop['imports_constant_U'] == 60
        assert loop['net_service_receipts_constant_U'] == 40
        assert loop['change_in_zone_net_foreign_assets_constant_U'] == 40
        for stage in loop['stages'].values():
            assert stage['global_financial_residual_constant_U'] == 0
            assert stage['aggregate_net_worth_constant_U'] == 120
            assert stage['net_worth_minus_real_assets_constant_U'] == 0
            assert stage['base_U_nom'] == stage['bank_deposits_U_nom'] == 120
        assert all(r['max_balance_sheet_residual_U'] == 0 for r in loop['journal'])
        assert loop['stages']['opening']['zone_net_foreign_assets_constant_U'] == -120
        end = loop['stages']['imports_paid']
        assert end['zone_net_foreign_assets_constant_U'] == -80
        assert end['actor_net_worth_constant_U']['zone_exporter'] == 40
        assert end['actor_net_worth_constant_U']['earth_customer'] == -100
        assert end['actor_net_worth_constant_U']['earth_supplier'] == 60


def test_complete_reconversion_restores_dealer_inventory_without_issuing_U():
    end = final(m.inputs(), 'complete_reconversion')
    assert end['base_U_nom'] == end['bank_deposits_U_nom'] == 120
    assert end['customer_U_nom'] == end['exporter_U_nom'] == 0
    assert end['dealer_U_nom'] == 120
    assert end['dealer_external_assets_S'] == 0
    assert end['exporter_external_assets_S'] == 40
    assert end['zone_net_foreign_assets_constant_U'] == -80


def test_held_receipts_and_foreign_buffer_have_matching_liabilities():
    config = m.inputs()
    domestic = final(config, 'exporter_holds_U')
    foreign = final(config, 'customer_and_exporter_hold_U')
    assert domestic['exporter_U_nom'] == foreign['exporter_U_nom'] == 40
    assert domestic['dealer_U_nom'] == 80
    assert foreign['dealer_U_nom'] == 60
    assert foreign['customer_U_nom'] == 20
    assert domestic['foreign_U_liability_constant_U'] == foreign['foreign_U_liability_constant_U'] == 80
    assert domestic['dealer_external_assets_S'] == 40
    assert foreign['dealer_external_assets_S'] == 60
    assert domestic['zone_net_foreign_assets_constant_U'] == foreign['zone_net_foreign_assets_constant_U']


def test_common_valuation_rescaling_changes_no_nominal_contract_or_real_comparison():
    config = m.inputs()
    scaled = copy.deepcopy(config)
    scaled['valuation_prices_constant_U'] = {k: v * 7 for k, v in config['valuation_prices_constant_U'].items()}
    for name in config['settlement_loops']['cases']:
        original, changed = final(config, name), final(scaled, name)
        for key in ('base_U_nom', 'bank_deposits_U_nom', 'customer_U_nom', 'exporter_U_nom', 'dealer_U_nom', 'dealer_external_assets_S', 'exporter_external_assets_S'):
            assert changed[key] == original[key]
        assert changed['zone_net_foreign_assets_constant_U'] == 7 * original['zone_net_foreign_assets_constant_U']
    assert m.transaction_balance(700, .6, 7, 365) == pytest.approx(7 * m.transaction_balance(100, .6, 7, 365))


def test_constituent_basket_external_value_and_no_external_currency_component():
    config = m.inputs()
    basket = m.currency_basket(config['currency_basket'])
    assert set(basket['constituent_units_per_U']) == {'Ua', 'Ub'}
    assert basket['reference_value_S_per_U'] == 1
    assert basket['scenario_value_S_per_U'] == pytest.approx(.92)
    assert basket['external_value_ratio'] == pytest.approx(.92)
    assert basket['external_value_change_percent'] == pytest.approx(-8)
    # Rescaling the external measurement unit preserves the relative change.
    changed = copy.deepcopy(config['currency_basket'])
    for field in ('reference_quotes_S_per_constituent', 'scenario_quotes_S_per_constituent'):
        changed[field] = {k: v * 7 for k, v in changed[field].items()}
    rescaled = m.currency_basket(changed)
    assert rescaled['scenario_value_S_per_U'] == pytest.approx(7 * .92)
    assert rescaled['external_value_ratio'] == pytest.approx(.92)


@pytest.mark.parametrize('quantities,quotes', [({}, {}), ({'Ua': -1}, {'Ua': 1}),
    ({'Ua': 0}, {'Ua': 1}), ({'Ua': 1}, {'Ua': 0}), ({'Ua': 1}, {'Ub': 1})])
def test_invalid_basket_cannot_silently_change_constituents(quantities, quotes):
    with pytest.raises(ValueError, match='Invalid constituent'):
        m.basket_value(quantities, quotes)


def test_partial_external_liquidity_allocates_each_source_once():
    spec = m.inputs()['partial_external_liquidity']
    result = m.partial_external_liquidity(spec)
    assert result['initial_reserve_allocation_S'] == 5
    assert result['annual_earned_cash_reserve_allocation_S'] == [2, 2, 2]
    assert result['accumulated_external_reserves_S'] == 11
    assert result['productive_capital_S_value'] == 95
    assert result['remaining_cash_for_other_uses_S'] == 24
    assert result['total_funding_S_value'] == 130
    assert result['funding_residual_S'] == 0
    assert result['coverage_of_liabilities'] == .275
    assert (result['productive_capital_S_value'] + result['accumulated_external_reserves_S'] +
            result['remaining_cash_for_other_uses_S']) == result['total_funding_S_value']
    # With no earnings retained, they remain available elsewhere; funding
    # never grows merely because the reserve allocation changes.
    changed = m.partial_external_liquidity(dict(spec, earned_cash_reserve_allocation_fraction=0))
    assert changed['accumulated_external_reserves_S'] == 5
    assert changed['remaining_cash_for_other_uses_S'] == 30
    assert changed['total_funding_S_value'] == result['total_funding_S_value']


def test_partial_reserves_expose_withdrawal_and_asset_value_shortfalls():
    result = m.partial_external_liquidity(m.inputs()['partial_external_liquidity'])
    rows = {(x['reserve_asset_haircut'], x['withdrawal_fraction']): x for x in result['stresses']}
    assert rows[0, .25]['payment_fully_covered']
    assert rows[0, .25]['remaining_reserves_S'] == 1
    assert rows[0, .4]['unmet_external_payment_S'] == 5
    assert not rows[.1, .25]['payment_fully_covered']
    assert rows[.1, .25]['unmet_external_payment_S'] == pytest.approx(.1)
    for row in rows.values():
        assert row['available_payment_S'] + row['remaining_reserves_S'] == pytest.approx(row['reserve_asset_value_after_haircut_S'])
        assert row['available_payment_S'] + row['unmet_external_payment_S'] == pytest.approx(row['requested_external_payment_S'])
        assert row['reserve_asset_value_after_haircut_S'] + row['reserve_valuation_loss_S'] == pytest.approx(11)


def test_transaction_stock_uses_holding_time_and_settlement_share():
    assert m.transaction_balance(365, .6, 7, 365) == pytest.approx(4.2)
    assert m.transaction_balance(365, 0, 30, 365) == 0
    assert m.transaction_balance(365, .85, 0, 365) == 0
    values = [m.transaction_balance(365, .6, days, 365) for days in (1, 7, 30)]
    assert values[0] < values[1] < values[2]
    assert values[2] < 365  # A stock at this assumed turnover, never another annual receipt.


def test_pinned_source_hash_rejects_substitution(tmp_path):
    config = copy.deepcopy(m.inputs())
    config['portfolio_source']['path'] = 'changed.json'
    (tmp_path / 'changed.json').write_text('{}')
    with pytest.raises(ValueError, match='Unreviewed portfolio'):
        m.portfolio_product(config, tmp_path)


def test_all_portfolios_and_failed_horizons_match_pinned_source():
    config = m.inputs()
    source = m.portfolio_product(config)
    output = m.portfolio_balances(source, config['transaction_balance'])
    assert len(output['rows']) == 9 * 3 * 4
    assert output['parameter_combinations_per_funded_row'] == 9
    summaries = {(r['strategy'], r['market'], r['horizon_elapsed_years']): r for r in source['horizons']}
    for row in output['rows']:
        summary = summaries[row['strategy'], row['market'], row['horizon']]
        if summary['failure'] and summary['failure']['year'] <= row['horizon']:
            assert row['status'] == 'unfunded_path_no_future_trade_estimate'
            assert 'balance_source_value_by_share_and_holding_days' not in row
        else:
            assert row['annual_gross_receipts_source_value'] == summary['final_gross_USD2041']
            for share, amounts in zip(output['matrix_rows_settlement_share'], row['balance_source_value_by_share_and_holding_days']):
                for days, amount in zip(output['matrix_columns_holding_days'], amounts):
                    assert amount / row['annual_gross_receipts_source_value'] == pytest.approx(share * days / 365)


def test_missing_funded_checkpoint_is_an_evidence_error():
    config = m.inputs()
    source = copy.deepcopy(m.portfolio_product(config))
    checkpoint = next(x for x in source['checkpoints'] if x['market'] == 'breakthrough')
    checkpoint['rows'] = [x for x in checkpoint['rows'] if x['year'] != 100]
    with pytest.raises(ValueError, match='Missing funded horizon'):
        m.portfolio_balances(source, config['transaction_balance'])


def test_saved_product_reproduces():
    saved = json.loads((m.HERE / 'results/trade_transition.json').read_text())
    assert m.results() == saved
