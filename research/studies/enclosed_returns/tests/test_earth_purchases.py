import math

import pytest

from research.studies.enclosed_returns import earth_purchases as model


@pytest.fixture(scope='module')
def account():
    return model.product()


def test_saved_product_reproduces_and_reconciles_all_financing_modes(account):
    assert account == model.read(model.HERE / 'results/earth_purchases.json.gz')
    assert account['reconciliation']['financial_checkpoints_compared'] == 486


def test_rejected_commodity_orders_are_not_reported_as_earth_purchases(account):
    for case in account['cases']:
        if case['portfolio'] in ('raw_resource_price', 'commodity_refined'):
            for point in case['checkpoints']:
                assert point['cumulative_mass_kg']['Earth_purchased_goods'] == 0
    commodity = next(c for c in account['cases'] if c['name'] == 'commodity_refined:breakthrough:wider')
    assert commodity['checkpoints'][1]['cumulative_mass_kg']['customer_propellant_delivered_in_space'] > 0


def test_propellant_returns_example_contains_no_earth_cargo(account):
    case = next(c for c in account['cases'] if c['name'] == 'propellant_logistics:breakthrough:wider')
    point = case['checkpoints'][1]
    assert point['cumulative_mass_kg']['Earth_purchased_goods'] == 0
    assert point['cumulative_mass_kg']['customer_propellant_delivered_in_space'] / 1e9 == pytest.approx(706.327812, rel=1e-7)
    assert point['cumulative_total_invoices_USD2041'] / 1e12 == pytest.approx(1297.9785, rel=1e-7)


def test_checkpoint_convention_and_disjoint_mass_categories(account):
    for case in account['cases']:
        for point in case['checkpoints']:
            assert point['first_sales_year'] == 8
            assert point['annual_rows_counted'] == point['horizon'] - 7
            m = point['cumulative_mass_kg']
            assert m['Earth_purchased_goods'] == pytest.approx(sum(m[k] for k in (
                'lunar_standard_products', 'lunar_qualified_premium_products', 'asteroid_refined_products')))
            assert m['gross_export_goods_before_losses'] == pytest.approx(
                m['Earth_purchased_goods'] + m['goods_not_sold_after_quality_return_and_qualification'])
            assert math.isclose(sum(point['cumulative_revenue_by_sector_USD2041'].values()),
                                point['cumulative_total_invoices_USD2041'], rel_tol=1e-12)


def test_invoice_inversion_matches_independent_acceptance_and_demand_limits():
    _, orders, cfg, q, profiles, _ = model.inputs()
    quality = cfg['delivered_quality_yield'] * (1 - cfg['return_loss_probability'])
    for case in orders['cases']:
        strategy = q['strategies'][case['portfolio']]
        if not strategy['premium_share']:
            continue
        market = q['markets'][case['market']]
        for row in case['annual']:
            result = model.annual(row, strategy, cfg, q, profiles['inflation_factor'])
            realization = row['industrial_realization']
            # Reconstruct the forward market limit before industrial-budget dispatch.
            attempts = row['goods'] / realization * quality * strategy['premium_share']
            acceptance = min(1., q['premium']['accepted_fraction'] +
                             (q['synergies']['premium_accepted_fraction_gain'] if row['service_kit_count'] else 0.))
            expected = min(attempts * acceptance if row['premium_qualified'] else 0.,
                           market['premium_demand_kg_year'] * row['demand_growth_factor']) * realization
            assert result['mass_kg']['lunar_qualified_premium_products'] == pytest.approx(expected, rel=1e-10, abs=1e-5)


def test_cumulative_goods_are_not_checkpoint_rate_times_elapsed_years(account):
    case = next(c for c in account['cases'] if c['name'] == 'diversified:breakthrough:wider')
    a, b = case['checkpoints']
    assert a['cumulative_mass_kg']['Earth_purchased_goods'] / 1e9 == pytest.approx(5.5158, abs=0.0001)
    assert b['cumulative_mass_kg']['Earth_purchased_goods'] / 1e9 == pytest.approx(273.812, abs=0.0001)
    assert b['cumulative_mass_kg']['Earth_purchased_goods'] < 100 * b['annual_at_checkpoint']['mass_kg']['Earth_purchased_goods']


def test_missing_checkpoint_is_an_error():
    with pytest.raises(ValueError, match='Missing purchase checkpoint'):
        model.checkpoint([], 100)
