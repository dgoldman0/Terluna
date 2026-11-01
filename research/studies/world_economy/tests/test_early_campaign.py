import pytest
from research.studies.world_economy.early_campaign import campaign


def run(**changes):
    x=dict(raw_nameplate_kg_year=100e9,years=25,eligible_orders_kg=15e12,
        buyer_budget_U_year=3.5e12,power_W=.5e12,outbound_capacity_kg_year=70e9,
        seconds_per_year=365.25*86400)
    x.update(changes);return campaign(**x)


def test_ramp_and_mass_conservation():
    result=run()
    assert result['delivered_qualified_mass_kg']==pytest.approx(63e9*23.5)
    assert result['delivered_qualified_mass_kg']+result['remaining_orders_kg']==pytest.approx(15e12)
    for row in result['annual']:
        assert row['qualified_delivery_kg']==pytest.approx(row['raw_processed_kg']*.9)
        assert row['local_feed_kg']==pytest.approx(row['raw_processed_kg']*.9)


def test_advance_changes_funding_not_earned_output_or_terminal_surplus():
    base=run();funded=run(customer_advance_fraction=.2)
    for key in ('delivered_qualified_mass_kg','total_earned_producer_revenue_U','total_gva_U','unlevered_project_npv_U','final_cash_after_all_capital_U'):
        assert funded[key]==base[key]
    assert funded['peak_funding_with_customer_advance_U']<base['peak_funding_without_customer_advance_U']
    assert funded['customer_finance_cash_closure_U']==pytest.approx(0,abs=.01)
    assert funded['terminal_customer_liability_U']==pytest.approx(0,abs=.01)


def test_all_production_caps_bind_when_tight():
    caps=[dict(buyer_budget_U_year=4.9e9),dict(outbound_capacity_kg_year=1e9),
          dict(power_W=1e9*100e6/.9/(365.25*86400))]
    for change in caps:
        result=run(**change)
        assert all(row['qualified_delivery_kg']<=1e9*(1+1e-12) for row in result['annual'])


def test_finite_orders_close_campaign_and_do_not_force_idle_replacement():
    result=run(years=50,eligible_orders_kg=50e9)
    assert result['effective_operating_years']==3
    assert result['delivered_qualified_mass_kg']==50e9
    assert not any(r['actual_replacement_U'] for r in result['annual'])


def test_replacement_and_full_capital_cash_close():
    result=run(years=50)
    assert [r['year'] for r in result['annual'] if r['actual_replacement_U']]==[40]
    expected=-result['initial_capital_U']+sum(r['operating_cash_before_capital_U']-r['actual_replacement_U'] for r in result['annual'])+result['terminal_recovery_U']-result['decommission_U']
    assert result['final_cash_after_all_capital_U']==pytest.approx(expected)


def test_joint_programme_receipts_do_not_create_independent_finance():
    result=run(lunar_cost_factor=.75,producer_price_U_kg=4.275)
    assert result['final_cash_after_all_capital_U']>0
    assert result['illustrative_independent_fraction_of_positive_final_surplus_U']==0


@pytest.mark.parametrize('change',[dict(years=0),dict(accepted_yield=0),dict(customer_advance_fraction=1.1),dict(raw_nameplate_kg_year=float('nan'))])
def test_bad_inputs(change):
    with pytest.raises(ValueError):run(**change)
