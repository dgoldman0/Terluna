"""Capacity, cash conservation, fixed debt and payment-liquidity controls."""
import copy
import json

import pytest

from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import constants_changed
from research.studies.world_economy import sports_finance as m


@pytest.fixture(scope='module')
def product():
    return m.results()


def test_complete_product_and_named_constants_reproduce(product):
    stored=json.loads((m.HERE/'results/sports_finance.json').read_text())
    assert product==stored
    assert not constants_changed(product['producer']['constants'])


def test_changed_upstream_product_requires_rebinding(tmp_path):
    spec,_=m.load_inputs()
    next(iter(spec['products'].values()))['sha256']='0'*64
    path=tmp_path/'changed.json';path.write_text(json.dumps(spec))
    with pytest.raises(ValueError,match='Input changed'):
        m.load_inputs(path)


def test_team_nights_reduce_guest_capacity_inside_occupancy_allowance():
    spec,old=m.load_inputs();venue=spec['venue']
    result=m.visitor_capacity(venue,old['stay_nights'])
    planned=venue['guest_beds']*JULIAN_YEAR_DAYS*venue['occupancy_ceiling']
    cap=result['caps_arrivals_year']['beds']
    assert cap*old['stay_nights']+venue['team_guest_nights']==pytest.approx(planned)
    assert cap==pytest.approx(1161.607142857143)
    unreserved=copy.deepcopy(venue);unreserved['team_guest_nights']=0
    extra=m.visitor_capacity(unreserved,old['stay_nights'])['caps_arrivals_year']['beds']-cap
    assert extra==pytest.approx(venue['team_guest_nights']/old['stay_nights'])
    assert result['limiting_arrivals_year']==cap


@pytest.mark.parametrize('constraint',['training','arena'])
def test_training_or_arena_can_bind_instead_of_beds(constraint):
    spec,old=m.load_inputs();spec=copy.deepcopy(spec)
    if constraint=='training':
        spec['venue']['training_completions_year']=500
    else:
        # 48,000 available participant-hours less 12,000 professional hours.
        spec['venue']['guest_arena_hours']=72
    result=m.operating_case(spec,old,spec['cases']['reference'])
    capacity=result['visitor_capacity']
    assert capacity['caps_arrivals_year'][constraint]==pytest.approx(500)
    assert capacity['caps_arrivals_year']['beds']>500
    assert result['realized_visitors_year']==pytest.approx(500)
    assert result['onsite_receipts_U']==pytest.approx(500*old['stay_nights']*old['destination_U_per_night'])


def test_operating_cash_closes_after_compensation_replacement_and_tax(product):
    for case in product['operating_scenarios'].values():
        row=case['operation']
        assert row['total_receipts_U']==pytest.approx(row['onsite_receipts_U']+row['rights_receipts_U']+row['sponsorship_receipts_U'])
        uses=sum(row[key] for key in ('variable_cash_cost_U','fixed_cash_cost_U','replacement_U','tax_allowance_U'))
        assert row['total_receipts_U']-uses==pytest.approx(row['cash_available_for_debt_service_U'])
        if row['pre_tax_cash_U']<=0:
            assert row['tax_allowance_U']==0
        else:
            assert row['tax_allowance_U']==pytest.approx(.2*row['pre_tax_cash_U'])
    base=product['operating_scenarios']['reference']['operation']
    assert base['onsite_receipts_U']==pytest.approx(14e6)
    assert base['replacement_U']==pytest.approx(4.5e6)
    assert base['cash_available_for_debt_service_U']==pytest.approx(7.12e6)


def test_capital_sensitivity_increases_renewal_without_scaling_sales(product):
    base=product['operating_scenarios']['reference']['operation']
    doubled=product['operating_scenarios']['capital_double']['operation']
    assert doubled['total_receipts_U']==base['total_receipts_U']
    assert doubled['installed_capital_U']==2*base['installed_capital_U']
    assert doubled['replacement_U']==2*base['replacement_U']
    assert doubled['fixed_cash_cost_U']-base['fixed_cash_cost_U']==pytest.approx(1.2e6)
    assert doubled['cash_available_for_debt_service_U']==pytest.approx(2.56e6)
    assert product['operating_scenarios']['capital_double']['additional_upfront_capital_U']==pytest.approx(120e6)
    assert product['operating_scenarios']['reference']['additional_upfront_capital_U']==0


def test_remote_audience_can_grow_while_visitors_remain_capacity_limited(product):
    base=product['operating_scenarios']['reference']['operation']
    growth=product['operating_scenarios']['audience_growth']['operation']
    assert growth['visitor_capacity']==base['visitor_capacity']
    assert growth['realized_visitors_year']==pytest.approx(growth['visitor_capacity']['limiting_arrivals_year'])
    assert growth['realized_visitors_year']<product['assumptions']['cases']['audience_growth']['visitor_demand_year']
    assert growth['rights_receipts_U']>base['rights_receipts_U']
    assert growth['installed_capital_U']==base['installed_capital_U']


def test_issue_caps_and_sources_uses_close(product):
    cases=[product['reference_bond'],*product['financing_sensitivity'].values(),
           *product['investor_demand_sensitivity'].values()]
    uses=product['operating_scenarios']['reference']['operation']['installed_capital_U']+product['assumptions']['working_cash_U']
    # Also exercise a fully funded project whose uses, rather than cash, cap debt.
    all_debt=m.bond_capacity(1e9,1000,.05,10,1.5,.5,.02,1e9)
    cases.append(all_debt)
    for bond in cases:
        principal=bond['gross_principal_U']
        assert principal==pytest.approx(min(bond['principal_limits_U'].values()))
        assert bond['sources_minus_uses_U']==pytest.approx(0,abs=1e-7)
        assert principal==pytest.approx(bond['net_proceeds_for_capital_and_working_cash_U']+bond['debt_service_reserve_U']+bond['issuance_cost_U'])
        required_uses=1000 if bond is all_debt else uses
        assert bond['net_proceeds_for_capital_and_working_cash_U']+bond['other_funding_required_U']==pytest.approx(required_uses)
    assert all_debt['other_funding_required_U']==pytest.approx(0,abs=1e-9)
    demand=product['investor_demand_sensitivity']['10000000.0']
    assert demand['gross_principal_U']==pytest.approx(10e6)
    assert demand['gross_principal_U']<product['reference_bond']['gross_principal_U']


def test_reference_debt_is_covered_and_amortizes_by_independent_cash_schedule(product):
    bond=product['reference_bond'];cash=product['operating_scenarios']['reference']['operation']['cash_available_for_debt_service_U']
    assert cash/bond['annual_service_U']==pytest.approx(bond['coverage_target'])
    balance=bond['gross_principal_U'];rate=bond['rate_real'];payment=bond['annual_service_U']
    for _ in range(bond['years']):
        balance=balance*(1+rate)-payment
    assert balance==pytest.approx(0,abs=1e-6)
    discounted=sum(payment/(1+rate)**year for year in range(1,bond['years']+1))
    assert discounted==pytest.approx(bond['gross_principal_U'])
    first=bond['first_year']
    assert first['interest']+first['principal_repaid']==pytest.approx(payment)
    assert bond['gross_principal_U']-first['principal_repaid']==pytest.approx(first['closing_debt'])


def test_operating_stresses_keep_reference_service_and_reserve_cash_closes(product):
    bond=product['reference_bond']
    for case in product['operating_scenarios'].values():
        cash=case['operation']['cash_available_for_debt_service_U'];stress=case['fixed_reference_bond']
        assert cash-stress['cash_after_debt_before_reserve_U']==pytest.approx(bond['annual_service_U'])
        assert stress['closing_debt_reserve_U']+stress['debt_reserve_draw_U']==pytest.approx(bond['debt_service_reserve_U'])
        supplied=cash+stress['funding_gap_after_debt_reserve_before_working_cash_U']+stress['debt_reserve_draw_U']
        assert supplied-bond['annual_service_U']==pytest.approx(stress['cash_after_service_positive_U'],abs=1e-8)
        assert stress['uncovered_debt_payment_U']+stress['debt_reserve_draw_U']==pytest.approx(stress['debt_gap_before_reserve_U'])


def test_debt_reserve_does_not_cover_a_shutdown_operating_deficit(product):
    case=product['operating_scenarios']['shared_venue_shutdown'];stress=case['fixed_reference_bond']
    assert stress['operating_funding_gap_U']==pytest.approx(11.7e6)
    assert stress['debt_reserve_draw_U']>0
    assert stress['funding_gap_after_debt_reserve_before_working_cash_U']>=stress['operating_funding_gap_U']
    assert stress['funding_gap_after_debt_reserve_before_working_cash_U']==pytest.approx(stress['operating_funding_gap_U']+stress['uncovered_debt_payment_U'])


def test_media_break_even_reaches_operating_and_debt_cash_targets(product):
    spec,old=m.load_inputs();case=copy.deepcopy(spec['cases']['reference'])
    for field,target in [('combined_media_and_sponsor_receipts_before_debt_U',0),
                         ('combined_media_and_sponsor_receipts_to_pay_reference_debt_U',product['reference_bond']['annual_service_U'])]:
        case['rights_receipts_U']=product['break_even'][field]
        case['sponsorship_receipts_U']=0
        row=m.operating_case(spec,old,case)
        assert row['cash_available_for_debt_service_U']==pytest.approx(target,abs=1e-8)


def test_fx_changes_only_local_net_cash_and_can_help_a_local_deficit():
    assert m.fx_cash(20e6,15e6,30e6,20e6,1)==pytest.approx(15e6)
    assert m.fx_cash(20e6,15e6,30e6,20e6,.7)==pytest.approx(12e6)
    assert m.fx_cash(20e6,15e6,20e6,20e6,.7)==pytest.approx(5e6)
    old=m.fx_cash(20e6,15e6,20e6,30e6,1)
    depreciated=m.fx_cash(20e6,15e6,20e6,30e6,.7)
    assert depreciated-old==pytest.approx(3e6)


def test_redemption_preserves_book_equity_while_exposing_illiquidity():
    for diversion in (0,6e6):
        row=m.redeemable_balances(10e6,7e6,diversion)
        assert row['redemption_paid_U']+row['outstanding_claims_U']==pytest.approx(10e6)
        assert row['redemption_paid_U']+row['unpaid_redemption_U']==pytest.approx(7e6)
        assert row['liquid_assets_U']+row['illiquid_asset_at_book_U']==pytest.approx(row['outstanding_claims_U'])
        assert row['net_assets_at_book_U']==pytest.approx(0)
        assert row['outstanding_claims_U']-row['liquid_assets_U']==pytest.approx(diversion)
    backed=m.redeemable_balances(10e6,7e6,0)
    diverted=m.redeemable_balances(10e6,7e6,6e6)
    assert backed['unpaid_redemption_U']==0
    assert diverted['unpaid_redemption_U']==pytest.approx(3e6)
    assert diverted['liquid_assets_U']==0


def test_negative_operating_cash_cannot_support_new_debt():
    bond=m.bond_capacity(-1e6,125e6,.05,10,1.5,.5,.02,60e6)
    assert bond['gross_principal_U']==bond['annual_service_U']==bond['debt_service_reserve_U']==0
    assert bond['other_funding_required_U']==pytest.approx(125e6)
    stress=m.fixed_bond_stress(-1e6,bond)
    assert stress['operating_coverage'] is None
    assert stress['funding_gap_after_debt_reserve_before_working_cash_U']==pytest.approx(1e6)
