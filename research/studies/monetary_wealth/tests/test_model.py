"""Ownership, monetary closure and immutable-source tests above the ledger kernel.

Outside sectors are the counterparties represented by this model. Their recorded
net worth is not a complete balance sheet of Earth households, firms or banks.
"""
from copy import deepcopy
import json
import math

import pytest
from research.studies.monetary_wealth import model, run


@pytest.fixture(scope='module')
def loaded():
    return run.load()


@pytest.fixture(scope='module')
def reference(loaded):
    cfg,finance,products=loaded
    case=next(c for c in finance['cases'] if c['name']==cfg['reference_case'])
    assets=next(c for c in products['asset_product']['cases'] if c['name']==case['name'])
    return cfg,case,assets


@pytest.fixture(scope='module')
def reference_paths(reference):
    cfg,case,assets=reference
    before=json.dumps((cfg,case,assets),sort_keys=True)
    result={name:model.path(case,assets,cfg,name) for name in cfg['ownership_scenarios']}
    assert json.dumps((cfg,case,assets),sort_keys=True)==before
    return result


def simple_flow(**changes):
    value=dict(cash=20.,reserves=10.,debt=0.,local_U_debt=0.,supplier_U_debt=0.,
               external_debt=0.,owner_custody=0.)
    value.update(changes)
    return value


def simple_asset():
    return dict(net_remaining_service_value_USD2041=100.,construction_WIP_value_USD2041=0.)


def settings(loaded,**changes):
    value=deepcopy(loaded[0]);value.update(changes)
    return value


def balance(flow,asset,cfg,**changes):
    arguments=dict(local_share=0.,household_cash=0.,fund_cash=0.,bank_equity=0.,
                   local_custody=0.,outside_U=0.,settings=cfg)
    arguments.update(changes)
    return model.balance(flow,asset,**arguments)


def transaction_money(result):
    return result['money']['resident_transaction_money']+result['money']['external_transaction_money']


def synthetic_path():
    flows=[];assets=[]
    for y in range(221):
        f=simple_flow()
        f.update(year=y,population=10.,outside_equity_contribution=130. if y==0 else 0.,
                 owner_payout=0.,dividend=0.,gross=0.,resident_payroll=0.,
                 local_U_draw=0.,local_U_principal=0.,supplier_U_draw=0.,
                 supplier_U_principal=0.,local_U_interest=0.)
        flows.append(f);assets.append({'year':y,'totals':simple_asset(),'by_functional_sector':{}})
    case=dict(name='synthetic:breakthrough:wider',portfolio='synthetic',growth='wider',
              seed_cost_scope='Synthetic accounting example',annual=flows)
    return case,{'name':case['name'],'annual':assets}


def test_purchase_and_gift_preserve_payouts_and_household_fund_cash(reference,reference_paths):
    cfg,case,_=reference
    for ownership,path in reference_paths.items():
        old_household=old_fund=0.
        for flow,row in zip(case['annual'],path['annual']):
            f=row['flows'];o=row['ownership_flow']
            assert f['local_payouts']+f['outside_payouts']==pytest.approx(flow['owner_payout'],rel=2e-15,abs=.1)
            assert f['gross_receipts']==flow['gross']
            assert f['outside_contributions']==flow['outside_equity_contribution']
            assert row['programme_cash']==flow['cash'] and row['programme_reserves']==flow['reserves']
            assert row['programme_debt']==flow['debt']
            assert row['household_retained_financial_assets']==pytest.approx(old_household+f['household_saving']-f['community_contributions'],rel=2e-15,abs=.1)
            assert row['community_uninvested_financial_assets']==pytest.approx(old_fund+f['community_contributions']+f['local_payouts']-f['purchase_spending'],rel=2e-15,abs=.1)
            assert row['local_economic_ownership_fraction']+row['outside_economic_ownership_fraction']==pytest.approx(1)
            assert row['local_economic_ownership_fraction']==pytest.approx(o['after_capital_dilution']+o['purchased_fraction']+o['gifted_fraction'])
            assert f['purchase_spending']==f['outside_buyout_receipts']
            assert min(f['purchase_spending'],f['outside_payouts'],f['local_payouts'])>=0
            old_household=row['household_retained_financial_assets'];old_fund=row['community_uninvested_financial_assets']
        if ownership=='external_equity':
            assert all(r['local_economic_ownership_fraction']==0 for r in path['annual'])
        elif ownership=='community_purchase':
            assert path['annual'][200]['local_economic_ownership_fraction']>0
            assert all(r['flows']['gift_book_value']==0 for r in path['annual'])
        else:
            assert path['annual'][100]['local_economic_ownership_fraction']==pytest.approx(cfg['vesting_target_fraction'])
            assert all(r['flows']['purchase_spending']==0 for r in path['annual'])
            for row in path['annual']:
                assert row['flows']['gift_book_value']==pytest.approx(row['ownership_flow']['gifted_fraction']*row['book_balance']['own_funds_equity'])


def test_fair_value_purchase_preserves_wealth_under_consistent_mark(loaded):
    cfg=settings(loaded);flow=simple_flow();asset=simple_asset();quote=200.;purchase=20.
    before=balance(flow,asset,cfg,fund_cash=purchase,equity_mark=quote)
    after=balance(flow,asset,cfg,local_share=purchase/quote,equity_mark=quote)
    assert after['consolidated_off_earth']['real_assets']==before['consolidated_off_earth']['real_assets']==100
    assert after['consolidated_off_earth']['net_worth']==pytest.approx(before['consolidated_off_earth']['net_worth'])
    assert after['sectors']['community_fund']['net_worth']==pytest.approx(before['sectors']['community_fund']['net_worth'])
    # At the different own-funds mark, the purchase premium reduces book wealth.
    old_book=balance(flow,asset,cfg,fund_cash=purchase)
    new_book=balance(flow,asset,cfg,local_share=purchase/quote)
    premium=purchase-(purchase/quote)*model.own_funds(flow,asset)
    assert new_book['consolidated_off_earth']['net_worth']-old_book['consolidated_off_earth']['net_worth']==pytest.approx(-premium)


def test_equity_gift_transfers_claim_value_without_creating_real_assets_or_money(loaded):
    cfg=settings(loaded);flow=simple_flow();asset=simple_asset();quote=200.
    before=balance(flow,asset,cfg,equity_mark=quote)
    after=balance(flow,asset,cfg,local_share=.25,equity_mark=quote)
    assert after['consolidated_off_earth']['net_worth']-before['consolidated_off_earth']['net_worth']==pytest.approx(50)
    assert after['sectors']['outside_owners_customers']['net_worth']-before['sectors']['outside_owners_customers']['net_worth']==pytest.approx(-50)
    assert after['consolidated_off_earth']['real_assets']==before['consolidated_off_earth']['real_assets']
    assert transaction_money(after)==transaction_money(before)


def test_outside_deposit_demand_changes_money_and_backing_without_local_wealth(loaded):
    cfg=settings(loaded);flow=simple_flow();asset=simple_asset()
    before=balance(flow,asset,cfg,outside_U=0)
    after=balance(flow,asset,cfg,outside_U=50)
    assert transaction_money(after)-transaction_money(before)==pytest.approx(50)
    assert after['money']['external_transaction_money']==50
    assert after['money']['resident_transaction_money']==before['money']['resident_transaction_money']
    assert after['monetary_external_assets']-before['monetary_external_assets']==pytest.approx(50)
    assert after['consolidated_off_earth']['external_financial_assets']-before['consolidated_off_earth']['external_financial_assets']==pytest.approx(50)
    assert after['consolidated_off_earth']['external_liabilities']-before['consolidated_off_earth']['external_liabilities']==pytest.approx(50)
    assert after['consolidated_off_earth']['net_worth']==pytest.approx(before['consolidated_off_earth']['net_worth'])
    assert after['sectors']['banks']['net_worth']==pytest.approx(0,abs=1e-12)
    assert after['sectors']['settlement_issuer']['net_worth']==pytest.approx(0,abs=1e-12)


def test_reserve_relocation_keeps_enterprise_claim_and_creates_no_public_gift(loaded):
    cfg=settings(loaded);asset=simple_asset()
    before=balance(simple_flow(cash=30,reserves=0),asset,cfg)
    after=balance(simple_flow(cash=15,reserves=15),asset,cfg)
    assert after['own_funds_equity']==before['own_funds_equity']==130
    assert after['consolidated_off_earth']['net_worth']==pytest.approx(before['consolidated_off_earth']['net_worth'])
    assert after['sectors']['community_fund']['net_worth']==0
    assert after['sectors']['settlement_issuer']['net_worth']==pytest.approx(0,abs=1e-12)
    restricted=[c for c in after['claims'] if c['instrument']=='restricted_fund']
    assert len(restricted)==1 and restricted[0]['holder']=='enterprises' and restricted[0]['value']==15
    assert restricted[0]['money_tags']==['excluded_instrument']
    assert after['money']['resident_transaction_money']<before['money']['resident_transaction_money']


def test_local_loan_draw_and_repayment_change_gross_claims_not_own_funds(loaded):
    cfg=settings(loaded,enterprise_U_cash_fraction=1.)
    asset=simple_asset();base_flow=simple_flow(cash=30)
    loan_flow=simple_flow(cash=50,debt=20,local_U_debt=20)
    before=balance(base_flow,asset,cfg);drawn=balance(loan_flow,asset,cfg);repaid=balance(base_flow,asset,cfg)
    assert drawn['own_funds_equity']==before['own_funds_equity']
    assert transaction_money(drawn)-transaction_money(before)==pytest.approx(20)
    assert drawn['sectors']['banks']['assets_by_instrument']['loan']==20
    assert drawn['sectors']['enterprises']['liabilities_by_instrument']['loan']==20
    assert drawn['consolidated_off_earth']['net_worth']==pytest.approx(before['consolidated_off_earth']['net_worth'])
    assert repaid==before
    # Full settlement reserves require explicit external funding for this draw.
    fully_backed=balance(loan_flow,asset,settings(loaded,enterprise_U_cash_fraction=1.,bank_settlement_reserve_fraction=1.))
    assert fully_backed['bank_external_funding_requirement']==20
    assert fully_backed['sectors']['banks']['liabilities_by_instrument']['loan']==20
    assert fully_backed['sectors']['banks']['net_worth']==0
    assert fully_backed['consolidated_off_earth']['net_worth']==pytest.approx(before['consolidated_off_earth']['net_worth'])


def test_supplier_term_credit_is_not_outside_transaction_money(loaded):
    cfg=settings(loaded);asset=simple_asset()
    before=balance(simple_flow(),asset,cfg)
    after=balance(simple_flow(debt=20,supplier_U_debt=20),asset,cfg)
    assert transaction_money(after)==transaction_money(before)
    term=[c for c in after['claims'] if c['instrument']=='term_credit']
    assert len(term)==1 and term[0]['value']==20 and term[0]['money_tags']==['excluded_instrument']
    assert after['consolidated_off_earth']['net_worth']==pytest.approx(before['consolidated_off_earth']['net_worth'])


def test_finite_quote_and_future_capital_call_purchase_gate():
    flows=[dict(cash=0,reserves=0,debt=0,outside_equity_contribution=0,owner_payout=10) for _ in range(4)]
    assets=[simple_asset() for _ in flows]
    rate=.05;q=model.pricing(flows,assets,rate)
    assert q[0]['quote']==pytest.approx(sum(10/(1+rate)**y for y in (1,2,3))+100/(1+rate)**3)
    assert q[-1]['quote']==100
    flows[2]['outside_equity_contribution']=20
    quoted=model.pricing(flows,assets,rate)
    assert quoted[0]['quote'] is None and quoted[1]['quote'] is None
    assert quoted[0]['remaining_capital_calls']==20
    assert quoted[2]['quote'] is not None # current call is not a hidden future call


def test_escrow_release_belongs_to_its_original_owners_not_later_gift_recipients(loaded):
    cfg=settings(loaded);case,assets=synthetic_path()
    case['annual'][35].update(gross=1100.,dividend=1000.,owner_payout=990.,owner_custody=10.)
    case['annual'][36].update(owner_payout=10.)
    result=model.path(case,assets,cfg,'negotiated_vesting')
    at35=result['annual'][35];at36=result['annual'][36]
    old_fraction=at35['ownership_flow']['after_capital_dilution']
    assert at35['local_custody_claim']==pytest.approx(10*old_fraction)
    assert at36['flows']['local_payouts']==pytest.approx(10*old_fraction)
    assert at36['local_custody_claim']==pytest.approx(0,abs=1e-12)
    assert at36['flows']['local_payouts']<10*at36['local_economic_ownership_fraction']
    assert at35['flows']['local_payouts']+at35['flows']['outside_payouts']==pytest.approx(990)
    assert at36['flows']['local_payouts']+at36['flows']['outside_payouts']==pytest.approx(10)


def test_annual_loan_creation_and_extinction_close_money_identity(loaded):
    cfg=settings(loaded,enterprise_U_cash_fraction=1.)
    case,assets=synthetic_path()
    for y in range(20,30):case['annual'][y].update(cash=40.,debt=20.,local_U_debt=20.)
    case['annual'][20]['local_U_draw']=20.
    case['annual'][30]['local_U_principal']=20.
    result=model.path(case,assets,cfg,'external_equity')
    for y,amount in ((20,20.),(30,-20.)):
        money=result['annual'][y]['money_stock_flow']
        assert money['net_local_loan_creation']==amount
        assert money['outstanding_domestic_issuer_nonbank_money']-result['annual'][y-1]['money_stock_flow']['outstanding_domestic_issuer_nonbank_money']==pytest.approx(amount)
        assert money['reconciliation_residual']==pytest.approx(0,abs=1e-12)
        assert result['annual'][y]['book_balance']['consolidated_off_earth']['net_worth']==pytest.approx(0,abs=1e-12)


def test_reporting_basket_and_money_flow_identity_add_no_money(reference_paths):
    for path in reference_paths.values():
        for row in path['annual']:
            f=row['money_stock_flow'];m=f['outstanding_domestic_issuer_nonbank_money']
            assert math.fsum(row['basket_constituent_money_values'].values())==pytest.approx(m,rel=2e-15,abs=.1)
            assert set(row['basket_constituent_money_values'])=={'u1','u2','u3'}
            assert abs(f['reconciliation_residual'])<=max(1.,m*1e-11)
            assert row['book_balance']['controls']['passed']
            assert row['book_balance']['sectors']['settlement_issuer']['net_worth']==pytest.approx(0,abs=1)


def test_external_retention_changes_money_not_modelled_zone_wealth(reference):
    cfg,case,assets=reference
    low=model.path(case,assets,{**cfg,'external_investor_U_retention_fraction':0.},'external_equity')
    high=model.path(case,assets,{**cfg,'external_investor_U_retention_fraction':.05},'external_equity')
    for y in (32,100,200,220):
        a=low['annual'][y];b=high['annual'][y]
        increment=.05*a['cumulative']['outside_payouts']
        assert b['money_stock_flow']['outstanding_domestic_issuer_nonbank_money']-a['money_stock_flow']['outstanding_domestic_issuer_nonbank_money']==pytest.approx(increment,rel=1e-12,abs=.1)
        assert b['book_balance']['consolidated_off_earth']['net_worth']==pytest.approx(a['book_balance']['consolidated_off_earth']['net_worth'],rel=1e-12,abs=.1)
        assert b['programme_cash']==a['programme_cash'] and b['physical_assets']==a['physical_assets']


def test_all_frozen_finance_checkpoints_and_seed_dates_remain_unchanged(loaded):
    cfg,finance,_=loaded
    for name,expected in cfg['finance_dependencies'].items():
        assert run.sha(run.HERE/name)==expected
    assert len(finance['cases'])==18
    for case in finance['cases']:
        flows=case['annual'];assert [r['year'] for r in flows]==list(range(221))
        assert all(r['gross']==r['owner_payout']==0 for r in flows[:8])
        assert math.fsum(r['outside_equity_contribution'] for r in flows[:8])==pytest.approx(case['seed']+case['initial_reserve'],rel=1e-12,abs=.1)
        assert flows[0]['reserves']==case['initial_reserve']
        assert all(r['cash']==0 for r in flows[:7])
        for checkpoint in case['baseline_checkpoints']:
            h=checkpoint['horizon'];last=flows[h]
            assert math.fsum(r['owner_payout'] for r in flows[:h+1])==pytest.approx(checkpoint['owner_payouts'],rel=1e-12,abs=.1)
            assert math.fsum(r['outside_equity_contribution'] for r in flows[8:h+1])==pytest.approx(checkpoint['additional_equity_required'],rel=1e-12,abs=.1)
            for actual,expected in [('cash','terminal_cash'),('reserves','terminal_reserves'),('debt','terminal_debt'),('owner_custody','terminal_custody')]:
                assert last[actual]==pytest.approx(checkpoint[expected],rel=1e-12,abs=.1)
        for prev,row in zip(flows,flows[1:]):
            assert row['local_U_debt']-prev['local_U_debt']==pytest.approx(row['local_U_draw']-row['local_U_principal'],abs=1)
            assert row['supplier_U_debt']-prev['supplier_U_debt']==pytest.approx(row['supplier_U_draw']-row['supplier_U_principal'],abs=1)
    case=next(c for c in finance['cases'] if c['name']==cfg['reference_case'])
    assert math.fsum(r['gross'] for r in case['annual'][:201])==pytest.approx(1297978529975457.2,rel=1e-12)
    assert math.fsum(r['owner_payout'] for r in case['annual'][:201])==pytest.approx(347654892021511.75,rel=1e-12)


@pytest.mark.parametrize('weights',[{}, {'u1':.8,'u2':.3}, {'u1':1.1,'u2':-.1}, {'u1':float('nan')}])
def test_invalid_constituent_allocations_are_rejected(loaded,weights):
    cfg=settings(loaded,constituent_reporting_value_shares=weights);case,assets=synthetic_path()
    with pytest.raises(ValueError,match='Constituent'):
        model.path(case,assets,cfg,'external_equity')


def test_new_outside_capital_dilutes_existing_claims_and_unpriced_gifts_stay_unknown(loaded):
    cfg=settings(loaded,vesting_end_year=45,vesting_target_fraction=.5)
    case,assets=synthetic_path()
    for y in range(50,221):case['annual'][y]['cash']=40.
    case['annual'][50]['outside_equity_contribution']=20.
    result=model.path(case,assets,cfg,'negotiated_vesting')
    before=result['annual'][49];after=result['annual'][50]
    assert before['local_economic_ownership_fraction']==pytest.approx(.5)
    assert after['local_economic_ownership_fraction']==pytest.approx(.5*130/150)
    assert after['book_balance']['sectors']['community_fund']['net_worth']==pytest.approx(before['book_balance']['sectors']['community_fund']['net_worth'])
    assert after['book_balance']['consolidated_off_earth']['net_worth']==pytest.approx(before['book_balance']['consolidated_off_earth']['net_worth'])
    gifted=[r for r in result['annual'] if r['ownership_flow']['gifted_fraction']>0]
    assert len(gifted)==13 and all(r['flows']['gift_cash_claim_value'] is None for r in gifted)
    assert result['annual'][45]['cumulative']['unpriced_gift_count']==13
    assert result['annual'][220]['cumulative']['gift_cash_claim_value'] is None
    assert result['annual'][45]['cumulative']['unpriced_gift_book_value']==pytest.approx(65)
    assert all(r['flows']['purchase_spending']==0 for r in result['annual'])


def test_constituent_financial_views_preserve_every_sector_claim(reference_paths):
    for path in reference_paths.values():
        for row in path['annual']:
            view=row['basket_constituent_financial_claim_values']
            for actor,sector in row['book_balance']['sectors'].items():
                for field,total_field in [('assets','assets_by_denomination'),('liabilities','liabilities_by_denomination'),('net_financial_claims','net_financial_wealth_by_denomination')]:
                    assert sum(c[actor][field] for c in view.values())==pytest.approx(sector[total_field].get(model.U,0.),rel=2e-15,abs=.1)
            for currency,sectors in view.items():
                assert abs(sum(s['net_financial_claims'] for s in sectors.values()))<.1
            # Constituent views classify the same claims; equity valuation and
            # outside currencies remain in their separate denomination columns.
            assert row['book_balance']['money']['classification_complete']
