"""Counterparty conservation, ownership transfers and monetary boundaries."""
from copy import deepcopy
import math
import pytest
from research.studies.monetary_wealth.ledger import snapshot


def actors():
    return {'enterprise':{'resident':True,'role':'enterprise'},
            'households':{'resident':True,'role':'household'},
            'community_fund':{'resident':True,'role':'community_fund'},
            'bank':{'resident':True,'role':'bank'},
            'issuer':{'resident':True,'role':'settlement_issuer'},
            'outside':{'resident':False,'role':'outside_nonbank'},
            'outside_bank':{'resident':False,'role':'outside_bank'}}


def claim(issuer,holder,instrument,value,denomination='U',liquid=False,**extra):
    return dict(issuer=issuer,holder=holder,instrument=instrument,value=value,
                denomination=denomination,liquid=liquid,**extra)


def check(result):
    assert result['controls']['passed']
    assert all(abs(result['controls'][k])<1e-9 for k in
               ('global_financial_net_zero','internal_cancellation_residual',
                'off_earth_net_worth_residual','global_net_worth_minus_real_assets',
                'maximum_claim_pair_residual'))


def test_domestic_credit_creates_deposit_and_loan_without_net_wealth():
    r=snapshot(actors(),{},[claim('enterprise','bank','loan',20),
                            claim('bank','enterprise','deposit',20,liquid=True)])
    check(r)
    assert r['money']['resident_transaction_money']==20
    assert r['money']['monetary_base']==0
    assert all(s['net_worth']==0 for s in r['sectors'].values())
    assert r['consolidated_off_earth']['net_worth']==0
    assert r['consolidated_off_earth']['internal_financial_assets_eliminated']==40


def test_bank_reserves_do_not_duplicate_nonbank_money():
    r=snapshot(actors(),{},[claim('enterprise','bank','loan',20),
                            claim('bank','households','deposit',20,liquid=True),
                            claim('issuer','bank','bank_reserve',10,liquid=True),
                            claim('issuer','outside','currency',2,liquid=True)])
    check(r)
    assert r['money']['resident_transaction_money']==20
    assert r['money']['external_transaction_money']==2
    assert r['money']['resident_bank_reserves']==10
    assert r['money']['monetary_base']==12
    assert r['claims'][2]['money_tags']==['resident_bank_reserves','monetary_base']


def test_bank_vault_currency_is_base_but_not_nonbank_money():
    r=snapshot(actors(),{},[claim('issuer','bank','currency',4,liquid=True)])
    assert r['money']['monetary_base']==4 and r['money']['resident_transaction_money']==0
    assert r['money']['resident_bank_currency']==4


def test_direct_issuer_deposit_and_reserve_alias_have_same_base_boundary():
    r=snapshot(actors(),{},[claim('issuer','households','deposit',6,liquid=True),
                            claim('issuer','bank','deposit',4,liquid=True)])
    assert r['money']['resident_transaction_money']==6
    assert r['money']['resident_bank_reserves']==4
    assert r['money']['monetary_base']==10


def test_external_equity_offsets_territorial_productive_asset():
    r=snapshot(actors(),{'enterprise':100},[claim('enterprise','outside','equity',100)])
    check(r)
    assert r['sectors']['enterprise']['net_worth']==0
    assert r['consolidated_off_earth']['real_assets']==100
    assert r['consolidated_off_earth']['external_liabilities_by_instrument']=={'equity':100}
    assert r['consolidated_off_earth']['net_international_investment_position']==-100
    assert r['consolidated_off_earth']['net_worth']==0
    assert r['sectors']['outside']['net_worth']==100


def test_fair_value_buyout_swaps_assets_without_creating_wealth():
    before=snapshot(actors(),{'enterprise':100},[
        claim('enterprise','outside','equity',100),
        claim('outside_bank','households','deposit',20,liquid=True)])
    after=snapshot(actors(),{'enterprise':100},[
        claim('enterprise','outside','equity',80),
        claim('enterprise','households','equity',20),
        claim('outside_bank','outside','deposit',20,liquid=True)])
    check(before);check(after)
    assert before['consolidated_off_earth']['net_worth']==after['consolidated_off_earth']['net_worth']==20
    assert before['sectors']['households']['net_worth']==after['sectors']['households']['net_worth']==20
    assert before['sectors']['outside']['net_worth']==after['sectors']['outside']['net_worth']==100
    assert before['money']['resident_transaction_money']==0
    assert before['money']['resident_external_issuer_liquid_assets']==20


def test_gift_transfers_net_wealth_without_creating_productive_assets():
    before=snapshot(actors(),{'enterprise':100},[claim('enterprise','outside','equity',100)])
    after=snapshot(actors(),{'enterprise':100},[
        claim('enterprise','outside','equity',75),claim('enterprise','community_fund','equity',25)])
    check(before);check(after)
    assert after['consolidated_off_earth']['net_worth']-before['consolidated_off_earth']['net_worth']==25
    assert after['sectors']['outside']['net_worth']-before['sectors']['outside']['net_worth']==-25
    assert before['consolidated_off_earth']['real_assets']==after['consolidated_off_earth']['real_assets']==100
    assert 'GDP' not in after and 'output' not in after


def test_local_fund_with_external_owners_does_not_create_local_wealth():
    r=snapshot(actors(),{'enterprise':100},[
        claim('enterprise','community_fund','equity',100),
        claim('community_fund','outside','equity',100)])
    check(r)
    assert r['consolidated_off_earth']['net_worth']==0
    assert r['sectors']['community_fund']['net_worth']==0


def test_signed_equity_is_retained_and_other_negative_claims_are_rejected():
    r=snapshot(actors(),{},[claim('enterprise','outside','loan',20),
                            claim('enterprise','outside','equity',-20)])
    check(r)
    assert r['sectors']['outside']['assets_by_instrument']['equity']==-20
    assert r['consolidated_off_earth']['external_liabilities']==0
    with pytest.raises(ValueError,match='Only equity'):
        snapshot(actors(),{},[claim('enterprise','outside','loan',-20)])


def test_restricted_custody_and_term_claims_are_not_transaction_money():
    r=snapshot(actors(),{},[
        claim('bank','enterprise','deposit',10,liquid=True,restricted=True),
        claim('bank','households','custody',20,liquid=True),
        claim('enterprise','outside','term_credit',30,liquid=True),
        claim('bank','community_fund','restricted_fund',40,liquid=True)])
    check(r)
    assert r['money']['resident_transaction_money']==r['money']['external_transaction_money']==0
    assert r['money']['excluded_restricted_or_illiquid']==10
    assert r['money']['excluded_instrument']==90
    assert all(x['money_reason'] for x in r['claims'])


def test_basket_denomination_does_not_create_an_issuer_or_claim():
    contracts=[claim('bank','households','deposit',25,denomination='basket_U',liquid=True)]
    a=snapshot(actors(),{},contracts)
    contracts[0]['denomination']='constituent_UA'
    b=snapshot(actors(),{},contracts)
    breakdowns={'assets_by_denomination','liabilities_by_denomination',
                'net_financial_wealth_by_denomination'}
    for actor in actors():
        assert {k:v for k,v in a['sectors'][actor].items() if k not in breakdowns}=={
            k:v for k,v in b['sectors'][actor].items() if k not in breakdowns}
    assert a['money']==b['money']
    assert a['consolidated_off_earth']==b['consolidated_off_earth']
    assert a['sectors']['households']['assets_by_denomination']=={'basket_U':25}
    assert b['sectors']['households']['assets_by_denomination']=={'constituent_UA':25}
    assert len(a['claims'])==1 and len(a['sectors'])==len(actors())


def test_boolean_actors_reconcile_without_inventing_money_classification():
    r=snapshot({'firm':True,'lender':True},{},[
        claim('firm','lender','loan',3),claim('lender','firm','deposit',3,liquid=True)])
    check(r)
    assert not r['money']['classification_complete']
    assert r['money']['unclassified_resident']==3


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),True,'10'])
def test_invalid_reference_values_are_rejected(bad):
    with pytest.raises(ValueError,match='finite reference'):
        snapshot(actors(),{},[claim('enterprise','outside','equity',bad)])


def test_unmatched_roles_self_claims_and_unknown_actors_are_rejected():
    with pytest.raises(ValueError,match='bank_reserve'):
        snapshot(actors(),{},[claim('enterprise','households','bank_reserve',2,liquid=True)])
    with pytest.raises(ValueError,match='Self-held'):
        snapshot(actors(),{},[claim('enterprise','enterprise','equity',1)])
    with pytest.raises(ValueError,match='Unknown claim'):
        snapshot(actors(),{},[claim('missing','outside','equity',1)])


def test_input_records_are_not_mutated_and_instrument_totals_remain_visible():
    contracts=[claim('enterprise','outside','loan',5),claim('enterprise','outside','equity',7)]
    original=deepcopy(contracts)
    r=snapshot(actors(),{'enterprise':12},contracts)
    assert contracts==original
    assert r['sectors']['enterprise']['liabilities_by_instrument']=={'equity':7.,'loan':5.}
    assert len(r['claims'])==2 and r['controls']['passed']


def test_restricted_reserve_custody_preserves_beneficial_ownership_once():
    r=snapshot(actors(),{},[
        claim('outside_bank','bank','deposit',10,liquid=True),
        claim('bank','enterprise','restricted_fund',10),
        claim('enterprise','community_fund','equity',10)])
    check(r)
    assert r['sectors']['bank']['net_worth']==0
    assert r['sectors']['enterprise']['net_worth']==0
    assert r['sectors']['community_fund']['net_worth']==10
    assert r['consolidated_off_earth']['external_financial_assets']==10
    assert r['consolidated_off_earth']['net_worth']==10
    assert r['money']['resident_transaction_money']==0
    assert r['money']['monetary_base']==0


def test_denomination_breakdowns_reconcile_by_sector_and_counterparty():
    r=snapshot(actors(),{'enterprise':12},[
        claim('bank','households','deposit',20,denomination='U_reference_value',liquid=True),
        claim('outside_bank','households','external_portfolio',30,denomination='external_reference_value'),
        claim('enterprise','households','equity',7,denomination='equity_reference_value'),
        claim('households','bank','loan',5,denomination='U_reference_value'),
        claim('enterprise','outside','equity',-2,denomination='equity_reference_value')])
    check(r)
    household=r['sectors']['households']
    assert household['assets_by_denomination']=={
        'U_reference_value':20,'external_reference_value':30,'equity_reference_value':7}
    assert household['liabilities_by_denomination']=={'U_reference_value':5}
    assert household['net_financial_wealth_by_denomination']=={
        'U_reference_value':15,'external_reference_value':30,'equity_reference_value':7}
    assert r['sectors']['outside']['assets_by_denomination']=={'equity_reference_value':-2}
    for sector in r['sectors'].values():
        for breakdown,total in [('assets_by_denomination','financial_assets'),
                                ('liabilities_by_denomination','financial_liabilities'),
                                ('net_financial_wealth_by_denomination','net_financial_wealth')]:
            assert math.fsum(sector[breakdown].values())==pytest.approx(sector[total])
    for denomination in ('U_reference_value','external_reference_value','equity_reference_value'):
        assert math.fsum(sector['net_financial_wealth_by_denomination'].get(denomination,0)
                         for sector in r['sectors'].values())==0
    assert r['money']['resident_transaction_money']==20


def test_equity_valuation_denomination_does_not_make_equity_money():
    r=snapshot(actors(),{'enterprise':100},[
        claim('enterprise','community_fund','equity',100,
              denomination='equity_reference_value',liquid=True)])
    check(r)
    assert r['sectors']['community_fund']['net_financial_wealth_by_denomination']=={
        'equity_reference_value':100}
    assert r['consolidated_off_earth']['net_worth']==100
    assert r['money']['resident_transaction_money']==0
    assert r['money']['monetary_base']==0
    assert len(r['claims'])==1
