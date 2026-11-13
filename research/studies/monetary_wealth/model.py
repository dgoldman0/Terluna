"""Conditional monetary holdings and ownership of unchanged programme cash flows.

All stocks are reference values, not nominal U faces or exchange-rate forecasts.
"""
from __future__ import annotations
import math
from .ledger import snapshot
from ..complex_returns import consolidated as valuation
ACTORS={
 'enterprises':{'resident':True,'role':'enterprise'},
 'households':{'resident':True,'role':'household'},
 'community_fund':{'resident':True,'role':'community_fund'},
 'banks':{'resident':True,'role':'bank'},
 'settlement_issuer':{'resident':True,'role':'settlement_issuer'},
 'outside_owners_customers':{'resident':False,'role':'outside_nonbank'},
 'outside_banks':{'resident':False,'role':'outside_bank'}}
U='U_system_reference_value';FX='external_currency_value';EQUITY='equity_reference_value'

def close(a,b,scale=1.):return abs(a-b)<=max(1.,1e-11*max(abs(a),abs(b),scale))

def own_funds(flow,asset):
    return asset['net_remaining_service_value_USD2041']+asset['construction_WIP_value_USD2041']+flow['cash']+flow['reserves']-flow['debt']

def pricing(flows,assets,rate):
    """Finite-horizon negotiated quote; unresolved if future capital calls remain.

    The terminal own-funds mark is a valuation convention, not actual sale cash.
    """
    n=len(flows);q=[0.]*n;future_calls=[0.]*n
    q[-1]=own_funds(flows[-1],assets[-1]);calls=0.
    for y in range(n-2,-1,-1):
        calls+=flows[y+1]['outside_equity_contribution'];future_calls[y]=calls
        q[y]=(q[y+1]+flows[y+1]['owner_payout']-flows[y+1]['outside_equity_contribution'])/(1+rate)
    return [{'quote':q[y] if future_calls[y]<=1 and q[y]>0 else None,'remaining_capital_calls':future_calls[y],
             'terminal_own_funds_PV':q[-1]/(1+rate)**(n-1-y)} for y in range(n)]

def balance(flow,asset,local_share,household_cash,fund_cash,bank_equity,local_custody,outside_U,settings,equity_mark=None):
    """One complete stock matrix; monetary institutions own no free reserve gift."""
    q=own_funds(flow,asset) if equity_mark is None else equity_mark
    eU=settings['enterprise_U_cash_fraction']*flow['cash'];hU=settings['household_U_liquid_fraction']*household_cash;fU=settings['community_U_liquid_fraction']*fund_cash
    deposits=eU+hU+fU+outside_U;bank_reserve=settings['bank_settlement_reserve_fraction']*deposits
    need=deposits+flow['owner_custody']+bank_equity-flow['local_U_debt']-bank_reserve
    bank_fx=max(0.,need);bank_external_funding=max(0.,-need)
    claims=[]
    def add(issuer,holder,instrument,value,denomination=U,liquid=False):
        if abs(value)<1e-8:return
        claims.append(dict(issuer=issuer,holder=holder,instrument=instrument,denomination=denomination,value=value,liquid=liquid))
    for holder,value in [('enterprises',eU),('households',hU),('community_fund',fU),('outside_owners_customers',outside_U)]:add('banks',holder,'deposit',value,liquid=True)
    add('outside_banks','enterprises','deposit',flow['cash']-eU,FX,True)
    add('outside_banks','households','external_portfolio',household_cash-hU,FX)
    add('outside_banks','community_fund','external_portfolio',fund_cash-fU,FX)
    add('outside_banks','banks','deposit',bank_fx,FX,True)
    add('outside_banks','settlement_issuer','deposit',flow['reserves']+bank_reserve,FX,True)
    add('settlement_issuer','enterprises','restricted_fund',flow['reserves'])
    add('settlement_issuer','banks','bank_reserve',bank_reserve,liquid=True)
    add('enterprises','banks','loan',flow['local_U_debt'])
    add('enterprises','outside_owners_customers','term_credit',flow['supplier_U_debt'])
    add('enterprises','outside_banks','loan',flow['external_debt'],FX)
    add('banks','outside_banks','loan',bank_external_funding,FX)
    add('banks','community_fund','custody',local_custody)
    add('banks','outside_owners_customers','custody',flow['owner_custody']-local_custody)
    add('banks','community_fund','equity',bank_equity,EQUITY)
    add('enterprises','community_fund','equity',local_share*q,EQUITY)
    add('enterprises','outside_owners_customers','equity',(1-local_share)*q,EQUITY)
    fixed=asset['net_remaining_service_value_USD2041']+asset['construction_WIP_value_USD2041']
    result=snapshot(ACTORS,{'enterprises':fixed},claims)
    if not result['controls']['passed']:raise ValueError('Counterparty balance failed')
    if equity_mark is None and abs(result['sectors']['enterprises']['net_worth'])>1:raise ValueError('Book equity double count')
    for actor in ('banks','settlement_issuer'):
        if abs(result['sectors'][actor]['net_worth'])>1:raise ValueError('Financial intermediary imbalance')
    for actor in ('outside_owners_customers','outside_banks'):
        result['sectors'][actor]['scope']='Partial off-Earth counterparty positions; not complete outside wealth.'
    result['own_funds_equity']=own_funds(flow,asset);result['equity_mark']=q
    result['bank_external_funding_requirement']=bank_external_funding
    result['monetary_external_assets']=bank_fx+flow['reserves']+bank_reserve
    return result

def compact_balance(b):
    return {k:b[k] for k in ('sectors','consolidated_off_earth','money','own_funds_equity','equity_mark','bank_external_funding_requirement','monetary_external_assets','controls')}

def path(case,asset_case,cfg,ownership,*,purchase_basis='cash_claim'):
    weights=cfg['constituent_reporting_value_shares']
    if not weights or any(not math.isfinite(w) or w<0 for w in weights.values()) or abs(sum(weights.values())-1)>1e-12:
        raise ValueError('Constituent value shares must be finite, nonnegative and sum to one')
    if ownership not in cfg['ownership_scenarios'] or purchase_basis not in ('cash_claim','own_funds'):
        raise ValueError('Unknown ownership or purchase convention')
    flows=case['annual'];assets=[a['totals'] for a in asset_case['annual']]
    if [f['year'] for f in flows]!=list(range(221)) or [a['year'] for a in asset_case['annual']]!=list(range(221)):raise ValueError('Annual boundary mismatch')
    price=pricing(flows,assets,cfg['purchase_discount_rate'])
    household=fund=share=bank_equity=local_custody=0.;outside_received=0.;rows=[]
    cumulative={k:0. for k in ('gross_receipts','outside_contributions','outside_payouts','local_payouts','outside_buyout_receipts','household_payroll','household_saving','community_contributions','purchase_spending','gift_book_value','gift_cash_claim_value','unpriced_gift_count','unpriced_gift_book_value','local_U_created','local_U_extinguished','supplier_U_issued','supplier_U_repaid','bank_retained_interest','outside_consumption')}
    old_money=old_mon_fx=old_reserve=old_custody=old_bank_equity=old_bank_external=0.
    old_nw=old_real=old_niip=0.;outside_cash=[]
    for y,(flow,a) in enumerate(zip(flows,assets)):
        q=own_funds(flow,a)
        if q<=0:raise ValueError('Nonpositive programme own funds requires separate insolvency treatment')
        opening_share=share;contribution=flow['outside_equity_contribution']
        # New outside funds buy the required recapitalization at the closing own-
        # funds proxy. Existing claims are diluted; no local capital call is hidden.
        if y>=8 and contribution>1:
            pre_money=q-contribution
            if pre_money<0:raise ValueError('Negative pre-money recapitalization requires a new loss-allocation case')
            share*=pre_money/q
        diluted_share=share
        intended=cfg['owner_custody_fraction']*flow['dividend'];opening_custody=flows[y-1]['owner_custody'] if y else 0.
        custody_pool=opening_custody+intended;local_pool=local_custody+intended*share
        release=custody_pool-flow['owner_custody']
        if release<-1 or release>custody_pool+1:raise ValueError('Invalid custody release')
        local_release=release*local_pool/custody_pool if custody_pool>0 else 0.
        local_custody=local_pool-local_release
        if local_custody < -1 or local_custody > flow['owner_custody']+1:
            raise ValueError('Local custody exceeds the outstanding custody account')
        # Sub-dollar cancellation at trillion-dollar scale is rounding only.
        local_custody=min(flow['owner_custody'],max(0.,local_custody))
        local_payout=(flow['dividend']-intended)*share+local_release
        outside_payout=flow['owner_payout']-local_payout
        if min(local_payout,outside_payout,local_custody)<-1:raise ValueError('Negative owner allocation')
        saving=cfg['household_net_saving_fraction']*flow['resident_payroll']
        donation=saving*cfg['household_savings_to_purchase_fund_fraction'] if ownership=='community_purchase' else 0.
        household+=saving-donation;fund+=donation+local_payout
        bank_profit=flow['local_U_interest']*(1-cfg['bank_operating_cost_fraction_of_local_interest']);bank_equity+=bank_profit
        purchase=acquired=gift=gift_book=gift_quote=0.
        executable=price[y]['quote'] is not None and y>=cfg['local_purchase_start_year'] and y<220
        quote=q if purchase_basis=='own_funds' else price[y]['quote']
        if ownership=='community_purchase' and executable:
            purchase=min(fund,(1-share)*quote);acquired=purchase/quote;share+=acquired;fund-=purchase
        if ownership=='negotiated_vesting' and cfg['local_purchase_start_year']<y<=cfg['vesting_end_year']:
            target=cfg['vesting_target_fraction']*(y-cfg['local_purchase_start_year'])/(cfg['vesting_end_year']-cfg['local_purchase_start_year'])
            gift=min(1-share,max(0.,target-share));share+=gift;gift_book=gift*q
            gift_quote=gift*price[y]['quote'] if price[y]['quote'] is not None else (None if gift>0 else 0.)
        if not -1e-12<=share<=1+1e-12:raise ValueError('Ownership outside unit interval')
        share=min(1.,max(0.,share));outside_received+=outside_payout+purchase
        elapsed=max(0,y-8);settlement_share=cfg['U_settlement_start_fraction']+(cfg['U_settlement_end_fraction']-cfg['U_settlement_start_fraction'])*min(1,elapsed/cfg['U_settlement_transition_years'])
        customer_U=settlement_share*flow['gross']*cfg['outside_customer_holding_days']/365
        investor_U=cfg['external_investor_U_retention_fraction']*outside_received
        b=balance(flow,a,share,household,fund,bank_equity,local_custody,customer_U+investor_U,cfg)
        m=b['money']['resident_transaction_money']+b['money']['external_transaction_money']
        dfx=b['monetary_external_assets']-old_mon_fx;dr=flow['reserves']-old_reserve;dc=flow['owner_custody']-old_custody;de=bank_equity-old_bank_equity;db=b['bank_external_funding_requirement']-old_bank_external
        loan_net=flow['local_U_draw']-flow['local_U_principal']
        money_residual=(m-old_money)-(loan_net+dfx-dr-dc-de-db)
        if not close(money_residual,0.,m):raise ValueError('Money stock-flow identity failed')
        old_money,old_mon_fx,old_reserve,old_custody,old_bank_equity,old_bank_external=m,b['monetary_external_assets'],flow['reserves'],flow['owner_custody'],bank_equity,b['bank_external_funding_requirement']
        adds=dict(gross_receipts=flow['gross'],outside_contributions=contribution,outside_payouts=outside_payout,local_payouts=local_payout,outside_buyout_receipts=purchase,household_payroll=flow['resident_payroll'],household_saving=saving,community_contributions=donation,purchase_spending=purchase,gift_book_value=gift_book,gift_cash_claim_value=gift_quote,unpriced_gift_count=int(gift_quote is None),unpriced_gift_book_value=gift_book if gift_quote is None else 0.,local_U_created=flow['local_U_draw'],local_U_extinguished=flow['local_U_principal'],supplier_U_issued=flow['supplier_U_draw'],supplier_U_repaid=flow['supplier_U_principal'],bank_retained_interest=bank_profit,outside_consumption=flow['resident_payroll']-saving)
        for k,v in adds.items():
            cumulative[k]=None if v is None or cumulative[k] is None else cumulative[k]+v
        marked=None
        if price[y]['quote'] is not None:
            alt=balance(flow,a,share,household,fund,bank_equity,local_custody,customer_U+investor_U,cfg,equity_mark=price[y]['quote'])
            marked={k:alt[k] for k in ('sectors','consolidated_off_earth','equity_mark')}
        row={'year':y,'population':flow.get('population',flow.get('population_convention')),'physical_assets':a,'physical_by_function':asset_case['annual'][y]['by_functional_sector'],'programme_cash':flow['cash'],'programme_reserves':flow['reserves'],'programme_debt':flow['debt'],'programme_local_U_debt':flow['local_U_debt'],'programme_supplier_U_debt':flow['supplier_U_debt'],'programme_external_debt':flow['external_debt'],'household_retained_financial_assets':household,'community_uninvested_financial_assets':fund,'local_economic_ownership_fraction':share,'outside_economic_ownership_fraction':1-share,'local_custody_claim':local_custody,'ownership_flow':{'opening_fraction':opening_share,'after_capital_dilution':diluted_share,'purchased_fraction':acquired,'gifted_fraction':gift,'purchase_price_basis':purchase_basis,'purchase_price':quote,'purchase_gate_satisfied':executable,'purchase_premium_above_own_funds':purchase-acquired*q,'gift_book_value':gift_book,'gift_cash_claim_value':gift_quote},'flows':adds,'cumulative':dict(cumulative),'book_balance':compact_balance(b),'cash_claim_balance':marked,'equity_quote':price[y],'money_stock_flow':{'outstanding_domestic_issuer_nonbank_money':m,'net_local_loan_creation':loan_net,'net_external_backing_acquisition':dfx,'restricted_reserve_funding':dr,'custody_funding':dc,'retained_bank_income':de,'new_external_bank_funding':db,'reconciliation_residual':money_residual},'outside_U_holders':{'customers':customer_U,'investors':investor_U},'basket_constituent_money_values':{k:v*m for k,v in cfg['constituent_reporting_value_shares'].items()}}
        row['basket_constituent_financial_claim_values']={currency:{actor:{
            'assets':weight*sector['assets_by_denomination'].get(U,0.),
            'liabilities':weight*sector['liabilities_by_denomination'].get(U,0.),
            'net_financial_claims':weight*sector['net_financial_wealth_by_denomination'].get(U,0.)}
            for actor,sector in b['sectors'].items()} for currency,weight in weights.items()}
        n=b['consolidated_off_earth'];real=n['real_assets'];niip=n['net_international_investment_position'];nw=n['net_worth']
        residual=(nw-old_nw)-(real-old_real)-(niip-old_niip)
        payout_residual=local_payout+outside_payout-flow['owner_payout']
        saving_residual=household+fund-(cumulative['household_saving']+cumulative['local_payouts']-cumulative['purchase_spending'])
        if not all(close(v,0.,max(abs(nw),m)) for v in (residual,payout_residual,saving_residual)):
            raise ValueError('Wealth/ownership flow conservation failed')
        row['wealth_change']={'real_asset_change':real-old_real,'net_external_position_change':niip-old_niip,'net_worth_change':nw-old_nw,'identity_residual':residual}
        row['allocation_controls']={'payout_residual':payout_residual,'saving_and_purchase_residual':saving_residual}
        old_nw,old_real,old_niip=nw,real,niip
        outside_cash.append(outside_payout+purchase-contribution)
        if y in cfg['checkpoints']:
            row['claims']=b['claims']
            row['outside_investor_cash_returns']={'dated_cash_NPV_at_3percent':valuation.pv(outside_cash,.03),'nonnegative_IRR_if_unique':valuation.nonnegative_irr(outside_cash),'durable_discounted_cash_recovery_year_at_3percent':valuation.durable_recovery(outside_cash,.03),'net_cash_received':sum(outside_cash),'remaining_equity_at_own_funds':(1-share)*q,'terminal_equity_added_to_cash_returns':False}

        rows.append(row)
    return {'name':case['name'],'portfolio':case['portfolio'],'growth':case['growth'],'ownership':ownership,'purchase_basis':purchase_basis,'annual':rows,'checkpoints':{str(y):rows[y] for y in cfg['checkpoints']},'baseline_scope':case['seed_cost_scope'],'physical_supply_qualified':False}
