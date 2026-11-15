"""Counterparty balance sheets for the enclosed-zone wealth comparison.

Values are already expressed in constant USD2041 reference valuation. Currency
labels are metadata: this kernel creates no exchange rate, basket issuer, output,
reserve-backing promise or ownership contract. Equity is an issuer liability.
An own-funds/book-equity proxy must be selected and supplied by the caller.
"""
from __future__ import annotations
from collections import defaultdict
from copy import deepcopy
import math
from numbers import Real

BANK_ROLES={'bank','banks','outside_bank'}
ISSUER_ROLES={'settlement_issuer','monetary_authority'}
NONBANK_ROLES={'enterprise','enterprises','household','households','community_fund',
               'nonprofit','government','outside_nonbank','nonbank'}
UNCLASSIFIED_ROLES={'unspecified','outside','external_aggregate'}
MONEY_INSTRUMENTS={'deposit','currency','bank_reserve'}


def _value(value,*,signed=False):
    if isinstance(value,bool) or not isinstance(value,Real) or not math.isfinite(value):
        raise ValueError('A finite reference value is required')
    if not signed and value<0:
        raise ValueError('Only equity claims can have negative values')
    return float(value)


def _actors(actors):
    if not isinstance(actors,dict) or not actors:
        raise ValueError('Named actors are required')
    result={}
    for name,item in actors.items():
        if not isinstance(name,str) or not name.strip():raise ValueError('Invalid actor name')
        if isinstance(item,bool):item={'resident':item,'role':'unspecified'}
        if not isinstance(item,dict) or not isinstance(item.get('resident'),bool):
            raise ValueError('Each actor needs a resident boolean')
        role=item.get('role','unspecified')
        if role not in BANK_ROLES|ISSUER_ROLES|NONBANK_ROLES|UNCLASSIFIED_ROLES:
            raise ValueError('Unknown institutional role: '+str(role))
        result[name]={'resident':item['resident'],'role':role}
    return result


def _money_tags(record,actors):
    """Return explicit instrument eligibility and separate base-money tags."""
    instrument=record['instrument'];issuer=actors[record['issuer']];holder=actors[record['holder']]
    residence='resident' if holder['resident'] else 'external'
    ir,hr=issuer['role'],holder['role'];tags=[]
    if instrument not in MONEY_INSTRUMENTS:
        return ['excluded_instrument'],'Instrument is outside the transaction-money and base-money definitions.'
    if not record['liquid'] or record.get('restricted',False):
        return ['excluded_restricted_or_illiquid'],'Restricted or illiquid balances are excluded from money aggregates.'
    if ir in UNCLASSIFIED_ROLES or hr in UNCLASSIFIED_ROLES:
        return ['unclassified_'+residence],'Issuer or holder role is unspecified; no complete monetary aggregate is inferred.'
    if instrument=='bank_reserve':
        if ir not in ISSUER_ROLES or hr not in BANK_ROLES:
            raise ValueError('A bank_reserve claim requires a settlement issuer and a bank holder')
        if issuer['resident']:
            tags.append(residence+'_bank_reserves')
            if holder['resident']:tags.append('monetary_base')
        else:tags.append(residence+'_bank_reserves_at_external_issuer')
        return tags,'Bank settlement reserves are separate from nonbank transaction money.'
    if instrument=='currency' and ir not in ISSUER_ROLES:
        raise ValueError('Currency requires a settlement-issuer role')
    if instrument=='deposit' and ir not in BANK_ROLES|ISSUER_ROLES:
        raise ValueError('Deposits require a bank or settlement-issuer role')
    if not issuer['resident']:
        return [residence+'_external_issuer_liquid_assets'],'Outside-issued liquid claims are reported separately from domestic-issuer money.'
    if hr in NONBANK_ROLES:
        tags.extend([residence+'_nonbank_'+instrument,residence+'_transaction_money'])
        if instrument=='deposit' and ir in ISSUER_ROLES and holder['resident']:
            tags.append('monetary_base')
    elif hr in BANK_ROLES:
        tags.append(residence+'_bank_'+instrument)
        # A bank deposit at the domestic settlement issuer is equivalent to a
        # settlement reserve, although the explicit bank_reserve label is preferred.
        if instrument=='deposit' and ir in ISSUER_ROLES:
            tags.append(residence+'_bank_reserves')
            if holder['resident']:tags.append('monetary_base')
    else:tags.append(residence+'_issuer_held_claims')
    if instrument=='currency':
        tags.append('currency_outside_issuer')
        tags.append('monetary_base')
    return tags,'Liquid domestic-issuer transaction claims, classified by holder; base-money tags are not additional wealth.'


def snapshot(actors,real_assets,claims):
    """Return sector and consolidated stocks with explicit counterpart controls.

    actors maps names either to resident booleans or to {resident, role}. Boolean
    actors reconcile wealth but leave monetary classification unqualified. Each
    claim has issuer, holder, instrument, denomination, value and liquid fields;
    optional restricted=True excludes it from money. Signed values are allowed
    only for instrument='equity'. Amounts are values, never nominal faces.
    """
    actors=_actors(actors)
    if not isinstance(real_assets,dict) or set(real_assets)-set(actors):
        raise ValueError('Real assets require known owners')
    if not isinstance(claims,(list,tuple)):raise ValueError('Claims must be a list of paired contracts')
    real={name:_value(real_assets.get(name,0)) for name in actors}
    asset_values={a:[] for a in actors};liability_values={a:[] for a in actors}
    assets_by={a:defaultdict(list) for a in actors};liabilities_by={a:defaultdict(list) for a in actors}
    denomination_assets={a:defaultdict(list) for a in actors}
    denomination_liabilities={a:defaultdict(list) for a in actors}
    external_assets=[];external_liabilities=[];internal_assets=[];internal_liabilities=[]
    monetary=defaultdict(list);records=[]
    for position,item in enumerate(claims):
        if not isinstance(item,dict):raise ValueError('A claim must be a record')
        missing={'issuer','holder','instrument','denomination','value','liquid'}-item.keys()
        if missing:raise ValueError('Missing claim fields: '+', '.join(sorted(missing)))
        issuer,holder=item['issuer'],item['holder']
        if issuer not in actors or holder not in actors:raise ValueError('Unknown claim counterparty')
        if issuer==holder:raise ValueError('Self-held claims require cancellation before entry')
        instrument,denomination=item['instrument'],item['denomination']
        if not isinstance(instrument,str) or not instrument.strip() or not isinstance(denomination,str) or not denomination.strip():
            raise ValueError('Instrument and denomination must be named')
        if not isinstance(item['liquid'],bool) or not isinstance(item.get('restricted',False),bool):
            raise ValueError('Liquid and restricted flags must be boolean')
        value=_value(item['value'],signed=instrument=='equity')
        record=deepcopy(item);record['value']=value
        tags,reason=_money_tags(record,actors)
        asset_values[holder].append(value);liability_values[issuer].append(value)
        assets_by[holder][instrument].append(value);liabilities_by[issuer][instrument].append(value)
        denomination_assets[holder][denomination].append(value)
        denomination_liabilities[issuer][denomination].append(value)
        ir,hr=actors[issuer]['resident'],actors[holder]['resident']
        if hr and not ir:external_assets.append(value)
        if ir and not hr:external_liabilities.append(value)
        if ir and hr:internal_assets.append(value);internal_liabilities.append(value)
        for tag in tags:monetary[tag].append(value)
        record.update(claim_index=position,issuer_resident=ir,holder_resident=hr,
                      asset_value=value,liability_value=value,money_tags=tags,money_reason=reason)
        records.append(record)
    sectors={}
    for name,metadata in actors.items():
        financial_assets=math.fsum(asset_values[name]);financial_liabilities=math.fsum(liability_values[name])
        nfw=math.fsum(asset_values[name]+[-x for x in liability_values[name]])
        da,dl=denomination_assets[name],denomination_liabilities[name]
        sectors[name]={**metadata,'financial_assets':financial_assets,'financial_liabilities':financial_liabilities,
                       'net_financial_wealth':nfw,'real_assets':real[name],
                       'net_worth':math.fsum([real[name],nfw]),
                       'assets_by_instrument':{k:math.fsum(v) for k,v in sorted(assets_by[name].items())},
                       'liabilities_by_instrument':{k:math.fsum(v) for k,v in sorted(liabilities_by[name].items())},
                       'assets_by_denomination':{k:math.fsum(v) for k,v in sorted(da.items())},
                       'liabilities_by_denomination':{k:math.fsum(v) for k,v in sorted(dl.items())},
                       'net_financial_wealth_by_denomination':{
                           k:math.fsum(da.get(k,[])+[-v for v in dl.get(k,[])])
                           for k in sorted(da.keys()|dl.keys())}}
    resident_real=math.fsum(real[a] for a in actors if actors[a]['resident'])
    ext_assets=math.fsum(external_assets);ext_liabilities=math.fsum(external_liabilities)
    niip=math.fsum(external_assets+[-x for x in external_liabilities])
    consolidated={'real_assets':resident_real,'external_financial_assets':ext_assets,
                  'external_liabilities':ext_liabilities,'net_international_investment_position':niip,
                  'net_worth':math.fsum([resident_real,niip]),
                  'internal_financial_assets_eliminated':math.fsum(internal_assets),
                  'internal_financial_liabilities_eliminated':math.fsum(internal_liabilities),
                  'external_liabilities_by_instrument':{k:math.fsum(r['value'] for r in records if r['issuer_resident'] and not r['holder_resident'] and r['instrument']==k) for k in sorted({r['instrument'] for r in records if r['issuer_resident'] and not r['holder_resident']})}}
    money={k:math.fsum(v) for k,v in sorted(monetary.items())}
    for key in ('resident_transaction_money','external_transaction_money','resident_bank_reserves','external_bank_reserves','currency_outside_issuer','monetary_base','unclassified_resident','unclassified_external'):
        money.setdefault(key,0.)
    # Non-negative monetary claims make these missing-role checks unambiguous.
    money['classification_complete']=not any(r['value'] and any(t.startswith('unclassified_') for t in r['money_tags']) for r in records)
    money['rule']='Domestic-issuer liquid currency/deposits held by nonbanks form resident/external transaction-money proxies. Monetary base includes domestic-issuer currency outside its issuer (including known external holdings), resident-bank settlement reserves and resident-nonbank direct issuer deposits. Outside-issued liquid assets, restricted/custody/term/equity claims are separate. Base and transaction money overlap and must not be added.'
    controls={'global_financial_net_zero':math.fsum(s['net_financial_wealth'] for s in sectors.values()),
              'internal_cancellation_residual':math.fsum(internal_assets+[-x for x in internal_liabilities]),
              'off_earth_net_worth_residual':math.fsum([s['net_worth'] for s in sectors.values() if s['resident']]+[-consolidated['net_worth']]),
              'global_net_worth_minus_real_assets':math.fsum([s['net_worth'] for s in sectors.values()]+[-x for x in real.values()]),
              'maximum_claim_pair_residual':max((abs(r['asset_value']-r['liability_value']) for r in records),default=0.)}
    scale=max(1.,math.fsum(abs(r['value']) for r in records),math.fsum(real.values()))
    controls['absolute_tolerance']=1e-10*scale
    controls['passed']=all(abs(v)<=controls['absolute_tolerance'] for k,v in controls.items() if k!='absolute_tolerance')
    return {'schema':'terluna.humanities.monetary-wealth-ledger/1',
            'valuation':'Constant USD2041 reference values; denominations are metadata, with no implicit U/USD peg or basket issuance.',
            'equity_rule':'Issued equity is a liability paired with each owner asset. A supplied own-funds/book-equity proxy is limited valuation evidence, not added a second time as corporate wealth. Negative equity is preserved; no legal loss cap is inferred.',
            'sectors':sectors,'consolidated_off_earth':consolidated,'money':money,'claims':records,'controls':controls,
            'evidence':'Counterparty stock accounting and consolidation only. No production, GDP, saving behaviour, cash-flow entitlement, market valuation or institutional mandate is inferred.'}
