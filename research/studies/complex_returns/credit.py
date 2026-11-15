"""Credit-supported investment with matched physical and low-price-pressure gates.

Production, customer budgets, capital costs and resource limits come from the
unchanged strategy producer. This new financing controller compares like-for-like
policies; its retained-earnings control is also compared with the historical path.
All cash/debt values use the source's constant USD2041 valuation, not a U/USD rate.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path

from . import strategies as s

HERE = s.HERE


def annuity(years, rate):
    return years if rate == 0 else (1-(1+rate)**(-years))/rate


def external_capital_coefficients(b, t, profile):
    c = t['source_structure_localization']['space_chemical']['cost']
    contingency = c['capital_with_contingency_USD2041']/c['capital_direct_USD2041']
    established = contingency*(sum(c['groups_USD2041'].values())+c['seed_USD2041']+c['additional_route_seed_USD2041'])
    stocks = contingency*(c['initial_stock_manufacture_USD2041']+c['initial_stock_Earth_freight_USD2041'])
    common = b['common_capital_fraction_of_establishment']*c['capital_with_contingency_USD2041']
    commercial = profile['commissioning']['additional_external_purchases_USD2041']
    g = profile['remaining_balanced_growth']
    return (max(0., established-common)+stocks+commercial,
            g['at_limit']['external']['USD2041']['expansion_external_purchases']/g['fraction_of_residents_year']+commercial)


def external_order_floor(b, coeff, committed, increment):
    first = min(increment,max(0.,1-committed))
    return b['marginal_capital_multiplier']*(first*coeff[0]+(increment-first)*coeff[1])


def capital_split(b, cap, coeff, committed, increment, kit_cost, monetary):
    cost = s.staged.expansion_cost(b,cap,committed,increment)
    external_floor = external_order_floor(b,coeff,committed,increment)
    if external_floor > cost+1:
        raise ValueError('Identified external capital exceeds modeled cash bill')
    eligible = monetary['unclassified_capital_local_fraction']*max(0.,cost-external_floor)+monetary['service_capital_local_fraction']*kit_cost
    outside_U=monetary.get('outside_supplier_U_fraction',0.)*external_floor
    return dict(cost=cost+kit_cost,identified_external_floor=external_floor,eligible_outside_U=outside_U,
                eligible_local=eligible,other_externally_funded=cost+kit_cost-eligible-outside_U)


def debt_due(loans, year):
    interest = principal = U_principal = external_interest = U_interest = supplier_principal = supplier_interest = 0.
    for loan in loans:
        if loan['remaining'] <= 0 or year <= loan['year']:
            continue
        charge = loan['remaining']*loan['rate']
        repay = min(loan['remaining'],loan['original']/loan['term']) if year > loan['year']+loan['grace'] else 0.
        interest += charge; principal += repay
        if loan['kind']=='U':
            U_principal += repay; U_interest += charge
            supplier_principal += repay*loan.get('supplier_share',0.)
            supplier_interest += charge*loan.get('supplier_share',0.)
        else: external_interest += charge
    return dict(interest=interest,principal=principal,U_principal=U_principal,
                external_interest=external_interest,U_interest=U_interest,supplier_principal=supplier_principal,supplier_interest=supplier_interest)


def repay(loans, year):
    due = debt_due(loans,year)
    for loan in loans:
        if year > loan['year']+loan['grace']:
            loan['remaining'] -= min(loan['remaining'],loan['original']/loan['term'])
    return due


def concession(draw, rate, term, grace, benchmark):
    balance = draw; pv = 0.
    for y in range(1,grace+term+1):
        principal = draw/term if y>grace else 0.
        pv += (rate*balance+principal)/(1+benchmark)**y
        balance -= principal
    return max(0.,draw-pv)


def price_screen(basket, local_spending, local_payout, gross, U_outstanding,
                 U_escrow, contractor_escrow, initial_money_need, age, monetary, outside_U=0.):
    planned = basket+local_spending+max(0.,local_payout)
    local_need = monetary['domestic_holding_days']/365*planned
    settlement = monetary['settlement_initial_share']+(monetary['settlement_final_share']-monetary['settlement_initial_share'])*min(1.,age/monetary['settlement_transition_years'])
    customer_need = settlement*gross*monetary['customer_holding_days']/365
    money_growth = max(0.,local_need+customer_need-initial_money_need)
    circulating_credit = max(0.,U_outstanding-U_escrow-contractor_escrow-outside_U)
    excess = max(0.,circulating_credit-money_growth)
    spending = planned+excess*365/monetary['domestic_holding_days']
    price = (spending/basket)**(1/(1+monetary['supply_elasticity']))
    return dict(price=price,planned_spending=planned,transaction_money_growth=money_growth,
                circulating_credit=circulating_credit,excess_credit_balances=excess,
                local_transaction_balances=local_need,customer_balances=customer_need,other_money_creation=max(0.,money_growth-circulating_credit),total_money_creation=max(money_growth,circulating_credit))


def stressed_inputs(q, scenario, year, config):
    if not scenario or not config['stress']['start_year'] <= year <= config['stress']['end_year']:
        return q
    v = deepcopy(q)
    for market in v['markets'].values():
        for key in ('industrial_budget_USD2041_year','mission_budget_USD2041_year',
                    'household_budget_USD2041_year','premium_demand_kg_year',
                    'institutional_services_budget_USD2041_year','funded_compatible_missions_year',
                    'compute_contracted_capacity_W'):
            market[key] *= config['stress']['buyer_budget_factor']
    return v


def investment_plan(q,d,b,cap,coeff,name,market,scale,committed,kits,total_kits,age,paid,licenses,maximum,monetary):
    desired = min(b['service_kit']['max_units'],math.ceil(max(
        b['scenarios'][market]['research_contracts']/b['service_kit']['research_contracts_year'],
        b['scenarios'][market]['professional_clients']/b['service_kit']['professional_clients_year'])*min(1.,b['service_demand_initial_fraction']*(1+b['service_demand_annual_growth'])**min(age,100))))
    commissioning = b['commissioning_years']*q['additional_commissioning_capacity_fraction']
    kit_admissible = total_kits<desired and s.resource_screen(q,d,b,name,scale,total_kits+1+commissioning)['feasible']
    kit_cost = b['service_kit']['capital_USD2041'] if kit_admissible else 0.
    reserved = kits+(total_kits-kits+int(kit_admissible))*(1+q['additional_commissioning_capacity_fraction'])
    resources = s.resource_screen(q,d,b,name,scale,reserved)
    slope = b['desired_growth_rate']*scale
    if committed>b['taper_fraction']*maximum:
        slope = b['desired_growth_rate']*b['taper_fraction']*maximum*math.sqrt(max(0.,(maximum-committed)/((1-b['taper_fraction'])*maximum)))
    limit = max(0.,min(slope,resources['growth_fraction_ceiling']*scale,maximum-committed))
    def op(stock,units):
        return s.operation(q,d,b,name,market,stock,stock,units,age,paid,licenses)
    now = op(committed,kits); plus = op(committed+1e-5,kits)
    joint = False; target = committed+limit; extra_kits = int(kit_admissible); carry = 0.
    if plus['available_USD2041']<=now['available_USD2041']:
        limit = 0.; target = committed
        if total_kits<desired and not kit_admissible:
            lo,hi=committed,maximum
            load=total_kits+1+commissioning
            if s.resource_screen(q,d,b,name,hi,load)['feasible']:
                for _ in range(35):
                    mid=(lo+hi)/2
                    if s.resource_screen(q,d,b,name,mid,load)['feasible']:hi=mid
                    else:lo=mid
                candidate=hi*(1+1e-8)
                if candidate<=maximum and op(candidate,total_kits+1)['available_USD2041']>op(committed,total_kits)['available_USD2041']:
                    joint=True; target=candidate; extra_kits=1
                    limit=max(0.,min(slope,resources['growth_fraction_ceiling']*scale,target-committed))
                    def burn(stock,units):
                        return s.operation(q,d,b,name,market,stock,0,units,0,synergies=False,service_demand=False)['total_cost_USD2041']
                    carry=max(0.,burn(target,total_kits+1)-burn(committed,total_kits))*(b['commissioning_years']+b['production_ramp_years']+b['liquidity_buffer_years'])
    package_kit = b['service_kit']['capital_USD2041'] if joint else kit_cost
    split = capital_split(b,cap,coeff,committed,max(0.,target-committed),package_kit,monetary)
    budget = split['cost']+carry
    margin = op(target,total_kits+extra_kits)['available_USD2041']-op(committed,total_kits)['available_USD2041'] if budget else 0.
    growth_delay = math.ceil(math.log(target/max(committed,1e-15))/math.log1p(min(b['desired_growth_rate'],resources['growth_fraction_ceiling']))) if joint and target>committed and resources['growth_fraction_ceiling']>0 else 0
    return dict(limit=limit,kit_cost=kit_cost,joint=joint,budget=budget,
        eligible_local=split['eligible_local'],eligible_U=split['eligible_local']+split['eligible_outside_U'],eligible_outside_U=split['eligible_outside_U'],carry=carry,target=target,margin=margin,
        appraisal_delay=growth_delay+b['commissioning_years']+b['production_ramp_years'],
        physical_growth_ceiling=resources['growth_fraction_ceiling'],desired_growth_ceiling=slope,
        desired_kits=desired,committed=committed)


def allocate(b,cap,coeff,plan,fx_cash,earned,fx_restricted,U_escrow,buffer,monetary,spending_cap=math.inf,distribution_lock=False,supplier_prepay=0.):
    spendable=max(0.,min(fx_cash-buffer,earned+fx_restricted))
    kit = starts = industrial = U_used = U_local_used = U_outside_used = FX_used = own_used = local = 0.
    U_kit_reserve=0.
    def pay(increment,kit_cost):
        nonlocal fx_cash,earned,fx_restricted,U_escrow,U_used,U_local_used,U_outside_used,FX_used,own_used,local,supplier_prepay
        split=capital_split(b,cap,coeff,plan['committed'],increment,kit_cost,monetary)
        available_U=max(0.,U_escrow-supplier_prepay-U_kit_reserve)
        u_local=min(available_U,split['eligible_local'])
        u_outside=min(supplier_prepay,split['eligible_outside_U'])
        u=u_local+u_outside;fx=split['cost']-u
        debt_fx=min(fx_restricted,fx); own=fx-debt_fx
        fx_cash-=fx;earned-=own;fx_restricted-=debt_fx;U_escrow-=u;supplier_prepay-=u_outside
        U_used+=u;U_local_used+=u_local;U_outside_used+=u_outside;FX_used+=debt_fx;own_used+=own;local+=split['eligible_local']
        return split['cost']
    if plan['kit_cost']:
        local_kit=monetary['service_capital_local_fraction']*plan['kit_cost']
        if spendable+min(max(0.,U_escrow-supplier_prepay),local_kit)>=plan['kit_cost'] and plan['kit_cost']<=spending_cap:
            kit=pay(0.,plan['kit_cost'])
    spendable=max(0.,min(fx_cash-buffer,earned+fx_restricted))
    if plan['kit_cost'] and not kit:
        U_kit_reserve=min(max(0.,U_escrow-supplier_prepay),monetary['service_capital_local_fraction']*plan['kit_cost'])
        spendable=max(0.,spendable-(plan['kit_cost']-U_kit_reserve))
    budget_left=max(0.,plan['budget']-kit)
    joint_ready=not plan['joint'] or spendable+min(U_escrow,plan.get('eligible_U',plan['eligible_local']))>=budget_left
    if joint_ready:
        lo,hi=0.,plan['limit']
        for _ in range(35):
            mid=(lo+hi)/2
            split=capital_split(b,cap,coeff,plan['committed'],mid,0.,monetary)
            if split['cost']<=min(spendable+min(max(0.,U_escrow-supplier_prepay-U_kit_reserve),split['eligible_local'])+min(supplier_prepay,split['eligible_outside_U']),spending_cap-kit):lo=mid
            else:hi=mid
        starts=lo if lo>1e-12 else 0.
        if starts:industrial=pay(starts,0.)
    reserve=max(0.,plan['budget']-kit-industrial) if plan['joint'] else (plan['kit_cost'] if not kit else 0.)
    own_reserve=max(0.,reserve-fx_restricted-U_escrow)
    excess=0. if distribution_lock else max(0.,min(fx_cash-buffer-fx_restricted,earned)-own_reserve)
    fx_cash-=excess;earned-=excess
    return dict(fx_cash=fx_cash,earned=earned,fx_restricted=fx_restricted,U_escrow=U_escrow,
                supplier_prepay=supplier_prepay,kit_order=kit,industrial_order=industrial,starts=starts,dividend=excess,
                local_capital=local,U_used=U_used,U_local_used=U_local_used,U_outside_used=U_outside_used,FX_used=FX_used,own_used=own_used,
                project_reserve=reserve,joint_ready=joint_ready)


def simulate(q,d,b,e,t,name,market,horizon,config,mode,control=None,stress=False):
    policy=config['financing'][mode];uw=config['underwriting'];mon=dict(config['monetary'],outside_supplier_U_fraction=policy.get('outside_supplier_U_fraction',0.))
    p=s.profile(q,d,name);cap=s.capital(q,d,b,e,t,name);coeff=external_capital_coefficients(b,t,p)
    seed=b['seed_scale'];start=b['bootstrap_years'];people=p['engineering']['population']
    maximum=e['growth']['0.01']['vacuum_architectural_reference_with_existing_cities_people']/people
    kits=1 if s.resource_screen(q,d,b,name,seed,1)['feasible'] else 0
    common=b['common_capital_fraction_of_establishment']*cap['establishment_USD2041']
    hardware=common+(cap['one_time_bill_USD2041']-common)*seed+cap['common_development_USD2041']+b['cultural_platform']['capital_USD2041']+kits*b['service_kit']['capital_USD2041']
    def burn(stock,units):
        return s.operation(q,d,b,name,market,stock,0,units,0,synergies=False,service_demand=False)['total_cost_USD2041']
    no_sales=burn(seed,kits);cash=no_sales*b['seed_no_receipt_runway_years'];carry=no_sales*b['bootstrap_carrying_fraction']*start
    external_seed=hardware+cash+carry;reserves=mon['initial_reserve_fraction_of_seed']*external_seed
    initial_reserve=reserves;earned=FX_restricted=U_escrow=paid=licenses=0.
    cohorts=[(start,seed)];pending=[];kit_pending=[];loans=[];rows=[];failure=None
    local_history=[];U_history=[];initial_money=None;previous_relative=1.;custody=0.;outside_U=0.;supplier_prepay=0.
    controls={r['year']:r for r in control['annual']} if control else {}
    for year in range(start,horizon+1):
        active_q=stressed_inputs(q,stress,year,config)
        cohorts.extend((due,n) for due,n in pending if due==year);pending=[x for x in pending if x[0]>year]
        kits+=sum(n for due,n in kit_pending if due==year);kit_pending=[x for x in kit_pending if x[0]>year]
        scale=sum(n for _,n in cohorts);productive=sum(n*min(1.,(year-due+1)/b['production_ramp_years']) for due,n in cohorts)
        op=s.operation(active_q,d,b,name,market,scale,productive,kits,year-start,paid,licenses)
        if not op['resource_screen']['feasible']:raise ValueError('Physical operating capacity exceeded')
        paid=op['paid_accounts'];licenses=op['licenses']
        pooled=min(0.,op['available_USD2041'])+max(0.,op['available_USD2041'])*b['collectible_fraction']
        opening_cash=cash;opening_U=U_escrow;opening_reserve=reserves
        cash+=pooled;earned+=pooled
        due=debt_due(loans,year);service=due['interest']+due['principal']
        dscr=max(0.,pooled)*uw['pledgeable_operating_cash_fraction']/service if service>0 else None
        # Restricted construction proceeds cannot pay existing debt service.
        if cash-FX_restricted<service:
            failure=dict(year=year,type='operating_or_debt_payment_shortfall',shortfall=service-(cash-FX_restricted),accrued_interest=due['interest'],U_accrued_interest=due['U_interest'],supplier_accrued_interest=due['supplier_interest'])
            rows.append(dict(year=year,scale=scale,productive_scale=productive,service_kits=kits,operation=op,
                failure=failure,debt=sum(x['remaining'] for x in loans),U_debt=sum(x['remaining'] for x in loans if x['kind']=='U'),
                cash=cash,U_escrow=U_escrow,reserves=reserves,custody=custody,outside_U=outside_U,supplier_prepay=supplier_prepay,DSCR=dscr))
            break
        supplier_principal=due['supplier_principal']
        supplier_interest=due['supplier_interest']
        repay(loans,year);cash-=service;earned-=service
        outside_U=sum(x['remaining']*x.get('supplier_share',0.) for x in loans if x['kind']=='U')
        covenant_lock=dscr is not None and dscr<uw['min_DSCR']
        buffer=burn(scale,kits)*max(b['liquidity_buffer_years'],b['seed_no_receipt_runway_years']-(year-start))
        U_debt=sum(x['remaining'] for x in loans if x['kind']=='U')
        reserve_target=mon['base_reserve_years_core_imports']*op['core_imports_USD2041']+(mon['U_credit_conversion_buffer_fraction']+mon['U_credit_loss_capital_fraction'])*U_debt
        contribution=min(max(0.,reserve_target-reserves),mon['reserve_contribution_fraction']*max(0.,pooled-service),max(0.,min(cash-buffer-FX_restricted,earned)))
        reserves+=contribution;cash-=contribution;earned-=contribution
        base_dividend=0. if covenant_lock else min(max(0.,min(cash-buffer-FX_restricted,earned)),max(0.,pooled-service)*(1-b['retained_fraction']))
        cash-=base_dividend;earned-=base_dividend
        committed=scale+sum(n for _,n in pending);total_kits=kits+sum(n for _,n in kit_pending)
        plan=investment_plan(active_q,d,b,cap,coeff,name,market,scale,committed,kits,total_kits,year-start,paid,licenses,maximum,mon)
        own_available=max(0.,min(cash-buffer,earned+FX_restricted))
        own_equity=max(0.,min(cash-buffer-FX_restricted,earned))
        gap=max(0.,plan['budget']-own_available-min(U_escrow,plan['eligible_U']))
        project_npv=plan['margin']*annuity(uw['project_appraisal_years'],uw['project_discount_rate'])/(1+uw['project_discount_rate'])**plan['appraisal_delay']-plan['budget']
        debt=sum(x['remaining'] for x in loans)
        annual_committed=sum(x['remaining']*x['rate']+x['original']/x['term'] for x in loans if x['remaining']>0)
        annual_headroom=max(0.,max(0.,pooled)*uw['pledgeable_operating_cash_fraction']/uw['min_DSCR']-annual_committed)
        max_credit=min(gap,uw['max_loan_fraction_of_budget']*plan['budget'],max(0.,uw['max_debt_fraction_of_seed']*external_seed-debt),annual_headroom/(policy['rate']+1/policy['amortization_years']))
        refusal=[]
        if gap>1 and policy['require_project_NPV'] and project_npv<=0:max_credit=0.;refusal.append('incremental_project_NPV')
        if annual_headroom<=0 and gap>1:refusal.append('debt_service_coverage')
        max_credit=min(max_credit,own_equity*uw['max_loan_fraction_of_budget']/(1-uw['max_loan_fraction_of_budget']))
        if mode=='retained_earnings' or covenant_lock:max_credit=0.
        if covenant_lock:plan=dict(plan,limit=0.,kit_cost=0.);refusal.append('coverage_covenant')
        reserve_free=max(0.,reserves-mon['base_reserve_years_core_imports']*op['core_imports_USD2041'])
        U_cap=max(0.,reserve_free/(mon['U_credit_conversion_buffer_fraction']+mon['U_credit_loss_capital_fraction'])-U_debt)
        U_cap=min(U_cap,max(0.,uw['max_loan_fraction_of_budget']*plan['eligible_U']-U_escrow)) if policy['U_credit'] else 0.
        if not policy['external_credit']:max_credit=min(max_credit,U_cap)
        basket=(p['domestic']['USD2041']['resident_consumption']+p['domestic']['USD2041']['domestic_capital_formation'])*scale
        if initial_money is None:
            initial_money=mon['domestic_holding_days']/365*basket+mon['settlement_initial_share']*op['gross_USD2041']*mon['customer_holding_days']/365
        period=b['commissioning_years']
        ceiling=min(1+mon['cumulative_relative_price_limit'],max(1.,previous_relative)*(1+mon['annual_relative_price_limit']))
        def proposal(draw,spending_cap=math.inf):
            u=min(draw,U_cap);fx=draw-u
            supplier_new=min(u,max(0.,uw['max_loan_fraction_of_budget']*plan['eligible_outside_U']-supplier_prepay))
            al=allocate(b,cap,coeff,plan,cash+fx,earned,FX_restricted+fx,U_escrow+u,buffer,mon,spending_cap,covenant_lock,supplier_prepay+supplier_new)
            local=(sum(local_history[-(period-1):])+al['local_capital'])/period if period>1 else al['local_capital']
            # The unspent fraction of construction prepayments remains restricted.
            prepaid=sum(value*(1-(year-y+1)/period) for y,value in U_history[-period:] if year-y+1<period)+al['U_local_used']*(1-1/period)
            planned_local=mon['incremental_distribution_spending_fraction']*(base_dividend+al['dividend'])
            available_local=custody+planned_local
            def priced(release):
                return price_screen(basket,local,release,op['gross_USD2041'],U_debt+u,al['U_escrow']-al['supplier_prepay'],prepaid,initial_money,year-start,mon,outside_U+supplier_new)
            pr=priced(available_local)
            reference=controls.get(year,{}).get('price',pr['price'])
            released=available_local
            if control and pr['price']>reference*ceiling:
                spending_ceiling=basket*(reference*ceiling)**(1+mon['supply_elasticity'])
                released=min(available_local,max(0.,spending_ceiling-basket-local))
                pr=priced(released)
            future_feasible=True
            previous_forecast_relative=pr['price']/reference
            # Test the unavoidable remainder of current construction prepayments.
            # Current sales are held fixed; no future demand growth underwrites it.
            for ahead in range(1,period):
                future_year=year+ahead
                if control and future_year in controls and 'price' in controls[future_year]:
                    future_local=sum((local_history+[al['local_capital']])[-(period-ahead):])/period
                    future_prepaid=sum(value*(1-(future_year-y+1)/period) for y,value in U_history[-period:]+[(year,al['U_local_used'])] if 0<future_year-y+1<period)
                    future_pr=price_screen(basket,future_local,0.,op['gross_USD2041'],U_debt+u,al['U_escrow']-al['supplier_prepay'],future_prepaid,initial_money,year-start+ahead,mon,outside_U+supplier_new)
                    future_relative=future_pr['price']/controls[future_year]['price']
                    future_limit=min(1+mon['cumulative_relative_price_limit'],max(1.,previous_forecast_relative)*(1+mon['annual_relative_price_limit']))
                    if future_relative>future_limit+1e-12:future_feasible=False
                    previous_forecast_relative=future_relative
            pr.update(local_payout=released,planned_local_payout=planned_local,custody=available_local-released,future_feasible=future_feasible,supplier_new=supplier_new)
            relative=pr['price']/reference
            return al,pr,relative,u,fx
        base=proposal(0.)
        candidate=proposal(max_credit) if max_credit>0 else base
        if control and (year not in controls or 'price' not in controls[year]):max_credit=0.;refusal.append('matched_control_unavailable')
        elif control and (candidate[2]>ceiling+1e-12 or not candidate[1]['future_feasible']) and max_credit>0:
            lo,hi=0.,max_credit
            for _ in range(35):
                mid=(lo+hi)/2
                trial=proposal(mid)
                if trial[2]<=ceiling+1e-12 and trial[1]['future_feasible']:lo=mid
                else:hi=mid
            max_credit=lo;refusal.append('price_pressure')
        candidate=proposal(max_credit)
        spending_cap=math.inf
        if control and (candidate[2]>ceiling+1e-12 or not candidate[1]['future_feasible']):
            lo,hi=0.,candidate[0]['industrial_order']+candidate[0]['kit_order']
            for _ in range(35):
                mid=(lo+hi)/2
                trial=proposal(max_credit,mid)
                if trial[2]<=ceiling+1e-12 and trial[1]['future_feasible']:lo=mid
                else:hi=mid
            spending_cap=lo;candidate=proposal(max_credit,spending_cap)
            refusal.append('investment_spending_schedule')
        if control and (year not in controls or 'price' not in controls[year]):spending_cap=0.;candidate=proposal(0.,0.)
        al,pr,relative,u,fx=candidate
        # Do not borrow into an unusable escrow when the whole package remains unfunded.
        if max_credit>0 and ((plan['joint'] and not al['joint_ready']) or al['industrial_order']+al['kit_order']<=base[0]['industrial_order']+base[0]['kit_order']+1):
            max_credit=u=fx=0.;al,pr,relative,_,_=proposal(0.,spending_cap);refusal.append('complete_package_or_use')
        for kind,draw in (('U',u),('external',fx)):
            if draw>1e-6:loans.append(dict(year=year,original=draw,remaining=draw,rate=policy['rate'],term=policy['amortization_years'],grace=uw['principal_grace_years'],kind=kind,supplier_share=pr['supplier_new']/draw if kind=='U' else 0.))
        cash=al['fx_cash'];earned=al['earned'];FX_restricted=al['fx_restricted'];U_escrow=al['U_escrow'];supplier_prepay=al['supplier_prepay']
        if al['kit_order']:kit_pending.append((year+b['commissioning_years'],1))
        if al['starts']:pending.append((year+b['commissioning_years'],al['starts']))
        local_history.append(al['local_capital']);U_history.append((year,al['U_local_used']));outside_U+=pr['supplier_new']
        dividend=base_dividend+al['dividend'];custody=pr['custody'];debt=sum(x['remaining'] for x in loans);U_debt=sum(x['remaining'] for x in loans if x['kind']=='U')
        relative_change=relative/previous_relative-1;premium_change=relative/max(1.,previous_relative)-1;previous_relative=relative
        cost=al['industrial_order']+al['kit_order']
        cash_residual=opening_cash+opening_U+pooled+u+fx-service-contribution-dividend-cost-cash-U_escrow
        if abs(cash_residual)>max(1.,1e-11*max(1.,abs(pooled),cost)):
            raise ValueError('Cash conservation failed')
        rows.append(dict(year=year,scale=scale,productive_scale=productive,service_kits=kits,operation=op,
            opening_cash=opening_cash,cash=cash,U_escrow=U_escrow,FX_restricted=FX_restricted,earned=earned,buffer=buffer,
            reserves=reserves,reserve_contribution=contribution,reserve_target=reserve_target,
            external_draw=fx,U_draw=u,debt=debt,U_debt=U_debt,interest=due['interest'],principal=due['principal'],
            U_principal=due['U_principal'],supplier_principal=supplier_principal,supplier_interest=supplier_interest,U_financed_external_capital=al['U_outside_used'],supplier_credit_origination=pr['supplier_new'],supplier_prepay=supplier_prepay,outside_U=outside_U,external_interest=due['external_interest'],U_interest=due['U_interest'],DSCR=dscr,
            dividend=dividend,custody=custody,local_payout=pr['local_payout'],cash_distribution_to_owners=dividend-pr['planned_local_payout']+pr['local_payout'],covenant_lock=covenant_lock,industrial_order=al['industrial_order'],service_order=al['kit_order'],starts=al['starts'],
            local_capital=al['local_capital'],local_construction_flow=pr['planned_spending']-basket-pr['local_payout'],
            project_budget=plan['budget'],project_NPV=project_npv,project_operating_increment=plan['margin'],
            project_reserve=al['project_reserve'],physical_growth_ceiling=plan['physical_growth_ceiling'],
            desired_growth_ceiling=plan['desired_growth_ceiling'],price=pr['price'],relative_price=relative,
            pressure_comparison_available=not control or (year in controls and 'price' in controls[year]),annual_relative_price_change=relative_change,annual_positive_premium_growth=premium_change,price_limit_breach=relative>1+mon['cumulative_relative_price_limit']+1e-9 or premium_change>mon['annual_relative_price_limit']+1e-9,
            transaction_money_growth=pr['transaction_money_growth'],circulating_credit=pr['circulating_credit'],excess_credit_balances=pr['excess_credit_balances'],
            concession=concession(u+fx,policy['rate'],policy['amortization_years'],uw['principal_grace_years'],uw['external_benchmark_rate']),
            cash_residual=cash_residual,refused_by=refusal))
    return dict(strategy=name,market=market,mode=mode,stress=stress,external_seed=external_seed,
                initial_reserve=initial_reserve,initial_hardware_and_carry=hardware+carry,
                initial_working_cash=no_sales*b['seed_no_receipt_runway_years'],annual=rows,failure=failure,loans=loans)


def summary(sim, horizon, config, control=None):
    rows=[r for r in sim['annual'] if r['year']<=horizon]
    completed=[r for r in rows if 'failure' not in r]
    last=rows[-1];failure=sim['failure'] if sim['failure'] and sim['failure']['year']<=horizon else None
    def total(key):return math.fsum(r.get(key,0.) for r in completed)
    # Match the inherited seed cash dates: construction spread over years0..7,
    # with working cash funded in year7. Reserves are funded in year0 separately.
    b,_,_=s.staged.inputs();start=b['bootstrap_years']
    def initial_pv(rate):
        return sum(-(sim['initial_hardware_and_carry']/start)/(1+rate)**y for y in range(start))-sim['initial_working_cash']/(1+rate)**(start-1)
    equity={str(rate): initial_pv(rate)+sum(r['cash_distribution_to_owners']/(1+rate)**r['year'] for r in completed) for rate in (.0,.01,.03,.07)}
    external_lender_npv=sum((r['external_interest']+r['principal']-r['U_principal']-r['external_draw'])/1.03**r['year'] for r in completed)
    dscr=[r['DSCR'] for r in rows if r.get('DSCR') is not None]
    debts=[r.get('debt',0.) for r in rows]
    relative=[r['relative_price'] for r in completed]
    claims=last.get('debt',0.)+(failure['accrued_interest'] if failure else 0.)
    U_claims=last.get('U_debt',0.)+(failure['U_accrued_interest'] if failure else 0.)
    recovery=min(claims,max(0.,last['cash']+last.get('U_escrow',0.)-last.get('supplier_prepay',0.))) if failure else 0.
    funding_loss=claims-recovery if failure else 0.
    U_loss=funding_loss*U_claims/claims if claims>0 else 0.
    if failure and claims>0:external_lender_npv+=recovery*(1-U_claims/claims)/1.03**failure['year']
    out=dict(strategy=sim['strategy'],market=sim['market'],mode=sim['mode'],stress=sim['stress'],horizon=horizon,
        status='default_or_operating_failure' if failure else 'continuing_conditional_path',failure=failure,
        original_seed=sim['external_seed'],initial_monetary_reserve=sim['initial_reserve'],
        final_scale=last['scale'],final_productive_scale=last['productive_scale'],final_service_kits=last['service_kits'],
        first_founding_scale_year=next((r['year'] for r in completed if r['scale']>=1-1e-8),None),
        first_second_service_kit_year=next((r['year'] for r in completed if r['service_kits']>=2),None),
        productive_scale_years=math.fsum(r['productive_scale'] for r in completed),
        service_kit_years=sum(r['service_kits'] for r in completed),
        cumulative_external_loan_draws=total('external_draw'),cumulative_U_credit_draws=total('U_draw'),
        cumulative_capital_orders=total('industrial_order')+total('service_order'),
        cumulative_interest=total('interest'),cumulative_principal=total('principal'),
        peak_debt=max(debts,default=0.),terminal_debt=last.get('debt',0.),
        terminal_cash_and_U_escrow=last['cash']+last.get('U_escrow',0.)-last.get('supplier_prepay',0.),terminal_supplier_capital_prepayments=last.get('supplier_prepay',0.),terminal_reserves=last['reserves'],terminal_supplier_U_claims=last['outside_U'],supplier_accrued_interest_on_default=failure['supplier_accrued_interest'] if failure else 0.,cumulative_supplier_U_finance=total('supplier_credit_origination'),cumulative_supplier_prepayments_used=total('U_financed_external_capital'),supplier_term_claim_financial_NPV_at_face_at_3_percent=sum((r['supplier_principal']+r['supplier_interest']-r['supplier_credit_origination'])/1.03**r['year'] for r in completed),reserve_capital_loss=min(U_loss,config['monetary']['U_credit_loss_capital_fraction']*last.get('U_debt',0.)),terminal_owner_stabilization_assets=last['custody'],
        minimum_DSCR=min(dscr) if dscr else None,
        max_relative_price_change_percent=100*(max(relative,default=1.)-1),
        max_annual_relative_price_change_percent=100*max((r['annual_relative_price_change'] for r in completed),default=0.),
        max_annual_positive_premium_growth_percent=100*max((r['annual_positive_premium_growth'] for r in completed),default=0.),
        price_limit_breach_years=[r['year'] for r in completed if r['price_limit_breach']],
        equity_distribution_NPV=equity,external_lender_cash_NPV_at_3_percent=external_lender_npv,U_bank_interest_income=total('U_interest'),U_bank_loan_face_before_default_resolution=last.get('U_debt',0.),U_created_deposit_liabilities_before_any_default_resolution=last.get('U_debt',0.),
        concession_PV_at_origination=total('concession'),
        lender_loss_under_cash_only_recovery=funding_loss,default_cash_recovery=recovery,U_credit_loss=U_loss,
        U_loss_above_allocated_capital=max(0.,U_loss-config['monetary']['U_credit_loss_capital_fraction']*last.get('U_debt',0.)),
        deposits_not_automatically_cancelled_on_default=completed[-1]['circulating_credit'] if failure and completed else 0.,
        U_bank_interest_after_supplier_interest=total('U_interest')-total('supplier_interest'),
        max_cash_identity_residual=max((abs(r['cash_residual']) for r in completed),default=0.),
        new_loan_refusal_counts={key:sum(key in r['refused_by'] for r in completed) for key in ('incremental_project_NPV','debt_service_coverage','price_pressure','complete_package_or_use','matched_control_unavailable','coverage_covenant','investment_spending_schedule')})
    for rate in (.0,.01,.03,.07):
        terminal_date=failure['year'] if failure else horizon
        terminal_equity=max(0.,last['cash']+last.get('U_escrow',0.)-last.get('supplier_prepay',0.)-claims) if failure else last['cash']+last.get('U_escrow',0.)-last.get('supplier_prepay',0.)-last.get('debt',0.)
        out.setdefault('equity_NPV_with_terminal_liquidation_claim',{})[str(rate)]=equity[str(rate)]+(terminal_equity+last['custody'])/(1+rate)**terminal_date
        out.setdefault('project_and_reserve_funders_cash_NPV',{})[str(rate)]=equity[str(rate)]-sim['initial_reserve']
        out.setdefault('project_and_reserve_funders_NPV_with_terminal_claims',{})[str(rate)]=out['equity_NPV_with_terminal_liquidation_claim'][str(rate)]-sim['initial_reserve']+(last['reserves']-min(U_loss,config['monetary']['U_credit_loss_capital_fraction']*last.get('U_debt',0.)))/(1+rate)**terminal_date
    if control:
        counterpart={r['year']:r for r in control['annual'] if r['year']<=horizon and 'failure' not in r}
        common=[r for r in completed if r['year'] in counterpart]
        valid=not failure and not (control['failure'] and control['failure']['year']<=horizon)
        out['matched_comparison_complete']=valid
        out['incremental_operating_less_capital_NPV_at_3_percent']=sum((r['operation']['available_USD2041']-counterpart[r['year']]['operation']['available_USD2041']-(r['industrial_order']+r['service_order']-counterpart[r['year']]['industrial_order']-counterpart[r['year']]['service_order']))/1.03**r['year'] for r in common) if valid else None
        out['incremental_productive_scale_years']=sum(r['productive_scale']-counterpart[r['year']]['productive_scale'] for r in common) if valid else None
        out['incremental_service_kit_years']=sum(r['service_kits']-counterpart[r['year']]['service_kits'] for r in common) if valid else None
    return out


def product():
    config=json.loads((HERE/'credit_inputs.json').read_text());q,d,b,e,t=s.inputs()
    saved_path=HERE/'results/strategies.json';saved=json.loads(saved_path.read_text())
    if s.run.digest(saved_path)!=config['strategy_result_sha256']:raise ValueError('Pinned strategy product changed')
    if saved['producer_sha256']!=s.run.digest(HERE/'strategies.py'):raise ValueError('Strategy producer changed')
    for path,sha in saved['dependencies'].items():
        if s.run.digest(HERE/path)!=sha:raise ValueError('Strategy dependency changed: '+path)
    expected={(r['strategy'],r['market'],r['horizon_elapsed_years']):r for r in saved['horizons']}
    results=[];checkpoints=[];stress_results=[];severe_results=[];sensitivity=[];historical=[]
    for name in q['strategies']:
        for market in q['markets']:
            original=s.simulate(q,d,b,e,t,name,market,max(config['horizons']))
            if original['summary']!=expected[name,market,max(config['horizons'])]:raise ValueError('Historical financial reproduction failed')
            control=simulate(q,d,b,e,t,name,market,max(config['horizons']),config,'retained_earnings')
            for mode in config['financing']:
                sim=control if mode=='retained_earnings' else simulate(q,d,b,e,t,name,market,max(config['horizons']),config,mode,control)
                results.extend(summary(sim,h,config,control if mode!='retained_earnings' else None) for h in config['horizons'])
                checkpoints.append(dict(strategy=name,market=market,mode=mode,rows=[{k:v for k,v in r.items() if k!='operation'}|dict(gross=r['operation']['gross_USD2041'],available=r['operation']['available_USD2041'],core_imports=r['operation']['core_imports_USD2041']) for r in sim['annual'] if r['year'] in config['checkpoint_years'] or r is sim['annual'][-1]]))
            for h in config['horizons']:
                ref=expected[name,market,h]
                matched=summary(control,h,config)
                historical.append(dict(strategy=name,market=market,horizon=h,historical_final_scale=ref['final_scale'],matched_final_scale=matched['final_scale'],historical_NPV_3_percent=ref['NPV_USD2041']['0.03'],matched_equity_distribution_NPV_3_percent=matched['equity_distribution_NPV']['0.03']))
            stress_control=simulate(q,d,b,e,t,name,market,100,config,'retained_earnings',stress=True)
            for mode in ('retained_earnings','funded_loans','mixed_credit','development_mix','international_U_mix'):
                sim=stress_control if mode=='retained_earnings' else simulate(q,d,b,e,t,name,market,100,config,mode,stress_control,stress=True)
                stress_results.append(summary(sim,100,config,stress_control if mode!='retained_earnings' else None))
        severe=deepcopy(config)
        severe['stress'].update(start_year=20,end_year=29,buyer_budget_factor=.25)
        severe_control=simulate(q,d,b,e,t,name,'breakthrough',100,severe,'retained_earnings',stress=True)
        for mode in ('retained_earnings','funded_loans','mixed_credit','development_mix','international_U_mix'):
            sim=severe_control if mode=='retained_earnings' else simulate(q,d,b,e,t,name,'breakthrough',100,severe,mode,severe_control,stress=True)
            severe_results.append(summary(sim,100,severe,severe_control if mode!='retained_earnings' else None))
        # Keep sensitivities bounded to the middle and strong markets at100years.
        for variant in ('tight_prices','looser_prices','no_local_eligibility','all_residual_local','higher_interest','higher_DSCR'):
            v=deepcopy(config)
            if variant=='tight_prices':v['monetary']['annual_relative_price_limit']=.0025
            if variant=='looser_prices':v['monetary']['annual_relative_price_limit']=.01
            if variant=='no_local_eligibility':v['monetary']['unclassified_capital_local_fraction']=0.;v['monetary']['service_capital_local_fraction']=0.
            if variant=='all_residual_local':v['monetary']['unclassified_capital_local_fraction']=1.;v['monetary']['service_capital_local_fraction']=.5
            if variant=='higher_interest':v['financing']['mixed_credit']['rate']=.06
            if variant=='higher_DSCR':v['underwriting']['min_DSCR']=2.
            for market in ('global','breakthrough'):
                control=simulate(q,d,b,e,t,name,market,100,v,'retained_earnings')
                sim=simulate(q,d,b,e,t,name,market,100,v,'mixed_credit',control)
                sensitivity.append(dict(variant=variant,**summary(sim,100,v,control)))
    paths=['credit_inputs.json','strategies.py','results/strategies.json',*saved['dependencies']]
    return dict(schema='terluna.research.credit-comparison/1',evidence=config['evidence'],boundary=config['boundary'],
        producer_sha256=s.run.digest(HERE/'credit.py'),dependencies={p:s.run.digest(HERE/p) for p in paths},
        reading_rule='Values are constant USD2041-equivalent source accounts, not a nominal U conversion. New matched policies share physical and customer constraints. Rate and price gates are scenarios; relative price pressure is not CPI inflation. Principal is not an operating expense; draws are not income; terminal debt and reserve assets remain separate. Historic financial products are unchanged.',
        horizons=results,checkpoints=checkpoints,stress=stress_results,severe_stress=severe_results,sensitivities=sensitivity,historical_comparison=historical)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=s.run.ROOT/'research/runs/complex_returns/credit.json')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(product(),indent=2,allow_nan=False)+'\n');print(args.output)
