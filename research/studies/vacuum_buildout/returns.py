"""Dated sponsor and socioeconomic accounts for the prescribed year-200 build.

The two financing arms share a physical schedule. Residual external equity is a
reported requirement; it is never an inferred contract or free project receipt.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from . import run as m
from ..complex_returns import strategies as s, credit as c, consolidated as a

HERE=m.HERE


def initial_provisions(b,t,eng):
    cost=t['source_structure_localization']['space_chemical']['cost']
    contingency=cost['capital_with_contingency_USD2041']/cost['capital_direct_USD2041']
    original=cost['preservation_USD2041']*contingency
    stocks=cost['initial_stock_restoration_USD2041']*contingency
    added=eng['added_initial_cost']
    additional=added['installed_with_preservation_and_contingency_USD2041']-added['inherited_installed_quote_USD2041']*contingency
    return original,stocks,additional


def inputs():
    cfg=json.loads((HERE/'returns_inputs.json').read_text());q,eng,inv=m.inputs()
    for path,sha in cfg['pinned_dependencies'].items():
        if m.digest(HERE/path)!=sha:raise ValueError('Changed returns dependency: '+path)
    sq,d,b,e,t=s.inputs()
    credit=json.loads((s.HERE/'credit_inputs.json').read_text())
    old,_,_=s.run.inputs()
    return cfg,q,eng,sq,d,b,e,t,credit,old


def screen(q,eng,sq,d,b,name,scale,g,kits):
    r=m.response(eng,name,g);ratios=dict(r['capacity_use']);p=s.profile(sq,d,name)
    premium=sq['strategies'][name]['premium_share'];fuel=bool(p['repricing']['propellant_kg_year'])
    power=b['cultural_platform']['power_W']+kits*b['service_kit']['power_W']+premium*p['engineering']['finished_export_kg_year']*scale*sq['premium']['additional_mean_W_per_attempted_kg_year']
    cap=eng['architecture']['expanded'];source=d['physical_capacity']
    ratios['ground_firm_power']+=power*source['capacity_margin']/source['lunar_process_availability']/(scale*cap['engineering']['total_ground_firm_W'])
    ratios['lunar_process_power']+=power/(scale*cap['peak_flows']['local_mean_W'])
    ratios['depot_power']+=fuel*sq['logistics']['additional_depot_use_fraction']
    jobs=(r['physical']['technical_resident_jobs']+q['additional_technical_staff_per_founding_scale']+premium*sq['premium']['staff_per_founding_scale']+fuel*sq['logistics']['support_staff_per_founding_scale'])*scale+b['cultural_platform']['staff']+kits*b['service_kit']['staff']
    ratios['staff']=jobs/(q['founding_population']*scale*q['technical_workforce_fraction'])
    return dict(feasible=max(ratios.values())<=1+1e-10,maximum=max(ratios.values()),ratios=ratios)


def operating(q,eng,sq,d,b,name,market,scale,productive,kits,age,paid=0,licenses=0,**kwargs):
    op=s.operation(sq,d,b,name,market,scale,productive,kits,age,paid,licenses,**kwargs)
    extra=eng['incremental_modular_cost'][name]
    imports=extra['added_steady_external_purchases_USD2041_year']*scale
    provision=extra['added_steady_resource_provision_USD2041_year']*scale
    wages=q['additional_technical_staff_per_founding_scale']*q['additional_staff_USD2041_year']*scale
    op['core_imports_USD2041']+=imports
    op['resource_provision_USD2041']+=provision
    op['total_cost_USD2041']+=imports+provision+wages
    op['available_USD2041']-=imports+provision+wages
    op['added_technical_payroll']=wages
    op.pop('resource_screen')
    return op


def plan(data,name,market):
    cfg,q,eng,sq,d,b,e,t,cc,old=data
    p=s.profile(sq,d,name);cap=s.capital(sq,d,b,e,t,name)
    coeff=c.external_capital_coefficients(b,t,p);extra=eng['incremental_modular_cost'][name]
    added=sum(extra[k] for k in ('added_modular_external_purchases_USD2041_per_scale','added_modular_resource_provision_USD2041_per_scale','added_modular_local_component_value_USD2041_per_scale'))
    seed=b['seed_scale'];start=b['bootstrap_years'];lead=b['commissioning_years']
    kits=1;pending=[];cohorts=[(start,seed)];paid=licenses=0.;rows=[]
    def op(scale,productive,units,age,**kwargs):
        return operating(q,eng,sq,d,b,name,market,scale,productive,units,age,**kwargs)
    common=b['common_capital_fraction_of_establishment']*cap['establishment_USD2041']
    hardware=common+(cap['one_time_bill_USD2041']-common)*seed+cap['common_development_USD2041']+b['cultural_platform']['capital_USD2041']+b['service_kit']['capital_USD2041']+seed*eng['added_initial_cost']['installed_with_preservation_and_contingency_USD2041']
    burn=op(seed,0,1,0,synergies=False,service_demand=False)['total_cost_USD2041']
    carry=burn*b['bootstrap_carrying_fraction']*start;working=burn*b['seed_no_receipt_runway_years']
    for year in range(start,q['completion_year']+2):
        scale=m.population(year,q)/q['founding_population'];previous=m.population(year-1,q)/q['founding_population']
        if year>start and scale>previous:cohorts.append((year,scale-previous))
        productive=sum(n*min(1.,(year-due+1)/b['production_ramp_years']) for due,n in cohorts)
        kits+=sum(n for due,n in pending if due==year);pending=[(due,n) for due,n in pending if due>year]
        starts=(m.population(year+lead,q)-m.population(year+lead-1,q))/q['founding_population']
        g=starts/scale;age=year-start
        demand=min(b['service_kit']['max_units'],math.ceil(max(b['scenarios'][market]['research_contracts']/b['service_kit']['research_contracts_year'],b['scenarios'][market]['professional_clients']/b['service_kit']['professional_clients_year'])*min(1.,b['service_demand_initial_fraction']*(1+b['service_demand_annual_growth'])**min(age,100))))
        total=kits+sum(n for _,n in pending)
        kit_order=int(year<=q['completion_year']-lead and total<demand and screen(q,eng,sq,d,b,name,scale,g,total+1+lead*sq['additional_commissioning_capacity_fraction'])['feasible'])
        if kit_order:pending.append((year+lead,1))
        resource=screen(q,eng,sq,d,b,name,scale,g,kits+(total-kits+kit_order)*(1+sq['additional_commissioning_capacity_fraction']))
        if not resource['feasible']:raise ValueError('Joint capacity exceeded '+str((name,market,year,resource)))
        operation=op(scale,productive,kits,age,paid=paid,licenses=licenses)
        paid=operation['paid_accounts'];licenses=operation['licenses']
        committed=m.population(year+lead-1,q)/q['founding_population']
        kitcost=kit_order*b['service_kit']['capital_USD2041']
        split=c.capital_split(b,cap,coeff,committed,starts,kitcost,cc['monetary'])
        cost=split['cost']+starts*added
        external=split['identified_external_floor']+starts*extra['added_modular_external_purchases_USD2041_per_scale']
        local=split['eligible_local']+starts*extra['added_modular_local_component_value_USD2041_per_scale']*cc['monetary']['unclassified_capital_local_fraction']
        baseline=op(committed,committed,total,age,paid=paid,licenses=licenses)
        enlarged=op(committed+starts,committed+starts,total+kit_order,age,paid=paid,licenses=licenses)
        margin=enlarged['available_USD2041']-baseline['available_USD2041']
        uw=cc['underwriting']
        appraisal=-cost+margin*c.annuity(uw['project_appraisal_years'],uw['project_discount_rate'])/(1+uw['project_discount_rate'])**(lead+b['production_ramp_years'])
        original_provision,stock_provision,_=initial_provisions(b,t,eng)
        first=min(starts,max(0.,1-committed))
        reserve=a.capital_reserve(p,b,committed,starts)+starts*extra['added_modular_resource_provision_USD2041_per_scale']+first*((1-b['common_capital_fraction_of_establishment'])*original_provision+stock_provision)
        rows.append(dict(year=year,population=scale*q['founding_population'],scale=scale,productive_scale=productive,service_kits=kits,starts=starts,operation=operation,capital=cost,industrial_capital=cost-kitcost,service_capital=kitcost,external_capital_floor=external,local_capital=local,capital_resource_provision=reserve,appraisal=appraisal,resource=resource))
    return dict(strategy=name,market=market,seed_hardware=hardware,seed_carry=carry,initial_working=working,seed=hardware+carry+working,annual=rows)


def finance(path,data,mode,control=None):
    cfg,q,eng,sq,d,b,e,t,cc,old=data;mon=cc['monetary'];uw=cc['underwriting']
    initial_reserve=path['seed']*mon['initial_reserve_fraction_of_seed'];reserves=initial_reserve
    cash=path['initial_working'];loans=[];rows=[];initial_money=None;previous_relative=1.;local_history=[];custody=0.
    controls={} if control is None else {r['year']:r for r in control['annual']}
    for index,row in enumerate(path['annual']):
        year=row['year'];op=row['operation'];opening=cash
        due=c.repay(loans,year);service=due['interest']+due['principal'];available=op['available_USD2041']
        debt=sum(x['remaining'] for x in loans);Udebt=sum(x['remaining'] for x in loans if x['kind']=='U')
        outside=sum(x['remaining']*x.get('supplier_share',0.) for x in loans if x['kind']=='U')
        reserve_target=mon['base_reserve_years_core_imports']*op['core_imports_USD2041']+(mon['U_credit_conversion_buffer_fraction']+mon['U_credit_loss_capital_fraction'])*Udebt
        contribution=min(max(0.,reserve_target-reserves),mon['reserve_contribution_fraction']*max(0.,available-service))
        reserves+=contribution;cash+=available-service-contribution
        runway=max(b['liquidity_buffer_years'],b['seed_no_receipt_runway_years']-(year-b['bootstrap_years']))
        buffer=op['total_cost_USD2041']*runway
        gap=max(0.,row['capital']-max(0.,cash-buffer))
        annual_committed=sum(x['remaining']*x['rate']+x['original']/x['term'] for x in loans if x['remaining']>0)
        annual_headroom=max(0.,max(0.,available)*uw['pledgeable_operating_cash_fraction']/uw['min_DSCR']-annual_committed)
        covenant_DSCR=available*uw['pledgeable_operating_cash_fraction']/service if service else None
        covenant_lock=covenant_DSCR is not None and covenant_DSCR<uw['min_DSCR']
        draw=min(gap,uw['max_loan_fraction_of_budget']*row['capital'],max(0.,uw['max_debt_fraction_of_seed']*path['seed']-debt),annual_headroom/(cfg['external_and_supplier_rate']+1/20)) if mode=='latest_U_completion' and row['appraisal']>0 and not covenant_lock else 0.
        reserve_free=max(0.,reserves-mon['base_reserve_years_core_imports']*op['core_imports_USD2041'])
        ucap=max(0.,reserve_free/(mon['U_credit_conversion_buffer_fraction']+mon['U_credit_loss_capital_fraction'])-Udebt)
        supplier=min(draw,ucap,uw['max_loan_fraction_of_budget']*cfg['outside_supplier_U_fraction']*row['external_capital_floor'])
        local=min(draw-supplier,max(0.,ucap-supplier),uw['max_loan_fraction_of_budget']*row['local_capital'])
        external_draw=draw-supplier-local
        p=s.profile(sq,d,path['strategy']);basket=(p['domestic']['USD2041']['resident_consumption']+p['domestic']['USD2041']['domestic_capital_formation'])*row['scale']
        if initial_money is None:initial_money=mon['domestic_holding_days']/365*basket+mon['settlement_initial_share']*op['gross_USD2041']*mon['customer_holding_days']/365
        local_history.append(row['local_capital']);local_flow=sum(local_history[-b['commissioning_years']:])/b['commissioning_years']
        def price(local_draw,payout=0.):
            return c.price_screen(basket,local_flow,payout,op['gross_USD2041'],Udebt+supplier+local_draw,0.,0.,initial_money,year-b['bootstrap_years'],mon,outside_U=outside+supplier)
        pr=price(local);reference=controls.get(year,{}).get('price',pr['price'])
        ceiling=min(1+mon['cumulative_relative_price_limit'],max(1.,previous_relative)*(1+mon['annual_relative_price_limit']))
        if pr['price']/reference>ceiling:
            lo,hi=0.,local
            for _ in range(45):
                mid=(lo+hi)/2
                if price(mid)['price']/reference<=ceiling:lo=mid
                else:hi=mid
            local=lo;pr=price(local)
        # Existing circulating credit can still create a breach as needs change.
        relative=pr['price']/reference;breach=relative>ceiling+1e-10
        draw=external_draw+supplier+local
        for value,kind,rate,share in ((local,'U',cfg['local_U_rate'],0.),(supplier,'U',cfg['external_and_supplier_rate'],1.),(external_draw,'external',cfg['external_and_supplier_rate'],0.)):
            if value>1e-6:loans.append(dict(year=year,original=value,remaining=value,rate=rate,term=20,grace=uw['principal_grace_years'],kind=kind,supplier_share=share))
        pre_topup=cash+draw-row['capital'];topup=max(0.,buffer-pre_topup)
        # Sequential purpose allocation: current operations/debt, then capital, then liquidity.
        free_opening=max(0.,opening-buffer)
        operating_gap=max(0.,-free_opening-available)
        debt_gap=max(0.,service-max(0.,free_opening+available))
        capital_gap=max(0.,row['capital']-max(0.,free_opening+available-service-contribution+draw))
        topup_purposes=dict(operating=min(topup,operating_gap),debt=min(max(0.,topup-operating_gap),debt_gap),capital=min(max(0.,topup-operating_gap-debt_gap),capital_gap))
        topup_purposes['liquidity']=topup-sum(topup_purposes.values())
        cash=pre_topup+topup
        nxt=path['annual'][index+1] if index+1<len(path['annual']) else row
        next_due=c.debt_due(loans,year+1)
        next_service=next_due['interest']+next_due['principal']
        next_U=sum(x['remaining'] for x in loans if x['kind']=='U')-next_due['U_principal']
        next_target=mon['base_reserve_years_core_imports']*nxt['operation']['core_imports_USD2041']+(mon['U_credit_conversion_buffer_fraction']+mon['U_credit_loss_capital_fraction'])*next_U
        next_contribution=min(max(0.,next_target-reserves),mon['reserve_contribution_fraction']*max(0.,nxt['operation']['available_USD2041']-next_service))
        next_runway=max(b['liquidity_buffer_years'],b['seed_no_receipt_runway_years']-(year+1-b['bootstrap_years']))
        next_buffer=nxt['operation']['total_cost_USD2041']*next_runway
        retained_requirement=max(buffer,next_buffer+nxt['capital']+next_service+next_contribution-nxt['operation']['available_USD2041'])
        dividend=0. if covenant_lock else max(0.,cash-retained_requirement);cash-=dividend
        intended_local=dividend*mon['incremental_distribution_spending_fraction']
        releasable=custody+intended_local;released=releasable
        pr=price(local,released)
        if control is None:reference=pr['price']
        elif pr['price']/reference>ceiling:
            lo,hi=0.,released
            for _ in range(45):
                mid=(lo+hi)/2
                if price(local,mid)['price']/reference<=ceiling:lo=mid
                else:hi=mid
            released=lo;pr=price(local,released)
        custody=releasable-released
        owner_payout=dividend-intended_local+released
        relative=pr['price']/reference;breach=relative>ceiling+1e-10;previous_relative=relative
        residual=opening+available+draw+topup-service-contribution-row['capital']-dividend-cash
        if abs(residual)>1.:raise ValueError('Cash identity failed')
        rows.append(dict(year=year,contribution=topup,topup_by_purpose=topup_purposes,pre_topup_cash=pre_topup,dividend=dividend,owner_payout=owner_payout,owner_custody=custody,local_payout=released,cash=cash,reserves=reserves,reserve_contribution=contribution,capital=row['capital'],operating_available=available,draw=draw,local_U_draw=local,supplier_U_draw=supplier,external_draw=external_draw,debt=sum(x['remaining'] for x in loans),U_debt=sum(x['remaining'] for x in loans if x['kind']=='U'),interest=due['interest'],principal=due['principal'],local_U_interest=due['U_interest']-due['supplier_interest'],supplier_interest=due['supplier_interest'],external_interest=due['external_interest'],DSCR=covenant_DSCR,gross_DSCR=available/service if service else None,covenant_lock=covenant_lock,price=pr['price'],relative_price=relative,price_breach=breach,excess_credit_balances=pr['excess_credit_balances'],cash_residual=residual))
    return dict(mode=mode,initial_reserve=initial_reserve,annual=rows)


def metrics(path,sim,data,horizon):
    cfg,q,eng,sq,d,b,e,t,cc,old=data;n=horizon+1;start=b['bootstrap_years'];p=s.profile(sq,d,path['strategy'])
    sponsor=[0.]*n;unlevered=[0.]*n;exports=[0.]*n;resident=[0.]*n;economic_cost=[0.]*n
    for year in range(start):
        sponsor[year]-=(path['seed_hardware']+path['seed_carry'])/start
        unlevered[year]=sponsor[year];economic_cost[year]=-sponsor[year]
    sponsor[start-1]-=path['initial_working'];sponsor[0]-=sim['initial_reserve']
    # Identified seed and carrying resource provisions are transfers, not consumed inputs.
    extra=eng['incremental_modular_cost'][path['strategy']]
    seed_resource=b['seed_scale']*p['commissioning']['additional_resource_provision_USD2041']
    original_provision,stock_provision,added_provision=initial_provisions(b,t,eng)
    fraction=b['common_capital_fraction_of_establishment']+b['seed_scale']*(1-b['common_capital_fraction_of_establishment'])
    seed_resource+=fraction*original_provision+b['seed_scale']*(stock_provision+added_provision)
    carry_resource=b['seed_scale']*(p['external']['USD2041']['resource_provision']+extra['added_steady_resource_provision_USD2041_year'])*b['bootstrap_carrying_fraction']*start
    for year in range(start):economic_cost[year]-=(seed_resource+carry_resource)/start
    for row,fin in zip(path['annual'],sim['annual']):
        year=row['year']
        if year>horizon:break
        op=row['operation'];sponsor[year]=fin['owner_payout']-fin['contribution']
        unlevered[year]=op['available_USD2041']-row['capital']
        exports[year]=op['gross_USD2041'];resident[year]=p['domestic']['USD2041']['resident_consumption']*row['scale']
        payroll=a.embedded_resident_payroll(sq,b,p,path['strategy'],row)+op['added_technical_payroll']
        opportunity=row['population']*old['resident_work_fraction']*old['alternative_USD2041_worker_year']
        economic_cost[year]=op['total_cost_USD2041']+row['capital']-op['resource_provision_USD2041']-row['capital_resource_provision']-payroll+opportunity
    rows=[r for r in sim['annual'] if r['year']<=horizon];last=rows[-1];end=next(r for r in path['annual'] if r['year']==horizon)
    terminal=last['cash']+last['reserves']+last['owner_custody']-last['debt'];discounted={}
    for rate in cfg['rates']:
        net=a.pv(sponsor,rate);raw=a.pv(unlevered,rate);E=a.pv(exports,rate);R=a.pv(resident,rate);C=a.pv(economic_cost,rate)
        factor=sum((1+rate)**(-year) for year in range(start,horizon+1))
        discounted[str(rate)]=dict(sponsor_cash_NPV=net,unlevered_operating_less_capital_NPV=raw,terminal_net_book_claim_PV=terminal/(1+rate)**horizon,sponsor_with_terminal_book_claim_NPV=net+terminal/(1+rate)**horizon,full_price_socioeconomic_NPV=E+R-C,half_price_socioeconomic_NPV=.5*(E+R)-C,full_price_benefit_cost_ratio=(E+R)/C,export_value_break_even_with_full_resident_value=(C-R)/E,additional_annual_net_social_benefit=max(0.,C-E-R)/factor,additional_annual_sponsor_cash=max(0.,-net)/factor,matched_local_interest_relief_PV=sum(r['local_U_interest']*(cfg['external_and_supplier_rate']-cfg['local_U_rate'])/cfg['local_U_rate']/(1+rate)**r['year'] for r in rows),discounted_exports=E,discounted_resident_value=R,discounted_economic_cost=C,durable_cash_recovery_year=a.durable_recovery(sponsor,rate))
    return dict(strategy=path['strategy'],market=path['market'],mode=sim['mode'],horizon=horizon,population=end['population'],service_kits=end['service_kits'],seed=path['seed'],initial_monetary_reserve=sim['initial_reserve'],additional_external_equity=sum(r['contribution'] for r in rows),external_equity_by_purpose={key:sum(r['topup_by_purpose'][key] for r in rows) for key in ('operating','debt','capital','liquidity')},peak_annual_external_equity=max(r['contribution'] for r in rows),first_negative_operating_year=next((r['year'] for r in rows if r['operating_available']<0),None),first_external_topup_year=next((r['year'] for r in rows if r['contribution']>1),None),cumulative_capital_orders=sum(r['capital'] for r in rows),cumulative_dividends=sum(r['owner_payout'] for r in rows),peak_debt=max(r['debt'] for r in rows),cumulative_local_U_credit=sum(r['local_U_draw'] for r in rows),cumulative_supplier_U_credit=sum(r['supplier_U_draw'] for r in rows),minimum_DSCR=min((r['DSCR'] for r in rows if r['DSCR'] is not None),default=None),price_breach_years=[r['year'] for r in rows if r['price_breach']],maximum_relative_price=max(r['relative_price'] for r in rows),terminal_cash=last['cash'],terminal_owner_custody=last['owner_custody'],terminal_reserves=last['reserves'],terminal_debt=last['debt'],annual_extra_service_gross_at_50_percent_margin=max(0.,-end['operation']['available_USD2041'])/.5,annual_gross=end['operation']['gross_USD2041'],annual_operating_available=end['operation']['available_USD2041'],revenue_by_sector=end['operation']['revenue_by_sector_USD2041'],discounted=discounted,nonnegative_single_root_IRR=a.nonnegative_irr(sponsor))


def product():
    data=inputs();cfg,q,eng,sq,d,b,e,t,cc,old=data;results=[];paths=[];finance_paths=[]
    for name in sq['strategies']:
        for market in sq['markets']:
            path=plan(data,name,market);control=finance(path,data,'equity_completion');latest=finance(path,data,'latest_U_completion',control)
            paths.append(path)
            for sim in (control,latest):
                finance_paths.append(dict(strategy=name,market=market,**sim))
                for horizon in cfg['horizons']:results.append(metrics(path,sim,data,horizon))
    return dict(schema='terluna.research.vacuum-buildout-returns/1',evidence=cfg['evidence'],producer_sha256=m.digest(__file__),dependencies={k:m.digest(HERE/k) for k in ['returns_inputs.json','run.py','buildout_inputs.json','inputs/engineering.json','results/buildout.json']},rules=cfg['rules'],horizons=results,physical_and_operating_paths=paths,financing_paths=finance_paths)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE.parents[1]/'runs/vacuum_buildout/returns.json');args=parser.parse_args()
    result=product();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(str(args.output));print('Return records:',len(result['horizons']))

if __name__=='__main__':main()
