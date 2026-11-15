"""Joint export portfolios, finite customer budgets and earned-cash expansion.

The financial consumer reads frozen provisioning responses. It imports no
provisioning implementation and makes no claim to qualify new material recipes.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from . import staged
from . import run

HERE = run.HERE


def inputs():
    q=json.loads((HERE/'strategy_inputs.json').read_text())
    p=HERE/q['profile_snapshot']
    if run.digest(p)!=q['profile_sha256']:raise ValueError('Unreviewed strategy responses')
    d=json.loads(p.read_text())
    if d['schema']!='terluna.research.strategy-profile-export/1':raise ValueError('Unexpected response schema')
    b,e,t=staged.inputs()
    return q,d,b,e,t


def profile(q,d,name):return d['profiles'][q['strategies'][name]['profile']]


def capital(q,d,b,e,t,name):
    p=profile(q,d,name);s=q['strategies'][name];a=staged.capital_account(e,t)
    commission=sum(p['commissioning'][k] for k in ('additional_external_purchases_USD2041',
        'local_civil_installation_allowance_USD2041','local_new_component_value_USD2041','additional_resource_provision_USD2041'))
    premium=s['premium_share'];logistics=bool(p['repricing']['propellant_kg_year'])
    a['one_time_bill_USD2041']+=commission-a['commercial_USD2041']+premium*q['premium']['installed_capital_USD2041_per_founding_scale']+logistics*q['logistics']['support_capital_USD2041_per_founding_scale']
    a['commercial_USD2041']=commission
    g=p['remaining_balanced_growth'];v=g['at_limit']['external']['USD2041'];base=p['external']['USD2041']
    a['marginal_post_founding_USD2041_per_scale']=(v['expansion_external_purchases']+v['resource_provision']-base['resource_provision'])/g['fraction_of_residents_year']+commission+premium*q['premium']['installed_capital_USD2041_per_founding_scale']+logistics*q['logistics']['support_capital_USD2041_per_founding_scale']
    a['common_development_USD2041']=premium*q['premium']['common_development_USD2041']
    return a


def resource_screen(q,d,b,name,scale,kits):
    p=profile(q,d,name);s=q['strategies'][name];capacity=d['physical_capacity']
    premium=s['premium_share'];logistics=bool(p['repricing']['propellant_kg_year'])
    extra_W=b['cultural_platform']['power_W']+kits*b['service_kit']['power_W']+premium*p['engineering']['finished_export_kg_year']*scale*q['premium']['additional_mean_W_per_attempted_kg_year']
    ratios=dict(p['engineering']['capacity_use'])
    ratios['ground_firm_power']+=extra_W*capacity['capacity_margin']/capacity['lunar_process_availability']/(scale*capacity['ground_firm_W'])
    ratios['lunar_process_power']+=extra_W/(scale*capacity['peak_flows']['local_mean_W'])
    ratios['depot_power']+=logistics*q['logistics']['additional_depot_use_fraction']
    people=p['engineering']['population']*scale
    staff=p['engineering']['physical']['technical_resident_jobs']*scale+b['cultural_platform']['staff']+kits*b['service_kit']['staff']+premium*q['premium']['staff_per_founding_scale']*scale+logistics*q['logistics']['support_staff_per_founding_scale']*scale
    staff_budget=people*(q['workforce']['total_working_fraction']-q['workforce']['civic_and_other_reserved_fraction'])
    ratios['staff']=staff/staff_budget
    g=p['joint_replication'];limit=g['fraction_of_residents_year']
    for key,used in p['engineering']['capacity_use'].items():
        increment=g['at_limit']['engineering']['capacity_use'][key]-used
        if increment>1e-12:limit=min(limit,g['fraction_of_residents_year']*max(0.,1-ratios[key])/increment)
    feasible=max(ratios.values())<=1+1e-8
    return dict(feasible=feasible,capacity_use=ratios,additional_mean_power_W=extra_W,
        technical_and_new_service_staff=staff,staff_budget=staff_budget,
        growth_fraction_ceiling=max(0.,limit) if feasible else 0.,binding=max(ratios,key=ratios.get))


def operation(q,d,b,name,market,scale,productive,kits,age,previous_paid=0,previous_licenses=0,*,synergies=True,service_demand=True):
    p=profile(q,d,name);strategy=q['strategies'][name];m=q['markets'][market];I=d['inflation_factor']
    premium=q['premium'];logistics=q['logistics'];syn=q['synergies'];competition=q['competition']
    service=staged.service_account(b,market,age,kits,previous_paid,previous_licenses)
    revenues=dict(service['revenue_by_activity_USD2041'])
    if not service_demand:
        revenues={k:0. for k in revenues};service['paid_accounts']=service['licenses']=0.
    maturation=min(1.,b['service_demand_initial_fraction']*(1+b['service_demand_annual_growth'])**min(age,100))
    internal_research=min(kits*b['service_kit']['research_contracts_year'],scale*strategy['premium_share']*syn['qualification_research_slots_per_founding_scale']) if synergies else 0.
    if service_demand:
        research_slots=min(b['scenarios'][market]['research_contracts']*maturation,max(0.,kits*b['service_kit']['research_contracts_year']-internal_research))
        revenues['research']=research_slots*b['scenarios'][market]['research_USD2041_contract']
    qualified=age>=max(0,premium['qualification_years_after_first_production']-(syn['premium_qualification_years_saved'] if synergies and kits else 0))
    accepted=min(1.,premium['accepted_fraction']+(syn['premium_accepted_fraction_gain'] if synergies and kits and premium['accepted_fraction']>0 else 0)) if qualified else 0.
    total_lunar=p['repricing']['delivered_lunar_goods_kg_year']*productive
    attempted_premium=p['engineering']['finished_export_kg_year']*strategy['premium_share']*scale
    premium_feed=total_lunar*strategy['premium_share']
    premium_qty=min(premium_feed*accepted,m['premium_demand_kg_year'])
    standard_qty=total_lunar-premium_feed+(premium_feed-premium_qty)*premium['fallback_recovery_fraction']
    premium_price=premium['equivalent_product_offer_USD2026_kg']*premium['willingness_to_pay_multiple']*competition['premium_competitor_price_factor']*I
    goods=standard_qty*strategy['lunar_price_USD2026_kg']*I+premium_qty*premium_price
    source_goods=p['repricing']['delivered_source_goods_kg_year']*productive*strategy['source_price_USD2026_kg']*I
    industrial_budget=m['industrial_budget_USD2041_year']*min(competition['industrial_budget_cap_multiple'],(1+competition['annual_industrial_budget_growth'])**min(age,1000))
    goods_realization=min(1.,industrial_budget/(goods+source_goods)) if goods+source_goods else 1.
    goods*=goods_realization;source_goods*=goods_realization;premium_qty*=goods_realization
    # Extra branded purchases are distinct transactions, capped by the SAME
    # combined revenue budget as the existing cultural/license programme.
    brand=syn['additional_separate_branded_license_fraction'] if synergies and premium_qty>0 and service_demand else 0.
    revenues['licenses']*=1+brand;service['licenses']*=1+brand
    household=min(1.,m['household_budget_USD2041_year']/max(1.,revenues['culture']+revenues['licenses']))
    revenues['culture']*=household;revenues['licenses']*=household
    service['paid_accounts']*=household;service['licenses']*=household
    share=competition['mission_share_of_research_invoices']
    requested_research=revenues['research']*share
    mission_research=min(requested_research,m['mission_budget_USD2041_year'])
    revenues['research']=revenues['research']-requested_research+mission_research
    # Research and external mission support draw from one budget. Internal fuel
    # consumption already belongs to the provisioning cost response.
    mr=d['market_reference'];launch=competition['Earth_launch_cost_factor']
    competitor=mr['Earth_propellant_competitor_USD2026_kg']*I*(competition['Earth_fuel_price_floor_share_of_delivered_alternative']+(1-competition['Earth_fuel_price_floor_share_of_delivered_alternative'])*launch)
    fuel_price=min(mr['propellant_USD2026_kg']*I,mr['Earth_propellant_price_fraction_ceiling']*competitor)
    fuel_kg=min(p['repricing']['propellant_kg_year']*productive,m['funded_compatible_missions_year']*logistics['delivered_fuel_kg_per_mission'])
    requested_fuel=fuel_kg*fuel_price
    requested_support=fuel_kg/logistics['delivered_fuel_kg_per_mission']*logistics['support_invoice_USD2041_per_mission']
    mission_realization=min(1.,max(0.,m['mission_budget_USD2041_year']-mission_research)/max(1.,requested_fuel+requested_support))
    fuel=requested_fuel*mission_realization;support=requested_support*mission_realization
    # Purchased digital work can be displaced only by capacity withdrawn from
    # outside computing contracts. Internal use carries no second sales invoice.
    compute_capacity=p['repricing']['compute_W']*productive
    needed=syn['internal_compute_W_for_culture']+kits*syn['internal_compute_W_per_service_kit']
    internal=min(compute_capacity,needed) if synergies and service_demand else 0.
    external_compute=min((compute_capacity-internal)*m['compute_capacity_offtake_fraction'],m['compute_contracted_capacity_W'])
    compute=external_compute*mr['compute_USD2026_per_installed_W_year']*I*m['compute_price_factor']
    nonmission_research=revenues['research']-mission_research
    institutional=min(1.,m['institutional_services_budget_USD2041_year']/max(1.,nonmission_research+revenues['professional']+compute))
    revenues['research']=mission_research+nonmission_research*institutional
    revenues['professional']*=institutional;compute*=institutional;external_compute*=institutional
    variable={key:revenues[key]*b['variable_cost_fractions'][key]*competition['service_cost_multiplier'] for key in revenues}
    digital_saving=(variable['culture']+variable['licenses'])*syn['avoided_external_digital_cost_fraction']*internal/needed
    acquisition=max(0.,service['paid_accounts']-previous_paid*(1-b['paid_churn']))*b['paid_acquisition_USD2041']+max(0.,service['licenses']-previous_licenses*(1-b['license_renewal_churn']))*b['license_acquisition_USD2041']
    service_cost=sum(variable.values())-digital_saving+acquisition+service['wages_USD2041']+service['fixed_purchases_USD2041']+service['renewal_USD2041']
    a=p['external'];v=a['USD2041'];old=a['sales_by_activity_USD2041']
    removed_cost=sum(old[key]*fraction for key,fraction in b['replaced_baseline_Earth_service_fractions'].items())
    imports=v['total_imports']*(scale+b['common_recurring_fraction']*max(0.,1-scale))-removed_cost*scale
    imports+=mr['compute_Earth_service_fraction']*(compute-old['batch_computing']*scale)
    # Linked launch competition changes the identified annual launch purchase,
    # not every equipment cost or the unchanged original establishment bill.
    imports+=p['engineering']['cash_USD2026']['Earth_launch_service']*I*scale*(launch-1)
    provision=v['resource_provision']*scale
    extra_premium_cost=strategy['premium_share']*scale*(premium['annual_development_USD2041_per_founding_scale']+premium['installed_capital_USD2041_per_founding_scale']/premium['life_years']+premium['staff_per_founding_scale']*b['service_staff_USD2041_year'])+attempted_premium*premium['incremental_cost_USD2041_attempted_kg']
    extra_support_cost=support*logistics['support_variable_cost_fraction']+bool(p['repricing']['propellant_kg_year'])*scale*(logistics['support_fixed_USD2041_per_founding_scale_year']+logistics['support_capital_USD2041_per_founding_scale']/logistics['life_years']+logistics['support_staff_per_founding_scale']*b['service_staff_USD2041_year'])
    destination=old.get('destination_services',0)*productive
    sectors=dict(materials=goods,source_products=source_goods,propellant=fuel,mission_support=support,computing=compute,destination=destination,
        sports=revenues['culture']*.45,entertainment=revenues['culture']*.55,licenses=revenues['licenses'],research=revenues['research'],professional=revenues['professional'])
    cost=imports+provision+service_cost+extra_premium_cost+extra_support_cost
    return dict(revenue_by_sector_USD2041=sectors,gross_USD2041=sum(sectors.values()),total_cost_USD2041=cost,
        available_USD2041=sum(sectors.values())-cost,core_imports_USD2041=imports,resource_provision_USD2041=provision,
        service_cost_USD2041=service_cost,premium_cost_USD2041=extra_premium_cost,logistics_cost_USD2041=extra_support_cost,
        audience=service['audience'] if service_demand else 0,paid_accounts=service['paid_accounts'],licenses=service['licenses'],
        premium_qualified=qualified,premium_accepted_fraction=accepted,premium_sold_kg=premium_qty,
        industrial_budget_USD2041=industrial_budget,industrial_realization=goods_realization,
        household_realization=household,mission_research_USD2041=mission_research,
        mission_budget_USD2041=m['mission_budget_USD2041_year'],mission_realization=mission_realization,
        delivered_propellant_kg=fuel_kg*mission_realization,propellant_price_USD2041_kg=fuel_price,
        internal_compute_W=internal,external_compute_W=external_compute,
        internal_research_slots=internal_research,institutional_realization=institutional,
        avoided_external_digital_cost_USD2041=digital_saving,branded_license_fraction=brand,
        resource_screen=resource_screen(q,d,b,name,scale,kits))


def simulate(q,d,b,e,t,name,market,horizon,*,synergies=True):
    p=profile(q,d,name);cap=capital(q,d,b,e,t,name);seed=b['seed_scale'];start=b['bootstrap_years']
    people=p['engineering']['population'];maximum=e['growth']['0.01']['vacuum_architectural_reference_with_existing_cities_people']/people
    if not p['engineering']['feasible'] or not p['commissioning']['feasible_overlap']:raise ValueError('Rejected source portfolio')
    if not resource_screen(q,d,b,name,seed,0)['feasible']:raise ValueError('Cultural/industrial seed exceeds shared capacity')
    kits=1 if resource_screen(q,d,b,name,seed,1)['feasible'] else 0
    common=b['common_capital_fraction_of_establishment']*cap['establishment_USD2041']
    hardware=common+(cap['one_time_bill_USD2041']-common)*seed+cap['common_development_USD2041']+b['cultural_platform']['capital_USD2041']+kits*b['service_kit']['capital_USD2041']
    def burn(scale,units):
        return operation(q,d,b,name,market,scale,0,units,0,synergies=False,service_demand=False)['total_cost_USD2041']
    no_sales=burn(seed,kits);working=no_sales*b['seed_no_receipt_runway_years'];carry=no_sales*b['bootstrap_carrying_fraction']*start
    flows=[-(hardware+carry)/start]*start;flows[-1]-=working
    cash=working;earned=0.;cohorts=[(start,seed)];pending=[];kit_pending=[];paid=licenses=0.
    rows=[];failure=None;first_full=reached_cap=None
    for year in range(start,horizon+1):
        cohorts.extend((due,n) for due,n in pending if due==year);pending=[a for a in pending if a[0]>year]
        kits+=sum(n for due,n in kit_pending if due==year);kit_pending=[a for a in kit_pending if a[0]>year]
        scale=sum(n for _,n in cohorts);productive=sum(n*min(1.,(year-due+1)/b['production_ramp_years']) for due,n in cohorts)
        a=operation(q,d,b,name,market,scale,productive,kits,year-start,paid,licenses,synergies=synergies)
        if not a['resource_screen']['feasible']:raise ValueError('Shared operating capacity exceeded')
        paid=a['paid_accounts'];licenses=a['licenses'];opening=cash
        pooled=min(0,a['available_USD2041'])+max(0,a['available_USD2041'])*b['collectible_fraction']
        cash+=pooled;earned+=pooled
        buffer=burn(scale,kits)*max(b['liquidity_buffer_years'],b['seed_no_receipt_runway_years']-(year-start))
        industrial=service_orders=dividend=starts=kit_saving=0.
        if cash<0:failure={'year':year,'shortfall_USD2041':-cash}
        else:
            dividend=min(max(0.,min(cash-buffer,earned)),max(0.,pooled)*(1-b['retained_fraction']))
            cash-=dividend;earned-=dividend
            total_kits=kits+sum(n for _,n in kit_pending)
            desired=min(b['service_kit']['max_units'],math.ceil(max(
                b['scenarios'][market]['research_contracts']/b['service_kit']['research_contracts_year'],
                b['scenarios'][market]['professional_clients']/b['service_kit']['professional_clients_year'])*min(1.,b['service_demand_initial_fraction']*(1+b['service_demand_annual_growth'])**min(year-start,100))))
            # Commissioning reserve is explicit; pending facilities do not earn.
            candidate=resource_screen(q,d,b,name,scale,total_kits+1+b['commissioning_years']*q['additional_commissioning_capacity_fraction'])
            if total_kits<desired and candidate['feasible']:
                kit_saving=b['service_kit']['capital_USD2041']
                if min(cash-buffer,earned)>=kit_saving:
                    service_orders=kit_saving;cash-=kit_saving;earned-=kit_saving;kit_saving=0.
                    kit_pending.append((year+b['commissioning_years'],1))
            reserved_kits=kits+sum(n for _,n in kit_pending)*(1+q['additional_commissioning_capacity_fraction'])
            available_resources=resource_screen(q,d,b,name,scale,reserved_kits)
            committed=scale+sum(n for _,n in pending)
            slope=b['desired_growth_rate']*scale
            if committed>b['taper_fraction']*maximum:
                slope=b['desired_growth_rate']*b['taper_fraction']*maximum*math.sqrt(max(0.,(maximum-committed)/((1-b['taper_fraction'])*maximum)))
            limit=max(0.,min(slope,available_resources['growth_fraction_ceiling']*scale,maximum-committed))
            # Marginal return uses current qualification and customer budgets,
            # keeping kits fixed. Anticipated future premium success cannot gate
            # today's industrial order.
            now=operation(q,d,b,name,market,committed,committed,kits,year-start,paid,licenses,synergies=synergies)
            plus=operation(q,d,b,name,market,committed+1e-5,committed+1e-5,kits,year-start,paid,licenses,synergies=synergies)
            joint_reserve=0.
            if plus['available_USD2041']<=now['available_USD2041']:
                limit=0.
                # A service business may justify its supporting infrastructure.
                # Fund the complete next host-plus-kit programme, including a
                # conservative no-receipt carrying allowance, before starting it.
                if total_kits<desired:
                    lo,hi=committed,maximum
                    kit_load=total_kits+1+b['commissioning_years']*q['additional_commissioning_capacity_fraction']
                    if resource_screen(q,d,b,name,hi,kit_load)['feasible']:
                        for _ in range(35):
                            mid=(lo+hi)/2
                            if resource_screen(q,d,b,name,mid,kit_load)['feasible']:hi=mid
                            else:lo=mid
                        target=hi*(1+1e-8)
                        future=operation(q,d,b,name,market,target,target,total_kits+1,year-start,paid,licenses,synergies=synergies)
                        current=operation(q,d,b,name,market,committed,committed,total_kits,year-start,paid,licenses,synergies=synergies)
                        if future['available_USD2041']>current['available_USD2041'] and target<=maximum:
                            incremental_burn=max(0.,burn(target,total_kits+1)-burn(committed,total_kits))
                            joint_reserve=staged.expansion_cost(b,cap,committed,max(0.,target-committed))+b['service_kit']['capital_USD2041']+incremental_burn*(b['commissioning_years']+b['production_ramp_years']+b['liquidity_buffer_years'])
                            if min(cash-buffer,earned)>=joint_reserve:
                                limit=max(0.,min(slope,available_resources['growth_fraction_ceiling']*scale,target-committed))
            investable=max(0.,min(cash-buffer,earned)-kit_saving)
            starts=staged.affordable_increment(b,cap,committed,limit,investable)
            if starts>1e-12:
                industrial=staged.expansion_cost(b,cap,committed,starts);cash-=industrial;earned-=industrial
                pending.append((year+b['commissioning_years'],starts))
            else:starts=0.
            reserve=max(kit_saving,max(0.,joint_reserve-industrial))
            excess=max(0.,min(cash-buffer,earned)-reserve);dividend+=excess;cash-=excess;earned-=excess
        flows.append(dividend)
        rows.append(dict(year=year,scale=scale,productive_scale=productive,service_kits=kits,
            opening_cash_USD2041=opening,pooled_USD2041=pooled,closing_cash_USD2041=cash,earned_USD2041=earned,
            industrial_orders_USD2041=industrial,service_orders_USD2041=service_orders,dividend_USD2041=dividend,
            starts_scale=starts,pipeline_scale=sum(n for _,n in pending),buffer_USD2041=buffer,
            joint_programme_reserve_USD2041=joint_reserve if not failure else 0,
            growth_fraction_ceiling=available_resources['growth_fraction_ceiling'] if not failure else 0,operation=a))
        if scale>=1-1e-8 and first_full is None:first_full=year
        if scale>=maximum-1e-8 and reached_cap is None:reached_cap=year
        if failure:break
    terminal=max(0.,cash-buffer) if not failure else 0.;flows[-1]+=terminal
    total=0.;balances=[]
    for flow in flows:total+=flow;balances.append(total)
    recovery=next((i for i,value in enumerate(balances) if value>=0 and min(balances[i:])>=0),None)
    summary=dict(strategy=name,market=market,horizon_elapsed_years=horizon,synergies=synergies,
        external_seed_USD2041=hardware+working+carry,seed_hardware_USD2041=hardware,
        bootstrap_carry_USD2041=carry,seed_liquidity_USD2041=working,
        external_funding_after_first_production_USD2041=0,failure=failure,
        NPV_USD2041={str(r):run.pv(flows,r) for r in q['rates']},durable_recovery_year=recovery,
        first_founding_scale_year=first_full,cap_year=reached_cap,final_scale=rows[-1]['scale'],
        final_available_USD2041=rows[-1]['operation']['available_USD2041'],
        final_gross_USD2041=rows[-1]['operation']['gross_USD2041'],
        total_industrial_orders_USD2041=sum(a['industrial_orders_USD2041'] for a in rows),
        total_service_orders_USD2041=sum(a['service_orders_USD2041'] for a in rows),
        terminal_liquid_USD2041=terminal,retained_operating_liquidity_USD2041=min(max(0.,cash),buffer))
    return {'summary':summary,'annual':rows,'capital_provider_flows_USD2041':flows}


def results():
    q,d,b,e,t=inputs();horizons=[];checkpoints=[];steady=[];cap_checks=[]
    for name in q['strategies']:
        for market in q['markets']:
            for h in q['horizons']:
                x=simulate(q,d,b,e,t,name,market,h);horizons.append(x['summary'])
                if h==max(q['horizons']):
                    checkpoints.append(dict(strategy=name,market=market,rows=[a for a in x['annual'] if a['year'] in (8,10,16,20,30,50,100,200,300,500,800) or a is x['annual'][-1]]))
            units=min(b['service_kit']['max_units'],math.ceil(max(b['scenarios'][market]['research_contracts']/b['service_kit']['research_contracts_year'],b['scenarios'][market]['professional_clients']/b['service_kit']['professional_clients_year'])))
            while units and not resource_screen(q,d,b,name,1,units)['feasible']:units-=1
            first=operation(q,d,b,name,market,1,1,units,30)
            a=operation(q,d,b,name,market,1,1,units,30,first['paid_accounts'],first['licenses'])
            alone=operation(q,d,b,name,market,1,1,0,30,synergies=False,service_demand=False)
            counterfactual=operation(q,d,b,name,market,1,1,units,30,synergies=False)
            without=operation(q,d,b,name,market,1,1,units,30,counterfactual['paid_accounts'],counterfactual['licenses'],synergies=False)
            steady.append(dict(strategy=name,market=market,service_kits=units,operation=a,
                without_worldwide_service_demand_USD2041=alone['available_USD2041'],
                without_synergies_USD2041=without['available_USD2041']))
            cap_scale=e['growth']['0.01']['vacuum_architectural_reference_with_existing_cities_people']/p_people(d)
            initial_cap=operation(q,d,b,name,market,cap_scale,cap_scale,units,792)
            plateau=operation(q,d,b,name,market,cap_scale,cap_scale,units,792,initial_cap['paid_accounts'],initial_cap['licenses'])
            cap_checks.append(dict(strategy=name,market=market,installed_scale=cap_scale,
                available_USD2041_year=plateau['available_USD2041'],gross_USD2041_year=plateau['gross_USD2041'],
                feasible_under_assumed_replication=plateau['resource_screen']['feasible']))
    sensitivity=[]
    for variant in ('no_synergies','premium_fails','Earth_competition','half_customer_budgets','higher_service_costs'):
        v=deepcopy(q)
        if variant=='premium_fails':v['premium']['accepted_fraction']=0
        if variant=='Earth_competition':
            v['competition']['Earth_launch_cost_factor']=.25;v['competition']['premium_competitor_price_factor']=.7
            for market in v['markets'].values():market['compute_price_factor']*=.5
        if variant=='half_customer_budgets':
            for market in v['markets'].values():
                for key in ('industrial_budget_USD2041_year','mission_budget_USD2041_year','household_budget_USD2041_year','premium_demand_kg_year','institutional_services_budget_USD2041_year','funded_compatible_missions_year','compute_contracted_capacity_W'):market[key]*=.5
        if variant=='higher_service_costs':v['competition']['service_cost_multiplier']=1.3
        for name in ('premium_lunar','propellant_logistics','diversified','compute_digital'):
            sensitivity.append(dict(variant=variant,**simulate(v,d,b,e,t,name,'breakthrough',100,synergies=variant!='no_synergies')['summary']))
    return dict(schema='terluna.research.export-strategy-returns/1',evidence=q['evidence'],
        producer_sha256=run.digest(HERE/'strategies.py'),input_sha256=run.digest(HERE/'strategy_inputs.json'),
        dependencies={p:run.digest(HERE/p) for p in ('staged.py','staged_inputs.json','run.py','returns_inputs.json')},
        horizons=horizons,checkpoints=checkpoints,steady_comparisons=steady,sensitivities=sensitivity,cap_operation_checks=cap_checks,
        rejected_independent_maxima={'binding':d['rejected_controls']['all_individual_maxima']['engineering']['binding_capacity'],
            'maximum_use':max(d['rejected_controls']['all_individual_maxima']['engineering']['capacity_use'].values())})


def p_people(d):return d['baseline']['engineering']['population']


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=run.ROOT/'research/runs/complex_returns/strategies.json')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n');print(args.output)
