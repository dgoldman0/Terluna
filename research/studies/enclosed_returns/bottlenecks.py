"""Stratify saved returns and compare independently solved interventions.

This consumer prices no machinery and never treats unused buyer budgets as
contracted sales. Common physical costs remain shared across export sectors.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import gzip
import json
import math
from pathlib import Path
from . import finance, markets

HERE = Path(__file__).resolve().parent


def pv_amounts(values, rate):
    return math.fsum(amount / (1 + rate)**year for year, amount in values)


def annual_costs(request, physical, row, cfg):
    cash = physical['cash']
    costs = {key: cash[key] for key in (
        'operating_external_replacement_USD2041', 'operating_external_logistics_materials_USD2041',
        'operating_external_services_USD2041', 'operating_local_payroll_USD2041',
        'operating_local_replacement_USD2041', 'operating_resource_provision_USD2041',
        'capital_external_USD2041', 'capital_local_USD2041', 'capital_resource_provision_USD2041')}
    costs.update(commercial_operating_USD2041=request['commercial_operating_USD2041'],
                 hotel_payroll_USD2041=request['visitor_beds'] * cfg['hotel_staff_per_bed'] * cfg['hotel_staff_USD2041_year'],
                 residual_service_renewal_USD2041=row['service_residual_renewal'],
                 residual_service_capital_USD2041=row['service_residual_capital'],
                 capital_allowance_USD2041=row['capital_cost_allowance'])
    if not math.isclose(sum(costs.values()), row['capital'] + row['operating_cost'], rel_tol=1e-12, abs_tol=.01):
        raise ValueError('Cost stratification does not reconcile')
    return costs


def constraint_row(request, physical, cfg, b):
    y = request['year']
    age = y-cfg['first_productive_year']
    maturity = min(1., b['service_demand_initial_fraction'] * (1+b['service_demand_annual_growth'])**min(age, 100))
    scenario = b['scenarios'][request['_market']]
    wanted = math.ceil(max(scenario['research_contracts']/b['service_kit']['research_contracts_year'],
                           scenario['professional_clients']/b['service_kit']['professional_clients_year'])
                       * maturity * request['demand_growth_factor'])
    staff = max(1, math.floor((cfg['service_staff_fraction'] * request['population'] - b['cultural_platform']['staff'])/b['service_kit']['staff']))
    rev = request['revenue_by_sector_USD2041']
    sold = dict(industrial=rev['materials']+rev['source_products'],
                mission=request['mission_research_USD2041']+rev['propellant']+rev['mission_support'],
                institutional=rev['research']-request['mission_research_USD2041']+rev['professional']+rev['computing'],
                household=rev['sports']+rev['entertainment']+rev['licenses']+rev['destination'])
    budgets = {k: request[k+'_budget_USD2041'] for k in sold}
    return dict(year=y, population=request['population'], installed_service_packages=request['service_kit_count'],
                demand_desired_packages_before_budget_competition=wanted, service_staff_planning_limit=staff,
                total_workforce_fraction_of_allowed_pool=physical['capacity_summary']['resident_workforce'],
                maximum_capacity_requirement_to_available_ratio=physical['maximum_capacity_use'],
                admitted_sales_by_customer_budget=sold, buyer_budget=budgets,
                buyer_budget_used_fraction={k:sold[k]/budgets[k] if budgets[k] else None for k in sold},
                offered_sales_budget_admission_fraction={k:request[k+'_realization'] for k in sold})


def funder_measures(path, sim, cfg, horizon):
    start=path['annual'][0]['year']
    contributions=[(y,path['seed_hardware_and_carry']/start) for y in range(start)]
    contributions += [(0,sim['initial_reserve']), (start-1,path['initial_working'])]
    fs=[f for f in sim['annual'] if f['year']<=horizon]
    contributions += [(f['year'],f['contribution']) for f in fs]
    payouts=[(f['year'],f['owner_payout']) for f in fs]
    cf=[0.]*(horizon+1)
    for y,v in contributions:cf[y]-=v
    for y,v in payouts:cf[y]+=v
    # The baseline uniqueness check can see sub-dollar floating top-ups as
    # extra sign changes. This diagnostic explicitly screens that roundoff.
    irr_cf=[0. if abs(v)<1. else v for v in cf]
    record=finance.metrics(path,sim,cfg,horizon)
    rates={}
    for r in cfg['rates']:
        pin=pv_amounts(contributions,r);pout=pv_amounts(payouts,r)
        rates[str(r)]=dict(PV_all_funder_contributions=pin,PV_distributions=pout,
                          PV_distribution_multiple=pout/pin,net_present_value=pout-pin,
                          baseline_record_residual=pout-pin-record['discounted'][str(r)]['sponsor_distributed_cash_NPV'])
    total=sum(v for _,v in contributions)
    return dict(horizon=horizon,initial_full_funding=path['seed']+sim['initial_reserve'],
                priced_seed_deployment=path['seed_direct_cost'],seed_capital_allowance=path['seed_capital_cost_allowance'],
                initial_working_cash=path['initial_working'],initial_monetary_reserve=sim['initial_reserve'],
                all_funder_contributions=total,total_distributions=sum(v for _,v in payouts),
                undiscounted_distribution_multiple=sum(v for _,v in payouts)/total,
                first_distribution_year=next((f['year'] for f in fs if f['owner_payout']>1.),None),
                real_IRR_roundoff_screened=finance.valuation.nonnegative_irr(irr_cf),
                IRR_roundoff_tolerance_USD2041=1.,discounted=rates,checkpoint=record)


def carry_bridge(path, sim, cfg, horizon, rate):
    """Reconcile PV of real cash to payouts plus terminal net financial claims."""
    start=path['annual'][0]['year']
    assets=0.;debt=0.;carry=interest_vs_discount=0.
    fs={f['year']:f for f in sim['annual']}
    for y in range(horizon+1):
        if y:
            carry += rate*assets/(1+rate)**y
            if y in fs:
                interest_vs_discount += (fs[y]['interest']-rate*debt)/(1+rate)**y
        if y==0:assets+=sim['initial_reserve']
        if y==start-1:assets+=path['initial_working']
        if y in fs:
            f=fs[y];assets=f['cash']+f['reserves']+f['owner_custody'];debt=f['debt']
    metric=finance.metrics(path,sim,cfg,horizon)['discounted'][str(rate)]
    wedge=metric['unlevered_cash_NPV']-metric['sponsor_cash_and_terminal_NPV']
    residual=wedge-carry-interest_vs_discount
    if not math.isclose(residual,0.,abs_tol=.25):
        raise ValueError(f'Financial carry bridge failed by {residual}')
    return dict(unlevered_cash_NPV=metric['unlevered_cash_NPV'],sponsor_distributed_cash_NPV=metric['sponsor_distributed_cash_NPV'],
                terminal_net_financial_claim_PV=metric['terminal_financial_claim_PV'],
                zero_real_yield_asset_carry_cost_PV=carry,interest_minus_discount_cost_on_debt_PV=interest_vs_discount,
                reconciliation_residual_USD2041=residual)


def stratify(case, physical, cfg, periods):
    path=finance.plan(case,physical,cfg)
    eq=finance.finance(path,cfg,'equity');sim=finance.finance(path,cfg,'U_credit',eq)
    mc,_,_,b,_,_=markets.inputs()
    rows=[]
    for request,pr,r,f in zip(case['annual'],physical['annual'],path['annual'],sim['annual']):
        costs=annual_costs(request,pr,r,cfg)
        rows.append(dict(year=r['year'],revenue=request['revenue_by_sector_USD2041'],costs=costs,
                         unlevered=r['operating_available']-r['capital'],funder_payout=f['owner_payout'],
                         extra_equity=f['contribution'],cash=f['cash'],reserve=f['reserves']))
    grouped=[]
    for first,last in periods:
        selected=[r for r in rows if first<=r['year']<=last]
        values={}
        for rate in (0.,.03):
            revenue={k:pv_amounts([(r['year'],r['revenue'][k]) for r in selected],rate) for k in rows[0]['revenue']}
            costs={k:pv_amounts([(r['year'],r['costs'][k]) for r in selected],rate) for k in rows[0]['costs']}
            seed=pv_amounts([(y,path['seed_hardware_and_carry']/mc['first_productive_year']) for y in range(mc['first_productive_year']) if first<=y<=last],rate)
            costs['seed_deployment_and_allowance_USD2041']=seed
            values[str(rate)]=dict(revenue_by_sector=revenue,costs_by_category=costs,
                                  net_real_cash=sum(revenue.values())-sum(costs.values()),
                                  distributions=pv_amounts([(r['year'],r['funder_payout']) for r in selected],rate),
                                  additional_equity=pv_amounts([(r['year'],r['extra_equity']) for r in selected],rate))
        grouped.append(dict(first_year=first,last_year=last,values=values))
    constraints=[constraint_row(dict(request,_market=case['market']),pr,mc,b)
                 for request,pr in zip(case['annual'],physical['annual']) if request['year'] in (8,20,50,100,150,200,220)]
    return dict(name=case['name'],portfolio=case['portfolio'],market=case['market'],growth=case['market_growth_name'],
                priced_physical_seed=physical['seed']['cost']['total_USD2041'],periods=grouped,constraints=constraints,
                checkpoints=[funder_measures(path,sim,cfg,h) for h in cfg['horizons']],
                financing_bridges={str(h):carry_bridge(path,sim,cfg,h,.03) for h in cfg['horizons']},
                peak_retained_cash=max((dict(year=r['year'],cash=r['cash']) for r in rows),key=lambda x:x['cash']))


def run_intervention(case, physical, cfg):
    if not physical['convergence']['converged'] or physical['first_capacity_shortfall'] is not None:
        return dict(name=case['name'],experiment=case['experiment'],status='physical_requirements_unresolved',
                    convergence=physical['convergence'],first_capacity_shortfall=physical['first_capacity_shortfall'],
                    portfolio=case['portfolio'],workforce_failures=[dict(year=r['year'],population=r['population'],required_jobs=r['technical_resident_jobs'],
                        available_jobs=.4*r['population'],workforce_ratio=r['capacity_summary']['resident_workforce'],core_jobs=r['core_resident_jobs'],commercial_jobs=r['commercial_resident_jobs'])
                        for r in physical['annual'] if r['capacity_summary']['resident_workforce']>1+1e-6])
    path=finance.plan(case,physical,cfg);eq=finance.finance(path,cfg,'equity');sim=finance.finance(path,cfg,'U_credit',eq)
    return dict(name=case['name'],portfolio=case['portfolio'],experiment=case['experiment'],status='conditional_account',
                peak_workforce_fraction_of_allowed_pool=max(r['capacity_summary']['resident_workforce'] for r in physical['annual']),
                checkpoints=[funder_measures(path,sim,cfg,h) for h in cfg['horizons']],
                physical_checkpoints=[dict(year=r['year'],service_packages=next(a['service_kit_count'] for a in case['annual'] if a['year']==r['year']),
                    total_physical_operation=r['cash']['operating_total_USD2041'],
                    imported_component_replacement=r['cash']['operating_external_replacement_USD2041'],
                    source_feed_kg_year=r['physical']['source_feed_kg_year']) for r in physical['annual'] if r['year'] in (20,50,100,200,220)])


def margin_target(case, physical, cfg, reduction=.05):
    """Price target on purchased commercial inputs; no invented physical output."""
    path=finance.plan(case,physical,cfg)
    _,q,_,commercial,_,_=markets.inputs()
    requests={r['year']:r for r in case['annual']};savings=[]
    for r in path['annual']:
        covered=sum(r['revenue_by_sector'][k] for k in ('sports','entertainment','licenses','research','professional'))
        saving=reduction*covered
        request=requests[r['year']];v=commercial['variable_cost_fractions'];rev=r['revenue_by_sector']
        culture=rev['sports']+rev['entertainment'];syn=q['synergies']
        internal_need=syn['internal_compute_W_for_culture']+request['service_kit_count']*syn['internal_compute_W_per_service_kit']
        digital_saving=(culture*v['culture']+rev['licenses']*v['licenses'])*syn['avoided_external_digital_cost_fraction']*request['internal_compute_W']/internal_need
        covered_inputs=culture*v['culture']+rev['licenses']*v['licenses']+rev['research']*v['research']+rev['professional']*v['professional']-digital_saving
        if saving>covered_inputs:
            raise ValueError('Margin target exceeds covered purchased commercial inputs')
        savings.append(dict(year=r['year'],covered_gross=covered,covered_variable_purchases=covered_inputs,cost_reduction=saving))
        r['operating_cost']-=saving;r['external_operating']-=saving;r['operating_available']+=saving
    for i,r in enumerate(path['annual']):
        future=path['annual'][i+1:i+1+cfg['underwriting']['project_appraisal_years']]
        r['programme_NPV_screen']=-r['capital']+math.fsum((v['operating_available']-v['capital'])/(1+cfg['underwriting']['project_discount_rate'])**(j+1) for j,v in enumerate(future))
    eq=finance.finance(path,cfg,'equity');sim=finance.finance(path,cfg,'U_credit',eq)
    financial=[funder_measures(path,sim,cfg,h) for h in cfg['horizons']]
    allowed={'sponsor_distributed_cash_NPV','sponsor_cash_and_terminal_NPV','terminal_financial_claim_PV','unlevered_cash_NPV','additional_initial_cost_at_unlevered_break_even','durable_distributed_cash_recovery_year','local_interest_relief_and_bank_margin_loss'}
    for record in financial:
        for rate,values in record['checkpoint']['discounted'].items():
            record['checkpoint']['discounted'][rate]={k:v for k,v in values.items() if k in allowed}
    return dict(name=case['name'],cost_reduction_fraction_of_covered_gross=reduction,annual_cost_reductions=savings,socioeconomic_change=None,
                evidence='Five percentage points of covered gross invoices are removed from purchased commercial input costs, at unchanged physical output and seed funding. This is an unpriced efficiency or bargaining target; its qualification costs and suppliers\' losses are not appraised. The result is a sponsor-cash sensitivity, not an estimate of additional social welfare.',
                checkpoints=financial)


def read_product(spec):
    p=HERE/spec['path']
    if markets.digest(p)!=spec['sha256']:raise ValueError('Unadmitted diagnostic product: '+str(p))
    raw=p.read_bytes();data=json.loads(gzip.decompress(raw) if p.suffix=='.gz' else raw)
    if data['schema']!=spec['schema']:raise ValueError('Unexpected diagnostic product schema')
    return data


def product():
    spec=json.loads((HERE/'bottleneck_inputs.json').read_text())
    for p,h in spec['baseline_dependencies'].items():
        if markets.digest(HERE/p)!=h:raise ValueError('Changed diagnostic baseline: '+p)
    cfg,orders,products=finance.inputs()
    by={p['name']:p for p in products['physical']['cases']}
    base=[stratify(c,by[c['name']],cfg,spec['periods']) for c in orders['cases']]
    admitted=json.loads((HERE/'bottleneck_products.json').read_text())
    experimental_orders=read_product(admitted['orders']);physical=read_product(admitted['physical'])
    mining=read_product(admitted['resources'])
    # Product coupling is checked explicitly, independently of filenames.
    if physical['orderbook']['sha256']!=admitted['orders']['sha256']:
        raise ValueError('Physical intervention used different orders')
    if not any(p['sha256']==admitted['physical']['sha256'] for p in mining['physical_products']):
        raise ValueError('Resource diagnostic used different physical schedules')
    by_physical={c['name']:c for c in physical['cases']}
    comparisons=[]
    for c in experimental_orders['cases']:
        record=run_intervention(c,by_physical[c['name']],cfg)
        pc=by_physical[c['name']]
        record['resources']=[dict(horizon=h,**finance.resource_condition(pc,mining,h)) for h in cfg['horizons']] if pc['annual'] else None
        comparisons.append(record)
    pairs=[]
    for portfolio in spec['experimental_portfolios']:
        cases={x['experiment']['variant']:x for x in comparisons if x.get('portfolio')==portfolio}
        for h in cfg['horizons']:
            metrics={k:next(r['checkpoint'] for r in v['checkpoints'] if r['horizon']==h) for k,v in cases.items() if v['status']=='conditional_account'}
            if 'baseline' not in metrics:continue
            values={k:r['discounted']['0.03']['sponsor_distributed_cash_NPV'] for k,r in metrics.items()}
            unlevered={k:r['discounted']['0.03']['unlevered_cash_NPV'] for k,r in metrics.items()}
            pairs.append(dict(portfolio=portfolio,horizon=h,NPV_delta_from_baseline={k:v-values['baseline'] for k,v in values.items()},
                              unlevered_NPV_delta_from_baseline={k:v-unlevered['baseline'] for k,v in unlevered.items()},
                              demand_rollout_interaction=(values['demand_and_rollout']-values['demand_25']-values['combined_rollout']+values['baseline']) if all(k in values for k in ('demand_and_rollout','demand_25','combined_rollout')) else None))
    targets=[margin_target(c,by[c['name']],cfg) for c in orders['cases'] if c['portfolio'] in spec['experimental_portfolios'] and c['market']==spec['market'] and c['market_growth_name']==spec['growth']]
    return dict(schema='terluna.humanities.enclosed-bottlenecks/1',producer=dict(path='research/studies/enclosed_returns/bottlenecks.py',sha256=markets.digest(__file__)),
                dependencies={p:markets.digest(HERE/p) for p in ('bottleneck_inputs.json','bottleneck_products.json','finance.py','markets.py',*[s['path'] for s in admitted.values()])},
                evidence=spec['evidence'],rules=spec['rules'],baseline_stratification=base,interventions=comparisons,
                intervention_differences=pairs,commercial_margin_targets=targets,
                discount_examples={str(y):(1.03)**(-y) for y in (20,50,100,150,200)},
                additional_annual_net_cash_for_10T_PV={str(start):1e13/math.fsum(1.03**(-y) for y in range(start,201)) for start in (20,50,100,150)},
                continuing_value_screen='No terminal enterprise value is included. Existing cash, reserves and debt are already reported separately. A terminal continuation sensitivity requires a qualified post-220 renewal and source plan; a revenue multiple would not establish it.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE.parents[1]/'runs/enclosed_returns/bottlenecks.json')
    args=parser.parse_args();result=product();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
    print('Stratified baselines:',len(result['baseline_stratification']),'interventions:',len(result['interventions']))
    print(args.output)


if __name__=='__main__':main()
