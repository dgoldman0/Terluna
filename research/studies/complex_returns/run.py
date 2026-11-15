"""Finite-horizon stakeholder appraisal of the founding off-Earth complex.

Reads frozen provisioning products, not physical model implementations. Monetary
flows use constant USD2041. Cash recovery is not investor entitlement or welfare.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def inputs():
    q = json.loads((HERE/'returns_inputs.json').read_text())
    products = {}
    for name, identity in q['products'].items():
        p = HERE/identity['snapshot']
        if digest(p) != identity['sha256']:
            raise ValueError(f'Unreviewed input change: {name}')
        products[name] = json.loads(p.read_text())
        if products[name]['schema'] != identity['schema']:
            raise ValueError('Unexpected schema')
    return q, products['economy'], products['industrial']


def pv(flows, rate):
    if rate < 0 or not math.isfinite(rate):
        raise ValueError('Rate must be finite and nonnegative')
    return sum(x/(1+rate)**t for t,x in enumerate(flows))


def target_population(t, initial, cap, rate, taper_fraction):
    """Exponential then quadratic approach matching slope, ending at fixed cap."""
    if not 0 < initial <= cap or not 0 < taper_fraction < 1 or rate < 0:
        raise ValueError('Invalid growth inputs')
    if rate == 0 or cap == initial:
        return initial
    start = max(initial, taper_fraction*cap)
    onset = math.log(start/initial)/rate
    if t <= onset:
        return initial*math.exp(rate*t)
    duration = 2*(cap-start)/(rate*start)
    u = min(1., (t-onset)/duration)
    return cap-(cap-start)*(1-u)**2


def growth_path(q, economy, rate, years):
    ref = economy['growth']['0.01']
    initial = ref['years'][0]['occupied_capacity_people']
    cap = ref['vacuum_architectural_reference_with_existing_cities_people']
    frontier = economy['high_value_capacity_with_growth']['engineering']['growth_people_year']/initial
    delay = q['commissioning_delay_years']
    population = initial
    pending = {}
    rows = []
    for t in range(1,years+1):
        population += pending.pop(t,0)
        desired = target_population(t,initial,cap,rate,q['taper_fraction'])-target_population(t-1,initial,cap,rate,q['taper_fraction'])
        starts = min(desired, frontier*population, max(0,cap-population-sum(pending.values())))
        pending[t+delay] = starts
        rows.append(dict(year=t,population=population,scale=population/initial,starts=starts,
                         pipeline=sum(pending.values()),capacity_fraction=starts/(frontier*population)))
    return rows


def capital_schedule(total, years):
    # A timing sensitivity, not a claimed dated procurement plan.
    return [total/years]*years


def evaluate(q, economy, portfolio, horizon, rate, growth=0, demand_cap=None,
             capture=1, service_value=1, counterfactual_scale=1, delay=0, export_value=1, spillover=0):
    if not 0 <= capture <= 1 or service_value < 0 or counterfactual_scale < 0:
        raise ValueError('Invalid valuation or collection fraction')
    row = (economy['high_value_capacity_case'] if portfolio=='large' else economy['portfolios'][portfolio])
    ext = row['external']['USD2041']; dom = row['domestic']['USD2041']
    commission = (economy['high_value_capacity_commissioning'] if portfolio=='large' else economy['portfolio_commissioning'])
    initial_capital = economy['initial_capital_context']['inherited_industrial_and_material_campaign_USD2041']
    initial_capital *= q['campaign_scope_multiplier']
    commercial = sum(commission[k] for k in ['additional_external_purchases_USD2041','local_civil_installation_allowance_USD2041','local_new_component_value_USD2041','additional_resource_provision_USD2041'])
    build_years = q['campaign_years']+delay
    cost = capital_schedule(initial_capital,build_years)
    for t in range(build_years-commission['years'],build_years):
        cost[t] += commercial/commission['years']
    cash = [-x for x in cost]
    resource_cost = cost[:]
    for t in range(build_years-commission['years'],build_years):
        resource_cost[t] -= commission['additional_resource_provision_USD2041']/commission['years']
    benefits = [0.]*build_years
    opportunity = [0.]*build_years
    cash_support = [0.]*build_years
    path = growth_path(q,economy,growth,horizon)
    initial_people = path[0]['population']
    grow = economy['high_value_capacity_with_growth']
    base = economy['high_value_capacity_case']
    n = grow['engineering']['growth_people_year']
    expansion_per_person = grow['external']['USD2041']['expansion_external_purchases']/n
    growth_provision_per_person = (grow['external']['USD2041']['resource_provision']-base['external']['USD2041']['resource_provision'])/n
    # Domestic capital formation is part of growth resource use, not Earth imports.
    local_expansion_per_person = (grow['domestic']['USD2041']['domestic_capital_formation']-base['domestic']['USD2041']['domestic_capital_formation'])/n
    operating = []
    for g in path:
        scale = g['scale']
        demand = scale if demand_cap is None else min(scale,demand_cap)
        sales = ext['gross_external_sales']*demand
        imports = ext['total_imports']*scale
        provision = ext['resource_provision']*scale
        expansion = (expansion_per_person+growth_provision_per_person)*g['starts']
        surplus = sales-imports-provision-expansion
        collected = capture*max(0,surplus)+min(0,surplus)
        # Global proxy: final service value (at purchase-price benchmarks) less
        # imported resource price proxies, local opportunity costs and growth.
        # Domestic GDP, supplier turnover and financing transfers are not added.
        resident_service = dom['resident_consumption']*scale*service_value
        labour = g['population']*q['resident_work_fraction']*q['alternative_USD2041_worker_year']*counterfactual_scale
        extra_local = q['unpriced_local_resource_USD2041_year']*scale
        costs = imports+expansion+q['local_growth_cost_fraction']*local_expansion_per_person*g['starts']+labour+extra_local
        # Reserve is a cash restriction, not itself ecological damage. It enters
        # expansion funding above; remove its growth part from resource costs.
        costs -= growth_provision_per_person*g['starts']
        benefit = sales*export_value+resident_service+spillover
        lost_exports = ext['gross_external_sales']*(scale-demand)
        cash.append(collected); resource_cost.append(costs); benefits.append(benefit)
        opportunity.append(labour); cash_support.append(max(0,-surplus))
        operating.append(dict(**g,exports_USD2041=sales,imports_USD2041=imports,
            resource_reserve_USD2041=provision,expansion_funding_USD2041=expansion,
            external_cash_USD2041=surplus,collectible_cash_USD2041=collected,
            final_service_value_USD2041=benefit,resource_cost_proxy_USD2041=costs,
            domestic_GVA_USD2041=dom['GDP_GVA']*scale-lost_exports+local_expansion_per_person*g['starts'],
            domestic_NDP_USD2041=dom['net_domestic_product']*scale-lost_exports+local_expansion_per_person*g['starts']))
    cash_npv = pv(cash,rate)
    social_net = [b-c for b,c in zip(benefits,resource_cost)]
    social_npv = pv(social_net,rate)
    recurring_factors = [0.]*build_years+[1.]*horizon
    annuity = pv(recurring_factors,rate)
    total=0.; peak=0.; balances=[]
    for v in cash:
        total+=v; peak=max(peak,-total); balances.append(total)
    # Durable recovery: no subsequent deficit within the evaluated window.
    recovery=next((t for t,v in enumerate(balances) if v>=0 and min(balances[t:])>=0),None)
    tail = operating[-1]
    return dict(portfolio=portfolio,horizon_operating_years=horizon,real_rate=rate,growth_rate=growth,
        demand_cap_multiple=demand_cap,capture_fraction=capture,service_value_fraction=service_value,
        counterfactual_earnings_scale=counterfactual_scale,export_counterfactual_value_fraction=export_value,extra_net_spillover_USD2041_year=spillover,
        campaign_years=build_years,initial_campaign_USD2041=initial_capital,commercial_commissioning_USD2041=commercial,
        cash_NPV_USD2041=cash_npv,undiscounted_cash_USD2041=sum(cash),peak_cumulative_funding_USD2041=peak,
        durable_recovery_year_from_start=recovery,required_extra_annual_cash_USD2041=max(0,-cash_npv)/annuity,
        proxy_economic_NPV_USD2041=social_npv,proxy_benefit_cost_ratio=pv(benefits,rate)/pv(resource_cost,rate),
        required_extra_annual_net_benefit_USD2041=max(0,-social_npv)/annuity,
        final_population=tail['population'],final_pipeline=tail['pipeline'],
        final_external_cash_USD2041=tail['external_cash_USD2041'],
        cumulative_source_feed_scale_years=sum(g['scale'] for g in path),
        cumulative_resident_years=sum(g['population'] for g in path),
        cash_flows_USD2041=cash,operating=operating)


def compact(row):
    return {k:v for k,v in row.items() if k not in ('cash_flows_USD2041','operating')}


def results():
    q,e,i=inputs()
    main=[]
    for name in ['specialty','qualified_high_value','large']:
        for horizon in q['horizons_operating_years']:
            for rate in q['real_rates']:
                main.append(compact(evaluate(q,e,name,horizon,rate)))
    growth={}
    for r in q['growth_rates']:
        for demand in [None,5]:
            name=f'{r}:{demand}'
            growth[name]=evaluate(q,e,'large',q['growth_horizon_years'],.01,growth=r,demand_cap=demand)
    growth_horizons=[compact(evaluate(q,e,'large',h,r,growth=g,demand_cap=cap))
                     for h in q['horizons_operating_years'] for r in q['real_rates']
                     for g in q['growth_rates'] for cap in [None,5]]
    valuation=[]
    for service in [0,0.5,1]:
        for counterfactual in [0.5,1,2]:
            valuation.append(compact(evaluate(q,e,'large',100,.03,service_value=service,counterfactual_scale=counterfactual)))
    stakeholder_values=[compact(evaluate(q,e,'large',100,.03,export_value=v,spillover=b))
                       for v in [0,.5,1] for b in [0,1e11,1e12]]
    collection=[compact(evaluate(q,e,'large',100,.03,capture=s)) for s in [0.25,.5,1]]
    sensitivities=[]
    for multiplier in [1,1.25,1.5]:
        for delay in [0,10]:
            x=deepcopy(q);x['campaign_scope_multiplier']=multiplier
            sensitivities.append(dict(scope_multiplier=multiplier,delay_years=delay,**compact(evaluate(x,e,'large',100,.03,delay=delay))))
    investor_base=evaluate(q,e,'large',100,.03)['cash_flows_USD2041']
    investors=[]
    for funding,claim in [(1,1),(.5,.5),(.25,.5)]:
        own=[min(0,x)*funding+max(0,x)*claim for x in investor_base]
        other=[x-y for x,y in zip(investor_base,own)]
        investors.append(dict(funding_fraction=funding,positive_cash_claim=claim,investor_NPV_USD2041=pv(own,.03),other_stakeholder_cash_NPV_USD2041=pv(other,.03),combined_NPV_USD2041=pv(investor_base,.03)))
    baseline=e['baseline_profiles']['local_structures']['external']['USD2041']['after_renewal_and_resource']
    high=e['high_value_capacity_case']['external']['USD2041']['after_renewal_and_resource']
    commission=e['high_value_capacity_commissioning']
    upgrade=sum(commission[k] for k in ['additional_external_purchases_USD2041','local_civil_installation_allowance_USD2041','local_new_component_value_USD2041','additional_resource_provision_USD2041'])
    incremental=[]
    for h in [10,30,100]:
        flow=[-upgrade/commission['years']]*commission['years']+[high-baseline]*h
        incremental.append(dict(operating_years=h,real_rate=.03,additional_funding_USD2041=upgrade,annual_external_balance_improvement_USD2041=high-baseline,incremental_cash_NPV_USD2041=pv(flow,.03),rule='Additional commercial programme against continuing the already funded modest-service system. Full additional programme cost charged; no whole-campaign capital recovery inferred.'))
    architecture={k:{'establishment_USD2041':v['cost']['capital_with_contingency_USD2041'],
                    'campaign_USD2041':v['cost']['industrial_and_initial_material_campaign_USD2041']}
                  for k,v in i['cases'].items() if k in ('reference','local_fleet_fabrication')}
    return dict(schema='terluna.research.complex-returns/1',evidence=q['evidence'],boundary=q['boundary'],inputs=q,
        producer={'path':'research/studies/complex_returns/run.py','sha256':digest(HERE/'run.py'),'inputs_sha256':digest(HERE/'returns_inputs.json')},
        historical_architectures=architecture,matched_initial_capital=e['initial_capital_context'],
        static_returns=main,growth_returns=growth,growth_horizons=growth_horizons,valuation_sensitivity=valuation,
        collection_sensitivity=collection,cost_and_delay_sensitivity=sensitivities,stakeholder_value_sensitivity=stakeholder_values,illustrative_capital_claims=investors,incremental_commercial_programme=incremental)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'research/runs/complex_returns/returns.json')
    args=p.parse_args();data=results();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');print(args.output)

if __name__=='__main__':main()
