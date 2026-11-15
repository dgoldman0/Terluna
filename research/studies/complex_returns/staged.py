"""Cash-constrained modular expansion; constant USD2041, no assumed bailouts.

Engineering small modules and service demand remain hypotheses. Pinned products
provide the matched account, not a qualified bootstrap deployment design.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
import math
from . import run

HERE = run.HERE


def inputs():
    q = json.loads((HERE / 'staged_inputs.json').read_text())
    _, e, _ = run.inputs()
    p = HERE / q['transport_snapshot']
    if run.digest(p) != q['transport_sha256']:
        raise ValueError('Unreviewed transport product')
    t = json.loads(p.read_text())
    if t['schema'] != 'terluna.research.industrial-transport/1':
        raise ValueError('Unexpected transport schema')
    return q, e, t


def capital_account(e, t):
    c = t['source_structure_localization']['space_chemical']['cost']
    contingency = c['capital_with_contingency_USD2041'] / c['capital_direct_USD2041']
    stock = sum(c[k] for k in ('initial_stock_manufacture_USD2041',
        'initial_stock_Earth_freight_USD2041', 'initial_stock_restoration_USD2041')) * contingency
    supplier = c['campaign_supplier_operation_and_renewal_USD2041'] * contingency
    makeup = c['construction_period_settlement_makeup_USD2041'] * contingency
    establishment = c['capital_with_contingency_USD2041']
    commercial = sum(e['high_value_capacity_commissioning'][k] for k in (
        'additional_external_purchases_USD2041', 'local_civil_installation_allowance_USD2041',
        'local_new_component_value_USD2041', 'additional_resource_provision_USD2041'))
    campaign = e['initial_capital_context']['inherited_industrial_and_material_campaign_USD2041']
    if not math.isclose(establishment + stock + supplier + makeup, campaign, rel_tol=1e-12):
        raise ValueError('Campaign reconciliation failed')
    b = e['high_value_capacity_case']; g = e['high_value_capacity_with_growth']
    n = g['engineering']['growth_people_year']
    people = b['engineering']['population']
    growth_cost = (g['external']['USD2041']['expansion_external_purchases'] +
        g['external']['USD2041']['resource_provision'] -
        b['external']['USD2041']['resource_provision']) / n * people
    return dict(establishment_USD2041=establishment, initial_stocks_USD2041=stock,
        original_supplier_operation_and_renewal_USD2041=supplier,
        original_construction_makeup_USD2041=makeup, original_campaign_USD2041=campaign,
        one_time_bill_USD2041=establishment + stock + commercial,
        commercial_USD2041=commercial, marginal_post_founding_USD2041_per_scale=growth_cost)


def kit_limit(q, scale, people):
    k = q['service_kit']; c = q['cultural_platform']
    staff = (people * scale * q['maximum_new_service_staff_fraction_of_residents'] - c['staff']) / k['staff']
    power = (q['reserved_service_power_W_per_founding_scale'] * scale - c['power_W']) / k['power_W']
    return max(0, min(k['max_units'], math.floor(min(staff, power) + 1e-9)))


def service_account(q, scenario, age, kits, previous_paid=0, previous_licenses=0):
    s = q['scenarios'][scenario]; k = q['service_kit']; c = q['cultural_platform']
    maturity = min(1., q['service_demand_initial_fraction'] *
                   (1 + q['service_demand_annual_growth']) ** min(age, 100))
    audience = s['audience'] * maturity
    paid = audience * s['paid_fraction']; licenses = s['licenses'] * maturity
    contracts = min(s['research_contracts'] * maturity, kits * k['research_contracts_year'])
    clients = min(s['professional_clients'] * maturity, kits * k['professional_clients_year'])
    revenue = dict(culture=paid * s['paid_USD2041_year'] + audience * s['advertising_and_rights_USD2041_fan_year'],
        licenses=licenses * s['license_USD2041_year'], research=contracts * s['research_USD2041_contract'],
        professional=clients * s['professional_USD2041_client_year'])
    variable = sum(revenue[n] * q['variable_cost_fractions'][n] for n in revenue)
    acquired_paid = max(0., paid - previous_paid * (1-q['paid_churn']))
    acquired_licenses = max(0., licenses - previous_licenses * (1-q['license_renewal_churn']))
    acquisition = acquired_paid * q['paid_acquisition_USD2041'] + acquired_licenses * q['license_acquisition_USD2041']
    staff = c['staff'] + kits*k['staff']
    wages = staff*q['service_staff_USD2041_year']
    fixed = c['fixed_purchases_USD2041_year'] + kits*k['fixed_purchases_USD2041_year']
    renewal = c['capital_USD2041']/c['life_years'] + kits*k['capital_USD2041']/k['life_years']
    costs = variable + acquisition + wages + fixed + renewal
    return dict(audience=audience, paid_accounts=paid, licenses=licenses,
        research_contracts=contracts, professional_clients=clients, staff=staff,
        power_W=c['power_W']+kits*k['power_W'], revenue_by_activity_USD2041=revenue,
        revenue_USD2041=sum(revenue.values()), variable_purchases_USD2041=variable,
        acquisition_USD2041=acquisition, wages_USD2041=wages, fixed_purchases_USD2041=fixed,
        renewal_USD2041=renewal, total_cost_USD2041=costs, available_USD2041=sum(revenue.values())-costs)


def core_account(q, e, installed, productive):
    a = e['high_value_capacity_case']['external']; v = a['USD2041']
    replaced = q['replaced_baseline_Earth_service_fractions']
    sales = a['sales_by_activity_USD2041']
    remaining = v['gross_external_sales'] - sum(sales[k] for k in replaced)
    sold_scale = min(productive, q['material_compute_demand_cap_multiple'])
    # Keep all physical service assets and logistics. Remove only identified
    # overlapping variable Earth service purchases at the installed scale.
    avoided = sum(sales[k]*fraction for k,fraction in replaced.items())
    # The full common floor operates even when only a small first module opens.
    common = q['common_recurring_fraction']*max(0.,1-installed) if installed else 0.
    imports = v['total_imports'] * (common + installed) - avoided*installed
    reserve = v['resource_provision']*installed
    return dict(revenue_USD2041=remaining*sold_scale, imports_USD2041=imports,
        resource_provision_USD2041=reserve, available_USD2041=remaining*sold_scale-imports-reserve)


def expansion_cost(q, capital, committed, increment):
    common = q['common_capital_fraction_of_establishment']*capital['establishment_USD2041']
    founding_part = min(increment, max(0., 1-committed))
    return q['marginal_capital_multiplier'] * (founding_part*(capital['one_time_bill_USD2041']-common) +
        (increment-founding_part)*capital['marginal_post_founding_USD2041_per_scale'])


def affordable_increment(q, capital, committed, limit, cash):
    lo, hi = 0., limit
    for _ in range(48):
        mid = (lo+hi)/2
        if expansion_cost(q,capital,committed,mid) <= cash: lo = mid
        else: hi = mid
    return lo


def simulate(q, e, t, scenario, horizon, *, grow=True):
    if horizon <= q['bootstrap_years'] or not 0 < q['seed_scale'] <= 1:
        raise ValueError('Horizon must include operations and seed must be in (0,1]')
    for key in ('collectible_fraction','retained_fraction','common_recurring_fraction',
                'common_capital_fraction_of_establishment'):
        if not 0 <= q[key] <= 1: raise ValueError('Invalid fraction: '+key)
    capital = capital_account(e,t)
    seed = q['seed_scale']; start = q['bootstrap_years']; k = q['service_kit']; c = q['cultural_platform']
    people = e['high_value_capacity_case']['engineering']['population']
    if seed*people*q['maximum_new_service_staff_fraction_of_residents'] < c['staff'] or seed*q['reserved_service_power_W_per_founding_scale'] < c['power_W']:
        raise ValueError('Seed cannot host even the cultural platform')
    kits = min(1,kit_limit(q,seed,people))
    common = capital['establishment_USD2041']*q['common_capital_fraction_of_establishment']
    hardware = common + seed*(capital['one_time_bill_USD2041']-common) + c['capital_USD2041']+kits*k['capital_USD2041']
    def fixed_burn(scale, units):
        a=core_account(q,e,scale,0)
        return a['imports_USD2041']+a['resource_provision_USD2041'] + c['fixed_purchases_USD2041_year'] + c['capital_USD2041']/c['life_years'] + units*(k['fixed_purchases_USD2041_year']+k['capital_USD2041']/k['life_years']) + (c['staff']+units*k['staff'])*q['service_staff_USD2041_year']
    seed_cash = fixed_burn(seed,kits)*q['seed_no_receipt_runway_years']
    bootstrap_carry = fixed_burn(seed,kits)*q['bootstrap_carrying_fraction']*start
    investor_flows = [-(hardware+bootstrap_carry)/start]*start
    investor_flows[-1] -= seed_cash
    cash = seed_cash; earned = 0.; cohorts = [(start, seed)]; pending = []; kit_pending = []
    rows=[]; paid=licenses=0.; failure=None; first_full=None; reached_cap=None
    cap=e['growth']['0.01']['vacuum_architectural_reference_with_existing_cities_people']/people
    frontier=e['high_value_capacity_with_growth']['engineering']['growth_people_year']/people
    for year in range(start,horizon+1):
        cash_open=cash
        for due,amount in pending:
            if due==year: cohorts.append((due,amount))
        pending=[x for x in pending if x[0]>year]
        kits += sum(n for due,n in kit_pending if due==year)
        kit_pending=[x for x in kit_pending if x[0]>year]
        installed=sum(n for _,n in cohorts)
        productive=sum(n*min(1.,(year-due+1)/q['production_ramp_years']) for due,n in cohorts)
        core=core_account(q,e,installed,productive)
        service=service_account(q,scenario,year-start,kits,paid,licenses)
        paid=service['paid_accounts']; licenses=service['licenses']
        available=core['available_USD2041']+service['available_USD2041']
        pooled=min(0.,available)+max(0.,available)*q['collectible_fraction']
        other=max(0.,available)*(1-q['collectible_fraction'])
        cash += pooled
        earned += pooled
        # Protect the unspent startup runway while production and audiences ramp.
        buffer=fixed_burn(installed,kits)*max(q['liquidity_buffer_years'],
            q['seed_no_receipt_runway_years']-(year-start))
        orders=kit_orders=dividend=0.; starts=0.; kit_saving=0.
        if cash < 0:
            failure=dict(year=year, shortfall_USD2041=-cash)
        else:
            room=max(0.,min(cash-buffer,earned))
            dividend=min(room,max(0.,pooled)*(1-q['retained_fraction']))
            cash-=dividend;earned-=dividend
            committed=installed+sum(n for _,n in pending)
            maturity=service['audience']/q['scenarios'][scenario]['audience']
            desired_kits=min(k['max_units'], math.ceil(maturity*max(
                q['scenarios'][scenario]['research_contracts']/k['research_contracts_year'],
                q['scenarios'][scenario]['professional_clients']/k['professional_clients_year'])))
            existing_kits=kits+sum(n for _,n in kit_pending)
            needs_kit=grow and existing_kits < min(kit_limit(q,installed,people),desired_kits)
            kit_saving=k['capital_USD2041'] if needs_kit else 0.
            # Kits must fit currently commissioned power/staff, not future promises.
            if needs_kit and min(cash-buffer,earned) >= k['capital_USD2041']:
                kit_orders=k['capital_USD2041'];cash-=kit_orders;earned-=kit_orders
                kit_pending.append((year+q['commissioning_years'],1))
                kit_saving=0.
            if grow:
                slope=q['desired_growth_rate']*installed
                if committed > q['taper_fraction']*cap:
                    slope=q['desired_growth_rate']*q['taper_fraction']*cap*math.sqrt(max(0.,(cap-committed)/((1-q['taper_fraction'])*cap)))
                limit=max(0.,min(slope,frontier*installed,cap-committed,
                    q['material_compute_demand_cap_multiple']-committed))
                # Demand alone does not make a lossmaking increment investable.
                eps=1e-5
                delta=core_account(q,e,committed+eps,committed+eps)['available_USD2041']-core_account(q,e,committed,committed)['available_USD2041']
                if delta <= 0: limit=0.
                starts=affordable_increment(q,capital,committed,limit,max(0.,min(cash-buffer,earned)-kit_saving))
                if starts>1e-12:
                    orders=expansion_cost(q,capital,committed,starts);cash-=orders;earned-=orders
                    pending.append((year+q['commissioning_years'],starts))
                else: starts=0.
            # Retained earnings with no funded use are released after the buffer.
            excess=max(0.,min(cash-buffer,earned)-kit_saving)
            dividend+=excess;cash-=excess;earned-=excess
        investor_flows.append(dividend)
        row=dict(year=year,installed_scale=installed,productive_scale=productive,
            supported_residents=installed*people,service_kits=kits,core=core,services=service,
            available_USD2041=available,pooled_USD2041=pooled,other_stakeholder_cash_USD2041=other,
            opening_cash_USD2041=cash_open,closing_cash_USD2041=cash,
            retained_earned_cash_USD2041=earned,
            industrial_orders_USD2041=orders,service_orders_USD2041=kit_orders,
            starts_scale=starts,pipeline_scale=sum(n for _,n in pending),
            pipeline_kits=sum(n for _,n in kit_pending),dividend_USD2041=dividend,
            planned_service_saving_USD2041=min(max(0.,earned),kit_saving),
            liquidity_buffer_USD2041=buffer)
        rows.append(row)
        if installed>=1-1e-9 and first_full is None: first_full=year
        if installed>=cap-1e-8 and reached_cap is None: reached_cap=year
        if failure: break
    # No asset sale or release of restricted resource provision. Cash can be
    # returned at a finite appraisal endpoint; operating obligations continue.
    terminal_liquid=max(0.,cash-buffer) if not failure else 0.
    investor_flows[-1]+=terminal_liquid
    totals=[];balance=0.
    for f in investor_flows: balance+=f;totals.append(balance)
    recovery=next((i for i,b in enumerate(totals) if b>=0 and min(totals[i:])>=0),None)
    summary=dict(scenario=scenario,horizon_elapsed_years=horizon,expansion_policy=grow,
        seed_hardware_USD2041=hardware,seed_liquidity_USD2041=seed_cash,
        bootstrap_carry_USD2041=bootstrap_carry,
        external_seed_USD2041=hardware+seed_cash+bootstrap_carry,external_funding_after_first_production_USD2041=0,
        failure=failure,first_founding_scale_year=first_full,cap_year=reached_cap,
        final_scale=rows[-1]['installed_scale'],final_supported_residents=rows[-1]['supported_residents'],
        final_service_revenue_USD2041=rows[-1]['services']['revenue_USD2041'],
        final_available_USD2041=rows[-1]['available_USD2041'],
        total_industrial_orders_USD2041=sum(x['industrial_orders_USD2041'] for x in rows),
        total_service_orders_USD2041=sum(x['service_orders_USD2041'] for x in rows),
        total_dividends_USD2041=sum(x['dividend_USD2041'] for x in rows),
        terminal_liquid_USD2041=terminal_liquid,
        retained_operating_liquidity_USD2041=min(max(0.,cash),buffer),durable_cash_recovery_year=recovery,
        NPV_USD2041={str(r):run.pv(investor_flows,r) for r in q['discount_rates']},
        additional_customer_surplus_PV_at_3pct_USD2041={str(v):sum(x['services']['audience']*v/(1.03)**x['year'] for x in rows) for v in (0,10,50)})
    return dict(summary=summary,annual=rows,capital_provider_flows_USD2041=investor_flows)


def results():
    q,e,t=inputs(); paths={name:simulate(q,e,t,name,max(q['horizons_elapsed_years'])) for name in q['scenarios']}
    horizons=[simulate(q,e,t,name,h)['summary'] for name in q['scenarios'] for h in q['horizons_elapsed_years']]
    variants=[('dividends_without_expansion',{},False), ('half_service_prices',{},True),
        ('no_additional_B2B_demand',{},True),('higher_variable_costs',{},True),
        ('lower_collection',{'collectible_fraction':.7},True),
        ('twice_marginal_capital',{'marginal_capital_multiplier':2},True),
        ('larger_common_floor',{'common_capital_fraction_of_establishment':.2,'common_recurring_fraction':.1},True),
        ('five_times_industrial_market',{'material_compute_demand_cap_multiple':5},True),
        ('one_percent_seed',{'seed_scale':.01},True),('ten_percent_seed',{'seed_scale':.1},True),
        ('slow_service_adoption',{'service_demand_annual_growth':.15},True),
        ('smaller_launch_audience',{'service_demand_initial_fraction':.05},True)]
    sensitivity=[]
    for name,overrides,grow in variants:
        v=deepcopy(q);v.update(overrides)
        if name=='half_service_prices':
            for case in v['scenarios'].values():
                for field in case:
                    if 'USD2041' in field: case[field]*=.5
        if name=='no_additional_B2B_demand':
            for case in v['scenarios'].values():
                case['research_contracts']=case['professional_clients']=0
        if name=='higher_variable_costs':
            v['variable_cost_fractions']={key:.8 for key in v['variable_cost_fractions']}
        for scenario in ('global','breakthrough'):
            x=simulate(v,e,t,scenario,800,grow=grow)['summary']
            sensitivity.append(dict(variant=name,**x))
    comparisons=[]
    capital=capital_account(e,t)
    people=e['high_value_capacity_case']['engineering']['population']
    for scenario in ('global','breakthrough'):
        s=q['scenarios'][scenario];k=q['service_kit'];c=q['cultural_platform']
        kits=min(kit_limit(q,1,people),math.ceil(max(s['research_contracts']/k['research_contracts_year'],s['professional_clients']/k['professional_clients_year'])))
        service_capital=c['capital_USD2041']+kits*k['capital_USD2041']
        for horizon in (100,300):
            flows=[-capital['original_campaign_USD2041']/33]*33
            for year in range(28,33):flows[year]-=(capital['commercial_USD2041']+service_capital)/5
            paid=licenses=0.
            for year in range(33,horizon+1):
                service=service_account(q,scenario,year-33,kits,paid,licenses)
                paid=service['paid_accounts'];licenses=service['licenses']
                flows.append(core_account(q,e,1,1)['available_USD2041']+service['available_USD2041'])
            comparisons.append(dict(scenario=scenario,horizon_elapsed_years=horizon,
                fixed_campaign_NPV_at_3pct_USD2041=run.pv(flows,.03),
                staged_NPV_at_3pct_USD2041=next(x['NPV_USD2041']['0.03'] for x in horizons if x['scenario']==scenario and x['horizon_elapsed_years']==horizon),
                fixed_supported_residents_after_year33=people,
                staged_supported_residents=next(x['final_supported_residents'] for x in horizons if x['scenario']==scenario and x['horizon_elapsed_years']==horizon)))
    return dict(schema='terluna.research.staged-returns/1',evidence=q['evidence'],
        input_sha256=run.digest(HERE/'staged_inputs.json'),producer_sha256=run.digest(HERE/'staged.py'),
        dependency_sha256={name:run.digest(HERE/name) for name in ('run.py','returns_inputs.json')},
        capital_account=capital,paths=paths,horizons=horizons,sensitivities=sensitivity,
        fixed_commitment_comparisons=comparisons)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=__import__('pathlib').Path,
        default=run.ROOT/'research/runs/complex_returns/staged.json')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(args.output)
