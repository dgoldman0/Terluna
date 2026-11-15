"""Reconciled domestic, external and growth accounts for the founding system.

Physical capacity is distinct from demand. Future prices and modular designs
are scenarios; no foreign-exchange peg or source entitlement is inferred.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
from pathlib import Path
from . import evolution as model
from . import run as base

HERE, ROOT = base.HERE, base.ROOT
SCHEMA = 'terluna.research.initial-system-economy/1'


def compute_price(q, electricity=.1, markup=1):
    """Competitive scenario using identical IT purchase and replacement terms.

Quoted per installed W-year, including billable duty, with a full annual IT
replacement bill. Network performance and customer acceptance remain gates.
"""
    e=q['exports']
    return markup*((e['compute_IT_USD2026_W']+e['compute_communications_USD2026_W'])/e['compute_life_years']+e['compute_Earth_overhead_USD2026_W_year']+electricity*base.HOURS*e['compute_duty']/1000)/(1-e['compute_common_Earth_service_share'])


def external(q,src,s,c,a,goods_price=750,source_price=30,compute_markup=1,inventory_price=10):
    args={k:a[k] for k in ['profile','population','plant_scale','service_portfolio']}
    args.update(goods=a['finished_export_kg_year'],source_goods=a['source_export_kg_year'],compute=a['compute_W'],propellant_sales=a['propellant_sales_kg_year'],extra_metal=a['stock_metal_kg_year'],refinery_capacity=a['refinery_capacity_kg_year'])
    steady=model.annual(q,src,s,c,**args)
    op=model.annual(q,src,s,c,**args,renewal_factor=0)
    gross={}
    survival=(1-src['revenue_inputs']['return']['loss_probability_flight'])*src['revenue_inputs']['manufacturing']['delivered_quality_yield']
    gross['lunar_goods_delivered']=args['goods']*survival*goods_price
    gross['source_goods_delivered']=args['source_goods']*survival*source_price
    gross['batch_computing']=args['compute']*compute_price(q,markup=compute_markup)
    gross['independent_propellant_customers']=args['propellant_sales']*src['revenue_inputs']['markets']['reference']['propellant_USD2026_kg']
    # All off-Earth productive operations share the boundary; co-products
    # retained there are domestic inventory/capital, even if Earth-financed.
    if a['service']:
        for key,value in a['service']['receipts_USD2026'].items():
            if key not in ('manufactured_goods_at_Earth','batch_computing','Earth_mission_propellant'):gross[key]=value
    # Earth communications/contract services retain the earlier 20% computing
    # sales allowance; imported IT and freight are additional physical costs.
    compute_services=gross['batch_computing']*src['revenue_inputs']['computing']['Earth_service_share_of_receipts']
    total=sum(gross.values());operation=op['cash_USD2026']['total_external_purchases']+compute_services
    renewal=steady['cash_USD2026']['total_external_purchases']-op['cash_USD2026']['total_external_purchases']
    expansion=a['cash_USD2026']['total_external_purchases']-steady['cash_USD2026']['total_external_purchases']
    reserve=a['cash_USD2026']['resource_provision'];I=c['cost']['inflation_factor']
    vals={'gross_external_sales':total,'operating_imports':operation,'operating_cash_surplus':total-operation,'external_renewal':renewal,'resource_provision':reserve,'after_renewal_and_resource':total-operation-renewal-reserve,'expansion_external_purchases':expansion,'after_expansion_and_resource':total-operation-renewal-expansion-reserve,'total_imports':operation+renewal+expansion}
    return {'USD2041':{k:v*I for k,v in vals.items()},'USD2026':vals,'sales_by_activity_USD2041':{k:v*I for k,v in gross.items()},'debt_service_USD2041':None,'price_scenarios':{'lunar_goods_USD2026_kg':goods_price,'source_goods_USD2026_kg':source_price,'compute_USD2026_W_year':compute_price(q,markup=compute_markup),'source_inventory_USD2026_kg':inventory_price},'definition':'Consolidated trade of all off-Earth operations. Internal sales cancel here. Net measures exclude debt, taxes and dividends; the resource provision is retained funding, not payment for resource rights.'}


def supplies(q,a):
    """Explicit transfer-price convention for the domestic supply-use example."""
    p=q['production'];v=q['valuation'];f=a['physical']
    return (f['lunar_bulk_kg_year']*p['local_bulk_value_USD2026_kg']+(f['lunar_metal_kg_year']+f['source_metal_kg_year'])*p['structural_metal_value_USD2026_kg']+f['source_water_kg_year']*p['source_water_value_USD2026_kg']+f['source_nitrogen_kg_year']*p['source_nitrogen_value_USD2026_kg']+(f['lunar_process_W']+f['resident_mean_W']+f['source_process_W']+f['depot_W'])*base.HOURS/1000*v['power_USD2026_kWh']+a['cash_USD2026']['total_external_purchases']*v['domestic_transport_service_fraction_of_external_cost'])


def domestic(q,src,s,c,a,trade):
    """Four-sector supply-use account with imported capital recorded separately.

All intermediate imports are consolidated into industrial supplies. This
preserves the GDP identity without pretending to possess a surveyed I/O table.
Food, household production, housing and civic services are one resident sector.
"""
    args={k:a[k] for k in ['profile','population','plant_scale']}
    basic=model.annual(q,src,s,c,**args)
    consumption=model.annual(q,src,s,c,**args,renewal_factor=0)
    steady=model.annual(q,src,s,c,**args,service_portfolio=a['service_portfolio'],goods=a['finished_export_kg_year'],source_goods=a['source_export_kg_year'],compute=a['compute_W'],propellant_sales=a['propellant_sales_kg_year'],extra_metal=a['stock_metal_kg_year'],refinery_capacity=a['refinery_capacity_kg_year'])
    # Consumer supplies, replacement inputs, trade inputs and expansion inputs
    # sum to the full year's industrial output, with every transfer once.
    Sc=supplies(q,consumption);Si=supplies(q,basic)-Sc+supplies(q,a)-supplies(q,steady);Sx=supplies(q,steady)-supplies(q,basic)
    stock_inputs=0.
    if a['stock_metal_kg_year']:
        no_stock=model.annual(q,src,s,c,**args,service_portfolio=a['service_portfolio'],goods=a['finished_export_kg_year'],source_goods=a['source_export_kg_year'],compute=a['compute_W'],propellant_sales=a['propellant_sales_kg_year'],refinery_capacity=a['refinery_capacity_kg_year'])
        stock_inputs=supplies(q,steady)-supplies(q,no_stock)
        Si+=stock_inputs;Sx-=stock_inputs
    inventory_value=a['stock_metal_kg_year']*trade['price_scenarios']['source_inventory_USD2026_kg']
    I=c['cost']['inflation_factor'];v=q['valuation'];ps=a['population']/c['demand']['population'];m=src['settlement']['cases']['reference']['maintenance'];op=m['variable_operations_U_year'];normal=v['inherited_U_service_normalization_USD2041_per_U']/I
    civic=(m['routine_upkeep_U_year']*v['routine_service_value_added_fraction']+op['production_operators']+op['basic_civic_services']+op['everyday_process_consumables_and_services']*v['routine_service_value_added_fraction'])*normal*ps
    # These are service/assembly valuations, not a conversion of the old budget.
    # Old air, water, power, feed, freight, renewal and contingency are excluded.
    parts=a['components'];cf=a['capital_formation'];p=q['production'];t=trade['USD2026']
    local_kg=sum(a['physical']['local_fabrication_kg_year'].values())
    finishing=local_kg*max(0,p['local_fabrication_value_USD2026_kg']-p['structural_metal_value_USD2026_kg'])
    commercial_growth=sum(r['import_growth_USD2026']+r['local_growth_value_USD2026'] for r in parts if r['owner']=='commerce')
    installation=(cf['imported_renewal_USD2026']+cf['imported_growth_USD2026']+cf['local_renewal_value_USD2026']+cf['local_growth_value_USD2026']-commercial_growth)*v['stock_installation_fraction_of_hardware_value']
    construction=a['growth_people_year']*v['civic_service_growth_USD2026_per_new_resident']+a['cash_USD2026']['commercial_civil_installation_allowance']
    local_I=Si-stock_inputs+finishing+installation+construction+inventory_value
    civil_CFC=(c['demand']['renewal_bulk']*ps*p['local_bulk_value_USD2026_kg']+c['demand']['renewal_metal']*ps*p['structural_metal_value_USD2026_kg'])
    city_CFC=sum(x['annual_CFC_USD2026'] for x in parts if x['owner']=='settlement')+civil_CFC
    # Residence and public-service prices include consumption of their capital.
    C=Sc+civic+city_CFC;X=t['gross_external_sales'];Mi=t['operating_imports'];Mk=t['external_renewal']+t['expansion_external_purchases']
    sectors=[{'sector':'Industrial materials, utilities and transport','gross_output':Sc+Si+Sx,'domestic_intermediates':0,'imported_intermediates':Mi},
        {'sector':'Food, ordinary goods, housing and civic services','gross_output':C,'domestic_intermediates':Sc,'imported_intermediates':0},
        {'sector':'Domestic replacement and expansion work','gross_output':local_I,'domestic_intermediates':Si,'imported_intermediates':0},
        {'sector':'Goods and services sold outside the system','gross_output':X,'domestic_intermediates':Sx,'imported_intermediates':0}]
    for row in sectors:row['value_added']=row['gross_output']-row['domestic_intermediates']-row['imported_intermediates']
    GDP=sum(row['value_added'] for row in sectors);investment=local_I+Mk;imports=Mi+Mk;expenditure=C+investment+X-imports
    # CFC uses current component replacement valuations plus delivered renewal
    # costs and civil renewal. Resource provision is separate from depreciation.
    delivery=max(0,t['external_renewal']-cf['imported_renewal_USD2026'])
    CFC=cf['component_CFC_USD2026']+civil_CFC+delivery
    total={'gross_domestic_output':sum(row['gross_output'] for row in sectors),'GDP_GVA':GDP,'resident_consumption':C,'gross_capital_formation':investment,'gross_fixed_capital_formation':investment-inventory_value,'inventory_change':inventory_value,'domestic_capital_formation':local_I,'capitalized_external_purchases':Mk,'exports':X,'imports':imports,'consumption_of_fixed_capital':CFC,'net_domestic_product':GDP-CFC,'GDP_from_final_expenditure':expenditure,'accounting_residual':GDP-expenditure,'resource_provision_separate':t['resource_provision']}
    return {'USD2041':{k:v*I for k,v in total.items()},'sectors_USD2041':[{k:v*I if isinstance(v,(int,float)) else v for k,v in row.items()} for row in sectors],
      'valuation_basis':{'civic_service_value_added_U_year_before_normalization':civic/normal,'USD2041_per_U_convention':v['inherited_U_service_normalization_USD2041_per_U'],'routine_service_share':v['routine_service_value_added_fraction']},
      'limits':'Illustrative cost and transfer-price account, not observed GDP or an independent currency calibration. Imported completed capital is absent from domestic gross output and enters investment and imports once. Sector grouping is consolidated, not a measured input-output survey; loss-making production can have negative value added. Changing assumed service prices changes the valuation, not physical output or Earth cash.'}


def brief(a):
    return {k:v for k,v in a.items() if k not in ('components','service')}


def evaluate(q,src,s,c,**args):
    prices={k:args.pop(k) for k in ['goods_price','source_price','compute_markup','inventory_price'] if k in args}
    a=model.annual(q,src,s,c,**args);trade=external(q,src,s,c,a,**prices)
    return {'engineering':brief(a),'external':trade,'domestic':domestic(q,src,s,c,a,trade)}


def commissioning(q,src,s,c,**args):
    """Full operations plus staged additions: conservative overlap requirement."""
    n=q['exports']['commercial_commission_years'];steady=model.annual(q,src,s,c,**args)
    trial=model.annual(q,src,s,c,**args,commission_fraction=1/n,qualification=args.get('profile')=='local_structures')
    I=c['cost']['inflation_factor']
    return {'years':n,'feasible_overlap':trial['feasible'],'binding_capacity':trial['binding_capacity'],'capacity_use':trial['capacity_use'],'additional_external_purchases_USD2041':(trial['cash_USD2026']['total_external_purchases']-steady['cash_USD2026']['total_external_purchases'])*n*I,'local_civil_installation_allowance_USD2041':trial['cash_USD2026']['commercial_civil_installation_allowance']*n*I,'local_new_component_value_USD2041':(trial['capital_formation']['local_growth_value_USD2026']-steady['capital_formation']['local_growth_value_USD2026'])*n*I,'additional_resource_provision_USD2041':(trial['cash_USD2026']['resource_provision']-steady['cash_USD2026']['resource_provision'])*n*I,'additional_source_feed_kg':(trial['physical']['source_feed_kg_year']-steady['physical']['source_feed_kg_year'])*n,'status':'Installed capacity limits still apply. This pays for extra hardware, qualified-design trials, air, civil stock and freight over five years alongside full operations. It does not establish development success or fund continuing base deficits.'}


def growth_path(q,src,s,c,rate=None,years=None):
    P=c['demand']['population'];K=1.;pending=[];rows=[];total_feed=0.;total_imports=0.;I=c['cost']['inflation_factor']
    crater=src['crater_programme'];factor=q['growth']['eventual_estate_factor']
    # This population reference bounds the exercise; it does not approve more
    # tube sites or turn a crater site quota into a demographic instruction.
    vacuum_cap=crater['initial_clustered']['architectural_population']['representative']+src['settlement']['cases']['reference']['city_population']
    eventual=crater['eventual_by_estate_factor'][factor]['architectural_population']['representative']
    delay=q['growth']['commissioning_duration_years']
    for year in range(1,(years or q['growth']['years'])+1):
        due=[x for x in pending if x[0]==year]
        for _,beds,plant in due:P+=beds;K+=plant
        pending=[x for x in pending if x[0]>year]
        frontier=model.frontier(q,src,s,c,'growth_people',population=P,plant_scale=K,balanced_growth=True)
        available=frontier['maximum'];room=max(0,vacuum_cap-P-sum(x[1] for x in pending))
        requested=available if rate is None else rate*P
        beds=min(available,requested,room);newK=beds/c['demand']['population']
        a=model.annual(q,src,s,c,profile='local_structures',population=P,plant_scale=K,growth_people=beds,new_plant_scale=newK)
        total_feed+=a['physical']['source_feed_kg_year'];total_imports+=a['cash_USD2026']['total_external_purchases']*I
        pending.append((year+delay,beds,newK))
        rows.append({'year_after_full_occupation':year,'occupied_capacity_people':P,'new_capacity_started_people':beds,'installed_industry_multiple':K,'new_industry_started_multiple':newK,'capacity_supported_maximum_people':available,'requested_capacity_people':requested,'clipped_by_capacity':requested>available+1e-5,'binding_capacity':a['binding_capacity'],'feasible':a['feasible'],'source_feed_kg':a['physical']['source_feed_kg_year'],'external_purchases_USD2041':a['cash_USD2026']['total_external_purchases']*I})
    return {'scenario':'capacity_limited' if rate is None else 'comparison_'+str(rate),'years':rows,'occupied_capacity_at_end_people':P,'pipeline_people':sum(x[1] for x in pending),'cumulative_source_feed_kg':total_feed,'cumulative_external_purchases_USD2041':total_imports,'vacuum_architectural_reference_with_existing_cities_people':vacuum_cap,'eventual_crater_architectural_reference_people':eventual,'rule':'Production-supported construction capacity, conditional on funding, additional approved sites, staffing and passenger access. Starts consume current capacity and enter service after two years. Founding material intensities are a proxy; the four-city/100-town ratio is not an expansion plan. No terraforming date or uncounted tube quota is inferred.'}


def run(spec=None):
    q,src,s,c=model.read_inputs(spec);result={'schema':SCHEMA,'date':q['date'],'evidence':q['evidence'],'inputs':q['inputs'],'producer':{'path':'research/studies/initial_industrial_system/economy.py','sha256':base.sha(__file__),'sources':{'path':'research/studies/initial_industrial_system/economy_sources.json','sha256':base.sha(HERE/'economy_sources.json')},'configuration':{'path':'research/studies/initial_industrial_system/economy_inputs.json','sha256':base.sha(HERE/'economy_inputs.json')},'model_dependencies':[{'path':'research/studies/initial_industrial_system/'+n,'sha256':base.sha(HERE/n)} for n in ['evolution.py','revenue.py','run.py','transport.py']]},'baseline_profiles':{},'frontiers':{},'portfolios':{},'rules':q['rules']}
    for profile in q['profiles']:
        result['baseline_profiles'][profile]=evaluate(q,src,s,c,profile=profile,service_portfolio=True)
    for variable in ['goods','source_goods','compute','growth_people','propellant_sales','extra_metal']:
        f=model.frontier(q,src,s,c,variable,service_portfolio=True)
        result['frontiers'][variable]={**{k:v for k,v in f.items() if k!='at_limit'},'at_limit':brief(f['at_limit'])}
    result['balanced_growth_frontier']={}
    for profile in q['profiles']:
        f=model.frontier(q,src,s,c,'growth_people',profile=profile,balanced_growth=True)
        result['balanced_growth_frontier'][profile]={'new_capacity_people_year':f['maximum'],'fraction_of_initial_population_year':f['maximum']/c['demand']['population'],'binding_capacity':f['at_limit']['binding_capacity'],'external_purchases_USD2041':f['at_limit']['cash_USD2026']['total_external_purchases']*c['cost']['inflation_factor']}
    # A shared, explicit portfolio lets price cases vary without granting extra
    # infrastructure, fuel or spare capacity to each activity independently.
    focus=dict(profile='local_structures',service_portfolio=True,goods=100e6,source_goods=50e6,compute=1e9,propellant_sales=25e6)
    for name,price,source_price,markup in [('commodity',3,3,1),('specialty',750,30,1),('qualified_high_value',2000,750,1.5)]:
        result['portfolios'][name]=evaluate(q,src,s,c,**focus,goods_price=price,source_price=source_price,compute_markup=markup)
    result['portfolio_commissioning']=commissioning(q,src,s,c,**focus)
    growth=model.frontier(q,src,s,c,'growth_people',balanced_growth=True,**focus)
    result['portfolio_remaining_balanced_growth']={'new_capacity_people_year':growth['maximum'],'binding_capacity':growth['at_limit']['binding_capacity'],'engineering':brief(growth['at_limit'])}
    result['portfolio_with_growth']=evaluate(q,src,s,c,**focus,growth_people=growth['maximum']*.5,new_plant_scale=growth['maximum']*.5/c['demand']['population'],goods_price=2000,source_price=750,compute_markup=1.5)
    no_trade=model.annual(q,src,s,c,profile='local_structures')
    result['coproducts']=model.coproducts(q,s,c,no_trade)
    result['coproduct_inventory_comparison']=evaluate(q,src,s,c,profile='local_structures',extra_metal=result['coproducts']['current_extra_metal_refining_capacity_kg_year']*.99,inventory_price=10)
    result['coproduct_inventory_comparison']['external']['demand_gate']='Refined inventory stays inside the consolidated off-Earth system. Its illustrative value enters domestic capital formation, with no external revenue. Future internal buyers must have a use and funding; no recurring unlimited stock accumulation is adopted.'
    result['coproduct_stock_incremental_external_cost_USD2041']=external(q,src,s,c,no_trade)['USD2041']['after_renewal_and_resource']-result['coproduct_inventory_comparison']['external']['USD2041']['after_renewal_and_resource']
    result['expanded_source_refinery']=evaluate(q,src,s,c,profile='local_structures',refinery_capacity=q['coproducts']['metal_refinery_expansion_kg_year'],extra_metal=q['coproducts']['metal_refinery_expansion_kg_year'],inventory_price=10)
    result['expanded_source_refinery_commissioning']=commissioning(q,src,s,c,profile='local_structures',refinery_capacity=q['coproducts']['metal_refinery_expansion_kg_year'],extra_metal=q['coproducts']['metal_refinery_expansion_kg_year'])
    result['expanded_source_refinery_incremental_external_cost_USD2041']=external(q,src,s,c,no_trade)['USD2041']['after_renewal_and_resource']-result['expanded_source_refinery']['external']['USD2041']['after_renewal_and_resource']
    feed=no_trade['physical']['source_feed_kg_year'];cp=q['coproducts'];city=src['settlement']['cases']['reference'];nutrient=city['food']['nutrient_element_makeup_kg_year']
    carbon=nutrient['carbon']+city['resources']['operating_import_kg_year']['everyday_carbon_feed'];required={**{k:nutrient[k] for k in ['N','P','K','trace']},'carbon':carbon}
    result['nutrient_coproduct_screen']={'annual_element_requirement_kg':required,'required_feed_mass_fraction_at_assumed_recovery':{k:v/feed/cp['nutrient_recovery_comparison'] for k,v in required.items()},'assumed_recovery':cp['nutrient_recovery_comparison'],'Bennu_analogue_contained_carbon_kg_year':feed*cp['carbon_analogue']['mass_fraction'],'Bennu_analogue_recovered_carbon_comparison_kg_year':feed*cp['carbon_analogue']['mass_fraction']*cp['carbon_analogue']['recovery_comparison'],'upper_bound_displaceable_Earth_feed_kg_year':city['resources']['operating_import_kg_year']['food_nutrients_and_carbon']+city['resources']['operating_import_kg_year']['everyday_carbon_feed'],'rule':'Element requirement is distinct from delivered salts or feed compounds. Carbon analogue and trial recovery do not certify production. Nitrogen food demand is additional to atmospheric nitrogen. Process chemistry, toxicity control, food-grade certification, power and delivery must be designed before crediting import substitution. Trace nutrients require element-specific assays; no aggregate trace-metal sales value is assigned.'}
    n=result['nutrient_coproduct_screen'];n['gross_import_avoidance_ceiling_USD2041']=n['upper_bound_displaceable_Earth_feed_kg_year']*(s['prices']['other_special_USD2026_kg']+s['transport']['LEO_USD2026_kg'])*c['cost']['inflation_factor']
    pgm=cp['platinum_sensitivity'];result['trace_metal_sensitivity']=[{'platinum_grade_ppm_trial':ppm,'recovery_trial':pgm['recovery_trial'],'recovered_kg':feed*ppm*1e-6*pgm['recovery_trial'],'share_of_2025_world_mine_output':feed*ppm*1e-6*pgm['recovery_trial']/pgm['world_mine_production_2025_kg'],'gross_content_value_at_frozen_2025_price_USD2025':feed*ppm*1e-6*pgm['recovery_trial']*pgm['USD2025_kg_rounded'],'credited_external_revenue_USD2041':0,'status':'Unassayed grade and recovery trial; contained value before refining, transport and price response. Trace metal is part of aggregate metals, never additive to that mass.'} for ppm in pgm['grade_ppm_trials']]
    result['refined_inventory_valuation_sensitivity']={str(price):{'current_spare_metal_stock_USD2041':result['coproduct_inventory_comparison']['engineering']['stock_metal_kg_year']*price*c['cost']['inflation_factor'],'expanded_refinery_annual_stock_USD2041':q['coproducts']['metal_refinery_expansion_kg_year']*price*c['cost']['inflation_factor']} for price in cp['refined_inventory_prices_USD2026_kg']}
    result['trade_capacity_grid']=[]
    for goods in q['exports']['finished_goods_kg_year_grid']:
        for W in q['exports']['compute_W_grid']:
            a=model.annual(q,src,s,c,**{**focus,'goods':goods,'source_goods':0,'compute':W})
            for price in q['exports']['prices_USD2026_kg']:
                t=external(q,src,s,c,a,goods_price=price)
                result['trade_capacity_grid'].append({'goods_kg_year':goods,'compute_W':W,'goods_USD2026_kg':price,'feasible':a['feasible'],'binding_capacity':a['binding_capacity'],**t['USD2041']})
    result['best_tested_trade_cases']={}
    for price in q['exports']['prices_USD2026_kg']:
        best=max((x for x in result['trade_capacity_grid'] if x['feasible'] and x['goods_USD2026_kg']==price),key=lambda x:x['after_renewal_and_resource'])
        result['best_tested_trade_cases'][str(price)]=best
    best=result['best_tested_trade_cases']['2000'];args=dict(profile='local_structures',service_portfolio=True,goods=best['goods_kg_year'],compute=best['compute_W'],propellant_sales=25e6)
    result['high_value_capacity_case']=evaluate(q,src,s,c,**args,goods_price=2000)
    result['high_value_capacity_commissioning']=commissioning(q,src,s,c,**args)
    f=model.frontier(q,src,s,c,'growth_people',balanced_growth=True,**args)
    result['high_value_capacity_with_growth']=evaluate(q,src,s,c,**args,growth_people=f['maximum'],new_plant_scale=f['maximum']/c['demand']['population'],goods_price=2000)
    a=model.annual(q,src,s,c,**args);zero=external(q,src,s,c,a,goods_price=0)['USD2026'];delivered=a['saleable_export_kg_year']
    result['high_value_break_even']={'saleable_kg_year':delivered,'whole_system_USD2026_per_saleable_kg':-zero['after_renewal_and_resource']/delivered,'whole_system_USD2041_per_saleable_kg':-zero['after_renewal_and_resource']/delivered*c['cost']['inflation_factor'],'rule':'Required average price at this physical portfolio and all other assumed service sales, before debt. Quantity and price must both find qualified buyers.'}
    result['growth']={str(rate):growth_path(q,src,s,c,rate=rate) for rate in [None,0,.01,.03]}
    result['sensitivities']={'active_module_price':{},'service_valuation':{},'inherited_upkeep_stress':{}}
    for mult in q['production']['component_price_sensitivity']:
        v=deepcopy(q);v['production']['special_component_price_multiplier']=mult
        result['sensitivities']['active_module_price'][str(mult)]=evaluate(v,src,s,c,**focus,goods_price=2000,source_price=750,compute_markup=1.5)['external']['USD2041']
    for share in q['valuation']['routine_service_fraction_sensitivity']:
        v=deepcopy(q);v['valuation']['routine_service_value_added_fraction']=share
        result['sensitivities']['service_valuation'][str(share)]=evaluate(v,src,s,c,**focus)['domestic']['USD2041']
    old=model.annual(q,src,s,c,profile='inherited');new=model.annual(q,src,s,c,profile='local_structures')
    penalty=max(0,old['cash_USD2026']['industrial_external_services']-new['cash_USD2026']['industrial_external_services'])*c['cost']['inflation_factor']
    result['sensitivities']['inherited_upkeep_stress']={'additional_USD2041_year':penalty,'qualified_high_value_after_renewal_and_resource_USD2041':result['portfolios']['qualified_high_value']['external']['USD2041']['after_renewal_and_resource']-penalty}
    # Continuous proportional retirement is a screening bound. It does not
    # certify a reduced discrete generator constellation or a new dispatch.
    lo=.05;hi=1.
    for _ in range(50):
        mid=(lo+hi)/2
        if model.annual(q,src,s,c,profile='local_structures',plant_scale=mid)['feasible']:hi=mid
        else:lo=mid
    r=model.annual(q,src,s,c,profile='local_structures',plant_scale=hi)
    result['proportional_rightsizing_bound']={'retained_fraction':hi,'binding_capacity':r['binding_capacity'],'external_purchases_USD2041':r['cash_USD2026']['total_external_purchases']*c['cost']['inflation_factor'],'status':'Continuous proportional-sizing bound only. Discrete fleets, ten-station power redundancy and safe retirement require separate design; retaining the construction system for growth is evaluated explicitly above.'}
    result['initial_capital_context']={'inherited_industrial_establishment_USD2041':c['cost']['capital_with_contingency_USD2041'],'inherited_industrial_and_material_campaign_USD2041':c['cost']['industrial_and_initial_material_campaign_USD2041'],'status':'Preserved space-fuel/local-tug procurement scenario, with its original campaign and replacement assumptions. The new annual modular cases do not silently rerate this sunk campaign; a new initial deployment design must solve the bootstrap and commissioning sequence.'}
    result['currency']={'working_price_convention_USD2041_per_U':q['valuation']['inherited_U_service_normalization_USD2041_per_U'],'adopted_exchange_rate':None,'independent_anchors':src['revenue']['U_calibration'],'rule':'A bookkeeping normalization, purchasing-power comparisons and a market exchange rate are distinct. Domestic production cannot settle USD invoices by choosing a different conversion. Capital recovery or sales per U of capital is a yield, not an exchange rate. A sustainable peg requires funded trade, convertibility and reserve policy; none is fixed here.'}
    result['resource_horizon']={'prior_inventory_screen':src['revenue']['inventory_summary'],'constant_local_structure_feed_kg_year':no_trade['physical']['source_feed_kg_year'],'one_million_year_constant_feed_kg':no_trade['physical']['source_feed_kg_year']*1e6,'rule':'No permitted accessible reserve has been proved. Joint products improve useful recovery but do not restore discharged propellant or eliminate feed. Growth paths include cumulative feed; exhaustion dates require permitted recoverable reserves and route-qualified sources, and local Earth-Moon/near-Earth industry has no demonstrated million-year entitlement.'}
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'research/runs/initial_economy.json');args=parser.parse_args()
    result=run();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    for name,case in result['portfolios'].items():print(name,case['engineering']['feasible'],{k:round(v/1e9,3) for k,v in case['external']['USD2041'].items() if k in ['gross_external_sales','operating_cash_surplus','after_renewal_and_resource']})


if __name__=='__main__':main()
