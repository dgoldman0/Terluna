"""Reappraise staged returns with pinned reclamation and resource distributions.

The existing financial engine remains unchanged. Physical schedules, commercial
orders, monetary claims and geological stock conditions retain separate accounts.
"""
from __future__ import annotations
from copy import deepcopy
import gzip,hashlib,json,math,time
from pathlib import Path
import numpy as np
from . import finance,bottlenecks
from ..complex_returns import credit,policy
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):
    raw=Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith('.gz') else raw)
def write(path,value):
    raw=(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_bytes(gzip.compress(raw,mtime=0) if str(path).endswith('.gz') else raw)
def quantiles(values):return np.quantile(values,[.05,.5,.95]).tolist()


def inputs():
    spec=read(HERE/'revised_inputs.json');products={}
    for path,digest in spec['dependencies'].items():
        if sha(HERE/path)!=digest:raise ValueError('Changed dependency: '+path)
    for name,item in spec['products'].items():
        path=HERE/item['path']
        if sha(path)!=item['sha256']:raise ValueError('Changed product: '+name)
        data=read(path)
        if data['schema']!=item['schema']:raise ValueError('Unexpected schema: '+name)
        products[name]=data
    cfg,orders,old=finance.inputs()
    physical=products['physical']
    for remote,local in [('provisioning/enclosed_industry/inputs/orders.json','results/orders.json'),('provisioning/enclosed_industry/inputs/revised_rollout_orders.json','results/revised_orders.json')]:
        if physical['dependencies'].get(remote)!=sha(HERE/local):raise ValueError('Physical orderbook mismatch: '+remote)
    order_cases={c['name']:c for product in (orders,products['rollout_orders']) for c in product['cases']}
    for case in physical['cases']:
        request=order_cases[case['order_case_name']]
        identity=hashlib.sha256((json.dumps(request,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()).hexdigest()
        if identity!=case['source_order_case_sha256']:raise ValueError('Physical case used different customer orders: '+case['name'])
    return spec,cfg,orders,old,products


def valid(physical):
    return physical['convergence']['converged'] and physical['first_capacity_shortfall'] is None


def sponsor_cash(path,sim,horizon):
    start=path['annual'][0]['year'];cash=[0.]*(horizon+1)
    for y in range(min(start,horizon+1)):cash[y]=-path['seed_hardware_and_carry']/start
    cash[0]-=sim['initial_reserve']
    if start-1<=horizon:cash[start-1]-=path['initial_working']
    for r in sim['annual']:
        if r['year']<=horizon:cash[r['year']]=r['owner_payout']-r['contribution']
    return cash


def measures(path,sim,cfg,horizon):
    m=finance.metrics(path,sim,cfg,horizon);rows=[r for r in sim['annual'] if r['year']<=horizon]
    initial=path['seed']+sim['initial_reserve'];contributions=initial+m['additional_equity_required'];cf=sponsor_cash(path,sim,horizon)
    m.update(initial_full_funding=initial,all_funder_contributions=contributions,
             cumulative_gross_receipts=sum(r['gross'] for r in rows),
             undiscounted_payout_multiple=m['owner_payouts']/contributions,
             first_payout_year=next((r['year'] for r in rows if r['owner_payout']>1),None),
             real_IRR_roundoff_screened=finance.valuation.nonnegative_irr([0. if abs(v)<1 else v for v in cf]),
             sustaining_years_with_negative_operating_cash=[r['year'] for r in rows if r['year']>=200 and r['operating_available']<0],
             evidence='Conditional funded-completion account; additional equity, if required, must be supplied. Cash distributions, retained financial claims and physical assets are separate valuations.')
    return m


def stratification(case,physical,path,sim,cfg,periods):
    rows=[]
    for request,pr,r,f in zip(case['annual'],physical['annual'],path['annual'],sim['annual']):
        costs=bottlenecks.annual_costs(request,pr,r,cfg)
        rows.append({'year':r['year'],'revenue':r['revenue_by_sector'],'costs':costs,'payout':f['owner_payout'],'contribution':f['contribution']})
    result=[];start=path['annual'][0]['year']
    for first,last in periods:
        selected=[r for r in rows if first<=r['year']<=last];values={}
        for rate in (0.,.03):
            revenues={k:math.fsum(r['revenue'][k]/(1+rate)**r['year'] for r in selected) for k in rows[0]['revenue']}
            costs={k:math.fsum(r['costs'][k]/(1+rate)**r['year'] for r in selected) for k in rows[0]['costs']}
            costs['seed_deployment_and_allowance_USD2041']=math.fsum(path['seed_hardware_and_carry']/start/(1+rate)**y for y in range(start) if first<=y<=last)
            values[str(rate)]={'revenue_by_sector':revenues,'costs_by_category':costs,'net_real_cash':sum(revenues.values())-sum(costs.values()),
                              'distributions':math.fsum(r['payout']/(1+rate)**r['year'] for r in selected),'additional_equity':math.fsum(r['contribution']/(1+rate)**r['year'] for r in selected)}
        result.append({'first_year':first,'last_year':last,'values':values})
    return {'periods':result,'financing_bridges':{str(h):bottlenecks.carry_bridge(path,sim,cfg,h,.03) for h in cfg['horizons']},
            'peak_retained_cash':max(({'year':r['year'],'cash':r['cash']} for r in sim['annual']),key=lambda r:r['cash'])}


def export_flow(case,path,sim,cfg,physical,variant='combined'):
    """Expose dated counterparties for the unchanged wealth-accounting engine."""
    rows=[];loans=[];start=path['annual'][0]['year']
    for y in range(start):
        funding=path['seed_hardware_and_carry']/start+(sim['initial_reserve'] if y==0 else 0)+(path['initial_working'] if y==start-1 else 0)
        rows.append({'year':y,'population_convention':path['seed_population']/2,'gross':0.,'owner_payout':0.,'dividend':0.,'owner_custody':0.,'resident_payroll':path['seed_resident_payroll']/start,'cash':path['initial_working'] if y==start-1 else 0.,'reserves':sim['initial_reserve'],'outside_equity_contribution':funding,'contribution':0.,'debt':0.,'local_U_debt':0.,'supplier_U_debt':0.,'external_debt':0.,'local_U_draw':0.,'supplier_U_draw':0.,'external_draw':0.,'local_U_principal':0.,'supplier_U_principal':0.,'external_principal':0.,'local_U_interest':0.,'supplier_interest':0.,'external_interest':0.,'resource_provision':physical['seed']['total_resource_provision_USD2041']/start,'external_operating':None,'capital':None,'service_residual_capital':None,'service_residual_renewal':None})
    for r,f in zip(path['annual'],sim['annual']):
        due=credit.repay(loans,f['year'])
        for field,kind,rate,share in [('local_U_draw','local',cfg['local_U_rate'],0),('supplier_U_draw','supplier',cfg['external_rate'],1),('external_draw','external',cfg['external_rate'],0)]:
            if f[field]>1e-6:loans.append({'year':f['year'],'original':f[field],'remaining':f[field],'rate':rate,'term':cfg['loan_term_years'],'grace':cfg['underwriting']['principal_grace_years'],'kind':'external' if kind=='external' else 'U','supplier_share':share,'sector_kind':kind})
        balances={k:sum(l['remaining'] for l in loans if l['sector_kind']==k) for k in ('local','supplier','external')}
        if not math.isclose(sum(balances.values()),f['debt'],abs_tol=1.,rel_tol=1e-12):raise ValueError('Debt cohort mismatch')
        rows.append({**r,**f,'outside_equity_contribution':f['contribution'],'local_U_debt':balances['local'],'supplier_U_debt':balances['supplier'],'external_debt':balances['external'],'local_U_principal':due['U_principal']-due['supplier_principal'],'supplier_U_principal':due['supplier_principal'],'external_principal':due['principal']-due['U_principal']})
    return {'name':case['name'],'baseline_case':case['name'],'variant':variant,'portfolio':path['portfolio'],'growth':path['growth'],'market':path['market'],'seed':path['seed'],'initial_reserve':sim['initial_reserve'],'seed_cost_scope':path['seed_cost_scope'],'unpriced_deployment_requirements':path['unpriced_deployment_requirements'],'annual':rows,'baseline_checkpoints':[measures(path,sim,cfg,h) for h in cfg['horizons']]}


def evaluate(case,physical,cfg,periods,variant):
    if not valid(physical):
        return {'name':case['name'],'variant':variant,'status':'physical_requirements_unresolved','first_capacity_shortfall':physical['first_capacity_shortfall'],'convergence':physical['convergence']}
    path=finance.plan(case,physical,cfg);eq=finance.finance(path,cfg,'equity');sims={'equity':eq}
    sims.update({mode:finance.finance(path,cfg,mode,eq) for mode in ('funded_credit','U_credit')})
    checkpoints=[{**measures(path,sim,cfg,h),'variant':variant} for sim in sims.values() for h in cfg['horizons']]
    result={'name':case['name'],'variant':variant,'portfolio':case['portfolio'],'market':case['market'],'growth':case['market_growth_name'],'status':'conditional_account','checkpoints':checkpoints,
            'physical_peak_workforce_ratio':max(r['capacity_summary']['resident_workforce'] for r in physical['annual'])}
    if case['market']=='breakthrough' and case['market_growth_name'] in ('reference','wider'):
        result['stratification']=stratification(case,physical,path,sims['U_credit'],cfg,periods)
    return result,path,sims


def cumulative_demand(physical,audit,h,seed_bounds=None):
    if audit is not None:
        rows=[r for r in audit['annual'] if r['year']<=h]
        if [r['year'] for r in rows]!=list(range(h+1)):raise ValueError('Missing dated physical audit')
        total=lambda key:math.fsum(r[key] for r in rows)
        return {'feed_kg':total('source_feed_kg'),'lunar_raw_kg':total('lunar_raw_kg'),'water_kg_bounds':[total('source_water_equivalent_kg')]*2,'nitrogen_kg':total('source_nitrogen_kg'),'virgin_metal_kg':total('source_virgin_metal_kg'),'seed_rule':'Exact inherited seed total counted once; seed years use its eight-year average.'}
    rows=[r for r in physical['annual'] if r['year']<=h];seed=physical['seed'];water=math.fsum(r['physical']['source_water_kg_year'] for r in rows)
    water_bounds=[water,water+seed['source_feed_kg']*.1*.7]
    seed_rule='Frozen compact baseline omits seed coproducts. Water lower bound excludes seed; upper bound adds seed feed times inherited 10 percent grade and 70 percent recovery.'
    if seed_bounds is not None:
        identity=hashlib.sha256((json.dumps(physical,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()).hexdigest()
        if identity!=seed_bounds['source_case_sha256']:raise ValueError('Baseline seed case mismatch')
        water_bounds=[water+seed_bounds['products']['water'][k] for k in ('total_lower_kg','total_upper_kg')]
        seed_rule=seed_bounds['evidence']
    return {'feed_kg':seed['source_feed_kg']+math.fsum(r['physical']['source_feed_kg_year'] for r in rows),'lunar_raw_kg':seed['lunar_raw_kg']+math.fsum(r['physical']['lunar_raw_kg_year'] for r in rows),
            'water_kg_bounds':water_bounds,
            'nitrogen_kg':None,'virgin_metal_kg':None,
            'seed_rule':seed_rule}


def resource_account(physical,audit,population,spec,seed_bounds=None):
    if not physical.get('annual') or (audit is not None and not audit.get('annual')):
        return {'checkpoints':[],'qualified_source_supply':None,'extraction_authorized':False,'rule':'Physical requirements unresolved; no completed delivery or stock account admitted.'}
    values=[]
    recovery=audit.get('source_coefficients',{}).get('water_recovery',spec['water_recovery']) if audit else spec['water_recovery']
    if not 0<recovery<=1:raise ValueError('Invalid source recovery')
    for h in spec['horizons']:
        d=cumulative_demand(physical,audit,h,seed_bounds);trials=[]
        for scenario,s in population['scenarios'].items():
            for pool in ('wider','all_no_don_quixote'):
                stock=s['pools'][pool]['total'];body=np.asarray(stock['body_mass_kg']);water=np.asarray(stock['water_equivalent_kg'])*recovery
                lo=np.maximum(d['feed_kg']/body,d['water_kg_bounds'][0]/water);hi=np.maximum(d['feed_kg']/body,d['water_kg_bounds'][1]/water)
                trials.append({'scenario':scenario,'pool':pool,'joint_required_share_p05_p50_p95_lower_demand':quantiles(lo),'joint_required_share_p05_p50_p95_upper_demand':quantiles(hi),
                               'availability':[{'share':f,'joint_stock_sufficiency_frequency_bounds':[float(np.mean(hi<=f)),float(np.mean(lo<=f))]} for f in spec['availability_shares']]})
        values.append({'horizon':h,'demand':d,'trials':trials})
    return {'checkpoints':values,'water_recovery_used':recovery,'qualified_source_supply':None,'extraction_authorized':False,
            'rule':'Conditional geological-stock tests, not delivery/permission probabilities or NPV risk weights. Access fractions do not grant extraction. Water recovery is applied once; all physical grades and route obligations remain with the producer. Nitrogen, alloys, archives, future reserves and site qualification remain separate.'}


def risk_cases(case,physical,cfg,spec):
    result=[]
    for shock in spec['receipt_shocks']:
        p=finance.plan(case,physical,cfg);eq=finance.finance(p,cfg,'equity',receipt_shock=shock);u=finance.finance(p,cfg,'U_credit',eq,receipt_shock=shock)
        result.append({'kind':'temporary_receipt_loss','shock':shock,'checkpoints':[measures(p,u,cfg,h) for h in cfg['horizons']]})
    for allowance in (0.,1.):
        alt=deepcopy(cfg);alt['capital_cost_allowance_fraction']=allowance;p=finance.plan(case,physical,alt);eq=finance.finance(p,alt,'equity');u=finance.finance(p,alt,'U_credit',eq)
        result.append({'kind':'capital_cost_allowance','fraction':allowance,'checkpoints':[measures(p,u,alt,h) for h in cfg['horizons']]})
    for multiple in (.5,2.):
        alt=deepcopy(cfg);alt['monetary_basket_multiplier']=multiple;p=finance.plan(case,physical,alt);eq=finance.finance(p,alt,'equity');u=finance.finance(p,alt,'U_credit',eq)
        result.append({'kind':'transaction_basket','multiplier':multiple,'checkpoints':[measures(p,u,alt,h) for h in cfg['horizons']]})
    alt=deepcopy(cfg);alt['local_U_rate']=cfg['external_rate'];p=finance.plan(case,physical,alt);eq=finance.finance(p,alt,'equity');u=finance.finance(p,alt,'U_credit',eq)
    result.append({'kind':'no_local_U_rate_relief','local_U_rate':alt['local_U_rate'],'checkpoints':[measures(p,u,alt,h) for h in cfg['horizons']]})
    alt=deepcopy(cfg);alt['distribution_policy']='retain_to_completion';p=finance.plan(case,physical,alt);eq=finance.finance(p,alt,'equity');u=finance.finance(p,alt,'U_credit',eq)
    result.append({'kind':'retain_surplus_through_completion','checkpoints':[measures(p,u,alt,h) for h in cfg['horizons']]})
    return result


def payment_episodes(case,path,sim,cfg,physical,spec):
    assumptions=read(HERE/'../complex_returns/policy_inputs.json');flow=export_flow(case,path,sim,cfg,physical)
    rows={r['year']:r for r in flow['annual']};result=[]
    for h in spec['horizons']:
        r=rows[h]
        snap={'revenue':r['gross'],'operating_cost':r['operating_cost'],'core_imports':r['external_operating'],
              'debt_service':r['interest']+r['principal'],'external_debt_service':r['external_interest']+r['external_principal'],
              'reserves':r['reserves'],'existing_U_loans':r['local_U_debt']+r['supplier_U_debt'],'existing_convertible_U':r['local_U_debt'],
              'basket':r['domestic_transaction_basket'],'price':r['price'],'relative_price':r['relative_price'],
              'holding_years':cfg['monetary']['domestic_holding_days']/365,'elasticity':cfg['monetary']['supply_elasticity'],
              'cumulative_limit':cfg['monetary']['cumulative_relative_price_limit'],'seed':path['seed'],
              'baseline_annual_payout':r['owner_payout'],'baseline_annual_new_orders':r['capital']}
        for available in (1.,.25):
            snap['cash']=r['cash']*available
            for months in (6,12):
                for name in ('cash_only','U_term_50','U_term_50_half_FX_access','external_bridge'):
                    result.append({'horizon':h,'cash_access_fraction':available,'months':months,'policy':name,
                                   'result':policy.bridge(snap,assumptions['liquidity'],assumptions['policies'][name],months)})
    return {'cases':result,'rule':'Independent monthly episodes at fixed annual activity. Ordinary payouts, new investment and reserve contributions pause; this is operating-payment resilience, not uninterrupted construction. Only the stated fraction of programme cash is available; custody is excluded. Delays collect later and are distinct from permanent losses. Existing funded-equity requirements remain obligations, not promised facilities.'}

def compare(old,new):
    keys=('initial_full_funding','additional_equity_required','owner_payouts','cumulative_gross_receipts','annual_operating_available','annual_post_capital_cash','terminal_cash','terminal_reserves','terminal_debt','cumulative_U_draws','maximum_relative_price')
    return {'name':new['name'],'mode':new['mode'],'horizon':new['horizon'],
            'delta':{k:new[k]-old[k] for k in keys},
            'discounted_delta':{r:{k:new['discounted'][r][k]-old['discounted'][r][k] for k in ('sponsor_distributed_cash_NPV','sponsor_cash_and_terminal_NPV','unlevered_cash_NPV','full_price_socioeconomic_NPV','half_price_socioeconomic_NPV')} for r in new['discounted']},
            'old_self_financing':old['no_additional_equity_to_checkpoint'],'revised_self_financing':new['no_additional_equity_to_checkpoint']}


def verify_original(result,reference):
    for row in result['checkpoints']:
        old=next(x for x in reference if x['name']==row['name'] and x['mode']==row['mode'] and x['horizon']==row['horizon'])
        for k in ('owner_payouts','additional_equity_required','terminal_cash','terminal_reserves','terminal_debt','annual_operating_available'):
            if row[k]!=old[k]:raise ValueError('Changed baseline: '+row['name']+' '+k)
        if row['discounted']!=old['discounted']:raise ValueError('Changed discounted baseline: '+row['name'])


def product():
    spec,cfg,orders,old,products=inputs();cfg=deepcopy(cfg);cfg['horizons']=spec['horizons']
    requests={c['name']:c for c in orders['cases']};rollout={c['name']:c for c in products['rollout_orders']['cases']}
    seed_bounds={c['name']:c for c in products['physical']['baseline']['seed_material_intervals']}
    original=read(HERE/'results/returns.json')['checkpoints'];old_physical={c['name']:c for c in old['physical']['cases'] if c['name'] in requests}
    primary={c['baseline_case']:c for c in products['physical']['cases'] if c['variant']=='combined' and c['baseline_case'] in requests}
    if set(primary)!=set(requests):raise ValueError('Full combined portfolio/market/growth coverage required')
    results=[];flows=[];conditions=[];differences=[];risks=[];sensitivity=[]
    cache=ROOT/'research/runs/enclosed_returns/revised';cache.mkdir(parents=True,exist_ok=True)
    def run(case,physical,variant):
        identity=hashlib.sha256(json.dumps({'producer':sha(__file__),'cfg':cfg,'case':case,'physical':physical,'periods':spec['periods'],'finance':sha(HERE/'finance.py'),'stratification':sha(HERE/'bottlenecks.py')},sort_keys=True).encode()).hexdigest()
        p=cache/(identity+'.json.gz')
        if p.exists():return read(p)
        value=evaluate(case,physical,cfg,spec['periods'],variant)
        if isinstance(value,dict):return value
        record,path,sims=value
        saved={'record':record,'path':path,'sims':sims};write(p,saved);return saved
    for name,case in requests.items():
        start=time.monotonic();baseline=run(case,old_physical[name],'baseline');revised=run(case,primary[name]['physical_case'],'combined')
        verify_original(baseline['record'],original)
        results.extend([baseline['record'],revised['record'] if 'record' in revised else revised])
        conditions.append({'name':name,'variant':'baseline',**resource_account(old_physical[name],None,products['population'],spec,seed_bounds[name])})
        conditions.append({'name':name,'variant':'combined',**resource_account(primary[name]['physical_case'],primary[name],products['population'],spec)})
        if 'record' in revised:
            a={(r['mode'],r['horizon']):r for r in baseline['record']['checkpoints']}
            differences.extend(compare(a[r['mode'],r['horizon']],r) for r in revised['record']['checkpoints'])
            if case['market']=='breakthrough' and case['market_growth_name'] in ('reference','wider'):
                flow=export_flow(case,revised['path'],revised['sims']['U_credit'],cfg,primary[name]['physical_case'])
                flow['source_case_sha256']=hashlib.sha256((json.dumps(primary[name]['physical_case'],sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()).hexdigest()
                flow['source_loop_case_sha256']=hashlib.sha256((json.dumps(primary[name],sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()).hexdigest()
                flows.append(flow)
                risks.append({'name':name,'comparisons':risk_cases(case,primary[name]['physical_case'],cfg,spec),'payment_episodes':payment_episodes(case,revised['path'],revised['sims']['U_credit'],cfg,primary[name]['physical_case'],spec)})
        print(name,round(time.monotonic()-start,2),'seconds',flush=True)
    # The earlier six-variant experiment supplies separately costed interventions.
    experiments=[c for c in products['reclamation']['cases'] if c['variant'] not in ('baseline','combined')]
    experiments += [c for c in products['physical']['cases'] if c['variant']!='combined']
    for ex in experiments:
        order_name=ex.get('order_case_name',ex['baseline_case']);case=deepcopy(requests[order_name] if order_name in requests else rollout[order_name]);case['name']=ex['name']
        value=run(case,ex['physical_case'],ex['variant']);record=value['record'] if 'record' in value else value
        record.update(baseline_case=ex['baseline_case'],order_case_name=order_name)
        record['resources']=resource_account(ex['physical_case'],ex,products['population'],spec)
        sensitivity.append(record)
    margins=[]
    for name,case in requests.items():
        if case['portfolio'] in spec['margin_portfolios'] and case['market']=='breakthrough' and case['market_growth_name']=='reference' and valid(primary[name]['physical_case']):
            margins.append(bottlenecks.margin_target(case,primary[name]['physical_case'],cfg))
    dependencies={path:sha(HERE/path) for path in ('revised.py','revised_inputs.json','revised_sources.json',*spec['dependencies'],*[p['path'] for p in spec['products'].values()])}
    result={'schema':'terluna.humanities.revised-enclosed-returns/1','producer':{'path':'research/studies/enclosed_returns/revised.py','sha256':sha(__file__)},'dependencies':dependencies,
            'evidence':spec['evidence'],'rules':spec['rules'],'primary_cases':results,'matched_differences':differences,'resource_conditions':conditions,'physical_sensitivities':sensitivity,'financial_risks':risks,'commercial_margin_targets':margins,
            'original_baseline_reproduction':{'cases':len(requests),'modes':3,'horizons':spec['horizons'],'exact':True}}
    flow={'schema':'terluna.humanities.wealth-finance-flows/1','producer':result['producer'],'dependencies':dependencies,'source_physical_sha256':spec['products']['physical']['sha256'],
          'evidence':'Revised combined-case financing with the unchanged financial rules; annual debt counterparties and cash exactly reconcile. Ownership changes are applied downstream, not added to programme receipts.','cases':flows}
    return result,flow


def main():
    result,flows=product();write(HERE/'results/revised.json.gz',result);write(HERE/'results/revised_flows.json.gz',flows)
    print('Primary cases',len(result['primary_cases']),'wealth flows',len(flows['cases']),flush=True)

if __name__=='__main__':main()
