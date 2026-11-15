"""Compare probabilistic geological stocks with dated virgin-material demands.

Availability fractions are explicit access/preservation scenarios. They are not
probabilities that scientific or governance release conditions will be met.
"""
from __future__ import annotations
import gzip,hashlib,json,math
from pathlib import Path
from . import costs
import numpy as np
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quantiles(a):return np.quantile(a,[.05,.5,.95]).tolist()

def fraction(v):
    if not math.isfinite(v) or not 0<=v<=1:raise ValueError('A fraction must lie in [0,1]')
    return v

def stock_test(stock,demand):
    """Conditional exceedance frequency, with Monte Carlo sampling error alone."""
    a=np.asarray(stock,dtype=float)
    if a.ndim!=1 or len(a)<2 or np.any(~np.isfinite(a)) or np.any(a<0) or not math.isfinite(demand) or demand<0:
        raise ValueError('Finite nonnegative stock draws and demand required')
    p=float(np.mean(a>=demand));n=len(a);z=1.959963984540054
    denominator=1+z*z/n;centre=(p+z*z/(2*n))/denominator
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/denominator
    return {'probability_sufficient_under_model':p,'monte_carlo_frequency_interval_95percent':[max(0,centre-half),min(1,centre+half)],
            'demand_kg':float(demand),'available_stock_kg_p05_p50_p95':quantiles(a),
            'interval_rule':'Wilson interval for simulation frequency only; geological model and engineering uncertainty remain in the stock model and separate scenario branches.'}

def minimum_share(stock,demand,recovery=1.):
    fraction(recovery)
    if recovery==0:return None
    a=np.asarray(stock)*recovery;low,median,high=quantiles(a)
    return {'share_at_median_stock':demand/median if median else None,
            'share_for_95percent_model_stock_exceedance':demand/low if low else None,
            'share_at_p95_stock':demand/high if high else None,
            'rule':'Shares above one exceed this inventory. These are required availability conditions, not predicted access or permission rates.'}

def cumulative_audit(case,horizon):
    rows=[r for r in case['annual'] if r['year']<=horizon]
    if [r['year'] for r in rows]!=list(range(horizon+1)):raise ValueError('Missing dated seed/operating row')
    keys=('source_feed_kg','lunar_raw_kg','source_water_equivalent_kg','source_nitrogen_kg','source_virgin_metal_kg')
    values={k:math.fsum(r[k] for r in rows) for k in keys}
    if any(not math.isfinite(v) or v<0 for v in values.values()):raise ValueError('Invalid demand')
    mature=[r for r in case['annual'] if 201<=r['year']<=220]
    if len(mature)!=20:raise ValueError('Mature demand window incomplete')
    values['mature_annual_mean']={k:math.fsum(r[k] for r in mature)/20 for k in keys}
    values.update(name=case['name'],baseline_case=case['baseline_case'],variant=case['variant'],horizon=horizon,
                  physical_capacity_qualified=case['physical_case']['first_capacity_shortfall'] is None)
    return values

def comparison(stock,demand,available_share,water_recovery,metal_recovery):
    f=fraction(available_share);wr=fraction(water_recovery);mr=fraction(metal_recovery)
    body=np.asarray(stock['body_mass_kg']);water=np.asarray(stock['water_equivalent_kg']);iron=np.asarray(stock['native_iron_kg'])
    if len(body)!=len(water) or len(body)!=len(iron):raise ValueError('Unpaired resource draws')
    if np.any(water+iron>body*(1+1e-10)):raise ValueError('Contained resource masses exceed whole-body mass')
    feed=stock_test(body*f,demand['source_feed_kg']);w=stock_test(water*f*wr,demand['source_water_equivalent_kg'])
    metal=stock_test(iron*f*mr,demand['source_virgin_metal_kg'])
    joint=float(np.mean((body*f>=demand['source_feed_kg'])&(water*f*wr>=demand['source_water_equivalent_kg'])))
    result={'name':demand['name'],'baseline_case':demand['baseline_case'],'variant':demand['variant'],'horizon':demand['horizon'],
        'available_share':f,'water_recovery':wr,'native_iron_recovery':mr,'feed_mass':feed,'water_product':w,
        'native_iron_mass_proxy':metal,'joint_feed_and_water_probability':joint,'physical_capacity_qualified':demand['physical_capacity_qualified'],
        'native_iron_proxy_rule':'Compares native-iron tonnes with generic virgin-metal tonnes as a mass screen only. It establishes no alloy recipe, grade suitability or completed material supply.',
        'nitrogen_rule':'No calibrated population nitrogen inventory is supplied; nitrogen demand remains explicit and is not declared satisfied by water or iron sufficiency.'}
    if demand['horizon']==220:
        annual=demand['mature_annual_mean']
        spans=[]
        for stock_values,used,rate in [(body*f,demand['source_feed_kg'],annual['source_feed_kg']),
            (water*f*wr,demand['source_water_equivalent_kg'],annual['source_water_equivalent_kg'])]:
            if rate<=0:raise ValueError('Positive sustaining demand required')
            spans.append(np.maximum(0,stock_values-used)/rate)
        result['remaining_feed_and_water_stock_years_p05_p50_p95']=quantiles(np.minimum(*spans))
        result['sustaining_rule']='Remaining stock divided by mean fresh demands in years 201–220, with the smaller feed/water span retained. No replenishment, further expansion, new resource discoveries or proven future transport capacity is implied.'
    return result

def baseline_ranges(product):
    result=[]
    for h in (100,200,220):
        records=[]
        for case in product['cases']:
            rows=[r for r in case['annual'] if r['year']<=h];seed=case['seed']
            feed=seed['source_feed_kg']+sum(r['physical']['source_feed_kg_year'] for r in rows)
            lunar=seed['lunar_raw_kg']+sum(r['physical']['lunar_raw_kg_year'] for r in rows)
            water=sum(r['physical']['source_water_kg_year'] for r in rows)
            records.append({'name':case['name'],'source_feed_kg':feed,'lunar_raw_kg':lunar,'operating_water_kg':water,
                'water_with_seed_upper_bound_kg':water+seed['source_feed_kg']*.1*.7})
        result.append({'horizon':h,'case_count':len(records),'ranges':{k:[min(r[k] for r in records),max(r[k] for r in records)] for k in records[0] if k!='name'},
            'seed_water_rule':'The old compact export omits seed product flows; operating water is a lower bound, with seed feed×the inherited 10% grade×70% recovery as an upper bound. New reclamation comparisons use exact seed product totals.'})
    return result

def representative_grade_cost(stock, coefficients, cases, rate, recovery):
    """Uniform-pool mining sensitivity, kept separate from selected ore mining."""
    grade=np.asarray(stock['water_equivalent_kg'])/np.asarray(stock['body_mass_kg'])
    if np.any(grade<=0) or np.any(grade>1):raise ValueError('Positive physical pool grade required')
    recovery=fraction(recovery)
    if recovery==0:raise ValueError('Positive recovery required')
    trials={name:costs.unit_feed_cost(coefficients,case,rate) for name,case in cases.items()}
    return {'water_recovery':recovery,'water_mass_fraction_p05_p50_p95':quantiles(grade),
        'feed_multiplier_relative_to_10percent_grade_at_same_recovery_p05_p50_p95':quantiles(.1/grade),
        'partial_source_water_USD2041_kg_by_cost_scenario':{name:quantiles(v['source_USD2041_kg_feed']/(grade*recovery)) for name,v in trials.items()},
        'rule':'A representative mixture of an entire orbital pool is compared with selectively mining the inherited 10 percent water grade. Named cases vary equipment, energy and availability while holding recovery fixed. This is a yield/cost sensitivity only: no fleet, freight, power-capacity or financing response is solved. Water equivalents include bound hydrogen, whose chemistry may change energy requirements. Partial source-process costs exclude delivery, common infrastructure, labour and unpriced preservation obligations.'}

def load():
    c=json.loads((HERE/'potential_inputs.json').read_text());p={}
    for name,s in c['products'].items():
        path=HERE/s['path']
        if sha(path)!=s['sha256']:raise ValueError('Changed '+name+' product')
        v=json.loads(gzip.decompress(path.read_bytes()))
        if v['schema']!=s['schema']:raise ValueError('Unexpected '+name+' schema')
        p[name]=v
    return c,p

def product():
    c,p=load();cost_cfg=c['partial_costs'];coeff_path=HERE/cost_cfg['coefficients_path']
    if sha(coeff_path)!=cost_cfg['sha256']:raise ValueError('Changed cost coefficients')
    coefficients=json.loads(coeff_path.read_text())
    demand=[cumulative_audit(case,h) for case in p['reclamation']['cases'] for h in c['checkpoints']]
    results=[];requirements=[];grade_costs=[]
    for scenario,s in p['population']['scenarios'].items():
        for pool,pool_data in s['pools'].items():
            stock=pool_data['total']
            grade_costs.append({'scenario':scenario,'pool':pool,**representative_grade_cost(stock,coefficients,cost_cfg['cases'],cost_cfg['capital_rate'],c['water_recovery'])})
            for d in demand:
                requirements.append({'scenario':scenario,'pool':pool,'name':d['name'],'horizon':d['horizon'],
                    'joint_required_availability_share_p05_p50_p95':quantiles(np.maximum(d['source_feed_kg']/np.asarray(stock['body_mass_kg']),d['source_water_equivalent_kg']/(np.asarray(stock['water_equivalent_kg'])*c['water_recovery']))),
                    'joint_rule':'Required availability is calculated jointly within each draw before taking quantiles. Marginal feed/water thresholds below are separate quantities.',
                    'feed':minimum_share(stock['body_mass_kg'],d['source_feed_kg']),
                    'water':minimum_share(stock['water_equivalent_kg'],d['source_water_equivalent_kg'],c['water_recovery'])})
                for f in c['available_resource_shares']:
                    results.append({'population_scenario':scenario,'pool':pool,**comparison(stock,d,f,c['water_recovery'],c['native_iron_recovery'])})
    return {'schema':'terluna.resources.potential-and-reclamation/1','producer':{'path':'resources/source_access/potential.py','sha256':sha(__file__)},
        'dependencies':{name:sha(HERE/name) for name in ['potential_inputs.json','potential_sources.json','costs.py',cost_cfg['coefficients_path'],*[s['path'] for s in c['products'].values()]]},
        'evidence':'Conditional probabilistic geological-stock assessment and stock-only sufficiency comparisons against separately solved physical reclamation scenarios. No probability is assigned to permits, preservation release, launch windows or future engineering qualification.',
        'population_source_schema':p['population']['schema'],'lunar_potential':p['lunar']['bulk_inventory'],
        'baseline_demand_ranges':baseline_ranges(p['baseline']),'demand_checkpoints':demand,'comparisons':results,'required_shares':requirements,
        'representative_pool_grade_costs':grade_costs,'rules':c['rules'],'scientific_protection_release':None,'extraction_authorized':False}

def main():
    x=product();path=HERE/'results/potential.json.gz';path.write_bytes(gzip.compress((json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0));print(path)
if __name__=='__main__':main()
