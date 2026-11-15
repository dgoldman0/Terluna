"""Read-only response export from provisioning commit146a5582; no repo writes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--producer-root',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args_cli=parser.parse_args()
ROOT=args_cli.producer_root.resolve()
expected_commit='146a558275f486c9fd9a5cf79908b7c6fdf21888'
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==expected_commit,'Use the pinned provisioning checkout'
subprocess.run(['git','diff','--exit-code','HEAD','--','research/studies/initial_industrial_system','shared'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT))
from research.studies.initial_industrial_system import economy as E
from research.studies.initial_industrial_system import evolution as M

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def identity(path):return {'path':path,'sha256':digest(ROOT/path)}

q,src,s,c=M.read_inputs()
baseargs={'profile':'local_structures','population':250000,'plant_scale':1,'service_portfolio':True}
# Name, physical arguments, quoted prices. Each row is an alternative, not additive.
cases=[
 ('domestic_services',{},{}),
 ('source_refined',{'source_goods':100e6},{'source_price':30}),
 ('source_refinery_inventory',{'refinery_capacity':1e9,'extra_metal':1e9},{'inventory_price':10}),
 ('commodity_mixed',{'goods':100e6,'source_goods':50e6,'compute':1e9,'propellant_sales':25e6},{'goods_price':3,'source_price':3}),
 ('specialty_mixed',{'goods':100e6,'source_goods':50e6,'compute':1e9,'propellant_sales':25e6},{'goods_price':750,'source_price':30}),
 ('premium_lunar',{'goods':100e6},{'goods_price':2000}),
 ('goods_heavy',{'goods':200e6},{'goods_price':2000}),
 ('current_large',{'goods':200e6,'compute':3e9,'propellant_sales':25e6},{'goods_price':2000}),
 ('logistics_250kt',{'propellant_sales':250e6},{}),
 ('logistics_1Mt',{'propellant_sales':1e9},{}),
 ('computing_5GW',{'compute':5e9},{}),
 ('balanced',{'goods':100e6,'source_goods':25e6,'compute':3e9,'propellant_sales':250e6},{'goods_price':750,'source_price':30}),
]
base=E.evaluate(q,src,s,c,**baseargs)
I=c['cost']['inflation_factor']

def export_case(name,physical,prices):
 args={**baseargs,**physical};pp={'goods_price':750,'source_price':30,'compute_markup':1,**prices}
 a=E.evaluate(q,src,s,c,**args,**pp)
 commissioning=E.commissioning(q,src,s,c,**args)
 growth=M.frontier(q,src,s,c,'growth_people',balanced_growth=True,**args)
 n=growth['maximum']
 growing=E.evaluate(q,src,s,c,**args,**pp,growth_people=n,new_plant_scale=n/250000)
 ratios=a['engineering']['capacity_use']
 growthratios=growing['engineering']['capacity_use']
 def headroom(r):
  fraction={k:1-v for k,v in r.items()}
  mean=min((1-r['ground_firm_power'])*c['engineering']['total_ground_firm_W']*s['lunar']['process_availability']/s['economics']['capacity_margin'],(1-r['lunar_process_power'])*c['peak_flows']['local_mean_W'])
  return {'fractions':fraction,'additional_lunar_mean_power_W_independent_screen':max(0,mean)}
 a['input_arguments']={'physical':args,'price':pp}
 a['commissioning']=commissioning
 a['remaining_balanced_growth']={'status':growth['status'],'people_year':n,'fraction_of_residents_year':n/250000,'at_limit':growing}
 samples=[]
 for f in (.25,.5,.75):
  z=E.evaluate(q,src,s,c,**args,**pp,growth_people=n*f,new_plant_scale=n*f/250000)
  predicted={k:ratios[k]+f*(growthratios[k]-ratios[k]) for k in ratios}
  max_excess=max(z['engineering']['capacity_use'][k]-predicted[k] for k in ratios)
  actual=z['external']['USD2041']['total_imports'];prediction=a['external']['USD2041']['total_imports']+f*(growing['external']['USD2041']['total_imports']-a['external']['USD2041']['total_imports'])
  samples.append({'fraction':f,'max_capacity_ratio_above_linear_interpolation':max_excess,'external_imports_actual_minus_linear_USD2041':actual-prediction})
 a['growth_interpolation_controls']=samples
 # Expansion must install the commercial stock as well as its host industry.
 lo,hi=0.,n
 for _ in range(50):
  mid=(lo+hi)/2
  trial=M.annual(q,src,s,c,**args,growth_people=mid,new_plant_scale=mid/250000,commission_fraction=mid/250000)
  if max(trial['capacity_use'].values())<=1:lo=mid
  else:hi=mid
 joint=E.evaluate(q,src,s,c,**args,**pp,growth_people=lo,new_plant_scale=lo/250000,commission_fraction=lo/250000)
 a['joint_replication']={'people_year':lo,'fraction_of_residents_year':lo/250000,'at_limit':joint,'rule':'Host and portfolio commercial assets expand together; qualification remains separately included in the conservative full commissioning allowance.'}
 a['headroom']=headroom(ratios)
 a['remaining_balanced_growth']['headroom_at_growth_limit']=headroom(growthratios)
 a['increment_relative_to_domestic_services_USD2041']={k:a['external']['USD2041'][k]-base['external']['USD2041'][k] for k in a['external']['USD2041']}
 # Imported IT and route costs are already priced in the engineering row. Only
 # commercial Earth service purchases change when compute price changes.
 a['repricing']={'delivered_lunar_goods_kg_year':args.get('goods',0)*(1-src['revenue_inputs']['return']['loss_probability_flight'])*src['revenue_inputs']['manufacturing']['delivered_quality_yield'],
  'delivered_source_goods_kg_year':args.get('source_goods',0)*(1-src['revenue_inputs']['return']['loss_probability_flight'])*src['revenue_inputs']['manufacturing']['delivered_quality_yield'],
  'compute_W':args.get('compute',0),'propellant_kg_year':args.get('propellant_sales',0),
  'compute_external_service_fraction_of_receipt_change':src['revenue_inputs']['computing']['Earth_service_share_of_receipts'],
  'price_rule':'Prices are USD2026/unit; multiply receipt changes once by inflation_factor to obtain constantUSD2041. All other physical costs held fixed; changing product identity or demandedquantity needs a new physical row. No extra total-cost percentage for goods or propellant is imposed on a repricing-only change.'}
 for field,sector,price in [('delivered_lunar_goods_kg_year','lunar_goods_delivered',pp['goods_price']),('delivered_source_goods_kg_year','source_goods_delivered',pp['source_price'])]:
  expected=a['repricing'][field]*price*I
  actual=a['external']['sales_by_activity_USD2041'][sector]
  assert abs(actual-expected)<=max(1e-3,abs(actual)*1e-12),(name,field,actual,expected)
 return a
rows={name:export_case(name,args,prices) for name,args,prices in cases}
old=json.loads((ROOT/'research/studies/initial_industrial_system/results/economy.json').read_text())
maxargs={k:old['frontiers'][k]['maximum'] for k in ('goods','source_goods','compute','propellant_sales')}
control=export_case('all_individual_maxima',maxargs,{})
assert not control['engineering']['feasible']
for name in ('commodity_mixed','specialty_mixed'):
 oldkey=name.split('_')[0]
 assert rows[name]['external']==old['portfolios'][oldkey]['external']
 assert rows[name]['domestic']==old['portfolios'][oldkey]['domestic']
allpaths=['research/studies/initial_industrial_system/'+x for x in ('economy.py','evolution.py','run.py','transport.py','revenue.py','economy_inputs.json','economy_sources.json')]+['shared/constants.json','shared/constants.py','shared/provenance.py']
result={'schema':'terluna.research.strategy-profile-export/1','source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
 'producer':{'export_recipe':'export_terluna_strategy_profiles.py','export_script_sha256':digest(__file__),'sources':[identity(p) for p in allpaths]},
 'source_inputs':q['inputs'],'matched_existing_product':identity('research/studies/initial_industrial_system/results/economy.json'),
 'evidence':'Read-only evaluation of the pinned provisioning producer. New portfolio combinations, not new engineering calibration. All markets, prices, source permission and equipment qualification remain conditional.',
 'units':{'money':'constantUSD2041 unless field explicitlyUSD2026','mass_flow':'kg/year','power':'W'},'inflation_factor':I,
 'baseline':base,'profiles':rows,'rejected_controls':{'all_individual_maxima':control},
 'physical_capacity':{'lunar_process_availability':s['lunar']['process_availability'],'capacity_margin':s['economics']['capacity_margin'],'ground_firm_W':c['engineering']['total_ground_firm_W'],'peak_flows':{k:c['peak_flows'][k] for k in ('source_mean_W','source_heat_W','source_metal_kg_year','source_ore_kg_year','local_mean_W','metal_kg_year','depot_mean_W')}},
 'market_reference':{'Earth_propellant_competitor_USD2026_kg':M.services(src,s,c,True)['Earth_fuel_alternative_USD2026_kg'],'Earth_propellant_price_fraction_ceiling':src['revenue_inputs']['propellant_customers']['Earth_alternative_price_fraction_ceiling'],'reference_shared_research_and_mission_budget_USD2026_year':src['revenue_inputs']['markets']['reference']['mission_customer_budget_USD2026_year'],'reference_research_receipts_USD2026_year':src['revenue_inputs']['markets']['reference']['research_receipts_USD2026_year'],'propellant_USD2026_kg':src['revenue_inputs']['markets']['reference']['propellant_USD2026_kg'],'compute_USD2026_per_installed_W_year':E.compute_price(q),'compute_Earth_service_fraction':src['revenue_inputs']['computing']['Earth_service_share_of_receipts']},
 'limits':['Profiles are alternatives on one retained founding system. Their standalone maxima and growth limits cannot be summed.',
 'One-at-a-time headroom is not jointly spendable by all services; combining additions requires a new shared-capacity calculation.',
 'source_goods denotes refined/manufactured source products, not raw residue or assayed precious-metal output.',
 'Premium product prices do not establish a low-gravity advantage; recipe, certification, costs and customer evidence remain independent.',
 'Propellant receipts are final Earth-mission purchases only; internal logistics sales cancel. Independent mission budgets and competitive delivered price must bound actual sales. The economy producer does not apply the earlier mission budget to its explicit propellant_sales argument. Source reference budget is30B USD2026 including5B research, leaving83.333million kg fuel at300USD2026/kg; larger logistics rows are physical/requested-demand sensitivities and require new independent budgets.',
 'Commissioning includes the source qualification programme for each alternative. Do not add whole-alternative commissioning totals when combining profiles. Shared qualification and base services count once.',
 'No minimum small-module construction or array settlement is inferred.'],
 'verification':{'existing_mixed_cases':'Commodity and specialty external/domestic outputs exactly reproduce saved economy product.','control':'Combined independent maxima rejected by the source physical model.','inputs':'M.read_inputs verified all recorded input hashes.'}}
args_cli.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print('Wrote',args_cli.output)
for name,a in rows.items():
 ext=a['external']['USD2041'];print(name,'feasible',a['engineering']['feasible'],'salesB',round(ext['gross_external_sales']/1e9,3),'netB',round(ext['after_renewal_and_resource']/1e9,3),'growth',round(a['remaining_balanced_growth']['people_year'],2),'binding',a['engineering']['binding_capacity'])
print('rejected maxima',control['engineering']['feasible'],control['engineering']['binding_capacity'],max(control['engineering']['capacity_use'].values()))
