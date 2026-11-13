"""Render integrated export-strategy comparisons from the reviewed product."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
SOURCE=ROOT/'research/studies/complex_returns/results/strategies.json'
OUTPUT=ROOT/'research/runs/complex_returns/figures';OUTPUT.mkdir(parents=True,exist_ok=True)
r=json.loads(SOURCE.read_text())
identity={'product_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
names=['raw_resource_price','commodity_refined','specialty','premium_lunar','propellant_logistics','compute_digital','services','diversified','prior_manufacturing_mix']
labels=['Raw-price route proxy','Commodity stock','Specialty products','Premium lunar goods','Propellant + support','Computing + digital','Services + domestic support','Diversified','Earlier manufacturing mix']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'strategy-returns-v1'})
fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
a=axes[0,0]
for market,offset,color,label in [('global',-.18,'#397c9f','Middle service market'),('breakthrough',.18,'#31926f','Strong service market')]:
    values=[next(x['NPV_USD2041']['0.03']/1e12 for x in r['horizons'] if x['strategy']==n and x['market']==market and x['horizon_elapsed_years']==100) for n in names]
    a.barh([i+offset for i in range(len(names))],values,height=.34,label=label,color=color)
a.set(yticks=range(len(names)),yticklabels=labels,xlabel='100-year NPV at 3%, trillion USD2041',title='Whole enterprise: services can finance weak material exports')
a.set_ylim(len(names)-.2,-1.6);a.axvline(0,color='gray',lw=.7);a.legend(frameon=False,fontsize=9,loc='upper right')
a=axes[0,1]
base=next(x['operation']['available_USD2041'] for x in r['steady_comparisons'] if x['strategy']=='services' and x['market']=='global')
values=[(next(x['operation']['available_USD2041'] for x in r['steady_comparisons'] if x['strategy']==n and x['market']=='global')-base)/1e9 for n in names]
a.barh(labels,values,color=['#31926f' if x>=0 else '#ad514d' for x in values]);a.invert_yaxis();a.axvline(0,color='gray',lw=.7)
a.set(xlabel='Annual improvement, billion USD2041',title='Same founding system and middle service market\nCompared with services and domestic support')
a=axes[1,0]
for name,label,color in [('premium_lunar','Premium lunar goods','#94624a'),('propellant_logistics','Propellant + support','#31926f'),('diversified','Diversified','#397c9f')]:
    rows=[x for x in r['horizons'] if x['strategy']==name and x['market']=='breakthrough']
    a.plot([x['horizon_elapsed_years'] for x in rows],[x['NPV_USD2041']['0.03']/1e12 for x in rows],marker='o',label=label,color=color)
a.set(xlabel='Elapsed years from first investment',ylabel='NPV at 3%, trillion USD2041',title='Strong service market across appraisal horizons');a.legend(frameon=False)
a=axes[1,1]
selected=['premium_lunar','propellant_logistics','diversified','prior_manufacturing_mix']
values=[next(x['available_USD2041_year']/1e12 for x in r['cap_operation_checks'] if x['strategy']==n and x['market']=='breakthrough') for n in selected]
a.barh(['Premium lunar goods','Propellant + support','Diversified','Earlier manufacturing mix'],values,color=['#31926f' if x>=0 else '#ad514d' for x in values]);a.invert_yaxis();a.axvline(0,color='gray',lw=.7)
a.set(xlabel='Annual available cash, trillion USD2041',title='Imposed enclosed-cap operation, strong market\nA separate stock test; expansion need not reach it')
fig.suptitle('Industrial exports and worldwide services share customers, capital and capacity\nConditional USD2041 scenarios; premium qualification, small modules and future markets remain unvalidated',fontsize=14)
fig.savefig(OUTPUT/'strategies.png',dpi=150,metadata={'Description':json.dumps(identity)})
fig.savefig(OUTPUT/'strategies.svg',metadata={'Date':None,'Description':json.dumps(identity)})
(OUTPUT/'strategies_manifest.json').write_text(json.dumps(identity,indent=2)+'\n')
