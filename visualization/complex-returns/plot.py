"""Standalone figures from the reviewed returns product."""
import json
import hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'research/studies/complex_returns/results/returns.json'
OUTPUT=ROOT/'research/runs/complex_returns/figures'
OUTPUT.mkdir(parents=True,exist_ok=True)
identity={'product_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'plot_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
p=json.loads(SOURCE.read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'complex-returns-v1'})
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
a=axes[0,0]
for rate in (0,.005,.01):
 rows=p['growth_returns'][f'{rate}:None']['operating']
 a.plot([r['year'] for r in rows],[r['population']/1e6 for r in rows],label=f'{rate:.1%} initial growth')
a.set(xlabel='Years after founding occupation',ylabel='Residents, million',title='Physical growth tapers to the enclosed reference')
a.legend(frameon=False)
a=axes[0,1]
for name,label in [('specialty','Mixed specialty'),('qualified_high_value','Same plant, higher prices'),('large','Larger manufacturing portfolio')]:
 rows=[r for r in p['static_returns'] if r['portfolio']==name and r['horizon_operating_years']==100]
 a.plot([r['real_rate']*100 for r in rows],[r['cash_NPV_USD2041']/1e12 for r in rows],marker='o',label=label)
a.axhline(0,color='gray',lw=.7);a.set(xlabel='Real discount rate, %',ylabel='External recovery NPV, trillion USD2041',title='100 operating years plus 33 construction draws')
a.legend(frameon=False,fontsize=8)
a=axes[1,0]
for demand,label in [('None','Demand scales with installed capacity'),('5','Earth sales capped at 5× founding sales')]:
 rows=p['growth_returns']['0.01:'+demand]['operating']
 a.plot([r['year'] for r in rows],[r['external_cash_USD2041']/1e12 for r in rows],label=label)
a.axhline(0,color='gray',lw=.7);a.set(xlabel='Years after founding occupation',ylabel='External balance, trillion USD2041/year',title='1% initial growth with independently constrained demand')
a.legend(frameon=False,fontsize=8)
a=axes[1,1]
rows=[r for r in p['stakeholder_value_sensitivity'] if r['extra_net_spillover_USD2041_year']==0]
a.bar([str(int(100*r['export_counterfactual_value_fraction']))+'%' for r in rows],[r['required_extra_annual_net_benefit_USD2041']/1e12 for r in rows],color='#547c96')
a.set(xlabel='Export purchase value credited economically',ylabel='Required extra net benefit, trillion USD2041/year',title='100-year socioeconomic switching values at 3%')
fig.suptitle('Founding off-Earth complex: conditional returns across methods\nNo array or Open Moon; constant 2041 dollars; demand and valuation are assumptions',fontsize=13)
fig.savefig(OUTPUT/'returns.png',dpi=160,metadata={'Description':json.dumps(identity)})
fig.savefig(OUTPUT/'returns.svg',metadata={'Date':None,'Description':json.dumps(identity)})
(OUTPUT/'manifest.json').write_text(json.dumps(identity,indent=2)+'\n')
