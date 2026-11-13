"""Render the reviewed staged-finance product; no model calculations."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'research/studies/complex_returns/results/staged.json'
OUTPUT=ROOT/'research/runs/complex_returns/figures'
OUTPUT.mkdir(parents=True,exist_ok=True)
identity={'product_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
p=json.loads(SOURCE.read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.hashsalt':'staged-complex-v1'})
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
colors={'limited':'#ad514d','global':'#397c9f','breakthrough':'#31926f'}
for name,path in p['paths'].items():
    rows=path['annual'];short=[r for r in rows if r['year']<=100]
    axes[0,0].plot([r['year'] for r in short],[r['services']['revenue_USD2041']/1e9 for r in short],label=name.capitalize(),color=colors[name])
    axes[0,1].plot([r['year'] for r in rows],[r['installed_scale'] for r in rows],label=name.capitalize(),color=colors[name])
axes[0,0].set(xlabel='Elapsed years from first investment',ylabel='New-service invoices, billion USD2041/year',title='Worldwide customers reach a small productive complex')
axes[0,0].legend(frameon=False)
axes[0,1].set(xlabel='Elapsed years from first investment',ylabel='Installed capacity / full founding system',yscale='log',title='Earned cash funds slow physical expansion')
axes[0,1].axhline(1,color='gray',lw=.8,ls='--')
axes[0,1].text(420,1.15,'Full founding capacity',color='#555555',fontsize=9)
axes[0,1].legend(frameon=False,loc='upper left')
a=axes[1,0]
for name in ('global','breakthrough'):
    rows=[r for r in p['horizons'] if r['scenario']==name]
    a.plot([r['horizon_elapsed_years'] for r in rows],[r['NPV_USD2041']['0.03']/1e12 for r in rows],marker='o',label=name.capitalize(),color=colors[name])
a.axhline(0,color='gray',lw=.8)
a.set(xlabel='Appraisal endpoint, elapsed year',ylabel='Capital recovery NPV, trillion USD2041',title='3% real; ongoing liquidity and renewal remain funded')
a.legend(frameon=False)
a=axes[1,1]
variants=['dividends_without_expansion','no_additional_B2B_demand','half_service_prices','higher_variable_costs','lower_collection']
labels=['No expansion','No new institutional demand','Half service prices','80% variable costs','70% cash collection']
values=[next(r['NPV_USD2041']['0.03']/1e12 for r in p['sensitivities'] if r['variant']==v and r['scenario']=='breakthrough') for v in variants]
a.barh(labels,values,color=[colors['breakthrough'] if x>=0 else colors['limited'] for x in values])
a.axvline(0,color='gray',lw=.8);a.invert_yaxis()
a.set(xlabel='800-year NPV at 3%, trillion USD2041',title='Strong demand still needs prices, margins and rights')
fig.suptitle('Self-financed expansion after a $2.789 trillion external seed\nConditional enclosed-complex scenarios; first receipts in year 8; small-module engineering unqualified',fontsize=13)
fig.savefig(OUTPUT/'staged.png',dpi=160,metadata={'Description':json.dumps(identity)})
fig.savefig(OUTPUT/'staged.svg',metadata={'Date':None,'Description':json.dumps(identity)})
(OUTPUT/'staged_manifest.json').write_text(json.dumps(identity,indent=2)+'\n')
