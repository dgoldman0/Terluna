"""Display the saved stratification and intervention accounts without modeling."""
from pathlib import Path
import argparse
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
PORTFOLIOS=['services','diversified','propellant_logistics','prior_manufacturing_mix']
NAMES=['Worldwide\nservices','Diversified','Propellant /\nlogistics','Earlier\nmanufacturing']
VARIANTS=['baseline','demand_25','capture_8','service_rollout','combined_rollout','demand_and_rollout','tug_life_30']
LABELS=['Baseline','25% stronger demand','Faster industrial capture','Faster service rollout','Both rollout changes','Demand + both rollout changes','Longer source-tug active life']


def plotted_data(d):
    c=next(c for c in d['baseline_stratification'] if c['name']=='propellant_logistics:breakthrough:wider')
    ps=[p for p in c['periods'] if p['last_year']<=200]
    costs={k:sum(p['values']['0.03']['costs_by_category'][k] for p in ps) for k in ps[0]['values']['0.03']['costs_by_category']}
    gross=sum(sum(p['values']['0.03']['revenue_by_sector'].values()) for p in ps)
    physical=sum(v for k,v in costs.items() if k.startswith('operating_'))+costs['hotel_payroll_USD2041']+costs['residual_service_renewal_USD2041']
    capital=sum(v for k,v in costs.items() if k.startswith('capital_'))+costs['residual_service_capital_USD2041']
    bridge=c['financing_bridges']['200']
    waterfall=[gross,-costs['commercial_operating_USD2041'],-physical,-capital,-costs['seed_deployment_and_allowance_USD2041'],
               -bridge['zero_real_yield_asset_carry_cost_PV']-bridge['interest_minus_discount_cost_on_debt_PV'],
               -bridge['terminal_net_financial_claim_PV'],bridge['sponsor_distributed_cash_NPV']]
    by={(x['portfolio'],x['experiment']['variant']):x for x in d['interventions']}
    matrix=[]
    for variant in VARIANTS:
        matrix.append([next(r['checkpoint']['discounted']['0.03']['sponsor_distributed_cash_NPV'] for r in by[p,variant]['checkpoints'] if r['horizon']==200)/1e12
                       if by[p,variant]['status']=='conditional_account' else None for p in PORTFOLIOS])
    op={variant:[next(r['total_physical_operation'] for r in by[p,variant]['physical_checkpoints'] if r['year']==200)/1e12 for p in PORTFOLIOS] for variant in ('baseline','tug_life_30')}
    return dict(waterfall=[v/1e12 for v in waterfall],period_labels=[f"{p['first_year']}–{p['last_year']}" for p in ps],
                period_gross=[sum(p['values']['0.03']['revenue_by_sector'].values())/1e12 for p in ps],
                period_net=[p['values']['0.03']['net_real_cash']/1e12 for p in ps],intervention_NPV=matrix,physical_operation=op)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=ROOT/'research/studies/enclosed_returns/results/bottlenecks.json');parser.add_argument('--output',type=Path,default=ROOT/'research/runs/enclosed_returns/bottlenecks')
    a=parser.parse_args();d=json.loads(a.input.read_text());assert d['schema']=='terluna.humanities.enclosed-bottlenecks/1';v=plotted_data(d)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(17,12),gridspec_kw={'height_ratios':[1,1.15]})
    ax=axes[0,0];running=0
    for i,value in enumerate(v['waterfall']):
        base=0 if i in (0,7) else running;ax.bar(i,value,bottom=base,color='#218096' if i in (0,7) else '#a95a41')
        end=base+value;ax.text(i,max(base,end)+.65,f'{value:.2f}',ha='center',fontsize=9)
        if i<7:running=end
    ax.set_xticks(range(8),['Receipts','Commercial\nexpense','Physical +\nother operation','Later\ncapital','Seed +\nallowance','Financial\ncarry','Terminal\ncash held','Sponsor\nNPV'],rotation=22,ha='right')
    ax.set_ylabel('Trillion constant USD2041, present value');ax.set_ylim(0,49);ax.set_title('How 43.84T of discounted receipts becomes 11.94T NPV\nPropellant / logistics portfolio; 1% real demand growth')
    ax=axes[0,1];x=np.arange(len(v['period_labels']))
    ax.bar(x-.18,v['period_gross'],width=.36,label='Gross receipts',color='#218096');ax.bar(x+.18,v['period_net'],width=.36,label='Cash after operations and capital',color='#d09c3e')
    ax.set_xticks(x,v['period_labels']);ax.set_xlabel('Programme years');ax.set_ylabel('Trillion constant USD2041, present value');ax.set_title('Earlier receipts carry more weight at 3% real\nSame portfolio and growth as the left panel');ax.legend(fontsize=9);ax.axhline(0,color='#555',lw=.7)
    ax=axes[1,0];arr=np.array([[np.nan if x is None else x for x in row] for row in v['intervention_NPV']]);cm=plt.get_cmap('RdYlGn').copy();cm.set_bad('#e5e5e5')
    ax.imshow(arr,cmap=cm,vmin=-5,vmax=20,aspect='auto');ax.set_yticks(range(7),LABELS);ax.set_xticks(range(4),NAMES)
    for i,row in enumerate(v['intervention_NPV']):
        for j,value in enumerate(row):ax.text(j,i,'Staff limit\nexceeded' if value is None else f'{value:.2f}',ha='center',va='center',fontsize=9)
    ax.set_title('Year-200 sponsor NPV under separate interventions\n0.5% demand growth; trillion USD2041 at 3% real')
    ax=axes[1,1];x=np.arange(4)
    ax.bar(x-.18,v['physical_operation']['baseline'],width=.36,color='#a95a41',label='20-year active source-tug life')
    ax.bar(x+.18,v['physical_operation']['tug_life_30'],width=.36,color='#218096',label='30-year life target')
    for j,p in enumerate(('baseline','tug_life_30')):
        for i,value in enumerate(v['physical_operation'][p]):ax.text(i+(-.18 if j==0 else .18),value+.12,f'{value:.2f}',ha='center',fontsize=9)
    ax.set_xticks(x,NAMES);ax.set_ylim(0,11.5);ax.set_ylabel('Trillion constant USD2041 per year');ax.set_title('Durability reduces the recurring physical bill\nYear 200; 0.5% demand growth; price unchanged');ax.legend(fontsize=9)
    for ax in (axes[0,0],axes[0,1],axes[1,1]):ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    fig.suptitle('What limits the returns from staged off-Earth development?',fontsize=19,y=.992)
    fig.text(.5,.95,'Highest market case: 2 billion potential audience; 35% capital allowance; U financing.\nDemand, rollout and component-life changes are conditional tests. Physical development and source access remain unqualified.',ha='center',va='top',fontsize=11)
    fig.text(.5,.015,'All four baseline controls reproduce exactly. Seven faster-rollout cases fail the workforce screen and receive no completed return estimate.\nThe longer-life case has no priced development premium. Values are financial returns; unpriced spillovers and future development options are separate.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.06,1,.925),h_pad=3,w_pad=3)
    a.output.mkdir(parents=True,exist_ok=True)
    for extension in ('png','svg'):fig.savefig(a.output/('bottlenecks.'+extension),dpi=150)
    (a.output/'manifest.json').write_text(json.dumps(dict(input=str(a.input.relative_to(ROOT)),input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),plotted_data=v),indent=2)+'\n')
    print(a.output)


if __name__=='__main__':main()
