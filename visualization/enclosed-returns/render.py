"""Render the admitted return accounts without recalculating their economics."""
from pathlib import Path
import argparse
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
LABELS={'raw_resource_price':'Raw-price offer','commodity_refined':'Commodity offer','specialty':'Specialty products','premium_lunar':'Premium lunar','propellant_logistics':'Propellant / logistics','compute_digital':'Computing / digital','services':'Worldwide services','diversified':'Diversified','prior_manufacturing_mix':'Prior manufacturing mix'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=ROOT/'research/studies/enclosed_returns/results/returns.json');p.add_argument('--output',type=Path,default=ROOT/'research/runs/enclosed_returns/figures');a=p.parse_args()
    d=json.loads(a.input.read_text());assert d['schema']=='terluna.humanities.enclosed-returns/1'
    rows=[r for r in d['checkpoints'] if r['mode']=='U_credit' and r['market']=='breakthrough' and len(r['name'].split(':'))==3]
    by={(r['portfolio'],r['growth'],r['horizon']):r for r in rows};names=list(LABELS);y=np.arange(len(names));series={}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(14,11));colors=['#2580a0','#a64b31','#7e9c38']
    for h,shift,color in [(100,-.16,colors[0]),(200,.16,colors[1])]:
        values=[by[n,'reference',h]['annual_operating_available']/1e12 for n in names];series['operating_'+str(h)]=values
        ax[0,0].barh(y+shift,values,height=.3,color=color,label=f'Year {h}')
        values=[by[n,'reference',h]['discounted']['0.03']['sponsor_distributed_cash_NPV']/1e12 for n in names];series['sponsor_NPV_'+str(h)]=values
        ax[0,1].plot(values,y+shift,'o',color=color,label=f'Through year {h}')
    for j,(growth,label) in enumerate([('static','0%'),('reference','0.5%'),('wider','1%')]):
        values=[by[n,growth,200]['additional_equity_required']/1e12 for n in names];series['extra_equity_'+growth]=values
        ax[1,0].barh(y+(j-1)*.23,values,height=.22,color=colors[j],label=label+' real demand growth')
    for axis in (ax[0,0],ax[0,1],ax[1,0]):
        axis.set_yticks(y,[LABELS[n] for n in names]);axis.invert_yaxis();axis.axvline(0,color='#444444',linewidth=.8);axis.grid(axis='x',alpha=.18);axis.set_axisbelow(True);axis.set_xlabel('Trillion constant USD2041');axis.legend(fontsize=9)
    ax[0,0].set_title('Annual operating cash after renewal');ax[0,1].set_title('Pooled funder cash NPV at 3% real\nDistributions less contributions; no terminal claims')
    ax[1,0].set_title('Additional equity required through year 200\nBeyond initial seed investment')
    selected={'services':'#2580a0','diversified':'#a64b31','propellant_logistics':'#7e9c38','premium_lunar':'#8064a2'}
    for portfolio,color in selected.items():
        path=next(p for p in d['reference_finance_paths'] if p['mode']=='U_credit' and p['name']==portfolio+':breakthrough:reference')
        years=[r['year'] for r in path['annual']];values=[r['operating_available']/1e12 for r in path['annual']]
        ax[1,1].plot(years,values,label=LABELS[portfolio],color=color);series['annual_'+portfolio]={'years':years,'cash':values}
    ax[1,1].axhline(0,color='#444444',linewidth=.8);ax[1,1].axvline(100,color='#aaaaaa',linestyle=':');ax[1,1].axvline(200,color='#aaaaaa',linestyle=':')
    ax[1,1].set_title('Operating cash during construction and at the cap');ax[1,1].set_xlabel('Elapsed programme year');ax[1,1].set_ylabel('Trillion constant USD2041 per year');ax[1,1].legend(fontsize=9);ax[1,1].grid(alpha=.18)
    fig.suptitle('Staged construction and worldwide exports',fontsize=18)
    fig.text(.5,.935,'Highest market case (2 billion potential audience); U financing; 35% capital-cost allowance\n0.5% real demand growth unless shown otherwise',ha='center',fontsize=10)
    fig.text(.5,.012,'Conditional appraisals. Source access, trajectories and complete deployment remain unqualified.\nRaw and commodity Earth-return material offers fail the price floor; other receipts remain in those portfolios.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.06,1,.925));a.output.mkdir(parents=True,exist_ok=True)
    for extension in ('png','svg'):fig.savefig(a.output/('returns.'+extension),dpi=160)
    manifest={'input':str(a.input.relative_to(ROOT)),'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'portfolios':names,'plotted_data':series,'evidence':d['evidence']}
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(a.output)

if __name__=='__main__':main()
