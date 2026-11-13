"""Display revised returns and ownership accounts from their pinned products."""
from pathlib import Path
import gzip,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
LABELS={'raw_resource_price':'Raw-price offer','commodity_refined':'Commodity offer','specialty':'Specialty products','premium_lunar':'Premium lunar','propellant_logistics':'Propellant / logistics','compute_digital':'Computing / digital','services':'Worldwide services','diversified':'Diversified','prior_manufacturing_mix':'Earlier manufacturing mix'}
def load(path):return json.loads(gzip.decompress(path.read_bytes()))
def main():
    source=ROOT/'research/studies/enclosed_returns/results/revised.json.gz'
    flows=ROOT/'research/studies/enclosed_returns/results/revised_flows.json.gz'
    d=load(source);f=load(flows)
    assert d['schema']=='terluna.humanities.revised-enclosed-returns/1'
    by={(r['portfolio'],r['growth'],r['variant'],v['horizon']):v for r in d['primary_cases'] if r.get('market')=='breakthrough' and r['status']=='conditional_account' for v in r['checkpoints'] if v['mode']=='U_credit'}
    names=list(LABELS);y=np.arange(len(names));series={}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(15,11))
    colors={'baseline':'#929aa0','combined':'#087f8c'}
    for variant,shift in [('baseline',-.16),('combined',.16)]:
        op=[by[n,'reference',variant,200]['annual_operating_available']/1e12 for n in names]
        npv=[by[n,'reference',variant,200]['discounted']['0.03']['sponsor_distributed_cash_NPV']/1e12 for n in names]
        series[variant]={'annual_operating':op,'funder_NPV':npv}
        axes[0,0].barh(y+shift,op,height=.3,color=colors[variant],label='Previous' if variant=='baseline' else 'Recovery + refurbishment')
        axes[0,1].plot(npv,y+shift,'o',color=colors[variant],label='Previous' if variant=='baseline' else 'Recovery + refurbishment')
    for ax in axes[0]:
        ax.set_yticks(y,[LABELS[n] for n in names]);ax.invert_yaxis();ax.axvline(0,color='#444',lw=.8);ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True);ax.legend(fontsize=9)
    axes[0,0].set_title('Year-200 annual operating cash\n0.5% annual real demand growth');axes[0,0].set_xlabel('Trillion USD2041 per year')
    axes[0,1].set_title('Year-200 pooled funder cash NPV at 3% real\n0.5% annual real demand growth');axes[0,1].set_xlabel('Trillion USD2041')
    selected={'services':'#936b32','propellant_logistics':'#087f8c','diversified':'#8756a2'}
    for portfolio,color in selected.items():
        for growth,style in [('reference','-'),('wider','--')]:
            case=next(c for c in f['cases'] if c['portfolio']==portfolio and c['growth']==growth)
            rows=[r for r in case['annual'] if r['year']>=8]
            years=[r['year'] for r in rows];values=[r['operating_available']/1e12 for r in rows]
            series[portfolio+'_'+growth]={'years':years,'annual_operating':values}
            axes[1,0].plot(years,values,style,color=color,label=LABELS[portfolio]+(' (0.5%)' if growth=='reference' else ' (1%)'))
    axes[1,0].set_title('Operating cash through buildout and sustaining years');axes[1,0].set_ylabel('Trillion USD2041 per year');axes[1,0].set_xlabel('Elapsed programme year');axes[1,0].axhline(0,color='#444',lw=.8);axes[1,0].axvline(200,color='#999',ls=':');axes[1,0].grid(alpha=.2);axes[1,0].legend(fontsize=8)
    horizons=[100,200,220];x=np.arange(3)
    for j,(portfolio,color) in enumerate(selected.items()):
        values=[by[portfolio,'wider','combined',h]['owner_payouts']/1e12 for h in horizons]
        series[portfolio+'_payouts_1percent']=values
        axes[1,1].bar(x+(j-1)*.24,values,width=.22,color=color,label=LABELS[portfolio])
    axes[1,1].set_title('Cumulative pooled funder payouts\n1% annual real demand growth; undiscounted');axes[1,1].set_xticks(x,[f'Year {h}' for h in horizons]);axes[1,1].set_ylabel('Trillion USD2041');axes[1,1].grid(axis='y',alpha=.2);axes[1,1].set_axisbelow(True);axes[1,1].legend(fontsize=9)
    fig.suptitle('Revised returns from recovery and equipment refurbishment',fontsize=17)
    fig.text(.5,.935,'Breakthrough worldwide demand; same 6.47-million-person capacity at year 200; U financing; 35% capital allowance',ha='center',fontsize=10)
    fig.text(.5,.014,'Conditional physical and commercial scenarios. Geological stocks do not establish qualified mines or delivery routes.\nOperating cash precedes new capital and debt service. Funder cash NPV and cumulative payouts are not system net wealth.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.065,1,.925));out=ROOT/'research/runs/enclosed_returns/revised_figures';out.mkdir(parents=True,exist_ok=True)
    for ext in ('png','svg'):fig.savefig(out/('revised_returns.'+ext),dpi=160)
    manifest={'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,flows)},'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'portfolios':names,'plotted_data':series,'evidence':d['evidence']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(out)
if __name__=='__main__':main()
