"""Display the admitted year-200 national-account distributions without refitting."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
LABELS={'raw_resource_price':'Raw-resource offer','commodity_refined':'Commodity refined','specialty':'Specialty products','premium_lunar':'Premium lunar products','propellant_logistics':'Propellant and logistics','compute_digital':'Computing and digital','services':'Worldwide services','diversified':'Diversified','prior_manufacturing_mix':'Earlier manufacturing mix'}

def main(source,out):
    source,out=Path(source),Path(out)
    raw=source.read_bytes()
    d=json.loads(gzip.decompress(raw))
    if d['schema']!='terluna.humanities.production-income/1':raise ValueError('Unexpected account product')
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    names=list(LABELS)
    fig,axs=plt.subplots(1,3,figsize=(14,6.2),sharey=True)
    drawn=[]
    for ax,metric,title in zip(axs,('gdp','outside_ownership','three_quarters_local'),('GDP','GNI: outside ownership','GNI: 75% locally owned')):
        for i,name in enumerate(names):
            for screen,offset,color,label in [('funded_completion',-.13,'#2166ac','Full build funded'),('previously_self_financing_paths',.13,'#b05c00','Previously self-financing paths')]:
                v=d['distributions']['balanced'][name][screen]
                z=v['gdp'] if metric=='gdp' else v['gni'][metric]
                qs=np.array([z['quantiles'][k] for k in ('0.1','0.5','0.9')])/1e12
                if not qs[0]<=qs[1]<=qs[2]:raise ValueError('Unordered interval')
                ax.plot([qs[0],qs[2]],[i+offset]*2,lw=3,color=color,alpha=.85,label=label if i==0 and metric=='gdp' else None)
                ax.plot(qs[1],i+offset,'o',color=color,ms=4)
                drawn.append(dict(portfolio=name,metric=metric,screen=screen,p10_p50_p90_trillion=qs.tolist()))
        ax.set_title(title,fontweight='bold')
        ax.set_xlabel('Trillion USD2041 equivalents / year')
        ax.grid(axis='x',alpha=.2)
        ax.set_xlim(0,24)
    axs[0].set_yticks(range(len(names)),[LABELS[n] for n in names]);axs[0].invert_yaxis()
    fig.suptitle('Year 200: output and income across the nine portfolios',fontweight='bold',y=.99)
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,bbox_to_anchor=(.59,.045),frameon=False)
    fig.text(.02,.015,'Dots: medians. Lines: central 80%. Subjective balanced prior; completion assumed. Source self-financing histories can incur losses under new draws.',fontsize=9)
    fig.tight_layout(rect=(0,.105,1,.95))
    fig.savefig(out/'intervals.png',dpi=170);fig.savefig(out/'intervals.svg');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,5.4))
    for ax,key,title in zip(axs,('gdp','three_quarters_local'),('Diversified portfolio: GDP','Diversified portfolio: GNI, 75% local')):
        for prior,color in [('restrained','#5f7184'),('balanced','#2166ac'),('expansive','#b05c00')]:
            r=d['distributions'][prior]['diversified']['funded_completion']
            z=r['gdp'] if key=='gdp' else r['gni'][key]
            curve=z['inverse_cdf']
            ax.plot(np.array(curve['value'])/1e12,curve['probability'],label=prior.capitalize(),color=color,lw=2)
        ax.set_title(title,fontweight='bold');ax.set_xlabel('Trillion USD2041 equivalents / year')
        ax.set_ylabel('Cumulative probability');ax.grid(alpha=.2);ax.set_ylim(0,1);ax.legend(frameon=False,loc='lower right')
    fig.suptitle('The distribution changes substantially with the market prior',fontweight='bold',y=.99)
    fig.text(.02,.02,'Completed, funded build; subjective probabilities within admitted market cases. No probability of completion or unmodelled upside is assigned.',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.94));fig.savefig(out/'prior_sensitivity.png',dpi=170);fig.savefig(out/'prior_sensitivity.svg');plt.close(fig)
    (out/'manifest.json').write_text(json.dumps(dict(source=str(source),source_sha256=hashlib.sha256(raw).hexdigest(),renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),displayed_intervals=drawn,evidence=d['evidence']),indent=2)+'\n')
    print(out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',default=str(ROOT/'research/studies/production_income/results/year200.json.gz'));p.add_argument('--out',default=str(ROOT/'research/runs/production_income/figures'));a=p.parse_args();main(a.input,a.out)
