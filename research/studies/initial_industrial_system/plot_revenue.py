"""Render scenario accounts from the recorded product; no fitted market data."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
HERE=Path(__file__).resolve().parent


def main():
    r=json.loads((HERE/'results/revenue.json').read_text())
    cases=[r['cases']['local_structures_space_chemical'][m] for m in ['weak','reference','strong']]
    cases.append(r['cases']['favourable_earth_chemical']['strong'])
    labels=['Weak sales\nLocal tug structures,\nspace-supplied fuel','Reference sales\nSame supply case','Strong sales\nSame supply case','Strong sales\nFavourable original\nindustrial case']
    fig,ax=plt.subplots(figsize=(11,6.8),layout='constrained')
    fig.patch.set_facecolor('#faf9f5');ax.set_facecolor('#faf9f5')
    x=np.arange(4);bottom=np.zeros(4)
    for key,label,color in [('external_operating_cost_USD2041','Earth operating purchases','#648a9a'),('external_renewal_USD2041','Earth replacement purchases','#b9ced6'),('resource_transition_provision_USD2041','Resource provision','#d5b46b')]:
        vals=np.array([c[key]/1e9 for c in cases]);ax.bar(x,vals,bottom=bottom,label=label,color=color,width=.6);bottom+=vals
    receipts=[c['gross_receipts_USD2041']/1e9 for c in cases]
    ax.scatter(x,receipts,color='#19683f',s=90,marker='D',zorder=3,label='Gross receipts from Earth')
    for i,c in enumerate(cases):
        net=c['balance_after_renewal_and_resource_provision_USD2041']/1e9
        ax.text(i,bottom[i]+17,f'Net {net:,.1f}',ha='center',color='#96372d',fontsize=11)
        ax.annotate(f'{receipts[i]:.1f}',(i,receipts[i]),xytext=(13,5),textcoords='offset points',fontsize=10,color='#19683f')
    ax.set_xticks(x,labels);ax.set_ylabel('Billion 2041 USD per stabilized year');ax.set_ylim(0,880)
    ax.set_title('Export receipts and recurring external requirements',loc='left',fontweight='bold',fontsize=16,pad=18)
    ax.spines[['right','top']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    ax.legend(loc='upper right',frameon=False,fontsize=10)
    fig.text(.02,-.03,'Conditional post-occupation portfolios for 250,000 residents and supporting industry.\nDebt is separate; internal transactions cancel. Construction capacity is retained. No market forecast or adopted route.',fontsize=10)
    fig.savefig(HERE/'results/revenue.png',dpi=160,bbox_inches='tight');plt.close(fig)
if __name__=='__main__':main()
