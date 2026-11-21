"""Standalone scientific overview of the capacity and accounting scenarios."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
HERE=Path(__file__).resolve().parent

def main():
    r=json.loads((HERE/'results/economy.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(14,9.6),layout='constrained')
    fig.suptitle('Founding economy: co-products, trade and capacity-supported growth',fontsize=18,fontweight='bold')
    ax=axes[0,0];co=r['coproducts'];vals=[co['products']['aggregate_recoverable_metals']['unused_recovery_potential_kg_year']/1e9,co['current_extra_metal_refining_capacity_kg_year']/1e9,1]
    bars=ax.barh(['Recoverable metal potential','Spare installed processing','Tested additional refinery'],vals,color=['#949aa2','#2e847c','#a66b2d'])
    ax.set_ylim(3.2,-.6);ax.set_xlim(0,7.3);ax.set_xlabel('Million tonnes per year');ax.set_title('Metal content exceeds present refining capacity',loc='left',fontweight='bold')
    for b,v in zip(bars,vals):ax.text(v+.07,b.get_y()+b.get_height()/2,f'{v:.3f}',va='center')
    ax.text(.01,.02,'No-growth local-structure case. Potential is unseparated content.\nAdditional refinery requires investment; it is not installed capacity.',transform=ax.transAxes,fontsize=9)
    ax=axes[0,1];names=['Mixed specialty','Mixed high prices','Larger high-value'];cases=[r['portfolios']['specialty'],r['portfolios']['qualified_high_value'],r['high_value_capacity_case']]
    x=np.arange(3);w=.25
    for off,key,label,color in [(-w,'gross_external_sales','Earth sales','#2e847c'),(0,'operating_cash_surplus','After operations','#6689ad'),(w,'after_renewal_and_resource','After renewal + resources','#a66b2d')]:
        ax.bar(x+off,[a['external']['USD2041'][key]/1e9 for a in cases],w,label=label,color=color)
    ax.axhline(0,color='black',lw=.7);ax.set_xticks(x,names,fontsize=9);ax.set_ylabel('Billion USD2041 per year');ax.set_title('Earth cash depends on products and buyers',loc='left',fontweight='bold');ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,.9));ax.set_ylim(-360,1040)
    ax.text(.02,.95,'Conditional prices and orders; debt excluded.',transform=ax.transAxes,fontsize=9)
    ax=axes[1,0]
    for key,label,color in [('None','Capacity-derived','#2e847c'),('0.03','3% comparison','#a66b2d'),('0.01','1% comparison','#6689ad'),('0','No growth','#777777')]:
        rows=r['growth'][key]['years'];ax.plot([x['year_after_full_occupation'] for x in rows],[x['occupied_capacity_people']/1000 for x in rows],label=label,color=color,lw=2,ls='--' if key=='0.03' else '-')
    ax.set_title('Growth using retained construction capacity',loc='left',fontweight='bold');ax.set_xlabel('Years after initial full occupation');ax.set_ylabel('Supported population capacity, thousands');ax.legend(fontsize=9)
    ax.text(.02,.97,'Two-year commissioning; no trade portfolio.\nFunding, approved sites and migration remain separate.',transform=ax.transAxes,va='top',fontsize=9)
    ax=axes[1,1];d=r['portfolios']['specialty']['domestic']['USD2041'];vals=[d[k]/1e9 for k in ['gross_domestic_output','GDP_GVA','exports','imports']]
    bars=ax.bar(['Gross output','Value added','Earth sales','Earth imports'],vals,color=['#6689ad','#2e847c','#a66b2d','#96969e']);ax.set_ylim(0,820);ax.set_ylabel('Billion USD2041 per year');ax.set_title('Domestic output, value added and Earth trade',loc='left',fontweight='bold');ax.tick_params(axis='x',labelsize=9)
    for b,v in zip(bars,vals):ax.text(b.get_x()+b.get_width()/2,v+10,f'{v:.1f}',ha='center',fontsize=10)
    ax.text(.02,.97,'Mixed specialty portfolio; illustrative domestic prices.\nGross output includes intermediate production.',transform=ax.transAxes,va='top',fontsize=9)
    for ax in axes.ravel():ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    fig.savefig(HERE/'results/economy.png',dpi=155)
    fig.savefig(HERE/'results/economy.pdf')
    plt.close(fig)

if __name__=='__main__':main()
