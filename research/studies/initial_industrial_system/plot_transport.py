"""Standalone comparison figure; all values come from the transport product."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent

def main():
    r=json.loads((HERE/'results/transport.json').read_text())
    names=['earth_chemical','split_chemical','space_chemical','electric_earth_leg','electric_both_legs']
    labels=['Earth chemical fuel\n(original case)','Depot fuels return\n(chemical Earth leg)','Space chemical fuel\n(Earth leg)','Electric Earth leg\n(space-derived water)','Both long-haul legs\nelectric, space water']
    fig,axes=plt.subplots(1,3,figsize=(14.5,5.7),layout='constrained')
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    y=list(range(len(names)))
    for shift,phase,label,color in [(-.18,'peak_flows','Construction','#245d80'),(.18,'steady_flows','Sustaining','#62a7a1')]:
        axes[0].barh([v+shift for v in y],[r['cases'][n]['launch_split'][phase]['launches_at_100t']/1000 for n in names],height=.33,label=label,color=color)
        axes[2].barh([v+shift for v in y],[r['cases'][n][phase]['source_ore_kg_year']/1e9 for n in names],height=.33,color=color)
    axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].set_xlabel('Thousand Earth launches/year\n100 tonnes delivered per launch');axes[0].legend(frameon=False)
    axes[1].barh(y,[r['cases'][n]['cost']['industrial_and_initial_material_campaign_USD2041']/1e12 for n in names],color='#997543',height=.6)
    axes[1].set_xlabel('Industrial/material campaign\nTrillion 2041 USD')
    axes[2].set_xlabel('Source material processed\nMillion tonnes/year')
    for ax in axes[1:]:ax.set_yticks(y,[]);ax.invert_yaxis()
    fig.suptitle('Earth launch demand changes with freight architecture',fontsize=16)
    fig.supxlabel('Conditional sizing: fixed settlement demand; source, trajectory and technology qualification remain open.\nElectric cases use 2,000 s total-water Isp, 30% bus-to-jet efficiency and 100 W/kg solar assemblies.',fontsize=10)
    fig.savefig(HERE/'results/transport.png',dpi=170)
    plt.close(fig)
if __name__=='__main__':main()
