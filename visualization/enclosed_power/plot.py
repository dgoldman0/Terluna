"""Render the enclosed surface-power product without adding a physical model."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT/'provisioning/enclosed_power/results/power_comparison.json'
NAMES = ['fission','beamed_legacy','solar_thermal','hybrid_50','solar_RFC','solar_battery']
LABELS = ['Fission','Beamed solar','Solar thermal + local store','Half fission + solar / fuel cells','Solar + regenerative fuel cells','Solar + batteries']


def selected_cost(case, imported):
    return next(x for x in case['costs'] if x['imported_mass_fraction']==imported and x.get('freight_case','reference')=='reference')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=ROOT/'research/runs/enclosed_power')
    args=ap.parse_args(); d=json.loads(PRODUCT.read_text())
    cases={x['architecture']:x for x in d['comparisons'] if x['technology_case']=='reference_assumptions'}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(15,10.8))
    colours=['#286f91','#c47438','#58986e']
    y=np.arange(len(NAMES))
    for imp,offset,colour,label in [(1,-.17,colours[0],'Imported equipment / fluid'),(.1,.17,colours[1],'90% local equipment / fluid')]:
        values=[selected_cost(cases[n],imp)['installed_supply_and_delivery_USD2026_per_mean_kW']/1000 for n in NAMES]
        axes[0,0].barh(y+offset,values,.32,color=colour,label=label)
        for yy,v in zip(y+offset,values):axes[0,0].text(v*1.06,yy,f'{v:,.1f}',va='center',fontsize=8)
    axes[0,0].set(yticks=y,yticklabels=LABELS,xscale='log',xlim=(10,7000),xlabel='Initial equipment + delivery (million USD2026 / MW mean)',title='Middle assumptions favor fission')
    axes[0,0].set_ylim(6.6,-.6);axes[0,0].legend(frameon=False,loc='lower left',fontsize=9)
    for j,(bundle,label) in enumerate([('development_target','Development target'),('reference_assumptions','Middle'),('high_burden','High burden')]):
        byname={x['architecture']:x for x in d['comparisons'] if x['technology_case']==bundle}
        vals=[selected_cost(byname[n],1)['lifecycle']['100']['levelized_USD2026_MWh'] for n in NAMES]
        axes[0,1].barh(y+(j-1)*.23,vals,.22,color=colours[j],label=label)
    axes[0,1].set(yticks=y,yticklabels=LABELS,xscale='log',xlabel='100-year levelized cost (USD2026 / MWh)',title='Development assumptions change the ordering')
    axes[0,1].set_ylim(6.3,-.6);axes[0,1].legend(frameon=False,loc='lower left',fontsize=8,ncol=3)
    axes[0,1].text(.97,.97,'Imported equipment / fluid\nThermal bulk material remains local',transform=axes[0,1].transAxes,ha='right',va='top',fontsize=9)
    for imp,colour,label in [(1,colours[0],'Imported equipment / fluid'),(.1,colours[1],'90% local equipment / fluid')]:
        rows=d['beam_cost_curves']
        axes[1,0].plot([r['useful_fraction']*100 for r in rows],[selected_cost(r,imp)['installed_supply_and_delivery_USD2026_per_mean_kW']/1000 for r in rows],'-',color=colour,label=label+' — beam')
        nf=selected_cost(cases['fission'],imp)['installed_supply_and_delivery_USD2026_per_mean_kW']/1000
        axes[1,0].axhline(nf,color=colour,ls='--',alpha=.8)
        point=next(x for x in d['initial_cost_crossovers'] if x['imported_mass_fraction']==imp)
        x=100*point['beam_useful_capacity_fraction_to_match_reference_fission']
        axes[1,0].scatter([x],[nf],color=colour,marker='x',s=75)
        axes[1,0].annotate(f'{x:.0f}% crossover',(x,nf),xytext=(4,9),textcoords='offset points',color=colour)
    axes[1,0].set(xlabel='Useful constellation capacity (%)',ylabel='Initial cost (million USD2026 / MW mean)',title='Better beam coverage could reverse the choice',ylim=(0,190))
    axes[1,0].legend(frameon=False,loc='upper right',fontsize=8)
    axes[1,0].text(.04,.04,'Dashed: middle fission cost. Curves use fixed hardware.\nCrossover marks are calculated; coverage is unqualified.',transform=axes[1,0].transAxes,fontsize=8)
    for name,colour,label in [('solar_RFC',colours[0],'Solar + regenerative fuel cells'),('solar_thermal',colours[1],'Solar thermal + local store')]:
        rows=[r for r in d['load_and_site_sensitivities'] if r['site']=='ordinary' and r['architecture']==name]
        axes[1,1].plot([100*r['night_load_shift'] for r in rows],[selected_cost(r,1)['installed_supply_and_delivery_USD2026_per_mean_kW']/1000 for r in rows],'o-',color=colour,label=label)
    axes[1,1].set(xlabel='Dark-time load transferred to daylight (%)',ylabel='Initial cost (million USD2026 / MW mean)',title='Scheduling reduces storage but adds factory demands')
    axes[1,1].legend(frameon=False,loc='upper right',fontsize=9)
    axes[1,1].text(.04,.05,'Imported hardware; extra factory capacity is unpriced.\nCurrent processing is only ~15% of year-200 surface demand.',transform=axes[1,1].transAxes,fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.15,axis='x' if ax in axes[0,:] else 'both')
    fig.suptitle('Power for the enclosed programme before Open Moon',fontsize=18,x=.04,ha='left')
    fig.text(.04,.025,'Conditional engineering comparisons, not forecasts or revised financial returns. Reference installation: 1 GW mean, 1.5 GW specified capacity.\nOrdinary lunar night: 354 hours. Surface / orbital delivery: 1,000 / 500 USD2026 per kg. No shield credit. Fixed development and site works unpriced.',fontsize=9,color='#444444')
    fig.tight_layout(rect=(.02,.085,.99,.95),w_pad=3,h_pad=3)
    args.out.mkdir(parents=True,exist_ok=True)
    for suffix in ['png','pdf']:fig.savefig(args.out/f'power-comparison.{suffix}',dpi=160)
    (args.out/'render-provenance.json').write_text(json.dumps({'schema':'terluna.visualization.enclosed-power/1','source_product':str(PRODUCT.relative_to(ROOT)),'source_sha256':hashlib.sha256(PRODUCT.read_bytes()).hexdigest(),'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
    print(args.out/'power-comparison.png')


if __name__=='__main__':main()
