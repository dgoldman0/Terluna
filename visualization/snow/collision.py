"""Plot the independently computed prescribed-bath collision distributions.

Reads frozen domain products only; no collision or geometry physics lives here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
PRODUCT='atmosphere/snow/results/collision_reference.json'
COLORS=['#bf5b1b','#824c9c','#247b83']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def curve(mass,counts,total,mean):
    m,n=np.asarray(mass,float),np.asarray(counts,float)
    if m.shape!=n.shape or not np.isfinite(m).all() or not np.isfinite(n).all() or (m<=0).any() or (n<0).any():
        raise ValueError('invalid stored distribution')
    if not np.isclose(n.sum(),total,rtol=2e-10,atol=0):raise ValueError('collector number normalization differs')
    if not np.isclose(m@n/total,mean,rtol=2e-10,atol=0):raise ValueError('stored collector mean differs')
    keep=n>0;m,n=m[keep],n[keep];order=np.argsort(m);m,n=m[order],n[order]
    return dict(mass_kg=m.tolist(),number_cdf=(np.cumsum(n)/total).tolist())


def prepare():
    data=json.loads((ROOT/PRODUCT).read_text())
    if data['schema']!='terluna.atmosphere.snow-collision-reference/1':raise ValueError('wrong schema')
    identities={PRODUCT:sha(ROOT/PRODUCT)}
    for group in ('files','inputs'):
        for path,expected in data['producer'][group].items():
            if sha(ROOT/path)!=expected:raise ValueError('producer/input changed: '+path)
            identities[path]=expected
    capture=data['capture'];ref=max(capture['continuous_mass_references'],key=lambda r:r['samples'])
    initial=capture['initial_collector_mass_kg'];panels=[];right=1.
    for index in (1,2):
        row=ref['snapshots'][index]
        reference=curve(row['masses_kg'],row['sample_counts'],ref['samples'],row['mean_mass_kg'])
        curves=[]
        for grid in capture['fixed_pivot']:
            r=grid['snapshots'][index]
            if r['time_s']!=row['time_s'] or abs(r['number_residual'])>1e-10 or r['upper_exit_number']>1e-10:
                raise ValueError('comparison has mismatched time or material number-domain loss')
            n=np.asarray(r['numbers_per_initial_collector'])
            c=curve(grid['mass_pivots_kg'],n,float(n.sum()),r['mass_kg_m3']/r['number_m3'])
            c['mass_bins']=grid['mass_bins'];curves.append(c)
        for c in [reference,*curves]:
            k=np.searchsorted(c['number_cdf'],.9999)
            right=max(right,c['mass_kg'][min(k,len(c['mass_kg'])-1)]/initial)
        panels.append(dict(time_s=row['time_s'],reference=reference,grids=curves))
    right=1.+(right-1.)*1.05
    for panel in panels:
        for c in [panel['reference'],*panel['grids']]:
            x=np.asarray(c['mass_kg'])/initial;cdf=np.asarray(c['number_cdf']);inside=np.flatnonzero(x<=right)
            c['number_fraction_above_display']=float(1-cdf[inside[-1]]) if len(inside) else 1.
            if c['number_fraction_above_display']>1e-4+1e-12:raise ValueError('display hides material tail')
    return data,identities,dict(initial_collector_mass_kg=initial,panels=panels,right_mass_ratio=right,
        reference_samples=ref['samples'],dkw_halfwidth=ref['number_cdf_uniform_95_percent_dkw_halfwidth'],
        interpretation='Number CDF per collector, normalized by all surviving collectors. Less than0.01% of each curve lies above the view; omitted tails remain in normalization. The band is the independent-reference95% DKW number-CDF band at each displayed time, not a joint two-time band or a mass-CDF uncertainty.')


def main():
    data,identities,plot=prepare();capture=data['capture']
    output=HERE/'generated';output.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.hashsalt':'terluna-snow-collision-v1'})
    fig,axes=plt.subplots(1,2,figsize=(12.4,6.0));fig.subplots_adjust(left=.08,right=.98,bottom=.22,top=.73,wspace=.24)
    fig.text(.08,.93,'Repeated small-particle capture',fontsize=21,weight='bold',color='#172b3a')
    fig.text(.08,.875,f'Imposed collector {capture["initial_collector_diameter_m"]*100:.2f} cm · captured particles {capture["captured_diameter_m"]*1000:.3f} mm · bath {capture["bath_number_m3"]/1000:g}/L',color='#52646d')
    for panel,ax in zip(plot['panels'],axes):
        for c,color in zip(panel['grids'],COLORS):
            x=(np.asarray(c['mass_kg'])/plot['initial_collector_mass_kg']-1)*100
            ax.step(np.r_[-1.,x],np.r_[0.,c['number_cdf']],where='post',color=color,lw=1.8,label=f'{c["mass_bins"]} mass pivots')
        c=panel['reference'];x=(np.asarray(c['mass_kg'])/plot['initial_collector_mass_kg']-1)*100;cdf=np.asarray(c['number_cdf'])
        x=np.r_[-1.,x,(plot['right_mass_ratio']-1)*100];cdf=np.r_[0.,cdf,1.]
        band=plot['dkw_halfwidth'];ax.fill_between(x,np.maximum(cdf-band,0),np.minimum(cdf+band,1),step='post',color='#172b3a',alpha=.13,label='95% reference CDF band')
        ax.step(x,cdf,where='post',color='#172b3a',lw=2,label=f'Exact events ({plot["reference_samples"]:,} histories)')
        ax.set(xlim=(-1.,(plot['right_mass_ratio']-1)*100),ylim=(-.015,1.015),xlabel='Collector mass increase (%)',ylabel='Cumulative fraction of collectors',title=f'{panel["time_s"]/60:.2f} minutes of capture')
        ax.yaxis.set_major_formatter(PercentFormatter(1.));ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#dce3e6',lw=.7);ax.set_axisbelow(True)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper left',bbox_to_anchor=(.073,.835),ncol=3,frameon=False,fontsize=9)
    fig.text(.08,.105,'Fixed −10°C air; replenished bath; no deposition, settling or breakup.',fontsize=10,color='#52646d')
    fig.text(.08,.063,'This isolates collision allocation around an imposed large collector. It does not predict a lunar flake size.',fontsize=10,color='#52646d')
    metadata={'input_sha256':identities[PRODUCT],'renderer_sha256':sha(__file__)}
    fig.savefig(output/'collision_comparison.png',dpi=160,metadata=metadata)
    fig.savefig(output/'collision_comparison.svg',metadata={'Date':None,'Description':json.dumps(metadata,sort_keys=True)})
    plt.close(fig)
    manifest=dict(schema='terluna.visualization.snow-collision/1',renderer_sha256=sha(__file__),inputs=identities,
                  upstream_evidence=data['evidence'],upstream_reading_rule=data['reading_rule'],plot=plot,
                  outputs={name:sha(output/name) for name in ('collision_comparison.png','collision_comparison.svg')})
    (output/'collision_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output),input_sha256=identities[PRODUCT],right_mass_ratio=plot['right_mass_ratio'])))


if __name__=='__main__':main()
