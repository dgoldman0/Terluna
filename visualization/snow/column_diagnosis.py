"""Display observed vapour losses and controlled process-removal results."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
PRODUCT='atmosphere/snow/results/column_diagnosis.json'
FOLLOWUP='atmosphere/snow/results/column_diagnosis_aggregation_only.json'
OUT=Path(__file__).resolve().parent/'generated'
FAMILIES=['full','no_deposition','no_aggregation','no_chipping','aggregation_only']
LABELS=['Complete\ncolumn','No vapour\ngrowth','No\ncollisions','No\nchipping','Sticking\nonly']


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(p,expected_cases=40):
    if p['schema']!='terluna.atmosphere.snow-column-diagnosis-controls/1':raise ValueError('Unknown schema')
    if not p['checks']['complete'] or p['completed_cases']!=expected_cases:raise ValueError('Incomplete diagnostic matrix')
    verified={}
    for category in ('files','inputs'):
        for name,sha in p['producer'][category].items():
            if digest(ROOT/name)!=sha:raise ValueError('Changed producer or input: '+name)
            verified[name]=sha
    if digest(ROOT/p['baseline']['path'])!=p['baseline']['sha256']:raise ValueError('Changed baseline')
    for name,g in p['groups'].items():
        rows=[r['result'] for r in p['cases'] if r['group']==name and r['status']=='complete']
        if len(rows)!=4:raise ValueError('Missing realization')
        traces=np.array([r['process_diagnostics']['blocks'] for r in rows]);gather=traces[:,1:].sum(axis=1)
        empty=gather[:,:,1]/gather[:,:,0]
        if not np.allclose(empty.mean(axis=0),g['mean_process_profiles']['empty_fraction'],rtol=1e-12,atol=1e-14):raise ValueError('Altered empty fractions')
        for j,key in [(5,'gather_deposition_mm_h'),(6,'gather_empty_unused_mm_h'),(7,'gather_occupied_unused_mm_h')]:
            values=np.array([gather[i,:,j].sum()/r['settings']['area_m2']/(r['settings']['gathering_time_s']/3600) for i,r in enumerate(rows)])
            if not np.isclose(values.mean(),g['mean_process_summary'][key],rtol=1e-12,atol=1e-14):raise ValueError('Altered water attribution')
        for r in rows:
            s=r['settings'];a=np.array(r['process_diagnostics']['blocks'])
            cap=s['snowfall_mm_h']/3600/s['depth_m']*s['layer_volume_m3'] if s['growth'] else 0
            if not np.allclose(a[:,:,5]+a[:,:,6]+a[:,:,7],cap*a[:,:,0],rtol=2e-12,atol=1e-20):raise ValueError('Water attribution does not close')
            if not np.isclose(a[:,:,5].sum()/s['area_m2'],r['mass_ledger_kg_m2']['deposited'],rtol=2e-11,atol=1e-11):raise ValueError('Water ledger mismatch')
        quant=np.array([r['summary']['mass_d99_m'] for r in rows])
        if not np.isclose(quant.mean(),g['mean_summary']['mass_d99_m'],rtol=1e-12):raise ValueError('Altered quantile mean')
        if not np.isclose(quant.std(ddof=1)/2,g['standard_error_summary']['mass_d99_m'],rtol=1e-12,atol=1e-15):raise ValueError('Altered quantile uncertainty')
        cdf=np.array([np.array(r['arriving_cumulative_moments'])[:,1]/r['arriving_cumulative_moments'][-1][1] for r in rows])
        if not np.allclose(cdf.mean(axis=0),g['mean_mass_cdf'],rtol=1e-12,atol=1e-14):raise ValueError('Altered CDF')
    for r in p['comparisons']:
        a,b=[p['groups'][r[k]] for k in ('left','right')]
        delta=np.max(np.abs(np.array(a['mean_mass_cdf'])-b['mean_mass_cdf']))
        changes=[abs(b['mean_summary'][k]/a['mean_summary'][k]-1) for k in ('snowfall_mm_h','mass_d50_m','mass_d90_m','mass_d99_m')]
        if not np.isclose(delta,r['maximum_mean_mass_cdf_difference'],rtol=1e-12,atol=1e-14):raise ValueError('Altered CDF comparison')
        if bool(delta<=.05 and max(changes)<=.05)!=r['passed']:raise ValueError('Altered refinement decision')
    return verified


def render(path=OUT):
    p=json.loads((ROOT/PRODUCT).read_text());extra=json.loads((ROOT/FOLLOWUP).read_text())
    verified=validate(p);verified.update(validate(extra,expected_cases=8))
    if extra['primary_diagnosis']['sha256']!=digest(ROOT/PRODUCT):raise ValueError('Changed primary diagnosis')
    allgroups=p['groups']|extra['groups'];comparisons=p['comparisons']+extra['comparisons']
    path.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
        'axes.spines.right':False,'svg.hashsalt':'terluna-snow-column-diagnosis-v1'})
    fig,axes=plt.subplots(2,2,figsize=(13.6,10.6));fig.subplots_adjust(left=.08,right=.98,top=.86,bottom=.17,hspace=.4,wspace=.26)
    fig.text(.08,.95,'Why the lunar snow column changes with volume',fontsize=21,weight='bold')
    fig.text(.08,.914,'Observed water allocation and controlled process removals · prescribed model',fontsize=12,color='#555555')
    colors=['#126e82','#c56a28','#80519b'];groups=[p['groups'][f'full_v{n}'] for n in (1,4,16)]
    volume=[g['settings']['layer_volume_m3'] for g in groups]
    row=next(r['result'] for r in p['cases'] if r['group']=='full_v1')
    height=(np.array(row['layer_height_m'])-row['layer_height_m'][0]+row['settings']['depth_m']/60)/1000
    ax=axes[0,0]
    for g,v,c in zip(groups,volume,colors):
        mean=np.array(g['mean_process_profiles']['empty_fraction'])*100
        se=np.array(g['standard_error_process_profiles']['empty_fraction'])*100
        ax.plot(mean,height,color=c,label=f'{v:.3g} m³ per layer')
        ax.fill_betweenx(height,np.maximum(0,mean-se),mean+se,color=c,alpha=.12)
    ax.set(xlabel='Sampled time with no particles (%)',ylabel='Height above freezing level (km)',title='A  Empty layers in the complete column')
    ax.legend(fontsize=9);ax.grid(alpha=.15)
    ax=axes[0,1];x=np.arange(3);bottom=np.zeros(3)
    for key,label,color in [('gather_deposition_mm_h','Deposited','#126e82'),('gather_empty_unused_mm_h','Unused: empty layer','#c56a28'),('gather_occupied_unused_mm_h','Unused: insufficient demand','#a7a7a7')]:
        val=np.array([g['mean_process_summary'][key] for g in groups]);ax.bar(x,val,bottom=bottom,label=label,color=color,width=.62);bottom+=val
    ax.set(xticks=x,xticklabels=[f'{v:.3g}' for v in volume],xlabel='Statistical volume per layer (m³)',ylabel='Water equivalent during gathering (mm/h)',title='B  Accounting for the 3 mm/h vapour supply',ylim=(0,3.2))
    ax.legend(fontsize=8,loc='lower right');ax.grid(axis='y',alpha=.15)
    ax=axes[1,0];x=np.arange(len(FAMILIES))
    for scale,offset,c in [(1,-.18,colors[0]),(4,.18,colors[1])]:
        g=[allgroups[f'{family}_v{scale}'] for family in FAMILIES]
        means=[100*r['mean_summary']['mass_d99_m'] for r in g]
        se=[100*r['standard_error_summary']['mass_d99_m'] for r in g]
        ax.bar(x+offset,means,yerr=se,width=.33,color=c,capsize=3,label=f'{volume[0]*scale:.3g} m³ per layer')
    ax.set(xticks=x,xticklabels=LABELS,yscale='log',ylabel='Mean mass 99th-percentile dimension (cm)',title='C  Upper sizes after process removal')
    ax.legend(fontsize=9);ax.grid(axis='y',alpha=.15)
    ax=axes[1,1];values=[]
    for family in FAMILIES:
        r=next(r for r in comparisons if r['left']==f'{family}_v1' and r['right']==f'{family}_v4')
        values.append(100*r['maximum_mean_mass_cdf_difference'])
    ax.bar(x,values,color=[colors[1] if v>5 else colors[0] for v in values],width=.6)
    ax.axhline(5,color='#555555',linestyle='--',linewidth=1,label='5-percentage-point screen')
    for i,v in enumerate(values):
        label_y=6.3 if 4.0<=v+.6<=5.5 else v+.6
        ax.text(i,label_y,f'{v:.2f}',ha='center',fontsize=9)
    ax.set(xticks=x,xticklabels=LABELS,ylabel='Maximum sampled mass-CDF change\n(percentage points)',title='D  Change under fourfold volume increase',ylim=(0,max(values)*1.17+1))
    ax.legend(fontsize=8);ax.grid(axis='y',alpha=.15)
    fig.text(.08,.092,'Study: 48 runs, four realizations per setting; bands and error bars show one whole-realization standard error.',fontsize=10)
    fig.text(.08,.069,'All runs cover 132.59 h. Water and occupancy statistics use the last 66.29 h and the growth-stage quadrature.',fontsize=10)
    fig.text(.08,.046,'Panel D tests the CDF alone; overall acceptance also requires flux and all three mass quantiles to change by at most 5%.',fontsize=9)
    fig.text(.08,.023,'Removing a process changes subsequent particle histories. These experiments do not validate giant-flake physics or storm optics.',fontsize=9)
    outputs={}
    for suffix in ('png','svg'):
        target=path/f'column_diagnosis.{suffix}'
        fig.savefig(target,dpi=180,metadata={'Date':None} if suffix=='svg' else {'Software':'Terluna'})
        outputs[target.name]=digest(target)
    plt.close(fig)
    manifest=dict(schema='terluna.visualization.snow-column-diagnosis/1',inputs={name:digest(ROOT/name) for name in (PRODUCT,FOLLOWUP)},
        producer={str(Path(__file__).resolve().relative_to(ROOT)):digest(Path(__file__))},verified_inputs=verified,
        outputs=outputs,physical_validation=False,reading_rule=p['evidence'])
    (path/'column_diagnosis_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args();print(json.dumps(render(args.output)['outputs']))

if __name__=='__main__':main()
