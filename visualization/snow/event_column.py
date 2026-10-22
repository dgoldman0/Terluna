"""Render source-bound full-column controls, including failed refinements."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
PRODUCT='atmosphere/snow/results/event_column.json'
COLORS=['#126b82','#c36a27','#844c8f']
SCHEMA='terluna.atmosphere.snow-event-column-controls/1'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(actual,expected,label,rtol=3e-12,atol=0.):
    if not np.allclose(actual,expected,rtol=rtol,atol=atol):raise ValueError(label+' differs from native quantities')


def validate(product,root=ROOT):
    if product.get('schema')!=SCHEMA:raise ValueError('Unsupported column schema')
    if product['completed_cases']!=product['expected_cases'] or len(product['cases'])!=product['expected_cases']:
        raise ValueError('Incomplete control matrix')
    verified={}
    for row in product['cases']:
        if row['status']!='complete':raise ValueError('Incomplete case '+row['id'])
        for group in ('files','inputs'):
            for path,digest in row['producer'][group].items():
                if path in verified and verified[path]!=digest:raise ValueError('Inconsistent producer identity')
                if path not in verified and sha(root/path)!=digest:raise ValueError('Stale producer or input '+path)
                verified[path]=digest
        r=row['result'];m=r['mass_ledger_kg_m2'];n=r['number_ledger_m2'];s=r['settings']
        source=m['injected']+m['deposited']
        residual=source-m['freezing_exit']-m['remaining']-m['fractional_fragment_remaining']
        close(residual/source,m['relative_residual'],'water residual',atol=3e-14)
        nsource=n['injected']+n['chipping_production']
        nres=nsource-n['freezing_exit']-n['remaining']-n['fractional_fragment_remaining']-n['coagulation_loss']
        close(nres/nsource,n['relative_residual'],'number residual',atol=3e-14)
        if abs(m['relative_residual'])>1e-10 or abs(n['relative_residual'])>1e-10:raise ValueError('Failed accounting')
        if m['deposited']>s['snowfall_mm_h']*s['total_time_s']/3600*(1+1e-12):raise ValueError('Vapour supply exceeded')
        if m['fractional_fragment_remaining']>r['diagnostics']['maximum_fractional_fragment_mass_bound_kg_m2']*(1+1e-12):
            raise ValueError('Fractional allocation bound exceeded')
        a=np.asarray(r['arriving_cumulative_moments'])
        if a.shape!=(len(r['query_diameter_m']),4) or not np.isfinite(a).all() or np.min(a)<0 or np.min(np.diff(a,axis=0))<0:
            raise ValueError('Invalid arrival moments')
        if r['query_diameter_m'][-1] is not None:raise ValueError('Missing all-arrival endpoint')
        close(a[-1],np.sum(r['arrival_time_blocks']['moments'],axis=0),'arrival block totals')
        close(a[-1,1]*3600/s['gathering_time_s'],r['summary']['snowfall_mm_h'],'arriving flux')
        close(a[-1,0],r['summary']['arriving_number_m2'],'arriving number')
        diameters=np.r_[r['query_diameter_m'][:-1],np.inf]
        mass_cdf=a[:,1]/a[-1,1]
        for key,probability in [('mass_d50_m',.5),('mass_d90_m',.9),('mass_d99_m',.99)]:
            index=int(np.searchsorted(diameters,r['summary'][key],side='left'))
            lower=mass_cdf[index-1] if index else 0.
            if lower>probability+1e-12 or mass_cdf[index]<probability-1e-12:
                raise ValueError('Quantile outside native CDF bracket')
        if a[-1,2]*a[-1,0]<a[-1,1]**2*(1-1e-12):raise ValueError('Negative particle-mass variance')
        close(r['diagnostics']['final_time_s'],s['total_time_s'],'completed model duration')
        if not np.isfinite(r['time_mean_layer_moments']).all() or np.min(r['time_mean_layer_moments'])<0:
            raise ValueError('Invalid height moments')
        if r['diagnostics']['build']['source_sha256']!=row['producer']['files']['atmosphere/snow/event_column.cpp']:
            raise ValueError('Compiled source differs')
    for name,g in product['groups'].items():
        runs=[r['result'] for r in product['cases'] if r['group']==name]
        if g['status']!='complete' or len(runs)!=g['realizations']:raise ValueError('Group incomplete')
        arrays=np.array([r['arriving_cumulative_moments'] for r in runs]);cdfs=arrays[:,:,1]/arrays[:,-1,1,None]
        close(g['mean_mass_cdf'],cdfs.mean(axis=0),'mean CDF')
        close(g['standard_error_mass_cdf'],cdfs.std(axis=0,ddof=1)/np.sqrt(len(runs)),'CDF standard error')
        for key in g['mean_summary']:
            values=[r['summary'][key] for r in runs]
            close(g['mean_summary'][key],np.mean(values),'mean '+key)
            close(g['standard_error_summary'][key],np.std(values,ddof=1)/np.sqrt(len(runs)),'SE '+key)
    for c in product['comparisons']:
        a,b=[product['groups'][c[k]] for k in ('left','right')]
        delta=float(np.max(np.abs(np.array(a['mean_mass_cdf'])-b['mean_mass_cdf'])))
        close(c['maximum_mean_mass_cdf_difference'],delta,'refinement CDF difference')
        for key,value in c['relative_summary_changes'].items():
            close(value,b['mean_summary'][key]/a['mean_summary'][key]-1,'refinement '+key,atol=1e-15)
        expected=delta<=.05 and all(abs(v)<=.05 for v in c['relative_summary_changes'].values())
        if c['passed']!=expected:raise ValueError('Incorrect refinement decision')
    expected=all(c['passed'] for c in product['comparisons'])
    if product['acceptance']['passed']!=expected:raise ValueError('Incorrect overall decision')
    return verified


def render(destination,root=ROOT):
    path=root/PRODUCT;product=json.loads(path.read_text());verified=validate(product,root)
    groups=product['groups'];query=np.array(product['cases'][0]['result']['query_diameter_m'][:-1])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
        'axes.labelsize':10,'svg.hashsalt':'terluna-event-column-v1','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(13.6,10.6))
    fig.subplots_adjust(left=.085,right=.97,bottom=.12,top=.875,wspace=.29,hspace=.34)
    accepted=product['acceptance']['passed']
    fig.suptitle('Full prescribed snow columns',fontsize=20,fontweight='bold',x=.085,ha='left',y=.965)
    fig.text(.085,.924,'Lunar volume refinement '+('passes' if all(c['passed'] for c in product['comparisons'] if c['axis']=='volume') else 'fails')+' · giant-flake sizes remain unvalidated',
             fontsize=13,color='#a13d28' if not accepted else '#126b82')
    def cdf(ax,key,color,label,scale=100,style='-'):
        g=groups[key];mean=np.array(g['mean_mass_cdf'])[:-1];se=np.array(g['standard_error_mass_cdf'])[:-1]
        ax.plot(query*scale,mean,color=color,lw=2,label=label,ls=style)
        ax.fill_between(query*scale,np.maximum(0,mean-se),np.minimum(1,mean+se),color=color,alpha=.13,lw=0)
        ax.set_xscale('log');ax.set_ylim(0,1.02);ax.grid(alpha=.2);ax.set_ylabel('Cumulative arriving mass fraction')
    ax=axes[0,0]
    for key,color in zip(('moon_v1','moon_v4','moon_v16'),COLORS):
        v=groups[key]['settings']['layer_volume_m3'];cdf(ax,key,color,f'{v:.3g} m³ per layer')
    ax.set_xlim(.2,200);ax.set_xlabel('Maximum dimension (cm)');ax.set_title('A  Lunar distribution changes with volume',loc='left')
    ax.legend(frameon=False,fontsize=9,loc='upper left')
    ax=axes[0,1];names=('moon_v1','moon_v4','moon_v16')
    x=[groups[k]['settings']['layer_volume_m3'] for k in names]
    for key,label,color,marker in [('mass_d50_m','Mass median',COLORS[0],'o'),('mass_d99_m','Mass 99th percentile',COLORS[1],'s')]:
        y=[100*groups[k]['mean_summary'][key] for k in names]
        se=[100*groups[k]['standard_error_summary'][key] for k in names]
        ax.errorbar(x,y,yerr=se,color=color,marker=marker,lw=2,capsize=4,label=label)
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xticks(x);ax.set_xticklabels([f'{v:.3g}' for v in x]);ax.set_yticks([10,20,40,80]);ax.set_yticklabels(['10','20','40','80']);ax.grid(alpha=.2);ax.set_xlabel('Statistical volume per layer (m³)')
    ax.set_ylabel('Mean realization quantile (cm)');ax.set_title('B  Volume-dependent lunar sizes',loc='left');ax.legend(frameon=False,fontsize=9)
    ax=axes[1,0]
    for key,color,label,style in [('moon_v4',COLORS[1],'30 layers, 60 s','-'),('moon_dt30',COLORS[0],'30 layers, 30 s','--'),
                                 ('moon_z60',COLORS[2],'60 layers, 60 s',':')]:cdf(ax,key,color,label,style=style)
    ax.set_xlim(.2,200);ax.set_xlabel('Maximum dimension (cm)');ax.set_title('C  Separate time and height controls',loc='left')
    ax.legend(frameon=False,fontsize=9,loc='upper left')
    ax.text(.97,.08,'Layer volume held at '+f"{groups['moon_v4']['settings']['layer_volume_m3']:.3g} m³",transform=ax.transAxes,ha='right',fontsize=9)
    ax=axes[1,1]
    for key,color,style in [('earth_v4',COLORS[0],'-'),('earth_v16',COLORS[1],'--')]:
        v=groups[key]['settings']['layer_volume_m3'];cdf(ax,key,color,f'{v:.3g} m³ per layer',scale=1000,style=style)
    ax.set_xlim(.1,30);ax.set_xlabel('Maximum dimension (mm)');ax.set_title('D  Generic Earth transfer control',loc='left')
    ax.legend(frameon=False,fontsize=9,loc='upper left')
    moon=groups['moon_v4']['settings'];res=max(g['max_absolute_relative_water_residual'] for g in groups.values())
    fig.text(.085,.064,f"Moon: −15 °C top · 3 mm/h supply · 100 seeds/L · {moon['total_time_s']/3600:.2f} h. Earth: 1 mm/h supply.",fontsize=9)
    fig.text(.085,.043,'Four independent realizations per setting; bands and bars show one realization standard error. Upwind settling on all panels.',fontsize=9)
    fig.text(.085,.022,f"Maximum relative water-ledger residual: {res:.2g}. Water conservation does not establish distribution convergence or physical validity.",fontsize=9)
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    metadata={'Creator':'Terluna scientific visualization','Description':PRODUCT+' sha256='+sha(path)}
    png=destination/'event_column_comparison.png';svg=destination/'event_column_comparison.svg'
    fig.savefig(png,dpi=180,metadata=metadata)
    fig.savefig(svg,metadata={**metadata,'Date':None});plt.close(fig)
    manifest=dict(schema='terluna.visualization.snow-event-column/1',inputs={PRODUCT:sha(path)},
        producer={str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__)},verified_producers_and_inputs=verified,
        environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
        outputs={p.name:sha(p) for p in (png,svg)},full_control_matrix=product['completed_cases'],
        numerical_acceptance=accepted,physical_validation=False,
        reading_rule='Panels display source-bound numerical controls. Failed volume or other refinements remain visible. '
                     'No rare-tail confidence bound, observed-storm validation or optical forecast is supplied.')
    (destination/'event_column_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=Path(__file__).parent/'generated')
    a=p.parse_args();print(json.dumps(render(a.output),indent=2))

if __name__=='__main__':main()
