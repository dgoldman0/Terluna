"""Render finite-time capture controls from versioned atmosphere products."""
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

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.snow-finite-capture-benchmarks/1'
REFERENCE = 'atmosphere/snow/results/collision_reference.json'
PATHS = {'boundary':'atmosphere/snow/results/collision_jump_boundary.json',
         'bath':'atmosphere/snow/results/collision_jump_bath.json', 'reference':REFERENCE}
COLORS = ['#17667a','#b34b32','#7b5d9e','#58883e','#b18525','#616970']
STYLES = ['-','--',':','-.',(0,(5,2,1,2)),(0,(2,2))]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_close(actual, expected, label, rtol=3e-12, atol=0.):
    if not np.allclose(actual,expected,rtol=rtol,atol=atol):
        raise ValueError(label+': stored and independently reconstructed quantities differ')


def validate_run(run,case,reference):
    """Reconstruct native moments, accounting and CDF errors without model imports."""
    capture=reference['capture'];target=capture['continuous_mass_references'][-1]
    unit=1. if case=='boundary' else capture['initial_collector_mass_kg']
    pivots=np.asarray(run['pivots_mass_ratio' if case=='boundary' else 'pivots_kg'],np.longdouble)
    widths=np.diff(np.asarray(run['bounds_mass_ratio' if case=='boundary' else 'bounds_kg'],np.longdouble))
    for j,row in enumerate(run['snapshots']):
        state=np.asarray(row['state'],np.longdouble)
        if state.shape!=(len(pivots),3) or not np.isfinite(state).all() or np.any(state[:,0]<0):
            raise ValueError('Invalid native moment state')
        n,u,u2=state.T
        moments=[n.sum(),np.sum(n*pivots+u*widths),np.sum(n*pivots**2+2*u*pivots*widths+u2*widths**2)]
        require_close(moments,row['moments'],'native moments')
        if row['time_s']:
            active=n>0;mean=moments[1]/moments[0];local_mean=u[active]/n[active]
            variance=np.sum(n[active]*((pivots[active]+widths[active]*local_mean-mean)**2+
                widths[active]**2*np.maximum(0.,u2[active]/n[active]-local_mean**2)))/moments[0]
            require_close(variance,row['variance'],'stable variance',rtol=2e-7)
        exits=np.asarray(row['ledgers']['upper_exit'])+row['ledgers']['lower_exit']
        require_close(moments[0]+exits[0],1.,'number ledger',atol=1e-10)
        require_close((moments[1]+exits[1]-row['ledgers']['supplied_mass'])/unit,1.,'mass ledger',atol=1e-10)
        query=np.asarray(row['query_mass_ratio' if case=='boundary' else 'number_cdf_query_mass_kg'])
        cdf=np.asarray(row['number_cdf'])
        if not np.isfinite(cdf).all() or np.any(np.diff(query)<=0) or np.min(np.diff(cdf)) < -1e-10 or cdf.min() < -1e-10 or cdf.max()>1+1e-10:
            raise ValueError('Invalid number CDF')
        if case=='boundary':
            exact=np.asarray(row['exact_poisson_number_cdf']);reported=row['number_cdf_supremum']
        else:
            ref=target['snapshots'][j]
            if ref['time_s']!=row['time_s']:raise ValueError('Reference times differ')
            mass=np.asarray(ref['masses_kg']);counts=np.asarray(ref['sample_counts'])
            cumulative=np.r_[0.,np.cumsum(counts)/counts.sum()]
            exact=cumulative[np.searchsorted(mass,query,side='right')];reported=row['empirical_cdf_supremum']
        require_close(np.max(abs(cdf-exact)),reported,'CDF difference',atol=2e-15)


def prepare(root=ROOT):
    products={};identities={};verified={}
    for key,name in PATHS.items():
        data=(root/name).read_bytes();products[key]=json.loads(data);identities[name]=hashlib.sha256(data).hexdigest()
        p=products[key]
        expected='terluna.atmosphere.snow-collision-reference/1' if key=='reference' else SCHEMA
        if p.get('schema')!=expected:raise ValueError('Unsupported product schema: '+name)
        for group in ('files','inputs'):
            for file,digest in p['producer'][group].items():
                if sha(root/file)!=digest:raise ValueError('Producer or input identity changed: '+file)
                verified[file]=digest
    reference=products['reference']
    for case in ('boundary','bath'):
        p=products[case]
        if p['case']!=('poisson_boundary' if case=='boundary' else 'prescribed_bath'):
            raise ValueError('Incorrect control case')
        if p['producer']['inputs'][REFERENCE]!=identities[REFERENCE]:raise ValueError('Reference binding differs')
        if p['execution_failures'] or len(p['runs'])!=len(p['requested_controls']) or not p['acceptance']['passed']:
            raise ValueError('Figure requires a completed, accepted control matrix')
        for run in p['runs']:validate_run(run,case,reference)
    capture=reference['capture']
    for key,value in products['bath']['settings'].items():
        if value!=capture[key]:raise ValueError('Bath setting differs from reference: '+key)
    ref=capture['continuous_mass_references'][-1];last=ref['snapshots'][-1]
    mass=np.asarray(last['masses_kg']);number=np.asarray(last['sample_counts']);total=number.sum()
    mean=float(number@mass/total)
    variance=float(number@(mass-mean)**2/(total-1))
    require_close(mean,last['mean_mass_kg'],'reference mean')
    require_close(variance,last['variance_mass_kg2'],'unbiased sample variance',rtol=2e-10)
    return products,identities,verified


def label(run,case):
    s=run['settings'];text=f"{run['mass_cells']} cells, q{s['order']}, {s['dt']:g} s"
    if case=='boundary' and s['phase']!=1.:text+=', shifted'
    if case=='bath' and s['capacity']!=256:text+=f", C{s['capacity']}"
    return text


def render(output,root=ROOT):
    products,inputs,verified=prepare(root)
    reference=products['reference']['capture'];target=reference['continuous_mass_references'][-1]
    final=target['snapshots'][-1];unit=reference['initial_collector_mass_kg']
    sd=final['variance_mass_kg2']**.5;mean=final['mean_mass_kg'];band=target['number_cdf_uniform_95_percent_dkw_halfwidth']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.labelcolor':'#172b3a','text.color':'#172b3a','svg.hashsalt':'terluna-capture-1'})
    fig,axes=plt.subplots(2,2,figsize=(13.6,10.2),gridspec_kw={'hspace':.52,'wspace':.28})
    fig.subplots_adjust(left=.085,right=.97,top=.88,bottom=.18)
    fig.suptitle('Finite-time capture: distribution controls',fontsize=19,x=.085,ha='left',y=.97)
    fig.text(.085,.925,'Six analytic and six snow-kernel controls • prescribed replenished bath • no lunar size prediction',fontsize=11)
    a,b,c,d=axes.ravel();plotted={'boundary':[],'bath':[]}
    for ax in axes.ravel():ax.grid(alpha=.22);ax.set_axisbelow(True)
    for i,run in enumerate(products['boundary']['runs']):
        rows=run['snapshots'][1:4];face=run['snapshots'][2]['poisson_mean']
        phase=[(row['poisson_mean']-face)/np.sqrt(face) for row in rows]
        error=[100*row['number_cdf_supremum'] for row in rows]
        a.plot(phase,error,color=COLORS[i],ls=STYLES[i],marker='o',ms=3,label=label(run,'boundary'))
        plotted['boundary'].append({'label':label(run,'boundary'),'phase':phase,'cdf_error_percentage_points':error,
                                    'all_saved_errors':[row['number_cdf_supremum'] for row in run['snapshots']],
                                    'wall_seconds':run['wall_seconds'],'macrosteps':run['diagnostics']['macrosteps']})
    a.axhline(.5,color='#a23c34',ls='--',lw=1)
    a.text(-1.95,.515,'Declared 0.5-point threshold',fontsize=9,color='#a23c34')
    a.set(xlim=(-2.15,2.15),ylim=(0,.57),xticks=[-2,0,2],xlabel='Expected captures relative to own cell-face crossing\n(units of crossing-count standard deviation)',ylabel='Maximum number-CDF error (percentage points)')
    a.set_title('A  Analytic Poisson crossings',loc='left',fontweight='bold')
    a.legend(fontsize=8,loc='lower right',ncol=2,frameon=False)
    sample=np.asarray(final['masses_kg']);counts=np.asarray(final['sample_counts']);cdf=np.cumsum(counts)/counts.sum()
    limits=mean+np.array([-6.,6.])*sd
    band_mass=np.r_[limits[0],sample,limits[1]];band_center=np.r_[0.,cdf,1.]
    sx=100*(band_mass-mean)/unit
    b.fill_between(sx,np.maximum(0,band_center-band),np.minimum(1,band_center+band),step='post',color='#9caeb7',alpha=.5,label='95% DKW band at this time')
    b.step(sx,band_center,where='post',color='#172b3a',lw=1.6,label='8,192 exact-event histories',zorder=4)
    b.text(.03,.91,'Six capture curves use the colours in panel C',transform=b.transAxes,fontsize=8)
    exclusions=[]
    for i,run in enumerate(products['bath']['runs']):
        rows=run['snapshots'];last=rows[-1]
        query=np.asarray(last['number_cdf_query_mass_kg']);curve=np.asarray(last['number_cdf'])
        x=100*(query-mean)/unit
        b.plot(x,curve,color=COLORS[i],ls=STYLES[i],lw=1.3,alpha=.9)
        errors=[100*row['empirical_cdf_supremum'] for row in rows[1:]]
        c.plot([0,1],errors,color=COLORS[i],ls=STYLES[i],marker='o',ms=4,label=label(run,'bath'))
        d.bar(i,run['wall_seconds'],color=COLORS[i],width=.65)
        d.text(i,run['wall_seconds']+4,f"{run['diagnostics']['macrosteps']} steps",ha='center',fontsize=8,rotation=0)
        excluded={'below':float(np.interp(limits[0],query,curve)),
                  'above':float(last['moments'][0]-np.interp(limits[1],query,curve))}
        exclusions.append(excluded)
        plotted['bath'].append({'label':label(run,'bath'),'final_query_mass_kg':query.tolist(),'final_number_cdf':curve.tolist(),
            'cdf_errors_percentage_points':errors,'variance_ratios':[row['variance_ratio'] for row in rows[1:]],
            'wall_seconds':run['wall_seconds'],'macrosteps':run['diagnostics']['macrosteps'],'display_exclusions':excluded})
    b.set(xlim=100*(limits-mean)/unit,ylim=(-.02,1.02),xlabel='Mass departure from empirical mean\n(% of initial collector mass)',ylabel='Cumulative fraction of initial collector number')
    b.set_title(f"B  Snow kernel after {final['time_s']/60:.2f} minutes",loc='left',fontweight='bold')
    b.legend(loc='lower right',fontsize=8,frameon=False)
    c.axhline(100*band,color='#a23c34',ls='--',lw=1)
    c.text(.02,100*band+.04,'95% DKW sampling half-width',fontsize=9,color='#a23c34')
    c.set(xticks=[0,1],xticklabels=[f"{row['time_s']/60:.2f} min" for row in target['snapshots'][1:]],ylim=(0,1.72),xlim=(-.08,1.08),ylabel='Maximum number-CDF difference (percentage points)')
    c.set_title('C  Every snow-kernel variation meets its CDF gate',loc='left',fontweight='bold',fontsize=11)
    c.legend(fontsize=8,loc='lower left',ncol=2,frameon=False)
    short=['191\nbase','191\nhalf step','191\ndouble q','96\ncells','381\ncells','191\n384 events']
    d.set(xticks=range(6),xticklabels=short,ylabel='Recorded wall time per control (s)',ylim=(0,max(r['wall_seconds'] for r in products['bath']['runs'])*1.17))
    d.set_title('D  Recorded runtime of each completed control',loc='left',fontweight='bold',fontsize=11)
    fig.text(.085,.065,'Each grid’s analytic samples use its own first face. Bath curves retain the original collector normalization.\nThe DKW band is uniform in mass at each time separately; it describes reference sampling, not numerical convergence.\nAll cases impose a 15.86 cm collector. Formation, donor depletion, aggregate survival and storm optics are outside these controls.\nInitial controls used two workers; the higher-quadrature repeat used one. Runtimes describe these executions.',fontsize=9,linespacing=1.45)
    output.mkdir(parents=True,exist_ok=True)
    identity={'inputs':inputs,'renderer':sha(Path(__file__))}
    encoded=json.dumps(identity,sort_keys=True)
    png=output/'capture_comparison.png';svg=output/'capture_comparison.svg'
    fig.savefig(png,dpi=180,metadata={'Description':encoded})
    fig.savefig(svg,metadata={'Description':encoded,'Date':None})
    plt.close(fig)
    for path,digest in {**inputs,**verified}.items():
        if sha(root/path)!=digest:raise ValueError('Input changed during rendering: '+path)
    manifest={'schema':'terluna.visualization.snow-finite-capture/1','identity':identity,'verified_upstream_files':verified,
              'environment':{'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__},
              'plotted':plotted,'reference':{'masses_kg':sample.tolist(),'number_cdf':cdf.tolist(),'mean_mass_kg':mean,
                'unbiased_variance_mass_kg2':final['variance_mass_kg2'],'sample_count':int(counts.sum()),'dkw_halfwidth':band,
                'time_s':final['time_s'],'window_kg':limits.tolist(),'band_plot_mass_kg':band_mass.tolist(),'band_plot_center':band_center.tolist()},
              'acceptance':{case:products[case]['acceptance'] for case in ('boundary','bath')},
              'evidence':{case:products[case]['evidence'] for case in ('boundary','bath')},
              'reading_rules':{case:products[case]['reading_rule'] for case in ('boundary','bath')},
              'display_rule':'Stored continuous CDF queries joined by lines; empirical CDF is right-continuous. Original number normalization retained, including tails outside the view. No domain-model imports or new physical calculation.',
              'output_hashes':{p.name:sha(p) for p in (png,svg)}}
    (output/'capture_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print('Rendered',png)
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'generated')
    args=parser.parse_args();render(args.output)


if __name__=='__main__':main()
