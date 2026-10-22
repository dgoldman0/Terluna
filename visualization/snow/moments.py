"""Show stored collision boundary controls and narrow/broad bath distributions.

    OPENBLAS_NUM_THREADS=1 python visualization/snow/moments.py

Reads domain products only. Within-cell reconstruction and collision physics remain in
atmosphere; this renderer checks exported moments and number-CDF comparisons.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter, FuncFormatter
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
REFERENCE='atmosphere/snow/results/collision_reference.json'
CANDIDATE='atmosphere/snow/results/collision_moment_bath.json'
SCHEMAS={REFERENCE:'terluna.atmosphere.snow-collision-reference/1',
         CANDIDATE:'terluna.atmosphere.snow-collision-moment-benchmarks/1'}
INK,MUTED,GRID='#172b3a','#52646d','#dce3e6'
BENCHMARK_SCHEMAS={
    'bounded_beta':'terluna.atmosphere.snow-collision-moment-benchmarks/1',
    'bounded_quadratic_maximum_entropy':'terluna.atmosphere.snow-collision-entropy-benchmarks/1'}
OLD_COLORS={96:'#bd602c',191:'#ba927c',381:'#d2b6a5'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_benchmark(product,case):
    if (product.get('measure_name') not in BENCHMARK_SCHEMAS
        or product.get('schema')!=BENCHMARK_SCHEMAS[product['measure_name']]
        or product.get('case')!=case):
        raise ValueError('unsupported benchmark schema, reconstruction or case')


def candidate_style(run,baseline_step):
    if run['quadrature_order']>2:return '-.'
    if run['max_step_s']<baseline_step:return (0,(5,1,1,1,1,1))
    return {1:'-',2:'--',4:':'}.get(run['grid_factor'],'-')


def check(actual,expected,label,*,rtol=2e-10,atol=1e-13):
    if not np.isclose(actual,expected,rtol=rtol,atol=atol):
        raise ValueError(f'{label}: reconstructed {actual} differs from stored {expected}')


def discrete_curve(mass,number,total,mean,variance,label,*,sample_variance=False):
    mass,number=np.asarray(mass,float),np.asarray(number,float)
    if mass.shape!=number.shape or mass.ndim!=1 or not np.isfinite(mass).all() or not np.isfinite(number).all():
        raise ValueError(label+': invalid distribution arrays')
    if np.any(mass<=0) or np.any(number<0):raise ValueError(label+': negative mass/count')
    check(number.sum(),total,label+'/number')
    actual_mean=np.dot(number,mass)/total
    cdf_variance=np.dot(number,(mass-actual_mean)**2)/total
    actual_variance=cdf_variance*total/(total-1) if sample_variance else cdf_variance
    check(actual_mean,mean,label+'/mean',atol=0)
    check(actual_variance,variance,label+'/variance',rtol=2e-7,atol=0)
    keep=number>0;order=np.argsort(mass[keep],kind='stable')
    mass,number=mass[keep][order],number[keep][order]
    return dict(label=label,mass_kg=mass.tolist(),number_cdf=(np.cumsum(number)/total).tolist(),
                represented_number_per_initial_collector=1.,mean_mass_kg=actual_mean,variance_mass_kg2=actual_variance,
                variance_rule='unbiased independent-sample variance (ddof=1)' if sample_variance else 'number-weighted population variance',
                cdf_population_variance_mass_kg2=cdf_variance,
                interpolation='right-continuous discrete empirical/pivot CDF')


def evaluate(curve,query):
    x,y=np.asarray(curve['mass_kg']),np.asarray(curve['number_cdf'])
    if curve['interpolation'].startswith('right-continuous'):
        index=np.searchsorted(x,query,side='right')-1
        return np.where(index<0,0.,y[np.maximum(index,0)])
    return np.interp(query,x,y,left=0.,right=y[-1])


def exclusions(curve,lower,upper):
    # Below lower is a strict inequality, above upper excludes its own CDF jump.
    low=float(evaluate(curve,np.nextafter(lower,-np.inf)))
    high=float(curve['represented_number_per_initial_collector']-evaluate(curve,upper))
    return dict(below=low,above=high,represented_total=curve['represented_number_per_initial_collector'])


def prepare(root=ROOT, candidate_relative=CANDIDATE):
    products,inputs,verified={},{},{}
    for relative,schema in {REFERENCE:SCHEMAS[REFERENCE],candidate_relative:SCHEMAS[CANDIDATE]}.items():
        raw=(root/relative).read_bytes();p=json.loads(raw)
        if relative==REFERENCE:
            if p.get('schema')!=schema:raise ValueError(relative+': unsupported reference schema')
        else:validate_benchmark(p,'prescribed_bath')
        if not p.get('evidence') or not p.get('reading_rule'):
            raise ValueError(relative+': missing evidence/reading rule')
        inputs[relative]=hashlib.sha256(raw).hexdigest();products[relative]=p
        for group in ('files','inputs'):
            for name,expected in p['producer'][group].items():
                if sha(root/name)!=expected:raise ValueError('producer/input changed: '+name)
                verified[name]=expected
    source,candidate=products[REFERENCE],products[candidate_relative]
    if candidate['producer']['inputs'].get(REFERENCE)!=inputs[REFERENCE]:
        raise ValueError('candidate does not bind the plotted frozen reference')
    for name,value in source['producer']['constants'].items():
        if candidate['producer']['constants'].get(name)!=value:raise ValueError('constant record differs: '+name)
    capture=source['capture'];settings=candidate['settings'];initial=capture['initial_collector_mass_kg']
    for name in ('initial_collector_mass_kg','captured_mass_kg','initial_collector_diameter_m',
                 'captured_diameter_m','bath_number_m3','times_s'):
        if settings[name]!=capture[name]:raise ValueError('candidate physical setting differs: '+name)
    ref=max(capture['continuous_mass_references'],key=lambda r:r['samples'])
    target=max(ref['snapshots'],key=lambda r:r['time_s']);duration=target['time_s']
    exact=discrete_curve(target['masses_kg'],target['sample_counts'],ref['samples'],
        target['mean_mass_kg'],target['variance_mass_kg2'],'Exact events',sample_variance=True)
    old=[]
    for grid in capture['fixed_pivot']:
        snap=next(s for s in grid['snapshots'] if s['time_s']==duration)
        n=np.asarray(snap['numbers_per_initial_collector']);total=float(n.sum())
        check(total+snap['upper_exit_number'],1.,'fixed-pivot number ledger')
        mean=snap['mass_kg_m3']/snap['number_m3']
        variance=snap['second_mass_moment_kg2_m3']/snap['number_m3']-mean**2
        c=discrete_curve(grid['mass_pivots_kg'],n,total,mean,variance,f'Old fixed pivots: {grid["mass_bins"]}')
        # The frozen reference has zero exits. Refuse silently renormalized curves.
        check(total,1.,'fixed-pivot retained normalization',atol=1e-12)
        c.update(mass_cells=grid['mass_bins'],color=OLD_COLORS[grid['mass_bins']],linestyle='--',
                 upper_exit_number=snap['upper_exit_number'])
        old.append(c)
    candidates=[]
    for run in candidate['runs']:
        row=next(s for s in run['snapshots'] if s['time_s']==duration)
        state=np.asarray(row['state'],np.longdouble);pivots=np.asarray(run['pivots_kg'],np.longdouble)
        bounds=np.asarray(run['bounds_kg'],np.longdouble);width=np.diff(bounds)
        if state.shape!=(len(pivots),3) or bounds.shape!=(len(pivots)+1,) or np.any(width<=0) or not np.isfinite(state).all() or np.any(state[:,0]<0):
            raise ValueError('invalid stored local-moment state or grid')
        numbers,first,second=state.T
        reconstructed=[numbers.sum(),np.sum(pivots*numbers+width*first),
                       np.sum(pivots*pivots*numbers+2*pivots*width*first+width*width*second)]
        for i,value in enumerate(reconstructed):check(float(value),row['moments'][i],f'candidate moment {i}',atol=0)
        n,m1,m2=map(float,row['moments']);mean=m1/n
        active=numbers>0;local_mean=first[active]/numbers[active]
        local_variance=second[active]/numbers[active]-local_mean**2
        stable=float(np.sum(numbers[active]*((pivots[active]+width[active]*local_mean-mean)**2
                                            +width[active]**2*local_variance))/n)
        variance=row['variance_mass_kg2']
        check(mean,row['mean_mass_kg'],'candidate mean',atol=0)
        check(stable,variance,'candidate stable variance',rtol=2e-7,atol=initial**2*3e-15)
        ledgers=row['ledgers'];exits=np.asarray(ledgers['lower_exit'])+np.asarray(ledgers['upper_exit'])
        check(n+exits[0]+ledgers['collision_number_loss']-1.,row['number_residual'],'candidate number ledger')
        check(m1+exits[1]-ledgers['supplied_mass']-initial,row['mass_residual'],'candidate mass ledger',atol=initial*1e-12)
        check(row['reference_mean_mass_kg'],exact['mean_mass_kg'],'candidate reference mean',atol=0)
        check(row['reference_variance_mass_kg2'],exact['variance_mass_kg2'],'candidate reference variance',rtol=2e-7,atol=0)
        check(row['variance_ratio'],variance/exact['variance_mass_kg2'],'candidate variance ratio',rtol=2e-7)
        query,cdf=np.asarray(row['number_cdf_query_mass_kg']),np.asarray(row['number_cdf'])
        if (query.ndim!=1 or cdf.shape!=query.shape or not np.isfinite(query).all() or not np.isfinite(cdf).all()
            or np.any(np.diff(query)<=0) or np.any(np.diff(cdf)<-1e-12) or cdf.min() < -1e-12 or cdf.max()>n+1e-10):
            raise ValueError('invalid exported candidate CDF')
        check(cdf[-1],n,'candidate CDF total')
        ref_mass=np.asarray(target['masses_kg']);counts=np.asarray(target['sample_counts'])
        if not np.isin(ref_mass,query).all() or not np.isin(np.nextafter(ref_mass,-np.inf),query).all():
            raise ValueError('CDF query omits an empirical jump or its left limit')
        empirical=np.cumsum(counts)/ref['samples']
        distance=float(max(np.max(np.abs(np.interp(ref_mass,query,cdf)-empirical)),
            np.max(np.abs(np.interp(np.nextafter(ref_mass,-np.inf),query,cdf)-np.r_[0.,empirical[:-1]]))))
        check(distance,row['empirical_cdf_supremum'],'candidate empirical CDF distance')
        order,dt=run['quadrature_order'],run['max_step_s']
        family='Beta' if candidate.get('measure_name')=='bounded_beta' else 'Entropy'
        label=f'Bath {family}: {run["mass_cells"]} cells · q{order} · {dt:g} s'
        candidates.append(dict(label=label,mass_cells=run['mass_cells'],grid_factor=run['grid_factor'],
            quadrature_order=order,max_step_s=dt,color='#c65b36' if family=='Beta' else '#246ca6',
            linestyle=candidate_style(run,1.),measure_name=candidate['measure_name'],mass_kg=query.tolist(),number_cdf=cdf.tolist(),
            represented_number_per_initial_collector=n,mean_mass_kg=mean,variance_mass_kg2=variance,
            variance_ratio=row['variance_ratio'],empirical_cdf_supremum=distance,ledgers=row['ledgers'],
            reconstructed_stable_variance_mass_kg2=stable,number_residual=row['number_residual'],mass_residual=row['mass_residual'],
            interpolation='linear display interpolation of producer-exported continuous CDF queries'))
    sd=exact['variance_mass_kg2']**.5;center=exact['mean_mass_kg']
    narrow=(center-6*sd,center+6*sd)
    right=initial
    for c in [exact,*old,*candidates]:
        j=np.searchsorted(c['number_cdf'],.9999)
        right=max(right,c['mass_kg'][min(j,len(c['mass_kg'])-1)])
    broad=(.99*initial,initial+1.05*(right-initial))
    for c in [exact,*old,*candidates]:
        c['view_exclusions']={name:exclusions(c,*bounds) for name,bounds in (('narrow',narrow),('broad',broad))}
        if c['view_exclusions']['broad']['above']>1e-4+1e-10:
            raise ValueError('broad view hides more than0.01% of represented collectors')
    return products,inputs,verified,dict(settings=settings,time_s=duration,reference_samples=ref['samples'],
        reference_dkw_halfwidth=ref['number_cdf_uniform_95_percent_dkw_halfwidth'],reference=exact,
        old_fixed_pivots=old,candidates=candidates,reference_standard_deviation_initial_percent=100*sd/initial,
        narrow_mass_window_kg=list(narrow),broad_mass_window_kg=list(broad),
        acceptance=candidate['acceptance'],
        interpretation='Number CDF per initial collector; no renormalization to surviving particles or display windows. The narrow view is centered on the exact-reference mean and spans six reference SD on either side. Its excluded shares are recorded for every curve. Broad-view exclusions remain in normalization. The reference SD uses its unbiased independent-sample variance (ddof=1); candidate variance is the reconstructed population moment. DKW bands describe exact-reference sampling at this one endpoint, not uncertainty in the candidate closure or a mass-weighted CDF. Candidate CDF physics is supplied by the domain producer; connecting its query points is a display interpolation.')


def prepare_boundary(relative, root=ROOT):
    raw=(root/relative).read_bytes();product=json.loads(raw)
    validate_benchmark(product,'poisson_boundary')
    if not product.get('evidence') or not product.get('reading_rule'):
        raise ValueError('boundary evidence/read rule missing')
    inputs={relative:hashlib.sha256(raw).hexdigest()};verified={}
    for group in ('files','inputs'):
        for name,expected in product['producer'][group].items():
            if sha(root/name)!=expected:raise ValueError('boundary source/input differs: '+name)
            verified[name]=expected
    runs=[]
    for run in product['runs']:
        rows=[]
        for row in run['snapshots']:
            if row['time_s']==0:continue
            x=np.asarray(row['query_mass_ratio']);cdf=np.asarray(row['number_cdf']);exact=np.asarray(row['exact_poisson_number_cdf'])
            if x.ndim!=1 or x.shape!=cdf.shape or x.shape!=exact.shape or not np.isfinite([x,cdf,exact]).all():
                raise ValueError('invalid stored boundary CDF arrays')
            if np.any(np.diff(x)<=0) or np.any(np.diff(cdf)<-1e-10) or np.any(np.diff(exact)<-1e-10):
                raise ValueError('nonmonotone boundary query/CDF')
            error=float(np.max(abs(cdf-exact)))
            check(error,row['number_cdf_supremum'],'boundary CDF supremum')
            rate=run['prescribed_capture_rate_s'];donor=run['donor_mass_ratio'];count=rate*row['time_s']
            check(count,row['poisson_mean'],'Poisson count mean')
            check(1+count*donor,row['expected_mean_mass_ratio'],'Poisson mass mean')
            check(count*donor*donor,row['expected_variance_mass_ratio'],'Poisson mass variance',atol=0)
            n,m1,m2=row['moments'];check(m1/n,row['mean_mass_ratio'],'boundary candidate mean')
            state=np.asarray(row['state'],np.longdouble);pivots=np.asarray(run['pivots_mass_ratio'],np.longdouble)
            width=np.diff(np.asarray(run['bounds_mass_ratio'],np.longdouble))
            if state.shape!=(len(pivots),3) or width.shape!=pivots.shape or np.any(width<=0) or not np.isfinite(state).all() or np.any(state[:,0]<0):
                raise ValueError('invalid boundary native state or grid')
            numbers,first,second=state.T
            native=[numbers.sum(),np.sum(pivots*numbers+width*first),
                    np.sum(pivots*pivots*numbers+2*pivots*width*first+width*width*second)]
            for i,value in enumerate(native):check(float(value),row['moments'][i],f'boundary native moment {i}',atol=0)
            active=numbers>0;local_mean=first[active]/numbers[active]
            stable=float(np.sum(numbers[active]*((pivots[active]+width[active]*local_mean-m1/n)**2
                               +width[active]**2*(second[active]/numbers[active]-local_mean**2)))/n)
            check(stable,row['variance_mass_ratio'],'boundary native stable variance',rtol=2e-7,atol=3e-15)
            ledger=row['ledgers'];exits=np.asarray(ledger['lower_exit'])+np.asarray(ledger['upper_exit'])
            check(n+exits[0]+ledger['collision_number_loss']-1,row['number_residual'],'boundary number ledger')
            check(m1+exits[1]-ledger['supplied_mass']-1,row['mass_residual'],'boundary mass ledger')
            if np.min(cdf)<-1e-10 or np.max(cdf)>n+1e-10 or np.min(exact)<0 or np.max(exact)>1:
                raise ValueError('boundary CDF leaves its number range')
            rows.append({key:row[key] for key in ('time_s','query_mass_ratio','number_cdf','exact_poisson_number_cdf',
                'expected_mean_mass_ratio','expected_variance_mass_ratio','mean_mass_ratio','variance_mass_ratio',
                'number_cdf_supremum','unqueried_left_probability','unqueried_right_probability','ledgers',
                'number_residual','mass_residual')})
        family='Beta' if product['measure_name']=='bounded_beta' else 'Entropy'
        runs.append(dict(label=f'Poisson {family}: {run["mass_cells"]} cells · q{run["quadrature_order"]} · {run["max_step_s"]:g} s',
            measure_name=product['measure_name'],mass_cells=run['mass_cells'],quadrature_order=run['quadrature_order'],
            max_step_s=run['max_step_s'],color='#c65b36' if family=='Beta' else '#246ca6',
            linestyle=candidate_style(run,.5),donor_mass_ratio=run['donor_mass_ratio'],prescribed_capture_rate_s=run['prescribed_capture_rate_s'],
            reference_cell_face_mass_ratio=run['reference_cell_face_mass_ratio'],snapshots=rows,
            acceptance=product['acceptance']))
    return product,inputs,verified,runs


def prepare_all(extra_bath=None,extra_boundary=None):
    products,inputs,verified,plot=prepare()
    if extra_bath:
        other,pins,checks,addition=prepare(candidate_relative=extra_bath)
        products.update(other);inputs.update(pins);verified.update(checks)
        plot['candidates'].extend(addition['candidates'])
        if addition['settings']!=plot['settings'] or addition['time_s']!=plot['time_s']:
            raise ValueError('bath conditions or observation differ between candidates')
        plot['broad_mass_window_kg']=[min(plot['broad_mass_window_kg'][0],addition['broad_mass_window_kg'][0]),
                                     max(plot['broad_mass_window_kg'][1],addition['broad_mass_window_kg'][1])]
    for curve in [plot['reference'],*plot['old_fixed_pivots'],*plot['candidates']]:
        curve['view_exclusions']={name:exclusions(curve,*plot[name+'_mass_window_kg']) for name in ('narrow','broad')}
        if sum(curve['view_exclusions']['broad'][key] for key in ('below','above'))>1e-4+1e-10:
            raise ValueError('common broad view hides more than 0.01% of represented collectors')
    boundaries=[]
    for path in ('atmosphere/snow/results/collision_moment_boundary.json',extra_boundary):
        if not path:continue
        product,pins,checks,runs=prepare_boundary(path)
        products[path]=product;inputs.update(pins);verified.update(checks);boundaries.extend(runs)
    first=boundaries[0]
    for run in boundaries[1:]:
        for name in ('donor_mass_ratio','prescribed_capture_rate_s','reference_cell_face_mass_ratio'):
            check(run[name],first[name],'matched boundary '+name)
        if [s['time_s'] for s in run['snapshots']]!=[s['time_s'] for s in first['snapshots']]:
            raise ValueError('boundary observations differ between candidates')
    face=first['reference_cell_face_mass_ratio']
    crossing=min(first['snapshots'],key=lambda s:abs(s['expected_mean_mass_ratio']-face))
    check(crossing['expected_mean_mass_ratio'],face,'analytic mean at cell face',atol=1e-12)
    plot['boundary_runs']=boundaries;plot['crossing_time_s']=crossing['time_s']
    plot['control_statuses']={name:{'completed_controls':len(p['runs']),'case':p.get('case'),'measure_name':p.get('measure_name'),
        **{key:p[key] for key in ('control_scope','requested_controls','execution_failure','acceptance') if key in p}}
        for name,p in products.items() if 'runs' in p}
    plot['boundary_interpretation']='Analytic constant-rate Poisson control, distinct from the snow-kernel bath below. Observation times are fixed at the 191-cell reference face for all grid comparisons; other grids need not have a face at this mass. CDF errors use saved queries at both sides of each Poisson jump within twelve count SD; unqueried analytic probabilities are retained. No Monte Carlo confidence band belongs to these exact Poisson CDFs. Agreement later in the same control does not clear a boundary failure.'
    return products,inputs,verified,plot


def render(plot,output,identity):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,
                         'axes.labelcolor':INK,'svg.hashsalt':'terluna-collision-moments-v2'})
    fig,axes=plt.subplots(2,2,figsize=(14.2,10.6))
    fig.subplots_adjust(left=.075,right=.98,bottom=.17,top=.72,wspace=.25,hspace=.48)
    settings=plot['settings'];initial=settings['initial_collector_mass_kg'];ref=plot['reference'];center=ref['mean_mass_kg']
    fig.text(.075,.955,'Collision closure at cell boundaries and late times',fontsize=21,weight='bold')
    beta=next(run for run in plot['boundary_runs'] if run['measure_name']=='bounded_beta')
    beta_crossing=next(row for row in beta['snapshots'] if row['time_s']==plot['crossing_time_s'])
    entropy_runs=[run for run in plot['boundary_runs'] if run['measure_name']=='bounded_quadratic_maximum_entropy']
    subtitle=f'Beta misses the cell-face number CDF by {100*beta_crossing["number_cdf_supremum"]:.1f} percentage points despite accurate moments.'
    if entropy_runs:
        entropy_crossing=next(row for row in entropy_runs[0]['snapshots'] if row['time_s']==plot['crossing_time_s'])
        subtitle=(f'Cell-face number-CDF error: {100*beta_crossing["number_cdf_supremum"]:.1f} percentage points for Beta; '
                  f'{100*entropy_crossing["number_cdf_supremum"]:.3f} for bounded maximum entropy.')
    fig.text(.075,.914,subtitle,fontsize=12,color=MUTED)
    handles=[Line2D([],[],color=INK,lw=2,label='Exact reference'),Patch(facecolor=INK,alpha=.12,label='95% empirical-CDF band (lower panels)')]
    labels=set()
    for c in [*plot['boundary_runs'],*plot['candidates'],*plot['old_fixed_pivots']]:
        if c['label'] in labels:continue
        labels.add(c['label']);handles.append(Line2D([],[],color=c['color'],ls=c.get('linestyle','-'),lw=1.7,label=c['label']))
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.069,.885),frameon=False,ncol=3,fontsize=8.8,columnspacing=1.8,handlelength=2.3)
    # Analytic boundary test: exact Poisson, with no empirical uncertainty band.
    face=plot['boundary_runs'][0]['reference_cell_face_mass_ratio'];rows=[]
    for run in plot['boundary_runs']:
        row=next(s for s in run['snapshots'] if s['time_s']==plot['crossing_time_s']);rows.append(row)
        sd=row['expected_variance_mass_ratio']**.5;x=(np.asarray(row['query_mass_ratio'])-face)/sd
        axes[0,0].plot(x,row['number_cdf'],color=run['color'],ls=run['linestyle'],lw=1.9)
    row=rows[0];sd=row['expected_variance_mass_ratio']**.5
    axes[0,0].step((np.asarray(row['query_mass_ratio'])-face)/sd,row['exact_poisson_number_cdf'],where='post',color=INK,lw=1.5)
    axes[0,0].axvline(0,color='#8c989e',lw=1,ls=':',zorder=0)
    axes[0,0].set(xlim=(-6,6),ylim=(-.015,1.015),xlabel='Mass relative to reference face (analytic SD units)',
                  ylabel='Fraction of initial collectors at or below mass')
    axes[0,0].set_title('a  Constant-rate control at the 191-cell reference face',loc='left',fontsize=11,weight='bold',pad=12)
    axes[0,0].yaxis.set_major_formatter(PercentFormatter(1.))
    gate=plot['boundary_runs'][0]['acceptance']['gates']['maximum_number_cdf_supremum']
    annotated=set()
    for run in plot['boundary_runs']:
        check(run['acceptance']['gates']['maximum_number_cdf_supremum'],gate,'common boundary shape gate')
        errors=[100*s['number_cdf_supremum'] for s in run['snapshots']]
        axes[0,1].plot(range(len(errors)),errors,marker='o',ls=run['linestyle'],color=run['color'],lw=1.6,ms=3)
        if run['measure_name'] not in annotated:
            for i,e in enumerate(errors):
                axes[0,1].annotate(f'{e:.3f} pp',(i,e),xytext=(0,-15 if i==len(errors)-1 and run['measure_name']!='bounded_beta' else 7),textcoords='offset points',ha='center',fontsize=8,color=run['color'],bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.3))
            annotated.add(run['measure_name'])
    axes[0,1].axhline(100*gate,color='#8c989e',lw=1,ls='--')
    axes[0,1].text(.02,100*gate*1.12,f'Shape gate: {100*gate:g} percentage points',fontsize=8.5,color=MUTED)
    minimum=min(s['number_cdf_supremum'] for r in plot['boundary_runs'] for s in r['snapshots'])
    axes[0,1].set(yscale='log',ylim=(max(100*minimum/3,1e-4),50),xlim=(-.2,3.2),
                  ylabel='Number-CDF error (percentage points)')
    axes[0,1].set_xticks(range(4),['Before face','At face','After face','Later'])
    axes[0,1].yaxis.set_major_formatter(FuncFormatter(lambda value,_:f'{value:g}'))
    axes[0,1].set_title('b  Shape error before and after the boundary crossing',loc='left',fontsize=11,weight='bold',pad=12)
    # Physical-kernel bath endpoints remain a separate comparison.
    for index,(axis,window) in enumerate(zip(axes[1],(plot['narrow_mass_window_kg'],plot['broad_mass_window_kg']))):
        offset=center if index==0 else initial
        transform=lambda mass:100*(np.asarray(mass)-offset)/initial
        lower,upper=transform(window)
        if index==1:
            for c in plot['old_fixed_pivots']:
                x=np.r_[min(window[0],c['mass_kg'][0])*.999,c['mass_kg'],max(window[1],c['mass_kg'][-1])*1.001]
                y=np.r_[0.,c['number_cdf'],c['number_cdf'][-1]]
                axis.step(transform(x),y,where='post',lw=1.5,color=c['color'],ls='--')
        x=np.r_[min(window[0],ref['mass_kg'][0])*.999,ref['mass_kg'],max(window[1],ref['mass_kg'][-1])*1.001]
        y=np.r_[0.,ref['number_cdf'],1.];band=plot['reference_dkw_halfwidth']
        axis.fill_between(transform(x),np.maximum(y-band,0),np.minimum(y+band,1),step='post',color=INK,alpha=.12,zorder=1)
        for c in plot['candidates']:
            axis.plot(transform(c['mass_kg']),c['number_cdf'],color=c['color'],ls=c['linestyle'],lw=1.8,zorder=3)
        axis.step(transform(x),y,where='post',color=INK,lw=1.4,zorder=4)
        axis.set_xlim(lower,upper);axis.set_ylim(-.015,1.015)
        axis.set_ylabel('Fraction of initial collectors at or below mass');axis.yaxis.set_major_formatter(PercentFormatter(1.))
    axes[1,0].set_xlabel('Departure from empirical mean (% of initial mass)')
    axes[1,0].set_title(f'c  Snow-kernel bath at {plot["time_s"]/60:.2f} min: narrow view',loc='left',fontsize=11,weight='bold',pad=12)
    axes[1,1].set_xlabel('Mass increase (% of initial collector mass)')
    axes[1,1].set_title('d  Same endpoint: old fixed-pivot spreading retained',loc='left',fontsize=11,weight='bold',pad=12)
    for axis in axes.flat:
        axis.spines[['top','right']].set_visible(False);axis.spines[['left','bottom']].set_color('#aab8bf')
        axis.tick_params(colors=MUTED,labelsize=9);axis.grid(axis='y',color=GRID,lw=.6);axis.set_axisbelow(True)
    fig.text(.075,.116,f'Upper panels: analytic Poisson capture, no sampling uncertainty. Lower panels: {plot["reference_samples"]:,} exact-event histories; '
             f'empirical SD = {plot["reference_standard_deviation_initial_percent"]:.4f}% of initial mass.',fontsize=9.4,color=MUTED)
    fig.text(.075,.084,'The lower close-up spans ±6 reference SD. Original collector normalization and all excluded tails remain in the manifest; q denotes quadrature order.',fontsize=9.4,color=MUTED)
    fig.text(.075,.047,'These imposed-collector tests omit deposition, settling and breakup. Beta fails the boundary-shape test; late endpoint agreement is not adoption or a lunar size prediction.',fontsize=9.4)
    incomplete=[status for status in plot['control_statuses'].values()
                if status.get('execution_failure') is not None
                or status.get('acceptance',{}).get('all_requested_controls_completed') is False]
    if incomplete:
        counts=[]
        for status in plot['control_statuses'].values():
            if status['measure_name']!='bounded_quadratic_maximum_entropy':continue
            name={'prescribed_bath':'physical bath','poisson_boundary':'Poisson'}[status['case']]
            counts.append(f'{name} {status["completed_controls"]}/{len(status.get("requested_controls",[]))}')
        fig.text(.075,.021,'Completed entropy controls: '+', '.join(counts)+'. Refinement remains incomplete; failure details are retained in the manifest.',fontsize=9.2,color='#a83b27')
    output.mkdir(parents=True,exist_ok=True);description=json.dumps(identity,sort_keys=True)
    paths=[output/'moment_comparison.png',output/'moment_comparison.svg']
    fig.savefig(paths[0],dpi=180,metadata={'Description':description})
    fig.savefig(paths[1],metadata={'Date':None,'Title':'Collision closure at boundaries and late times','Description':description})
    plt.close(fig)
    return {p.name:sha(p) for p in paths}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'generated')
    parser.add_argument('--entropy-bath',help='Explicit additional saved product path, relative to repository root')
    parser.add_argument('--entropy-boundary',help='Explicit additional saved product path, relative to repository root')
    args=parser.parse_args();products,inputs,verified,plot=prepare_all(args.entropy_bath,args.entropy_boundary)
    identity=dict(renderer={'visualization/snow/moments.py':sha(__file__)},inputs=inputs)
    outputs=render(plot,args.output,identity)
    manifest=dict(schema='terluna.visualization.snow-collision-moments/2',**identity,
        verified_upstream_files_and_inputs=verified,environment=dict(matplotlib=matplotlib.__version__,numpy=np.__version__),
        upstream={name:{k:p[k] for k in ('schema','producer','evidence','reading_rule','measure_name','control_scope',
                                        'requested_controls','execution_failure','acceptance') if k in p}
                  for name,p in products.items()},
        plotted=plot,outputs=outputs)
    path=args.output/'moment_manifest.json';path.write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(manifest=str(path),outputs=outputs,
        reference_SD_initial_percent=plot['reference_standard_deviation_initial_percent'],
        boundary_maximum_errors={c['label']:max(s['number_cdf_supremum'] for s in c['snapshots']) for c in plot['boundary_runs']}),indent=2))


if __name__=='__main__':main()
