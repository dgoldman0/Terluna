"""Conditional kinematic histories in an archived lunar CM1 storm window.

No climate calculation or growing-particle model is run. Fixed settling drifts
sample resolved, interpolated air. The archive and release selection limit all
exposure statements; neither trajectories nor selected cells estimate weather
occurrence or complete flake lifetimes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from atmosphere.snow import turbulence as t
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT/'climate/crm/runs'
OUT = ROOT/'atmosphere/snow/results/moist_history.json'
CASE = 'box_0e_elec_corrected_storms_second'
SCHEMA = 'terluna.atmosphere.snow-moist-history/1'
CODE = ('atmosphere/snow/moist_history.py', 'atmosphere/snow/turbulence.py',
        'atmosphere/snow/flakes.py', 'atmosphere/snow/moist_history_sources.json')
INVENTORY = ('box_0e_small', 'box_0e_elec_corrected',
             'box_0e_elec_corrected_storms_first', CASE,
             'box_0e_elec_fine', 'ring_45n_corrected', 'ring_equator')
DRIFTS = (0., .3, .6, 1.)
TEMPERATURE_C = (-20., 0.)
SEED = 20261021
SNOW_THRESHOLD = LIQUID_THRESHOLD = 1e-5
SEED_RELEASE_HOURS = (0., 6., 12., 18., 24., 30.)
CONTROL_FILES = ('case.json', 'input_grid_z', 'namelist.input',
                 'cm1out_s.ctl', 'cm1out_metadata.ctl', 'cm1out_metadata.dat')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inventory(name):
    path = RUNS/name
    record = json.loads((path/'case.json').read_text())
    shape, fields, length = t.layout(path/'cm1out_s.ctl')
    _, mf, ml = t.layout(path/'cm1out_metadata.ctl')
    metadata = np.fromfile(path/'cm1out_metadata.dat', dtype='<f4').reshape(-1, ml)
    lookup = {int(row[mf['nwrite'][0]]): float(row[mf['mtime'][0]]) for row in metadata}
    numbers = sorted(int(p.name[8:14]) for p in path.glob('cm1out_t*_s.dat'))
    times = np.array([lookup[n] for n in numbers])
    if len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError('need unique increasing archived times')
    interfaces = np.loadtxt(path/'input_grid_z')
    z = (interfaces[:-1]+interfaces[1:])/2
    levels = np.flatnonzero((z >= 20000.) & (z <= 60000.))
    return dict(case=name, shape=shape, fields=fields, record_length=length,
                numbers=numbers, times=times, z=z, levels=levels,
                dx=float(record['grid']['dx_m']), record=record,
                control_hashes={f: sha(path/f) for f in CONTROL_FILES})


def inventory_record(item):
    r = item['record']
    return dict(case=item['case'], dimensions=item['shape'], dx_m=item['dx'],
                winter_layer_thickness_m=np.unique(np.diff(np.loadtxt(RUNS/item['case']/'input_grid_z'))[item['levels']]).tolist(),
                fields_available=list(item['fields']), saved_frames=len(item['numbers']),
                model_day_range=(item['times'][[0,-1]]/86400).tolist(),
                cadence_s=np.unique(np.diff(item['times'])).tolist(),
                selected_winter_center_range_m=item['z'][item['levels']][[0,-1]].tolist(),
                forcing_run=r.get('reference', {}).get('run'),
                purpose=r['configuration'].get('purpose'),
                build_identity=r.get('build'), controls=item['control_hashes'])


def air(fields):
    """CM1 theta/pressure conversion; humidity is recomputed after interpolation."""
    a = np.asarray(fields)
    th, p, qv, qc, qi, qs, u, v, w = np.moveaxis(a, -1, 0)
    T = th*(p/t.CM1_P00)**(CM1_DRY_AIR_GAS_CONSTANT/t.CM1_CP)
    pv = p*qv/(CM1_DRY_AIR_GAS_CONSTANT/CM1_VAPOUR_GAS_CONSTANT+qv)
    esi = np.exp(9.550426-5723.265/T+3.53068*np.log(T)-.00728332*T)
    return dict(T=T, rh=pv/esi, qc=qc, qs=qs, u=u, v=v, w=w)


def masks(a):
    winter = (a['T'] >= t.sf.T0+TEMPERATURE_C[0]) & (a['T'] <= t.sf.T0)
    return dict(winter=winter, moist=winter & (a['rh'] >= 1.),
                near_saturated=winter & (a['rh'] >= .95),
                dry=winter & (a['rh'] < .95),
                liquid=winter & (a['qc'] >= LIQUID_THRESHOLD),
                snow=winter & (a['qs'] >= SNOW_THRESHOLD))


def stats(values):
    values = np.asarray(values, float)
    if not values.size:
        return dict(count=0, min=None, p10=None, p50=None, p90=None, p99=None, max=None)
    quantiles = np.percentile(values, [0,10,50,90,99,100])
    return dict(count=int(values.size), **dict(zip(('min','p10','p50','p90','p99','max'),map(float,quantiles))))


def load_window(item, end_index=210):
    """Freeze selected raw float32 fields in memory; hash exactly what was read."""
    path = RUNS/item['case']; levels = item['levels']; shape = item['shape']
    numbers = item['numbers'][:end_index+1]
    if len(numbers) != end_index+1 or not np.allclose(np.diff(item['times'][:end_index+1]),600.,rtol=0,atol=0):
        raise ValueError('the declared comparison requires consecutive ten-minute fields')
    data = np.empty((len(numbers),len(t.FIELDS),len(levels),shape['y'],shape['x']),dtype=np.float32)
    sources = []; samples = {g: {k: [] for k in ('horizontal_frame_cells','vertical_frame_layers','rh','qc_g_kg')}
                            for g in ('snow','moist_snow','liquid_snow')}
    dz = np.diff(np.loadtxt(path/'input_grid_z'))[levels,None,None]
    running = np.zeros(data.shape[2:],dtype=int); spans=[]; censored=[]; began_at_start=np.zeros_like(running,bool)
    total_moist=both=0; previous=None
    for j,n in enumerate(numbers):
        digest=hashlib.sha256();digest.update(np.asarray(levels,dtype='<i4').tobytes())
        file=path/f'cm1out_t{n:06d}_s.dat'; st=file.stat()
        fields=t.read_selected(file,shape,item['fields'],item['record_length'],levels,digest)
        data[j]=np.stack([fields[k] for k in t.FIELDS])
        sources.append(dict(number=n,model_time_s=float(item['times'][j]),selected_sha256=digest.hexdigest(),
                            bytes=st.st_size,mtime_ns=st.st_mtime_ns))
        a=air(np.moveaxis(data[j].astype(float),0,-1));m=masks(a)
        groups={'snow':m['snow'],'moist_snow':m['snow']&m['moist'],'liquid_snow':m['snow']&m['liquid']}
        q=dict(horizontal_frame_cells=np.hypot(a['u'],a['v'])*600/item['dx'],
               vertical_frame_layers=np.abs(a['w'])*600/dz,rh=a['rh'],qc_g_kg=a['qc']*1000)
        for g,mask in groups.items():
            for k,v in q.items():samples[g][k].append(v[mask])
        wet=m['moist']
        if j==0:began_at_start=wet.copy()
        if previous is not None:
            total_moist+=int(previous.sum());both+=int((previous&wet).sum())
        finished=(running>0)&~wet
        if finished.any():
            spans.extend(((running[finished]-1)*600).tolist())
            censored.extend(began_at_start[finished].tolist())
        began_at_start=np.where(wet,began_at_start,False)
        running=np.where(wet,running+1,0);previous=wet
        if j%30==0:print(f'read moist-history frame {j+1}/{len(numbers)}',flush=True)
    unfinished=running>0
    spans.extend(((running[unfinished]-1)*600).tolist());censored.extend([True]*int(unfinished.sum()))
    data.setflags(write=False)
    resolution={g:{k:stats(np.concatenate(v)) for k,v in qs.items()} for g,qs in samples.items()}
    for g in samples:
        v=np.concatenate(samples[g]['vertical_frame_layers'])
        resolution[g]['fraction_vertical_displacement_over_one_layer']=float(np.mean(v>1)) if len(v) else None
    return data,sources,dict(conditioned_cell_time_diagnostics=resolution,
        eulerian=dict(moist_spell_observed_span_s=stats(spans),edge_censored_spell_count=int(np.sum(censored)),
                      consecutive_frame_same_cell_moist_fraction=both/total_moist,
                      reading_rule='Spells count consecutive moist samples at fixed grid cells. Their first-to-last sampled span is not a continuous-exposure lower bound or particle residence. With no sub-frame temporal regularity assumption, positive exposure has no nonzero rigorous lower bound; unsampled dry gaps can intervene.'))


def releases(data,z,dx,per_release=64):
    rng=np.random.default_rng(SEED); rows=[]; manifest=[]
    for hour in SEED_RELEASE_HOURS:
        index=int(hour*6)
        a=air(np.moveaxis(data[index].astype(float),0,-1));m=masks(a)
        cells=np.flatnonzero(m['moist']&m['snow'])
        selected=np.sort(rng.choice(cells,min(per_release,len(cells)),replace=False)) if len(cells) else []
        manifest.append(dict(hour=hour,qualifying_cells=int(len(cells)),selected_cells=len(selected)))
        for flat in selected:
            iz,iy,ix=np.unravel_index(flat,data.shape[2:])
            for drift in DRIFTS:
                rows.append([hour*3600,(ix+.5)*dx,(iy+.5)*dx,z[iz],drift])
    if not rows:raise ValueError('no qualifying moist snowy release cells')
    return np.asarray(rows,float),manifest


def spatial(frame,positions,z,dx):
    """Trilinear scalar-point interpolation, periodic horizontally, no z extrapolation."""
    pos=np.asarray(positions,float)
    if ((pos[:,2]<z[0])|(pos[:,2]>z[-1])).any():raise ValueError('vertical interpolation outside selected centers')
    ny,nx=frame.shape[-2:]
    x=(pos[:,0]/dx-.5)%nx;y=(pos[:,1]/dx-.5)%ny
    ix=np.floor(x).astype(int);iy=np.floor(y).astype(int)
    iz=np.clip(np.searchsorted(z,pos[:,2],side='right')-1,0,len(z)-2)
    fx=x-ix;fy=y-iy;fz=(pos[:,2]-z[iz])/(z[iz+1]-z[iz])
    result=np.zeros((len(pos),frame.shape[0]))
    for kz in (0,1):
        for ky in (0,1):
            for kx in (0,1):
                weight=(fz if kz else 1-fz)*(fy if ky else 1-fy)*(fx if kx else 1-fx)
                result+=frame[:,iz+kz,(iy+ky)%ny,(ix+kx)%nx].T*weight[:,None]
    return result


def sample(data,positions,z,dx,seconds,stride):
    interval=600.*stride; last=(len(data)-1)//stride
    left=min(int(seconds//interval),last-1);fraction=(seconds-left*interval)/interval
    return ((1-fraction)*spatial(data[left*stride],positions,z,dx)+
            fraction*spatial(data[(left+1)*stride],positions,z,dx))


def boundary_fraction(start,trial,low,high):
    """First boundary on a linear stage segment, or one if it stays in range."""
    fraction=np.ones(len(start));delta=trial-start
    for boundary,mask in ((low,trial<low),(high,trial>high)):
        fraction[mask]=(boundary-start[mask])/delta[mask]
    return np.clip(fraction,0.,1.)


def integrate(data,z,dx,seeds,*,step_s=60.,stride=1,histories=False):
    started=time.monotonic();n=len(seeds);end=(len(data)-1)*600.
    if end%(600*stride) or 600%step_s or any((seeds[:,0]%step_s)!=0):
        raise ValueError('controls must share exact release and end times')
    position=seeds[:,1:4].copy();live=np.ones(n,bool);status=np.full(n,'time_window',dtype='U20')
    duration=np.zeros(n);exposure={k:np.zeros(n) for k in ('winter','moist','near_saturated','dry','liquid','snow')}
    integral={k:np.zeros(n) for k in ('ice_supersaturation_h','ice_subsaturation_h','cloud_liquid_g_kg_h')}
    temperature_min=np.full(n,np.inf);temperature_max=np.full(n,-np.inf);rh_min=np.full(n,np.inf);rh_max=np.full(n,-np.inf)
    longest=np.zeros(n);run=np.zeros(n);first=np.zeros(n);initial_open=np.ones(n,bool)
    history_ids=np.unique(np.linspace(0,n-1,min(n,24)).round().astype(int)) if histories else []
    history={int(i):[] for i in history_ids}
    for it in range(int(round(end/step_s))):
        seconds=it*step_s;ids=np.flatnonzero(live&(seeds[:,0]<=seconds))
        if not len(ids):continue
        old=position[ids];drift=seeds[ids,4]
        a=air(sample(data,old,z,dx,seconds,stride));v=np.column_stack((a['u'],a['v'],a['w']-drift))
        middle=old+step_s/2*v
        fmid=boundary_fraction(old[:,2],middle[:,2],z[0],z[-1])
        middle[:,2]=np.clip(middle[:,2],z[0],z[-1])
        amid=air(sample(data,middle,z,dx,seconds+step_s/2,stride))
        vmid=np.column_stack((amid['u'],amid['v'],amid['w']-drift))
        trial=old+step_s*vmid
        fraction=boundary_fraction(old[:,2],trial[:,2],z[0],z[-1])
        midexit=fmid<1
        trial[midexit]=old[midexit]+step_s/2*v[midexit]
        fraction[midexit]=fmid[midexit]
        effective=np.where(midexit,step_s/2,step_s)*fraction
        final=old+(trial-old)*fraction[:,None]
        final[:,2]=np.clip(final[:,2],z[0],z[-1])
        # Midpoint exposure on each accepted path segment; partial boundary steps
        # use their own time and position, without extrapolated atmosphere.
        environment=[]
        for elapsed in np.unique(effective):
            selected=effective==elapsed
            environment.append((selected,air(sample(data,(old[selected]+final[selected])/2,z,dx,
                                                    seconds+elapsed/2,stride))))
        for selected,env in environment:
            ii=ids[selected];dt=effective[selected];flags=masks(env)
            duration[ii]+=dt
            for key,flag in flags.items():exposure[key][ii]+=dt*flag
            wet=flags['moist'];run[ii]=np.where(wet,run[ii]+dt,0);longest[ii]=np.maximum(longest[ii],run[ii])
            first[ii]+=dt*(initial_open[ii]&wet);initial_open[ii]&=wet
            integral['ice_supersaturation_h'][ii]+=np.maximum(env['rh']-1,0)*dt/3600*flags['winter']
            integral['ice_subsaturation_h'][ii]+=np.maximum(1-env['rh'],0)*dt/3600*flags['winter']
            integral['cloud_liquid_g_kg_h'][ii]+=env['qc']*1000*dt/3600*flags['winter']
            temperature_min[ii]=np.minimum(temperature_min[ii],env['T']-t.sf.T0)
            temperature_max[ii]=np.maximum(temperature_max[ii],env['T']-t.sf.T0)
            rh_min[ii]=np.minimum(rh_min[ii],env['rh']);rh_max[ii]=np.maximum(rh_max[ii],env['rh'])
        position[ids]=final
        exited=effective<step_s
        status[ids[exited]]=np.where(final[exited,2]<=z[0],'lower_height','upper_height');live[ids[exited]]=False
        if histories and ((it+1)*step_s)%600==0:
            for i in history_ids:
                if seeds[i,0]>(it+1)*step_s or not live[i]:continue
                env=air(sample(data,position[i:i+1],z,dx,(it+1)*step_s,stride))
                history[int(i)].append([float((it+1)*step_s),float(position[i,2]),float(env['T'][0]-t.sf.T0),
                                        float(env['rh'][0]),float(env['qc'][0]*1000),float(env['qs'][0]*1000)])
    arrays=dict(tracked_h=duration/3600,first_moist_episode_h=first/3600,longest_moist_episode_h=longest/3600,
                temperature_min_c=temperature_min,temperature_max_c=temperature_max,rh_min=rh_min,rh_max=rh_max,
                **{k+'_h':v/3600 for k,v in exposure.items()},**integral)
    summary={str(drift):dict(count=int(np.sum(seeds[:,4]==drift)),
                    censoring={s:int(np.sum((seeds[:,4]==drift)&(status==s))) for s in ('time_window','lower_height','upper_height')},
                    quantities={k:stats(v[seeds[:,4]==drift]) for k,v in arrays.items()}) for drift in DRIFTS}
    return dict(settings=dict(step_s=step_s,field_cadence_s=600*stride),wall_seconds=time.monotonic()-started,
                summary=summary,per_tracer={k:v.tolist() for k,v in arrays.items()},
                censoring=status.tolist(),first_moist_episode_right_censored=initial_open.tolist(),
                final_moist_episode_right_censored=(run>0).tolist(),end_height_m=position[:,2].tolist(),
                histories=dict(columns=['window_time_s','height_m','temperature_c','relative_humidity_ice','cloud_liquid_g_kg','snow_g_kg'],
                               selection='24 evenly spaced tracer indices, diagnostic examples rather than representative weights',rows=history))


def comparison(base,other):
    return dict(absolute_paired_difference={key:stats(np.abs(np.asarray(base['per_tracer'][key])-np.asarray(other['per_tracer'][key])))
                                           for key in base['per_tracer']},
                changed_censoring_count=sum(a!=b for a,b in zip(base['censoring'],other['censoring'])))


def analyse(per_release=64):
    before={p:sha(ROOT/p) for p in CODE};constant_record=constants_used(CODE)
    inventories=[inventory(name) for name in INVENTORY];item=next(i for i in inventories if i['case']==CASE)
    namelist=(RUNS/CASE/'namelist.input').read_text()
    import re
    if any(int(re.search(r'^\s*'+k+r'\s*=\s*(-?\d+)',namelist,re.M).group(1))!=value
           for k,value in [('wbc',1),('ebc',1),('sbc',1),('nbc',1),('imove',0)]):
        raise ValueError('selected calculation assumes periodic, stationary horizontal box')
    data,sources,eulerian=load_window(item);z=item['z'][item['levels']]
    seeds,release_manifest=releases(data,z,item['dx'],per_release)
    cases={}
    for name,step,stride in [('base',60.,1),('step30',30.,1),('fields20',60.,2),('fields30',60.,3)]:
        print(f'moist histories {name}: {len(seeds)} tracers',flush=True)
        cases[name]=integrate(data,z,item['dx'],seeds,step_s=step,stride=stride,histories=name=='base')
    if before!={p:sha(ROOT/p) for p in CODE} or constant_record!=constants_used(CODE):
        raise RuntimeError('producer changed during calculation')
    for inv in inventories:
        if inv['control_hashes']!={f:sha(RUNS/inv['case']/f) for f in CONTROL_FILES}:
            raise RuntimeError('archive control changed during calculation')
    for source in sources:
        st=(RUNS/CASE/f'cm1out_t{source["number"]:06d}_s.dat').stat()
        if (st.st_size,st.st_mtime_ns)!=(source['bytes'],source['mtime_ns']):
            raise RuntimeError('archive snapshot metadata changed during calculation')
    return dict(schema=SCHEMA,producer=dict(domain='atmosphere',model=CODE[0],files=before,constants=constant_record),
        evidence='Conditional kinematic tracer histories through archived, coarsely resolved CM1 air. Fixed settling drifts, no flake growth/loss or turbulent dispersion, no empirical validation.',
        reading_rule='All controls share releases and a 35 h analysis window within a selected active episode. Exposure integrals follow interpolated fields; field subsampling is a sensitivity, not evidence that 10 min data resolve missing variability. Time-window and height exits censor histories. Releases are selected positions, not inferred particle births or cloud onsets. Episode durations start at release; no earlier exposure is recovered. These are not complete flake lifetimes or weather frequencies. Raw per-tracer durations are retained for matched comparisons. First moist episodes retain a separate right-censor flag; minima and maxima are sampled values, not bounds on unresolved fluctuations.',
        units=dict(position='m',velocity='m/s',exposure='h',temperature='degrees Celsius',humidity='vapour pressure / ice saturation pressure',cloud_water='g per kg dry air',supersaturation_integrals='dimensionless humidity excess or deficit times hours'),
        source_register=json.loads((ROOT/'atmosphere/snow/moist_history_sources.json').read_text()),
        inventory=[inventory_record(i) for i in inventories],
        source=dict(case=CASE,field_order=list(t.FIELDS),selected_levels=item['levels'].tolist(),frames=sources,
                    selected_byte_hash_rule='SHA256 of little-endian int32 selected levels, then each field name plus NUL and selected C-order little-endian float32 bytes',
                    analysis_model_day_range=(item['times'][[0,210]]/86400).tolist(),available_but_unused_final_frames=5,
                    memory='Selected source fields copied once into immutable memory; source stat and control hashes rechecked afterward'),
        configuration=dict(temperature_c=TEMPERATURE_C,height_center_domain_m=z[[0,-1]].tolist(),drifts_m_s=DRIFTS,
                    snow_threshold_kg_kg=SNOW_THRESHOLD,liquid_threshold_kg_kg=LIQUID_THRESHOLD,
                    release_condition='cold, qs>=1e-5 kg/kg and RHi>=1 at original 10 min fields; sample up to the declared maximum without replacement per release, then repeat each cell for every drift. Equal sampled-cell weights are conditional on release; they are not particle number or mass weights.',
                    seed=SEED,maximum_cells_per_release=per_release,releases=release_manifest,
                    seeds_columns=['release_window_s','x_m','y_m','height_m','settling_drift_m_s'],seeds=seeds.tolist(),
                    interpolation='trilinear scalar-point fields, periodic x/y, linear in time; recompute T and ice humidity from interpolated theta/pressure/vapour',
                    integration='explicit midpoint; midpoint exposure quadrature; height crossing within a step uses linear segment fraction',
                    absent_processes=['growth, sublimation and melting','collisions and breakup','unresolved turbulent dispersion','material tracking of nudging and imposed large-scale advection tendencies']),
        eulerian=eulerian,cases=cases,controls={name:comparison(cases['base'],row) for name,row in cases.items() if name!='base'},
        required_for_particle_prediction=['mass/shape evolution and fall-speed closure with validated geometry','higher-cadence resolved velocity and thermodynamics or demonstrated temporal convergence','spatial refinement under matched forcing and domain','subgrid transport closure where required','longer and representative weather sampling with particle source histories'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--per-release',type=int,default=64)
    args=parser.parse_args()
    if args.per_release<1:parser.error('positive per-release count required')
    result=analyse(args.per_release)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
