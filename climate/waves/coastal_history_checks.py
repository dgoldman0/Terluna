"""Occurrence, timing and energy checks for the evolving Smythii coast."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np

from climate.waves.coastal import load_terrain
from climate.waves.coastal_checks import iter_nesting_spectra
from climate.waves.cycle import CYCLE_HOURS, exceedance, window
from climate.waves.model import ORIGIN, ROOT, sha256
from climate.waves.shore import SITES
from climate.waves.shore_checks import shore_normal, spectral_transport
from climate.waves.shore_history import RUNS
from shared.constants import MOON_SURFACE_GRAVITY

OUTPUT=Path(__file__).parent/'results/coastal_history.json'


def load_case(path):
    path=Path(path).absolute()
    record=json.loads(path.read_text())
    if record['schema']!='terluna.climate.coastal-history-case/1':
        raise ValueError('Unsupported coastal history')
    for filename,digest in {**record['run']['identity']['inputs'],**record['run']['output_sha256']}.items():
        if sha256(path.parent/filename)!=digest:
            raise ValueError(f'Coastal history bytes changed: {filename}')
    if sha256(path.parent/'producer_source.py')!=record['producer_sha256']:
        raise ValueError('Coastal history producer bytes changed')
    if sha256(ROOT/record['terrain_path'])!=record['terrain_sha256']:
        raise ValueError('Coastal history terrain changed')
    upstream=record['boundary']['identity']
    source=ROOT/upstream['source']
    if (sha256(source)!=upstream['source_sha256'] or
            sha256(source.parent/'product.json')!=upstream['source_record_sha256']):
        raise ValueError('Upstream coastal-history product changed')
    data=np.loadtxt(path.parent/'sites.tbl').reshape(-1,len(SITES),12)
    _,unique=np.unique(data[:,0,0],return_index=True)
    data=data[np.sort(unique)]
    if np.any(np.diff(data[:,0,0])<=0) or not np.isfinite(data).all():
        raise ValueError('Coastal station history must be finite and ordered')
    np.testing.assert_allclose(data[:,0,0],
        np.arange(record['start_hour'],record['end_hour']+.01,.25)*3600,atol=.01,rtol=0)
    np.testing.assert_allclose(data[:,:,1:3],np.broadcast_to(np.array(list(SITES.values())),data[:,:,1:3].shape),
                               atol=1e-5,rtol=0)
    return record,data


def power_history(path,record,data):
    terrain=load_terrain(ROOT/record['terrain_path'])
    normals=[None if name=='offshore' else shore_normal(terrain,xy) for name,xy in SITES.items()]
    settings=next(l for l in (path.parent/'INPUT').read_text().splitlines() if l.startswith('SET '))
    rho=float(re.search(r'RHO\s+([\d.]+)',settings)[1])
    gravity=float(re.search(r'GRAV\s+([\d.]+)',settings)[1])
    np.testing.assert_allclose(gravity,MOON_SURFACE_GRAVITY,rtol=1e-10,atol=0)
    times,values,edges,errors=[],[],[],[]
    columns=['incoming_w_m','outgoing_w_m','net_vector_magnitude_w_m','eastward_w_m','northward_w_m','resolved_hs_m']
    for spectrum in iter_nesting_spectra(path.parent/'sites.spc'):
        np.testing.assert_allclose(spectrum['locations'],np.array(list(SITES.values())),atol=1e-5,rtol=0)
        hour=(spectrum['dates'][0]-ORIGIN).total_seconds()/3600
        position=np.searchsorted(data[:,0,0],hour*3600)
        np.testing.assert_allclose(data[position,0,0],hour*3600,atol=.01,rtol=0)
        rows=[]; edge=[]; error=[]
        for i,name in enumerate(SITES):
            if not spectrum['available'][0,i]:
                raise ValueError(f'A coastal spectrum is missing at {hour}: {name}')
            f=spectrum['frequency']; s=spectrum['spectra'][0,i]
            power=spectral_transport(f,spectrum['direction'],s,data[position,i,6],rho,MOON_SURFACE_GRAVITY,normals[i])
            variance=s.sum(axis=1)*360/len(spectrum['direction'])
            m0=float(np.trapezoid(variance,f))
            if name=='offshore':
                incoming,outgoing=power['magnitude_w_m'],0.
            else:
                incoming,outgoing=power['shoreward_w_m'],power['seaward_w_m']
            rows.append([incoming,outgoing,power['magnitude_w_m'],*power['vector_w_m'],4*np.sqrt(m0)])
            edge.append([float(np.trapezoid(variance[:3],f[:3])/max(m0,1e-30)),
                         float(np.trapezoid(variance[-3:],f[-3:])/max(m0,1e-30))])
            table=data[position,i,9:11]*rho*MOON_SURFACE_GRAVITY
            error.append(float(np.linalg.norm(table-power['vector_w_m'])/max(np.linalg.norm(table),1.)))
        times.append(hour); values.append(rows); edges.append(edge); errors.append(error)
    np.testing.assert_allclose(times,np.arange(record['start_hour'],record['end_hour']+1),atol=1e-9,rtol=0)
    return np.array(times),np.array(values),np.array(edges),np.array(errors),columns


def analyse(path,output=OUTPUT):
    path=Path(path).absolute()
    record,data=load_case(path)
    if record['level']!='coast' or record['start_hour']!=0 or record['end_hour']<2*CYCLE_HOURS:
        raise ValueError('Cycle statistics require the complete coast and its full first-cycle spin-up')
    t,h=window(data[:,0,0]/3600,data,CYCLE_HOURS,2*CYCLE_HOURS)
    pt,power,edge,error,pcolumns=power_history(path,record,data)
    qt,q=window(pt,power,CYCLE_HOURS,2*CYCLE_HOURS)
    stations={}
    for i,name in enumerate(SITES):
        peak=int(np.argmax(h[:,i,3])); ppeak=int(np.argmax(q[:,i,0]))
        # Saved spectra, used for follow-on individual-wave runs, stay at integer hours.
        available=(pt>=CYCLE_HOURS)&(pt<=2*CYCLE_HOURS)
        saved=np.flatnonzero(available)[np.argmax(power[available,i,0])]
        eligible=available&(power[:,i,5]>=.1)
        band_failure=eligible&(edge[:,i,1]>=.01)
        metre_scale=available&(power[:,i,5]>=1.)
        stations[name]=dict(longitude_latitude=list(SITES[name]),depth_m=float(h[0,i,6]),
            maximum_hs_m=float(h[peak,i,3]),peak_hs_hour=float(t[peak]),
            peak_hs_day_in_second_cycle=float((t[peak]-CYCLE_HOURS)/24),
            mean_period_at_peak_s=float(h[peak,i,5]),
            time_mean_hs_m=float(np.trapezoid(h[:,i,3],t)/CYCLE_HOURS),
            height_thresholds={str(level):exceedance(t,h[:,i,3],level) for level in (.5,1.,1.5,2.)},
            maximum_incoming_w_m=float(q[ppeak,i,0]),peak_incoming_hour=float(qt[ppeak]),
            selected_spectrum_hour=float(pt[saved]),
            mean_incoming_w_m=float(np.trapezoid(q[:,i,0],qt)/CYCLE_HOURS),
            incoming_energy_mj_m=float(np.trapezoid(q[:,i,0],qt)*3600/1e6),
            source_band_checks=dict(eligible_spectra=int(eligible.sum()),
                high_edge_failures=int(band_failure.sum()),
                maximum_high_edge_fraction=float(edge[eligible,i,1].max()),
                failed_resolved_hs_range_m=[float(power[band_failure,i,5].min()),
                    float(power[band_failure,i,5].max())] if band_failure.any() else None,
                metre_scale_eligible_spectra=int(metre_scale.sum()),
                metre_scale_high_edge_failures=int((band_failure&metre_scale).sum())),
            power_thresholds={str(level):exceedance(qt,q[:,i,0],level) for level in (10.,100.,300.)},
            power_definition='Net transport magnitude at offshore reference' if name=='offshore' else 'Incoming power across the local depth contour toward shallower water')
    selected=(pt>=CYCLE_HOURS)&(pt<=2*CYCLE_HOURS)
    energetic=power[selected,:,5]>=.1
    terrain=load_terrain(ROOT/record['terrain_path'])
    stride=record['stride']
    nx=len(terrain['longitude_deg'][::stride]);ny=len(terrain['latitude_deg'][::stride])
    maps=np.loadtxt(path.parent/'maps.tbl').reshape(-1,ny,nx,12)
    final=np.loadtxt(path.parent/'final.tbl').reshape(1,ny,nx,12)
    maps=np.concatenate((maps,final))
    _,indices=np.unique(maps[:,0,0,0],return_index=True)
    maps=maps[np.sort(indices)]
    mt,maps_window=window(maps[:,0,0,0]/3600,maps,CYCLE_HOURS,2*CYCLE_HOURS)
    wet=terrain['wet'][::stride,::stride] & (terrain['depth_m'][::stride,::stride]>.05)
    maximum=np.max(maps_window[:,:,:,3],axis=0)
    mean=np.trapezoid(maps_window[:,:,:,3],mt,axis=0)/CYCLE_HOURS
    arrays=path.parent/'analysis.npz'
    np.savez_compressed(arrays,time_hours=data[:,0,0]/3600,stations=data,
        spectrum_time_hours=pt,power=power,spectral_edges=edge,transport_error=error,
        longitude_deg=terrain['longitude_deg'][::stride],latitude_deg=terrain['latitude_deg'][::stride],
        wet=wet,map_max_hs_m=maximum,map_mean_hs_m=mean)
    result=dict(schema='terluna.climate.coastal-history/1',
        evidence='A complete evolving coastal second cycle after a full cycle of spin-up at basin, regional and coastal scales.',
        reading_rule='Occurrence describes this single simulated cycle. Threshold crossings use piecewise-linear histories. Coast-to-offshore onset differences combine propagation, local growth and direction changes; they require those qualifications when interpreted as arrival delays.',
        producer=dict(path='climate/waves/coastal_history_checks.py',sha256=sha256(Path(__file__))),
        case=dict(path=str(path.relative_to(ROOT)),sha256=sha256(path),stride=record['stride'],step_s=record['step_s']),
        arrays=dict(path=str(arrays.relative_to(ROOT)),sha256=sha256(arrays),station_columns=record['columns'],power_columns=pcolumns),
        clock=dict(spinup_hours=[0,CYCLE_HOURS],report_hours=[CYCLE_HOURS,2*CYCLE_HOURS],
                   station_interval_hours=.25,spectrum_interval_hours=1.),
        units=dict(height='m',period='s',time='Earth hours from basin origin',power='W per metre of contour'),
        stations=stations,
        maps=dict(sampling_hours=6,reading_rule='Spatial maps summarize six-hourly fields, with the final state supplying the reporting endpoint. Station histories have their separately stated sampling.'),
        spectral_checks=dict(maximum_transport_relative_error=float(error[selected].max()),
            transport_pass=bool(error[selected].max()<.01),
            high_edge_failure_fraction=float(np.mean(edge[selected,:,1][energetic]>=.01)),
            maximum_high_edge_fraction=float(edge[selected,:,1][energetic].max())),
        limitations=['One simulated reporting cycle; variability across cycles and other seas remains open.',
                     'The 474 m propagation grid retains the local spatial limits measured in the shore study.',
                     'Temporal and spectral refinements accompany the record; weak waves retain their frequency-band question.',
                     'Currents, changing water levels, diffraction, reflection and sediment evolution are absent from this spectral history.',
                     'Lunar growth and dissipation parameterizations require physical calibration.'])
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({name:{k:s[k] for k in ('maximum_hs_m','maximum_incoming_w_m','selected_spectrum_hour')} for name,s in stations.items()}),flush=True)
    return result


def compare(reference_path,test_path,*,start,end):
    reference,a=load_case(reference_path); test,b=load_case(test_path)
    t,aa=window(a[:,0,0]/3600,a,start,end)
    tt,bb=window(b[:,0,0]/3600,b,start,end)
    np.testing.assert_allclose(t,tt,atol=.001,rtol=0)
    results={}
    for i,name in enumerate(SITES):
        selected=(aa[:,i,3]>=.1)&(bb[:,i,3]>=.1)
        height=abs(aa[selected,i,3]-bb[selected,i,3])/bb[selected,i,3]
        period=abs(aa[selected,i,5]-bb[selected,i,5])/np.maximum(bb[selected,i,5],1.)
        av,bv=aa[selected,i,9:11],bb[selected,i,9:11]
        flux=np.linalg.norm(av-bv,axis=1)/np.maximum(np.linalg.norm(bv,axis=1),1e-3)
        bands={}
        for minimum in (.5,1.):
            energetic=(aa[selected,i,3]>=minimum)&(bb[selected,i,3]>=minimum)
            bands[str(minimum)]=dict(samples=int(energetic.sum()))
            if energetic.any():
                bands[str(minimum)].update(
                    hs_relative_p95=float(np.percentile(height[energetic],95)),
                    hs_relative_max=float(height[energetic].max()),
                    period_relative_p95=float(np.percentile(period[energetic],95)),
                    flux_relative_p95=float(np.percentile(flux[energetic],95)))
        threshold={}
        for level in (.5,1.):
            one=exceedance(t,aa[:,i,3],level);two=exceedance(t,bb[:,i,3],level)
            threshold[str(level)]=dict(reference_hours=one['hours'],test_hours=two['hours'],
                duration_difference_hours=two['hours']-one['hours'],
                reference_intervals_hours=one['intervals_hours'],test_intervals_hours=two['intervals_hours'])
        results[name]=dict(samples=int(selected.sum()),
            hs_relative_p95=float(np.percentile(height,95)),hs_relative_max=float(height.max()),
            period_relative_p95=float(np.percentile(period,95)),flux_relative_p95=float(np.percentile(flux,95)),
            paired_height_bands=bands,
            threshold_comparison=threshold,
            passed=bool(np.percentile(height,95)<.05 and height.max()<.15 and
                        np.percentile(period,95)<.05 and np.percentile(flux,95)<.1))
    return dict(reference=str(Path(reference_path).absolute().relative_to(ROOT)),reference_sha256=sha256(reference_path),
                test=str(Path(test_path).absolute().relative_to(ROOT)),test_sha256=sha256(test_path),
                numerical_bounds=dict(hs_relative_p95=.05,hs_relative_max=.15,
                    period_relative_p95=.05,flux_relative_p95=.1),
                flux_definition='Relative difference of the net east/north energy-transport vectors, with the refined vector as denominator.',
                interval_hours=[start,end],stations=results,all_stations_pass=all(v['passed'] for v in results.values()))


def comparison_product(reference,test,start,end,output,*,append=False):
    comparison=compare(reference,test,start=start,end=end)
    comparisons=[]
    if append:
        previous=json.loads(output.read_text())
        if previous['schema']!='terluna.climate.coastal-history-checks/1':
            raise ValueError('Expected a coastal comparison product to extend')
        for item in previous['comparisons']:
            a,b=ROOT/item['reference'],ROOT/item['test']
            if sha256(a)!=item['reference_sha256'] or sha256(b)!=item['test_sha256']:
                raise ValueError('A previously compared coastal record changed')
            if item['reference']==comparison['reference'] and item['test']==comparison['test']:
                continue
            comparisons.append(compare(a,b,start=item['interval_hours'][0],end=item['interval_hours'][1]))
    comparisons.append(comparison)
    pulse_path=RUNS/'timing_controls/controls.json'
    pulses=json.loads(pulse_path.read_text())
    for step,case in pulses['cases'].items():
        directory=pulse_path.parent/f'dt{step}'
        for name,digest in {**case['run']['identity']['inputs'],**case['run']['output_sha256']}.items():
            if sha256(directory/name)!=digest:
                raise ValueError('Analytic pulse-control bytes changed')
    result=dict(schema='terluna.climate.coastal-history-checks/1',
        evidence='A selected-window timestep refinement with unchanged upstream forcing, and analytic propagation controls.',
        reading_rule='A 24-hour adjustment after the restart precedes the selected comparison window. Other times and upstream refinements retain their stated limits.',
        producer=dict(path='climate/waves/coastal_history_checks.py',sha256=sha256(Path(__file__))),
        comparisons=comparisons,pulse_control=dict(path=str(pulse_path.relative_to(ROOT)),
            sha256=sha256(pulse_path),record=pulses))
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,default=RUNS/'coast_s4_dt900_0_1419/product.json')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--compare',type=Path)
    parser.add_argument('--start',type=float)
    parser.add_argument('--end',type=float)
    parser.add_argument('--append',action='store_true',help='Retain and recheck earlier comparisons in the output product')
    args=parser.parse_args()
    if args.append and not args.compare:
        parser.error('--append requires --compare')
    output=args.output or (OUTPUT.with_name('coastal_history_checks.json') if args.compare else OUTPUT)
    if args.compare:
        if args.start is None or args.end is None:
            parser.error('--compare requires --start and --end')
        comparison_product(args.case,args.compare,args.start,args.end,output,append=args.append)
    else:
        analyse(args.case,output)
