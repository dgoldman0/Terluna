"""Account for nonlinear allocation and particle-number changes in saved runs.

The mean-demand calculation smooths each saved gathering block before applying
its supply limit. It is a frozen-history arithmetic comparison, not a rerun with
redistributed water or a predicted correction to the snow distribution.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from atmosphere.snow import column_diagnosis_runs as runs

ROOT=runs.ROOT
INPUT=runs.OUT
OUT=ROOT/'atmosphere/snow/results/column_diagnosis_analysis.json'
FOLLOWUP=ROOT/'atmosphere/snow/results/column_diagnosis_aggregation_only.json'


def allocation_gap(demand,admitted,supply):
    demand,admitted,supply=np.broadcast_arrays(np.asarray(demand,float),np.asarray(admitted,float),np.asarray(supply,float))
    if not all(np.isfinite(x).all() and (x>=0).all() for x in (demand,admitted,supply)):
        raise ValueError('Nonnegative finite water integrals required')
    upper=np.minimum(demand,supply);tol=2e-12*np.maximum(supply,1e-30)
    if (admitted>upper+tol).any():raise ValueError('Admitted water exceeds demand or supply')
    return np.maximum(0.,upper-admitted)


def replay(raw):
    runs.validate_trace(raw)
    a=np.array(raw['process_diagnostics']['blocks'])[1:];s=raw['settings']
    c=s['snowfall_mm_h']/3600/s['depth_m']*s['layer_volume_m3'] if s['growth'] else 0
    supply=c*a[:,:,0]
    gap=allocation_gap(a[:,:,4],a[:,:,5],supply)
    norm=s['area_m2']*s['gathering_time_s']/3600
    by_layer=gap.sum(axis=0)/norm
    return dict(gap_mm_h=float(by_layer.sum()),gap_by_layer_mm_h=by_layer,
        observed_deposition_mm_h=float(a[:,:,5].sum()/norm),
        mean_demand_admission_mm_h=float(np.minimum(a[:,:,4],supply).sum()/norm),
        gathering_blocks=4,block_duration_s=s['gathering_time_s']/4)


def quantile_bracket(cdf,diameters,probability):
    cdf=np.asarray(cdf,float)
    if not 0<probability<1 or len(cdf)!=len(diameters) or (np.diff(cdf)<-1e-12).any() or cdf[-1]<probability:
        raise ValueError('Invalid CDF or probability')
    index=int(np.searchsorted(cdf,probability))
    return [0. if index==0 else diameters[index-1],diameters[index]]


def analyse(product,additional=None):
    if not product['checks']['complete']:raise ValueError('Complete matrix required')
    if additional is not None:
        if not additional['checks']['complete']:raise ValueError('Complete follow-up required')
        if set(product['groups']) & set(additional['groups']):raise ValueError('Overlapping experiment groups')
        product=dict(cases=product['cases']+additional['cases'],groups=product['groups']|additional['groups'])
    groups={}
    for name,group in product['groups'].items():
        rows=[r['result'] for r in product['cases'] if r['group']==name]
        values=[replay(r) for r in rows]
        gap=np.array([r['gap_mm_h'] for r in values]);profiles=np.array([r['gap_by_layer_mm_h'] for r in values])
        averaged=np.array([r['mean_demand_admission_mm_h'] for r in values])
        moments=np.array([r['arriving_cumulative_moments'][-1] for r in rows]);duration=rows[0]['settings']['gathering_time_s']
        number_flux=moments[:,0]*3600/duration;mass_flux=moments[:,1]*3600/duration
        area=rows[0]['settings']['area_m2']
        actual_counts=area*moments[:,0];effective_counts=area*moments[:,1]**2/moments[:,2]
        if (effective_counts>actual_counts*(1+1e-12)).any():raise ValueError('Invalid arrival second moment')
        groups[name]=dict(mean_allocation_gap_mm_h=float(gap.mean()),standard_error_allocation_gap_mm_h=float(gap.std(ddof=1)/np.sqrt(len(gap))),
            mean_admission_after_demand_averaging_mm_h=float(averaged.mean()),standard_error_admission_after_demand_averaging_mm_h=float(averaged.std(ddof=1)/np.sqrt(len(gap))),
            mean_gap_by_layer_mm_h=profiles.mean(axis=0),standard_error_gap_by_layer_mm_h=profiles.std(axis=0,ddof=1)/np.sqrt(len(gap)),
            mean_number_flux_m2_h=float(number_flux.mean()),standard_error_number_flux_m2_h=float(number_flux.std(ddof=1)/np.sqrt(len(gap))),
            ratio_of_mean_mass_to_number_flux_kg=float(mass_flux.mean()/number_flux.mean()),
            individual_represented_arrival_counts=actual_counts,individual_effective_arrival_mass_counts=effective_counts,
            mean_mass_cdf_99th_bracket_m=quantile_bracket(group['mean_mass_cdf'],rows[0]['query_diameter_m'],.99),
            individual_allocation_replays=values)
    decomposition=[]
    for left,right in [('full_v1','full_v4'),('full_v4','full_v16')]:
        a,b=[product['groups'][name]['mean_process_summary'] for name in (left,right)]
        gained=a['gather_empty_unused_mm_h']-b['gather_empty_unused_mm_h']
        offset=b['gather_occupied_unused_mm_h']-a['gather_occupied_unused_mm_h']
        actual=b['gather_deposition_mm_h']-a['gather_deposition_mm_h']
        if not np.isclose(gained-offset,actual,rtol=1e-12,atol=1e-12):raise ValueError('Change in water allocation does not close')
        mass_ratio=product['groups'][right]['mean_summary']['snowfall_mm_h']/product['groups'][left]['mean_summary']['snowfall_mm_h']
        number_ratio=groups[right]['mean_number_flux_m2_h']/groups[left]['mean_number_flux_m2_h']
        decomposition.append(dict(left=left,right=right,less_unused_supply_in_empty_layers_mm_h=gained,
            more_unused_supply_in_occupied_layers_mm_h=offset,actual_deposition_increase_mm_h=actual,
            arriving_mass_flux_ratio=mass_ratio,arriving_number_flux_ratio=number_ratio,
            mean_arriving_particle_mass_ratio=mass_ratio/number_ratio))
    return dict(groups=groups,full_volume_decomposition=decomposition)


def main():
    product=json.loads(INPUT.read_text());additional=json.loads(FOLLOWUP.read_text())
    if additional['primary_diagnosis']['sha256']!=hashlib.sha256(INPUT.read_bytes()).hexdigest():raise RuntimeError('Changed primary matrix')
    for source in (product,additional):
        for category in ('files','inputs'):
            for path,sha in source['producer'][category].items():
                if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=sha:raise RuntimeError('Changed input: '+path)
    result=analyse(product,additional)
    result.update(schema='terluna.atmosphere.snow-column-allocation-analysis/1',
        inputs={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (INPUT,FOLLOWUP)},
        producer={str(Path(__file__).resolve().relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        evidence='Arithmetic on observed histories. min(mean demand, supply) is at least mean min(demand, supply); '
            'the recorded gap includes all demand variability within each 16.57-hour gathering block. It does not isolate '
            'stochastic variability from temporal evolution and supplies no feedback-corrected size distribution.',
        reading_rule='Allocation gaps have units kg/m2/h, numerically mm/h of liquid-water equivalent. They retain four '
            'separate gathering blocks, then sum over layers. Mean-CDF quantiles are brackets on the native query grid; '
            'they differ from means of individual realization quantiles and are not confidence intervals.')
    runs.old.write(OUT,result)
    print(json.dumps({k:v['mean_allocation_gap_mm_h'] for k,v in result['groups'].items()}))

if __name__=='__main__':main()
