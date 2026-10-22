"""Verify native growth traces, morphology parity and recovered source data."""
from pathlib import Path
import json
import math
import numpy as np
from atmosphere.snow import scale_growth_benchmarks as benchmark
from atmosphere.snow import scale_growth_controls as controls
from atmosphere.snow import takahashi_measurements as measurements

CACHE=benchmark.CACHE
ROOT=benchmark.lunar.ROOT


def trace_check(case):
    d=CACHE/'runs'/case
    result=json.loads((d/'benchmark_analysis.json').read_text());s=result['settings']
    path=d/'growth_trace.txt'
    header=path.read_text().splitlines()[0].split()[1:]
    if header!=list(controls.TRACE_FIELDS):raise ValueError('Unidentified growth trace columns')
    data=np.atleast_2d(np.loadtxt(path))
    if not np.isfinite(data).all() or data.shape!=(round(s['duration_s']/s['dt_s']),len(header)):
        raise ValueError('Incomplete or nonfinite trace')
    v={name:data[:,i] for i,name in enumerate(header)}
    n=len(data);dt=s['dt_s'];mode=s['growth_control']
    for key,wanted in [('call',np.arange(1,n+1)),('elapsed_after_s',dt*np.arange(1,n+1)),
                       ('elapsed_before_s',dt*np.arange(n)),('mode',mode),('dt_s',dt)]:
        if not np.array_equal(v[key],np.broadcast_to(wanted,(n,))):raise ValueError('Trace timing or mode mismatch')
    if not np.allclose(v['temperature_K'],s['temperature_K'],rtol=0,atol=1e-9):raise ValueError('Wrong trace temperature')
    if not np.allclose(v['pressure_Pa'],s['pressure_Pa'],rtol=0,atol=1e-6):raise ValueError('Wrong trace pressure')
    if np.any(v['mean_ventilation']<1):raise ValueError('Mean ventilation was disabled')
    expected=v['native_axis_ventilation_ratio'] if mode==0 else np.ones(n)
    if not np.array_equal(v['axis_ventilation_ratio'],expected):raise ValueError('Wrong controlled growth ratio')
    for row in result['snapshots']:
        if row['elapsed_s']==0:continue
        i=round(row['elapsed_s']/dt)-1
        values={'mass_kg':v['new_mass_kg'][i],
                'equatorial_dimension_m':2*v['new_equatorial_radius_m'][i],
                'polar_dimension_m':2*v['new_polar_radius_m'][i],
                'density_kg_m3':v['new_density_kg_m3'][i]}
        if any(not math.isclose(row[k],x,rel_tol=1e-11) for k,x in values.items()):
            raise ValueError('Native growth trace differs from saved particle state')
    # The identical mean-ventilation expression is evaluated using each
    # evolving shape. Its values need not stay equal between the two runs.
    q=['intrinsic_growth_ratio','effective_growth_ratio','mean_ventilation',
       'native_axis_ventilation_ratio','axis_ventilation_ratio','deposition_density_kg_m3',
       'capacitance_m','terminal_speed_m_s','ice_supersaturation']
    ranges={k:[float(v[k].min()),float(v[k].max())] for k in q}
    first=result['snapshots'][0];last=result['snapshots'][-1]
    ratio=math.log(last['polar_dimension_m']/first['polar_dimension_m'])/math.log(last['equatorial_dimension_m']/first['equatorial_dimension_m'])
    return dict(case=case,trace_rows=n,trace_sha256=benchmark.lunar.base.digest(path),
        snapshot_matches=len(result['snapshots'])-1,ranges=ranges,
        integrated_axis_growth_ratio=ratio,endpoint=last)


def run():
    build=controls.check_build()
    products=[json.loads(p.read_text()) for p in (benchmark.OUT,benchmark.OUT.with_name('scale_growth_controls.json'))]
    checks=[]
    for product in products:
        for name,result in product['cases'].items():
            d=CACHE/'runs'/name;s=result['settings']
            spec=dict(temperature_C=s['temperature_C'],dt=s['dt_s'])
            if 'growth_control' in s:spec['control']=s['growth_control']
            repeated=benchmark.analyze(d)
            # Compare actual saved product to independently reread raw output.
            for key in ('snapshots','water','input_sha256','execution','air_ranges'):
                if repeated[key]!=result[key]:raise ValueError('Product differs from native outputs')
            checks.append(dict(case=name,physics_input_matches=benchmark.verify_inputs(spec,d),
                               analysis_sha256=benchmark.lunar.base.digest(d/'benchmark_analysis.json')))
    for dt in (5,1):
        native=CACHE/'runs'/benchmark.case_name(-14.4,dt)
        legacy=CACHE/'runs'/f'growth_earth_observed_dt{dt}'
        instrumented=CACHE/'runs'/benchmark.case_name(-14.4,dt,0)
        if benchmark.native_parity(legacy,native)!=products[0]['baseline_ledger_parity'][str(dt)]:
            raise ValueError('Recorded legacy parity changed')
        if benchmark.native_parity(native,instrumented)!=products[1]['baseline_ledger_parity'][str(dt)]:
            raise ValueError('Recorded mode-zero parity changed')
    traces=[trace_check(name) for name in products[1]['cases']]
    workbook=CACHE/'growth_environment_sources/1991_takahashi.xlsx'
    rows=measurements.read_rows(workbook)
    public=json.loads((ROOT/'atmosphere/snow/results/takahashi_measurements.json').read_text())
    if rows!=public['rows'] or measurements.diagnostic_windows(rows)!=public['diagnostic_windows']:
        raise ValueError('Measurement product differs from identified workbook')
    record=dict(schema='terluna.atmosphere.scale-growth-verification/1',
        evidence='Independent raw-output, trace, input and data-source checks; empirical adequacy is reported separately.',
        verifier_sha256=benchmark.lunar.base.digest(Path(__file__)),native_runs=checks,traces=traces,
        build_producer_sha256=build['producer_sha256'],
        measurements=dict(count=len(rows),missing_polar_axis=sum(r['c_full_mm'] is None for r in rows),
            md5_and_sha256_verified=True,unit_evidence='Author plotting script and original paper, recorded in source register',
            workbook_sha256=benchmark.lunar.base.digest(workbook)),
        mode_zero_parity=products[1]['baseline_ledger_parity'],legacy_case_parity=products[0]['baseline_ledger_parity'])
    benchmark.lunar.boxes.save(ROOT/'atmosphere/snow/results/scale_growth_verification.json',record)
    print('Verified',len(checks),'native runs,',len(traces),'step traces,',len(rows),'measurements',flush=True)
    return record


if __name__=='__main__':run()
