"""Independent input, native-output and refinement checks for SCALE sensitivities."""
from pathlib import Path
import json
import math
import shutil
import tempfile
import numpy as np
from atmosphere.snow import scale_crystal_growth as growth
from atmosphere.snow import scale_separation as separation
from atmosphere.snow.scale_size_verify import verify as verify_column

lunar=growth.lunar


def growth_checks(spec):
    from netCDF4 import Dataset
    d=lunar.CACHE/'runs'/spec['name']
    r=json.loads((d/'growth_analysis.json').read_text());s=r['settings']
    hashes=r['input_sha256']
    if any(lunar.base.digest(d/n)!=h for n,h in hashes.items()):
        raise ValueError('Changed laboratory input')
    execution=r['execution'];expected=s['build']['build']['executable_sha256']
    if len(execution)!=2 or any(e['exit_code']!=0 or e['sha256']!=expected[e['executable']] for e in execution):
        raise ValueError('Incomplete or unidentified native growth execution')
    with tempfile.TemporaryDirectory(prefix='terluna-growth-reproduce-') as tmp:
        cache=Path(tmp);(cache/'lunar').mkdir();(cache/'limitations_sources').mkdir()
        shutil.copyfile(lunar.CACHE/'lunar/manifest.json',cache/'lunar/manifest.json')
        (cache/growth.SOURCE_PDF).symlink_to(lunar.CACHE/growth.SOURCE_PDF)
        copy=growth.configure(cache=cache,**spec)
        names=('init.conf','run.conf','env.txt','particles.pe000000','lunar_forcing.txt')
        matches={n:lunar.base.digest(d/n)==lunar.base.digest(copy/n) for n in names}
        if not all(matches.values()):raise ValueError(f'Laboratory input reproduction failed: {matches}')
    with Dataset(d/'history.pe000000.nc') as nc:
        pressure=np.asarray(nc['PRES'][:])[1:]
        temp=np.asarray(nc['T'][:])[1:]
        rh=np.asarray(nc['RH'][:])[1:]
    if np.max(np.abs(pressure-growth.PRESSURE))>1e-6 or np.max(np.abs(temp-growth.TEMPERATURE))>1e-9:
        raise ValueError('Laboratory air state was not maintained')
    if np.max(np.abs(rh-100))>.1:
        raise ValueError('Liquid saturation differs by more than 0.1 percentage point')
    w=r['water']
    if w['ice_deposition_relative_residual']>1e-11 or w['maximum_used_fraction_of_growth_ceiling']>1e-4:
        raise ValueError('Water accounting or an active growth ceiling affects the laboratory test')
    if len(r['snapshots'])!=28 or [v['age_min'] for v in r['snapshots']]!=list(range(3,31)):
        raise ValueError('Incomplete measured-age comparison')
    if max(v['maximum_relative_between_particle_difference'] for v in r['snapshots'])>1e-12:
        raise ValueError('Identical isolated crystals differ')
    initial=r['snapshots'][0];final=r['snapshots'][-1]
    observed0=initial['observed_fit'];observed1=final['observed_fit']
    slope=lambda a,b:math.log(b['thickness_m']/a['thickness_m'])/math.log(b['diameter_m']/a['diameter_m'])
    return dict(case=d.name,physics_input_matches=matches,execution_identity_checked=True,
                snapshot_count=len(r['snapshots']),water=w,
                air_ranges_after_first_step=r['air_ranges_after_first_step'],
                simulated_integrated_axis_growth_ratio=slope(initial,final),
                observed_integrated_axis_growth_ratio=slope(observed0,observed1),
                history_sha256=lunar.base.digest(d/'history.pe000000.nc'),
                analysis_sha256=lunar.base.digest(d/'growth_analysis.json'))


def zero_parity():
    from netCDF4 import Dataset
    old=lunar.CACHE/'runs/size_moon_baseline';new=lunar.CACHE/'runs/separation_zero_moon'
    names=('lunar_arrivals.bin','lunar_ledger.txt','collision_observer.txt')
    files={n:dict(sha256=lunar.base.digest(old/n),identical=lunar.base.digest(old/n)==lunar.base.digest(new/n)) for n in names}
    arrays={}
    with Dataset(old/'history.pe000000.nc') as a,Dataset(new/'history.pe000000.nc') as b:
        if set(a.variables)!=set(b.variables):raise ValueError('Native history fields differ')
        for k in a.variables:
            arrays[k]=bool(np.array_equal(np.asarray(a[k][:]),np.asarray(b[k][:]),equal_nan=True))
    if not all(v['identical'] for v in files.values()) or not all(arrays.values()):
        raise ValueError('Zero-separation numerical parity failed')
    return dict(files=files,numerical_history_arrays=arrays,
                note='NetCDF global build metadata differ after the progress commit; numerical arrays, arrivals and ledgers match exactly.')


def run():
    columns=[verify_column({k:v for k,v in spec.items() if k!='separation'}) for spec in separation.CASES]
    for check in columns:
        d=lunar.CACHE/'runs'/check['case']
        settings=json.loads((d/'settings.json').read_text())
        execution=json.loads((d/'execution.json').read_text())
        expected=settings['build']['build']['executable_sha256']
        if len(execution)!=2 or any(e['exit_code']!=0 or e['sha256']!=expected[e['executable']] for e in execution):
            raise ValueError('Incomplete or unidentified native column execution')
        check['execution_identity_checked']=True
    crystals=[growth_checks(s) for s in growth.CASES]
    g=json.loads(growth.OUT.read_text())['cases']
    pairs=[]
    for a,b in [('growth_earth_observed_dt10','growth_earth_observed_dt5'),
                ('growth_earth_observed_dt5','growth_earth_observed_dt1')]:
        x=g[a]['snapshots'];y=g[b]['snapshots']
        final={k:y[-1][k]/x[-1][k]-1 for k in ('mass_kg','diameter_m','thickness_m')}
        trajectory={k:max(abs(v[k]/u[k]-1) for u,v in zip(x,y)) for k in final}
        pairs.append(dict(left=a,right=b,endpoint_relative_change=final,
                          maximum_relative_change_over_sampled_trajectory=trajectory,
                          evidence='Measured numerical sensitivity; no asymptotic convergence or observational acceptance inferred.'))
    result=dict(schema='terluna.atmosphere.scale-limitations-verification/1',
                verifier_sha256=lunar.base.digest(Path(__file__)),columns=columns,crystals=crystals,
                zero_separation_parity=zero_parity(),growth_refinements=pairs)
    lunar.boxes.save(lunar.CACHE/'limitations_verification.json',result)
    print('Verified',len(columns),'columns and',len(crystals),'growth cases',flush=True)
    return result


if __name__=='__main__':run()
