"""Independent input and native-history checks for completed lunar cases.

Requires netCDF4 only when invoked. Native time labels precede the microphysics
update by one nominal step, including a shortened final physical step.
"""
from pathlib import Path
import json
import tempfile
import shutil
import subprocess
import numpy as np
from atmosphere.snow import scale_lunar as l, scale_lunar_comparison as comparison


def verify(case):
    d=l.CACHE/'runs'/case['name'];settings=json.loads((d/'settings.json').read_text())
    expected=json.loads((d/'inputs.json').read_text())
    mismatched=[name for name,sha in expected.items() if l.base.digest(d/name)!=sha]
    if mismatched:raise ValueError(f'Input identity mismatch in {d.name}: {mismatched}')
    rows=np.fromfile(d/'lunar_arrivals.bin',dtype='>f8').reshape(-1,len(l.ARRIVAL_FIELDS))
    from netCDF4 import Dataset
    with Dataset(d/'history.pe000000.nc') as nc:
        times=np.array(nc['time'][:]);precip=np.array(nc['PREC_sd'][:]).mean(axis=(1,2))*1000
    actual_times=np.minimum(times+settings['dt_s'],settings['duration_s'])
    order=np.argsort(rows[:,0]);rows=rows[order]
    cumulative=np.r_[0,np.cumsum(rows[:,2]*rows[:,3])]/(settings['horizontal_cells']*settings['dx_m']**2)
    expected_precip=cumulative[np.searchsorted(rows[:,0],actual_times+1e-7,side='right')]
    error=float(np.max(np.abs(precip-expected_precip)))
    relative=error/max(float(precip.max()),1e-15)
    if relative>1e-9:raise ValueError(f'Native precipitation mismatch: {d.name}, {relative}')
    # Recreate native inputs with the final adapter. Config identity is stronger
    # than assuming a long-lived Python process reloaded later analysis edits.
    scratch=Path(tempfile.mkdtemp(prefix='terluna-lunar-reproduce-'))
    try:
        (scratch/'runs').mkdir()
        build_dir='lunar' if case.get('collection','native')=='native' else 'lunar_sticking'
        (scratch/build_dir).mkdir()
        shutil.copyfile(l.CACHE/build_dir/'manifest.json',scratch/build_dir/'manifest.json')
        generator=shutil.which('python')
        code="import json,sys,numpy;from atmosphere.snow import scale_lunar as l;p=l.configure(**json.loads(sys.argv[1]),cache=sys.argv[2]);print(json.dumps(dict(path=str(p),python=sys.executable,numpy=numpy.__version__)))"
        native=subprocess.run([generator,'-B','-c',code,json.dumps({**case,'dt':settings['dt_s']}),str(scratch)],cwd=l.ROOT,check=True,capture_output=True,text=True)
        generation=json.loads(native.stdout);q=Path(generation['path'])
        names=('init.conf','run.conf','env.txt','particles.pe000000','lunar_forcing.txt')
        reproduced={name:l.base.digest(q/name)==l.base.digest(d/name) for name in names}
        if not all(reproduced.values()):raise ValueError(f'Native input reproduction failed: {d.name}: {reproduced}')
    finally:shutil.rmtree(scratch)
    return dict(case=d.name,input_hashes_checked=len(expected),native_history_samples=len(times),
       maximum_precipitation_difference_kg_m2=error,maximum_relative_precipitation_difference=relative,
       reproduced_native_inputs=reproduced,input_reproduction_environment=generation,history_sha256=l.base.digest(d/'history.pe000000.nc'))


if __name__=='__main__':
    rows=[]
    for spec in comparison.NATIVE_CASES:
        path=l.CACHE/'runs'/spec['name']/'analysis.json'
        if path.exists():rows.append(verify(spec));print('verified',spec['name'],flush=True)
    l.boxes.save(l.CACHE/'lunar_verification.json',rows)
