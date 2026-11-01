"""Verify diagnostic forcing bytes and precipitation against native SCALE history."""
from pathlib import Path
import json
import shutil
import subprocess
import tempfile
import numpy as np
from atmosphere.snow import scale_size_diagnosis as diagnosis


def verify(spec):
    l=diagnosis.lunar
    d=l.CACHE/'runs'/spec['name']
    settings=json.loads((d/'settings.json').read_text())
    expected=json.loads((d/'inputs.json').read_text())
    mismatched=[name for name,sha in expected.items() if l.base.digest(d/name)!=sha]
    if mismatched:raise ValueError(f'Input identity mismatch: {d.name} {mismatched}')
    from netCDF4 import Dataset
    with Dataset(d/'history.pe000000.nc') as nc:
        times=np.array(nc['time'][:]);native=np.array(nc['PREC_sd'][:]).mean(axis=(1,2))*1000
    arrivals=np.fromfile(d/'lunar_arrivals.bin',dtype='>f8').reshape(-1,len(l.ARRIVAL_FIELDS))
    arrivals=arrivals[np.argsort(arrivals[:,0])]
    cumulative=np.r_[0,np.cumsum(arrivals[:,2]*arrivals[:,3])]/(settings['horizontal_cells']*settings['dx_m']**2)
    indices=np.searchsorted(arrivals[:,0],np.minimum(times+settings['dt_s'],settings['duration_s'])+1e-7,side='right')
    error=float(np.max(np.abs(native-cumulative[indices])))
    relative=error/max(float(native.max()),1e-15)
    if relative>1e-9:raise ValueError(f'Native precipitation differs: {d.name} {relative}')
    with tempfile.TemporaryDirectory(prefix='terluna-size-reproduce-') as tmp:
        scratch=Path(tmp);(scratch/'lunar').mkdir()
        shutil.copyfile(l.CACHE/'lunar/manifest.json',scratch/'lunar/manifest.json')
        code="import json,sys,numpy;from atmosphere.snow.scale_size_diagnosis import configure;p=configure(**json.loads(sys.argv[1]),cache=sys.argv[2]);print(json.dumps(dict(path=str(p),python=sys.executable,numpy=numpy.__version__)))"
        result=subprocess.run([shutil.which('python'),'-B','-c',code,json.dumps(spec),str(scratch)],cwd=l.ROOT,check=True,capture_output=True,text=True)
        generation=json.loads(result.stdout);recreated=Path(generation['path'])
        names=('init.conf','run.conf','env.txt','particles.pe000000','lunar_forcing.txt')
        matches={name:l.base.digest(d/name)==l.base.digest(recreated/name) for name in names}
        if not all(matches.values()):raise ValueError(f'Physics input recreation failed: {d.name} {matches}')
    return dict(case=d.name,native_history_samples=len(times),input_hash_count=len(expected),
        maximum_precipitation_difference_kg_m2=error,maximum_relative_precipitation_difference=relative,
        reproduced_physics_inputs=matches,generation=generation,
        history_sha256=l.base.digest(d/'history.pe000000.nc'))


if __name__=='__main__':
    checks=[]
    for spec in diagnosis.CASES:
        checks.append(verify(spec));print('verified',spec['name'],flush=True)
    diagnosis.lunar.boxes.save(diagnosis.lunar.CACHE/'size_verification.json',checks)
