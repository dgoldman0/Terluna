"""Independent numerical and artifact-consistency checks; not physical validation."""
from pathlib import Path
import hashlib,json
import numpy as np
import model as m

root=Path(__file__).resolve().parent
j=json.loads((root/'results/results.json').read_text())
w=np.linspace(.4,2.5,1000);n=1.5;d=1.2
r,t=m.stack_rt(w,[np.full(w.shape,n,dtype=complex)],[d])
R=((n-1)/(n+1))**2
analytic=1/(1+4*R/(1-R)**2*np.sin(2*np.pi*n*d/w)**2)
err=float(np.max(abs(t-analytic)))
assert err<1e-12
_,_,coarse=m.holding_acceleration(78000e3,count=2048)
_,_,fine=m.holding_acceleration(78000e3,count=8192)
force_error=float(abs(coarse.max()/fine.max()-1))
assert force_error<1e-5
cases=[j['orbits']['reference']]+j['orbits']['energy_storage_sensitivity']
closure=[]
for z in cases:
    if not z['closed']:continue
    components=j['orbits']['reference']['area_m2']*.05+z['power_system_mass_kg']+z['fuel_buffer_mass_kg']+z.get('energy_storage_mass_kg',0)
    e=abs(components/z['total_mass_kg']-1);closure.append(e);assert e<1e-12
for name,h in j['input_hashes'].items():assert hashlib.sha256((root/'sources'/name).read_bytes()).hexdigest()==h
output={'scope':'Numerical and data consistency only; no physical validation of future systems.',
        'analytic_single_slab_T_max_error':err,
        'holding_peak_relative_grid_error_2048_vs8192':force_error,
        'maximum_component_mass_closure_relative_error':max(closure),
        'verified_input_hashes':len(j['input_hashes'])}
(root/'results/validation.json').write_text(json.dumps(output,indent=2))
print(json.dumps(output,indent=2))
