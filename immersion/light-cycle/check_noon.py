from solver import *
from render_data import evaluate_file
from reference_mc import check
import json, numpy as np, hashlib
if __name__=='__main__':
 records=[]
 for atm in [EARTH,MOON]:
  p=ROOT/'work'/f'{atm.name}_standard_20.npz';d=np.load(p)
  det=evaluate_file(p,[90.],[30.],[0.],nshell=128)[0,0,0]
  for wave in [440.,560.,680.]:
   r=check(atm,wave,90.,30.,0.,photons=100000,seed=84+int(wave))
   val=float(det[np.argmin(abs(d['lam']-wave))]);r.update(deterministic=val,relative_difference=(val-r['mean'])/r['mean'],sigma_difference=(val-r['mean'])/r['standard_error'],source_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
   print(json.dumps(r),flush=True);records.append(r)
  (ROOT/'validation'/'noon_monte_carlo.json').write_text(json.dumps(records,indent=2))
