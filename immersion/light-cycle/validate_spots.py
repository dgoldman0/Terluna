"""Run independent transport spot checks and record, rather than hide, errors."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from solver import ROOT,EARTH,MOON,MOON_ZERO
from render_data import evaluate_file
from reference_mc import check

def run(photons=500000):
    records=[];paths={a.name:ROOT/'work'/f'{a.name}_standard_20.npz' for a in [EARTH,MOON,MOON_ZERO]}
    # Four viewing directions through three twilight states, safely away from limb.
    targets=[(3.,90.,0.),(-3.,30.,0.),(-9.,10.,0.),(-3.,10.,180.)]
    for atm in [EARTH,MOON]:
        p=paths[atm.name];d=np.load(p);lam=d['lam'];m=json.loads(str(d['meta']))
        for sun,elev,az in targets:
            det=evaluate_file(p,[sun],[elev],[az],nshell=100)[0,0,0]
            for wave in [440.,560.,680.]:
                rec=check(atm,wave,sun,elev,az,photons=photons,seed=83+int(wave))
                j=int(np.argmin(abs(lam-wave)));v=float(det[j]);se=rec['standard_error'];mean=rec['mean']
                rec.update(deterministic=v,relative_difference=(v-mean)/max(mean,1e-30),
                  difference_in_sampling_standard_errors=(v-mean)/max(se,1e-30),
                  source_orders=len(m['history']),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                records.append(rec);print(json.dumps(rec),flush=True)
    out=ROOT/'validation';out.mkdir(exist_ok=True)
    (out/'monte_carlo_spots.json').write_text(json.dumps(records,indent=2))
    # More ray segments with the identical source field, all 22 wavelengths.
    resolution=[]
    for atm in [EARTH,MOON,MOON_ZERO]:
        p=paths[atm.name]
        s=[-12.,-6.,-3.,0.,3.];e=[1.,10.,30.,90.];a=[0.,90.,180.]
        low=evaluate_file(p,s,e,a,nshell=64);high=evaluate_file(p,s,e,a,nshell=128)
        bright=high>1e-7
        rel=np.abs(low-high)/np.maximum(high,1e-30)
        resolution.append(dict(world=atm.name,suns=s,elevations=e,azimuths=a,wavelengths_nm=np.load(p)['lam'].tolist(),
          max_relative_above_1e_7=float(np.max(rel[bright])),p95_relative_above_1e_7=float(np.quantile(rel[bright],.95)),
          max_absolute_W_m2_sr_nm=float(np.max(abs(low-high)))))
    (out/'ray_resolution.json').write_text(json.dumps(resolution,indent=2))
    return records

if __name__=='__main__':run()
