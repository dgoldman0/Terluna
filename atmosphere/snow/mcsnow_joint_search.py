"""Local joint calibration; only the separate pointwise audit can accept a result.

Every sampled point contributes excess above its individual limit. A bounded
local optimizer reaching its stopping criterion is not physical acceptance.
"""
import copy,json,math,time
import numpy as np
from scipy.optimize import least_squares
from atmosphere.snow import mcsnow_joint_growth as g
from atmosphere.snow import mcsnow_joint_acceptance as a


def run(temp,label,initial,setup,max_nfev=80):
 g.build();path=g.CACHE/'trajectory_fits'/label;path.mkdir(parents=True,exist_ok=True)
 spec=dict(temperature_C=temp,setup=setup,initial=list(initial),max_nfev=max_nfev,dt_s=1/16,
  baseline_sha256=a.BASELINE_SHA256,build_sha256=g.native.sha(g.CACHE/'build/build.json'),
  objective='Sum of squared pointwise error excess beyond min(10%, baseline error), including all four trajectory quantities, seed shapes and supplied-state rates. Search ages are a subset; full acceptance follows.')
 if (path/'specification.json').exists() and json.loads((path/'specification.json').read_text())!=spec:raise ValueError('Search changed')
 g.model.save(path/'specification.json',spec)
 if (path/'result.json').exists():return json.loads((path/'result.json').read_text())
 baseline=[p for p in a.baseline_inventory()['points'] if p['temperature_C']==temp];index={p['id']:p for p in baseline}
 ages=(3.,4.,5.,6-1/16/60,6.,7.,10.,15.,20.,25.,30.)
 entries=g.entries_for(temp)
 for e in entries:
  e['sample_ages_min']=[age for age in e['sample_ages_min'] if age in ages]
  e['acceptance_family']='trajectory';e['acceptance_start']=e['initialization']
 for old in g.baseline_fit(temp)['seed_shape_sensitivity']:
  phi=old['initial_aspect_ratio']
  if phi==1:continue
  e=copy.deepcopy(entries[0]);d=e['configuration']['initial_equatorial_m'];e['sample_ages_min']=[r['age_min'] for r in old['assessment']['rows']]
  e['configuration'].update(initial_equatorial_m=d*phi**(-1/3),initial_polar_m=d*phi**(2/3))
  e['acceptance_family']='seed_sensitivity';e['acceptance_start']='seed_phi='+format(phi,'.17g');entries.append(e)
 def mat(ps):return np.array([[*p['rates'],p['effective_axis_ratio'],p['deposition_density_kg_m3']] for p in ps])
 oldp=g.baseline_fit(temp)['fixed_state_probes'];pt=np.array([[*p['target_rates'],p['target_effective_axis_ratio'],p['target_deposition_density_kg_m3']] for p in oldp]);pa=np.minimum(.1,abs(mat(oldp)/pt-1))
 initial=np.array(initial);free=[0,1,3,4,5,6,7,8,9,10,11,12]
 def expand(z):
  x=initial.copy();x[free]=z;return x
 with g.batch(entries,setup,dt=1/16) as b,g.batch(g.probe_entries(temp),setup) as pb:
  allowed=[];target=g.tuning.targets(b)
  for i,age in b.rows:
   e=entries[i];limits=[]
   for q in a.QUANTITIES:
    p=index[a.point_id(dict(family=e['acceptance_family'],temperature_C=temp,start=e['acceptance_start'],age_min=age,quantity=q))]
    limits.append(min(.1,abs(p['value']/p['target']-1)))
   allowed.append(limits)
  allowed=np.array(allowed);n=target.size+pt.size;started=time.monotonic();calls=0;best=math.inf;best_x=initial.copy();last=started
  def residual(z):
   nonlocal calls,best,best_x,last
   x=expand(z);r=b.evaluate(x);pr=pb.evaluate(x);calls+=1
   if r['error'] or pr['error']:return np.full(n,100.)
   rel=g.tuning.reported(b,r)/target-1;prel=mat(g.probe_values(pr,temp))/pt-1
   vec=np.r_[(np.sign(rel)*np.maximum(abs(rel)-allowed,0)).ravel(),(np.sign(prel)*np.maximum(abs(prel)-pa,0)).ravel()]
   loss=float(vec@vec)
   if loss<best:best=loss;best_x=x.copy()
   now=time.monotonic()
   if now-last>30:
    print(label,'calls',calls,'best_violation_loss',best,'elapsed_s',round(now-started),flush=True);last=now
   return vec
  fit=least_squares(residual,np.clip(initial[free],1e-9,1-1e-9),bounds=(np.zeros(len(free)),np.ones(len(free))),diff_step=.005,x_scale='jac',max_nfev=max_nfev,ftol=1e-8,xtol=1e-7,gtol=1e-8)
  result=dict(label=label,specification=spec,coordinates=best_x.tolist(),physical_controls=g.decode(best_x,setup).tolist(),objective=best,
   message=str(fit.message),optimizer_success=bool(fit.success),calls=calls,elapsed_s=time.monotonic()-started,
   evaluation_directories=[str(b.path),str(pb.path)])
 g.model.save(path/'result.json',result)
 full=g.evaluate_fit(best_x,temp,setup);other=g.baseline_fit(-14.4 if temp==-5.3 else -5.3)
 inv=a.inventory(full,other) if temp==-5.3 else a.inventory(other,full)
 result['full_fit']=full;result['acceptance_with_other_temperature_unchanged']=a.compare(inv)
 g.model.save(path/'result.json',result)
 print(label,'done',result['message'],{k:result['acceptance_with_other_temperature_unchanged'][k] for k in ['accepted','above_ten_percent','regression_count']},flush=True)
 return result

if __name__=='__main__':
 for t in [-5.3,-14.4]:
  p=g.CACHE/'probe_fits'/('state_'+g.diagnosis.key(t))/'result.json';x=np.array(json.loads(p.read_text())['coordinates'])
  if t==-14.4:x[6]=20/150;x[7]=0.;x[10]=0.
  run(t,'joint_state_'+g.diagnosis.key(t),x,g.STATE)

 setup=dict(g.NEUTRAL,normal_compact=True,original_best=.85,prism_width_fraction=3/math.pi)
 run(-5.3,'joint_preserved_density_m5p3',g.baseline_coordinates(-5.3),setup,max_nfev=160)
