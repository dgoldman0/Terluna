import concurrent.futures,json,math,time
from pathlib import Path
import numpy as np
import coupled_transport as t
import deterministic as dt

def run():
 c=t.load_inputs();checks=[];o=np.array([0.,2.,0.]);up=np.array([0.,1.,0.]);sun=np.array([2**-.5,2**-.5,0.]);k=10;start=time.monotonic()
 def add(name,passed,**values):checks.append(dict(name=name,passed=bool(passed),**values))
 for label,g,cl in [('zero cloud',True,False),('zero gas',False,True)]:
  a=t.estimate(*t.args(c,o,up,sun,k,g,cl,1),8192,45670)
  b=t.estimate(*t.args(c,o,up,sun,k,g,cl,3),8192,45670)
  add('clear-source continuation reduces exactly to full transport with '+label,np.array_equal(a[0],b[0]) and np.array_equal(a[1],b[1]) and a[2]+b[2]==0)
 def estimate(seed):return t.estimate(*t.args(c,o,up,sun,k,True,True,1),1024,seed)
 serial=[estimate(s) for s in [404,606,808,1010]]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:parallel=list(pool.map(estimate,[404,606,808,1010]))
 add('independent-seed batches reproduce exactly across serial and threaded dispatch',all(np.array_equal(a[0],b[0]) and np.array_equal(a[1],b[1]) and a[2]==b[2] for a,b in zip(serial,parallel)))
 # Normalized isotropic moment field should produce the same radiance in every
 # outgoing direction, including the solar pole where the tangent basis collapses.
 moment=np.zeros((2,2,1,6));I=2.;moment[:,:,:,0]=4*math.pi*I;moment[:,:,:,1:4]=4*math.pi*I/3
 rg=np.array([c['R'],c['R']+c['top']]);ag=np.array([-math.pi/2,math.pi/2]);values=[]
 for ds,ss in [(up,up),(np.array([1.,0.,0.]),up),(up,sun),(np.array([1.,0.,0.]),sun)]:values.append(t.molecular_source(o,ds,ss,c['R'],rg,ag,moment,0))
 add('A2 moment-source adapter preserves isotropic radiance including pole cases',max(abs(v-I) for v in values)<1e-12,values=values,expected=I)
 A2=dt.a2_clear(c,o,up,sun)
 for wave in [400,460,580,700]:
  k=int(np.where(c['lam']==wave)[0][0]);m,se,f=t.estimate(*t.args(c,o,up,sun,k,True,False,2),131072,62338+k)
  add(f'A2 source adapter recovers saved clear field at {wave} nm',abs(m[0]-A2[k])<=6*se[0]+.002*A2[k] and f==0,mean=float(m[0]),standard_error=float(se[0]),A2=float(A2[k]))
 n,w=np.polynomial.legendre.leggauss(24)
 for ri,d in enumerate([up,np.array([.25,1.,-.15])]):
  d=d/np.linalg.norm(d);air,oz,blocked=dt.gas_columns(o,d,c['R'],c['top'],c['rows'],c['oz'],c['shells'],n,w);cl=dt.cloud_depth(o,d,c['lo'],c['hi'],c['dims'],c['ext'],c['ice'])
  for wave in [400,580,700]:
   k=int(np.where(c['lam']==wave)[0][0]);expected=math.exp(-air*c['beta'][k]-oz*t.DU*c['sigma'][k]-cl)
   m,se,fail=t.estimate_transmission(o,d,sun,c['R'],c['top'],c['rows'],c['oz'],c['beta'][k],c['sigma'][k],c['edges'],c['air_bounds'],c['oz_bounds'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['cmax'],True,True,131072,400042+ri*100+k)
   add(f'combined transmission equals independent gas times cloud attenuation ray {ri}, {wave} nm',abs(m-expected)<=6*se+1e-6 and fail==0,mean=m,standard_error=se,expected=expected,unresolved=int(fail))
 out=dict(schema='open-moon-a3-additional-checks/1',checks=checks,passed=all(v['passed'] for v in checks),elapsed_s=time.monotonic()-start)
 (t.ROOT/'validation/a3/extension-tests.json').write_text(json.dumps(out,indent=2)+'\n')
 for ch in checks:print(json.dumps(ch),flush=True)
 return out
if __name__=='__main__':
 out=run();raise SystemExit(0 if out['passed'] else 1)
