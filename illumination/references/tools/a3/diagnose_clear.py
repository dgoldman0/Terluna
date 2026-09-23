"""Independent clear-medium samples and roulette sensitivity at three wavelengths."""
import numpy as np,time,json,math,concurrent.futures
import coupled_transport as t
import deterministic as dt
c=t.load_inputs();o=np.array([0.,2.,0.]);d=np.array([0.,1.,0.]);sun=np.array([2**-.5,2**-.5,0.]);A2=dt.a2_clear(c,o,d,sun)
def run(job):
 k,n,seed,mode,rr=job;m,se,f=t.estimate(*t.args(c,o,d,sun,k,True,False,mode),n,seed,0.,rr)
 return dict(wavelength_nm=float(c['lam'][k]),samples=n,seed=seed,mode=mode,roulette=rr,mean=float(m[0]),SE=float(se[0]),A2=float(A2[k]),unresolved=f)
jobs=[]
for w in [400,460,580]:
 k=int(np.where(c['lam']==w)[0][0])
 for mode,rr,seed in [(1,.95,11003),(1,.99,88327),(2,.95,33211)]:jobs.append((k,1048576,seed+k,mode,rr))
start=time.monotonic();out=[]
# Warm cached code before threads.
run((10,2,1,1,.95))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for r in pool.map(run,jobs):out.append(r);print(json.dumps(r),flush=True)
p=t.ROOT/'validation/a3/clear-diagnosis.json';p.write_text(json.dumps(dict(schema='open-moon-a3-clear-diagnosis/1',profile_sha256=c['profile_sha256'],source_sha256=t.source_hashes(),results=out,elapsed_s=time.monotonic()-start),indent=2)+'\n')

if any(r["unresolved"] for r in out):raise SystemExit("Unresolved histories require review")
