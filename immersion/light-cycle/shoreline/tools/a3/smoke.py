import coupled_transport as t
import numpy as np,time,json
c=t.load_inputs();print('inputs',c['dims'],c['profile_sha256']);u=np.array([2**-.5,2**-.5,0.]);o=np.array([0.,2.,0.]);d=np.array([0.,1.,0.]);start=time.monotonic()
for w in [460,580,700]:
 k=int(np.where(c['lam']==w)[0][0])
 for name,gas,cl,mode in [('clear',True,False,1),('coupled',True,True,1),('frozen',True,True,2)]:
  mean,se,fail=t.estimate(*t.args(c,o,d,u,k,gas,cl,mode),4096,384838+k)
  print(json.dumps({'wave':w,'mode':name,'mean':mean.tolist(),'SE':se.tolist(),'fail':fail,'elapsed':time.monotonic()-start}),flush=True)
