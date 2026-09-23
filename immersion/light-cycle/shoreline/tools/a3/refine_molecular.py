"""Targeted independent A2 discretization diagnostic; original solver unchanged."""
import numpy as np,time,json,math,argparse
import coupled_transport as t
m=t.molecular
ap=argparse.ArgumentParser();ap.add_argument('--quality',default='dense');a=ap.parse_args();c=t.load_inputs();ks=np.array([4,10]);R=c['R'];top=c['top'];start=time.monotonic()
if a.quality=='dense':nr,nmu,nphi,nshell,npath,na=64,32,32,192,512,2.
else:nr,nmu,nphi,nshell,npath,na=40,24,24,128,384,3.
rg=R+top*np.linspace(0,1,nr)**2
ag=np.deg2rad(np.unique(np.r_[np.arange(-90,-30,4.),np.arange(-30,30.001,na),np.arange(32,90.001,na),90.]))
srg=R+top*np.linspace(0,1,100)**2;sag=np.deg2rad(np.unique(np.r_[np.arange(-90,-30,2.),np.arange(-30,30.001,.25),np.arange(30.5,90.001,.5)]))
beta,sigma,solar=c['beta'][ks],c['sigma'][ks],c['solar'][ks]
beam=m.solar_table(srg,sag,R,R+top,c['rows'],c['oz'],beta,sigma,solar,npath)
mu,mw,cp,sp,dp=m.angular_quadrature(nmu,nphi);paths=m.build_paths(R,top,rg,mu,beta,sigma,c['rows'],c['oz'],nshell,c['profile']['sounding']['end_m'])
prev=np.zeros((len(rg),len(ag),2,6));total=np.zeros_like(prev);hist=[]
for order in range(1,101):
 cur=m.transport(prev,beam,rg,ag,srg,sag,mu,mw,cp,sp,dp,*paths,order==1);total+=cur
 frac=float(np.max(cur[:,:,:,0])/np.max(total[:,:,:,0]));hist.append(frac);prev=cur
 if order%5==0:print(order,frac,'time',time.monotonic()-start,flush=True)
 if order>=10 and frac<1e-5:break
shells=m.shell_altitudes(top,384,c['oz'],c['profile']['sounding']['end_m']);out=m.eval_ray(math.pi/2,0,math.pi/4,R,top,c['rows'],c['oz'],beta,sigma,beam,rg,ag,srg,sag,total-prev,shells)
j=dict(quality=a.quality,grid=dict(nr=nr,nmu=nmu,nphi=nphi,nshell=nshell,npath=npath),lam=c['lam'][ks].tolist(),zenith=out.tolist(),orders=order,history=hist,elapsed_s=time.monotonic()-start)
print(json.dumps(j),flush=True);(t.ROOT/f'validation/a3/clear-refine-{a.quality}.json').write_text(json.dumps(j,indent=2)+'\n')
