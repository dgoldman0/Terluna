"""Independent B1 geometry/boundary checks; retains every measured failure."""
from __future__ import annotations
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from numba import njit
import geometry as g
import scene_transport as b
import render_scene as render
import primary

@njit(cache=True)
def brute(p,d,R,V,limit=1e30):
    hit=min(limit,g.sphere_hit(p,d,R));n=np.zeros(3)
    for z in range(V.shape[0]-1):
        for x in range(V.shape[1]-1):
            a=V[z,x];v=V[z,x+1];c=V[z+1,x];e=V[z+1,x+1]
            h,k=g.triangle(p,d,a,v,c,hit)
            if h<hit:hit=h;n=k
            h,k=g.triangle(p,d,v,e,c,hit)
            if h<hit:hit=h;n=k
    return hit,n

def path_args(c,V,glo,ghi,step,origin,d,sun,lam=550.,radius=0.,gas=False,cloud=False,terrain=False,scale=1.,water=False):
    return (np.array(origin,float),np.array(d,float),np.array(sun,float),radius,lam,c['R'],c['top'],c['rows'],c['oz'],1.24062e-6*(lam/1000)**-4,float(np.interp(lam,c['lam'],c['sigma'])),c['edges'],c['air_bounds'],c['oz_bounds'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['cmax'],V,glo,ghi,step,gas,cloud,terrain,scale,water)

@njit(cache=True)
def samples(args,n,seed,boundary=0.):
    np.random.seed(seed);s=s2=0.;f=np.zeros(8,np.int64)
    for j in range(n):
        v,sv,fail=b.path_sample(*args,.95,boundary,4096)
        s+=v;s2+=v*v;f[fail]+=1
    return s/n,math.sqrt(max(0.,(s2-s*s/n)/max(1,n-1)/n)),f

@njit(cache=True)
def cosine_test(n):
    np.random.seed(571);v=np.array([0.,1.,0.]);s=0.;mn=1.
    for j in range(n):
        d=b.cosine_direction(v);s+=d[1];mn=min(mn,d[1])
    return s/n,mn

@njit(cache=True)
def compare_a3(a,ar,n):
    mx=0.
    for i in range(n):
        np.random.seed(2319+i);v,_,f=b.path_sample(*a)
        np.random.seed(2319+i);w,gf=b.t.path_sample(*ar)
        mx=max(mx,abs(v-w.sum()))
        if f or gf:return mx,False
    return mx,True

def main():
    out=[]
    def record(name,passed,**details):out.append(dict(name=name,passed=bool(passed),**details));print(name,passed,details,flush=True)
    c=b.t.load_inputs();V,glo,ghi,step=g.make_coast(c['R'],25)
    rng=np.random.default_rng(190);errors=[]
    for k in range(150):
        o=np.array([rng.uniform(-25000,25000),rng.uniform(30,13000),rng.uniform(-25000,25000)])
        target=np.array([rng.uniform(-18000,18000),-300.,rng.uniform(-18000,18000)]);d=target-o;d/=np.linalg.norm(d)
        h,normal,mat=g.surface(o,d,c['R'],V,glo,ghi,step);ref,_=brute(o,d,c['R'],V)
        errors.append(abs(h-ref))
    record('DDA versus brute-force triangle intersection',max(errors)<1e-6,max_distance_error_m=max(errors),rays=150)
    hit=g.sphere_hit(np.array([0.,2.,0.]),np.array([0.,-1.,0.]),c['R']);record('Spherical sea intercept',abs(hit-2)<1e-7,distance_m=hit)
    record('No self intersection for outward sea ray',g.sphere_hit(np.array([0.,.02,0.]),np.array([0.,1.,0.]),c['R'])>1e20)
    F0=b.fresnel(1.);F1=b.fresnel(0.);exact=((1.333-1)/(1.333+1))**2
    record('Fresnel normal and grazing limits',abs(F0-exact)<1e-15 and F1==1.,normal=F0,grazing=F1)
    vals=np.array([b.fresnel(x) for x in np.linspace(0,1,1000)]);record('Passive Fresnel partition',np.all((vals>=0)&(vals<=1)))
    mu,mn=cosine_test(100000);record('Cosine hemisphere sampler',abs(mu-2/3)<.003 and mn>=0,mean_cosine=mu)
    for mat in [2,3,4,5]:
        rr=np.array([b.reflectance(mat,l,1.) for l in np.arange(360,831,10)]);record('Passive spectral material '+str(mat),rr.min()>=0 and rr.max()<=1,minimum=float(rr.min()),maximum=float(rr.max()))
    a=path_args(c,V,glo,ghi,step,[0,2,0],[0,-1,0],[0,1,0])
    mean,se,f=samples(a,16,33,boundary=1.);rho=b.reflectance(1,550,1.)
    record('Lambertian surface under unit isotropic sky',abs(mean-rho)<1e-12 and not f[1:].any(),mean=mean,expected=rho)
    for elev in [90,45,6]:
        sun=[math.cos(math.radians(elev)),math.sin(math.radians(elev)),0]
        a=path_args(c,V,glo,ghi,step,[0,2,0],[0,-1,0],sun)
        mean,se,f=samples(a,32,313);expected=rho/math.pi*math.sin(math.radians(elev))
        record('Vacuum direct Lambertian '+str(elev),abs(mean-expected)<1e-12 and not f[1:].any(),mean=mean,expected=expected)
    a=path_args(c,V,glo,ghi,step,[0,2,0],[0,-1,0],[0,1,0],water=True)
    mean,se,f=samples(a,16,42,boundary=1.)
    record('Deep water under unit isotropic sky',abs(mean-F0)<1e-12 and not f[1:].any(),mean=mean,expected=F0)
    radius=math.radians(.26667);solid=2*math.pi*(1-math.cos(radius));source=1/(math.pi*math.sin(radius)**2)
    record('Finite solar disk normal irradiance',abs(source*solid*(1+math.cos(radius))*.5-1)<1e-11)
    V,glo,ghi,step=g.make_coast(c['R']);o=[0,80,-6500];sun=np.array([.3,.7,.64]);sun/=np.linalg.norm(sun)
    for direction in [[0,1,0],[.3,.5,.8124]]:
        d=np.array(direction,float);d/=np.linalg.norm(d);band=19
        a=path_args(c,V,glo,ghi,step,o,d,sun,lam=c['lam'][band],gas=True,cloud=True,scale=0.)
        ar=list(b.t.args(c,o,d,sun,band));ar[3]=1.
        err,done=compare_a3(a,tuple(ar),1000)
        record('Unchanged A3 black-boundary path reduction '+str(direction),err<1e-12 and done,maximum_absolute_path_difference=err)
    a=path_args(c,V,glo,ghi,step,o,[0,1,0],sun,lam=550.,radius=radius,gas=True,cloud=True,terrain=True,water=True)
    mean,se,f=samples(a,20000,51631)
    record('Mixed medium surface paths complete',not f[1:].any(),mean=mean,standard_error=se,completion_codes=f.tolist())
    # Re-bin measured A3 full spectra to quantify spectral grouping alone.
    benchmark=json.loads(b.locate('data/a3/generated/benchmark.json').read_text())
    data=json.loads(b.locate('data/a1/spectral-inputs.json').read_text());dw=np.full(48,10.);dw[[0,-1]]=5.;W=683*np.array(data['cie1931_xyz'])*dw[:,None]
    group_errors=[]
    for case in benchmark['cases']:
        spec=np.array(case['treatments']['coupled_mc']['mean'])[:,0];transfer=spec/c['solar'];lam=np.array([np.mean(c['lam'][ix]) for ix in np.array_split(np.arange(48),12)])
        collapsed=np.array([(W[ix]*c['solar'][ix,None]).sum(axis=0) for ix in np.array_split(np.arange(48),12)])
        approximate=np.interp(lam,c['lam'],transfer)@collapsed;exact=spec@W
        group_errors.append(dict(case=case['id'],Y_relative=float((approximate[1]-exact[1])/exact[1]),XYZ_L1=float(abs(approximate-exact).sum()/exact.sum())))
    record('12-group spectral reduction on saved A3 rays',max(abs(x['Y_relative']) for x in group_errors)<.03,comparisons=group_errors,scope='selected A3 rays; Monte Carlo noise and changed B1 surfaces are separate')
    sea=np.array([[g.sea_y(v[0],v[2],c['R']) for v in row] for row in V])
    heights=V[:,:,1]-sea;edge_heights=np.r_[heights[0],heights[-1],heights[:,0],heights[:,-1]]
    record('Closed terrain boundary below sea',edge_heights.max()<0.,maximum_edge_height_m=float(edge_heights.max()))
    sun=np.array([0.,.5,math.sqrt(.75)])
    aa=path_args(c,V,glo,ghi,step,[0,5000,0],sun,sun,radius=radius,scale=0.)
    np.random.seed(32);with_disk=b.path_sample(*aa)[0]
    np.random.seed(32);without_disk=b.path_sample(*aa,.95,0.,4096,True,True)[0]
    record('Primary solar term separated without double counting',abs(with_disk-1/(math.pi*math.sin(radius)**2))<1e-9 and without_disk==0.,with_disk=with_disk,without_disk=without_disk)
    zero=np.zeros_like(c['ext']);d1,_=primary.primary_frame(np.array([0.,5000.,0.]),64,40,60.,0.,30.,0,sun,radius,np.array([0.]),np.array([0.]),c['R'],c['top'],c['rows'],c['oz'],c['shells'],c['lo'],c['hi'],c['dims'],zero,c['ice'],V,glo,ghi,step,False,64)
    omega=(2*math.tan(math.radians(30))/64)**2;flux=float(d1.sum()*omega)
    record('Resolved solar disk vacuum flux',abs(flux-1.)<.03,normal_flux=flux,subpixel_grid=64,scope='selected centred disk quadrature')
    path=b.ROOT/'validation/b1/boundary-tests.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(dict(schema='open-moon-b1-tests/1',tests=out,passed=sum(x['passed'] for x in out),failed=sum(not x['passed'] for x in out)),indent=2)+'\n')
    if any(not x['passed'] for x in out):raise SystemExit(1)
if __name__=='__main__':main()
