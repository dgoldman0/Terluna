"""Null-collision multiple-scattering reference for one immutable cloud voxel field.

Analog free-flight sampling, ratio-tracked sunlight, reciprocal two-HG phase
sampling, next-event estimation and unbiased Russian roulette. No exposure or
empirical lighting boost. Black or uniform-radiance external boundaries. The
finite crop, gray cloud optics and excluded molecular atmosphere are explicit.
"""
from __future__ import annotations
import argparse, base64, hashlib, json, math, time
from pathlib import Path
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[2]

@njit(cache=True)
def box(p,d,lo,hi):
    entry,exit=0.,1e99
    for k in range(3):
        if abs(d[k])<1e-14:
            if p[k]<lo[k] or p[k]>hi[k]:return 1.,0.
        else:
            a,b=(lo[k]-p[k])/d[k],(hi[k]-p[k])/d[k]
            entry=max(entry,min(a,b));exit=min(exit,max(a,b))
    return entry,exit

@njit(cache=True)
def voxel(p,lo,hi,dims,ext,ice):
    q=np.empty(3);a=np.empty(3,dtype=np.int64);f=np.empty(3)
    for k in range(3):
        if p[k]<lo[k] or p[k]>hi[k]:return 0.,0.
        q[k]=(p[k]-lo[k])/(hi[k]-lo[k])*(dims[k]-1)
        a[k]=min(dims[k]-2,int(math.floor(q[k])));f[k]=q[k]-a[k]
    s,g=0.,0.
    for z in range(2):
        wz=f[2] if z else 1-f[2]
        for y in range(2):
            wy=f[1] if y else 1-f[1]
            for x in range(2):
                w=wz*wy*(f[0] if x else 1-f[0]);i=((a[2]+z)*dims[1]+a[1]+y)*dims[0]+a[0]+x
                s+=w*ext[i];g+=w*ice[i]
    return s,g

@njit(cache=True)
def phase(mu,ice):
    g=.78+.04*ice;back=-.25
    return .85*(1-g*g)/(4*math.pi*(1+g*g-2*g*mu)**1.5)+.15*(1-back*back)/(4*math.pi*(1+back*back-2*back*mu)**1.5)

@njit(cache=True)
def scatter(d,ice):
    g=.78+.04*ice if np.random.random()<.85 else -.25
    u=np.random.random();c=(1+g*g-((1-g*g)/(1-g+2*g*u))**2)/(2*g);c=max(-1.,min(1.,c))
    phi=2*math.pi*np.random.random();s=math.sqrt(max(0.,1-c*c))
    axis=np.array([0.,0.,1.]) if abs(d[2])<.9 else np.array([0.,1.,0.])
    a=np.cross(axis,d);a/=math.sqrt(np.dot(a,a));b=np.cross(d,a)
    v=d*c+s*(math.cos(phi)*a+math.sin(phi)*b)
    return v/math.sqrt(np.dot(v,v))

@njit(cache=True)
def ratio_transmission(p,d,lo,hi,dims,ext,ice,majorant):
    entry,end=box(p,d,lo,hi)
    if end<=entry or majorant==0:return 1.,0
    t,tr=entry,1.
    for _ in range(100000):
        t+=-math.log(max(1e-16,np.random.random()))/majorant
        if t>=end:return tr,0
        sigma,_ice=voxel(p+d*t,lo,hi,dims,ext,ice)
        if sigma>majorant*(1+1e-10):return 0.,1
        tr*=max(0.,1-sigma/majorant)
        if tr<.001:
            if np.random.random()>=.1:return 0.,0
            tr/=.1
    return 0.,1

@njit(cache=True)
def path_sample(origin,direction,sun,E,background,albedo,lo,hi,dims,ext,ice,majorant,multiple):
    p,d=origin.copy(),direction.copy();weight,radiance=1.,0.
    if majorant==0:return background,0
    for bounce in range(4096):
        entry,end=box(p,d,lo,hi)
        if end<=entry:return radiance+weight*background,0
        t=entry;found=False
        for _ in range(100000):
            t+=-math.log(max(1e-16,np.random.random()))/majorant
            if t>=end:return radiance+weight*background,0
            sigma,ice_fraction=voxel(p+d*t,lo,hi,dims,ext,ice)
            if sigma>majorant*(1+1e-10):return radiance,1
            if np.random.random()*majorant<sigma:
                found=True;break
        if not found:return radiance,1
        p=p+d*t;weight*=albedo
        if weight==0:return radiance,0
        if E>0:
            trans,fail=ratio_transmission(p,sun,lo,hi,dims,ext,ice,majorant)
            if fail:return radiance,1
            radiance+=weight*E*phase(np.dot(d,sun),ice_fraction)*trans
        if not multiple:return radiance,0
        d=scatter(d,ice_fraction)
        if bounce>=7:
            survival=min(.95,max(.05,weight))
            if np.random.random()>=survival:return radiance,0
            weight/=survival
    return radiance,1

@njit(cache=True)
def estimate(origin,direction,sun,E,background,albedo,lo,hi,dims,ext,ice,majorant,multiple,n,seed):
    np.random.seed(seed);total,total2,fail=0.,0.,0
    for _ in range(n):
        v,f=path_sample(origin,direction,sun,E,background,albedo,lo,hi,dims,ext,ice,majorant,multiple)
        total+=v;total2+=v*v;fail+=f
    mean=total/n;se=math.sqrt(max(0.,(total2-total*total/n)/(n-1))/n)
    return mean,se,fail

@njit(cache=True)
def deterministic_transmission(origin,d,lo,hi,dims,ext,ice,steps):
    a,b=box(origin,d,lo,hi)
    if b<=a:return 1.
    tau=0.;ds=(b-a)/steps
    for j in range(steps):
        sigma,_=voxel(origin+d*(a+(j+.5)*ds),lo,hi,dims,ext,ice);tau+=sigma*ds
    return math.exp(-tau)

def load_field(path):
    raw=path.read_bytes();s=json.loads(raw)
    if s['schema']!='open-moon-frozen-cloud/1':raise ValueError('Unsupported frozen field')
    dims=np.array(s['dimensions'],dtype=np.int64);lo=np.array(s['bounds_min_m'],dtype=float);hi=np.array(s['bounds_max_m'],dtype=float)
    ext=np.frombuffer(base64.b64decode(s['extinction_float32le']),dtype='<f4').astype(np.float64)
    ice=np.frombuffer(base64.b64decode(s['ice_fraction_float32le']),dtype='<f4').astype(np.float64)
    if dims.shape!=(3,) or np.any(dims<2) or lo.shape!=(3,) or hi.shape!=(3,) or not np.all(hi>lo) or len(ext)!=np.prod(dims) or len(ice)!=len(ext):raise ValueError('Grid shape mismatch')
    if not (np.isfinite(ext).all() and np.isfinite(ice).all() and np.all(ext>=0) and np.all((ice>=0)&(ice<=1))):raise ValueError('Invalid optical field')
    albedo=float(s.get('single_scattering_albedo',1.))
    if not 0<=albedo<=1:raise ValueError('Invalid single-scattering albedo')
    return s,lo,hi,dims,ext,ice,float(ext.max())*(1+1e-6),albedo,hashlib.sha256(raw).hexdigest()

@njit(cache=True)
def sampled_phase_check(n=131072):
    np.random.seed(314159)
    d=np.array([0.,0.,1.]);total,total2,norm_error=0.,0.,0.
    for _ in range(n):
        v=scatter(d,.5);c=np.dot(d,v);total+=c;total2+=c*c
        norm_error=max(norm_error,abs(np.dot(v,v)-1.))
    mean=total/n;se=math.sqrt((total2-total*total/n)/(n-1)/n)
    return mean,se,norm_error

def analytic_tests():
    lo=np.array([0.,0.,0.]);hi=np.array([10.,1.,10.]);dims=np.array([2,2,2]);ice=np.zeros(8);origin=np.array([5.,-1.,5.]);d=np.array([0.,1.,0.]);sun=d.copy();results=[]
    for name,tau,E,bg,albedo,multiple,expected in [
        ('vacuum',0.,0.,1.,1.,True,1.),
        ('pure absorption Beer',.7,0.,1.,0.,True,math.exp(-.7)),
        ('homogeneous single scattering',.4,1.,0.,1.,False,phase(1.,0.)*.4*math.exp(-.4)),
        ('conservative isotropic boundary',2.,0.,1.,1.,True,1.)]:
        ext=np.full(8,tau);mean,se,fail=estimate(origin,d,sun,E,bg,albedo,lo,hi,dims,ext,ice,tau*(1+1e-6),multiple,65536,271828)
        passed=abs(mean-expected)<=6*se+1e-8 and fail==0
        results.append({'name':name,'mean':mean,'standard_error':se,'expected':expected,'failures':fail,'passed':bool(passed)})
    # Numerical phase normalization and moment check for the mixture itself.
    mu,w=np.polynomial.legendre.leggauss(256)
    values=np.array([phase(x,.5) for x in mu]);integral=float(2*math.pi*np.sum(w*values));moment=float(2*math.pi*np.sum(w*mu*values))
    results.append({'name':'phase normalization and first moment','integral':integral,'moment':moment,'expectedMoment':.85*.8+.15*(-.25),'passed':abs(integral-1)<1e-10 and abs(moment-(.85*.8-.15*.25))<1e-10})
    mean,se,norm_error=sampled_phase_check()
    results.append({'name':'sampled phase moment and direction normalization','mean':mean,'standard_error':se,'expected':.85*.8-.15*.25,'max_norm_error':norm_error,'passed':abs(mean-(.85*.8-.15*.25))<=6*se and norm_error<1e-12})
    return results

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--field',type=Path,default=ROOT/'data/a1/generated/moon-fair.cloud.json')
    ap.add_argument('--gpu',type=Path,default=ROOT/'data/a1/generated/moon-fair.gpu.json')
    ap.add_argument('--samples',type=int,default=16384)
    ap.add_argument('--seed',type=int,default=73471)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    if a.samples<100:ap.error('At least 100 independent paths per ray required')
    s,lo,hi,dims,ext,ice,majorant,albedo,sha=load_field(a.field)
    gpu=json.loads(a.gpu.read_text())
    if not gpu['passed'] or gpu['field_sha256']!=sha:raise ValueError('GPU and reference field identities differ')
    start=time.monotonic();tests=analytic_tests()
    if not all(t['passed'] for t in tests):raise AssertionError(tests)
    sun=np.array(gpu['transport']['sun_direction'],dtype=float);sun/=np.linalg.norm(sun)
    result={'schema':'open-moon-a1-cloud-reference/1','field_sha256':sha,'field_file':a.field.name,'profile_sha256':s['provenance']['shared_atmospheric_profile_sha256'],
        'solver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'samples_per_ray_per_order':a.samples,'seed':a.seed,'majorant_m1':majorant,
        'algorithm':'analog delta-tracked free flight; ratio-tracked solar visibility; two-HG phase; next-event estimation; Russian roulette from bounce 8',
        'model':gpu['transport']['model'].replace('first order','first order versus multiple scattering'),
        'limitations':['Finite voxelized crop; sampling errors against continuous production density are separate.','Gray conservative cloud with selected phase function.','Unit solar boundary, zero molecular atmosphere, black external boundaries.','Production empirical radiance closure and whole-image cache remain unbenchmarked.'],
        'analytic_tests':tests,'rays':[],'unresolved_paths':0,'passed':True}
    for i,ray in enumerate(gpu['transport']['rays']):
        origin=np.array(ray['origin'],dtype=float);d=np.array(ray['direction'],dtype=float);d/=np.linalg.norm(d)
        low=estimate(origin,d,sun,1.,0.,albedo,lo,hi,dims,ext,ice,majorant,False,a.samples,a.seed+i*2)
        high=estimate(origin,d,sun,1.,0.,albedo,lo,hi,dims,ext,ice,majorant,True,a.samples,a.seed+i*2+1)
        reference=gpu['transport']['readbacks'][-1]['values'][i]
        agreement=abs(low[0]-reference['single_scatter_radiance'])<=6*low[1]+1e-5
        tr=deterministic_transmission(origin,d,lo,hi,dims,ext,ice,16384)
        entry={'ray':ray,'gpu_single_scatter':reference['single_scatter_radiance'],'gpu_transmittance':reference['transmittance'],'cpu_transmittance_16384':tr,
            'single':{'mean':low[0],'standard_error':low[1]},'multiple':{'mean':high[0],'standard_error':high[1]},'single_order_agrees_with_GPU_within_6SE_plus_1e_minus5':bool(agreement),
            'multiple_minus_single':high[0]-low[0]}
        result['rays'].append(entry);result['unresolved_paths']+=int(low[2]+high[2]);result['passed'] &= bool(agreement and abs(tr-reference['transmittance'])<.001)
        print(json.dumps({'ray':i,'single':low[:2],'multiple':high[:2],'agreement':bool(agreement),'elapsed_s':time.monotonic()-start}),flush=True)
    result['passed'] &= result['unresolved_paths']==0
    result['elapsed_s']=time.monotonic()-start
    output=a.output or a.field.with_name(a.field.stem.replace('.cloud','')+'.mc.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
