"""Analytic, identity, geometry and independent-reference gates for A3."""
from __future__ import annotations
import json,math,time
from pathlib import Path
import numpy as np
from numba import njit
import coupled_transport as t
import deterministic as dt

@njit(cache=True)
def phase_samples(n):
    np.random.seed(1707);d=np.array([0.,1.,0.]);s=s2=nerr=0.
    for _ in range(n):
        v=t.rayleigh_scatter(d);mu=v[1];s+=mu;s2+=mu*mu;nerr=max(nerr,abs(np.dot(v,v)-1.))
    return s/n,s2/n,nerr

def run(output=None):
    start=time.monotonic();c=t.load_inputs();checks=[]
    def add(name,passed,**data):checks.append(dict(name=name,passed=bool(passed),**data))
    # Bounds tested at every knot plus dense samples, including inversions.
    z=np.unique(np.r_[c['rows'][:,0],np.linspace(0,c['top'],10001)])
    air=np.array([t.molecular.sample_air(float(h),c['rows']) for h in z]);oz=np.array([t.molecular.ozone_du_m(float(h),c['oz']) for h in z]);idx=np.minimum(len(c['edges'])-2,np.searchsorted(c['edges'],z,side='right')-1)
    add('molecular and ozone majorants bound all profile knots and 10001 uniform probes',np.all(air<=c['air_bounds'][idx]) and np.all(oz<=c['oz_bounds'][idx]),max_air_fraction=float(np.max(air/c['air_bounds'][idx])))
    mu,w=np.polynomial.legendre.leggauss(128);P=np.array([t.rayleigh_phase(v) for v in mu]);integral=2*math.pi*np.sum(w*P);moment2=2*math.pi*np.sum(w*P*mu*mu)
    add('Rayleigh phase integrates to unity and has second moment 0.4',abs(integral-1)<1e-12 and abs(moment2-.4)<1e-12,integral=integral,second_moment=moment2)
    m,m2,err=phase_samples(262144)
    add('Rayleigh sampler has correct first/second moments and unit directions',abs(m)<.006 and abs(m2-.4)<.004 and err<1e-12,mean=m,second_moment=m2,max_norm_error=err)
    o=np.array([0.,2.,0.]);up=np.array([0.,1.,0.]);down=-up
    end,g=t.domain_end(o,up,c['R'],c['top']);end2,g2=t.domain_end(o,down,c['R'],c['top'])
    add('spherical top and solid-ground intersections',abs(end-(c['top']-2))<1e-6 and abs(end2-2)<1e-6 and not g and g2,top_distance=end,ground_distance=end2)
    # A homogeneous mixed-medium test uses the same spherical collision algorithm.
    a=c.copy();a.update(R=1e7,top=20.,rows=np.array([[0.,math.log(101325.),288.15],[20.,math.log(101325.),288.15]]),oz=np.array([1.,2.,3.,0.]),edges=np.array([0.,20.]),lo=np.array([-100.,0.,-100.]),hi=np.array([100.,20.,100.]),dims=np.array([2,2,2]),ice=np.zeros(8),beta=np.array([.01]),sigma=np.array([0.]),solar=np.array([1.]))
    a['air_bounds'],a['oz_bounds']=t.shell_majorants(a['rows'],a['oz'],a['edges'])
    for j,(name,gb,cl,albedo,mode,bg,expect) in enumerate([
        ('vacuum radiance',0.,0.,1.,1,0.,0.),
        ('vacuum constant boundary',0.,0.,1.,1,1.,1.),
        ('pure absorption Beer transmission',0.,.04,0.,1,1.,math.exp(-.04*18)),
        ('mixed homogeneous first-order',.01,.02,1.,0,0.,(.01*t.rayleigh_phase(1)+.02*t.cloud.phase(1,0))*18*math.exp(-.03*18)),
        ('conservative mixed medium with isotropic boundaries',.01,.03,1.,1,1.,1.),
    ]):
        a['beta']=np.array([gb]);a['ext']=np.full(8,cl);a['cmax']=cl*(1+1e-6);a['calbedo']=albedo;a['solar']=np.array([1. if 'first-order' in name else 0.])
        mean,se,fail=t.estimate(*t.args(a,o,up,up,0,gb>0,cl>0,mode),65536,9133+j,bg)
        add(name,abs(mean[0]-expect)<=6*se[0]+1e-10 and fail==0,mean=float(mean[0]),standard_error=float(se[0]),expected=expect,unresolved=int(fail))
    # Cubic along-ray cloud integral: exact 2-point Gauss per trilinear cell.
    lo=np.array([0.,0.,0.]);hi=np.array([1.,1.,1.]);dims=np.array([2,2,2]);ext=np.array([1+x+2*y+3*z+4*x*y*z for z in [0,1] for y in [0,1] for x in [0,1]],float);ice=np.zeros(8)
    d=np.ones(3)/math.sqrt(3);v=dt.cloud_depth(np.zeros(3),d,lo,hi,dims,ext,ice)
    # Along x=y=z=u: sigma=1+6u+4u^3; ds=sqrt(3)du.
    expected=5*math.sqrt(3)
    add('cell-exact frozen trilinear cloud line integral',abs(v-expected)<1e-12,computed=v,expected=expected)
    # Exact transmittance additivity, finite crop, reciprocal traversal.
    o=np.array([0.,2.,0.]);d=np.array([.2,1.,.1]);d/=np.linalg.norm(d)
    total=dt.cloud_depth(o,d,c['lo'],c['hi'],c['dims'],c['ext'],c['ice']);a,b=t.cloud.box(o,d,c['lo'],c['hi']);mid=(a+b)/2
    split=dt.cloud_depth(o,d,c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],mid)+dt.cloud_depth(o+d*mid,d,c['lo'],c['hi'],c['dims'],c['ext'],c['ice']);rev=dt.cloud_depth(o+d*(b+1),-d,c['lo'],c['hi'],c['dims'],c['ext'],c['ice'])
    add('cloud depth additivity and reversal',abs(total-split)<1e-10 and abs(total-rev)<1e-10,total=total,split=split,reverse=rev)
    # Clear reduction directly compares an independent all-order sampler with A2.
    sun=np.array([2**-.5,2**-.5,0.]);o=np.array([0.,2.,0.]);d=up
    A2=dt.a2_clear(c,o,d,sun)
    for wave in [400.,460.,520.,580.,640.,700.]:
        k=int(np.where(c['lam']==wave)[0][0]);m,se,fail=t.estimate(*t.args(c,o,d,sun,k,True,False,1),65536,813700+k)
        add(f'zero-cloud full MC reproduces A2 at {wave:g} nm',abs(m[0]-A2[k])<=6*se[0]+.01*A2[k] and fail==0,mean=float(m[0]),standard_error=float(se[0]),A2=float(A2[k]),relative_difference=float((m[0]-A2[k])/A2[k]),unresolved=int(fail))
    # Independently compute coupled single scattering using cell-exact columns.
    S4=dt.independent_single(c,o,d,sun,4,8);S8=dt.independent_single(c,o,d,sun,8,16)
    err=float(np.sum(abs(S4-S8))/np.sum(S8));add('independent single-scattering quadrature refinement',err<.005,spectral_L1=err)
    for wave in [460.,580.,700.]:
        k=int(np.where(c['lam']==wave)[0][0]);m,se,fail=t.estimate(*t.args(c,o,d,sun,k,True,True,0),65536,71818+k)
        add(f'coupled first-order MC agrees with deterministic transfer at {wave:g} nm',abs(m[0]-S8[k])<=6*se[0]+.002*S8[k] and fail==0,mean=float(m[0]),standard_error=float(se[0]),deterministic=float(S8[k]),unresolved=int(fail))
    # Check the interface to the inherited cloud-only solver with the exact field.
    k=int(np.where(c['lam']==580)[0][0]);m,se,fail=t.estimate(*t.args(c,o,d,sun,k,False,True,1),131072,18991)
    cm,cs,cf=t.cloud.estimate(o,d,sun,c['solar'][k],0.,c['calbedo'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['cmax'],True,131072,17871)
    add('zero-gas reduction agrees with inherited cloud-only MC',abs(m[0]-cm)<=6*math.sqrt(se[0]**2+cs**2) and fail+cf==0,mean=float(m[0]),standard_error=float(se[0]),inherited_mean=cm,inherited_standard_error=cs,unresolved=int(fail+cf))
    out=dict(schema='open-moon-a3-transport-tests/1',checks=checks,passed=all(v['passed'] for v in checks),elapsed_s=time.monotonic()-start)
    if output:Path(output).write_text(json.dumps(out,indent=2)+'\n')
    for x in checks:print(('PASS ' if x['passed'] else 'FAIL ')+json.dumps(x),flush=True)
    return out
if __name__=='__main__':
    r=run(t.ROOT/'validation/a3/transport-tests.json');raise SystemExit(0 if r['passed'] else 1)
