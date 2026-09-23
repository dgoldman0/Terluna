"""A3 matched spectral gas+cloud transport.

Full-domain backward null-collision Monte Carlo in the A2 spherical atmosphere
and the exact frozen A2 cloud. Both media participate at every path vertex and
on sunlight shadow rays. A second mode freezes only the diffuse molecular source
to A2's saved angular moments, isolating cloud feedback into molecular transport.

Coordinates are the production local (x, y, z), planet centre (0, -R, 0).
Radiance excludes the solar disk. The production boundary is black ground/space,
with collimated solar normal irradiance. No A2 sky is added to a full Monte Carlo
estimate: doing so would count the atmospheric photons twice.
"""
from __future__ import annotations
import hashlib, json, math, sys
from pathlib import Path
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[2]
for sub in ['a1','a2']:
    p=str(ROOT/'tools'/sub)
    if p not in sys.path: sys.path.insert(0,p)
import cloud_reference as cloud
import reference_optics as optics
import molecular_multiple_scattering as molecular
KB=1.380649e-23; NREF=101325/(KB*288.15); DU=2.687e20

@njit(cache=True)
def radial(p,R):
    r=math.sqrt(p[0]*p[0]+(p[1]+R)**2+p[2]*p[2])
    h=(p[0]*p[0]+p[2]*p[2]+p[1]*(2*R+p[1]))/(r+R)
    return r,h

@njit(cache=True)
def domain_end(p,d,R,top):
    q=np.array([p[0],p[1]+R,p[2]])
    r2=np.dot(q,q);pd=np.dot(q,d)
    end=max(0.,-pd+math.sqrt(max(0.,pd*pd+(R+top)**2-r2)))
    disc=pd*pd+R*R-r2
    if pd<0. and disc>0.:
        ground=-pd-math.sqrt(disc)
        if ground>=-1e-6 and ground<end: return max(0.,ground),True
    return end,False

@njit(cache=True)
def segments(p,d,R,top,edges,lo,hi,cloud_on):
    """Partition into spherical majorant shells and the finite cloud crop."""
    end,ground=domain_end(p,d,R,top)
    ts=np.empty(2*len(edges)+4);ts[0]=0.;ts[1]=end;n=2
    q=np.array([p[0],p[1]+R,p[2]]);pd=np.dot(q,d);r2=np.dot(q,q)
    for h in edges[1:-1]:
        disc=pd*pd+(R+h)**2-r2
        if disc>0:
            rt=math.sqrt(disc)
            for t in (-pd-rt,-pd+rt):
                if 1e-7<t<end-1e-7:ts[n]=t;n+=1
    if cloud_on:
        a,b=cloud.box(p,d,lo,hi)
        if b>a:
            if 1e-7<a<end-1e-7:ts[n]=a;n+=1
            if 1e-7<b<end-1e-7:ts[n]=b;n+=1
    return np.sort(ts[:n]),ground

@njit(cache=True)
def coefficients(p,R,rows,oz,beta,sigma,lo,hi,dims,ext,ice,gas_on,cloud_on):
    r,h=radial(p,R);gscat=gabs=c=ic=0.
    if gas_on:
        gscat=beta*molecular.sample_air(max(0.,h),rows)
        gabs=sigma*DU*molecular.ozone_du_m(max(0.,h),oz)
    if cloud_on: c,ic=cloud.voxel(p,lo,hi,dims,ext,ice)
    return gscat,gabs,c,ic

@njit(cache=True)
def segment_majorant(p,R,edges,air_bounds,oz_bounds,beta,sigma,lo,hi,cmax,gas_on,cloud_on):
    r,h=radial(p,R);i=min(len(edges)-2,max(0,np.searchsorted(edges,h,side='right')-1))
    m=(beta*air_bounds[i]+sigma*DU*oz_bounds[i]) if gas_on else 0.
    if cloud_on and lo[0]<=p[0]<=hi[0] and lo[1]<=p[1]<=hi[1] and lo[2]<=p[2]<=hi[2]:m+=cmax
    return m

@njit(cache=True)
def next_collision(p,d,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,calbedo,gas_on,cloud_on):
    ts,ground=segments(p,d,R,top,edges,lo,hi,cloud_on);trials=0
    for j in range(1,len(ts)):
        a,b=ts[j-1],ts[j]
        if b-a<1e-7:continue
        mid=p+d*((a+b)/2)
        mg=segment_majorant(mid,R,edges,air_bounds,oz_bounds,beta,sigma,lo,hi,0.,gas_on,False)
        mc=cmax if cloud_on and lo[0]<=mid[0]<=hi[0] and lo[1]<=mid[1]<=hi[1] and lo[2]<=mid[2]<=hi[2] else 0.
        M=mg+mc
        if M<=0:continue
        flight=a
        while trials<100000:
            trials+=1;flight+=-math.log(max(1e-300,np.random.random()))/M
            if flight>=b:break
            q=p+d*flight;u=np.random.random()*M
            # Superposed independent Poisson proposals. Evaluate only the
            # selected constituent: gas profile lookup is unnecessary at cloud
            # proposals, and vice versa. This changes cost, not the estimator.
            if u<mg:
                _,h=radial(q,R);gs=beta*molecular.sample_air(max(0.,h),rows)
                ga=sigma*DU*molecular.ozone_du_m(max(0.,h),oz)
                if gs+ga>mg*(1+1e-10):return q,-1,0.,ground,trials
                if u<gs:return q,1,0.,ground,trials
                if u<gs+ga:return q,3,0.,ground,trials
            else:
                u-=mg;c,ic=cloud.voxel(q,lo,hi,dims,ext,ice)
                if c>mc*(1+1e-10):return q,-1,ic,ground,trials
                if u<c*calbedo:return q,2,ic,ground,trials
                if u<c:return q,3,ic,ground,trials
        if trials>=100000:return p,-1,0.,ground,trials
    return p+d*ts[-1],0,0.,ground,trials

@njit(cache=True)
def cloud_transmission(p,d,limit,lo,hi,dims,ext,ice,cmax):
    a,b=cloud.box(p,d,lo,hi);b=min(b,limit)
    if b<=a or cmax<=0:return 1.,0
    tr=1.;flight=a
    for _ in range(100000):
        flight+=-math.log(max(1e-300,np.random.random()))/cmax
        if flight>=b:return tr,0
        sigma,_=cloud.voxel(p+d*flight,lo,hi,dims,ext,ice)
        if sigma>cmax*(1+1e-10):return 0.,1
        tr*=max(0.,1-sigma/cmax)
        if tr<1e-5:
            if np.random.random()>=.1:return 0.,0
            tr/=.1
    return 0.,1

@njit(cache=True)
def transmission(p,d,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,gas_on,cloud_on,ground_blocks=True):
    """Unbiased product of independent gas/cloud ratio-tracking estimates.

    Splitting constituent Poisson processes avoids evaluating the slowly varying
    gas at every dense-cloud null event. Both factors stop at the same spherical
    boundary; planetary shadowing is evaluated geometrically.
    """
    limit,ground=domain_end(p,d,R,top)
    if ground and ground_blocks:return 0.,0
    tr=1.;trials=0
    if gas_on:
        ts,_=segments(p,d,R,top,edges,lo,hi,False)
        for j in range(1,len(ts)):
            a,b=ts[j-1],ts[j]
            if b-a<1e-7:continue
            M=segment_majorant(p+d*((a+b)/2),R,edges,air_bounds,oz_bounds,beta,sigma,lo,hi,0.,True,False)
            if M<=0:continue
            flight=a
            while trials<100000:
                trials+=1;flight+=-math.log(max(1e-300,np.random.random()))/M
                if flight>=b:break
                _,h=radial(p+d*flight,R)
                gs=beta*molecular.sample_air(max(0.,h),rows);ga=sigma*DU*molecular.ozone_du_m(max(0.,h),oz)
                if gs+ga>M*(1+1e-10):return 0.,1
                tr*=max(0.,1-(gs+ga)/M)
                if tr<1e-5:
                    if np.random.random()>=.1:return 0.,0
                    tr/=.1
            if trials>=100000:return 0.,1
    if cloud_on:
        ct,fail=cloud_transmission(p,d,limit,lo,hi,dims,ext,ice,cmax)
        return tr*ct,fail
    return tr,0

@njit(cache=True)
def rotate(d,c,phi):
    axis=np.array([0.,0.,1.]) if abs(d[2])<.9 else np.array([0.,1.,0.])
    a=np.cross(axis,d);a/=math.sqrt(np.dot(a,a));b=np.cross(d,a)
    v=c*d+math.sqrt(max(0.,1-c*c))*(math.cos(phi)*a+math.sin(phi)*b)
    return v/math.sqrt(np.dot(v,v))

@njit(cache=True)
def rayleigh_phase(mu):return 3*(1+mu*mu)/(16*math.pi)

@njit(cache=True)
def rayleigh_scatter(d):
    # Inversion of F(mu)=1/2 + 3mu/8 + mu^3/8; exact scalar Rayleigh.
    c=2*math.sinh(math.asinh(4*np.random.random()-2)/3)
    return rotate(d,c,2*math.pi*np.random.random())

@njit(cache=True)
def molecular_source(p,d,sun,R,rg,ag,mom,band):
    """A2 saved diffuse moment field convolved with the Rayleigh phase."""
    q=np.array([p[0],p[1]+R,p[2]]);r=math.sqrt(np.dot(q,q));n=q/r
    ms=min(1.,max(-1.,np.dot(n,sun)));qr=np.dot(n,d);nu=np.dot(d,sun)
    qt=(nu-ms*qr)/math.sqrt(max(1e-16,1-ms*ms));qp2=max(0.,1-qr*qr-qt*qt)
    i,u=molecular.bracket(rg,r);j,v=molecular.bracket(ag,math.asin(ms))
    M=np.empty(5)
    for c in range(5):M[c]=(1-u)*((1-v)*mom[i,j,band,c]+v*mom[i,j+1,band,c])+u*((1-v)*mom[i+1,j,band,c]+v*mom[i+1,j+1,band,c])
    return max(0.,3/(16*math.pi)*(M[0]+qr*qr*M[1]+qt*qt*M[2]+qp2*M[3]+2*qr*qt*M[4]))

@njit(cache=True)
def path_sample(origin,direction,sun,E,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,calbedo,gas_on,cloud_on,mode,rg,ag,mom,band,boundary=0.,roulette=.95):
    """mode 0: first order; 1: fully coupled; 2: frozen A2 diffuse gas source.
    mode 3: same frozen-source closure, evaluated by explicit clear-gas
    continuation after the first gas event, without the A2 moment grid.

    Components 0..2 are gas-only, cloud-only and mixed scattering histories in
    modes 0/1. In mode 2 they instead label explicit direct-gas, direct-cloud and
    terminated A2 diffuse-source contributions; these labels are kept separate.
    """
    p=origin.copy();d=direction.copy();weight=1.;out=np.zeros(3);flags=0;active_cloud=cloud_on
    for bounce in range(4096):
        p,kind,ic,ground,trials=next_collision(p,d,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,calbedo,gas_on,active_cloud)
        if kind<0:return out,1
        if kind==0:
            out[0]+=weight*boundary
            return out,0
        if kind==3:return out,0
        flags|=kind
        if E>0.:
            tr,failed=transmission(p,sun,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,gas_on,active_cloud)
            if failed:return out,1
            ph=rayleigh_phase(np.dot(d,sun)) if kind==1 else cloud.phase(np.dot(d,sun),ic)
            idx=kind-1 if mode==2 else flags-1
            out[idx]+=weight*E*ph*tr
        if mode==0:return out,0
        if mode==2 and kind==1:
            out[2]+=weight*molecular_source(p,d,sun,R,rg,ag,mom,band)
            return out,0
        if mode==3 and kind==1:active_cloud=False
        d=rayleigh_scatter(d) if kind==1 else cloud.scatter(d,ic)
        if bounce>=11:
            if np.random.random()>=roulette:return out,0
            weight/=roulette
    return out,1

@njit(cache=True,nogil=True)
def estimate(origin,direction,sun,E,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,calbedo,gas_on,cloud_on,mode,rg,ag,mom,band,n,seed,boundary=0.,roulette=.95):
    np.random.seed(seed);s=np.zeros(4);s2=np.zeros(4);failed=0
    for _ in range(n):
        c,f=path_sample(origin,direction,sun,E,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,calbedo,gas_on,cloud_on,mode,rg,ag,mom,band,boundary,roulette)
        failed+=f;vals=np.array([c[0]+c[1]+c[2],c[0],c[1],c[2]])
        s+=vals;s2+=vals*vals
    return s/n,np.sqrt(np.maximum(0.,(s2-s*s/n)/(n-1)/n)),failed

@njit(cache=True,nogil=True)
def estimate_transmission(origin,d,sun_unused,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,gas_on,cloud_on,n,seed):
    np.random.seed(seed);s=s2=0.;failed=0
    for _ in range(n):
        tr,f=transmission(origin,d,R,top,rows,oz,beta,sigma,edges,air_bounds,oz_bounds,lo,hi,dims,ext,ice,cmax,gas_on,cloud_on)
        s+=tr;s2+=tr*tr;failed+=f
    return s/n,math.sqrt(max(0.,(s2-s*s/n)/(n-1)/n)),failed

def shell_majorants(rows,oz,edges):
    """Certified bounds under log-linear pressure, linear temperature/ozone.

    p_max / T_min bounds density on each shell, even through inversions.
    Include all interior profile nodes and interpolated interval endpoints.
    """
    air=[];o3=[]
    for a,b in zip(edges[:-1],edges[1:]):
        z=np.unique(np.r_[a,b,rows[(rows[:,0]>a)&(rows[:,0]<b),0]])
        lp=np.interp(z,rows[:,0],rows[:,1]);T=np.interp(z,rows[:,0],rows[:,2])
        air.append(math.exp(float(lp.max()))/(KB*float(T.min())*NREF)*(1+1e-9))
        q=np.unique(np.r_[a,b,oz[1] if a<oz[1]<b else a])
        o3.append(max(molecular.ozone_du_m(float(h),oz) for h in q)*(1+1e-9))
    return np.array(air),np.array(o3)

def load_inputs(field_path=None):
    profile=ROOT/'data/a1/generated/moon-fair.profile.json';p,rows,oz,shells,sha=optics.unpack_profile(profile)
    field_path=field_path or ROOT/'data/a2/generated/moon-fair-grid-fine.cloud.json'
    f,lo,hi,dims,ext,ice,cmax,albedo,fsha=cloud.load_field(field_path)
    molecular_path=ROOT/'data/a2/generated/moon-fair-molecular-standard-full.npz';m=np.load(molecular_path)
    if f['provenance']['shared_atmospheric_profile_sha256']!=sha or str(m['profile_sha256'])!=sha:raise ValueError('Profile identities differ')
    R=float(p['planet']['radius']);top=float(p['upper']['top_m'])
    edges=np.array([0.,10000.,25000.,40000.,80000.,160000.,320000.,top])
    airb,ozb=shell_majorants(rows,oz,edges)
    return dict(profile=p,profile_sha256=sha,field=f,field_sha256=fsha,field_file=field_path.name,molecular_field_sha256=hashlib.sha256(molecular_path.read_bytes()).hexdigest(),R=R,top=top,rows=rows,oz=oz,shells=shells,edges=edges,air_bounds=airb,oz_bounds=ozb,lo=lo,hi=hi,dims=dims,ext=ext,ice=ice,cmax=cmax,calbedo=albedo,rg=m['rg'],ag=m['ag'],mom=m['source_mom'],beam=m['beam'],srg=m['srg'],sag=m['sag'],lam=m['lam'],beta=m['beta'],sigma=m['sigma'],solar=m['solar'])

def args(c,origin,direction,sun,band,gas=True,clouds=True,mode=1):
    o=np.asarray(origin,float);d=np.asarray(direction,float);u=np.asarray(sun,float)
    if not all(np.isfinite(x).all() for x in (o,d,u)) or np.linalg.norm(d)==0 or np.linalg.norm(u)==0:raise ValueError('Invalid ray')
    r,h=radial(o,c['R'])
    if not 0<h<c['top']:raise ValueError('Origin must lie inside atmospheric domain')
    return (o,d/np.linalg.norm(d),u/np.linalg.norm(u),float(c['solar'][band]),c['R'],c['top'],c['rows'],c['oz'],float(c['beta'][band]),float(c['sigma'][band]),c['edges'],c['air_bounds'],c['oz_bounds'],c['lo'],c['hi'],c['dims'],c['ext'],c['ice'],c['cmax'],c['calbedo'],gas,clouds,mode,c['rg'],c['ag'],c['mom'],band)

def source_hashes():
    files=['tools/a3/coupled_transport.py','tools/a1/reference_optics.py','tools/a1/cloud_reference.py','tools/a2/molecular_multiple_scattering.py','data/a1/spectral-inputs.json']
    return {f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
