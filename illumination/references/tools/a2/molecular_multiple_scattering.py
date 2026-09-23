"""Profile-consistent molecular multiple-scattering reference for A2.

Successive scattering orders use the same immutable A1 atmospheric profile as
its first-order reference. The solver is deliberately offline and additive: it
produces numerical evidence and does not alter the shoreline renderer.

The angular-source formulation follows the existing Terluna spectral transport
prototype's exact scalar-Rayleigh moment reduction, generalized here from an
exponential atmosphere to the A1 pressure/temperature/ozone profile.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, sys, time
from pathlib import Path
import numpy as np
from numba import njit, prange
ROOT=Path(__file__).resolve().parents[2]
A1_PATH=ROOT/'tools/a1/reference_optics.py'
_spec=importlib.util.spec_from_file_location('a1_reference_optics',A1_PATH);a1=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(a1)
KB=1.380649e-23;NREF=101325/(KB*288.15);DU=2.687e20;PI=math.pi

@njit(cache=True)
def bracket(grid,x):
    j=np.searchsorted(grid,x)-1;j=max(0,min(j,len(grid)-2));return j,max(0.,min(1.,(x-grid[j])/(grid[j+1]-grid[j])))
@njit(cache=True)
def sample_air(z,rows):
    lo,hi=0,rows.shape[0]-1
    while hi-lo>1:
        m=(lo+hi)//2
        if rows[m,0]<=z:lo=m
        else:hi=m
    t=min(1.,max(0.,(z-rows[lo,0])/(rows[hi,0]-rows[lo,0])));lp=rows[lo,1]+t*(rows[hi,1]-rows[lo,1]);T=rows[lo,2]+t*(rows[hi,2]-rows[lo,2])
    return math.exp(lp)/(KB*T*NREF)
@njit(cache=True)
def ozone_du_m(z,oz):
    b,p,t,amount=oz
    if z<=b or z>=t:return 0.
    f=(z-b)/(p-b) if z<p else (t-z)/(t-p)
    return 2*f/(t-b)*amount
@njit(cache=True)
def ray_path(r,mu,R,top_radius,shell_radii,rows,oz):
    disc=r*r*(mu*mu-1.)+R*R;ground=mu<0. and disc>0.
    if ground:end=-r*mu-math.sqrt(disc)
    else:end=-r*mu+math.sqrt(max(0.,r*r*(mu*mu-1.)+top_radius*top_radius))
    end=max(0.,end);ts=np.empty(2*len(shell_radii)+2);n=1;ts[0]=0.
    for rad in shell_radii:
        dd=r*r*(mu*mu-1.)+rad*rad
        if dd>0.:
            rt=math.sqrt(dd);t1=-r*mu-rt;t2=-r*mu+rt
            if t1>1e-7 and t1<end-1e-7:ts[n]=t1;n+=1
            if t2>1e-7 and t2<end-1e-7:ts[n]=t2;n+=1
    ts[n]=end;n+=1;ts=np.sort(ts[:n]);count=n-1;out=np.zeros((count,8));st=math.sqrt(max(0.,1-mu*mu))
    for k in range(count):
        tt=(ts[k]+ts[k+1])*.5;rad=math.sqrt(r*r+2*r*mu*tt+tt*tt);z=max(0.,rad-R)
        out[k,0]=rad;out[k,1]=(r*mu+tt)/rad;out[k,2]=(r+mu*tt)/rad;out[k,3]=tt*st/rad;out[k,4]=sample_air(z,rows);out[k,5]=ozone_du_m(z,oz);out[k,6]=ts[k+1]-ts[k];out[k,7]=z
    den=R if ground else top_radius;cend=(r+mu*end)/den;send=end*st/den
    return out,ground,cend,send
@njit(cache=True)
def weights_path(path,beta,sigma):
    n=len(path);nl=len(beta);w=np.zeros((n,nl));T=np.ones(nl)
    for k in range(n):
        for l in range(nl):
            bs=path[k,4]*beta[l];chi=bs+path[k,5]*DU*sigma[l];dt=chi*path[k,6];v=-math.expm1(-dt)
            w[k,l]=T[l]*v*(bs/chi if chi>0 else 0.);T[l]*=math.exp(-dt)
    return w,T
@njit(cache=True)
def solar_columns(r,mu,R,top_radius,rows,oz,n=256):
    if mu<0 and r*r*(1-mu*mu)<R*R-1e-4:return -1.,-1.
    d=-r*mu+math.sqrt(max(0.,top_radius*top_radius-r*r*(1-mu*mu)));ca=0.;co=0.
    for j in range(n):
        t0=d*(j/n)**2;t1=d*((j+1)/n)**2;t=(t0+t1)/2;rad=math.sqrt(r*r+2*r*mu*t+t*t);z=max(0.,rad-R)
        ca+=sample_air(z,rows)*(t1-t0);co+=ozone_du_m(z,oz)*(t1-t0)
    return ca,co
@njit(parallel=True,cache=True)
def solar_table(rgrid,agrid,R,top_radius,rows,oz,beta,sigma,solar,npath):
    nr,na,nl=len(rgrid),len(agrid),len(beta);out=np.zeros((nr,na,nl,2))
    for q in prange(nr*na):
        ir=q//na;ia=q%na;mu=math.sin(agrid[ia]);ca,co=solar_columns(rgrid[ir],mu,R,top_radius,rows,oz,npath)
        if ca<0:continue
        for l in range(nl):
            F=solar[l]*math.exp(-beta[l]*ca-sigma[l]*DU*co);out[ir,ia,l,0]=F;out[ir,ia,l,1]=F*max(0.,mu)
    return out
@njit(cache=True)
def interpolate(tab,rg,ag,r,alpha,l,c):
    i,u=bracket(rg,r);j,v=bracket(ag,alpha);return (1-u)*((1-v)*tab[i,j,l,c]+v*tab[i,j+1,l,c])+u*((1-v)*tab[i+1,j,l,c]+v*tab[i+1,j+1,l,c])
def angular_quadrature(nmu=12,nphi=12):
    if nmu%4:raise ValueError('nmu must be divisible by 4')
    g,w=np.polynomial.legendre.leggauss(nmu//4);mus=[];ws=[]
    for lo,hi in [(-1.,-.15),(-.15,0.),(0.,.15),(.15,1.)]:mus.extend((g+1)*(hi-lo)/2+lo);ws.extend(w*(hi-lo)/2)
    phi=(np.arange(nphi)+.5)*2*PI/nphi;return np.array(mus),np.array(ws),np.cos(phi),np.sin(phi),2*PI/nphi
def shell_altitudes(top,nshell,oz,sounding_end):
    h=top*np.linspace(0,1,nshell+1)**2
    extras=np.array([100.,500.,1000.,5000.,10000.,25000.,40000.,80000.,160000.,320000.,sounding_end,*oz[:3]])
    x=np.unique(np.r_[h,extras]);return x[(x>0)&(x<top)]
def build_paths(R,top,rgrid,mus,beta,sigma,rows,oz,nshell,sounding_end):
    shells=R+shell_altitudes(top,nshell,oz,sounding_end);items=[];K=0
    for r in rgrid:
        row=[]
        for mu in mus:
            p,g,c,s=ray_path(r,mu,R,R+top,shells,rows,oz);w,T=weights_path(p,beta,sigma);row.append((p,w,T,g,c,s));K=max(K,len(p))
        items.append(row)
    nr,nm,nl=len(rgrid),len(mus),len(beta);geom=np.zeros((nr,nm,K,4));weights=np.zeros((nr,nm,K,nl));ends=np.zeros((nr,nm,nl));nseg=np.zeros((nr,nm),np.int32);gb=np.zeros((nr,nm,3))
    for i in range(nr):
        for j in range(nm):
            p,w,T,g,c,s=items[i][j];n=len(p);geom[i,j,:n]=p[:,:4];weights[i,j,:n]=w;ends[i,j]=T;nseg[i,j]=n;gb[i,j]=[g,c,s]
    return geom,weights,ends,nseg,gb
@njit(parallel=True,cache=True)
def transport(prev,beam,rg,ag,srg,sag,mus,muw,cp,sp,dphi,geom,weights,ends,nseg,gb,include_beam):
    nr,na,nl=len(rg),len(ag),prev.shape[2];out=np.zeros_like(prev);fac=3./(16*math.pi)
    for idx in prange(nr*na):
        ir=idx//na;ia=idx%na;ca=math.cos(ag[ia]);sa=math.sin(ag[ia]);acc=np.zeros((nl,6))
        for im in range(len(mus)):
            mu=mus[im];st=math.sqrt(max(0.,1-mu*mu))
            for ip in range(len(cp)):
                vx=st*cp[ip];vy=st*sp[ip];nu=sa*mu+ca*vx;I=np.zeros(nl)
                # Black ground in A2: no reflected boundary term.
                for k in range(nseg[ir,im]):
                    r=geom[ir,im,k,0];qr=geom[ir,im,k,1];ss=sa*geom[ir,im,k,2]+ca*geom[ir,im,k,3]*cp[ip];ss=max(-1.,min(1.,ss));aa=math.asin(ss)
                    qt=(nu-ss*qr)/math.sqrt(max(1e-16,1-ss*ss));qp2=max(0.,1-qr*qr-qt*qt);ii,u=bracket(rg,r);jj,v=bracket(ag,aa)
                    if include_beam:bi,bu=bracket(srg,r);bj,bv=bracket(sag,aa)
                    for l in range(nl):
                        M=np.empty(5)
                        for c in range(5):M[c]=(1-u)*((1-v)*prev[ii,jj,l,c]+v*prev[ii,jj+1,l,c])+u*((1-v)*prev[ii+1,jj,l,c]+v*prev[ii+1,jj+1,l,c])
                        src=M[0]+qr*qr*M[1]+qt*qt*M[2]+qp2*M[3]+2*qr*qt*M[4]
                        if include_beam:src+=((1-bu)*((1-bv)*beam[bi,bj,l,0]+bv*beam[bi,bj+1,l,0])+bu*((1-bv)*beam[bi+1,bj,l,0]+bv*beam[bi+1,bj+1,l,0]))*(1+nu*nu)
                        I[l]+=weights[ir,im,k,l]*fac*max(0.,src)
                dw=muw[im]*dphi
                for l in range(nl):
                    x=I[l]*dw;acc[l,0]+=x;acc[l,1]+=x*mu*mu;acc[l,2]+=x*vx*vx;acc[l,3]+=x*vy*vy;acc[l,4]+=x*mu*vx;acc[l,5]+=x*max(mu,0.)
        out[ir,ia]=acc
    return out
def grids(R,top,quality):
    if quality=='draft':nr,nmu,nphi,nshell,npath=16,8,8,32,128;core=3.
    elif quality=='standard':nr,nmu,nphi,nshell,npath=22,12,12,48,192;core=2.
    else:nr,nmu,nphi,nshell,npath=28,16,16,64,256;core=1.5
    r=R+top*np.linspace(0,1,nr)**2
    deg=np.unique(np.r_[np.arange(-90,-30,6.),np.arange(-30,30.001,core),np.arange(36,91,6.)]);a=np.deg2rad(deg)
    sr=R+top*np.linspace(0,1,max(31,nr*2))**2;sdeg=np.unique(np.r_[np.arange(-90,-30,4.),np.arange(-30,30.001,core/2),np.arange(34,91,4.)]);sa=np.deg2rad(sdeg)
    return r,a,sr,sa,nmu,nphi,nshell,npath
@njit(cache=True)
def eval_ray(elevation,azimuth,sun_elevation,R,top,rows,oz,beta,sigma,beam,rg,ag,srg,sag,mom,shells):
    mu=math.sin(elevation);st=math.cos(elevation);cp=math.cos(azimuth);ss=math.sin(sun_elevation);cs=math.cos(sun_elevation);nu=ss*mu+cs*st*cp
    p,g,ce,se=ray_path(R+2.,mu,R,R+top,R+shells,rows,oz);w,T=weights_path(p,beta,sigma);out=np.zeros(len(beta));fac=3./(16*math.pi)
    if g:return out
    for k in range(len(p)):
        r=p[k,0];qr=p[k,1];ms=max(-1.,min(1.,ss*p[k,2]+cs*p[k,3]*cp));aa=math.asin(ms);qt=(nu-ms*qr)/math.sqrt(max(1e-16,1-ms*ms));qp2=max(0.,1-qr*qr-qt*qt)
        i,u=bracket(rg,r);j,v=bracket(ag,aa);bi,bu=bracket(srg,r);bj,bv=bracket(sag,aa)
        for l in range(len(beta)):
            F=(1-bu)*((1-bv)*beam[bi,bj,l,0]+bv*beam[bi,bj+1,l,0])+bu*((1-bv)*beam[bi+1,bj,l,0]+bv*beam[bi+1,bj+1,l,0])
            M=np.empty(5)
            for c in range(5):M[c]=(1-u)*((1-v)*mom[i,j,l,c]+v*mom[i,j+1,l,c])+u*((1-v)*mom[i+1,j,l,c]+v*mom[i+1,j+1,l,c])
            src=F*(1+nu*nu)+M[0]+qr*qr*M[1]+qt*qt*M[2]+qp2*M[3]+2*qr*qt*M[4];out[l]+=w[k,l]*fac*max(0.,src)
    return out
def solve_profile(profile_path,spec_path,quality,wavelengths,orders):
    state,rows,oz,a1shells,sha=a1.unpack_profile(profile_path);R=state['planet']['radius'];top=state['upper']['top_m'];sounding_end=state['sounding']['end_m'];spec=json.loads(spec_path.read_text());allw=np.array(spec['wavelength_nm'],float);idx=[]
    for w in wavelengths:
        q=np.where(np.isclose(allw,w))[0]
        if len(q)!=1:
            raise ValueError(f'wavelength {w} not in spectral input')
        idx.append(int(q[0]))
    lam=allw[idx];beta=1.24062e-6*(lam/1000.)**-4;sigma=np.array(spec['ozone_m2'],float)[idx];solar=np.array(spec['solar_W_m2_nm'],float)[idx]
    rg,ag,srg,sag,nmu,nphi,nshell,npath=grids(R,top,quality);t=time.monotonic();beam=solar_table(srg,sag,R,R+top,rows,oz,beta,sigma,solar,npath);mus,muw,cp,sp,dphi=angular_quadrature(nmu,nphi);pp=build_paths(R,top,rg,mus,beta,sigma,rows,oz,nshell,sounding_end)
    prev=np.zeros((len(rg),len(ag),len(lam),6));total=np.zeros_like(prev);history=[];first=None
    for order in range(1,orders+1):
        cur=transport(prev,beam,rg,ag,srg,sag,mus,muw,cp,sp,dphi,*pp,order==1);cur[np.abs(cur)<1e-30]=0.;total+=cur
        if order==1:first=cur.copy()
        roi=(ag>=math.radians(-20))&(ag<=math.radians(90));rel=float(np.max(cur[:,roi,:,0])/max(np.max(total[:,roi,:,0]),1e-30));history.append({'order':order,'max_radiance_moment_increment_fraction':rel,'elapsed_s':time.monotonic()-t});print(json.dumps(history[-1]),flush=True);prev=cur
        if order>=5 and rel<5e-4:break
    return state,rows,oz,sha,lam,beta,sigma,solar,rg,ag,srg,sag,beam,total,prev,history,dict(nmu=nmu,nphi=nphi,nshell=nshell,npath=npath,quality=quality,numerical_zero_floor=1e-30,elapsed_s=time.monotonic()-t)
def a1_single_reference(state,rows,oz,lam,beta,sigma,solar,elev,az,sun):
    # Independent A1 first-order implementation, with high quadrature settings.
    _,_,_,shells,_=a1.unpack_profile(Path(state['_profile_path'])) if False else (None,None,None,None,None)
    inputs=state['sounding']['inputs'];xs=sorted(set(x[0] for k in ('temperature','humidity') for x in inputs[k] if x[0]<=2.2));heights=np.interp(xs,np.log(state['planet']['ps'])-rows[:,1],rows[:,0]);shells=np.unique(np.concatenate((heights,[100,1000,5000,10000,25000,40000,80000,160000,320000,state['sounding']['end_m']],oz[:3])));shells=shells[(shells>0)&(shells<state['upper']['top_m'])];nodes,weights=np.polynomial.legendre.leggauss(24)
    return a1.sky_ray(state['planet']['radius'],state['upper']['top_m'],2.,math.radians(elev),math.radians(az),math.radians(sun),rows,oz,shells,nodes,weights,beta,sigma,solar,384)
def run_case(profile_path,spec_path,quality,wavelengths,orders,rays,field_output=None):
    state,rows,oz,sha,lam,beta,sigma,solar,rg,ag,srg,sag,beam,total,last,history,meta=solve_profile(profile_path,spec_path,quality,wavelengths,orders)
    source_mom=total-last
    if field_output is not None:
        field_output.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(field_output,rg=rg,ag=ag,srg=srg,sag=sag,beam=beam,source_mom=source_mom,lam=lam,beta=beta,sigma=sigma,solar=solar,profile_sha256=np.array(sha),quality=np.array(quality),orders=np.array(len(history),dtype=np.int32),history_json=np.array(json.dumps(history)),solver_meta_json=np.array(json.dumps(meta)))
    eval_shells=shell_altitudes(state['upper']['top_m'],128,oz,state['sounding']['end_m'])
    records=[];zero=np.zeros_like(source_mom)
    validation_targets=np.array([400.,460.,520.,580.,640.,700.])
    validation_idx=np.array([i for i,w in enumerate(lam) if np.any(np.isclose(w,validation_targets))],dtype=np.int64)
    if len(validation_idx)==0: validation_idx=np.arange(min(len(lam),6),dtype=np.int64)
    for ray in rays:
        e,a,sun_elev=ray
        single=eval_ray(math.radians(e),math.radians(a),math.radians(sun_elev),state['planet']['radius'],state['upper']['top_m'],rows,oz,beta,sigma,beam,rg,ag,srg,sag,zero,eval_shells)
        multi=eval_ray(math.radians(e),math.radians(a),math.radians(sun_elev),state['planet']['radius'],state['upper']['top_m'],rows,oz,beta,sigma,beam,rg,ag,srg,sag,source_mom,eval_shells)
        vi=validation_idx;ref=a1_single_reference(state,rows,oz,lam[vi],beta[vi],sigma[vi],solar[vi],e,a,sun_elev)
        error=float(np.sum(np.abs(single[vi]-ref))/max(1e-30,np.sum(ref)))
        records.append({'sun_deg':sun_elev,'view_elevation_deg':e,'view_azimuth_deg':a,'single':single.tolist(),'multiple_through_final_order':multi.tolist(),
            'a1_validation_wavelength_nm':lam[vi].tolist(),'a1_independent_single_validation':ref.tolist(),'single_relative_spectral_L1_vs_a1_validation_bands':error,
            'multiple_to_single_ratio':[float(multi[i]/single[i]) if single[i]>1e-30 else None for i in range(len(lam))]})
    return {'profile_sha256':sha,'wavelength_nm':lam.tolist(),'quality':quality,'solver_meta':meta,'order_history':history,'orders_computed':len(history),'rays':records}
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--profile',type=Path,default=ROOT/'data/a1/generated/moon-fair.profile.json');ap.add_argument('--spectral-inputs',type=Path,default=ROOT/'data/a1/spectral-inputs.json');ap.add_argument('--quality',choices=['draft','standard','fine'],default='draft');ap.add_argument('--wavelengths',nargs='+',type=float,default=[400,460,520,580,640,700]);ap.add_argument('--full-spectrum',action='store_true');ap.add_argument('--orders',type=int,default=8);ap.add_argument('--field-output',type=Path);ap.add_argument('--output',type=Path,default=ROOT/'data/a2/generated/moon-fair-molecular-multiple.json');a=ap.parse_args();wavelengths=json.loads(a.spectral_inputs.read_text())['wavelength_nm'] if a.full_spectrum else a.wavelengths;rays=[(90,0,45),(45,0,45),(15,0,45),(5,0,6),(5,90,6),(1,0,0),(15,180,45)];start=time.monotonic();case=run_case(a.profile,a.spectral_inputs,a.quality,wavelengths,a.orders,rays,a.field_output);case.update({'schema':'open-moon-a2-molecular-multiple/1','solver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'profile_file':a.profile.name,'spectral_inputs_sha256':hashlib.sha256(a.spectral_inputs.read_bytes()).hexdigest(),'model':'Spherical scalar Rayleigh transport by successive scattering orders; A1 molecular/ozone profile; collimated Sun; black surface; no aerosol, refraction, clouds or water-vapour absorption.','limitations':['Selected spectral bands only unless the caller requests the full input grid.','Angular/radial convergence must be checked by comparing quality levels.','Specified atmospheric profile is conditional, not a climate prediction.'],'elapsed_s':time.monotonic()-start});a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(case,indent=2)+'\n');mx=max(r['single_relative_spectral_L1_vs_a1_validation_bands'] for r in case['rays']);print(json.dumps({'output':str(a.output),'max_single_L1_vs_A1_validation_bands':mx,'orders':case['orders_computed'],'elapsed_s':case['elapsed_s']},indent=2))
if __name__=='__main__':main()
