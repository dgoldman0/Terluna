"""Whole-sphere energy accounting for a solved spectral atmosphere.

At each wavelength: incoming at atmosphere top = outgoing direct + outgoing
scattered + ground absorption + ozone absorption. Spatial quadrature and a
finite scattering-order truncation give a measurable residual.
"""
import json,sys
import numpy as np
from solver import Atmosphere,grids,angular_quadrature,build_paths,transport,interpolate,densities,ROOT

def audit(fn):
    d=np.load(fn);meta=json.loads(str(d['meta']));atm=Atmosphere(**meta['atmosphere'])
    # Use three actual solved wavelengths for this whole-sphere diagnostic.
    inds=np.array([int(np.argmin(abs(d['lam']-w))) for w in (440,560,680)])
    sr=d['sr'];sa=d['sa'];r=d['r'];a=d['a'];lam=d['lam'][inds];mu=np.sin(a);solar=d['solar'][inds];beta=d['beta'][inds];absorb=d['absorb'][inds]
    mm=np.zeros((len(r),len(a),len(inds),7));mm[...,:6]=(d['moments']-d['last_order'])[:,:,inds,:]
    beam=d['beam'][:,:,inds,:]
    _,_,_,_,nm,np_,ns=grids(atm,meta['quality']);mus,muw,cp,sp,dp=angular_quadrature(nm,np_)
    pp=build_paths(atm,r,mus,beta,absorb,ns)
    ans=transport(mm,beam,r,a,sr,sa,mus,muw,cp,sp,dp,*pp,atm.ground_albedo,True)
    solar_on_grid=np.zeros((len(r),len(a),len(inds),2))
    for i,rr in enumerate(r):
        for j,aa in enumerate(a):
            for l in range(len(inds)):
                for c in range(2):solar_on_grid[i,j,l,c]=interpolate(beam,sr,sa,rr,aa,l,c)
    outgoing_direct=2*np.trapezoid(solar_on_grid[-1,:,:,0]*np.maximum(-mu,0)[:,None],mu,axis=0)/solar
    outgoing_diffuse=2*np.trapezoid(ans[-1,:,:,6],mu,axis=0)/solar
    absorbed_ground=2*(atm.radius_m/atm.top)**2*(1-atm.ground_albedo)*np.trapezoid(ans[0,:,:,5]+solar_on_grid[0,:,:,1],mu,axis=0)/solar
    oz=np.array([densities(rr-atm.radius_m,atm.scale_height_m,atm.ozone_scale)[1] for rr in r])
    volume_integrand=np.trapezoid(d['moments'][:,:,inds,0]+solar_on_grid[:,:,:,0],mu,axis=1)*oz[:,None]*r[:,None]**2
    absorbed_ozone=2*np.trapezoid(volume_integrand,r,axis=0)*absorb/(atm.top**2*solar)
    residual=1-outgoing_direct-outgoing_diffuse-absorbed_ground-absorbed_ozone
    result=dict(world=atm.name,quality=meta['quality'],wavelength_nm=lam.tolist(),outgoing_direct_fraction=outgoing_direct.tolist(),outgoing_diffuse_fraction=outgoing_diffuse.tolist(),ground_absorption_fraction=absorbed_ground.tolist(),ozone_absorption_fraction=absorbed_ozone.tolist(),energy_residual_fraction=residual.tolist())
    print(json.dumps(result),flush=True);return result
if __name__=='__main__':
    out=[audit(x) for x in sys.argv[1:]];(ROOT/'work'/'energy_checks.json').write_text(json.dumps(out,indent=2))
