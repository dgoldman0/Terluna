"""Independent 48-bin versus 12-group scene-ray checks with sampling errors."""
import json,math,time,concurrent.futures
import numpy as np
import render_scene as r
import scene_transport as b
import geometry as g
from test_b1 import samples,path_args

def run(case,sun_deg,uv,N=4096):
    c,_,_,_,_,M,V,glo,ghi,step=r.prepare(case)
    origin=np.array([-2200.,5000.,-6500.]);h,n,_=g.surface(origin,np.array([0.,-1.,0.]),c['R'],V,glo,ghi,step);origin[1]-=h;origin+=n*2
    d=g.ray_direction(uv[0],uv[1],160/96,104,10,16);a=math.radians(sun_deg);az=math.radians(25);sun=np.array([math.sin(az)*math.cos(a),math.sin(a),math.cos(az)*math.cos(a)])
    s=json.loads((b.ROOT/'data/a1/spectral-inputs.json').read_text());dw=np.full(48,10.);dw[[0,-1]]=5.;fullW=683*np.array(s['cie1931_xyz'])*dw[:,None]*c['solar'][:,None]
    records={}
    for groups in [48,12]:
        blocks=np.array_split(np.arange(48),groups);ls=np.array([np.mean(c['lam'][ix]) for ix in blocks]);W=np.array([fullW[ix].sum(axis=0) for ix in blocks]);means=[];ses=[];fail=0
        for k,lam in enumerate(ls):
            ar=path_args(c,V,glo,ghi,step,origin,d,sun,lam=lam,radius=math.radians(.26667),gas=True,cloud=case!='clear',terrain=True,water=True)
            mean,se,ff=samples(ar,N,95101+k*73+groups*101);means.append(mean);ses.append(se);fail+=int(ff[1:].sum())
        means=np.array(means);ses=np.array(ses);X=means@W;E=np.sqrt((ses*ses)@(W*W))
        records[str(groups)]=dict(wavelength_nm=ls.tolist(),unit_transfer=means.tolist(),unit_transfer_se=ses.tolist(),XYZ=X.tolist(),XYZ_standard_error=E.tolist(),unresolved=fail)
    full=np.array(records['48']['XYZ']);coarse=np.array(records['12']['XYZ']);e=np.sqrt(np.array(records['48']['XYZ_standard_error'])**2+np.array(records['12']['XYZ_standard_error'])**2)
    delta=coarse-full
    return dict(scenario=case,sun_deg=sun_deg,pixel_fraction=uv,samples_per_band=N,treatments=records,Y_relative_difference=float(delta[1]/full[1]),Y_difference_in_standard_errors=float(delta[1]/e[1]),Y_relative_two_SE=float(2*e[1]/full[1]),provisional_gate=bool(abs(delta[1])<.03*full[1]+3*e[1]),gate_definition='absolute Y difference <=3% of full-spectrum Y plus 3 combined MC standard errors; selected-ray implementation diagnostic, not a global colour bound')

if __name__=='__main__':
    t0=time.monotonic();results=[]
    for args in [('fair',12,(.2,.86)),('fair',12,(.6,.2)),('high',6,(.5,.04))]:
        q=run(*args);results.append(q);print(json.dumps({k:q[k] for k in ['scenario','sun_deg','Y_relative_difference','Y_relative_two_SE','provisional_gate']}),flush=True)
    out=dict(schema='open-moon-b1-spectral-check/1',cases=results,total_paths=sum(60*q['samples_per_band'] for q in results),elapsed_s=time.monotonic()-t0,qualified_scope='Independent finite-disk, surface-enabled reference rays. Monte Carlo and spectral discretization are separated only to the reported sampling uncertainty. Image sampling, other spectra/angles/materials remain unqualified.')
    p=b.ROOT/'validation/b1/scene-spectral-check.json';p.write_text(json.dumps(out,indent=2)+'\n')
