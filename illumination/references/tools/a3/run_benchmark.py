"""Execute reproducible matched 48-band gas+cloud experiments and controls.

The saved A2 moment field is a comparator, not an added background in the full
solver. All treatments use one atmosphere, cloud grid, Sun and camera ray.
Each wavelength receives independent random histories. Two independent seed
batches are retained; uncertainty is pooled before spectral integration.
"""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, json, math, platform, time
from pathlib import Path
import numpy as np
import numba
import coupled_transport as t
import deterministic as dt

def cases(c):
    E=float(c['field']['provenance']['a2']['baseExtent_m']);out=[]
    def add(name,x,z,sd,d=(0.,1.,0.)):
        R=c['R'];y=((R+2)**2-R**2-x*x-z*z)/(math.sqrt((R+2)**2-x*x-z*z)+R)
        dr=np.array(d,float);dr/=np.linalg.norm(dr)
        out.append(dict(id=name,sun_elevation_at_anchor_deg=sd,origin_m=[x,y,z],direction=dr.tolist(),sun_direction=[math.cos(math.radians(sd)),math.sin(math.radians(sd)),0.]))
    for name,i in [('dense',0),('moderate',1),('clear-gap',3),('edge',4),('broken',6)]:
        ix=i%4;iz=i//4;x=-E+(ix+.5)/4*(2*E);z=-E+(iz+.5)/3*(2*E);add(name+'-sun45',x,z,45.)
    add('oblique-sun45',0.,0.,45.,(.5,1.,-.3))
    add('dense-sun6',-E+(.5)/4*(2*E),-E+(.5)/3*(2*E),6.)
    add('zenith-sun6',0.,0.,6.)
    return out

def pooled(batch):
    n=np.array([r['n'] for r in batch]);mu=np.array([r['mean'] for r in batch]);se=np.array([r['standard_error'] for r in batch]);N=int(n.sum())
    M=np.sum(mu*n[:,None],axis=0)/N
    ss=np.sum((n*(n-1))[:,None]*se*se+n[:,None]*(mu-M)**2,axis=0)
    return dict(mean=M.tolist(),standard_error=np.sqrt(np.maximum(0.,ss/(N-1)/N)).tolist(),samples=N,unresolved=sum(r['unresolved'] for r in batch),batches=batch)

def color_weights():
    spec=json.loads((t.ROOT/'data/a1/spectral-inputs.json').read_text());cie=np.array(spec['cie1931_xyz']);dw=np.full(len(cie),10.);dw[[0,-1]]=5.;W=683*cie*dw[:,None]
    return W,W@np.array(spec['xyz_to_linear_srgb']).T

def summarize_spectrum(s,xyz,rgb):
    means=np.array(s['mean']);ses=np.array(s['standard_error']);L=means[:,0];err=ses[:,0]
    return dict(XYZ=(L@xyz).tolist(),linear_srgb=(L@rgb).tolist(),XYZ_standard_error=np.sqrt((err*err)@(xyz*xyz)).tolist(),components_Y=(means[:,1:].T@xyz[:,1]).tolist(),components_Y_standard_error=np.sqrt((ses[:,1:].T**2)@(xyz[:,1]**2)).tolist())

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--samples',type=int,default=65536);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--output',type=Path,default=t.ROOT/'data/a3/generated/benchmark.json');ap.add_argument('--case-ids',nargs='*');a=ap.parse_args()
    if a.samples<1024 or a.samples%2:ap.error('Use an even sample count >=1024')
    c=t.load_inputs();C=cases(c);C=[r for r in C if not a.case_ids or r['id'] in a.case_ids]
    modes=[('clear_mc',True,False,1),('coupled_mc',True,True,1),('a2_frozen_source',True,True,2),('clear_source_mc',True,True,3)]
    n=a.samples//2;cache=t.ROOT/'data/a3/generated/jobs';cache.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    identity=dict(profile=c['profile_sha256'],field=c['field_sha256'],molecular=c['molecular_field_sha256'],sources=t.source_hashes(),n=n)
    runhash=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest();jobs=[]
    for ci,r in enumerate(C):
        for mi,(name,g,cl,mode) in enumerate(modes):
            for k in range(48):
                for b in range(2):jobs.append((ci,mi,k,b))
    def work(job):
        ci,mi,k,b=job;r=C[ci];name,g,cl,mode=modes[mi];key=f'{runhash[:12]}-{r["id"]}-{name}-{k}-{b}';path=cache/(key+'.json')
        if path.exists():return job,json.loads(path.read_text())
        seed=772399+ci*100000+mi*10000+k*31+b*7
        mean,se,fail=t.estimate(*t.args(c,r['origin_m'],r['direction'],r['sun_direction'],k,g,cl,mode),n,seed)
        record=dict(n=n,seed=seed,mean=mean.tolist(),standard_error=se.tolist(),unresolved=int(fail))
        path.write_text(json.dumps(record)+'\n');return job,record
    # Compile before parallel dispatch; subsequent calls release the GIL.
    t.estimate(*t.args(c,C[0]['origin_m'],C[0]['direction'],C[0]['sun_direction'],10,True,True,1),2,101)
    records={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for i,(job,result) in enumerate(pool.map(work,jobs)):
            records[job]=result
            if (i+1)%96==0:print(json.dumps(dict(done=i+1,total=len(jobs),elapsed_s=time.monotonic()-start)),flush=True)
    xyz,rgb=color_weights();results=[]
    for ci,case in enumerate(C):
        r=dict(case);r['treatments']={}
        for mi,(name,g,cl,mode) in enumerate(modes):
            values=[pooled([records[(ci,mi,k,b)] for b in range(2)]) for k in range(48)]
            s=dict(mean=[v['mean'] for v in values],standard_error=[v['standard_error'] for v in values],batches_by_wavelength=[v['batches'] for v in values],unresolved=sum(v['unresolved'] for v in values),samples_per_wavelength=a.samples)
            s['color']=summarize_spectrum(s,xyz,rgb);r['treatments'][name]=s
        p=np.array(case['origin_m']);d=np.array(case['direction']);sun=np.array(case['sun_direction'])
        # Gray cloud-only transfer needs one estimate per ray, not 48 identical
        # spectral simulations. Its spectral uncertainty is perfectly correlated.
        ka=22;cloud_batches=[]
        for b in range(2):
            ar=list(t.args(c,p,d,sun,ka,False,True,1));ar[3]=1.
            mean,se,fail=t.estimate(*ar,a.samples*2,220003+ci*100+b)
            cloud_batches.append(dict(n=a.samples*2,seed=220003+ci*100+b,mean=mean.tolist(),standard_error=se.tolist(),unresolved=int(fail)))
        gray=pooled(cloud_batches);L=c['solar']*gray['mean'][0];SE=c['solar']*gray['standard_error'][0]
        r['cloud_only']=dict(gray=gray,mean=L.tolist(),standard_error=SE.tolist(),XYZ=(L@xyz).tolist(),linear_srgb=(L@rgb).tolist(),XYZ_standard_error=(SE@xyz).tolist(),spectral_uncertainty='perfectly correlated; gray transfer times the fixed solar spectrum')
        r['a2_clear']=dt.a2_clear(c,p,d,sun).tolist()
        # Independent first-order references and their quadrature convergence.
        first4=dt.independent_single(c,p,d,sun,4,8);first8=dt.independent_single(c,p,d,sun,8,16)
        r['coupled_single_scattering']=first8.tolist();r['single_quadrature_spectral_L1']=float(np.sum(abs(first8-first4))/max(1e-30,np.sum(first8)))
        replay6,diffuse,H=dt.production_replay(c,p,d,sun,np.array(r['a2_clear']),6);replay10,_,_=dt.production_replay(c,p,d,sun,np.array(r['a2_clear']),10)
        r['matched_production_expression']=replay10.tolist();r['replay_quadrature_spectral_L1']=float(np.sum(abs(replay6-replay10))/max(1e-30,np.sum(replay10)));r['replay_diffuse_falloff_scale_m']=H
        r['cloud_optical_depth']=float(dt.cloud_depth(p,d,c['lo'],c['hi'],c['dims'],c['ext'],c['ice']))
        n24,w24=np.polynomial.legendre.leggauss(24);air,oz,block=dt.gas_columns(p,d,c['R'],c['top'],c['rows'],c['oz'],c['shells'],n24,w24)
        r['direct_view_transmittance']=np.exp(-c['beta']*air-c['sigma']*t.DU*oz-r['cloud_optical_depth']).tolist()
        r['deterministic_color']={name:dict(XYZ=(np.array(r[name])@xyz).tolist(),linear_srgb=(np.array(r[name])@rgb).tolist()) for name in ['a2_clear','coupled_single_scattering','matched_production_expression']}
        results.append(r);print(json.dumps(dict(completed_case=r['id'],tau=r['cloud_optical_depth'],coupled_Y=r['treatments']['coupled_mc']['color']['XYZ'][1],elapsed_s=time.monotonic()-start)),flush=True)
    out=dict(schema='open-moon-a3-coupled-benchmark/1',identity=identity,run_sha256=runhash,wavelength_nm=c['lam'].tolist(),samples_per_wavelength_per_treatment=a.samples,seed_batches=2,software=dict(python=platform.python_version(),numpy=np.__version__,numba=numba.__version__),cloud_grid=dict(dimensions=c['dims'].tolist(),bounds_min_m=c['lo'].tolist(),bounds_max_m=c['hi'].tolist(),file=c['field_file']),model=dict(atmospheric_top_m=c['top'],radius_m=c['R'],surface='black absorbing sphere',external_boundary='black; collimated solar beam enters at top',cloud_phase='0.85 HG(g=.78+.04 ice) + 0.15 HG(g=-.25)',cloud_single_scattering_albedo=c['calbedo'],molecular_phase='scalar Rayleigh',ozone='shared A1/A2 profile; absorption only'),treatments=dict(clear_mc='Full spherical molecular MC; cloud absent.',coupled_mc='Full gas+cloud MC; every real and shadow path includes both media.',a2_frozen_source='Cloud transport with diffuse gas source frozen to the saved A2 moments; shadowed direct beam remains coupled.',clear_source_mc='The same frozen-source closure evaluated without A2 interpolation: after first gas scattering, continue through clear gas; retain coupled direct visibility at that event.',matched_production_expression='Spectral, matched-input replay of revision-04 diffuse/internal-light/composite equations; black ground. Not the production GPU/cache output.'),units='spectral radiance W m^-2 sr^-1 nm^-1; integrated XYZ photopic equivalents, Y cd m^-2',uncertainty='Independent wavelength and seed batches. Pooled ordinary Monte Carlo standard errors; no physical-model or discretization uncertainty included.',cases=results,unresolved_paths=sum(r['treatments'][name]['unresolved'] for r in results for name,_,_,_ in modes)+sum(r['cloud_only']['gray']['unresolved'] for r in results),elapsed_s=time.monotonic()-start)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print('Saved',a.output,'elapsed',out['elapsed_s'],flush=True)
if __name__=='__main__':main()
