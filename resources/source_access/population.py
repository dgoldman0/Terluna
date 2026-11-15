"""Conditional NEA material potential beyond the diameter-known catalogue.

This is a population/characterization model, not an extraction release, delivery
schedule or empirical confidence interval. Quantities are contained kg before
recovery. Known objects and the inferred residual population are disjoint.
"""
from __future__ import annotations
import argparse, gzip, hashlib, io, json, math
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from . import composition, neo, protection

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
MATERIALS=('body_mass_kg','water_equivalent_kg','native_iron_kg')
POOLS=('narrow','wider','all_no_don_quixote','all_nea')


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked(path, expected):
    raw=Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('Source identity mismatch: '+str(path))
    return raw


def size_spline(header):
    """Cumulative log10 N(D), independently reconstructed from supplied knots."""
    x=np.asarray(header['log10_diameter_knots_km']);g=np.asarray(header['slopes'])
    ref=math.log10(header['reference_diameter_km']);anchor=np.searchsorted(x,ref)-1
    y=np.zeros(len(x));y[anchor]=math.log10(header['reference_number'])+(x[anchor]-ref)*g[anchor]
    for j in range(anchor+1,len(x)):y[j]=y[j-1]+(x[j]-x[j-1])*g[j-1]
    for j in range(anchor-1,-1,-1):y[j]=y[j+1]-(x[j+1]-x[j])*g[j]
    return CubicSpline(x,y,bc_type='natural')


def orbit_table(raw,geometry):
    """Integrate the published piecewise-constant H/orbit grid over nested cuts.

    Sixteen-point Gaussian quadrature resolves q<1.3 AU intersections inside
    a/e cells. This integration is a grid approximation, not a trajectory fit.
    """
    lines=raw.decode().splitlines();first=lambda i:lines[i].split('!')[0].split()
    nh,hmin,hmax=map(float,first(1));na,amin,amax=map(float,first(6))
    ne,emin,emax=map(float,first(7));ni,imin,imax=map(float,first(8))
    values=np.loadtxt(io.StringIO('\n'.join(lines[13:])),usecols=range(5))
    ih,ia,ie,ii=(values[:,k].astype(int)-1 for k in range(4));weight=values[:,4]
    if np.any(weight<0):raise ValueError('Negative orbit density')
    da,de,di=(amax-amin)/na,(emax-emin)/ne,(imax-imin)/ni
    alo=amin+ia*da;elo=emin+ie*de;ilo=imin+ii*di
    nodes,qw=np.polynomial.legendre.leggauss(16)
    def fractional(g):
        left=np.maximum(alo,g.get('a_min',amin));right=np.minimum(alo+da,g.get('a_max',amax))
        width=np.maximum(0,right-left);area=np.zeros(len(alo))
        for node,w in zip(nodes,qw):
            a=left+(node+1)*width/2
            lower=np.maximum(elo,1-1.3/np.maximum(a,1e-12))
            upper=np.minimum(elo+de,g.get('e_max',emax))
            area+=w*np.maximum(0,upper-lower)
        inclination=np.maximum(0,np.minimum(ilo+di,g.get('i_max',imax))-ilo)/di
        return area*width/(2*da*de)*inclination
    total=np.bincount(ih,weights=weight*fractional({}),minlength=int(nh))
    result=[]
    cuts={name:np.bincount(ih,weights=weight*fractional(g),minlength=int(nh))/total for name,g in geometry.items()}
    for j in range(int(nh)):
        p=[float(cuts['narrow'][j]),float(cuts['wider'][j]-cuts['narrow'][j]),float(1-cuts['wider'][j])]
        if min(p)<-1e-10 or abs(sum(p)-1)>1e-10:raise ValueError('Non-nested orbital probabilities')
        result.append({'H_low':hmin+j*(hmax-hmin)/nh,'H_high':hmin+(j+1)*(hmax-hmin)/nh,'category_probabilities':p})
    header={'model_name':lines[0],'reference_diameter_km':float(first(2)[0]),'reference_number':float(first(2)[1]),'log10_diameter_knots_km':list(map(float,first(4))),'slopes':list(map(float,first(5)))}
    return header,result


def rayleigh(p,sigma):return p/sigma**2*np.exp(-p*p/(2*sigma**2))


def albedo_parameters(d):
    """Published simulator's size-dependent two-Rayleigh coefficients."""
    x=np.log10(d);knots=np.log10(np.sqrt([.1*.3,.3*1.,1.*3.]))
    return tuple(np.interp(x,knots,v) for v in ([.183,.212,.306],[.0566,.0367,.0251],[.182,.182,.162]))


class Characterization:
    def __init__(self,cfg,header):
        self.cfg=cfg;self.spline=size_spline(header);self.cache={}
        self.pgrid=np.geomspace(cfg['albedo_bounds'][0],cfg['albedo_bounds'][1],1024)

    def albedo(self,H,group,rng):
        """p(albedo|H) uses the size prior and the H/diameter Jacobian.

        H is held at its catalogued value. Quarter-mag lookup bins are an
        explicit numerical approximation, while final diameters use exact H.
        """
        H=np.asarray(H);group=np.asarray(group);result=np.empty(H.shape)
        bins=np.round(H*4).astype(int)
        for hb,g in set(zip(bins.ravel().tolist(),group.ravel().tolist())):
            key=(hb,g)
            if key not in self.cache:
                h=hb/4;p=self.pgrid;D=1329*10**(-h/5)/np.sqrt(p);ld=np.log10(D)
                fd,d,b=albedo_parameters(D)
                phi=rayleigh(p,.03) if g==1 else rayleigh(p,.17) if g==0 else fd*rayleigh(p,d)+(1-fd)*rayleigh(p,b)
                slope=self.spline(ld,1);logn=self.spline(ld)
                prior=np.maximum(0,-slope)*10**np.clip(logn,-30,30)
                weight=phi*prior*np.gradient(p)
                if not np.any(weight>0):raise ValueError('Empty conditional-albedo distribution')
                self.cache[key]=np.cumsum(weight)/sum(weight)
            mask=(bins==hb)&(group==g);result[mask]=np.interp(rng.random(np.count_nonzero(mask)),self.cache[key],self.pgrid)
        return result

    def materials(self,diameter_km,albedo,group,rng):
        """Return actual material fractions, never numerical wt.% as fractions."""
        p=np.asarray(albedo);g=np.array(group,copy=True);shape=p.shape
        prior=self.cfg['composition_prior'];pc=prior['carbonaceous']*rayleigh(p,.03)
        ps=prior['stony']*rayleigh(p,.17);px=prior['ambiguous_other']*rayleigh(p,.17)
        den=pc+ps+px;u=rng.random(shape)
        inferred=np.where(u<pc/den,1,np.where(u<(pc+ps)/den,0,2))
        g=np.where(g<0,inferred,g)
        # X/other labels are not identified iron bodies. Nonmetal outcomes retain
        # albedo-conditioned C/stony ambiguity instead of becoming pure metal.
        x=g==2;metal=rng.random(shape)<self.cfg['ambiguous_metal_probability']
        fallback=np.where(rng.random(shape)<pc/(pc+ps),1,0)
        g=np.where(x,np.where(metal,2,fallback),g)
        bounds=np.asarray(self.cfg['density_bounds_kg_m3']);density=bounds[g,0]+rng.random(shape)*(bounds[g,1]-bounds[g,0])
        mass=math.pi/6*(1000*np.asarray(diameter_km))**3*density
        uw=np.maximum(rng.random(shape),np.finfo(float).tiny);uf=np.maximum(rng.random(shape),np.finfo(float).tiny)
        water=np.where(g==0,-1.098*np.log(720*uw)+6.9235,np.where(g==1,-7.037*np.log(72*uw)+29.05,-1.469*np.log(15*uw)+3.6966))/100
        xiron=55*uf
        iron=np.where(g==0,-5.3999*np.log(615*uf)+37.8162,np.where(g==1,-3.265*np.log(45*uf)+11.344,-.001276*xiron**3+.09227*xiron**2-2.56*xiron+69.64))/100
        water=np.clip(water,0,1);iron=np.minimum(np.clip(iron,0,1),1-water)
        return np.stack([mass,mass*water,mass*iron],axis=-1)

    def synthetic(self,H,rng):
        group=np.full(np.shape(H),-1,dtype=int);p=self.albedo(H,group,rng)
        return self.materials(1329*10**(-np.asarray(H)/5)/np.sqrt(p),p,group,rng)


def category(row,geometry):
    def inside(g):return all(neo.number(row,k) is not None for k in ('a','e','i')) and g['a_min']<=float(row['a'])<=g['a_max'] and 0<=float(row['e'])<=g['e_max'] and 0<=float(row['i'])<=g['i_max']
    return 0 if inside(geometry['narrow']) else 1 if inside(geometry['wider']) else 2


def material_group(state):
    g=state['family']
    if g=='C_B' or g in ('D','T'):return 1
    if g=='S_Q' or g in ('A','K','L','O','R','V'):return 0
    return 2 if g=='X' else -1


def catalogued(records,matched,cfg,char,rng):
    draws=cfg['draws'];out={p:np.zeros((draws,3)) for p in POOLS};counts={p:0 for p in POOLS}
    known=[];excluded=[];unquantified=[];protection_cfg=protection.load_inputs()
    for r in records:
        if protection.exclusion_reasons(r,protection_cfg):excluded.append(r['spkid']);continue
        if neo.number(r,'H') is None and neo.number(r,'diameter') is None:unquantified.append(r['spkid']);continue
        known.append(r)
    # Object blocks keep the single-worker memory footprint bounded.
    for start in range(0,len(known),cfg['catalogue_block_size']):
        block=known[start:start+cfg['catalogue_block_size']];n=len(block)
        H=np.array([neo.number(r,'H') if neo.number(r,'H') is not None else 18 for r in block])
        group=np.array([material_group(composition.object_state(composition.entries_for(r,matched))) for r in block])
        h=np.broadcast_to(H,(draws,n));g=np.broadcast_to(group,(draws,n));p=char.albedo(h,g,rng)
        measured=np.array([neo.number(r,'diameter') is not None for r in block])
        D=1329*10**(-h/5)/np.sqrt(p)
        for j,r in enumerate(block):
            albedo=neo.number(r,'albedo')
            if albedo is not None and albedo>0:
                sigma=cfg['reported_albedo_log_sigma'];p[:,j]=np.clip(albedo*np.exp(sigma*rng.normal(size=draws)),*cfg['albedo_bounds'])
                if not measured[j]:D[:,j]=1329*10**(-H[j]/5)/np.sqrt(p[:,j])
            if measured[j]:
                d=neo.number(r,'diameter');sigma=neo.number(r,'diameter_sigma')
                if sigma is None:sigma=d*cfg['missing_diameter_sigma_fraction']
                if sigma<0 or d<=0:raise ValueError('Invalid reported diameter/sigma')
                D[:,j]=np.maximum(d+sigma*rng.normal(size=draws),d*.01)
                if (albedo is None or albedo<=0) and neo.number(r,'H') is not None:
                    p[:,j]=np.clip((1329*10**(-H[j]/5)/D[:,j])**2,*cfg['albedo_bounds'])
        values=char.materials(D,p,g,rng)
        for j,r in enumerate(block):
            cat=category(r,cfg['geometry']);selected=['all_nea']
            if not r['full_name'].lstrip().startswith('3552 '):selected.append('all_no_don_quixote')
            if cat<2:selected.append('wider')
            if cat==0:selected.append('narrow')
            for pool in selected:out[pool]+=values[:,j,:];counts[pool]+=1
    return out,{'quantified_objects_by_pool':counts,'protected_ids_excluded':excluded,'unquantified_no_H_or_diameter_ids':unquantified,'known_measured_diameter_count':sum(neo.number(r,'diameter') is not None for r in known),'known_inferred_diameter_count':sum(neo.number(r,'diameter') is None for r in known)}


def residual_counts(records,counts,orbits,geometry):
    """Subtract all known objects, including protected ones, before inference."""
    result=[]
    for b in counts:
        lo,hi=b['H_low'],b['H_high'];selected=[r for r in records if neo.number(r,'H') is not None and lo<=float(r['H'])<hi and 0<=float(r['a'])<4.2 and 0<=float(r['i'])<88]
        observed=np.bincount([category(r,geometry) for r in selected],minlength=3)
        covered=[o for o in orbits if o['H_low']>=lo-1e-9 and o['H_high']<=hi+1e-9]
        if not covered:raise ValueError('Missing orbital bins')
        prob=np.mean([o['category_probabilities'] for o in covered],axis=0)
        raw=np.array([b[k] for k in ('count_low','count_central','count_high')])[:,None]*prob[None,:]-observed
        central=np.maximum(0,raw[1]);low=np.minimum(central,np.maximum(0,raw[0]));high=np.maximum(central,np.maximum(0,raw[2]))
        result.append({**b,'known_count_by_category':observed.tolist(),'intrinsic_orbital_probabilities':prob.tolist(),'missing_mean_by_branch':[low.tolist(),central.tolist(),high.tolist()],'catalogue_excess_over_model_by_branch':np.maximum(0,-raw).tolist()})
    return result


def aggregate_population(mean_count,moments,char,lo,hi,cfg,rng):
    """Compound-Poisson sums; explicit bodies for sparse, moment sums for dense.

    Gaussian compound sums apply only above the declared count threshold. A
    projection onto nonnegative contained fractions preserves physical bounds.
    Their input moments are independently sampled and their MC error reported.
    """
    draws=cfg['draws'];n=rng.poisson(mean_count,size=draws);out=np.zeros((draws,3))
    if mean_count<cfg['explicit_count_threshold']:
        for j in range(int(n.max(initial=0))):
            mask=n>j
            H=rng.uniform(lo,hi,np.count_nonzero(mask));out[mask]+=char.synthetic(H,rng)
    else:
        mu,cov=moments
        noise=rng.multivariate_normal(np.zeros(3),cov,size=draws)
        out=n[:,None]*mu+np.sqrt(n)[:,None]*noise
        out=np.maximum(0,out);contained=out[:,1:].sum(axis=1)
        factor=np.minimum(1,np.divide(out[:,0],contained,out=np.ones(draws),where=contained>0));out[:,1:]*=factor[:,None]
    return out,n


def undiscovered(bins,cfg,char,rng):
    totals=[{p:np.zeros((cfg['draws'],3)) for p in POOLS} for _ in range(3)]
    numbers=[{p:np.zeros(cfg['draws'],dtype=int) for p in POOLS} for _ in range(3)];diagnostic=[]
    for b in bins:
        lo,hi=b['H_low'],b['H_high'];sample=char.synthetic(rng.uniform(lo,hi,cfg['moment_samples']),rng)
        mu=sample.mean(axis=0);cov=np.cov(sample,rowvar=False)
        diagnostic.append({'H_low':lo,'H_high':hi,'single_body_mean_kg':mu.tolist(),'moment_sample_count':len(sample),'mean_relative_MC_standard_error':(sample.std(axis=0)/np.sqrt(len(sample))/np.maximum(mu,1e-30)).tolist()})
        means=np.array(b['missing_mean_by_branch']);increments=np.vstack([means[0],np.diff(means,axis=0)])
        for cat in range(3):
            for branch in range(3):
                if increments[branch,cat]<=0:continue
                sums,n=aggregate_population(increments[branch,cat],(mu,cov),char,lo,hi,cfg,rng)
                pools=['all_nea','all_no_don_quixote']+(['wider'] if cat<2 else [])+(['narrow'] if cat==0 else [])
                for k in range(branch,3):
                    for p in pools:totals[k][p]+=sums;numbers[k][p]+=n
    return totals,numbers,diagnostic


def pack(values):
    return {m:values[:,j].tolist() for j,m in enumerate(MATERIALS)}


def quantiles(values):
    return {m:dict(zip(('p05','p50','p95'),np.quantile(values[:,j],[.05,.5,.95]).tolist())) for j,m in enumerate(MATERIALS)}


def product(cfg=None):
    cfg=json.loads((HERE/'population_inputs.json').read_text()) if cfg is None else cfg
    source=cfg['source_files'];archive=ROOT/cfg['archive_directory']
    raw=checked(ROOT/source['catalogue']['path'],source['catalogue']['sha256']);records=composition.jpl_rows(raw)
    survey=checked(ROOT/source['taxonomy']['path'],source['taxonomy']['sha256']);matched,_=composition.match_observations(records,composition.parse_table(survey))
    grid=checked(archive/source['neomod_grid']['filename'],source['neomod_grid']['sha256']);header,orbits=orbit_table(grid,cfg['geometry'])
    rng=np.random.default_rng(cfg['random_seed']);char=Characterization(cfg,header)
    known,audit=catalogued(records,matched,cfg,char,rng);bins=residual_counts(records,cfg['population_count_bins'],orbits,cfg['geometry'])
    missing,numbers,moments=undiscovered(bins,cfg,char,rng);scenarios={}
    for j,name in enumerate(('population_low','central','population_high')):
        scenarios[name]={'probability_weight':None,'pools':{p:{'catalogued':pack(known[p]),'undiscovered':pack(missing[j][p]),'total':pack(known[p]+missing[j][p]),'total_quantiles_kg':quantiles(known[p]+missing[j][p]),'undiscovered_count_draws':numbers[j][p].tolist()} for p in POOLS}}
    # Common calibration and grade sensitivities use identical realization IDs.
    # They are outer scenarios, not extra probability in the reference prior.
    for label,diameter_scale,grade_scale in cfg['sensitivity_branches']:
        pools={}
        for p in POOLS:
            parts={}
            for key,v in (('catalogued',known[p]),('undiscovered',missing[1][p]),('total',known[p]+missing[1][p])):
                z=v.copy()*diameter_scale**3;z[:,1:]*=grade_scale
                factor=np.minimum(1,np.divide(z[:,0],z[:,1:].sum(axis=1),out=np.ones(cfg['draws']),where=z[:,1:].sum(axis=1)>0));z[:,1:]*=factor[:,None]
                parts[key]=pack(z)
            parts['total_quantiles_kg']={m:dict(zip(('p05','p50','p95'),np.quantile(parts['total'][m],[.05,.5,.95]).tolist())) for m in MATERIALS}
            pools[p]=parts
        scenarios[label]={'probability_weight':None,'diameter_scale':diameter_scale,'grade_scale':grade_scale,'pools':pools}
    deps={str(p.relative_to(ROOT)):digest(p) for p in (Path(__file__),HERE/'population_inputs.json',HERE/'population_sources.json',HERE/'composition.py',HERE/'protection.py',HERE/'protection_inputs.json',HERE/'neo.py')}
    return {'schema':'terluna.resources.nea-population/1','producer':{'path':str(Path(__file__).relative_to(ROOT)),'sha256':digest(__file__)},'dependencies':deps,'source_files':source,'draws':cfg['draws'],'random_seed':cfg['random_seed'],'draw_id_rule':'Same ordinal draw identifies a coupled realization across nested pools and population/scenario branches; pools and branches must not be added.','evidence':cfg['evidence'],'scope':cfg['scope'],'priors':{k:v for k,v in cfg.items() if k not in ('source_files','population_count_bins')},'catalogue_audit':audit,'population_bins':bins,'orbital_model':{'header':header,'probabilities':orbits,'integration':'Piecewise-constant supplied cells; 16-point quadrature in a/e, exact inclination and cut intersections.'},'moment_diagnostics':moments,'scenarios':scenarios,'release_status':'All inferred supply remains unqualified; named catalogue science holds are excluded and unknown future targets require every existing physical/data gate.','qualified_reserve_kg':None}


def main():
    x=product();path=HERE/'results/population.json.gz';path.write_bytes(gzip.compress((json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0));print(path,flush=True)

if __name__=='__main__':main()
