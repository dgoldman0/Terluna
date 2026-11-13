"""Regional predictive resource distributions and explicit lunar depth scenarios."""
from __future__ import annotations
import gzip,hashlib,json,math
from pathlib import Path
import numpy as np
from . import lunar,costs
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quantiles(values):return np.quantile(values,[.05,.5,.95]).tolist()

def column_mass(depth,surface,deep,scale):
    """Analytic integral of rho(z)=deep-(deep-surface)*exp(-z/scale)."""
    if depth<=0 or np.any(np.asarray(scale)<=0) or np.any(np.asarray(deep)<surface):
        raise ValueError('Positive depth/scale and ordered densities required')
    return deep*depth-(deep-surface)*scale*(-np.expm1(-depth/scale))

def ranked_pool(pixels,selection):
    if selection is None:return pixels
    ordered=sorted(pixels,key=lambda p:p['abundances'][selection['oxide']],reverse=True)
    target=sum(p['solid_angle_sr'] for p in ordered)*selection['top_area_fraction']
    selected=[];remaining=target
    for pixel in ordered:
        weight=min(remaining,pixel['solid_angle_sr'])
        if weight>0:selected.append({**pixel,'solid_angle_sr':weight})
        remaining-=weight
        if remaining<=1e-15:break
    return selected

def regional_draws(pixels,n,rng):
    weights=np.array([p['solid_angle_sr'] for p in pixels],dtype=float);weights/=weights.sum()
    ids=rng.choice(len(pixels),n,p=weights)
    means=np.array([[p['abundances'][k] for k in lunar.OXIDES] for p in pixels])[ids]
    cov=np.array([p['covariance'] for p in pixels])[:,:6,:6]
    factors=np.linalg.cholesky(cov)[ids]
    draws=means+np.einsum('nij,nj->ni',factors,rng.normal(size=(n,6)))
    rejected=0
    for _ in range(100):
        bad=np.any(draws<0,axis=1);count=int(bad.sum())
        if not count:break
        rejected+=count
        draws[bad]=means[bad]+np.einsum('nij,nj->ni',factors[bad],rng.normal(size=(count,6)))
    else:raise ValueError('Regional nonnegative assay sampler did not converge')
    return draws,{'draws':n,'rejected_proposals':rejected,'selected_pixel_count':len(set(ids.tolist()))}

def product(config=None):
    c=config or json.loads((HERE/'lunar_potential_inputs.json').read_text())
    for key in ('lunar_product','cost_coefficients'):
        if sha(HERE/c[key]['path'])!=c[key]['sha256']:raise ValueError('Changed '+key)
    p=json.loads(gzip.decompress((HERE/c['lunar_product']['path']).read_bytes()))
    coeff=json.loads((HERE/c['cost_coefficients']['path']).read_text())
    constants=json.loads((ROOT/'shared/constants.json').read_text());radius=constants['moon']['radius_m']
    pixels=lunar.latitude_band(p['pixels'],*c['nonpolar_latitude_bounds'])
    n=c['draws'];rng=np.random.default_rng(c['random_seed'])
    density=rng.triangular(*c['deep_density_kg_m3_triangular'],n)
    scale=rng.triangular(*c['density_scale_height_m_triangular'],n)
    recovery=rng.triangular(*c['oxide_recovery_triangular'],n)
    process={key:rng.triangular(*value,n) for key,value in c['cost_priors'].items()}
    handling=[];bulk=[]
    for i in range(n):
        inputs={key:float(v[i]) for key,v in process.items()}
        inputs.update(water_grade=.0684,water_recovery=.7,ore_energy_multiplier=1.)
        cost=costs.unit_feed_cost(coeff,inputs,c['capital_rate'])
        handling.append(cost['lunar_handling_USD2041_kg_feed']);bulk.append(cost['lunar_bulk_USD2041_kg_product'])
    handling=np.array(handling);bulk=np.array(bulk)
    pools={};bulk_inventory=[]
    nonpolar_area=sum(p['solid_angle_sr'] for p in pixels)*radius**2
    for depth in c['depth_scenarios_m']:
        mass_column=column_mass(depth,c['surface_density_kg_m3'],density,scale)
        bulk_inventory.append({'depth_m':depth,'nonpolar_bulk_mass_kg_p05_p50_p95':quantiles(nonpolar_area*mass_column),
            'column_kg_m2_p05_p50_p95':quantiles(mass_column),
            'footprints':[{'area_km2':area,'bulk_mass_kg_p05_p50_p95':quantiles(area*1e6*mass_column)} for area in c['prospective_footprints_km2']]})
    for name,selection in c['selection'].items():
        selected=ranked_pool(pixels,selection);draws,diag=regional_draws(selected,n,rng);area=sum(p['solid_angle_sr'] for p in selected)*radius**2
        grades={}
        for key in ('FeO','Al2O3','TiO2'):
            grade=draws[:,lunar.ORDER.index(key)]
            grades[key]={'grade_p05_p50_p95':quantiles(grade),'recovered_oxide_equivalent_per_feed_p05_p50_p95':quantiles(grade*recovery),
                'feed_kg_per_kg_oxide_equivalent_p05_p50_p95':quantiles(1/(grade*recovery)),
                'handling_USD2041_per_kg_oxide_equivalent_p05_p50_p95':quantiles(handling/(grade*recovery)),
                'depth_transfer_scenarios':[{'grade_depth_factor':factor,'at_5m_per_1000km2_recovered_oxide_equivalent_kg_p05_p50_p95':quantiles(1000*1e6*column_mass(5,c['surface_density_kg_m3'],density,scale)*grade*factor*recovery)} for factor in c['grade_depth_factors']]}
        pools[name]={'area_km2':area/1e6,'fraction_of_lunar_area':area/(4*math.pi*radius**2),'region_count':len(selected),'minimum_region_area_km2':min(v['solid_angle_sr'] for v in selected)*radius**2/1e6,'sampling':diag,'composition':grades,'six_oxide_sum_p05_p50_p95':quantiles(draws.sum(axis=1))}
    return {'schema':'terluna.resources.lunar-potential/1','producer':{'path':'resources/source_access/lunar_potential.py','sha256':sha(__file__)},
        'dependencies':{name:sha(HERE/name) for name in ('lunar_potential_inputs.json','lunar_potential_sources.json','results/lunar.json.gz','inputs/cost_coefficients.json','lunar.py','costs.py','../../shared/constants.json')},
        'evidence':'Conditional Monte Carlo regional predictive intervals from mapped composition/covariance and disclosed density, process and recovery priors; separate fixed depth-transfer scenarios. Not a permitted reserve or calibrated global sampling interval.',
        'quantiles':[.05,.5,.95],'nonpolar_area_km2':nonpolar_area/1e6,'bulk_inventory':bulk_inventory,'regional_pools':pools,
        'partial_costs':{'lunar_handling_USD2041_kg_feed_p05_p50_p95':quantiles(handling),'lunar_bulk_USD2041_kg_product_p05_p50_p95':quantiles(bulk)},
        'rules':c['rules'],'extraction_authorized':False,'site_specific_exclusion_area_km2':None,
        'uncertainty_not_sampled':['Unknown systematic composition calibration','Continuity of surface grades with depth (separate factor scenarios)','Site-specific protected fractions and future-resource allocations','Finished-metal separability and quality','Lunar polar ice tonnage and yield']}

def main():
    x=product();p=HERE/'results/lunar_potential.json.gz';p.write_bytes(gzip.compress((json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0));print(p)
if __name__=='__main__':main()
