"""Bounded source-access assessment for the planned enclosed phase."""
from __future__ import annotations
import argparse,gzip,hashlib,json
from pathlib import Path
from . import neo,lunar,protection,costs
HERE=Path(__file__).resolve().parent

def read_checked(item):
    raw=(HERE/item['path']).read_bytes()
    if neo.sha(raw)!=item['sha256']:raise ValueError('Product hash mismatch: '+item['path'])
    result=json.loads(gzip.decompress(raw) if item['path'].endswith('.gz') else raw)
    if item.get('schema') and result.get('schema')!=item['schema']:raise ValueError('Unexpected product schema')
    return result

def demand_summary(demand,years,groups):
    # The frozen producer contains 81 primary portfolio/market cases.
    cases=[c for c in demand['cases'] if c['physical_scenario']=='base']
    if len(cases)!=81:
        raise ValueError('Expected 81 primary physical demand cases; inspect producer scenarios')
    out=[]
    for year in years:
        points=[next(a for a in c['annual'] if a['year']==year) for c in cases]
        feed=[p['cumulative_source_feed_kg']+c['seed']['source_feed_kg'] for p,c in zip(points,cases)]
        lunar_feed=[p['cumulative_lunar_raw_kg']+c['seed']['lunar_raw_kg'] for p,c in zip(points,cases)]
        fractions={}
        for name in ('narrow','wider'):
            group=groups[name];bounds=group['structural_envelopes']['0.2']['mass_kg']
            fractions[name]={'reference_mass_required_fraction':[min(feed)/group['mass_reference_kg'],max(feed)/group['mass_reference_kg']], 'joint_demand_mass_envelope_required_fraction':[min(feed)/bounds[1],max(feed)/bounds[0]],'rule':'Minimum fraction of hypothetical homogeneous whole-body feed required, before unknown source suitability, retained stocks and legal/site exclusions. Not a permitted extraction fraction.'}
        out.append({'year':year,'case_count':len(points),'cumulative_source_feed_kg_range':[min(feed),max(feed)],'cumulative_lunar_feed_kg_range':[min(lunar_feed),max(lunar_feed)],'required_feed_fractions':fractions})
    return out

def survey_products():
    products={}
    for name,filename in [('composition','composition_inputs.json'),('lightcurves','lightcurve_inputs.json')]:
        path=HERE/'results'/(name+'.json.gz');x=json.loads(gzip.decompress(path.read_bytes()))
        if x['producer']['sha256']!=neo.sha((HERE/(name+'.py')).read_bytes()):raise ValueError('Stale '+name+' producer')
        if name=='composition':
            for relative,expected in x['protection_dependencies'].items():
                if neo.sha((HERE/relative).read_bytes())!=expected:raise ValueError('Changed composition protection dependency')
            config=json.loads((HERE/filename).read_text())
            if x['inputs_sha256']!=neo.sha(json.dumps(config,sort_keys=True,separators=(',',':')).encode()):raise ValueError('Changed composition inputs')
        else:
            for relative,expected in x['dependencies'].items():
                if neo.sha((neo.ROOT/relative).read_bytes())!=expected:raise ValueError('Changed lightcurve dependency: '+relative)
            if x['inputs']['canonical_sha256']!=neo.sha((json.dumps(json.loads((HERE/filename).read_text()),sort_keys=True,separators=(',',':'))+'\n').encode()):raise ValueError('Changed lightcurve inputs')
        products[name]=x
    return products

def integrated_candidates(n,surveys):
    spectra={r['spkid']:r for r in surveys['composition']['wider_objects']}
    curves={r['source_id']:r for r in surveys['lightcurves']['candidates']}
    objects=n['groups']['wider']['objects']
    if set(spectra)!={r['spkid'] for r in objects}:raise ValueError('Composition pool differs from current source screen')
    result=[]
    for r in objects:
        spec=spectra[r['spkid']];lc=curves[r['spkid']]['lcdb']
        result.append({'source_id':r['spkid'],'name':r['full_name'].strip(),'spectral_family':spec['family'],'spectral_uncertainty_flag':spec['has_uncertain_class'],'has_period_record':lc is not None,'period_quality_2minus_or_better':bool(lc and lc['period_quality_meets_LCDB_statistical_threshold']),'has_nhats_trajectory_in_query':r['nhats'] is not None,'mass_reference_kg':r['mass_reference_kg'],'protection_status':protection.evaluate_release({'region':'near_earth','source_id':r['spkid'],'name':r['full_name']})['release_status']})
    return {'objects':result,'spectral_and_period_count':sum(r['spectral_family']!='no_observation' and r['has_period_record'] for r in result),'spectral_and_period_and_nhats_count':sum(r['spectral_family']!='no_observation' and r['has_period_record'] and r['has_nhats_trajectory_in_query'] for r in result),'qualification_rule':'Evidence overlap identifies follow-up candidates; no spectrum, period or human-trajectory record supplies an ore assay, cargo route or protection release.'}

def product():
    config=json.loads((HERE/'access_inputs.json').read_text())
    lpath=HERE/'results/lunar.json.gz';lp=json.loads(gzip.decompress(lpath.read_bytes()))
    if lp['producer']['sha256']!=neo.sha((HERE/'lunar.py').read_bytes()) or lp['source_register']['sha256']!=neo.sha((HERE/'lunar_sources.json').read_bytes()):raise ValueError('Stale lunar product')
    if lp['inputs_canonical_sha256']!=neo.sha((json.dumps(json.loads((HERE/'lunar_inputs.json').read_text()),sort_keys=True,separators=(',',':'))+'\n').encode()):raise ValueError('Changed lunar inputs')
    n=neo.product(config);coeff=read_checked(config['cost_coefficients']);d=read_checked(config['demand'])
    checks=demand_summary(d,config['checkpoints'],n['groups'])
    surveys=survey_products()
    integrated=integrated_candidates(n,surveys)
    pc=protection.load_inputs()
    release_examples={region:protection.evaluate_release({'region':region,'name':'Unspecified candidate; no attestations supplied'},config=pc) for region in ('lunar','near_earth')}
    return {'schema':'terluna.resources.source-access/1','producer':{'path':'resources/source_access/run.py','sha256':neo.sha(Path(__file__).read_bytes())},'model_files':{f:neo.sha((HERE/f).read_bytes()) for f in ('neo.py','costs.py','protection.py','lunar.py')},'input_files':{f:neo.sha((HERE/f).read_bytes()) for f in ('access_inputs.json','protection_inputs.json','lunar_inputs.json','lunar_sources.json')},'dependencies':{'lunar_product_sha256':neo.sha(lpath.read_bytes()),'physical_demand':config['demand'],'engineering_coefficients':config['cost_coefficients']},'evidence':'Survey-derived composition and geometry; conditional cost and structural uncertainty scenarios. No empirical economic validation, joint probabilistic reserve confidence, released extraction quantity or new financial return estimate.', 'boundary':'Off-Earth enclosed build through year 200 and sustained capacity to 220. Near-Earth asteroids are heliocentric sources supplying cislunar facilities; no Open Moon, array or outer Solar System expansion.', 'units':{'mass':'kg','currency':'constant USD2041 common valuation, not U exchange rate'},'composition_stratified_estimate':surveys['composition']['wider_stratified_completion'],'composition_pools':surveys['composition']['pools'],'lightcurve_groups':surveys['lightcurves']['groups'],'integrated_candidate_evidence':integrated,'survey_product_hashes':{name:neo.sha((HERE/'results'/(name+'.json.gz')).read_bytes()) for name in surveys},'neo':n,'lunar_summary':lp['global_summary'],'costs':costs.cost_product(coeff,config['costs'],lp['pixels']),'demand_checkpoints':checks,'protection_examples':release_examples,'qualified_supply_kg':None,'qualified_supply_note':'No site-specific release records supplied. Unknown qualification is not zero physical resource.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=HERE/'results/access.json.gz');a=parser.parse_args()
    result=product();raw=(json.dumps(result,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_bytes(gzip.compress(raw,mtime=0));print(a.output)

if __name__=='__main__':main()
