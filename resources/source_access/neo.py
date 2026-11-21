"""Catalogue geometry, reported diameter errors and mission-evidence coverage.

A geometric screen is neither a composition assay nor a qualified reserve.
"""
from __future__ import annotations
import hashlib,json,math
from pathlib import Path
from . import protection
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def sha(raw): return hashlib.sha256(raw).hexdigest()

def load_archive(config):
    folder=ROOT/config['archive_directory']; result={}
    for item in config['external_files']:
        raw=(folder/item['path']).read_bytes()
        if len(raw)!=item['bytes'] or sha(raw)!=item['sha256']:
            raise ValueError('External source identity mismatch: '+item['path'])
        blob=json.loads(raw)
        if item.get('signature') and blob.get('signature')!=item['signature']:
            raise ValueError('API signature mismatch')
        result[item['path']]=blob
    return result

def rows(catalogue):
    result=[dict(zip(catalogue['fields'],r)) for r in catalogue['data']]
    for r in result:r['spkid']=str(r['spkid'])
    if int(catalogue['count'])!=len(result) or len({r['spkid'] for r in result})!=len(result):
        raise ValueError('Incomplete or duplicate catalogue')
    return result

def number(row,key):
    value=row.get(key)
    if value in (None,''):return None
    value=float(value)
    if not math.isfinite(value):raise ValueError('Nonfinite '+key)
    return value

def mass(row,density):
    d=number(row,'diameter')
    if d is None or d<=0 or density<=0:raise ValueError('Positive diameter/density required')
    return math.pi/6*(1000*d)**3*density

def pools(records,config):
    dark=[r for r in records if number(r,'diameter') is not None and number(r,'diameter')>0 and number(r,'albedo') is not None and 0<=number(r,'albedo')<=config['maximum_albedo']]
    result={'dark_catalogue':dark}
    for name,g in config['geometry'].items():
        result[name]=[r for r in dark if all(number(r,k) is not None for k in ('a','e','i')) and g['a_min']<=number(r,'a')<=g['a_max'] and 0<=number(r,'e')<=g['e_max'] and 0<=number(r,'i')<=g['i_max']]
    return protection.filter_source_sets(result)

def name_key(s):return ' '.join(s.split())

def summarize(selected,config,nhats):
    density=config['density_reference_kg_m3']; ordered=sorted(selected,key=lambda r:mass(r,density),reverse=True)
    total=math.fsum(mass(r,density) for r in ordered)
    have=[r for r in ordered if number(r,'diameter_sigma') is not None]
    for r in have:
        if number(r,'diameter_sigma')<0:raise ValueError('Negative diameter uncertainty')
    # First-order propagation of reported independent diameter errors, fixed density.
    sigma=math.sqrt(math.fsum((3*mass(r,density)*number(r,'diameter_sigma')/number(r,'diameter'))**2 for r in have))
    nhats_names={name_key(r['fullname']):r for r in nhats['data']}
    if len(nhats_names)!=len(nhats['data']) or int(nhats['count'])!=len(nhats['data']):raise ValueError('Ambiguous/incomplete NHATS data')
    joined=[(r,nhats_names[name_key(r['full_name'])]) for r in ordered if name_key(r['full_name']) in nhats_names]
    envelopes={}
    for delta in config['common_diameter_scale_scenarios']:
        envelopes[str(delta)]={'mass_kg':[total*(1-delta)**3*config['density_scenario_kg_m3'][0]/density,total*(1+delta)**3*config['density_scenario_kg_m3'][1]/density], 'probability':None}
    return {'count':len(ordered),'mass_reference_kg':total,'known_taxonomy_count':sum(bool(r.get('spec_B') or r.get('spec_T')) for r in ordered),
        'largest_fraction':mass(ordered[0],density)/total if total else None,'top_five_fraction':sum(mass(r,density) for r in ordered[:5])/total if total else None,
        'reported_diameter_error':{'count':len(have),'covered_reference_mass_kg':sum(mass(r,density) for r in have),'missing_sigma_count':len(ordered)-len(have),'missing_sigma_reference_mass_kg':sum(mass(r,density) for r in ordered if number(r,'diameter_sigma') is None),'linearized_independent_sigma_kg':sigma,'interpretation':'Reported-error propagation for covered objects only, fixed density and membership; no systematic-error or reserve confidence claim.'},
        'structural_envelopes':envelopes,
        'nhats':{'matched_count':len(joined),'matched_mass_kg':sum(mass(r,density) for r,h in joined),'join':'Unique normalized full names, checked for duplicates; candidate SPK IDs retained. No approximate-name matching.', 'dv_bins':[{'maximum_human_round_trip_dv_km_s':dv,'count':sum(float(h['min_dv']['dv'])<=dv for r,h in joined),'mass_kg':sum(mass(r,density) for r,h in joined if float(h['min_dv']['dv'])<=dv)} for dv in (6,8,12)],'interpretation':'NHATS Earth human round trips in 2020–2045, <=450 days and >=8 days stay. Absence is not impossibility; presence does not qualify a cislunar cargo route.'},
        'objects':[{**r,'mass_reference_kg':mass(r,density),'nhats':nhats_names.get(name_key(r['full_name']))} for r in ordered], 'qualified_reserve_kg':None}

def product(config,archive=None):
    data=load_archive(config) if archive is None else archive
    records=rows(data['nea_catalogue.json']); screened=pools(records,config)
    details=[]
    for filename,d in data.items():
        if filename.startswith('sbdb_'):
            details.append({'spkid':d['object']['spkid'],'name':d['object']['fullname'],'condition_code':d['orbit']['condition_code'],'orbit_id':d['orbit']['orbit_id'],'physical_parameters':d.get('phys_par',[]),'source_file':filename})
    return {'schema':'terluna.resources.source-catalogue/1','catalogue_count':len(records),'diameter_count':sum(number(r,'diameter') is not None for r in records),'dark_diameter_count_before_holds':len(screened['source_sets']['dark_catalogue'])+len(screened['excluded']['dark_catalogue']), 'groups':{k:summarize(v,config,data['nhats_summary.json']) for k,v in screened['source_sets'].items()},'exclusions':screened['excluded'],'detail_audit':details,'qualified_reserve_kg':None,'evidence':'Observed catalogue screen; formal errors and assumed common-scale/density envelopes separated. No unobserved-object extrapolation or dark-to-water conversion.'}
