"""Join pinned LCDB rotational evidence to the shared, protected source pools.

Raw archive bytes stay in the ignored research input cache. The derived product
preserves solution quality and qualifiers, and grants no extraction approval.
"""
from __future__ import annotations
import argparse
from collections import Counter
from copy import deepcopy
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import urllib.request
import zipfile
from . import neo, protection

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def load_inputs(path=None):
    config=json.loads(Path(path or HERE/'lightcurve_inputs.json').read_text())
    if config.get('schema')!='terluna.resources.source-lightcurve-inputs/1':
        raise ValueError('Unexpected lightcurve input schema')
    return config


def checked_bytes(raw,identity,label):
    if len(raw)!=identity['bytes'] or sha(raw)!=identity['sha256']:
        raise ValueError('Lightcurve source identity mismatch: '+label)
    return raw


def read_members(config,archive_directory=None):
    """Verify the ZIP and selected members without extracting arbitrary paths."""
    folder=Path(archive_directory) if archive_directory is not None else ROOT/config['archive_directory']
    archive=config['archive']
    raw=checked_bytes((folder/archive['filename']).read_bytes(),archive,archive['filename'])
    members={}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for item in config['members']:
            if z.namelist().count(item['name'])!=1:
                raise ValueError('Missing or duplicate LCDB archive member: '+item['name'])
            members[item['name']]=checked_bytes(z.read(item['name']),item,item['name'])
    return members


def fetch_archive(config,archive_directory=None):
    """Explicit retrieval of the public snapshot; changed remote bytes fail closed."""
    folder=Path(archive_directory) if archive_directory is not None else ROOT/config['archive_directory']
    item=config['archive'];folder.mkdir(parents=True,exist_ok=True)
    target=folder/item['filename']
    if target.exists():
        checked_bytes(target.read_bytes(),item,item['filename'])
        return target
    with urllib.request.urlopen(item['url'],timeout=60) as response:
        raw=response.read(item['bytes']+1)
    checked_bytes(raw,item,item['filename'])
    temporary=target.with_suffix(target.suffix+'.part')
    temporary.write_bytes(raw);temporary.replace(target)
    return target


def _numeric(text):
    if not text.strip():return None
    value=float(text)
    if not math.isfinite(value) or value<0:
        raise ValueError('Invalid LCDB physical value')
    return value


def parse_summary(raw,quality_codes):
    """Read the pinned README section 4.1.3 fixed-width column map.

    LCDB's number zero is an unnumbered object, not a shared catalogue identity.
    Size, albedo and taxonomic fields may be assumed and are not imported.
    """
    by_number={};by_designation={};records=[]
    lines=raw.decode('ascii').splitlines()
    if not lines or not lines[0].startswith('ASTEROID LIGHTCURVE DATABASE (LCDB) SUMMARY'):
        raise ValueError('Unexpected LCDB summary header')
    for line in lines[5:]:
        if not line.strip():continue
        if len(line)<190:raise ValueError('Truncated LCDB summary row')
        number=line[:7].strip()
        if number and not number.isdigit():raise ValueError('Invalid LCDB object number')
        number=int(number) if number else 0
        name=line[10:40].strip();designation=line[41:61].strip()
        record={'number':number or None,'name':name,'designation':designation,
                'period_flag':line[143:144].strip(),'period_hours':_numeric(line[145:158]),
                'period_description':line[159:174].strip(),'amplitude_flag':line[175:176].strip(),
                'amplitude_min_mag':_numeric(line[177:181]),'amplitude_max_mag':_numeric(line[182:186]),
                'period_quality':line[187:189].strip(),'notes':line[190:195].strip(),
                'binary_flag':line[196:199].strip(),'pole_reported':line[200:203].strip(),
                'survey_flag':line[204:209].strip(),'extended_note':line[210:213].strip(),
                'private_unpublished':line[214:217].strip(),
                'period_quality_meets_LCDB_statistical_threshold':line[187:189].strip() in quality_codes}
        if not name:raise ValueError('LCDB row lacks object name/designation')
        if number:
            if number in by_number:raise ValueError('Duplicate LCDB numbered object')
            by_number[number]=record
        for key in set((name,designation))-{''}:
            by_designation.setdefault(key,[]).append(record)
        records.append(record)
    return {'by_number':by_number,'by_designation':by_designation,'record_count':len(records)}


def match(row,index):
    found=re.match(r'^\s*(\d+)\s',row['full_name'])
    record=index['by_number'].get(int(found.group(1))) if found else None
    if record:return deepcopy(record),'MPC_number'
    designation=re.search(r'\(([^()]*)\)\s*$',row['full_name'])
    if designation:
        possible=index['by_designation'].get(designation.group(1),[])
        if len(possible)>1:raise ValueError('Ambiguous exact LCDB designation join')
        if possible:return deepcopy(possible[0]),'exact_designation'
    return None,None


def join(selected,index,density):
    """Join already screened shared pools; retain their original memberships."""
    membership={};objects={}
    for name,rows in selected.items():
        for row in rows:
            key=str(row['spkid'])
            if key in objects and objects[key]!=row:raise ValueError('Conflicting source rows')
            objects[key]=row;membership.setdefault(key,[]).append(name)
    result=[]
    for key,row in sorted(objects.items()):
        evidence,method=match(row,index)
        result.append({'source_id':key,'name':row['full_name'].strip(),
                       'groups':membership[key],'diameter_km_from_JPL':neo.number(row,'diameter'),
                       'scenario_mass_kg':neo.mass(row,density),'lcdb_match_method':method,'lcdb':evidence,
                       'extraction_authorized':False})
    return result


def summarize(candidates,group):
    rows=[r for r in candidates if group in r['groups']]
    matched=[r for r in rows if r['lcdb']]
    rated=[r for r in matched if r['lcdb']['period_quality_meets_LCDB_statistical_threshold']]
    return {'candidates':len(rows),'lcdb_matches':len(matched),
            'numeric_period_entries_including_limits':sum(r['lcdb']['period_hours'] is not None for r in matched),
            'quality_2minus_or_better':len(rated),
            'quality_3minus_or_3':sum(r['lcdb']['period_quality'] in ('3-','3') for r in matched),
            'quality_counts':dict(sorted(Counter(r['lcdb']['period_quality'] or 'unrated' for r in matched).items())),
            'binary_flags':dict(sorted(Counter(r['lcdb']['binary_flag'] or 'blank' for r in matched).items())),
            'tumbling_note_source_ids':[r['source_id'] for r in matched if 'T' in r['lcdb']['notes']],
            'catalogue_scenario_mass_kg':math.fsum(r['scenario_mass_kg'] for r in rows),
            'mass_with_LCDB_kg':math.fsum(r['scenario_mass_kg'] for r in matched)}


def product(config=None,archive_directory=None):
    config=load_inputs() if config is None else config
    input_path=HERE/config['neo_inputs']
    source_config=json.loads(input_path.read_text())
    # Only the retained JPL catalogue is needed; NHATS and detailed API calls
    # belong to the separate source-access producer.
    jpl_items=[i for i in source_config['external_files'] if i['path']=='nea_catalogue.json']
    if len(jpl_items)!=1:raise ValueError('Exactly one retained JPL catalogue identity required')
    catalogue=neo.load_archive({**source_config,'external_files':jpl_items})['nea_catalogue.json']
    members=read_members(config,archive_directory)
    index=parse_summary(members['lc_summary_pub.txt'],config['quality_codes_for_statistical_studies'])
    pc=protection.load_inputs()
    screened=neo.pools(neo.rows(catalogue),source_config)
    candidates=join(screened['source_sets'],index,source_config['density_reference_kg_m3'])
    dependencies={str(path.relative_to(ROOT)):sha(path.read_bytes()) for path in
                  (HERE/'neo.py',HERE/'protection.py',HERE/'protection_inputs.json',input_path)}
    return {'schema':'terluna.resources.source-lightcurves/1',
            'producer':{'path':'resources/source_access/lightcurves.py','sha256':sha(Path(__file__).read_bytes())},
            'inputs':{'path':'resources/source_access/lightcurve_inputs.json',
                      'canonical_sha256':sha((json.dumps(config,sort_keys=True,separators=(',',':'))+'\n').encode())},
            'dependencies':dependencies,'protection_sources':pc['sources'],
            'sources':{'lcdb_archive':config['archive'],'lcdb_members':config['members'],
                       'jpl_catalogue':jpl_items[0],'auxiliary_metadata_only':config['auxiliary_sources']},
            'source_release_date':config['source_release_date'],'retrieved':config['retrieved'],
            'evidence':'Factual join of curated pro-am lightcurve summaries to the shared protected candidate pools. Not an assay, qualified spin/anchoring model, reserve, route qualification or extraction authorization.',
            'reading_rule':config['reading_rule'],'rights':config['rights'],
            'units':{'period':'hours','amplitude':'magnitudes','diameter':'km','mass':'kg'},
            'density_assumption_kg_m3':source_config['density_reference_kg_m3'],
            'lcdb_summary_records':index['record_count'],
            'groups':{name:summarize(candidates,name) for name in screened['source_sets']},
            'exclusions':screened['excluded'],'candidates':candidates,'extraction_authorized':False}


def write_product(path,result):
    raw=(json.dumps(result,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(gzip.compress(raw,mtime=0))
    return sha(raw)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch',action='store_true',help='Fetch the pinned LCDB archive into the ignored cache')
    parser.add_argument('--archive-directory',type=Path)
    parser.add_argument('--output',type=Path,default=HERE/'results/lightcurves.json.gz')
    args=parser.parse_args();config=load_inputs()
    if args.fetch:fetch_archive(config,args.archive_directory)
    result=product(config,args.archive_directory)
    digest=write_product(args.output,result)
    print(json.dumps({'path':str(args.output),'decompressed_sha256':digest,'groups':result['groups']},sort_keys=True))


if __name__=='__main__':main()
