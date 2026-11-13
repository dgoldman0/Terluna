"""Join published asteroid spectral classifications without inferring ore grades.

The published MITHNEOS table and the frozen JPL catalogue are distinct source
products, not independent random samples. Each JPL object contributes once to a
pool's distribution. Raw classes, uncertainty flags and disagreement survive the
join. No spectrum is reclassified and no object receives an inferred label;
aggregate completion scenarios carry explicit assumptions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import urllib.request
from . import protection

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
# This grouping retains the near-infrared C/X degeneracy. It groups spectral
# labels only: X does not mean metal-rich; C/B does not assign a water assay.
FAMILIES = {
    **dict.fromkeys(('S','Sa','Sk','Sl','Sq','Sr','Sv','S_comp','S(IV)','Q'), 'S_Q'),
    **dict.fromkeys(('B','C','Cb','Cg','Cgh','Ch','F','G'), 'C_B'),
    **dict.fromkeys(('X','Xc','Xe','Xk','Xn','E','M','P'), 'X'),
    **{k:k for k in ('A','D','K','L','O','R','T','V')}, 'Ld':'L',
}
UNRESOLVED = {'no_observation', 'unclassified', 'conflicting_families'}


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical_designation(value: str) -> str:
    """Case/whitespace normalization only; no fuzzy or orbital matching."""
    return re.sub(r'\s+', '', value).upper()


def parse_table(raw: bytes) -> list[dict]:
    """Read identifier and taxonomy columns from the journal's original Table 7."""
    text = raw.decode('ascii')
    if 'ajac532ft7_mrt.txt' not in text or '124-129' not in text:
        raise ValueError('Unexpected journal table header/layout')
    observations = []
    for line in text.splitlines():
        if not re.match(r'a\S+\.txt\s', line):
            continue
        if len(line) < 129:
            raise ValueError('Truncated journal observation')
        filename = line[:32].strip()
        number = int(line[33:39]) if line[33:39].strip() else None
        designation = line[40:50].strip()
        if not number and not designation:
            raise ValueError('Observation without an object identifier')
        observations.append({'spectrum': filename, 'number': number,
                             'designation': designation, 'taxonomy': line[123:129].strip()})
    if not observations or len({o['spectrum'] for o in observations}) != len(observations):
        raise ValueError('Missing or duplicate spectrum records')
    return observations


def jpl_rows(raw: bytes) -> list[dict]:
    data = json.loads(raw)
    required = {'spkid','full_name','diameter','albedo','a','e','i','spec_B','spec_T'}
    if not required.issubset(data['fields']):
        raise ValueError('Missing JPL fields')
    rows = [dict(zip(data['fields'], values)) for values in data['data']]
    if len(rows) != int(data['count']) or len({str(r['spkid']) for r in rows}) != len(rows):
        raise ValueError('JPL count or unique SPK identity mismatch')
    for r in rows:
        r['spkid'] = str(r['spkid'])
        n = re.match(r'^\s*(\d+)\s', r['full_name'])
        d = re.search(r'\(([^)]+)\)', r['full_name'])
        r['number'] = int(n[1]) if n else None
        r['designation'] = d[1] if d else ''
    return rows


def match_observations(rows: list[dict], observations: list[dict]) -> tuple[dict, list[dict]]:
    """Require consistent exact number/designation matches; preserve unmatched."""
    by_number, by_designation = {}, {}
    for r in rows:
        for key, value, index in (('number', r['number'], by_number),
                                  ('designation', canonical_designation(r['designation']), by_designation)):
            if value:
                if value in index and index[value] != r['spkid']:
                    raise ValueError(f'Ambiguous JPL {key}')
                index[value] = r['spkid']
    matched, unmatched = defaultdict(list), []
    for o in observations:
        n = by_number.get(o['number'])
        d = by_designation.get(canonical_designation(o['designation']))
        if n and d and n != d:
            raise ValueError('Number and designation identify different JPL objects')
        identity = n or d
        if identity:
            matched[identity].append({**o, 'matched_by': 'both' if n and d else 'number' if n else 'designation'})
        else:
            unmatched.append(o)
    return dict(matched), unmatched


def taxonomy(raw: str | None, system: str = 'mithneos') -> dict:
    """Interpret flags and explicit spectral families, without material assays."""
    raw = (raw or '').strip()
    flag = 'low_snr' if '::' in raw else 'uncertain' if ':' in raw else 'unflagged'
    clean = raw.replace(':', '').strip()
    if not clean or clean == 'U' or flag == 'low_snr':
        return {'raw': raw, 'flag': flag, 'usable': False, 'family': None}
    if clean in FAMILIES:
        group = FAMILIES[clean]
    elif clean in ('C;X', 'C,X', 'CX', 'XC'):
        group = 'C_or_X'
    else:
        parts = re.split(r'[;,/]', clean)
        # Tholen compound labels are retained as ambiguous when they include U
        # or span families. Their letter ordering is not converted to a vote.
        if len(parts) == 1 and system == 'jpl_spec_T' and clean.isupper():
            parts = list(clean)
        groups = {FAMILIES.get(p) for p in parts}
        group = next(iter(groups)) if len(groups) == 1 and None not in groups else None
    return {'raw': raw, 'flag': flag, 'usable': group is not None, 'family': group}


def object_state(entries: list[dict]) -> dict:
    """Count an object once; repeated observations cannot amplify its weight."""
    labels = [taxonomy(e['taxonomy'], e['source']) for e in entries]
    groups = sorted({c['family'] for c in labels if c['usable']})
    group = ('no_observation' if not entries else 'unclassified' if not groups
             else groups[0] if len(groups) == 1 else 'conflicting_families')
    return {'family': group, 'candidate_families': groups,
            'has_uncertain_class': any(c['flag'] == 'uncertain' for c in labels),
            'has_excluded_low_snr': any(c['flag'] == 'low_snr' for c in labels),
            'source_classifications': entries}


def entries_for(row: dict, matched: dict) -> list[dict]:
    entries = [{'source': f'jpl_{field}', 'taxonomy': row[field]}
               for field in ('spec_B','spec_T') if row[field]]
    entries.extend({'source': 'mithneos', **o} for o in matched.get(row['spkid'], []))
    return entries


def assumed_mass(row: dict, density: float) -> float | None:
    """Uniform-density spherical scenario; missing diameters retain unknown mass."""
    if not math.isfinite(density) or density <= 0:
        raise ValueError('Positive finite assumed density required')
    try:
        diameter = float(row['diameter'])
    except (TypeError, ValueError):
        return None
    return math.pi/6*(1000*diameter)**3*density if math.isfinite(diameter) and diameter > 0 else None


def measure_distribution(states: list[dict], weights: list[float | None]) -> dict:
    """Finite-catalogue missing-label bounds and an explicit completion scenario.

    Bounds are identification ranges conditional on accepting observed labels,
    not sampling confidence intervals. Missing measure is outside the denominator.
    The completion scenario gives unknown labels the observed weighted mixture.
    It does not claim that observing selection is ignorable in the actual survey.
    """
    if len(states) != len(weights):
        raise ValueError('State/weight lengths differ')
    measures = defaultdict(list)
    for state, weight in zip(states, weights):
        if weight is not None:
            if not math.isfinite(weight) or weight < 0:
                raise ValueError('Invalid distribution weight')
            measures[state['family']].append(weight)
    totals = {k:math.fsum(v) for k,v in sorted(measures.items())}
    total = math.fsum(totals.values())
    unknown = math.fsum(v for k,v in totals.items() if k in UNRESOLVED)
    known = math.fsum(v for k,v in totals.items() if k not in UNRESOLVED)
    observed = {k:v/known for k,v in totals.items() if k not in UNRESOLVED} if known else {}
    return {'total_known_measure':total, 'missing_measure_objects':sum(v is None for v in weights),
            'by_observed_family':totals, 'unresolved_measure':unknown,
            'fraction_with_usable_group':known/total if total else None,
            'fractions_among_usable_groups':observed,
            'identification_bounds':{k:[v/total,(v+unknown)/total]
                for k,v in totals.items() if k not in UNRESOLVED} if total else {},
            'unseen_family_upper_fraction':unknown/total if total else None,
            'conditional_observed_mixture_completion':{
                'assumption':'Unresolved labels follow the observed weighted spectral-family mixture, including C_or_X as its own category; selection bias is not corrected.',
                'family_measure':{k:totals[k]+unknown*f for k,f in observed.items()},
                'family_fractions':observed} if known else None}


def summary(rows: list[dict], matched: dict, density: float) -> dict:
    def distribution(source=None):
        states = [object_state([e for e in entries_for(r,matched)
                                if source is None or (e['source']=='mithneos') == (source=='mithneos')])
                  for r in rows]
        counts = dict(sorted(Counter(s['family'] for s in states).items()))
        usable = sum(v for k,v in counts.items() if k not in UNRESOLVED)
        return {'object_count': len(rows), 'counts': counts, 'usable_group_count': usable,
                'unambiguous_family_count': usable-counts.get('C_or_X',0),
                'fractions_of_pool': {k:v/len(rows) for k,v in counts.items()} if rows else {},
                'fractions_among_usable_groups': {k:v/usable for k,v in counts.items() if k not in UNRESOLVED} if usable else {},
                'objects_with_uncertain_class': sum(s['has_uncertain_class'] for s in states),
                'objects_with_excluded_low_snr': sum(s['has_excluded_low_snr'] for s in states),
                'count_inference': measure_distribution(states,[1.0]*len(rows)),
                'known_size_mass_kg': measure_distribution(states,[assumed_mass(r,density) for r in rows])}
    jpl = {r['spkid'] for r in rows if r['spec_B'] or r['spec_T']}
    survey = {r['spkid'] for r in rows if r['spkid'] in matched}
    return {'coverage': {'pool': len(rows), 'jpl': len(jpl), 'mithneos': len(survey),
                         'both': len(jpl & survey), 'either': len(jpl | survey),
                         'newly_observed_in_mithneos': len(survey-jpl),
                         'neither': len(rows)-len(jpl | survey)},
            'jpl_only_classifications': distribution('jpl'),
            'mithneos_only_classifications': distribution('mithneos'),
            'integrated_classifications': distribution(),
            'population_estimate': None, 'statistical_confidence_interval': None}


def stratified_completion(rows, matched, density):
    """Conditional allocation within albedo strata; never add overlapping pools."""
    strata={'dark':[], 'brighter':[], 'albedo_missing':[]}
    for row in rows:
        key='albedo_missing' if row['albedo'] is None else 'dark' if float(row['albedo'])<=0.1 else 'brighter'
        strata[key].append(row)
    summaries={k:summary(v,matched,density)['integrated_classifications'] for k,v in strata.items()}
    estimates={}
    for measure in ('count_inference','known_size_mass_kg'):
        families=defaultdict(float);unallocated=0.;total=0.
        for result in summaries.values():
            d=result[measure];total+=d['total_known_measure']
            completion=d['conditional_observed_mixture_completion']
            if completion is None:unallocated+=d['total_known_measure']
            else:
                for family,value in completion['family_measure'].items():families[family]+=value
        estimates[measure]={'total_known_measure':total,'unallocated_measure':unallocated,'family_measure':dict(families),'family_fractions':{k:v/total for k,v in families.items()} if total else {}}
    return {'strata':summaries,'conditional_estimates':estimates,'assumption':'Within each of dark (albedo <=0.10), brighter and missing-albedo strata, unresolved classes follow the observed weighted family mixture. Empty/unclassified strata retain unallocated measure. Observing selection within strata remains uncorrected; no statistical confidence or ore grade follows. This disjoint allocation preserves subset consistency; separate overlapping-pool completion scenarios must not be added.'}

def selected(row: dict, criteria: dict, dark: bool = True) -> bool:
    """Select the existing geometric screen; no mission accessibility inference."""
    try:
        values = [float(row[k]) for k in ('a','e','i','diameter')]
        if not all(math.isfinite(v) for v in values):
            return False
        a,e,i,diameter = values
        if not (criteria.get('a_min_au',0) <= a and 0 < a <= criteria['a_max_au'] and 0 <= e <= criteria['e_max']
                and 0 <= i <= criteria['i_max_deg'] and diameter > 0):
            return False
        if dark:
            albedo = float(row['albedo'])
            if not math.isfinite(albedo) or not 0 <= albedo <= criteria['albedo_max']:
                return False
        return row['number'] not in criteria['excluded_numbers']
    except (TypeError, ValueError):
        return False


def checked_bytes(path: Path, specification: dict) -> bytes:
    raw = path.read_bytes()
    if len(raw) != specification['bytes'] or digest(raw) != specification['sha256']:
        raise ValueError(f'Source hash or length mismatch: {path}')
    return raw


def product(config: dict | None = None, *, archive: Path | None = None, jpl_path: Path | None = None) -> dict:
    config = config or json.loads((HERE/'composition_inputs.json').read_text())
    archive = Path(archive) if archive else ROOT/config['archive_directory']
    rows = jpl_rows(checked_bytes(Path(jpl_path or ROOT/config['jpl'].get('cache_path',config['jpl']['path'])), config['jpl']))
    observations = parse_table(checked_bytes(archive/config['survey']['filename'], config['survey']))
    identities = {('number',o['number']) if o['number'] else ('designation',canonical_designation(o['designation'])) for o in observations}
    if (len(rows) != config['jpl']['expected_rows'] or len(observations) != config['survey']['expected_observations']
            or len(identities) != config['survey']['expected_identifier_groups']):
        raise ValueError('Unexpected source row or identity count')
    matched, unmatched = match_observations(rows, observations)
    wider = [r for r in rows if selected(r,config['wider_screen'])]
    if len(wider) != config['wider_screen']['expected_count']:
        raise ValueError('Wider pool no longer matches frozen reference')
    geometry = [r for r in rows if selected(r,config['wider_screen'],dark=False)]
    narrow = [r for r in rows if selected(r,config['narrow_screen'])]
    if len(narrow) != config['narrow_screen']['expected_count']:
        raise ValueError('Narrow pool no longer matches frozen reference')
    screened=protection.filter_source_sets({'wider':wider,'narrow':narrow,'geometry':geometry})['source_sets']
    wider,narrow,geometry=(screened[k] for k in ('wider','narrow','geometry'))
    density = config['density_assumption_kg_m3']
    survey_rows = [r for r in rows if r['spkid'] in matched]
    return {'schema': 'terluna.resources.asteroid-composition/1',
            'producer': {'path':'resources/source_access/composition.py','sha256':digest(Path(__file__).read_bytes())},
            'inputs_sha256': digest(json.dumps(config,sort_keys=True,separators=(',',':')).encode()),
            'protection_dependencies':{f: digest((HERE/f).read_bytes()) for f in ('protection.py','protection_inputs.json')},
            'sources': {'jpl':config['jpl'],'mithneos':config['survey']},
            'grouping': FAMILIES, 'limits': config['limits'], 'rights_scope':config['rights_scope'],
            'source_reading': config['source_reading'], 'aggregation':config['aggregation'],
            'survey_account': {'observations':len(observations), 'distinct_number_designation_identities':len(identities),
                'distinct_filename_prefixes':len({o['spectrum'].split('.')[0] for o in observations}),
                'published_object_count':config['survey']['published_object_count'],
                'matched_jpl_objects':len(matched), 'matched_observations':sum(map(len,matched.values())),
                'unmatched_observations':unmatched,
                'observation_class_counts':dict(sorted(Counter(o['taxonomy'] for o in observations).items())),
                'matched_observed_sample': summary(survey_rows,matched,density)['mithneos_only_classifications']},
            'pools': {'all_frozen_nea':summary(rows,matched,density),
                      'wider_geometry_all_albedos':summary(geometry,matched,density),
                      'wider_dark_screen':summary(wider,matched,density),
                      'narrow_dark_screen':summary(narrow,matched,density)},
            'wider_stratified_completion':stratified_completion(geometry,matched,density),
            'wider_screen':config['wider_screen'], 'narrow_screen':config['narrow_screen'],
            'density_assumption_kg_m3':density,
            'mass_scope':'Uniform-density spherical mass for objects with positive catalogue diameters; the mass of objects without measured size is unbounded here. Taxonomy assigns no density.',
            'narrow_object_ids':[r['spkid'] for r in narrow],
            'wider_objects': [{'spkid':r['spkid'],'full_name':r['full_name'],
                'diameter_km':float(r['diameter']),'assumed_mass_kg':assumed_mass(r,density),'geometric_albedo':float(r['albedo']),
                **object_state(entries_for(r,matched))} for r in wider],
            'grade_inference':None, 'qualified_reserve':None}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--download',action='store_true',help='Restore the exact publisher table to ignored inputs')
    p.add_argument('--output',type=Path,default=HERE/'results/composition.json.gz')
    args=p.parse_args()
    config=json.loads((HERE/'composition_inputs.json').read_text())
    if args.download:
        location=ROOT/config['archive_directory']/config['survey']['filename']
        if not location.exists():
            raw=urllib.request.urlopen(config['survey']['url'],timeout=60).read()
            if digest(raw)!=config['survey']['sha256'] or len(raw)!=config['survey']['bytes']:
                raise ValueError('Retrieved source differs from the pinned publisher table')
            location.parent.mkdir(parents=True,exist_ok=True)
            location.write_bytes(raw)
    result=product(config)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_bytes(gzip.compress(json.dumps(result,sort_keys=True,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    print(json.dumps({'output':str(args.output),'sha256':digest(args.output.read_bytes()),
                      'survey':result['survey_account']['matched_jpl_objects'],
                      'wider_coverage':result['pools']['wider_dark_screen']['coverage']}))


if __name__ == '__main__':
    main()
