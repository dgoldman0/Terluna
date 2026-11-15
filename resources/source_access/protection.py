"""Evaluate recorded project protection conditions without authorizing extraction.

Scientific samples and records are separate from the places and functions they
represent. Missing evidence remains unknown; money and catalogue mass cannot
substitute for a completed preservation condition.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATUSES = {'passed', 'failed', 'unknown', 'hold'}


def load_inputs(path=None):
    path = Path(path or HERE/'protection_inputs.json')
    config = json.loads(path.read_text())
    if config.get('schema') != 'terluna.resources.source-protection-inputs/1':
        raise ValueError('Unexpected source protection schema')
    for name, source in config['sources'].items():
        actual = hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()
        if actual != source['sha256']:
            raise ValueError('Changed protection protocol source: '+name)
    identifiers = [gate['id'] for gate in config['gates']]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError('Duplicate protection gate')
    return config


def _identifiers(subject):
    return {str(subject[key]).strip() for key in ('source_id', 'spkid', 'id')
            if subject.get(key) is not None and str(subject[key]).strip()}


def _identifier(subject):
    identifiers = _identifiers(subject)
    if len(identifiers) > 1:
        raise ValueError('Conflicting source IDs require explicit reconciliation')
    return next(iter(identifiers), None)


def exclusion_reasons(subject, config=None):
    """Retained illustrative science holds and adopted protected categories."""
    config = config or load_inputs()
    name = str(subject.get('name', subject.get('full_name', '')))
    identifiers = _identifiers(subject)
    reasons = []
    for held in config['named_science_holds']:
        matches_name = re.search(r'\b'+re.escape(held['name'])+r'\b', name, flags=re.IGNORECASE)
        if held['source_id'] in identifiers or matches_name:
            reasons.append(dict(kind='named_science_hold',source_id=held['source_id'],name=held['name'],
                                status=config['named_science_hold_status']))
    category = subject.get('protected_category')
    if category in config['protected_categories']:
        reasons.append(dict(kind='protected_category',category=category,
                            basis=deepcopy(config['protected_category_basis'])))
    if subject.get('protected') is True:
        reasons.append(dict(kind='recorded_protection_hold',record=subject.get('protection_record'),
                            status='Protection asserted in the source record; this evaluator cannot remove it.'))
    return reasons


def _attestation(value):
    if value is None:
        return dict(status='unknown',evidence=[])
    if not isinstance(value, dict) or value.get('status') not in STATUSES:
        raise ValueError('Attestation needs an explicit passed/failed/unknown/hold status')
    evidence = value.get('evidence', [])
    if not isinstance(evidence, list) or any(not isinstance(item, str) or not item.strip() for item in evidence):
        raise ValueError('Evidence must be a list of nonempty record references')
    if value['status'] in ('passed','failed','hold') and not evidence:
        raise ValueError('An asserted condition needs a supporting record reference')
    return dict(status=value['status'],evidence=deepcopy(evidence),
                note=value.get('note'),supplied_attestation=True)


def evaluate_release(subject, attestations=None, config=None):
    """Return completeness/status of supplied records, never a mining permit.

    subject.region is lunar or near_earth. Attestations map gate IDs to records
    containing status and evidence references. Their truth is not verified here.
    """
    config = config or load_inputs()
    region = subject.get('region')
    if region not in ('lunar','near_earth'):
        raise ValueError('Source region must be lunar or near_earth')
    attestations = {} if attestations is None else attestations
    if not isinstance(attestations, dict):
        raise ValueError('Attestations must be keyed by protection gate')
    known = {gate['id'] for gate in config['gates']}
    unexpected = set(attestations)-known
    if unexpected:
        raise ValueError('Unknown protection attestations: '+', '.join(sorted(unexpected)))
    conditions=[]
    for gate in config['gates']:
        applicable = gate['scope'] in ('both',region)
        status = _attestation(attestations.get(gate['id'])) if applicable else dict(status='not_applicable',evidence=[])
        conditions.append({**deepcopy(gate),**status, 'applicable':applicable})
    exclusions = exclusion_reasons(subject,config)
    failed = [g['id'] for g in conditions if g['status']=='failed']
    held = [g['id'] for g in conditions if g['status']=='hold']
    unknown = [g['id'] for g in conditions if g['status']=='unknown']
    blocked = bool(exclusions or failed or held)
    release = 'blocked' if blocked else ('unresolved' if unknown else 'conditions_satisfied')
    return dict(schema='terluna.resources.source-protection-status/1',
                subject=deepcopy(subject),release_status=release,extraction_authorized=False,
                conditions_satisfied=release=='conditions_satisfied',
                failed_conditions=failed,held_conditions=held,unknown_conditions=unknown,
                source_holds=exclusions,conditions=conditions,
                evidence='Evaluation of supplied project-condition attestations. Record references are not independently verified; passing them grants no permission, reserve status or extraction quantity.',
                authority_rule='Decided requirements and conservative conditions of proposed origin retain their individual authority_status. No approval is inferred from a missing objection.',
                protection_system_scope=deepcopy(config['protection_system_scope']))


def filter_source_sets(source_sets, config=None):
    """Apply identical holds independently to every candidate source pool.

    Filtering is not qualification. Every retained source still needs the
    release evaluation and a comprehensive scientific/heritage review.
    """
    config = config or load_inputs()
    clean={}; excluded={}
    for name, subjects in source_sets.items():
        clean[name]=[];excluded[name]=[]
        for subject in subjects:
            reasons=exclusion_reasons(subject,config)
            if reasons:excluded[name].append(dict(subject=deepcopy(subject),reasons=reasons))
            else:clean[name].append(deepcopy(subject))
    return dict(source_sets=clean,excluded=excluded,
                rule='The retained illustrative exclusions apply to every pool before set operations; other candidates remain unapproved.')


def source_difference(broad, narrow, config=None):
    """Return screened broad-minus-narrow rows without reintroducing a hold."""
    screened=filter_source_sets({'broad':broad,'narrow':narrow},config)
    broad=screened['source_sets']['broad'];narrow=screened['source_sets']['narrow']
    if any(_identifier(row) is None for row in broad+narrow):
        raise ValueError('A stable source ID is required for set subtraction')
    remove={_identifier(row) for row in narrow};seen=set();result=[]
    for row in broad:
        key=_identifier(row)
        if key not in seen and key not in remove:
            seen.add(key);result.append(row)
    return dict(sources=result,excluded=screened['excluded'],
                rule=screened['rule'],extraction_authorized=False)
