"""Attach the numerical audit to raw snow products without changing their arrays.

Run after population, growth_depth and numerical_audit, before optics. The
assessment is a separate, hash-identified producer stage; the original numerical
producer identities and raw-product hash are retained. Re-running a raw producer
removes its assessment, so repeat this stage before downstream use.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
SELF = 'atmosphere/snow/assessment.py'
AUDIT = 'atmosphere/snow/results/numerical_audit.json'
TARGETS = {
    'atmosphere/snow/results/population.json': 'terluna.atmosphere.snow-population/1',
    'atmosphere/snow/results/growth_depth.json': 'terluna.atmosphere.snow-growth-depth/1',
}
SCHEMA = 'terluna.atmosphere.snow-assessment/1'
REPLAY_CHECKS = ('summary_exactly_equal', 'mass_ledger_exactly_equal',
                 'maximum_arriving_diameter_exactly_equal')
AUDITED_MODEL = (
    'atmosphere/snow/aggregation.py', 'atmosphere/snow/habits.py',
    'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py',
    'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py',
)
WARNING = (
    'NUMERICAL ADEQUACY UNRESOLVED: the inherited fractional-weight collision '
    'algorithm mixes real weights with integer splitting and an absolute fragment '
    'threshold. Its heavy-case replay exhibits large all-or-none collision batches '
    'already at decimetre growth and a physically impossible extreme tail. '
    'These are exploratory numerical populations, not quantitative snowfall '
    'predictions. Removing extreme bins does not validate the smaller population. '
    'The audit does not quantify the bias in every case or rule out physical '
    'decimetre aggregates.'
)
RULE = (
    'Read the attached numerical_assessment and its hash-identified audit before '
    'using sizes or deriving storm optics. The assessment changes evidence '
    'metadata only; all original numerical results remain intact. '
)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical_digest(value):
    return digest(json.dumps(value, sort_keys=True, separators=(',', ':'),
                             allow_nan=False).encode())


def numerical_digest(product):
    return canonical_digest(product['runs'])


def verify_current(product):
    producer = product['producer']
    for group in ('files', 'inputs'):
        for path, expected in producer[group].items():
            if digest((ROOT/path).read_bytes()) != expected:
                raise ValueError(f'Stale {group} identity: {path}')
    if constants_changed(producer['constants']):
        raise ValueError('Shared constants changed')


def assess(raw, audit, audit_identity, assessor_identity):
    """Return assessed bytes, or identical bytes for an already current assessment."""
    product = json.loads(raw)
    if product['schema'] not in TARGETS.values():
        raise ValueError('Unsupported population schema')
    if audit['schema'] != 'terluna.atmosphere.snow-numerical-audit/1':
        raise ValueError('Unsupported numerical audit schema')
    heavy = audit.get('heavy_replay') or {}
    replay = heavy.get('unchanged_replay_checks', {})
    if set(replay) != set(REPLAY_CHECKS) or not all(replay[key] is True for key in REPLAY_CHECKS):
        raise ValueError('Assessment requires the successful unchanged heavy replay')
    for path in AUDITED_MODEL:
        expected = audit['producer']['files'].get(path)
        if expected is None or product['producer']['files'].get(path) != expected:
            raise ValueError(f'Population differs from audited model: {path}')
    numerical_hash = numerical_digest(product)
    existing = product.get('numerical_assessment')
    if existing is not None:
        if (existing.get('schema') != SCHEMA or existing.get('audit') != audit_identity
                or existing.get('assessor') != assessor_identity
                or existing.get('numerical_payload_sha256') != numerical_hash):
            raise ValueError('Stale assessment; reproduce the raw product before reassessing')
        if (product['evidence'] != WARNING+' '+existing['original_evidence']
                or product['reading_rule'] != RULE+existing['original_reading_rule']):
            raise ValueError('Assessment text contradicts its recorded original')
        restored = dict(product)
        del restored['numerical_assessment']
        restored['evidence'] = existing['original_evidence']
        restored['reading_rule'] = existing['original_reading_rule']
        if canonical_digest(restored) != existing.get('original_canonical_sha256'):
            raise ValueError('Original product metadata changed after assessment')
        return raw
    original_evidence, original_rule = product['evidence'], product['reading_rule']
    product['numerical_assessment'] = dict(
        schema=SCHEMA, status='exploratory_numerically_unresolved',
        audit=audit_identity, assessor=assessor_identity,
        original_product_sha256=digest(raw), original_canonical_sha256=canonical_digest(product),
        numerical_payload_sha256=numerical_hash,
        original_evidence=original_evidence, original_reading_rule=original_rule,
        scope='Evidence metadata only. Original producer/input identities and every run are preserved.',
    )
    product['evidence'] = WARNING+' '+original_evidence
    product['reading_rule'] = RULE+original_rule
    result = (json.dumps(product, indent=1, allow_nan=False)+'\n').encode()
    if numerical_digest(json.loads(result)) != numerical_hash:
        raise ValueError('Assessment changed numerical payload')
    return result


def main():
    audit_raw = (ROOT/AUDIT).read_bytes()
    audit = json.loads(audit_raw)
    verify_current(audit)
    audit_identity = dict(path=AUDIT, sha256=digest(audit_raw))
    assessor_identity = dict(path=SELF, sha256=digest((ROOT/SELF).read_bytes()))
    # Prepare and verify both products before writing either one.
    prepared = []
    for relative, schema in TARGETS.items():
        path = ROOT/relative
        raw = path.read_bytes()
        product = json.loads(raw)
        if product['schema'] != schema:
            raise ValueError(f'Wrong schema: {relative}')
        verify_current(product)
        prepared.append((path, assess(raw, audit, audit_identity, assessor_identity)))
    for path, raw in prepared:
        if path.read_bytes() != raw:
            temporary = path.with_suffix('.assessed.tmp')
            temporary.write_bytes(raw)
            temporary.replace(path)
        print(f'assessed {path.relative_to(ROOT)}', flush=True)


if __name__ == '__main__':
    main()
