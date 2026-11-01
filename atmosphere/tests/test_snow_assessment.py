"""Evidence qualifiers preserve numerical results and reject stale model identity."""
import copy
import json

import pytest

from atmosphere.snow import assessment as a
from shared.baselines import recorded_file


def fixture():
    files = {path: f'hash-{i}' for i,path in enumerate(a.AUDITED_MODEL)}
    product = dict(schema=next(iter(a.TARGETS.values())), producer=dict(files=files, inputs={'column':'input'}),
                   evidence='Original conditional calculation.', reading_rule='Original units.',
                   runs=[dict(id='case', numbers=[0., .2, 1e-25], mass=[2.3, 4.5], edges=[.01, 10., None])])
    audit = dict(schema='terluna.atmosphere.snow-numerical-audit/1', producer=dict(files=files.copy()),
                 heavy_replay=dict(unchanged_replay_checks={key:True for key in a.REPLAY_CHECKS}))
    return (json.dumps(product, indent=1)+'\n').encode(), audit


def test_assessment_preserves_payload_identity_original_text_and_provenance():
    raw,audit=fixture(); before=json.loads(raw)
    result=a.assess(raw,audit,{'path':'audit','sha256':'audit-hash'},{'path':'assessor','sha256':'code-hash'})
    after=json.loads(result); review=after['numerical_assessment']
    assert after['runs']==before['runs']
    assert after['producer']==before['producer']
    assert review['original_product_sha256']==a.digest(raw)
    assert review['original_canonical_sha256']==a.canonical_digest(before)
    assert review['numerical_payload_sha256']==a.numerical_digest(before)
    assert review['original_evidence']==before['evidence']
    assert review['original_reading_rule']==before['reading_rule']
    assert after['evidence'].startswith(a.WARNING)
    assert a.assess(result,audit,review['audit'],review['assessor'])==result


def test_assessment_rejects_different_model_incomplete_replay_and_changed_numerical_payload():
    raw,audit=fixture()
    wrong=copy.deepcopy(audit);wrong['producer']['files'][a.AUDITED_MODEL[0]]='changed'
    with pytest.raises(ValueError,match='differs from audited'):
        a.assess(raw,wrong,{}, {})
    wrong=copy.deepcopy(audit);wrong['heavy_replay']['unchanged_replay_checks']['mass_ledger_exactly_equal']=False
    with pytest.raises(ValueError,match='successful unchanged'):
        a.assess(raw,wrong,{}, {})
    result=a.assess(raw,audit,{}, {});changed=json.loads(result);changed['runs'][0]['mass'][0]+=1
    with pytest.raises(ValueError,match='Stale assessment'):
        a.assess(json.dumps(changed).encode(),audit,{}, {})
    with pytest.raises(ValueError,match='Stale assessment'):
        a.assess(result,audit,{'sha256':'new-audit'}, {})


def test_reassessment_rejects_changed_original_metadata_even_when_runs_are_unchanged():
    raw,audit=fixture()
    result=a.assess(raw,audit,{}, {})
    changed=json.loads(result)
    changed['producer']['inputs'].clear()
    with pytest.raises(ValueError,match='Original product metadata changed'):
        a.assess(json.dumps(changed).encode(),audit,{}, {})


def current(product):
    """assessment.verify_current, reading a baselined file's kept copy where the product recorded it
    (shared/baselines.py): the audited products were built with the rain solver of 24 October."""
    from shared.provenance import constants_changed
    for group in ('files','inputs'):
        for path,expected in product['producer'][group].items():
            assert a.digest(recorded_file(path,expected).read_bytes())==expected,path
    assert not constants_changed(product['producer']['constants'])


def test_published_assessments_bind_current_audit_and_unchanged_payloads():
    for relative in a.TARGETS:
        path=a.ROOT/relative
        if not path.exists() or 'numerical_assessment' not in json.loads(path.read_bytes()):
            pytest.skip('Both final population assessments have not been produced')
    audit_raw=(a.ROOT/a.AUDIT).read_bytes();audit=json.loads(audit_raw)
    current(audit)
    for relative in a.TARGETS:
        raw=(a.ROOT/relative).read_bytes();product=json.loads(raw)
        current(product)
        review=product['numerical_assessment']
        assert review['audit']==dict(path=a.AUDIT,sha256=a.digest(audit_raw))
        assert review['assessor']==dict(path=a.SELF,sha256=a.digest((a.ROOT/a.SELF).read_bytes()))
        assert review['numerical_payload_sha256']==a.numerical_digest(product)
        assert a.assess(raw,audit,review['audit'],review['assessor'])==raw
