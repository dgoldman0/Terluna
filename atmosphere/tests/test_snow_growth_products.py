"""Final metadata qualification preserves frozen runs and rejects stale inputs."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from atmosphere.snow import growth_products as gp


def raw_row(identity, *, reference=False):
    row = dict(
        id=identity, specification={'maximum_diameter_m': 10.},
        producer={'files': {}, 'inputs': {}},
        summary={'snowfall_mm_h': 8., 'arriving_mass_share_above': {'1.0': 1., '10.0': 0.}},
        domain_adequacy={'passed': not identity.startswith('moon')},
        arriving_numbers_m2=[2., 3.], gathering_seconds=3600.,
        evidence='Original evidence', interpretation='Original interpretation',
    )
    masses = [1., 2.]
    diameters = [1., np.nextafter(10., 0.)]
    if reference:
        row.update(arriving_mass_pivots_kg=masses, arriving_diameter_pivots_m=diameters,
                   reading_rule='Original Qext=2 homogeneous direct-beam rule',
                   numerical_status='domain_passed_refinement_pending')
        row['summary'].update(extinction_m_inv=.1, beam_5_percent_distance_m=29.9573227355)
    else:
        row.update(masses_kg=masses, diameters_m=diameters)
    return row


@pytest.fixture
def assembly(monkeypatch, tmp_path):
    monkeypatch.setattr(gp, 'ROOT', tmp_path)
    for relative in gp.FILES:
        p=tmp_path/relative; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(relative)
    limited=tmp_path/'limited'; reference=tmp_path/'reference'
    limited.mkdir(); reference.mkdir()
    monkeypatch.setattr(gp.frozen, 'CACHE', limited)
    monkeypatch.setattr(gp.repair, 'CACHE', reference)
    monkeypatch.setattr(gp.frozen, 'cases', lambda: [{'id': 'earth_limited_m144'}])
    monkeypatch.setattr(gp, 'REFERENCE_IDS', ('earth_moderate_m144',))
    destination=tmp_path/'public.json'; destination.write_text('previous public product\n')
    monkeypatch.setattr(gp, 'OUT', destination)
    monkeypatch.setattr(gp.frozen, 'OUT', destination)
    case=raw_row('earth_limited_m144')
    ref=raw_row('earth_moderate_m144', reference=True)
    (limited/(case['id']+'.json')).write_text(json.dumps(case))
    (reference/(ref['id']+'.json')).write_text(json.dumps(ref))
    reduced={key:ref[key] for key in ('id','specification','producer','summary','domain_adequacy')}
    product=dict(schema='terluna.atmosphere.snow-growth-repair/1', cases=[case],
                 references=[reduced], reading_rule='Original product reading rule.', unperformed=[])
    calls=[]
    def admitted():
        calls.append(gp.frozen.OUT)
        gp.frozen.OUT.write_text(json.dumps(product))
        return copy.deepcopy(product)
    monkeypatch.setattr(gp.frozen, 'assemble', admitted)
    return dict(product=product, case=case, reference=ref, destination=destination,
                limited=limited, reference_cache=reference, calls=calls)


def test_final_export_recovers_qualifications_preserves_summaries_and_binds_caches(assembly):
    result=gp.export()
    assert json.loads(assembly['destination'].read_text())==result
    assert gp.frozen.OUT==assembly['destination']
    assert assembly['calls'][0]!=assembly['destination']
    for group in ('cases','references'):
        row=result[group][0]
        assert row['summary']==assembly['product'][group][0]['summary']
        assert row['producer']==assembly['product'][group][0]['producer']
        assert not row['optical_prediction_eligible']
        assert 'Generic Earth' in row['optical_prediction_reason']
        diagnostic=row['nominal_endpoint_diagnostic']
        assert diagnostic['endpoint_matches_within_roundoff']
        assert diagnostic['top_resolved_pivot_arriving_mass_share']==.75
        assert diagnostic['arriving_mass_share_at_or_above_nominal_endpoint']==.75
        assert row['summary']['arriving_mass_share_above']['10.0']==0.
    ref=result['references'][0]
    for key in ('reading_rule','evidence','numerical_status'):
        assert ref[key]==assembly['reference'][key]
    assert result['assembly']['producer']==gp.identity()
    assert result['assembly']['identity_checked_before_and_after']


def test_build_does_not_write_public_product(assembly):
    gp.build_product()
    assert assembly['destination'].read_text()=='previous public product\n'


def test_frozen_admission_failure_is_propagated_without_overwriting(assembly, monkeypatch):
    def stale(): raise RuntimeError('Stale checkpoint reference')
    monkeypatch.setattr(gp.frozen,'assemble',stale)
    with pytest.raises(RuntimeError,match='Stale checkpoint'):
        gp.export()
    assert gp.frozen.OUT==assembly['destination']
    assert assembly['destination'].read_text()=='previous public product\n'


def test_cache_change_between_admission_and_qualification_rejects(assembly, monkeypatch):
    admitted=gp.frozen.assemble
    def changed():
        result=admitted()
        with (assembly['limited']/'earth_limited_m144.json').open('a') as handle:
            handle.write(' ')
        return result
    monkeypatch.setattr(gp.frozen,'assemble',changed)
    with pytest.raises(RuntimeError,match='Raw cache changed'):
        gp.export()
    assert assembly['destination'].read_text()=='previous public product\n'


def test_postprocessor_source_change_rejects_without_overwriting(assembly, monkeypatch):
    qualify=gp._qualify
    def changed(*args,**kwargs):
        qualify(*args,**kwargs)
        (gp.ROOT/gp.FILES[0]).write_text('changed postprocessor')
    monkeypatch.setattr(gp,'_qualify',changed)
    with pytest.raises(RuntimeError,match='changed during final export'):
        gp.export()
    assert assembly['destination'].read_text()=='previous public product\n'


def test_nominal_match_does_not_snap_a_physically_different_endpoint():
    row=raw_row('case'); row['diameters_m'][-1]=9.99
    diagnostic=gp.endpoint_diagnostic(row)
    assert not diagnostic['endpoint_matches_within_roundoff']
    assert diagnostic['arriving_mass_share_at_or_above_nominal_endpoint']==0.
    assert diagnostic['top_resolved_pivot_arriving_mass_share']==.75


def test_native_flux_disagreement_rejects():
    row=raw_row('case'); row['arriving_numbers_m2'][0]+=1.
    with pytest.raises(ValueError,match='disagrees with saved flux'):
        gp.endpoint_diagnostic(row)


def test_missing_reference_qualification_rejects(assembly):
    p=assembly['reference_cache']/'earth_moderate_m144.json'
    raw=json.loads(p.read_text()); del raw['reading_rule']; p.write_text(json.dumps(raw))
    with pytest.raises(ValueError,match='Missing reference qualification'):
        gp.export()
    assert assembly['destination'].read_text()=='previous public product\n'
