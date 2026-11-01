"""Admission checks for diagnostic products reject corrupted physical accounts."""
import copy
import shutil
import numpy as np
import pytest
from atmosphere.snow import column_diagnosis as model, column_diagnosis_runs as runs

pytestmark=pytest.mark.skipif(shutil.which('g++') is None,reason='C++17 compiler required.')

@pytest.fixture(scope='module')
def raw():
    return model.column(model.sf.moon_column()[0],top_c=-1.,layers=3,area_m2=.0001,fill=.2,gather=.2,seed=3144)


def test_trace_admits_consistent_observations(raw):
    assert runs.validate_trace(raw)
    summary,profiles=runs.metrics(raw)
    assert sum(summary[k] for k in ('gather_deposition_mm_h','gather_empty_unused_mm_h','gather_occupied_unused_mm_h'))==pytest.approx(3.,rel=2e-12)
    assert summary['whole_deposition_kg_m2']==pytest.approx(raw['mass_ledger_kg_m2']['deposited'],rel=2e-12)
    assert len(profiles['empty_fraction'])==3


@pytest.mark.parametrize('field',[0,1,5,6,7,12,13,14])
def test_trace_rejects_altered_accounting(raw,field):
    changed=copy.deepcopy(raw)
    changed['process_diagnostics']['blocks'][1][0][field]+=1.
    with pytest.raises(ValueError):runs.validate_trace(changed)


def test_statistics_use_independent_realizations(raw):
    rows=[dict(status='complete',result=copy.deepcopy(raw)) for _ in range(4)]
    result=runs.summarize(rows)
    assert result['realizations']==4
    assert result['standard_error_process_summary']['bottom_empty_fraction']==0
    np.testing.assert_array_equal(result['standard_error_mass_cdf'],0)


def test_full_baselines_are_bound_to_the_previous_cases():
    specs=runs.cases()
    assert len(specs)==40 and len({s['id'] for s in specs})==40
    for spec in specs:
        ref=runs.original_case(spec)
        if ref:
            assert all(ref['settings'][k]==v for k,v in spec['parameters'].items())
